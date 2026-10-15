"""Altair charts for the results view. Data prep is split out so it can be tested without a UI."""

from __future__ import annotations

from collections.abc import Sequence

import altair as alt
import pandas as pd

from src import config

ANCHOR_LABEL = "Anchor (P2 / C-terminus)"
TCR_LABEL = "T-cell-facing"


# --- Data prep ------------------------------------------------------------------
def funnel_data(stages: Sequence[tuple[str, int]]) -> pd.DataFrame:
    """One row per funnel stage, in pipeline order, with the unit each count is in."""
    df = pd.DataFrame(stages, columns=["stage", "count"])
    df["order"] = range(len(df))
    df["unit"] = [config.FUNNEL_STAGE_UNITS.get(s, "") for s in df["stage"]]
    return df


def scatter_data(table: pd.DataFrame, max_points: int = config.VIZ_MAX_POINTS) -> pd.DataFrame:
    """Mutant binders (peptide x allele rows) with both affinities, capped to the top rows.

    The table arrives sorted by SORT_ORDER, so the cap keeps the best-ranked rows.
    """
    binders = table[(table["binder_rank"] < 2) & table["wt_affinity_nM"].notna()]
    df = binders.head(max_points).copy()
    df["position_class"] = df["anchor_position"].map({True: ANCHOR_LABEL, False: TCR_LABEL})
    df["variant"] = df["gene"] + " " + df["hgvsp"].fillna("")
    cols = [
        "variant", "allele", "mut_peptide", "wt_peptide", "offset", "position_class",
        "mut_affinity_nM", "wt_affinity_nM", "median_tumor_vaf", "n_patients",
    ]
    return df[cols].reset_index(drop=True)


def recurrence_data(table: pd.DataFrame, top_n: int = config.VIZ_TOP_RECURRENT) -> pd.DataFrame:
    """Variants carried by >= 2 patients with at least one predicted binder, most shared first."""
    binders = table[(table["binder_rank"] < 2) & (table["n_patients"] >= 2)]
    if binders.empty:
        return pd.DataFrame(
            columns=["variant", "n_patients", "best_mut_affinity_nM", "n_alleles", "n_peptides"]
        )
    df = (
        binders.assign(variant=binders["gene"] + " " + binders["hgvsp"].fillna(""))
        .groupby("variant", as_index=False)
        .agg(
            n_patients=("n_patients", "max"),
            best_mut_affinity_nM=("mut_affinity_nM", "min"),
            n_alleles=("allele", "nunique"),
            n_peptides=("mut_peptide", "nunique"),
        )
        .sort_values(["n_patients", "best_mut_affinity_nM"], ascending=[False, True])
    )
    return df.head(top_n).reset_index(drop=True)


# --- Charts -----------------------------------------------------------------------
def funnel_chart(stages: Sequence[tuple[str, int]]) -> alt.VConcatChart:
    """Two panels with their own x scales, because the stages are counted in different units."""
    df = funnel_data(stages)
    is_row = df["unit"] == "MAF rows"
    return alt.vconcat(
        _funnel_panel(df[is_row], "Variant rows (MAF rows)"),
        _funnel_panel(df[~is_row], "Peptides and predictions"),
    ).resolve_scale(x="independent")


def _funnel_panel(df: pd.DataFrame, title: str) -> alt.LayerChart:
    base = alt.Chart(df, title=alt.TitleParams(title, anchor="start", fontSize=13)).encode(
        y=alt.Y("stage:N", sort=alt.SortField("order"), title=None),
        x=alt.X("count:Q", title="Count", axis=alt.Axis(format=",d")),
        tooltip=[
            alt.Tooltip("stage:N", title="Stage"),
            alt.Tooltip("count:Q", title="Count", format=","),
            alt.Tooltip("unit:N", title="Counted as"),
        ],
    )
    bars = base.mark_bar(color=config.VIZ_COLOR_PRIMARY, cornerRadiusEnd=4, height={"band": 0.7})
    labels = base.mark_text(align="left", dx=4).encode(text=alt.Text("count:Q", format=","))
    return (bars + labels).properties(height=34 * len(df), width="container")


def affinity_scatter(table: pd.DataFrame) -> alt.LayerChart:
    df = scatter_data(table)
    lo = max(min(df["mut_affinity_nM"].min(), df["wt_affinity_nM"].min()) / 1.5, 0.1)
    x_scale = alt.Scale(type="log", domain=[lo, df["wt_affinity_nM"].max() * 1.5])
    # Every point is a mutant binder, so the y range stops just above the weak cut-off.
    y_scale = alt.Scale(type="log", domain=[lo, config.WEAK_BINDER_NM * 1.5])
    points = (
        alt.Chart(df)
        .mark_point(filled=True, size=64, opacity=0.55, clip=True)
        .encode(
            x=alt.X("wt_affinity_nM:Q", scale=x_scale, title="Wild-type affinity (nM, log)"),
            y=alt.Y("mut_affinity_nM:Q", scale=y_scale, title="Mutant affinity (nM, log)"),
            color=alt.Color(
                "position_class:N",
                title="Mutation position",
                scale=alt.Scale(
                    domain=[TCR_LABEL, ANCHOR_LABEL],
                    range=[config.VIZ_COLOR_PRIMARY, config.VIZ_COLOR_SECONDARY],
                ),
                legend=alt.Legend(orient="top"),
            ),
            shape=alt.Shape(
                "position_class:N",
                scale=alt.Scale(domain=[TCR_LABEL, ANCHOR_LABEL], range=["circle", "triangle-up"]),
                legend=None,
            ),
            tooltip=[
                alt.Tooltip("variant:N", title="Variant"),
                alt.Tooltip("allele:N", title="Allele"),
                alt.Tooltip("mut_peptide:N", title="Mutant"),
                alt.Tooltip("wt_peptide:N", title="Wild type"),
                alt.Tooltip("offset:Q", title="Mutation at position"),
                alt.Tooltip("mut_affinity_nM:Q", title="Mutant nM", format=".1f"),
                alt.Tooltip("wt_affinity_nM:Q", title="Wild-type nM", format=".1f"),
                alt.Tooltip("median_tumor_vaf:Q", title="Median tumour VAF", format=".2f"),
                alt.Tooltip("n_patients:Q", title="Patients"),
            ],
        )
    )
    top = config.WEAK_BINDER_NM * 1.5
    diagonal = (
        alt.Chart(pd.DataFrame({"x": [lo, top], "y": [lo, top]}))
        .mark_line(color="gray", strokeDash=[4, 4], strokeWidth=1, clip=True)
        .encode(x=alt.X("x:Q", scale=x_scale), y=alt.Y("y:Q", scale=y_scale))
    )
    thresholds = pd.DataFrame(
        {
            "nM": [config.STRONG_BINDER_NM, config.WEAK_BINDER_NM],
            "label": [
                f"strong < {config.STRONG_BINDER_NM:g} nM",
                f"weak < {config.WEAK_BINDER_NM:g} nM",
            ],
        }
    )
    rules = alt.Chart(thresholds).mark_rule(color="gray", strokeWidth=1, opacity=0.6).encode(
        y=alt.Y("nM:Q", scale=y_scale)
    )
    rule_labels = (
        alt.Chart(thresholds)
        .mark_text(align="left", baseline="bottom", dx=4, dy=-2, color="gray")
        .encode(x=alt.value(0), y=alt.Y("nM:Q", scale=y_scale), text="label:N")
    )
    return (
        (diagonal + rules + rule_labels + points)
        .properties(height=480)
        .configure_axis(gridOpacity=0.35)
        .interactive()
    )


def recurrence_chart(table: pd.DataFrame) -> alt.Chart:
    df = recurrence_data(table)
    base = alt.Chart(df).encode(
        y=alt.Y("variant:N", sort=None, title=None),
        x=alt.X(
            "n_patients:Q",
            title="Patients carrying the variant (in the loaded files)",
            axis=alt.Axis(values=list(range(0, int(df["n_patients"].max()) + 1)), format="d"),
            scale=alt.Scale(domain=[0, int(df["n_patients"].max())]),
        ),
        tooltip=[
            alt.Tooltip("variant:N", title="Variant"),
            alt.Tooltip("n_patients:Q", title="Patients"),
            alt.Tooltip("best_mut_affinity_nM:Q", title="Best mutant nM", format=".1f"),
            alt.Tooltip("n_alleles:Q", title="Alleles with a binder"),
            alt.Tooltip("n_peptides:Q", title="Binder peptides"),
        ],
    )
    bars = base.mark_bar(color=config.VIZ_COLOR_PRIMARY, cornerRadiusEnd=4, height={"band": 0.7})
    labels = base.mark_text(align="left", dx=4).encode(text="n_patients:Q")
    return (bars + labels).properties(height=max(22 * len(df), 120))
