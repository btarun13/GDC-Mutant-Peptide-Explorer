"""Join peptides, variant evidence and predictions; add flags; sort by a config-driven tuple."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd

from src import config

VARIANT_KEY = ["ensp", "position", "ref", "alt"]
PRED_VALUES = ["affinity_nM", "percentile", "presentation_score", "presentation_percentile"]

OUTPUT_COLUMNS = [
    "gene",
    "hgvsp",
    "transcript_id",
    "ensp",
    "allele",
    "length",
    "mut_peptide",
    "wt_peptide",
    "offset",
    "anchor_position",
    "mut_affinity_nM",
    "wt_affinity_nM",
    "delta_affinity_nM",
    "mut_percentile",
    "wt_percentile",
    "delta_percentile",
    "mut_presentation_score",
    "wt_presentation_score",
    "binder_strong",
    "binder_weak",
    "better_than_wt",
    "median_tumor_vaf",
    "low_vaf",
    "max_t_alt_count",
    "n_patients",
    "n_samples",
    "binder_rank",
    "gdc_file_ids",
]


def build_table(
    pairs: pd.DataFrame,
    variants: pd.DataFrame,
    predictions: pd.DataFrame,
    sort_order: Sequence[tuple[str, bool]] = config.SORT_ORDER,
) -> pd.DataFrame:
    if pairs.empty or predictions.empty:
        return pd.DataFrame(columns=OUTPUT_COLUMNS)
    preds = predictions.drop_duplicates(["peptide", "allele"])
    mut = preds.rename(columns={c: f"mut_{c}" for c in PRED_VALUES})
    wt = preds.rename(columns={c: f"wt_{c}" for c in PRED_VALUES})
    df = pairs.merge(mut, left_on="mut_peptide", right_on="peptide").drop(columns="peptide")
    df = df.merge(wt, left_on=["wt_peptide", "allele"], right_on=["peptide", "allele"], how="left")
    df = df.drop(columns="peptide")
    df = df.merge(variants, on=VARIANT_KEY, how="left")

    df["delta_affinity_nM"] = df["wt_affinity_nM"] - df["mut_affinity_nM"]
    df["delta_percentile"] = df["wt_percentile"] - df["mut_percentile"]
    df["binder_strong"] = df["mut_affinity_nM"] < config.STRONG_BINDER_NM
    df["binder_weak"] = ~df["binder_strong"] & (df["mut_affinity_nM"] < config.WEAK_BINDER_NM)
    df["better_than_wt"] = df["mut_affinity_nM"] < df["wt_affinity_nM"]
    df["low_vaf"] = df["median_tumor_vaf"] < config.LOW_VAF_THRESHOLD
    df["binder_rank"] = np.select([df["binder_strong"], df["binder_weak"]], [0, 1], default=2)

    cols = [c for c, _ in sort_order]
    asc = [a for _, a in sort_order]
    df = df.sort_values(cols, ascending=asc, na_position="last", kind="mergesort")
    return df[OUTPUT_COLUMNS].reset_index(drop=True)
