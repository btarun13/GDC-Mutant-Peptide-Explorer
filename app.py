"""Streamlit entry point. UI wiring only; all logic lives in src/."""

from __future__ import annotations

import logging
import re

import pandas as pd
import streamlit as st

from src import config, pipeline
from src.gdc.client import GdcClient, MafFile
from src.gdc.disease_map import match_projects
from src.predict.base import PredictionCache
from src.protein.source import ProteinSource
from src.provenance import Provenance
from src.viz import charts

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")

st.set_page_config(page_title="GDC mutant peptide binding", layout="wide")


# --- Cached resources --------------------------------------------------------
@st.cache_data(ttl=config.GDC_LISTING_TTL_S, show_spinner="Loading GDC projects...")
def cached_projects() -> list[dict]:
    return GdcClient().list_projects()


@st.cache_data(ttl=config.GDC_LISTING_TTL_S, show_spinner="Counting MAF workflows...")
def cached_workflow_types(project_ids: tuple[str, ...]) -> dict[str, int]:
    return GdcClient().workflow_types(list(project_ids))


@st.cache_data(ttl=config.GDC_LISTING_TTL_S, show_spinner="Listing open MAF files...")
def cached_mafs(project_ids: tuple[str, ...], workflow_type: str) -> list[MafFile]:
    return GdcClient().find_open_mafs(list(project_ids), workflow_type)


@st.cache_resource(show_spinner="Loading MHCflurry models...")
def class1_predictor():
    from src.predict.mhcflurry_backend import MhcflurryClass1

    return MhcflurryClass1()


def parse_alleles(text: str) -> list[str]:
    return [a for a in re.split(r"[\s,;]+", text.strip()) if a]


def progress_cb(bar, label: str):
    def _cb(done: int, total: int) -> None:
        bar.progress(min(done / total, 1.0) if total else 1.0, text=f"{label}: {done}/{total}")

    return _cb


# --- Header ------------------------------------------------------------------
st.title("GDC somatic mutations -> mutant peptides -> HLA binding")
st.error(config.RESEARCH_DISCLAIMER)
st.info(config.MASKED_MAF_NOTE)

# --- 1. Disease -> projects ----------------------------------------------------
st.header("1. Disease")
query = st.text_input("Disease, site or project", key="np_query", placeholder="e.g. melanoma")
projects = cached_projects() if query else []
matches = match_projects(query, projects) if query else []
selected_projects: list[str] = []
if query and not matches:
    st.warning("No GDC project matched.")
elif matches:
    if matches[0].fuzzy:
        st.warning("No exact match. Closest projects:")
    st.dataframe(
        pd.DataFrame(
            [
                {"project": m.project_id, "name": m.name, "matched on": "; ".join(m.matched_on)}
                for m in matches
            ]
        ),
        hide_index=True,
        width="stretch",
    )
    default = [m.project_id for m in matches if m.program == "TCGA" and not m.fuzzy][:1]
    selected_projects = st.multiselect(
        "Projects to use", [m.project_id for m in matches], default=default, key="np_projects"
    )

# --- 2. Open MAF files -------------------------------------------------------------
chosen_mafs: list[MafFile] = []
workflow = ""
if selected_projects:
    st.header("2. Open-access MAF files")
    wf_counts = cached_workflow_types(tuple(selected_projects))
    if not wf_counts:
        st.warning(
            "No open Masked Somatic Mutation MAFs for these projects. Filter used: "
            f"data_type={config.MAF_DATA_TYPE}, data_format={config.MAF_DATA_FORMAT}, "
            f"access={config.MAF_ACCESS}."
        )
    else:
        options = sorted(wf_counts, key=lambda w: w != config.DEFAULT_WORKFLOW_TYPE)
        workflow = st.selectbox(
            "Variant-calling workflow (analysis.workflow_type)",
            options,
            format_func=lambda w: f"{w} ({wf_counts[w]} files)",
            key="np_workflow",
        )
        mafs = cached_mafs(tuple(selected_projects), workflow)
        total_mb = sum(m.file_size for m in mafs) / 1e6
        st.write(
            f"**{len(mafs)}** open MAF files ({total_mb:.1f} MB). GDC open MAFs are one file "
            "per tumour aliquot, so each file is roughly one patient sample."
        )
        n_files = st.number_input(
            "Number of MAF files to use (taken in file-name order)",
            min_value=1,
            max_value=min(len(mafs), config.MAX_MAF_FILES_LIMIT),
            value=min(len(mafs), config.DEFAULT_MAX_MAF_FILES),
            key="np_n_files",
        )
        chosen_mafs = mafs[: int(n_files)]
        st.caption(
            f"Selected {len(chosen_mafs)} files, "
            f"{sum(m.file_size for m in chosen_mafs) / 1e6:.1f} MB to download (cached after)."
        )

# --- 3. Prediction settings ----------------------------------------------------------
if chosen_mafs:
    st.header("3. Peptides and HLA alleles (MHC class I, MHCflurry)")
    st.caption(
        "The default allele list is a demo panel of common alleles, **not a patient HLA "
        "type**. MAFs carry no HLA typing; replace them with the alleles you want to study."
    )
    lengths = sorted(
        st.multiselect(
            "Peptide lengths",
            list(config.CLASS1_LENGTHS),
            default=[config.CLASS1_DEFAULT_LENGTH],
            key="np_c1_len",
        )
    )
    alleles = parse_alleles(
        st.text_area(
            "HLA class I alleles (HLA-A*02:01 style)",
            "\n".join(config.DEMO_CLASS1_ALLELES),
            key="np_c1_alleles",
        )
    )

    if st.button("Run", type="primary", disabled=not (lengths and alleles)):
        funnel = pipeline.Funnel()
        table = None
        rejected: list[str] = []
        versions: dict[str, str] = {}
        protein_source = ProteinSource()
        try:
            with st.status("Running pipeline...", expanded=True) as status:
                st.write(f"Downloading {len(chosen_mafs)} MAF files...")
                bar = st.progress(0.0)
                paths, errors = GdcClient().download_many(
                    chosen_mafs, progress=progress_cb(bar, "MAF files")
                )
                if errors:
                    st.warning(f"{len(errors)} files failed to download and were skipped.")
                used_mafs = [m for m in chosen_mafs if m.file_id in paths]

                st.write("Parsing MAFs...")
                variants = pipeline.load_variants(paths, funnel)

                st.write(f"Fetching {len({v.ensp for v in variants})} proteins from Ensembl...")
                bar = st.progress(0.0)
                proteins = pipeline.fetch_proteins(
                    variants, protein_source, funnel, progress_cb(bar, "Ensembl batches")
                )
                if funnel.proteins_fetch_failed:
                    st.warning(
                        f"Ensembl failed for {funnel.proteins_fetch_failed:,} proteins after "
                        "retries; their variants are counted as 'no protein'. Run again later "
                        "to fill the gaps (fetched proteins are cached)."
                    )

                st.write("Building mutant / wild-type peptide pairs...")
                peptides = pipeline.make_peptides(variants, proteins, lengths, funnel)

                predictor = class1_predictor()
                versions[predictor.name] = predictor.version
                ok, rejected = predictor.validate_alleles(alleles)
                size = pipeline.job_size(peptides, len(ok))
                if size > config.MAX_CLASS1_JOB:
                    st.error(
                        f"Job is {size:,} peptide x allele predictions, over the cap of "
                        f"{config.MAX_CLASS1_JOB:,}. Use fewer files, lengths or alleles."
                    )
                else:
                    st.write(f"Scoring {size:,} peptide x allele pairs "
                             f"with {predictor.name} {predictor.version}...")
                    bar = st.progress(0.0)
                    table = pipeline.score_peptides(
                        predictor, peptides, ok, PredictionCache(), funnel,
                        progress_cb(bar, "Predictions"),
                    )
                status.update(label="Done", state="complete")
        except Exception as exc:  # show, do not freeze
            logging.exception("Pipeline failed")
            st.error(f"Pipeline failed: {exc}")
        else:
            st.session_state["np_results"] = {
                "table": table,
                "funnel": funnel,
                "rejected": rejected,
                "provenance": Provenance(
                    disease_query=query,
                    project_ids=selected_projects,
                    workflow_type=workflow,
                    gdc_files=Provenance.files_from(used_mafs),
                    ensembl_release=protein_source.release(),
                    predictors=versions,
                    parameters={
                        "lengths": lengths,
                        "alleles": alleles,
                        "strong_binder_nM": config.STRONG_BINDER_NM,
                        "weak_binder_nM": config.WEAK_BINDER_NM,
                        "low_vaf_threshold": config.LOW_VAF_THRESHOLD,
                        "sort_order": config.SORT_ORDER,
                        "download_errors": errors,
                    },
                    funnel=funnel.to_dict(),
                ),
            }

# --- 4. Results ---------------------------------------------------------------------
results = st.session_state.get("np_results")
if results:
    st.header("4. Results")
    st.error(config.RESEARCH_DISCLAIMER)
    funnel: pipeline.Funnel = results["funnel"]
    prov: Provenance = results["provenance"]

    table = results["table"]

    st.subheader("MHC class I (MHCflurry)")
    if results["rejected"]:
        st.warning(f"Rejected alleles: {', '.join(results['rejected'])}")

    st.markdown("**Funnel: where variants drop out**")
    st.altair_chart(charts.funnel_chart(funnel.stages()), width="stretch")
    st.caption(
        "Top panel: MAF rows surviving each filter. Bottom panel: unique mutant / wild-type "
        "peptide pairs, and peptide x allele predictions under "
        f"{config.WEAK_BINDER_NM:g} nM. The panels have separate scales."
    )

    if table is not None:
        n_binders = int((table["binder_rank"] < 2).sum())
        n_patients_loaded = len({f["case_id"] for f in prov.gdc_files if f.get("case_id")})
        tab_scatter, tab_recur, tab_table = st.tabs(
            ["Mutant vs wild type", "Recurrence", "Table"]
        )
        with tab_scatter:
            if n_binders:
                st.altair_chart(charts.affinity_scatter(table), width="stretch")
                shown = min(n_binders, config.VIZ_MAX_POINTS)
                st.caption(
                    f"Each point is a mutant binder (peptide x allele); showing the top "
                    f"{shown:,} of {n_binders:,} by the table's sort order. Below the dashed "
                    "diagonal the mutant binds more tightly than wild type. Points far below "
                    "it at an anchor position usually mean the mutation creates a binder "
                    "while the T-cell-facing surface still looks like self. T-cell-facing "
                    "changes where both peptides bind present an altered surface on an "
                    "already-presented peptide. Neither is evidence of immunogenicity."
                )
            else:
                st.info("No predicted binders to plot.")
        with tab_recur:
            recur = charts.recurrence_data(table)
            if recur.empty:
                st.info(
                    "No variant with a predicted binder is shared by two or more of the "
                    f"{n_patients_loaded:,} patients loaded. Load more MAF files to look for "
                    "recurrent (shared) variants."
                )
            else:
                st.altair_chart(charts.recurrence_chart(table), width="stretch")
                st.caption(
                    f"{len(recur)} variants carried by two or more patients that give at least "
                    "one predicted binder, ranked by distinct patients carrying the exact "
                    "variant among the "
                    f"{n_patients_loaded:,} patients loaded, then by best mutant affinity. "
                    "Shared variants are candidates for shared (off-the-shelf) targets, but "
                    "patients' real HLA types are unknown here, so a shared variant is not "
                    "necessarily presented in every carrier."
                )
        with tab_table:
            only_binders = st.checkbox("Show binders only", value=True, key="np_binders_only")
            view = table[table["binder_rank"] < 2] if only_binders else table
            st.dataframe(view, hide_index=True, width="stretch")
            st.download_button(
                "Download CSV",
                table.to_csv(index=False),
                file_name="binding_classI.csv",
                mime="text/csv",
                key="np_dl_table",
            )
            st.caption(
                "Sorted by binder class, then mutant affinity, then wild-type minus mutant "
                "affinity, then patient count (see SORT_ORDER in src/config.py)."
            )
    st.caption(
        "Binding is the first gate only: expression, processing, T-cell repertoire and HLA "
        "loss are not scored."
    )

    st.subheader("Quality checks")
    st.write(
        f"Reference-check mismatches: {funnel.ref_mismatch_rows:,} rows "
        f"({100 * funnel.ref_mismatch_rate:.2f}%). Ensembl serves release "
        f"{prov.ensembl_release}; the GDC annotation follows {config.GDC_GENE_MODEL}."
    )
    with st.expander("Funnel details and skip reasons"):
        st.json(funnel.to_dict())
    with st.expander("Provenance"):
        st.json(prov.to_json())
    st.download_button(
        "Download provenance (JSON)",
        prov.to_json(),
        file_name="provenance.json",
        mime="application/json",
        key="np_dl_prov",
    )
