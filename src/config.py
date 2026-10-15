"""Project-wide settings. UI code reads these; it never hardcodes its own values."""

from __future__ import annotations

import os
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
CACHE_DIR = Path(os.environ.get("PREDICT_CACHE_DIR", PROJECT_ROOT / "data" / "cache"))
MAF_CACHE_DIR = CACHE_DIR / "maf"
PROTEIN_CACHE_DIR = CACHE_DIR / "protein"
PREDICTION_CACHE_DIR = CACHE_DIR / "predictions"

# --- Network -----------------------------------------------------------------
GDC_BASE_URL = "https://api.gdc.cancer.gov"
ENSEMBL_BASE_URL = "https://rest.ensembl.org"
HTTP_TIMEOUT_S = 60
HTTP_DOWNLOAD_TIMEOUT_S = 300
HTTP_RETRIES = 5
HTTP_BACKOFF_FACTOR = 1.0
GDC_PAGE_SIZE = 500
GDC_DOWNLOAD_WORKERS = 4
ENSEMBL_BATCH_SIZE = 50  # documented POST /sequence/id limit
# Stop calling Ensembl after this many failed batches in a row (it is likely down).
ENSEMBL_MAX_CONSECUTIVE_FAILURES = 3
GDC_LISTING_TTL_S = 6 * 60 * 60

# --- GDC data ------------------------------------------------------------------
MAF_DATA_TYPE = "Masked Somatic Mutation"
MAF_DATA_FORMAT = "MAF"
MAF_ACCESS = "open"
DEFAULT_WORKFLOW_TYPE = "Aliquot Ensemble Somatic Variant Merging and Masking"
GENOME_BUILD = "GRCh38"
GDC_GENE_MODEL = "GENCODE v36"
# Open MAFs are one file per tumour aliquot, so a cohort can be hundreds of files.
DEFAULT_MAX_MAF_FILES = 50
MAX_MAF_FILES_LIMIT = 2000

# --- Variants ----------------------------------------------------------------
SUPPORTED_VARIANT_CLASSES = ("Missense_Mutation",)  # v1 scope: missense only
STANDARD_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")
LOW_VAF_THRESHOLD = 0.10

# --- Peptide windows -----------------------------------------------------------
CLASS1_DEFAULT_LENGTH = 9
CLASS1_LENGTHS = (8, 9, 10, 11)
# Class I anchors: position 2 and the C-terminal residue (1-based, "-1" = last).
CLASS1_ANCHOR_POSITIONS = (2, -1)

# --- Alleles -------------------------------------------------------------------
# Demo convenience only. This is NOT a patient HLA type; the UI says so.
DEMO_CLASS1_ALLELES = (
    "HLA-A*01:01",
    "HLA-A*02:01",
    "HLA-A*03:01",
    "HLA-A*24:02",
    "HLA-B*07:02",
    "HLA-B*08:01",
)

# --- Binding thresholds ------------------------------------------------------
# Affinity cut-offs are the commonly used ones. Flags use affinity only, because
# MHCflurry percentile cut-offs are not confirmed.
STRONG_BINDER_NM = 50.0
WEAK_BINDER_NM = 500.0

# --- Job size caps (peptides x alleles) ---------------------------------------
MAX_CLASS1_JOB = 200_000

# --- Ranking -------------------------------------------------------------------
# Transparent tuple sort: (column, ascending). No composite score.
SORT_ORDER = (
    ("binder_rank", True),  # 0 strong, 1 weak, 2 non-binder
    ("mut_affinity_nM", True),
    ("delta_affinity_nM", False),  # wild type minus mutant; larger = mutant binds better
    ("n_patients", False),
)

# --- Charts ----------------------------------------------------------------------
VIZ_MAX_POINTS = 5000  # scatter shows the top-ranked binder rows; Altair's default row limit
VIZ_TOP_RECURRENT = 25
# Categorical slots 1-2 of the reference palette; validated as a pair for scatter use.
VIZ_COLOR_PRIMARY = "#2a78d6"
VIZ_COLOR_SECONDARY = "#eb6834"
# Funnel stages are counted in different units; the chart says which.
FUNNEL_STAGE_UNITS = {
    "Rows in MAFs": "MAF rows",
    "Missense rows": "MAF rows",
    "HGVSp consistent": "MAF rows",
    "Protein found": "MAF rows",
    "Reference check passed": "MAF rows",
    "Peptide pairs": "unique mutant / wild-type peptide pairs",
    "Binder rows": "peptide x allele rows under the weak-binder cut-off",
}

RESEARCH_DISCLAIMER = (
    "Research use only. Predicted MHC binding is only the first gate: it is not "
    "immunogenicity, and nothing here supports clinical decisions."
)
MASKED_MAF_NOTE = (
    "Open-access GDC MAFs are germline-masked with deliberately aggressive criteria, "
    "so some true somatic mutations are missing. Normal-sample read counts are blanked, "
    "so normal VAF cannot be computed."
)
