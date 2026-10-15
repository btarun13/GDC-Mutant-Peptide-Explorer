"""Read GDC .maf.gz files and validate the columns the pipeline needs."""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from src import config

log = logging.getLogger(__name__)

REQUIRED_COLUMNS = (
    "Hugo_Symbol",
    "Chromosome",
    "Start_Position",
    "Reference_Allele",
    "Tumor_Seq_Allele2",
    "Variant_Classification",
    "Tumor_Sample_Barcode",
    "case_id",
    "HGVSp_Short",
    "Transcript_ID",
    "ENSP",
    "Protein_position",
    "Amino_acids",
    "t_depth",
    "t_alt_count",
    "NCBI_Build",
)


class MafFormatError(ValueError):
    pass


def read_maf(path: Path) -> pd.DataFrame:
    df = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        compression="gzip" if str(path).endswith(".gz") else None,
        usecols=lambda c: c in REQUIRED_COLUMNS,
        dtype=str,
        keep_default_na=False,
    )
    missing = sorted(set(REQUIRED_COLUMNS) - set(df.columns))
    if missing:
        raise MafFormatError(f"{path.name} is missing columns: {missing}")
    builds = set(df["NCBI_Build"].unique()) - {""}
    if builds and builds != {config.GENOME_BUILD}:
        raise MafFormatError(f"{path.name} has genome build {sorted(builds)}, expected GRCh38")
    return df


def read_mafs(paths: dict[str, Path]) -> pd.DataFrame:
    """Concatenate MAFs, tagging each row with the GDC file_id it came from."""
    frames = []
    for file_id, path in paths.items():
        df = read_maf(path)
        df["gdc_file_id"] = file_id
        frames.append(df)
    if not frames:
        return pd.DataFrame(columns=[*REQUIRED_COLUMNS, "gdc_file_id"])
    return pd.concat(frames, ignore_index=True)


def class_counts(df: pd.DataFrame) -> dict[str, int]:
    return df["Variant_Classification"].value_counts().to_dict()
