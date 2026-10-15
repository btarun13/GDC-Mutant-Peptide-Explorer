"""Turn MAF rows into Variant objects, with a reason recorded for every skipped row."""

from __future__ import annotations

import logging
import re
from collections import Counter
from dataclasses import dataclass

import pandas as pd

from src import config

log = logging.getLogger(__name__)

HGVSP_MISSENSE = re.compile(r"^p\.([A-Z])(\d+)([A-Z])$")


@dataclass(frozen=True)
class Variant:
    gene: str
    transcript_id: str
    ensp: str
    position: int  # 1-based, on this transcript's protein
    ref: str
    alt: str
    hgvsp: str
    case_id: str
    tumor_sample_barcode: str
    t_depth: int
    t_alt_count: int
    gdc_file_id: str

    @property
    def tumor_vaf(self) -> float | None:
        return self.t_alt_count / self.t_depth if self.t_depth > 0 else None

    @property
    def key(self) -> tuple[str, int, str, str]:
        return (self.ensp, self.position, self.ref, self.alt)


def _to_int(value: str) -> int | None:
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def parse_missense_row(row: pd.Series) -> tuple[Variant | None, str | None]:
    """Returns (variant, None) or (None, skip_reason)."""
    m = HGVSP_MISSENSE.match(row["HGVSp_Short"])
    if not m:
        return None, "hgvsp_unparsed"
    ref, pos, alt = m.group(1), int(m.group(2)), m.group(3)
    if ref not in config.STANDARD_AMINO_ACIDS or alt not in config.STANDARD_AMINO_ACIDS:
        return None, "non_standard_residue"
    if _to_int(row["Protein_position"].split("/")[0]) != pos:
        return None, "position_disagrees"
    if row["Amino_acids"] != f"{ref}/{alt}":
        return None, "amino_acids_disagree"
    if not row["ENSP"]:
        return None, "no_ensp"
    t_depth = _to_int(row["t_depth"])
    t_alt = _to_int(row["t_alt_count"])
    if not t_depth or t_alt is None:
        return None, "no_tumour_depth"
    return (
        Variant(
            gene=row["Hugo_Symbol"],
            transcript_id=row["Transcript_ID"],
            ensp=row["ENSP"].split(".")[0],
            position=pos,
            ref=ref,
            alt=alt,
            hgvsp=row["HGVSp_Short"],
            case_id=row["case_id"],
            tumor_sample_barcode=row["Tumor_Sample_Barcode"],
            t_depth=t_depth,
            t_alt_count=t_alt,
            gdc_file_id=row.get("gdc_file_id", ""),
        ),
        None,
    )


def extract_missense(df: pd.DataFrame) -> tuple[list[Variant], Counter]:
    """Keep supported classes and return variants plus a Counter of skip reasons."""
    skipped: Counter = Counter()
    keep = df["Variant_Classification"].isin(config.SUPPORTED_VARIANT_CLASSES)
    skipped["not_missense"] = int((~keep).sum())
    variants: list[Variant] = []
    for _, row in df[keep].iterrows():
        v, reason = parse_missense_row(row)
        if v is None:
            skipped[reason] += 1
            log.debug("Skipped %s %s: %s", row["Hugo_Symbol"], row["HGVSp_Short"], reason)
        else:
            variants.append(v)
    return variants, skipped
