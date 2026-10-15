"""Mutant and wild-type peptide windows around missense variants. Pure functions, no I/O."""

from __future__ import annotations

import statistics
from collections import Counter, defaultdict
from collections.abc import Iterator, Sequence
from dataclasses import asdict, dataclass

import pandas as pd

from src import config
from src.maf.variants import Variant


class RefMismatch(ValueError):
    pass


def missense_windows(
    protein: str, pos: int, ref: str, alt: str, k: int
) -> Iterator[tuple[int, str, str, int]]:
    """pos is 1-based. Yields (start, mutant_peptide, wildtype_peptide, offset_in_peptide),
    with start and offset 1-based."""
    if not (1 <= pos <= len(protein)) or protein[pos - 1] != ref:
        raise RefMismatch(f"{ref}{pos} not found in protein of length {len(protein)}")
    mutant = protein[: pos - 1] + alt + protein[pos:]
    lo = max(1, pos - k + 1)
    hi = min(pos, len(protein) - k + 1)
    for s in range(lo, hi + 1):
        yield s, mutant[s - 1 : s - 1 + k], protein[s - 1 : s - 1 + k], pos - s + 1


def is_standard(peptide: str) -> bool:
    return set(peptide) <= config.STANDARD_AMINO_ACIDS


def is_anchor(offset: int, k: int, anchor_positions: Sequence[int] | None) -> bool | None:
    if anchor_positions is None:
        return None
    resolved = {k + 1 + a if a < 0 else a for a in anchor_positions}
    return offset in resolved


@dataclass(frozen=True)
class PeptidePair:
    ensp: str
    position: int
    ref: str
    alt: str
    length: int
    start: int
    offset: int
    mut_peptide: str
    wt_peptide: str
    anchor_position: bool | None


@dataclass
class PeptideResult:
    pairs: pd.DataFrame  # one row per unique (ensp, position, ref, alt, mut_peptide)
    variants: pd.DataFrame  # one row per unique (ensp, position, ref, alt)
    counts: Counter  # rows_no_protein, rows_ref_mismatch, rows_ref_ok, pairs_*


def summarise_variants(variants: Sequence[Variant]) -> pd.DataFrame:
    groups: dict[tuple, list[Variant]] = defaultdict(list)
    for v in variants:
        groups[v.key].append(v)
    rows = []
    for (ensp, pos, ref, alt), vs in groups.items():
        vafs = [v.tumor_vaf for v in vs if v.tumor_vaf is not None]
        rows.append(
            {
                "ensp": ensp,
                "position": pos,
                "ref": ref,
                "alt": alt,
                "gene": vs[0].gene,
                "transcript_id": vs[0].transcript_id,
                "hgvsp": vs[0].hgvsp,
                "n_patients": len({v.case_id for v in vs}),
                "n_samples": len({v.tumor_sample_barcode for v in vs}),
                "median_tumor_vaf": statistics.median(vafs) if vafs else None,
                "max_t_alt_count": max(v.t_alt_count for v in vs),
                "gdc_file_ids": ";".join(sorted({v.gdc_file_id for v in vs})),
            }
        )
    return pd.DataFrame(rows)


def generate_pairs(
    variants: Sequence[Variant],
    proteins: dict[str, str],
    lengths: Sequence[int],
    anchor_positions: Sequence[int] | None = None,
) -> PeptideResult:
    counts: Counter = Counter()
    ok_variants: list[Variant] = []
    pairs: dict[tuple, PeptidePair] = {}
    checked: dict[tuple, bool] = {}
    for v in variants:
        protein = proteins.get(v.ensp)
        if protein is None:
            counts["rows_no_protein"] += 1
            continue
        if v.key not in checked:
            ref_ok = 1 <= v.position <= len(protein) and protein[v.position - 1] == v.ref
            checked[v.key] = ref_ok
            if ref_ok:
                for k in lengths:
                    for start, mut, wt, offset in missense_windows(
                        protein, v.position, v.ref, v.alt, k
                    ):
                        if mut == wt:
                            counts["pairs_identical"] += 1
                            continue
                        if not (is_standard(mut) and is_standard(wt)):
                            counts["pairs_non_standard"] += 1
                            continue
                        dedup = (v.ensp, v.position, v.ref, v.alt, mut)
                        pairs.setdefault(
                            dedup,
                            PeptidePair(
                                v.ensp, v.position, v.ref, v.alt, k, start, offset, mut, wt,
                                is_anchor(offset, k, anchor_positions),
                            ),
                        )
        if checked[v.key]:
            counts["rows_ref_ok"] += 1
            ok_variants.append(v)
        else:
            counts["rows_ref_mismatch"] += 1
    counts["pairs"] = len(pairs)
    pair_df = pd.DataFrame([asdict(p) for p in pairs.values()])
    if pair_df.empty:
        pair_df = pd.DataFrame(columns=list(PeptidePair.__dataclass_fields__))
    return PeptideResult(pairs=pair_df, variants=summarise_variants(ok_variants), counts=counts)
