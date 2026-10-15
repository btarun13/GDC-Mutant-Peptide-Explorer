"""Pipeline stages the UI calls in order. No Streamlit here; progress goes through callbacks."""

from __future__ import annotations

import logging
from collections import Counter
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from src import config
from src.maf.io import class_counts, read_mafs
from src.maf.variants import Variant, extract_missense
from src.peptides.generate import PeptideResult, generate_pairs
from src.predict.base import PredictionCache, Predictor, Progress, predict_with_cache
from src.protein.source import ProteinSource
from src.rank.score import build_table

log = logging.getLogger(__name__)


@dataclass
class Funnel:
    maf_rows: int = 0
    variant_classes: dict[str, int] = field(default_factory=dict)
    missense_rows: int = 0
    skipped_parse: dict[str, int] = field(default_factory=dict)
    parsed_rows: int = 0
    no_protein_rows: int = 0
    proteins_fetch_failed: int = 0  # Ensembl errors after retries; re-run to fill from cache
    ref_mismatch_rows: int = 0
    ref_ok_rows: int = 0
    unique_variants: int = 0
    peptide_pairs: int = 0
    scored_rows: int = 0
    binders: int = 0

    @property
    def ref_mismatch_rate(self) -> float:
        checked = self.ref_ok_rows + self.ref_mismatch_rows
        return self.ref_mismatch_rows / checked if checked else 0.0

    def stages(self) -> list[tuple[str, int]]:
        return [
            ("Rows in MAFs", self.maf_rows),
            ("Missense rows", self.missense_rows),
            ("HGVSp consistent", self.parsed_rows),
            ("Protein found", self.parsed_rows - self.no_protein_rows),
            ("Reference check passed", self.ref_ok_rows),
            ("Peptide pairs", self.peptide_pairs),
            ("Binder rows", self.binders),
        ]

    def to_dict(self) -> dict:
        d = dict(self.__dict__)
        d["ref_mismatch_rate"] = round(self.ref_mismatch_rate, 4)
        return d


def load_variants(paths: dict[str, Path], funnel: Funnel) -> list[Variant]:
    df = read_mafs(paths)
    funnel.maf_rows = len(df)
    funnel.variant_classes = class_counts(df) if len(df) else {}
    variants, skipped = extract_missense(df)
    funnel.missense_rows = funnel.maf_rows - skipped.pop("not_missense", 0)
    funnel.skipped_parse = dict(skipped)
    funnel.parsed_rows = len(variants)
    return variants


def fetch_proteins(
    variants: Sequence[Variant],
    source: ProteinSource,
    funnel: Funnel,
    progress: Progress | None = None,
) -> dict[str, str]:
    proteins = source.fetch_many({v.ensp for v in variants}, progress=progress)
    funnel.proteins_fetch_failed = len(source.failed)
    return proteins


def make_peptides(
    variants: Sequence[Variant],
    proteins: dict[str, str],
    lengths: Sequence[int],
    funnel: Funnel,
) -> PeptideResult:
    res = generate_pairs(
        variants, proteins, lengths, anchor_positions=config.CLASS1_ANCHOR_POSITIONS
    )
    c: Counter = res.counts
    funnel.peptide_pairs = int(c["pairs"])
    funnel.no_protein_rows = int(c["rows_no_protein"])
    funnel.ref_mismatch_rows = int(c["rows_ref_mismatch"])
    funnel.ref_ok_rows = int(c["rows_ref_ok"])
    funnel.unique_variants = len(res.variants)
    if funnel.ref_mismatch_rate > 0.05:
        log.warning("Reference mismatch rate is %.1f%%", 100 * funnel.ref_mismatch_rate)
    return res


def job_size(peptides: PeptideResult, n_alleles: int) -> int:
    pairs = peptides.pairs
    if pairs.empty:
        return 0
    unique = set(pairs["mut_peptide"]) | set(pairs["wt_peptide"])
    return len(unique) * n_alleles


def score_peptides(
    predictor: Predictor,
    peptides: PeptideResult,
    alleles: Sequence[str],
    cache: PredictionCache,
    funnel: Funnel,
    progress: Progress | None = None,
) -> pd.DataFrame:
    pairs = peptides.pairs
    if pairs.empty or not alleles:
        return build_table(pairs, peptides.variants, pd.DataFrame())
    # Mutant and wild-type go through the same call, so both use the same model.
    to_score = sorted(set(pairs["mut_peptide"]) | set(pairs["wt_peptide"]))
    preds, n_new = predict_with_cache(predictor, to_score, alleles, cache, progress)
    log.info("%s: %d predictions, %d newly computed", predictor.name, len(preds), n_new)
    table = build_table(pairs, peptides.variants, preds)
    funnel.scored_rows = len(table)
    funnel.binders = int((table["binder_rank"] < 2).sum())
    return table
