"""MHC class I predictions with MHCflurry's Class1PresentationPredictor."""

from __future__ import annotations

import logging
from collections.abc import Sequence
from importlib.metadata import version as pkg_version

import pandas as pd

from src.predict.base import RESULT_COLUMNS, Progress

log = logging.getLogger(__name__)

BATCH_SIZE = 5000

# MHCflurry output column -> our schema. Checked against mhcflurry 2.3.14.
COLUMN_MAP = {
    "peptide": "peptide",
    "affinity": "affinity_nM",
    "affinity_percentile": "percentile",
    "presentation_score": "presentation_score",
    "presentation_percentile": "presentation_percentile",
}


class MhcflurryClass1:
    name = "mhcflurry"

    def __init__(self, predictor=None):
        if predictor is None:
            from mhcflurry import Class1PresentationPredictor

            predictor = Class1PresentationPredictor.load()
        self._predictor = predictor
        self.version = pkg_version("mhcflurry")
        self._supported = set(predictor.supported_alleles)

    def validate_alleles(self, alleles: Sequence[str]) -> tuple[list[str], list[str]]:
        accepted, rejected = [], []
        for a in alleles:
            if a in self._supported:
                accepted.append(a)
            else:
                rejected.append(a)
        return accepted, rejected

    def predict(
        self, peptides: Sequence[str], alleles: Sequence[str], progress: Progress | None = None
    ) -> pd.DataFrame:
        accepted, rejected = self.validate_alleles(alleles)
        if rejected:
            raise ValueError(f"Not supported by MHCflurry class I: {rejected}")
        peptides = list(peptides)
        frames = []
        for i in range(0, len(peptides), BATCH_SIZE):
            batch = peptides[i : i + BATCH_SIZE]
            df = self._predictor.predict(
                peptides=batch,
                alleles={a: [a] for a in accepted},  # one sample per allele -> row per allele
                include_affinity_percentile=True,
                verbose=0,
                throw=True,
            )
            frames.append(self._to_schema(df))
            if progress:
                progress(min(i + BATCH_SIZE, len(peptides)), len(peptides))
        if not frames:
            return pd.DataFrame(columns=list(RESULT_COLUMNS))
        return pd.concat(frames, ignore_index=True)

    @staticmethod
    def _to_schema(df: pd.DataFrame) -> pd.DataFrame:
        missing = set(COLUMN_MAP) - set(df.columns)
        if missing:
            raise RuntimeError(f"MHCflurry output lacks expected columns: {sorted(missing)}")
        out = df.rename(columns=COLUMN_MAP)
        out["allele"] = df["sample_name"]
        return out[list(RESULT_COLUMNS)]
