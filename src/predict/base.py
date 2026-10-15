"""Predictor protocol, shared result schema, and an on-disk prediction cache."""

from __future__ import annotations

import hashlib
import logging
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Protocol

import pandas as pd

from src import config

log = logging.getLogger(__name__)

# Every backend returns exactly these columns. Missing values are NaN, never invented.
RESULT_COLUMNS = (
    "peptide",
    "allele",  # as requested, "HLA-A*02:01" style
    "affinity_nM",
    "percentile",
    "presentation_score",
    "presentation_percentile",
)

Progress = Callable[[int, int], None]


class Predictor(Protocol):
    name: str
    version: str

    def validate_alleles(self, alleles: Sequence[str]) -> tuple[list[str], list[str]]:
        """Return (accepted, rejected)."""
        ...

    def predict(
        self, peptides: Sequence[str], alleles: Sequence[str], progress: Progress | None = None
    ) -> pd.DataFrame:
        """Return a DataFrame with RESULT_COLUMNS, one row per (peptide, allele)."""
        ...


class PredictionCache:
    """One parquet file per (predictor, version, allele), keyed by peptide."""

    def __init__(self, cache_dir: Path | None = None):
        self.cache_dir = cache_dir or config.PREDICTION_CACHE_DIR
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def _path(self, predictor: str, version: str, allele: str) -> Path:
        key = hashlib.sha1(f"{predictor}|{version}|{allele}".encode()).hexdigest()[:16]
        safe = allele.replace("*", "").replace(":", "").replace("/", "_")
        return self.cache_dir / f"{predictor}_{safe}_{key}.parquet"

    def get(self, predictor: str, version: str, allele: str) -> pd.DataFrame:
        p = self._path(predictor, version, allele)
        if not p.exists():
            return pd.DataFrame(columns=list(RESULT_COLUMNS))
        return pd.read_parquet(p)

    def put(self, predictor: str, version: str, allele: str, df: pd.DataFrame) -> None:
        if df.empty:
            return
        old = self.get(predictor, version, allele)
        merged = pd.concat([old, df[list(RESULT_COLUMNS)]], ignore_index=True)
        merged = merged.drop_duplicates("peptide", keep="last")
        p = self._path(predictor, version, allele)
        tmp = p.with_suffix(".tmp")
        merged.to_parquet(tmp, index=False)
        tmp.replace(p)


def predict_with_cache(
    predictor: Predictor,
    peptides: Sequence[str],
    alleles: Sequence[str],
    cache: PredictionCache,
    progress: Progress | None = None,
) -> tuple[pd.DataFrame, int]:
    """Score peptides x alleles, using the cache. Returns (results, n_newly_computed)."""
    unique = sorted(set(peptides))
    frames = []
    todo: dict[str, list[str]] = {}
    for allele in alleles:
        cached = cache.get(predictor.name, predictor.version, allele)
        hit = cached[cached["peptide"].isin(unique)]
        frames.append(hit)
        missing = sorted(set(unique) - set(hit["peptide"]))
        if missing:
            todo[allele] = missing
    total = sum(len(v) for v in todo.values())
    done = 0
    for allele, missing in todo.items():

        def _step(n: int, _of: int, base: int = done) -> None:
            if progress:
                progress(base + n, total)

        fresh = predictor.predict(missing, [allele], progress=_step)
        cache.put(predictor.name, predictor.version, allele, fresh)
        frames.append(fresh)
        done += len(missing)
    if progress and total == 0:
        progress(1, 1)
    frames = [f for f in frames if not f.empty]
    if not frames:
        return pd.DataFrame(columns=list(RESULT_COLUMNS)), total
    return pd.concat(frames, ignore_index=True)[list(RESULT_COLUMNS)], total
