import math

import pandas as pd

from src.maf.variants import Variant
from src.peptides.generate import generate_pairs
from src.predict.base import RESULT_COLUMNS, PredictionCache, predict_with_cache
from src.rank.score import build_table

PROTEIN = "A" * 20 + "V" + "A" * 20


class FakePredictor:
    """Affinity depends on whether the peptide contains E, so mutants bind better."""

    name, version = "fake", "1.0"

    def __init__(self):
        self.calls = []

    def validate_alleles(self, alleles):
        return list(alleles), []

    def predict(self, peptides, alleles, progress=None):
        self.calls.append((list(peptides), list(alleles)))
        rows = [
            {
                "peptide": p, "allele": a,
                "affinity_nM": 20.0 if "E" in p else 5000.0,
                "percentile": 0.1 if "E" in p else 50.0,
                "presentation_score": math.nan, "presentation_percentile": math.nan,
            }
            for a in alleles
            for p in peptides
        ]
        return pd.DataFrame(rows, columns=list(RESULT_COLUMNS))


def _pairs():
    v = Variant("G", "ENST1", "ENSP1", 21, "V", "E", "p.V21E", "c1", "s1", 100, 5, "f1")
    return generate_pairs([v], {"ENSP1": PROTEIN}, [9], anchor_positions=(2, -1))


def test_cache_avoids_recompute(tmp_path):
    cache = PredictionCache(tmp_path)
    pred = FakePredictor()
    res = _pairs()
    peps = sorted(set(res.pairs["mut_peptide"]) | set(res.pairs["wt_peptide"]))
    df1, n1 = predict_with_cache(pred, peps, ["HLA-A*02:01"], cache)
    df2, n2 = predict_with_cache(pred, peps, ["HLA-A*02:01"], cache)
    assert n1 == len(peps) and n2 == 0
    assert len(pred.calls) == 1
    assert len(df1) == len(df2) == len(peps)


def test_build_table_flags_and_sort(tmp_path):
    res = _pairs()
    peps = sorted(set(res.pairs["mut_peptide"]) | set(res.pairs["wt_peptide"]))
    preds, _ = predict_with_cache(FakePredictor(), peps, ["HLA-A*02:01"], PredictionCache(tmp_path))
    table = build_table(res.pairs, res.variants, preds)
    assert len(table) == 9
    assert table["binder_strong"].all()
    assert table["better_than_wt"].all()
    assert (table["delta_affinity_nM"] == 4980.0).all()
    assert table["low_vaf"].all()  # VAF 0.05
    assert table["wt_peptide"].str.contains("V").all()

