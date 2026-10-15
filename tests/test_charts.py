from src import pipeline
from src.predict.base import PredictionCache, predict_with_cache
from src.rank.score import build_table
from src.viz import charts
from tests.test_rank_and_cache import FakePredictor, _pairs


def _table(tmp_path):
    res = _pairs()
    peps = sorted(set(res.pairs["mut_peptide"]) | set(res.pairs["wt_peptide"]))
    preds, _ = predict_with_cache(
        FakePredictor(), peps, ["HLA-A*02:01", "HLA-B*07:02"], PredictionCache(tmp_path)
    )
    return build_table(res.pairs, res.variants, preds)


def test_funnel_data_keeps_order_and_units():
    df = charts.funnel_data(pipeline.Funnel(maf_rows=10, missense_rows=4).stages())
    assert list(df["stage"])[:2] == ["Rows in MAFs", "Missense rows"]
    assert list(df["order"]) == list(range(len(df)))
    assert (df["unit"] != "").all()  # every stage has a unit in config


def test_scatter_data_splits_anchor_and_caps(tmp_path):
    table = _table(tmp_path)
    df = charts.scatter_data(table)
    assert len(df) == 18  # 9 windows x 2 alleles, all binders
    # 9-mer anchors are P2 and P9: 2 of the 9 offsets, for each allele.
    assert (df["position_class"] == charts.ANCHOR_LABEL).sum() == 4
    assert (df["position_class"] == charts.TCR_LABEL).sum() == 14
    assert len(charts.scatter_data(table, max_points=5)) == 5


def test_recurrence_data_only_shared_variants(tmp_path):
    table = _table(tmp_path)
    assert charts.recurrence_data(table).empty  # single-patient variant is not recurrent
    table["n_patients"] = 3
    df = charts.recurrence_data(table)
    assert len(df) == 1
    row = df.iloc[0]
    assert row["variant"] == "G p.V21E"
    assert row["n_patients"] == 3 and row["n_alleles"] == 2 and row["n_peptides"] == 9


def test_charts_build(tmp_path):
    table = _table(tmp_path)
    for chart in (
        charts.funnel_chart(pipeline.Funnel().stages()),
        charts.affinity_scatter(table),
        charts.recurrence_chart(table.assign(n_patients=2)),
    ):
        assert chart.to_dict()  # valid Vega-Lite spec
