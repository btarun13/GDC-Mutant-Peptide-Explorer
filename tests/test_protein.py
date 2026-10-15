import requests

from src import config, pipeline
from src.protein.source import ProteinSource


class FlakySession:
    """Batch POSTs fail for any batch containing ENSP1; everything else returns 'MSEQ'."""

    def __init__(self):
        self.posts = 0

    def post(self, url, params=None, json=None, **kwargs):
        self.posts += 1
        if "ENSP1" in json["ids"]:
            raise requests.exceptions.RetryError("too many 503 error responses")
        resp = requests.Response()
        resp.status_code = 200
        resp._content = (
            "[" + ",".join(f'{{"query": "{i}", "seq": "MSEQ"}}' for i in json["ids"]) + "]"
        ).encode()
        return resp


def test_failed_batch_keeps_partial_results(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "ENSEMBL_BATCH_SIZE", 1)
    src = ProteinSource(cache_dir=tmp_path, session=FlakySession())
    out = src.fetch_many(["ENSP1.4", "ENSP2.1", "ENSP3"])
    assert out == {"ENSP2": "MSEQ", "ENSP3": "MSEQ"}
    assert src.failed == {"ENSP1"}
    assert not (tmp_path / "ENSP1.txt").exists()  # failures are not cached, so a re-run retries


def test_funnel_counts_fetch_failures(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "ENSEMBL_BATCH_SIZE", 1)
    from src.maf.variants import Variant

    variants = [
        Variant("G", "ENST1", e, 2, "S", "A", "p.S2A", "c1", "s1", 10, 5, "f1")
        for e in ("ENSP1", "ENSP2")
    ]
    funnel = pipeline.Funnel()
    src = ProteinSource(cache_dir=tmp_path, session=FlakySession())
    proteins = pipeline.fetch_proteins(variants, src, funnel)
    assert list(proteins) == ["ENSP2"]
    assert funnel.proteins_fetch_failed == 1


def test_stops_calling_after_consecutive_failures(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "ENSEMBL_BATCH_SIZE", 1)
    monkeypatch.setattr(config, "ENSEMBL_MAX_CONSECUTIVE_FAILURES", 2)

    class DownSession(FlakySession):
        def post(self, url, params=None, json=None, **kwargs):
            self.posts += 1
            raise requests.exceptions.RetryError("too many 503 error responses")

    session = DownSession()
    src = ProteinSource(cache_dir=tmp_path, session=session)
    assert src.fetch_many([f"ENSP{i}" for i in range(10)]) == {}
    assert session.posts == 2
    assert len(src.failed) == 10
