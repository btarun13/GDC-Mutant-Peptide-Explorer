"""Live tests against GDC, Ensembl and UniProt. Run with: pytest -m network"""

import pytest
import requests

from src.gdc.client import GdcClient
from src.peptides.generate import missense_windows
from src.protein.source import ProteinSource

pytestmark = pytest.mark.network

BRAF_ENSP = "ENSP00000493543"  # BRAF MANE Select protein (ENST00000646891), used by the GDC


def _uniprot(accession: str) -> str:
    r = requests.get(f"https://rest.uniprot.org/uniprotkb/{accession}.fasta", timeout=30)
    r.raise_for_status()
    return "".join(line.strip() for line in r.text.splitlines() if not line.startswith(">"))


def test_braf_v600e_golden(tmp_path):
    seq = ProteinSource(cache_dir=tmp_path).fetch_many([BRAF_ENSP])[BRAF_ENSP]
    assert seq == _uniprot("P15056")
    assert seq[599] == "V"
    windows = list(missense_windows(seq, 600, "V", "E", 9))
    assert len(windows) == 9
    for _start, mut, wt, offset in windows:
        assert mut[offset - 1] == "E" and wt[offset - 1] == "V"
        assert sum(a != b for a, b in zip(mut, wt, strict=True)) == 1


def test_protein_cache_hit_makes_no_call(tmp_path):
    src = ProteinSource(cache_dir=tmp_path)
    src.fetch_many([BRAF_ENSP])

    class Boom:
        def post(self, *a, **k):
            raise AssertionError("network call on warm cache")

        get = post

    assert BRAF_ENSP in ProteinSource(cache_dir=tmp_path, session=Boom()).fetch_many([BRAF_ENSP])


def test_unknown_ensp_is_absent(tmp_path):
    out = ProteinSource(cache_dir=tmp_path).fetch_many([BRAF_ENSP, "ENSP99999999999"])
    assert list(out) == [BRAF_ENSP]


def test_gdc_skcm_has_open_mafs():
    client = GdcClient()
    projects = {p["project_id"] for p in client.list_projects()}
    assert "TCGA-SKCM" in projects
    wf = client.workflow_types(["TCGA-SKCM"])
    assert wf and all(n > 0 for n in wf.values())
