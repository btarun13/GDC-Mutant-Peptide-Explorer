import json
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from src.gdc.client import GdcClient, MafFile, md5_of
from src.gdc.disease_map import match_projects

FIX = Path(__file__).parent / "fixtures"
PROJECTS = json.loads((FIX / "gdc_projects.json").read_text())["data"]["hits"]
FILES_PAGE = json.loads((FIX / "gdc_files_skcm_page.json").read_text())


def test_melanoma_matches_skcm_first_with_reason():
    matches = match_projects("melanoma", PROJECTS)
    assert matches and matches[0].project_id == "TCGA-SKCM"
    assert not matches[0].fuzzy
    assert any("melanoma" in m.lower() for m in matches[0].matched_on)


def test_typo_falls_back_to_fuzzy():
    matches = match_projects("melanomma", PROJECTS)
    assert matches and all(m.fuzzy for m in matches)
    assert "TCGA-SKCM" in [m.project_id for m in matches]


def test_empty_query():
    assert match_projects("  ", PROJECTS) == []


def _response(payload):
    r = MagicMock()
    r.json.return_value = payload
    r.raise_for_status.return_value = None
    return r


def test_find_open_mafs_pages_until_total():
    hit = FILES_PAGE["data"]["hits"][0]
    page = lambda n, total: {"data": {"hits": [hit] * n, "pagination": {"total": total}}}  # noqa: E731
    session = MagicMock()
    session.post.side_effect = [_response(page(2, 3)), _response(page(1, 3))]
    files = GdcClient(session=session).find_open_mafs(["TCGA-SKCM"])
    assert len(files) == 3
    assert session.post.call_count == 2
    body = session.post.call_args_list[1].kwargs["json"]
    assert body["from"] == "2"
    fields = {c["content"]["field"]: c["content"]["value"] for c in body["filters"]["content"]}
    assert fields["access"] == "open"
    assert fields["data_type"] == "Masked Somatic Mutation"
    f = files[0]
    assert f.project_id == "TCGA-SKCM" and f.workflow_type and f.data_release and f.case_id


def _maf(md5: str) -> MafFile:
    return MafFile("fid", "x.maf.gz", 3, md5, "46.0", "wf", "TCGA-SKCM", "case")


def test_download_uses_cache(tmp_path):
    cached = tmp_path / "fid.maf.gz"
    cached.write_bytes(b"abc")
    session = MagicMock()
    path = GdcClient(session=session).download_maf(_maf(md5_of(cached)), tmp_path)
    assert path == cached
    session.get.assert_not_called()


def test_download_md5_mismatch_removes_file(tmp_path):
    resp = MagicMock()
    resp.iter_content.return_value = [b"abc"]
    resp.__enter__.return_value = resp
    session = MagicMock()
    session.get.return_value = resp
    with pytest.raises(ValueError):
        GdcClient(session=session).download_maf(_maf("0" * 32), tmp_path)
    assert list(tmp_path.iterdir()) == []
