from pathlib import Path

import pandas as pd
import pytest

from src.maf.io import MafFormatError, read_maf, read_mafs
from src.maf.variants import extract_missense, parse_missense_row

FIXTURE = Path(__file__).parent / "fixtures" / "skcm_slice.maf.gz"


def test_read_real_slice():
    df = read_maf(FIXTURE)
    assert len(df) == 32
    assert set(df["NCBI_Build"]) == {"GRCh38"}


def test_extract_missense_from_real_slice():
    df = read_mafs({"file-1": FIXTURE})
    variants, skipped = extract_missense(df)
    assert skipped["not_missense"] == 9
    assert len(variants) + sum(v for k, v in skipped.items() if k != "not_missense") == 23
    v = variants[0]
    assert v.gene == "CHD5" and v.position == 1886 and (v.ref, v.alt) == ("P", "S")
    assert v.gdc_file_id == "file-1"
    assert v.tumor_vaf == pytest.approx(6 / 27)


def _row(**over):
    base = {
        "Hugo_Symbol": "BRAF", "HGVSp_Short": "p.V600E", "Protein_position": "600/766",
        "Amino_acids": "V/E", "ENSP": "ENSP00000493543.1", "Transcript_ID": "ENST1",
        "case_id": "c", "Tumor_Sample_Barcode": "b", "t_depth": "50", "t_alt_count": "20",
    }
    base.update(over)
    return pd.Series(base)


def test_parse_row_ok_strips_version():
    v, reason = parse_missense_row(_row())
    assert reason is None and v.ensp == "ENSP00000493543" and v.position == 600


@pytest.mark.parametrize(
    "over,reason",
    [
        ({"HGVSp_Short": "p.V600_K601delinsE"}, "hgvsp_unparsed"),
        ({"Protein_position": "601/766"}, "position_disagrees"),
        ({"Amino_acids": "V/D"}, "amino_acids_disagree"),
        ({"ENSP": ""}, "no_ensp"),
        ({"t_depth": "0"}, "no_tumour_depth"),
        ({"HGVSp_Short": "p.V600X", "Amino_acids": "V/X"}, "non_standard_residue"),
    ],
)
def test_parse_row_rejects(over, reason):
    v, got = parse_missense_row(_row(**over))
    assert v is None and got == reason


def test_wrong_build_refused(tmp_path):
    df = read_maf(FIXTURE)
    df["NCBI_Build"] = "GRCh37"
    p = tmp_path / "x.maf"
    df.to_csv(p, sep="\t", index=False)
    with pytest.raises(MafFormatError):
        read_maf(p)
