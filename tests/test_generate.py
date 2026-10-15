import pytest

from src.maf.variants import Variant
from src.peptides.generate import RefMismatch, generate_pairs, is_anchor, missense_windows

PROTEIN = "A" * 20 + "V" + "A" * 20  # synthetic, 41 residues, V at position 21


def _variant(pos=21, ref="V", alt="E", case="c1", ensp="ENSP1", depth=100, alt_count=30):
    return Variant(
        gene="G", transcript_id="ENST1", ensp=ensp, position=pos, ref=ref, alt=alt,
        hgvsp=f"p.{ref}{pos}{alt}", case_id=case, tumor_sample_barcode=f"S-{case}",
        t_depth=depth, t_alt_count=alt_count, gdc_file_id="f1",
    )


def test_mid_protein_gives_k_windows():
    windows = list(missense_windows(PROTEIN, 21, "V", "E", 9))
    assert len(windows) == 9
    for start, mut, wt, offset in windows:
        assert len(mut) == len(wt) == 9
        assert "E" in mut and "V" in wt
        diffs = [i for i, (a, b) in enumerate(zip(mut, wt, strict=True)) if a != b]
        assert diffs == [offset - 1]
        assert PROTEIN[start - 1 : start + 8] == wt


def test_offsets_cover_every_position():
    offsets = sorted(o for *_, o in missense_windows(PROTEIN, 21, "V", "E", 9))
    assert offsets == list(range(1, 10))


@pytest.mark.parametrize("pos,expected", [(1, 1), (3, 3), (41, 1), (39, 3)])
def test_termini_give_fewer_windows(pos, expected):
    protein = "M" + "A" * 39 + "K"
    windows = list(missense_windows(protein, pos, protein[pos - 1], "W", 9))
    assert len(windows) == expected


def test_reference_mismatch_raises():
    with pytest.raises(RefMismatch):
        list(missense_windows(PROTEIN, 21, "L", "E", 9))
    with pytest.raises(RefMismatch):
        list(missense_windows(PROTEIN, 99, "V", "E", 9))


def test_protein_shorter_than_k_gives_nothing():
    assert list(missense_windows("AVA", 2, "V", "E", 9)) == []


def test_anchor_resolution():
    assert is_anchor(2, 9, (2, -1)) is True
    assert is_anchor(9, 9, (2, -1)) is True
    assert is_anchor(5, 9, (2, -1)) is False
    assert is_anchor(10, 10, (2, -1)) is True
    assert is_anchor(5, 15, None) is None


def test_generate_pairs_counts_and_dedup():
    variants = [
        _variant(case="c1"),
        _variant(case="c2", alt_count=10),  # same variant, second patient
        _variant(pos=5, ref="L", alt="E"),  # reference mismatch (position 5 is A)
        _variant(ensp="ENSP_MISSING"),
    ]
    res = generate_pairs(variants, {"ENSP1": PROTEIN}, [9], anchor_positions=(2, -1))
    assert res.counts["rows_no_protein"] == 1
    assert res.counts["rows_ref_mismatch"] == 1
    assert res.counts["rows_ref_ok"] == 2
    assert len(res.pairs) == 9  # deduplicated across the two patients
    assert res.pairs["anchor_position"].sum() == 2
    (row,) = res.variants.to_dict("records")
    assert row["n_patients"] == 2
    assert row["median_tumor_vaf"] == pytest.approx(0.2)
    assert row["max_t_alt_count"] == 30


def test_non_standard_and_identical_pairs_dropped():
    protein = "A" * 10 + "V" + "A" * 3 + "X" + "A" * 10
    res = generate_pairs([_variant(pos=11)], {"ENSP1": protein}, [9])
    assert res.counts["pairs_non_standard"] > 0
    assert all("X" not in p for p in res.pairs["mut_peptide"])
    res = generate_pairs([_variant(alt="V")], {"ENSP1": PROTEIN}, [9])
    assert res.counts["pairs_identical"] == 9
    assert res.pairs.empty


def test_multiple_lengths():
    res = generate_pairs([_variant()], {"ENSP1": PROTEIN}, [9, 15])
    assert sorted(res.pairs["length"].value_counts().items()) == [(9, 9), (15, 15)]
