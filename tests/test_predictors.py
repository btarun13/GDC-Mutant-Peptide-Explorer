"""Backend checks against the real tools. Run with: pytest -m models
Needs `mhcflurry-downloads fetch` to have been run."""

import csv
import io
import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.models

CLASS1_PEPTIDES = ["SIINFEKL", "GILGFVFTL", "NLVPMVATV", "KTWGQYWQV", "YLQPRTFLL",
                   "LLFGYPVYV", "ELAGIGILTV", "AAAWYLWEV", "RMFPNAPYL", "FLPSDFFPSV"]


@pytest.fixture(scope="module")
def class1():
    from src.predict.mhcflurry_backend import MhcflurryClass1

    return MhcflurryClass1()


def test_class1_rejects_unsupported_alleles(class1):
    ok, bad = class1.validate_alleles(["HLA-A*02:01", "HLA-A*99:99"])
    assert ok == ["HLA-A*02:01"] and bad == ["HLA-A*99:99"]


def test_class1_matches_cli(class1):
    df = class1.predict(CLASS1_PEPTIDES, ["HLA-A*02:01", "HLA-B*07:02"])
    assert len(df) == 20
    assert set(df["allele"]) == {"HLA-A*02:01", "HLA-B*07:02"}
    cli = subprocess.run(
        [str(Path(sys.executable).parent / "mhcflurry-predict"), "--alleles", "HLA-A*02:01",
         "--peptides", *CLASS1_PEPTIDES, "--out", "/dev/stdout"],
        capture_output=True, text=True, check=True,
    ).stdout
    cli_rows = {r["peptide"]: float(r["mhcflurry_affinity"]) for r in csv.DictReader(
        io.StringIO(cli[cli.index("allele,"):]))}
    ours = df[df["allele"] == "HLA-A*02:01"].set_index("peptide")["affinity_nM"]
    for pep, aff in cli_rows.items():
        assert ours[pep] == pytest.approx(aff, rel=1e-3)

