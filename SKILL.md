---
name: gdc-maf-peptide-binding
description: Use when building or changing the Streamlit app that takes a disease, pulls open-access somatic MAF files from the NCI GDC API, rebuilds mutant 8-11-mer peptides from protein-changing mutations, and predicts HLA class I binding with MHCflurry. MHC class II is out of scope. Covers GDC queries, MAF parsing, peptide windowing, prediction, ranking and verification.
---

# GDC MAF -> mutant peptides -> binding prediction

Follow the steps in order. Each step lists what to build, a reference snippet, and how to check it. Snippets are starting points: confirm field names and return columns against the live API or `df.columns` before relying on them. See CLAUDE.md for project rules and caveats.

## Step 1. Disease text -> GDC projects

The GDC has a small number of projects (78 when the docs example was captured), so fetch them all once, cache them, and match locally. This is simpler and more robust than guessing server-side filters.

```python
import requests

GDC = "https://api.gdc.cancer.gov"

def list_projects() -> list[dict]:
    r = requests.get(
        f"{GDC}/projects",
        params={"size": 1000, "fields": "project_id,name,disease_type,primary_site,program.name"},
        timeout=30,
    )
    r.raise_for_status()
    return r.json()["data"]["hits"]
```

Match the user's text (case-insensitive substring, then fuzzy) against `name`, each entry of `disease_type`, and each entry of `primary_site`. Show matches in a multiselect so the user confirms which `project_id`s to use (for example `TCGA-SKCM` for melanoma). If nothing matches, show the closest few and say so. Check available fields with `GET /projects/_mapping`.

Check: searching "melanoma" returns at least one project and the user can see why it matched.

## Step 2. Projects -> open MAF files

```python
def find_open_mafs(project_ids: list[str]) -> list[dict]:
    payload = {
        "filters": {"op": "and", "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": project_ids}},
            {"op": "=",  "content": {"field": "data_type",   "value": "Masked Somatic Mutation"}},
            {"op": "=",  "content": {"field": "data_format", "value": "MAF"}},
            {"op": "=",  "content": {"field": "access",      "value": "open"}},
        ]},
        "fields": ",".join([
            "file_id", "file_name", "file_size", "md5sum", "data_release",
            "analysis.workflow_type", "cases.project.project_id",
        ]),
        "format": "JSON",
        "size": "100",
    }
    r = requests.post(f"{GDC}/files", json=payload, timeout=60)
    r.raise_for_status()
    return r.json()["data"]["hits"]
```

- Open MAFs are one file per tumour aliquot (TCGA-SKCM had 472 on 2026-10-13), so a cohort is hundreds of small files. Cap how many are used in the UI, download in parallel, and cache. Show `analysis.workflow_type` (facet query) and default to "Aliquot Ensemble Somatic Variant Merging and Masking" when present.
- Read `pagination.total`; page with `from` if it exceeds `size`.
- Show file sizes before downloading.

## Step 3. Download and cache

```python
import hashlib
from pathlib import Path

def download_maf(file_id: str, md5: str, cache_dir: Path) -> Path:
    dest = cache_dir / f"{file_id}.maf.gz"
    if dest.exists() and _md5(dest) == md5:
        return dest
    with requests.get(f"{GDC}/data/{file_id}", stream=True, timeout=120) as r:
        r.raise_for_status()
        with open(dest, "wb") as fh:
            for chunk in r.iter_content(1 << 20):
                fh.write(chunk)
    if _md5(dest) != md5:
        dest.unlink()
        raise ValueError(f"md5 mismatch for {file_id}")
    return dest

def _md5(p: Path) -> str:
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()
```

Check: a second run with a warm cache makes no network call for the file.

## Step 4. Parse the MAF

```python
import pandas as pd

USECOLS = [
    "Hugo_Symbol", "Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2",
    "Variant_Classification", "Tumor_Sample_Barcode", "case_id", "HGVSp_Short",
    "Transcript_ID", "ENSP", "Protein_position", "Amino_acids", "t_depth", "t_alt_count",
    "NCBI_Build",
]

def read_maf(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", comment="#", compression="gzip",
                     usecols=lambda c: c in USECOLS, low_memory=False)
    missing = set(USECOLS) - set(df.columns)
    if missing:
        raise ValueError(f"MAF is missing columns: {sorted(missing)}")
    return df
```

Then:

- Keep `Variant_Classification == "Missense_Mutation"` (v1). Count every other class for the funnel.
- Parse `HGVSp_Short` with `^p\.([A-Z])(\d+)([A-Z])$` to get ref, position, alt.
- Cross-check `Protein_position` (text before the `/`) and `Amino_acids` (`V/E`). Skip and log rows that disagree.
- Tumour VAF = `t_alt_count / t_depth`; skip depth 0. Open MAFs have no normal counts.
- Confirm `NCBI_Build` is `GRCh38`; refuse otherwise.

Check: for a known project, the missense count is plausible against the class counts you print.

## Step 5. Protein sequences

```python
def fetch_protein(ensp: str, cache_dir: Path) -> str:
    cached = cache_dir / f"{ensp}.txt"
    if cached.exists():
        return cached.read_text().strip()
    r = requests.get(
        f"https://rest.ensembl.org/sequence/id/{ensp}",
        params={"type": "protein"},
        headers={"Content-Type": "text/plain"},
        timeout=30,
    )
    r.raise_for_status()
    seq = r.text.strip()
    cached.write_text(seq)
    return seq
```

- Strip any version suffix (`.3`) if the service rejects it; try both.
- Batch and throttle requests; keep a disk cache so each ENSP is fetched once.
- The service serves its current Ensembl release, while the GDC annotation follows GENCODE v36. The reference check in Step 6 catches disagreements; report how many rows it dropped.

## Step 6. Mutant and wild-type windows

```python
class RefMismatch(ValueError):
    pass

def missense_windows(protein: str, pos: int, ref: str, alt: str, k: int):
    """pos is 1-based. Yields (start, mutant_peptide, wildtype_peptide, offset_in_peptide)."""
    if not (1 <= pos <= len(protein)) or protein[pos - 1] != ref:
        raise RefMismatch(f"{ref}{pos} not found in protein")
    mutant = protein[: pos - 1] + alt + protein[pos:]
    lo = max(1, pos - k + 1)
    hi = min(pos, len(protein) - k + 1)
    for s in range(lo, hi + 1):
        yield s, mutant[s - 1 : s - 1 + k], protein[s - 1 : s - 1 + k], pos - s + 1
```

Window size:

| Pathway | Length k | Notes |
|---|---|---|
| MHC class I (CD8 T cells) | 9 by default; allow 8-11 | MHCflurry |

Rules:

- A mid-protein missense gives exactly `k` windows. Near the termini it gives fewer.
- Drop a pair if mutant == wild type, and any peptide with letters outside the 20 standard amino acids.
- Record `offset_in_peptide`. For a 9-mer, positions 2 and 9 are the usual anchors and positions 3-8 mostly face the T cell receptor; keep this as a flag, with the positions in config.
- Later phases: in-frame indels need the junction rebuilt from the `HGVSp_Short` range; frameshifts need the CDS sequence and the genomic change, translated until the new stop.

Check (offline test, synthetic protein): for `protein = "A" * 20 + "V" + "A" * 20`, `pos=21`, `ref="V"`, `alt="E"`, `k=9` gives 9 windows, each mutant peptide contains `E`, each wild-type partner contains `V`, and they differ only at the mutated offset.

Check (live, `network` marker): BRAF V600E. Fetch the BRAF protein, assert position 600 is `V`, assert nine 9-mers. Compare against UniProt P15056. Do not hardcode the sequence.

## Step 7. Binding prediction (class I)

```python
from mhcflurry import Class1PresentationPredictor

predictor = Class1PresentationPredictor.load()      # cache once (st.cache_resource)

def predict_class1(peptides: list[str], alleles: list[str]) -> pd.DataFrame:
    unsupported = [a for a in alleles if a not in predictor.supported_alleles]
    if unsupported:
        raise ValueError(f"Unsupported alleles: {unsupported}")
    df = predictor.predict(
        peptides=peptides,
        alleles={a: [a] for a in alleles},   # one result row per allele
        verbose=0,
        throw=False,
    )
    return df   # inspect df.columns; map names in one place
```

- Models must be downloaded once: `mhcflurry-downloads fetch`.
- `predict` returns a DataFrame with `peptide`, `sample_name`, `best_allele` and `presentation_score`, plus affinity and percentile columns. Confirm the exact names with `df.columns` on the installed version and translate them to your own schema (`affinity_nM`, `percentile`, `presentation_score`) in one adapter function.
- Optional flanking context: `n_flanks` and `c_flanks` take the residues around each peptide from the mutant protein; check the API reference for the expected flank length. Alternatively `predict_sequences` scans a mutant protein segment and adds flanks itself; filter the output to peptides that include the mutated position.
- MHC class II is out of scope: no class II alleles, peptide lengths or predictors.
- Run mutant and wild-type peptides in the same call and the same model version.
- Cache results on disk keyed by (predictor, version, allele, peptide). Cap peptides x alleles in the UI and show progress.

Check: 10 peptides through the app match `mhcflurry-predict --alleles ... --peptides ...` on the command line.

## Step 8. Combine, flag, rank

Build one row per (variant, peptide, allele) with:

| Column | Meaning |
|---|---|
| gene, protein change, transcript | where it came from |
| mutant / wild-type peptide, offset | the pair and the mutation position inside it |
| mutant / wild-type affinity (nM), percentile, presentation score | binding for both |
| delta percentile | wild type minus mutant; positive means the mutant binds better |
| tumour VAF, mutant read count | evidence the variant is real and common in the tumour |
| patients carrying variant | distinct `case_id` count across the loaded MAFs |

Flags: `binder_strong`, `binder_weak` (thresholds from config), `better_than_wt`, `anchor_position`, `low_vaf`.

Sort with a documented tuple (for example binder class, then percentile, then delta), configurable in `config.py`. Do not invent a composite score. Say plainly that binding is the first gate only: expression, processing, T cell availability and HLA loss all matter and are not scored here.

## Step 9. Verify before calling it done

- Funnel counts add up at every stage and are shown in the UI.
- Reference-mismatch rate is logged and shown; investigate if it is more than a few percent.
- The BRAF V600E trace passes (Step 6).
- A cached second run makes no GDC, Ensembl or model calls for unchanged inputs.
- Export includes provenance: GDC file UUIDs and names, `data_release`, genome build, predictor name and version, parameters, timestamp.
- The UI states: open-access MAFs are masked; normal read counts are unavailable; HLA alleles are user-chosen; research use only.

## Failure modes and what to show

| Symptom | Likely cause | What to do |
|---|---|---|
| Many reference mismatches | Ensembl release differs from GENCODE v36, or wrong transcript | Report the rate; try the unversioned ID; consider an archive release matching the GDC annotation |
| GDC returns zero MAFs | Filters too narrow or project has no open MAF | Show the filter used; list what the project does have |
| MHCflurry `ValueError` on an allele | Unsupported allele | Show which alleles were rejected |
| Very slow predictions | Too many peptides x alleles | Enforce the cap; batch; use the cache |
| 429 or 5xx from GDC or Ensembl | Rate limiting or outage | Back off and retry; keep partial results |

## Later phases

1. In-frame insertions and deletions.
2. Frameshifts (CDS-level reconstruction).
3. Expression filter from GDC open gene-expression (STAR-Counts) files, matching patients on barcode prefix and marking missing RNA as missing, not zero.
4. Patient mode from VCFs, which keeps phasing and nearby-variant context that a MAF loses.
