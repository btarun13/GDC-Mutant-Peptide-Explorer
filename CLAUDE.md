# CLAUDE.md

Project memory for Claude Code. Read this first; load the `gdc-maf-peptide-binding` skill (SKILL.md) before touching the pipeline code.

Last verified against live docs: 2026-10-13.

## What this project is

A Streamlit web app. The user picks a disease; the app finds matching projects on the NCI Genomic Data Commons (GDC), downloads their open-access somatic MAF files through the GDC API, rebuilds mutant peptides around each protein-changing mutation (9-mers by default, 8-11 allowed), and predicts HLA class I binding with MHCflurry for the mutant peptide and its normal (wild-type) counterpart. Output is a ranked, downloadable table.

Research use only. Binding is not immunogenicity, and nothing here supports clinical decisions. Say so in the UI.

## Pipeline at a glance

1. Disease text -> GDC project IDs (`/projects`, matched on `name`, `disease_type`, `primary_site`).
2. Project IDs -> open MAF file list (`/files`, filtered to Masked Somatic Mutation, access open).
3. Download `.maf.gz` (`/data/<file_id>`), verify md5, cache on disk.
4. Parse MAF -> keep protein-changing rows -> `Variant` objects.
5. Fetch the affected protein sequence (Ensembl REST by `ENSP`), cache on disk.
6. Check the reference amino acid at the reported position, substitute, slide a window of length k.
7. Predict binding for mutant and wild-type peptides against user-chosen HLA alleles.
8. Rank, show a funnel (variants -> peptides -> binders), export CSV plus a provenance sidecar.

## Layout

```
app.py                       # Streamlit entry point. UI only, no business logic.
src/pipeline.py              # stage functions the UI calls in order; funnel counts
src/provenance.py            # provenance record for views and exports
src/gdc/http.py              # requests session with timeouts and retry/backoff
src/gdc/client.py            # GDC API calls: projects, file search, download
src/gdc/disease_map.py       # free-text disease -> project_ids
src/maf/io.py                # read .maf.gz, validate columns, filter classes
src/maf/variants.py          # Variant dataclass built from MAF rows
src/protein/source.py        # ENSP -> protein sequence, disk cache
src/peptides/generate.py     # mutant and wild-type windows (pure functions)
src/predict/base.py          # Predictor protocol
src/predict/mhcflurry_backend.py
src/rank/score.py            # flags and sorting, thresholds from config
src/viz/charts.py            # Altair charts: funnel, mutant vs wild-type scatter, recurrence
src/config.py                # thresholds, window sizes, default alleles
tests/
data/cache/                  # gitignored: MAFs, protein FASTA, predictions
```

## Commands

```
python3.11 -m venv .venv && source .venv/bin/activate
python -m pip install -r requirements.txt
mhcflurry-downloads fetch          # one-time model download, large
streamlit run app.py
pytest -q -m "not network and not models"   # fast, offline
pytest -q -m network               # live GDC/Ensembl/UniProt calls
pytest -q -m models                # real MHCflurry checks (needs models)
ruff check .
```

## GDC rules

- Base URL `https://api.gdc.cancer.gov`. Endpoints used: `/projects`, `/files`, `/data/<uuid>`.
- Use only open-access data. Never request controlled files, never ask users for a GDC token in v1.
- Open MAF filter: `data_type = "Masked Somatic Mutation"`, `access = "open"`, `data_format = "MAF"`, `cases.project.project_id in [...]`. Show `analysis.workflow_type` to the user; the TCGA default is "Aliquot Ensemble Somatic Variant Merging and Masking".
- Use POST with a JSON body for file queries (GET has URL length limits). Page with `size` and `from`; read `pagination.total` instead of assuming a count.
- Never hardcode project lists, counts or file IDs. They change with each GDC data release. Record the `data_release` of each file in the output.
- Open MAFs are one file per tumour aliquot (checked 2026-10-13: TCGA-SKCM has 472 files of ~180 KB, all "Aliquot Ensemble Somatic Variant Merging and Masking"). A disease cohort is hundreds of small files, so the UI caps how many are used (`DEFAULT_MAX_MAF_FILES`) and downloads them in parallel with a disk cache.
- Do not scrape the portal website. Cache every download on disk; re-use before re-fetching.

## Open-access MAF caveats (these shape what the app can claim)

- Masked, not complete. The GDC states the germline-masking criteria are deliberately aggressive and recommends protected MAFs if missing true somatic mutations is a concern. Show this note in the UI.
- Normal-sample read counts and matched-normal alleles are blanked in open MAFs (`n_ref_count`, `n_alt_count`, `Match_Norm_Seq_Allele1/2`). Normal VAF cannot be computed. Tumour VAF (`t_alt_count / t_depth`) can.
- Coordinates are GRCh38. Gene models follow GENCODE v36 (listed in the GDC reference files). Never mix in hg19 data.
- MAFs carry no HLA type. The user must choose alleles (see Prediction rules).

## MAF -> peptide rules

- Use the row's own transcript: `Transcript_ID` and `ENSP`. Position numbers are only valid for that transcript.
- v1 handles missense only (`Variant_Classification == "Missense_Mutation"`). Phase 2: in-frame indels. Phase 3: frameshifts (need CDS plus the genomic change; `HGVSp_Short` alone is not enough). Nonsense, silent, intron, UTR and intergenic rows yield no new peptide: skip, but count them in the funnel.
- Parse `HGVSp_Short` (`p.V600E`) and cross-check against `Protein_position` ("600/766" style) and `Amino_acids` ("V/E"). If they disagree, skip the row and log it.
- Reference check is mandatory: the residue at the reported position in the fetched protein must equal the reference amino acid. On mismatch, skip and count it. The Ensembl REST service serves its current release, which can differ from the GENCODE v36 annotation the GDC used, so expect a few mismatches and report the rate.
- Window math for a 1-based position `pos`, length `k`: starts run from `max(1, pos-k+1)` to `min(pos, len(protein)-k+1)`. A mid-protein missense gives exactly `k` windows (nine 9-mers for k=9).
- Always keep the wild-type peptide from the same window next to the mutant one. Drop pairs where they are identical, and any peptide containing non-standard letters (`X`, `*`, `U`, `B`, `Z`).
- Deduplicate on (protein, position, ref, alt, peptide) and keep a count of distinct patients (`case_id`) per variant.

## Prediction rules

- Predictor choice is fixed: MHCflurry for MHC class I. Do not add or swap in other predictors (NetMHCpan, IEDB API, etc.) without discussing it in an issue first.
- MHC class II is out of scope: no class II alleles, peptide lengths, predictors (MHCnuggets was removed on 2026-10-13) or UI. Do not add it back without discussing it in an issue first.
- Class I (default, 9-mers; allow 8-11): MHCflurry, `Class1PresentationPredictor`.
- Allele names in the form `HLA-A*02:01`. Validate alleles against MHCflurry `predictor.supported_alleles` and show the user which were rejected.
- The default allele panel in `config.py` is a convenience for demos, not a patient HLA type. Label it that way in the UI and let the user replace it.
- Per-allele results: MHCflurry reports the best allele per sample, so pass `alleles={a: [a] for a in alleles}` to get one row per allele.
- Score the mutant and wild-type peptide in the same call so both come from the same model.
- Thresholds (affinity nM, percentile rank) live in `config.py` as named constants. Commonly used cut-offs are under 50 nM strong, under 500 nM weak; confirm MHCflurry percentile cut-offs before hard-coding.
- Cache predictions on disk keyed by (predictor name, predictor version, allele, peptide).
- Cap job size in the UI (peptides x alleles) and show progress. Load the MHCflurry predictor once with `st.cache_resource`.

## Ranking

Report, do not over-claim. Per mutant peptide show: mutant and wild-type binding, rank difference, tumour VAF, mutant read count, mutation position inside the peptide (anchors versus T-cell-facing positions), and number of patients carrying the variant. Sorting is a transparent tuple sort driven by config, not an invented composite score. Gene expression (GDC open STAR-Counts files) is a planned phase 2 filter, not part of v1.

## Streamlit conventions

- `app.py` only wires widgets to functions in `src/`. Everything in `src/` is importable and testable without Streamlit.
- `@st.cache_data(ttl=...)` for GDC listings; `@st.cache_resource` for the predictor; widget state in `st.session_state` with prefixed keys.
- Long operations run inside `st.status` or a progress bar. Never leave the user with a frozen page.
- Every result view and every export carries provenance: GDC file UUIDs and names, `data_release`, genome build, predictor name and version, parameters, timestamp.
- Show the funnel: rows in MAF -> missense kept -> reference check passed -> peptides generated -> binders.

## Code conventions

- Python 3.10+, type hints, dataclasses, `pathlib`, `logging` (no `print`), pandas for tables.
- Network code is isolated in `src/gdc` and `src/protein` so tests can mock it. Set timeouts and retry with backoff on 5xx and 429.
- No secrets in the repo. Cache and downloads are gitignored.
- Never paste protein sequences from memory into code or tests. Fetch and cache them.
- `mhcflurry` is pinned to 2.3.14 (stable on PyPI 2026-10-13; it runs on PyTorch, not TensorFlow). Check PyPI before bumping.
- On macOS, if the checkout sits on an exFAT or other non-APFS drive, macOS writes `._*` AppleDouble files; a `._*.pth` in site-packages crashes Python at startup. Delete `._*` under `.venv` after installing (`find .venv -name '._*' -delete`), or keep the virtual environment on an APFS disk.

## Testing

- Offline unit tests use a short synthetic protein: window counts, edge positions near the N- and C-terminus, reference mismatch handling, wild-type pairing.
- A `network`-marked golden test traces BRAF V600E end to end: fetch the sequence, assert position 600 is V, assert nine mutant 9-mers each containing E, assert each wild-type partner differs only at that position. Check the sequence against UniProt P15056. The BRAF windows in the earlier planning notes were written from memory and must not be copied.
- Mock the GDC API with saved JSON responses. Keep one tiny real `.maf.gz` slice as a fixture.
- Cross-check a handful of predictions against the `mhcflurry-predict` command line.

## Do not

- Do not call binding predictions "neoantigens" or "immunogenic" without the caveat that binding is only the first gate.
- Do not claim normal-sample evidence from an open MAF.
- Do not add controlled-access downloads, patient-identifying data, or credentials.
- Do not hardcode thresholds, counts, release numbers or project IDs inside UI code.

## Settled scope decisions (open an issue before changing them)

1. MHC class II: decided, out of scope (removed 2026-10-13). Class I only, MHCflurry.
2. HLA input: decided, manual allele picker pre-filled with the labelled demo panel from `config.py`.
3. Hosting: decided, local only (cache in `data/cache/`).
4. v1 scope: decided, missense only.

## References

- GDC API search and retrieval: https://docs.gdc.cancer.gov/API/Users_Guide/Search_and_Retrieval/
- GDC API downloading files: https://docs.gdc.cancer.gov/API/Users_Guide/Downloading_Files/
- GDC MAF format and masking rules: https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/
- GDC reference files (GENCODE versions): https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files
- MHCflurry Python tutorial: https://openvax.github.io/mhcflurry/python_tutorial.html
- Ensembl REST sequence endpoint: https://rest.ensembl.org/documentation/info/sequence_id
