# Contributing

Thanks for your interest. Bug reports, questions and pull requests are welcome.

## Getting set up

```bash
python3.11 -m venv .venv && source .venv/bin/activate
python -m pip install -r requirements.txt
mhcflurry-downloads fetch          # only needed to run the app or the `models` tests
```

## Before opening a pull request

```bash
ruff check .
pytest -q -m "not network and not models"   # fast, offline; this is what CI runs
pytest -q -m network                        # optional: live GDC / Ensembl / UniProt calls
pytest -q -m models                         # optional: real MHCflurry predictions
```

Add or update tests for any change to parsing, peptide windows, prediction or ranking.

## Project rules

The design rules live in [CLAUDE.md](CLAUDE.md) (project rules) and [SKILL.md](SKILL.md)
(step-by-step pipeline notes). The ones most likely to matter for a pull request:

- `app.py` only wires widgets to functions in `src/`. Logic goes in `src/` and must be
  testable without Streamlit.
- Network code stays in `src/gdc` and `src/protein`, with timeouts and retries.
- Open-access GDC data only. No controlled-access downloads, tokens or credentials.
- Thresholds, limits and defaults go in `src/config.py`, not in UI code.
- Never paste protein sequences from memory into code or tests; fetch and cache them.
- Wording: predicted binding is not immunogenicity. Do not describe results as
  "neoantigens" or "immunogenic" without that caveat.

## Scope

MHC class I with MHCflurry, missense mutations. Planned next: in-frame indels, then
frameshifts, and a tumour-expression filter. Please open an issue to discuss before adding a
new predictor, MHC class II, or anything that changes the scope above.

## Reporting bugs

Include the disease or project, the number of MAF files, peptide lengths, alleles, and the
provenance JSON from the run if you have one.
