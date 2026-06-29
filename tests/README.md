# Test Suite

This directory contains the pytest-based (for robust automated execution) test suite for the RPA Primer Specificity Evaluation Tool.
Tests are validated against published biological ground truth from peer-reviewed literature.

## Structure

```
tests/
  fixtures/          # Saved BLAST output CSVs from validated pipeline runs
    MBV/             # Monkey B virus validation (Chen et al. 2023)
  conftest.py        # Shared fixtures and pipeline setup
  test_oligo_parsing.py       # Unit tests: CSV parsing and header matching
  test_thermodynamics.py      # Unit tests: NN thermodynamics, theta, alignment parsing
  test_specificity_pipeline.py # Integration tests: MBV ground truth validation
  README.md
```

## Testing Architecture

```text
tests/
├── conftest.py                  # Global pytest configurations, command-line arguments (--api-key),
│                                # and shared mock pipelines for downstream verification.
├── test_lookup_gene.py          # Live integration suite: validates coordinate mapping, intronic 
│                                # boundaries, and throttling guidelines against NCBI E-utilities.
├── test_oligo_parsing.py        # Unit suite: evaluates candidate input formatting, malformed CSV records,
│                                # and header verification pipelines.
└── test_thermodynamics.py       # Unit suite: evaluates nearest-neighbor Gibbs free energy changes (ΔG),
                                 # theta (θ) off-target risks, and primer-alignment parsers.
```
## Running the Tests

### Install dependencies:
```bash
pip install pytest
```

### Run the full suite:
```bash
python -m pytest tests/ -v
```

### Run a specific file:
```bash
python -m pytest tests/test_specificity_pipeline.py -v
```

### Run unit tests only (no pipeline execution):
```bash
python -m pytest tests/test_oligo_parsing.py tests/test_thermodynamics.py -v
```

## Biological Ground Truth

### Chen et al. 2023 — Monkey B Virus (MBV)
**Citation:** Chen et al. (2023). *Zoonoses*. DOI: 10.15212/ZOONOSES-2022-0040

Two RPA assays targeting Macacine alphaherpesvirus 1 were evaluated against
Human alphaherpesvirus 1 (HSV-1, NC_001806.2) as the background organism:

| Assay | Target Gene | Expected Result | Confirmed by Pipeline |
|-------|------------|-----------------|----------------------|
| UL27  | Glycoprotein B | Cross-reactive with HSV-1 | ❌ FAIL (Tier 1) |
| US6   | Glycoprotein D | Specific — no HSV-1 cross-reactivity | ✅ PASS |

## Adding New Validation Cases

1. Run the pipeline with the new primer set and correct background organisms
2. Save the BLAST fixture CSVs to `tests/fixtures/<STUDY_NAME>/`
3. Add a new test class to `test_specificity_pipeline.py` with a docstring citing the paper
4. Assert the known TP/TN outcomes from the publication