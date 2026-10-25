# specificity-eval

A command-line tool for assessing RPA (Recombinase Polymerase Amplification) primer cross-reactivity by combining sequence, thermodynamic and geometric evidence.

Developed as part of an MSc Bioinformatics research project at Queen Mary University of London. Full methodology: [link to thesis].

## Project rationale

RPA runs isothermally at around 37–42 °C, and sequence similarity alone is a poor predictor of whether an off-target interaction will actually amplify. This tool separates **productive off-target amplification** from **isolated binding events** (primer sequestration or "primer sinks") using three layers of evidence:

1. **Sequence specificity**: BLAST identifies candidate off-target binding sites.
2. **Thermodynamic binding**: nearest-neighbour calculations give the fractional occupancy (`θ`) of each primer–template interaction under RPA conditions.
3. **Amplification geometry**: the two primers must bind the same sequence in a convergent orientation, give an amplicon within the permitted size range, and have 3′ ends compatible with extension.

Interactions are assigned a risk tier based on all three layers, not on any single threshold. The same logic is run against the target genome to confirm that each assay should detect what it is meant to.

## Installation

Requires a Unix-like system (macOS, Linux, or Windows via WSL) and internet access to NCBI. Native Windows is not supported because the workflow depends on BLAST+ and samtools.

```bash
git clone https://github.com/madebyemmanuel-r/specificity-eval.git
cd specificity-eval
conda env create -f env.yml
conda activate env
```

`env.yml` pins Python 3.11, BLAST+ 2.17, samtools 1.23 and the Python dependencies. Check the external tools are on your path:

```bash
blastn -version
samtools --version
```

### NCBI API key (optional)

An [NCBI API key](https://www.ncbi.nlm.nih.gov/account/settings/) raises the lookup rate limit from 5 to 10 requests per second. Pass it with `--api-key YOUR_KEY`.

## Usage

```bash
python main.py \
  --oligos validation/MBV_Primers.csv \
  --target "Macacine alphaherpesvirus 1" \
  --background 10298 \
  --sensitivity Advanced
```

- **`--oligos`**: CSV with one row per primer and the columns `OligoID`, `Sequence`, `PrimerPair`, `Orientation` (`Forward` / `Reverse`), plus an optional `AssayName`. Common header variants (`pair_id`, `type`, `setid`, …) are matched automatically, and IUPAC ambiguity codes are accepted.
- **`--target`**: scientific name or NCBI taxon ID. The best RefSeq assembly is downloaded.
- **`--background`**: one or more scientific names, taxon IDs or local FASTA files (`.fna`, `.fa`, `.fasta`).

Run `python main.py --help` for all options.

### Defaults

| Option | Default | Meaning |
| --- | ---: | --- |
| `--operating-temp-c` | 37.0 | Reaction temperature (°C) |
| `--divalent-salt-mM` | 14.0 | Mg²⁺ (MgAc) |
| `--monovalent-salt-mM` | 100.0 | KAc |
| `--dntp-mM` | 1.8 | dNTPs |
| `--oligo-conc-uM` | 0.5 | Primer concentration |
| `--template-conc-uM` | 1e-6 | Template concentration (1 pM) |
| `--theta-threshold` | 0.05 | Minimum `θ` for a site to count as bound |
| `--max-amplicon-size` | 500 | Maximum off-target amplicon (bp) |
| `--three-prime-window` | 5 | 3′ bases checked for mismatches |
| `--max-3prime-mismatches` | 3 | Maximum 3′ mismatches for extension competence |
| `--sensitivity` | Standard | BLAST preset: Fast, Standard, Advanced, Very Advanced |
| `--threads` | 4 | BLAST threads |

## Graphical Interface (Prototype)

A prototype GUI wraps the command-line tool. `prototype-qt` launches `main.py` as a subprocess through `runner.py`, so run it from the repository root. 

It lets you:
* Upload an oligo CSV
* Set the target and background organisms
* Adjust the BLAST sensitivity, e-value, maximum amplicon size, and thread count
* Optionally set an analysis ID, BLAST database folder, and NCBI API key

It also shows a progress bar and live console output, lets you cancel a run, and displays the report and output files when it finishes.

### Interface Reference

| Interface | Run with | Notes |
| :--- | :--- | :--- |
| **Desktop (PySide6)** | `python prototype-qt.py` | Background organisms are added to a list; includes Open Report and Open Output Folder buttons. |

### Installation

```bash
pip install PySide6   
```

> **Note:** The GUI does not yet expose the thermodynamic and reaction-condition options (temperature, salts, concentrations, θ threshold, 3′ settings), so those run at their defaults. Use the command line if changes are needed. 

## Output

Each run writes to `outputs/<analysis_id>/`:

| File | Contents |
| --- | --- |
| `specificity_report.md` | Human-readable report with per-assay verdicts |
| `dual_primer_cross_amplification_risks.csv` | Off-target primer pairs classed Tier 1 or Tier 2 |
| `single_primer_sequestration_sinks.csv` | Off-target sites that bind a primer but cannot amplify |
| `on_target_amplification.csv` | Convergent pairs found on the target genome |
| `run_metadata.json` | Parameters, genomes used and Python version for the run |

A high `θ` shows favourable binding, **not** confirmed amplification. Interactions are classified as follows:

| Class | Criteria |
| --- | --- |
| **Tier 1: Severe Cross-Reactivity Risk** | Both primers bind the same background sequence in convergent orientation, with an amplicon within `--max-amplicon-size`, `θ` ≥ threshold for both, and at most 2 mismatches in each primer's 3′ window |
| **Tier 2: Plausible Amplification Risk** | As Tier 1, but at least one primer has more than 2 mismatches in its 3′ window (still within `--max-3prime-mismatches`) |
| **Sequestration site** | Appreciable `θ` with no convergent partner, so no amplification but possible loss of primer and reduced sensitivity |

Each assay then gets a verdict: **FAIL** if it has any Tier 1 event, **WARN** if it has only Tier 2 events, otherwise **PASS**. On-target detection is reported as Confirmed, Weak or Not detected.

## How `θ` is calculated

`θ` is the fraction of template sites occupied by primer at equilibrium, calculated from nearest-neighbour ΔH and ΔS for the primer–template duplex, including mismatches, and the primer and template concentrations (primer in excess). Mismatch parameters follow Oliveira Martins & Weber (2024) (AOP-MM-60 and AOP-OW04-69; see thesis).

## Validation and case studies

- **Thermodynamics:** checked against [UFMG CompareTm / reference] for ΔH, ΔS, ΔG, Tm and occupancy. Differences in Tm come from salt correction.
- **Case studies:** an example primer sets from published RPA assays are included in `validation/`. The Monkey B virus assays are covered by automated tests using saved BLAST output in `tests/fixtures/MBV/`.

| Study | Target | Primer file | Expected result |
| --- | --- | --- | --- | --- |
| Chen et al. (2023), *Zoonoses* | Monkey B virus | `MBV_Primers.csv` | UL27 cross-reacts with HSV-1 (Tier 1); US6 is specific |


## Testing

```bash
python -m pytest tests/ -v
```

Tests that call NCBI (gene annotation of Tier 1 hits and the live gene-lookup suite) need internet access. For an offline run of the unit tests:

```bash
python -m pytest tests/test_oligo_parsing.py tests/test_thermodynamics.py -v
```

See `tests/README.md` for the test layout.

## Limitations

This is a prioritisation tool, not a substitute for experimental validation. It does not model the full RPA enzymatic mechanism (recombinase/polymerase kinetics, inhibitors), results depend on genome assembly quality, each organism is represented by a single reference genome, classification depends on the chosen thresholds, and real reaction conditions may differ from the defaults.

## Future work

Integration with PrimedRPA v2, a structured RPA assay knowledge base with semi-automated literature curation, expanded RPA-specific thermodynamic parameters, HPC-scale execution, and a graphical interface (early prototypes in `prototype.py` and `prototype-qt.py`).

## Citation

## Licence

MIT 

## Acknowledgements

Thanks to Matt Higgins [MatthewHiggins2017].
