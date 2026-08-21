import argparse

def parse_args():
    parser = argparse.ArgumentParser(
        description="Oligo specificity analysis tool for isothermal RPA design validation."
    )

    # 1. Input Oligos
    parser.add_argument(
        "--oligos",
        required=True,
        help="Path to the input CSV file containing primer metadata (OligoID, Sequence, PairID, Type)."
    )
    
    # 2. Genome Specifications
    parser.add_argument(
        "--target",
        required=True,
        help="Target organism scientific name (e.g., 'Homo sapiens') or NCBI Taxon ID (e.g., 9606)."
    )
    
    parser.add_argument(
        "--background",
        nargs="+",
        required=True,
        help="Space-separated list of background organism names, Taxon IDs, or local FASTA paths."
    )
    
    # 3. Execution Infrastructure & Isolation Environment
    parser.add_argument(
        "--analysis-id",
        default=None,
        help="Unique custom identifier for the run folder. Automatically generated if omitted."
    )
    
    parser.add_argument(
        "--db-output",
        default="outputs/blast_dbs/",
        help="Root directory where the shared, dynamically hashed BLAST databases are cached."
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Number of computational CPU threads to allocate to external BLAST operations."
    )

    # 4. Alignment & Specificity Search Tuning
    parser.add_argument(
        "--sensitivity",
        default="Standard",
        choices=["Fast", "Standard", "Advanced", "Very Advanced"],
        help="BLAST search sensitivity preset dictating internal word sizes and seed match penalties."
    )

    parser.add_argument(
        "--evalue",
        type=float,
        default=1000.0,
        help="Expectation value threshold for filtering out statistically weak background hits."
    )

    parser.add_argument(
        "--max-amplicon-size",
        type=int,
        default=500,
        help="Maximum allowable distance in base pairs between opposite hits to flag an off-target risk."
    )

    # 5. External Network Infrastructure
    parser.add_argument(
        "--api-key",
        type=str,
        default=None,
        help="Optional NCBI API key to increase the gene lookup rate limit from 5 to 10 requests per second."
    )

    # 6. Thermodynamics & specificity tuning
    parser.add_argument(
        "--operating-temp-c",
        type=float,
        default=37.0,
        help="Operating temperature (°C) used for thermodynamic calculations."
    )

    parser.add_argument(
        "--theta-threshold",
        type=float,
        default=0.05,
        help="Fractional occupancy (theta) cutoff used to flag elevated sequestration/amplification risk."
    )

    parser.add_argument(
        "--max-3prime-mismatches",
        type=int,
        default=3,
        help="Maximum allowed mismatches in the primer 3' window for extension competency checks."
    )

    parser.add_argument(
        "--three-prime-window",
        type=int,
        default=5,
        help="Length (bp) of the primer 3' window used when counting mismatches (n_bases)."
    )

    # 7. Reaction condition overrides (maps to AssayConditions)
    parser.add_argument(
        "--monovalent-salt-mM",
        type=float,
        default=100.0,
        help="Monovalent salt concentration in mM (used for Tm/thermo calculations)."
    )

    parser.add_argument(
        "--divalent-salt-mM",
        type=float,
        default=14.0,
        help="Divalent cation (Mg2+) concentration in mM (used for Tm/thermo calculations)."
    )

    parser.add_argument(
        "--dntp-mM",
        type=float,
        default=1.8,
        help="dNTP concentration in mM (used for Tm/thermo calculations)."
    )

    parser.add_argument(
        "--oligo-conc-uM",
        type=float,
        default=0.5,
        help="Oligonucleotide (primer) concentration in μM used for occupancy calculations."
    )

    parser.add_argument(
        "--template-conc-uM",
        type=float,
        default=1e-6,
        help="Template concentration in μM used for occupancy calculations."
    )

    # 8. Plot control
    parser.add_argument(
        "--no-plot",
        dest="write_plot",
        action="store_false",
        default=True,
        help="Disable generation of graphical plots (turn off writing plots)."
    )

    return parser.parse_args()