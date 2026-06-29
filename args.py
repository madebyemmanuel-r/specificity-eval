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

    return parser.parse_args()