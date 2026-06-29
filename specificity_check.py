import logging
from pathlib import Path
from typing import Optional, Any
import pandas as pd

# Import mechanisms from the components file
from specificity_check_components import (
    AssayConditions,
    load_blast_csv,
    check_amplification_risk,
    evaluate_primer_sequestration_risks,
    plot_specificity_landscape
)

logger = logging.getLogger(__name__)

def run_specificity_pipeline(
    fwd_csv: str | Path,
    rev_csv: str | Path,
    output_dir: str | Path,
    alignment_col: str = "Advanced Cross Reactivity Percentage String",
    n_bases: int = 5,
    max_amplicon_size: int = 500,
    max_3prime_mismatches: int = 3,
    theta_threshold: float = 0.5,
    operating_temp_c: float = 37.0,
    write_plot: bool = True,
) -> dict[str, Any]:
    """
    Executes the end-to-end primer specificity assessment workflow.
    
    Workflow Steps:
      1. Loads and parses forward/reverse BLAST data blocks via core utilities.
      2. Computes duplex nearest-neighbor thermodynamics and fractional occupancy (theta).
      3. Performs an intersection analysis to trap convergent dual-primer amplification events.
      4. Queries NCBI Entrez annotations ONLY for verified Tier 1 structural risks.
      5. Isolates independent high-theta binder events showing active primer-sequestration profiles.
      6. Writes organized data metrics and landscape visualization diagnostic assets to disk.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    conditions = AssayConditions()

    # Step 1 & 2: Parse BLAST files and derive thermodynamic values in-memory
    logger.info("Scoring thermodynamic occupancy profiles for individual primer sets...")
    fwd_hits = load_blast_csv(fwd_csv, "forward", alignment_col, conditions, operating_temp_c, n_bases, max_3prime_mismatches)
    rev_hits = load_blast_csv(rev_csv, "reverse", alignment_col, conditions, operating_temp_c, n_bases, max_3prime_mismatches)

    # Step 3 & 4: Evaluate co-binding constraints across shared accessions
    shared_accessions = set(fwd_hits.keys()) & set(rev_hits.keys())
    pairwise_rows = []
    for accession in shared_accessions:
        for fwd_hit in fwd_hits[accession]:
            for rev_hit in rev_hits[accession]:
                pair = check_amplification_risk(fwd_hit, rev_hit, max_amplicon_size, theta_threshold, ncbi_email)
                if pair is not None:
                    # Convert data class to dictionary natively
                    from dataclasses import asdict
                    pairwise_rows.append(asdict(pair))

    pairwise_df = pd.DataFrame(pairwise_rows) if pairwise_rows else pd.DataFrame()
    pairwise_df["pair_id"] = pairwise_df["OligoID"].map(
    pd.Series(oligos_df.PairID.values, index=oligos_df.OligoID).to_dict()
    )

    # Step 5: Isolate independent off-target high-theta duplexes, which would lead to primer depletion
    sequestration_df = pd.concat([
        evaluate_primer_sequestration_risks(fwd_hits, theta_threshold),
        evaluate_primer_sequestration_risks(rev_hits, theta_threshold)
    ], ignore_index=True, sort=False).drop_duplicates(
        subset=["accession", "primer_direction", "back_seq_start", "back_seq_end", "theta"], errors="ignore"
    )

    # Step 6: Flush structured output matrices and diagnostic charts to disk
    logger.info(f"Writing diagnostic summaries and plots to: {output_dir.resolve()}")
    pairwise_df.to_csv(output_dir / "pairwise_risks.csv", index=False)
    sequestration_df.to_csv(output_dir / "sequestration_risks.csv", index=False)

    if write_plot:
        plot_specificity_landscape(pairwise_df, output_dir / "specificity_landscape.png")

    return {
        "pairwise_df": pairwise_df,
        "sequestration_df": sequestration_df,
        "plot_path": output_dir / "specificity_landscape.png" if write_plot else None
    }