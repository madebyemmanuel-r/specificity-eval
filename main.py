import argparse
import logging
import sys
from pathlib import Path
import datetime
import json

from NCBI import resolve_identity, fetch_reference_genome, extract_accession_labels
from build_blast_database import build_blast_db, generate_db_prefix
from log_setup import setup_logging
from args import parse_args
from QC import verify_ncbi_connectivity
from oligo_handling import parse_and_validate_oligos, parse_oligos
from analysis_env import initialise_analysis_environment, get_target_name_for_run
from BLAST import run_blast_searches
from specificity_analysis import run_specificity_pipeline
from gui_events import emit_gui_event


def main():
    # 1. Parse CLI arguments
    args = parse_args()
    setup_logging()
    
    logger = logging.getLogger(__name__) 

    logger.info("Initialising primer specificity analysis tool...")
    emit_gui_event("STAGE", "INITIALISING_PIPELINE")
    logger.info("================ RUN PARAMETERS ================")
    for key, value in vars(args).items():
        logger.info(f"  {key}: {value}")
    logger.info("================================================")

    # 2. Fail-fast: confirm NCBI is reachable before any analysis begins
    try:
        verify_ncbi_connectivity()
    except ConnectionError as exc:
        logger.error(f"Critical: {exc}")
        sys.exit(1)
    logger.info("NCBI connectivity verified successfully.")
    emit_gui_event("STAGE", "VERIFYING_NCBI_CONNECTIVITY")

    # 3. Parse and QC the input oligo CSV
    logger.info(f"Parsing and validating oligo profiles: {args.oligos}")
    emit_gui_event("STAGE", "VALIDATING_OLIGOS")
    oligos = parse_oligos(args.oligos)
    oligos_df = oligos

    # 4. Resolve target organism taxonomy and obtain best RefSeq assembly
    logger.info(f"Resolving target organism's identity: '{args.target}'")
    target_meta = resolve_identity(args.target)
    emit_gui_event("STAGE", "RESOLVING_TARGET_ORGANISM")

    # 5. Initialise an isolated output directory for this run
    logger.info("Initialising analysis environment...")
    emit_gui_event("STAGE", "INITIALISING_ANALYSIS_ENVIRONMENT")
    target_name_for_run = get_target_name_for_run(target_meta)
    analysis_id, run_dir = initialise_analysis_environment(target_name_for_run, args.analysis_id)
    emit_gui_event("ANALYSIS_ID", analysis_id)
    emit_gui_event("RUN_DIR", str(run_dir))
    setup_logging(run_dir) #creates the run-specific log file, in the dedicated output directory

    # 6. Fetch the target genome — used for on-target amplification confirmation
    # Maps accession IDs to FASTA header descriptions for cleaner reporting downstream
    logger.info("Fetching target genome...")
    emit_gui_event("STAGE", "FETCHING_TARGET_GENOME")
    ref_path = fetch_reference_genome(target_meta, output_dir="data/genomes", api_key = args.api_key)
    accession_labels = {}
    accession_labels.update(extract_accession_labels(ref_path["fasta_path"]))

    # 7. Build the on-target BLAST DB from the target genome only
    # This is kept separate from the background DB so the two analyses never collide
    target_prefix, target_db_title = generate_db_prefix(
        [args.target], base_dir="outputs/blast_dbs/target"
    )
    target_db_path = build_blast_db(
        fasta_inputs=[ref_path["fasta_path"]],
        output_prefix=target_prefix,
        db_title=target_db_title
    )
    logger.info(f"Target BLAST DB ready: {target_db_path}")
    emit_gui_event("STAGE", "TARGET_BLAST_DB_READY")

    # 8. Fetch background genomes and build the background BLAST DB
    # Background organisms are what cross-reactivity is checked against
    logger.info(f"Resolving background genomes: {args.background}")
    emit_gui_event("STAGE", "RESOLVING_BACKGROUND_GENOMES")
    background_fastas = []
    background_genome_metas = [] 
    
    for bg_entry in args.background:
        # Check if the input is a valid, existing local genomic file
        if Path(bg_entry).exists() and Path(bg_entry).suffix in {".fna", ".fa", ".fasta"}:
            logger.info(f"  Using local FASTA: {bg_entry}")
            background_fastas.append(bg_entry)
            background_genome_metas.append({
                "organism": Path(bg_entry).stem,
                "fasta_path": str(bg_entry)
            })
            accession_labels.update(extract_accession_labels(bg_entry))
        else:
            bg_meta   = resolve_identity(bg_entry)
            bg_genome = fetch_reference_genome(bg_meta, output_dir="data/genomes", api_key = args.api_key)
            background_fastas.append(bg_genome["fasta_path"])
            background_genome_metas.append(bg_genome) 
            logger.info(f"  → {bg_genome['organism']} | {bg_genome['common_name']} ({bg_genome['accession']})")
            accession_labels.update(extract_accession_labels(bg_genome["fasta_path"]))

    if not background_fastas:
        logger.error("No background genomes resolved. Cannot build background BLAST DB.")
        sys.exit(1)

    bg_prefix, bg_db_title = generate_db_prefix(args.background, base_dir="outputs/blast_dbs")
    bg_db_path = build_blast_db(
        fasta_inputs=background_fastas,
        output_prefix=bg_prefix,
        db_title=bg_db_title
    )
    logger.info(f"Background BLAST DB ready: {bg_db_path}")
    emit_gui_event("STAGE", "BACKGROUND_BLAST_DB_READY")

    # 9. BLAST primers against the background — detects off-target cross-reactivity
    logger.info("Running background BLASTn (off-target screening)...")
    emit_gui_event("STAGE", "RUNNING_BACKGROUND_BLAST")
    fwd_bg_path, rev_bg_path = run_blast_searches(
        oligos_df, args, run_dir, bg_db_path, hit_prefix="bg"
    )

    # 10. BLAST primers against the target — confirms productive on-target amplification
    logger.info("Running target BLASTn (on-target confirmation)...")
    emit_gui_event("STAGE", "RUNNING_TARGET_BLAST")

    # Backup the original permissive e-value used for off-target screening
    original_evalue = args.evalue
    
    # Temporarily set a stricter e-value for on-target validation 
    args.evalue = 10 
    
    try:
        fwd_target_path, rev_target_path = run_blast_searches(
            oligos_df, args, run_dir, target_db_path, hit_prefix="target"
        )
    finally:
        # Guarantee the original value is restored even if the BLAST run crashes
        args.evalue = original_evalue

    # 11. Run the full specificity and confirmation analysis
    logger.info("Initiating downstream specificity evaluation pipeline...")
    emit_gui_event("STAGE", "INITIATING_SPECIFICITY_PIPELINE")
    summary_report = run_specificity_pipeline(
        fwd_csv=fwd_bg_path,
        rev_csv=rev_bg_path,
        fwd_target_csv=fwd_target_path,
        rev_target_csv=rev_target_path,
        output_dir=run_dir,
        alignment_col="Advanced Cross Reactivity Percentage String",
        n_bases=args.three_prime_window,
        max_amplicon_size=args.max_amplicon_size,
        max_3prime_mismatches=args.max_3prime_mismatches,
        theta_threshold=args.theta_threshold,
        operating_temp_c=args.operating_temp_c,
        write_plot=args.write_plot,
        api_key=args.api_key,
        target_organism=args.target,
        background_organisms=args.background,
        sensitivity=args.sensitivity,
        accession_labels=accession_labels,
        monovalent_salt_mM=args.monovalent_salt_mM,
        divalent_salt_mM=args.divalent_salt_mM,
        dntp_mM=args.dntp_mM,
        oligo_conc_uM=args.oligo_conc_uM,
        template_conc_uM=args.template_conc_uM,
    )

    #12. Write a single consolidated analysis metadata file covering run-specific data
    logger.info("Analysis complete. Writing metadata JSON file for this run...")
    emit_gui_event("STAGE", "WRITING_METADATA")
    run_metadata = {
        "pipeline_execution": {
            "analysis_id": analysis_id,
            "timestamp": datetime.datetime.now().isoformat(),
            "python_version": sys.version.split()[0],
            "command_line_args": {k: str(v) for k, v in vars(args).items()},
            "workspace_dir": str(run_dir)
        },
        "assay_parameters": {
            "operating_temperature_c": args.operating_temp_c,
            "max_amplicon_size_bp": args.max_amplicon_size,
            "max_3prime_mismatches_allowed": args.max_3prime_mismatches,
            "theta_risk_threshold": args.theta_threshold,
            "three_prime_window_bp": args.three_prime_window,
            "reaction_conditions": {
                "monovalent_salt_mM": args.monovalent_salt_mM,
                "divalent_salt_mM": args.divalent_salt_mM,
                "dntp_mM": args.dntp_mM,
                "oligo_conc_uM": args.oligo_conc_uM,
                "template_conc_uM": args.template_conc_uM,
            },
            "thermodynamic_units": {
                "enthalpy": "kcal/mol",
                "entropy": "cal/(mol·K)"
            },
            "sensitivity_mode": args.sensitivity
        },
        "target_organism": {
            "organism": ref_path["organism"],
            "common_name": ref_path.get("common_name", "N/A"),
            "tax_id": ref_path.get("tax_id"),
            "accession": ref_path["accession"],
            "assembly_level": ref_path["assembly_level"],
            "seq_length": ref_path["seq_length"],
            "local_fasta_path": ref_path["fasta_path"],
            "blast_db_path": target_db_path
        },
        "background_organisms": [
            {
                "organism": bg.get("organism", "Unknown"),
                "common_name": bg.get("common_name", "N/A"),
                "tax_id": bg.get("tax_id"),
                "accession": bg.get("accession", "Local Override"),
                "assembly_level": bg.get("assembly_level", "N/A"),
                "seq_length": bg.get("seq_length", 0),
                "local_fasta_path": bg.get("fasta_path")
            }
            for bg in background_genome_metas
        ]
    }
    
    output_json_path = Path(run_dir) / "run_metadata.json"
    with open(output_json_path, "w") as f:
        json.dump(run_metadata, f, indent=4)
        
    logger.info(f"Comprehensive run metadata written to: {output_json_path}")
    logger.info(f"Specificity pipeline complete. Report: {summary_report['report_path']}")
    emit_gui_event("REPORT_PATH", str(summary_report['report_path']))
    emit_gui_event("METADATA_PATH", str(output_json_path))
    emit_gui_event("STAGE", "PIPELINE_COMPLETE")
    logger.info("Primer specificity evaluation complete. Results saved to workspace.")

if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        logging.basicConfig(
            level=logging.ERROR,
            format="%(asctime)s - %(levelname)s - %(message)s",
            stream=sys.stdout,
        )
        logging.error("ERROR::Unhandled exception occurred during pipeline execution")
        logging.exception(exc)
        print(f"ERROR::{exc}")
        sys.exit(1)
