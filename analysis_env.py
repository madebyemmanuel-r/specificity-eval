import os
import logging
from datetime import datetime
from pathlib import Path

logger = logging.getLogger(__name__)

def _normalise_analysis_id(raw_id: str) -> str:
    """Convert an arbitrary analysis identifier into a safe filesystem name.

    The input may contain spaces, punctuation, or other characters that are
    not safe or convenient for directory names. This helper normalises that
    string by replacing spaces with underscores and removing any character
    that is not a standard alphanumeric, hyphen, underscore, or dot.

    Returns a fallback name if the cleaned string would otherwise be empty.
    """
    normalised = raw_id.replace(" ", "_")
    allowed = set("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-.")
    normalised = "".join(c for c in normalised if c in allowed)
    return normalised or "analysis_run"


def get_target_name_for_run(target_meta) -> str:
    """Return a short filesystem-safe name for the selected target organism."""
    if isinstance(target_meta, dict):
        organism_info = target_meta.get("organism", {})
        name = (
            organism_info.get("organism_name")
            or organism_info.get("sci_name")
            or target_meta.get("accession")
            or target_meta.get("current_accession")
            or target_meta.get("common_name")
            or str(target_meta)
        )
    else:
        name = str(target_meta)
    return "_".join(name.strip().split())


def initialise_analysis_environment(target_name: str, explicit_id: str = None, base_output_dir: str = "outputs") -> tuple:
    """Generates an Analysis ID and sets up an isolated output folder and run-specific log file."""
    # 1. Determine the Analysis ID
    if explicit_id:
        analysis_id = _normalise_analysis_id(explicit_id.strip())
    else:
        # Auto-generate: e.g., "E_coli_run_29May_1430"
        timestamp = datetime.now().strftime("%d%b_%H%M") # e.g., 29May_1430
        analysis_id = _normalise_analysis_id(f"{target_name}_run_{timestamp}")

    # 2. Establish directory paths
    run_dir = Path(base_output_dir) / analysis_id
    run_dir.mkdir(parents=True, exist_ok=True)

    # 3. Dynamically point a dedicated log file to this run's folder
    os.makedirs("logs", exist_ok=True)
    root_logger = logging.getLogger()
    
    # Add a file handler specifically for this run's output directory
    run_log_path = run_dir / f"{analysis_id}.log"
    file_handler = logging.FileHandler(run_log_path)
    file_handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
    root_logger.addHandler(file_handler)

    logging.info(f"Initialised environment for Analysis ID: {analysis_id}")
    logging.info(f"All run-specific assets and duplicate logs will save to: {run_dir}")

    return analysis_id, run_dir