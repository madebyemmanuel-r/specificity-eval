import os
import logging
from datetime import datetime
from pathlib import Path

logger = logging.getLogger(__name__)

def initialise_analysis_environment(target_name: str, explicit_id: str = None, base_output_dir: str = "outputs") -> tuple:
    """Generates an Analysis ID and sets up an isolated output folder and run-specific log file."""
    # 1. Determine the Analysis ID
    if explicit_id:
        analysis_id = explicit_id.strip().replace(" ", "_")
    else:
        # Auto-generate: e.g., "E_coli_run_29May_1430"
        timestamp = datetime.now().strftime("%d%b_%H%M") # e.g., 29May_1430
        analysis_id = f"{target_name}_run_{timestamp}"

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