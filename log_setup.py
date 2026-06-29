import logging
import os
import sys
from datetime import datetime

logger = logging.getLogger(__name__)

def setup_logging(output_dir: str = None):
    """
    Configure logging for a single isolated run.
    - Clears any handlers attached by previous calls or early module imports
    - Writes to logs/pipeline.log (overwrite mode, not append)
    - Optionally mirrors the log into output_dir if supplied
    - Always streams to stdout
    Handles two distinct phases of pipeline execution safely:
    Phase 1: Initial call sets up root 'logs/pipeline.log' and stdout.
    Phase 2: Later call attaches a dynamic run-specific output directory and 
             flushes historical initialization steps to prevent partial logging.
    """
    root = logging.getLogger()
    central_log_path = os.path.abspath("logs/pipeline.log")
    
    # =========================================================================
    # PHASE 1: Canonical Central Log File & Terminal Handlers
    # =========================================================================
    if not root.handlers:
        root.setLevel(logging.INFO)
        fmt = logging.Formatter('%(asctime)s - %(levelname)s - %(name)s - %(message)s')
        
        # 1. Primary central workspace ledger
        os.makedirs("logs", exist_ok=True)
        fh = logging.FileHandler(central_log_path, mode="w", encoding="utf-8")
        fh.setFormatter(fmt)
        root.addHandler(fh)
        
        # 2. Synchronous terminal readout
        sh = logging.StreamHandler(sys.stdout)
        sh.setFormatter(fmt)
        root.addHandler(sh)

    # =========================================================================
    # PHASE 2: Dynamic Run Directory Injection & Cache Catch-Up
    # =========================================================================
    if output_dir:
        fmt = root.handlers[0].formatter  # Re-use the existing pipeline layout
        os.makedirs(output_dir, exist_ok=True)
        run_log_path = os.path.abspath(os.path.join(output_dir, "pipeline.log"))
        
        # Guard clause: Prevent duplicate handler registration if called iteratively
        already_attached = any(
            isinstance(h, logging.FileHandler) and h.baseFilename == run_log_path
            for h in root.handlers
        )
        
        if not already_attached:
            # 💡 THE CATCH-UP FIX: If historical records exist, flush them to the run log
            if os.path.exists(central_log_path):
                try:
                    with open(central_log_path, "r", encoding="utf-8") as central_file:
                        history = central_file.read()
                    with open(run_log_path, "w", encoding="utf-8") as run_file:
                        run_file.write(history)
                except (IOError, OSError) as e:
                    root.warning(f"Unable to seed historical log entries to run workspace: {e}")

            # Open the file handler in append ('a') mode to seamlessly pick up from here
            run_fh = logging.FileHandler(run_log_path, mode="a", encoding="utf-8")
            run_fh.setFormatter(fmt)
            root.addHandler(run_fh)
            
            logging.getLogger(__name__).info(
                f"Dynamic run-specific logging handler attached cleanly to: {output_dir}"
            )