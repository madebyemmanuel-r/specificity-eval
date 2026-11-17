"""
smoke_test_runner.py

Run this once, standalone, before wiring PipelineRunner into either GUI
prototype. It exercises the exact same code path a GUI would use — launch,
live log drain, progress inference, result retrieval — so any problems
(wrong flag names, buffering issues, missing RUN_DIR/REPORT_PATH markers)
surface here rather than inside a half-built UI.

Edit the CONFIG dict below to point at real inputs, then:

    python smoke_test_runner.py

Expect: log lines printing live as the pipeline runs, a progress percentage
increasing as stage markers are matched, and a final summary block.
"""

import sys
import time
from pathlib import Path

ROOT_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT_DIR))

from runner import PipelineRunner

CONFIG = {
    "oligos": str(ROOT_DIR / "MBV_Primers.csv"),  # sample CSV already in the repo
    "target": "Homo sapiens",                     # point at a real target
    "background": ["Escherichia coli"],           # point at real background(s)
    "sensitivity": "Fast",                        # fast preset for a quick smoke test
}

MAIN_PY_PATH = ROOT_DIR / "main.py"


def main():
    runner = PipelineRunner(main_py_path=MAIN_PY_PATH)
    runner.start(CONFIG)

    print(f"Launched: {' '.join([runner.python_executable, '-u', str(MAIN_PY_PATH)])}")
    print("Waiting for output...\n")

    saw_run_dir_marker = False
    saw_report_path_marker = False
    saw_metadata_path_marker = False
    saw_analysis_id_marker = False
    saw_error_marker = False
    stage_events_seen = []

    while runner.poll() is None:
        for line in runner.drain_logs():
            print(f"  [log] {line}")
            saw_run_dir_marker |= "GUI::RUN_DIR=" in line
            saw_report_path_marker |= "GUI::REPORT_PATH=" in line
            saw_metadata_path_marker |= "GUI::METADATA_PATH=" in line
            saw_analysis_id_marker |= "GUI::ANALYSIS_ID=" in line
            saw_error_marker |= line.startswith("ERROR::")
            
            # Track stage events in order
            if "GUI::STAGE=" in line:
                import re
                match = re.search(r"GUI::STAGE=(\S+)", line)
                if match:
                    stage_events_seen.append(match.group(1))
        
        print(f"  [progress] {runner.progress:.0%} — {runner.current_stage_label}", end="\r")
        time.sleep(0.3)

    # Drain anything left in the queue after the process exited
    for line in runner.drain_logs():
        print(f"  [log] {line}")
        saw_run_dir_marker |= "GUI::RUN_DIR=" in line
        saw_report_path_marker |= "GUI::REPORT_PATH=" in line
        saw_metadata_path_marker |= "GUI::METADATA_PATH=" in line
        saw_analysis_id_marker |= "GUI::ANALYSIS_ID=" in line
        saw_error_marker |= line.startswith("ERROR::")
        
        if "GUI::STAGE=" in line:
            import re
            match = re.search(r"GUI::STAGE=(\S+)", line)
            if match:
                stage_events_seen.append(match.group(1))

    result = runner.get_result()
    print("\n\n--- RESULT ---")
    print(f"status:        {result.status.value}")
    print(f"return_code:   {result.return_code}")
    print(f"run_dir:       {result.run_dir}")
    print(f"report_path:   {result.report_path}")
    print(f"metadata_path: {result.metadata_path}")
    print(f"analysis_id:   {result.analysis_id}")
    if result.metadata:
        print(f"metadata keys: {list(result.metadata.keys())}")

    print("\n--- MARKER VALIDATION ---")
    print(f"RUN_DIR marker seen:      {saw_run_dir_marker}")
    print(f"REPORT_PATH marker seen:  {saw_report_path_marker}")
    print(f"METADATA_PATH marker seen:{saw_metadata_path_marker}")
    print(f"ANALYSIS_ID marker seen:  {saw_analysis_id_marker}")
    print(f"ERROR marker seen:        {saw_error_marker}")

    print("\n--- PATH EXISTENCE CHECKS ---")
    paths_ok = True
    if result.run_dir and result.run_dir.exists():
        print(f"✓ run_dir exists:      {result.run_dir}")
    else:
        print(f"✗ run_dir missing:     {result.run_dir}")
        paths_ok = False
    
    if result.report_path and result.report_path.exists():
        print(f"✓ report_path exists:  {result.report_path}")
    else:
        print(f"✗ report_path missing: {result.report_path}")
        paths_ok = False
    
    if result.metadata_path and result.metadata_path.exists():
        print(f"✓ metadata_path exists:{result.metadata_path}")
    else:
        print(f"✗ metadata_path missing:{result.metadata_path}")
        paths_ok = False

    print("\n--- STAGE EVENT ORDER ---")
    print(f"Stage events captured: {len(stage_events_seen)}")
    for i, stage in enumerate(stage_events_seen[:5], 1):
        print(f"  {i}. {stage}")
    if len(stage_events_seen) > 5:
        print(f"  ... and {len(stage_events_seen) - 5} more")
    
    # Verify the expected sequence starts correctly
    expected_start = ["INITIALISING_PIPELINE", "VERIFYING_NCBI_CONNECTIVITY", "VALIDATING_OLIGOS"]
    stage_order_ok = stage_events_seen[:3] == expected_start if len(stage_events_seen) >= 3 else False
    if stage_order_ok:
        print("✓ Stage order begins correctly")
    else:
        print("✗ Stage order mismatch")

    print("\n--- SMOKE VALIDATION ---")
    all_markers_ok = saw_run_dir_marker and saw_report_path_marker and saw_metadata_path_marker and saw_analysis_id_marker
    if result.status.value == "success":
        if all_markers_ok and paths_ok and runner.progress > 0:
            print("PASS: All checks passed!")
        else:
            print("FAIL: Success but missing markers or paths")
            if not all_markers_ok:
                print("  - Missing marker(s)")
            if not paths_ok:
                print("  - Extracted paths don't exist")
    else:
        print("FAIL: Pipeline did not succeed")



if __name__ == "__main__":
    main()