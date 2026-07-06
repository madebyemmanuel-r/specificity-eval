"""Integration tests for PipelineRunner and GUI event API.

This test suite validates that:
  - The runner correctly launches main.py as a subprocess
  - GUI event markers (RUN_DIR, REPORT_PATH, METADATA_PATH, ANALYSIS_ID, STAGE) are emitted
  - Extracted paths and metadata are accessible after execution
  - Pipeline output can be parsed for progress inference
  - The runner handles file path resolution for GUI uploads
"""

import re
import time
from pathlib import Path

import pytest

from runner import PipelineRunner, RunStatus


@pytest.fixture
def runner_config():
    """Minimal pipeline config for testing."""
    root_dir = Path(__file__).resolve().parent.parent.parent
    return {
        "oligos": str(root_dir / "MBV_Primers.csv"),
        "target": "Homo sapiens",
        "background": ["Escherichia coli"],
        "sensitivity": "Fast",
    }


@pytest.fixture
def main_py_path():
    """Path to main.py relative to the runner."""
    return Path(__file__).resolve().parent.parent.parent / "main.py"


class TestPipelineRunner:
    """Integration tests for the PipelineRunner class."""

    def test_runner_starts_and_completes(self, main_py_path, runner_config):
        """Verify runner can start main.py and wait for completion."""
        runner = PipelineRunner(main_py_path=main_py_path)
        runner.start(runner_config)
        
        # Wait for process to complete
        while runner.poll() is None:
            time.sleep(0.5)
        
        result = runner.get_result()
        assert result.status == RunStatus.SUCCESS
        assert result.return_code == 0

    def test_gui_event_markers_emitted(self, main_py_path, runner_config):
        """Verify all required GUI event markers are emitted."""
        runner = PipelineRunner(main_py_path=main_py_path)
        runner.start(runner_config)
        
        saw_run_dir = False
        saw_report_path = False
        saw_metadata_path = False
        saw_analysis_id = False
        
        while runner.poll() is None:
            for line in runner.drain_logs():
                saw_run_dir |= "GUI::RUN_DIR=" in line
                saw_report_path |= "GUI::REPORT_PATH=" in line
                saw_metadata_path |= "GUI::METADATA_PATH=" in line
                saw_analysis_id |= "GUI::ANALYSIS_ID=" in line
            time.sleep(0.3)
        
        # Drain final logs
        for line in runner.drain_logs():
            saw_run_dir |= "GUI::RUN_DIR=" in line
            saw_report_path |= "GUI::REPORT_PATH=" in line
            saw_metadata_path |= "GUI::METADATA_PATH=" in line
            saw_analysis_id |= "GUI::ANALYSIS_ID=" in line
        
        runner.get_result()
        
        assert saw_run_dir, "GUI::RUN_DIR marker not emitted"
        assert saw_report_path, "GUI::REPORT_PATH marker not emitted"
        assert saw_metadata_path, "GUI::METADATA_PATH marker not emitted"
        assert saw_analysis_id, "GUI::ANALYSIS_ID marker not emitted"

    def test_extracted_paths_exist(self, main_py_path, runner_config):
        """Verify extracted paths from GUI markers actually exist on disk."""
        runner = PipelineRunner(main_py_path=main_py_path)
        runner.start(runner_config)
        
        while runner.poll() is None:
            runner.drain_logs()
            time.sleep(0.3)
        
        runner.drain_logs()
        result = runner.get_result()
        
        assert result.run_dir is not None, "run_dir not extracted"
        assert result.run_dir.exists(), f"run_dir does not exist: {result.run_dir}"
        
        assert result.report_path is not None, "report_path not extracted"
        assert result.report_path.exists(), f"report_path does not exist: {result.report_path}"
        
        assert result.metadata_path is not None, "metadata_path not extracted"
        assert result.metadata_path.exists(), f"metadata_path does not exist: {result.metadata_path}"

    def test_analysis_id_captured(self, main_py_path, runner_config):
        """Verify analysis_id is extracted from GUI marker and stored in result."""
        runner = PipelineRunner(main_py_path=main_py_path)
        runner.start(runner_config)
        
        while runner.poll() is None:
            runner.drain_logs()
            time.sleep(0.3)
        
        runner.drain_logs()
        result = runner.get_result()
        
        assert result.analysis_id is not None, "analysis_id not extracted"
        assert len(result.analysis_id) > 0, "analysis_id is empty"
        assert result.analysis_id.replace("_", "").replace("-", "").isalnum(), \
            f"analysis_id contains unexpected characters: {result.analysis_id}"

    def test_stage_events_in_correct_order(self, main_py_path, runner_config):
        """Verify stage events are emitted in the expected sequence."""
        runner = PipelineRunner(main_py_path=main_py_path)
        runner.start(runner_config)
        
        stage_events = []
        
        while runner.poll() is None:
            for line in runner.drain_logs():
                if "GUI::STAGE=" in line:
                    match = re.search(r"GUI::STAGE=(\S+)", line)
                    if match:
                        stage_events.append(match.group(1))
            time.sleep(0.3)
        
        for line in runner.drain_logs():
            if "GUI::STAGE=" in line:
                match = re.search(r"GUI::STAGE=(\S+)", line)
                if match:
                    stage_events.append(match.group(1))
        
        runner.get_result()
        
        # Verify expected start sequence
        expected_start = ["INITIALISING_PIPELINE", "VERIFYING_NCBI_CONNECTIVITY", "VALIDATING_OLIGOS"]
        assert len(stage_events) >= len(expected_start), \
            f"Expected at least {len(expected_start)} stage events, got {len(stage_events)}"
        
        assert stage_events[:len(expected_start)] == expected_start, \
            f"Stage sequence does not begin correctly. Got: {stage_events[:3]}"
        
        # Verify pipeline completes with PIPELINE_COMPLETE
        assert "PIPELINE_COMPLETE" in stage_events, \
            "PIPELINE_COMPLETE event not emitted"

    def test_metadata_valid_json(self, main_py_path, runner_config):
        """Verify extracted metadata is valid JSON with expected keys."""
        runner = PipelineRunner(main_py_path=main_py_path)
        runner.start(runner_config)
        
        while runner.poll() is None:
            runner.drain_logs()
            time.sleep(0.3)
        
        runner.drain_logs()
        result = runner.get_result()
        
        assert result.metadata is not None, "metadata not loaded"
        assert isinstance(result.metadata, dict), "metadata is not a dict"
        
        # Check for expected top-level keys
        expected_keys = {"pipeline_execution", "assay_parameters", "target_organism"}
        missing_keys = expected_keys - set(result.metadata.keys())
        assert not missing_keys, f"metadata missing keys: {missing_keys}"

    def test_file_path_resolution(self, main_py_path):
        """Verify the runner resolves GUI upload paths correctly."""
        runner = PipelineRunner(main_py_path=main_py_path)
        
        root_dir = main_py_path.parent
        existing_file = root_dir / "MBV_Primers.csv"
        
        assert existing_file.exists(), f"Test file does not exist: {existing_file}"
        
        # Test resolution of a GUI upload path
        resolved = runner._resolve_file_path("gui_uploads/MBV_Primers.csv")
        resolved_path = Path(resolved)
        
        assert resolved_path.exists(), \
            f"Path resolver did not find file: {resolved}"

    def test_progress_inference(self, main_py_path, runner_config):
        """Verify progress is inferred from stage markers."""
        runner = PipelineRunner(main_py_path=main_py_path)
        runner.start(runner_config)
        
        progresses_seen = []
        
        while runner.poll() is None:
            runner.drain_logs()
            progresses_seen.append(runner.progress)
            time.sleep(0.5)
        
        runner.drain_logs()
        runner.get_result()
        
        # Progress should increase over time (or stay at final value)
        assert len(progresses_seen) > 0, "No progress samples collected"
        assert max(progresses_seen) > 0.0, "Progress never increased from 0"
        assert max(progresses_seen) <= 1.0, "Progress exceeded 1.0"
