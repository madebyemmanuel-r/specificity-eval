"""
Wraps main.py as a subprocess so the GUI can start a pipeline run,
stream live logs, estimate progress from log markers, cancel cleanly,
and read the final result.
"""

import json
import os
import queue
import re
import signal
import subprocess
import sys
import threading
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Optional


class RunStatus(Enum):
    RUNNING = "running"
    SUCCESS = "success"
    FAILED = "failed"
    CANCELLED = "cancelled"


@dataclass
class RunResult:
    status: RunStatus
    return_code: Optional[int]
    run_dir: Optional[Path]
    report_path: Optional[Path]
    metadata_path: Optional[Path]
    analysis_id: Optional[str] = None
    error_message: Optional[str] = None
    metadata: Optional[dict] = None


class PipelineRunner:
    """
    Launches main.py as a subprocess, streams live logs, estimates progress
    from known log markers, and supports clean cancellation.
    """

    # Substring markers used to infer progress from main.py log output.
    STAGE_MARKERS = [
        ("Initialising primer specificity analysis tool", "Starting"),
        ("NCBI connectivity verified", "Verifying NCBI connection"),
        ("Parsing and validating oligo profiles", "Validating input CSV"),
        ("Resolving target organism's identity", "Resolving target organism"),
        ("Initialising analysis environment", "Initialising analysis environment"),
        ("Fetching target genome", "Fetching target genome"),
        ("Target BLAST DB ready", "Building target BLAST DB"),
        ("Resolving background genomes", "Resolving background genomes"),
        ("Background BLAST DB ready", "Building background BLAST DB"),
        ("Running background BLASTn", "Screening off-target (background BLAST)"),
        ("Running target BLASTn", "Confirming on-target amplification"),
        ("Initiating downstream specificity evaluation", "Running specificity + thermodynamic analysis"),
        ("Writing metadata JSON", "Writing run metadata"),
        ("Specificity pipeline complete", "Finalising report"),
    ]

    RUN_DIR_PATTERN = re.compile(r"GUI::RUN_DIR=(.+)$")
    REPORT_PATH_PATTERN = re.compile(r"GUI::REPORT_PATH=(.+)$")
    METADATA_PATH_PATTERN = re.compile(r"GUI::METADATA_PATH=(.+)$")
    ANALYSIS_ID_PATTERN = re.compile(r"GUI::ANALYSIS_ID=(.+)$")
    STAGE_EVENT_PATTERN = re.compile(r"GUI::STAGE=(.+)$")
    ERROR_PATTERN = re.compile(r"ERROR::(.+)$")
    STAGE_EVENT_ORDER = [
        "INITIALISING_PIPELINE",
        "VERIFYING_NCBI_CONNECTIVITY",
        "VALIDATING_OLIGOS",
        "RESOLVING_TARGET_ORGANISM",
        "INITIALISING_ANALYSIS_ENVIRONMENT",
        "FETCHING_TARGET_GENOME",
        "TARGET_BLAST_DB_READY",
        "RESOLVING_BACKGROUND_GENOMES",
        "BACKGROUND_BLAST_DB_READY",
        "RUNNING_BACKGROUND_BLAST",
        "RUNNING_TARGET_BLAST",
        "INITIATING_SPECIFICITY_PIPELINE",
        "WRITING_METADATA",
        "PIPELINE_COMPLETE",
    ]

    def __init__(self, main_py_path: Path, python_executable: Optional[str] = None):
        self.main_py_path = Path(main_py_path)
        self.python_executable = python_executable or sys.executable
        self.log_queue: "queue.Queue[str]" = queue.Queue()
        self.process: Optional[subprocess.Popen] = None
        self._reader_thread: Optional[threading.Thread] = None
        self._run_dir: Optional[Path] = None
        self._report_path: Optional[Path] = None
        self._metadata_path: Optional[Path] = None
        self._analysis_id: Optional[str] = None
        self._current_stage_index: int = -1
        self._cancelled = False

    @property
    def total_stages(self) -> int:
        return len(self.STAGE_MARKERS)

    @property
    def progress(self) -> float:
        """Rough 0-1 estimate based on the furthest stage marker seen so far."""
        if self._current_stage_index < 0:
            return 0.0
        return (self._current_stage_index + 1) / self.total_stages

    @property
    def current_stage_label(self) -> str:
        if self._current_stage_index < 0:
            return "Starting..."
        return self.STAGE_MARKERS[self._current_stage_index][1]

    def start(self, config: dict) -> None:
        """
        config maps CLI flag names (without leading --) to values, e.g.:
            {"oligos": "primers.csv", "target": "Homo sapiens",
             "background": ["Escherichia coli", "Staphylococcus aureus"]}
        """
        # Resolve common GUI upload and fallback locations for file paths.
        # If the provided path doesn't exist, try: gui_uploads/, repo root, etc.
        config = dict(config)  # shallow copy to avoid mutating caller's config
        if "oligos" in config and config["oligos"]:
            config["oligos"] = self._resolve_file_path(config["oligos"])

        argv = self._build_argv(config)
        cmd = [self.python_executable, "-u", str(self.main_py_path), *argv]

        self.process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            cwd=self.main_py_path.parent,
            start_new_session=True,  # own process group -> enables clean group kill
        )
        self._cancelled = False
        self._run_dir = None
        self._report_path = None
        self._metadata_path = None
        self._error_message = None
        self._current_stage_index = -1

        self._reader_thread = threading.Thread(target=self._read_stdout, daemon=True)
        self._reader_thread.start()

    def _build_argv(self, config: dict) -> list[str]:
        argv: list[str] = []
        for key, value in config.items():
            # Convert pythonic underscore keys to CLI-style hyphen flags
            flag_name = str(key).replace("_", "-")
            flag = f"--{flag_name}"
            if isinstance(value, bool):
                if value:
                    argv.append(flag)
            elif isinstance(value, (list, tuple)):
                # list/tuple values become repeated flag args for nargs="+"
                argv.append(flag)
                argv.extend(str(v) for v in value)
            elif value is not None:
                argv.append(flag)
                argv.append(str(value))
        return argv

    def _resolve_file_path(self, provided_path: str) -> str:
        """Resolve file path with fallback locations.

        Try the provided path first, then check common upload/repo locations.
        Returns the first existing path, or the provided path if none are found.
        """
        candidates = [
            provided_path,
            str(Path("gui_uploads") / Path(provided_path).name),
            str(Path(self.main_py_path.parent) / Path(provided_path).name),
            str(Path(self.main_py_path.parent) / provided_path),
        ]

        for candidate in candidates:
            p = Path(candidate)
            if p.exists():
                return candidate

        return provided_path  # return original; let main.py raise the error

    def _read_stdout(self) -> None:
        assert self.process is not None and self.process.stdout is not None
        for line in self.process.stdout:
            line = line.rstrip("\n")
            self.log_queue.put(line)

            run_dir_match = self.RUN_DIR_PATTERN.search(line)
            if run_dir_match:
                self._run_dir = Path(run_dir_match.group(1).strip())

            report_match = self.REPORT_PATH_PATTERN.search(line)
            if report_match:
                self._report_path = Path(report_match.group(1).strip())

            metadata_match = self.METADATA_PATH_PATTERN.search(line)
            if metadata_match:
                self._metadata_path = Path(metadata_match.group(1).strip())

            analysis_id_match = self.ANALYSIS_ID_PATTERN.search(line)
            if analysis_id_match:
                self._analysis_id = analysis_id_match.group(1).strip()

            error_match = self.ERROR_PATTERN.search(line)
            if error_match:
                self._error_message = error_match.group(1).strip()

            stage_event_match = self.STAGE_EVENT_PATTERN.search(line)
            if stage_event_match:
                stage_name = stage_event_match.group(1).strip()
                if stage_name in self.STAGE_EVENT_ORDER:
                    idx = self.STAGE_EVENT_ORDER.index(stage_name)
                    if idx > self._current_stage_index:
                        self._current_stage_index = idx

            for i, (marker, _label) in enumerate(self.STAGE_MARKERS):
                if marker in line and i > self._current_stage_index:
                    self._current_stage_index = i
                    break

        self.process.stdout.close()

    def poll(self) -> Optional[int]:
        """Non-blocking check. Returns the return code once finished, else None."""
        if self.process is None:
            return None
        return self.process.poll()

    def drain_logs(self) -> list[str]:
        """Call from a UI timer to pull all lines queued since the last drain."""
        lines = []
        while True:
            try:
                lines.append(self.log_queue.get_nowait())
            except queue.Empty:
                break
        return lines

    def cancel(self) -> None:
        """Terminate the pipeline process group, including BLAST children."""
        if self.process and self.process.poll() is None:
            self._cancelled = True
            try:
                os.killpg(os.getpgid(self.process.pid), signal.SIGTERM)
            except ProcessLookupError:
                pass  # already exited

    def get_result(self) -> RunResult:
        """Return the final run result after the subprocess has finished."""
        if self.process is None:
            raise RuntimeError("Pipeline has not been started.")

        return_code = self.process.wait()
        if self._reader_thread is not None:
            self._reader_thread.join(timeout=5)

        if self._cancelled:
            status = RunStatus.CANCELLED
        elif return_code == 0:
            status = RunStatus.SUCCESS
        else:
            status = RunStatus.FAILED

        metadata = None
        metadata_path = self._metadata_path
        if metadata_path is None and self._run_dir and self._run_dir.exists():
            candidate = self._run_dir / "run_metadata.json"
            if candidate.exists():
                metadata_path = candidate

        if metadata_path and metadata_path.exists():
            try:
                metadata = json.loads(metadata_path.read_text())
            except json.JSONDecodeError:
                metadata = None

        return RunResult(
            status=status,
            return_code=return_code,
            run_dir=self._run_dir,
            report_path=self._report_path,
            metadata_path=metadata_path,
            analysis_id=self._analysis_id,
            error_message=self._error_message,
            metadata=metadata,
        )