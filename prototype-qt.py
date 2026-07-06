import sys
from pathlib import Path

from PySide6.QtCore import QTimer, QUrl
from PySide6.QtGui import QDesktopServices
from PySide6.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QHBoxLayout, QFormLayout,
    QLabel, QLineEdit, QComboBox, QSpinBox, QDoubleSpinBox, QPushButton,
    QProgressBar, QPlainTextEdit, QFileDialog, QListWidget, QGroupBox,
    QMessageBox, QTextBrowser, QFrame, QToolBox, QWidget as QtWidget,
)

from runner import PipelineRunner, RunStatus

MAIN_PY_PATH = Path("main.py")  # adjust if main.py lives elsewhere relative to this app


class PipelineWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("RPA Pipeline Orchestrator")
        self.resize(980, 820)

        self.runner: PipelineRunner | None = None
        self.uploaded_csv_path: Path | None = None

        self.poll_timer = QTimer(self)
        self.poll_timer.setInterval(200)
        self.poll_timer.timeout.connect(self.poll_runner)

        root = QVBoxLayout(self)
        root.setContentsMargins(18, 18, 18, 18)
        root.setSpacing(16)

        self.setStyleSheet('''
            QWidget { background-color: #f8fafc; color: #111827; font-family: 'Segoe UI', Arial, sans-serif; }
            QGroupBox { border: 1px solid #e5e7eb; border-radius: 10px; margin-top: 12px; padding: 12px; background: #ffffff; }
            QGroupBox:title { subcontrol-origin: margin; left: 14px; padding: 0 8px 0 8px; color: #111827; font-weight: 600; }
            QLabel { color: #111827; }
            QLineEdit, QComboBox, QSpinBox, QDoubleSpinBox, QPlainTextEdit, QTextBrowser, QListWidget { border: 1px solid #d1d5db; border-radius: 8px; padding: 6px; background: #ffffff; }
            QLineEdit:focus, QComboBox:focus, QSpinBox:focus, QDoubleSpinBox:focus { border-color: #2563eb; }
            QPushButton { background-color: #2563eb; color: white; border-radius: 10px; padding: 10px 16px; font-weight: 600; }
            QPushButton:hover { background-color: #1e40af; }
            QPushButton:disabled { background-color: #9ca3af; color: #f3f4f6; }
            QProgressBar { border: 1px solid #d1d5db; border-radius: 10px; text-align: center; }
            QProgressBar::chunk { background-color: #2563eb; border-radius: 10px; }
            QPlainTextEdit { background-color: #111827; color: #f8fafc; border: 1px solid #374151; }
            QTextBrowser { background-color: #ffffff; }
            QListWidget { background-color: #ffffff; }
        ''')

        title = QLabel("RPA Specificity Pipeline")
        title.setStyleSheet("font-size: 24px; font-weight: 800; color: #111827;")
        subtitle = QLabel("Launch the pipeline, monitor progress, and inspect results from one place.")
        subtitle.setStyleSheet("font-size: 12px; color: #4b5563; margin-bottom: 10px;")
        root.addWidget(title)
        root.addWidget(subtitle)

        separator = QFrame()
        separator.setFrameShape(QFrame.HLine)
        separator.setFrameShadow(QFrame.Sunken)
        separator.setStyleSheet("color: #d1d5db;")
        root.addWidget(separator)

        # --- Config group ---
        config_group = QGroupBox("Run Configuration")
        config_layout = QVBoxLayout()
        toolbox = QToolBox()

        # Basic settings
        basic_widget = QtWidget()
        basic_form = QFormLayout()

        upload_row = QHBoxLayout()
        self.upload_button = QPushButton("Choose Oligo CSV...")
        self.upload_button.clicked.connect(self.choose_csv)
        self.upload_status = QLabel("No file selected.")
        upload_row.addWidget(self.upload_button)
        upload_row.addWidget(self.upload_status)
        basic_form.addRow(upload_row)

        self.target_input = QLineEdit("Homo sapiens")
        basic_form.addRow("Target organism:", self.target_input)

        background_row = QHBoxLayout()
        self.background_input = QLineEdit()
        self.background_input.setPlaceholderText("Enter organism and click Add")
        self.background_add_button = QPushButton("Add")
        self.background_add_button.clicked.connect(self.add_background)
        background_row.addWidget(self.background_input)
        background_row.addWidget(self.background_add_button)
        basic_form.addRow("Background organisms:", background_row)

        self.background_list = QListWidget()
        self.background_list.addItem("Escherichia coli")
        self.background_list.setMaximumHeight(120)
        basic_form.addRow(self.background_list)

        self.background_remove_button = QPushButton("Remove selected")
        self.background_remove_button.clicked.connect(self.remove_background)
        basic_form.addRow(self.background_remove_button)

        self.sensitivity_select = QComboBox()
        self.sensitivity_select.addItems(["Fast", "Standard", "Advanced", "Very Advanced"])
        self.sensitivity_select.setCurrentText("Standard")
        basic_form.addRow("Sensitivity:", self.sensitivity_select)

        basic_widget.setLayout(basic_form)
        toolbox.addItem(basic_widget, "Basic")

        # Advanced settings
        advanced_widget = QtWidget()
        advanced_form = QFormLayout()

        self.evalue_input = QDoubleSpinBox()
        self.evalue_input.setRange(0, 1_000_000)
        self.evalue_input.setValue(1000.0)
        self.evalue_input.setDecimals(2)
        advanced_form.addRow("E-value:", self.evalue_input)

        self.max_amplicon_input = QSpinBox()
        self.max_amplicon_input.setRange(1, 100_000)
        self.max_amplicon_input.setValue(500)
        advanced_form.addRow("Max amplicon size (bp):", self.max_amplicon_input)

        self.threads_input = QSpinBox()
        self.threads_input.setRange(1, 128)
        self.threads_input.setValue(4)
        advanced_form.addRow("Threads:", self.threads_input)

        self.analysis_id_input = QLineEdit()
        self.analysis_id_input.setPlaceholderText("optional")
        advanced_form.addRow("Analysis ID:", self.analysis_id_input)

        self.db_output_input = QLineEdit("outputs/blast_dbs/")
        advanced_form.addRow("BLAST DB output dir:", self.db_output_input)

        self.api_key_input = QLineEdit()
        self.api_key_input.setEchoMode(QLineEdit.Password)
        self.api_key_input.setPlaceholderText("optional")
        advanced_form.addRow("NCBI API key:", self.api_key_input)

        advanced_widget.setLayout(advanced_form)
        toolbox.addItem(advanced_widget, "Advanced")

        config_layout.addWidget(toolbox)
        config_group.setLayout(config_layout)
        root.addWidget(config_group)

        # --- Controls ---
        controls_row = QHBoxLayout()
        self.run_button = QPushButton("Start Pipeline")
        self.run_button.clicked.connect(self.start_pipeline)
        self.cancel_button = QPushButton("Cancel")
        self.cancel_button.clicked.connect(self.cancel_pipeline)
        self.cancel_button.setEnabled(False)
        controls_row.addWidget(self.run_button)
        controls_row.addWidget(self.cancel_button)
        controls_row.addStretch()
        root.addLayout(controls_row)

        # --- Progress ---
        progress_row = QHBoxLayout()
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setTextVisible(True)
        self.stage_label = QLabel("Idle")
        self.stage_label.setStyleSheet("font-weight: 600; color: #2563eb;")
        progress_row.addWidget(self.progress_bar, 3)
        progress_row.addWidget(self.stage_label, 1)
        root.addLayout(progress_row)

        # --- Console ---
        root.addWidget(QLabel("Console Output:"))
        self.console = QPlainTextEdit()
        self.console.setReadOnly(True)
        self.console.setStyleSheet(
            "background-color: #111827; color: #d1fae5; font-family: 'Courier New', monospace; font-size: 12px; border-radius: 10px;"
        )
        self.console.setMinimumHeight(260)
        root.addWidget(self.console)

        # --- Results ---
        self.results_group = QGroupBox("Results")
        self.results_group.setVisible(False)
        results_layout = QVBoxLayout()
        results_layout.setSpacing(10)
        results_layout.addWidget(QLabel("Report:"))
        self.report_viewer = QTextBrowser()
        self.report_viewer.setStyleSheet("border: 1px solid #d1d5db; border-radius: 10px; background: #f9fafb;")
        results_layout.addWidget(self.report_viewer)

        summary_row = QHBoxLayout()
        self.summary_label = QLabel("Runtime: N/A  |  Files: N/A  |  Status: N/A")
        self.summary_label.setStyleSheet("font-size: 12px; color: #374151;")
        summary_row.addWidget(self.summary_label)
        summary_row.addStretch()
        results_layout.addLayout(summary_row)

        actions_row = QHBoxLayout()
        self.open_report_button = QPushButton("Open Report")
        self.open_report_button.setEnabled(False)
        self.open_report_button.clicked.connect(self.open_report)
        self.open_folder_button = QPushButton("Open Output Folder")
        self.open_folder_button.setEnabled(False)
        self.open_folder_button.clicked.connect(self.open_output_folder)
        actions_row.addWidget(self.open_report_button)
        actions_row.addWidget(self.open_folder_button)
        actions_row.addStretch()
        results_layout.addLayout(actions_row)

        results_layout.addWidget(QLabel("Output files:"))
        self.files_list = QListWidget()
        self.files_list.setMaximumHeight(140)
        results_layout.addWidget(self.files_list)
        self.results_group.setLayout(results_layout)
        root.addWidget(self.results_group)

    # --- Handlers ---

    def choose_csv(self):
        path, _ = QFileDialog.getOpenFileName(self, "Select Oligo CSV", "", "CSV Files (*.csv)")
        if path:
            self.uploaded_csv_path = Path(path)
            self.upload_status.setText(f"Selected: {self.uploaded_csv_path.name}")
            self.upload_status.setStyleSheet("color: #059669; font-weight: 600;")

    def add_background(self):
        value = self.background_input.text().strip()
        if not value:
            return
        existing = [self.background_list.item(i).text() for i in range(self.background_list.count())]
        if value in existing:
            self.background_input.clear()
            return
        self.background_list.addItem(value)
        self.background_input.clear()

    def remove_background(self):
        for item in self.background_list.selectedItems():
            self.background_list.takeItem(self.background_list.row(item))

    def open_report(self):
        if getattr(self, 'report_path', None) and self.report_path.exists():
            QDesktopServices.openUrl(QUrl.fromLocalFile(str(self.report_path)))

    def open_output_folder(self):
        if getattr(self, 'run_dir', None) and self.run_dir.exists():
            QDesktopServices.openUrl(QUrl.fromLocalFile(str(self.run_dir)))

    def start_pipeline(self):
        if self.uploaded_csv_path is None:
            QMessageBox.warning(self, "Missing input", "Please select an oligo CSV first.")
            return
        if not self.target_input.text().strip():
            QMessageBox.warning(self, "Missing input", "Target organism is required.")
            return

        background_list = [self.background_list.item(i).text() for i in range(self.background_list.count())]
        if not background_list:
            QMessageBox.warning(self, "Missing input", "At least one background organism is required.")
            return

        config = {
            "oligos": str(self.uploaded_csv_path),
            "target": self.target_input.text().strip(),
            "background": background_list,
            "analysis_id": self.analysis_id_input.text().strip() or None,
            "db_output": self.db_output_input.text().strip() or "outputs/blast_dbs/",
            "threads": self.threads_input.value(),
            "sensitivity": self.sensitivity_select.currentText(),
            "evalue": self.evalue_input.value(),
            "max_amplicon_size": self.max_amplicon_input.value(),
            "api_key": self.api_key_input.text().strip() or None,
        }

        self.results_group.setVisible(False)
        self.console.clear()
        self.console.appendPlainText("--- Initialising pipeline execution ---")
        self.run_button.setEnabled(False)
        self.cancel_button.setEnabled(True)
        self.progress_bar.setValue(0)
        self.stage_label.setText("Starting...")

        self.runner = PipelineRunner(main_py_path=MAIN_PY_PATH)
        self.runner.start(config)
        self.poll_timer.start()

    def cancel_pipeline(self):
        if self.runner:
            self.runner.cancel()
            self.console.appendPlainText("\n[SYSTEM] Cancellation requested...")
            self.cancel_button.setEnabled(False)

    def poll_runner(self):
        if not self.runner:
            return

        for line in self.runner.drain_logs():
            self.console.appendPlainText(line)

        self.progress_bar.setValue(int(self.runner.progress * 100))
        self.stage_label.setText(self.runner.current_stage_label)

        # poll() is the authoritative completion signal — progress is only
        # an estimate and may not hit 100% exactly before the process exits.
        if self.runner.poll() is not None:
            self.poll_timer.stop()
            self.run_button.setEnabled(True)
            self.cancel_button.setEnabled(False)

            result = self.runner.get_result()

            if result.status == RunStatus.SUCCESS:
                self.console.appendPlainText(f"\n[SYSTEM] Pipeline finished. Report: {result.report_path}")
                QMessageBox.information(self, "Done", "Pipeline execution finished!")
                self.show_results(result)
            elif result.status == RunStatus.CANCELLED:
                self.console.appendPlainText("\n[SYSTEM] Pipeline was cancelled.")
            else:
                error_detail = f"\nError: {result.error_message}" if result.error_message else ""
                self.console.appendPlainText(
                    f"\n[SYSTEM] Pipeline failed (exit code {result.return_code}).{error_detail}"
                )
                QMessageBox.critical(self, "Failed", "Pipeline failed — see console output.")

    def show_results(self, result):
        self.files_list.clear()
        self.run_dir = result.run_dir
        self.report_path = result.report_path

        if result.report_path and result.report_path.exists():
            self.report_viewer.setMarkdown(result.report_path.read_text())
            self.open_report_button.setEnabled(True)
        else:
            self.report_viewer.setMarkdown("_No report file found at the expected path._")
            self.open_report_button.setEnabled(False)

        file_count = 0
        if result.run_dir and result.run_dir.exists():
            for path in sorted(result.run_dir.iterdir()):
                if path.is_file():
                    size_kb = path.stat().st_size / 1024
                    self.files_list.addItem(f"{path.name}  ({size_kb:.1f} KB)")
                    file_count += 1
            self.open_folder_button.setEnabled(True)
        else:
            self.open_folder_button.setEnabled(False)

        self.summary_label.setText(
            f"Runtime: {result.runtime or 'N/A'}  |  Files: {file_count}  |  Status: {result.status.name}"
        )
        self.results_group.setVisible(True)


def main():
    app = QApplication(sys.argv)
    window = PipelineWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()