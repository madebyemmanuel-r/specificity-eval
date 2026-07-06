from pathlib import Path

from nicegui import ui, events

from runner import PipelineRunner, RunStatus

MAIN_PY_PATH = Path("main.py")  # adjust if main.py lives elsewhere relative to this app
UPLOAD_DIR = Path("gui_uploads")
UPLOAD_DIR.mkdir(exist_ok=True)


class PipelineAppScaffold:
    def __init__(self):
        self.runner: PipelineRunner | None = None
        self.polling_timer = None
        self.uploaded_csv_path: Path | None = None

        ui.add_head_html('''
            <style>
                body { background: #f8fafc; }
                .nicegui-footer, .nicegui-topbar { display: none !important; }
            </style>
        ''')

        with ui.column().classes('items-center'):
            with ui.card().classes('w-full max-w-5xl m-4 p-6 shadow-lg border border-slate-200 items-center'):
                with ui.column().classes('gap-2 items-center'):
                    ui.label('RPA Pipeline Orchestrator').classes('text-4xl font-extrabold text-slate-900 text-center w-full')
                    ui.label('Launch the primer specificity pipeline, monitor progress, and inspect results from one place.').classes('text-sm text-slate-600 text-center max-w-3xl w-full')

        ui.separator().classes('mx-auto w-11/12')

        # --- Config form ---
        with ui.card().classes('w-full max-w-5xl mx-auto m-4 p-6 gap-4 shadow-lg border border-slate-200'):
            ui.label('Run Configuration').classes('text-xl font-semibold text-slate-900')

            ui.upload(on_upload=self.handle_upload, auto_upload=True, label='Oligo CSV').classes('w-full')
            self.upload_status = ui.label('No file uploaded yet.').classes('text-sm text-slate-500')

            self.target_input = ui.input('Target organism', value='Homo sapiens').classes('w-full')
            self.background_input = ui.input(
                'Background organisms (comma-separated)', value='Escherichia coli'
            ).classes('w-full')

            with ui.row().classes('w-full gap-4'):
                self.sensitivity_select = ui.select(
                    ["Fast", "Standard", "Advanced", "Very Advanced"],
                    value="Standard",
                    label="Sensitivity",
                )
                self.evalue_input = ui.number('E-value', value=1000.0, format='%.2f')
                self.max_amplicon_input = ui.number('Max amplicon size (bp)', value=500, format='%.0f')
                self.threads_input = ui.number('Threads', value=4, format='%.0f')

            with ui.expansion('Advanced options', icon='settings').classes('w-full'):
                self.analysis_id_input = ui.input('Analysis ID (optional)')
                self.db_output_input = ui.input('BLAST DB output dir', value='outputs/blast_dbs/')
                self.api_key_input = ui.input('NCBI API key (optional)').props('type=password')

        # --- Controls ---
        with ui.row().classes('justify-center m-4 gap-3'):
            self.run_button = ui.button('Start Pipeline', on_click=self.start_pipeline, color='primary').classes('w-48 h-12 text-base font-semibold')
            self.cancel_button = ui.button('Cancel', on_click=self.cancel_pipeline, color='negative').classes('w-28 h-12 text-base font-semibold')
            self.cancel_button.disable()

        # --- Progress ---
        with ui.card().classes('w-full max-w-5xl mx-auto m-4 p-4 shadow-lg border border-slate-200'):
            with ui.row().classes('items-center gap-4'):
                self.progress_bar = ui.linear_progress(value=0).classes('flex-1')
                self.stage_label = ui.label('Idle').classes('text-sm font-semibold text-slate-700')

        # --- Console ---
        with ui.card().classes('w-full max-w-5xl mx-auto m-4 p-4 shadow-lg border border-slate-200 bg-slate-950'):
            ui.label('Console Output').classes('text-lg font-semibold text-white')
            self.console = ui.log(max_lines=2000).classes(
                'w-full h-72 p-3 bg-slate-950 text-green-300 font-mono text-sm'
            )

        # --- Results (hidden until a run completes successfully) ---
        self.results_card = ui.card().classes('w-full max-w-5xl mx-auto m-4 p-5 shadow-lg border border-slate-200')
        self.results_card.visible = False
        with self.results_card:
            ui.label('Results').classes('text-xl font-semibold text-slate-900')
            self.report_viewer = ui.markdown('').classes('bg-slate-50 p-4 rounded-xl border border-slate-200')
            ui.label('Output files:').classes('font-semibold mt-4 text-slate-800')
            self.files_column = ui.column().classes('gap-2')

    # --- Handlers ---

    def handle_upload(self, e: events.UploadEventArguments):
        # NiceGUI 3.x: UploadEventArguments has a single `.file` (FileUpload)
        # with `.name`, `.content_type`, and a `.save(path)` method.
        try:
            name = e.file.name
            dest = UPLOAD_DIR / name
            e.file.save(str(dest))
        except Exception as exc:
            print(f"[handle_upload] failed: {exc!r}")  # visible in server terminal
            ui.notify(f'Upload failed: {exc}', type='negative')
            return

        self.uploaded_csv_path = dest
        self.upload_status.text = f'Uploaded: {name} ({dest.stat().st_size} bytes)'
        ui.notify(f'File uploaded: {name}', type='positive')

    def start_pipeline(self):
        if self.uploaded_csv_path is None:
            ui.notify('Please upload an oligo CSV first.', type='warning')
            return
        if not self.target_input.value.strip():
            ui.notify('Target organism is required.', type='warning')
            return

        background_list = [b.strip() for b in self.background_input.value.split(',') if b.strip()]
        if not background_list:
            ui.notify('At least one background organism is required.', type='warning')
            return

        config = {
            "oligos": str(self.uploaded_csv_path),
            "target": self.target_input.value.strip(),
            "background": background_list,
            "analysis_id": self.analysis_id_input.value.strip() or None,
            "db_output": self.db_output_input.value.strip() or "outputs/blast_dbs/",
            "threads": int(self.threads_input.value),
            "sensitivity": self.sensitivity_select.value,
            "evalue": float(self.evalue_input.value),
            "max_amplicon_size": int(self.max_amplicon_input.value),
            "api_key": self.api_key_input.value.strip() or None,
        }

        self.results_card.visible = False
        self.console.clear()
        self.console.push("--- Initialising pipeline execution ---")
        self.run_button.disable()
        self.cancel_button.enable()
        self.progress_bar.value = 0
        self.stage_label.text = 'Starting...'

        self.runner = PipelineRunner(main_py_path=MAIN_PY_PATH)
        self.runner.start(config)

        self.polling_timer = ui.timer(0.2, self.poll_runner)

    def cancel_pipeline(self):
        if self.runner:
            self.runner.cancel()
            self.console.push("\n[SYSTEM] Cancellation requested...")
            self.cancel_button.disable()

    def poll_runner(self):
        if not self.runner:
            return

        for line in self.runner.drain_logs():
            self.console.push(line)

        self.progress_bar.value = self.runner.progress
        self.stage_label.text = self.runner.current_stage_label

        # poll() is the authoritative completion signal — progress is only
        # an estimate and may not hit 1.0 exactly before the process exits.
        if self.runner.poll() is not None:
            self.polling_timer.deactivate()
            self.run_button.enable()
            self.cancel_button.disable()

            result = self.runner.get_result()

            if result.status == RunStatus.SUCCESS:
                self.console.push(f"\n[SYSTEM] Pipeline finished. Report: {result.report_path}")
                ui.notify('Pipeline execution finished!', type='positive')
                self.show_results(result)
            elif result.status == RunStatus.CANCELLED:
                self.console.push("\n[SYSTEM] Pipeline was cancelled.")
                ui.notify('Pipeline cancelled.', type='warning')
            else:
                self.console.push(f"\n[SYSTEM] Pipeline failed (exit code {result.return_code}). Check the log above.")
                ui.notify('Pipeline failed — see console output.', type='negative')

    def show_results(self, result):
        self.files_column.clear()

        if result.report_path and result.report_path.exists():
            self.report_viewer.set_content(result.report_path.read_text())
        else:
            self.report_viewer.set_content('_No report file found at the expected path._')

        if result.run_dir and result.run_dir.exists():
            with self.files_column:
                for path in sorted(result.run_dir.iterdir()):
                    if path.is_file():
                        size_kb = path.stat().st_size / 1024
                        ui.label(f'{path.name}  ({size_kb:.1f} KB)').classes('text-xs font-mono')

        self.results_card.visible = True


app = PipelineAppScaffold()
ui.run(title="Pipeline Orchestrator", host="0.0.0.0", port=8080)