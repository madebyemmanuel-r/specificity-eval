"""Provide a reusable helper for emitting machine-readable pipeline events.

This module centralises GUI event formatting so the pipeline can publish
structured status and result markers without mixing that logic into the
main application flow.
"""


def emit_gui_event(key: str, value: str) -> None:
    """Emit a dedicated machine-readable event for GUI/parser consumption."""
    print(f"GUI::{key}={value}", flush=True)
