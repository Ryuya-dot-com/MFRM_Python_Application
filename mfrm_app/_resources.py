"""Shared reporting text in installed wheels, source checkouts and desktop bundles."""
from importlib.resources import files
import json
from pathlib import Path


def metric_note(name: str) -> str:
    resource = files("mfrm_app").joinpath("data", "en.json")
    if not resource.is_file():
        # Source checkouts and existing desktop builds retain the shared locale here.
        resource = Path(__file__).resolve().parents[1] / "locales" / "en.json"
    return json.loads(resource.read_text(encoding="utf-8"))["metric_notes"][name]
