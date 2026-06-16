"""Lightweight command-line diagnostics that do not require Streamlit."""

from __future__ import annotations

import importlib.metadata as importlib_metadata
import json
import sys
from pathlib import Path


RUNTIME_PACKAGE_FLOORS = {
    "numpy": "1.24",
    "pandas": "2.0",
    "scipy": "1.10",
    "plotly": "6.1.1",
    "kaleido": "1.0",
    "matplotlib": "3.8",
    "streamlit": "1.54",
    "openpyxl": "3.1",
    "python-docx": "1.0",
    "reportlab": "4.0",
    "arviz": "0.17",
    "netCDF4": "1.6",
    "pyarrow": "15.0",
    "networkx": "3.0",
}

BUNDLED_ANCHOR_ASSETS = [
    "anchor_table_blank.csv",
    "anchor_table_example.csv",
    "group_anchor_table_blank.csv",
    "group_anchor_table_example.csv",
    "anchor_user_guidelines.md",
]


def _version_tuple(version_text: str) -> tuple[int, ...]:
    parts = []
    for token in str(version_text).replace("-", ".").split("."):
        if token.isdigit():
            parts.append(int(token))
        else:
            digits = "".join(ch for ch in token if ch.isdigit())
            if digits:
                parts.append(int(digits))
            break
    return tuple(parts or [0])


def _version_meets_floor(version_text: str, floor_text: str) -> bool:
    observed = _version_tuple(version_text)
    floor = _version_tuple(floor_text)
    max_len = max(len(observed), len(floor))
    observed = observed + (0,) * (max_len - len(observed))
    floor = floor + (0,) * (max_len - len(floor))
    return observed >= floor


def _app_root() -> Path:
    return Path(__file__).resolve().parents[1]


def run_lightweight_doctor(json_output: bool = False) -> int:
    """Run dependency/file checks before Streamlit can be imported."""
    root = _app_root()
    rows: list[dict[str, str]] = []

    def record(check: str, status: str, detail: str) -> None:
        rows.append({"Check": check, "Status": status, "Detail": detail})

    record(
        "python",
        "ok" if sys.version_info >= (3, 11) else "warn",
        sys.version.split()[0],
    )
    for package_name, floor in RUNTIME_PACKAGE_FLOORS.items():
        try:
            version_text = importlib_metadata.version(package_name)
            status = "ok" if _version_meets_floor(version_text, floor) else "warn"
            record(f"package:{package_name}", status, f"{version_text} (required >= {floor})")
        except importlib_metadata.PackageNotFoundError:
            record(f"package:{package_name}", "error", f"missing (required >= {floor})")

    for asset_name in BUNDLED_ANCHOR_ASSETS:
        path = root / "anchor_templates_and_guideline" / asset_name
        try:
            size = path.stat().st_size
            record(f"asset:{asset_name}", "ok" if size else "error", f"{size} bytes at {path.name}")
        except OSError as exc:
            record(f"asset:{asset_name}", "error", str(exc))

    for locale_name in ("en.json", "ja.json"):
        path = root / "locales" / locale_name
        try:
            with path.open(encoding="utf-8") as fh:
                json.load(fh)
            record(f"locale:{locale_name}", "ok", "valid JSON")
        except Exception as exc:
            record(f"locale:{locale_name}", "error", str(exc))

    record(
        "doctor_mode",
        "warn",
        "lightweight doctor ran before Streamlit import; install missing packages for full app checks",
    )

    if json_output:
        print(json.dumps(rows, indent=2, ensure_ascii=False))
    else:
        width = max(len(row["Check"]) for row in rows) if rows else 5
        for row in rows:
            print(f"[{row['Status'].upper():5}] {row['Check']:<{width}}  {row['Detail']}")

    return 1 if any(row["Status"] == "error" for row in rows) else 0

