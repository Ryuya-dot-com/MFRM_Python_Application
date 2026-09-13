"""Run the frozen official-example qualification of the BATCH=NO launcher."""

from __future__ import annotations

import json
from pathlib import Path
import shutil
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from validation.facets_visible_launcher import invoke_facets_batch_no
from validation.operating_characteristics_facets import sha256_file


PLAN_PATH = ROOT / "validation" / "facets_visible_launcher_qualification_plan_20260811.json"
OUTPUT_DIR = ROOT / "validation" / "facets_visible_launcher_qualification_20260811"


def main() -> None:
    if OUTPUT_DIR.exists():
        raise FileExistsError(f"Refusing to overwrite launcher qualification: {OUTPUT_DIR}")
    plan = json.loads(PLAN_PATH.read_text(encoding="utf-8"))
    dependencies = (
        ("trigger report", ROOT / plan["trigger_control_report"], plan["trigger_control_report_sha256"]),
        ("launcher", ROOT / plan["launcher"], plan["launcher_sha256"]),
        ("FACETS executable", Path(plan["facets_executable"]), plan["facets_executable_sha256"]),
        ("official example", Path(plan["official_example"]), plan["official_example_sha256"]),
    )
    for name, path, expected in dependencies:
        if sha256_file(path).lower() != expected.lower():
            raise ValueError(f"{name} hash mismatch")

    temporary = Path(tempfile.mkdtemp(prefix="facets_visible_qualification_"))
    spec = temporary / "Kct.txt"
    report = temporary / "report.txt"
    shutil.copy2(Path(plan["official_example"]), spec)
    invocation = invoke_facets_batch_no(
        [str(Path(plan["facets_executable"])), "BATCH=NO", str(spec), str(report)],
        cwd=temporary,
        readiness_paths=(report,),
        timeout_seconds=float(plan["execution"]["timeout_seconds"]),
        stable_seconds=float(plan["execution"]["stable_seconds"]),
        hide_window=True,
    )
    completed = invocation.completed
    report_text = report.read_text(encoding="utf-8", errors="replace") if report.is_file() else ""
    expected_sections = {
        "table_7": "Table 7.1" in report_text,
        "table_8": "Table 8.1" in report_text,
        "table_4": "Table 4.1" in report_text,
        "title": "Knox Cube Test" in report_text,
    }
    passed = bool(
        completed.returncode == 0
        and report.is_file()
        and report.stat().st_size > 0
        and all(expected_sections.values())
        and invocation.windows_closed >= 1
        and not invocation.forced_termination
    )
    shutil.copytree(temporary, OUTPUT_DIR)
    result = {
        "schema_version": "facets_visible_launcher_qualification_result_v1",
        "pass": passed,
        "return_code": int(completed.returncode),
        "report_created": report.is_file(),
        "report_bytes": report.stat().st_size if report.is_file() else 0,
        "expected_sections": expected_sections,
        "stable_seconds": invocation.stable_seconds,
        "windows_closed": invocation.windows_closed,
        "forced_termination": invocation.forced_termination,
        "shared_helper_change_authorized": passed,
        "scientific_calibration_authorized": False,
        "scientific_endpoint_read": False,
        "claim_boundary": plan["claim_boundary"],
    }
    (OUTPUT_DIR / "result.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, ensure_ascii=False))
    if not passed:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
