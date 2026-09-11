"""Capture one frozen FACETS 4.5.0 crash using Microsoft ProcDump."""

from __future__ import annotations

import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from validation.operating_characteristics_facets import sha256_file


PLAN_PATH = ROOT / "validation" / "known_assignment_facets_crash_capture_plan_20260811.json"
OUTPUT_DIR = ROOT / "validation" / "kafp_crash_capture_20260811"


def main() -> None:
    if OUTPUT_DIR.exists():
        raise FileExistsError(f"Refusing to overwrite FACETS crash capture: {OUTPUT_DIR}")
    plan = json.loads(PLAN_PATH.read_text(encoding="utf-8"))
    trigger = ROOT / plan["trigger_result"]
    source = Path(plan["official_example"])
    executable = Path(plan["facets_executable"])
    capture_tool = ROOT / plan["capture_tool"]
    for name, path, expected in (
        ("trigger result", trigger, plan["trigger_result_sha256"]),
        ("official example", source, plan["official_example_sha256"]),
        ("FACETS executable", executable, plan["facets_executable_sha256"]),
        ("capture tool", capture_tool, plan["capture_tool_sha256"]),
    ):
        if sha256_file(path).lower() != expected.lower():
            raise ValueError(f"{name} hash mismatch")

    temporary = Path(tempfile.mkdtemp(prefix="kafp_crash_capture_"))
    process_temp = temporary / "process_temp"
    dump_dir = temporary / "dumps"
    process_temp.mkdir()
    dump_dir.mkdir()
    spec = temporary / source.name
    report = temporary / "report.txt"
    shutil.copy2(source, spec)
    environment = os.environ.copy()
    environment["TEMP"] = str(process_temp.resolve())
    environment["TMP"] = str(process_temp.resolve())
    command = [
        str(capture_tool),
        "-accepteula",
        "-ma",
        "-e",
        "-x",
        str(dump_dir.resolve()),
        str(executable),
        "BATCH=YES",
        str(spec.resolve()),
        str(report.resolve()),
    ]
    completed = subprocess.run(
        command,
        cwd=temporary,
        capture_output=True,
        text=True,
        timeout=float(plan["execution"]["timeout_seconds"]),
        check=False,
        env=environment,
    )
    (temporary / "procdump_stdout.txt").write_text(
        completed.stdout or "", encoding="utf-8"
    )
    (temporary / "procdump_stderr.txt").write_text(
        completed.stderr or "", encoding="utf-8"
    )
    dump_files = [path for path in sorted(dump_dir.glob("*.dmp")) if path.stat().st_size]
    captured = bool(dump_files)
    shutil.copytree(temporary, OUTPUT_DIR)
    result = {
        "schema_version": "known_assignment_facets_crash_capture_result_v1",
        "diagnostic_capture_pass": captured,
        "procdump_return_code": int(completed.returncode),
        "dump_files": [
            {
                "name": path.name,
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
            for path in dump_files
        ],
        "report_created": report.is_file(),
        "report_bytes": report.stat().st_size if report.is_file() else 0,
        "existing_temp_files_deleted": False,
        "temporary_execution_root": str(temporary),
        "supplemental_calibration_authorized": False,
        "scientific_endpoint_read": False,
        "claim_boundary": plan["claim_boundary"],
    }
    (OUTPUT_DIR / "result.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    artifacts = [path for path in sorted(OUTPUT_DIR.rglob("*")) if path.is_file()]
    manifest = {
        "schema_version": "known_assignment_facets_crash_capture_manifest_v1",
        "dependency_sha256": {
            "plan": sha256_file(PLAN_PATH),
            "runner": sha256_file(Path(__file__).resolve()),
            "trigger_result": sha256_file(trigger),
            "official_example": sha256_file(source),
            "facets_executable": sha256_file(executable),
            "capture_tool": sha256_file(capture_tool),
        },
        "artifact_sha256": {
            path.relative_to(OUTPUT_DIR).as_posix(): sha256_file(path)
            for path in artifacts
        },
    }
    (OUTPUT_DIR / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, ensure_ascii=False))
    if not captured:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
