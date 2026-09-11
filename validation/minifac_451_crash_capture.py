"""Capture one MINIFAC 4.5.1 mini dump for cross-version crash localization."""

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


PLAN_PATH = ROOT / "validation" / "minifac_451_crash_capture_plan_20260811.json"
OUTPUT_DIR = ROOT / "validation" / "minifac_451_crash_capture_20260811"


def main() -> None:
    if OUTPUT_DIR.exists():
        raise FileExistsError(f"Refusing to overwrite MINIFAC crash capture: {OUTPUT_DIR}")
    plan = json.loads(PLAN_PATH.read_text(encoding="utf-8"))
    dependencies = (
        ("trigger control", ROOT / plan["trigger_control"], plan["trigger_control_sha256"]),
        ("MINIFAC executable", Path(plan["executable"]), plan["executable_sha256"]),
        ("official example", Path(plan["official_example"]), plan["official_example_sha256"]),
        ("Xojo framework", Path(plan["shared_xojo_framework"]), plan["shared_xojo_framework_sha256"]),
        ("capture tool", ROOT / plan["capture_tool"], plan["capture_tool_sha256"]),
    )
    for name, path, expected in dependencies:
        if sha256_file(path).lower() != expected.lower():
            raise ValueError(f"{name} hash mismatch")

    executable = Path(plan["executable"])
    source = Path(plan["official_example"])
    capture_tool = ROOT / plan["capture_tool"]
    temporary = Path(tempfile.mkdtemp(prefix="minifac_451_crash_"))
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
    completed = subprocess.run(
        [
            str(capture_tool), "-accepteula", "-mm", "-e", "-x",
            str(dump_dir.resolve()), str(executable), "BATCH=YES",
            str(spec.resolve()), str(report.resolve()),
        ],
        cwd=temporary,
        capture_output=True,
        timeout=float(plan["execution"]["timeout_seconds"]),
        check=False,
        env=environment,
    )
    (temporary / "procdump_stdout.bin").write_bytes(completed.stdout or b"")
    (temporary / "procdump_stderr.bin").write_bytes(completed.stderr or b"")
    dump_files = [path for path in sorted(dump_dir.glob("*.dmp")) if path.stat().st_size]
    captured = bool(dump_files)
    shutil.copytree(temporary, OUTPUT_DIR)
    result = {
        "schema_version": "minifac_451_crash_capture_result_v1",
        "diagnostic_capture_pass": captured,
        "procdump_return_code": int(completed.returncode),
        "dump_files": [
            {"name": path.name, "bytes": path.stat().st_size, "sha256": sha256_file(path)}
            for path in dump_files
        ],
        "report_created": report.is_file(),
        "supplemental_calibration_authorized": False,
        "scientific_endpoint_read": False,
        "claim_boundary": plan["claim_boundary"],
    }
    (OUTPUT_DIR / "result.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    artifacts = [path for path in sorted(OUTPUT_DIR.rglob("*")) if path.is_file()]
    manifest = {
        "schema_version": "minifac_451_crash_capture_manifest_v1",
        "dependency_sha256": {
            "plan": sha256_file(PLAN_PATH),
            **{name: sha256_file(path) for name, path, _expected in dependencies},
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
