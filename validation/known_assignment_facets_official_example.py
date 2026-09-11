"""Execute the frozen bundled-example FACETS 4.5.0 diagnostic."""

from __future__ import annotations

import json
from pathlib import Path
import shutil
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from validation.operating_characteristics_facets import invoke_facets, sha256_file


PLAN_PATH = ROOT / "validation" / "known_assignment_facets_official_example_plan_20260811.json"
OUTPUT_DIR = ROOT / "validation" / "kafp_official_example_20260811"


def main() -> None:
    if OUTPUT_DIR.exists():
        raise FileExistsError(f"Refusing to overwrite official-example diagnostic: {OUTPUT_DIR}")
    plan = json.loads(PLAN_PATH.read_text(encoding="utf-8"))
    trigger = ROOT / plan["trigger_result"]
    source = Path(plan["official_example"])
    executable = Path(plan["facets_executable"])
    for name, path, expected in (
        ("trigger result", trigger, plan["trigger_result_sha256"]),
        ("official example", source, plan["official_example_sha256"]),
        ("FACETS executable", executable, plan["facets_executable_sha256"]),
    ):
        if sha256_file(path).lower() != expected.lower():
            raise ValueError(f"{name} hash mismatch")

    temporary = Path(tempfile.mkdtemp(prefix="kafp_official_example_"))
    spec = temporary / source.name
    report = temporary / "report.txt"
    shutil.copy2(source, spec)
    completed = invoke_facets(
        executable,
        spec,
        report,
        timeout_seconds=float(plan["execution"]["timeout_seconds"]),
    )
    (temporary / "stdout.txt").write_text(completed.stdout or "", encoding="utf-8")
    (temporary / "stderr.txt").write_text(completed.stderr or "", encoding="utf-8")
    passed = bool(
        completed.returncode == 0
        and report.is_file()
        and report.stat().st_size > 0
    )
    shutil.copytree(temporary, OUTPUT_DIR)
    result = {
        "schema_version": "known_assignment_facets_official_example_result_v1",
        "pass": passed,
        "return_code": int(completed.returncode),
        "return_code_hex_unsigned": f"0x{int(completed.returncode) & 0xFFFFFFFF:08X}",
        "report_created": report.is_file(),
        "report_bytes": report.stat().st_size if report.is_file() else 0,
        "temporary_execution_root": str(temporary),
        "failure_scope": None if passed else "FACETS executable or runtime, independent of project input",
        "supplemental_calibration_authorized": passed,
        "scientific_endpoint_read": False,
        "claim_boundary": plan["claim_boundary"],
    }
    result_path = OUTPUT_DIR / "result.json"
    result_path.write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    artifacts = [path for path in sorted(OUTPUT_DIR.rglob("*")) if path.is_file()]
    manifest = {
        "schema_version": "known_assignment_facets_official_example_manifest_v1",
        "dependency_sha256": {
            "plan": sha256_file(PLAN_PATH),
            "runner": sha256_file(Path(__file__).resolve()),
            "trigger_result": sha256_file(trigger),
            "official_example": sha256_file(source),
            "facets_executable": sha256_file(executable),
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
    if not passed:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
