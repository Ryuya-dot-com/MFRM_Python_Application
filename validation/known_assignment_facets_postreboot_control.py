"""Execute the frozen post-reboot FACETS availability control."""

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


PLAN_PATH = ROOT / "validation" / "known_assignment_facets_postreboot_control_plan_20260811.json"
OUTPUT_DIR = ROOT / "validation" / "kafp_postreboot_control_20260811"


def main() -> None:
    if OUTPUT_DIR.exists():
        raise FileExistsError(f"Refusing to overwrite post-reboot control: {OUTPUT_DIR}")
    plan = json.loads(PLAN_PATH.read_text(encoding="utf-8"))
    trigger = ROOT / plan["trigger_audit"]
    source = ROOT / plan["control_spec"]
    executable = Path(plan["facets_executable"])
    for name, path, expected in (
        ("trigger audit", trigger, plan["trigger_audit_sha256"]),
        ("control spec", source, plan["control_spec_sha256"]),
        ("FACETS executable", executable, plan["facets_executable_sha256"]),
    ):
        if sha256_file(path).lower() != expected.lower():
            raise ValueError(f"{name} hash mismatch")

    temporary = Path(tempfile.mkdtemp(prefix="kafp_postreboot_"))
    score_base = (temporary / "scores.txt").resolve()
    lines = source.read_text(encoding="utf-8").splitlines()
    replaced = [
        f"Scorefile={score_base}" if line.startswith("Scorefile=") else line
        for line in lines
    ]
    if replaced == lines:
        raise ValueError("Scorefile specification was not replaced")
    spec = temporary / "analysis.txt"
    report = temporary / "report.txt"
    spec.write_text("\n".join(replaced) + "\n", encoding="utf-8", newline="\n")
    completed = invoke_facets(
        executable,
        spec,
        report,
        timeout_seconds=float(plan["execution"]["timeout_seconds"]),
    )
    (temporary / "stdout.txt").write_text(completed.stdout or "", encoding="utf-8")
    (temporary / "stderr.txt").write_text(completed.stderr or "", encoding="utf-8")
    score_files = sorted(temporary.glob("scores.*.txt"))
    passed = bool(
        completed.returncode == 0
        and report.is_file()
        and report.stat().st_size > 0
        and len(score_files) == 4
    )
    shutil.copytree(temporary, OUTPUT_DIR)
    result = {
        "schema_version": "known_assignment_facets_postreboot_control_result_v1",
        "pass": passed,
        "return_code": int(completed.returncode),
        "report_created": report.is_file(),
        "report_bytes": report.stat().st_size if report.is_file() else 0,
        "score_files": len(score_files),
        "temporary_execution_root": str(temporary),
        "supplemental_calibration_authorized": passed,
        "scientific_endpoint_read": False,
        "claim_boundary": plan["claim_boundary"],
    }
    result_path = OUTPUT_DIR / "result.json"
    result_path.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    artifacts = [path for path in sorted(OUTPUT_DIR.rglob("*")) if path.is_file()]
    manifest = {
        "schema_version": "known_assignment_facets_postreboot_control_manifest_v1",
        "dependency_sha256": {
            "plan": sha256_file(PLAN_PATH),
            "runner": sha256_file(Path(__file__).resolve()),
            "trigger_audit": sha256_file(trigger),
            "control_spec": sha256_file(source),
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
