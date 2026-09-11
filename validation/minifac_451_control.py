"""Execute the frozen MINIFAC 4.5.1 bundled-example version control."""

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


PLAN_PATH = ROOT / "validation" / "minifac_451_control_plan_20260811.json"
OUTPUT_DIR = ROOT / "validation" / "minifac_451_version_control_20260811"


def main() -> None:
    if OUTPUT_DIR.exists():
        raise FileExistsError(f"Refusing to overwrite MINIFAC control: {OUTPUT_DIR}")
    plan = json.loads(PLAN_PATH.read_text(encoding="utf-8"))
    dependencies = (
        (
            "version sensitivity plan",
            ROOT / plan["version_sensitivity_plan"],
            plan["version_sensitivity_plan_sha256"],
        ),
        ("install result", ROOT / plan["install_result"], plan["install_result_sha256"]),
        ("MINIFAC executable", Path(plan["executable"]), plan["executable_sha256"]),
        ("official example", Path(plan["official_example"]), plan["official_example_sha256"]),
        (
            "shared Xojo framework",
            Path(plan["shared_xojo_framework"]),
            plan["shared_xojo_framework_sha256"],
        ),
    )
    for name, path, expected in dependencies:
        if sha256_file(path).lower() != expected.lower():
            raise ValueError(f"{name} hash mismatch")

    executable = Path(plan["executable"])
    source = Path(plan["official_example"])
    temporary = Path(tempfile.mkdtemp(prefix="minifac_451_control_"))
    process_temp = temporary / "process_temp"
    process_temp.mkdir()
    spec = temporary / source.name
    report = temporary / "report.txt"
    shutil.copy2(source, spec)
    environment = os.environ.copy()
    environment["TEMP"] = str(process_temp.resolve())
    environment["TMP"] = str(process_temp.resolve())
    completed = subprocess.run(
        [str(executable), "BATCH=YES", str(spec), str(report)],
        cwd=temporary,
        capture_output=True,
        text=True,
        timeout=float(plan["execution"]["timeout_seconds"]),
        check=False,
        creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0),
        env=environment,
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
        "schema_version": "minifac_451_control_result_v1",
        "pass": passed,
        "return_code": int(completed.returncode),
        "return_code_hex_unsigned": f"0x{int(completed.returncode) & 0xFFFFFFFF:08X}",
        "report_created": report.is_file(),
        "report_bytes": report.stat().st_size if report.is_file() else 0,
        "version_sensitivity": "4.5.0-specific main-executable path" if passed else "not resolved by 4.5.1",
        "facets_450_replaced": False,
        "supplemental_calibration_authorized": False,
        "scientific_endpoint_read": False,
        "claim_boundary": plan["claim_boundary"],
    }
    (OUTPUT_DIR / "result.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    artifacts = [path for path in sorted(OUTPUT_DIR.rglob("*")) if path.is_file()]
    manifest = {
        "schema_version": "minifac_451_control_manifest_v1",
        "dependency_sha256": {
            "plan": sha256_file(PLAN_PATH),
            **{
                name: sha256_file(path)
                for name, path, _expected in dependencies
            },
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
