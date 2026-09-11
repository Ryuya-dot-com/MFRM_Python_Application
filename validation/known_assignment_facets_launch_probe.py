"""Run the no-CREATE_NO_WINDOW FACETS localization probe."""

from __future__ import annotations

import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from validation.operating_characteristics_facets import sha256_file

PLAN_PATH = ROOT / "validation" / "known_assignment_facets_launch_probe_plan_20260811.json"
OUTPUT_DIR = ROOT / "validation" / "kafp_launch_20260811"


def main() -> None:
    if OUTPUT_DIR.exists():
        raise FileExistsError(f"Refusing to overwrite: {OUTPUT_DIR}")
    plan = json.loads(PLAN_PATH.read_text(encoding="utf-8"))
    executable = Path(plan["facets_executable"])
    if sha256_file(executable).lower() != plan["facets_executable_sha256"].lower():
        raise ValueError("FACETS executable hash mismatch")
    temporary = Path(tempfile.mkdtemp(prefix="kafp_launch_"))
    rows = []
    for probe_id, relative, expected_hash in plan["probes"]:
        source = ROOT / relative
        if sha256_file(source).lower() != expected_hash.lower():
            raise ValueError(f"Source hash mismatch: {probe_id}")
        work = temporary / probe_id
        work.mkdir()
        score = (work / "scores.txt").resolve()
        lines = source.read_text(encoding="utf-8").splitlines()
        lines = [f"Scorefile={score}" if line.startswith("Scorefile=") else line for line in lines]
        spec = work / "analysis.txt"
        report = work / "report.txt"
        spec.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
        completed = subprocess.run(
            [str(executable), "BATCH=YES", str(spec), str(report)],
            cwd=work,
            capture_output=True,
            text=True,
            timeout=float(plan["execution"]["timeout_seconds"]),
            check=False,
            creationflags=0,
        )
        (work / "stdout.txt").write_text(completed.stdout or "", encoding="utf-8")
        (work / "stderr.txt").write_text(completed.stderr or "", encoding="utf-8")
        rows.append({
            "ProbeId": probe_id,
            "ReturnCode": int(completed.returncode),
            "ReportCreated": report.is_file(),
            "ReportBytes": report.stat().st_size if report.is_file() else 0,
            "ScoreFiles": len(list(work.glob("scores.*.txt"))),
            "Succeeded": bool(completed.returncode == 0 and report.is_file()),
        })
    shutil.copytree(temporary, OUTPUT_DIR)
    result = {
        "schema_version": "known_assignment_facets_launch_probe_result_v1",
        "probes": rows,
        "both_pass": bool(all(row["Succeeded"] for row in rows)),
        "scientific_endpoint_read": False,
        "registered_calibration_replaced": False,
        "claim_boundary": plan["claim_boundary"],
    }
    (OUTPUT_DIR / "result.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, ensure_ascii=False))


if __name__ == "__main__":
    main()
