"""Execute the frozen two-control qualification of the shared FACETS helper."""

from __future__ import annotations

import json
from pathlib import Path
import re
import shutil
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from validation.operating_characteristics_facets import invoke_facets, sha256_file


PLAN_PATH = ROOT / "validation" / "facets_shared_helper_qualification_plan_20260811.json"
OUTPUT_DIR = ROOT / "validation" / "facets_shared_helper_qualification_20260811"


def _launcher_evidence(stderr: str) -> dict[str, object]:
    lines = [line for line in stderr.splitlines() if line.startswith("MFRM_FACETS_LAUNCH=")]
    if len(lines) != 1:
        raise ValueError(f"Expected one launcher evidence line, found {len(lines)}")
    return json.loads(lines[0].split("=", 1)[1])


def main() -> None:
    if OUTPUT_DIR.exists():
        raise FileExistsError(f"Refusing to overwrite helper qualification: {OUTPUT_DIR}")
    plan = json.loads(PLAN_PATH.read_text(encoding="utf-8"))
    dependencies = (
        ("trigger result", ROOT / plan["trigger_result"], plan["trigger_result_sha256"]),
        ("FACETS executable", Path(plan["facets_executable"]), plan["facets_executable_sha256"]),
        ("shared helper", ROOT / plan["shared_helper"], plan["shared_helper_sha256"]),
        ("visible launcher", ROOT / plan["visible_launcher"], plan["visible_launcher_sha256"]),
    )
    for name, path, expected in dependencies:
        if sha256_file(path).lower() != expected.lower():
            raise ValueError(f"{name} hash mismatch")

    temporary = Path(tempfile.mkdtemp(prefix="facets_shared_helper_"))
    control_results: list[dict[str, object]] = []
    for control in plan["controls"]:
        control_dir = temporary / control["id"]
        control_dir.mkdir()
        source = Path(control["specification"])
        if not source.is_absolute():
            source = ROOT / source
        if sha256_file(source).lower() != control["specification_sha256"].lower():
            raise ValueError(f"Control specification hash mismatch: {control['id']}")
        spec = control_dir / "analysis.txt"
        text = source.read_text(encoding="utf-8")
        score_base = control_dir / "scores.txt"
        if int(control["expected_score_files"]):
            pattern = re.compile(r"^\s*Scorefile\s*=.*$", re.IGNORECASE | re.MULTILINE)
            replacement = f"Scorefile={score_base.resolve()}"
            text, replacements = pattern.subn(lambda _match: replacement, text)
            if replacements != 1:
                raise ValueError(f"Expected one Scorefile replacement, found {replacements}")
        spec.write_text(text, encoding="utf-8", newline="\n")
        report = control_dir / "report.txt"
        completed = invoke_facets(
            Path(plan["facets_executable"]),
            spec,
            report,
            timeout_seconds=float(plan["execution"]["timeout_seconds_each"]),
        )
        (control_dir / "stdout.txt").write_text(completed.stdout or "", encoding="utf-8")
        (control_dir / "stderr.txt").write_text(completed.stderr or "", encoding="utf-8")
        evidence = _launcher_evidence(completed.stderr or "")
        scores = sorted(control_dir.glob("scores.*.txt"))
        expected_scores = int(control["expected_score_files"])
        passed = bool(
            completed.returncode == 0
            and report.is_file()
            and report.stat().st_size > 0
            and len(scores) == expected_scores
            and all(path.stat().st_size > 0 for path in scores)
            and int(evidence["windows_closed"]) >= 1
            and evidence["forced_termination"] is False
        )
        control_results.append({
            "control_id": control["id"],
            "pass": passed,
            "return_code": int(completed.returncode),
            "report_bytes": report.stat().st_size if report.is_file() else 0,
            "score_files": len(scores),
            "score_bytes": [path.stat().st_size for path in scores],
            "launcher_evidence": evidence,
        })

    passed = bool(all(row["pass"] for row in control_results))
    shutil.copytree(temporary, OUTPUT_DIR)
    result = {
        "schema_version": "facets_shared_helper_qualification_result_v1",
        "pass": passed,
        "controls": control_results,
        "supplemental_calibration_plan_authorized": passed,
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
