"""Run the frozen two-input FACETS exit-code localization probe."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from validation.operating_characteristics_facets import invoke_facets, sha256_file


PLAN_PATH = ROOT / "validation" / "known_assignment_facets_exit_probe_plan_20260811.json"
OUTPUT_DIR = ROOT / "validation" / "kafp_20260811"


def main() -> None:
    if OUTPUT_DIR.exists():
        raise FileExistsError(f"Refusing to overwrite probe: {OUTPUT_DIR}")
    plan = json.loads(PLAN_PATH.read_text(encoding="utf-8"))
    facets_exe = Path(plan["facets_executable"])
    if sha256_file(facets_exe).lower() != plan["facets_executable_sha256"].lower():
        raise ValueError("FACETS executable hash mismatch")
    for probe in plan["probes"]:
        source = ROOT / probe["source_spec"]
        if sha256_file(source).lower() != probe["source_spec_sha256"].lower():
            raise ValueError(f"Source specification hash mismatch: {probe['id']}")

    OUTPUT_DIR.mkdir(parents=False, exist_ok=False)
    rows: list[dict[str, object]] = []
    artifact_paths: list[Path] = []
    for probe in plan["probes"]:
        probe_dir = OUTPUT_DIR / probe["id"]
        probe_dir.mkdir()
        source = ROOT / probe["source_spec"]
        score_base = (probe_dir / "scores.txt").resolve()
        lines = source.read_text(encoding="utf-8").splitlines()
        replaced = [
            f"Scorefile={score_base}" if line.startswith("Scorefile=") else line
            for line in lines
        ]
        if replaced == lines:
            raise ValueError(f"Scorefile was not replaced: {probe['id']}")
        spec_path = probe_dir / "analysis.txt"
        report_path = probe_dir / "report.txt"
        spec_path.write_text("\n".join(replaced) + "\n", encoding="utf-8", newline="\n")
        completed = invoke_facets(
            facets_exe,
            spec_path,
            report_path,
            timeout_seconds=float(plan["execution"]["timeout_seconds"]),
        )
        (probe_dir / "stdout.txt").write_text(completed.stdout or "", encoding="utf-8")
        (probe_dir / "stderr.txt").write_text(completed.stderr or "", encoding="utf-8")
        score_files = sorted(probe_dir.glob("scores.*.txt"))
        rows.append(
            {
                "ProbeId": probe["id"],
                "ReturnCode": int(completed.returncode),
                "ReportCreated": report_path.is_file(),
                "ReportBytes": report_path.stat().st_size if report_path.is_file() else 0,
                "ScoreFiles": len(score_files),
                "Succeeded": bool(completed.returncode == 0 and report_path.is_file()),
            }
        )
        artifact_paths.extend(path for path in probe_dir.iterdir() if path.is_file())
    control, new = (bool(row["Succeeded"]) for row in rows)
    if not control and not new:
        classification = "both_fail"
    elif control and not new:
        classification = "control_passes_new_fails"
    elif control and new:
        classification = "both_pass"
    else:
        classification = "control_fails_new_passes"
    result = {
        "schema_version": "known_assignment_facets_exit_probe_result_v1",
        "classification": classification,
        "interpretation": plan["interpretation"][classification],
        "probes": rows,
        "scientific_endpoint_read": False,
        "registered_calibration_replaced": False,
        "claim_boundary": plan["claim_boundary"],
    }
    result_path = OUTPUT_DIR / "result.json"
    result_path.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    artifact_paths.append(result_path)
    manifest = {
        "schema_version": "known_assignment_facets_exit_probe_manifest_v1",
        "dependency_sha256": {
            "plan": sha256_file(PLAN_PATH),
            "runner": sha256_file(Path(__file__).resolve()),
            "facets_executable": sha256_file(facets_exe),
        },
        "artifact_sha256": {
            path.relative_to(OUTPUT_DIR).as_posix(): sha256_file(path)
            for path in artifact_paths
        },
    }
    (OUTPUT_DIR / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, ensure_ascii=False))


if __name__ == "__main__":
    main()
