"""Audit preflight completion and external calibration before endpoint access."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import pandas as pd

from validation import known_assignment_multivector_preflight as preflight
from validation.operating_characteristics_facets import sha256_file


OUTPUT_PATH = preflight.STUDY_DIR / "preendpoint_operational_audit.json"
IDENTITY_PATH = preflight.STUDY_DIR / "preendpoint_operational_audit_identity.json"
PROBE_PATHS = (
    ROOT / "validation" / "kafp_20260811" / "result.json",
    ROOT / "validation" / "kafp_temp_20260811" / "result.json",
    ROOT / "validation" / "kafp_launch_20260811" / "result.json",
)


def main() -> None:
    if OUTPUT_PATH.exists() or IDENTITY_PATH.exists():
        raise FileExistsError("Refusing to overwrite pre-endpoint operational audit")
    identity = preflight.validate_study_identity()
    attempts = pd.read_csv(
        preflight.STUDY_DIR / "retained_input" / "attempt_manifest.csv"
    )
    rows = []
    marker_hashes = {}
    for _, attempt in attempts.iterrows():
        marker = preflight._completion_path(attempt)
        if not marker.is_file():
            raise FileNotFoundError(f"Missing completion marker: {attempt['AttemptId']}")
        completion = json.loads(marker.read_text(encoding="utf-8"))
        preflight._validate_completion(completion, attempt, identity)
        marker_hashes[str(attempt["AttemptId"])] = sha256_file(marker)
        facets_failure = ""
        if str(attempt["AttemptType"]) == "RESILIENT_FACETS_PYTHON_JMLE_PCM":
            metrics_path = (
                preflight.STUDY_DIR
                / str(completion["artifact_root"])
                / "resilient_pair_metrics.json"
            )
            metrics = json.loads(metrics_path.read_text(encoding="utf-8"))
            facets_failure = str(metrics.get("facets_final_error", ""))
        rows.append(
            {
                "AttemptType": str(attempt["AttemptType"]),
                "ExecutionCompleted": bool(completion["execution_completed"]),
                "StatisticalEvidenceReady": bool(
                    completion["statistical_evidence_ready"]
                ),
                "FACETSCalibrationReady": completion.get(
                    "facets_calibration_ready"
                ),
                "FailureReason": str(completion.get("failure_reason", "")),
                "FACETSFailureReason": facets_failure,
            }
        )
    outcomes = pd.DataFrame(rows)
    jmle = outcomes.loc[
        outcomes["AttemptType"].eq("RESILIENT_FACETS_PYTHON_JMLE_PCM")
    ]
    probes = [json.loads(path.read_text(encoding="utf-8")) for path in PROBE_PATHS]
    repeated_code = "3221225477"
    facets_failure_reasons = jmle["FACETSFailureReason"].astype(str)
    audit = {
        "schema_version": "known_assignment_multivector_preendpoint_audit_v1",
        "audit_timing": "After all 48 completion markers, before aggregate or registered endpoint access.",
        "completion_markers": len(marker_hashes),
        "execution_completed": int(outcomes["ExecutionCompleted"].sum()),
        "statistical_evidence_ready": int(
            outcomes["StatisticalEvidenceReady"].sum()
        ),
        "python_jmle_evidence_ready": int(
            jmle["StatisticalEvidenceReady"].sum()
        ),
        "facets_calibration_ready": int(
            jmle["FACETSCalibrationReady"].fillna(False).astype(bool).sum()
        ),
        "facets_pairs_planned": len(jmle),
        "facets_failures_all_same_access_violation": bool(
            len(jmle) == 12
            and facets_failure_reasons.str.contains(repeated_code, regex=False).all()
        ),
        "localization_probes": [
            {
                "path": path.relative_to(ROOT).as_posix(),
                "sha256": sha256_file(path),
                "classification": probe.get("classification"),
                "both_pass": probe.get("both_pass"),
            }
            for path, probe in zip(PROBE_PATHS, probes, strict=True)
        ],
        "operational_gate_pass": False,
        "blocking_gate": "all_12_facets_python_jmle_pairs_calibration_ready",
        "failure_localization": "Current FACETS 4.5.0 executable/runtime session: retained-success and new specs both fail in workspace and OS-temp paths, with and without CREATE_NO_WINDOW.",
        "scientific_endpoint_read": False,
        "aggregate_created": False,
        "python_evidence_invalidated_by_facets_failure": False,
        "next_action": "Restart Windows as recommended by the official FACETS/Winsteps problems page, rerun a retained-success control, and register any supplemental calibration execution before use.",
        "claim_boundary": "No preflight sign diagnostic, estimator ranking, or confirmatory claim is available while the external JMLE calibration gate is false.",
    }
    OUTPUT_PATH.write_text(
        json.dumps(audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    marker_digest = hashlib.sha256(
        json.dumps(marker_hashes, sort_keys=True, separators=(",", ":")).encode(
            "utf-8"
        )
    ).hexdigest()
    identity_record = {
        "schema_version": "known_assignment_multivector_preendpoint_audit_identity_v1",
        "study_identity_sha256": sha256_file(
            preflight.STUDY_DIR / "study_identity.json"
        ),
        "completion_marker_set_sha256": marker_digest,
        "audit_runner_sha256": sha256_file(Path(__file__).resolve()),
        "audit_sha256": sha256_file(OUTPUT_PATH),
    }
    IDENTITY_PATH.write_text(
        json.dumps(identity_record, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(audit, ensure_ascii=False))


if __name__ == "__main__":
    main()
