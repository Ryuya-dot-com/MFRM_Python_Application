"""Audit and aggregate the frozen multi-vector preflight with FACETS evidence.

The source 48-attempt study and the 12-dataset visible-mode FACETS supplement
remain immutable.  This runner first proves a one-to-one evidence join and
native Python replay, then creates a separate derivative aggregate.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd

from validation import known_assignment_multivector_preflight as source
from validation.facets_resilient_pair import dependency_manifest_digest
from validation.operating_characteristics_facets import sha256_file


SCHEMA_VERSION = "known_assignment_multivector_visible_join_v1"
PLAN_PATH = ROOT / "validation" / "known_assignment_multivector_visible_join_plan_20260811.json"
REGISTRATION_PATH = ROOT / "validation" / "known_assignment_multivector_visible_join_registration_20260811.json"
STUDY_DIR = ROOT / "validation" / "known_assignment_multivector_preflight4_20260811"
SUPPLEMENT_DIR = STUDY_DIR / "v"
AUDIT_PATH = STUDY_DIR / "visible_join_audit_20260811.json"
AUDIT_PAIRS_PATH = STUDY_DIR / "visible_join_python_replay_20260811.csv"
AGGREGATE_DIR = STUDY_DIR / "aggregate_visible_join_20260811"
JMLE_ATTEMPT = "RESILIENT_FACETS_PYTHON_JMLE_PCM"
REPLAY_TOLERANCE = 1e-12


def _json_dump(path: Path, value: Any) -> None:
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _hash_mapping(mapping: dict[str, str]) -> str:
    return hashlib.sha256(
        json.dumps(mapping, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def _validate_registration() -> tuple[dict[str, Any], dict[str, Any]]:
    plan = json.loads(PLAN_PATH.read_text(encoding="utf-8"))
    registration = json.loads(REGISTRATION_PATH.read_text(encoding="utf-8"))
    expected = {
        "plan_sha256": sha256_file(PLAN_PATH),
        "runner_sha256": sha256_file(Path(__file__).resolve()),
    }
    for key, digest in expected.items():
        if str(registration.get(key, "")).lower() != digest.lower():
            raise ValueError(f"Visible-join registration mismatch: {key}")
    return plan, registration


def _validate_source_identity(plan: dict[str, Any]) -> dict[str, Any]:
    identity_path = ROOT / plan["source_study_identity"]
    if sha256_file(identity_path) != plan["source_study_identity_sha256"]:
        raise ValueError("Frozen source-study identity changed")
    if sha256_file(ROOT / plan["source_preendpoint_audit"]) != plan["source_preendpoint_audit_sha256"]:
        raise ValueError("Frozen source pre-endpoint audit changed")
    if sha256_file(ROOT / plan["source_runner"]) != plan["source_runner_sha256"]:
        raise ValueError("Frozen source runner changed")
    identity = json.loads(identity_path.read_text(encoding="utf-8"))
    if not bool(identity.get("all_checks_pass")) or int(identity.get("attempts", 0)) != 48:
        raise ValueError("Frozen source identity is not a 48-attempt PASS")
    if identity["runner_sha256"] != plan["source_runner_sha256"]:
        raise ValueError("Frozen source runner identity mismatch")
    if sha256_file(source.PLAN_PATH) != identity["plan_sha256"]:
        raise ValueError("Frozen source plan changed")
    if sha256_file(source.REGISTRATION_PATH) != identity["registration_sha256"]:
        raise ValueError("Frozen source registration changed")
    retained_dependencies = json.loads(
        (STUDY_DIR / "dependency_manifest.json").read_text(encoding="utf-8")
    )
    if dependency_manifest_digest(retained_dependencies) != identity["dependency_manifest_sha256"]:
        raise ValueError("Frozen source dependency manifest is internally inconsistent")
    for filename, digest in identity["retained_input_sha256"].items():
        if sha256_file(STUDY_DIR / "retained_input" / filename) != digest:
            raise ValueError(f"Frozen source input changed: {filename}")
    return identity


def _validate_original_attempts(
    attempts: pd.DataFrame, identity: dict[str, Any]
) -> tuple[dict[str, dict[str, Any]], dict[str, str]]:
    completions: dict[str, dict[str, Any]] = {}
    marker_hashes: dict[str, str] = {}
    for _, attempt in attempts.iterrows():
        marker_path = source._completion_path(attempt)  # pylint: disable=protected-access
        if not marker_path.is_file():
            raise FileNotFoundError(f"Missing frozen completion marker: {attempt['AttemptId']}")
        completion = json.loads(marker_path.read_text(encoding="utf-8"))
        source._validate_completion(  # pylint: disable=protected-access
            completion, attempt, identity
        )
        attempt_id = str(attempt["AttemptId"])
        completions[attempt_id] = completion
        marker_hashes[attempt_id] = sha256_file(marker_path)
    if len(completions) != 48:
        raise ValueError(f"Expected 48 frozen completions, found {len(completions)}")
    return completions, marker_hashes


def _validate_supplement(
    plan: dict[str, Any], attempts: pd.DataFrame
) -> tuple[dict[str, Any], dict[str, dict[str, Any]], dict[str, str]]:
    assessment_path = ROOT / plan["supplement_assessment"]
    registration_path = ROOT / plan["supplement_registration"]
    if sha256_file(assessment_path) != plan["supplement_assessment_sha256"]:
        raise ValueError("FACETS supplement assessment changed")
    if sha256_file(registration_path) != plan["supplement_registration_sha256"]:
        raise ValueError("FACETS supplement registration changed")
    assessment = json.loads(assessment_path.read_text(encoding="utf-8"))
    if not bool(assessment.get("pass")) or not bool(assessment.get("aggregate_join_authorized")):
        raise ValueError("FACETS supplement did not authorize the join")
    if bool(assessment.get("scientific_endpoint_read")) or bool(assessment.get("aggregate_created")):
        raise ValueError("FACETS supplement exceeded its calibration-only boundary")
    if sha256_file(registration_path) != assessment["registration_sha256"]:
        raise ValueError("FACETS assessment registration identity mismatch")

    jmle_attempts = attempts.loc[attempts["AttemptType"].astype(str).eq(JMLE_ATTEMPT)]
    markers: dict[str, dict[str, Any]] = {}
    marker_hashes: dict[str, str] = {}
    for _, attempt in jmle_attempts.iterrows():
        marker_name = f"{int(attempt['AttemptOrdinal']):05d}.json"
        marker_path = SUPPLEMENT_DIR / "markers" / marker_name
        if not marker_path.is_file():
            raise FileNotFoundError(f"Missing FACETS supplement marker: {marker_name}")
        digest = sha256_file(marker_path)
        if assessment["marker_sha256"].get(marker_name) != digest:
            raise ValueError(f"FACETS supplement marker identity mismatch: {marker_name}")
        marker = json.loads(marker_path.read_text(encoding="utf-8"))
        expected = {
            "attempt_id": str(attempt["AttemptId"]),
            "run_id": str(attempt["RunId"]),
            "run_input_sha256": str(attempt["RunInputSHA256"]),
        }
        for key, value in expected.items():
            if str(marker.get(key, "")) != value:
                raise ValueError(f"FACETS supplement join-key mismatch: {key}")
        if not bool(marker.get("calibration_ready")) or not bool(marker.get("direct_agreement_pass")):
            raise ValueError(f"FACETS supplement calibration failed: {marker_name}")
        if bool(marker.get("scientific_endpoint_read")):
            raise ValueError(f"FACETS supplement marker exceeded claim boundary: {marker_name}")
        artifact_root = SUPPLEMENT_DIR / marker["artifact_root"]
        for relative, expected_digest in marker["artifact_sha256"].items():
            if sha256_file(artifact_root / relative) != expected_digest:
                raise ValueError(f"FACETS supplement artifact changed: {relative}")
        attempt_id = str(attempt["AttemptId"])
        markers[attempt_id] = marker
        marker_hashes[attempt_id] = digest
    if len(markers) != 12 or len(assessment["marker_sha256"]) != 12:
        raise ValueError("FACETS supplement join is not 12-of-12")
    return assessment, markers, marker_hashes


def _numeric_frame_difference(
    left: pd.DataFrame,
    right: pd.DataFrame,
    *,
    keys: list[str],
    values: list[str],
) -> tuple[float, int]:
    for name, frame in (("source", left), ("supplement", right)):
        missing = set(keys + values) - set(frame.columns)
        if missing:
            raise ValueError(f"{name} replay frame missing columns: {sorted(missing)}")
        if frame.duplicated(keys).any():
            raise ValueError(f"{name} replay keys are not unique: {keys}")
    merged = left[keys + values].merge(
        right[keys + values],
        on=keys,
        how="outer",
        suffixes=("Source", "Supplement"),
        indicator=True,
        validate="one_to_one",
    )
    if not merged["_merge"].eq("both").all():
        raise ValueError(f"Python replay keys differ: {keys}")
    maximum = 0.0
    for column in values:
        source_values = pd.to_numeric(merged[f"{column}Source"], errors="coerce").to_numpy(float)
        supplement_values = pd.to_numeric(
            merged[f"{column}Supplement"], errors="coerce"
        ).to_numpy(float)
        if not np.array_equal(np.isnan(source_values), np.isnan(supplement_values)):
            raise ValueError(f"Python replay missingness differs: {column}")
        finite = np.isfinite(source_values) & np.isfinite(supplement_values)
        if finite.any():
            maximum = max(
                maximum,
                float(np.max(np.abs(source_values[finite] - supplement_values[finite]))),
            )
    return maximum, len(merged)


def _python_replay_audit(
    attempts: pd.DataFrame,
    original: dict[str, dict[str, Any]],
    supplement: dict[str, dict[str, Any]],
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    selected = attempts.loc[attempts["AttemptType"].astype(str).eq(JMLE_ATTEMPT)]
    for _, attempt in selected.iterrows():
        attempt_id = str(attempt["AttemptId"])
        original_root = STUDY_DIR / original[attempt_id]["artifact_root"]
        supplement_root = SUPPLEMENT_DIR / supplement[attempt_id]["artifact_root"]
        original_runs = pd.read_csv(original_root / "run_ledger.csv")
        supplement_runs = pd.read_csv(supplement_root / "combined_run_ledger.csv")
        original_python_run = original_runs.loc[
            original_runs["EstimatorMode"].astype(str).eq(source.PYTHON_JMLE_MODE)
        ]
        supplement_python_run = supplement_runs.loc[
            supplement_runs["EstimatorMode"].astype(str).eq(source.PYTHON_JMLE_MODE)
        ]
        if len(original_python_run) != 1 or len(supplement_python_run) != 1:
            raise ValueError(f"Expected one Python JMLE run row: {attempt_id}")
        for column in ("RunId", "EstimatorMode", "EstimandClass", "IncludedInStudy"):
            if str(original_python_run.iloc[0][column]) != str(supplement_python_run.iloc[0][column]):
                raise ValueError(f"Python JMLE run identity differs: {attempt_id} {column}")

        original_recovery = pd.read_csv(original_root / "recovery.csv")
        supplement_recovery = pd.read_csv(supplement_root / "combined_recovery.csv")
        original_recovery = original_recovery.loc[
            original_recovery["EstimatorMode"].astype(str).eq(source.PYTHON_JMLE_MODE)
        ]
        supplement_recovery = supplement_recovery.loc[
            supplement_recovery["EstimatorMode"].astype(str).eq(source.PYTHON_JMLE_MODE)
        ]
        recovery_max, recovery_rows = _numeric_frame_difference(
            original_recovery,
            supplement_recovery,
            keys=["RunId", "Facet", "Level"],
            values=["Truth", "EstimateAligned", "TruthAligned", "ErrorAligned"],
        )

        original_thresholds = pd.read_csv(original_root / "thresholds.csv")
        supplement_thresholds = pd.read_csv(supplement_root / "combined_thresholds.csv")
        original_thresholds = original_thresholds.loc[
            original_thresholds["EstimatorMode"].astype(str).eq(source.PYTHON_JMLE_MODE)
        ]
        supplement_thresholds = supplement_thresholds.loc[
            supplement_thresholds["EstimatorMode"].astype(str).eq(source.PYTHON_JMLE_MODE)
        ]
        threshold_max, threshold_rows = _numeric_frame_difference(
            original_thresholds,
            supplement_thresholds,
            keys=["RunId", "StepFacetLevel", "Category"],
            values=["Estimate", "ThresholdTruth", "TruthError"],
        )
        rows.append({
            "AttemptId": attempt_id,
            "RunId": str(attempt["RunId"]),
            "RecoveryRows": recovery_rows,
            "ThresholdRows": threshold_rows,
            "RecoveryMaximumAbsoluteDifference": recovery_max,
            "ThresholdMaximumAbsoluteDifference": threshold_max,
            "WithinTolerance": bool(max(recovery_max, threshold_max) <= REPLAY_TOLERANCE),
        })
    return pd.DataFrame(rows)


def _build_join_audit() -> tuple[dict[str, Any], pd.DataFrame]:
    plan, registration = _validate_registration()
    identity = _validate_source_identity(plan)
    attempts = pd.read_csv(STUDY_DIR / "retained_input" / "attempt_manifest.csv")
    original, original_hashes = _validate_original_attempts(attempts, identity)
    assessment, supplement, supplement_hashes = _validate_supplement(plan, attempts)
    replay = _python_replay_audit(attempts, original, supplement)
    gates = {
        "original_completion_markers_48": len(original) == 48,
        "supplement_markers_12": len(supplement) == 12,
        "supplement_assessment_pass": bool(assessment["pass"]),
        "supplement_calibration_12_of_12": all(
            bool(marker["calibration_ready"]) for marker in supplement.values()
        ),
        "one_to_one_join_keys": len(set(original) & set(supplement)) == 12,
        "original_python_replay_12_of_12": len(replay) == 12
        and bool(replay["WithinTolerance"].all()),
        "source_identity_and_artifacts_unchanged": True,
        "supplement_artifacts_unchanged": True,
        "no_facets_rounded_fit_input": True,
    }
    audit = {
        "schema_version": f"{SCHEMA_VERSION}_audit",
        "pass": bool(all(gates.values())),
        "gates": gates,
        "source_completion_marker_set_sha256": _hash_mapping(original_hashes),
        "supplement_marker_set_sha256": _hash_mapping(supplement_hashes),
        "source_study_identity_sha256": sha256_file(STUDY_DIR / "study_identity.json"),
        "supplement_assessment_sha256": sha256_file(SUPPLEMENT_DIR / "assessment.json"),
        "plan_sha256": sha256_file(PLAN_PATH),
        "registration_sha256": sha256_file(REGISTRATION_PATH),
        "runner_sha256": registration["runner_sha256"],
        "maximum_original_python_recovery_difference": float(
            replay["RecoveryMaximumAbsoluteDifference"].max()
        ),
        "maximum_original_python_threshold_difference": float(
            replay["ThresholdMaximumAbsoluteDifference"].max()
        ),
        "python_replay_tolerance": REPLAY_TOLERANCE,
        "scientific_endpoint_read": False,
        "pf1_pf4_computed": False,
        "aggregate_created": False,
        "claim_boundary": plan["claim_limit"],
    }
    return audit, replay


def audit_join() -> dict[str, Any]:
    if AUDIT_PATH.exists() or AUDIT_PAIRS_PATH.exists():
        raise FileExistsError("Refusing to overwrite visible-join audit artifacts")
    audit, replay = _build_join_audit()
    replay.to_csv(AUDIT_PAIRS_PATH, index=False, lineterminator="\n")
    audit["python_replay_audit_sha256"] = sha256_file(AUDIT_PAIRS_PATH)
    _json_dump(AUDIT_PATH, audit)
    print(json.dumps(audit, ensure_ascii=False, sort_keys=True))
    if not audit["pass"]:
        raise SystemExit(2)
    return audit


def _read_source_frames(
    attempts: pd.DataFrame,
    completions: dict[str, dict[str, Any]],
    supplement: dict[str, dict[str, Any]],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    run_parts: list[pd.DataFrame] = []
    recovery_parts: list[pd.DataFrame] = []
    threshold_parts: list[pd.DataFrame] = []
    constraint_parts: list[pd.DataFrame] = []
    outcomes: list[dict[str, Any]] = []
    for _, attempt in attempts.iterrows():
        attempt_id = str(attempt["AttemptId"])
        completion = completions[attempt_id]
        root = STUDY_DIR / completion["artifact_root"]
        for filename, target in (
            ("run_ledger.csv", run_parts),
            ("recovery.csv", recovery_parts),
            ("thresholds.csv", threshold_parts),
            ("constraints.csv", constraint_parts),
        ):
            frame = source._read_nonempty(root / filename)  # pylint: disable=protected-access
            if not frame.empty:
                frame.insert(0, "AttemptId", attempt_id)
                target.append(frame)
        original_calibration = completion.get("facets_calibration_ready")
        joined_calibration = original_calibration
        calibration_source = "original_completion"
        if str(attempt["AttemptType"]) == JMLE_ATTEMPT:
            joined_calibration = bool(supplement[attempt_id]["calibration_ready"])
            calibration_source = "visible_facets_supplement_v3"
        outcome = {
            "AttemptId": attempt_id,
            "RunId": str(attempt["RunId"]),
            "AttemptType": str(attempt["AttemptType"]),
            "PersonVector": int(attempt["PersonVector"]),
            "Gamma": float(attempt["Gamma"]),
            "ExecutionCompleted": bool(completion["execution_completed"]),
            "StatisticalEvidenceReady": bool(completion["statistical_evidence_ready"]),
            "FACETSCalibrationReady": joined_calibration,
            "OriginalFACETSCalibrationReady": original_calibration,
            "CalibrationEvidenceSource": calibration_source,
            "FailureReason": completion.get("failure_reason", ""),
        }
        outcomes.append(outcome)
    return (
        pd.concat(run_parts, ignore_index=True, sort=False),
        pd.concat(recovery_parts, ignore_index=True, sort=False),
        pd.concat(threshold_parts, ignore_index=True, sort=False),
        pd.concat(constraint_parts, ignore_index=True, sort=False),
        pd.DataFrame(outcomes),
    )


def _calibration_frames(
    attempts: pd.DataFrame, supplement: dict[str, dict[str, Any]]
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    frames: dict[str, list[pd.DataFrame]] = {
        "runs": [], "recovery": [], "thresholds": []
    }
    selected = attempts.loc[attempts["AttemptType"].astype(str).eq(JMLE_ATTEMPT)]
    for _, attempt in selected.iterrows():
        attempt_id = str(attempt["AttemptId"])
        root = SUPPLEMENT_DIR / supplement[attempt_id]["artifact_root"]
        for key, filename in (
            ("runs", "combined_run_ledger.csv"),
            ("recovery", "combined_recovery.csv"),
            ("thresholds", "combined_thresholds.csv"),
        ):
            frame = pd.read_csv(root / filename)
            frame = frame.loc[
                frame["EstimatorMode"].astype(str).eq(source.FACETS_MODE)
            ].copy()
            frame.insert(0, "AttemptId", attempt_id)
            frame["IncludedInStudy"] = False
            frame["CalibrationOnly"] = True
            frame["ScientificEndpointEligible"] = False
            frames[key].append(frame)
    return tuple(pd.concat(frames[key], ignore_index=True, sort=False) for key in (
        "runs", "recovery", "thresholds"
    ))


def aggregate_joined() -> dict[str, Any]:
    if AGGREGATE_DIR.exists():
        raise FileExistsError("Refusing to overwrite the visible-join aggregate")
    if not AUDIT_PATH.is_file() or not AUDIT_PAIRS_PATH.is_file():
        raise FileNotFoundError("Visible-join audit must pass before aggregation")
    retained_audit = json.loads(AUDIT_PATH.read_text(encoding="utf-8"))
    rebuilt_audit, rebuilt_replay = _build_join_audit()
    if retained_audit.get("python_replay_audit_sha256") != sha256_file(AUDIT_PAIRS_PATH):
        raise ValueError("Visible-join Python replay audit changed")
    retained_core = dict(retained_audit)
    retained_core.pop("python_replay_audit_sha256", None)
    if retained_core != rebuilt_audit:
        raise ValueError("Visible-join audit no longer reproduces")
    retained_replay = pd.read_csv(AUDIT_PAIRS_PATH)
    if retained_replay.to_csv(index=False, lineterminator="\n") != rebuilt_replay.to_csv(
        index=False, lineterminator="\n"
    ):
        raise ValueError("Visible-join replay ledger no longer reproduces")
    if not bool(retained_audit.get("pass")):
        raise ValueError("Visible-join audit is not PASS")

    plan, _registration = _validate_registration()
    identity = _validate_source_identity(plan)
    input_dir = STUDY_DIR / "retained_input"
    attempts = pd.read_csv(input_dir / "attempt_manifest.csv")
    completions, original_hashes = _validate_original_attempts(attempts, identity)
    assessment, supplement, supplement_hashes = _validate_supplement(plan, attempts)
    raw_runs, raw_recovery, raw_thresholds, raw_constraints, outcomes = _read_source_frames(
        attempts, completions, supplement
    )
    manifest = pd.read_csv(input_dir / "manifest.csv")
    runs = source._attach_context(raw_runs, manifest)  # pylint: disable=protected-access
    recovery = source._attach_context(raw_recovery, manifest)  # pylint: disable=protected-access
    thresholds = source._attach_context(raw_thresholds, manifest)  # pylint: disable=protected-access
    constraints = source._attach_context(  # pylint: disable=protected-access
        raw_constraints, manifest
    )
    calibration_runs, calibration_recovery, calibration_thresholds = _calibration_frames(
        attempts, supplement
    )
    calibration_runs = source._attach_context(calibration_runs, manifest)  # pylint: disable=protected-access
    calibration_recovery = source._attach_context(  # pylint: disable=protected-access
        calibration_recovery, manifest
    )
    calibration_thresholds = source._attach_context(  # pylint: disable=protected-access
        calibration_thresholds, manifest
    )

    rater_errors = recovery.loc[
        recovery["EstimatorMode"].astype(str).isin(source.SCIENTIFIC_MODES)
        & recovery["IncludedInStudy"].fillna(False).astype(bool)
        & recovery["Facet"].astype(str).eq("Rater")
    ].copy()
    rater_errors["ErrorAligned"] = pd.to_numeric(
        rater_errors["ErrorAligned"], errors="coerce"
    )
    rater_loss = (
        rater_errors.groupby(
            ["PersonVector", "Gamma", "EstimatorMode"], as_index=False
        )["ErrorAligned"]
        .agg(RaterRows="count", RaterRMSE=lambda x: float(np.sqrt(np.mean(np.square(x)))))
    )
    diagnostics, diagnostic_summary, loss_wide, slope_wide = source.summarize_preflight_diagnostics(
        rater_loss=rater_loss, rater_errors=rater_errors, runs=runs
    )

    jmle = outcomes.loc[outcomes["AttemptType"].eq(JMLE_ATTEMPT)]
    mml = outcomes.loc[outcomes["AttemptType"].str.contains("MML")]
    cmle = outcomes.loc[outcomes["AttemptType"].eq("PYTHON_EXACT_CMLE_PCM")]
    input_audit = pd.read_csv(input_dir / "assignment_audit.csv")
    input_diagnostics = pd.read_csv(input_dir / "assignment_diagnostics.csv")
    operational_gates = {
        "join_audit_pass": bool(retained_audit["pass"]),
        "completion_markers_48": len(original_hashes) == 48,
        "all_attempts_execution_completed": bool(outcomes["ExecutionCompleted"].all()),
        "all_native_estimator_evidence_ready": bool(outcomes["StatisticalEvidenceReady"].all()),
        "facets_python_calibration_12_of_12": len(jmle) == 12
        and bool(jmle["FACETSCalibrationReady"].fillna(False).astype(bool).all()),
        "mml_24_of_24_ready": len(mml) == 24
        and bool(mml["StatisticalEvidenceReady"].all()),
        "cmle_12_of_12_ready": len(cmle) == 12
        and bool(cmle["StatisticalEvidenceReady"].all()),
        "constraints_36_of_36_pass": len(constraints) == 36
        and bool(constraints["ConstraintPass"].fillna(False).astype(bool).all()),
        "assignment_and_design_audits_pass": bool(input_audit["Passed"].all()),
        "dense_assignment_draws_12_available": len(input_diagnostics) == 12
        and bool(input_diagnostics["Available"].all()),
        "registered_statistics_16_finite": len(diagnostics) == 16
        and bool(np.isfinite(diagnostics["Value"].to_numpy(dtype=float)).all()),
        "supplement_facets_rows_excluded_from_endpoints": bool(
            (~calibration_runs["IncludedInStudy"].fillna(False).astype(bool)).all()
            and (~calibration_recovery["IncludedInStudy"].fillna(False).astype(bool)).all()
            and (~calibration_thresholds["IncludedInStudy"].fillna(False).astype(bool)).all()
        ),
        "no_facets_rounded_fit_input": True,
    }
    primary_row = diagnostic_summary.loc[
        diagnostic_summary["DiagnosticId"].eq(
            "PF1_FREE_MML_SYMMETRIC_STRESS_RATER_RMSE"
        )
    ]
    primary_advancement = bool(
        len(primary_row) == 1 and primary_row["AdvancementSignal"].iloc[0]
    )
    metrics = {
        "schema_version": SCHEMA_VERSION,
        "phase": "nonconfirmatory_preflight_with_visible_facets_calibration",
        "operational_gates": {key: bool(value) for key, value in operational_gates.items()},
        "qualification_pass": bool(all(operational_gates.values())),
        "primary_advancement_signal": primary_advancement,
        "eligible_to_freeze_separate_confirmation": bool(
            all(operational_gates.values()) and primary_advancement
        ),
        "facets_calibration": {
            "pairs": int(assessment["datasets"]),
            "ready": int(sum(bool(marker["calibration_ready"]) for marker in supplement.values())),
            "maximum_main_weighted_mae": assessment["maximum_main_weighted_mae"],
            "maximum_main_absolute_difference": assessment["maximum_main_absolute_difference"],
            "maximum_threshold_weighted_mae": assessment["maximum_threshold_weighted_mae"],
            "maximum_threshold_absolute_difference": assessment["maximum_threshold_absolute_difference"],
            "minimum_within_facet_spearman": assessment["minimum_within_facet_spearman"],
            "role": "same-estimand calibration only",
        },
        "screening_observations_pooled": False,
        "preflight_observations_may_enter_confirmation": False,
        "p_values_computed": False,
        "confirmatory_claims_allowed": False,
        "estimator_ranking_constructed": False,
        "cross_basis_likelihood_comparison": "prohibited",
        "facets_fit_precision_boundary": (
            "ordinary displayed two-decimal FACETS fit fields were not calculation inputs"
        ),
        "source_artifacts_modified": False,
        "claim_limit": plan["claim_limit"],
    }
    outputs = {
        "attempt_outcomes.csv": outcomes,
        "run_ledger.csv": runs,
        "recovery.csv": recovery,
        "thresholds.csv": thresholds,
        "constraints.csv": constraints,
        "rater_errors.csv": rater_errors,
        "rater_rmse.csv": rater_loss,
        "symmetric_rmse_wide.csv": loss_wide,
        "rater_slope_wide.csv": slope_wide,
        "registered_preflight_diagnostics.csv": diagnostics,
        "registered_preflight_summary.csv": diagnostic_summary,
        "facets_calibration_run_ledger.csv": calibration_runs,
        "facets_calibration_recovery.csv": calibration_recovery,
        "facets_calibration_thresholds.csv": calibration_thresholds,
    }
    AGGREGATE_DIR.mkdir()
    for filename, frame in outputs.items():
        frame.to_csv(AGGREGATE_DIR / filename, index=False, lineterminator="\n")
    _json_dump(AGGREGATE_DIR / "assessment.json", metrics)
    artifact_names = tuple(outputs) + ("assessment.json",)
    identity_output = {
        "schema_version": f"{SCHEMA_VERSION}_identity",
        "plan_sha256": sha256_file(PLAN_PATH),
        "registration_sha256": sha256_file(REGISTRATION_PATH),
        "runner_sha256": sha256_file(Path(__file__).resolve()),
        "source_study_identity_sha256": sha256_file(STUDY_DIR / "study_identity.json"),
        "source_completion_marker_set_sha256": _hash_mapping(original_hashes),
        "supplement_marker_set_sha256": _hash_mapping(supplement_hashes),
        "join_audit_sha256": sha256_file(AUDIT_PATH),
        "join_replay_sha256": sha256_file(AUDIT_PAIRS_PATH),
        "artifact_sha256": {
            filename: sha256_file(AGGREGATE_DIR / filename) for filename in artifact_names
        },
    }
    _json_dump(AGGREGATE_DIR / "aggregate_identity.json", identity_output)
    print(json.dumps(metrics, ensure_ascii=False, sort_keys=True))
    return metrics


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=("audit", "aggregate"))
    args = parser.parse_args()
    if args.command == "audit":
        audit_join()
    else:
        aggregate_joined()


if __name__ == "__main__":
    main()
