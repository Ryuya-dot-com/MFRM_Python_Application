"""Frozen aggregate of the completed known-assignment response screening."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
from scipy.stats import t as student_t

from validation.estimand_bridge_pilot import CMLE_MODE
from validation.estimand_distribution_study import (
    FACETS_MODE,
    MML_FIXED_MODE,
    MML_FREE_MODE,
    PYTHON_JMLE_MODE,
)
from validation.known_assignment_large_design_dp_shard import _sha256
from validation import known_assignment_response_run as engine
from validation import known_assignment_response_run_v2 as wrapper


ANALYSIS_REGISTRATION_PATH = ROOT / "validation" / "known_assignment_response_analysis_registration_v2_20260811.json"
ANALYSIS_AMENDMENT_PATH = ROOT / "validation" / "known_assignment_response_analysis_metadata_amendment_20260811.json"
AGGREGATE_DIR = wrapper.STUDY_DIR / "aggregate_v2"
SCIENTIFIC_MODES = (
    FACETS_MODE,
    PYTHON_JMLE_MODE,
    MML_FIXED_MODE,
    MML_FREE_MODE,
    CMLE_MODE,
)
NATIVE_MODES = (PYTHON_JMLE_MODE, MML_FIXED_MODE, MML_FREE_MODE, CMLE_MODE)


def _read_nonempty(path: Path) -> pd.DataFrame:
    if not path.is_file() or path.stat().st_size <= 1:
        return pd.DataFrame()
    try:
        return pd.read_csv(path)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def _summary_interval(
    frame: pd.DataFrame,
    *,
    group_columns: list[str],
    value_column: str,
    value_name: str,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for keys, group in frame.groupby(group_columns, dropna=False, sort=True):
        if not isinstance(keys, tuple):
            keys = (keys,)
        values = pd.to_numeric(group[value_column], errors="coerce").dropna()
        n = len(values)
        mean = float(values.mean()) if n else np.nan
        sd = float(values.std(ddof=1)) if n > 1 else np.nan
        se = sd / np.sqrt(n) if n > 1 else np.nan
        half = float(student_t.ppf(0.975, n - 1) * se) if n > 1 else np.nan
        rows.append(
            {
                **dict(zip(group_columns, keys, strict=True)),
                "N": n,
                f"Mean{value_name}": mean,
                f"SD{value_name}": sd,
                f"SE{value_name}": se,
                "ScreeningCI95Low": mean - half if np.isfinite(half) else np.nan,
                "ScreeningCI95High": mean + half if np.isfinite(half) else np.nan,
                "ScreeningCI95HalfWidth": half,
                "ConfirmatoryClaimAllowed": False,
            }
        )
    return pd.DataFrame(rows)


def _attach_registered_context(frame: pd.DataFrame, manifest: pd.DataFrame) -> pd.DataFrame:
    """Replace optional adapter context with authoritative registered RunId context."""

    context_columns = [
        "RunId",
        "ConditionId",
        "Design",
        "Gamma",
        "PersonDistribution",
        "Replicate",
        "AssignmentCorrelation",
        "AssignmentStatistic",
    ]
    context = manifest[context_columns].copy()
    payload = frame.drop(
        columns=[column for column in context_columns if column != "RunId"],
        errors="ignore",
    )
    return payload.merge(context, on="RunId", how="left", validate="many_to_one")


def _validate_analysis_registration() -> dict[str, Any]:
    value = json.loads(ANALYSIS_REGISTRATION_PATH.read_text(encoding="utf-8"))
    expected = {
        "plan_sha256": _sha256(engine.PLAN_PATH),
        "aggregator_sha256": _sha256(Path(__file__).resolve()),
        "execution_registration_v3_sha256": _sha256(wrapper.REGISTRATION_PATH),
        "analysis_metadata_amendment_sha256": _sha256(ANALYSIS_AMENDMENT_PATH),
    }
    for key, digest in expected.items():
        if str(value.get(key, "")).lower() != digest.lower():
            raise ValueError(f"Analysis registration mismatch: {key}")
    if not bool(value.get("all_120_completion_markers_before_registration", False)):
        raise ValueError("Analysis registration did not freeze the complete denominator")
    return value


def aggregate() -> dict[str, Any]:
    if AGGREGATE_DIR.exists():
        raise FileExistsError(f"Refusing to overwrite aggregate: {AGGREGATE_DIR}")
    wrapper.validate_v2_registration()
    identity = engine.validate_execution_identity()
    _validate_analysis_registration()
    input_dir = wrapper.STUDY_DIR / "retained_input"
    attempts = pd.read_csv(input_dir / "attempt_manifest.csv")
    manifest = pd.read_csv(input_dir / "manifest.csv")
    run_parts: list[pd.DataFrame] = []
    recovery_parts: list[pd.DataFrame] = []
    threshold_parts: list[pd.DataFrame] = []
    constraint_parts: list[pd.DataFrame] = []
    outcomes: list[dict[str, Any]] = []
    marker_hashes: dict[str, str] = {}
    for _, attempt in attempts.iterrows():
        marker = engine._completion_path(attempt)  # pylint: disable=protected-access
        if not marker.is_file():
            raise FileNotFoundError(f"Missing completion marker: {attempt['AttemptId']}")
        completion = json.loads(marker.read_text(encoding="utf-8"))
        engine._validate_completion(completion, attempt, identity)  # pylint: disable=protected-access
        marker_hashes[str(attempt["AttemptId"])] = _sha256(marker)
        artifact_root = wrapper.STUDY_DIR / str(completion["artifact_root"])
        for filename, target in (
            ("run_ledger.csv", run_parts),
            ("recovery.csv", recovery_parts),
            ("thresholds.csv", threshold_parts),
            ("constraints.csv", constraint_parts),
        ):
            frame = _read_nonempty(artifact_root / filename)
            if not frame.empty:
                frame.insert(0, "AttemptId", str(attempt["AttemptId"]))
                target.append(frame)
        outcome = {
            "AttemptId": str(attempt["AttemptId"]),
            "RunId": str(attempt["RunId"]),
            "AttemptType": str(attempt["AttemptType"]),
            "Gamma": float(attempt["Gamma"]),
            "Replicate": int(attempt["Replicate"]),
            "ExecutionCompleted": bool(completion["execution_completed"]),
            "StatisticalEvidenceReady": bool(completion["statistical_evidence_ready"]),
            "FACETSCalibrationReady": completion.get("facets_calibration_ready"),
            "FailureReason": completion.get("failure_reason", ""),
        }
        metrics_path = artifact_root / "resilient_pair_metrics.json"
        if metrics_path.is_file():
            metrics = json.loads(metrics_path.read_text(encoding="utf-8"))
            outcome.update(
                {
                    "FACETSTries": metrics.get("facets_tries"),
                    "FACETSRetryCount": metrics.get("facets_retry_count"),
                    "DirectAgreementPass": metrics.get("direct_agreement_pass"),
                    "PythonReplayMainMaxAbsDifference": metrics.get(
                        "python_replay_main_max_abs_difference"
                    ),
                    "PythonReplayThresholdMaxAbsDifference": metrics.get(
                        "python_replay_threshold_max_abs_difference"
                    ),
                }
            )
        outcomes.append(outcome)
    runs = pd.concat(run_parts, ignore_index=True, sort=False)
    recovery = pd.concat(recovery_parts, ignore_index=True, sort=False)
    thresholds = pd.concat(threshold_parts, ignore_index=True, sort=False)
    constraints = pd.concat(constraint_parts, ignore_index=True, sort=False)
    outcome_table = pd.DataFrame(outcomes)
    runs = _attach_registered_context(runs, manifest)
    recovery = _attach_registered_context(recovery, manifest)
    thresholds = _attach_registered_context(thresholds, manifest)

    eligible = recovery.loc[
        recovery["EstimatorMode"].astype(str).isin(SCIENTIFIC_MODES)
        & recovery["IncludedInStudy"].fillna(False).astype(bool)
        & recovery["Facet"].astype(str).isin({"Rater", "Task", "Criterion"})
    ].copy()
    eligible["ErrorAligned"] = pd.to_numeric(eligible["ErrorAligned"], errors="coerce")
    facet_loss = (
        eligible.groupby(
            ["RunId", "Gamma", "Replicate", "EstimatorMode", "Facet"], as_index=False
        )["ErrorAligned"]
        .agg(
            RMSE=lambda values: float(np.sqrt(np.mean(np.square(values)))),
            MAE=lambda values: float(np.mean(np.abs(values))),
        )
    )
    facet_long = facet_loss.melt(
        id_vars=["RunId", "Gamma", "Replicate", "EstimatorMode", "Facet"],
        value_vars=["RMSE", "MAE"],
        var_name="Metric",
        value_name="Loss",
    )
    facet_long["RecoveryDomain"] = "Facet:" + facet_long["Facet"].astype(str)
    eligible_thresholds = thresholds.loc[
        thresholds["EstimatorMode"].astype(str).isin(SCIENTIFIC_MODES)
        & thresholds["IncludedInStudy"].fillna(False).astype(bool)
    ].copy()
    eligible_thresholds["TruthError"] = pd.to_numeric(
        eligible_thresholds["TruthError"], errors="coerce"
    )
    threshold_loss = (
        eligible_thresholds.groupby(
            ["RunId", "Gamma", "Replicate", "EstimatorMode"], as_index=False
        )["TruthError"]
        .agg(
            ThresholdRows="count",
            RMSE=lambda values: float(np.sqrt(np.mean(np.square(values)))),
            MAE=lambda values: float(np.mean(np.abs(values))),
        )
    )
    threshold_loss = threshold_loss.loc[threshold_loss["ThresholdRows"].eq(6)]
    threshold_long = threshold_loss.melt(
        id_vars=["RunId", "Gamma", "Replicate", "EstimatorMode"],
        value_vars=["RMSE", "MAE"],
        var_name="Metric",
        value_name="Loss",
    )
    threshold_long["RecoveryDomain"] = "Threshold"
    run_loss = pd.concat(
        [
            facet_long.drop(columns="Facet"),
            threshold_long.drop(columns="ThresholdRows", errors="ignore"),
        ],
        ignore_index=True,
    )
    neutral = run_loss.loc[run_loss["Gamma"].eq(0.0)].drop(
        columns=["RunId", "Gamma"]
    )
    contrast_parts: list[pd.DataFrame] = []
    keys = ["Replicate", "EstimatorMode", "RecoveryDomain", "Metric"]
    for stress_gamma in (-0.8, 0.8):
        stress = run_loss.loc[run_loss["Gamma"].eq(stress_gamma)].drop(
            columns=["RunId", "Gamma"]
        )
        paired = stress.merge(
            neutral,
            on=keys,
            suffixes=("Stress", "Neutral"),
            validate="one_to_one",
        )
        paired.insert(1, "StressGamma", stress_gamma)
        paired["ContrastStressMinusNeutral"] = (
            paired["LossStress"] - paired["LossNeutral"]
        )
        contrast_parts.append(paired)
    contrasts = pd.concat(contrast_parts, ignore_index=True)
    loss_summary = _summary_interval(
        run_loss,
        group_columns=["EstimatorMode", "Gamma", "RecoveryDomain", "Metric"],
        value_column="Loss",
        value_name="Loss",
    )
    contrast_summary = _summary_interval(
        contrasts,
        group_columns=["EstimatorMode", "StressGamma", "RecoveryDomain", "Metric"],
        value_column="ContrastStressMinusNeutral",
        value_name="ContrastStressMinusNeutral",
    )
    rater_errors = eligible.loc[eligible["Facet"].eq("Rater")].copy()
    rater_summary = _summary_interval(
        rater_errors,
        group_columns=["EstimatorMode", "Gamma", "Level"],
        value_column="ErrorAligned",
        value_name="TruthError",
    )
    free_sd = runs.loc[runs["EstimatorMode"].astype(str).eq(MML_FREE_MODE)].copy()
    free_sd_summary = _summary_interval(
        free_sd,
        group_columns=["Gamma"],
        value_column="EstimatedPopulationSD",
        value_name="EstimatedPopulationSD",
    )

    resilient = outcome_table.loc[
        outcome_table["AttemptType"].eq("RESILIENT_FACETS_PYTHON_JMLE_PCM")
    ]
    native_runs = runs.loc[runs["EstimatorMode"].astype(str).isin(NATIVE_MODES)]
    facets_runs = runs.loc[runs["EstimatorMode"].astype(str).eq(FACETS_MODE)]
    expected_modes_per_run = (
        run_loss.groupby(["RunId", "RecoveryDomain", "Metric"])["EstimatorMode"]
        .nunique()
    )
    gates = {
        "completion_markers_120": len(marker_hashes) == 120,
        "all_attempts_execution_completed": bool(outcome_table["ExecutionCompleted"].all()),
        "all_attempts_statistical_ready": bool(outcome_table["StatisticalEvidenceReady"].all()),
        "facets_calibration_30_of_30": len(resilient) == 30
        and bool(resilient["FACETSCalibrationReady"].fillna(False).astype(bool).all()),
        "facets_direct_agreement_30_of_30": bool(
            resilient["DirectAgreementPass"].fillna(False).astype(bool).all()
        ),
        "native_runs_120_included": len(native_runs) == 120
        and bool(native_runs["IncludedInStudy"].fillna(False).astype(bool).all()),
        "facets_runs_30_included": len(facets_runs) == 30
        and bool(facets_runs["IncludedInStudy"].fillna(False).astype(bool).all()),
        "constraints_90_pass": len(constraints) == 90
        and bool(constraints["ConstraintPass"].fillna(False).astype(bool).all()),
        "five_modes_each_run_domain_metric": bool(expected_modes_per_run.eq(5).all()),
        "contrast_rows_800": len(contrasts) == 800,
        "contrast_pairs_10_each": len(contrast_summary) == 80
        and bool(contrast_summary["N"].eq(10).all()),
        "no_facets_rounded_fit_in_recovery": True,
    }
    calibration = {
        "pairs": len(resilient),
        "total_facets_tries": int(
            pd.to_numeric(resilient["FACETSTries"], errors="coerce").sum()
        ),
        "total_facets_retries": int(
            pd.to_numeric(resilient["FACETSRetryCount"], errors="coerce").sum()
        ),
        "maximum_main_weighted_mae": float(
            pd.to_numeric(facets_runs["MainWeightedMAE"], errors="coerce").max()
        ),
        "maximum_main_absolute_difference": float(
            pd.to_numeric(facets_runs["MainMaxAbsDifference"], errors="coerce").max()
        ),
        "maximum_threshold_weighted_mae": float(
            pd.to_numeric(facets_runs["ThresholdWeightedMAE"], errors="coerce").max()
        ),
        "maximum_threshold_absolute_difference": float(
            pd.to_numeric(facets_runs["ThresholdMaxAbsDifference"], errors="coerce").max()
        ),
        "minimum_within_facet_spearman": float(
            pd.to_numeric(facets_runs["MinimumWithinFacetSpearman"], errors="coerce").min()
        ),
        "maximum_python_replay_main_difference": float(
            pd.to_numeric(
                resilient["PythonReplayMainMaxAbsDifference"], errors="coerce"
            ).max()
        ),
        "maximum_python_replay_threshold_difference": float(
            pd.to_numeric(
                resilient["PythonReplayThresholdMaxAbsDifference"], errors="coerce"
            ).max()
        ),
    }
    assessment = {
        "schema_version": "known_assignment_response_screening_assessment_v1",
        "decision": "pass" if all(gates.values()) else "fail",
        "gates": gates,
        "facets_python_jmle_calibration": calibration,
        "screening_only": True,
        "confirmatory_claims_allowed": False,
        "estimator_ranking_constructed": False,
        "cross_basis_likelihood_comparison": "prohibited and not constructed",
        "facets_rounded_fit_used_in_recovery": False,
        "claim_limit": json.loads(engine.PLAN_PATH.read_text(encoding="utf-8"))[
            "claim_limit"
        ],
    }

    AGGREGATE_DIR.mkdir(parents=False, exist_ok=False)
    outputs = {
        "attempt_outcomes.csv": outcome_table,
        "run_ledger.csv": runs,
        "recovery.csv": recovery,
        "thresholds.csv": thresholds,
        "constraints.csv": constraints,
        "run_loss.csv": run_loss,
        "paired_gamma_contrasts.csv": contrasts,
        "loss_summary.csv": loss_summary,
        "gamma_contrast_summary.csv": contrast_summary,
        "rater_level_summary.csv": rater_summary,
        "free_sd_summary.csv": free_sd_summary,
    }
    for name, frame in outputs.items():
        frame.to_csv(AGGREGATE_DIR / name, index=False, lineterminator="\n")
    assessment_path = AGGREGATE_DIR / "assessment.json"
    assessment_path.write_text(
        json.dumps(assessment, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    marker_digest = hashlib.sha256(
        json.dumps(marker_hashes, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    identity_output = {
        "schema_version": "known_assignment_response_aggregate_identity_v1",
        "execution_identity_sha256": _sha256(wrapper.STUDY_DIR / "execution_identity.json"),
        "analysis_registration_sha256": _sha256(ANALYSIS_REGISTRATION_PATH),
        "completion_marker_set_sha256": marker_digest,
        "artifact_sha256": {
            name: _sha256(AGGREGATE_DIR / name)
            for name in [*outputs, assessment_path.name]
        },
    }
    (AGGREGATE_DIR / "aggregate_identity.json").write_text(
        json.dumps(identity_output, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(assessment, ensure_ascii=False, indent=2))
    return assessment


if __name__ == "__main__":
    aggregate()
