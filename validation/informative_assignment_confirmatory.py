#!/usr/bin/env python3
"""Prepare, execute, resume, and aggregate the registered assignment confirmation."""

from __future__ import annotations

import argparse
from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path
import platform
import sys
from typing import Any, Iterable, Iterator

import numpy as np
import pandas as pd
from scipy.stats import t as student_t


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from validation import informative_assignment_screening as screening  # noqa: E402
from validation.estimand_bridge_pilot import CMLE_MODE  # noqa: E402
from validation.estimand_distribution_study import (  # noqa: E402
    FACETS_MODE,
    MML_FIXED_MODE,
    MML_FREE_MODE,
    PYTHON_JMLE_MODE,
    enrich_registered_metadata,
)
from validation.facets_resilient_pair import (  # noqa: E402
    dependency_manifest,
    dependency_manifest_digest,
)
from validation.informative_assignment_screening_aggregate_corrected import (  # noqa: E402
    build_corrected_outputs,
)
from validation.operating_characteristics_facets import sha256_file  # noqa: E402


SCHEMA_VERSION = "mfrm-informative-assignment-confirmatory-v1"
PLAN_PATH = REPO_ROOT / "validation" / "informative_assignment_confirmatory_plan_20260811.json"
REGISTRATION_PATH = (
    REPO_ROOT
    / "validation"
    / "informative_assignment_confirmatory_execution_registration_20260811.json"
)
FACETS_EXE_DEFAULT = Path(r"C:\Facets\Facets.exe")
CONFIRMATORY_DESIGNS = (
    "planned_connected",
    "ability_severity_aligned_connected",
)
SCIENTIFIC_MODES = (PYTHON_JMLE_MODE, MML_FIXED_MODE, MML_FREE_MODE, CMLE_MODE)
PRIMARY_ID = "IA1_FREE_MML_RATER_RMSE_ALIGNED_GT_PLANNED"
SECONDARY_IDS = (
    "IA2_FIXED_MML_RATER_RMSE_ALIGNED_GT_PLANNED",
    "IA3_FREE_MML_POPULATION_SD_ALIGNED_LT_PLANNED",
    "IA4_FREE_MML_RATER_COMPRESSION_SLOPE_LT_ZERO",
    "IA5_FIXED_MML_RATER_COMPRESSION_SLOPE_LT_ZERO",
)
RATER_TRUTH = {"R01": -0.45, "R02": -0.15, "R03": 0.15, "R04": 0.45}


def _json_dump(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _atomic_json(path: Path, value: Any) -> None:
    temporary = path.with_name(f"{path.name}.tmp.{os.getpid()}")
    _json_dump(temporary, value)
    os.replace(temporary, path)


def _portable_name(path: Path) -> str:
    resolved = path.resolve()
    try:
        return resolved.relative_to(REPO_ROOT).as_posix()
    except ValueError:
        return str(resolved)


def _resolve_recorded_path(value: str) -> Path:
    path = Path(value)
    return path.resolve() if path.is_absolute() else (REPO_ROOT / path).resolve()


def registered_range(plan: dict[str, Any]) -> tuple[int, int, int]:
    evidence = plan["evidence_status"]
    start, end = (int(value) for value in evidence["replicate_range"])
    count = int(evidence["replicates"])
    if (start, end, count) != (241, 340, 100) or end - start + 1 != count:
        raise ValueError("Confirmatory replicate range must be the frozen 241..340")
    if bool(evidence.get("pool_replicates_1_through_240", True)):
        raise ValueError("Confirmation must not pool replicates 1 through 240")
    if not bool(evidence.get("fixed_sample_size", False)):
        raise ValueError("Confirmation must use a fixed sample size")
    if bool(evidence.get("optional_stopping", True)):
        raise ValueError("Confirmation must prohibit optional stopping")
    if int(evidence.get("interim_scientific_looks", -1)) != 0:
        raise ValueError("Confirmation must prohibit interim scientific looks")
    return start, end, count


def build_dependency_manifest(
    *, plan_path: Path = PLAN_PATH, registration_path: Path = REGISTRATION_PATH
) -> dict[str, str]:
    manifest = dependency_manifest()
    for path in (
        Path(__file__).resolve(),
        REPO_ROOT / "validation" / "informative_assignment_design.py",
        REPO_ROOT / "validation" / "informative_assignment_screening.py",
        REPO_ROOT / "validation" / "informative_assignment_screening_aggregate_corrected.py",
        plan_path.resolve(),
        registration_path.resolve(),
    ):
        if not path.is_file():
            raise FileNotFoundError(f"Confirmatory dependency missing: {path}")
        manifest[_portable_name(path)] = sha256_file(path)
    return dict(sorted(manifest.items()))


def validate_registration(
    plan_path: Path = PLAN_PATH,
    registration_path: Path = REGISTRATION_PATH,
) -> dict[str, Any]:
    value = json.loads(registration_path.read_text(encoding="utf-8"))
    expected = {
        "plan_sha256": sha256_file(plan_path),
        "runner_sha256": sha256_file(Path(__file__).resolve()),
        "design_component_sha256": sha256_file(
            REPO_ROOT / "validation" / "informative_assignment_design.py"
        ),
        "resilient_pair_component_sha256": sha256_file(
            REPO_ROOT / "validation" / "facets_resilient_pair.py"
        ),
    }
    for key, digest in expected.items():
        if str(value.get(key, "")).lower() != digest.lower():
            raise ValueError(f"Execution registration mismatch: {key}")
    if not bool(value.get("tests_passed_before_registration", False)):
        raise ValueError("Execution registration does not assert passing tests")
    return value


@contextmanager
def _screening_engine_context() -> Iterator[None]:
    """Temporarily parameterize the frozen execution engine without mutating it."""

    replacements = {
        "DESIGNS": CONFIRMATORY_DESIGNS,
        "_registered_range": registered_range,
        "_validate_registration": validate_registration,
        "build_dependency_manifest": build_dependency_manifest,
        "validate_study_identity": validate_study_identity,
    }
    originals = {name: getattr(screening, name) for name in replacements}
    try:
        for name, value in replacements.items():
            setattr(screening, name, value)
        yield
    finally:
        for name, value in originals.items():
            setattr(screening, name, value)


def prepare_study(
    study_dir: Path,
    *,
    facets_exe: Path,
    plan_path: Path = PLAN_PATH,
    registration_path: Path = REGISTRATION_PATH,
) -> dict[str, Any]:
    study_dir = study_dir.resolve()
    if study_dir.exists():
        raise FileExistsError(f"Study directory already exists: {study_dir}")
    plan_path = plan_path.resolve()
    registration_path = registration_path.resolve()
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    start, end, count = registered_range(plan)
    facets_exe = facets_exe.resolve()
    if not facets_exe.is_file():
        raise FileNotFoundError(f"FACETS executable missing: {facets_exe}")
    validate_registration(plan_path, registration_path)
    dependencies = build_dependency_manifest(
        plan_path=plan_path, registration_path=registration_path
    )
    with _screening_engine_context():
        bundle = screening.generate_screening_bundle(
            replicates=count,
            replicate_start=start,
            base_seed=int(plan["data_generating_process"]["base_seed"]),
            plan_path=plan_path,
            registration_path=registration_path,
        )
    study_dir.mkdir(parents=True)
    input_dir = study_dir / "retained_input"
    input_dir.mkdir()
    for filename, frame in bundle.items():
        frame.to_csv(input_dir / filename, index=False, lineterminator="\n")
    input_hashes = {
        filename: sha256_file(input_dir / filename) for filename in screening.INPUT_FILES
    }
    _json_dump(input_dir / "bundle_hashes.json", input_hashes)
    _json_dump(study_dir / "dependency_manifest.json", dependencies)
    identity = {
        "schema_version": SCHEMA_VERSION,
        "phase": "registered_confirmatory",
        "replicate_start": start,
        "replicate_end": end,
        "replicates": count,
        "datasets": int(len(bundle["manifest.csv"])),
        "attempts": int(len(bundle["attempt_manifest.csv"])),
        "base_seed": int(plan["data_generating_process"]["base_seed"]),
        "plan_file": _portable_name(plan_path),
        "plan_sha256": sha256_file(plan_path),
        "registration_file": _portable_name(registration_path),
        "registration_sha256": sha256_file(registration_path),
        "runner_sha256": sha256_file(Path(__file__).resolve()),
        "dependency_manifest_sha256": dependency_manifest_digest(dependencies),
        "facets_executable": str(facets_exe),
        "facets_executable_sha256": sha256_file(facets_exe),
        "retained_input_sha256": input_hashes,
        "python_version": sys.version,
        "platform": platform.platform(),
        "prior_replicates_pooled": False,
        "confirmatory_claims_allowed_after_qualification": True,
        "claim_limit": plan["claim_limit"],
    }
    _json_dump(study_dir / "study_identity.json", identity)
    return identity


def validate_study_identity(study_dir: Path) -> dict[str, Any]:
    study_dir = study_dir.resolve()
    identity = json.loads((study_dir / "study_identity.json").read_text(encoding="utf-8"))
    plan_path = _resolve_recorded_path(identity["plan_file"])
    registration_path = _resolve_recorded_path(identity["registration_file"])
    for key, actual in {
        "plan_sha256": sha256_file(plan_path),
        "registration_sha256": sha256_file(registration_path),
        "runner_sha256": sha256_file(Path(__file__).resolve()),
    }.items():
        if str(identity.get(key, "")).lower() != actual.lower():
            raise ValueError(f"Study identity mismatch: {key}")
    validate_registration(plan_path, registration_path)
    dependencies = build_dependency_manifest(
        plan_path=plan_path, registration_path=registration_path
    )
    if identity["dependency_manifest_sha256"] != dependency_manifest_digest(dependencies):
        raise ValueError("Dependency manifest changed after preparation")
    retained = json.loads(
        (study_dir / "dependency_manifest.json").read_text(encoding="utf-8")
    )
    if retained != dependencies:
        raise ValueError("Retained dependency manifest differs from current dependencies")
    if identity["facets_executable_sha256"] != sha256_file(
        Path(identity["facets_executable"])
    ):
        raise ValueError("FACETS executable changed after preparation")
    for filename, expected in identity["retained_input_sha256"].items():
        if sha256_file(study_dir / "retained_input" / filename) != expected:
            raise ValueError(f"Retained input changed: {filename}")
    return identity


def run_shard(
    study_dir: Path,
    *,
    shard_index: int,
    shard_count: int,
    resume: bool,
    timeout_seconds: float,
) -> dict[str, int]:
    with _screening_engine_context():
        return screening.run_shard(
            study_dir,
            shard_index=shard_index,
            shard_count=shard_count,
            resume=resume,
            timeout_seconds=timeout_seconds,
        )


def preexecution_audit(study_dir: Path) -> dict[str, Any]:
    study_dir = study_dir.resolve()
    identity = validate_study_identity(study_dir)
    input_dir = study_dir / "retained_input"
    manifest = pd.read_csv(input_dir / "manifest.csv")
    ratings = pd.read_csv(input_dir / "generated_ratings.csv")
    attempts = pd.read_csv(input_dir / "attempt_manifest.csv")
    assignment = pd.read_csv(input_dir / "assignment_audit.csv")
    row_counts = ratings.groupby("RunId").size()
    expected_rows = manifest.set_index("RunId")["ExpectedRows"].astype(int)
    longest = study_dir / (
        "w/031f/g00/facets_attempts/pair_try_03/facets_runs/"
        "aligned__normal_rep-00340__pcm/report_u6.txt"
    )
    aligned = manifest[manifest["Design"].eq("ability_severity_aligned_connected")]
    checks = {
        "identity_valid": True,
        "registered_replicates_100": manifest["Replicate"].nunique() == 100,
        "replicate_range_241_340": [
            int(manifest["Replicate"].min()),
            int(manifest["Replicate"].max()),
        ] == [241, 340],
        "datasets_200": len(manifest) == 200,
        "attempts_800": len(attempts) == 800,
        "designs_exact": set(manifest["Design"].astype(str))
        == set(CONFIRMATORY_DESIGNS),
        "row_counts_match": row_counts.reindex(expected_rows.index).astype(int).equals(
            expected_rows
        ),
        "all_rows_960": manifest["ExpectedRows"].eq(960).all(),
        "all_nullity_zero": manifest["ExpectedStructuralNullity"].eq(0).all(),
        "all_one_component": manifest["PersonRaterComponents"].eq(1).all(),
        "aligned_assignment_strong": aligned["ThetaAssignedSeveritySpearman"].gt(0.9).all(),
        "complete_category_support": manifest["CompleteCriterionCategorySupport"].all(),
        "positive_category_counts": manifest["MinimumCriterionCategoryCount"].gt(0).all(),
        "shared_complete_response_within_replicate": manifest.groupby("Replicate")[
            "CompleteResponseSHA256"
        ].nunique().eq(1).all(),
        "assignment_invariants_all_pass": assignment["AuditValue"]
        .fillna(False)
        .astype(bool)
        .all(),
        "realized_mean_zero": manifest["PersonMeanRealized"].abs().max() <= 1e-12,
        "realized_sd_08": (manifest["PersonSDRealized"] - 0.8).abs().max() <= 1e-12,
        "no_screening_replicates": int(manifest["Replicate"].min()) > 240,
        "worst_case_facets_path_below_240": len(str(longest)) < 240,
        "facets_hash_bound": bool(identity["facets_executable_sha256"]),
    }
    audit = {
        "schema_version": f"{SCHEMA_VERSION}-preexecution-audit",
        "checks": {key: bool(value) for key, value in checks.items()},
        "all_checks_pass": all(bool(value) for value in checks.values()),
        "worst_case_facets_report_path_length": len(str(longest)),
        "study_identity_sha256": sha256_file(study_dir / "study_identity.json"),
    }
    _json_dump(study_dir / "preexecution_audit.json", audit)
    return audit


def _read_attempt_artifacts(
    study_dir: Path,
    identity: dict[str, Any],
) -> tuple[
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    dict[str, str],
]:
    attempts = pd.read_csv(study_dir / "retained_input" / "attempt_manifest.csv")
    parts: dict[str, list[pd.DataFrame]] = {
        "run_ledger.csv": [],
        "recovery.csv": [],
        "thresholds.csv": [],
        "constraints.csv": [],
    }
    outcomes: list[dict[str, Any]] = []
    marker_hashes: dict[str, str] = {}
    for _, attempt in attempts.iterrows():
        marker = screening._attempt_dir(study_dir, attempt) / "completion.json"
        if not marker.is_file():
            raise FileNotFoundError(f"Missing completion marker: {attempt['AttemptId']}")
        completion = json.loads(marker.read_text(encoding="utf-8"))
        artifact_root = study_dir / str(completion["retained_artifact_root"])
        screening._validate_completion(completion, attempt, artifact_root, identity)
        marker_hashes[str(attempt["AttemptId"])] = sha256_file(marker)
        for filename in parts:
            frame = screening._read_csv_if_nonempty(artifact_root / filename)
            if len(frame):
                frame.insert(0, "AttemptId", str(attempt["AttemptId"]))
                parts[filename].append(frame)
        outcome = {
            "AttemptId": str(attempt["AttemptId"]),
            "RunId": str(attempt["RunId"]),
            "AttemptType": str(attempt["AttemptType"]),
            "ExecutionCompleted": bool(completion["attempt_execution_completed"]),
            "StatisticalEvidenceReady": bool(completion["statistical_evidence_ready"]),
            "FACETSCalibrationReady": completion.get("facets_calibration_ready"),
            "FailureReason": completion.get("failure_reason", ""),
        }
        resilient_metrics = artifact_root / "resilient_pair_metrics.json"
        if resilient_metrics.is_file():
            outcome.update(json.loads(resilient_metrics.read_text(encoding="utf-8")))
        outcomes.append(outcome)
    combined = {
        filename: pd.concat(frames, ignore_index=True, sort=False)
        if frames
        else pd.DataFrame()
        for filename, frames in parts.items()
    }
    return (
        combined["run_ledger.csv"],
        combined["recovery.csv"],
        combined["thresholds.csv"],
        combined["constraints.csv"],
        pd.DataFrame(outcomes),
        marker_hashes,
    )


def _enrich_ledgers(
    runs: pd.DataFrame,
    recovery: pd.DataFrame,
    thresholds: pd.DataFrame,
    constraints: pd.DataFrame,
    manifest: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    source = {
        "runs": runs,
        "recovery": recovery,
        "thresholds": thresholds,
        "constraints": constraints,
    }
    enriched: dict[str, pd.DataFrame] = {}
    missing_before: dict[str, dict[str, int]] = {}
    for label, frame in source.items():
        missing_before[label] = {
            column: int(frame[column].isna().sum()) if column in frame else int(len(frame))
            for column in ("ConditionId", "Design", "PersonDistribution", "Replicate")
        }
        enriched[label] = enrich_registered_metadata(frame, manifest, label=label)
    audit = {
        "rows": {label: int(len(frame)) for label, frame in enriched.items()},
        "missing_before": missing_before,
        "missing_after": {
            label: {
                column: int(frame[column].isna().sum())
                for column in ("ConditionId", "Design", "PersonDistribution", "Replicate")
            }
            for label, frame in enriched.items()
        },
        "all_registered_metadata_complete": all(
            frame[["ConditionId", "Design", "PersonDistribution", "Replicate"]]
            .notna()
            .all()
            .all()
            for frame in enriched.values()
            if len(frame)
        ),
    }
    return (
        enriched["runs"],
        enriched["recovery"],
        enriched["thresholds"],
        enriched["constraints"],
        audit,
    )


def holm_adjust(p_values: Iterable[float]) -> np.ndarray:
    values = np.asarray(list(p_values), dtype=float)
    if len(values) == 0 or not np.isfinite(values).all():
        raise ValueError("Holm adjustment requires finite p values")
    order = np.argsort(values, kind="mergesort")
    adjusted_sorted = np.empty(len(values), dtype=float)
    running = 0.0
    for rank, index in enumerate(order):
        running = max(running, (len(values) - rank) * values[index])
        adjusted_sorted[rank] = min(1.0, running)
    adjusted = np.empty(len(values), dtype=float)
    adjusted[order] = adjusted_sorted
    return adjusted


def _endpoint_summary(
    endpoint_id: str,
    values: pd.Series,
    *,
    alternative: str,
    precision_target: float | None = None,
) -> dict[str, Any]:
    numeric = pd.to_numeric(values, errors="coerce").dropna().to_numpy(dtype=float)
    n = len(numeric)
    mean = float(np.mean(numeric)) if n else np.nan
    sd = float(np.std(numeric, ddof=1)) if n > 1 else np.nan
    se = sd / np.sqrt(n) if n > 1 else np.nan
    statistic = mean / se if np.isfinite(se) and se > 0 else np.nan
    df = n - 1
    raw_p = (
        float(student_t.cdf(statistic, df))
        if alternative == "less" and np.isfinite(statistic)
        else float(student_t.sf(statistic, df))
        if alternative == "greater" and np.isfinite(statistic)
        else np.nan
    )
    half = float(student_t.ppf(0.975, df)) * se if np.isfinite(se) else np.nan
    direction_pass = bool(
        np.isfinite(mean)
        and ((alternative == "greater" and mean > 0) or (alternative == "less" and mean < 0))
    )
    return {
        "EndpointId": endpoint_id,
        "Alternative": alternative,
        "FinitePairedReplicates": n,
        "RequiredPairedReplicates": 100,
        "FullPairGate": n == 100,
        "MeanContrast": mean,
        "MonteCarloSD": sd,
        "MonteCarloSE": se,
        "TStatistic": statistic,
        "DegreesOfFreedom": df,
        "RawOneSidedP": raw_p,
        "Lower95": mean - half if np.isfinite(half) else np.nan,
        "Upper95": mean + half if np.isfinite(half) else np.nan,
        "TwoSided95HalfWidth": half,
        "PrecisionTarget": precision_target,
        "PrecisionPass": bool(
            precision_target is not None
            and np.isfinite(half)
            and half <= precision_target
        ),
        "DirectionPass": direction_pass,
    }


def build_endpoint_contrasts(
    paired: pd.DataFrame,
    rater_level: pd.DataFrame,
    free_sd: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    for endpoint_id, mode in (
        (PRIMARY_ID, MML_FREE_MODE),
        (SECONDARY_IDS[0], MML_FIXED_MODE),
    ):
        selected = paired[
            paired["EstimatorMode"].eq(mode)
            & paired["RecoveryDomain"].eq("Facet:Rater")
            & paired["Metric"].eq("RMSE")
        ][["Replicate", "ContrastAlignedMinusPlanned"]].copy()
        selected.insert(0, "EndpointId", endpoint_id)
        selected = selected.rename(columns={"ContrastAlignedMinusPlanned": "Contrast"})
        rows.append(selected)
    sd_selected = free_sd[["Replicate", "ContrastAlignedMinusPlanned"]].copy()
    sd_selected.insert(0, "EndpointId", SECONDARY_IDS[1])
    sd_selected = sd_selected.rename(columns={"ContrastAlignedMinusPlanned": "Contrast"})
    rows.append(sd_selected)
    truth = pd.Series(RATER_TRUTH, dtype=float)
    denominator = float(np.square(truth.to_numpy()).sum())
    for endpoint_id, mode in (
        (SECONDARY_IDS[2], MML_FREE_MODE),
        (SECONDARY_IDS[3], MML_FIXED_MODE),
    ):
        selected = rater_level[rater_level["EstimatorMode"].eq(mode)].copy()
        selected["TrueSeverity"] = selected["Level"].map(truth)
        slope = (
            selected.assign(
                Product=lambda frame: frame["TrueSeverity"]
                * frame["ErrorContrastAlignedMinusPlanned"]
            )
            .groupby("Replicate", as_index=False)["Product"]
            .sum()
        )
        slope["Contrast"] = slope.pop("Product") / denominator
        slope.insert(0, "EndpointId", endpoint_id)
        rows.append(slope)
    output = pd.concat(rows, ignore_index=True)
    observed = tuple(output["EndpointId"].drop_duplicates())
    if observed != (PRIMARY_ID, *SECONDARY_IDS):
        raise RuntimeError("Confirmatory endpoint construction order changed")
    return output


def summarize_endpoints(contrasts: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    primary_values = contrasts.loc[contrasts["EndpointId"].eq(PRIMARY_ID), "Contrast"]
    primary = pd.DataFrame([
        _endpoint_summary(
            PRIMARY_ID,
            primary_values,
            alternative="greater",
            precision_target=0.015,
        )
    ])
    primary["AdjustedP"] = primary["RawOneSidedP"]
    primary["DirectionConfirmed"] = (
        primary["FullPairGate"].astype(bool)
        & primary["DirectionPass"].astype(bool)
        & primary["RawOneSidedP"].le(0.05)
    )
    alternatives = {
        SECONDARY_IDS[0]: "greater",
        SECONDARY_IDS[1]: "less",
        SECONDARY_IDS[2]: "less",
        SECONDARY_IDS[3]: "less",
    }
    secondary = pd.DataFrame([
        _endpoint_summary(
            endpoint_id,
            contrasts.loc[contrasts["EndpointId"].eq(endpoint_id), "Contrast"],
            alternative=alternatives[endpoint_id],
        )
        for endpoint_id in SECONDARY_IDS
    ])
    secondary["HolmAdjustedP"] = holm_adjust(secondary["RawOneSidedP"])
    primary_pass = bool(primary.iloc[0]["DirectionConfirmed"])
    secondary["PrimaryGatePass"] = primary_pass
    secondary["DirectionConfirmed"] = (
        secondary["PrimaryGatePass"].astype(bool)
        & secondary["FullPairGate"].astype(bool)
        & secondary["DirectionPass"].astype(bool)
        & secondary["HolmAdjustedP"].le(0.05)
    )
    return primary, secondary


def aggregate_study(study_dir: Path) -> dict[str, Any]:
    study_dir = study_dir.resolve()
    identity = validate_study_identity(study_dir)
    manifest = pd.read_csv(study_dir / "retained_input" / "manifest.csv")
    runs, recovery, thresholds, constraints, outcomes, marker_hashes = (
        _read_attempt_artifacts(study_dir, identity)
    )
    runs, recovery, thresholds, constraints, join_audit = _enrich_ledgers(
        runs, recovery, thresholds, constraints, manifest
    )
    jmle = outcomes[outcomes["AttemptType"].eq("RESILIENT_FACETS_PYTHON_JMLE_PCM")]
    mml = outcomes[outcomes["AttemptType"].str.contains("MML")]
    cmle = outcomes[outcomes["AttemptType"].eq("PYTHON_EXACT_CMLE_PCM")]
    calibration_ready = jmle["FACETSCalibrationReady"].fillna(False).astype(bool)
    pre_gates = {
        "completion_markers_800": len(marker_hashes) == 800,
        "all_attempts_execution_completed": outcomes["ExecutionCompleted"].all(),
        "all_statistical_evidence_ready": outcomes["StatisticalEvidenceReady"].all(),
        "facets_calibration_200_of_200": len(jmle) == 200 and calibration_ready.all(),
        "mml_400_of_400": len(mml) == 400 and mml["StatisticalEvidenceReady"].all(),
        "cmle_200_of_200": len(cmle) == 200 and cmle["StatisticalEvidenceReady"].all(),
        "constraints_600_of_600": len(constraints) == 600
        and constraints["ConstraintPass"].fillna(False).astype(bool).all(),
        "metadata_join_complete": bool(join_audit["all_registered_metadata_complete"]),
    }
    if not all(bool(value) for value in pre_gates.values()):
        raise RuntimeError(f"Operational qualification failed before endpoint analysis: {pre_gates}")
    scientific = build_corrected_outputs(runs, recovery, thresholds)
    contrasts = build_endpoint_contrasts(
        scientific["paired_design_contrasts.csv"],
        scientific["rater_level_contrasts.csv"],
        scientific["free_sd_paired.csv"],
    )
    primary, secondary = summarize_endpoints(contrasts)
    endpoint_gate = bool(
        len(contrasts) == 500
        and contrasts.groupby("EndpointId").size().eq(100).all()
        and primary["FullPairGate"].all()
        and secondary["FullPairGate"].all()
    )
    gates = {**pre_gates, "five_endpoints_100_pairs_each": endpoint_gate}
    numeric = lambda column: pd.to_numeric(jmle[column], errors="coerce")
    facets_runs = runs[runs["EstimatorMode"].astype(str).eq(FACETS_MODE)].copy()
    facets_numeric = lambda column: pd.to_numeric(facets_runs[column], errors="coerce")
    calibration = {
        "pairs": int(len(jmle)),
        "ready": int(calibration_ready.sum()),
        "facets_total_tries": int(numeric("facets_tries").fillna(0).sum()),
        "facets_total_retries": int(numeric("facets_retry_count").fillna(0).sum()),
        "maximum_main_weighted_mae": float(facets_numeric("MainWeightedMAE").max()),
        "maximum_main_absolute_difference": float(
            facets_numeric("MainMaxAbsDifference").max()
        ),
        "maximum_threshold_weighted_mae": float(
            facets_numeric("ThresholdWeightedMAE").max()
        ),
        "maximum_threshold_absolute_difference": float(
            facets_numeric("ThresholdMaxAbsDifference").max()
        ),
        "minimum_within_facet_spearman": float(
            facets_numeric("MinimumWithinFacetSpearman").min()
        ),
    }
    marker_digest = hashlib.sha256(
        json.dumps(marker_hashes, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    aggregate_dir = study_dir / "aggregate"
    if aggregate_dir.exists():
        raise FileExistsError(f"Aggregate directory already exists: {aggregate_dir}")
    aggregate_dir.mkdir()
    outputs: dict[str, pd.DataFrame] = {
        "attempt_outcomes.csv": outcomes,
        "run_ledger.csv": runs,
        "recovery.csv": recovery,
        "thresholds.csv": thresholds,
        "constraints.csv": constraints,
        **scientific,
        "endpoint_contrasts.csv": contrasts,
        "primary_result.csv": primary,
        "secondary_results.csv": secondary,
    }
    for filename, frame in outputs.items():
        frame.to_csv(aggregate_dir / filename, index=False, lineterminator="\n")
    _json_dump(aggregate_dir / "metadata_join_audit.json", join_audit)
    metrics = {
        "schema_version": SCHEMA_VERSION,
        "phase": "confirmatory",
        "datasets": int(identity["datasets"]),
        "attempts": int(identity["attempts"]),
        "gates": {key: bool(value) for key, value in gates.items()},
        "qualification_pass": all(bool(value) for value in gates.values()),
        "facets_calibration": calibration,
        "completion_marker_set_sha256": marker_digest,
        "primary_endpoint_id": PRIMARY_ID,
        "primary_direction_confirmed": bool(primary.iloc[0]["DirectionConfirmed"]),
        "primary_precision_pass": bool(primary.iloc[0]["PrecisionPass"]),
        "secondary_family_method": "Holm",
        "estimator_ranking_constructed": False,
        "cross_basis_likelihood_comparison": "prohibited",
        "claim_limit": identity["claim_limit"],
    }
    _json_dump(aggregate_dir / "confirmatory_metrics.json", metrics)
    artifact_names = tuple(outputs) + (
        "metadata_join_audit.json",
        "confirmatory_metrics.json",
    )
    _atomic_json(aggregate_dir / "aggregate_identity.json", {
        "schema_version": f"{SCHEMA_VERSION}-aggregate-identity",
        "study_identity_sha256": sha256_file(study_dir / "study_identity.json"),
        "completion_marker_set_sha256": marker_digest,
        "artifact_sha256": {
            filename: sha256_file(aggregate_dir / filename) for filename in artifact_names
        },
    })
    return metrics


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    prepare = subparsers.add_parser("prepare")
    prepare.add_argument("--study-dir", type=Path, required=True)
    prepare.add_argument("--facets-exe", type=Path, default=FACETS_EXE_DEFAULT)
    audit = subparsers.add_parser("audit")
    audit.add_argument("--study-dir", type=Path, required=True)
    run = subparsers.add_parser("run")
    run.add_argument("--study-dir", type=Path, required=True)
    run.add_argument("--shard-index", type=int, required=True)
    run.add_argument("--shard-count", type=int, required=True)
    run.add_argument("--resume", action="store_true")
    run.add_argument("--timeout-seconds", type=float, default=120.0)
    aggregate = subparsers.add_parser("aggregate")
    aggregate.add_argument("--study-dir", type=Path, required=True)
    return parser.parse_args(argv)


def main() -> None:
    args = parse_args()
    if args.command == "prepare":
        result = prepare_study(args.study_dir, facets_exe=args.facets_exe)
    elif args.command == "audit":
        result = preexecution_audit(args.study_dir)
    elif args.command == "run":
        result = run_shard(
            args.study_dir,
            shard_index=args.shard_index,
            shard_count=args.shard_count,
            resume=args.resume,
            timeout_seconds=args.timeout_seconds,
        )
    else:
        result = aggregate_study(args.study_dir)
    print(json.dumps(result, indent=2, sort_keys=True, default=str))


if __name__ == "__main__":
    main()
