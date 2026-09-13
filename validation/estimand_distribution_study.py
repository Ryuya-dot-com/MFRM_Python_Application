#!/usr/bin/env python3
"""Prepare, shard, resume, and aggregate the Person-shape estimand study.

The study deliberately fits the correct heterogeneous-threshold PCM only.
FACETS/Python JMLE parity, normal-population MML sensitivity, and exact CMLE
behavior are separate evidence layers; cross-likelihood comparisons are never
constructed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import sys
import time
from typing import Any, Iterable
import warnings

import numpy as np
import pandas as pd
from scipy.stats import kurtosis, skew


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from mfrm_app.operating_characteristics import deterministic_replicate_seed  # noqa: E402
from validation.estimand_bridge_pilot import (  # noqa: E402
    CMLE_MODE,
    CONSTRAINT_TOLERANCE,
    _normalize_recovery,
    _normalize_steps,
    fit_exact_cmle,
)
from validation.facets_pcm_boundary_pilot import (  # noqa: E402
    _fit_facets_and_python,
    apply_observation_design,
    constrained_adjacent_design_audit,
)
from validation.facets_pcm_known_truth_smoke import (  # noqa: E402
    CATEGORIES,
    CRITERION_TRUTH,
    N_PERSONS,
    RATER_TRUTH,
    TASK_TRUTH,
    THRESHOLD_CONDITIONS,
    _frame_sha256,
    apply_threshold_condition,
)
from validation.operating_characteristics_facets import (  # noqa: E402
    sha256_file,
    validate_bundle,
)
from validation.operating_characteristics_python_mml import build_mml_kwargs  # noqa: E402


SCHEMA_VERSION = "mfrm-estimand-distribution-study-v1"
PLAN_PATH = REPO_ROOT / "validation" / "estimand_distribution_screening_plan_20260811.json"
AMENDMENT_PATH = REPO_ROOT / "validation" / "estimand_distribution_preflight_infrastructure_amendment2_20260811.json"
AGGREGATION_AMENDMENT_PATH = REPO_ROOT / "validation" / "estimand_distribution_preflight_aggregation_amendment_20260811.json"
BASE_SEED = 2026081107
PERSON_SD = 0.8
FIT_MODEL = "PCM"
THRESHOLD_CONDITION = "heterogeneous"
DISTRIBUTIONS = ("normal", "right_skew", "symmetric_mixture", "heavy_tail_t3")
DESIGNS = ("complete", "planned_connected")
ATTEMPT_TYPES = (
    "FACETS_PYTHON_JMLE_PCM",
    "PYTHON_MML_FIXED_SD08_Q31_PCM",
    "PYTHON_MML_FREE_SD_Q31_PCM",
    "PYTHON_EXACT_CMLE_PCM",
)
MML_FIXED_MODE = "PYTHON_MML_FIXED_SD08_Q31"
MML_FREE_MODE = "PYTHON_MML_FREE_SD_Q31"
FACETS_MODE = "FACETS_4_5_JMLE"
PYTHON_JMLE_MODE = "PYTHON_JMLE"
INPUT_FILES = (
    "manifest.csv",
    "generated_ratings.csv",
    "generated_facet_truth.csv",
    "generated_anchors.csv",
    "generated_pcm_threshold_truth.csv",
    "attempt_manifest.csv",
)


def _json_dump(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _frame_digest(frame: pd.DataFrame) -> str:
    return hashlib.sha256(frame.to_csv(index=False, lineterminator="\n").encode("utf-8")).hexdigest()


def _identity_file_name(path: Path) -> str:
    """Store repository files portably while allowing explicit external plans."""

    resolved = path.resolve()
    try:
        return resolved.relative_to(REPO_ROOT).as_posix()
    except ValueError:
        return str(resolved)


def _identity_file_path(value: Any, *, fallback: Path) -> Path:
    """Resolve a file recorded in a study identity, including legacy identities."""

    if value in (None, ""):
        return fallback.resolve()
    path = Path(str(value))
    return path.resolve() if path.is_absolute() else (REPO_ROOT / path).resolve()


def _validate_confirmatory_prepare_contract(
    plan_path: Path,
    *,
    replicates: int,
    replicate_start: int,
) -> None:
    """Prevent a confirmatory label from being attached to an unregistered range."""

    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    evidence = plan.get("evidence_status", {})
    registered_range = evidence.get("confirmatory_replicates")
    if not isinstance(registered_range, list) or len(registered_range) != 2:
        raise ValueError("Confirmatory plan must register a two-element replicate range")
    registered_start, registered_end = (int(value) for value in registered_range)
    registered_count = int(evidence.get("confirmatory_replicate_count", -1))
    actual_end = int(replicate_start) + int(replicates) - 1
    if (
        int(replicates) != registered_count
        or int(replicate_start) != registered_start
        or actual_end != registered_end
    ):
        raise ValueError("Prepared replicate range does not match confirmatory plan")
    if not bool(evidence.get("fresh_data_required", False)):
        raise ValueError("Confirmatory plan must require fresh data")
    screening_range = evidence.get("screening_replicates", [])
    if isinstance(screening_range, list) and len(screening_range) == 2:
        screening_start, screening_end = (int(value) for value in screening_range)
        if max(screening_start, registered_start) <= min(screening_end, registered_end):
            raise ValueError("Screening and confirmatory replicate ranges overlap")


def _standardize_persons(values: np.ndarray) -> np.ndarray:
    centered = np.asarray(values, dtype=float) - float(np.mean(values))
    scale = float(np.std(centered, ddof=0))
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("Person draw has no finite dispersion")
    output = centered * (PERSON_SD / scale)
    if abs(float(output.mean())) > 1e-12 or abs(float(output.std(ddof=0)) - PERSON_SD) > 1e-12:
        raise RuntimeError("Person standardization contract failed")
    return output


def generate_person_shapes(seed: int, n_persons: int = N_PERSONS) -> dict[str, np.ndarray]:
    """Generate paired shapes and force identical realized mean/SD."""

    rng = np.random.default_rng(int(seed))
    z1 = rng.normal(size=n_persons)
    z2 = rng.normal(size=n_persons)
    component = rng.random(n_persons) < 0.5
    draws = {
        "normal": z1,
        "right_skew": np.exp(0.75 * z1),
        "symmetric_mixture": np.where(component, -1.3 + 0.45 * z1, 1.3 + 0.45 * z2),
        "heavy_tail_t3": rng.standard_t(df=3, size=n_persons),
    }
    return {name: _standardize_persons(draws[name]) for name in DISTRIBUTIONS}


def _latent_rows(person_values: np.ndarray, row_uniforms: np.ndarray) -> pd.DataFrame:
    expected = len(person_values) * len(RATER_TRUTH) * len(TASK_TRUTH) * len(CRITERION_TRUTH)
    if len(row_uniforms) != expected:
        raise ValueError(f"Expected {expected} row uniforms, got {len(row_uniforms)}")
    rows: list[dict[str, Any]] = []
    uniform_index = 0
    for person_index, theta in enumerate(person_values, start=1):
        for rater, rater_value in RATER_TRUTH.items():
            for task, task_value in TASK_TRUTH.items():
                for criterion, criterion_value in CRITERION_TRUTH.items():
                    rows.append({
                        "Person": f"P{person_index:03d}",
                        "Rater": rater,
                        "Task": task,
                        "Criterion": criterion,
                        "Theta": float(theta),
                        "Eta": float(theta - rater_value - task_value - criterion_value),
                        "Uniform": float(row_uniforms[uniform_index]),
                    })
                    uniform_index += 1
    return pd.DataFrame(rows)


def generate_study_bundle(
    *,
    replicates: int,
    base_seed: int = BASE_SEED,
    replicate_start: int = 1,
    plan_path: Path = PLAN_PATH,
    execution_amendment_path: Path = AMENDMENT_PATH,
) -> dict[str, pd.DataFrame]:
    if int(replicates) < 1:
        raise ValueError("replicates must be >= 1")
    if int(replicate_start) < 1:
        raise ValueError("replicate_start must be >= 1")
    plan_path = plan_path.resolve()
    execution_amendment_path = execution_amendment_path.resolve()
    plan_sha256 = sha256_file(plan_path)
    amendment_sha256 = sha256_file(execution_amendment_path)
    script_sha256 = sha256_file(Path(__file__).resolve())
    manifests: list[dict[str, Any]] = []
    ratings_parts: list[pd.DataFrame] = []
    truth_parts: list[pd.DataFrame] = []
    threshold_parts: list[dict[str, Any]] = []
    replicate_stop = int(replicate_start) + int(replicates)
    for replicate in range(int(replicate_start), replicate_stop):
        seed = deterministic_replicate_seed(int(base_seed), "person-shape-study", replicate)
        shapes = generate_person_shapes(seed)
        uniform_seed = deterministic_replicate_seed(int(base_seed), "row-uniforms", replicate)
        uniform_rng = np.random.default_rng(uniform_seed)
        row_count = N_PERSONS * len(RATER_TRUTH) * len(TASK_TRUTH) * len(CRITERION_TRUTH)
        row_uniforms = uniform_rng.random(row_count)
        uniform_hash = hashlib.sha256(row_uniforms.tobytes()).hexdigest()
        for distribution in DISTRIBUTIONS:
            persons = shapes[distribution]
            latent = _latent_rows(persons, row_uniforms)
            complete = apply_threshold_condition(
                latent, THRESHOLD_CONDITIONS[THRESHOLD_CONDITION]
            )
            complete_hash = _frame_sha256(
                complete, ["Person", "Rater", "Task", "Criterion", "Score"]
            )
            person_truth = {
                f"P{index:03d}": float(value)
                for index, value in enumerate(persons, start=1)
            }
            for design in DESIGNS:
                run_id = f"{design}__{distribution}::rep-{replicate:05d}"
                run_ratings = apply_observation_design(complete, design)
                design_audit = constrained_adjacent_design_audit(run_ratings, FIT_MODEL)
                if design_audit["Nullity"] != 0:
                    raise RuntimeError(f"Prepared rank-deficient registered run: {run_id}")
                manifests.append({
                    "RunId": run_id,
                    "ConditionId": f"{design}__{distribution}",
                    "Design": design,
                    "PersonDistribution": distribution,
                    "TruthBias": 0.0,
                    "Replicate": replicate,
                    "Seed": int(seed),
                    "UniformSeed": int(uniform_seed),
                    "Categories": CATEGORIES,
                    "ThresholdCondition": THRESHOLD_CONDITION,
                    "PersonMeanRealized": float(np.mean(persons)),
                    "PersonSDRealized": float(np.std(persons, ddof=0)),
                    "PersonSkewRealized": float(skew(persons, bias=False)),
                    "PersonExcessKurtosisRealized": float(kurtosis(persons, fisher=True, bias=False)),
                    "ExpectedRows": int(len(run_ratings)),
                    "ExpectedStructuralNullity": int(design_audit["Nullity"]),
                    "PersonRaterComponents": int(design_audit["PersonRaterComponents"]),
                    "SharedUniformSHA256": uniform_hash,
                    "CompleteResponseSHA256": complete_hash,
                })
                run_ratings.insert(0, "RunId", run_id)
                ratings_parts.append(run_ratings)
                facet_rows = [
                    {"RunId": run_id, "Facet": "Person", "Level": level, "Truth": value}
                    for level, value in person_truth.items()
                ]
                for facet, values in (
                    ("Rater", RATER_TRUTH),
                    ("Task", TASK_TRUTH),
                    ("Criterion", CRITERION_TRUTH),
                ):
                    facet_rows.extend(
                        {"RunId": run_id, "Facet": facet, "Level": level, "Truth": value}
                        for level, value in values.items()
                    )
                truth_parts.append(pd.DataFrame(facet_rows))
                for criterion, vector in THRESHOLD_CONDITIONS[THRESHOLD_CONDITION].items():
                    for category, value in enumerate(vector, start=1):
                        threshold_parts.append({
                            "RunId": run_id,
                            "ConditionId": f"{design}__{distribution}",
                            "Replicate": replicate,
                            "Seed": int(seed),
                            "StepFacetLevel": criterion,
                            "Category": category,
                            "ThresholdTruth": float(value),
                        })
    manifest = pd.DataFrame(manifests)
    ratings = pd.concat(ratings_parts, ignore_index=True)
    truth = pd.concat(truth_parts, ignore_index=True)
    thresholds = pd.DataFrame(threshold_parts)
    attempts: list[dict[str, Any]] = []
    ordinal = 0
    for manifest_row in manifest.itertuples(index=False):
        run_ratings = ratings[ratings["RunId"].astype(str).eq(str(manifest_row.RunId))]
        run_truth = truth[truth["RunId"].astype(str).eq(str(manifest_row.RunId))]
        run_thresholds = thresholds[
            thresholds["RunId"].astype(str).eq(str(manifest_row.RunId))
        ]
        input_hash = hashlib.sha256(
            (
                _frame_digest(run_ratings)
                + _frame_digest(run_truth)
                + _frame_digest(run_thresholds)
            ).encode("ascii")
        ).hexdigest()
        for attempt_type in ATTEMPT_TYPES:
            attempt_id = f"{manifest_row.RunId}::{attempt_type}"
            fingerprint = hashlib.sha256(
                (
                    f"{attempt_id}|{input_hash}|{plan_sha256}|"
                    f"{amendment_sha256}|{script_sha256}"
                ).encode("utf-8")
            ).hexdigest()
            attempts.append({
                "AttemptOrdinal": ordinal,
                "AttemptId": attempt_id,
                "RunId": manifest_row.RunId,
                "AttemptType": attempt_type,
                "Replicate": manifest_row.Replicate,
                "Design": manifest_row.Design,
                "PersonDistribution": manifest_row.PersonDistribution,
                "RunInputSHA256": input_hash,
                "AttemptFingerprint": fingerprint,
            })
            ordinal += 1
    return {
        "manifest.csv": manifest,
        "generated_ratings.csv": ratings,
        "generated_facet_truth.csv": truth,
        "generated_anchors.csv": pd.DataFrame(columns=["RunId", "Facet", "Level", "Anchor"]),
        "generated_pcm_threshold_truth.csv": thresholds,
        "attempt_manifest.csv": pd.DataFrame(attempts),
    }


def prepare_study(
    study_dir: Path,
    *,
    replicates: int,
    base_seed: int,
    replicate_start: int = 1,
    phase: str | None = None,
    plan_path: Path = PLAN_PATH,
    execution_amendment_path: Path = AMENDMENT_PATH,
) -> None:
    study_dir = study_dir.resolve()
    if study_dir.exists():
        raise FileExistsError(f"Study directory already exists: {study_dir}")
    if phase is None:
        phase = "preflight" if int(replicates) == 1 else "screening"
    if phase not in {"preflight", "screening", "confirmatory"}:
        raise ValueError("phase must be preflight, screening, or confirmatory")
    plan_path = plan_path.resolve()
    execution_amendment_path = execution_amendment_path.resolve()
    if phase == "confirmatory":
        _validate_confirmatory_prepare_contract(
            plan_path,
            replicates=replicates,
            replicate_start=replicate_start,
        )
    study_dir.mkdir(parents=True)
    input_dir = study_dir / "retained_input"
    input_dir.mkdir()
    bundle = generate_study_bundle(
        replicates=replicates,
        base_seed=base_seed,
        replicate_start=replicate_start,
        plan_path=plan_path,
        execution_amendment_path=execution_amendment_path,
    )
    for filename, frame in bundle.items():
        frame.to_csv(input_dir / filename, index=False, lineterminator="\n")
    hashes = {filename: sha256_file(input_dir / filename) for filename in INPUT_FILES}
    _json_dump(input_dir / "bundle_hashes.json", hashes)
    manifest = bundle["manifest.csv"]
    attempts = bundle["attempt_manifest.csv"]
    identity = {
        "schema_version": SCHEMA_VERSION,
        "phase": phase,
        "replicates": int(replicates),
        "replicate_start": int(replicate_start),
        "replicate_end": int(replicate_start) + int(replicates) - 1,
        "base_seed": int(base_seed),
        "datasets": int(len(manifest)),
        "attempts": int(len(attempts)),
        "plan_file": _identity_file_name(plan_path),
        "plan_sha256": sha256_file(plan_path),
        "infrastructure_amendment_file": _identity_file_name(execution_amendment_path),
        "infrastructure_amendment_sha256": sha256_file(execution_amendment_path),
        "script_sha256_at_prepare": sha256_file(Path(__file__).resolve()),
        "retained_input_sha256": hashes,
        "python_version": sys.version,
        "platform": platform.platform(),
        "claim_limit": (
            "Operational preflight only; no estimator ranking or bias claim."
            if phase == "preflight" else
            "Screening only; below confirmatory study depth and no estimator ranking."
            if phase == "screening" else
            "Confirmatory only for endpoints registered in the study plan; no estimator "
            "ranking or cross-likelihood comparison."
        ),
    }
    _json_dump(study_dir / "study_identity.json", identity)
    print(f"Prepared {len(manifest)} datasets and {len(attempts)} attempts in {study_dir}")


def validate_study_identity(
    study_dir: Path,
    *,
    require_current_execution_script: bool = True,
) -> dict[str, Any]:
    identity_path = study_dir / "study_identity.json"
    if not identity_path.is_file():
        raise FileNotFoundError(f"Study identity missing: {identity_path}")
    identity = json.loads(identity_path.read_text(encoding="utf-8"))
    plan_path = _identity_file_path(identity.get("plan_file"), fallback=PLAN_PATH)
    amendment_path = _identity_file_path(
        identity.get("infrastructure_amendment_file"), fallback=AMENDMENT_PATH
    )
    if identity.get("plan_sha256") != sha256_file(plan_path):
        raise ValueError("Registered plan hash changed")
    if identity.get("infrastructure_amendment_sha256") != sha256_file(amendment_path):
        raise ValueError("Registered infrastructure amendment hash changed")
    if (
        require_current_execution_script
        and identity.get("script_sha256_at_prepare") != sha256_file(Path(__file__).resolve())
    ):
        raise ValueError("Execution script changed after study preparation")
    input_dir = study_dir / "retained_input"
    for filename, expected in identity.get("retained_input_sha256", {}).items():
        actual = sha256_file(input_dir / filename)
        if actual.lower() != str(expected).lower():
            raise ValueError(f"Retained input hash mismatch for {filename}")
    return identity


def _mml_constraint(facets: pd.DataFrame, thresholds: pd.DataFrame) -> dict[str, Any]:
    facet_sums = facets.groupby("Facet")["Estimate"].sum().to_dict()
    constrained = {key: float(value) for key, value in facet_sums.items() if key in {"Rater", "Task"}}
    step_sums = thresholds.groupby("StepFacetLevel")["Estimate"].sum().to_dict()
    residuals = [abs(value) for value in constrained.values()]
    residuals.extend(abs(float(value)) for value in step_sums.values())
    maximum = max(residuals) if residuals else np.nan
    return {
        "ConstrainedFacetSums": constrained,
        "StepSums": {key: float(value) for key, value in step_sums.items()},
        "MaxAbsConstraintResidual": float(maximum),
        "ConstraintPass": bool(np.isfinite(maximum) and maximum <= CONSTRAINT_TOLERANCE),
    }


def _fit_mml_shape(
    app: Any,
    ratings: pd.DataFrame,
    truth: pd.DataFrame,
    threshold_truth: pd.DataFrame,
    manifest_row: pd.Series,
    *,
    fixed_sd: bool,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    mode = MML_FIXED_MODE if fixed_sd else MML_FREE_MODE
    registered_mode = "PYTHON_MML_FIXED_SD1_Q31" if fixed_sd else MML_FREE_MODE
    kwargs = build_mml_kwargs(registered_mode, int(manifest_row["Categories"]))
    kwargs.update({
        "model": FIT_MODEL,
        "step_facet": "Criterion",
        "estimate_population_sd": not fixed_sd,
        "population_prior_sd": PERSON_SD,
    })
    started = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = app.mfrm_estimate(
            data=ratings[["Person", "Rater", "Task", "Criterion", "Score"]].copy(),
            **kwargs,
        )
    summary = result["summary"].iloc[0]
    facets = result["facets"]["others"].copy()
    thresholds = _normalize_steps(
        result,
        fit_model=FIT_MODEL,
        estimator_mode=registered_mode,
        run_id=str(manifest_row["RunId"]),
        threshold_condition=THRESHOLD_CONDITION,
        threshold_truth=threshold_truth,
    )
    thresholds["EstimatorMode"] = mode
    constraint = _mml_constraint(facets, thresholds)
    finite_main = pd.to_numeric(facets["Estimate"], errors="coerce").notna().all()
    finite_steps = pd.to_numeric(thresholds["Estimate"], errors="coerce").notna().all()
    included = bool(
        summary.get("Converged", False)
        and summary.get("InferenceReady", False)
        and finite_main
        and finite_steps
        and constraint["ConstraintPass"]
    )
    recovery = _normalize_recovery(
        facets,
        truth,
        manifest_row,
        estimator_family="MML",
        estimator_mode=registered_mode,
        fit_model=FIT_MODEL,
        included=included,
    )
    recovery["EstimatorMode"] = mode
    recovery["IncludedInStudy"] = recovery.pop("IncludedInBridge")
    thresholds["IncludedInStudy"] = bool(included) & thresholds["TruthError"].notna()
    run = pd.DataFrame([{
        "RunId": manifest_row["RunId"],
        "ConditionId": manifest_row["ConditionId"],
        "Design": manifest_row["Design"],
        "PersonDistribution": manifest_row["PersonDistribution"],
        "Replicate": int(manifest_row["Replicate"]),
        "EstimatorFamily": "MML",
        "EstimatorMode": mode,
        "EstimandClass": "normal_person_population_marginal_likelihood",
        "LikelihoodBasis": "marginal",
        "Rows": int(len(ratings)),
        "FitReturned": True,
        "Converged": bool(summary.get("Converged", False)),
        "InferenceReady": bool(summary.get("InferenceReady", False)),
        "IncludedInStudy": included,
        "ConstraintPass": constraint["ConstraintPass"],
        "MaxAbsConstraintResidual": constraint["MaxAbsConstraintResidual"],
        "PopulationSDMode": "fixed" if fixed_sd else "estimated",
        "PopulationSDInput": PERSON_SD,
        "EstimatedPopulationSD": summary.get("EstimatedPopulationSD"),
        "PersonSDRealized": manifest_row["PersonSDRealized"],
        "LogLik": summary.get("LogLik"),
        "AIC": summary.get("AIC"),
        "BIC": summary.get("BIC"),
        "Warnings": " | ".join(str(item.message) for item in caught),
        "ElapsedSeconds": time.perf_counter() - started,
        "FailureReason": "" if included else "MML readiness or constraint gate failed",
    }])
    return run, recovery, thresholds, constraint


def _jmle_outputs(
    manifest_row: pd.Series,
    metrics: dict[str, Any],
    main: pd.DataFrame,
    thresholds: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    base = {
        "RunId": manifest_row["RunId"],
        "ConditionId": manifest_row["ConditionId"],
        "Design": manifest_row["Design"],
        "PersonDistribution": manifest_row["PersonDistribution"],
        "Replicate": int(manifest_row["Replicate"]),
        "EstimatorFamily": "JMLE",
        "EstimandClass": "fixed_person_joint_likelihood",
        "LikelihoodBasis": "joint",
        "Rows": int(metrics["Rows"]),
        "ComparisonEligible": bool(metrics["ComparisonEligible"]),
        "DirectAgreementPass": bool(metrics["DirectAgreementPass"]),
        "MainWeightedMAE": metrics["MainWeightedMAE"],
        "MainMaxAbsDifference": metrics["MainMaxAbsDifference"],
        "ThresholdWeightedMAE": metrics["ThresholdWeightedMAE"],
        "ThresholdMaxAbsDifference": metrics["ThresholdMaxAbsDifference"],
        "MinimumWithinFacetSpearman": metrics["MinimumWithinFacetSpearman"],
        "ElapsedSeconds": metrics["ElapsedSeconds"],
    }
    runs = []
    for mode, engine, converged, ready in (
        (FACETS_MODE, "FACETS", metrics["FACETSConverged"], metrics["ComparisonEligible"]),
        (PYTHON_JMLE_MODE, "PythonApp", metrics["PythonConverged"], metrics["PythonInferenceReady"]),
    ):
        runs.append({
            **base,
            "EstimatorMode": mode,
            "Engine": engine,
            "FitReturned": True,
            "Converged": bool(converged),
            "InferenceReady": bool(ready),
            "IncludedInStudy": bool(metrics["ComparisonEligible"]),
            "FailureReason": "" if metrics["ComparisonEligible"] else "JMLE parity eligibility failed",
        })
    recovery_parts = []
    threshold_parts = []
    for mode, engine, estimate_col, error_col in (
        (FACETS_MODE, "FACETS", "FACETSEstimateAligned", "FACETSTruthErrorAligned"),
        (PYTHON_JMLE_MODE, "PythonApp", "PythonEstimateAligned", "PythonTruthErrorAligned"),
    ):
        recovery_parts.append(pd.DataFrame({
            "RunId": main["RunId"],
            "ConditionId": main["ConditionId"],
            "Design": manifest_row["Design"],
            "PersonDistribution": manifest_row["PersonDistribution"],
            "Replicate": manifest_row["Replicate"],
            "EstimatorFamily": "JMLE",
            "EstimatorMode": mode,
            "EstimandClass": "fixed_person_joint_likelihood",
            "Engine": engine,
            "FitModel": FIT_MODEL,
            "Facet": main["Facet"],
            "Level": main["Level"],
            "Truth": main["Truth"],
            "EstimateAligned": main[estimate_col],
            "TruthAligned": main["TruthAligned"],
            "ErrorAligned": main[error_col],
            "SE": main["SE"] if mode == FACETS_MODE else np.nan,
            "IncludedInStudy": main["ComparisonEligible"],
        }))
        threshold_parts.append(pd.DataFrame({
            "RunId": manifest_row["RunId"],
            "Design": manifest_row["Design"],
            "PersonDistribution": manifest_row["PersonDistribution"],
            "Replicate": manifest_row["Replicate"],
            "EstimatorFamily": "JMLE",
            "EstimatorMode": mode,
            "EstimandClass": "fixed_person_joint_likelihood",
            "Engine": engine,
            "FitModel": FIT_MODEL,
            "StepFacetLevel": thresholds["StepFacetLevel"],
            "Category": thresholds["Category"],
            "Estimate": (
                thresholds["ThresholdMeasureDisplayed"]
                if mode == FACETS_MODE else thresholds["PythonThreshold"]
            ),
            "ThresholdTruth": thresholds["ThresholdTruth"],
            "TruthError": (
                thresholds["FACETSTruthError"]
                if mode == FACETS_MODE else thresholds["PythonTruthError"]
            ),
            "TruthTargetStatus": thresholds["TruthTargetStatus"],
            "IncludedInStudy": thresholds["ComparisonEligible"],
            "DisplayPrecisionContract": (
                "FACETS displayed measure; Umean=6 primary, no invented raw precision"
                if mode == FACETS_MODE else "native numeric estimate"
            ),
        }))
    return pd.DataFrame(runs), pd.concat(recovery_parts), pd.concat(threshold_parts)


def _run_one_attempt(
    attempt: pd.Series,
    manifest_row: pd.Series,
    ratings: pd.DataFrame,
    truth: pd.DataFrame,
    threshold_truth: pd.DataFrame,
    *,
    facets_exe: Path,
    stage_dir: Path,
    timeout_seconds: float,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list[dict[str, Any]]]:
    attempt_type = str(attempt["AttemptType"])
    if attempt_type == "FACETS_PYTHON_JMLE_PCM":
        work_root = stage_dir / "facets_runs"
        work_root.mkdir()
        metrics, main, thresholds, _, _ = _fit_facets_and_python(
            manifest_row,
            ratings,
            truth,
            threshold_truth,
            fit_model=FIT_MODEL,
            facets_exe=facets_exe,
            work_root=work_root,
            timeout_seconds=timeout_seconds,
        )
        runs, recovery, threshold_rows = _jmle_outputs(manifest_row, metrics, main, thresholds)
        return runs, recovery, threshold_rows, []
    import streamlit_app as app  # pylint: disable=import-outside-toplevel
    if attempt_type in {"PYTHON_MML_FIXED_SD08_Q31_PCM", "PYTHON_MML_FREE_SD_Q31_PCM"}:
        run, recovery, thresholds, constraint = _fit_mml_shape(
            app,
            ratings,
            truth,
            threshold_truth,
            manifest_row,
            fixed_sd=attempt_type.startswith("PYTHON_MML_FIXED"),
        )
        constraint_row = {
            "RunId": manifest_row["RunId"],
            "EstimatorMode": run.iloc[0]["EstimatorMode"],
            **constraint,
        }
        return run, recovery, thresholds, [constraint_row]
    if attempt_type == "PYTHON_EXACT_CMLE_PCM":
        run_dict, recovery, thresholds, constraint = fit_exact_cmle(
            ratings,
            truth,
            threshold_truth,
            manifest_row,
            fit_model=FIT_MODEL,
        )
        run_dict["ConditionId"] = manifest_row["ConditionId"]
        run_dict["PersonDistribution"] = manifest_row["PersonDistribution"]
        run_dict["EstimandClass"] = "person_total_conditional_likelihood"
        run_dict["LikelihoodBasis"] = "conditional"
        run_dict["IncludedInStudy"] = run_dict.pop("IncludedInBridge")
        run = pd.DataFrame([run_dict])
        recovery["PersonDistribution"] = manifest_row["PersonDistribution"]
        recovery["IncludedInStudy"] = recovery.pop("IncludedInBridge")
        thresholds["PersonDistribution"] = manifest_row["PersonDistribution"]
        thresholds["IncludedInStudy"] = run_dict["IncludedInStudy"] & thresholds["TruthError"].notna()
        constraint_row = {
            "RunId": manifest_row["RunId"],
            "EstimatorMode": CMLE_MODE,
            **constraint,
        }
        return run, recovery, thresholds, [constraint_row]
    raise ValueError(f"Unknown attempt type: {attempt_type}")


def _attempt_dir(study_dir: Path, attempt: pd.Series) -> Path:
    ordinal = int(attempt["AttemptOrdinal"])
    return study_dir / "attempts" / f"{ordinal:05d}"


def _validate_completion(
    completion: dict[str, Any],
    attempt: pd.Series,
    *,
    expected_plan_sha256: str | None = None,
    expected_script_sha256: str | None = None,
    expected_amendment_sha256: str | None = None,
) -> None:
    expected = str(attempt["AttemptFingerprint"])
    if completion.get("attempt_fingerprint") != expected:
        raise ValueError(f"Completed attempt fingerprint mismatch: {attempt['AttemptId']}")
    expected_plan = expected_plan_sha256 or sha256_file(PLAN_PATH)
    if completion.get("plan_sha256") != expected_plan:
        raise ValueError(f"Completed attempt plan mismatch: {attempt['AttemptId']}")
    expected_amendment = expected_amendment_sha256 or sha256_file(AMENDMENT_PATH)
    expected_script = expected_script_sha256 or sha256_file(Path(__file__).resolve())
    if completion.get("infrastructure_amendment_sha256") != expected_amendment:
        raise ValueError(f"Completed attempt amendment mismatch: {attempt['AttemptId']}")
    if completion.get("script_sha256") != expected_script:
        raise ValueError(f"Completed attempt script mismatch: {attempt['AttemptId']}")


def select_shard_attempts(
    attempts: pd.DataFrame,
    *,
    shard_index: int,
    shard_count: int,
) -> pd.DataFrame:
    """Return one deterministic, disjoint modulo shard."""

    if shard_count < 1 or shard_index < 0 or shard_index >= shard_count:
        raise ValueError("Require 0 <= shard-index < shard-count")
    if "AttemptOrdinal" not in attempts:
        raise ValueError("attempt manifest is missing AttemptOrdinal")
    ordinals = pd.to_numeric(attempts["AttemptOrdinal"], errors="raise").astype(int)
    if ordinals.duplicated().any():
        raise ValueError("attempt manifest contains duplicate ordinals")
    return attempts[ordinals.mod(shard_count).eq(shard_index)].copy()


def run_shard(
    study_dir: Path,
    *,
    facets_exe: Path,
    shard_index: int,
    shard_count: int,
    resume: bool,
    timeout_seconds: float,
) -> dict[str, int]:
    study_dir = study_dir.resolve()
    identity = validate_study_identity(study_dir)
    if not facets_exe.is_file():
        raise FileNotFoundError(f"FACETS executable missing: {facets_exe}")
    input_dir = study_dir / "retained_input"
    tables = validate_bundle(input_dir)
    attempts = pd.read_csv(input_dir / "attempt_manifest.csv")
    selected = select_shard_attempts(
        attempts, shard_index=shard_index, shard_count=shard_count
    )
    manifest = tables["manifest.csv"].set_index("RunId", drop=False)
    ratings_all = tables["generated_ratings.csv"]
    truth_all = tables["generated_facet_truth.csv"]
    threshold_all = pd.read_csv(input_dir / "generated_pcm_threshold_truth.csv")
    (study_dir / "attempts").mkdir(exist_ok=True)
    completed_count = skipped_count = failed_count = 0
    for sequence, (_, attempt) in enumerate(selected.iterrows(), start=1):
        final_dir = _attempt_dir(study_dir, attempt)
        completion_path = final_dir / "completion.json"
        if completion_path.is_file():
            completion = json.loads(completion_path.read_text(encoding="utf-8"))
            _validate_completion(
                completion,
                attempt,
                expected_plan_sha256=str(identity["plan_sha256"]),
                expected_script_sha256=str(identity["script_sha256_at_prepare"]),
                expected_amendment_sha256=str(identity["infrastructure_amendment_sha256"]),
            )
            if not resume:
                raise FileExistsError(f"Attempt already completed: {attempt['AttemptId']}")
            skipped_count += 1
            print(f"[{sequence}/{len(selected)}] skip {attempt['AttemptId']}", flush=True)
            continue
        final_dir.mkdir(parents=True, exist_ok=True)
        stage_dir = final_dir / f"_work__{os.getpid()}__{time.time_ns()}"
        stage_dir.mkdir()
        run_id = str(attempt["RunId"])
        manifest_row = manifest.loc[run_id]
        ratings = ratings_all[ratings_all["RunId"].astype(str).eq(run_id)].drop(columns="RunId")
        truth = truth_all[truth_all["RunId"].astype(str).eq(run_id)]
        threshold_truth = threshold_all[threshold_all["RunId"].astype(str).eq(run_id)]
        succeeded = True
        error = ""
        try:
            runs, recovery, thresholds, constraints = _run_one_attempt(
                attempt,
                manifest_row,
                ratings,
                truth,
                threshold_truth,
                facets_exe=facets_exe,
                stage_dir=stage_dir,
                timeout_seconds=timeout_seconds,
            )
        except Exception as exc:  # retain the planned denominator and traceback class
            succeeded = False
            failed_count += 1
            error = f"{type(exc).__name__}: {exc}"
            runs = pd.DataFrame([{
                "RunId": run_id,
                "ConditionId": manifest_row["ConditionId"],
                "Design": manifest_row["Design"],
                "PersonDistribution": manifest_row["PersonDistribution"],
                "Replicate": int(manifest_row["Replicate"]),
                "EstimatorMode": attempt["AttemptType"],
                "FitReturned": False,
                "Converged": False,
                "InferenceReady": False,
                "IncludedInStudy": False,
                "FailureReason": error,
            }])
            recovery = pd.DataFrame()
            thresholds = pd.DataFrame()
            constraints = []
        runs.insert(0, "AttemptId", attempt["AttemptId"])
        runs.to_csv(stage_dir / "run_ledger.csv", index=False, lineterminator="\n")
        recovery.to_csv(stage_dir / "recovery.csv", index=False, lineterminator="\n")
        thresholds.to_csv(stage_dir / "thresholds.csv", index=False, lineterminator="\n")
        pd.DataFrame(constraints).to_csv(
            stage_dir / "constraints.csv", index=False, lineterminator="\n"
        )
        artifact_files = ("run_ledger.csv", "recovery.csv", "thresholds.csv", "constraints.csv")
        for filename in artifact_files:
            shutil.copy2(stage_dir / filename, final_dir / filename)
        completion = {
            "schema_version": SCHEMA_VERSION,
            "attempt_id": attempt["AttemptId"],
            "attempt_fingerprint": attempt["AttemptFingerprint"],
            "run_input_sha256": attempt["RunInputSHA256"],
            "plan_sha256": identity["plan_sha256"],
            "infrastructure_amendment_sha256": identity["infrastructure_amendment_sha256"],
            "script_sha256": identity["script_sha256_at_prepare"],
            "attempt_succeeded": succeeded,
            "failure_reason": error,
            "artifact_sha256": {
                filename: sha256_file(final_dir / filename)
                for filename in artifact_files
            },
            "retained_work_generation": stage_dir.name,
        }
        marker_temp = final_dir / f"completion.json.tmp.{os.getpid()}.{time.time_ns()}"
        _json_dump(marker_temp, completion)
        os.replace(marker_temp, completion_path)
        completed_count += 1
        print(
            f"[{sequence}/{len(selected)}] complete {attempt['AttemptId']} succeeded={succeeded}",
            flush=True,
        )
    summary = {
        "assigned": int(len(selected)),
        "completed_now": completed_count,
        "skipped_valid_completed": skipped_count,
        "failed_now": failed_count,
    }
    print(json.dumps(summary, sort_keys=True), flush=True)
    return summary


def _rmse(values: pd.Series) -> float:
    numeric = pd.to_numeric(values, errors="coerce").dropna().to_numpy(dtype=float)
    return float(np.sqrt(np.mean(np.square(numeric)))) if len(numeric) else np.nan


def enrich_registered_metadata(
    frame: pd.DataFrame,
    registered_manifest: pd.DataFrame,
    *,
    label: str,
) -> pd.DataFrame:
    """Attach immutable condition metadata by RunId and reject conflicts."""

    if frame.empty:
        return frame.copy()
    metadata_columns = ["RunId", "ConditionId", "Design", "PersonDistribution", "Replicate"]
    missing = set(metadata_columns).difference(registered_manifest.columns)
    if missing:
        raise ValueError(f"Registered manifest is missing metadata: {sorted(missing)}")
    registered_metadata = registered_manifest[metadata_columns].copy()
    if registered_metadata["RunId"].duplicated().any():
        raise ValueError("Registered manifest contains duplicate RunId values")
    if "RunId" not in frame:
        raise ValueError(f"{label} is missing RunId")
    output = frame.copy()
    for column in metadata_columns[1:]:
        if column not in output:
            continue
        check = output[["RunId", column]].merge(
            registered_metadata[["RunId", column]],
            on="RunId",
            how="left",
            suffixes=("Artifact", "Registered"),
            validate="many_to_one",
        )
        artifact = check[f"{column}Artifact"]
        registered = check[f"{column}Registered"]
        if column == "Replicate":
            mismatch = (
                artifact.notna()
                & pd.to_numeric(artifact, errors="coerce").ne(
                    pd.to_numeric(registered, errors="coerce")
                )
            )
        else:
            mismatch = artifact.notna() & artifact.astype(str).ne(registered.astype(str))
        if mismatch.any():
            raise ValueError(f"{label} conflicts with registered {column}")
        output = output.drop(columns=column)
    return output.merge(
        registered_metadata,
        on="RunId",
        how="left",
        validate="many_to_one",
    )


def aggregate_study(
    study_dir: Path,
    *,
    aggregate_name: str = "aggregate",
) -> dict[str, Any]:
    study_dir = study_dir.resolve()
    if not aggregate_name or any(character not in "abcdefghijklmnopqrstuvwxyz0123456789_" for character in aggregate_name):
        raise ValueError("aggregate-name must contain only lowercase letters, digits, and underscore")
    identity = validate_study_identity(
        study_dir, require_current_execution_script=False
    )
    input_dir = study_dir / "retained_input"
    attempts = pd.read_csv(input_dir / "attempt_manifest.csv")
    registered_manifest = pd.read_csv(input_dir / "manifest.csv")
    completion_rows = []
    run_parts = []
    recovery_parts = []
    threshold_parts = []
    constraint_parts = []
    marker_hashes: dict[str, str] = {}
    for _, attempt in attempts.iterrows():
        attempt_dir = _attempt_dir(study_dir, attempt)
        marker = attempt_dir / "completion.json"
        if not marker.is_file():
            completion_rows.append({
                **attempt.to_dict(), "Completed": False, "AttemptSucceeded": False,
                "FailureReason": "missing completion marker",
            })
            continue
        completion = json.loads(marker.read_text(encoding="utf-8"))
        _validate_completion(
            completion,
            attempt,
            expected_plan_sha256=str(identity["plan_sha256"]),
            expected_script_sha256=str(identity["script_sha256_at_prepare"]),
            expected_amendment_sha256=str(identity["infrastructure_amendment_sha256"]),
        )
        for filename, expected in completion["artifact_sha256"].items():
            actual = sha256_file(attempt_dir / filename)
            if actual.lower() != str(expected).lower():
                raise ValueError(f"Attempt artifact hash mismatch: {attempt['AttemptId']} {filename}")
        marker_hashes[str(attempt["AttemptId"])] = sha256_file(marker)
        completion_rows.append({
            **attempt.to_dict(),
            "Completed": True,
            "AttemptSucceeded": bool(completion["attempt_succeeded"]),
            "FailureReason": completion.get("failure_reason", ""),
        })
        for filename, target in (
            ("run_ledger.csv", run_parts),
            ("recovery.csv", recovery_parts),
            ("thresholds.csv", threshold_parts),
            ("constraints.csv", constraint_parts),
        ):
            path = attempt_dir / filename
            if path.stat().st_size > 1:
                try:
                    frame = pd.read_csv(path)
                except pd.errors.EmptyDataError:
                    continue
                if len(frame):
                    target.append(frame)
    completion_ledger = pd.DataFrame(completion_rows)
    runs = pd.concat(run_parts, ignore_index=True, sort=False) if run_parts else pd.DataFrame()
    recovery = pd.concat(recovery_parts, ignore_index=True, sort=False) if recovery_parts else pd.DataFrame()
    thresholds = pd.concat(threshold_parts, ignore_index=True, sort=False) if threshold_parts else pd.DataFrame()
    constraints = pd.concat(constraint_parts, ignore_index=True, sort=False) if constraint_parts else pd.DataFrame()

    runs = enrich_registered_metadata(runs, registered_manifest, label="run ledger")
    recovery = enrich_registered_metadata(recovery, registered_manifest, label="recovery")
    thresholds = enrich_registered_metadata(thresholds, registered_manifest, label="thresholds")
    constraints = enrich_registered_metadata(constraints, registered_manifest, label="constraints")
    aggregate_dir = study_dir / aggregate_name
    if aggregate_dir.exists():
        raise FileExistsError(f"Aggregate directory already exists: {aggregate_dir}")
    aggregate_dir.mkdir()
    completion_ledger.to_csv(aggregate_dir / "attempt_ledger.csv", index=False)
    runs.to_csv(aggregate_dir / "study_run_ledger.csv", index=False)
    recovery.to_csv(aggregate_dir / "study_recovery.csv", index=False)
    thresholds.to_csv(aggregate_dir / "study_thresholds.csv", index=False)
    constraints.to_csv(aggregate_dir / "study_constraints.csv", index=False)

    included_recovery = recovery[
        recovery.get("IncludedInStudy", pd.Series(False, index=recovery.index)).fillna(False).astype(bool)
    ].copy()
    recovery_summary = pd.DataFrame()
    run_facet = pd.DataFrame()
    contrasts = pd.DataFrame()
    if len(included_recovery):
        recovery_summary = (
            included_recovery.groupby(
                ["EstimatorMode", "PersonDistribution", "Design", "Facet"], dropna=False
            )["ErrorAligned"]
            .agg(Parameters="size", MeanError="mean", MAE=lambda x: float(x.abs().mean()), RMSE=_rmse)
            .reset_index()
        )
        run_facet = (
            included_recovery.groupby(
                ["RunId", "Replicate", "EstimatorMode", "PersonDistribution", "Design", "Facet"],
                dropna=False,
            )["ErrorAligned"]
            .agg(MeanError="mean", RMSE=_rmse, MAE=lambda x: float(x.abs().mean()))
            .reset_index()
        )
        normal = run_facet[run_facet["PersonDistribution"].eq("normal")].drop(
            columns=["RunId", "PersonDistribution"]
        )
        nonnormal = run_facet[~run_facet["PersonDistribution"].eq("normal")]
        paired = nonnormal.merge(
            normal,
            on=["Replicate", "EstimatorMode", "Design", "Facet"],
            suffixes=("", "Normal"),
            validate="many_to_one",
        )
        if len(paired):
            for metric in ("MeanError", "RMSE", "MAE"):
                paired[f"{metric}ContrastVsNormal"] = paired[metric] - paired[f"{metric}Normal"]
            contrast_rows = []
            for keys, group in paired.groupby(
                ["EstimatorMode", "PersonDistribution", "Design", "Facet"], dropna=False
            ):
                for metric in ("MeanError", "RMSE", "MAE"):
                    values = group[f"{metric}ContrastVsNormal"].dropna()
                    n = len(values)
                    mean = float(values.mean()) if n else np.nan
                    sd = float(values.std(ddof=1)) if n > 1 else np.nan
                    se = sd / np.sqrt(n) if n > 1 else np.nan
                    contrast_rows.append({
                        "EstimatorMode": keys[0], "PersonDistribution": keys[1],
                        "Design": keys[2], "Facet": keys[3], "Metric": metric,
                        "PairedReplicates": n, "MeanContrastVsNormal": mean,
                        "MonteCarloSD": sd, "MonteCarloSE": se,
                        "ScreeningNormalApproxLower95": mean - 1.96 * se if np.isfinite(se) else np.nan,
                        "ScreeningNormalApproxUpper95": mean + 1.96 * se if np.isfinite(se) else np.nan,
                    })
            contrasts = pd.DataFrame(contrast_rows)
    recovery_summary.to_csv(aggregate_dir / "recovery_summary.csv", index=False)
    run_facet.to_csv(aggregate_dir / "run_facet_recovery.csv", index=False)
    contrasts.to_csv(aggregate_dir / "paired_shape_contrasts.csv", index=False)

    expected_attempts = len(attempts)
    completed = int(completion_ledger["Completed"].sum())
    succeeded = int(completion_ledger["AttemptSucceeded"].sum())
    jmle = runs[runs.get("EstimatorFamily", pd.Series(index=runs.index, dtype=object)).eq("JMLE")]
    mml = runs[runs.get("EstimatorFamily", pd.Series(index=runs.index, dtype=object)).eq("MML")]
    cmle = runs[runs.get("EstimatorFamily", pd.Series(index=runs.index, dtype=object)).eq("CMLE")]
    jmle_pairs = jmle[jmle.get("EstimatorMode", pd.Series(index=jmle.index, dtype=object)).eq(FACETS_MODE)]
    constraint_pass = bool(
        len(constraints) == int(identity["datasets"]) * 3
        and constraints["ConstraintPass"].fillna(False).astype(bool).all()
        and pd.to_numeric(constraints["MaxAbsConstraintResidual"], errors="coerce")
        .le(CONSTRAINT_TOLERANCE).all()
    )
    metadata_contract_pass = bool(
        all(
            frame.empty
            or (
                frame["PersonDistribution"].notna().all()
                and set(frame["PersonDistribution"].astype(str)).issubset(set(DISTRIBUTIONS))
            )
            for frame in (runs, recovery, thresholds)
        )
    )
    gates = {
        "attempts_complete": completed == expected_attempts,
        "attempts_succeeded": succeeded == expected_attempts,
        "facets_python_parity": bool(
            len(jmle_pairs) == int(identity["datasets"])
            and jmle_pairs["ComparisonEligible"].fillna(False).astype(bool).all()
            and jmle_pairs["DirectAgreementPass"].fillna(False).astype(bool).all()
        ),
        "mml_operational": bool(
            len(mml) == int(identity["datasets"]) * 2
            and mml["IncludedInStudy"].fillna(False).astype(bool).all()
        ),
        "cmle_operational": bool(
            len(cmle) == int(identity["datasets"])
            and cmle["IncludedInStudy"].fillna(False).astype(bool).all()
            and cmle["ConditionalNullity"].fillna(-1).astype(int).eq(0).all()
            and cmle["FiniteMLEExistenceQualified"].fillna(False).astype(bool).all()
        ),
        "constraints": constraint_pass,
        "metadata_complete": metadata_contract_pass,
    }
    metrics = {
        "schema_version": SCHEMA_VERSION,
        "phase": identity["phase"],
        "datasets": int(identity["datasets"]),
        "expected_attempts": int(expected_attempts),
        "completed_attempts": completed,
        "successful_attempts": succeeded,
        "run_ledger_rows": int(len(runs)),
        "recovery_rows": int(len(recovery)),
        "threshold_rows": int(len(thresholds)),
        "facets_python_pairs": int(len(jmle_pairs)),
        "mml_included": int(mml.get("IncludedInStudy", pd.Series(False, index=mml.index)).fillna(False).sum()),
        "cmle_included": int(cmle.get("IncludedInStudy", pd.Series(False, index=cmle.index)).fillna(False).sum()),
        "free_sd_range": [
            float(pd.to_numeric(mml.loc[mml["EstimatorMode"].eq(MML_FREE_MODE), "EstimatedPopulationSD"], errors="coerce").min()),
            float(pd.to_numeric(mml.loc[mml["EstimatorMode"].eq(MML_FREE_MODE), "EstimatedPopulationSD"], errors="coerce").max()),
        ] if len(mml) else [np.nan, np.nan],
        "cmle_extreme_person_range": [
            int(pd.to_numeric(cmle["PersonsExtreme"], errors="coerce").min()),
            int(pd.to_numeric(cmle["PersonsExtreme"], errors="coerce").max()),
        ] if len(cmle) else [0, 0],
        "gates": gates,
        "qualification_pass": bool(all(gates.values())),
        "claim_limit": identity["claim_limit"],
        "cross_basis_likelihood_comparison": "prohibited",
        "aggregate_name": aggregate_name,
    }
    _json_dump(aggregate_dir / "study_metrics.json", metrics)
    _json_dump(aggregate_dir / "completion_marker_hashes.json", marker_hashes)
    aggregate_identity = {
        "schema_version": SCHEMA_VERSION,
        "plan_sha256": identity["plan_sha256"],
        "frozen_execution_amendment_sha256": identity["infrastructure_amendment_sha256"],
        "frozen_execution_script_sha256": identity["script_sha256_at_prepare"],
        "aggregation_amendment_sha256": sha256_file(AGGREGATION_AMENDMENT_PATH),
        "aggregation_script_sha256": sha256_file(Path(__file__).resolve()),
        "study_identity_sha256": sha256_file(study_dir / "study_identity.json"),
        "completion_markers_sha256": hashlib.sha256(
            json.dumps(marker_hashes, sort_keys=True).encode("utf-8")
        ).hexdigest(),
        "metrics_sha256": sha256_file(aggregate_dir / "study_metrics.json"),
    }
    _json_dump(aggregate_dir / "aggregate_identity.json", aggregate_identity)
    print(json.dumps(metrics, indent=2, sort_keys=True))
    return metrics


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    prepare = subparsers.add_parser("prepare")
    prepare.add_argument("--study-dir", type=Path, required=True)
    prepare.add_argument("--replicates", type=int, required=True)
    prepare.add_argument("--base-seed", type=int, default=BASE_SEED)
    prepare.add_argument("--replicate-start", type=int, default=1)
    prepare.add_argument(
        "--phase", choices=("preflight", "screening", "confirmatory"), default=None
    )
    prepare.add_argument("--plan-path", type=Path, default=PLAN_PATH)
    prepare.add_argument(
        "--execution-amendment-path", type=Path, default=AMENDMENT_PATH
    )
    run = subparsers.add_parser("run")
    run.add_argument("--study-dir", type=Path, required=True)
    run.add_argument("--facets-exe", type=Path, required=True)
    run.add_argument("--shard-index", type=int, default=0)
    run.add_argument("--shard-count", type=int, default=1)
    run.add_argument("--resume", action="store_true")
    run.add_argument("--timeout-seconds", type=float, default=120.0)
    aggregate = subparsers.add_parser("aggregate")
    aggregate.add_argument("--study-dir", type=Path, required=True)
    aggregate.add_argument("--aggregate-name", default="aggregate")
    return parser.parse_args(argv)


def main() -> None:
    args = parse_args()
    if args.command == "prepare":
        prepare_study(
            args.study_dir,
            replicates=args.replicates,
            base_seed=args.base_seed,
            replicate_start=args.replicate_start,
            phase=args.phase,
            plan_path=args.plan_path,
            execution_amendment_path=args.execution_amendment_path,
        )
    elif args.command == "run":
        run_shard(
            args.study_dir,
            facets_exe=args.facets_exe.resolve(),
            shard_index=args.shard_index,
            shard_count=args.shard_count,
            resume=args.resume,
            timeout_seconds=args.timeout_seconds,
        )
    else:
        aggregate_study(args.study_dir, aggregate_name=args.aggregate_name)


if __name__ == "__main__":
    main()
