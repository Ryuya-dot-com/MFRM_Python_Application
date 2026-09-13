#!/usr/bin/env python3
"""Independent fixed-N replication of the registered H6 design effect.

The primary evidence is the independently completion-marked, native-precision
Python JMLE threshold recovery.  FACETS 4.5.0 is a separately gated calibration
layer and can never replace or invalidate an otherwise ready Python result.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import sys
import time
from typing import Any, Iterable

import numpy as np
import pandas as pd
from scipy.stats import kurtosis, skew, t as student_t


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from mfrm_app.operating_characteristics import deterministic_replicate_seed  # noqa: E402
from validation.estimand_distribution_study import (  # noqa: E402
    BASE_SEED,
    DESIGNS,
    FIT_MODEL,
    PERSON_SD,
    PYTHON_JMLE_MODE,
    THRESHOLD_CONDITION,
    _frame_digest,
    _latent_rows,
    generate_person_shapes,
)
from validation.facets_pcm_boundary_pilot import (  # noqa: E402
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
from validation.facets_resilient_pair import (  # noqa: E402
    dependency_manifest,
    dependency_manifest_digest,
    fit_resilient_pair,
)
from validation.operating_characteristics_facets import sha256_file  # noqa: E402


SCHEMA_VERSION = "mfrm-h6-design-replication-v1"
PLAN_PATH = REPO_ROOT / "validation" / "h6_design_replication_plan_20260811.json"
REGISTRATION_PATH = (
    REPO_ROOT / "validation" / "h6_design_replication_execution_registration_20260811.json"
)
FACETS_EXE_DEFAULT = Path(r"C:\Facets\Facets.exe")
INPUT_FILES = (
    "manifest.csv",
    "generated_ratings.csv",
    "generated_facet_truth.csv",
    "generated_anchors.csv",
    "generated_pcm_threshold_truth.csv",
    "attempt_manifest.csv",
)
PRIMARY_ID = "H6R_JMLE_NORMAL_THRESHOLD_RMSE_PLANNED_GT_COMPLETE"


def _json_dump(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _atomic_json(path: Path, value: Any) -> None:
    temporary = path.with_name(f"{path.name}.tmp.{os.getpid()}.{time.time_ns()}")
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


def _path_key(path: Path) -> str:
    return _portable_name(path)


def build_dependency_manifest(
    *,
    plan_path: Path = PLAN_PATH,
    registration_path: Path = REGISTRATION_PATH,
) -> dict[str, str]:
    """Extend the resilient-pair manifest with this registered execution layer."""

    manifest = dependency_manifest()
    additions = (Path(__file__).resolve(), plan_path.resolve(), registration_path.resolve())
    for path in additions:
        if not path.is_file():
            raise FileNotFoundError(f"Replication dependency missing: {path}")
        manifest[_path_key(path)] = sha256_file(path)
    return dict(sorted(manifest.items()))


def _validate_registration(plan_path: Path, registration_path: Path) -> dict[str, Any]:
    registration = json.loads(registration_path.read_text(encoding="utf-8"))
    expected = {
        "plan_sha256": sha256_file(plan_path),
        "runner_sha256": sha256_file(Path(__file__).resolve()),
        "resilient_pair_component_sha256": sha256_file(
            REPO_ROOT / "validation" / "facets_resilient_pair.py"
        ),
    }
    for key, value in expected.items():
        if str(registration.get(key, "")).lower() != value.lower():
            raise ValueError(f"Execution registration mismatch: {key}")
    if not bool(registration.get("tests_passed_before_registration", False)):
        raise ValueError("Execution registration does not assert passing tests")
    return registration


def _registered_range(plan: dict[str, Any]) -> tuple[int, int, int]:
    evidence = plan.get("evidence_separation", {})
    registered = evidence.get("replication_replicates")
    if not isinstance(registered, list) or len(registered) != 2:
        raise ValueError("Plan must register a two-element replication range")
    start, end = (int(value) for value in registered)
    count = int(evidence.get("replication_count", -1))
    if end - start + 1 != count:
        raise ValueError("Registered replication range/count mismatch")
    if bool(evidence.get("pool_previous_replicates_into_primary", True)):
        raise ValueError("Replication plan must prohibit pooling prior replicates")
    if not bool(evidence.get("fixed_sample_size", False)):
        raise ValueError("Replication plan must use a fixed sample size")
    return start, end, count


def generate_replication_bundle(
    *,
    replicates: int,
    replicate_start: int,
    base_seed: int,
    plan_path: Path = PLAN_PATH,
    registration_path: Path = REGISTRATION_PATH,
) -> dict[str, pd.DataFrame]:
    """Generate normal-person, paired complete/planned-connected inputs only."""

    plan_path = plan_path.resolve()
    registration_path = registration_path.resolve()
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    registered_start, registered_end, registered_count = _registered_range(plan)
    if (
        int(replicates) != registered_count
        or int(replicate_start) != registered_start
        or int(replicate_start) + int(replicates) - 1 != registered_end
    ):
        raise ValueError("Prepared range does not match the registered replication range")
    _validate_registration(plan_path, registration_path)
    dependencies = build_dependency_manifest(
        plan_path=plan_path, registration_path=registration_path
    )
    dependency_digest = dependency_manifest_digest(dependencies)

    manifests: list[dict[str, Any]] = []
    ratings_parts: list[pd.DataFrame] = []
    truth_parts: list[pd.DataFrame] = []
    threshold_parts: list[dict[str, Any]] = []
    for replicate in range(int(replicate_start), int(replicate_start) + int(replicates)):
        seed = deterministic_replicate_seed(int(base_seed), "person-shape-study", replicate)
        persons = generate_person_shapes(seed)["normal"]
        uniform_seed = deterministic_replicate_seed(int(base_seed), "row-uniforms", replicate)
        row_count = N_PERSONS * len(RATER_TRUTH) * len(TASK_TRUTH) * len(CRITERION_TRUTH)
        row_uniforms = np.random.default_rng(uniform_seed).random(row_count)
        latent = _latent_rows(persons, row_uniforms)
        complete = apply_threshold_condition(
            latent, THRESHOLD_CONDITIONS[THRESHOLD_CONDITION]
        )
        complete_hash = _frame_sha256(
            complete, ["Person", "Rater", "Task", "Criterion", "Score"]
        )
        uniform_hash = hashlib.sha256(row_uniforms.tobytes()).hexdigest()
        person_truth = {
            f"P{index:03d}": float(value)
            for index, value in enumerate(persons, start=1)
        }
        for design in DESIGNS:
            run_id = f"{design}__normal::rep-{replicate:05d}"
            run_ratings = apply_observation_design(complete, design)
            audit = constrained_adjacent_design_audit(run_ratings, FIT_MODEL)
            if int(audit["Nullity"]) != 0:
                raise RuntimeError(f"Prepared rank-deficient registered run: {run_id}")
            manifests.append({
                "RunId": run_id,
                "ConditionId": f"{design}__normal",
                "Design": design,
                "PersonDistribution": "normal",
                "TruthBias": 0.0,
                "Replicate": replicate,
                "Seed": int(seed),
                "UniformSeed": int(uniform_seed),
                "Categories": CATEGORIES,
                "ThresholdCondition": THRESHOLD_CONDITION,
                "PersonMeanRealized": float(np.mean(persons)),
                "PersonSDRealized": float(np.std(persons, ddof=0)),
                "PersonSkewRealized": float(skew(persons, bias=False)),
                "PersonExcessKurtosisRealized": float(
                    kurtosis(persons, fisher=True, bias=False)
                ),
                "ExpectedRows": int(len(run_ratings)),
                "ExpectedStructuralNullity": int(audit["Nullity"]),
                "PersonRaterComponents": int(audit["PersonRaterComponents"]),
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
                        "ConditionId": f"{design}__normal",
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
    for ordinal, row in enumerate(manifest.itertuples(index=False)):
        run_ratings = ratings[ratings["RunId"].astype(str).eq(str(row.RunId))]
        run_truth = truth[truth["RunId"].astype(str).eq(str(row.RunId))]
        run_thresholds = thresholds[thresholds["RunId"].astype(str).eq(str(row.RunId))]
        input_hash = hashlib.sha256(
            (
                _frame_digest(run_ratings)
                + _frame_digest(run_truth)
                + _frame_digest(run_thresholds)
            ).encode("ascii")
        ).hexdigest()
        attempt_id = f"{row.RunId}::RESILIENT_FACETS_PYTHON_JMLE_PCM"
        fingerprint = hashlib.sha256(
            f"{attempt_id}|{input_hash}|{dependency_digest}".encode("utf-8")
        ).hexdigest()
        attempts.append({
            "AttemptOrdinal": ordinal,
            "AttemptId": attempt_id,
            "RunId": row.RunId,
            "AttemptType": "RESILIENT_FACETS_PYTHON_JMLE_PCM",
            "Replicate": int(row.Replicate),
            "Design": row.Design,
            "PersonDistribution": "normal",
            "RunInputSHA256": input_hash,
            "DependencyManifestSHA256": dependency_digest,
            "AttemptFingerprint": fingerprint,
        })
    return {
        "manifest.csv": manifest,
        "generated_ratings.csv": ratings,
        "generated_facet_truth.csv": truth,
        "generated_anchors.csv": pd.DataFrame(
            columns=["RunId", "Facet", "Level", "Anchor"]
        ),
        "generated_pcm_threshold_truth.csv": thresholds,
        "attempt_manifest.csv": pd.DataFrame(attempts),
    }


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
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    start, end, count = _registered_range(plan)
    facets_exe = facets_exe.resolve()
    if not facets_exe.is_file():
        raise FileNotFoundError(f"FACETS executable missing: {facets_exe}")
    _validate_registration(plan_path.resolve(), registration_path.resolve())
    dependencies = build_dependency_manifest(
        plan_path=plan_path, registration_path=registration_path
    )
    bundle = generate_replication_bundle(
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
    input_hashes = {filename: sha256_file(input_dir / filename) for filename in INPUT_FILES}
    _json_dump(input_dir / "bundle_hashes.json", input_hashes)
    _json_dump(study_dir / "dependency_manifest.json", dependencies)
    identity = {
        "schema_version": SCHEMA_VERSION,
        "phase": "independent_registered_replication",
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
        "original_confirmatory_evidence_replaced": False,
    }
    _json_dump(study_dir / "study_identity.json", identity)
    return identity


def validate_study_identity(study_dir: Path) -> dict[str, Any]:
    study_dir = study_dir.resolve()
    identity = json.loads((study_dir / "study_identity.json").read_text(encoding="utf-8"))
    plan_path = _resolve_recorded_path(identity["plan_file"])
    registration_path = _resolve_recorded_path(identity["registration_file"])
    checks = {
        "plan_sha256": sha256_file(plan_path),
        "registration_sha256": sha256_file(registration_path),
        "runner_sha256": sha256_file(Path(__file__).resolve()),
    }
    for key, actual in checks.items():
        if str(identity.get(key, "")).lower() != actual.lower():
            raise ValueError(f"Study identity mismatch: {key}")
    _validate_registration(plan_path, registration_path)
    dependencies = build_dependency_manifest(
        plan_path=plan_path, registration_path=registration_path
    )
    if identity["dependency_manifest_sha256"] != dependency_manifest_digest(dependencies):
        raise ValueError("Dependency manifest changed after preparation")
    retained_manifest = json.loads(
        (study_dir / "dependency_manifest.json").read_text(encoding="utf-8")
    )
    if retained_manifest != dependencies:
        raise ValueError("Retained dependency manifest differs from current dependencies")
    if identity["facets_executable_sha256"] != sha256_file(Path(identity["facets_executable"])):
        raise ValueError("FACETS executable changed after preparation")
    for filename, expected in identity["retained_input_sha256"].items():
        if sha256_file(study_dir / "retained_input" / filename) != expected:
            raise ValueError(f"Retained input changed: {filename}")
    return identity


def select_shard_attempts(
    attempts: pd.DataFrame, *, shard_index: int, shard_count: int
) -> pd.DataFrame:
    if shard_count < 1 or shard_index < 0 or shard_index >= shard_count:
        raise ValueError("Require 0 <= shard-index < shard-count")
    ordinals = pd.to_numeric(attempts["AttemptOrdinal"], errors="raise").astype(int)
    if ordinals.duplicated().any():
        raise ValueError("Attempt ordinals must be unique")
    return attempts.loc[ordinals.mod(shard_count).eq(shard_index)].copy()


def _generation_hashes(generation_dir: Path) -> dict[str, str]:
    return {
        path.relative_to(generation_dir).as_posix(): sha256_file(path)
        for path in sorted(generation_dir.rglob("*"))
        if path.is_file() and not path.name.startswith("completion.json")
    }


def _validate_completion(
    completion: dict[str, Any], attempt: pd.Series, generation_dir: Path, identity: dict[str, Any]
) -> None:
    expected = {
        "attempt_id": str(attempt["AttemptId"]),
        "attempt_fingerprint": str(attempt["AttemptFingerprint"]),
        "run_input_sha256": str(attempt["RunInputSHA256"]),
        "dependency_manifest_sha256": identity["dependency_manifest_sha256"],
        "runner_sha256": identity["runner_sha256"],
        "facets_executable_sha256": identity["facets_executable_sha256"],
    }
    for key, value in expected.items():
        if str(completion.get(key, "")) != str(value):
            raise ValueError(f"Completion mismatch for {attempt['AttemptId']}: {key}")
    for relative, expected_hash in completion.get("artifact_sha256", {}).items():
        path = generation_dir / relative
        if not path.is_file() or sha256_file(path) != expected_hash:
            raise ValueError(f"Completion artifact mismatch: {path}")


def _attempt_dir(study_dir: Path, attempt: pd.Series) -> Path:
    return study_dir / "attempts" / f"{int(attempt['AttemptOrdinal']):05d}"


def _next_short_artifact_root(study_dir: Path, attempt: pd.Series) -> Path:
    """Keep legacy FACETS report paths safely below the Windows MAX_PATH boundary."""

    parent = study_dir / "w" / f"{int(attempt['AttemptOrdinal']):04x}"
    parent.mkdir(parents=True, exist_ok=True)
    generation = 0
    while (parent / f"g{generation:02x}").exists():
        generation += 1
    return parent / f"g{generation:02x}"


def run_shard(
    study_dir: Path,
    *,
    shard_index: int,
    shard_count: int,
    resume: bool,
    timeout_seconds: float,
) -> dict[str, int]:
    study_dir = study_dir.resolve()
    identity = validate_study_identity(study_dir)
    input_dir = study_dir / "retained_input"
    attempts = pd.read_csv(input_dir / "attempt_manifest.csv")
    selected = select_shard_attempts(
        attempts, shard_index=shard_index, shard_count=shard_count
    )
    manifest = pd.read_csv(input_dir / "manifest.csv").set_index("RunId", drop=False)
    ratings_all = pd.read_csv(input_dir / "generated_ratings.csv")
    truth_all = pd.read_csv(input_dir / "generated_facet_truth.csv")
    threshold_all = pd.read_csv(input_dir / "generated_pcm_threshold_truth.csv")
    (study_dir / "attempts").mkdir(exist_ok=True)
    (study_dir / "locks").mkdir(exist_ok=True)
    summary = {"assigned": int(len(selected)), "completed_now": 0, "skipped": 0, "failed": 0}
    for sequence, (_, attempt) in enumerate(selected.iterrows(), start=1):
        attempt_dir = _attempt_dir(study_dir, attempt)
        completion_path = attempt_dir / "completion.json"
        if completion_path.is_file():
            completion = json.loads(completion_path.read_text(encoding="utf-8"))
            artifact_root = study_dir / str(completion["retained_artifact_root"])
            _validate_completion(completion, attempt, artifact_root, identity)
            if not resume:
                raise FileExistsError(f"Attempt already completed: {attempt['AttemptId']}")
            summary["skipped"] += 1
            print(f"[{sequence}/{len(selected)}] skip {attempt['AttemptId']}", flush=True)
            continue
        attempt_dir.mkdir(parents=True, exist_ok=True)
        artifact_root = _next_short_artifact_root(study_dir, attempt)
        run_id = str(attempt["RunId"])
        row = manifest.loc[run_id]
        ratings = ratings_all[ratings_all["RunId"].astype(str).eq(run_id)].drop(columns="RunId")
        truth = truth_all[truth_all["RunId"].astype(str).eq(run_id)]
        thresholds = threshold_all[threshold_all["RunId"].astype(str).eq(run_id)]
        succeeded = False
        failure_reason = ""
        outcome: dict[str, Any] = {}
        try:
            outcome = fit_resilient_pair(
                row,
                ratings,
                truth,
                thresholds,
                facets_exe=Path(identity["facets_executable"]),
                output_dir=artifact_root,
                lock_path=study_dir / "locks" / "facets_global.lock",
                timeout_seconds=timeout_seconds,
            )
            succeeded = True
        except Exception as exc:  # retain the planned denominator and failed generation
            failure_reason = f"{type(exc).__name__}: {exc}"
            artifact_root.mkdir(parents=True, exist_ok=True)
            _json_dump(artifact_root / "unhandled_failure.json", {
                "attempt_id": str(attempt["AttemptId"]),
                "failure_reason": failure_reason,
            })
            summary["failed"] += 1
        artifact_hashes = _generation_hashes(artifact_root)
        completion = {
            "schema_version": SCHEMA_VERSION,
            "attempt_id": str(attempt["AttemptId"]),
            "attempt_fingerprint": str(attempt["AttemptFingerprint"]),
            "run_input_sha256": str(attempt["RunInputSHA256"]),
            "dependency_manifest_sha256": identity["dependency_manifest_sha256"],
            "runner_sha256": identity["runner_sha256"],
            "facets_executable_sha256": identity["facets_executable_sha256"],
            "attempt_execution_completed": succeeded,
            "statistical_evidence_ready": bool(outcome.get("statistical_evidence_ready", False)),
            "calibration_ready": bool(outcome.get("calibration_ready", False)),
            "failure_reason": failure_reason,
            "retained_artifact_root": artifact_root.relative_to(study_dir).as_posix(),
            "artifact_sha256": artifact_hashes,
        }
        _atomic_json(completion_path, completion)
        summary["completed_now"] += 1
        print(
            f"[{sequence}/{len(selected)}] complete {attempt['AttemptId']} "
            f"python={completion['statistical_evidence_ready']} "
            f"facets={completion['calibration_ready']}",
            flush=True,
        )
    print(json.dumps(summary, sort_keys=True), flush=True)
    return summary


def summarize_primary(
    python_thresholds: pd.DataFrame,
    *,
    required_pairs: int = 100,
    precision_target: float = 0.015,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Calculate the single registered paired RMSE contrast and descriptives."""

    eligible = python_thresholds[
        python_thresholds["EstimatorMode"].astype(str).eq(PYTHON_JMLE_MODE)
        & python_thresholds["IncludedInStudy"].astype(bool)
    ].copy()
    eligible["TruthError"] = pd.to_numeric(eligible["TruthError"], errors="coerce")
    counts = eligible.groupby(["RunId", "Replicate", "Design"], as_index=False).agg(
        ThresholdRows=("TruthError", "count"),
        ThresholdRMSE=("TruthError", lambda values: float(np.sqrt(np.mean(np.square(values))))),
    )
    valid_runs = counts[counts["ThresholdRows"].eq(6)]
    pivot = valid_runs.pivot(index="Replicate", columns="Design", values="ThresholdRMSE")
    if {"complete", "planned_connected"}.issubset(pivot.columns):
        paired = pivot[["complete", "planned_connected"]].dropna().reset_index()
        paired["ContrastPlannedMinusComplete"] = (
            paired["planned_connected"] - paired["complete"]
        )
    else:
        paired = pd.DataFrame(columns=[
            "Replicate", "complete", "planned_connected", "ContrastPlannedMinusComplete"
        ])
    values = pd.to_numeric(paired["ContrastPlannedMinusComplete"], errors="coerce").dropna()
    n = int(len(values))
    mean = float(values.mean()) if n else np.nan
    sd = float(values.std(ddof=1)) if n > 1 else np.nan
    se = sd / np.sqrt(n) if n > 1 else np.nan
    critical = float(student_t.ppf(0.975, n - 1)) if n > 1 else np.nan
    half_width = critical * se if n > 1 else np.nan
    if n > 1 and np.isfinite(se) and se > 0:
        statistic = mean / se
        p_one_sided = float(student_t.sf(statistic, n - 1))
    elif n > 1 and mean > 0:
        statistic, p_one_sided = np.inf, 0.0
    else:
        statistic, p_one_sided = np.nan, np.nan
    replicated = bool(
        n == int(required_pairs)
        and np.isfinite(mean)
        and mean > 0
        and np.isfinite(p_one_sided)
        and p_one_sided <= 0.05
    )
    precision_pass = bool(
        n == int(required_pairs) and np.isfinite(half_width) and half_width <= precision_target
    )
    result = pd.DataFrame([{
        "PrimaryId": PRIMARY_ID,
        "FinitePairs": n,
        "RequiredPairs": int(required_pairs),
        "MeanContrast": mean,
        "SDContrast": sd,
        "SEContrast": se,
        "TStatistic": statistic,
        "DegreesFreedom": n - 1 if n else 0,
        "OneSidedP": p_one_sided,
        "CI95Low": mean - half_width if np.isfinite(half_width) else np.nan,
        "CI95High": mean + half_width if np.isfinite(half_width) else np.nan,
        "CI95HalfWidth": half_width,
        "PrecisionTarget": float(precision_target),
        "ReplicationDecision": replicated,
        "PrecisionQualification": precision_pass,
        "MissingRuleSatisfied": n == int(required_pairs),
    }])
    element = eligible.pivot_table(
        index=["Replicate", "StepFacetLevel", "Category"],
        columns="Design",
        values="TruthError",
        aggfunc="first",
    ).reset_index()
    if {"complete", "planned_connected"}.issubset(element.columns):
        element["TruthErrorContrastPlannedMinusComplete"] = (
            element["planned_connected"] - element["complete"]
        )
        element_summary = element.groupby(
            ["StepFacetLevel", "Category"], as_index=False
        )["TruthErrorContrastPlannedMinusComplete"].agg(
            N="count", Mean="mean", SD="std"
        )
    else:
        element_summary = pd.DataFrame()
    return result, paired, element_summary


def aggregate_study(study_dir: Path) -> dict[str, Any]:
    study_dir = study_dir.resolve()
    identity = validate_study_identity(study_dir)
    attempts = pd.read_csv(study_dir / "retained_input" / "attempt_manifest.csv")
    run_parts: list[pd.DataFrame] = []
    threshold_parts: list[pd.DataFrame] = []
    outcomes: list[dict[str, Any]] = []
    marker_hashes: dict[str, str] = {}
    for _, attempt in attempts.iterrows():
        attempt_dir = _attempt_dir(study_dir, attempt)
        completion_path = attempt_dir / "completion.json"
        if not completion_path.is_file():
            raise FileNotFoundError(f"Missing completion marker: {attempt['AttemptId']}")
        completion = json.loads(completion_path.read_text(encoding="utf-8"))
        artifact_root = study_dir / str(completion["retained_artifact_root"])
        _validate_completion(completion, attempt, artifact_root, identity)
        marker_hashes[str(attempt["AttemptId"])] = sha256_file(completion_path)
        pair_dir = artifact_root
        if (pair_dir / "combined_run_ledger.csv").is_file():
            frame = pd.read_csv(pair_dir / "combined_run_ledger.csv")
            frame.insert(0, "AttemptId", str(attempt["AttemptId"]))
            run_parts.append(frame)
            frame = pd.read_csv(pair_dir / "combined_thresholds.csv")
            frame.insert(0, "AttemptId", str(attempt["AttemptId"]))
            threshold_parts.append(frame)
            outcome = json.loads(
                (pair_dir / "resilient_pair_metrics.json").read_text(encoding="utf-8")
            )
        else:
            outcome = {
                "run_id": str(attempt["RunId"]),
                "statistical_evidence_ready": False,
                "calibration_ready": False,
                "pair_fully_qualified": False,
                "facets_tries": 0,
                "facets_retry_count": 0,
                "facets_final_error": completion.get("failure_reason", ""),
            }
        outcome["AttemptId"] = str(attempt["AttemptId"])
        outcomes.append(outcome)
    runs = pd.concat(run_parts, ignore_index=True, sort=False) if run_parts else pd.DataFrame()
    thresholds = (
        pd.concat(threshold_parts, ignore_index=True, sort=False)
        if threshold_parts else pd.DataFrame()
    )
    outcome_frame = pd.DataFrame(outcomes)
    plan = json.loads(_resolve_recorded_path(identity["plan_file"]).read_text(encoding="utf-8"))
    primary_plan = plan["primary"]
    primary, paired, element_summary = summarize_primary(
        thresholds,
        required_pairs=int(primary_plan["required_finite_pairs"]),
        precision_target=float(primary_plan["precision_target_two_sided_95_half_width"]),
    )
    expected = int(identity["attempts"])
    python_ready = int(outcome_frame["statistical_evidence_ready"].astype(bool).sum())
    calibration_ready = int(outcome_frame["calibration_ready"].astype(bool).sum())
    facets_runs = runs[runs.get("EstimatorMode", pd.Series(dtype=str)).astype(str).eq("FACETS_4_5_JMLE")]
    numeric_maxima = {}
    for column in (
        "MainWeightedMAE",
        "MainMaxAbsDifference",
        "ThresholdWeightedMAE",
        "ThresholdMaxAbsDifference",
    ):
        values = pd.to_numeric(facets_runs.get(column, pd.Series(dtype=float)), errors="coerce")
        numeric_maxima[column] = float(values.max()) if values.notna().any() else np.nan
    min_spearman_values = pd.to_numeric(
        facets_runs.get("MinimumWithinFacetSpearman", pd.Series(dtype=float)), errors="coerce"
    )
    calibration = {
        "expected_pairs": expected,
        "python_statistical_evidence_ready": python_ready,
        "facets_calibration_ready": calibration_ready,
        "all_python_ready": python_ready == expected,
        "all_facets_calibration_ready": calibration_ready == expected,
        "pair_fully_qualified": int(outcome_frame["pair_fully_qualified"].astype(bool).sum()),
        "total_facets_tries": int(pd.to_numeric(outcome_frame["facets_tries"]).fillna(0).sum()),
        "total_facets_retries": int(
            pd.to_numeric(outcome_frame["facets_retry_count"]).fillna(0).sum()
        ),
        "minimum_within_facet_spearman": (
            float(min_spearman_values.min()) if min_spearman_values.notna().any() else np.nan
        ),
        **{f"maximum_{key}": value for key, value in numeric_maxima.items()},
    }
    calibration["strict_workbench_gate"] = bool(
        calibration["all_python_ready"] and calibration["all_facets_calibration_ready"]
    )
    aggregate_dir = study_dir / "aggregate"
    aggregate_dir.mkdir(exist_ok=True)
    runs.to_csv(aggregate_dir / "run_ledger.csv", index=False, lineterminator="\n")
    thresholds.to_csv(aggregate_dir / "thresholds.csv", index=False, lineterminator="\n")
    outcome_frame.to_csv(aggregate_dir / "attempt_outcomes.csv", index=False, lineterminator="\n")
    primary.to_csv(aggregate_dir / "primary_results.csv", index=False, lineterminator="\n")
    paired.to_csv(aggregate_dir / "paired_threshold_rmse.csv", index=False, lineterminator="\n")
    element_summary.to_csv(
        aggregate_dir / "threshold_element_descriptives.csv", index=False, lineterminator="\n"
    )
    marker_digest = hashlib.sha256(
        json.dumps(marker_hashes, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    metrics = {
        "schema_version": SCHEMA_VERSION,
        "primary": primary.iloc[0].to_dict(),
        "facets_calibration_secondary": calibration,
        "completion_markers": len(marker_hashes),
        "completion_marker_set_sha256": marker_digest,
        "prior_replicates_pooled": False,
        "original_confirmatory_evidence_replaced": False,
        "claim_limit": plan["claim_limit"],
    }
    _json_dump(aggregate_dir / "replication_metrics.json", metrics)
    aggregate_files = (
        "run_ledger.csv",
        "thresholds.csv",
        "attempt_outcomes.csv",
        "primary_results.csv",
        "paired_threshold_rmse.csv",
        "threshold_element_descriptives.csv",
        "replication_metrics.json",
    )
    aggregate_identity = {
        "schema_version": f"{SCHEMA_VERSION}-aggregate-identity",
        "study_identity_sha256": sha256_file(study_dir / "study_identity.json"),
        "completion_marker_set_sha256": marker_digest,
        "artifact_sha256": {
            filename: sha256_file(aggregate_dir / filename) for filename in aggregate_files
        },
    }
    _atomic_json(aggregate_dir / "aggregate_identity.json", aggregate_identity)
    return metrics


def preexecution_audit(study_dir: Path) -> dict[str, Any]:
    study_dir = study_dir.resolve()
    identity = validate_study_identity(study_dir)
    input_dir = study_dir / "retained_input"
    manifest = pd.read_csv(input_dir / "manifest.csv")
    ratings = pd.read_csv(input_dir / "generated_ratings.csv")
    attempts = pd.read_csv(input_dir / "attempt_manifest.csv")
    row_counts = ratings.groupby("RunId").size()
    expected_rows = manifest.set_index("RunId")["ExpectedRows"].astype(int)
    longest_report = study_dir / (
        "w/00c7/g00/facets_attempts/pair_try_03/facets_runs/"
        "planned_connected__normal_rep-00220__pcm/report_u6.txt"
    )
    checks = {
        "identity_valid": True,
        "registered_replicates_exact": int(manifest["Replicate"].nunique()) == 100,
        "replicate_range_exact": [int(manifest["Replicate"].min()), int(manifest["Replicate"].max())]
        == [121, 220],
        "normal_only": set(manifest["PersonDistribution"].astype(str)) == {"normal"},
        "two_registered_designs_only": set(manifest["Design"].astype(str)) == set(DESIGNS),
        "dataset_count_200": len(manifest) == 200,
        "attempt_count_200": len(attempts) == 200,
        "row_counts_match": row_counts.reindex(expected_rows.index).astype(int).equals(
            expected_rows
        ),
        "all_structurally_full_rank": manifest["ExpectedStructuralNullity"].eq(0).all(),
        "realized_mean_contract": manifest["PersonMeanRealized"].abs().max() <= 1e-12,
        "realized_sd_contract": (manifest["PersonSDRealized"] - PERSON_SD).abs().max() <= 1e-12,
        "no_previous_replicates": manifest["Replicate"].min() > 120,
        "facets_executable_hash_bound": bool(identity["facets_executable_sha256"]),
        "facets_worst_case_report_path_below_240": len(str(longest_report)) < 240,
    }
    audit = {
        "schema_version": f"{SCHEMA_VERSION}-preexecution-audit",
        "checks": {key: bool(value) for key, value in checks.items()},
        "all_checks_pass": all(bool(value) for value in checks.values()),
        "study_identity_sha256": sha256_file(study_dir / "study_identity.json"),
        "facets_worst_case_report_path_length": len(str(longest_report)),
    }
    _json_dump(study_dir / "preexecution_audit.json", audit)
    return audit


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
