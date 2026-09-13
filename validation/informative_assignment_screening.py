#!/usr/bin/env python3
"""Prepare, execute, resume, and aggregate the registered assignment screen."""

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
from validation.estimand_bridge_pilot import CMLE_MODE  # noqa: E402
from validation.estimand_distribution_study import (  # noqa: E402
    BASE_SEED,
    FACETS_MODE,
    FIT_MODEL,
    MML_FIXED_MODE,
    MML_FREE_MODE,
    PERSON_SD,
    PYTHON_JMLE_MODE,
    THRESHOLD_CONDITION,
    _frame_digest,
    _latent_rows,
    _run_one_attempt,
    generate_person_shapes,
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
from validation.informative_assignment_design import (  # noqa: E402
    DESIGNS,
    apply_informative_assignment_design,
    informative_assignment_audit,
)
from validation.operating_characteristics_facets import sha256_file  # noqa: E402


SCHEMA_VERSION = "mfrm-informative-assignment-screening-v1"
PLAN_PATH = REPO_ROOT / "validation" / "informative_assignment_screening_plan_20260811.json"
REGISTRATION_PATH = (
    REPO_ROOT
    / "validation"
    / "informative_assignment_screening_execution_registration_20260811.json"
)
FACETS_EXE_DEFAULT = Path(r"C:\Facets\Facets.exe")
RUN_DESIGN_SLUG = {
    "complete": "complete",
    "planned_connected": "planned",
    "ability_severity_aligned_connected": "aligned",
}
ATTEMPT_TYPES = (
    "RESILIENT_FACETS_PYTHON_JMLE_PCM",
    "PYTHON_MML_FIXED_SD08_Q31_PCM",
    "PYTHON_MML_FREE_SD_Q31_PCM",
    "PYTHON_EXACT_CMLE_PCM",
)
SCIENTIFIC_MODES = (PYTHON_JMLE_MODE, MML_FIXED_MODE, MML_FREE_MODE, CMLE_MODE)
INPUT_FILES = (
    "manifest.csv",
    "generated_ratings.csv",
    "generated_facet_truth.csv",
    "generated_anchors.csv",
    "generated_pcm_threshold_truth.csv",
    "assignment_audit.csv",
    "attempt_manifest.csv",
)


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


def build_dependency_manifest(
    *, plan_path: Path = PLAN_PATH, registration_path: Path = REGISTRATION_PATH
) -> dict[str, str]:
    manifest = dependency_manifest()
    additions = (
        Path(__file__).resolve(),
        REPO_ROOT / "validation" / "informative_assignment_design.py",
        plan_path.resolve(),
        registration_path.resolve(),
    )
    for path in additions:
        if not path.is_file():
            raise FileNotFoundError(f"Screening dependency missing: {path}")
        manifest[_portable_name(path)] = sha256_file(path)
    return dict(sorted(manifest.items()))


def _validate_registration(plan_path: Path, registration_path: Path) -> dict[str, Any]:
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


def _registered_range(plan: dict[str, Any]) -> tuple[int, int, int]:
    evidence = plan["evidence_status"]
    start, end = (int(value) for value in evidence["replicate_range"])
    count = int(evidence["replicates"])
    if end - start + 1 != count:
        raise ValueError("Registered replicate range/count mismatch")
    if bool(evidence.get("pool_replicates_1_through_220", True)):
        raise ValueError("Screening plan must prohibit pooling previous replicates")
    if not bool(evidence.get("fixed_sample_size", False)):
        raise ValueError("Screening plan must use fixed sample size")
    if bool(evidence.get("confirmatory_claims_allowed", True)):
        raise ValueError("Screening plan must prohibit confirmatory claims")
    return start, end, count


def generate_screening_bundle(
    *,
    replicates: int,
    replicate_start: int,
    base_seed: int,
    plan_path: Path = PLAN_PATH,
    registration_path: Path = REGISTRATION_PATH,
) -> dict[str, pd.DataFrame]:
    plan_path = plan_path.resolve()
    registration_path = registration_path.resolve()
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    start, end, count = _registered_range(plan)
    if (
        int(replicates) != count
        or int(replicate_start) != start
        or int(replicate_start) + int(replicates) - 1 != end
    ):
        raise ValueError("Prepared range does not match registered screening range")
    _validate_registration(plan_path, registration_path)
    dependencies = build_dependency_manifest(
        plan_path=plan_path, registration_path=registration_path
    )
    dependency_digest = dependency_manifest_digest(dependencies)

    manifests: list[dict[str, Any]] = []
    ratings_parts: list[pd.DataFrame] = []
    truth_parts: list[pd.DataFrame] = []
    threshold_parts: list[dict[str, Any]] = []
    audit_parts: list[dict[str, Any]] = []
    for replicate in range(int(replicate_start), int(replicate_start) + int(replicates)):
        seed = deterministic_replicate_seed(int(base_seed), "person-shape-study", replicate)
        persons = generate_person_shapes(seed)["normal"]
        uniform_seed = deterministic_replicate_seed(int(base_seed), "row-uniforms", replicate)
        row_count = N_PERSONS * len(RATER_TRUTH) * len(TASK_TRUTH) * len(CRITERION_TRUTH)
        row_uniforms = np.random.default_rng(uniform_seed).random(row_count)
        latent = _latent_rows(persons, row_uniforms)
        complete_scores = apply_threshold_condition(
            latent, THRESHOLD_CONDITIONS[THRESHOLD_CONDITION]
        )
        theta = latent[["Person", "Theta"]].drop_duplicates("Person")
        complete_for_design = complete_scores.merge(
            theta, on="Person", how="left", validate="many_to_one"
        )
        complete_hash = _frame_sha256(
            complete_scores, ["Person", "Rater", "Task", "Criterion", "Score"]
        )
        uniform_hash = hashlib.sha256(row_uniforms.tobytes()).hexdigest()
        person_truth = {
            f"P{index:03d}": float(value)
            for index, value in enumerate(persons, start=1)
        }
        for design in DESIGNS:
            run_id = f"{RUN_DESIGN_SLUG[design]}__normal::rep-{replicate:05d}"
            audit = informative_assignment_audit(complete_for_design, design)
            if not audit["AllInvariantsPass"]:
                raise RuntimeError(f"Prepared design invariants failed: {run_id}")
            run_ratings = apply_informative_assignment_design(
                complete_for_design, design
            ).drop(columns="Theta")
            category_counts = (
                run_ratings.groupby(["Criterion", "Score"])
                .size()
                .reindex(
                    pd.MultiIndex.from_product(
                        [list(CRITERION_TRUTH), range(CATEGORIES)],
                        names=["Criterion", "Score"],
                    ),
                    fill_value=0,
                )
            )
            complete_category_support = bool(category_counts.gt(0).all())
            if not complete_category_support:
                raise RuntimeError(f"Prepared category support failed: {run_id}")
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
                "ExpectedStructuralNullity": int(audit["ConstrainedPCMNullity"]),
                "PersonRaterComponents": int(audit["PersonRaterComponents"]),
                "ThetaAssignedSeveritySpearman": audit["ThetaAssignedSeveritySpearman"],
                "MinimumCriterionCategoryCount": int(category_counts.min()),
                "CompleteCriterionCategorySupport": complete_category_support,
                "SharedUniformSHA256": uniform_hash,
                "CompleteResponseSHA256": complete_hash,
            })
            for key, passed in audit["Invariants"].items():
                audit_parts.append({
                    "RunId": run_id,
                    "Replicate": replicate,
                    "Design": design,
                    "AuditType": "Invariant",
                    "AuditKey": key,
                    "AuditValue": bool(passed),
                    "AuditJSON": "",
                })
            audit_parts.append({
                "RunId": run_id,
                "Replicate": replicate,
                "Design": design,
                "AuditType": "FullAudit",
                "AuditKey": "InformativeAssignmentAudit",
                "AuditValue": True,
                "AuditJSON": json.dumps(audit, sort_keys=True),
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
    ordinal = 0
    for row in manifest.itertuples(index=False):
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
        for attempt_type in ATTEMPT_TYPES:
            attempt_id = f"{row.RunId}::{attempt_type}"
            fingerprint = hashlib.sha256(
                f"{attempt_id}|{input_hash}|{dependency_digest}".encode("utf-8")
            ).hexdigest()
            attempts.append({
                "AttemptOrdinal": ordinal,
                "AttemptId": attempt_id,
                "RunId": row.RunId,
                "AttemptType": attempt_type,
                "Replicate": int(row.Replicate),
                "Design": row.Design,
                "PersonDistribution": "normal",
                "RunInputSHA256": input_hash,
                "DependencyManifestSHA256": dependency_digest,
                "AttemptFingerprint": fingerprint,
            })
            ordinal += 1
    return {
        "manifest.csv": manifest,
        "generated_ratings.csv": ratings,
        "generated_facet_truth.csv": truth,
        "generated_anchors.csv": pd.DataFrame(
            columns=["RunId", "Facet", "Level", "Anchor"]
        ),
        "generated_pcm_threshold_truth.csv": thresholds,
        "assignment_audit.csv": pd.DataFrame(audit_parts),
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
    bundle = generate_screening_bundle(
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
        "phase": "registered_screening",
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
        "confirmatory_claims_allowed": False,
        "claim_limit": plan["claim_limit"],
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
    retained = json.loads((study_dir / "dependency_manifest.json").read_text(encoding="utf-8"))
    if retained != dependencies:
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


def _next_artifact_root(study_dir: Path, attempt: pd.Series) -> Path:
    parent = study_dir / "w" / f"{int(attempt['AttemptOrdinal']):04x}"
    parent.mkdir(parents=True, exist_ok=True)
    generation = 0
    while (parent / f"g{generation:02x}").exists():
        generation += 1
    return parent / f"g{generation:02x}"


def _artifact_hashes(root: Path) -> dict[str, str]:
    return {
        path.relative_to(root).as_posix(): sha256_file(path)
        for path in sorted(root.rglob("*"))
        if path.is_file()
    }


def _attempt_dir(study_dir: Path, attempt: pd.Series) -> Path:
    return study_dir / "attempts" / f"{int(attempt['AttemptOrdinal']):05d}"


def _validate_completion(
    completion: dict[str, Any], attempt: pd.Series, artifact_root: Path, identity: dict[str, Any]
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
        path = artifact_root / relative
        if not path.is_file() or sha256_file(path) != expected_hash:
            raise ValueError(f"Completion artifact mismatch: {path}")


def _write_standard_artifacts(
    root: Path,
    runs: pd.DataFrame,
    recovery: pd.DataFrame,
    thresholds: pd.DataFrame,
    constraints: list[dict[str, Any]] | pd.DataFrame,
) -> None:
    runs.to_csv(root / "run_ledger.csv", index=False, lineterminator="\n")
    recovery.to_csv(root / "recovery.csv", index=False, lineterminator="\n")
    thresholds.to_csv(root / "thresholds.csv", index=False, lineterminator="\n")
    pd.DataFrame(constraints).to_csv(root / "constraints.csv", index=False, lineterminator="\n")


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
        artifact_root = _next_artifact_root(study_dir, attempt)
        run_id = str(attempt["RunId"])
        row = manifest.loc[run_id]
        ratings = ratings_all[ratings_all["RunId"].astype(str).eq(run_id)].drop(columns="RunId")
        truth = truth_all[truth_all["RunId"].astype(str).eq(run_id)]
        thresholds = threshold_all[threshold_all["RunId"].astype(str).eq(run_id)]
        execution_completed = False
        statistical_ready = False
        calibration_ready: bool | None = None
        failure_reason = ""
        try:
            if str(attempt["AttemptType"]) == "RESILIENT_FACETS_PYTHON_JMLE_PCM":
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
                runs = pd.read_csv(artifact_root / "combined_run_ledger.csv")
                recovery = pd.read_csv(artifact_root / "combined_recovery.csv")
                threshold_rows = pd.read_csv(artifact_root / "combined_thresholds.csv")
                _write_standard_artifacts(
                    artifact_root, runs, recovery, threshold_rows, []
                )
                statistical_ready = bool(outcome["statistical_evidence_ready"])
                calibration_ready = bool(outcome["calibration_ready"])
            else:
                artifact_root.mkdir(parents=True)
                runs, recovery, threshold_rows, constraints = _run_one_attempt(
                    attempt,
                    row,
                    ratings,
                    truth,
                    thresholds,
                    facets_exe=Path(identity["facets_executable"]),
                    stage_dir=artifact_root,
                    timeout_seconds=timeout_seconds,
                )
                _write_standard_artifacts(
                    artifact_root, runs, recovery, threshold_rows, constraints
                )
                statistical_ready = bool(runs["IncludedInStudy"].fillna(False).astype(bool).all())
            execution_completed = True
        except Exception as exc:  # retain planned denominator and work generation
            failure_reason = f"{type(exc).__name__}: {exc}"
            artifact_root.mkdir(parents=True, exist_ok=True)
            _json_dump(artifact_root / "unhandled_failure.json", {
                "attempt_id": str(attempt["AttemptId"]),
                "failure_reason": failure_reason,
            })
            summary["failed"] += 1
        completion = {
            "schema_version": SCHEMA_VERSION,
            "attempt_id": str(attempt["AttemptId"]),
            "attempt_type": str(attempt["AttemptType"]),
            "attempt_fingerprint": str(attempt["AttemptFingerprint"]),
            "run_input_sha256": str(attempt["RunInputSHA256"]),
            "dependency_manifest_sha256": identity["dependency_manifest_sha256"],
            "runner_sha256": identity["runner_sha256"],
            "facets_executable_sha256": identity["facets_executable_sha256"],
            "attempt_execution_completed": execution_completed,
            "statistical_evidence_ready": statistical_ready,
            "facets_calibration_ready": calibration_ready,
            "failure_reason": failure_reason,
            "retained_artifact_root": artifact_root.relative_to(study_dir).as_posix(),
            "artifact_sha256": _artifact_hashes(artifact_root),
        }
        _atomic_json(completion_path, completion)
        summary["completed_now"] += 1
        print(
            f"[{sequence}/{len(selected)}] complete {attempt['AttemptId']} "
            f"ready={statistical_ready} calibration={calibration_ready}",
            flush=True,
        )
    print(json.dumps(summary, sort_keys=True), flush=True)
    return summary


def screening_contrast_summary(contrasts: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    keys = ["EstimatorMode", "RecoveryDomain", "Metric"]
    for group_keys, group in contrasts.groupby(keys, dropna=False):
        values = pd.to_numeric(group["ContrastAlignedMinusPlanned"], errors="coerce").dropna()
        n = int(len(values))
        mean = float(values.mean()) if n else np.nan
        sd = float(values.std(ddof=1)) if n > 1 else np.nan
        se = sd / np.sqrt(n) if n > 1 else np.nan
        half = float(student_t.ppf(0.975, n - 1)) * se if n > 1 else np.nan
        rows.append({
            "EstimatorMode": group_keys[0],
            "RecoveryDomain": group_keys[1],
            "Metric": group_keys[2],
            "FinitePairs": n,
            "MeanContrastAlignedMinusPlanned": mean,
            "SDContrast": sd,
            "SEContrast": se,
            "ScreeningCI95Low": mean - half if np.isfinite(half) else np.nan,
            "ScreeningCI95High": mean + half if np.isfinite(half) else np.nan,
            "ScreeningCI95HalfWidth": half,
            "IntervalExcludesZero": bool(
                np.isfinite(half) and (mean - half > 0 or mean + half < 0)
            ),
            "ConfirmatoryClaimAllowed": False,
        })
    return pd.DataFrame(rows)


def _read_csv_if_nonempty(path: Path) -> pd.DataFrame:
    if not path.is_file() or path.stat().st_size <= 1:
        return pd.DataFrame()
    try:
        return pd.read_csv(path)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def aggregate_study(study_dir: Path) -> dict[str, Any]:
    study_dir = study_dir.resolve()
    identity = validate_study_identity(study_dir)
    input_dir = study_dir / "retained_input"
    attempts = pd.read_csv(input_dir / "attempt_manifest.csv")
    run_parts: list[pd.DataFrame] = []
    recovery_parts: list[pd.DataFrame] = []
    threshold_parts: list[pd.DataFrame] = []
    constraint_parts: list[pd.DataFrame] = []
    outcomes: list[dict[str, Any]] = []
    marker_hashes: dict[str, str] = {}
    for _, attempt in attempts.iterrows():
        attempt_dir = _attempt_dir(study_dir, attempt)
        marker = attempt_dir / "completion.json"
        if not marker.is_file():
            raise FileNotFoundError(f"Missing completion marker: {attempt['AttemptId']}")
        completion = json.loads(marker.read_text(encoding="utf-8"))
        artifact_root = study_dir / str(completion["retained_artifact_root"])
        _validate_completion(completion, attempt, artifact_root, identity)
        marker_hashes[str(attempt["AttemptId"])] = sha256_file(marker)
        for filename, target in (
            ("run_ledger.csv", run_parts),
            ("recovery.csv", recovery_parts),
            ("thresholds.csv", threshold_parts),
            ("constraints.csv", constraint_parts),
        ):
            frame = _read_csv_if_nonempty(artifact_root / filename)
            if len(frame):
                frame.insert(0, "AttemptId", str(attempt["AttemptId"]))
                target.append(frame)
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
    runs = pd.concat(run_parts, ignore_index=True, sort=False)
    recovery = pd.concat(recovery_parts, ignore_index=True, sort=False)
    thresholds = pd.concat(threshold_parts, ignore_index=True, sort=False)
    constraints = (
        pd.concat(constraint_parts, ignore_index=True, sort=False)
        if constraint_parts else pd.DataFrame()
    )
    outcome_frame = pd.DataFrame(outcomes)

    eligible_recovery = recovery[
        recovery["EstimatorMode"].astype(str).isin(SCIENTIFIC_MODES)
        & recovery["IncludedInStudy"].fillna(False).astype(bool)
        & recovery["Facet"].astype(str).isin({"Rater", "Task", "Criterion"})
    ].copy()
    eligible_recovery["ErrorAligned"] = pd.to_numeric(
        eligible_recovery["ErrorAligned"], errors="coerce"
    )
    facet_loss = eligible_recovery.groupby(
        ["RunId", "Replicate", "Design", "EstimatorMode", "Facet"], as_index=False
    )["ErrorAligned"].agg(
        RMSE=lambda values: float(np.sqrt(np.mean(np.square(values)))),
        MAE=lambda values: float(np.mean(np.abs(values))),
    )
    facet_long = facet_loss.melt(
        id_vars=["RunId", "Replicate", "Design", "EstimatorMode", "Facet"],
        value_vars=["RMSE", "MAE"],
        var_name="Metric",
        value_name="Loss",
    )
    facet_long["RecoveryDomain"] = "Facet:" + facet_long["Facet"].astype(str)

    eligible_thresholds = thresholds[
        thresholds["EstimatorMode"].astype(str).isin(SCIENTIFIC_MODES)
        & thresholds["IncludedInStudy"].fillna(False).astype(bool)
    ].copy()
    eligible_thresholds["TruthError"] = pd.to_numeric(
        eligible_thresholds["TruthError"], errors="coerce"
    )
    threshold_loss = eligible_thresholds.groupby(
        ["RunId", "Replicate", "Design", "EstimatorMode"], as_index=False
    )["TruthError"].agg(
        ThresholdRows="count",
        RMSE=lambda values: float(np.sqrt(np.mean(np.square(values)))),
        MAE=lambda values: float(np.mean(np.abs(values))),
    )
    threshold_loss = threshold_loss[threshold_loss["ThresholdRows"].eq(6)]
    threshold_long = threshold_loss.melt(
        id_vars=["RunId", "Replicate", "Design", "EstimatorMode"],
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
        sort=False,
    )
    planned = run_loss[run_loss["Design"].eq("planned_connected")].drop(
        columns=["RunId", "Design"]
    )
    aligned = run_loss[
        run_loss["Design"].eq("ability_severity_aligned_connected")
    ].drop(columns=["RunId", "Design"])
    contrasts = aligned.merge(
        planned,
        on=["Replicate", "EstimatorMode", "RecoveryDomain", "Metric"],
        suffixes=("Aligned", "Planned"),
        validate="one_to_one",
    )
    contrasts["ContrastAlignedMinusPlanned"] = (
        contrasts["LossAligned"] - contrasts["LossPlanned"]
    )
    screening = screening_contrast_summary(contrasts)

    rater_errors = eligible_recovery[eligible_recovery["Facet"].eq("Rater")]
    rater_planned = rater_errors[rater_errors["Design"].eq("planned_connected")][
        ["Replicate", "EstimatorMode", "Level", "ErrorAligned"]
    ]
    rater_aligned = rater_errors[
        rater_errors["Design"].eq("ability_severity_aligned_connected")
    ][["Replicate", "EstimatorMode", "Level", "ErrorAligned"]]
    rater_contrasts = rater_aligned.merge(
        rater_planned,
        on=["Replicate", "EstimatorMode", "Level"],
        suffixes=("Aligned", "Planned"),
        validate="one_to_one",
    )
    rater_contrasts["ErrorContrastAlignedMinusPlanned"] = (
        rater_contrasts["ErrorAlignedAligned"] - rater_contrasts["ErrorAlignedPlanned"]
    )
    rater_summary = rater_contrasts.groupby(
        ["EstimatorMode", "Level"], as_index=False
    )["ErrorContrastAlignedMinusPlanned"].agg(N="count", Mean="mean", SD="std")

    run_context = run_loss.groupby(
        ["EstimatorMode", "Design", "RecoveryDomain", "Metric"], as_index=False
    )["Loss"].agg(N="count", Mean="mean", SD="std")
    free_sd = runs[runs["EstimatorMode"].astype(str).eq(MML_FREE_MODE)].copy()
    free_sd_pivot = free_sd.pivot(
        index="Replicate", columns="Design", values="EstimatedPopulationSD"
    ).reset_index()
    if {"planned_connected", "ability_severity_aligned_connected"}.issubset(
        free_sd_pivot.columns
    ):
        free_sd_pivot["ContrastAlignedMinusPlanned"] = (
            free_sd_pivot["ability_severity_aligned_connected"]
            - free_sd_pivot["planned_connected"]
        )

    aggregate_dir = study_dir / "aggregate"
    if aggregate_dir.exists():
        raise FileExistsError(f"Aggregate directory already exists: {aggregate_dir}")
    aggregate_dir.mkdir()
    outputs = {
        "attempt_outcomes.csv": outcome_frame,
        "run_ledger.csv": runs,
        "recovery.csv": recovery,
        "thresholds.csv": thresholds,
        "constraints.csv": constraints,
        "run_loss.csv": run_loss,
        "paired_design_contrasts.csv": contrasts,
        "screening_contrast_summary.csv": screening,
        "rater_level_contrasts.csv": rater_contrasts,
        "rater_level_summary.csv": rater_summary,
        "design_context_summary.csv": run_context,
        "free_sd_paired.csv": free_sd_pivot,
    }
    for filename, frame in outputs.items():
        frame.to_csv(aggregate_dir / filename, index=False, lineterminator="\n")

    jmle_outcomes = outcome_frame[
        outcome_frame["AttemptType"].eq("RESILIENT_FACETS_PYTHON_JMLE_PCM")
    ]
    mml_outcomes = outcome_frame[outcome_frame["AttemptType"].str.contains("MML")]
    cmle_outcomes = outcome_frame[outcome_frame["AttemptType"].eq("PYTHON_EXACT_CMLE_PCM")]
    calibration_ready = jmle_outcomes["FACETSCalibrationReady"].fillna(False).astype(bool)
    constraints_ready = bool(
        len(constraints) == 180
        and constraints["ConstraintPass"].fillna(False).astype(bool).all()
    )
    gates = {
        "completion_markers_240": len(marker_hashes) == 240,
        "all_attempts_execution_completed": outcome_frame["ExecutionCompleted"].all(),
        "all_statistical_evidence_ready": outcome_frame["StatisticalEvidenceReady"].all(),
        "facets_calibration_60_of_60": len(jmle_outcomes) == 60 and calibration_ready.all(),
        "mml_120_of_120": len(mml_outcomes) == 120
        and mml_outcomes["StatisticalEvidenceReady"].all(),
        "cmle_60_of_60": len(cmle_outcomes) == 60
        and cmle_outcomes["StatisticalEvidenceReady"].all(),
        "constraints_180_of_180": constraints_ready,
        "screening_pairs_20_each": len(screening) == 32
        and screening["FinitePairs"].eq(20).all(),
    }
    numeric = lambda column: pd.to_numeric(jmle_outcomes[column], errors="coerce")
    facets_runs = runs[runs["EstimatorMode"].astype(str).eq(FACETS_MODE)].copy()
    facets_numeric = lambda column: pd.to_numeric(facets_runs[column], errors="coerce")
    calibration = {
        "pairs": int(len(jmle_outcomes)),
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
    metrics = {
        "schema_version": SCHEMA_VERSION,
        "phase": "screening",
        "datasets": int(identity["datasets"]),
        "attempts": int(identity["attempts"]),
        "gates": {key: bool(value) for key, value in gates.items()},
        "qualification_pass": all(bool(value) for value in gates.values()),
        "facets_calibration": calibration,
        "completion_marker_set_sha256": marker_digest,
        "confirmatory_claims_allowed": False,
        "estimator_ranking_constructed": False,
        "cross_basis_likelihood_comparison": "prohibited",
        "claim_limit": identity["claim_limit"],
    }
    _json_dump(aggregate_dir / "screening_metrics.json", metrics)
    artifact_names = tuple(outputs) + ("screening_metrics.json",)
    _atomic_json(aggregate_dir / "aggregate_identity.json", {
        "schema_version": f"{SCHEMA_VERSION}-aggregate-identity",
        "study_identity_sha256": sha256_file(study_dir / "study_identity.json"),
        "completion_marker_set_sha256": marker_digest,
        "artifact_sha256": {
            filename: sha256_file(aggregate_dir / filename) for filename in artifact_names
        },
    })
    return metrics


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
        "w/00ef/g00/facets_attempts/pair_try_03/facets_runs/"
        "aligned__normal_rep-00240__pcm/report_u6.txt"
    )
    sparse = manifest[manifest["Design"].ne("complete")]
    aligned = manifest[manifest["Design"].eq("ability_severity_aligned_connected")]
    checks = {
        "identity_valid": True,
        "registered_replicates_20": manifest["Replicate"].nunique() == 20,
        "replicate_range_221_240": [manifest["Replicate"].min(), manifest["Replicate"].max()]
        == [221, 240],
        "datasets_60": len(manifest) == 60,
        "attempts_240": len(attempts) == 240,
        "designs_exact": set(manifest["Design"].astype(str)) == set(DESIGNS),
        "row_counts_match": row_counts.reindex(expected_rows.index).astype(int).equals(
            expected_rows
        ),
        "sparse_rows_960": sparse["ExpectedRows"].eq(960).all(),
        "all_nullity_zero": manifest["ExpectedStructuralNullity"].eq(0).all(),
        "all_one_component": manifest["PersonRaterComponents"].eq(1).all(),
        "aligned_assignment_strong": aligned["ThetaAssignedSeveritySpearman"].gt(0.9).all(),
        "complete_category_support": manifest["CompleteCriterionCategorySupport"].all(),
        "positive_category_counts": manifest["MinimumCriterionCategoryCount"].gt(0).all(),
        "shared_complete_response_within_replicate": manifest.groupby("Replicate")[
            "CompleteResponseSHA256"
        ].nunique().eq(1).all(),
        "assignment_invariants_all_pass": assignment["AuditValue"].fillna(False).astype(bool).all(),
        "realized_mean_zero": manifest["PersonMeanRealized"].abs().max() <= 1e-12,
        "realized_sd_08": (manifest["PersonSDRealized"] - PERSON_SD).abs().max() <= 1e-12,
        "no_prior_replicates": manifest["Replicate"].min() > 220,
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
