#!/usr/bin/env python3
"""Prepare, execute, and aggregate the frozen known-assignment confirmation.

The frozen preflight engine is parameterized in a context manager rather than
edited.  This preserves the hashes of the retained preflight while reusing its
data-generation and attempt-execution machinery.
"""

from __future__ import annotations

import argparse
from contextlib import contextmanager
import copy
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
from typing import Any, Iterable, Iterator

import numpy as np
import pandas as pd
from scipy.stats import chi2, t as student_t


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from validation import known_assignment_multivector_preflight as engine  # noqa: E402
from validation import known_assignment_mml_crossfit as mml_crossfit  # noqa: E402
from validation.estimand_bridge_pilot import CMLE_MODE  # noqa: E402
from validation.estimand_distribution_study import (  # noqa: E402
    FACETS_MODE,
    MML_FIXED_MODE,
    MML_FREE_MODE,
    PYTHON_JMLE_MODE,
)
from validation.facets_resilient_pair import (  # noqa: E402
    dependency_manifest,
    dependency_manifest_digest,
)
from validation.operating_characteristics_facets import sha256_file  # noqa: E402


SCHEMA_VERSION = "known_assignment_confirmatory_v1"
PLAN_PATH = ROOT / "validation" / "known_assignment_confirmatory_plan_20260811.json"
REGISTRATION_PATH = (
    ROOT / "validation" / "known_assignment_confirmatory_execution_registration_v2_20260811.json"
)
DEFAULT_STUDY_DIR = ROOT / "validation" / "kac200_20260811"
FACETS_EXE_DEFAULT = Path(r"C:\Facets\Facets.exe")
PRIMARY_ID = "KA1_FREE_MML_EXPECTED_WITHIN_DATASET_FIXED_RATER_RMSE_STRESS_CONTRAST"
SECONDARY_IDS = (
    "KA2_FIXED_MML_EXPECTED_WITHIN_DATASET_FIXED_RATER_RMSE_STRESS_CONTRAST",
    "KA3_FREE_MML_SYMMETRIC_STRESS_POPULATION_SD_SHIFT",
    "KA4_FREE_MML_DIRECTION_ALIGNED_RATER_ERROR_SLOPE",
)
MML_ATTEMPT_TYPES = (
    "PYTHON_MML_FIXED_SD08_Q31_PCM",
    "PYTHON_MML_FREE_SD_Q31_PCM",
)
PAIR_ATTEMPT_TYPE = "RESILIENT_FACETS_PYTHON_JMLE_PCM"
CMLE_ATTEMPT_TYPE = "PYTHON_EXACT_CMLE_PCM"
SCIENTIFIC_MODES = (MML_FIXED_MODE, MML_FREE_MODE)
REQUIRED_TRIPLETS = 200
REQUIRED_DATASETS = 600
CALIBRATION_TRIPLETS = (5, 22, 41, 57, 72, 99, 116, 128, 143, 161, 176, 200)
CALIBRATION_DATASETS = 36
REQUIRED_ATTEMPTS = 1272
INPUT_FILES = engine.INPUT_FILES


def _json_dump(path: Path, value: Any) -> None:
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _atomic_json(path: Path, value: Any) -> None:
    temporary = path.with_name(f"{path.name}.tmp.{os.getpid()}")
    _json_dump(temporary, value)
    os.replace(temporary, path)


def _portable_name(path: Path) -> str:
    resolved = path.resolve()
    try:
        return resolved.relative_to(ROOT).as_posix()
    except ValueError:
        return str(resolved)


def _canonical_subset(values: Iterable[int]) -> str:
    return json.dumps([int(value) for value in values], separators=(",", ":"))


def validate_plan(plan: dict[str, Any] | None = None) -> dict[str, Any]:
    """Fail closed on scientific drift in the human-readable frozen plan."""

    value = plan or json.loads(PLAN_PATH.read_text(encoding="utf-8"))
    evidence = value["evidence_status"]
    if value.get("status") != (
        "scientific_design_frozen_before_any_confirmatory_person_assignment_or_response_generation"
    ):
        raise ValueError("Confirmatory plan is not frozen before data generation")
    exact_evidence = {
        "independent_person_vector_triplets": REQUIRED_TRIPLETS,
        "datasets": REQUIRED_DATASETS,
        "gammas_per_triplet": 3,
        "fixed_sample_size": True,
        "optional_stopping": False,
        "interim_scientific_looks": 0,
        "replacement_triplets": False,
        "screening_observations_pooled": False,
        "preflight_observations_pooled": False,
        "extension_after_precision_failure": False,
    }
    for key, expected in exact_evidence.items():
        if evidence.get(key) != expected:
            raise ValueError(f"Confirmatory evidence contract drift: {key}")

    primary = value["primary_endpoint"]
    if primary.get("id") != PRIMARY_ID or float(primary.get("alpha", np.nan)) != 0.05:
        raise ValueError("Primary endpoint or alpha changed")
    if "not E(estimate-truth)" not in str(primary.get("not_bias", "")):
        raise ValueError("Primary endpoint must retain its non-bias boundary")
    secondary = tuple(row["id"] for row in value["secondary_family"]["endpoints"])
    if secondary != SECONDARY_IDS:
        raise ValueError("Secondary endpoint family changed")

    lanes = value["execution_lanes"]
    subset = tuple(int(item) for item in lanes["facets_external_calibration"]["triplet_subset"])
    if subset != CALIBRATION_TRIPLETS:
        raise ValueError("FACETS calibration subset changed")
    if tuple(lanes["independent_r_free_mml"]["triplet_subset"]) != subset:
        raise ValueError("R MML subset must equal the FACETS subset")
    if tuple(lanes["exact_cmle_descriptive"]["triplet_subset"]) != subset:
        raise ValueError("CMLE subset must equal the frozen calibration subset")
    canonical = _canonical_subset(subset)
    subset_hash = hashlib.sha256(canonical.encode("utf-8")).hexdigest()
    if canonical != lanes["facets_external_calibration"]["subset_json_canonical"]:
        raise ValueError("Calibration subset canonical serialization changed")
    if subset_hash != lanes["facets_external_calibration"]["subset_sha256"]:
        raise ValueError("Calibration subset hash changed")
    if int(lanes["registered_python_attempt_units_total"]) != REQUIRED_ATTEMPTS:
        raise ValueError("Registered attempt denominator changed")

    sample = value["sample_size_contract"]
    planning_sd = float(sample["planning_sd"])
    df = int(sample["planning_sd_df"])
    sigma_upper = planning_sd * np.sqrt(df / chi2.ppf(0.025, df))
    half_width = float(student_t.ppf(0.975, REQUIRED_TRIPLETS - 1)) * sigma_upper / np.sqrt(
        REQUIRED_TRIPLETS
    )
    if not np.isclose(sigma_upper, sample["planning_sd_two_sided_95_percent_upper"], atol=1e-15):
        raise ValueError("Registered planning SD upper bound is not reproducible")
    if not np.isclose(half_width, sample["expected_two_sided_95_percent_half_width"], atol=1e-15):
        raise ValueError("Registered precision calculation is not reproducible")
    if float(sample["target_half_width"]) != 0.015:
        raise ValueError("Primary precision target changed")

    dgm = value["data_generating_process"]
    assignment = value["assignment"]
    person_seeds = set(range(int(dgm["person_seed_range"][0]), int(dgm["person_seed_range"][1]) + 1))
    response_seeds = set(
        range(int(dgm["response_uniform_seed_range"][0]), int(dgm["response_uniform_seed_range"][1]) + 1)
    )
    assignment_seeds = set(
        range(int(assignment["assignment_seed_range"][0]), int(assignment["assignment_seed_range"][1]) + 1)
    )
    if len(person_seeds) != 200 or len(response_seeds) != 200 or len(assignment_seeds) != 600:
        raise ValueError("Seed ranges do not match frozen denominators")
    if person_seeds & response_seeds or person_seeds & assignment_seeds or response_seeds & assignment_seeds:
        raise ValueError("Confirmatory seed streams overlap")

    for path_key, hash_key in (
        ("preflight_plan", "preflight_plan_sha256"),
        ("preflight_visible_join_assessment", "preflight_visible_join_assessment_sha256"),
        ("preflight_r_aggregation_assessment", "preflight_r_aggregation_assessment_sha256"),
    ):
        path = ROOT / value["prior_evidence"][path_key]
        if sha256_file(path).lower() != str(value["prior_evidence"][hash_key]).lower():
            raise ValueError(f"Prior evidence hash mismatch: {path_key}")
    return value


def runtime_plan(plan: dict[str, Any]) -> dict[str, Any]:
    """Expand frozen seed formulae into the shape expected by the preflight engine."""

    validate_plan(plan)
    value = copy.deepcopy(plan)
    vectors = list(range(1, REQUIRED_TRIPLETS + 1))
    value["data_generating_process"]["persons"] = int(
        value["data_generating_process"]["persons_per_triplet"]
    )
    value["data_generating_process"]["person_vector_ids"] = vectors
    value["data_generating_process"]["person_vector_seeds"] = [26100000 + item for item in vectors]
    value["data_generating_process"]["response_uniform_seeds"] = [
        26101000 + item for item in vectors
    ]
    gamma_values = tuple(float(item) for item in value["assignment"]["gammas"])
    value["assignment"]["assignment_seeds_by_person_vector"] = {
        str(vector): {
            str(gamma): 26102000 + 3 * (vector - 1) + ordinal
            for ordinal, gamma in enumerate(gamma_values, start=1)
        }
        for vector in vectors
    }
    return value


def _verify_prior_evidence(plan: dict[str, Any]) -> None:
    validate_plan(plan)
    visible = json.loads(
        (ROOT / plan["prior_evidence"]["preflight_visible_join_assessment"]).read_text(
            encoding="utf-8"
        )
    )
    if not bool(visible.get("qualification_pass", False)):
        raise ValueError("Known-assignment preflight did not qualify")
    if not bool(visible.get("eligible_to_freeze_separate_confirmation", False)):
        raise ValueError("Known-assignment preflight did not authorize confirmation")


def build_dependency_manifest(*, include_registration: bool = True) -> dict[str, str]:
    manifest = dependency_manifest()
    additions = [
        Path(__file__).resolve(),
        PLAN_PATH,
        ROOT / "validation" / "known_assignment_multivector_preflight.py",
        ROOT / "validation" / "known_assignment_mml_crossfit.py",
        ROOT / "validation" / "known_assignment_mml_crossfit.R",
        ROOT / "mfrm_app" / "assignment_mechanism.py",
    ]
    if include_registration:
        additions.append(REGISTRATION_PATH)
    for path in additions:
        if not path.is_file():
            raise FileNotFoundError(f"Confirmatory dependency missing: {path}")
        manifest[_portable_name(path)] = sha256_file(path)
    return dict(sorted(manifest.items()))


def validate_registration() -> dict[str, Any]:
    value = json.loads(REGISTRATION_PATH.read_text(encoding="utf-8"))
    expected = {
        "plan_sha256": sha256_file(PLAN_PATH),
        "runner_sha256": sha256_file(Path(__file__).resolve()),
        "preflight_engine_sha256": sha256_file(
            ROOT / "validation" / "known_assignment_multivector_preflight.py"
        ),
        "assignment_mechanism_sha256": sha256_file(
            ROOT / "mfrm_app" / "assignment_mechanism.py"
        ),
        "facets_resilient_pair_sha256": sha256_file(
            ROOT / "validation" / "facets_resilient_pair.py"
        ),
        "mml_crossfit_python_sha256": sha256_file(
            ROOT / "validation" / "known_assignment_mml_crossfit.py"
        ),
        "mml_crossfit_r_sha256": sha256_file(
            ROOT / "validation" / "known_assignment_mml_crossfit.R"
        ),
    }
    for key, digest in expected.items():
        if str(value.get(key, "")).lower() != digest.lower():
            raise ValueError(f"Execution registration mismatch: {key}")
    if not bool(value.get("tests_passed_before_confirmatory_generation", False)):
        raise ValueError("Registration does not assert passing pre-generation tests")
    if int(value.get("person_vectors_generated_before_registration", -1)) != 0:
        raise ValueError("Registration does not establish zero prior Person generation")
    if value.get("r_mml_tolerance") != mml_crossfit.default_tolerance():
        raise ValueError("Registered R MML tolerances changed")
    qualification = ROOT / str(value.get("r_mml_preflight_qualification", ""))
    if not qualification.is_file() or sha256_file(qualification) != value.get(
        "r_mml_preflight_qualification_sha256"
    ):
        raise ValueError("R MML preflight qualification is missing or changed")
    qualification_value = json.loads(qualification.read_text(encoding="utf-8"))
    if not bool(qualification_value.get("pass", False)):
        raise ValueError("R MML preflight qualification is not PASS")
    rscript = Path(str(value.get("rscript_executable", "")))
    if not rscript.is_file() or sha256_file(rscript) != value.get("rscript_executable_sha256"):
        raise ValueError("Registered Rscript executable is missing or changed")
    return value


@contextmanager
def _engine_context(study_dir: Path) -> Iterator[None]:
    """Temporarily bind the immutable preflight engine to this confirmation."""

    replacements = {
        "SCHEMA_VERSION": SCHEMA_VERSION,
        "PLAN_PATH": PLAN_PATH,
        "REGISTRATION_PATH": REGISTRATION_PATH,
        "STUDY_DIR": study_dir.resolve(),
        "_verify_prior_evidence": _verify_prior_evidence,
        "build_dependency_manifest": build_dependency_manifest,
        "validate_registration": validate_registration,
        "validate_study_identity": lambda: validate_study_identity(study_dir),
    }
    originals = {name: getattr(engine, name) for name in replacements}
    try:
        for name, replacement in replacements.items():
            setattr(engine, name, replacement)
        yield
    finally:
        for name, original in originals.items():
            setattr(engine, name, original)


def filter_attempt_manifest(attempts: pd.DataFrame) -> pd.DataFrame:
    attempt_type = attempts["AttemptType"].astype(str)
    vector = pd.to_numeric(attempts["PersonVector"], errors="raise").astype(int)
    keep = attempt_type.isin(MML_ATTEMPT_TYPES) | (
        vector.isin(CALIBRATION_TRIPLETS)
        & attempt_type.isin((PAIR_ATTEMPT_TYPE, CMLE_ATTEMPT_TYPE))
    )
    output = attempts.loc[keep].copy().reset_index(drop=True)
    counts = output["AttemptType"].value_counts().to_dict()
    expected = {
        MML_ATTEMPT_TYPES[0]: 600,
        MML_ATTEMPT_TYPES[1]: 600,
        PAIR_ATTEMPT_TYPE: 36,
        CMLE_ATTEMPT_TYPE: 36,
    }
    if len(output) != REQUIRED_ATTEMPTS or counts != expected:
        raise ValueError(f"Attempt filtering changed the registered denominator: {counts}")
    if output["AttemptId"].duplicated().any():
        raise ValueError("Attempt identifiers are not unique")
    return output


def _worst_case_facets_paths(study_dir: Path) -> dict[str, int]:
    base = (
        study_dir.resolve()
        / "work"
        / "02391"
        / "facets_attempts"
        / "pair_try_03"
        / "facets_runs"
        / "gamma_neg_0p8__vector-200__pcm"
    )
    return {
        name: len(str(base / name))
        for name in ("analysis.txt", "report_u6.txt", "scores_u6.txt")
    }


def prepare_study(study_dir: Path, *, facets_exe: Path) -> dict[str, Any]:
    study_dir = study_dir.resolve()
    if study_dir.exists():
        raise FileExistsError(f"Refusing to overwrite confirmatory study: {study_dir}")
    plan = validate_plan()
    validate_registration()
    facets_exe = facets_exe.resolve()
    if not facets_exe.is_file():
        raise FileNotFoundError(f"FACETS executable missing: {facets_exe}")
    dependencies = build_dependency_manifest()
    dependency_digest = dependency_manifest_digest(dependencies)
    with _engine_context(study_dir):
        bundle = engine.generate_bundle(runtime_plan(plan), dependency_digest=dependency_digest)
    bundle["attempt_manifest.csv"] = filter_attempt_manifest(bundle["attempt_manifest.csv"])

    study_dir.mkdir(parents=True, exist_ok=False)
    input_dir = study_dir / "retained_input"
    input_dir.mkdir()
    for filename, frame in bundle.items():
        frame.to_csv(input_dir / filename, index=False, lineterminator="\n")
    input_hashes = {filename: sha256_file(input_dir / filename) for filename in INPUT_FILES}
    manifest = bundle["manifest.csv"]
    attempts = bundle["attempt_manifest.csv"]
    audit = bundle["assignment_audit.csv"]
    diagnostics = bundle["assignment_diagnostics.csv"]
    path_lengths = _worst_case_facets_paths(study_dir)
    checks = {
        "triplets_200": manifest["PersonVector"].nunique() == REQUIRED_TRIPLETS,
        "datasets_600": len(manifest) == REQUIRED_DATASETS,
        "attempts_1272": len(attempts) == REQUIRED_ATTEMPTS,
        "ratings_576000": len(bundle["generated_ratings.csv"]) == 576_000,
        "all_assignment_and_design_audits": bool(audit["Passed"].all()),
        "all_category_support": bool(manifest["CompleteCriterionCategorySupport"].all()),
        "all_nullity_zero": bool(manifest["ExpectedStructuralNullity"].eq(0).all()),
        "all_connected": bool(manifest["PersonRaterComponents"].eq(1).all()),
        "shared_complete_response_within_triplet": bool(
            manifest.groupby("PersonVector")["CompleteResponseSHA256"].nunique().eq(1).all()
        ),
        "unique_person_vectors_200": manifest["PersonVectorSHA256"].nunique() == 200,
        "three_gamma_conditions_per_triplet": bool(
            manifest.groupby("PersonVector")["Gamma"].nunique().eq(3).all()
        ),
        "all_dense_draws_available": len(diagnostics) == 600
        and bool(diagnostics["Available"].all()),
        "dense_conditional_residual": bool(
            diagnostics["MaximumConditionalLogResidual"].le(1e-10).all()
        ),
        "dense_statistic_residual": bool(
            diagnostics["MaximumStatisticUpdateResidual"].le(1e-10).all()
        ),
        "facets_path_budget_220": max(path_lengths.values()) <= 220,
        "no_attempts_or_endpoints_before_prepare": True,
    }
    identity = {
        "schema_version": SCHEMA_VERSION,
        "phase": "registered_confirmatory_input_prepared",
        "plan_file": _portable_name(PLAN_PATH),
        "plan_sha256": sha256_file(PLAN_PATH),
        "registration_file": _portable_name(REGISTRATION_PATH),
        "registration_sha256": sha256_file(REGISTRATION_PATH),
        "runner_sha256": sha256_file(Path(__file__).resolve()),
        "dependency_manifest_sha256": dependency_digest,
        "facets_executable": str(facets_exe),
        "facets_executable_sha256": sha256_file(facets_exe),
        "rscript_executable": str(Path(validate_registration()["rscript_executable"]).resolve()),
        "rscript_executable_sha256": validate_registration()["rscript_executable_sha256"],
        "r_mml_tolerance": validate_registration()["r_mml_tolerance"],
        "retained_input_sha256": input_hashes,
        "triplets": REQUIRED_TRIPLETS,
        "datasets": len(manifest),
        "attempts": len(attempts),
        "calibration_triplets": list(CALIBRATION_TRIPLETS),
        "worst_case_facets_path_lengths": path_lengths,
        "python": sys.version,
        "platform": platform.platform(),
        "checks": {key: bool(value) for key, value in checks.items()},
        "all_checks_pass": bool(all(checks.values())),
        "screening_or_preflight_observations_pooled": False,
        "endpoint_summary_computed": False,
        "claim_limit": plan["claim_limit"],
    }
    _json_dump(study_dir / "dependency_manifest.json", dependencies)
    _json_dump(study_dir / "study_identity.json", identity)
    if not identity["all_checks_pass"]:
        raise RuntimeError(f"Confirmatory input audit failed: {checks}")
    return identity


def validate_study_identity(study_dir: Path) -> dict[str, Any]:
    study_dir = study_dir.resolve()
    identity = json.loads((study_dir / "study_identity.json").read_text(encoding="utf-8"))
    validate_plan()
    validate_registration()
    expected = {
        "plan_sha256": sha256_file(PLAN_PATH),
        "registration_sha256": sha256_file(REGISTRATION_PATH),
        "runner_sha256": sha256_file(Path(__file__).resolve()),
    }
    for key, digest in expected.items():
        if str(identity.get(key, "")).lower() != digest.lower():
            raise ValueError(f"Study identity mismatch: {key}")
    dependencies = build_dependency_manifest()
    if identity["dependency_manifest_sha256"] != dependency_manifest_digest(dependencies):
        raise ValueError("Execution dependency manifest changed")
    retained = json.loads((study_dir / "dependency_manifest.json").read_text(encoding="utf-8"))
    if retained != dependencies:
        raise ValueError("Retained dependency manifest differs from current dependencies")
    if sha256_file(Path(identity["facets_executable"])) != identity["facets_executable_sha256"]:
        raise ValueError("FACETS executable changed")
    if sha256_file(Path(identity["rscript_executable"])) != identity["rscript_executable_sha256"]:
        raise ValueError("Rscript executable changed")
    for filename, digest in identity["retained_input_sha256"].items():
        if sha256_file(study_dir / "retained_input" / filename) != digest:
            raise ValueError(f"Retained input changed: {filename}")
    if not bool(identity.get("all_checks_pass", False)):
        raise ValueError("Prepared confirmatory input audit was not PASS")
    return identity


def run_shard(
    study_dir: Path,
    *,
    shard_index: int,
    shard_count: int,
    resume: bool,
    timeout_seconds: float,
) -> dict[str, int]:
    with _engine_context(study_dir):
        return engine.run_shard(
            shard_index=shard_index,
            shard_count=shard_count,
            resume=resume,
            timeout_seconds=timeout_seconds,
        )


def run_r_crossfit(study_dir: Path, *, maxit: int = 200) -> dict[str, Any]:
    """Run the frozen independent R Q31/Q61 cross-fit before endpoint analysis."""

    study_dir = study_dir.resolve()
    identity = validate_study_identity(study_dir)
    output_dir = study_dir / "r_mml"
    if output_dir.exists():
        raise FileExistsError(f"Refusing to overwrite R MML evidence: {output_dir}")
    attempts = pd.read_csv(study_dir / "retained_input" / "attempt_manifest.csv")
    free = attempts.loc[attempts["AttemptType"].astype(str).eq(MML_ATTEMPT_TYPES[1])]
    selected = free.loc[
        pd.to_numeric(free["PersonVector"], errors="raise")
        .astype(int)
        .isin(CALIBRATION_TRIPLETS)
    ]
    if len(selected) != 36:
        raise ValueError("Frozen R cross-fit subset is not 36 datasets")
    with _engine_context(study_dir):
        for _, attempt in selected.iterrows():
            marker = engine._completion_path(attempt)  # pylint: disable=protected-access
            if not marker.is_file():
                raise FileNotFoundError(f"R cross-fit requires completed MML marker: {attempt['AttemptId']}")
            completion = json.loads(marker.read_text(encoding="utf-8"))
            engine._validate_completion(completion, attempt, identity)  # pylint: disable=protected-access
            if not bool(completion.get("statistical_evidence_ready", False)):
                raise RuntimeError(f"R cross-fit source is not evidence-ready: {attempt['AttemptId']}")
    output_dir.mkdir()
    command = [
        identity["rscript_executable"],
        str(ROOT / "validation" / "known_assignment_mml_crossfit.R"),
        "--study",
        str(study_dir),
        "--output",
        str(output_dir),
        "--vectors",
        ",".join(map(str, CALIBRATION_TRIPLETS)),
        "--maxit",
        str(int(maxit)),
    ]
    completed = subprocess.run(
        command,
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
        timeout=1800,
    )
    _json_dump(
        output_dir / "r_process.json",
        {
            "command": command,
            "returncode": completed.returncode,
            "stdout": completed.stdout,
            "stderr": completed.stderr,
        },
    )
    if completed.returncode != 0:
        raise RuntimeError(f"Independent R MML cross-fit failed: {completed.stderr}")
    assessment = mml_crossfit.assess_crossfit(
        study_dir=study_dir,
        output_dir=output_dir,
        expected_datasets=36,
        tolerance=dict(identity["r_mml_tolerance"]),
    )
    if not bool(assessment["pass"]):
        raise RuntimeError("Independent R MML cross-fit did not pass frozen tolerances")
    return assessment


def preexecution_audit(study_dir: Path) -> dict[str, Any]:
    identity = validate_study_identity(study_dir)
    input_dir = study_dir.resolve() / "retained_input"
    manifest = pd.read_csv(input_dir / "manifest.csv")
    attempts = pd.read_csv(input_dir / "attempt_manifest.csv")
    ratings = pd.read_csv(input_dir / "generated_ratings.csv")
    audit = pd.read_csv(input_dir / "assignment_audit.csv")
    counts = attempts["AttemptType"].value_counts().to_dict()
    checks = {
        "identity_valid": True,
        "triplets_200": manifest["PersonVector"].nunique() == 200,
        "datasets_600": len(manifest) == 600,
        "attempts_1272": len(attempts) == 1272,
        "mml_free_600": counts.get(MML_ATTEMPT_TYPES[1], 0) == 600,
        "mml_fixed_600": counts.get(MML_ATTEMPT_TYPES[0], 0) == 600,
        "facets_pairs_36": counts.get(PAIR_ATTEMPT_TYPE, 0) == 36,
        "cmle_36": counts.get(CMLE_ATTEMPT_TYPE, 0) == 36,
        "ratings_576000": len(ratings) == 576_000,
        "all_input_audits_pass": bool(audit["Passed"].all()),
        "all_facets_paths_at_most_220": max(
            identity["worst_case_facets_path_lengths"].values()
        )
        <= 220,
        "no_completion_markers": not (study_dir / "attempts").exists(),
        "no_work_artifacts": not (study_dir / "work").exists(),
        "no_aggregate": not (study_dir / "aggregate").exists(),
    }
    result = {
        "schema_version": f"{SCHEMA_VERSION}_preexecution_audit_v1",
        "checks": {key: bool(value) for key, value in checks.items()},
        "all_checks_pass": bool(all(checks.values())),
        "endpoint_summary_computed": False,
        "study_identity_sha256": sha256_file(study_dir / "study_identity.json"),
    }
    _json_dump(study_dir / "preexecution_audit.json", result)
    return result


def holm_adjust(p_values: Iterable[float]) -> np.ndarray:
    values = np.asarray(list(p_values), dtype=float)
    if not len(values) or not np.isfinite(values).all():
        raise ValueError("Holm adjustment requires finite p values")
    order = np.argsort(values, kind="mergesort")
    sorted_adjusted = np.empty(len(values), dtype=float)
    running = 0.0
    for rank, index in enumerate(order):
        running = max(running, (len(values) - rank) * values[index])
        sorted_adjusted[rank] = min(1.0, running)
    adjusted = np.empty(len(values), dtype=float)
    adjusted[order] = sorted_adjusted
    return adjusted


def build_endpoint_contrasts(
    *, rater_errors: pd.DataFrame, runs: pd.DataFrame
) -> pd.DataFrame:
    errors = rater_errors.copy()
    errors["ErrorAligned"] = pd.to_numeric(errors["ErrorAligned"], errors="coerce")
    loss = (
        errors.groupby(["PersonVector", "Gamma", "EstimatorMode"], as_index=False)[
            "ErrorAligned"
        ]
        .agg(RaterRows="count", RaterRMSE=lambda x: float(np.sqrt(np.mean(np.square(x)))))
    )
    loss_wide = loss.pivot(
        index=["PersonVector", "EstimatorMode"], columns="Gamma", values="RaterRMSE"
    ).reset_index()
    for gamma in (-0.8, 0.0, 0.8):
        if gamma not in loss_wide:
            loss_wide[gamma] = np.nan
    loss_wide["Contrast"] = 0.5 * (loss_wide[-0.8] + loss_wide[0.8]) - loss_wide[0.0]

    rows: list[pd.DataFrame] = []
    for endpoint_id, mode in ((PRIMARY_ID, MML_FREE_MODE), (SECONDARY_IDS[0], MML_FIXED_MODE)):
        selected = loss_wide.loc[
            loss_wide["EstimatorMode"].astype(str).eq(mode), ["PersonVector", "Contrast"]
        ].copy()
        selected.insert(0, "EndpointId", endpoint_id)
        rows.append(selected)

    free = runs.loc[runs["EstimatorMode"].astype(str).eq(MML_FREE_MODE)].copy()
    free["EstimatedPopulationSD"] = pd.to_numeric(free["EstimatedPopulationSD"], errors="coerce")
    sd_wide = free.pivot(index="PersonVector", columns="Gamma", values="EstimatedPopulationSD").reset_index()
    for gamma in (-0.8, 0.0, 0.8):
        if gamma not in sd_wide:
            sd_wide[gamma] = np.nan
    sd_endpoint = sd_wide[["PersonVector"]].copy()
    sd_endpoint["Contrast"] = 0.5 * (sd_wide[-0.8] + sd_wide[0.8]) - sd_wide[0.0]
    sd_endpoint.insert(0, "EndpointId", SECONDARY_IDS[1])
    rows.append(sd_endpoint)

    truth = pd.Series({str(key): float(value) for key, value in engine.RATER_TRUTH.items()})
    free_errors = errors.loc[errors["EstimatorMode"].astype(str).eq(MML_FREE_MODE)].copy()
    free_errors["TrueSeverity"] = free_errors["Level"].astype(str).map(truth)
    denominator = float(np.square(truth.to_numpy(dtype=float)).sum())
    slopes = (
        free_errors.assign(Product=lambda frame: frame["TrueSeverity"] * frame["ErrorAligned"])
        .groupby(["PersonVector", "Gamma"], as_index=False)["Product"]
        .sum()
    )
    slopes["Slope"] = slopes.pop("Product") / denominator
    slope_wide = slopes.pivot(index="PersonVector", columns="Gamma", values="Slope").reset_index()
    for gamma in (-0.8, 0.8):
        if gamma not in slope_wide:
            slope_wide[gamma] = np.nan
    slope_endpoint = slope_wide[["PersonVector"]].copy()
    slope_endpoint["Contrast"] = 0.5 * (slope_wide[-0.8] - slope_wide[0.8])
    slope_endpoint.insert(0, "EndpointId", SECONDARY_IDS[2])
    rows.append(slope_endpoint)

    output = pd.concat(rows, ignore_index=True)
    observed = tuple(output["EndpointId"].drop_duplicates())
    if observed != (PRIMARY_ID, *SECONDARY_IDS):
        raise RuntimeError("Endpoint construction order changed")
    return output


def _endpoint_summary(endpoint_id: str, values: pd.Series, *, alternative: str) -> dict[str, Any]:
    numeric = pd.to_numeric(values, errors="coerce").dropna().to_numpy(dtype=float)
    n = len(numeric)
    mean = float(np.mean(numeric)) if n else np.nan
    sd = float(np.std(numeric, ddof=1)) if n > 1 else np.nan
    se = sd / np.sqrt(n) if n > 1 else np.nan
    if np.isfinite(se) and se > 0:
        statistic = mean / se
    elif np.isfinite(se) and se == 0 and np.isfinite(mean) and mean != 0:
        statistic = np.copysign(np.inf, mean)
    else:
        statistic = np.nan
    degrees = n - 1
    p_value = (
        float(student_t.sf(statistic, degrees))
        if alternative == "greater" and not np.isnan(statistic)
        else float(student_t.cdf(statistic, degrees))
        if alternative == "less" and not np.isnan(statistic)
        else np.nan
    )
    half = float(student_t.ppf(0.975, degrees)) * se if np.isfinite(se) else np.nan
    direction = bool(
        np.isfinite(mean)
        and ((alternative == "greater" and mean > 0) or (alternative == "less" and mean < 0))
    )
    return {
        "EndpointId": endpoint_id,
        "Alternative": alternative,
        "FiniteTriplets": n,
        "RequiredTriplets": REQUIRED_TRIPLETS,
        "FullTripletGate": n == REQUIRED_TRIPLETS,
        "MeanContrast": mean,
        "MonteCarloSD": sd,
        "MonteCarloSE": se,
        "TStatistic": statistic,
        "DegreesOfFreedom": degrees,
        "RawOneSidedP": p_value,
        "Lower95": mean - half if np.isfinite(half) else np.nan,
        "Upper95": mean + half if np.isfinite(half) else np.nan,
        "TwoSided95HalfWidth": half,
        "DirectionPass": direction,
    }


def summarize_endpoints(contrasts: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    primary = pd.DataFrame(
        [
            _endpoint_summary(
                PRIMARY_ID,
                contrasts.loc[contrasts["EndpointId"].eq(PRIMARY_ID), "Contrast"],
                alternative="greater",
            )
        ]
    )
    primary["AdjustedP"] = primary["RawOneSidedP"]
    primary["DirectionConfirmed"] = (
        primary["FullTripletGate"].astype(bool)
        & primary["DirectionPass"].astype(bool)
        & primary["RawOneSidedP"].le(0.05)
    )
    primary["PrecisionTarget"] = 0.015
    primary["PrecisionQualified"] = (
        primary["FullTripletGate"].astype(bool)
        & primary["TwoSided95HalfWidth"].le(0.015)
    )

    alternatives = {
        SECONDARY_IDS[0]: "greater",
        SECONDARY_IDS[1]: "less",
        SECONDARY_IDS[2]: "greater",
    }
    secondary = pd.DataFrame(
        [
            _endpoint_summary(
                endpoint_id,
                contrasts.loc[contrasts["EndpointId"].eq(endpoint_id), "Contrast"],
                alternative=alternatives[endpoint_id],
            )
            for endpoint_id in SECONDARY_IDS
        ]
    )
    if secondary["RawOneSidedP"].notna().all():
        secondary["HolmAdjustedP"] = holm_adjust(secondary["RawOneSidedP"])
    else:
        secondary["HolmAdjustedP"] = np.nan
    primary_gate = bool(primary.iloc[0]["DirectionConfirmed"])
    secondary["PrimaryGatePass"] = primary_gate
    secondary["DirectionConfirmed"] = (
        secondary["PrimaryGatePass"].astype(bool)
        & secondary["FullTripletGate"].astype(bool)
        & secondary["DirectionPass"].astype(bool)
        & secondary["HolmAdjustedP"].le(0.05)
    )
    secondary["IntervalMultiplicityStatus"] = "unadjusted_descriptive_95_percent"
    return primary, secondary


def bias_mse_decomposition(rater_errors: pd.DataFrame) -> pd.DataFrame:
    selected = rater_errors.loc[
        rater_errors["EstimatorMode"].astype(str).isin(SCIENTIFIC_MODES)
    ].copy()
    selected["ErrorAligned"] = pd.to_numeric(selected["ErrorAligned"], errors="coerce")
    rows: list[dict[str, Any]] = []
    for keys, group in selected.groupby(["EstimatorMode", "Gamma", "Level"], sort=True):
        finite = group["ErrorAligned"].dropna().to_numpy(dtype=float)
        n = len(finite)
        bias = float(np.mean(finite)) if n else np.nan
        variance = float(np.var(finite, ddof=1)) if n > 1 else np.nan
        mse = float(np.mean(np.square(finite))) if n else np.nan
        identity = (
            mse - bias**2 - ((n - 1) / n) * variance
            if n > 1 and np.isfinite(bias) and np.isfinite(variance) and np.isfinite(mse)
            else np.nan
        )
        rows.append(
            {
                "EstimatorMode": str(keys[0]),
                "Gamma": float(keys[1]),
                "Rater": str(keys[2]),
                "FiniteTriplets": n,
                "RequiredTriplets": REQUIRED_TRIPLETS,
                "FullTripletGate": n == REQUIRED_TRIPLETS,
                "MonteCarloBias": bias,
                "ErrorVarianceSample": variance,
                "MSE": mse,
                "RMSEAcrossTriplets": float(np.sqrt(mse)) if np.isfinite(mse) else np.nan,
                "MSEIdentityResidual": identity,
                "InferenceStatus": "descriptive_registered_decomposition",
            }
        )
    return pd.DataFrame(rows)


def _read_completed_artifacts(
    study_dir: Path, identity: dict[str, Any]
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, str]]:
    attempts = pd.read_csv(study_dir / "retained_input" / "attempt_manifest.csv")
    parts: dict[str, list[pd.DataFrame]] = {
        "run_ledger.csv": [],
        "recovery.csv": [],
        "thresholds.csv": [],
        "constraints.csv": [],
    }
    outcomes: list[dict[str, Any]] = []
    marker_hashes: dict[str, str] = {}
    with _engine_context(study_dir):
        for _, attempt in attempts.iterrows():
            marker = engine._completion_path(attempt)  # pylint: disable=protected-access
            if not marker.is_file():
                raise FileNotFoundError(f"Missing completion marker: {attempt['AttemptId']}")
            completion = json.loads(marker.read_text(encoding="utf-8"))
            engine._validate_completion(completion, attempt, identity)  # pylint: disable=protected-access
            marker_hashes[str(attempt["AttemptId"])] = sha256_file(marker)
            artifact_root = study_dir / str(completion["artifact_root"])
            for filename in parts:
                frame = engine._read_nonempty(artifact_root / filename)  # pylint: disable=protected-access
                if len(frame):
                    frame.insert(0, "AttemptId", str(attempt["AttemptId"]))
                    parts[filename].append(frame)
            outcome = {
                "AttemptId": str(attempt["AttemptId"]),
                "RunId": str(attempt["RunId"]),
                "AttemptType": str(attempt["AttemptType"]),
                "PersonVector": int(attempt["PersonVector"]),
                "Gamma": float(attempt["Gamma"]),
                "ExecutionCompleted": bool(completion["execution_completed"]),
                "StatisticalEvidenceReady": bool(completion["statistical_evidence_ready"]),
                "FACETSCalibrationReady": completion.get("facets_calibration_ready"),
                "FailureReason": completion.get("failure_reason", ""),
            }
            metrics_path = artifact_root / "resilient_pair_metrics.json"
            if metrics_path.is_file():
                outcome.update(json.loads(metrics_path.read_text(encoding="utf-8")))
            outcomes.append(outcome)
    combined = {
        filename: pd.concat(frames, ignore_index=True, sort=False) if frames else pd.DataFrame()
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


def aggregate_study(study_dir: Path) -> dict[str, Any]:
    study_dir = study_dir.resolve()
    aggregate_dir = study_dir / "aggregate"
    if aggregate_dir.exists():
        raise FileExistsError(f"Refusing to overwrite aggregate: {aggregate_dir}")
    identity = validate_study_identity(study_dir)
    manifest = pd.read_csv(study_dir / "retained_input" / "manifest.csv")
    runs, recovery, thresholds, constraints, outcomes, marker_hashes = _read_completed_artifacts(
        study_dir, identity
    )
    runs = engine._attach_context(runs, manifest)  # pylint: disable=protected-access
    recovery = engine._attach_context(recovery, manifest)  # pylint: disable=protected-access
    thresholds = engine._attach_context(thresholds, manifest)  # pylint: disable=protected-access
    constraints = (
        engine._attach_context(constraints, manifest)  # pylint: disable=protected-access
        if len(constraints)
        else constraints
    )

    mml = outcomes[outcomes["AttemptType"].isin(MML_ATTEMPT_TYPES)].copy()
    free_outcomes = outcomes[outcomes["AttemptType"].eq(MML_ATTEMPT_TYPES[1])]
    fixed_outcomes = outcomes[outcomes["AttemptType"].eq(MML_ATTEMPT_TYPES[0])]
    jmle = outcomes[outcomes["AttemptType"].eq(PAIR_ATTEMPT_TYPE)].copy()
    cmle = outcomes[outcomes["AttemptType"].eq(CMLE_ATTEMPT_TYPE)].copy()
    mml_constraints = constraints.loc[
        constraints["EstimatorMode"].astype(str).isin(SCIENTIFIC_MODES)
    ]
    rater_errors = recovery.loc[
        recovery["EstimatorMode"].astype(str).isin(SCIENTIFIC_MODES)
        & recovery["IncludedInStudy"].fillna(False).astype(bool)
        & recovery["Facet"].astype(str).eq("Rater")
    ].copy()
    contrasts = build_endpoint_contrasts(rater_errors=rater_errors, runs=runs)
    primary, secondary = summarize_endpoints(contrasts)
    decomposition = bias_mse_decomposition(rater_errors)

    scientific_gates = {
        "completion_markers_1272": len(marker_hashes) == REQUIRED_ATTEMPTS,
        "mml_attempts_1200": len(mml) == 1200,
        "mml_all_execution_completed": bool(mml["ExecutionCompleted"].all()),
        "mml_all_statistical_evidence_ready": bool(mml["StatisticalEvidenceReady"].all()),
        "free_mml_600_ready": len(free_outcomes) == 600
        and bool(free_outcomes["StatisticalEvidenceReady"].all()),
        "fixed_mml_600_ready": len(fixed_outcomes) == 600
        and bool(fixed_outcomes["StatisticalEvidenceReady"].all()),
        "mml_constraints_1200_pass": len(mml_constraints) == 1200
        and bool(mml_constraints["ConstraintPass"].fillna(False).astype(bool).all()),
        "four_endpoints_200_triplets_each": len(contrasts) == 800
        and bool(contrasts.groupby("EndpointId").size().eq(200).all()),
        "primary_complete_triplets_200": bool(primary["FullTripletGate"].all()),
        "bias_decomposition_24_cells": len(decomposition) == 24
        and bool(decomposition["FullTripletGate"].all()),
    }
    r_assessment_path = study_dir / "r_mml" / "assessment.json"
    r_assessment = (
        json.loads(r_assessment_path.read_text(encoding="utf-8"))
        if r_assessment_path.is_file()
        else {"pass": False, "status": "missing"}
    )
    r_gate = bool(r_assessment.get("pass", False))
    scientific_endpoint_computed = bool(all(scientific_gates.values()))
    scientific_inference_validated = bool(scientific_endpoint_computed and r_gate)

    calibration_ready = jmle["FACETSCalibrationReady"].fillna(False).astype(bool)
    facets_gates = {
        "facets_pair_attempts_36": len(jmle) == 36,
        "facets_pair_execution_completed_36": len(jmle) == 36
        and bool(jmle["ExecutionCompleted"].all()),
        "facets_python_calibration_ready_36": len(jmle) == 36
        and bool(calibration_ready.all()),
    }
    facets_complement_qualification = bool(all(facets_gates.values()))
    cmle_gates = {
        "cmle_attempts_36": len(cmle) == 36,
        "cmle_ready_36": len(cmle) == 36 and bool(cmle["StatisticalEvidenceReady"].all()),
    }

    facets_runs = runs.loc[runs["EstimatorMode"].astype(str).eq(FACETS_MODE)].copy()
    numeric = lambda column: pd.to_numeric(facets_runs.get(column), errors="coerce")
    calibration = {
        "planned_pairs": 36,
        "attempted_pairs": int(len(jmle)),
        "ready_pairs": int(calibration_ready.sum()),
        "maximum_main_weighted_mae": float(numeric("MainWeightedMAE").max())
        if len(facets_runs)
        else None,
        "maximum_main_absolute_difference": float(numeric("MainMaxAbsDifference").max())
        if len(facets_runs)
        else None,
        "maximum_threshold_weighted_mae": float(numeric("ThresholdWeightedMAE").max())
        if len(facets_runs)
        else None,
        "maximum_threshold_absolute_difference": float(numeric("ThresholdMaxAbsDifference").max())
        if len(facets_runs)
        else None,
        "minimum_within_facet_spearman": float(numeric("MinimumWithinFacetSpearman").min())
        if len(facets_runs)
        else None,
        "ordinary_two_decimal_fit_used_as_raw_input": False,
    }

    aggregate_dir.mkdir()
    outputs = {
        "attempt_outcomes.csv": outcomes,
        "run_ledger.csv": runs,
        "recovery.csv": recovery,
        "thresholds.csv": thresholds,
        "constraints.csv": constraints,
        "registered_endpoint_contrasts.csv": contrasts,
        "primary_result.csv": primary,
        "secondary_results.csv": secondary,
        "registered_bias_mse_decomposition.csv": decomposition,
    }
    for filename, frame in outputs.items():
        frame.to_csv(aggregate_dir / filename, index=False, lineterminator="\n")
    marker_digest = hashlib.sha256(
        json.dumps(marker_hashes, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    assessment = {
        "schema_version": SCHEMA_VERSION,
        "phase": "confirmatory",
        "scientific_gates": {key: bool(value) for key, value in scientific_gates.items()},
        "independent_r_mml_gate": r_gate,
        "facets_gates": {key: bool(value) for key, value in facets_gates.items()},
        "cmle_descriptive_gates": {key: bool(value) for key, value in cmle_gates.items()},
        "ScientificEndpointComputed": scientific_endpoint_computed,
        "ScientificInferenceValidated": scientific_inference_validated,
        "FACETSComplementQualification": facets_complement_qualification,
        "FACETSValidatedComplementaryWorkbench": bool(
            scientific_inference_validated and facets_complement_qualification
        ),
        "primary_direction_confirmed": bool(primary.iloc[0]["DirectionConfirmed"]),
        "primary_precision_qualified": bool(primary.iloc[0]["PrecisionQualified"]),
        "facets_calibration": calibration,
        "completion_marker_set_sha256": marker_digest,
        "facets_failure_does_not_change_mml_denominator": True,
        "ordinary_facets_fit_used_as_raw_input": False,
        "estimator_ranking_constructed": False,
        "cross_basis_likelihood_comparison": "prohibited",
        "claim_limit": identity["claim_limit"],
    }
    _json_dump(aggregate_dir / "assessment.json", assessment)
    artifact_names = tuple(outputs) + ("assessment.json",)
    _atomic_json(
        aggregate_dir / "aggregate_identity.json",
        {
            "schema_version": f"{SCHEMA_VERSION}_aggregate_identity_v1",
            "study_identity_sha256": sha256_file(study_dir / "study_identity.json"),
            "completion_marker_set_sha256": marker_digest,
            "r_mml_assessment_sha256": sha256_file(r_assessment_path)
            if r_assessment_path.is_file()
            else None,
            "artifact_sha256": {
                filename: sha256_file(aggregate_dir / filename) for filename in artifact_names
            },
        },
    )
    return assessment


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    prepare = subparsers.add_parser("prepare")
    prepare.add_argument("--study-dir", type=Path, default=DEFAULT_STUDY_DIR)
    prepare.add_argument("--facets-exe", type=Path, default=FACETS_EXE_DEFAULT)
    audit = subparsers.add_parser("audit")
    audit.add_argument("--study-dir", type=Path, default=DEFAULT_STUDY_DIR)
    run = subparsers.add_parser("run")
    run.add_argument("--study-dir", type=Path, default=DEFAULT_STUDY_DIR)
    run.add_argument("--shard-index", type=int, required=True)
    run.add_argument("--shard-count", type=int, required=True)
    run.add_argument("--resume", action="store_true")
    run.add_argument("--timeout-seconds", type=float, default=120.0)
    aggregate = subparsers.add_parser("aggregate")
    aggregate.add_argument("--study-dir", type=Path, default=DEFAULT_STUDY_DIR)
    r_crossfit = subparsers.add_parser("r-crossfit")
    r_crossfit.add_argument("--study-dir", type=Path, default=DEFAULT_STUDY_DIR)
    r_crossfit.add_argument("--maxit", type=int, default=200)
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
    elif args.command == "r-crossfit":
        result = run_r_crossfit(args.study_dir, maxit=args.maxit)
    else:
        result = aggregate_study(args.study_dir)
    print(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True, default=str))


if __name__ == "__main__":
    main()
