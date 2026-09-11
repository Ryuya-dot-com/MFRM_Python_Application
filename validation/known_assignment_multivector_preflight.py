#!/usr/bin/env python3
"""Prepare, execute, and aggregate the frozen multi-vector preflight."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import sys
import time
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd

from mfrm_app.assignment_mechanism import sample_degree_conditioned_assignment_dense_dp
from validation.estimand_bridge_pilot import CMLE_MODE
from validation.estimand_distribution_study import (
    FACETS_MODE,
    FIT_MODEL,
    MML_FIXED_MODE,
    MML_FREE_MODE,
    PYTHON_JMLE_MODE,
    THRESHOLD_CONDITION,
    _frame_digest,
    _latent_rows,
    _run_one_attempt,
    generate_person_shapes,
)
from validation.facets_pcm_boundary_pilot import (
    constrained_adjacent_design_audit,
    person_rater_components,
)
from validation.facets_pcm_known_truth_smoke import (
    CATEGORIES,
    CRITERION_TRUTH,
    RATER_TRUTH,
    TASK_TRUTH,
    THRESHOLD_CONDITIONS,
    _frame_sha256,
    apply_threshold_condition,
)
from validation.facets_resilient_pair import (
    dependency_manifest,
    dependency_manifest_digest,
    fit_resilient_pair,
)
from validation.operating_characteristics_facets import sha256_file


SCHEMA_VERSION = "known_assignment_multivector_preflight_v1"
PLAN_PATH = ROOT / "validation" / "known_assignment_multivector_preflight_plan_20260811.json"
REGISTRATION_PATH = ROOT / "validation" / "known_assignment_multivector_preflight_execution_registration_v2_20260811.json"
AMENDMENT_PATH = ROOT / "validation" / "known_assignment_multivector_preflight_evidence_key_amendment_20260811.json"
STUDY_DIR = ROOT / "validation" / "known_assignment_multivector_preflight4_20260811"
FACETS_EXE_DEFAULT = Path(r"C:\Facets\Facets.exe")
ATTEMPT_TYPES = (
    "RESILIENT_FACETS_PYTHON_JMLE_PCM",
    "PYTHON_MML_FIXED_SD08_Q31_PCM",
    "PYTHON_MML_FREE_SD_Q31_PCM",
    "PYTHON_EXACT_CMLE_PCM",
)
SCIENTIFIC_MODES = (
    FACETS_MODE,
    PYTHON_JMLE_MODE,
    MML_FIXED_MODE,
    MML_FREE_MODE,
    CMLE_MODE,
)
INPUT_FILES = (
    "manifest.csv",
    "generated_ratings.csv",
    "generated_facet_truth.csv",
    "generated_anchors.csv",
    "generated_pcm_threshold_truth.csv",
    "assignment_audit.csv",
    "assignment_draws.csv",
    "assignment_diagnostics.csv",
    "attempt_manifest.csv",
)


def _json_dump(path: Path, value: Any) -> None:
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _atomic_json(path: Path, value: Any) -> None:
    temporary = path.with_name(f"{path.name}.tmp.{os.getpid()}.{time.time_ns()}")
    _json_dump(temporary, value)
    os.replace(temporary, path)


def _parse_edges(signature: str) -> set[tuple[str, str]]:
    return {
        tuple(value.split("::", maxsplit=1))
        for value in str(signature).split("|")
        if value
    }


def _gamma_design(gamma: float) -> str:
    return f"gamma_{gamma:+.1f}".replace("+", "pos_").replace("-", "neg_").replace(".", "p")


def _verify_prior_evidence(plan: dict[str, Any]) -> None:
    prior = plan["prior_evidence"]
    pairs = (
        ("dense_dp_assessment", "dense_dp_assessment_sha256"),
        ("fixed_person_screening_assessment", "fixed_person_screening_assessment_sha256"),
        (
            "fixed_person_screening_contrast_summary",
            "fixed_person_screening_contrast_summary_sha256",
        ),
    )
    for path_key, hash_key in pairs:
        path = ROOT / prior[path_key]
        if sha256_file(path).lower() != str(prior[hash_key]).lower():
            raise ValueError(f"Prior-evidence hash mismatch: {path_key}")
    dense = json.loads((ROOT / prior["dense_dp_assessment"]).read_text(encoding="utf-8"))
    screening = json.loads(
        (ROOT / prior["fixed_person_screening_assessment"]).read_text(encoding="utf-8")
    )
    if not bool(dense.get("qualification_pass", False)):
        raise ValueError("Dense-DP qualification is not PASS.")
    screening_gates = screening.get("gates", {})
    screening_pass = bool(
        str(screening.get("decision", "")).lower() == "pass"
        and isinstance(screening_gates, dict)
        and screening_gates
        and all(bool(value) for value in screening_gates.values())
    )
    if not screening_pass:
        raise ValueError("Fixed-Person screening qualification is not PASS.")


def build_dependency_manifest() -> dict[str, str]:
    manifest = dependency_manifest()
    additions = {
        "validation/known_assignment_multivector_preflight.py": Path(__file__).resolve(),
        "validation/known_assignment_multivector_preflight_plan_20260811.json": PLAN_PATH,
        "validation/known_assignment_multivector_preflight_execution_registration_v2_20260811.json": REGISTRATION_PATH,
        "validation/known_assignment_multivector_preflight_evidence_key_amendment_20260811.json": AMENDMENT_PATH,
        "validation/known_assignment_dense_dp_20260811/assessment.json": ROOT
        / "validation"
        / "known_assignment_dense_dp_20260811"
        / "assessment.json",
    }
    for name, path in additions.items():
        if not path.is_file():
            raise FileNotFoundError(f"Preflight dependency missing: {path}")
        manifest[name] = sha256_file(path)
    return dict(sorted(manifest.items()))


def validate_registration() -> dict[str, Any]:
    value = json.loads(REGISTRATION_PATH.read_text(encoding="utf-8"))
    expected = {
        "plan_sha256": sha256_file(PLAN_PATH),
        "runner_sha256": sha256_file(Path(__file__).resolve()),
        "assignment_mechanism_sha256": sha256_file(
            ROOT / "mfrm_app" / "assignment_mechanism.py"
        ),
        "estimand_runner_sha256": sha256_file(
            ROOT / "validation" / "estimand_distribution_study.py"
        ),
        "resilient_pair_sha256": sha256_file(
            ROOT / "validation" / "facets_resilient_pair.py"
        ),
        "evidence_key_amendment_sha256": sha256_file(AMENDMENT_PATH),
    }
    for key, digest in expected.items():
        if str(value.get(key, "")).lower() != digest.lower():
            raise ValueError(f"Execution registration mismatch: {key}")
    if not bool(value.get("tests_passed_before_person_generation", False)):
        raise ValueError("Registration does not assert passing pre-generation tests.")
    return value


def generate_bundle(
    plan: dict[str, Any], *, dependency_digest: str
) -> dict[str, pd.DataFrame]:
    _verify_prior_evidence(plan)
    dgm = plan["data_generating_process"]
    assignment = plan["assignment"]
    vector_ids = [int(value) for value in dgm["person_vector_ids"]]
    person_seeds = [int(value) for value in dgm["person_vector_seeds"]]
    uniform_seeds = [int(value) for value in dgm["response_uniform_seeds"]]
    if not (len(vector_ids) == len(person_seeds) == len(uniform_seeds)):
        raise ValueError("Person-vector identifiers and seed lists differ in length.")

    manifests: list[dict[str, Any]] = []
    rating_parts: list[pd.DataFrame] = []
    truth_parts: list[pd.DataFrame] = []
    threshold_rows: list[dict[str, Any]] = []
    audit_rows: list[dict[str, Any]] = []
    draw_rows: list[dict[str, Any]] = []
    diagnostic_parts: list[pd.DataFrame] = []
    plan_hash = sha256_file(PLAN_PATH)
    row_count = (
        int(dgm["persons"])
        * len(RATER_TRUTH)
        * len(TASK_TRUTH)
        * len(CRITERION_TRUTH)
    )

    for vector_id, person_seed, uniform_seed in zip(
        vector_ids, person_seeds, uniform_seeds, strict=True
    ):
        abilities = generate_person_shapes(person_seed, int(dgm["persons"]))["normal"]
        people = [f"P{index:03d}" for index in range(1, len(abilities) + 1)]
        person_coordinates = dict(zip(people, map(float, abilities), strict=True))
        person_frame = pd.DataFrame(
            {"Person": people, "Theta": np.asarray(abilities, dtype=float)}
        )
        person_hash = _frame_digest(person_frame)
        row_uniforms = np.random.default_rng(uniform_seed).random(row_count)
        uniform_hash = hashlib.sha256(row_uniforms.tobytes()).hexdigest()
        latent = _latent_rows(abilities, row_uniforms)
        complete = apply_threshold_condition(
            latent, THRESHOLD_CONDITIONS[THRESHOLD_CONDITION]
        )
        complete_hash = _frame_sha256(
            complete, ["Person", "Rater", "Task", "Criterion", "Score"]
        )
        for gamma_value in assignment["gammas"]:
            gamma = float(gamma_value)
            assignment_seed = int(
                assignment["assignment_seeds_by_person_vector"][str(vector_id)][
                    str(gamma)
                ]
            )
            sampled = sample_degree_conditioned_assignment_dense_dp(
                person_degrees={
                    person: int(assignment["person_degree"]) for person in people
                },
                rater_degrees={
                    rater: int(assignment["rater_degree"]) for rater in RATER_TRUTH
                },
                person_coordinates=person_coordinates,
                rater_coordinates={
                    str(key): float(value) for key, value in dgm["raters"].items()
                },
                gamma=gamma,
                n_samples=1,
                seed=assignment_seed,
                require_connected=bool(assignment["connected"]),
                max_dense_dp_cells=int(assignment["max_dense_dp_cells"]),
                retain_edge_signatures=True,
            )
            draw = sampled["samples"].iloc[0]
            edges = _parse_edges(str(draw["EdgeSignature"]))
            edge_frame = pd.DataFrame(sorted(edges), columns=["Person", "Rater"])
            ratings = edge_frame.merge(
                complete,
                on=["Person", "Rater"],
                how="left",
                validate="one_to_many",
            )[["Person", "Rater", "Task", "Criterion", "Score"]]
            design = _gamma_design(gamma)
            run_id = f"{design}__vector-{vector_id:02d}"
            components = person_rater_components(ratings)
            rank = constrained_adjacent_design_audit(ratings, FIT_MODEL)
            person_degree = (
                ratings[["Person", "Rater"]]
                .drop_duplicates()
                .groupby("Person")["Rater"]
                .nunique()
            )
            rater_persons = (
                ratings[["Person", "Rater"]]
                .drop_duplicates()
                .groupby("Rater")["Person"]
                .nunique()
            )
            rater_rows = ratings.groupby("Rater").size()
            context_rows = ratings.groupby(["Rater", "Task", "Criterion"]).size()
            category_counts = (
                ratings.groupby(["Criterion", "Score"])
                .size()
                .reindex(
                    pd.MultiIndex.from_product(
                        [list(CRITERION_TRUTH), range(CATEGORIES)],
                        names=["Criterion", "Score"],
                    ),
                    fill_value=0,
                )
            )
            invariants = {
                "Rows": len(ratings)
                == int(plan["locked_design_invariants"]["rows_per_dataset"]),
                "Persons": ratings["Person"].nunique() == int(dgm["persons"]),
                "Raters": ratings["Rater"].nunique() == len(RATER_TRUTH),
                "RatersPerPerson": bool(
                    person_degree.eq(int(assignment["person_degree"])).all()
                ),
                "PersonsPerRater": bool(
                    rater_persons.eq(int(assignment["rater_degree"])).all()
                ),
                "RowsPerRater": bool(
                    rater_rows.eq(
                        int(plan["locked_design_invariants"]["rating_rows_per_rater"])
                    ).all()
                ),
                "ContextRowsPerRater": bool(
                    context_rows.eq(int(assignment["rater_degree"])).all()
                ),
                "PersonRaterComponents": len(components) == 1,
                "ConstrainedPCMNullity": int(rank["Nullity"]) == 0,
                "CriterionCategorySupport": bool(category_counts.gt(0).all()),
                "DenseAssignmentAvailable": bool(sampled["available"]),
            }
            manifests.append(
                {
                    "RunId": run_id,
                    "ConditionId": f"{design}__fresh_normal_vector",
                    "Design": design,
                    "Gamma": gamma,
                    "PersonVector": vector_id,
                    "PersonDistribution": "fresh_standardized_normal_vector",
                    "TruthBias": 0.0,
                    "Replicate": vector_id,
                    "Seed": person_seed,
                    "UniformSeed": uniform_seed,
                    "AssignmentSeed": assignment_seed,
                    "Categories": CATEGORIES,
                    "ThresholdCondition": THRESHOLD_CONDITION,
                    "PersonMeanRealized": float(np.mean(abilities)),
                    "PersonSDRealized": float(np.std(abilities, ddof=0)),
                    "PersonVectorSHA256": person_hash,
                    "SharedUniformSHA256": uniform_hash,
                    "CompleteResponseSHA256": complete_hash,
                    "AssignmentEdgeSHA256": _frame_digest(
                        edge_frame.sort_values(["Person", "Rater"])
                    ),
                    "AssignmentStatistic": float(draw["Statistic"]),
                    "AssignmentCorrelation": float(draw["AssignmentCorrelation"]),
                    "ExpectedRows": len(ratings),
                    "ExpectedStructuralNullity": int(rank["Nullity"]),
                    "PersonRaterComponents": len(components),
                    "MinimumCriterionCategoryCount": int(category_counts.min()),
                    "CompleteCriterionCategorySupport": bool(category_counts.gt(0).all()),
                }
            )
            for key, passed in invariants.items():
                audit_rows.append(
                    {
                        "RunId": run_id,
                        "PersonVector": vector_id,
                        "Gamma": gamma,
                        "AuditKey": key,
                        "Passed": bool(passed),
                    }
                )
            draw_rows.append(
                {
                    "RunId": run_id,
                    "PersonVector": vector_id,
                    "Gamma": gamma,
                    "Sample": int(draw["Sample"]),
                    "Statistic": float(draw["Statistic"]),
                    "AssignmentCorrelation": float(draw["AssignmentCorrelation"]),
                    "EdgeSignature": str(draw["EdgeSignature"]),
                }
            )
            diagnostic = sampled["diagnostics"].copy()
            diagnostic.insert(0, "RunId", run_id)
            diagnostic.insert(1, "PersonVector", vector_id)
            diagnostic["Available"] = bool(sampled["available"])
            diagnostic_parts.append(diagnostic)
            run_ratings = ratings.copy()
            run_ratings.insert(0, "RunId", run_id)
            rating_parts.append(run_ratings)
            truth_rows = [
                {
                    "RunId": run_id,
                    "Facet": "Person",
                    "Level": person,
                    "Truth": theta,
                }
                for person, theta in person_coordinates.items()
            ]
            for facet, values in (
                ("Rater", RATER_TRUTH),
                ("Task", TASK_TRUTH),
                ("Criterion", CRITERION_TRUTH),
            ):
                truth_rows.extend(
                    {
                        "RunId": run_id,
                        "Facet": facet,
                        "Level": level,
                        "Truth": float(value),
                    }
                    for level, value in values.items()
                )
            truth_parts.append(pd.DataFrame(truth_rows))
            for criterion, vector in THRESHOLD_CONDITIONS[THRESHOLD_CONDITION].items():
                for category, value in enumerate(vector, start=1):
                    threshold_rows.append(
                        {
                            "RunId": run_id,
                            "ConditionId": f"{design}__fresh_normal_vector",
                            "Design": design,
                            "Gamma": gamma,
                            "PersonVector": vector_id,
                            "PersonDistribution": "fresh_standardized_normal_vector",
                            "Replicate": vector_id,
                            "StepFacetLevel": criterion,
                            "Category": category,
                            "ThresholdTruth": float(value),
                        }
                    )

    manifest = pd.DataFrame(manifests)
    ratings = pd.concat(rating_parts, ignore_index=True)
    truth = pd.concat(truth_parts, ignore_index=True)
    thresholds = pd.DataFrame(threshold_rows)
    attempts: list[dict[str, Any]] = []
    ordinal = 0
    for row in manifest.itertuples(index=False):
        run_ratings = ratings.loc[ratings["RunId"].eq(row.RunId)]
        run_truth = truth.loc[truth["RunId"].eq(row.RunId)]
        run_thresholds = thresholds.loc[thresholds["RunId"].eq(row.RunId)]
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
                f"{attempt_id}|{input_hash}|{plan_hash}|{dependency_digest}".encode(
                    "utf-8"
                )
            ).hexdigest()
            attempts.append(
                {
                    "AttemptOrdinal": ordinal,
                    "AttemptId": attempt_id,
                    "RunId": row.RunId,
                    "AttemptType": attempt_type,
                    "PersonVector": int(row.PersonVector),
                    "Replicate": int(row.Replicate),
                    "Gamma": float(row.Gamma),
                    "Design": row.Design,
                    "PersonDistribution": row.PersonDistribution,
                    "RunInputSHA256": input_hash,
                    "DependencyManifestSHA256": dependency_digest,
                    "AttemptFingerprint": fingerprint,
                }
            )
            ordinal += 1
    return {
        "manifest.csv": manifest,
        "generated_ratings.csv": ratings,
        "generated_facet_truth.csv": truth,
        "generated_anchors.csv": pd.DataFrame(
            columns=["RunId", "Facet", "Level", "Anchor"]
        ),
        "generated_pcm_threshold_truth.csv": thresholds,
        "assignment_audit.csv": pd.DataFrame(audit_rows),
        "assignment_draws.csv": pd.DataFrame(draw_rows),
        "assignment_diagnostics.csv": pd.concat(diagnostic_parts, ignore_index=True),
        "attempt_manifest.csv": pd.DataFrame(attempts),
    }


def prepare_study(*, facets_exe: Path) -> dict[str, Any]:
    if STUDY_DIR.exists():
        raise FileExistsError(f"Refusing to overwrite frozen study: {STUDY_DIR}")
    plan = json.loads(PLAN_PATH.read_text(encoding="utf-8"))
    validate_registration()
    _verify_prior_evidence(plan)
    facets_exe = facets_exe.resolve()
    if not facets_exe.is_file():
        raise FileNotFoundError(f"FACETS executable missing: {facets_exe}")
    dependencies = build_dependency_manifest()
    dependency_digest = dependency_manifest_digest(dependencies)
    bundle = generate_bundle(plan, dependency_digest=dependency_digest)
    STUDY_DIR.mkdir(parents=False, exist_ok=False)
    input_dir = STUDY_DIR / "retained_input"
    input_dir.mkdir()
    for filename, frame in bundle.items():
        frame.to_csv(input_dir / filename, index=False, lineterminator="\n")
    input_hashes = {filename: sha256_file(input_dir / filename) for filename in INPUT_FILES}
    manifest = bundle["manifest.csv"]
    audit = bundle["assignment_audit.csv"]
    diagnostics = bundle["assignment_diagnostics.csv"]
    checks = {
        "datasets_12": len(manifest) == 12,
        "attempts_48": len(bundle["attempt_manifest.csv"]) == 48,
        "ratings_11520": len(bundle["generated_ratings.csv"]) == 11_520,
        "all_assignment_and_design_audits": bool(audit["Passed"].all()),
        "all_category_support": bool(manifest["CompleteCriterionCategorySupport"].all()),
        "all_nullity_zero": bool(manifest["ExpectedStructuralNullity"].eq(0).all()),
        "all_connected": bool(manifest["PersonRaterComponents"].eq(1).all()),
        "same_complete_response_within_vector": bool(
            manifest.groupby("PersonVector")["CompleteResponseSHA256"].nunique().eq(1).all()
        ),
        "four_unique_person_vectors": manifest["PersonVectorSHA256"].nunique() == 4,
        "three_gamma_conditions_each_vector": bool(
            manifest.groupby("PersonVector")["Gamma"].nunique().eq(3).all()
        ),
        "all_dense_draws_available": len(diagnostics) == 12
        and bool(diagnostics["Available"].all()),
        "dense_conditional_residual": bool(
            diagnostics["MaximumConditionalLogResidual"].le(1e-10).all()
        ),
        "dense_statistic_residual": bool(
            diagnostics["MaximumStatisticUpdateResidual"].le(1e-10).all()
        ),
        "assignment_correlations_finite": bool(
            np.isfinite(manifest["AssignmentCorrelation"].to_numpy(dtype=float)).all()
        ),
    }
    identity = {
        "schema_version": SCHEMA_VERSION,
        "phase": "registered_nonconfirmatory_preflight",
        "plan_sha256": sha256_file(PLAN_PATH),
        "registration_sha256": sha256_file(REGISTRATION_PATH),
        "runner_sha256": sha256_file(Path(__file__).resolve()),
        "dependency_manifest_sha256": dependency_digest,
        "facets_executable": str(facets_exe),
        "facets_executable_sha256": sha256_file(facets_exe),
        "retained_input_sha256": input_hashes,
        "datasets": len(manifest),
        "attempts": len(bundle["attempt_manifest.csv"]),
        "python": sys.version,
        "platform": platform.platform(),
        "checks": checks,
        "all_checks_pass": bool(all(checks.values())),
        "screening_observations_pooled": False,
        "confirmatory_claims_allowed": False,
        "claim_limit": plan["claim_limit"],
    }
    _json_dump(STUDY_DIR / "dependency_manifest.json", dependencies)
    _json_dump(STUDY_DIR / "study_identity.json", identity)
    if not identity["all_checks_pass"]:
        raise RuntimeError(f"Frozen preflight input audit failed: {checks}")
    print(json.dumps(identity, ensure_ascii=False, sort_keys=True))
    return identity


def validate_study_identity() -> dict[str, Any]:
    identity = json.loads((STUDY_DIR / "study_identity.json").read_text(encoding="utf-8"))
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
        raise ValueError("Execution dependency manifest changed.")
    retained = json.loads((STUDY_DIR / "dependency_manifest.json").read_text(encoding="utf-8"))
    if retained != dependencies:
        raise ValueError("Retained dependency manifest differs from current dependencies.")
    facets_exe = Path(identity["facets_executable"])
    if sha256_file(facets_exe) != identity["facets_executable_sha256"]:
        raise ValueError("FACETS executable changed.")
    for filename, digest in identity["retained_input_sha256"].items():
        if sha256_file(STUDY_DIR / "retained_input" / filename) != digest:
            raise ValueError(f"Retained input changed: {filename}")
    if not bool(identity.get("all_checks_pass", False)):
        raise ValueError("Prepared input audit was not PASS.")
    return identity


def _artifact_root(attempt: pd.Series) -> Path:
    return STUDY_DIR / "work" / f"{int(attempt['AttemptOrdinal']):05d}"


def _completion_path(attempt: pd.Series) -> Path:
    return STUDY_DIR / "attempts" / f"{int(attempt['AttemptOrdinal']):05d}" / "completion.json"


def _artifact_hashes(root: Path) -> dict[str, str]:
    return {
        path.relative_to(root).as_posix(): sha256_file(path)
        for path in sorted(root.rglob("*"))
        if path.is_file()
    }


def _validate_completion(
    completion: dict[str, Any], attempt: pd.Series, identity: dict[str, Any]
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
    root = STUDY_DIR / str(completion["artifact_root"])
    for relative, digest in completion.get("artifact_sha256", {}).items():
        path = root / relative
        if not path.is_file() or sha256_file(path) != digest:
            raise ValueError(f"Completion artifact mismatch: {path}")


def _write_standard(
    root: Path,
    runs: pd.DataFrame,
    recovery: pd.DataFrame,
    thresholds: pd.DataFrame,
    constraints: Any,
) -> None:
    runs.to_csv(root / "run_ledger.csv", index=False, lineterminator="\n")
    recovery.to_csv(root / "recovery.csv", index=False, lineterminator="\n")
    thresholds.to_csv(root / "thresholds.csv", index=False, lineterminator="\n")
    pd.DataFrame(constraints).to_csv(root / "constraints.csv", index=False, lineterminator="\n")


def run_shard(
    *, shard_index: int, shard_count: int, resume: bool, timeout_seconds: float
) -> dict[str, int]:
    identity = validate_study_identity()
    if shard_count < 1 or shard_index < 0 or shard_index >= shard_count:
        raise ValueError("Require 0 <= shard-index < shard-count")
    input_dir = STUDY_DIR / "retained_input"
    attempts = pd.read_csv(input_dir / "attempt_manifest.csv")
    ordinals = pd.to_numeric(attempts["AttemptOrdinal"], errors="raise").astype(int)
    selected = attempts.loc[ordinals.mod(shard_count).eq(shard_index)].copy()
    manifest = pd.read_csv(input_dir / "manifest.csv").set_index("RunId", drop=False)
    ratings_all = pd.read_csv(input_dir / "generated_ratings.csv")
    truth_all = pd.read_csv(input_dir / "generated_facet_truth.csv")
    threshold_all = pd.read_csv(input_dir / "generated_pcm_threshold_truth.csv")
    (STUDY_DIR / "attempts").mkdir(exist_ok=True)
    (STUDY_DIR / "locks").mkdir(exist_ok=True)
    summary = {"assigned": len(selected), "completed_now": 0, "skipped": 0, "failed": 0}
    for sequence, (_, attempt) in enumerate(selected.iterrows(), start=1):
        completion_path = _completion_path(attempt)
        if completion_path.is_file():
            completion = json.loads(completion_path.read_text(encoding="utf-8"))
            _validate_completion(completion, attempt, identity)
            if not resume:
                raise FileExistsError(f"Attempt already complete: {attempt['AttemptId']}")
            summary["skipped"] += 1
            print(f"[{sequence}/{len(selected)}] skip {attempt['AttemptId']}", flush=True)
            continue
        artifact_root = _artifact_root(attempt)
        if artifact_root.exists():
            raise FileExistsError(f"Unmarked artifact root requires adjudication: {artifact_root}")
        completion_path.parent.mkdir(parents=True, exist_ok=True)
        run_id = str(attempt["RunId"])
        manifest_row = manifest.loc[run_id]
        ratings = ratings_all.loc[ratings_all["RunId"].eq(run_id)].drop(columns="RunId")
        truth = truth_all.loc[truth_all["RunId"].eq(run_id)]
        thresholds = threshold_all.loc[threshold_all["RunId"].eq(run_id)]
        execution_completed = False
        statistical_ready = False
        calibration_ready: bool | None = None
        failure_reason = ""
        try:
            if str(attempt["AttemptType"]) == "RESILIENT_FACETS_PYTHON_JMLE_PCM":
                outcome = fit_resilient_pair(
                    manifest_row,
                    ratings,
                    truth,
                    thresholds,
                    facets_exe=Path(identity["facets_executable"]),
                    output_dir=artifact_root,
                    lock_path=STUDY_DIR / "locks" / "facets_global.lock",
                    timeout_seconds=timeout_seconds,
                )
                runs = pd.read_csv(artifact_root / "combined_run_ledger.csv")
                recovery = pd.read_csv(artifact_root / "combined_recovery.csv")
                threshold_rows = pd.read_csv(artifact_root / "combined_thresholds.csv")
                _write_standard(artifact_root, runs, recovery, threshold_rows, [])
                statistical_ready = bool(outcome["statistical_evidence_ready"])
                calibration_ready = bool(outcome["calibration_ready"])
            else:
                artifact_root.mkdir(parents=True)
                runs, recovery, threshold_rows, constraints = _run_one_attempt(
                    attempt,
                    manifest_row,
                    ratings,
                    truth,
                    thresholds,
                    facets_exe=Path(identity["facets_executable"]),
                    stage_dir=artifact_root,
                    timeout_seconds=timeout_seconds,
                )
                _write_standard(artifact_root, runs, recovery, threshold_rows, constraints)
                statistical_ready = bool(
                    runs["IncludedInStudy"].fillna(False).astype(bool).all()
                )
            execution_completed = True
        except Exception as exc:  # planned denominator is retained
            failure_reason = f"{type(exc).__name__}: {exc}"
            artifact_root.mkdir(parents=True, exist_ok=True)
            _json_dump(
                artifact_root / "unhandled_failure.json",
                {"attempt_id": str(attempt["AttemptId"]), "failure_reason": failure_reason},
            )
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
            "execution_completed": execution_completed,
            "statistical_evidence_ready": statistical_ready,
            "facets_calibration_ready": calibration_ready,
            "failure_reason": failure_reason,
            "artifact_root": artifact_root.relative_to(STUDY_DIR).as_posix(),
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


def _read_nonempty(path: Path) -> pd.DataFrame:
    if not path.is_file() or path.stat().st_size <= 1:
        return pd.DataFrame()
    try:
        return pd.read_csv(path)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def _attach_context(frame: pd.DataFrame, manifest: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "RunId",
        "PersonVector",
        "Gamma",
        "Design",
        "Replicate",
        "AssignmentStatistic",
        "AssignmentCorrelation",
    ]
    payload = frame.drop(columns=[column for column in columns if column != "RunId"], errors="ignore")
    return payload.merge(manifest[columns], on="RunId", how="left", validate="many_to_one")


def summarize_preflight_diagnostics(
    *, rater_loss: pd.DataFrame, rater_errors: pd.DataFrame, runs: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Construct the four registered n=4 diagnostics without inferential tests."""

    loss_wide = rater_loss.pivot(
        index=["PersonVector", "EstimatorMode"], columns="Gamma", values="RaterRMSE"
    ).reset_index()
    for gamma in (-0.8, 0.0, 0.8):
        if gamma not in loss_wide.columns:
            loss_wide[gamma] = np.nan
    loss_wide["SymmetricStressRaterRMSE"] = (
        0.5 * (loss_wide[-0.8] + loss_wide[0.8]) - loss_wide[0.0]
    )

    severity = {str(key): float(value) for key, value in RATER_TRUTH.items()}
    slope_rows: list[dict[str, Any]] = []
    for (vector, mode, gamma), group in rater_errors.groupby(
        ["PersonVector", "EstimatorMode", "Gamma"], sort=True
    ):
        x = group["Level"].astype(str).map(severity).to_numpy(dtype=float)
        y = pd.to_numeric(group["ErrorAligned"], errors="coerce").to_numpy(dtype=float)
        slope_rows.append(
            {
                "PersonVector": int(vector),
                "EstimatorMode": str(mode),
                "Gamma": float(gamma),
                "RaterErrorSeveritySlope": float(np.polyfit(x, y, 1)[0])
                if len(x) == 4 and np.isfinite(x).all() and np.isfinite(y).all()
                else np.nan,
            }
        )
    slopes = pd.DataFrame(slope_rows)
    slope_wide = slopes.pivot(
        index=["PersonVector", "EstimatorMode"],
        columns="Gamma",
        values="RaterErrorSeveritySlope",
    ).reset_index()
    for gamma in (-0.8, 0.8):
        if gamma not in slope_wide.columns:
            slope_wide[gamma] = np.nan
    slope_wide["DirectionAlignedSlopeHalfDifference"] = 0.5 * (
        slope_wide[-0.8] - slope_wide[0.8]
    )

    free = runs.loc[runs["EstimatorMode"].astype(str).eq(MML_FREE_MODE)].copy()
    free["EstimatedPopulationSD"] = pd.to_numeric(
        free["EstimatedPopulationSD"], errors="coerce"
    )
    free_wide = free.pivot(
        index="PersonVector", columns="Gamma", values="EstimatedPopulationSD"
    ).reset_index()
    for gamma in (-0.8, 0.0, 0.8):
        if gamma not in free_wide.columns:
            free_wide[gamma] = np.nan
    free_wide["SymmetricStressSDShift"] = (
        0.5 * (free_wide[-0.8] + free_wide[0.8]) - free_wide[0.0]
    )

    diagnostics: list[dict[str, Any]] = []
    specifications = (
        (
            "PF1_FREE_MML_SYMMETRIC_STRESS_RATER_RMSE",
            loss_wide.loc[loss_wide["EstimatorMode"].eq(MML_FREE_MODE)],
            "SymmetricStressRaterRMSE",
            "positive",
            3,
        ),
        (
            "PF2_FIXED_MML_SYMMETRIC_STRESS_RATER_RMSE",
            loss_wide.loc[loss_wide["EstimatorMode"].eq(MML_FIXED_MODE)],
            "SymmetricStressRaterRMSE",
            "positive",
            3,
        ),
        (
            "PF3_FREE_MML_SYMMETRIC_STRESS_SD_SHIFT",
            free_wide,
            "SymmetricStressSDShift",
            "negative",
            3,
        ),
        (
            "PF4_FREE_MML_DIRECTION_ALIGNED_RATER_SLOPE",
            slope_wide.loc[slope_wide["EstimatorMode"].eq(MML_FREE_MODE)],
            "DirectionAlignedSlopeHalfDifference",
            "positive",
            4,
        ),
    )
    for diagnostic_id, frame, value_column, direction, required_count in specifications:
        for row in frame[["PersonVector", value_column]].itertuples(index=False):
            diagnostics.append(
                {
                    "DiagnosticId": diagnostic_id,
                    "PersonVector": int(row[0]),
                    "Value": float(row[1]),
                    "RegisteredDirection": direction,
                    "RequiredDirectionalCount": required_count,
                    "ConfirmatoryClaimAllowed": False,
                }
            )
    diagnostic_table = pd.DataFrame(diagnostics)
    summary_rows: list[dict[str, Any]] = []
    for diagnostic_id, group in diagnostic_table.groupby("DiagnosticId", sort=True):
        values = group["Value"].to_numpy(dtype=float)
        direction = str(group["RegisteredDirection"].iloc[0])
        required = int(group["RequiredDirectionalCount"].iloc[0])
        directional_count = int(np.sum(values > 0)) if direction == "positive" else int(np.sum(values < 0))
        mean = float(np.mean(values)) if np.isfinite(values).all() else np.nan
        mean_direction = bool(mean > 0) if direction == "positive" else bool(mean < 0)
        summary_rows.append(
            {
                "DiagnosticId": diagnostic_id,
                "N": len(values),
                "Mean": mean,
                "Minimum": float(np.min(values)) if np.isfinite(values).all() else np.nan,
                "Maximum": float(np.max(values)) if np.isfinite(values).all() else np.nan,
                "RegisteredDirection": direction,
                "DirectionalCount": directional_count,
                "RequiredDirectionalCount": required,
                "AdvancementSignal": bool(
                    len(values) == 4
                    and np.isfinite(values).all()
                    and mean_direction
                    and directional_count >= required
                ),
                "PValueComputed": False,
                "ConfirmatoryClaimAllowed": False,
            }
        )
    return diagnostic_table, pd.DataFrame(summary_rows), loss_wide, slope_wide


def aggregate_study() -> dict[str, Any]:
    if (STUDY_DIR / "aggregate").exists():
        raise FileExistsError("Refusing to overwrite frozen aggregate.")
    identity = validate_study_identity()
    input_dir = STUDY_DIR / "retained_input"
    attempts = pd.read_csv(input_dir / "attempt_manifest.csv")
    manifest = pd.read_csv(input_dir / "manifest.csv")
    run_parts: list[pd.DataFrame] = []
    recovery_parts: list[pd.DataFrame] = []
    threshold_parts: list[pd.DataFrame] = []
    constraint_parts: list[pd.DataFrame] = []
    outcomes: list[dict[str, Any]] = []
    marker_hashes: dict[str, str] = {}
    for _, attempt in attempts.iterrows():
        marker = _completion_path(attempt)
        if not marker.is_file():
            raise FileNotFoundError(f"Missing completion marker: {attempt['AttemptId']}")
        completion = json.loads(marker.read_text(encoding="utf-8"))
        _validate_completion(completion, attempt, identity)
        marker_hashes[str(attempt["AttemptId"])] = sha256_file(marker)
        root = STUDY_DIR / str(completion["artifact_root"])
        for filename, target in (
            ("run_ledger.csv", run_parts),
            ("recovery.csv", recovery_parts),
            ("thresholds.csv", threshold_parts),
            ("constraints.csv", constraint_parts),
        ):
            frame = _read_nonempty(root / filename)
            if not frame.empty:
                frame.insert(0, "AttemptId", str(attempt["AttemptId"]))
                target.append(frame)
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
        metrics_path = root / "resilient_pair_metrics.json"
        if metrics_path.is_file():
            outcome.update(json.loads(metrics_path.read_text(encoding="utf-8")))
        outcomes.append(outcome)
    runs = _attach_context(pd.concat(run_parts, ignore_index=True, sort=False), manifest)
    recovery = _attach_context(
        pd.concat(recovery_parts, ignore_index=True, sort=False), manifest
    )
    thresholds = _attach_context(
        pd.concat(threshold_parts, ignore_index=True, sort=False), manifest
    )
    constraints = (
        _attach_context(pd.concat(constraint_parts, ignore_index=True, sort=False), manifest)
        if constraint_parts
        else pd.DataFrame()
    )
    outcome_table = pd.DataFrame(outcomes)

    rater_errors = recovery.loc[
        recovery["EstimatorMode"].astype(str).isin(SCIENTIFIC_MODES)
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
    diagnostics, diagnostic_summary, loss_wide, slope_wide = summarize_preflight_diagnostics(
        rater_loss=rater_loss, rater_errors=rater_errors, runs=runs
    )

    jmle = outcome_table.loc[
        outcome_table["AttemptType"].eq("RESILIENT_FACETS_PYTHON_JMLE_PCM")
    ]
    mml = outcome_table.loc[outcome_table["AttemptType"].str.contains("MML")]
    cmle = outcome_table.loc[
        outcome_table["AttemptType"].eq("PYTHON_EXACT_CMLE_PCM")
    ]
    input_audit = pd.read_csv(input_dir / "assignment_audit.csv")
    input_diagnostics = pd.read_csv(input_dir / "assignment_diagnostics.csv")
    operational_gates = {
        "completion_markers_48": len(marker_hashes) == 48,
        "all_attempts_execution_completed": bool(outcome_table["ExecutionCompleted"].all()),
        "all_native_estimator_evidence_ready": bool(
            outcome_table["StatisticalEvidenceReady"].all()
        ),
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
    facets_runs = runs.loc[runs["EstimatorMode"].astype(str).eq(FACETS_MODE)]
    calibration = {
        "pairs": len(jmle),
        "ready": int(jmle["FACETSCalibrationReady"].fillna(False).astype(bool).sum()),
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
    }
    metrics = {
        "schema_version": SCHEMA_VERSION,
        "phase": "nonconfirmatory_preflight",
        "operational_gates": {
            key: bool(value) for key, value in operational_gates.items()
        },
        "qualification_pass": bool(all(operational_gates.values())),
        "primary_advancement_signal": primary_advancement,
        "eligible_to_freeze_separate_confirmation": bool(
            all(operational_gates.values()) and primary_advancement
        ),
        "facets_calibration": calibration,
        "screening_observations_pooled": False,
        "preflight_observations_may_enter_confirmation": False,
        "p_values_computed": False,
        "confirmatory_claims_allowed": False,
        "estimator_ranking_constructed": False,
        "cross_basis_likelihood_comparison": "prohibited",
        "facets_fit_precision_boundary": "ordinary displayed two-decimal fit fields were not calculation inputs",
        "claim_limit": identity["claim_limit"],
    }

    aggregate_dir = STUDY_DIR / "aggregate"
    aggregate_dir.mkdir()
    outputs = {
        "attempt_outcomes.csv": outcome_table,
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
    }
    for filename, frame in outputs.items():
        frame.to_csv(aggregate_dir / filename, index=False, lineterminator="\n")
    _json_dump(aggregate_dir / "assessment.json", metrics)
    marker_digest = hashlib.sha256(
        json.dumps(marker_hashes, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    artifacts = tuple(outputs) + ("assessment.json",)
    _atomic_json(
        aggregate_dir / "aggregate_identity.json",
        {
            "schema_version": f"{SCHEMA_VERSION}_aggregate_identity_v1",
            "study_identity_sha256": sha256_file(STUDY_DIR / "study_identity.json"),
            "completion_marker_set_sha256": marker_digest,
            "artifact_sha256": {
                filename: sha256_file(aggregate_dir / filename)
                for filename in artifacts
            },
        },
    )
    print(json.dumps(metrics, ensure_ascii=False, sort_keys=True))
    return metrics


def main() -> None:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    prepare_parser = subparsers.add_parser("prepare")
    prepare_parser.add_argument("--facets-exe", type=Path, default=FACETS_EXE_DEFAULT)
    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("--shard-index", type=int, default=0)
    run_parser.add_argument("--shard-count", type=int, default=1)
    run_parser.add_argument("--resume", action="store_true")
    run_parser.add_argument("--timeout-seconds", type=float, default=120.0)
    subparsers.add_parser("aggregate")
    args = parser.parse_args()
    if args.command == "prepare":
        prepare_study(facets_exe=args.facets_exe)
    elif args.command == "run":
        run_shard(
            shard_index=args.shard_index,
            shard_count=args.shard_count,
            resume=args.resume,
            timeout_seconds=args.timeout_seconds,
        )
    else:
        aggregate_study()


if __name__ == "__main__":
    main()
