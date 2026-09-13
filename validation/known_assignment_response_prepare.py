"""Prepare and audit the frozen known-assignment response-screening inputs."""

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

from mfrm_app.operating_characteristics import deterministic_replicate_seed
from validation.estimand_distribution_study import (
    FIT_MODEL,
    THRESHOLD_CONDITION,
    _frame_digest,
    _latent_rows,
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
from validation.known_assignment_large_design_dp_shard import _sha256, _slug


PLAN_PATH = ROOT / "validation" / "known_assignment_response_screening_plan_20260811.json"
ASSIGNMENT_ROOT = ROOT / "validation" / "known_assignment_response_assignment_shards_20260811"
STUDY_DIR = ROOT / "validation" / "known_assignment_response_screening_20260811"
ATTEMPT_TYPES = (
    "RESILIENT_FACETS_PYTHON_JMLE_PCM",
    "PYTHON_MML_FIXED_SD08_Q31_PCM",
    "PYTHON_MML_FREE_SD_Q31_PCM",
    "PYTHON_EXACT_CMLE_PCM",
)
INPUT_FILES = (
    "manifest.csv",
    "generated_ratings.csv",
    "generated_facet_truth.csv",
    "generated_anchors.csv",
    "generated_pcm_threshold_truth.csv",
    "assignment_audit.csv",
    "attempt_manifest.csv",
)


def _gamma_design(gamma: float) -> str:
    return f"gamma_{gamma:+.1f}".replace("+", "pos_").replace("-", "neg_").replace(".", "p")


def _load_assignment_shards(plan: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, str]]:
    edges: list[pd.DataFrame] = []
    draws: list[pd.DataFrame] = []
    hashes: dict[str, str] = {}
    for gamma_value in plan["assignment"]["gammas"]:
        gamma = float(gamma_value)
        root = ASSIGNMENT_ROOT / _slug(gamma)
        manifest_path = root / "manifest.json"
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if float(manifest["gamma"]) != gamma:
            raise ValueError(f"Assignment shard gamma mismatch: {root}")
        for name, expected in manifest["artifact_sha256"].items():
            if _sha256(root / name) != expected:
                raise ValueError(f"Assignment shard artifact mismatch: {root / name}")
        edges.append(pd.read_csv(root / "assignment_edges.csv"))
        draw = pd.read_csv(root / "assignment_draws.csv")
        draw["Gamma"] = gamma
        draws.append(draw)
        audit = pd.read_csv(root / "materialization_audit.csv")
        if not audit["Passed"].astype(bool).all():
            raise ValueError(f"Assignment shard materialization failed: {root}")
        hashes[_slug(gamma)] = _sha256(manifest_path)
    return pd.concat(edges, ignore_index=True), pd.concat(draws, ignore_index=True), hashes


def generate_bundle(plan: dict[str, Any]) -> tuple[dict[str, pd.DataFrame], dict[str, Any]]:
    assignment_edges, assignment_draws, shard_hashes = _load_assignment_shards(plan)
    coordinate_path = ROOT / plan["data_generating_process"]["fixed_person_coordinate_source"]
    coordinates = pd.read_csv(coordinate_path).sort_values("Person").reset_index(drop=True)
    abilities = coordinates["Theta"].to_numpy(dtype=float)
    manifests: list[dict[str, Any]] = []
    ratings_parts: list[pd.DataFrame] = []
    truth_parts: list[pd.DataFrame] = []
    threshold_rows: list[dict[str, Any]] = []
    audit_rows: list[dict[str, Any]] = []
    plan_hash = _sha256(PLAN_PATH)
    row_count = len(abilities) * len(RATER_TRUTH) * len(TASK_TRUTH) * len(CRITERION_TRUTH)
    for replicate in range(1, int(plan["evidence_status"]["replicates"]) + 1):
        uniform_seed = deterministic_replicate_seed(
            int(plan["data_generating_process"]["response_uniform_seed_base"]),
            "known-assignment-response-uniforms",
            replicate,
        )
        row_uniforms = np.random.default_rng(uniform_seed).random(row_count)
        latent = _latent_rows(abilities, row_uniforms)
        complete = apply_threshold_condition(latent, THRESHOLD_CONDITIONS[THRESHOLD_CONDITION])
        complete_hash = _frame_sha256(
            complete, ["Person", "Rater", "Task", "Criterion", "Score"]
        )
        uniform_hash = hashlib.sha256(row_uniforms.tobytes()).hexdigest()
        for gamma_value in plan["assignment"]["gammas"]:
            gamma = float(gamma_value)
            design = _gamma_design(gamma)
            run_id = f"{design}__fixedtheta::rep-{replicate:05d}"
            edges = assignment_edges.loc[
                assignment_edges["Gamma"].eq(gamma)
                & assignment_edges["Replicate"].eq(replicate),
                ["Person", "Rater"],
            ].copy()
            if len(edges) != int(plan["assignment"]["edges"]) or edges.duplicated().any():
                raise RuntimeError(f"Assignment edge contract failed: {run_id}")
            ratings = edges.merge(
                complete,
                on=["Person", "Rater"],
                how="left",
                validate="one_to_many",
            )
            ratings = ratings[["Person", "Rater", "Task", "Criterion", "Score"]]
            components = person_rater_components(ratings)
            rank = constrained_adjacent_design_audit(ratings, FIT_MODEL)
            person_degree = ratings[["Person", "Rater"]].drop_duplicates().groupby("Person")["Rater"].nunique()
            rater_persons = ratings[["Person", "Rater"]].drop_duplicates().groupby("Rater")["Person"].nunique()
            rater_rows = ratings.groupby("Rater").size()
            context_rows = ratings.groupby(["Rater", "Task", "Criterion"]).size()
            category_counts = (
                ratings.groupby(["Criterion", "Score"]).size().reindex(
                    pd.MultiIndex.from_product(
                        [list(CRITERION_TRUTH), range(CATEGORIES)],
                        names=["Criterion", "Score"],
                    ),
                    fill_value=0,
                )
            )
            invariants = {
                "Rows": len(ratings) == int(plan["locked_design_invariants"]["rows_per_dataset"]),
                "Persons": ratings["Person"].nunique() == len(abilities),
                "Raters": ratings["Rater"].nunique() == len(RATER_TRUTH),
                "RatersPerPerson": person_degree.eq(int(plan["assignment"]["person_degree"])).all(),
                "PersonsPerRater": rater_persons.eq(int(plan["assignment"]["rater_degree"])).all(),
                "RowsPerRater": rater_rows.eq(int(plan["locked_design_invariants"]["rating_rows_per_rater"])).all(),
                "ContextRowsPerRater": context_rows.eq(40).all(),
                "PersonRaterComponents": len(components) == 1,
                "ConstrainedPCMNullity": int(rank["Nullity"]) == 0,
                "CriterionCategorySupport": category_counts.gt(0).all(),
            }
            if not all(bool(value) for value in invariants.values()):
                raise RuntimeError(f"Prepared response design invariant failed: {run_id}: {invariants}")
            draw = assignment_draws.loc[
                assignment_draws["Gamma"].eq(gamma)
                & assignment_draws["Sample"].eq(replicate)
            ].iloc[0]
            manifests.append(
                {
                    "RunId": run_id,
                    "ConditionId": f"{design}__fixedtheta",
                    "Design": design,
                    "Gamma": gamma,
                    "PersonDistribution": "fixed_normal_vector",
                    "TruthBias": 0.0,
                    "Replicate": replicate,
                    "UniformSeed": int(uniform_seed),
                    "Categories": CATEGORIES,
                    "ThresholdCondition": THRESHOLD_CONDITION,
                    "PersonMeanRealized": float(np.mean(abilities)),
                    "PersonSDRealized": float(np.std(abilities, ddof=0)),
                    "ExpectedRows": len(ratings),
                    "ExpectedStructuralNullity": int(rank["Nullity"]),
                    "PersonRaterComponents": len(components),
                    "AssignmentStatistic": float(draw["Statistic"]),
                    "AssignmentCorrelation": float(draw["AssignmentCorrelation"]),
                    "MinimumCriterionCategoryCount": int(category_counts.min()),
                    "CompleteCriterionCategorySupport": bool(category_counts.gt(0).all()),
                    "SharedUniformSHA256": uniform_hash,
                    "CompleteResponseSHA256": complete_hash,
                    "AssignmentEdgeSHA256": _frame_digest(edges.sort_values(["Person", "Rater"])),
                }
            )
            for key, passed in invariants.items():
                audit_rows.append(
                    {
                        "RunId": run_id,
                        "Gamma": gamma,
                        "Replicate": replicate,
                        "AuditKey": key,
                        "Passed": bool(passed),
                    }
                )
            run_ratings = ratings.copy()
            run_ratings.insert(0, "RunId", run_id)
            ratings_parts.append(run_ratings)
            truth_rows = [
                {"RunId": run_id, "Facet": "Person", "Level": row.Person, "Truth": float(row.Theta)}
                for row in coordinates.itertuples(index=False)
            ]
            for facet, values in (
                ("Rater", RATER_TRUTH),
                ("Task", TASK_TRUTH),
                ("Criterion", CRITERION_TRUTH),
            ):
                truth_rows.extend(
                    {"RunId": run_id, "Facet": facet, "Level": level, "Truth": float(value)}
                    for level, value in values.items()
                )
            truth_parts.append(pd.DataFrame(truth_rows))
            for criterion, vector in THRESHOLD_CONDITIONS[THRESHOLD_CONDITION].items():
                for category, value in enumerate(vector, start=1):
                    threshold_rows.append(
                        {
                            "RunId": run_id,
                            "ConditionId": f"{design}__fixedtheta",
                            "Design": design,
                            "Gamma": gamma,
                            "PersonDistribution": "fixed_normal_vector",
                            "Replicate": replicate,
                            "StepFacetLevel": criterion,
                            "Category": category,
                            "ThresholdTruth": float(value),
                        }
                    )
    manifest = pd.DataFrame(manifests)
    ratings = pd.concat(ratings_parts, ignore_index=True)
    truth = pd.concat(truth_parts, ignore_index=True)
    thresholds = pd.DataFrame(threshold_rows)
    attempts: list[dict[str, Any]] = []
    ordinal = 0
    for row in manifest.itertuples(index=False):
        run_ratings = ratings[ratings["RunId"].eq(row.RunId)]
        run_truth = truth[truth["RunId"].eq(row.RunId)]
        run_thresholds = thresholds[thresholds["RunId"].eq(row.RunId)]
        input_hash = hashlib.sha256(
            (_frame_digest(run_ratings) + _frame_digest(run_truth) + _frame_digest(run_thresholds)).encode("ascii")
        ).hexdigest()
        for attempt_type in ATTEMPT_TYPES:
            attempt_id = f"{row.RunId}::{attempt_type}"
            fingerprint = hashlib.sha256(
                f"{attempt_id}|{input_hash}|{plan_hash}".encode("utf-8")
            ).hexdigest()
            attempts.append(
                {
                    "AttemptOrdinal": ordinal,
                    "AttemptId": attempt_id,
                    "RunId": row.RunId,
                    "AttemptType": attempt_type,
                    "Replicate": int(row.Replicate),
                    "Gamma": float(row.Gamma),
                    "Design": row.Design,
                    "PersonDistribution": row.PersonDistribution,
                    "RunInputSHA256": input_hash,
                    "AttemptFingerprint": fingerprint,
                }
            )
            ordinal += 1
    bundle = {
        "manifest.csv": manifest,
        "generated_ratings.csv": ratings,
        "generated_facet_truth.csv": truth,
        "generated_anchors.csv": pd.DataFrame(columns=["RunId", "Facet", "Level", "Anchor"]),
        "generated_pcm_threshold_truth.csv": thresholds,
        "assignment_audit.csv": pd.DataFrame(audit_rows),
        "attempt_manifest.csv": pd.DataFrame(attempts),
    }
    metadata = {"assignment_shard_manifest_sha256": shard_hashes}
    return bundle, metadata


def prepare() -> dict[str, Any]:
    if STUDY_DIR.exists():
        raise FileExistsError(f"Refusing to overwrite study directory: {STUDY_DIR}")
    plan = json.loads(PLAN_PATH.read_text(encoding="utf-8"))
    bundle, metadata = generate_bundle(plan)
    STUDY_DIR.mkdir(parents=False, exist_ok=False)
    input_dir = STUDY_DIR / "retained_input"
    input_dir.mkdir()
    for name, table in bundle.items():
        table.to_csv(input_dir / name, index=False, lineterminator="\n")
    input_hashes = {name: _sha256(input_dir / name) for name in INPUT_FILES}
    manifest = bundle["manifest.csv"]
    ratings = bundle["generated_ratings.csv"]
    audit = bundle["assignment_audit.csv"]
    checks = {
        "datasets_30": len(manifest) == 30,
        "attempts_120": len(bundle["attempt_manifest.csv"]) == 120,
        "ratings_28800": len(ratings) == 28_800,
        "all_assignment_audits": bool(audit["Passed"].all()),
        "all_category_support": bool(manifest["CompleteCriterionCategorySupport"].all()),
        "all_nullity_zero": bool(manifest["ExpectedStructuralNullity"].eq(0).all()),
        "all_connected": bool(manifest["PersonRaterComponents"].eq(1).all()),
        "shared_complete_response_within_replicate": bool(
            manifest.groupby("Replicate")["CompleteResponseSHA256"].nunique().eq(1).all()
        ),
        "fixed_person_mean": abs(float(manifest["PersonMeanRealized"].iloc[0])) <= 1e-12,
        "fixed_person_sd": abs(float(manifest["PersonSDRealized"].iloc[0]) - 0.8) <= 1e-12,
        "gamma_correlation_order": bool(
            manifest.groupby("Gamma")["AssignmentCorrelation"].mean().sort_index().is_monotonic_increasing
        ),
    }
    identity = {
        "schema_version": "known_assignment_response_screening_input_identity_v1",
        "plan": str(PLAN_PATH.relative_to(ROOT)),
        "plan_sha256": _sha256(PLAN_PATH),
        "prepare_runner_sha256": _sha256(Path(__file__).resolve()),
        "retained_input_sha256": input_hashes,
        **metadata,
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }
    (STUDY_DIR / "input_identity.json").write_text(
        json.dumps(identity, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    if not identity["all_checks_pass"]:
        raise RuntimeError(f"Prepared input audit failed: {checks}")
    print(json.dumps(identity, ensure_ascii=False, indent=2))
    return identity


if __name__ == "__main__":
    prepare()
