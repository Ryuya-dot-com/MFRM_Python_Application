from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from validation.informative_assignment_screening import (
    BASE_SEED,
    REPO_ROOT,
    generate_screening_bundle,
    screening_contrast_summary,
    select_shard_attempts,
)
from validation.operating_characteristics_facets import sha256_file


def _registered_files(tmp_path: Path, *, start: int = 221, count: int = 1) -> tuple[Path, Path]:
    plan = json.loads(
        (
            REPO_ROOT
            / "validation"
            / "informative_assignment_screening_plan_20260811.json"
        ).read_text(encoding="utf-8")
    )
    plan["evidence_status"]["replicate_range"] = [start, start + count - 1]
    plan["evidence_status"]["replicates"] = count
    plan["operational_contract"]["datasets"] = count * 3
    plan["operational_contract"]["attempt_units"] = count * 12
    plan_path = tmp_path / "plan.json"
    plan_path.write_text(json.dumps(plan, indent=2) + "\n", encoding="utf-8")
    paths = {
        "runner_sha256": REPO_ROOT / "validation" / "informative_assignment_screening.py",
        "design_component_sha256": REPO_ROOT
        / "validation"
        / "informative_assignment_design.py",
        "resilient_pair_component_sha256": REPO_ROOT
        / "validation"
        / "facets_resilient_pair.py",
    }
    registration = {
        "plan_sha256": sha256_file(plan_path),
        **{key: sha256_file(path) for key, path in paths.items()},
        "tests_passed_before_registration": True,
    }
    registration_path = tmp_path / "registration.json"
    registration_path.write_text(
        json.dumps(registration, indent=2) + "\n", encoding="utf-8"
    )
    return plan_path, registration_path


def test_one_replicate_bundle_locks_three_designs_and_four_attempt_types(
    tmp_path: Path,
) -> None:
    plan, registration = _registered_files(tmp_path)
    bundle = generate_screening_bundle(
        replicates=1,
        replicate_start=221,
        base_seed=BASE_SEED,
        plan_path=plan,
        registration_path=registration,
    )
    manifest = bundle["manifest.csv"]
    ratings = bundle["generated_ratings.csv"]
    attempts = bundle["attempt_manifest.csv"]
    assert len(manifest) == 3
    assert manifest.set_index("Design")["ExpectedRows"].to_dict() == {
        "complete": 1920,
        "planned_connected": 960,
        "ability_severity_aligned_connected": 960,
    }
    assert manifest["ExpectedStructuralNullity"].eq(0).all()
    assert manifest["PersonRaterComponents"].eq(1).all()
    assert manifest["CompleteCriterionCategorySupport"].all()
    assert manifest["MinimumCriterionCategoryCount"].gt(0).all()
    assert manifest["CompleteResponseSHA256"].nunique() == 1
    assert manifest["SharedUniformSHA256"].nunique() == 1
    assert (
        manifest.loc[
            manifest["Design"].eq("ability_severity_aligned_connected"),
            "ThetaAssignedSeveritySpearman",
        ].iloc[0]
        > 0.9
    )
    assert ratings.groupby("RunId").size().sort_values().tolist() == [960, 960, 1920]
    assert "Theta" not in ratings.columns
    assert len(bundle["generated_facet_truth.csv"]) == 267
    assert len(bundle["generated_pcm_threshold_truth.csv"]) == 18
    assert len(attempts) == 12
    assert attempts["AttemptType"].nunique() == 4
    assert attempts["AttemptFingerprint"].nunique() == 12


def test_bundle_rejects_unregistered_range(tmp_path: Path) -> None:
    plan, registration = _registered_files(tmp_path)
    with pytest.raises(ValueError, match="does not match"):
        generate_screening_bundle(
            replicates=1,
            replicate_start=220,
            base_seed=BASE_SEED,
            plan_path=plan,
            registration_path=registration,
        )


def test_modulo_shards_are_disjoint_and_exhaustive() -> None:
    attempts = pd.DataFrame({"AttemptOrdinal": range(240), "AttemptId": range(240)})
    shards = [
        select_shard_attempts(attempts, shard_index=index, shard_count=8)
        for index in range(8)
    ]
    observed = [int(value) for frame in shards for value in frame["AttemptOrdinal"]]
    assert len(observed) == 240
    assert len(set(observed)) == 240
    assert sorted(observed) == list(range(240))


def test_screening_summary_has_32_descriptive_cells_and_no_confirmation() -> None:
    rows = []
    modes = ["PYTHON_JMLE", "MML_FIXED", "MML_FREE", "CMLE"]
    domains = ["Facet:Rater", "Facet:Task", "Facet:Criterion", "Threshold"]
    for replicate in range(221, 241):
        for mode in modes:
            for domain in domains:
                for metric in ("RMSE", "MAE"):
                    rows.append({
                        "Replicate": replicate,
                        "EstimatorMode": mode,
                        "RecoveryDomain": domain,
                        "Metric": metric,
                        "ContrastAlignedMinusPlanned": 0.04
                        + (replicate % 5) * 0.001,
                    })
    summary = screening_contrast_summary(pd.DataFrame(rows))
    assert len(summary) == 32
    assert summary["FinitePairs"].eq(20).all()
    assert summary["IntervalExcludesZero"].all()
    assert not summary["ConfirmatoryClaimAllowed"].any()
