from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from validation.h6_design_replication import (
    BASE_SEED,
    REPO_ROOT,
    _next_short_artifact_root,
    generate_replication_bundle,
    select_shard_attempts,
    summarize_primary,
)
from validation.operating_characteristics_facets import sha256_file


def _registered_files(tmp_path: Path, *, start: int = 121, count: int = 1) -> tuple[Path, Path]:
    plan = json.loads(
        (REPO_ROOT / "validation" / "h6_design_replication_plan_20260811.json").read_text(
            encoding="utf-8"
        )
    )
    plan["evidence_separation"]["replication_replicates"] = [start, start + count - 1]
    plan["evidence_separation"]["replication_count"] = count
    plan_path = tmp_path / "plan.json"
    plan_path.write_text(json.dumps(plan, indent=2) + "\n", encoding="utf-8")
    runner = REPO_ROOT / "validation" / "h6_design_replication.py"
    component = REPO_ROOT / "validation" / "facets_resilient_pair.py"
    registration = {
        "plan_sha256": sha256_file(plan_path),
        "runner_sha256": sha256_file(runner),
        "resilient_pair_component_sha256": sha256_file(component),
        "tests_passed_before_registration": True,
    }
    registration_path = tmp_path / "registration.json"
    registration_path.write_text(json.dumps(registration, indent=2) + "\n", encoding="utf-8")
    return plan_path, registration_path


def test_one_replicate_bundle_preserves_registered_dgm_contract(tmp_path: Path) -> None:
    plan_path, registration_path = _registered_files(tmp_path)
    bundle = generate_replication_bundle(
        replicates=1,
        replicate_start=121,
        base_seed=BASE_SEED,
        plan_path=plan_path,
        registration_path=registration_path,
    )

    manifest = bundle["manifest.csv"]
    ratings = bundle["generated_ratings.csv"]
    assert len(manifest) == 2
    assert set(manifest["PersonDistribution"]) == {"normal"}
    assert set(manifest["Design"]) == {"complete", "planned_connected"}
    assert manifest.set_index("Design")["ExpectedRows"].to_dict() == {
        "complete": 1920,
        "planned_connected": 960,
    }
    assert manifest["ExpectedStructuralNullity"].eq(0).all()
    assert manifest["PersonRaterComponents"].eq(1).all()
    assert manifest["PersonMeanRealized"].abs().max() <= 1e-12
    assert (manifest["PersonSDRealized"] - 0.8).abs().max() <= 1e-12
    assert ratings.groupby("RunId").size().sort_values().tolist() == [960, 1920]
    assert len(bundle["generated_facet_truth.csv"]) == 178
    assert len(bundle["generated_pcm_threshold_truth.csv"]) == 12
    assert len(bundle["attempt_manifest.csv"]) == 2
    assert bundle["attempt_manifest.csv"]["Replicate"].eq(121).all()


def test_bundle_rejects_range_not_matching_registration(tmp_path: Path) -> None:
    plan_path, registration_path = _registered_files(tmp_path)
    with pytest.raises(ValueError, match="does not match"):
        generate_replication_bundle(
            replicates=1,
            replicate_start=120,
            base_seed=BASE_SEED,
            plan_path=plan_path,
            registration_path=registration_path,
        )


def test_registration_detects_runner_hash_tampering(tmp_path: Path) -> None:
    plan_path, registration_path = _registered_files(tmp_path)
    value = json.loads(registration_path.read_text(encoding="utf-8"))
    value["runner_sha256"] = "0" * 64
    registration_path.write_text(json.dumps(value), encoding="utf-8")
    with pytest.raises(ValueError, match="runner_sha256"):
        generate_replication_bundle(
            replicates=1,
            replicate_start=121,
            base_seed=BASE_SEED,
            plan_path=plan_path,
            registration_path=registration_path,
        )


def test_modulo_shards_are_disjoint_and_exhaustive() -> None:
    attempts = pd.DataFrame({"AttemptOrdinal": range(200), "AttemptId": range(200)})
    shards = [
        select_shard_attempts(attempts, shard_index=index, shard_count=8)
        for index in range(8)
    ]
    observed = [int(value) for frame in shards for value in frame["AttemptOrdinal"]]
    assert len(observed) == 200
    assert len(set(observed)) == 200
    assert sorted(observed) == list(range(200))


def test_row_count_contract_is_invariant_to_groupby_index_order() -> None:
    expected = pd.Series(
        [1920, 960], index=["complete__normal", "planned_connected__normal"]
    )
    observed = pd.Series(
        [960, 1920], index=["planned_connected__normal", "complete__normal"]
    )
    assert observed.reindex(expected.index).astype(int).equals(expected)


def test_short_artifact_root_is_retained_and_generation_safe(tmp_path: Path) -> None:
    attempt = pd.Series({"AttemptOrdinal": 199})
    first = _next_short_artifact_root(tmp_path, attempt)
    assert first.relative_to(tmp_path).as_posix() == "w/00c7/g00"
    first.mkdir()
    second = _next_short_artifact_root(tmp_path, attempt)
    assert second.relative_to(tmp_path).as_posix() == "w/00c7/g01"


def test_absolute_worst_case_facets_path_stays_below_registered_gate() -> None:
    study_dir = (REPO_ROOT / "validation" / "h6_design_replication100_v4_20260811").resolve()
    longest = study_dir / (
        "w/00c7/g00/facets_attempts/pair_try_03/facets_runs/"
        "planned_connected__normal_rep-00220__pcm/report_u6.txt"
    )
    assert longest.is_absolute()
    assert len(str(longest)) < 240


def test_primary_summary_uses_exactly_six_native_thresholds_per_run() -> None:
    rows = []
    for replicate in range(121, 221):
        complete_error = 0.10 + (replicate - 170) * 0.0002
        planned_error = complete_error + 0.04 + (replicate % 7) * 0.0003
        for design, error in (
            ("complete", complete_error),
            ("planned_connected", planned_error),
        ):
            for criterion in ("C01", "C02"):
                for category in (1, 2, 3):
                    rows.append({
                        "RunId": f"{design}__normal::rep-{replicate:05d}",
                        "Replicate": replicate,
                        "Design": design,
                        "EstimatorMode": "PYTHON_JMLE",
                        "IncludedInStudy": True,
                        "StepFacetLevel": criterion,
                        "Category": category,
                        "TruthError": error,
                    })
    result, paired, elements = summarize_primary(pd.DataFrame(rows))
    primary = result.iloc[0]
    assert int(primary["FinitePairs"]) == 100
    assert bool(primary["MissingRuleSatisfied"])
    assert bool(primary["ReplicationDecision"])
    assert float(primary["MeanContrast"]) > 0.04
    assert len(paired) == 100
    assert len(elements) == 6
    assert elements["N"].eq(100).all()


def test_primary_summary_is_inconclusive_at_99_pairs() -> None:
    rows = []
    for replicate in range(121, 220):
        for design, error in (("complete", 0.1), ("planned_connected", 0.2)):
            for criterion in ("C01", "C02"):
                for category in (1, 2, 3):
                    rows.append({
                        "RunId": f"{design}-{replicate}",
                        "Replicate": replicate,
                        "Design": design,
                        "EstimatorMode": "PYTHON_JMLE",
                        "IncludedInStudy": True,
                        "StepFacetLevel": criterion,
                        "Category": category,
                        "TruthError": error,
                    })
    result, _, _ = summarize_primary(pd.DataFrame(rows))
    primary = result.iloc[0]
    assert int(primary["FinitePairs"]) == 99
    assert not bool(primary["MissingRuleSatisfied"])
    assert not bool(primary["ReplicationDecision"])
