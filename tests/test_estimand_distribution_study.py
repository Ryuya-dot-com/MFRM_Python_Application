import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from validation.estimand_distribution_study import (
    AMENDMENT_PATH,
    ATTEMPT_TYPES,
    DESIGNS,
    DISTRIBUTIONS,
    PERSON_SD,
    PLAN_PATH,
    _validate_completion,
    aggregate_study,
    enrich_registered_metadata,
    generate_person_shapes,
    generate_study_bundle,
    prepare_study,
    select_shard_attempts,
    sha256_file,
    validate_study_identity,
)


def test_person_shapes_have_locked_first_two_moments_and_distinct_shapes():
    shapes = generate_person_shapes(812345)
    assert tuple(shapes) == DISTRIBUTIONS
    for values in shapes.values():
        assert len(values) == 80
        assert np.mean(values) == pytest.approx(0.0, abs=1e-12)
        assert np.std(values, ddof=0) == pytest.approx(PERSON_SD, abs=1e-12)
    assert not np.allclose(shapes["normal"], shapes["right_skew"])
    assert not np.allclose(shapes["normal"], shapes["heavy_tail_t3"])


def test_one_replicate_bundle_is_paired_rank_full_and_attempt_complete():
    bundle = generate_study_bundle(replicates=1)
    manifest = bundle["manifest.csv"]
    attempts = bundle["attempt_manifest.csv"]
    ratings = bundle["generated_ratings.csv"]
    assert len(manifest) == len(DISTRIBUTIONS) * len(DESIGNS) == 8
    assert len(attempts) == len(manifest) * len(ATTEMPT_TYPES) == 32
    assert attempts["AttemptOrdinal"].tolist() == list(range(32))
    assert attempts["AttemptId"].is_unique
    assert attempts["AttemptFingerprint"].is_unique
    assert manifest["ExpectedStructuralNullity"].eq(0).all()
    assert manifest["PersonRaterComponents"].eq(1).all()
    assert manifest["SharedUniformSHA256"].nunique() == 1
    assert manifest.groupby("PersonDistribution")["CompleteResponseSHA256"].nunique().eq(1).all()
    observed_rows = ratings.groupby("RunId").size()
    expected_rows = manifest.set_index("RunId")["ExpectedRows"]
    pd.testing.assert_series_equal(
        observed_rows.sort_index(), expected_rows.sort_index(), check_names=False
    )
    assert set(manifest.loc[manifest["Design"].eq("complete"), "ExpectedRows"]) == {1920}
    assert set(manifest.loc[manifest["Design"].eq("planned_connected"), "ExpectedRows"]) == {960}


def test_shards_are_disjoint_exhaustive_and_stable():
    attempts = generate_study_bundle(replicates=1)["attempt_manifest.csv"]
    shards = [
        select_shard_attempts(attempts, shard_index=index, shard_count=3)
        for index in range(3)
    ]
    observed = pd.concat(shards)["AttemptOrdinal"]
    assert observed.is_unique
    assert sorted(observed) == list(range(len(attempts)))
    replay = select_shard_attempts(attempts, shard_index=1, shard_count=3)
    pd.testing.assert_frame_equal(shards[1], replay)
    with pytest.raises(ValueError, match="0 <= shard-index"):
        select_shard_attempts(attempts, shard_index=3, shard_count=3)


def test_prepare_hash_validation_and_incomplete_aggregate_fail_closed(tmp_path):
    study_dir = tmp_path / "preflight"
    prepare_study(study_dir, replicates=1, base_seed=2026081107)
    identity = validate_study_identity(study_dir)
    assert identity["phase"] == "preflight"
    assert identity["datasets"] == 8
    assert identity["attempts"] == 32
    metrics = aggregate_study(study_dir)
    assert metrics["completed_attempts"] == 0
    assert metrics["qualification_pass"] is False
    assert metrics["gates"]["attempts_complete"] is False


def test_confirmatory_prepare_uses_fresh_replicate_range_and_custom_plan(tmp_path):
    custom_plan = tmp_path / "confirmatory_plan.json"
    custom_plan.write_text(json.dumps({
        "phase": "confirmatory",
        "evidence_status": {
            "screening_replicates": [1, 20],
            "confirmatory_replicates": [21, 21],
            "confirmatory_replicate_count": 1,
            "fresh_data_required": True,
        },
    }) + "\n", encoding="utf-8")
    study_dir = tmp_path / "confirmatory"
    prepare_study(
        study_dir,
        replicates=1,
        base_seed=2026081107,
        replicate_start=21,
        phase="confirmatory",
        plan_path=custom_plan,
    )
    identity = validate_study_identity(study_dir)
    assert identity["phase"] == "confirmatory"
    assert identity["replicate_start"] == 21
    assert identity["replicate_end"] == 21
    assert identity["plan_sha256"] == sha256_file(custom_plan)
    assert "Confirmatory only" in identity["claim_limit"]
    manifest = pd.read_csv(study_dir / "retained_input" / "manifest.csv")
    assert manifest["Replicate"].unique().tolist() == [21]
    assert manifest["RunId"].str.endswith("rep-00021").all()
    screening = generate_study_bundle(replicates=1)["manifest.csv"]
    assert set(manifest["Seed"]).isdisjoint(set(screening["Seed"]))
    with pytest.raises(ValueError, match="does not match confirmatory plan"):
        prepare_study(
            tmp_path / "wrong_range",
            replicates=1,
            base_seed=2026081107,
            replicate_start=22,
            phase="confirmatory",
            plan_path=custom_plan,
        )


def test_custom_plan_completion_validation_is_explicit(tmp_path):
    custom_plan = tmp_path / "confirmatory_plan.json"
    custom_plan.write_text('{"phase": "confirmatory"}\n', encoding="utf-8")
    attempt = generate_study_bundle(
        replicates=1,
        replicate_start=21,
        plan_path=custom_plan,
    )["attempt_manifest.csv"].iloc[0]
    completion = {
        "attempt_fingerprint": attempt["AttemptFingerprint"],
        "plan_sha256": sha256_file(custom_plan),
        "infrastructure_amendment_sha256": sha256_file(AMENDMENT_PATH),
        "script_sha256": sha256_file(Path(
            __import__(
                "validation.estimand_distribution_study", fromlist=["__file__"]
            ).__file__
        )),
    }
    _validate_completion(
        completion,
        attempt,
        expected_plan_sha256=sha256_file(custom_plan),
    )
    with pytest.raises(ValueError, match="plan mismatch"):
        _validate_completion(completion, attempt)


def test_completion_validation_rejects_plan_or_attempt_mismatch():
    attempt = generate_study_bundle(replicates=1)["attempt_manifest.csv"].iloc[0]
    completion = {
        "attempt_fingerprint": attempt["AttemptFingerprint"],
        "plan_sha256": sha256_file(PLAN_PATH),
        "infrastructure_amendment_sha256": sha256_file(AMENDMENT_PATH),
        "script_sha256": sha256_file(Path(
            __import__(
                "validation.estimand_distribution_study", fromlist=["__file__"]
            ).__file__
        )),
    }
    _validate_completion(completion, attempt)
    bad = json.loads(json.dumps(completion))
    bad["attempt_fingerprint"] = "0" * 64
    with pytest.raises(ValueError, match="fingerprint mismatch"):
        _validate_completion(bad, attempt)


def test_registered_metadata_repairs_missing_values_and_rejects_conflicts():
    manifest = pd.DataFrame([
        {
            "RunId": "run-1",
            "ConditionId": "complete__normal",
            "Design": "complete",
            "PersonDistribution": "normal",
            "Replicate": 1,
        }
    ])
    artifact = pd.DataFrame([
        {"RunId": "run-1", "Design": "complete", "PersonDistribution": np.nan, "Estimate": 0.2}
    ])
    repaired = enrich_registered_metadata(artifact, manifest, label="test artifact")
    assert repaired.loc[0, "PersonDistribution"] == "normal"
    assert repaired.loc[0, "ConditionId"] == "complete__normal"
    assert repaired.loc[0, "Replicate"] == 1
    conflicting = artifact.assign(Design="planned_connected")
    with pytest.raises(ValueError, match="conflicts with registered Design"):
        enrich_registered_metadata(conflicting, manifest, label="test artifact")
