from __future__ import annotations

import numpy as np

from validation import facets_pcm_known_truth_smoke as smoke


def test_pcm_probabilities_are_normalized_and_shift_toward_high_scores():
    thresholds = np.array([-1.1, 0.0, 1.1])
    low = smoke.pcm_probabilities(-1.5, thresholds)
    high = smoke.pcm_probabilities(1.5, thresholds)

    np.testing.assert_allclose(low.sum(), 1.0)
    np.testing.assert_allclose(high.sum(), 1.0)
    assert np.dot(high, np.arange(4)) > np.dot(low, np.arange(4))


def test_locked_threshold_vectors_are_sum_to_zero():
    for condition in smoke.THRESHOLD_CONDITIONS.values():
        for vector in condition.values():
            assert abs(float(np.sum(vector))) < 1e-12


def test_generated_conditions_share_latent_design_and_uniforms():
    bundle = smoke.generate_bundle()
    manifest = bundle["manifest.csv"]
    uniforms = bundle["paired_uniforms.csv"]

    for replicate, group in manifest.groupby("Replicate"):
        assert group["PairedLatentSHA256"].nunique() == 1
        selected = uniforms[uniforms["RunId"].isin(group["RunId"])]
        pivot = selected.pivot(
            index=["Person", "Rater", "Task", "Criterion"],
            columns="RunId",
            values="Uniform",
        )
        assert pivot.shape[1] == 2
        np.testing.assert_allclose(pivot.iloc[:, 0], pivot.iloc[:, 1])


def test_generated_bundle_has_locked_rows_and_full_support():
    bundle = smoke.generate_bundle()
    ratings = bundle["generated_ratings.csv"]

    assert len(bundle["manifest.csv"]) == 4
    assert ratings.groupby("RunId").size().eq(1920).all()
    assert ratings.groupby("RunId")["Score"].nunique().eq(4).all()


def test_larger_bundle_is_an_exact_parent_extension(tmp_path):
    parent = tmp_path / "parent" / "retained_input"
    child = tmp_path / "child" / "retained_input"
    smoke.write_bundle(smoke.generate_bundle(seeds=(731101, 731102)), parent)
    smoke.write_bundle(smoke.generate_bundle(seeds=(731101, 731102, 731103)), child)

    audit = smoke.audit_parent_extension(child, parent)

    assert len(audit) == 4
    assert audit["AllMatch"].all()


def test_expanded_pcm_design_is_full_column_rank():
    bundle = smoke.generate_bundle()
    first_run = bundle["manifest.csv"].iloc[0]["RunId"]
    ratings = bundle["generated_ratings.csv"]
    ratings = ratings[ratings["RunId"].eq(first_run)].drop(columns="RunId")

    audit = smoke.expanded_pcm_design_audit(ratings)

    assert audit["Nullity"] == 0
    assert audit["Rank"] == audit["Columns"]
