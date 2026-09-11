from __future__ import annotations

import pandas as pd

from validation import facets_pcm_boundary_pilot as pilot
from validation.facets_pcm_known_truth_smoke import (
    THRESHOLD_CONDITIONS,
    apply_threshold_condition,
    generate_latent_replicate,
)


def _complete() -> pd.DataFrame:
    latent, _ = generate_latent_replicate(pilot.SEEDS[0])
    return apply_threshold_condition(latent, THRESHOLD_CONDITIONS["heterogeneous"])


def test_observation_designs_have_locked_row_counts_and_components():
    complete = _complete()
    expected = {
        "complete": (1920, 1),
        "planned_connected": (960, 1),
        "disconnected_negative_control": (960, 2),
    }
    for design, (rows, components) in expected.items():
        observed = pilot.apply_observation_design(complete, design)
        assert len(observed) == rows
        assert len(pilot.person_rater_components(observed)) == components


def test_constrained_design_audit_separates_rank_full_and_disconnected():
    complete = _complete()
    for model in pilot.FIT_MODELS:
        full = pilot.constrained_adjacent_design_audit(
            pilot.apply_observation_design(complete, "complete"), model
        )
        connected = pilot.constrained_adjacent_design_audit(
            pilot.apply_observation_design(complete, "planned_connected"), model
        )
        disconnected = pilot.constrained_adjacent_design_audit(
            pilot.apply_observation_design(complete, "disconnected_negative_control"), model
        )
        assert full["Nullity"] == 0
        assert connected["Nullity"] == 0
        assert disconnected["Nullity"] == 1


def test_boundary_bundle_contains_paired_conditions_and_expected_rows():
    bundle = pilot.generate_boundary_bundle()
    manifest = bundle["manifest.csv"]
    ratings = bundle["generated_ratings.csv"]

    assert len(manifest) == 12
    assert manifest.groupby(["Design", "ThresholdCondition"]).size().eq(2).all()
    observed_rows = ratings.groupby("RunId").size()
    expected_rows = manifest.set_index("RunId")["ExpectedRows"]
    assert observed_rows.eq(expected_rows).all()
    assert ratings.groupby("RunId")["Score"].nunique().eq(4).all()
    for (_, replicate), group in manifest.groupby(["Design", "Replicate"]):
        assert group["PairedLatentSHA256"].nunique() == 1


def test_rsm_spec_has_one_common_scale_and_no_bias_model(tmp_path):
    bundle = pilot.generate_boundary_bundle()
    manifest_row = bundle["manifest.csv"].iloc[0]
    run_id = manifest_row["RunId"]
    ratings = bundle["generated_ratings.csv"]
    ratings = ratings[ratings["RunId"].eq(run_id)].drop(columns="RunId")

    spec, _ = pilot.build_rsm_spec(
        manifest_row, ratings, score_base=tmp_path / "scores.txt"
    )

    assert "Models=\n?,?,?,?,R3\n*" in spec
    assert "#" not in spec
    assert "?B" not in spec


def test_rsm_count_audit_matches_independent_counts():
    ratings = pd.DataFrame({"Score": [0, 0, 1, 2, 3, 3]})
    categories = pd.DataFrame({
        "TableNumber": ["8.1"] * 4,
        "Category": [0, 1, 2, 3],
        "TotalCount": [2, 1, 1, 2],
    })

    audit = pilot.rsm_count_audit(categories, ratings)

    assert audit["CountMatches"].all()


def test_direction_check_is_paired_within_design_and_replicate():
    model_selection = pd.DataFrame({
        "Design": ["complete", "complete"],
        "Replicate": [1, 1],
        "ThresholdCondition": ["shared", "heterogeneous"],
        "ModelComparisonEligible": [True, True],
        "PCMMinusRSMLogLikPerObs": [0.001, 0.08],
    })

    checks = pilot.build_direction_checks(model_selection)

    assert len(checks) == 1
    assert bool(checks.iloc[0]["DirectionPass"])
