from __future__ import annotations

import numpy as np
import pandas.testing as pdt

from mfrm_app.assignment_generator import (
    MML_RANK_PRESERVING_MODE,
    SOURCE_FITTED_MODE,
    build_assignment_person_generation_plan,
    draw_assignment_person_coordinates,
)


PERSON_SCORES = {f"p{index}": float(index) / 3.0 for index in range(8)}


def test_source_fitted_generator_is_available_for_jmle_and_is_not_known_truth():
    plan = build_assignment_person_generation_plan(
        method="JMLE",
        person_scores=PERSON_SCORES,
        mode=SOURCE_FITTED_MODE,
    )
    draw = draw_assignment_person_coordinates(
        plan,
        rng=np.random.default_rng(1),
        replicate=1,
    )

    assert plan["available"] is True
    assert plan["contract"].iloc[0]["SourceFittedPersonIsKnownTruth"] is False or not bool(
        plan["contract"].iloc[0]["SourceFittedPersonIsKnownTruth"]
    )
    assert draw["coordinates"] == PERSON_SCORES
    assert not draw["summary"].iloc[0]["RandomPopulationDraw"]


def test_rank_preserving_mml_draw_uses_population_sd_and_common_strict_order():
    plan = build_assignment_person_generation_plan(
        method="MML",
        person_scores=PERSON_SCORES,
        mode=MML_RANK_PRESERVING_MODE,
        population_sd=1.25,
    )
    first = draw_assignment_person_coordinates(
        plan,
        rng=np.random.default_rng(260811),
        replicate=1,
    )
    second = draw_assignment_person_coordinates(
        plan,
        rng=np.random.default_rng(260811),
        replicate=1,
    )

    assert plan["available"] is True
    assert plan["gates"]["Passed"].all()
    assert plan["contract"].iloc[0]["GeneratedPersonTruthKnownWithinReplicate"]
    assert plan["contract"].iloc[0]["UnconditionalNewSample"] is False or not bool(
        plan["contract"].iloc[0]["UnconditionalNewSample"]
    )
    assert first["summary"].iloc[0]["RankPreservedExactly"]
    assert first["summary"].iloc[0]["GeneratorPopulationSD"] == 1.25
    assert first["draws"]["GeneratedCoordinate"].is_monotonic_increasing
    pdt.assert_frame_equal(first["draws"], second["draws"])


def test_population_mode_blocks_jmle_ties_latent_regression_and_invalid_sd():
    jmle = build_assignment_person_generation_plan(
        method="JMLE",
        person_scores=PERSON_SCORES,
        mode=MML_RANK_PRESERVING_MODE,
        population_sd=1.0,
    )
    tied = build_assignment_person_generation_plan(
        method="MML",
        person_scores={"p0": 0.0, "p1": 0.0, "p2": 1.0},
        mode=MML_RANK_PRESERVING_MODE,
        population_sd=1.0,
    )
    latent = build_assignment_person_generation_plan(
        method="MML",
        person_scores=PERSON_SCORES,
        mode=MML_RANK_PRESERVING_MODE,
        population_sd=1.0,
        population_model_enabled=True,
    )
    invalid_sd = build_assignment_person_generation_plan(
        method="MML",
        person_scores=PERSON_SCORES,
        mode=MML_RANK_PRESERVING_MODE,
        population_sd=0.0,
    )

    assert jmle["available"] is False
    assert "Generator-estimator compatibility" in jmle["reason"]
    assert tied["available"] is False
    assert "Strict source Person ordering" in tied["reason"]
    assert latent["available"] is False
    assert "No latent-regression population model" in latent["reason"]
    assert invalid_sd["available"] is False
    assert "Finite positive MML population SD" in invalid_sd["reason"]


def test_population_draws_change_across_replicates_but_reproduce_by_seed():
    plan = build_assignment_person_generation_plan(
        method="MML",
        person_scores=PERSON_SCORES,
        mode=MML_RANK_PRESERVING_MODE,
        population_sd=0.8,
    )
    rng = np.random.default_rng(99)
    first = draw_assignment_person_coordinates(plan, rng=rng, replicate=1)
    second = draw_assignment_person_coordinates(plan, rng=rng, replicate=2)
    replay_rng = np.random.default_rng(99)
    replay = draw_assignment_person_coordinates(plan, rng=replay_rng, replicate=1)

    assert first["coordinates"] != second["coordinates"]
    assert first["coordinates"] == replay["coordinates"]
