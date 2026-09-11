"""Contracts for no-human-data confirmatory sample-size sensitivity."""

from __future__ import annotations

import math

from mfrm_app.cmle_one_click_comprehension import HUMAN_STUDY_STATUS
from mfrm_app.cmle_one_click_confirmatory_gate import maximum_passing_errors
from mfrm_app.cmle_one_click_confirmatory_planning import (
    binomial_cdf_direct,
    binomial_cdf_logspace,
    build_attrition_assurance_table,
    build_blocked_mechanism_exposure_table,
    build_cell_probability_surface,
    build_cluster_sensitivity_table,
    build_minimum_n_sensitivity,
    cell_gate_pass_probability,
    cluster_design_effect,
    confirmatory_planning_human_gate_status,
    joint_pass_projections,
    maximum_passing_errors_fast,
    minimum_planned_slots_for_valid_target,
)


def test_logspace_binomial_matches_independent_direct_sums_and_edges() -> None:
    for n in (5, 10, 25, 50):
        for probability in (0.01, 0.05, 0.10, 0.50):
            for maximum in (0, 1, min(3, n), n):
                observed = binomial_cdf_logspace(maximum, n, probability)
                expected = binomial_cdf_direct(maximum, n, probability)
                assert math.isclose(observed, expected, rel_tol=2e-13, abs_tol=2e-15)
    assert binomial_cdf_logspace(-1, 10, 0.1) == 0.0
    assert binomial_cdf_logspace(0, 10, 0.0) == 1.0
    assert binomial_cdf_logspace(9, 10, 1.0) == 0.0
    assert binomial_cdf_logspace(10, 10, 1.0) == 1.0


def test_fast_passing_error_search_matches_registered_exhaustive_function() -> None:
    for n in (24, 25, 30, 50, 75, 100, 200, 500):
        assert maximum_passing_errors_fast(n) == maximum_passing_errors(n)
        assert maximum_passing_errors_fast(
            n, confidence=1.0 - 0.05 / 12.0
        ) == maximum_passing_errors(n, confidence=1.0 - 0.05 / 12.0)


def test_cell_surface_retains_all_assumptions_without_selecting_n() -> None:
    surface = build_cell_probability_surface()
    assert len(surface) == 6 * 11 * 2
    assert surface["TrueDangerousErrorProbability"].nunique() == 6
    assert surface["ValidEligibleN"].nunique() == 11
    assert surface["ConfidenceScheme"].nunique() == 2
    assert surface["CellPassProbabilityRaw"].between(0, 1).all()
    assert not surface["SampleSizeSelected"].any()
    zero = surface.loc[
        surface["TrueDangerousErrorProbability"].eq(0)
        & surface["ConfidenceScheme"].eq("primary_per_cell")
    ].set_index("ValidEligibleN")
    assert zero.loc[24, "CellPassProbabilityRaw"] == 0.0
    assert zero.loc[25, "CellPassProbabilityRaw"] == 1.0


def test_minimum_searches_meet_target_and_predecessor_fails() -> None:
    cell, joint = build_minimum_n_sensitivity(
        true_error_probabilities=(0.01, 0.05, 0.10),
        probability_targets=(0.80, 0.90),
        search_max_n=1500,
    )
    for _, row in cell.loc[cell["FirstCrossingAvailable"]].iterrows():
        assert row["ProbabilityAtFirstCrossing"] >= row["CellPassProbabilityTarget"]
        assert row["ProbabilityAtFirstCrossingPredecessor"] < row["CellPassProbabilityTarget"]
    for _, row in cell.loc[cell["SustainedThroughSearchAvailable"]].iterrows():
        assert row["MinimumProbabilityFromSustainedThroughSearchMax"] >= row[
            "CellPassProbabilityTarget"
        ]
    for _, row in joint.loc[joint["FirstCrossingAvailable"]].iterrows():
        assert row["ProbabilityAtFirstCrossing"] >= row["JointUnionBoundTarget"]
        assert row["ProbabilityAtFirstCrossingPredecessor"] < row["JointUnionBoundTarget"]
        assert row["ProbabilityAtFirstCrossing"] <= row[
            "IndependenceProjectionAtFirstCrossing"
        ] + 1e-15
    for _, row in joint.loc[joint["SustainedThroughSearchAvailable"]].iterrows():
        assert row["MinimumProbabilityFromSustainedThroughSearchMax"] >= row[
            "JointUnionBoundTarget"
        ]
    assert cell["DownwardAdjacentStepCount"].gt(0).any()
    assert not cell["SampleSizeSelected"].any()
    assert not joint["SampleSizeSelected"].any()


def test_joint_projection_separates_independence_and_union_bound() -> None:
    for probability in (0.0, 0.5, 0.9, 0.99, 1.0):
        result = joint_pass_projections(probability)
        assert 0 <= result["UnionBoundLowerRaw"] <= 1
        assert 0 <= result["IndependenceProjectionRaw"] <= 1
        assert result["UnionBoundLowerRaw"] <= result["IndependenceProjectionRaw"] + 1e-15


def test_attrition_assurance_is_minimal_and_not_a_selection() -> None:
    result = minimum_planned_slots_for_valid_target(50, 0.90, 0.95)
    assert result["AssuranceAtMinimumRaw"] >= 0.95
    assert result["AssuranceAtPredecessorRaw"] < 0.95
    table = build_attrition_assurance_table(
        target_valid_n=(25, 50),
        retention_probabilities=(0.9, 0.8),
        assurance_targets=(0.90, 0.95),
    )
    assert len(table) == 8
    assert not table["PlannedRecruitmentNSelected"].any()
    assert (table["MinimumPlannedSlots"] >= table["TargetValidN"]).all()
    assert (table["AssuranceAtMinimumRaw"] >= table["AssuranceTarget"]).all()
    assert (table["AssuranceAtPredecessorRaw"] < table["AssuranceTarget"]).all()


def test_cluster_design_effect_is_explicitly_heuristic() -> None:
    assert cluster_design_effect(1, 0.10) == 1.0
    assert cluster_design_effect(20, 0.0) == 1.0
    assert cluster_design_effect(10, 0.05) == 1.45
    table = build_cluster_sensitivity_table(
        effective_targets=(25, 100),
        average_cluster_sizes=(1, 10),
        intraclass_correlations=(0.0, 0.05),
    )
    assert len(table) == 8
    assert table["Interpretation"].str.startswith("heuristic_").all()
    assert not table["SampleSizeSelected"].any()
    assert (table["HeuristicNominalValidN"] >= table["EffectiveTargetN"]).all()


def test_blocked_mechanism_rotation_and_human_gate_remain_unselected() -> None:
    exposure = build_blocked_mechanism_exposure_table()
    assert len(exposure) == 8 * 2 * 3
    spread = exposure.groupby(
        ["CandidateSlotsPerLanguage", "Language"]
    )["MechanismExposureN"].agg(lambda values: int(values.max() - values.min()))
    assert spread.le(1).all()
    assert not exposure["MechanismSpecificPrimaryGateActivated"].any()
    status = confirmatory_planning_human_gate_status().iloc[0]
    assert status["HumanStudyStatus"] == HUMAN_STUDY_STATUS
    assert status["HumanParticipants"] == 0
    assert not bool(status["MinimumValidPerCellRegistered"])
    assert not bool(status["PlannedRecruitmentNSelected"])
    assert not bool(status["PilotResultAvailable"])
    assert not bool(status["ConfirmatoryResultAvailable"])
    assert not bool(status["PublicSurfaceEnabled"])
