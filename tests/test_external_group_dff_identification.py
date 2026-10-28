"""Structural checks, not effect recovery, precision or MML qualification."""
from copy import deepcopy

import numpy as np
import pytest

from validation import external_group_dff_identification_check as check


@pytest.mark.parametrize("interaction,case,rank,nullity,cells,estimable_count", [
    ("reference_relative", "crossed", 20, 0, 8, 10),
    ("reference_relative", "nested_references", 16, 4, 4, 2),
    ("reference_relative", "missing_group_cell", 19, 1, 7, 4),
    ("symmetric", "crossed", 17, 0, 8, 10),
    ("symmetric", "nested_references", 16, 1, 4, 4),
    ("symmetric", "missing_group_cell", 17, 0, 7, 10),
])
def test_free_coordinate_svd_agrees_with_independent_rational_constraints(
    interaction, case, rank, nullity, cells, estimable_count,
):
    summary, contrasts, _ = check.audit_case(case, interaction)
    assert summary["Rank"] == summary["ExactRank"] == rank
    assert summary["FreeParameters"] == rank + nullity
    assert summary["Nullity"] == nullity
    assert summary["ObservedStudentCells"] == cells
    assert sum(row["Estimable"] for row in contrasts) == estimable_count
    assert all(row["Estimable"] == row["ExactEstimable"] for row in contrasts)
    if case == "nested_references" and interaction == "reference_relative":
        assert [row["Estimable"] for row in contrasts] == [False] * 8 + [True, True]
    if case == "missing_group_cell":
        # Full rank does not turn the absent A_R1/B cell into observed evidence.
        assert [row["AllComparedCellsObserved"] for row in contrasts] == [False, True, True, True, False, True, False, True, True, True]


def _nested_null_shift(columns, groups, interaction):
    shift = np.zeros(len(columns))
    for person, group in groups.items():
        if group == "A" or (interaction == "symmetric" and group == "B"):
            shift[columns["Person", person]] = -1
    for rater, q in zip(check.RATERS, (1, 1, -1, -1)):
        shift[columns["Gamma", (rater, "A")]] = q
        if interaction == "symmetric":
            shift[columns["Gamma", (rater, "B")]] = -q
    return shift


@pytest.mark.parametrize("full_lesson", [False, True])
@pytest.mark.parametrize("interaction,nullity,n_gamma", [("reference_relative", 4, 6), ("symmetric", 1, 3)])
def test_shared_references_link_baseline_but_leave_a_dff_null_direction(full_lesson, interaction, nullity, n_gamma):
    data, groups = check.teaching_design("nested_references", full_lesson=full_lesson)
    raw, constraints, expansion, columns = check.model_matrices(data, groups, interaction)
    design = raw @ expansion
    # This includes the free, centered RSM steps, not only additive eta.
    assert np.linalg.matrix_rank(design[:, :-n_gamma]) == design.shape[1] - n_gamma
    assert np.linalg.matrix_rank(design) == design.shape[1] - nullity
    shift = _nested_null_shift(columns, groups, interaction)
    np.testing.assert_array_equal(constraints @ shift, 0)
    np.testing.assert_array_equal(raw @ shift, 0)
    contrasts = check.planned_contrasts(columns)
    multiplier = 2 if interaction == "symmetric" else 1
    np.testing.assert_array_equal([row @ shift for _, row, _, _ in contrasts],
                                  multiplier * np.array([1, 1, -1, -1, 0, 0, 2, 2, 0, 0]))

    # A different parameter vector gives identical five-category probabilities.
    beta = expansion @ np.random.default_rng(42).normal(size=expansion.shape[1])
    def probabilities(parameters):
        adjacent = (raw @ parameters).reshape(-1, 4)
        logits = np.column_stack([np.zeros(len(adjacent)), adjacent.cumsum(axis=1)])
        mass = np.exp(logits - logits.max(axis=1, keepdims=True))
        return mass / mass.sum(axis=1, keepdims=True)
    np.testing.assert_allclose(probabilities(beta), probabilities(beta + 0.5 * shift), atol=1e-14)


@pytest.mark.parametrize("interaction,nullity", [("reference_relative", 3), ("symmetric", 0)])
def test_one_cross_rated_student_is_not_a_general_design_remedy(interaction, nullity):
    data, groups = check.teaching_design("one_cross_rating")
    raw, _, expansion, columns = check.model_matrices(data, groups, interaction)
    assert np.linalg.matrix_rank(raw @ expansion) == expansion.shape[1] - nullity
    shift = _nested_null_shift(columns, groups, interaction)
    assert np.count_nonzero(raw @ shift) == 3 * 4
    assert set((raw @ shift)[raw @ shift != 0]) == {-2}


def test_four_constructive_nested_ambiguities_leave_all_reference_values_unchanged():
    data, groups = check.teaching_design("nested_references")
    raw, constraints, _, columns = check.model_matrices(data, groups)
    shifts = np.zeros((4, len(columns)), dtype=int)
    for i, (group, q, person_shift) in enumerate([
        ("A", (1, 1, 0, -2), -1), ("A", (0, 0, 1, -1), 0),
        ("B", (0, -2, 1, 1), -1), ("B", (1, -1, 0, 0), 0),
    ]):
        for person in groups.index[groups.eq(group)]:
            shifts[i, columns["Person", person]] = person_shift
        for rater, value in zip(check.RATERS, q):
            shifts[i, columns["Gamma", (rater, group)]] = value
    assert len(check.exact_row_basis(shifts)) == 4
    np.testing.assert_array_equal(raw @ shifts.T, 0)
    np.testing.assert_array_equal(constraints @ shifts.T, 0)
    # Neither externally fixing reference persons nor rater severities removes these shifts.
    for person in groups.index[groups.isna()]:
        np.testing.assert_array_equal(shifts[:, columns["Person", person]], 0)
    for rater in check.RATERS:
        np.testing.assert_array_equal(shifts[:, columns["Rater", rater]], 0)


def test_rational_check_detects_dependence_without_svd_tolerances():
    basis = check.exact_row_basis(np.array([[2, 4, 1], [4, 8, 2], [0, 1, 3]]))
    assert len(basis) == 2
    assert check.exact_estimable([2, 5, 4], basis)
    assert not check.exact_estimable([1, 0, 0], basis)


@pytest.fixture(scope="module")
def audited_designs():
    return {case: check.audit_case(case)[:2] for case in check.CASES}


@pytest.mark.parametrize("case,state,count", [
    ("crossed", "model_identified", 8),
    ("nested_references", "target_unidentified", 0),
    ("missing_group_cell", "partially_identified", 2),
])
def test_three_output_states_use_only_declared_between_group_targets(audited_designs, case, state, count):
    summary, rows = audited_designs[case]
    before = deepcopy((summary, rows))
    plan = check.build_output_plan(summary, rows)
    assert (plan["State"], plan["EstimableTargetCount"], plan["TargetCount"]) == (state, count, 8)
    assert not any(plan[k] for k in ("PointEstimatesReleased", "IntervalsReleased", "PValuesReleased", "GroupCurvesReleased"))
    assert plan["Stage"] == "design_only"
    assert all("Within" not in row["Contrast"] for row in plan["Comparisons"])
    assert all(row["OutputScope"] == "Design review only; no effect estimate" for row in plan["Comparisons"])
    assert all(not {"Estimate", "SE", "LowerCI", "UpperCI", "p_value"}.intersection(row) for row in plan["Comparisons"])
    assert (summary, rows) == before


def test_identified_targets_do_not_turn_a_deficient_model_into_a_full_rank_fit(audited_designs):
    summary, rows = audited_designs["missing_group_cell"]
    targets = [r["Contrast"] for r in rows if r["ComparisonScope"] == "between_groups" and r["Estimable"]]
    plan = check.build_output_plan(summary, rows, targets=targets)
    assert (plan["EstimableTargetCount"], plan["TargetCount"]) == (2, 2)
    assert plan["State"] == "partially_identified" and plan["ModelNullity"] == 1
    # The same design cannot answer a request for individual rater A/B effects.
    individual = [r["Contrast"] for r in rows[:4]]
    blocked = check.build_output_plan(summary, rows, targets=individual)
    assert blocked["State"] == "target_unidentified"
    assert blocked["EstimableTargetCount"] == 0


def test_added_symmetry_is_recorded_and_does_not_release_unobserved_effects():
    summary, rows, _ = check.audit_case("missing_group_cell", "symmetric")
    plan = check.build_output_plan(summary, rows)
    assert plan["Model"] == "symmetric" and plan["State"] == "model_identified"
    assert not plan["Comparisons"][0]["AllComparedCellsObserved"]
    assert not plan["PointEstimatesReleased"] and not plan["GroupCurvesReleased"]


@pytest.mark.parametrize("field,value", [
    ("Estimable", None), ("Estimable", float("nan")), ("Estimable", "False"),
    ("ExactEstimable", False), ("AllComparedCellsObserved", None),
    ("Case", "nested_references"), ("Model", "symmetric"), ("ComparisonScope", None),
])
def test_missing_or_inconsistent_contrast_evidence_cannot_certify_a_design(audited_designs, field, value):
    summary, rows = deepcopy(audited_designs["crossed"])
    rows[0][field] = value
    with pytest.raises(ValueError, match="evidence"):
        check.build_output_plan(summary, rows)


@pytest.mark.parametrize("field,value", [("Nullity", None), ("Rank", 21), ("ExactRank", 19), ("FreeParameters", 0)])
def test_invalid_rank_audit_is_not_an_unidentified_or_identified_result(audited_designs, field, value):
    summary, rows = deepcopy(audited_designs["crossed"])
    summary[field] = value
    with pytest.raises(ValueError, match="rank audit"):
        check.build_output_plan(summary, rows)


@pytest.mark.parametrize("targets", [[], ["unknown"], ["Within A: A_R1 minus A_R2"], ["A_R1: A minus B"] * 2])
def test_no_target_or_a_within_school_target_cannot_be_counted_as_school_dff(audited_designs, targets):
    summary, rows = audited_designs["nested_references"]
    with pytest.raises(ValueError, match="between-group"):
        check.build_output_plan(summary, rows, targets=targets)


def test_observed_counts_cannot_substitute_for_a_design_audit():
    with pytest.raises(ValueError, match="declared design fixtures"):
        check.build_output_plan({"levels_in_all_groups": 4, "zero_cells": 0}, [])


@pytest.mark.parametrize("field,value", [("Estimator", "MML"), ("RatingModel", "PCM"), ("RatingModel", "GPCM")])
def test_other_estimators_and_models_do_not_inherit_this_output_specification(audited_designs, field, value):
    summary, rows = deepcopy(audited_designs["crossed"])
    summary[field] = value
    with pytest.raises(ValueError, match="declared design fixtures"):
        check.build_output_plan(summary, rows)
