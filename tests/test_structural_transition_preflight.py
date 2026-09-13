"""Identifiability boundary for the non-zero K-1 repair overlay."""

from __future__ import annotations

from dataclasses import replace
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

from mfrm_app.simulation import (
    BASE_JMLE_BLOCKED,
    BASE_MML_PRIOR_CONDITIONAL,
    CLAIM_NOT_ESTIMABLE_COMMON_SCALE,
    CLAIM_PLANNED_CONNECTIVITY_PLUS_PRIOR_PENDING_FIT,
    CLAIM_PRIOR_CONDITIONAL_EXPLORATORY,
    COMPARISON_ESTIMABILITY_TRANSITION,
    DISCONNECTED_GROUPS,
    FULLY_CROSSED,
    IDENTIFICATION_FIXED_POPULATION_PRIOR,
    IDENTIFICATION_NONE,
    IDENTIFICATION_PLANNED_PERSON_RATER_PLUS_FIXED_PRIOR,
    NESTED_BRIDGE,
    PERSON_ROLE_STUDY,
    PRECISION_CHANGE_NOT_COMPARABLE,
    PRECISION_NOT_EVALUATED,
    RANDOM_STREAM_RESPONSE,
    STRUCTURE_NOT_EVALUATED_LIMIT,
    ClassificationRuleV1,
    CostPolicyV1,
    DesignSpecV1,
    MonteCarloSpecV1,
    RandomizationSpecV1,
    RsmTruthSpecV1,
    SimulationScenarioSpecV1,
    SimulationConditionValidationError,
    StructuralTransitionValidationError,
    build_formula_plan_summary,
    build_structural_transition_preflight,
    compile_rating_design,
    evaluate_plan_structure,
    jmle_estimator_spec,
    keyed_uniform01,
    mml_estimator_spec,
    normalize_structural_transition_preflight,
    propose_structural_bridge_repair,
    symmetric_adjacent_thresholds,
)


def _disconnected_spec() -> DesignSpecV1:
    return DesignSpecV1(
        plan=NESTED_BRIDGE,
        n_persons=4,
        n_raters=4,
        score_records_per_session=3,
        n_bridge_artifacts=0,
        bridge_extra_raters_per_artifact=0,
        assignment_seed=20260722,
    )


def _scenario(method: str) -> SimulationScenarioSpecV1:
    truth = RsmTruthSpecV1(
        rating_min=0,
        n_categories=5,
        adjacent_thresholds=symmetric_adjacent_thresholds(5),
    )
    estimator = (
        jmle_estimator_spec(max_iterations=100)
        if method == "JMLE"
        else mml_estimator_spec(max_iterations=100, population_prior_sd=1.0)
    )
    return SimulationScenarioSpecV1(
        truth=truth,
        estimator=estimator,
        monte_carlo=MonteCarloSpecV1(
            requested_replicates=100,
            randomization=RandomizationSpecV1(master_seed=41),
        ),
        classification_rule=ClassificationRuleV1(cut_scores=(0.0,)),
    )


def test_disconnected_jmle_is_blocked_and_improvement_claims_are_none():
    structural = evaluate_plan_structure(build_formula_plan_summary(_disconnected_spec()))
    preflight = build_structural_transition_preflight(structural, _scenario("JMLE"))

    assert preflight.comparison_kind == COMPARISON_ESTIMABILITY_TRANSITION
    assert preflight.planned_person_rater_graph_connected_before is False
    assert preflight.planned_person_rater_graph_connected_after is True
    assert preflight.base_fit_preflight_status == BASE_JMLE_BLOCKED
    assert preflight.base_identification_basis == IDENTIFICATION_NONE
    assert preflight.base_claim_tier == CLAIM_NOT_ESTIMABLE_COMMON_SCALE
    assert preflight.precision_change_status == PRECISION_CHANGE_NOT_COMPARABLE
    assert preflight.primary_claim_eligible is False
    assert preflight.paired_precision_evaluation_applicable is False
    assert preflight.precision_evaluation_status == PRECISION_NOT_EVALUATED
    assert preflight.expected_se_improvement_percent is None
    assert preflight.recovery_improvement_percent is None
    assert preflight.classification_improvement_percent is None
    assert preflight.sample_size_recommendation is None


def test_disconnected_mml_is_separate_prior_conditional_claim_tier():
    structural = evaluate_plan_structure(build_formula_plan_summary(_disconnected_spec()))
    preflight = build_structural_transition_preflight(structural, _scenario("MML"))

    assert preflight.base_fit_preflight_status == BASE_MML_PRIOR_CONDITIONAL
    assert preflight.base_identification_basis == IDENTIFICATION_FIXED_POPULATION_PRIOR
    assert preflight.base_claim_tier == CLAIM_PRIOR_CONDITIONAL_EXPLORATORY
    assert preflight.repaired_identification_basis == (
        IDENTIFICATION_PLANNED_PERSON_RATER_PLUS_FIXED_PRIOR
    )
    assert (
        preflight.repaired_claim_tier
        == CLAIM_PLANNED_CONNECTIVITY_PLUS_PRIOR_PENDING_FIT
    )
    assert preflight.primary_claim_eligible is False
    assert preflight.expected_se_improvement_percent is None


def test_nonclaim_fields_cannot_be_forged_into_a_percentage():
    structural = evaluate_plan_structure(build_formula_plan_summary(_disconnected_spec()))
    preflight = build_structural_transition_preflight(structural, _scenario("JMLE"))

    with pytest.raises(StructuralTransitionValidationError, match="must remain None"):
        replace(preflight, expected_se_improvement_percent=12.5)
    with pytest.raises(StructuralTransitionValidationError, match="base estimator status"):
        replace(preflight, base_fit_preflight_status="converged")
    with pytest.raises(StructuralTransitionValidationError, match="components_before - 1"):
        replace(preflight, components_before=5)
    with pytest.raises(StructuralTransitionValidationError, match="formula_evaluation_id"):
        replace(preflight, formula_evaluation_id="0" * 16)
    with pytest.raises(StructuralTransitionValidationError, match="not applicable"):
        replace(preflight, paired_precision_evaluation_applicable=True)


def test_preflight_payload_is_compact_json_and_does_not_retain_overlay_rows():
    structural = evaluate_plan_structure(build_formula_plan_summary(_disconnected_spec()))
    preflight = build_structural_transition_preflight(structural, _scenario("JMLE"))
    payload = preflight.to_dict()

    assert json.loads(json.dumps(payload, allow_nan=False)) == payload
    encoded = json.dumps(payload, sort_keys=True)
    assert "overlay_sessions" not in encoded
    assert "assignments" not in encoded
    assert payload["base_rating_sessions"] == 4
    assert payload["repaired_rating_sessions"] == 7
    assert payload["base_score_records"] == 12
    assert payload["repaired_score_records"] == 21
    assert payload["score_records_per_session"] == 3
    assert payload["base_calibration_events"] == 0
    assert payload["repaired_calibration_events"] == 0
    assert payload["cost_policy"] == preflight.cost_policy.to_dict()
    restored = normalize_structural_transition_preflight(
        json.loads(json.dumps(payload, allow_nan=False))
    )
    assert restored == preflight


def test_saved_preflight_rejects_missing_unknown_and_nested_tampering():
    preflight = build_structural_transition_preflight(
        evaluate_plan_structure(build_formula_plan_summary(_disconnected_spec())),
        _scenario("JMLE"),
    )
    payload = preflight.to_dict()
    missing = dict(payload)
    missing.pop("repair_fingerprint")
    with pytest.raises(StructuralTransitionValidationError, match="missing fields"):
        normalize_structural_transition_preflight(missing)

    unknown = dict(payload, display_label="recommended")
    with pytest.raises(StructuralTransitionValidationError, match="unknown fields"):
        normalize_structural_transition_preflight(unknown)

    nested = json.loads(json.dumps(payload))
    nested["scenario"]["truth"]["missingness"] = "random"
    with pytest.raises(SimulationConditionValidationError, match="missingness"):
        normalize_structural_transition_preflight(nested)

    repriced_without_identity = json.loads(json.dumps(payload))
    repriced_without_identity["cost_policy"]["cost_unit"] = "undeclared-unit"
    with pytest.raises(
        StructuralTransitionValidationError,
        match="cost_policy_fingerprint",
    ):
        normalize_structural_transition_preflight(repriced_without_identity)

    huge_cost = dict(payload, base_total_cost=10**10_000)
    with pytest.raises(StructuralTransitionValidationError, match="finite"):
        normalize_structural_transition_preflight(huge_cost)

    huge_count = dict(payload, components_before=10**5_000)
    with pytest.raises(StructuralTransitionValidationError, match="<="):
        normalize_structural_transition_preflight(huge_count)

    added_calibration = dict(payload, repaired_calibration_events=1)
    with pytest.raises(StructuralTransitionValidationError, match="cannot add"):
        normalize_structural_transition_preflight(added_calibration)


def test_repair_union_reuses_every_base_response_draw_and_only_adds_new_keys():
    bundle = compile_rating_design(_disconnected_spec())
    repair = propose_structural_bridge_repair(bundle)
    criterion_ids = tuple(
        f"C{index:03d}"
        for index in range(1, bundle.design_spec.score_records_per_session + 1)
    )
    base_keys = {
        (
            row.person_role,
            row.person_id,
            row.artifact_id,
            row.rater_id,
            criterion_id,
        )
        for row in bundle.assignments
        for criterion_id in criterion_ids
    }
    added_keys = {
        (
            PERSON_ROLE_STUDY,
            row.person_id,
            row.artifact_id,
            row.additional_rater_id,
            criterion_id,
        )
        for row in repair.overlay_sessions
        for criterion_id in criterion_ids
    }
    assert base_keys.isdisjoint(added_keys)
    assert len(added_keys) == repair.added_score_records

    randomization = RandomizationSpecV1(master_seed=41)
    base_draws = {
        key: keyed_uniform01(
            randomization,
            replicate_index=3,
            stream=RANDOM_STREAM_RESPONSE,
            key_parts=key,
        )
        for key in base_keys
    }
    repaired_draws = {
        key: keyed_uniform01(
            randomization,
            replicate_index=3,
            stream=RANDOM_STREAM_RESPONSE,
            key_parts=key,
        )
        for key in reversed(sorted(base_keys | added_keys))
    }
    assert all(repaired_draws[key] == draw for key, draw in base_draws.items())
    assert set(repaired_draws) - set(base_draws) == added_keys


def test_repricing_changes_only_economic_not_scientific_transition_identity():
    spec = _disconnected_spec()
    scenario = _scenario("JMLE")
    minutes = evaluate_plan_structure(build_formula_plan_summary(
        spec,
        CostPolicyV1(
            rating_session_unit_cost=2.0,
            score_record_unit_cost=0.0,
            calibration_event_unit_cost=0.0,
            cost_unit="staff-minutes",
        ),
    ))
    currency = evaluate_plan_structure(build_formula_plan_summary(
        spec,
        CostPolicyV1(
            rating_session_unit_cost=10.0,
            score_record_unit_cost=0.5,
            calibration_event_unit_cost=0.0,
            cost_unit="JPY-equivalent",
        ),
    ))

    first = build_structural_transition_preflight(minutes, scenario)
    second = build_structural_transition_preflight(currency, scenario)
    assert first.transition_id == second.transition_id
    assert first.economic_context_id != second.economic_context_id
    assert first.formula_evaluation_id != second.formula_evaluation_id
    assert first.cost_policy_fingerprint != second.cost_policy_fingerprint


def test_preflight_retains_calibration_count_for_standalone_cost_audit():
    spec = replace(_disconnected_spec(), calibration_events=5)
    policy = CostPolicyV1(
        rating_session_unit_cost=1.25,
        score_record_unit_cost=0.2,
        calibration_event_unit_cost=3.5,
        cost_unit="staff-minutes",
    )
    preflight = build_structural_transition_preflight(
        evaluate_plan_structure(build_formula_plan_summary(spec, policy)),
        _scenario("JMLE"),
    )

    assert preflight.base_calibration_events == 5
    assert preflight.repaired_calibration_events == 5
    assert preflight.base_total_cost == (
        preflight.base_rating_sessions * policy.rating_session_unit_cost
        + preflight.base_score_records * policy.score_record_unit_cost
        + preflight.base_calibration_events * policy.calibration_event_unit_cost
    )


def test_connected_zero_repair_requires_future_precision_overlay():
    connected = evaluate_plan_structure(build_formula_plan_summary(DesignSpecV1(
        plan=FULLY_CROSSED,
        n_persons=4,
        n_raters=2,
        raters_per_artifact=2,
    )))
    with pytest.raises(StructuralTransitionValidationError, match="precision-enhancement"):
        build_structural_transition_preflight(connected, _scenario("JMLE"))


def test_not_evaluated_limit_cannot_enter_transition_preflight():
    limited = evaluate_plan_structure(
        build_formula_plan_summary(_disconnected_spec()),
        max_rating_sessions=1,
    )
    assert limited.structure_evaluation_status == STRUCTURE_NOT_EVALUATED_LIMIT
    with pytest.raises(StructuralTransitionValidationError, match="explicitly evaluated"):
        build_structural_transition_preflight(limited, _scenario("JMLE"))


@pytest.mark.parametrize(
    "capacity_override",
    (
        {"n_artifacts_per_person": 1_000},
        {"score_records_per_session": 1_000},
    ),
)
def test_transition_rejects_designs_beyond_fixed_width_simulation_ids(
    capacity_override,
):
    spec = DesignSpecV1(
        plan=NESTED_BRIDGE,
        n_persons=2,
        n_raters=2,
        n_bridge_artifacts=0,
        bridge_extra_raters_per_artifact=0,
        **capacity_override,
    )
    structural = evaluate_plan_structure(build_formula_plan_summary(spec))
    with pytest.raises(StructuralTransitionValidationError, match="ID capacity"):
        build_structural_transition_preflight(structural, _scenario("JMLE"))


def test_connected_common_anchor_plan_is_not_mislabeled_as_k_minus_one_transition():
    anchored = evaluate_plan_structure(build_formula_plan_summary(DesignSpecV1(
        plan=DISCONNECTED_GROUPS,
        n_persons=4,
        n_raters=4,
        n_groups=2,
        n_common_anchor_artifacts=1,
        common_anchor_raters_per_group=2,
    )))
    assert anchored.primary_study is not None
    assert anchored.anchor_augmented is not None
    assert anchored.primary_study.connected is False
    assert anchored.anchor_augmented.connected is True
    with pytest.raises(StructuralTransitionValidationError, match="precision-enhancement"):
        build_structural_transition_preflight(anchored, _scenario("MML"))


def test_transition_import_stays_free_of_heavy_and_ui_dependencies():
    project_root = Path(__file__).resolve().parents[1]
    code = "\n".join((
        "import sys",
        "import mfrm_app.simulation.transition",
        "for name in ('streamlit', 'pandas', 'numpy', 'networkx', 'scipy', 'plotly'):",
        "    assert name not in sys.modules, name",
    ))
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
    subprocess.run(
        [sys.executable, "-c", code],
        cwd=project_root,
        env=env,
        check=True,
        capture_output=True,
        text=True,
    )
