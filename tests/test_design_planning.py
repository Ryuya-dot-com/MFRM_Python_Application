"""Compact two-stage prospective design-planning contracts."""

from __future__ import annotations

from dataclasses import replace
import json
from pathlib import Path
import subprocess
import sys

import pytest

import mfrm_app.simulation.planning as planning_module
from mfrm_app.simulation import (
    DISCONNECTED_GROUPS,
    FULLY_CROSSED,
    NESTED_BRIDGE,
    PRECISION_NOT_EVALUATED,
    REPAIR_NOT_NEEDED,
    REPAIR_PROPOSED,
    SPIRAL,
    STRUCTURE_EVALUATED,
    STRUCTURE_NOT_EVALUATED,
    STRUCTURE_NOT_EVALUATED_LIMIT,
    CostPolicyV1,
    DesignSpecV1,
    apply_cost_policy,
    build_formula_plan_summary,
    compute_design_workload,
    compute_repair_cost_delta,
    cost_policy_fingerprint,
    design_spec_fingerprint,
    design_spec_to_dict,
    evaluate_plan_structure,
)


FOUR_PLAN_SPECS = (
    DesignSpecV1(
        plan=FULLY_CROSSED,
        n_persons=3,
        n_raters=2,
        n_artifacts_per_person=2,
        score_records_per_session=2,
        raters_per_artifact=2,
        calibration_events=1,
        assignment_seed=11,
    ),
    DesignSpecV1(
        plan=SPIRAL,
        n_persons=4,
        n_raters=3,
        n_artifacts_per_person=2,
        score_records_per_session=3,
        raters_per_artifact=2,
        assignment_seed=12,
    ),
    DesignSpecV1(
        plan=NESTED_BRIDGE,
        n_persons=4,
        n_raters=3,
        n_artifacts_per_person=2,
        score_records_per_session=2,
        n_bridge_artifacts=1,
        bridge_extra_raters_per_artifact=1,
        assignment_seed=13,
    ),
    DesignSpecV1(
        plan=DISCONNECTED_GROUPS,
        n_persons=4,
        n_raters=4,
        n_groups=2,
        score_records_per_session=4,
        n_common_anchor_artifacts=1,
        common_anchor_raters_per_group=2,
        assignment_seed=14,
    ),
)


def _all_mapping_keys(value):
    if isinstance(value, dict):
        for key, child in value.items():
            yield key
            yield from _all_mapping_keys(child)
    elif isinstance(value, list):
        for child in value:
            yield from _all_mapping_keys(child)


@pytest.mark.parametrize("spec", FOUR_PLAN_SPECS, ids=lambda spec: spec.plan)
def test_four_plans_match_direct_formula_and_lazy_structure_contracts(spec):
    policy = CostPolicyV1(
        rating_session_unit_cost=2,
        score_record_unit_cost=0.25,
        calibration_event_unit_cost=5,
        cost_unit="staff-minutes",
    )

    formula = build_formula_plan_summary(spec, policy)
    direct_workload = compute_design_workload(spec)

    assert formula.workload == direct_workload
    assert formula.weighted_cost == apply_cost_policy(direct_workload, policy)
    assert formula.design_spec_fingerprint == design_spec_fingerprint(spec)
    assert formula.cost_policy_fingerprint == cost_policy_fingerprint(policy)
    assert formula.structure_evaluation_status == STRUCTURE_NOT_EVALUATED

    structural = evaluate_plan_structure(
        formula,
        max_rating_sessions=direct_workload.total_rating_sessions,
    )
    assert structural.structure_evaluation_status == STRUCTURE_EVALUATED
    assert structural.required_rating_sessions == direct_workload.total_rating_sessions
    assert structural.formula_summary is formula
    assert structural.primary_study is not None
    assert structural.anchor_augmented is not None
    assert structural.repair is not None
    assert structural.repair_cost_delta is not None


def test_formula_and_structure_payloads_are_json_safe_and_compact():
    formula = build_formula_plan_summary(FOUR_PLAN_SPECS[2])
    structural = evaluate_plan_structure(formula)

    formula_payload = formula.to_dict()
    structural_payload = structural.to_dict()
    assert json.loads(json.dumps(formula_payload, allow_nan=False)) == formula_payload
    assert json.loads(json.dumps(structural_payload, allow_nan=False)) == structural_payload

    keys = set(_all_mapping_keys(structural_payload))
    assert "assignments" not in keys
    assert "components" not in keys
    assert "incidence_edges" not in keys
    assert "rater_loads" not in keys
    assert structural_payload["primary_study"]["rater_load_min"] >= 0
    assert structural_payload["anchor_augmented"]["rater_load_max"] >= 1


def test_evaluation_id_is_deterministic_across_mapping_and_fresh_process():
    spec = FOUR_PLAN_SPECS[1]
    expected = build_formula_plan_summary(spec).evaluation_id
    from_mapping = build_formula_plan_summary(
        json.loads(json.dumps(design_spec_to_dict(spec)))
    )
    assert from_mapping.evaluation_id == expected

    project_root = Path(__file__).resolve().parents[1]
    code = "\n".join(
        (
            "from mfrm_app.simulation import DesignSpecV1, SPIRAL, build_formula_plan_summary",
            "spec = DesignSpecV1(plan=SPIRAL, n_persons=4, n_raters=3, "
            "n_artifacts_per_person=2, score_records_per_session=3, "
            "raters_per_artifact=2, assignment_seed=12)",
            "print(build_formula_plan_summary(spec).evaluation_id)",
        )
    )
    completed = subprocess.run(
        [sys.executable, "-c", code],
        cwd=project_root,
        check=True,
        capture_output=True,
        text=True,
    )
    assert completed.stdout.strip() == expected


def test_seed_and_policy_identities_change_only_their_owned_fingerprints():
    base = DesignSpecV1(
        plan=FULLY_CROSSED,
        n_persons=4,
        n_raters=2,
        raters_per_artifact=2,
        assignment_seed=1,
    )
    changed_seed = replace(base, assignment_seed=2)
    default = build_formula_plan_summary(base)
    seeded = build_formula_plan_summary(changed_seed)
    repriced = build_formula_plan_summary(
        base,
        CostPolicyV1(rating_session_unit_cost=2),
    )

    assert default.workload.total_rating_sessions == seeded.workload.total_rating_sessions
    assert default.design_spec_fingerprint != seeded.design_spec_fingerprint
    assert default.cost_policy_fingerprint == seeded.cost_policy_fingerprint
    assert default.evaluation_id != seeded.evaluation_id

    assert default.design_spec_fingerprint == repriced.design_spec_fingerprint
    assert default.cost_policy_fingerprint != repriced.cost_policy_fingerprint
    assert default.evaluation_id != repriced.evaluation_id


def test_repair_cost_and_projected_totals_use_exact_base_plus_delta_arithmetic():
    spec = DesignSpecV1(
        plan=DISCONNECTED_GROUPS,
        n_persons=4,
        n_raters=4,
        n_groups=2,
        score_records_per_session=4,
        calibration_events=2,
        assignment_seed=17,
    )
    policy = CostPolicyV1(
        rating_session_unit_cost=2,
        score_record_unit_cost=0.25,
        calibration_event_unit_cost=10,
        cost_unit="staff-minutes",
    )
    formula = build_formula_plan_summary(spec, policy)
    structural = evaluate_plan_structure(formula)
    repair = structural.repair
    delta = structural.repair_cost_delta

    assert repair is not None and repair.status == REPAIR_PROPOSED
    assert delta is not None
    assert repair.added_rating_sessions == 3
    assert repair.added_score_records == 12
    assert repair.added_calibration_events == 0
    assert delta.design_spec_fingerprint == formula.design_spec_fingerprint
    assert delta.assignment_fingerprint == structural.assignment_fingerprint
    assert delta.repair_fingerprint == repair.repair_fingerprint
    assert delta.cost_policy_fingerprint == formula.cost_policy_fingerprint
    assert delta.rating_session_cost_delta == pytest.approx(6)
    assert delta.score_record_cost_delta == pytest.approx(3)
    assert delta.calibration_event_cost_delta == pytest.approx(0)
    assert delta.total_cost_delta == pytest.approx(9)
    assert structural.projected_total_rating_sessions == (
        formula.workload.total_rating_sessions + 3
    )
    assert structural.projected_total_score_records == (
        formula.workload.total_score_records + 12
    )
    assert structural.projected_calibration_events == 2
    assert structural.projected_total_cost == pytest.approx(
        formula.weighted_cost.total_cost + 9
    )


def test_connected_augmented_design_has_a_zero_repair_delta():
    formula = build_formula_plan_summary(FOUR_PLAN_SPECS[3])
    structural = evaluate_plan_structure(formula)

    assert structural.primary_study is not None
    assert structural.anchor_augmented is not None
    assert structural.primary_study.connected is False
    assert structural.anchor_augmented.connected is True
    assert structural.repair is not None
    assert structural.repair.status == REPAIR_NOT_NEEDED
    assert structural.repair_cost_delta is not None
    assert structural.repair_cost_delta.total_cost_delta == 0
    assert structural.projected_total_rating_sessions == (
        formula.workload.total_rating_sessions
    )
    assert structural.projected_total_cost == formula.weighted_cost.total_cost


def test_over_cap_returns_formula_evidence_without_calling_the_compiler(monkeypatch):
    spec = DesignSpecV1(
        plan=FULLY_CROSSED,
        n_persons=20,
        n_raters=3,
        raters_per_artifact=3,
        score_records_per_session=4,
    )
    formula = build_formula_plan_summary(spec)

    def forbidden_materialization(*args, **kwargs):
        raise AssertionError("over-cap planning must stop before assignment rows")

    monkeypatch.setattr(
        planning_module,
        "compile_rating_design",
        forbidden_materialization,
    )
    structural = evaluate_plan_structure(formula, max_rating_sessions=59)

    assert formula.workload.total_rating_sessions == 60
    assert formula.weighted_cost.total_cost == 60
    assert structural.structure_evaluation_status == STRUCTURE_NOT_EVALUATED_LIMIT
    assert structural.materialization_reason == "rating_session_limit"
    assert structural.required_rating_sessions == 60
    assert structural.max_rating_sessions == 59
    for name in (
        "assignment_bundle_schema_version",
        "assignment_fingerprint",
        "topology_audit_schema_version",
        "primary_study",
        "anchor_augmented",
        "repair",
        "repair_cost_delta",
        "projected_total_rating_sessions",
        "projected_total_score_records",
        "projected_calibration_events",
        "projected_total_cost",
    ):
        assert getattr(structural, name) is None


@pytest.mark.parametrize("stage", ["formula", "structural"])
def test_precision_and_sample_size_nonclaims_are_fixed(stage):
    formula = build_formula_plan_summary(FOUR_PLAN_SPECS[0])
    result = formula if stage == "formula" else evaluate_plan_structure(formula)

    assert result.precision_evaluation_status == PRECISION_NOT_EVALUATED
    assert result.expected_se is None
    assert result.recovery_accuracy is None
    assert result.classification_agreement is None
    assert result.sample_size_recommendation is None

    with pytest.raises(ValueError, match="expected_se"):
        replace(result, expected_se=0.2)
    with pytest.raises(ValueError, match="precision_evaluation_status"):
        replace(result, precision_evaluation_status="evaluated")


def test_evaluation_rejects_conflicting_policy_and_invalid_limits():
    policy = CostPolicyV1(rating_session_unit_cost=2, cost_unit="minutes")
    formula = build_formula_plan_summary(FOUR_PLAN_SPECS[0], policy)

    with pytest.raises(ValueError, match="policy conflicts"):
        evaluate_plan_structure(
            formula,
            CostPolicyV1(rating_session_unit_cost=3, cost_unit="minutes"),
        )

    for invalid in (True, False, 0, -1, 1.5, "100"):
        with pytest.raises(ValueError, match="positive integer"):
            evaluate_plan_structure(formula, max_rating_sessions=invalid)

    incomplete_saved_policy = {"rating_session_unit_cost": 2.0}
    with pytest.raises(ValueError, match="missing fields"):
        build_formula_plan_summary(FOUR_PLAN_SPECS[0], incomplete_saved_policy)


def test_repair_cost_delta_requires_the_exact_repair_type():
    with pytest.raises(TypeError, match="StructuralBridgeRepairV1"):
        compute_repair_cost_delta(object())


@pytest.mark.parametrize(
    ("field", "message"),
    (
        ("design_spec_fingerprint", "repair design"),
        ("assignment_fingerprint", "repair assignment"),
    ),
)
def test_structural_summary_rejects_repair_cost_with_relabelled_source(
    field,
    message,
):
    result = evaluate_plan_structure(
        build_formula_plan_summary(
            DesignSpecV1(
                plan=NESTED_BRIDGE,
                n_persons=2,
                n_raters=2,
            )
        )
    )
    assert result.repair_cost_delta is not None
    replacement = "0" * 16
    if getattr(result.repair_cost_delta, field) == replacement:
        replacement = "1" * 16
    relabelled = replace(result.repair_cost_delta, **{field: replacement})

    with pytest.raises(ValueError, match=message):
        replace(result, repair_cost_delta=relabelled)


def test_repair_cost_delta_rejects_noncanonical_cost_unit():
    result = evaluate_plan_structure(
        build_formula_plan_summary(
            DesignSpecV1(
                plan=NESTED_BRIDGE,
                n_persons=2,
                n_raters=2,
            ),
            CostPolicyV1(rating_session_unit_cost=2, cost_unit="staff-minutes"),
        )
    )
    assert result.repair_cost_delta is not None

    with pytest.raises(ValueError, match="canonical without outer whitespace"):
        replace(result.repair_cost_delta, cost_unit=" staff-minutes ")


def test_structural_summary_recomputes_the_complete_repair_cost_delta():
    result = evaluate_plan_structure(
        build_formula_plan_summary(
            DesignSpecV1(
                plan=DISCONNECTED_GROUPS,
                n_persons=4,
                n_raters=4,
                n_groups=2,
                score_records_per_session=4,
            ),
            CostPolicyV1(
                rating_session_unit_cost=2,
                score_record_unit_cost=0.25,
                cost_unit="staff-minutes",
            ),
        )
    )
    delta = result.repair_cost_delta
    assert delta is not None
    forged = replace(
        delta,
        added_rating_sessions=delta.added_rating_sessions + 1,
        rating_session_cost_delta=delta.rating_session_cost_delta + 2,
        total_cost_delta=delta.total_cost_delta + 2,
    )

    with pytest.raises(ValueError, match="does not exactly match"):
        replace(result, repair_cost_delta=forged)


def test_structural_summary_binds_repair_component_count_to_augmented_scope():
    result = evaluate_plan_structure(
        build_formula_plan_summary(
            DesignSpecV1(
                plan=DISCONNECTED_GROUPS,
                n_persons=4,
                n_raters=4,
                n_groups=2,
            )
        )
    )
    assert result.anchor_augmented is not None
    assert result.anchor_augmented.connected is False
    forged_scope = replace(
        result.anchor_augmented,
        n_components=result.anchor_augmented.n_components + 1,
    )

    with pytest.raises(ValueError, match="components_before"):
        replace(result, anchor_augmented=forged_scope)


@pytest.mark.parametrize(
    ("scope_name", "message"),
    (
        ("primary_study", "primary-study sessions"),
        ("anchor_augmented", "anchor-augmented sessions"),
    ),
)
def test_structural_summary_binds_compact_scope_sessions_to_formula_workload(
    scope_name,
    message,
):
    result = evaluate_plan_structure(build_formula_plan_summary(FOUR_PLAN_SPECS[2]))
    scope = getattr(result, scope_name)
    assert scope is not None
    forged_scope = replace(scope, rating_sessions=scope.rating_sessions + 1)

    with pytest.raises(ValueError, match=message):
        replace(result, **{scope_name: forged_scope})


def test_planning_import_stays_free_of_heavy_or_ui_dependencies():
    project_root = Path(__file__).resolve().parents[1]
    code = "\n".join(
        (
            "import sys",
            "import mfrm_app.simulation.planning",
            "for name in ('streamlit', 'pandas', 'numpy', 'networkx', 'scipy', 'plotly'):",
            "    assert name not in sys.modules, name",
        )
    )
    subprocess.run(
        [sys.executable, "-c", code],
        cwd=project_root,
        check=True,
        capture_output=True,
        text=True,
    )
