"""Prospective Person-Rater topology and structural repair contracts."""

from __future__ import annotations

from dataclasses import replace
import json
import time

import pytest

from mfrm_app.simulation import (
    ANCHOR_AUGMENTED_SCOPE,
    DISCONNECTED_GROUPS,
    FULLY_CROSSED,
    NESTED_BRIDGE,
    PRECISION_NOT_EVALUATED,
    PRIMARY_STUDY_SCOPE,
    REPAIR_NOT_NEEDED,
    REPAIR_PROPOSED,
    SPIRAL,
    AssignmentLimitError,
    DesignSpecV1,
    TopologyComponentV1,
    TopologyScopeAuditV1,
    TopologyValidationError,
    assignment_bundle_to_dict,
    audit_rating_design_topology,
    compile_rating_design,
    propose_structural_bridge_repair,
    structural_repair_to_dict,
    topology_audit_to_dict,
)


def _component_count_after_overlay(scope, overlay_sessions) -> int:
    nodes = set()
    parent = {}

    def add(node):
        nodes.add(node)
        parent.setdefault(node, node)

    def find(node):
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    def union(left, right):
        add(left)
        add(right)
        root_left = find(left)
        root_right = find(right)
        if root_left != root_right:
            parent[root_right] = root_left

    for component in scope.components:
        for person_id in component.study_person_ids:
            add(("study", person_id))
        for person_id in component.common_anchor_person_ids:
            add(("common_anchor", person_id))
        for rater_id in component.rater_ids:
            add(("rater", rater_id))
    for edge in scope.incidence_edges:
        union((edge.person_role, edge.person_id), ("rater", edge.rater_id))
    for session in overlay_sessions:
        union(("study", session.person_id), ("rater", session.additional_rater_id))
    return len({find(node) for node in nodes})


def test_fully_crossed_topology_is_connected_and_collapses_parallel_artifacts():
    bundle = compile_rating_design(DesignSpecV1(
        plan=FULLY_CROSSED,
        n_persons=2,
        n_raters=2,
        n_artifacts_per_person=3,
        raters_per_artifact=2,
    ))

    audit = audit_rating_design_topology(bundle)
    primary = audit.primary_study

    assert primary.scope == PRIMARY_STUDY_SCOPE
    assert audit.anchor_augmented.scope == ANCHOR_AUGMENTED_SCOPE
    assert primary.connected is True
    assert primary.n_components == 1
    assert primary.rating_sessions == 12
    assert primary.unique_incidence_edges == 4
    assert primary.repeated_sessions_collapsed == 8
    assert [edge.rating_sessions for edge in primary.incidence_edges] == [3, 3, 3, 3]
    assert sum(edge.rating_sessions for edge in primary.incidence_edges) == 12
    assert primary.rater_loads[0].rating_sessions == 6
    assert primary.rater_loads[0].unique_persons == 2


@pytest.mark.parametrize(
    ("anchors", "raters_per_group", "expected_augmented_components"),
    [
        (1, 1, 3),
        (1, 2, 1),
        (2, 1, 2),
    ],
)
def test_common_anchor_scope_is_derived_from_actual_assignments(
    anchors,
    raters_per_group,
    expected_augmented_components,
):
    """Anchor counts and plan labels alone must never imply connectivity."""
    bundle = compile_rating_design(DesignSpecV1(
        plan=DISCONNECTED_GROUPS,
        n_persons=2,
        n_raters=4,
        n_groups=2,
        n_common_anchor_artifacts=anchors,
        common_anchor_raters_per_group=raters_per_group,
        assignment_seed=0,
    ))

    audit = audit_rating_design_topology(bundle)

    assert audit.primary_study.n_components == 4
    assert audit.primary_study.isolated_raters == 2
    assert audit.anchor_augmented.n_components == expected_augmented_components
    assert audit.anchor_augmented.declared_common_anchor_persons == anchors
    assert audit.primary_study.declared_common_anchor_persons == 0
    if raters_per_group == 2:
        assert audit.primary_study.connected is False
        assert audit.anchor_augmented.connected is True


def test_spiral_person_node_can_link_artifacts_rated_by_different_raters():
    bundle = compile_rating_design(DesignSpecV1(
        plan=SPIRAL,
        n_persons=1,
        n_raters=2,
        n_artifacts_per_person=2,
        raters_per_artifact=1,
        assignment_seed=0,
    ))

    audit = audit_rating_design_topology(bundle)

    assert audit.primary_study.connected is True
    assert audit.primary_study.n_components == 1
    assert audit.primary_study.unique_incidence_edges == 2
    assert audit.primary_study.rating_sessions == 2


def test_nested_repeated_artifacts_are_edge_multiplicity_not_extra_components():
    bundle = compile_rating_design(DesignSpecV1(
        plan=NESTED_BRIDGE,
        n_persons=2,
        n_raters=2,
        n_artifacts_per_person=3,
    ))

    audit = audit_rating_design_topology(bundle)
    repair = propose_structural_bridge_repair(bundle, topology_audit=audit)

    assert audit.primary_study.n_components == 2
    assert audit.primary_study.unique_incidence_edges == 2
    assert audit.primary_study.rating_sessions == 6
    assert [edge.rating_sessions for edge in audit.primary_study.incidence_edges] == [3, 3]
    assert repair.added_rating_sessions == 1
    assert repair.added_score_records == 1


def test_declared_bridge_sessions_are_part_of_the_primary_study_graph():
    bundle = compile_rating_design(DesignSpecV1(
        plan=NESTED_BRIDGE,
        n_persons=2,
        n_raters=3,
        n_artifacts_per_person=2,
        n_bridge_artifacts=2,
        bridge_extra_raters_per_artifact=1,
        assignment_seed=0,
    ))

    audit = audit_rating_design_topology(bundle)

    assert audit.primary_study.connected is True
    assert audit.primary_study.n_components == 1
    assert audit.primary_study.rating_sessions == 6
    assert audit.primary_study.unique_incidence_edges == 4


def test_structural_repair_is_a_separate_exact_k_minus_one_overlay():
    bundle = compile_rating_design(DesignSpecV1(
        plan=DISCONNECTED_GROUPS,
        n_persons=4,
        n_raters=4,
        n_groups=2,
        score_records_per_session=4,
        assignment_seed=17,
    ))
    audit = audit_rating_design_topology(bundle)

    repair = propose_structural_bridge_repair(bundle, topology_audit=audit)

    assert repair.status == REPAIR_PROPOSED
    assert repair.components_before == 4
    assert repair.structural_minimum_rating_sessions == 3
    assert repair.added_rating_sessions == 3
    assert repair.added_score_records == 12
    assert repair.added_calibration_events == 0
    assert repair.projected_components_after == 1
    assert len(bundle.assignments) == bundle.workload.total_rating_sessions
    assert _component_count_after_overlay(
        audit.anchor_augmented,
        repair.overlay_sessions,
    ) == 1

    base_session_keys = {row.session_key for row in bundle.assignments}
    component_by_study_person = {
        person_id: component.component_id
        for component in audit.anchor_augmented.components
        for person_id in component.study_person_ids
    }
    component_by_rater = {
        rater_id: component.component_id
        for component in audit.anchor_augmented.components
        for rater_id in component.rater_ids
    }
    for index, session in enumerate(repair.overlay_sessions, start=1):
        assert session.repair_index == index
        assert (
            session.person_id,
            session.artifact_id,
            session.additional_rater_id,
        ) not in base_session_keys
        assert component_by_study_person[session.person_id] == session.source_component_id
        assert component_by_rater[session.additional_rater_id] == session.target_component_id
        assert session.source_component_id != session.target_component_id
        assert session.added_score_records == 4
    assert len({row.target_component_id for row in repair.overlay_sessions}) == 3
    assert len({row.source_component_id for row in repair.overlay_sessions}) == 1


def test_connected_augmented_schedule_needs_no_structural_overlay():
    bundle = compile_rating_design(DesignSpecV1(
        plan=DISCONNECTED_GROUPS,
        n_persons=2,
        n_raters=4,
        n_groups=2,
        n_common_anchor_artifacts=1,
        common_anchor_raters_per_group=2,
        score_records_per_session=3,
        assignment_seed=0,
    ))

    audit = audit_rating_design_topology(bundle)
    repair = propose_structural_bridge_repair(bundle, topology_audit=audit)

    assert audit.primary_study.connected is False
    assert audit.anchor_augmented.connected is True
    assert repair.status == REPAIR_NOT_NEEDED
    assert repair.overlay_sessions == ()
    assert repair.added_rating_sessions == 0
    assert repair.added_score_records == 0
    assert repair.added_calibration_events == 0


def test_repair_contract_hard_codes_precision_nonclaims():
    bundle = compile_rating_design(DesignSpecV1(
        plan=NESTED_BRIDGE,
        n_persons=2,
        n_raters=2,
    ))
    repair = propose_structural_bridge_repair(bundle)

    assert repair.structural_connectivity_only is True
    assert repair.precision_evaluation_status == PRECISION_NOT_EVALUATED
    assert repair.expected_se_improvement_percent is None
    assert repair.classification_agreement_improvement_percent is None

    with pytest.raises(TopologyValidationError, match="expected SE improvement"):
        replace(repair, expected_se_improvement_percent=10.0)


def test_audit_and_repair_are_json_safe_and_bound_to_the_base_schedule():
    bundle = compile_rating_design(DesignSpecV1(
        plan=SPIRAL,
        n_persons=4,
        n_raters=3,
        raters_per_artifact=2,
        score_records_per_session=2,
        assignment_seed=41,
    ))
    mapping = assignment_bundle_to_dict(bundle)

    audit = audit_rating_design_topology(mapping)
    repair = propose_structural_bridge_repair(mapping, topology_audit=audit)
    audit_payload = topology_audit_to_dict(audit)
    repair_payload = structural_repair_to_dict(repair)

    assert json.loads(json.dumps(audit_payload, allow_nan=False)) == audit_payload
    assert json.loads(json.dumps(repair_payload, allow_nan=False)) == repair_payload
    assert audit.cache_identity[2:] == (
        bundle.design_spec_fingerprint,
        bundle.assignment_fingerprint,
    )
    assert repair.design_spec_fingerprint == bundle.design_spec_fingerprint
    assert repair.assignment_fingerprint == bundle.assignment_fingerprint

    other_bundle = compile_rating_design(replace(bundle.design_spec, assignment_seed=42))
    with pytest.raises(TopologyValidationError, match="does not belong"):
        propose_structural_bridge_repair(other_bundle, topology_audit=audit)
    with pytest.raises(TopologyValidationError, match="repair_fingerprint"):
        replace(repair, repair_fingerprint="0" * 16)


def test_repair_recomputes_topology_instead_of_trusting_relabelled_fingerprints():
    """Matching 16-hex labels cannot legitimize scopes from another schedule."""
    base_spec = DesignSpecV1(
        plan=SPIRAL,
        n_persons=2,
        n_raters=4,
        n_artifacts_per_person=2,
        raters_per_artifact=1,
        assignment_seed=0,
    )
    first_bundle = compile_rating_design(base_spec)
    second_bundle = compile_rating_design(replace(base_spec, assignment_seed=1))
    first_audit = audit_rating_design_topology(first_bundle)
    relabelled_audit = replace(
        first_audit,
        design_spec_fingerprint=second_bundle.design_spec_fingerprint,
        assignment_fingerprint=second_bundle.assignment_fingerprint,
    )

    with pytest.raises(TopologyValidationError, match="exactly reproduce"):
        propose_structural_bridge_repair(
            second_bundle,
            topology_audit=relabelled_audit,
        )


def test_scope_contract_rejects_disconnected_edges_claimed_as_one_component():
    bundle = compile_rating_design(DesignSpecV1(
        plan=NESTED_BRIDGE,
        n_persons=2,
        n_raters=2,
    ))
    scope = audit_rating_design_topology(bundle).primary_study
    collapsed_component = TopologyComponentV1(
        schema_version=scope.components[0].schema_version,
        component_id="C001",
        study_person_ids=tuple(sorted(
            person_id
            for component in scope.components
            for person_id in component.study_person_ids
        )),
        common_anchor_person_ids=(),
        rater_ids=tuple(sorted(
            rater_id
            for component in scope.components
            for rater_id in component.rater_ids
        )),
        unique_incidence_edges=scope.unique_incidence_edges,
        rating_sessions=scope.rating_sessions,
    )
    relabelled_edges = tuple(
        replace(edge, component_id="C001") for edge in scope.incidence_edges
    )

    with pytest.raises(TopologyValidationError, match="not connected"):
        replace(
            scope,
            components=(collapsed_component,),
            incidence_edges=relabelled_edges,
            n_components=1,
            largest_component_nodes=4,
            largest_component_study_persons=2,
            largest_component_raters=2,
            connected=True,
        )


def test_scope_contract_cannot_represent_an_empty_compiled_design():
    with pytest.raises(TopologyValidationError, match="declared_study_persons"):
        TopologyScopeAuditV1(
            schema_version="mfrm_rating_topology_scope_v1",
            scope=PRIMARY_STUDY_SCOPE,
            declared_study_persons=0,
            declared_common_anchor_persons=0,
            declared_raters=0,
            nodes=0,
            unique_incidence_edges=0,
            rating_sessions=0,
            repeated_sessions_collapsed=0,
            n_components=0,
            largest_component_nodes=0,
            largest_component_study_persons=0,
            largest_component_raters=0,
            isolated_raters=0,
            connected=False,
            components=(),
            incidence_edges=(),
            rater_loads=(),
        )


def test_topology_and_repair_thread_the_explicit_materialization_limit():
    bundle = compile_rating_design(DesignSpecV1(
        plan=FULLY_CROSSED,
        n_persons=3,
        n_raters=2,
        raters_per_artifact=2,
    ))
    with pytest.raises(AssignmentLimitError):
        audit_rating_design_topology(bundle, max_rating_sessions=5)
    audit = audit_rating_design_topology(bundle, max_rating_sessions=6)
    with pytest.raises(AssignmentLimitError):
        propose_structural_bridge_repair(
            bundle,
            topology_audit=audit,
            max_rating_sessions=5,
        )
    repair = propose_structural_bridge_repair(
        bundle,
        topology_audit=audit,
        max_rating_sessions=6,
    )
    assert repair.status == REPAIR_NOT_NEEDED


def test_many_disconnected_components_stay_in_the_linear_preview_layer():
    bundle = compile_rating_design(DesignSpecV1(
        plan=NESTED_BRIDGE,
        n_persons=5_000,
        n_raters=5_000,
    ))

    started = time.perf_counter()
    audit = audit_rating_design_topology(bundle)
    elapsed = time.perf_counter() - started

    assert audit.primary_study.n_components == 5_000
    assert audit.primary_study.unique_incidence_edges == 5_000
    assert elapsed < 2.0


def test_repair_component_ids_keep_numeric_order_beyond_c999():
    bundle = compile_rating_design(DesignSpecV1(
        plan=NESTED_BRIDGE,
        n_persons=1_000,
        n_raters=1_000,
    ))

    repair = propose_structural_bridge_repair(bundle)

    assert repair.components_before == 1_000
    assert repair.base_component_ids[998:] == ("C999", "C1000")
    assert repair.added_rating_sessions == 999


def test_repair_serializes_the_common_existing_rater_set_only_once():
    n_persons = 200
    bundle = compile_rating_design(DesignSpecV1(
        plan=NESTED_BRIDGE,
        n_persons=n_persons,
        n_raters=2 * n_persons,
        n_bridge_artifacts=1,
        bridge_extra_raters_per_artifact=n_persons,
    ))

    repair = propose_structural_bridge_repair(bundle)
    payload = structural_repair_to_dict(repair)

    assert len(repair.target_existing_rater_ids) == n_persons + 1
    assert len(repair.overlay_sessions) == n_persons - 1
    assert all(
        "existing_rater_ids" not in session_payload
        for session_payload in payload["overlay_sessions"]
    )
    assert json.dumps(payload).count('"target_existing_rater_ids"') == 1


def test_standard_target_stays_inside_immediate_topology_layer():
    bundle = compile_rating_design(DesignSpecV1(
        plan=FULLY_CROSSED,
        n_persons=200,
        n_raters=6,
        raters_per_artifact=6,
        score_records_per_session=4,
    ))

    audit = audit_rating_design_topology(bundle)

    assert audit.primary_study.rating_sessions == 1_200
    assert audit.primary_study.unique_incidence_edges == 1_200
    assert audit.primary_study.n_components == 1
    assert len(audit.primary_study.rater_loads) == 6
