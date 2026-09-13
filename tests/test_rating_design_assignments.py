"""Deterministic prospective rating-session compiler contracts."""

from __future__ import annotations

from dataclasses import FrozenInstanceError, replace
import json
from pathlib import Path
import subprocess
import sys
import time

import pytest

import mfrm_app.simulation.assignments as assignment_module
from mfrm_app.simulation import (
    ASSIGNMENT_BUNDLE_VERSION,
    DISCONNECTED_GROUPS,
    FULLY_CROSSED,
    NESTED_BRIDGE,
    PERSON_ROLE_COMMON_ANCHOR,
    PERSON_ROLE_STUDY,
    SESSION_BRIDGE,
    SESSION_COMMON_ANCHOR,
    SESSION_PRIMARY,
    SPIRAL,
    AssignmentBundleV1,
    AssignmentLimitError,
    AssignmentValidationError,
    DesignSpecV1,
    RatingAssignmentV1,
    assignment_bundle_to_dict,
    assignment_fingerprint,
    compile_rating_design,
    normalize_assignment_bundle,
)


def _row_tuple(row: RatingAssignmentV1) -> tuple[str, ...]:
    return (
        row.person_id,
        row.artifact_id,
        row.rater_id,
        row.group_id,
        row.person_role,
        row.session_type,
    )


class _ForgedAssignmentBundle(AssignmentBundleV1):
    """Adversarial subclass that deliberately bypasses base validation."""

    def __post_init__(self) -> None:
        pass


GOLDEN_CASES = (
    (
        DesignSpecV1(
            plan=FULLY_CROSSED,
            n_persons=2,
            n_raters=2,
            raters_per_artifact=2,
            assignment_seed=0,
        ),
        "b9aaa93a08232673",
        (
            ("P000001", "P000001-A001", "R000001", "G001", "study", "primary"),
            ("P000001", "P000001-A001", "R000002", "G001", "study", "primary"),
            ("P000002", "P000002-A001", "R000001", "G001", "study", "primary"),
            ("P000002", "P000002-A001", "R000002", "G001", "study", "primary"),
        ),
    ),
    (
        DesignSpecV1(
            plan=SPIRAL,
            n_persons=2,
            n_raters=4,
            raters_per_artifact=2,
            assignment_seed=0,
        ),
        "ae21ec28dda3078f",
        (
            ("P000001", "P000001-A001", "R000001", "G001", "study", "primary"),
            ("P000001", "P000001-A001", "R000004", "G001", "study", "primary"),
            ("P000002", "P000002-A001", "R000002", "G001", "study", "primary"),
            ("P000002", "P000002-A001", "R000003", "G001", "study", "primary"),
        ),
    ),
    (
        DesignSpecV1(
            plan=NESTED_BRIDGE,
            n_persons=2,
            n_raters=3,
            n_artifacts_per_person=2,
            n_bridge_artifacts=1,
            bridge_extra_raters_per_artifact=1,
            assignment_seed=0,
        ),
        "5bcf0e8139afe80c",
        (
            ("P000001", "P000001-A001", "R000001", "G001", "study", "primary"),
            ("P000001", "P000001-A002", "R000001", "G001", "study", "primary"),
            ("P000001", "P000001-A002", "R000002", "G001", "study", "bridge"),
            ("P000002", "P000002-A001", "R000003", "G001", "study", "primary"),
            ("P000002", "P000002-A002", "R000003", "G001", "study", "primary"),
        ),
    ),
    (
        DesignSpecV1(
            plan=DISCONNECTED_GROUPS,
            n_persons=2,
            n_raters=4,
            n_groups=2,
            n_common_anchor_artifacts=1,
            common_anchor_raters_per_group=2,
            assignment_seed=0,
        ),
        "1066dbf7ba3eaa80",
        (
            (
                "CA000001",
                "CA000001-A001",
                "R000001",
                "G001",
                "common_anchor",
                "common_anchor",
            ),
            (
                "CA000001",
                "CA000001-A001",
                "R000002",
                "G002",
                "common_anchor",
                "common_anchor",
            ),
            (
                "CA000001",
                "CA000001-A001",
                "R000003",
                "G002",
                "common_anchor",
                "common_anchor",
            ),
            (
                "CA000001",
                "CA000001-A001",
                "R000004",
                "G001",
                "common_anchor",
                "common_anchor",
            ),
            ("P000001", "P000001-A001", "R000004", "G001", "study", "primary"),
            ("P000002", "P000002-A001", "R000003", "G002", "study", "primary"),
        ),
    ),
)


@pytest.mark.parametrize(("spec", "fingerprint", "expected_rows"), GOLDEN_CASES)
def test_four_plans_have_cross_version_golden_rows(spec, fingerprint, expected_rows):
    bundle = compile_rating_design(spec)
    assert bundle.schema_version == ASSIGNMENT_BUNDLE_VERSION
    assert bundle.assignment_fingerprint == fingerprint
    assert tuple(_row_tuple(row) for row in bundle.assignments) == expected_rows
    assert len(bundle.assignments) == bundle.workload.total_rating_sessions
    assert {row.rater_id for row in bundle.assignments} == {
        f"R{index:06d}" for index in range(1, spec.n_raters + 1)
    }


def test_golden_compiler_is_identical_in_a_fresh_process():
    project_root = Path(__file__).resolve().parents[1]
    code = "\n".join(
        (
            "from mfrm_app.simulation import DesignSpecV1, SPIRAL, compile_rating_design",
            "spec = DesignSpecV1(plan=SPIRAL, n_persons=2, n_raters=4, "
            "raters_per_artifact=2, assignment_seed=0)",
            "print(compile_rating_design(spec).assignment_fingerprint)",
        )
    )
    completed = subprocess.run(
        [sys.executable, "-c", code],
        cwd=project_root,
        check=True,
        capture_output=True,
        text=True,
    )
    assert completed.stdout.strip() == "ae21ec28dda3078f"


def test_spiral_consumes_flattened_slots_and_activates_every_rater():
    spec = DesignSpecV1(
        plan=SPIRAL,
        n_persons=2,
        n_raters=6,
        raters_per_artifact=3,
        assignment_seed=0,
    )
    bundle = compile_rating_design(spec)
    by_artifact = {}
    for row in bundle.assignments:
        by_artifact.setdefault(row.artifact_id, set()).add(row.rater_id)
    assert [len(raters) for raters in by_artifact.values()] == [3, 3]
    assert len(set().union(*by_artifact.values())) == 6


def test_nested_primary_assignment_is_at_person_level_and_bridges_are_distinct():
    spec = DesignSpecV1(
        plan=NESTED_BRIDGE,
        n_persons=5,
        n_raters=6,
        n_artifacts_per_person=3,
        n_bridge_artifacts=2,
        bridge_extra_raters_per_artifact=1,
        assignment_seed=17,
    )
    rows = compile_rating_design(spec).assignments
    primary = [row for row in rows if row.session_type == SESSION_PRIMARY]
    bridge = [row for row in rows if row.session_type == SESSION_BRIDGE]
    for person_id in {row.person_id for row in primary}:
        assert len({row.rater_id for row in primary if row.person_id == person_id}) == 1
    assert len({(row.person_id, row.artifact_id) for row in bridge}) == 2
    assert len({row.rater_id for row in rows}) == 6
    assert len({row.session_key for row in rows}) == len(rows)


def test_common_anchor_ids_are_shared_but_group_membership_is_preserved():
    spec = DesignSpecV1(
        plan=DISCONNECTED_GROUPS,
        n_persons=5,
        n_raters=7,
        n_groups=3,
        n_common_anchor_artifacts=1,
        common_anchor_raters_per_group=1,
        assignment_seed=0,
    )
    rows = compile_rating_design(spec).assignments
    anchor = [row for row in rows if row.session_type == SESSION_COMMON_ANCHOR]
    assert {row.person_id for row in anchor} == {"CA000001"}
    assert {row.artifact_id for row in anchor} == {"CA000001-A001"}
    assert {row.group_id for row in anchor} == {"G001", "G002", "G003"}
    primary = [row for row in rows if row.session_type == SESSION_PRIMARY]
    person_groups = {}
    rater_groups = {}
    for row in primary:
        person_groups.setdefault(row.person_id, set()).add(row.group_id)
        rater_groups.setdefault(row.rater_id, set()).add(row.group_id)
    assert all(len(groups) == 1 for groups in person_groups.values())
    assert all(len(groups) == 1 for groups in rater_groups.values())


def test_seed_changes_sparse_schedules_but_not_fully_crossed_edges():
    spiral = DesignSpecV1(
        plan=SPIRAL,
        n_persons=10,
        n_raters=4,
        raters_per_artifact=2,
        assignment_seed=1,
    )
    assert compile_rating_design(spiral).assignment_fingerprint != (
        compile_rating_design(replace(spiral, assignment_seed=2)).assignment_fingerprint
    )
    fully = DesignSpecV1(
        plan=FULLY_CROSSED,
        n_persons=5,
        n_raters=3,
        raters_per_artifact=3,
        assignment_seed=1,
    )
    first = compile_rating_design(fully)
    second = compile_rating_design(replace(fully, assignment_seed=2))
    assert first.design_spec_fingerprint != second.design_spec_fingerprint
    assert first.assignment_fingerprint == second.assignment_fingerprint
    assert first.assignments == second.assignments


def test_score_row_multiplier_and_calibration_do_not_materialize_sessions():
    base = DesignSpecV1(
        plan=FULLY_CROSSED,
        n_persons=10,
        n_raters=2,
        raters_per_artifact=2,
    )
    expanded = replace(
        base,
        score_records_per_session=8,
        calibration_events=12,
    )
    first = compile_rating_design(base)
    second = compile_rating_design(expanded)
    assert first.assignments == second.assignments
    assert first.assignment_fingerprint == second.assignment_fingerprint
    assert second.workload.total_score_records == 8 * first.workload.total_score_records
    assert second.workload.calibration_events == 12


def test_materialization_limit_is_inclusive_and_rejects_before_id_generation(
    monkeypatch,
):
    spec = DesignSpecV1(
        plan=FULLY_CROSSED,
        n_persons=20,
        n_raters=3,
        raters_per_artifact=3,
    )
    assert len(compile_rating_design(spec, max_rating_sessions=60).assignments) == 60

    def forbidden_id_generation(*args, **kwargs):
        raise AssertionError("ID generation must not run after a failed preflight")

    monkeypatch.setattr(assignment_module, "_canonical_ids", forbidden_id_generation)
    with pytest.raises(AssignmentLimitError) as exc_info:
        compile_rating_design(spec, max_rating_sessions=59)
    assert exc_info.value.required_rating_sessions == 60
    assert exc_info.value.max_rating_sessions == 59


@pytest.mark.parametrize("limit", [True, False, 0, -1, 1.5, "100"])
def test_materialization_limit_requires_a_positive_builtin_integer(limit):
    spec = DesignSpecV1(
        plan=FULLY_CROSSED,
        n_persons=2,
        n_raters=2,
        raters_per_artifact=2,
    )
    with pytest.raises(ValueError, match="positive integer"):
        compile_rating_design(spec, max_rating_sessions=limit)


def test_billion_group_formula_spec_fails_cap_without_group_materialization():
    spec = DesignSpecV1(
        plan=DISCONNECTED_GROUPS,
        n_persons=10**9,
        n_raters=10**9,
        n_groups=10**9,
    )
    with pytest.raises(AssignmentLimitError) as exc_info:
        compile_rating_design(spec, max_rating_sessions=100_000)
    assert exc_info.value.required_rating_sessions == 10**9


def test_bundle_json_round_trip_is_exact_and_strict():
    spec = DesignSpecV1(
        plan=NESTED_BRIDGE,
        n_persons=6,
        n_raters=4,
        n_artifacts_per_person=2,
        n_bridge_artifacts=3,
        bridge_extra_raters_per_artifact=1,
        assignment_seed=99,
    )
    bundle = compile_rating_design(spec)
    payload = assignment_bundle_to_dict(bundle)
    restored = normalize_assignment_bundle(
        json.loads(json.dumps(payload, allow_nan=False))
    )
    assert restored == bundle
    assert isinstance(restored.assignments, tuple)

    missing = dict(payload)
    missing.pop("assignment_algorithm")
    with pytest.raises(AssignmentValidationError, match="missing fields"):
        normalize_assignment_bundle(missing)
    with pytest.raises(AssignmentValidationError, match="unknown fields"):
        normalize_assignment_bundle({**payload, "future_field": 1})
    bad_row = dict(payload)
    bad_row["assignments"] = [dict(row) for row in payload["assignments"]]
    bad_row["assignments"][0]["future_field"] = 1
    with pytest.raises(AssignmentValidationError, match="unknown fields"):
        normalize_assignment_bundle(bad_row)


def test_normalizer_enforces_its_own_session_limit_before_row_construction():
    bundle = compile_rating_design(
        DesignSpecV1(
            plan=FULLY_CROSSED,
            n_persons=3,
            n_raters=2,
            raters_per_artifact=2,
        )
    )
    payload = bundle.to_dict()
    with pytest.raises(AssignmentLimitError) as exc_info:
        normalize_assignment_bundle(payload, max_rating_sessions=5)
    assert exc_info.value.required_rating_sessions == 6


def test_normalizer_rejects_assignment_bundle_subclasses_that_bypass_validation():
    bundle = compile_rating_design(GOLDEN_CASES[0][0])
    forged = _ForgedAssignmentBundle(
        schema_version=bundle.schema_version,
        design_spec=bundle.design_spec,
        design_spec_fingerprint=bundle.design_spec_fingerprint,
        assignment_algorithm=bundle.assignment_algorithm,
        assignment_fingerprint=bundle.assignment_fingerprint,
        workload=bundle.workload,
        assignments=bundle.assignments,
    )

    with pytest.raises(AssignmentValidationError, match="subclasses"):
        normalize_assignment_bundle(forged)


def test_assignment_fingerprint_is_order_independent_and_content_sensitive():
    bundle = compile_rating_design(GOLDEN_CASES[0][0])
    assert assignment_fingerprint(reversed(bundle.assignments)) == (
        bundle.assignment_fingerprint
    )
    changed = replace(bundle.assignments[0], group_id="G002")
    assert assignment_fingerprint((changed, *bundle.assignments[1:])) != (
        bundle.assignment_fingerprint
    )
    assert len(assignment_fingerprint(bundle.assignments, length=64)) == 64
    for invalid in (True, 7, 65, 8.0):
        with pytest.raises(ValueError, match="between 8 and 64"):
            assignment_fingerprint(bundle.assignments, length=invalid)


def test_bundle_rejects_duplicate_wrong_artifact_and_inconsistent_fingerprints():
    bundle = compile_rating_design(GOLDEN_CASES[0][0])
    duplicated = tuple(
        sorted(
            (bundle.assignments[0], *bundle.assignments[:-1]),
            key=assignment_module._assignment_sort_key,
        )
    )
    with pytest.raises(AssignmentValidationError, match="duplicate"):
        replace(
            bundle,
            assignments=duplicated,
            assignment_fingerprint=assignment_fingerprint(duplicated),
        )

    wrong_artifact_row = replace(
        bundle.assignments[0],
        artifact_id="P000001-A999",
    )
    wrong_artifact = tuple(
        sorted(
            (wrong_artifact_row, *bundle.assignments[1:]),
            key=assignment_module._assignment_sort_key,
        )
    )
    with pytest.raises(AssignmentValidationError, match="study artifacts"):
        replace(
            bundle,
            assignments=wrong_artifact,
            assignment_fingerprint=assignment_fingerprint(wrong_artifact),
        )
    with pytest.raises(AssignmentValidationError, match="does not match"):
        replace(bundle, assignment_fingerprint="0" * 16)


def test_bundle_and_json_normalizer_reject_a_relabelled_noncanonical_schedule():
    bundle = compile_rating_design(GOLDEN_CASES[1][0])
    forged = list(bundle.assignments)
    forged[1] = replace(forged[1], rater_id="R000002")
    forged[2] = replace(forged[2], rater_id="R000004")
    forged_rows = tuple(sorted(forged, key=assignment_module._assignment_sort_key))
    forged_fingerprint = assignment_fingerprint(forged_rows)
    with pytest.raises(AssignmentValidationError, match="declared assignment algorithm"):
        replace(
            bundle,
            assignments=forged_rows,
            assignment_fingerprint=forged_fingerprint,
        )

    payload = bundle.to_dict()
    payload["assignments"] = [row.to_dict() for row in forged_rows]
    payload["assignment_fingerprint"] = forged_fingerprint
    with pytest.raises(AssignmentValidationError, match="declared assignment algorithm"):
        normalize_assignment_bundle(payload)


def test_nested_bundle_rejects_a_bridge_to_an_undeclared_artifact():
    bundle = compile_rating_design(GOLDEN_CASES[2][0])
    forged = tuple(
        sorted(
            (
                *bundle.assignments[:2],
                replace(bundle.assignments[2], artifact_id="P000001-A999"),
                *bundle.assignments[3:],
            ),
            key=assignment_module._assignment_sort_key,
        )
    )
    with pytest.raises(AssignmentValidationError, match="declared primary artifacts"):
        replace(
            bundle,
            assignments=forged,
            assignment_fingerprint=assignment_fingerprint(forged),
        )


@pytest.mark.parametrize(
    "kwargs",
    [
        {"person_id": "bad"},
        {"person_id": "P000000", "artifact_id": "P000000-A001"},
        {
            "person_id": "CA000000",
            "artifact_id": "CA000000-A001",
            "person_role": PERSON_ROLE_COMMON_ANCHOR,
            "session_type": SESSION_COMMON_ANCHOR,
        },
        {"artifact_id": "P000002-A001"},
        {"artifact_id": "P000001-A000"},
        {"rater_id": "R1"},
        {"rater_id": "R000000"},
        {"group_id": "group-1"},
        {"group_id": "G000"},
        {"person_role": "anchor-ish"},
        {"session_type": "repair"},
        {
            "person_role": PERSON_ROLE_COMMON_ANCHOR,
            "session_type": SESSION_PRIMARY,
            "person_id": "CA000001",
            "artifact_id": "CA000001-A001",
        },
    ],
)
def test_rating_assignment_rejects_ambiguous_ids_roles_and_types(kwargs):
    base = dict(
        person_id="P000001",
        artifact_id="P000001-A001",
        rater_id="R000001",
        group_id="G001",
        person_role=PERSON_ROLE_STUDY,
        session_type=SESSION_PRIMARY,
    )
    with pytest.raises(AssignmentValidationError):
        RatingAssignmentV1(**{**base, **kwargs})


@pytest.mark.parametrize(
    "kwargs",
    [
        {"person_id": 1},
        {"artifact_id": 1},
        {"rater_id": 1},
        {"group_id": 1},
    ],
)
def test_rating_assignment_rejects_non_string_identifiers_with_domain_error(kwargs):
    base = dict(
        person_id="P000001",
        artifact_id="P000001-A001",
        rater_id="R000001",
        group_id="G001",
        person_role=PERSON_ROLE_STUDY,
        session_type=SESSION_PRIMARY,
    )
    with pytest.raises(AssignmentValidationError):
        RatingAssignmentV1(**{**base, **kwargs})


def test_rating_assignment_handles_a_very_long_numeric_suffix_without_int():
    person_id = "P" + "0" * 5_000 + "1"
    row = RatingAssignmentV1(
        person_id=person_id,
        artifact_id=f"{person_id}-A001",
        rater_id="R" + "0" * 5_000 + "1",
        group_id="G" + "0" * 5_000 + "1",
        person_role=PERSON_ROLE_STUDY,
        session_type=SESSION_PRIMARY,
    )
    assert row.person_id == person_id


def test_assignment_rows_and_bundle_are_immutable():
    bundle = compile_rating_design(GOLDEN_CASES[0][0])
    with pytest.raises(FrozenInstanceError):
        bundle.assignments[0].rater_id = "R000002"
    with pytest.raises(FrozenInstanceError):
        bundle.assignment_fingerprint = "0" * 16
    assert bundle.cache_identity == (
        bundle.schema_version,
        bundle.design_spec_fingerprint,
        bundle.assignment_fingerprint,
    )


def test_target_interactive_preview_compiles_well_below_two_seconds():
    spec = DesignSpecV1(
        plan=FULLY_CROSSED,
        n_persons=200,
        n_raters=6,
        raters_per_artifact=6,
    )
    started = time.perf_counter()
    bundle = compile_rating_design(spec)
    elapsed = time.perf_counter() - started
    assert len(bundle.assignments) == 1_200
    assert elapsed < 2.0


def test_worst_shape_extra_rater_selection_is_linear_enough_for_preview():
    n_raters = 8_000
    spec = DesignSpecV1(
        plan=NESTED_BRIDGE,
        n_persons=1,
        n_raters=n_raters,
        n_bridge_artifacts=1,
        bridge_extra_raters_per_artifact=n_raters - 1,
        assignment_seed=0,
    )
    started = time.perf_counter()
    bundle = compile_rating_design(spec)
    elapsed = time.perf_counter() - started
    assert len(bundle.assignments) == n_raters
    assert elapsed < 2.0
