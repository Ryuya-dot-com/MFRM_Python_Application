"""Prospective rating-design identity and cost-currency contracts."""

from __future__ import annotations

from dataclasses import FrozenInstanceError, replace
import json
from pathlib import Path
import subprocess
import sys

import pytest

import mfrm_app.simulation as simulation
from mfrm_app.simulation import (
    ASSIGNMENT_ALGORITHM_V1,
    DISCONNECTED_GROUPS,
    FULLY_CROSSED,
    NESTED_BRIDGE,
    SPIRAL,
    CostPolicyV1,
    DesignSpecV1,
    DesignValidationError,
    WeightedDesignCostV1,
    apply_cost_policy,
    compute_design_workload,
    cost_policy_fingerprint,
    cost_policy_to_dict,
    design_spec_fingerprint,
    design_spec_to_dict,
    normalize_design_spec,
    normalize_cost_policy,
)


def test_simulation_contract_import_stays_streamlit_and_pandas_free():
    project_root = Path(__file__).resolve().parents[1]
    code = "\n".join(
        (
            "import sys",
            "import mfrm_app.simulation",
            "for name in ('streamlit', 'pandas', 'numpy', 'scipy', "
            "'networkx', 'plotly'):",
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


def test_simulation_package_exports_are_unique_and_resolvable():
    assert len(simulation.__all__) == len(set(simulation.__all__))
    for name in simulation.__all__:
        assert getattr(simulation, name) is not None


@pytest.mark.parametrize(
    ("spec", "expected"),
    [
        (
            DesignSpecV1(
                plan=FULLY_CROSSED,
                n_persons=20,
                n_raters=3,
                raters_per_artifact=3,
                score_records_per_session=4,
            ),
            (60, 0, 0, 60, 240),
        ),
        (
            DesignSpecV1(
                plan=SPIRAL,
                n_persons=30,
                n_raters=4,
                n_artifacts_per_person=2,
                raters_per_artifact=2,
                score_records_per_session=3,
            ),
            (120, 0, 0, 120, 360),
        ),
        (
            DesignSpecV1(
                plan=NESTED_BRIDGE,
                n_persons=20,
                n_raters=4,
                score_records_per_session=4,
                n_bridge_artifacts=5,
                bridge_extra_raters_per_artifact=1,
            ),
            (20, 5, 0, 25, 100),
        ),
        (
            DesignSpecV1(
                plan=DISCONNECTED_GROUPS,
                n_persons=20,
                n_raters=4,
                n_groups=2,
                score_records_per_session=4,
                n_common_anchor_artifacts=5,
                common_anchor_raters_per_group=1,
                calibration_events=3,
            ),
            (20, 0, 10, 30, 120),
        ),
    ],
)
def test_four_design_plans_have_auditable_workload_arithmetic(spec, expected):
    workload = compute_design_workload(spec)
    assert (
        workload.primary_rating_sessions,
        workload.bridge_rating_sessions,
        workload.common_anchor_rating_sessions,
        workload.total_rating_sessions,
        workload.total_score_records,
    ) == expected
    assert workload.total_rating_sessions == (
        workload.primary_rating_sessions
        + workload.bridge_rating_sessions
        + workload.common_anchor_rating_sessions
    )
    assert workload.total_score_records == (
        workload.primary_score_records
        + workload.bridge_score_records
        + workload.common_anchor_score_records
    )


def test_default_cost_currency_is_rating_sessions_not_long_rows():
    base = DesignSpecV1(
        plan=SPIRAL,
        n_persons=30,
        n_raters=4,
        raters_per_artifact=2,
        score_records_per_session=1,
    )
    expanded_rubric = replace(base, score_records_per_session=8)
    base_workload = compute_design_workload(base)
    expanded_workload = compute_design_workload(expanded_rubric)

    assert base_workload.total_rating_sessions == expanded_workload.total_rating_sessions
    assert expanded_workload.total_score_records == 8 * base_workload.total_score_records
    assert apply_cost_policy(base_workload).total_cost == pytest.approx(60.0)
    assert apply_cost_policy(expanded_workload).total_cost == pytest.approx(60.0)


def test_score_rows_and_calibration_are_counted_only_by_explicit_cost_policy():
    spec = DesignSpecV1(
        plan=DISCONNECTED_GROUPS,
        n_persons=20,
        n_raters=4,
        n_groups=2,
        score_records_per_session=4,
        n_common_anchor_artifacts=5,
        common_anchor_raters_per_group=1,
        calibration_events=3,
    )
    workload = compute_design_workload(spec)
    weighted = apply_cost_policy(
        workload,
        CostPolicyV1(
            rating_session_unit_cost=2.0,
            score_record_unit_cost=0.25,
            calibration_event_unit_cost=10.0,
            cost_unit="  minutes  ",
        ),
    )
    assert weighted.rating_session_cost == pytest.approx(60.0)
    assert weighted.score_record_cost == pytest.approx(30.0)
    assert weighted.calibration_event_cost == pytest.approx(30.0)
    assert weighted.total_cost == pytest.approx(120.0)
    assert weighted.cost_unit == "minutes"


def test_calibration_events_are_plan_independent_and_create_no_rating_sessions():
    without_calibration = DesignSpecV1(
        plan=DISCONNECTED_GROUPS,
        n_persons=20,
        n_raters=4,
        n_groups=2,
    )
    with_calibration = replace(without_calibration, calibration_events=7)
    base = compute_design_workload(without_calibration)
    calibrated = compute_design_workload(with_calibration)
    assert calibrated.total_rating_sessions == base.total_rating_sessions
    assert calibrated.total_score_records == base.total_score_records
    assert calibrated.calibration_events == 7

    fully_crossed = DesignSpecV1(
        plan=FULLY_CROSSED,
        n_persons=20,
        n_raters=2,
        raters_per_artifact=2,
        calibration_events=4,
    )
    assert compute_design_workload(fully_crossed).calibration_events == 4


def test_design_spec_json_round_trip_and_fingerprint_are_stable():
    spec = DesignSpecV1(
        plan=NESTED_BRIDGE,
        n_persons=50,
        n_raters=6,
        n_artifacts_per_person=2,
        score_records_per_session=3,
        n_bridge_artifacts=12,
        bridge_extra_raters_per_artifact=1,
        assignment_seed=42,
    )
    payload = design_spec_to_dict(spec)
    assert payload["assignment_algorithm"] == ASSIGNMENT_ALGORITHM_V1
    encoded = json.dumps(payload, allow_nan=False, sort_keys=True)
    restored = normalize_design_spec(json.loads(encoded))
    assert restored == spec
    assert design_spec_fingerprint(restored) == design_spec_fingerprint(spec)
    assert design_spec_fingerprint(replace(spec, assignment_seed=43)) != (
        design_spec_fingerprint(spec)
    )


def test_design_spec_is_immutable():
    spec = DesignSpecV1(
        plan=FULLY_CROSSED,
        n_persons=20,
        n_raters=2,
        raters_per_artifact=2,
    )
    with pytest.raises(FrozenInstanceError):
        spec.n_persons = 30


def test_saved_mapping_requires_exact_versioned_shape():
    spec = DesignSpecV1(
        plan=FULLY_CROSSED,
        n_persons=20,
        n_raters=2,
        raters_per_artifact=2,
    )
    payload = design_spec_to_dict(spec)
    missing = dict(payload)
    missing.pop("assignment_seed")
    with pytest.raises(DesignValidationError, match="missing fields"):
        normalize_design_spec(missing)

    with pytest.raises(DesignValidationError, match="unknown fields"):
        normalize_design_spec({**payload, "future_field": 1})
    with pytest.raises(DesignValidationError, match="field names must be strings"):
        normalize_design_spec({**payload, 1: "not-a-field"})
    with pytest.raises(DesignValidationError, match="Unsupported design spec version"):
        normalize_design_spec({**payload, "schema_version": "future"})


@pytest.mark.parametrize(
    "kwargs",
    [
        {"n_persons": True},
        {"n_persons": "20"},
        {"n_persons": 0},
        {"assignment_seed": -1},
        {"assignment_seed": 2**63},
    ],
)
def test_non_json_integer_or_out_of_range_design_values_are_rejected(kwargs):
    base = dict(
        plan=FULLY_CROSSED,
        n_persons=20,
        n_raters=2,
        raters_per_artifact=2,
    )
    with pytest.raises(DesignValidationError):
        DesignSpecV1(**{**base, **kwargs})


def test_plan_specific_fields_cannot_silently_change_another_plan():
    with pytest.raises(DesignValidationError, match="bridge-specific"):
        DesignSpecV1(
            plan=SPIRAL,
            n_persons=20,
            n_raters=3,
            raters_per_artifact=2,
            n_bridge_artifacts=1,
            bridge_extra_raters_per_artifact=1,
        )
    with pytest.raises(DesignValidationError, match="anchor-specific"):
        DesignSpecV1(
            plan=NESTED_BRIDGE,
            n_persons=20,
            n_raters=3,
            n_common_anchor_artifacts=1,
            common_anchor_raters_per_group=1,
        )
    with pytest.raises(DesignValidationError, match="observational bridge"):
        DesignSpecV1(
            plan=DISCONNECTED_GROUPS,
            n_persons=20,
            n_raters=4,
            n_groups=2,
            n_bridge_artifacts=1,
            bridge_extra_raters_per_artifact=1,
        )


def test_plan_bounds_prevent_ambiguous_or_impossible_schedules():
    with pytest.raises(DesignValidationError, match="fully_crossed requires"):
        DesignSpecV1(
            plan=FULLY_CROSSED,
            n_persons=20,
            n_raters=3,
            raters_per_artifact=2,
        )
    with pytest.raises(DesignValidationError, match="too few primary rating slots"):
        DesignSpecV1(
            plan=SPIRAL,
            n_persons=1,
            n_raters=4,
            raters_per_artifact=2,
        )
    # Coverage is determined by available rating slots, not artifact count.
    DesignSpecV1(
        plan=SPIRAL,
        n_persons=2,
        n_raters=4,
        raters_per_artifact=2,
    )
    with pytest.raises(DesignValidationError, match="primary artifact count"):
        DesignSpecV1(
            plan=NESTED_BRIDGE,
            n_persons=20,
            n_raters=3,
            n_bridge_artifacts=21,
            bridge_extra_raters_per_artifact=1,
        )
    with pytest.raises(DesignValidationError, match="smallest rater group"):
        DesignSpecV1(
            plan=DISCONNECTED_GROUPS,
            n_persons=20,
            n_raters=5,
            n_groups=2,
            n_common_anchor_artifacts=1,
            common_anchor_raters_per_group=3,
        )


def test_fingerprint_length_is_bounded_by_sha256_width():
    spec = DesignSpecV1(
        plan=FULLY_CROSSED,
        n_persons=20,
        n_raters=2,
        raters_per_artifact=2,
    )
    assert len(design_spec_fingerprint(spec, length=64)) == 64
    for invalid in (True, 7, 65, 8.0):
        with pytest.raises(ValueError, match="between 8 and 64"):
            design_spec_fingerprint(spec, length=invalid)


def test_huge_workload_uses_count_arithmetic_without_materializing_rows():
    spec = DesignSpecV1(
        plan=FULLY_CROSSED,
        n_persons=10**9,
        n_raters=6,
        n_artifacts_per_person=10,
        score_records_per_session=20,
        raters_per_artifact=6,
    )
    workload = compute_design_workload(spec)
    assert workload.total_rating_sessions == 60_000_000_000
    assert workload.total_score_records == 1_200_000_000_000


def test_declared_raters_are_active_not_an_unused_available_pool():
    with pytest.raises(DesignValidationError, match="activate every rater"):
        DesignSpecV1(
            plan=NESTED_BRIDGE,
            n_persons=1,
            n_raters=100,
        )
    with pytest.raises(DesignValidationError, match="activate every rater"):
        DesignSpecV1(
            plan=DISCONNECTED_GROUPS,
            n_persons=2,
            n_raters=10,
            n_groups=2,
        )
    # Common-anchor sessions may supply the otherwise missing active slots.
    DesignSpecV1(
        plan=DISCONNECTED_GROUPS,
        n_persons=2,
        n_raters=10,
        n_groups=2,
        n_common_anchor_artifacts=2,
        common_anchor_raters_per_group=2,
    )


def test_cost_policy_round_trip_fingerprint_and_provenance_are_stable():
    policy = CostPolicyV1(
        rating_session_unit_cost=2,
        score_record_unit_cost=0.25,
        calibration_event_unit_cost=10,
        cost_unit="  staff-minutes ",
    )
    payload = cost_policy_to_dict(policy)
    assert payload["rating_session_unit_cost"] == 2.0
    assert payload["cost_unit"] == "staff-minutes"
    restored = normalize_cost_policy(json.loads(json.dumps(payload)))
    assert restored == policy
    assert cost_policy_fingerprint(restored) == cost_policy_fingerprint(policy)

    spec = DesignSpecV1(
        plan=FULLY_CROSSED,
        n_persons=20,
        n_raters=2,
        raters_per_artifact=2,
    )
    weighted = apply_cost_policy(compute_design_workload(spec), policy)
    assert isinstance(weighted, WeightedDesignCostV1)
    assert weighted.cost_policy_fingerprint == cost_policy_fingerprint(policy)
    assert weighted.cost_policy_schema_version == policy.schema_version
    assert weighted.rating_session_unit_cost == 2.0
    assert weighted.score_record_unit_cost == 0.25
    assert weighted.calibration_event_unit_cost == 10.0
    assert weighted.total_rating_sessions == 40
    assert weighted.total_score_records == 40
    assert weighted.calibration_events == 0

    workload_round_trip = type(compute_design_workload(spec))(
        **json.loads(json.dumps(compute_design_workload(spec).to_dict()))
    )
    weighted_round_trip = type(weighted)(
        **json.loads(json.dumps(weighted.to_dict()))
    )
    assert workload_round_trip == compute_design_workload(spec)
    assert weighted_round_trip == weighted


def test_saved_cost_policy_requires_exact_versioned_shape():
    payload = cost_policy_to_dict(CostPolicyV1())
    missing = dict(payload)
    missing.pop("cost_unit")
    with pytest.raises(ValueError, match="missing fields"):
        normalize_cost_policy(missing)
    with pytest.raises(ValueError, match="unknown fields"):
        normalize_cost_policy({**payload, "future_field": 1})
    with pytest.raises(ValueError, match="Unsupported cost policy version"):
        normalize_cost_policy({**payload, "schema_version": "future"})


def test_cost_policy_fingerprint_has_bounded_sha256_prefix():
    policy = CostPolicyV1()
    assert len(cost_policy_fingerprint(policy, length=8)) == 8
    assert len(cost_policy_fingerprint(policy, length=64)) == 64
    for invalid in (True, 7, 65, 8.0):
        with pytest.raises(ValueError, match="between 8 and 64"):
            cost_policy_fingerprint(policy, length=invalid)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"rating_session_unit_cost": 0},
        {"rating_session_unit_cost": True},
        {"score_record_unit_cost": -0.1},
        {"score_record_unit_cost": float("nan")},
        {"calibration_event_unit_cost": float("inf")},
        {"calibration_event_unit_cost": 10**400},
        {"cost_unit": "  "},
    ],
)
def test_cost_policy_rejects_nonfinite_negative_or_ambiguous_values(kwargs):
    with pytest.raises(ValueError):
        CostPolicyV1(**kwargs)


def test_cost_policy_can_explicitly_use_score_rows_as_its_only_currency():
    policy = CostPolicyV1(
        rating_session_unit_cost=0,
        score_record_unit_cost=1,
    )
    assert policy.rating_session_unit_cost == 0.0
    assert policy.score_record_unit_cost == 1.0


def test_equivalent_positive_and_negative_zero_policy_values_are_canonical():
    positive = CostPolicyV1(score_record_unit_cost=0.0)
    negative = CostPolicyV1(score_record_unit_cost=-0.0)
    assert positive == negative
    assert cost_policy_fingerprint(positive) == cost_policy_fingerprint(negative)


def test_derived_workload_and_weighted_cost_reject_invalid_direct_construction():
    spec = DesignSpecV1(
        plan=FULLY_CROSSED,
        n_persons=20,
        n_raters=2,
        raters_per_artifact=2,
    )
    workload = compute_design_workload(spec)
    with pytest.raises(ValueError, match="component sum"):
        replace(workload, total_rating_sessions=999)
    with pytest.raises(ValueError, match="times score_records_per_session"):
        replace(workload, primary_score_records=999)
    with pytest.raises(ValueError, match="nonnegative integer"):
        replace(workload, calibration_events=-1)
    with pytest.raises(ValueError, match="must be positive"):
        replace(
            workload,
            primary_rating_sessions=0,
            total_rating_sessions=0,
            primary_score_records=0,
            total_score_records=0,
        )
    with pytest.raises(ValueError, match="Unsupported design workload version"):
        replace(workload, schema_version="future")

    weighted = apply_cost_policy(workload)
    with pytest.raises(ValueError, match="finite"):
        replace(weighted, total_cost=float("nan"))
    with pytest.raises(ValueError, match="component sum"):
        replace(weighted, total_cost=999.0)
    with pytest.raises(ValueError, match="count times unit cost"):
        replace(weighted, total_rating_sessions=999)
    with pytest.raises(ValueError, match="fingerprint"):
        replace(weighted, cost_policy_fingerprint="not-a-fingerprint")
    other_policy = CostPolicyV1(rating_session_unit_cost=2.0)
    with pytest.raises(ValueError, match="does not match the embedded"):
        replace(
            weighted,
            cost_policy_fingerprint=cost_policy_fingerprint(other_policy),
        )
    with pytest.raises(ValueError, match="Unsupported weighted design cost version"):
        replace(weighted, schema_version="future")


def test_weighted_cost_rejects_large_relative_error_at_tiny_magnitude():
    spec = DesignSpecV1(
        plan=FULLY_CROSSED,
        n_persons=1,
        n_raters=1,
        raters_per_artifact=1,
    )
    policy = CostPolicyV1(rating_session_unit_cost=1e-300)
    weighted = apply_cost_policy(compute_design_workload(spec), policy)
    with pytest.raises(ValueError, match="count times unit cost"):
        replace(weighted, rating_session_cost=1e-13, total_cost=1e-13)


def test_weighted_cost_fails_clearly_when_float_range_is_exceeded():
    spec = DesignSpecV1(
        plan=FULLY_CROSSED,
        n_persons=10**400,
        n_raters=2,
        raters_per_artifact=2,
    )
    workload = compute_design_workload(spec)
    with pytest.raises(ValueError, match="finite weighted-cost range"):
        apply_cost_policy(workload)
