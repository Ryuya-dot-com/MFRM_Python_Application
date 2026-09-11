"""Identifiability preflight for the existing K-1 structural repair.

The current repair overlay is non-empty only when the anchor-augmented
Person-Rater graph is disconnected.  It therefore represents a transition
from no planned-incidence common-scale connection to a connected graph, not a
connected-to-connected precision enhancement.  This module makes that
distinction machine-enforced before any simulation or estimator can run.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, fields
import hashlib
import json
import math

from .conditions import (
    ESTIMATOR_JMLE,
    ESTIMATOR_MML,
    MAX_SIMULATION_ARTIFACT_INDEX,
    MAX_SIMULATION_CRITERION_INDEX,
    MAX_SIMULATION_PERSON_INDEX,
    MAX_SIMULATION_RATER_INDEX,
    SimulationScenarioSpecV1,
    normalize_simulation_scenario_spec,
    simulation_scenario_spec_fingerprint,
)
from .cost import CostPolicyV1, cost_policy_fingerprint, normalize_cost_policy
from .planning import (
    FORMULA_PLAN_SUMMARY_VERSION,
    PLANNING_ALGORITHM_V1,
    STRUCTURE_EVALUATED,
    StructuralPlanSummaryV1,
)
from .topology import PRECISION_NOT_EVALUATED, REPAIR_PROPOSED


STRUCTURAL_TRANSITION_PREFLIGHT_VERSION = (
    "mfrm_structural_transition_preflight_v1"
)
STRUCTURAL_TRANSITION_ALGORITHM_V1 = "structural_estimability_transition_v1"
COMPARISON_ESTIMABILITY_TRANSITION = "estimability_transition"
PRECISION_CHANGE_NOT_COMPARABLE = "not_comparable_estimability_transition"
EXECUTION_NOT_STARTED = "not_executed"

BASE_JMLE_BLOCKED = "blocked_nonidentified_topology"
BASE_MML_PRIOR_CONDITIONAL = "eligible_prior_conditional_only"
REPAIRED_PENDING_PREFLIGHT = "eligible_pending_full_execution_preflight"

IDENTIFICATION_NONE = "none"
IDENTIFICATION_FIXED_POPULATION_PRIOR = "fixed_population_prior"
IDENTIFICATION_PLANNED_PERSON_RATER = "planned_person_rater_incidence_connectivity"
IDENTIFICATION_PLANNED_PERSON_RATER_PLUS_FIXED_PRIOR = (
    "planned_person_rater_incidence_connectivity_plus_fixed_population_prior"
)

CLAIM_NOT_ESTIMABLE_COMMON_SCALE = "not_estimable_common_scale"
CLAIM_PRIOR_CONDITIONAL_EXPLORATORY = "prior_conditional_exploratory"
CLAIM_STRUCTURAL_ONLY_PENDING_FIT = "structural_connectivity_only_pending_fit"
CLAIM_PLANNED_CONNECTIVITY_PLUS_PRIOR_PENDING_FIT = (
    "planned_connectivity_plus_prior_pending_fit"
)
MAX_STRUCTURAL_TRANSITION_COUNT = 2**63 - 1


class StructuralTransitionValidationError(ValueError):
    """Raised when a structural transition is misrepresented as precision."""


def _strict_mapping(value: Mapping) -> dict:
    expected = tuple(field.name for field in fields(StructuralTransitionPreflightV1))
    expected_set = frozenset(expected)
    supplied = set(value)
    if any(not isinstance(name, str) for name in supplied):
        raise StructuralTransitionValidationError(
            "Saved transition preflight field names must be strings"
        )
    missing = expected_set - supplied
    unknown = supplied - expected_set
    if missing:
        raise StructuralTransitionValidationError(
            f"Saved transition preflight is missing fields: {sorted(missing)}"
        )
    if unknown:
        raise StructuralTransitionValidationError(
            f"Saved transition preflight contains unknown fields: {sorted(unknown)}"
        )
    return {name: value[name] for name in expected}


def _fingerprint(payload: Mapping) -> str:
    encoded = json.dumps(
        dict(payload),
        allow_nan=False,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()[:16]


def _require_fingerprint(name: str, value: object) -> str:
    if (
        not isinstance(value, str)
        or len(value) != 16
        or value != value.lower()
        or any(character not in "0123456789abcdef" for character in value)
    ):
        raise StructuralTransitionValidationError(
            f"{name} must be a 16-character lowercase hex fingerprint"
        )
    return value


def _require_nonnegative_int(name: str, value: object) -> int:
    if (
        isinstance(value, bool)
        or not isinstance(value, int)
        or value < 0
        or value > MAX_STRUCTURAL_TRANSITION_COUNT
    ):
        raise StructuralTransitionValidationError(
            f"{name} must be a nonnegative integer <= "
            f"{MAX_STRUCTURAL_TRANSITION_COUNT}"
        )
    return value


def _require_cost(name: str, value: object) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise StructuralTransitionValidationError(f"{name} must be numeric")
    try:
        result = float(value)
    except (OverflowError, ValueError) as exc:
        raise StructuralTransitionValidationError(
            f"{name} must be finite and nonnegative"
        ) from exc
    if not math.isfinite(result) or result < 0:
        raise StructuralTransitionValidationError(
            f"{name} must be finite and nonnegative"
        )
    return result


def _transition_id(
    *,
    design_spec_fingerprint: str,
    assignment_fingerprint: str,
    repair_fingerprint: str,
    scenario_fingerprint: str,
    components_before: int,
    score_records_per_session: int,
    base_rating_sessions: int,
    repaired_rating_sessions: int,
    base_score_records: int,
    repaired_score_records: int,
) -> str:
    return _fingerprint({
        "algorithm": STRUCTURAL_TRANSITION_ALGORITHM_V1,
        "comparison_kind": COMPARISON_ESTIMABILITY_TRANSITION,
        "design_spec_fingerprint": design_spec_fingerprint,
        "assignment_fingerprint": assignment_fingerprint,
        "repair_fingerprint": repair_fingerprint,
        "scenario_fingerprint": scenario_fingerprint,
        "components_before": components_before,
        "score_records_per_session": score_records_per_session,
        "base_rating_sessions": base_rating_sessions,
        "repaired_rating_sessions": repaired_rating_sessions,
        "base_score_records": base_score_records,
        "repaired_score_records": repaired_score_records,
    })


def _economic_context_id(
    *,
    transition_id: str,
    formula_evaluation_id: str,
    cost_policy_fingerprint: str,
    base_calibration_events: int,
    repaired_calibration_events: int,
    base_total_cost: float,
    repaired_total_cost: float,
) -> str:
    return _fingerprint({
        "transition_id": transition_id,
        "formula_evaluation_id": formula_evaluation_id,
        "cost_policy_fingerprint": cost_policy_fingerprint,
        "base_calibration_events": base_calibration_events,
        "repaired_calibration_events": repaired_calibration_events,
        "base_total_cost": float(base_total_cost),
        "repaired_total_cost": float(repaired_total_cost),
    })


def _weighted_total_cost(
    *,
    rating_sessions: int,
    score_records: int,
    calibration_events: int,
    policy: CostPolicyV1,
) -> float:
    components = (
        rating_sessions * policy.rating_session_unit_cost,
        score_records * policy.score_record_unit_cost,
        calibration_events * policy.calibration_event_unit_cost,
    )
    total = float(components[0] + components[1] + components[2])
    if not all(math.isfinite(value) for value in (*components, total)):
        raise StructuralTransitionValidationError(
            "workload and cost policy exceed the finite weighted-cost range"
        )
    return total


def _formula_evaluation_id(
    *,
    design_spec_fingerprint: str,
    cost_policy_fingerprint: str,
) -> str:
    return _fingerprint({
        "schema_version": FORMULA_PLAN_SUMMARY_VERSION,
        "planning_algorithm": PLANNING_ALGORITHM_V1,
        "design_spec_fingerprint": design_spec_fingerprint,
        "cost_policy_fingerprint": cost_policy_fingerprint,
    })


def _validate_simulation_id_capacity(structural: StructuralPlanSummaryV1) -> None:
    """Reject plans whose compiler IDs exceed the fixed-width v1 RNG contract."""
    spec = structural.formula_summary.design_spec
    limits = (
        ("n_persons", spec.n_persons, MAX_SIMULATION_PERSON_INDEX),
        ("n_raters", spec.n_raters, MAX_SIMULATION_RATER_INDEX),
        (
            "n_common_anchor_artifacts",
            spec.n_common_anchor_artifacts,
            MAX_SIMULATION_PERSON_INDEX,
        ),
        (
            "n_artifacts_per_person",
            spec.n_artifacts_per_person,
            MAX_SIMULATION_ARTIFACT_INDEX,
        ),
        (
            "score_records_per_session",
            spec.score_records_per_session,
            MAX_SIMULATION_CRITERION_INDEX,
        ),
    )
    exceeded = tuple(
        f"{name}={value} > {maximum}"
        for name, value, maximum in limits
        if value > maximum
    )
    if exceeded:
        raise StructuralTransitionValidationError(
            "design exceeds the v1 simulation ID capacity: " + ", ".join(exceeded)
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class StructuralTransitionPreflightV1:
    """Compact planned-incidence claim boundary for a non-zero K-1 repair."""

    schema_version: str
    algorithm: str
    transition_id: str
    economic_context_id: str
    comparison_kind: str
    scenario: SimulationScenarioSpecV1
    scenario_fingerprint: str
    formula_evaluation_id: str
    design_spec_fingerprint: str
    cost_policy: CostPolicyV1
    cost_policy_fingerprint: str
    assignment_fingerprint: str
    repair_fingerprint: str
    components_before: int
    projected_components_after: int
    score_records_per_session: int
    base_rating_sessions: int
    repaired_rating_sessions: int
    base_score_records: int
    repaired_score_records: int
    base_calibration_events: int
    repaired_calibration_events: int
    base_total_cost: float
    repaired_total_cost: float
    planned_person_rater_graph_connected_before: bool
    planned_person_rater_graph_connected_after: bool
    base_fit_preflight_status: str
    base_identification_basis: str
    base_claim_tier: str
    repaired_fit_preflight_status: str
    repaired_identification_basis: str
    repaired_claim_tier: str
    precision_change_status: str
    primary_claim_eligible: bool
    precision_enhancement_overlay_required: bool
    paired_precision_evaluation_applicable: bool
    execution_status: str
    precision_evaluation_status: str
    expected_se_improvement_percent: None
    recovery_improvement_percent: None
    classification_improvement_percent: None
    sample_size_recommendation: None

    def __post_init__(self) -> None:
        if self.schema_version != STRUCTURAL_TRANSITION_PREFLIGHT_VERSION:
            raise StructuralTransitionValidationError(
                f"Unsupported transition preflight version: {self.schema_version!r}"
            )
        if self.algorithm != STRUCTURAL_TRANSITION_ALGORITHM_V1:
            raise StructuralTransitionValidationError(
                f"Unsupported transition algorithm: {self.algorithm!r}"
            )
        for name in (
            "transition_id",
            "economic_context_id",
            "scenario_fingerprint",
            "formula_evaluation_id",
            "design_spec_fingerprint",
            "cost_policy_fingerprint",
            "assignment_fingerprint",
            "repair_fingerprint",
        ):
            _require_fingerprint(name, getattr(self, name))
        if self.comparison_kind != COMPARISON_ESTIMABILITY_TRANSITION:
            raise StructuralTransitionValidationError(
                "comparison_kind must remain estimability_transition"
            )
        if type(self.scenario) is not SimulationScenarioSpecV1:
            raise TypeError("scenario must be a SimulationScenarioSpecV1")
        expected_scenario_fingerprint = simulation_scenario_spec_fingerprint(
            self.scenario
        )
        if self.scenario_fingerprint != expected_scenario_fingerprint:
            raise StructuralTransitionValidationError(
                "scenario_fingerprint does not match scenario"
            )
        if type(self.cost_policy) is not CostPolicyV1:
            raise TypeError("cost_policy must be a CostPolicyV1")
        if self.cost_policy_fingerprint != cost_policy_fingerprint(self.cost_policy):
            raise StructuralTransitionValidationError(
                "cost_policy_fingerprint does not match cost_policy"
            )
        for name in (
            "components_before",
            "projected_components_after",
            "score_records_per_session",
            "base_rating_sessions",
            "repaired_rating_sessions",
            "base_score_records",
            "repaired_score_records",
            "base_calibration_events",
            "repaired_calibration_events",
        ):
            _require_nonnegative_int(name, getattr(self, name))
        if self.components_before <= 1:
            raise StructuralTransitionValidationError(
                "a non-zero K-1 transition requires more than one base component"
            )
        if self.projected_components_after != 1:
            raise StructuralTransitionValidationError(
                "the structural repair must project one connected component"
            )
        if self.score_records_per_session < 1:
            raise StructuralTransitionValidationError(
                "score_records_per_session must be positive"
            )
        if (
            self.repaired_rating_sessions - self.base_rating_sessions
            != self.components_before - 1
        ):
            raise StructuralTransitionValidationError(
                "added rating sessions must equal components_before - 1"
            )
        if self.repaired_rating_sessions <= self.base_rating_sessions:
            raise StructuralTransitionValidationError(
                "repaired rating sessions must exceed the disconnected base"
            )
        if self.repaired_score_records <= self.base_score_records:
            raise StructuralTransitionValidationError(
                "repaired score records must exceed the disconnected base"
            )
        if self.base_score_records != (
            self.base_rating_sessions * self.score_records_per_session
        ):
            raise StructuralTransitionValidationError(
                "base score records do not match rating sessions times rubric width"
            )
        if self.repaired_score_records != (
            self.repaired_rating_sessions * self.score_records_per_session
        ):
            raise StructuralTransitionValidationError(
                "repaired score records do not match rating sessions times rubric width"
            )
        if self.repaired_calibration_events != self.base_calibration_events:
            raise StructuralTransitionValidationError(
                "a v1 structural bridge cannot add calibration events"
            )
        base_cost = _require_cost("base_total_cost", self.base_total_cost)
        repaired_cost = _require_cost("repaired_total_cost", self.repaired_total_cost)
        if repaired_cost < base_cost:
            raise StructuralTransitionValidationError(
                "repaired_total_cost cannot be below base_total_cost"
            )
        expected_base_cost = _weighted_total_cost(
            rating_sessions=self.base_rating_sessions,
            score_records=self.base_score_records,
            calibration_events=self.base_calibration_events,
            policy=self.cost_policy,
        )
        expected_delta_cost = _weighted_total_cost(
            rating_sessions=self.repaired_rating_sessions - self.base_rating_sessions,
            score_records=self.repaired_score_records - self.base_score_records,
            calibration_events=(
                self.repaired_calibration_events - self.base_calibration_events
            ),
            policy=self.cost_policy,
        )
        if base_cost != expected_base_cost:
            raise StructuralTransitionValidationError(
                "base_total_cost does not match workload times unit prices"
            )
        if repaired_cost != expected_base_cost + expected_delta_cost:
            raise StructuralTransitionValidationError(
                "repaired_total_cost does not match base plus repair cost"
            )
        if self.planned_person_rater_graph_connected_before is not False:
            raise StructuralTransitionValidationError(
                "planned_person_rater_graph_connected_before must remain false"
            )
        if self.planned_person_rater_graph_connected_after is not True:
            raise StructuralTransitionValidationError(
                "planned_person_rater_graph_connected_after must remain true"
            )

        if self.scenario.estimator.method == ESTIMATOR_JMLE:
            expected_base = (
                BASE_JMLE_BLOCKED,
                IDENTIFICATION_NONE,
                CLAIM_NOT_ESTIMABLE_COMMON_SCALE,
            )
        elif self.scenario.estimator.method == ESTIMATOR_MML:
            expected_base = (
                BASE_MML_PRIOR_CONDITIONAL,
                IDENTIFICATION_FIXED_POPULATION_PRIOR,
                CLAIM_PRIOR_CONDITIONAL_EXPLORATORY,
            )
        else:  # The estimator contract already prevents this branch.
            raise StructuralTransitionValidationError("unsupported estimator method")
        actual_base = (
            self.base_fit_preflight_status,
            self.base_identification_basis,
            self.base_claim_tier,
        )
        if actual_base != expected_base:
            raise StructuralTransitionValidationError(
                "base estimator status does not match the disconnected design"
            )
        if self.scenario.estimator.method == ESTIMATOR_JMLE:
            expected_repaired = (
                REPAIRED_PENDING_PREFLIGHT,
                IDENTIFICATION_PLANNED_PERSON_RATER,
                CLAIM_STRUCTURAL_ONLY_PENDING_FIT,
            )
        else:
            expected_repaired = (
                REPAIRED_PENDING_PREFLIGHT,
                IDENTIFICATION_PLANNED_PERSON_RATER_PLUS_FIXED_PRIOR,
                CLAIM_PLANNED_CONNECTIVITY_PLUS_PRIOR_PENDING_FIT,
            )
        actual_repaired = (
            self.repaired_fit_preflight_status,
            self.repaired_identification_basis,
            self.repaired_claim_tier,
        )
        if actual_repaired != expected_repaired:
            raise StructuralTransitionValidationError(
                "repaired estimator status does not match the connected projection"
            )
        if self.precision_change_status != PRECISION_CHANGE_NOT_COMPARABLE:
            raise StructuralTransitionValidationError(
                "precision change must remain not comparable for K-1 repair"
            )
        if self.primary_claim_eligible is not False:
            raise StructuralTransitionValidationError(
                "a structural estimability transition is not primary-claim eligible"
            )
        if self.precision_enhancement_overlay_required is not True:
            raise StructuralTransitionValidationError(
                "precision enhancement requires a separate connected-design overlay"
            )
        if self.paired_precision_evaluation_applicable is not False:
            raise StructuralTransitionValidationError(
                "paired precision evaluation is not applicable to a K-1 transition"
            )
        if self.execution_status != EXECUTION_NOT_STARTED:
            raise StructuralTransitionValidationError(
                "the preflight cannot claim simulation execution"
            )
        if self.precision_evaluation_status != PRECISION_NOT_EVALUATED:
            raise StructuralTransitionValidationError(
                "the preflight cannot claim evaluated precision"
            )
        for name in (
            "expected_se_improvement_percent",
            "recovery_improvement_percent",
            "classification_improvement_percent",
            "sample_size_recommendation",
        ):
            if getattr(self, name) is not None:
                raise StructuralTransitionValidationError(
                    f"{name} must remain None for an estimability transition"
                )

        expected_formula_id = _formula_evaluation_id(
            design_spec_fingerprint=self.design_spec_fingerprint,
            cost_policy_fingerprint=self.cost_policy_fingerprint,
        )
        if self.formula_evaluation_id != expected_formula_id:
            raise StructuralTransitionValidationError(
                "formula_evaluation_id does not match design and cost identities"
            )

        expected_transition_id = _transition_id(
            design_spec_fingerprint=self.design_spec_fingerprint,
            assignment_fingerprint=self.assignment_fingerprint,
            repair_fingerprint=self.repair_fingerprint,
            scenario_fingerprint=self.scenario_fingerprint,
            components_before=self.components_before,
            score_records_per_session=self.score_records_per_session,
            base_rating_sessions=self.base_rating_sessions,
            repaired_rating_sessions=self.repaired_rating_sessions,
            base_score_records=self.base_score_records,
            repaired_score_records=self.repaired_score_records,
        )
        if self.transition_id != expected_transition_id:
            raise StructuralTransitionValidationError(
                "transition_id does not match design, repair, and scenario"
            )
        expected_economic_id = _economic_context_id(
            transition_id=self.transition_id,
            formula_evaluation_id=self.formula_evaluation_id,
            cost_policy_fingerprint=self.cost_policy_fingerprint,
            base_calibration_events=self.base_calibration_events,
            repaired_calibration_events=self.repaired_calibration_events,
            base_total_cost=base_cost,
            repaired_total_cost=repaired_cost,
        )
        if self.economic_context_id != expected_economic_id:
            raise StructuralTransitionValidationError(
                "economic_context_id does not match transition and cost context"
            )

    def to_dict(self) -> dict:
        payload = {
            "schema_version": self.schema_version,
            "algorithm": self.algorithm,
            "transition_id": self.transition_id,
            "economic_context_id": self.economic_context_id,
            "comparison_kind": self.comparison_kind,
            "scenario": self.scenario.to_dict(),
            "scenario_fingerprint": self.scenario_fingerprint,
            "formula_evaluation_id": self.formula_evaluation_id,
            "design_spec_fingerprint": self.design_spec_fingerprint,
            "cost_policy": self.cost_policy.to_dict(),
            "cost_policy_fingerprint": self.cost_policy_fingerprint,
            "assignment_fingerprint": self.assignment_fingerprint,
            "repair_fingerprint": self.repair_fingerprint,
            "components_before": self.components_before,
            "projected_components_after": self.projected_components_after,
            "score_records_per_session": self.score_records_per_session,
            "base_rating_sessions": self.base_rating_sessions,
            "repaired_rating_sessions": self.repaired_rating_sessions,
            "base_score_records": self.base_score_records,
            "repaired_score_records": self.repaired_score_records,
            "base_calibration_events": self.base_calibration_events,
            "repaired_calibration_events": self.repaired_calibration_events,
            "base_total_cost": float(self.base_total_cost),
            "repaired_total_cost": float(self.repaired_total_cost),
            "planned_person_rater_graph_connected_before": (
                self.planned_person_rater_graph_connected_before
            ),
            "planned_person_rater_graph_connected_after": (
                self.planned_person_rater_graph_connected_after
            ),
            "base_fit_preflight_status": self.base_fit_preflight_status,
            "base_identification_basis": self.base_identification_basis,
            "base_claim_tier": self.base_claim_tier,
            "repaired_fit_preflight_status": self.repaired_fit_preflight_status,
            "repaired_identification_basis": self.repaired_identification_basis,
            "repaired_claim_tier": self.repaired_claim_tier,
            "precision_change_status": self.precision_change_status,
            "primary_claim_eligible": self.primary_claim_eligible,
            "precision_enhancement_overlay_required": (
                self.precision_enhancement_overlay_required
            ),
            "paired_precision_evaluation_applicable": (
                self.paired_precision_evaluation_applicable
            ),
            "execution_status": self.execution_status,
            "precision_evaluation_status": self.precision_evaluation_status,
            "expected_se_improvement_percent": self.expected_se_improvement_percent,
            "recovery_improvement_percent": self.recovery_improvement_percent,
            "classification_improvement_percent": (
                self.classification_improvement_percent
            ),
            "sample_size_recommendation": self.sample_size_recommendation,
        }
        json.dumps(payload, allow_nan=False, ensure_ascii=False)
        return payload


def normalize_structural_transition_preflight(
    value: StructuralTransitionPreflightV1 | Mapping,
) -> StructuralTransitionPreflightV1:
    """Strictly restore a saved preflight and re-run all cross-field checks."""
    if type(value) is StructuralTransitionPreflightV1:
        return value
    if not isinstance(value, Mapping):
        raise TypeError(
            "transition preflight must be a StructuralTransitionPreflightV1 or mapping"
        )
    payload = _strict_mapping(value)
    payload["scenario"] = normalize_simulation_scenario_spec(payload["scenario"])
    payload["cost_policy"] = normalize_cost_policy(payload["cost_policy"])
    return StructuralTransitionPreflightV1(**payload)


def build_structural_transition_preflight(
    structural: StructuralPlanSummaryV1,
    scenario: SimulationScenarioSpecV1,
) -> StructuralTransitionPreflightV1:
    """Bind a scientific scenario to a non-zero K-1 repair without fitting.

    A connected base with a zero repair is rejected: it needs a future
    connected-to-connected precision-enhancement overlay, not this API.
    """
    if type(structural) is not StructuralPlanSummaryV1:
        raise TypeError("structural must be a StructuralPlanSummaryV1")
    if type(scenario) is not SimulationScenarioSpecV1:
        raise TypeError("scenario must be a SimulationScenarioSpecV1")
    if structural.structure_evaluation_status != STRUCTURE_EVALUATED:
        raise StructuralTransitionValidationError(
            "structural topology must be explicitly evaluated before simulation"
        )
    if structural.anchor_augmented is None or structural.repair is None:
        raise StructuralTransitionValidationError(
            "evaluated topology and repair evidence are required"
        )
    repair = structural.repair
    if (
        structural.anchor_augmented.connected
        or repair.status != REPAIR_PROPOSED
        or repair.added_rating_sessions <= 0
    ):
        raise StructuralTransitionValidationError(
            "a connected or zero-repair design requires a separate "
            "precision-enhancement overlay"
        )
    if structural.assignment_fingerprint is None:
        raise StructuralTransitionValidationError(
            "evaluated structure requires an assignment fingerprint"
        )
    if (
        structural.projected_total_rating_sessions is None
        or structural.projected_total_score_records is None
        or structural.projected_calibration_events is None
        or structural.projected_total_cost is None
    ):
        raise StructuralTransitionValidationError(
            "evaluated structure requires projected workload and cost"
        )

    formula = structural.formula_summary
    _validate_simulation_id_capacity(structural)
    scenario_fingerprint = simulation_scenario_spec_fingerprint(scenario)
    transition_id = _transition_id(
        design_spec_fingerprint=formula.design_spec_fingerprint,
        assignment_fingerprint=structural.assignment_fingerprint,
        repair_fingerprint=repair.repair_fingerprint,
        scenario_fingerprint=scenario_fingerprint,
        components_before=repair.components_before,
        score_records_per_session=formula.design_spec.score_records_per_session,
        base_rating_sessions=formula.workload.total_rating_sessions,
        repaired_rating_sessions=structural.projected_total_rating_sessions,
        base_score_records=formula.workload.total_score_records,
        repaired_score_records=structural.projected_total_score_records,
    )
    economic_id = _economic_context_id(
        transition_id=transition_id,
        formula_evaluation_id=formula.evaluation_id,
        cost_policy_fingerprint=formula.cost_policy_fingerprint,
        base_calibration_events=formula.workload.calibration_events,
        repaired_calibration_events=structural.projected_calibration_events,
        base_total_cost=formula.weighted_cost.total_cost,
        repaired_total_cost=structural.projected_total_cost,
    )
    if scenario.estimator.method == ESTIMATOR_JMLE:
        base_status = BASE_JMLE_BLOCKED
        base_basis = IDENTIFICATION_NONE
        base_tier = CLAIM_NOT_ESTIMABLE_COMMON_SCALE
        repaired_basis = IDENTIFICATION_PLANNED_PERSON_RATER
        repaired_tier = CLAIM_STRUCTURAL_ONLY_PENDING_FIT
    else:
        base_status = BASE_MML_PRIOR_CONDITIONAL
        base_basis = IDENTIFICATION_FIXED_POPULATION_PRIOR
        base_tier = CLAIM_PRIOR_CONDITIONAL_EXPLORATORY
        repaired_basis = IDENTIFICATION_PLANNED_PERSON_RATER_PLUS_FIXED_PRIOR
        repaired_tier = CLAIM_PLANNED_CONNECTIVITY_PLUS_PRIOR_PENDING_FIT

    return StructuralTransitionPreflightV1(
        schema_version=STRUCTURAL_TRANSITION_PREFLIGHT_VERSION,
        algorithm=STRUCTURAL_TRANSITION_ALGORITHM_V1,
        transition_id=transition_id,
        economic_context_id=economic_id,
        comparison_kind=COMPARISON_ESTIMABILITY_TRANSITION,
        scenario=scenario,
        scenario_fingerprint=scenario_fingerprint,
        formula_evaluation_id=formula.evaluation_id,
        design_spec_fingerprint=formula.design_spec_fingerprint,
        cost_policy=formula.cost_policy,
        cost_policy_fingerprint=formula.cost_policy_fingerprint,
        assignment_fingerprint=structural.assignment_fingerprint,
        repair_fingerprint=repair.repair_fingerprint,
        components_before=repair.components_before,
        projected_components_after=repair.projected_components_after,
        score_records_per_session=formula.design_spec.score_records_per_session,
        base_rating_sessions=formula.workload.total_rating_sessions,
        repaired_rating_sessions=structural.projected_total_rating_sessions,
        base_score_records=formula.workload.total_score_records,
        repaired_score_records=structural.projected_total_score_records,
        base_calibration_events=formula.workload.calibration_events,
        repaired_calibration_events=structural.projected_calibration_events,
        base_total_cost=formula.weighted_cost.total_cost,
        repaired_total_cost=structural.projected_total_cost,
        planned_person_rater_graph_connected_before=False,
        planned_person_rater_graph_connected_after=True,
        base_fit_preflight_status=base_status,
        base_identification_basis=base_basis,
        base_claim_tier=base_tier,
        repaired_fit_preflight_status=REPAIRED_PENDING_PREFLIGHT,
        repaired_identification_basis=repaired_basis,
        repaired_claim_tier=repaired_tier,
        precision_change_status=PRECISION_CHANGE_NOT_COMPARABLE,
        primary_claim_eligible=False,
        precision_enhancement_overlay_required=True,
        paired_precision_evaluation_applicable=False,
        execution_status=EXECUTION_NOT_STARTED,
        precision_evaluation_status=PRECISION_NOT_EVALUATED,
        expected_se_improvement_percent=None,
        recovery_improvement_percent=None,
        classification_improvement_percent=None,
        sample_size_recommendation=None,
    )


__all__ = [
    "BASE_JMLE_BLOCKED",
    "BASE_MML_PRIOR_CONDITIONAL",
    "CLAIM_NOT_ESTIMABLE_COMMON_SCALE",
    "CLAIM_PLANNED_CONNECTIVITY_PLUS_PRIOR_PENDING_FIT",
    "CLAIM_PRIOR_CONDITIONAL_EXPLORATORY",
    "CLAIM_STRUCTURAL_ONLY_PENDING_FIT",
    "COMPARISON_ESTIMABILITY_TRANSITION",
    "EXECUTION_NOT_STARTED",
    "IDENTIFICATION_FIXED_POPULATION_PRIOR",
    "IDENTIFICATION_NONE",
    "IDENTIFICATION_PLANNED_PERSON_RATER",
    "IDENTIFICATION_PLANNED_PERSON_RATER_PLUS_FIXED_PRIOR",
    "MAX_STRUCTURAL_TRANSITION_COUNT",
    "PRECISION_CHANGE_NOT_COMPARABLE",
    "REPAIRED_PENDING_PREFLIGHT",
    "STRUCTURAL_TRANSITION_ALGORITHM_V1",
    "STRUCTURAL_TRANSITION_PREFLIGHT_VERSION",
    "StructuralTransitionPreflightV1",
    "StructuralTransitionValidationError",
    "build_structural_transition_preflight",
    "normalize_structural_transition_preflight",
]
