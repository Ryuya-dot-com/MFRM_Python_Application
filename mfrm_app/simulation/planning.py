"""Compact two-stage planning summaries for the prospective-design UI.

The formula stage is O(1) in the number of scheduled ratings and is safe to
recompute while widgets change.  The structural stage is explicit and lazy:
it materializes a deterministic assignment only below a caller-supplied
synchronous limit, audits Person-Rater topology, and keeps only compact
evidence plus the small K-1 repair overlay.

Neither stage generates scores, fits an estimator, or evaluates precision,
recovery, classification agreement, anchor invariance, or sample-size
adequacy.  Those fixed nonclaims are part of the returned contract so a UI
cannot accidentally turn an unevaluated field into a blank success state.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
import hashlib
import json
import math
from statistics import median

from .assignments import ASSIGNMENT_BUNDLE_VERSION, compile_rating_design
from .cost import (
    COST_POLICY_VERSION,
    CostPolicyV1,
    DesignWorkloadV1,
    WeightedDesignCostV1,
    apply_cost_policy,
    compute_design_workload,
    cost_policy_fingerprint,
    normalize_cost_policy,
)
from .schema import (
    DesignSpecV1,
    design_spec_fingerprint,
    design_spec_to_dict,
    normalize_design_spec,
)
from .topology import (
    ANCHOR_AUGMENTED_SCOPE,
    PRECISION_NOT_EVALUATED,
    PRIMARY_STUDY_SCOPE,
    STRUCTURAL_REPAIR_VERSION,
    TOPOLOGY_AUDIT_VERSION,
    StructuralBridgeRepairV1,
    TopologyScopeAuditV1,
    audit_rating_design_topology,
    propose_structural_bridge_repair,
)


FORMULA_PLAN_SUMMARY_VERSION = "mfrm_formula_plan_summary_v1"
STRUCTURAL_PLAN_SUMMARY_VERSION = "mfrm_structural_plan_summary_v1"
TOPOLOGY_SCOPE_PLAN_SUMMARY_VERSION = "mfrm_topology_scope_plan_summary_v1"
REPAIR_COST_DELTA_VERSION = "mfrm_repair_cost_delta_v1"
PLANNING_ALGORITHM_V1 = "prospective_rating_design_planner_v1"

STRUCTURE_NOT_EVALUATED = "not_evaluated"
STRUCTURE_EVALUATED = "evaluated"
STRUCTURE_NOT_EVALUATED_LIMIT = "not_evaluated_limit"
STRUCTURE_EVALUATION_STATUS_CHOICES = (
    STRUCTURE_NOT_EVALUATED,
    STRUCTURE_EVALUATED,
    STRUCTURE_NOT_EVALUATED_LIMIT,
)
MATERIALIZATION_EXPLICIT_REQUEST_REQUIRED = "explicit_request_required"
MATERIALIZATION_WITHIN_LIMIT = "within_limit"
MATERIALIZATION_RATING_SESSION_LIMIT = "rating_session_limit"
DEFAULT_SYNCHRONOUS_RATING_SESSION_LIMIT = 10_000


def _fingerprint(payload: Mapping, *, length: int = 16) -> str:
    encoded = json.dumps(
        dict(payload),
        allow_nan=False,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()[:length]


def _require_fingerprint(name: str, value: object) -> str:
    if (
        not isinstance(value, str)
        or len(value) != 16
        or value != value.lower()
        or any(character not in "0123456789abcdef" for character in value)
    ):
        raise ValueError(f"{name} must be a 16-character lowercase hex fingerprint")
    return value


def _require_positive_int(name: str, value: object) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f"{name} must be a positive integer")
    return value


def _require_nonnegative_int(name: str, value: object) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


def _json_safe(payload: dict) -> dict:
    json.dumps(payload, allow_nan=False, ensure_ascii=False)
    return payload


@dataclass(frozen=True, slots=True, kw_only=True)
class FormulaPlanSummaryV1:
    """Formula-only workload and cost evidence for one exact design/policy."""

    schema_version: str
    planning_algorithm: str
    evaluation_id: str
    design_spec: DesignSpecV1
    design_spec_fingerprint: str
    cost_policy: CostPolicyV1
    cost_policy_fingerprint: str
    workload: DesignWorkloadV1
    weighted_cost: WeightedDesignCostV1
    structure_evaluation_status: str
    precision_evaluation_status: str
    expected_se: None
    recovery_accuracy: None
    classification_agreement: None
    sample_size_recommendation: None

    def __post_init__(self) -> None:
        if self.schema_version != FORMULA_PLAN_SUMMARY_VERSION:
            raise ValueError(f"Unsupported formula summary version: {self.schema_version!r}")
        if self.planning_algorithm != PLANNING_ALGORITHM_V1:
            raise ValueError(f"Unsupported planning algorithm: {self.planning_algorithm!r}")
        _require_fingerprint("evaluation_id", self.evaluation_id)
        if type(self.design_spec) is not DesignSpecV1:
            raise TypeError("design_spec must be a DesignSpecV1")
        expected_design_fingerprint = design_spec_fingerprint(self.design_spec)
        if self.design_spec_fingerprint != expected_design_fingerprint:
            raise ValueError("design_spec_fingerprint does not match design_spec")
        if type(self.cost_policy) is not CostPolicyV1:
            raise TypeError("cost_policy must be a CostPolicyV1")
        expected_policy_fingerprint = cost_policy_fingerprint(self.cost_policy)
        if self.cost_policy_fingerprint != expected_policy_fingerprint:
            raise ValueError("cost_policy_fingerprint does not match cost_policy")
        if type(self.workload) is not DesignWorkloadV1:
            raise TypeError("workload must be a DesignWorkloadV1")
        if self.workload != compute_design_workload(self.design_spec):
            raise ValueError("workload does not match design_spec")
        if type(self.weighted_cost) is not WeightedDesignCostV1:
            raise TypeError("weighted_cost must be a WeightedDesignCostV1")
        if self.weighted_cost != apply_cost_policy(self.workload, self.cost_policy):
            raise ValueError("weighted_cost does not match workload and cost_policy")
        expected_id = _formula_evaluation_id(
            self.design_spec_fingerprint,
            self.cost_policy_fingerprint,
        )
        if self.evaluation_id != expected_id:
            raise ValueError("evaluation_id does not match design and cost identities")
        if self.structure_evaluation_status != STRUCTURE_NOT_EVALUATED:
            raise ValueError("a formula summary cannot claim evaluated topology")
        _validate_fixed_nonclaims(self)

    def to_dict(self) -> dict:
        return _json_safe({
            "schema_version": self.schema_version,
            "planning_algorithm": self.planning_algorithm,
            "evaluation_id": self.evaluation_id,
            "design_spec": design_spec_to_dict(self.design_spec),
            "design_spec_fingerprint": self.design_spec_fingerprint,
            "cost_policy": self.cost_policy.to_dict(),
            "cost_policy_fingerprint": self.cost_policy_fingerprint,
            "workload": self.workload.to_dict(),
            "weighted_cost": self.weighted_cost.to_dict(),
            "structure_evaluation_status": self.structure_evaluation_status,
            "precision_evaluation_status": self.precision_evaluation_status,
            "expected_se": self.expected_se,
            "recovery_accuracy": self.recovery_accuracy,
            "classification_agreement": self.classification_agreement,
            "sample_size_recommendation": self.sample_size_recommendation,
        })


def _formula_evaluation_id(
    design_fingerprint: str,
    policy_fingerprint: str,
) -> str:
    return _fingerprint({
        "schema_version": FORMULA_PLAN_SUMMARY_VERSION,
        "planning_algorithm": PLANNING_ALGORITHM_V1,
        "design_spec_fingerprint": design_fingerprint,
        "cost_policy_fingerprint": policy_fingerprint,
    })


def _validate_fixed_nonclaims(value: object) -> None:
    if getattr(value, "precision_evaluation_status") != PRECISION_NOT_EVALUATED:
        raise ValueError("precision_evaluation_status must remain not_evaluated")
    for name in (
        "expected_se",
        "recovery_accuracy",
        "classification_agreement",
        "sample_size_recommendation",
    ):
        if getattr(value, name) is not None:
            raise ValueError(f"{name} is not evaluated by the planning contract")


def _resolve_cost_policy(value: CostPolicyV1 | Mapping | None) -> CostPolicyV1:
    if value is None:
        return CostPolicyV1()
    return normalize_cost_policy(value)


def build_formula_plan_summary(
    spec: DesignSpecV1 | Mapping,
    policy: CostPolicyV1 | Mapping | None = None,
) -> FormulaPlanSummaryV1:
    """Return exact O(1) workload/cost arithmetic without compiling rows."""
    resolved_spec = normalize_design_spec(spec)
    resolved_policy = _resolve_cost_policy(policy)
    workload = compute_design_workload(resolved_spec)
    weighted_cost = apply_cost_policy(workload, resolved_policy)
    design_fingerprint = design_spec_fingerprint(resolved_spec)
    policy_fingerprint = cost_policy_fingerprint(resolved_policy)
    return FormulaPlanSummaryV1(
        schema_version=FORMULA_PLAN_SUMMARY_VERSION,
        planning_algorithm=PLANNING_ALGORITHM_V1,
        evaluation_id=_formula_evaluation_id(design_fingerprint, policy_fingerprint),
        design_spec=resolved_spec,
        design_spec_fingerprint=design_fingerprint,
        cost_policy=resolved_policy,
        cost_policy_fingerprint=policy_fingerprint,
        workload=workload,
        weighted_cost=weighted_cost,
        structure_evaluation_status=STRUCTURE_NOT_EVALUATED,
        precision_evaluation_status=PRECISION_NOT_EVALUATED,
        expected_se=None,
        recovery_accuracy=None,
        classification_agreement=None,
        sample_size_recommendation=None,
    )


@dataclass(frozen=True, slots=True, kw_only=True)
class TopologyScopePlanSummaryV1:
    """Compact topology evidence without component membership or edge rows."""

    schema_version: str
    scope: str
    connected: bool
    n_components: int
    isolated_raters: int
    unique_incidence_edges: int
    rating_sessions: int
    rater_load_min: int
    rater_load_median: float
    rater_load_max: int
    rater_load_spread: int

    def __post_init__(self) -> None:
        if self.schema_version != TOPOLOGY_SCOPE_PLAN_SUMMARY_VERSION:
            raise ValueError(f"Unsupported topology summary version: {self.schema_version!r}")
        if self.scope not in (PRIMARY_STUDY_SCOPE, ANCHOR_AUGMENTED_SCOPE):
            raise ValueError(f"Unknown topology scope: {self.scope!r}")
        if not isinstance(self.connected, bool):
            raise TypeError("connected must be boolean")
        for name in (
            "n_components",
            "isolated_raters",
            "unique_incidence_edges",
            "rating_sessions",
            "rater_load_min",
            "rater_load_max",
            "rater_load_spread",
        ):
            _require_nonnegative_int(name, getattr(self, name))
        if self.n_components < 1:
            raise ValueError("n_components must be positive")
        if (
            isinstance(self.rater_load_median, bool)
            or not isinstance(self.rater_load_median, (int, float))
            or not math.isfinite(float(self.rater_load_median))
            or self.rater_load_median < 0
        ):
            raise ValueError("rater_load_median must be finite and nonnegative")
        if self.rater_load_min > self.rater_load_max:
            raise ValueError("rater load minimum cannot exceed maximum")
        if self.rater_load_spread != self.rater_load_max - self.rater_load_min:
            raise ValueError("rater_load_spread must equal max minus min")
        if self.connected != (self.n_components == 1):
            raise ValueError("connected must agree with n_components")

    def to_dict(self) -> dict:
        return _json_safe({
            "schema_version": self.schema_version,
            "scope": self.scope,
            "connected": self.connected,
            "n_components": self.n_components,
            "isolated_raters": self.isolated_raters,
            "unique_incidence_edges": self.unique_incidence_edges,
            "rating_sessions": self.rating_sessions,
            "rater_load_min": self.rater_load_min,
            "rater_load_median": float(self.rater_load_median),
            "rater_load_max": self.rater_load_max,
            "rater_load_spread": self.rater_load_spread,
        })


def _compact_scope(scope: TopologyScopeAuditV1) -> TopologyScopePlanSummaryV1:
    if type(scope) is not TopologyScopeAuditV1:
        raise TypeError("scope must be a TopologyScopeAuditV1")
    loads = [load.rating_sessions for load in scope.rater_loads]
    if not loads:  # DesignSpecV1 requires at least one active rater.
        raise ValueError("topology scope does not contain declared rater loads")
    minimum = min(loads)
    maximum = max(loads)
    return TopologyScopePlanSummaryV1(
        schema_version=TOPOLOGY_SCOPE_PLAN_SUMMARY_VERSION,
        scope=scope.scope,
        connected=scope.connected,
        n_components=scope.n_components,
        isolated_raters=scope.isolated_raters,
        unique_incidence_edges=scope.unique_incidence_edges,
        rating_sessions=scope.rating_sessions,
        rater_load_min=minimum,
        rater_load_median=float(median(loads)),
        rater_load_max=maximum,
        rater_load_spread=maximum - minimum,
    )


@dataclass(frozen=True, slots=True, kw_only=True)
class RepairCostDeltaV1:
    """Exact cost of the structural overlay in the three workload currencies."""

    schema_version: str
    design_spec_fingerprint: str
    assignment_fingerprint: str
    repair_fingerprint: str
    cost_policy_schema_version: str
    cost_policy_fingerprint: str
    cost_unit: str
    added_rating_sessions: int
    added_score_records: int
    added_calibration_events: int
    rating_session_unit_cost: float
    score_record_unit_cost: float
    calibration_event_unit_cost: float
    rating_session_cost_delta: float
    score_record_cost_delta: float
    calibration_event_cost_delta: float
    total_cost_delta: float

    def __post_init__(self) -> None:
        if self.schema_version != REPAIR_COST_DELTA_VERSION:
            raise ValueError(f"Unsupported repair cost version: {self.schema_version!r}")
        for name in (
            "design_spec_fingerprint",
            "assignment_fingerprint",
            "repair_fingerprint",
            "cost_policy_fingerprint",
        ):
            _require_fingerprint(name, getattr(self, name))
        if self.cost_policy_schema_version != COST_POLICY_VERSION:
            raise ValueError("repair cost requires the current cost-policy version")
        if not isinstance(self.cost_unit, str) or not self.cost_unit.strip():
            raise ValueError("cost_unit must be non-empty")
        if self.cost_unit != self.cost_unit.strip():
            raise ValueError("cost_unit must be canonical without outer whitespace")
        for name in (
            "added_rating_sessions",
            "added_score_records",
            "added_calibration_events",
        ):
            _require_nonnegative_int(name, getattr(self, name))
        numeric_names = (
            "rating_session_unit_cost",
            "score_record_unit_cost",
            "calibration_event_unit_cost",
            "rating_session_cost_delta",
            "score_record_cost_delta",
            "calibration_event_cost_delta",
            "total_cost_delta",
        )
        for name in numeric_names:
            value = getattr(self, name)
            if (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(float(value))
                or value < 0
            ):
                raise ValueError(f"{name} must be finite and nonnegative")
        expected_policy = CostPolicyV1(
            rating_session_unit_cost=self.rating_session_unit_cost,
            score_record_unit_cost=self.score_record_unit_cost,
            calibration_event_unit_cost=self.calibration_event_unit_cost,
            cost_unit=self.cost_unit,
            schema_version=self.cost_policy_schema_version,
        )
        if cost_policy_fingerprint(expected_policy) != self.cost_policy_fingerprint:
            raise ValueError("cost_policy_fingerprint does not match embedded prices")
        expected = (
            self.added_rating_sessions * self.rating_session_unit_cost,
            self.added_score_records * self.score_record_unit_cost,
            self.added_calibration_events * self.calibration_event_unit_cost,
        )
        if self.rating_session_cost_delta != expected[0]:
            raise ValueError("rating_session_cost_delta does not match count times price")
        if self.score_record_cost_delta != expected[1]:
            raise ValueError("score_record_cost_delta does not match count times price")
        if self.calibration_event_cost_delta != expected[2]:
            raise ValueError("calibration_event_cost_delta does not match count times price")
        if self.total_cost_delta != sum(expected):
            raise ValueError("total_cost_delta must equal component cost deltas")

    def to_dict(self) -> dict:
        return _json_safe({name: getattr(self, name) for name in self.__slots__})


def _finite_product(count: int, unit_cost: float, *, name: str) -> float:
    try:
        value = float(count * unit_cost)
    except OverflowError as exc:
        raise ValueError(f"{name} exceeds the finite weighted-cost range") from exc
    if not math.isfinite(value):
        raise ValueError(f"{name} exceeds the finite weighted-cost range")
    return value


def compute_repair_cost_delta(
    repair: StructuralBridgeRepairV1,
    policy: CostPolicyV1 | Mapping | None = None,
) -> RepairCostDeltaV1:
    """Price a deterministic repair overlay without mutating its base design."""
    if type(repair) is not StructuralBridgeRepairV1:
        raise TypeError("repair must be a StructuralBridgeRepairV1")
    resolved = _resolve_cost_policy(policy)
    session_cost = _finite_product(
        repair.added_rating_sessions,
        resolved.rating_session_unit_cost,
        name="rating_session_cost_delta",
    )
    record_cost = _finite_product(
        repair.added_score_records,
        resolved.score_record_unit_cost,
        name="score_record_cost_delta",
    )
    calibration_cost = _finite_product(
        repair.added_calibration_events,
        resolved.calibration_event_unit_cost,
        name="calibration_event_cost_delta",
    )
    total = session_cost + record_cost + calibration_cost
    if not math.isfinite(total):
        raise ValueError("total_cost_delta exceeds the finite weighted-cost range")
    return RepairCostDeltaV1(
        schema_version=REPAIR_COST_DELTA_VERSION,
        design_spec_fingerprint=repair.design_spec_fingerprint,
        assignment_fingerprint=repair.assignment_fingerprint,
        repair_fingerprint=repair.repair_fingerprint,
        cost_policy_schema_version=resolved.schema_version,
        cost_policy_fingerprint=cost_policy_fingerprint(resolved),
        cost_unit=resolved.cost_unit,
        added_rating_sessions=repair.added_rating_sessions,
        added_score_records=repair.added_score_records,
        added_calibration_events=repair.added_calibration_events,
        rating_session_unit_cost=resolved.rating_session_unit_cost,
        score_record_unit_cost=resolved.score_record_unit_cost,
        calibration_event_unit_cost=resolved.calibration_event_unit_cost,
        rating_session_cost_delta=session_cost,
        score_record_cost_delta=record_cost,
        calibration_event_cost_delta=calibration_cost,
        total_cost_delta=total,
    )


@dataclass(frozen=True, slots=True, kw_only=True)
class StructuralPlanSummaryV1:
    """Compact result of one explicit assignment/topology evaluation."""

    schema_version: str
    planning_algorithm: str
    formula_summary: FormulaPlanSummaryV1
    structure_evaluation_status: str
    materialization_reason: str
    required_rating_sessions: int
    max_rating_sessions: int
    assignment_bundle_schema_version: str | None
    assignment_fingerprint: str | None
    topology_audit_schema_version: str | None
    primary_study: TopologyScopePlanSummaryV1 | None
    anchor_augmented: TopologyScopePlanSummaryV1 | None
    repair: StructuralBridgeRepairV1 | None
    repair_cost_delta: RepairCostDeltaV1 | None
    projected_total_rating_sessions: int | None
    projected_total_score_records: int | None
    projected_calibration_events: int | None
    projected_total_cost: float | None
    precision_evaluation_status: str
    expected_se: None
    recovery_accuracy: None
    classification_agreement: None
    sample_size_recommendation: None

    def __post_init__(self) -> None:
        if self.schema_version != STRUCTURAL_PLAN_SUMMARY_VERSION:
            raise ValueError(f"Unsupported structural summary version: {self.schema_version!r}")
        if self.planning_algorithm != PLANNING_ALGORITHM_V1:
            raise ValueError(f"Unsupported planning algorithm: {self.planning_algorithm!r}")
        if type(self.formula_summary) is not FormulaPlanSummaryV1:
            raise TypeError("formula_summary must be a FormulaPlanSummaryV1")
        _require_positive_int("required_rating_sessions", self.required_rating_sessions)
        _require_positive_int("max_rating_sessions", self.max_rating_sessions)
        if self.required_rating_sessions != self.formula_summary.workload.total_rating_sessions:
            raise ValueError("required_rating_sessions must match formula workload")
        if self.structure_evaluation_status not in (
            STRUCTURE_EVALUATED,
            STRUCTURE_NOT_EVALUATED_LIMIT,
        ):
            raise ValueError("structural summary must be evaluated or explicitly limited")
        _validate_fixed_nonclaims(self)

        evidence = (
            self.assignment_bundle_schema_version,
            self.assignment_fingerprint,
            self.topology_audit_schema_version,
            self.primary_study,
            self.anchor_augmented,
            self.repair,
            self.repair_cost_delta,
            self.projected_total_rating_sessions,
            self.projected_total_score_records,
            self.projected_calibration_events,
            self.projected_total_cost,
        )
        if self.structure_evaluation_status == STRUCTURE_NOT_EVALUATED_LIMIT:
            if self.materialization_reason != MATERIALIZATION_RATING_SESSION_LIMIT:
                raise ValueError("limited structure evaluation requires rating_session_limit reason")
            if self.required_rating_sessions <= self.max_rating_sessions:
                raise ValueError("limited status requires required sessions above the limit")
            if any(item is not None for item in evidence):
                raise ValueError("limited structure evaluation cannot contain topology evidence")
            return

        if self.materialization_reason != MATERIALIZATION_WITHIN_LIMIT:
            raise ValueError("evaluated structure requires within_limit reason")
        if self.required_rating_sessions > self.max_rating_sessions:
            raise ValueError("evaluated structure exceeds max_rating_sessions")
        if self.assignment_bundle_schema_version != ASSIGNMENT_BUNDLE_VERSION:
            raise ValueError("evaluated structure requires the current assignment bundle")
        if self.topology_audit_schema_version != TOPOLOGY_AUDIT_VERSION:
            raise ValueError("evaluated structure requires the current topology audit")
        _require_fingerprint("assignment_fingerprint", self.assignment_fingerprint)
        if type(self.primary_study) is not TopologyScopePlanSummaryV1 or (
            self.primary_study.scope != PRIMARY_STUDY_SCOPE
        ):
            raise ValueError("primary_study must contain compact primary-study evidence")
        if type(self.anchor_augmented) is not TopologyScopePlanSummaryV1 or (
            self.anchor_augmented.scope != ANCHOR_AUGMENTED_SCOPE
        ):
            raise ValueError("anchor_augmented must contain compact augmented evidence")
        expected_primary_sessions = (
            self.formula_summary.workload.primary_rating_sessions
            + self.formula_summary.workload.bridge_rating_sessions
        )
        if self.primary_study.rating_sessions != expected_primary_sessions:
            raise ValueError(
                "primary-study sessions do not match formula primary plus bridge workload"
            )
        if (
            self.anchor_augmented.rating_sessions
            != self.formula_summary.workload.total_rating_sessions
        ):
            raise ValueError(
                "anchor-augmented sessions do not match formula total workload"
            )
        if type(self.repair) is not StructuralBridgeRepairV1:
            raise ValueError("evaluated structure requires a structural repair contract")
        if self.repair.schema_version != STRUCTURAL_REPAIR_VERSION:
            raise ValueError("repair contract version mismatch")
        if self.repair.assignment_fingerprint != self.assignment_fingerprint:
            raise ValueError("repair is not bound to the evaluated assignment")
        if self.repair.design_spec_fingerprint != self.formula_summary.design_spec_fingerprint:
            raise ValueError("repair is not bound to the formula design")
        if self.repair.components_before != self.anchor_augmented.n_components:
            raise ValueError(
                "repair components_before does not match anchor-augmented topology"
            )
        if type(self.repair_cost_delta) is not RepairCostDeltaV1:
            raise ValueError("evaluated structure requires a repair cost delta")
        if self.repair_cost_delta.design_spec_fingerprint != self.repair.design_spec_fingerprint:
            raise ValueError("repair cost is not bound to the repair design")
        if self.repair_cost_delta.assignment_fingerprint != self.repair.assignment_fingerprint:
            raise ValueError("repair cost is not bound to the repair assignment")
        if self.repair_cost_delta.repair_fingerprint != self.repair.repair_fingerprint:
            raise ValueError("repair cost is not bound to the repair overlay")
        if self.repair_cost_delta.cost_policy_fingerprint != self.formula_summary.cost_policy_fingerprint:
            raise ValueError("repair cost is not bound to the formula cost policy")
        expected_repair_cost = compute_repair_cost_delta(
            self.repair,
            self.formula_summary.cost_policy,
        )
        if self.repair_cost_delta != expected_repair_cost:
            raise ValueError(
                "repair cost delta does not exactly match the repair and cost policy"
            )
        expected_counts = (
            self.required_rating_sessions + self.repair.added_rating_sessions,
            self.formula_summary.workload.total_score_records + self.repair.added_score_records,
            self.formula_summary.workload.calibration_events + self.repair.added_calibration_events,
        )
        if (
            self.projected_total_rating_sessions,
            self.projected_total_score_records,
            self.projected_calibration_events,
        ) != expected_counts:
            raise ValueError("projected workload does not equal base plus repair")
        expected_cost = self.formula_summary.weighted_cost.total_cost + self.repair_cost_delta.total_cost_delta
        if self.projected_total_cost != expected_cost or not math.isfinite(expected_cost):
            raise ValueError("projected_total_cost does not equal base plus repair")

    def to_dict(self) -> dict:
        return _json_safe({
            "schema_version": self.schema_version,
            "planning_algorithm": self.planning_algorithm,
            "formula_summary": self.formula_summary.to_dict(),
            "structure_evaluation_status": self.structure_evaluation_status,
            "materialization_reason": self.materialization_reason,
            "required_rating_sessions": self.required_rating_sessions,
            "max_rating_sessions": self.max_rating_sessions,
            "assignment_bundle_schema_version": self.assignment_bundle_schema_version,
            "assignment_fingerprint": self.assignment_fingerprint,
            "topology_audit_schema_version": self.topology_audit_schema_version,
            "primary_study": self.primary_study.to_dict() if self.primary_study else None,
            "anchor_augmented": self.anchor_augmented.to_dict() if self.anchor_augmented else None,
            "repair": self.repair.to_dict() if self.repair else None,
            "repair_cost_delta": self.repair_cost_delta.to_dict() if self.repair_cost_delta else None,
            "projected_total_rating_sessions": self.projected_total_rating_sessions,
            "projected_total_score_records": self.projected_total_score_records,
            "projected_calibration_events": self.projected_calibration_events,
            "projected_total_cost": self.projected_total_cost,
            "precision_evaluation_status": self.precision_evaluation_status,
            "expected_se": self.expected_se,
            "recovery_accuracy": self.recovery_accuracy,
            "classification_agreement": self.classification_agreement,
            "sample_size_recommendation": self.sample_size_recommendation,
        })


def evaluate_plan_structure(
    summary_or_spec: FormulaPlanSummaryV1 | DesignSpecV1 | Mapping,
    policy: CostPolicyV1 | Mapping | None = None,
    *,
    max_rating_sessions: int = DEFAULT_SYNCHRONOUS_RATING_SESSION_LIMIT,
) -> StructuralPlanSummaryV1:
    """Explicitly compile and audit one plan, or return a bounded non-result."""
    limit = _require_positive_int("max_rating_sessions", max_rating_sessions)
    if isinstance(summary_or_spec, FormulaPlanSummaryV1):
        formula = summary_or_spec
        if policy is not None and _resolve_cost_policy(policy) != formula.cost_policy:
            raise ValueError("policy conflicts with the formula summary cost policy")
    else:
        formula = build_formula_plan_summary(summary_or_spec, policy)
    required = formula.workload.total_rating_sessions
    common = dict(
        schema_version=STRUCTURAL_PLAN_SUMMARY_VERSION,
        planning_algorithm=PLANNING_ALGORITHM_V1,
        formula_summary=formula,
        required_rating_sessions=required,
        max_rating_sessions=limit,
        precision_evaluation_status=PRECISION_NOT_EVALUATED,
        expected_se=None,
        recovery_accuracy=None,
        classification_agreement=None,
        sample_size_recommendation=None,
    )
    if required > limit:
        return StructuralPlanSummaryV1(
            **common,
            structure_evaluation_status=STRUCTURE_NOT_EVALUATED_LIMIT,
            materialization_reason=MATERIALIZATION_RATING_SESSION_LIMIT,
            assignment_bundle_schema_version=None,
            assignment_fingerprint=None,
            topology_audit_schema_version=None,
            primary_study=None,
            anchor_augmented=None,
            repair=None,
            repair_cost_delta=None,
            projected_total_rating_sessions=None,
            projected_total_score_records=None,
            projected_calibration_events=None,
            projected_total_cost=None,
        )

    bundle = compile_rating_design(formula.design_spec, max_rating_sessions=limit)
    audit = audit_rating_design_topology(bundle, max_rating_sessions=limit)
    repair = propose_structural_bridge_repair(
        bundle,
        topology_audit=audit,
        max_rating_sessions=limit,
    )
    repair_cost = compute_repair_cost_delta(repair, formula.cost_policy)
    return StructuralPlanSummaryV1(
        **common,
        structure_evaluation_status=STRUCTURE_EVALUATED,
        materialization_reason=MATERIALIZATION_WITHIN_LIMIT,
        assignment_bundle_schema_version=bundle.schema_version,
        assignment_fingerprint=bundle.assignment_fingerprint,
        topology_audit_schema_version=audit.schema_version,
        primary_study=_compact_scope(audit.primary_study),
        anchor_augmented=_compact_scope(audit.anchor_augmented),
        repair=repair,
        repair_cost_delta=repair_cost,
        projected_total_rating_sessions=required + repair.added_rating_sessions,
        projected_total_score_records=(
            formula.workload.total_score_records + repair.added_score_records
        ),
        projected_calibration_events=(
            formula.workload.calibration_events + repair.added_calibration_events
        ),
        projected_total_cost=formula.weighted_cost.total_cost + repair_cost.total_cost_delta,
    )
