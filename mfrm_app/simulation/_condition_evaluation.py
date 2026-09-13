"""Monte Carlo, classification, and evaluation-policy contracts."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

from ._condition_randomization import (
    RandomizationSpecV1,
    normalize_randomization_spec,
)
from ._condition_shared import (
    FORMAL_MONTE_CARLO_MIN_REPLICATES,
    MAX_MONTE_CARLO_REPLICATES,
    SimulationConditionValidationError,
    _fingerprint,
    _json_safe,
    _require_int,
    _require_literal,
    _require_real,
    _strict_mapping,
    _tuple_payload_field,
)


MONTE_CARLO_SPEC_VERSION = "mfrm_monte_carlo_spec_v1"
CLASSIFICATION_RULE_VERSION = "mfrm_classification_rule_v1"
EVALUATION_POLICY_VERSION = "mfrm_evaluation_policy_v1"

PAIR_COMPARABILITY_REQUIREMENT_V1 = (
    "both_arms_common_scale_identified_same_claim_tier"
)

EVALUATION_FAILURE_CODES = (
    "generator_failure",
    "preflight_nonidentified",
    "preflight_unsupported",
    "fit_exception",
    "optimizer_nonconverged",
    "nonfinite_objective",
    "gradient_threshold",
    "active_parameter_bound",
    "invalid_person_output",
    "invalid_se_output",
)


@dataclass(frozen=True, slots=True, kw_only=True)
class MonteCarloSpecV1:
    """Fixed-N formal Monte Carlo policy, separate from execution chunking."""

    requested_replicates: int
    randomization: RandomizationSpecV1 = RandomizationSpecV1()
    aggregation: str = "replicate_first_equal_weight"
    quantile_probs: tuple[float, ...] = (0.05, 0.50, 0.95)
    quantile_method: str = "linear_type7"
    stopping_rule: str = "fixed_n"
    partial_run_policy: str = "no_final_aggregate"
    schema_version: str = MONTE_CARLO_SPEC_VERSION

    def __post_init__(self) -> None:
        if self.schema_version != MONTE_CARLO_SPEC_VERSION:
            raise SimulationConditionValidationError(
                f"Unsupported Monte Carlo version: {self.schema_version!r}"
            )
        _require_int(
            "requested_replicates",
            self.requested_replicates,
            minimum=FORMAL_MONTE_CARLO_MIN_REPLICATES,
            maximum=MAX_MONTE_CARLO_REPLICATES,
        )
        if type(self.randomization) is not RandomizationSpecV1:
            raise TypeError("randomization must be a RandomizationSpecV1")
        _require_literal(
            "aggregation", self.aggregation, "replicate_first_equal_weight"
        )
        if self.quantile_probs != (0.05, 0.50, 0.95):
            raise SimulationConditionValidationError(
                "quantile_probs must remain (0.05, 0.50, 0.95) in v1"
            )
        _require_literal("quantile_method", self.quantile_method, "linear_type7")
        _require_literal("stopping_rule", self.stopping_rule, "fixed_n")
        _require_literal(
            "partial_run_policy", self.partial_run_policy, "no_final_aggregate"
        )

    def to_dict(self) -> dict:
        return _json_safe({
            "requested_replicates": self.requested_replicates,
            "randomization": self.randomization.to_dict(),
            "aggregation": self.aggregation,
            "quantile_probs": [float(value) for value in self.quantile_probs],
            "quantile_method": self.quantile_method,
            "stopping_rule": self.stopping_rule,
            "partial_run_policy": self.partial_run_policy,
            "schema_version": self.schema_version,
        })


def normalize_monte_carlo_spec(value: MonteCarloSpecV1 | Mapping) -> MonteCarloSpecV1:
    if type(value) is MonteCarloSpecV1:
        return value
    payload = _strict_mapping(value, MonteCarloSpecV1, label="Monte Carlo spec")
    payload["randomization"] = normalize_randomization_spec(payload["randomization"])
    payload["quantile_probs"] = _tuple_payload_field(
        payload, "quantile_probs", label="Monte Carlo spec"
    )
    return MonteCarloSpecV1(**payload)


def monte_carlo_spec_fingerprint(
    value: MonteCarloSpecV1 | Mapping,
    *,
    length: int = 16,
) -> str:
    return _fingerprint(normalize_monte_carlo_spec(value).to_dict(), length=length)


@dataclass(frozen=True, slots=True, kw_only=True)
class ClassificationRuleV1:
    """Truth-versus-estimate latent-scale classification rule."""

    cut_scores: tuple[float, ...]
    target_population: str = "study_person"
    tie_rule: str = "estimate_equal_cut_goes_upper"
    truth_tie_rule: str = "truth_equal_cut_goes_upper"
    metric: str = "classification_exact_accuracy"
    comparison_scale_alignment: str = (
        "declared_estimator_identification_to_generating_logit_v1"
    )
    schema_version: str = CLASSIFICATION_RULE_VERSION

    def __post_init__(self) -> None:
        if self.schema_version != CLASSIFICATION_RULE_VERSION:
            raise SimulationConditionValidationError(
                f"Unsupported classification rule version: {self.schema_version!r}"
            )
        if not isinstance(self.cut_scores, tuple) or not self.cut_scores:
            raise SimulationConditionValidationError(
                "cut_scores must be a non-empty tuple"
            )
        cuts = tuple(
            _require_real(f"cut_scores[{index}]", value)
            for index, value in enumerate(self.cut_scores)
        )
        if any(right <= left for left, right in zip(cuts, cuts[1:])):
            raise SimulationConditionValidationError(
                "cut_scores must be strictly increasing"
            )
        _require_literal("target_population", self.target_population, "study_person")
        _require_literal(
            "tie_rule", self.tie_rule, "estimate_equal_cut_goes_upper"
        )
        _require_literal(
            "truth_tie_rule", self.truth_tie_rule, "truth_equal_cut_goes_upper"
        )
        _require_literal("metric", self.metric, "classification_exact_accuracy")
        _require_literal(
            "comparison_scale_alignment",
            self.comparison_scale_alignment,
            "declared_estimator_identification_to_generating_logit_v1",
        )

    def to_dict(self) -> dict:
        return _json_safe({
            "cut_scores": [float(value) for value in self.cut_scores],
            "target_population": self.target_population,
            "tie_rule": self.tie_rule,
            "truth_tie_rule": self.truth_tie_rule,
            "metric": self.metric,
            "comparison_scale_alignment": self.comparison_scale_alignment,
            "schema_version": self.schema_version,
        })


def normalize_classification_rule(
    value: ClassificationRuleV1 | Mapping,
) -> ClassificationRuleV1:
    if type(value) is ClassificationRuleV1:
        return value
    payload = _strict_mapping(value, ClassificationRuleV1, label="classification rule")
    payload["cut_scores"] = _tuple_payload_field(
        payload, "cut_scores", label="classification rule"
    )
    return ClassificationRuleV1(**payload)


def classification_rule_fingerprint(
    value: ClassificationRuleV1 | Mapping,
    *,
    length: int = 16,
) -> str:
    return _fingerprint(normalize_classification_rule(value).to_dict(), length=length)


@dataclass(frozen=True, slots=True, kw_only=True)
class EvaluationPolicyV1:
    """Fixed success, pair eligibility/inclusion, and partial-run rules."""

    minimum_complete_pairs: int = FORMAL_MONTE_CARLO_MIN_REPLICATES
    minimum_complete_pair_rate: float = 0.80
    require_optimizer_success: bool = True
    require_finite_objective: bool = True
    require_finite_all_study_person_estimates: bool = True
    require_finite_all_study_person_se: bool = True
    max_gradient_norm: float | None = None
    gradient_norm_basis: str = "l2_unconstrained_free_parameter"
    active_bound_abs_tolerance: float = 1e-8
    parameter_bound_policy: str = "fail_if_active"
    retry_policy: str = "none"
    pair_inclusion: str = "both_arms_success_same_replicate"
    pair_comparability_requirement: str = PAIR_COMPARABILITY_REQUIREMENT_V1
    failure_denominator: str = "all_requested_replicates"
    partial_run_policy: str = "no_final_aggregate"
    schema_version: str = EVALUATION_POLICY_VERSION

    def __post_init__(self) -> None:
        if self.schema_version != EVALUATION_POLICY_VERSION:
            raise SimulationConditionValidationError(
                f"Unsupported evaluation policy version: {self.schema_version!r}"
            )
        _require_int(
            "minimum_complete_pairs",
            self.minimum_complete_pairs,
            minimum=FORMAL_MONTE_CARLO_MIN_REPLICATES,
            maximum=MAX_MONTE_CARLO_REPLICATES,
        )
        rate = _require_real(
            "minimum_complete_pair_rate", self.minimum_complete_pair_rate, minimum=0.0
        )
        if rate > 1.0:
            raise SimulationConditionValidationError(
                "minimum_complete_pair_rate must be <= 1"
            )
        for name in (
            "require_optimizer_success",
            "require_finite_objective",
            "require_finite_all_study_person_estimates",
            "require_finite_all_study_person_se",
        ):
            if getattr(self, name) is not True:
                raise SimulationConditionValidationError(
                    f"{name} must remain true in v1"
                )
        if self.max_gradient_norm is not None:
            _require_real(
                "max_gradient_norm", self.max_gradient_norm, strictly_positive=True
            )
        _require_literal(
            "gradient_norm_basis",
            self.gradient_norm_basis,
            "l2_unconstrained_free_parameter",
        )
        _require_real(
            "active_bound_abs_tolerance",
            self.active_bound_abs_tolerance,
            strictly_positive=True,
        )
        _require_literal(
            "parameter_bound_policy", self.parameter_bound_policy, "fail_if_active"
        )
        _require_literal("retry_policy", self.retry_policy, "none")
        _require_literal(
            "pair_inclusion",
            self.pair_inclusion,
            "both_arms_success_same_replicate",
        )
        _require_literal(
            "pair_comparability_requirement",
            self.pair_comparability_requirement,
            PAIR_COMPARABILITY_REQUIREMENT_V1,
        )
        _require_literal(
            "failure_denominator",
            self.failure_denominator,
            "all_requested_replicates",
        )
        _require_literal(
            "partial_run_policy", self.partial_run_policy, "no_final_aggregate"
        )

    def to_dict(self) -> dict:
        return _json_safe({
            "minimum_complete_pairs": self.minimum_complete_pairs,
            "minimum_complete_pair_rate": float(self.minimum_complete_pair_rate),
            "require_optimizer_success": self.require_optimizer_success,
            "require_finite_objective": self.require_finite_objective,
            "require_finite_all_study_person_estimates": (
                self.require_finite_all_study_person_estimates
            ),
            "require_finite_all_study_person_se": (
                self.require_finite_all_study_person_se
            ),
            "max_gradient_norm": (
                float(self.max_gradient_norm)
                if self.max_gradient_norm is not None
                else None
            ),
            "gradient_norm_basis": self.gradient_norm_basis,
            "active_bound_abs_tolerance": float(self.active_bound_abs_tolerance),
            "parameter_bound_policy": self.parameter_bound_policy,
            "retry_policy": self.retry_policy,
            "pair_inclusion": self.pair_inclusion,
            "pair_comparability_requirement": self.pair_comparability_requirement,
            "failure_denominator": self.failure_denominator,
            "partial_run_policy": self.partial_run_policy,
            "schema_version": self.schema_version,
        })


def normalize_evaluation_policy(
    value: EvaluationPolicyV1 | Mapping,
) -> EvaluationPolicyV1:
    if type(value) is EvaluationPolicyV1:
        return value
    return EvaluationPolicyV1(
        **_strict_mapping(value, EvaluationPolicyV1, label="evaluation policy")
    )


def evaluation_policy_fingerprint(
    value: EvaluationPolicyV1 | Mapping,
    *,
    length: int = 16,
) -> str:
    return _fingerprint(normalize_evaluation_policy(value).to_dict(), length=length)


__all__ = [
    "CLASSIFICATION_RULE_VERSION",
    "EVALUATION_FAILURE_CODES",
    "EVALUATION_POLICY_VERSION",
    "MONTE_CARLO_SPEC_VERSION",
    "PAIR_COMPARABILITY_REQUIREMENT_V1",
    "ClassificationRuleV1",
    "EvaluationPolicyV1",
    "MonteCarloSpecV1",
    "classification_rule_fingerprint",
    "evaluation_policy_fingerprint",
    "monte_carlo_spec_fingerprint",
    "normalize_classification_rule",
    "normalize_evaluation_policy",
    "normalize_monte_carlo_spec",
]
