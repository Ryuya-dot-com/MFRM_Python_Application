"""Composite scientific-scenario contract."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
import math

from ._condition_estimator import (
    ESTIMATOR_MML,
    EstimatorSpecV1,
    normalize_estimator_spec,
)
from ._condition_evaluation import (
    ClassificationRuleV1,
    EvaluationPolicyV1,
    MonteCarloSpecV1,
    normalize_classification_rule,
    normalize_evaluation_policy,
    normalize_monte_carlo_spec,
)
from ._condition_shared import (
    SimulationConditionValidationError,
    _fingerprint,
    _json_safe,
    _strict_mapping,
)
from ._condition_truth import RsmTruthSpecV1, normalize_rsm_truth_spec


SIMULATION_SCENARIO_SPEC_VERSION = "mfrm_simulation_scenario_spec_v1"


@dataclass(frozen=True, slots=True, kw_only=True)
class SimulationScenarioSpecV1:
    """One estimator-specific scientific scenario, excluding design and cost."""

    truth: RsmTruthSpecV1
    estimator: EstimatorSpecV1
    monte_carlo: MonteCarloSpecV1
    evaluation_policy: EvaluationPolicyV1 = EvaluationPolicyV1()
    classification_rule: ClassificationRuleV1 | None = None
    schema_version: str = SIMULATION_SCENARIO_SPEC_VERSION

    def __post_init__(self) -> None:
        if self.schema_version != SIMULATION_SCENARIO_SPEC_VERSION:
            raise SimulationConditionValidationError(
                f"Unsupported simulation scenario version: {self.schema_version!r}"
            )
        if type(self.truth) is not RsmTruthSpecV1:
            raise TypeError("truth must be a RsmTruthSpecV1")
        if type(self.estimator) is not EstimatorSpecV1:
            raise TypeError("estimator must be an EstimatorSpecV1")
        if type(self.monte_carlo) is not MonteCarloSpecV1:
            raise TypeError("monte_carlo must be a MonteCarloSpecV1")
        if type(self.evaluation_policy) is not EvaluationPolicyV1:
            raise TypeError("evaluation_policy must be an EvaluationPolicyV1")
        if self.classification_rule is not None and (
            type(self.classification_rule) is not ClassificationRuleV1
        ):
            raise TypeError(
                "classification_rule must be a ClassificationRuleV1 or None"
            )
        if (
            self.evaluation_policy.minimum_complete_pairs
            > self.monte_carlo.requested_replicates
        ):
            raise SimulationConditionValidationError(
                "minimum_complete_pairs cannot exceed requested_replicates"
            )
        if self.estimator.method == ESTIMATOR_MML:
            if not math.isclose(
                float(self.estimator.population_prior_mean),
                float(self.truth.person_mean),
                rel_tol=0.0,
                abs_tol=1e-12,
            ):
                raise SimulationConditionValidationError(
                    "v1 MML population prior mean must match the truth person mean"
                )
            if not math.isclose(
                float(self.estimator.population_prior_sd),
                float(self.truth.person_sd),
                rel_tol=0.0,
                abs_tol=1e-12,
            ):
                raise SimulationConditionValidationError(
                    "v1 MML population prior SD must match the truth person SD"
                )

    def to_dict(self) -> dict:
        return _json_safe({
            "truth": self.truth.to_dict(),
            "estimator": self.estimator.to_dict(),
            "monte_carlo": self.monte_carlo.to_dict(),
            "evaluation_policy": self.evaluation_policy.to_dict(),
            "classification_rule": (
                self.classification_rule.to_dict()
                if self.classification_rule is not None
                else None
            ),
            "schema_version": self.schema_version,
        })


def normalize_simulation_scenario_spec(
    value: SimulationScenarioSpecV1 | Mapping,
) -> SimulationScenarioSpecV1:
    if type(value) is SimulationScenarioSpecV1:
        return value
    payload = _strict_mapping(
        value, SimulationScenarioSpecV1, label="simulation scenario spec"
    )
    payload["truth"] = normalize_rsm_truth_spec(payload["truth"])
    payload["estimator"] = normalize_estimator_spec(payload["estimator"])
    payload["monte_carlo"] = normalize_monte_carlo_spec(payload["monte_carlo"])
    payload["evaluation_policy"] = normalize_evaluation_policy(
        payload["evaluation_policy"]
    )
    if payload["classification_rule"] is not None:
        payload["classification_rule"] = normalize_classification_rule(
            payload["classification_rule"]
        )
    return SimulationScenarioSpecV1(**payload)


def simulation_scenario_spec_fingerprint(
    value: SimulationScenarioSpecV1 | Mapping,
    *,
    length: int = 16,
) -> str:
    return _fingerprint(
        normalize_simulation_scenario_spec(value).to_dict(), length=length
    )


__all__ = [
    "SIMULATION_SCENARIO_SPEC_VERSION",
    "SimulationScenarioSpecV1",
    "normalize_simulation_scenario_spec",
    "simulation_scenario_spec_fingerprint",
]
