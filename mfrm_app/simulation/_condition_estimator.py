"""Versioned estimator conditions for prospective MFRM simulation."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

from ._condition_shared import (
    MAX_ESTIMATOR_ITERATIONS,
    MAX_GENERATING_STANDARD_DEVIATION,
    MODEL_RSM,
    SimulationConditionValidationError,
    _fingerprint,
    _json_safe,
    _require_int,
    _require_literal,
    _require_real,
    _strict_mapping,
)


ESTIMATOR_SPEC_VERSION = "mfrm_estimator_spec_v1"

ESTIMATOR_JMLE = "JMLE"
ESTIMATOR_MML = "MML"
ESTIMATOR_METHOD_CHOICES = (ESTIMATOR_JMLE, ESTIMATOR_MML)


@dataclass(frozen=True, slots=True, kw_only=True)
class EstimatorSpecV1:
    """No-fallback adapter settings for one Python JMLE or MML fit."""

    method: str
    max_iterations: int = 400
    relative_tolerance: float = 1e-6
    model: str = MODEL_RSM
    person_score: str = "MLE"
    reported_person_uncertainty: str = (
        "conditional_element_information_se_other_parameters_fixed"
    )
    identification_policy: str = "center_rater_and_criterion_person_free"
    optimizer: str = "L-BFGS-B"
    mml_engine: str | None = None
    quadrature_nodes: int | None = None
    population_prior_mean: float | None = None
    population_prior_sd: float | None = None
    population_prior_distribution: str | None = None
    estimate_population_sd: bool = False
    rater_constraint: str = "sum_to_zero"
    criterion_constraint: str = "sum_to_zero"
    threshold_constraint: str = "sum_to_zero"
    regularization: str = "none"
    fallback: str = "none"
    fit_adapter_version: str = "python_mfrm_sim_adapter_v1"
    schema_version: str = ESTIMATOR_SPEC_VERSION

    def __post_init__(self) -> None:
        if self.schema_version != ESTIMATOR_SPEC_VERSION:
            raise SimulationConditionValidationError(
                f"Unsupported estimator spec version: {self.schema_version!r}"
            )
        if self.method not in ESTIMATOR_METHOD_CHOICES:
            raise SimulationConditionValidationError(
                f"method must be one of {ESTIMATOR_METHOD_CHOICES}"
            )
        _require_int(
            "max_iterations",
            self.max_iterations,
            minimum=1,
            maximum=MAX_ESTIMATOR_ITERATIONS,
        )
        tolerance = _require_real(
            "relative_tolerance", self.relative_tolerance, strictly_positive=True
        )
        if tolerance > 1.0:
            raise SimulationConditionValidationError(
                "relative_tolerance must be <= 1"
            )
        if not isinstance(self.estimate_population_sd, bool):
            raise SimulationConditionValidationError(
                "estimate_population_sd must be boolean"
            )
        if self.estimate_population_sd:
            raise SimulationConditionValidationError(
                "v1 simulation estimators require a fixed population SD"
            )
        _require_literal("model", self.model, MODEL_RSM)
        _require_literal("rater_constraint", self.rater_constraint, "sum_to_zero")
        _require_literal(
            "criterion_constraint", self.criterion_constraint, "sum_to_zero"
        )
        _require_literal(
            "threshold_constraint", self.threshold_constraint, "sum_to_zero"
        )
        _require_literal("regularization", self.regularization, "none")
        _require_literal("fallback", self.fallback, "none")
        _require_literal(
            "fit_adapter_version",
            self.fit_adapter_version,
            "python_mfrm_sim_adapter_v1",
        )

        if self.method == ESTIMATOR_JMLE:
            _require_literal("person_score", self.person_score, "MLE")
            _require_literal(
                "reported_person_uncertainty",
                self.reported_person_uncertainty,
                "conditional_element_information_se_other_parameters_fixed",
            )
            _require_literal(
                "identification_policy",
                self.identification_policy,
                "center_rater_and_criterion_person_free",
            )
            _require_literal("optimizer", self.optimizer, "L-BFGS-B")
            for name in (
                "mml_engine",
                "quadrature_nodes",
                "population_prior_mean",
                "population_prior_sd",
                "population_prior_distribution",
            ):
                if getattr(self, name) is not None:
                    raise SimulationConditionValidationError(
                        f"{name} must be None for JMLE"
                    )
            return

        _require_literal("person_score", self.person_score, "EAP")
        _require_literal(
            "reported_person_uncertainty",
            self.reported_person_uncertainty,
            "eap_posterior_standard_deviation",
        )
        _require_literal(
            "identification_policy",
            self.identification_policy,
            "fixed_zero_person_population_mean_center_rater_and_criterion",
        )
        _require_literal("optimizer", self.optimizer, "EM")
        _require_literal("mml_engine", self.mml_engine, "EM")
        nodes = _require_int(
            "quadrature_nodes", self.quadrature_nodes, minimum=7, maximum=101
        )
        if nodes % 2 == 0:
            raise SimulationConditionValidationError(
                "quadrature_nodes must be odd for the v1 MML adapter"
            )
        prior_mean = _require_real("population_prior_mean", self.population_prior_mean)
        if prior_mean != 0.0:
            raise SimulationConditionValidationError(
                "population_prior_mean must remain 0.0 in v1"
            )
        _require_real(
            "population_prior_sd",
            self.population_prior_sd,
            maximum=MAX_GENERATING_STANDARD_DEVIATION,
            strictly_positive=True,
        )
        _require_literal(
            "population_prior_distribution",
            self.population_prior_distribution,
            "normal",
        )

    def to_dict(self) -> dict:
        return _json_safe({
            "method": self.method,
            "max_iterations": self.max_iterations,
            "relative_tolerance": float(self.relative_tolerance),
            "model": self.model,
            "person_score": self.person_score,
            "reported_person_uncertainty": self.reported_person_uncertainty,
            "identification_policy": self.identification_policy,
            "optimizer": self.optimizer,
            "mml_engine": self.mml_engine,
            "quadrature_nodes": self.quadrature_nodes,
            "population_prior_mean": (
                float(self.population_prior_mean)
                if self.population_prior_mean is not None
                else None
            ),
            "population_prior_sd": (
                float(self.population_prior_sd)
                if self.population_prior_sd is not None
                else None
            ),
            "population_prior_distribution": self.population_prior_distribution,
            "estimate_population_sd": self.estimate_population_sd,
            "rater_constraint": self.rater_constraint,
            "criterion_constraint": self.criterion_constraint,
            "threshold_constraint": self.threshold_constraint,
            "regularization": self.regularization,
            "fallback": self.fallback,
            "fit_adapter_version": self.fit_adapter_version,
            "schema_version": self.schema_version,
        })


def jmle_estimator_spec(
    *,
    max_iterations: int = 400,
    relative_tolerance: float = 1e-6,
) -> EstimatorSpecV1:
    return EstimatorSpecV1(
        method=ESTIMATOR_JMLE,
        max_iterations=max_iterations,
        relative_tolerance=relative_tolerance,
    )


def mml_estimator_spec(
    *,
    max_iterations: int = 400,
    relative_tolerance: float = 1e-6,
    quadrature_nodes: int = 15,
    population_prior_sd: float = 1.0,
) -> EstimatorSpecV1:
    return EstimatorSpecV1(
        method=ESTIMATOR_MML,
        max_iterations=max_iterations,
        relative_tolerance=relative_tolerance,
        person_score="EAP",
        reported_person_uncertainty="eap_posterior_standard_deviation",
        identification_policy=(
            "fixed_zero_person_population_mean_center_rater_and_criterion"
        ),
        optimizer="EM",
        mml_engine="EM",
        quadrature_nodes=quadrature_nodes,
        population_prior_mean=0.0,
        population_prior_sd=population_prior_sd,
        population_prior_distribution="normal",
    )


def normalize_estimator_spec(value: EstimatorSpecV1 | Mapping) -> EstimatorSpecV1:
    if type(value) is EstimatorSpecV1:
        return value
    return EstimatorSpecV1(
        **_strict_mapping(value, EstimatorSpecV1, label="estimator spec")
    )


def estimator_spec_fingerprint(
    value: EstimatorSpecV1 | Mapping,
    *,
    length: int = 16,
) -> str:
    return _fingerprint(normalize_estimator_spec(value).to_dict(), length=length)


__all__ = [
    "ESTIMATOR_JMLE",
    "ESTIMATOR_METHOD_CHOICES",
    "ESTIMATOR_MML",
    "ESTIMATOR_SPEC_VERSION",
    "EstimatorSpecV1",
    "estimator_spec_fingerprint",
    "jmle_estimator_spec",
    "mml_estimator_spec",
    "normalize_estimator_spec",
]
