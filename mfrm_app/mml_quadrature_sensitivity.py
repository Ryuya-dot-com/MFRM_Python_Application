"""Typed Q-fit/Q-check numerical sensitivity for free-SD MML.

This module never accepts a quadrature PASS boolean and never emits scientific
``InferenceReady``.  It recomputes both stationarity assessments from their raw
runs, evaluates each solution under both finite-quadrature objectives, and
applies an explicit numerical sensitivity contract.  Artifact registration and
scientific inference remain separate future study-level responsibilities.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import operator
import re

import numpy as np

from mfrm_app.mml_engine_v2 import (
    FreeSdStationarityRun,
    StationarityAssessment,
    StationarityContract,
    assess_free_sd_stationarity,
)
from mfrm_app.mml_stationarity import ValueFunction


SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")


@dataclass(frozen=True)
class QuadratureSensitivityContract:
    """Unregistered numerical tolerances; not a scientific contract."""

    primary_quadrature_points: int
    sensitivity_quadrature_points: int
    max_structural_parameter_difference: float
    max_log_sigma_difference: float
    max_sensitivity_optimization_gain_per_observation: float
    max_negative_sensitivity_gain_per_observation: float
    max_abs_primary_back_evaluation_change_per_observation: float
    max_objective_reconstruction_disagreement: float

    def validate(self) -> None:
        primary = _quadrature_points(self.primary_quadrature_points, "primary")
        sensitivity = _quadrature_points(
            self.sensitivity_quadrature_points,
            "sensitivity",
        )
        if sensitivity <= primary:
            raise ValueError("sensitivity quadrature must exceed primary quadrature")
        for name in (
            "max_structural_parameter_difference",
            "max_log_sigma_difference",
            "max_sensitivity_optimization_gain_per_observation",
            "max_negative_sensitivity_gain_per_observation",
            "max_abs_primary_back_evaluation_change_per_observation",
            "max_objective_reconstruction_disagreement",
        ):
            value = float(getattr(self, name))
            if not np.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be positive and finite")


@dataclass(frozen=True)
class QuadratureSensitivityAssessment:
    problem_digest: str
    primary_quadrature_points: int
    sensitivity_quadrature_points: int
    observations: float
    primary_run: FreeSdStationarityRun
    sensitivity_run: FreeSdStationarityRun
    primary_stationarity: StationarityAssessment
    sensitivity_stationarity: StationarityAssessment
    primary_objective_at_primary_solution: float
    primary_objective_at_sensitivity_solution: float
    sensitivity_objective_at_primary_solution: float
    sensitivity_objective_at_sensitivity_solution: float
    primary_objective_reconstruction_difference: float
    sensitivity_objective_reconstruction_difference: float
    maximum_structural_parameter_difference: float
    absolute_log_sigma_difference: float
    sensitivity_optimization_gain_per_observation: float
    primary_back_evaluation_change_per_observation: float
    finite_cross_evaluation: bool
    primary_stationarity_pass: bool
    sensitivity_stationarity_pass: bool
    structural_difference_pass: bool
    log_sigma_difference_pass: bool
    sensitivity_gain_pass: bool
    primary_back_evaluation_pass: bool
    objective_reconstruction_pass: bool
    numerical_sensitivity_pass: bool
    scientific_inference_ready: bool
    status: str
    contract: dict[str, object]

    def to_dict(self) -> dict[str, object]:
        result = asdict(self)
        result.update(
            {
                "NumericalSensitivityPass": self.numerical_sensitivity_pass,
                "ScientificInferenceReady": self.scientific_inference_ready,
            }
        )
        return result


def _quadrature_points(value: object, label: str) -> int:
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{label} quadrature points must be an integer")
    try:
        result = operator.index(value)
    except TypeError as exc:
        raise ValueError(f"{label} quadrature points must be an integer") from exc
    if result < 3:
        raise ValueError(f"{label} quadrature points must be at least 3")
    return int(result)


def _problem_digest(value: str) -> str:
    result = str(value).strip().lower()
    if not SHA256_PATTERN.fullmatch(result):
        raise ValueError("problem_digest must be a lowercase SHA256 digest")
    return result


def _observations(primary: FreeSdStationarityRun, sensitivity: FreeSdStationarityRun) -> float:
    left = float(primary.observations)
    right = float(sensitivity.observations)
    if (
        not np.isfinite(left)
        or not np.isfinite(right)
        or left <= 0
        or right <= 0
        or abs(left - right) > 1e-12 * max(1.0, abs(left), abs(right))
    ):
        raise ValueError("primary and sensitivity observation counts must agree")
    return left


def _solution(run: FreeSdStationarityRun) -> tuple[np.ndarray, float]:
    structural = np.asarray(run.restart_polish.structural_parameters, dtype=float)
    sigma = float(run.restart_polish.sigma)
    if (
        structural.ndim != 1
        or structural.size == 0
        or not np.isfinite(structural).all()
        or not np.isfinite(sigma)
        or sigma <= 0
    ):
        raise ValueError("quadrature solution is nonfinite or malformed")
    return structural, sigma


def assess_quadrature_sensitivity(
    *,
    problem_digest: str,
    primary_quadrature_points: int,
    sensitivity_quadrature_points: int,
    primary_run: FreeSdStationarityRun,
    sensitivity_run: FreeSdStationarityRun,
    primary_value_function: ValueFunction,
    sensitivity_value_function: ValueFunction,
    stationarity_contract: StationarityContract,
    sensitivity_contract: QuadratureSensitivityContract,
) -> QuadratureSensitivityAssessment:
    """Recompute Q-fit/Q-check evidence without a boolean promotion path."""

    digest = _problem_digest(problem_digest)
    primary_points = _quadrature_points(primary_quadrature_points, "primary")
    sensitivity_points = _quadrature_points(
        sensitivity_quadrature_points,
        "sensitivity",
    )
    sensitivity_contract.validate()
    if (
        primary_points != sensitivity_contract.primary_quadrature_points
        or sensitivity_points != sensitivity_contract.sensitivity_quadrature_points
    ):
        raise ValueError("quadrature evidence differs from the numerical contract")
    observations = _observations(primary_run, sensitivity_run)
    primary_stationarity = assess_free_sd_stationarity(
        primary_run,
        stationarity_contract,
    )
    sensitivity_stationarity = assess_free_sd_stationarity(
        sensitivity_run,
        stationarity_contract,
    )
    primary_parameters, primary_sigma = _solution(primary_run)
    sensitivity_parameters, sensitivity_sigma = _solution(sensitivity_run)
    if primary_parameters.shape != sensitivity_parameters.shape:
        raise ValueError("primary and sensitivity structural dimensions differ")

    primary_at_primary = float(
        primary_value_function(primary_parameters, primary_sigma)
    )
    primary_at_sensitivity = float(
        primary_value_function(sensitivity_parameters, sensitivity_sigma)
    )
    sensitivity_at_primary = float(
        sensitivity_value_function(primary_parameters, primary_sigma)
    )
    sensitivity_at_sensitivity = float(
        sensitivity_value_function(sensitivity_parameters, sensitivity_sigma)
    )
    values = np.array(
        [
            primary_at_primary,
            primary_at_sensitivity,
            sensitivity_at_primary,
            sensitivity_at_sensitivity,
        ],
        dtype=float,
    )
    finite = bool(np.isfinite(values).all())
    primary_reconstruction = abs(
        primary_at_primary - float(primary_run.final_objective)
    )
    sensitivity_reconstruction = abs(
        sensitivity_at_sensitivity - float(sensitivity_run.final_objective)
    )
    structural_difference = float(
        np.max(np.abs(primary_parameters - sensitivity_parameters))
    )
    log_sigma_difference = abs(float(np.log(primary_sigma) - np.log(sensitivity_sigma)))
    sensitivity_gain = float(
        (sensitivity_at_primary - sensitivity_at_sensitivity) / observations
    )
    primary_back_change = float(
        (primary_at_sensitivity - primary_at_primary) / observations
    )
    reconstruction_pass = bool(
        finite
        and primary_reconstruction
        <= sensitivity_contract.max_objective_reconstruction_disagreement
        and sensitivity_reconstruction
        <= sensitivity_contract.max_objective_reconstruction_disagreement
    )
    structural_pass = bool(
        np.isfinite(structural_difference)
        and structural_difference
        <= sensitivity_contract.max_structural_parameter_difference
    )
    sigma_pass = bool(
        np.isfinite(log_sigma_difference)
        and log_sigma_difference <= sensitivity_contract.max_log_sigma_difference
    )
    gain_pass = bool(
        np.isfinite(sensitivity_gain)
        and sensitivity_gain
        >= -sensitivity_contract.max_negative_sensitivity_gain_per_observation
        and sensitivity_gain
        <= sensitivity_contract.max_sensitivity_optimization_gain_per_observation
    )
    back_pass = bool(
        np.isfinite(primary_back_change)
        and abs(primary_back_change)
        <= sensitivity_contract.max_abs_primary_back_evaluation_change_per_observation
    )
    numerical_pass = all(
        (
            finite,
            primary_stationarity.stationarity_pass,
            sensitivity_stationarity.stationarity_pass,
            structural_pass,
            sigma_pass,
            gain_pass,
            back_pass,
            reconstruction_pass,
        )
    )
    return QuadratureSensitivityAssessment(
        problem_digest=digest,
        primary_quadrature_points=primary_points,
        sensitivity_quadrature_points=sensitivity_points,
        observations=observations,
        primary_run=primary_run,
        sensitivity_run=sensitivity_run,
        primary_stationarity=primary_stationarity,
        sensitivity_stationarity=sensitivity_stationarity,
        primary_objective_at_primary_solution=primary_at_primary,
        primary_objective_at_sensitivity_solution=primary_at_sensitivity,
        sensitivity_objective_at_primary_solution=sensitivity_at_primary,
        sensitivity_objective_at_sensitivity_solution=sensitivity_at_sensitivity,
        primary_objective_reconstruction_difference=primary_reconstruction,
        sensitivity_objective_reconstruction_difference=sensitivity_reconstruction,
        maximum_structural_parameter_difference=structural_difference,
        absolute_log_sigma_difference=log_sigma_difference,
        sensitivity_optimization_gain_per_observation=sensitivity_gain,
        primary_back_evaluation_change_per_observation=primary_back_change,
        finite_cross_evaluation=finite,
        primary_stationarity_pass=primary_stationarity.stationarity_pass,
        sensitivity_stationarity_pass=sensitivity_stationarity.stationarity_pass,
        structural_difference_pass=structural_pass,
        log_sigma_difference_pass=sigma_pass,
        sensitivity_gain_pass=gain_pass,
        primary_back_evaluation_pass=back_pass,
        objective_reconstruction_pass=reconstruction_pass,
        numerical_sensitivity_pass=numerical_pass,
        scientific_inference_ready=False,
        status=(
            "NUMERICAL_SENSITIVITY_PASS_REGISTRATION_REQUIRED"
            if numerical_pass
            else "NUMERICAL_SENSITIVITY_NOT_QUALIFIED"
        ),
        contract=asdict(sensitivity_contract),
    )


__all__ = [
    "QuadratureSensitivityAssessment",
    "QuadratureSensitivityContract",
    "assess_quadrature_sensitivity",
]
