from __future__ import annotations

from dataclasses import replace
import hashlib

import numpy as np
import pytest

from mfrm_app.mml_engine_v2 import StationarityContract, run_free_sd_stationarity_v2
from mfrm_app.mml_quadrature_sensitivity import (
    QuadratureSensitivityContract,
    assess_quadrature_sensitivity,
)


PRIMARY_TARGET = np.array([0.5, -0.25])
SENSITIVITY_TARGET = np.array([0.5002, -0.2501])
PRIMARY_LOG_SIGMA = np.log(1.2)
SENSITIVITY_LOG_SIGMA = np.log(1.2001)
PROBLEM_DIGEST = hashlib.sha256(b"typed-qfit-qcheck-test-problem").hexdigest()


def objective(target, log_sigma_target):
    def value(parameters: np.ndarray, sigma: float) -> float:
        return 0.5 * float(
            np.sum(np.square(parameters - target))
            + (np.log(sigma) - log_sigma_target) ** 2
        )

    def value_gradient(parameters: np.ndarray, sigma: float):
        return value(parameters, sigma), parameters - target

    return value, value_gradient


PRIMARY_VALUE, PRIMARY_GRADIENT = objective(PRIMARY_TARGET, PRIMARY_LOG_SIGMA)
SENSITIVITY_VALUE, SENSITIVITY_GRADIENT = objective(
    SENSITIVITY_TARGET,
    SENSITIVITY_LOG_SIGMA,
)


def stationarity_contract() -> StationarityContract:
    return StationarityContract(
        max_projected_gradient_supnorm=1e-6,
        max_standardized_score_supnorm=1e-6,
        max_newton_correction_supnorm=1e-6,
        max_restart_improvement_per_observation=1e-10,
        max_restart_displacement=1e-5,
        max_gradient_fd_disagreement=1e-7,
        max_objective_value_disagreement=1e-10,
        max_information_relative_symmetry_residual=1e-6,
        max_information_condition_number=1e6,
        max_constraint_residual=1e-10,
        objective_worsening_tolerance_total=1e-12,
        objective_worsening_tolerance_per_observation=1e-14,
        sigma_boundary_tolerance=1e-6,
    )


def sensitivity_contract(**overrides) -> QuadratureSensitivityContract:
    values = {
        "primary_quadrature_points": 31,
        "sensitivity_quadrature_points": 61,
        "max_structural_parameter_difference": 1e-3,
        "max_log_sigma_difference": 1e-3,
        "max_sensitivity_optimization_gain_per_observation": 1e-6,
        "max_negative_sensitivity_gain_per_observation": 1e-10,
        "max_abs_primary_back_evaluation_change_per_observation": 1e-6,
        "max_objective_reconstruction_disagreement": 1e-9,
    }
    values.update(overrides)
    return QuadratureSensitivityContract(**values)


@pytest.fixture(scope="module")
def runs():
    primary = run_free_sd_stationarity_v2(
        np.array([-1.0, 1.0]),
        0.8,
        PRIMARY_VALUE,
        PRIMARY_GRADIENT,
        observations=100,
        constraint_residual_function=lambda _parameters: 0.0,
    )
    sensitivity = run_free_sd_stationarity_v2(
        primary.restart_polish.structural_parameters,
        primary.restart_polish.sigma,
        SENSITIVITY_VALUE,
        SENSITIVITY_GRADIENT,
        observations=100,
        constraint_residual_function=lambda _parameters: 0.0,
    )
    return primary, sensitivity


def assess(runs, **kwargs):
    primary, sensitivity = runs
    values = {
        "problem_digest": PROBLEM_DIGEST,
        "primary_quadrature_points": 31,
        "sensitivity_quadrature_points": 61,
        "primary_run": primary,
        "sensitivity_run": sensitivity,
        "primary_value_function": PRIMARY_VALUE,
        "sensitivity_value_function": SENSITIVITY_VALUE,
        "stationarity_contract": stationarity_contract(),
        "sensitivity_contract": sensitivity_contract(),
    }
    values.update(kwargs)
    return assess_quadrature_sensitivity(**values)


def test_qfit_qcheck_is_recomputed_but_cannot_emit_scientific_readiness(runs) -> None:
    result = assess(runs)
    assert result.primary_stationarity_pass is True
    assert result.sensitivity_stationarity_pass is True
    assert result.numerical_sensitivity_pass is True
    assert result.scientific_inference_ready is False
    assert result.status == "NUMERICAL_SENSITIVITY_PASS_REGISTRATION_REQUIRED"
    assert result.to_dict()["ScientificInferenceReady"] is False
    assert result.maximum_structural_parameter_difference == pytest.approx(2e-4, abs=1e-7)


def test_tight_sensitivity_threshold_fails_without_hiding_stationarity(runs) -> None:
    result = assess(
        runs,
        sensitivity_contract=sensitivity_contract(
            max_structural_parameter_difference=1e-8
        ),
    )
    assert result.primary_stationarity_pass is True
    assert result.sensitivity_stationarity_pass is True
    assert result.structural_difference_pass is False
    assert result.numerical_sensitivity_pass is False
    assert result.scientific_inference_ready is False


@pytest.mark.parametrize("bad_points", [True, 31.5, "31"])
def test_noninteger_quadrature_evidence_is_rejected(runs, bad_points) -> None:
    with pytest.raises(ValueError, match="must be an integer"):
        assess(runs, primary_quadrature_points=bad_points)


def test_contract_and_evidence_quadrature_mismatch_is_rejected(runs) -> None:
    with pytest.raises(ValueError, match="differs from the numerical contract"):
        assess(runs, sensitivity_quadrature_points=81)


def test_problem_digest_is_not_a_free_form_label(runs) -> None:
    with pytest.raises(ValueError, match="lowercase SHA256"):
        assess(runs, problem_digest="same-problem-trust-me")


def test_stationarity_failure_cannot_be_replaced_by_quadrature_agreement(runs) -> None:
    primary, sensitivity = runs
    broken_information = replace(
        sensitivity.information,
        positive_definite=False,
        condition_number=np.inf,
        newton_correction_supnorm=np.inf,
    )
    broken = replace(sensitivity, information=broken_information)
    result = assess(runs, sensitivity_run=broken)
    assert result.sensitivity_stationarity_pass is False
    assert result.numerical_sensitivity_pass is False
    assert result.scientific_inference_ready is False
