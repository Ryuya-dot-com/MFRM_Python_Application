from __future__ import annotations

import numpy as np
import pytest

from mfrm_app import mml_stationarity as stationarity


TARGET = np.array([1.25, -0.75])
LOG_SIGMA_TARGET = np.log(1.7)


def quadratic_value(parameters: np.ndarray, sigma: float) -> float:
    structural = 0.5 * float(np.sum(np.square(parameters - TARGET)))
    scale = 0.5 * float((np.log(sigma) - LOG_SIGMA_TARGET) ** 2)
    return structural + scale


def quadratic_value_gradient(
    parameters: np.ndarray,
    sigma: float,
) -> tuple[float, np.ndarray]:
    return quadratic_value(parameters, sigma), parameters - TARGET


def test_joint_gradient_matches_central_difference_in_log_sigma_coordinates() -> None:
    objective, value_gradient = stationarity.make_joint_free_sd_functions(
        quadratic_value,
        quadratic_value_gradient,
    )
    point = np.array([0.2, 0.4, np.log(0.9)])
    audit = stationarity.audit_joint_gradient(objective, value_gradient, point)
    assert audit.coordinates == 3
    assert audit.maximum_absolute_difference < 1e-9
    assert audit.analytical_gradient == pytest.approx(
        [point[0] - TARGET[0], point[1] - TARGET[1], point[2] - LOG_SIGMA_TARGET],
        abs=1e-9,
    )


def test_joint_polish_recovers_structural_and_population_scale_optimum() -> None:
    result = stationarity.polish_joint_free_sd(
        structural_start=np.array([-2.0, 2.5]),
        sigma_start=0.4,
        value_function=quadratic_value,
        structural_value_gradient=quadratic_value_gradient,
        sigma_bounds=(0.05, 10.0),
    )
    assert result.optimizer_success is True
    assert result.structural_parameters == pytest.approx(TARGET, abs=1e-7)
    assert result.sigma == pytest.approx(np.exp(LOG_SIGMA_TARGET), abs=1e-7)
    assert result.final_projected_gradient_supnorm < 1e-7
    assert result.objective_improvement > 0
    assert result.objective_scale == "sum_negative_log_likelihood"


def test_supplied_joint_gradient_avoids_scalar_differences_and_validates_output() -> None:
    def unexpected(*_args):
        raise AssertionError("analytic evaluation must not call finite-difference helpers")

    point = np.array([0.2, 0.4, np.log(0.9)])
    expected = point - np.r_[TARGET, LOG_SIGMA_TARGET]
    _, gradient = stationarity.make_joint_free_sd_functions(
        unexpected, unexpected,
        joint_value_gradient=lambda x: (quadratic_value(x[:-1], np.exp(x[-1])), expected),
    )
    assert gradient(point)[1] == pytest.approx(expected, abs=1e-14)
    for value, bad, error in (
        (0.0, expected[:-1], ValueError),
        (0.0, np.full(3, np.nan), ValueError),
        (np.inf, expected, FloatingPointError),
    ):
        _, invalid = stationarity.make_joint_free_sd_functions(
            unexpected, unexpected, joint_value_gradient=lambda _x: (value, bad),
        )
        with pytest.raises(error):
            invalid(point)


def test_projected_gradient_uses_kkt_direction_at_active_bounds() -> None:
    projected = stationarity.projected_gradient(
        coordinates=np.array([0.0, 1.0, 0.5]),
        gradient=np.array([2.0, -3.0, 4.0]),
        bounds=[(0.0, None), (None, 1.0), (0.0, 1.0)],
    )
    assert projected.tolist() == [0.0, 0.0, 4.0]


def test_information_diagnostics_reconstruct_quadratic_curvature() -> None:
    _, value_gradient = stationarity.make_joint_free_sd_functions(
        quadratic_value,
        quadratic_value_gradient,
    )
    point = np.array([1.3, -0.7, np.log(1.8)])
    diagnostics = stationarity.information_diagnostics(value_gradient, point)
    expected_score = np.array(
        [point[0] - TARGET[0], point[1] - TARGET[1], point[2] - LOG_SIGMA_TARGET]
    )
    assert diagnostics.finite is True
    assert diagnostics.positive_definite is True
    assert diagnostics.minimum_eigenvalue == pytest.approx(1.0, abs=1e-6)
    assert diagnostics.maximum_eigenvalue == pytest.approx(1.0, abs=1e-6)
    assert diagnostics.condition_number == pytest.approx(1.0, abs=1e-6)
    assert diagnostics.standardized_score_supnorm == pytest.approx(
        np.max(np.abs(expected_score)), abs=1e-6
    )
    assert diagnostics.newton_correction_supnorm == pytest.approx(
        np.max(np.abs(expected_score)), abs=1e-6
    )


@pytest.mark.parametrize(
    "options",
    [
        stationarity.JointPolishOptions(maxiter=0),
        stationarity.JointPolishOptions(maxls=True),
        stationarity.JointPolishOptions(gtol=0),
        stationarity.JointPolishOptions(ftol=np.nan),
        stationarity.JointPolishOptions(log_sigma_relative_step=-1),
    ],
)
def test_invalid_polish_options_fail_closed(
    options: stationarity.JointPolishOptions,
) -> None:
    with pytest.raises(ValueError):
        options.validate()
