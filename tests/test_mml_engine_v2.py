from __future__ import annotations

from dataclasses import replace
from types import SimpleNamespace

import numpy as np
import pytest

from mfrm_app.mml_engine_v2 import (
    StationarityContract,
    assess_free_sd_stationarity,
    run_free_sd_stationarity_v2,
    run_free_sd_two_stage,
)
from mfrm_app import mml_stationarity


TARGET = np.array([1.25, -0.75])
LOG_SIGMA_TARGET = np.log(1.7)


def value(parameters: np.ndarray, sigma: float) -> float:
    return 0.5 * float(
        np.sum(np.square(parameters - TARGET))
        + (np.log(sigma) - LOG_SIGMA_TARGET) ** 2
    )


def value_gradient(parameters: np.ndarray, sigma: float) -> tuple[float, np.ndarray]:
    return value(parameters, sigma), parameters - TARGET


def contract(**overrides: float) -> StationarityContract:
    values = {
        "max_projected_gradient_supnorm": 1e-6,
        "max_standardized_score_supnorm": 1e-6,
        "max_newton_correction_supnorm": 1e-6,
        "max_restart_improvement_per_observation": 1e-10,
        "max_restart_displacement": 1e-5,
        "max_gradient_fd_disagreement": 1e-7,
        "max_objective_value_disagreement": 1e-10,
        "max_information_relative_symmetry_residual": 1e-6,
        "max_information_condition_number": 1e6,
        "max_constraint_residual": 1e-10,
        "objective_worsening_tolerance_total": 1e-12,
        "objective_worsening_tolerance_per_observation": 1e-14,
        "sigma_boundary_tolerance": 1e-6,
    }
    values.update(overrides)
    return StationarityContract(**values)


@pytest.fixture(scope="module")
def completed_run():
    return run_free_sd_stationarity_v2(
        structural_start=np.array([-2.0, 2.5]),
        sigma_start=0.4,
        value_function=value,
        structural_value_gradient=value_gradient,
        observations=100,
        constraint_residual_function=lambda _parameters: 0.0,
    )


def test_two_polish_run_retains_threshold_free_diagnostics(completed_run) -> None:
    run = completed_run
    assert run.algorithm_terminated is True
    assert run.finite_solution is True
    assert run.objective_nonworsening is True
    assert run.final_sigma == pytest.approx(1.7, abs=1e-7)
    assert run.final_projected_gradient_supnorm < 1e-7
    assert run.information.positive_definite is True
    assert run.restart_improvement_per_observation <= 1e-12


def test_stationarity_cannot_self_promote_to_inference_ready(completed_run) -> None:
    incomplete = assess_free_sd_stationarity(
        completed_run,
        contract(),
    )
    assert incomplete.stationarity_pass is True
    assert incomplete.inference_ready is False
    assert incomplete.quadrature_sensitivity_pass is None
    assert incomplete.status == "STATIONARITY_PASS_QUADRATURE_EVIDENCE_REQUIRED"
    assert incomplete.to_dict()["InferenceReady"] is False
    with pytest.raises(TypeError):
        assess_free_sd_stationarity(  # type: ignore[call-arg]
            completed_run,
            contract(),
            quadrature_sensitivity_pass=True,
        )


def test_optimizer_success_cannot_hide_a_large_terminal_score(completed_run) -> None:
    assessment = assess_free_sd_stationarity(
        completed_run,
        contract(max_projected_gradient_supnorm=1e-20),
    )
    assert assessment.algorithm_terminated is True
    assert assessment.projected_gradient_pass is False
    assert assessment.stationarity_pass is False
    assert assessment.inference_ready is False


def test_nonpositive_information_is_fail_closed(completed_run) -> None:
    broken_information = replace(
        completed_run.information,
        positive_definite=False,
        condition_number=np.inf,
        newton_correction_supnorm=np.inf,
    )
    broken = replace(completed_run, information=broken_information)
    assessment = assess_free_sd_stationarity(
        broken,
        contract(),
    )
    assert assessment.information_positive_definite is False
    assert assessment.stationarity_pass is False
    assert assessment.inference_ready is False


def test_asymmetric_information_is_fail_closed(completed_run) -> None:
    asymmetric = replace(
        completed_run.information,
        relative_symmetry_residual=0.1,
    )
    assessment = assess_free_sd_stationarity(
        replace(completed_run, information=asymmetric),
        contract(),
    )
    assert assessment.information_symmetry_pass is False
    assert assessment.stationarity_pass is False


@pytest.mark.parametrize(
    "field,value",
    [
        ("max_projected_gradient_supnorm", 0.0),
        ("max_information_condition_number", np.inf),
        ("max_constraint_residual", 0.0),
        ("objective_worsening_tolerance_total", -1.0),
        ("sigma_boundary_tolerance", np.nan),
    ],
)
def test_invalid_contract_fails_before_assessment(
    completed_run,
    field: str,
    value: float,
) -> None:
    invalid = replace(contract(), **{field: value})
    with pytest.raises(ValueError):
        assess_free_sd_stationarity(
            completed_run,
            invalid,
        )


def test_missing_constraint_evidence_is_fail_closed(completed_run) -> None:
    missing = replace(completed_run, constraint_residual=None)
    assessment = assess_free_sd_stationarity(missing, contract())
    assert assessment.constraint_pass is False
    assert assessment.stationarity_pass is False
    assert assessment.inference_ready is False


def test_tampered_redundant_run_fields_are_rejected(completed_run) -> None:
    tampered = replace(completed_run, algorithm_terminated=False)
    with pytest.raises(ValueError, match="do not reconstruct"):
        assess_free_sd_stationarity(tampered, contract())

    tampered_score = replace(completed_run, final_projected_gradient_supnorm=0.13)
    with pytest.raises(ValueError, match="do not reconstruct"):
        assess_free_sd_stationarity(tampered_score, contract())

    bad_restart = replace(completed_run.restart_polish, sigma=2.0)
    coordinated_tamper = replace(
        completed_run,
        restart_polish=bad_restart,
        final_sigma=2.0,
    )
    with pytest.raises(ValueError, match="restart polish derived fields"):
        assess_free_sd_stationarity(coordinated_tamper, contract())

    bad_structural = replace(
        completed_run.restart_polish,
        structural_parameters=(99.0, *completed_run.restart_polish.structural_parameters[1:]),
    )
    with pytest.raises(ValueError, match="restart polish derived fields"):
        assess_free_sd_stationarity(
            replace(completed_run, restart_polish=bad_structural),
            contract(),
        )


def test_structural_bounds_are_out_of_scope_for_phase_one() -> None:
    with pytest.raises(ValueError, match="bounded structural"):
        run_free_sd_stationarity_v2(
            structural_start=np.array([-2.0, 2.5]),
            sigma_start=0.4,
            value_function=value,
            structural_value_gradient=value_gradient,
            observations=100,
            structural_bounds=[(-3.0, 3.0), (None, None)],
            constraint_residual_function=lambda _parameters: 0.0,
        )


def test_scalar_and_value_gradient_objectives_must_agree() -> None:
    inconsistent = run_free_sd_stationarity_v2(
        structural_start=np.array([-2.0, 2.5]),
        sigma_start=0.4,
        value_function=lambda parameters, sigma: value(parameters, sigma) + 1.0,
        structural_value_gradient=value_gradient,
        observations=100,
        constraint_residual_function=lambda _parameters: 0.0,
    )
    assessment = assess_free_sd_stationarity(inconsistent, contract())
    assert assessment.objective_value_consistency_pass is False
    assert assessment.stationarity_pass is False
    assert assessment.inference_ready is False


def test_optimizer_failure_is_not_hidden_by_finite_coordinates(monkeypatch) -> None:
    def failed_minimize(fun, x0, **_kwargs):
        fun(np.asarray(x0, dtype=float))
        return SimpleNamespace(
            x=np.asarray(x0, dtype=float),
            success=False,
            status=2,
            message="line search failed",
            nit=1,
            nfev=1,
            njev=1,
        )

    monkeypatch.setattr(mml_stationarity, "minimize", failed_minimize)
    run = run_free_sd_stationarity_v2(
        structural_start=np.array([-2.0, 2.5]),
        sigma_start=0.4,
        value_function=value,
        structural_value_gradient=value_gradient,
        observations=100,
        constraint_residual_function=lambda _parameters: 0.0,
    )
    assessment = assess_free_sd_stationarity(run, contract())
    assert assessment.algorithm_terminated is False
    assert assessment.stationarity_pass is False
    assert assessment.inference_ready is False


def test_nonfinite_objective_fails_before_a_run_is_created() -> None:
    with pytest.raises(FloatingPointError):
        run_free_sd_stationarity_v2(
            structural_start=np.array([-2.0, 2.5]),
            sigma_start=0.4,
            value_function=lambda _parameters, _sigma: float("nan"),
            structural_value_gradient=lambda parameters, sigma: (
                value(parameters, sigma),
                parameters - TARGET,
            ),
            observations=100,
            constraint_residual_function=lambda _parameters: 0.0,
        )


def shifted_problem(anchor):
    """A raw objective whose small local improvements round to zero."""
    target = np.r_[TARGET, LOG_SIGMA_TARGET]
    offset = 1e8
    raw = lambda parameters, sigma: offset + value(parameters, sigma)
    joint = lambda point: (raw(point[:-1], np.exp(point[-1])), point - target)
    structural = lambda parameters, sigma: (raw(parameters, sigma), parameters - TARGET)
    def difference(point):
        change = point - anchor
        return float(change @ (anchor - target) + 0.5 * (change @ change))
    return raw, structural, joint, difference


def test_shifted_optimization_preserves_raw_nll_and_quadrature_reconstruction() -> None:
    from mfrm_app.mml_quadrature_sensitivity import assess_quadrature_sensitivity
    from tests.test_mml_quadrature_sensitivity import sensitivity_contract

    anchor = np.r_[TARGET, LOG_SIGMA_TARGET] + 1e-5
    raw, structural, joint, difference = shifted_problem(anchor)
    run = run_free_sd_stationarity_v2(
        anchor[:-1], np.exp(anchor[-1]), raw, structural, observations=100,
        joint_value_gradient=joint, objective_difference=difference, objective_anchor=anchor,
        options=mml_stationarity.JointPolishOptions(ftol=0., gtol=1e-10),
        constraint_residual_function=lambda p: 0.,
    )
    assert assess_free_sd_stationarity(run, contract()).stationarity_pass
    primary = run.primary_polish
    assert run.final_objective == primary.initial_objective == 1e8
    assert primary.objective_improvement == 0.
    assert primary.optimization_improvement > 1e-10
    assert primary.objective_shift.final_difference < 0
    assert primary.objective_shift.anchor_coordinates == run.restart_polish.objective_shift.anchor_coordinates
    assert primary.last_query_raw_objective == 1e8
    assert abs(primary.last_query_optimizer_objective) < 1e-8
    assert abs(primary.optimizer_reported_objective) < 1e-8
    args = dict(problem_digest='0'*64, primary_quadrature_points=31, sensitivity_quadrature_points=61,
                primary_run=run, sensitivity_run=run, stationarity_contract=contract(),
                sensitivity_contract=sensitivity_contract())
    sensitivity = assess_quadrature_sensitivity(
        primary_value_function=raw, sensitivity_value_function=raw, **args)
    assert sensitivity.objective_reconstruction_pass and sensitivity.numerical_sensitivity_pass
    assert not sensitivity.scientific_inference_ready
    wrong = lambda parameters, sigma: difference(np.r_[parameters, np.log(sigma)])
    rejected = assess_quadrature_sensitivity(
        primary_value_function=wrong, sensitivity_value_function=wrong, **args)
    assert not rejected.objective_reconstruction_pass and not rejected.numerical_sensitivity_pass
    broken = replace(primary.objective_shift, anchor_objective=1e8+1)
    with pytest.raises(ValueError, match='objective origin'):
        assess_free_sd_stationarity(replace(run, primary_polish=replace(primary, objective_shift=broken)), contract())


def test_two_stage_retains_failed_preliminary_and_refuses_unusable_anchor(monkeypatch) -> None:
    anchor = np.r_[TARGET, LOG_SIGMA_TARGET] + 1e-5
    raw, structural, joint, _ = shifted_problem(anchor)
    original = mml_stationarity.minimize
    calls = []
    def preliminary_failure(fun, x0, **kwargs):
        calls.append(np.asarray(x0).copy())
        if len(calls) <= 2:
            f, _ = fun(x0)
            return SimpleNamespace(x=np.asarray(x0).copy(), fun=f, success=False,
                                   status=2, message='retained preliminary failure', nit=0, nfev=1, njev=1)
        return original(fun, x0, **kwargs)
    monkeypatch.setattr(mml_stationarity, 'minimize', preliminary_failure)
    factories = []
    def factory(a):
        factories.append(a.copy())
        return shifted_problem(a)[-1]
    args = dict(observations=100, joint_value_gradient=joint, difference_factory=factory,
                refinement_options=mml_stationarity.JointPolishOptions(ftol=0., gtol=1e-10),
                constraint_residual_function=lambda p: 0.)
    result = run_free_sd_two_stage(anchor[:-1], np.exp(anchor[-1]), raw, structural,
                                   anchor_gradient_limit=1e-3, **args)
    assert result.anchor_admitted and len(calls) == 4 and len(factories) == 1
    assert not result.preliminary.algorithm_terminated
    assert result.preliminary.primary_polish.optimizer_status == 2
    assert result.refinement.algorithm_terminated
    assert assess_free_sd_stationarity(result.refinement, contract()).stationarity_pass
    assert np.array_equal(factories[0], result.preliminary.restart_polish.joint_coordinates)
    calls.clear(); factories.clear()
    refused = run_free_sd_two_stage(anchor[:-1], np.exp(anchor[-1]), raw, structural,
                                    anchor_gradient_limit=1e-7, **args)
    assert not refused.anchor_admitted and refused.refinement is None
    assert len(calls) == 2 and not factories


def test_wrong_difference_cannot_pass_via_raw_nll_or_zero_terminal_gradient() -> None:
    anchor = np.r_[TARGET, LOG_SIGMA_TARGET] + 0.2
    raw, structural, joint, difference = shifted_problem(anchor)
    args = dict(observations=100, joint_value_gradient=joint, objective_anchor=anchor,
                options=mml_stationarity.JointPolishOptions(ftol=0.),
                constraint_residual_function=lambda p: 0.)
    bad = run_free_sd_stationarity_v2(anchor[:-1], np.exp(anchor[-1]), raw, structural,
                                      objective_difference=lambda p: 2*difference(p), **args)
    assessment = assess_free_sd_stationarity(bad, contract())
    assert not assessment.objective_value_consistency_pass and not assessment.stationarity_pass
    with pytest.raises(ValueError, match='zero at its anchor'):
        run_free_sd_stationarity_v2(anchor[:-1], np.exp(anchor[-1]), raw, structural,
                                    objective_difference=lambda p: difference(p)+1, **args)
    with pytest.raises(ValueError, match='supplied together'):
        run_free_sd_stationarity_v2(anchor[:-1], np.exp(anchor[-1]), raw, structural,
                                    objective_difference=None, **args)
