from __future__ import annotations

import copy
from dataclasses import replace

import numpy as np
import pytest

import streamlit_app as app
from mfrm_app.mml_engine_v2 import StationarityContract
from mfrm_app.mml_quadrature_sensitivity import QuadratureSensitivityContract
from mfrm_app.mml_stationarity import JointPolishOptions
from validation.mml_free_sd_quadrature_adapter import (
    app_free_sd_evaluator_implementation_sha256,
    app_free_sd_problem_digest,
    resolve_app_free_sd_objective_evaluator,
    run_app_free_sd_quadrature_sensitivity,
)
from validation.mml_free_sd_stationarity_adapter import prepare_app_free_sd_problem


@pytest.fixture(scope="module")
def native_result():
    n_person = 24
    data = app._generate_mfrm_rsm_from_params(
        {
            "persons": [f"P{i:03d}" for i in range(n_person)],
            "raters": ["R1", "R2"],
            "tasks": ["T1", "T2"],
            "criteria": ["C1", "C2"],
            "theta_sd": 1.2,
            "rater_severities": np.array([-0.2, 0.2]),
            "task_difficulties": np.array([-0.15, 0.15]),
            "criterion_difficulties": np.array([-0.1, 0.1]),
            "tau": np.array([-1.0, 0.0, 1.0]),
        },
        seed=20260811,
    )
    return app.mfrm_estimate(
        data=data,
        person_col="Person",
        facet_cols=["Rater", "Task", "Criterion"],
        score_col="Score",
        rating_min=0,
        rating_max=3,
        model="RSM",
        method="MML",
        noncenter_facet="Criterion",
        mml_engine="EM",
        quad_points=9,
        estimate_population_sd=True,
        population_prior_sd=1.0,
        maxit=300,
        reltol=1e-6,
        min_obs_per_element=1,
        min_obs_per_category=1,
    )


def development_stationarity_contract() -> StationarityContract:
    return StationarityContract(
        max_projected_gradient_supnorm=1e-3,
        max_standardized_score_supnorm=1e-3,
        max_newton_correction_supnorm=1e-3,
        max_restart_improvement_per_observation=1e-8,
        max_restart_displacement=1e-3,
        max_gradient_fd_disagreement=1e-4,
        max_objective_value_disagreement=1e-8,
        max_information_relative_symmetry_residual=1e-3,
        max_information_condition_number=1e8,
        max_constraint_residual=1e-8,
        objective_worsening_tolerance_total=1e-8,
        objective_worsening_tolerance_per_observation=1e-10,
        sigma_boundary_tolerance=1e-6,
    )


def development_sensitivity_contract() -> QuadratureSensitivityContract:
    return QuadratureSensitivityContract(
        primary_quadrature_points=9,
        sensitivity_quadrature_points=15,
        max_structural_parameter_difference=0.1,
        max_log_sigma_difference=0.1,
        max_sensitivity_optimization_gain_per_observation=1e-4,
        max_negative_sensitivity_gain_per_observation=1e-8,
        max_abs_primary_back_evaluation_change_per_observation=1e-4,
        max_objective_reconstruction_disagreement=1e-8,
    )


def test_problem_digest_is_q_independent_but_data_dependent(native_result) -> None:
    q9 = prepare_app_free_sd_problem(native_result, quadrature_points=9)
    q15 = prepare_app_free_sd_problem(native_result, quadrature_points=15)
    digest = app_free_sd_problem_digest(native_result, q9)
    assert digest == app_free_sd_problem_digest(native_result, q15)
    assert len(digest) == 64

    changed = dict(native_result)
    changed["prep"] = dict(native_result["prep"])
    changed["prep"]["data"] = native_result["prep"]["data"].copy()
    changed["prep"]["data"].loc[0, "Score"] += 1
    assert app_free_sd_problem_digest(changed, q9) != digest


def test_problem_digest_binds_effective_signs_indices_bounds_and_shape(native_result) -> None:
    problem = prepare_app_free_sd_problem(native_result, quadrature_points=9)
    baseline = app_free_sd_problem_digest(native_result, problem)

    sign_config = copy.deepcopy(problem.config)
    sign_config["facet_signs"]["Rater"] = 1
    sign_config["positive_facets"] = ["Rater"]
    assert (
        app_free_sd_problem_digest(
            native_result,
            replace(problem, config=sign_config),
        )
        != baseline
    )

    inconsistent_config = copy.deepcopy(problem.config)
    inconsistent_config["facet_signs"]["Rater"] = 1
    with pytest.raises(ValueError, match="positive_facets.*inconsistent"):
        app_free_sd_problem_digest(
            native_result,
            replace(problem, config=inconsistent_config),
        )

    reshaped_idx = dict(problem.idx)
    reshaped_idx["score_k"] = np.asarray(problem.idx["score_k"]).reshape(-1, 1)
    assert (
        app_free_sd_problem_digest(
            native_result,
            replace(problem, idx=reshaped_idx),
        )
        != baseline
    )

    changed_bounds = replace(problem, sigma_bounds=(0.1, problem.sigma_bounds[1]))
    assert app_free_sd_problem_digest(native_result, changed_bounds) != baseline


def test_problem_digest_rejects_nonfinite_unmodelled_inputs(native_result) -> None:
    problem = prepare_app_free_sd_problem(native_result, quadrature_points=9)
    bad_config = copy.deepcopy(problem.config)
    population = dict(bad_config.get("population_model", {}))
    population["X"] = np.array([[np.nan]])
    bad_config["population_model"] = population
    with pytest.raises(ValueError, match="population design must be a finite"):
        app_free_sd_problem_digest(
            native_result,
            replace(problem, config=bad_config),
        )


def test_real_app_qfit_qcheck_remains_numerical_only(native_result) -> None:
    bundle = run_app_free_sd_quadrature_sensitivity(
        native_result,
        primary_quadrature_points=9,
        sensitivity_quadrature_points=15,
        stationarity_contract=development_stationarity_contract(),
        sensitivity_contract=development_sensitivity_contract(),
        options=JointPolishOptions(maxiter=250, gtol=1e-8, ftol=1e-15),
    )
    assessment = bundle.assessment
    assert assessment.problem_digest == bundle.problem_digest
    assert assessment.primary_quadrature_points == 9
    assert assessment.sensitivity_quadrature_points == 15
    assert assessment.finite_cross_evaluation is True
    assert assessment.scientific_inference_ready is False
    assert bundle.objective_evaluator.problem_digest == bundle.problem_digest
    bundle.objective_evaluator.validate()
    assert (
        bundle.objective_evaluator_implementation_sha256
        == app_free_sd_evaluator_implementation_sha256()
    )
    replay = resolve_app_free_sd_objective_evaluator(
        native_result,
        primary_quadrature_points=9,
        sensitivity_quadrature_points=15,
        expected_problem_digest=bundle.problem_digest,
        expected_implementation_sha256=(
            bundle.objective_evaluator_implementation_sha256
        ),
    )
    assert replay.problem_digest == bundle.problem_digest
    assert replay.primary_quadrature_points == 9
    assert replay.sensitivity_quadrature_points == 15


def test_invalid_q_contract_fails_before_numerical_run(native_result, monkeypatch) -> None:
    import validation.mml_free_sd_quadrature_adapter as adapter

    called = False

    def unexpected_run(*_args, **_kwargs):
        nonlocal called
        called = True
        raise AssertionError("optimizer should not run")

    monkeypatch.setattr(adapter, "run_free_sd_stationarity_v2", unexpected_run)
    with pytest.raises(ValueError, match="differ from the sensitivity contract"):
        adapter.run_app_free_sd_quadrature_sensitivity(
            native_result,
            primary_quadrature_points=9,
            sensitivity_quadrature_points=13,
            stationarity_contract=development_stationarity_contract(),
            sensitivity_contract=development_sensitivity_contract(),
        )
    assert called is False


def test_evaluator_identity_rejects_source_drift(monkeypatch) -> None:
    import validation.mml_free_sd_quadrature_adapter as adapter

    monkeypatch.setattr(adapter, "_source_sha256", lambda _path: "0" * 64)
    with pytest.raises(ValueError, match="source changed after module import"):
        adapter.app_free_sd_evaluator_implementation_sha256()
