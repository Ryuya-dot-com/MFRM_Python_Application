from __future__ import annotations

import copy
from dataclasses import replace

import numpy as np
import pytest

import streamlit_app as app
from mfrm_app.mml_stationarity import JointPolishOptions, audit_joint_gradient, make_joint_free_sd_functions
from validation.mml_free_sd_stationarity_adapter import (
    prepare_app_free_sd_problem,
    run_app_free_sd_stationarity_v2,
)


@pytest.fixture(scope="module")
def native_free_sd_result():
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


def test_adapter_uses_same_native_objective_without_mutating_result(
    native_free_sd_result,
) -> None:
    before = float(native_free_sd_result["config"]["estimated_population_sd"])
    problem = prepare_app_free_sd_problem(native_free_sd_result)
    value, gradient = problem.value_gradient(
        np.asarray(problem.structural_start),
        problem.sigma_start,
    )
    assert np.isfinite(value)
    assert np.isfinite(gradient).all()
    assert gradient.size == len(problem.structural_start)
    assert problem.value(np.asarray(problem.structural_start), problem.sigma_start) == pytest.approx(
        value, abs=1e-10
    )
    assert problem.constraint_residual(np.asarray(problem.structural_start)) <= 1e-10
    assert native_free_sd_result["config"]["estimated_population_sd"] == before
    joint_value, joint_gradient = problem.joint_value_gradient(
        np.r_[problem.structural_start, np.log(problem.sigma_start)]
    )
    assert joint_value == value
    assert np.array_equal(joint_gradient[:-1], gradient)


def test_native_free_sd_joint_polish_reaches_finite_stationary_solution(
    native_free_sd_result,
) -> None:
    problem, run = run_app_free_sd_stationarity_v2(
        native_free_sd_result,
        options=JointPolishOptions(maxiter=250, gtol=1e-8, ftol=1e-15),
    )
    assert run.algorithm_terminated is True
    assert run.primary_polish.gradient_method == "provided_joint_gradient"
    assert run.restart_polish.gradient_method == "provided_joint_gradient"
    assert run.gradient_audit_h.maximum_absolute_difference < 1e-5
    assert run.gradient_audit_h_over_2.maximum_absolute_difference < 1e-5
    assert run.finite_solution is True
    assert run.objective_nonworsening is True
    assert problem.sigma_bounds[0] < run.final_sigma < problem.sigma_bounds[1]
    assert run.final_projected_gradient_supnorm < 1e-3
    assert run.gradient_step_agreement_supnorm < 1e-5
    assert run.constraint_residual is not None
    assert run.constraint_residual <= 1e-10


def test_adapter_rejects_fixed_sd_and_gpcm(native_free_sd_result) -> None:
    fixed = dict(native_free_sd_result)
    fixed["config"] = dict(native_free_sd_result["config"], estimate_population_sd=False)
    with pytest.raises(ValueError, match="free population SD"):
        prepare_app_free_sd_problem(fixed)

    gpcm = dict(native_free_sd_result)
    gpcm["config"] = dict(native_free_sd_result["config"], model="GPCM")
    with pytest.raises(ValueError, match="RSM and PCM only"):
        prepare_app_free_sd_problem(gpcm)
    with pytest.raises(ValueError, match="RSM and PCM only"):
        app.mfrm_loglik_mml_value_grad(
            np.array([]), {}, gpcm["config"], {}, {}, include_log_sigma=True,
        )


@pytest.mark.parametrize("model", ["RSM", "PCM"])
def test_joint_gradient_handles_weights_and_nonzero_person_means(native_free_sd_result, model):
    base = prepare_app_free_sd_problem(native_free_sd_result, quadrature_points=31)
    config = copy.deepcopy(base.config)
    config.update(model=model, step_facet="Criterion")
    config["population_model"] = {
        "enabled": True, "n_params": 1,
        "X": np.linspace(-1.0, 2.0, config["n_person"])[:, None],
    }
    sizes = app.build_param_sizes(config)
    idx = app.build_indices(native_free_sd_result["prep"], step_facet="Criterion")
    idx["weight"] = np.resize([0.0, 0.5, 1.5, 2.0], len(idx["score_k"]))
    problem = replace(base, config=config, sizes=sizes, idx=idx)
    point = np.r_[np.linspace(0.8, -0.3, sum(sizes.values())), np.log(1.7)]
    objective, analytic = make_joint_free_sd_functions(
        problem.value, problem.value_gradient, joint_value_gradient=problem.joint_value_gradient,
    )
    assert audit_joint_gradient(objective, analytic, point).maximum_absolute_difference < 1e-6
    assert audit_joint_gradient(objective, analytic, point, relative_step=5e-6).maximum_absolute_difference < 1e-6


def test_adapter_rejects_sigma_disagreement_and_fractional_quadrature(
    native_free_sd_result,
) -> None:
    mismatched = dict(native_free_sd_result)
    mismatched["config"] = dict(native_free_sd_result["config"])
    mismatched["config"]["estimated_population_sd"] += 0.01
    with pytest.raises(ValueError, match="population SD evidence differs"):
        prepare_app_free_sd_problem(mismatched)

    with pytest.raises(ValueError, match="integer of at least 3"):
        prepare_app_free_sd_problem(native_free_sd_result, quadrature_points=9.7)
