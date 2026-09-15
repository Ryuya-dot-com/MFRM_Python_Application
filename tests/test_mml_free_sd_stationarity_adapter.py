from __future__ import annotations

import copy
import json
from dataclasses import replace
from decimal import Decimal, localcontext

import numpy as np
import pytest

import streamlit_app as app
from mfrm_app.mml_stationarity import JointPolishOptions, audit_joint_gradient, make_joint_free_sd_functions
from validation.mml_free_sd_stationarity_adapter import (
    prepare_app_free_sd_problem,
    run_app_free_sd_stationarity_v2,
    run_app_free_sd_two_stage,
    _log_average_exp,
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


def _observed_problem(native, model):
    base = prepare_app_free_sd_problem(native, quadrature_points=31)
    config = copy.deepcopy(base.config)
    config.update(model=model, step_facet="Criterion", facet_signs={"Rater": 1, "Task": -1, "Criterion": -1})
    config["facet_specs"]["Rater"] = app.build_facet_constraint(["R1", "R2"], anchors={"R1": 0.35})
    config["facet_specs"]["Task"] = app.build_facet_constraint(
        ["T1", "T2"], groups={"T1": "g", "T2": "g"}, group_values={"g": 0.2},
    )
    config["population_model"] = {
        "enabled": True, "n_params": 1,
        "X": np.linspace(-1.0, 2.0, config["n_person"])[:, None],
    }
    prep = dict(native["prep"])
    data = prep["data"]
    # One wholly missing person, scattered missing rows, one zero-weight person.
    keep = (data["Person"].cat.codes != 0) & (np.arange(len(data)) % 5 != 0)
    prep["data"] = data.loc[keep].copy()
    prep["data"]["Weight"] = np.resize([0.0, 0.5, 1.5, 2.0], keep.sum())
    prep["data"].loc[prep["data"]["Person"].cat.codes == 1, "Weight"] = 0.0
    idx = app.build_indices(prep, step_facet="Criterion")
    sizes = app.build_param_sizes(config)
    return replace(base, config=config, sizes=sizes, idx=idx, observations=float(idx["weight"].sum()))


def _decimal_observed_nll(problem, point, precision):
    """Independent scalar formula for this explicitly specified two-level design.

    No app expansion, probability, posterior or difference helper is used.
    Binary GH nodes/weights define the target; Decimal evaluates that same rule.
    """
    with localcontext() as context:
        context.prec = precision
        d = lambda x: Decimal.from_float(float(x))
        p = list(map(d, point))
        sigma = p[-1].exp()
        rater, task, criterion = [d(0.35), p[1]], [p[2], 2*d(0.2)-p[2]], p[3:5]
        step_free = [p[5:7]] if problem.config["model"] == "RSM" else [p[5:7], p[7:9]]
        steps = [[Decimal(0), s[0], s[0]+s[1], Decimal(0)] for s in step_free]
        quad = app.gauss_hermite_normal(problem.quadrature_points)
        nll = Decimal(0)
        idx = problem.idx
        for person in range(problem.config["n_person"]):
            rows = np.flatnonzero((idx["person"] == person) & (idx["weight"] > 0))
            if not len(rows):
                continue
            mu = d(problem.config["population_model"]["X"][person, 0]) * p[0]
            mass = Decimal(0)
            for z, w in zip(quad["nodes"], quad["weights"]):
                ll = Decimal(0)
                for row in rows:
                    r, t, c = (idx["facets"][name][row] for name in ("Rater", "Task", "Criterion"))
                    eta = sigma*d(z) + mu + rater[r] - task[t] - criterion[c]
                    cumulative = steps[0 if problem.config["model"] == "RSM" else c]
                    logits = [k*eta-cumulative[k] for k in range(4)]
                    normalizer = sum(v.exp() for v in logits).ln()
                    ll += d(idx["weight"][row]) * (logits[idx["score_k"][row]] - normalizer)
                mass += d(w)*ll.exp()
            nll -= mass.ln()
        return nll


@pytest.mark.parametrize("model", ["RSM", "PCM"])
def test_observed_row_difference_matches_independent_decimal_and_gradient(native_free_sd_result, model, record_property):
    problem = _observed_problem(native_free_sd_result, model)
    anchor = np.r_[np.linspace(0.8, -0.3, sum(problem.sizes.values())), np.log(1.7)]
    difference = problem.likelihood_difference(anchor)
    assert difference(anchor) == 0.0
    assert problem.constraint_residual(anchor[:-1]) < 1e-12
    direction = np.resize([-1.0, 0.0, 1.0], len(anchor))
    comparisons = []
    for magnitude in (1e-11, 0.25):
        point = anchor + magnitude*direction
        values = []
        for precision in (60, 90):
            with localcontext() as context:
                context.prec = precision
                values.append(_decimal_observed_nll(problem, point, precision)
                              - _decimal_observed_nll(problem, anchor, precision))
        assert abs(values[0] - values[1]) < Decimal("1e-45")
        expected = float(values[1])
        comparisons.append(dict(magnitude=magnitude, point=point.tolist(), difference=difference(point),
                                decimal60=str(values[0]), decimal90=str(values[1]),
                                absolute_error=abs(difference(point)-expected)))
        assert abs(difference(point) - expected) <= 1e-18 + 1e-11*abs(expected)
        assert abs(problem.value(anchor[:-1], np.exp(anchor[-1])) + difference(point)
                   - problem.value(point[:-1], np.exp(point[-1]))) < 1e-9
    gradients = []
    for step in (1e-5, 1e-6):
        audit = audit_joint_gradient(difference, problem.joint_value_gradient, anchor, relative_step=step)
        gradients.append(dict(step=step, maximum_absolute_difference=audit.maximum_absolute_difference))
        assert audit.maximum_absolute_difference < 1e-7
    record_property("numerical_comparisons", json.dumps(dict(model=model, anchor=anchor.tolist(),
                    rows=len(problem.idx["score_k"]), weight_sum=problem.observations,
                    comparisons=comparisons, gradients=gradients)))
    # A pure SD move must include moving nodes, even with nonzero person means.
    moved = anchor.copy()
    moved[-1] += 0.01
    assert abs(difference(moved)) > 1e-5
    # Zero-weight observations carry no likelihood or gradient information.
    reduced = copy.deepcopy(problem.idx)
    keep = reduced["weight"] > 0
    for key in ("person", "score_k", "weight", "step_idx"):
        reduced[key] = reduced[key][keep]
    reduced["facets"] = {key: value[keep] for key, value in reduced["facets"].items()}
    reduced["rows_by_person"] = [np.flatnonzero(reduced["person"] == p) for p in range(problem.config["n_person"])]
    other = replace(problem, idx=reduced)
    assert abs(other.likelihood_difference(anchor)(moved) - difference(moved)) < 1e-12
    assert np.allclose(other.joint_value_gradient(moved)[1], problem.joint_value_gradient(moved)[1], rtol=0, atol=1e-12)


def test_observed_difference_rejects_invalid_inputs_and_retains_log_tails(native_free_sd_result):
    problem = _observed_problem(native_free_sd_result, "PCM")
    anchor = np.zeros(sum(problem.sizes.values()) + 1)
    difference = problem.likelihood_difference(anchor)
    for point in (anchor[:-1], anchor + np.nan, np.r_[anchor[:-1], np.log(20)]):
        with pytest.raises(ValueError):
            difference(point)
    for invalid in (-0.5, np.nan, np.inf):
        idx = copy.deepcopy(problem.idx)
        idx["weight"][0] = invalid
        with pytest.raises(ValueError, match="weights"):
            replace(problem, idx=idx).likelihood_difference(anchor)
    idx = copy.deepcopy(problem.idx)
    idx["rows_by_person"][0] = np.array([0])
    with pytest.raises(ValueError, match="row groups"):
        replace(problem, idx=idx).likelihood_difference(anchor)
    assert abs(_log_average_exp(np.array([0., -1000.]), np.array([0., 1000.])) - np.log(2)) < 1e-14
    for change in (0., 1e-20, .5, np.nextafter(.5, 1.), -.5, np.nextafter(-.5, -1.)):
        assert abs(_log_average_exp(np.log([.2, .8]), np.full(2, change)) - change) < 1e-15


def test_native_two_stage_adapter_retains_raw_objective(native_free_sd_result, record_property):
    problem, run = run_app_free_sd_two_stage(
        native_free_sd_result, anchor_gradient_limit=1e-3,
        preliminary_options=JointPolishOptions(maxiter=250, gtol=1e-8, ftol=1e-15),
        refinement_options=JointPolishOptions(maxiter=250, gtol=1e-8, ftol=0.0),
    )
    assert run.anchor_admitted and run.refinement is not None
    record_property("two_stage_run", json.dumps(run.to_dict()))
    assert run.refinement.algorithm_terminated
    assert run.refinement.final_projected_gradient_supnorm <= 1e-8
    for polish in (run.refinement.primary_polish, run.refinement.restart_polish):
        assert polish.objective_shift.reconstruction_max_abs_difference < 1e-9
        assert abs(polish.final_objective - problem.value(np.array(polish.structural_parameters), polish.sigma)) < 1e-9
