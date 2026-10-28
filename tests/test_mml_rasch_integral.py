"""Checks of the independent finite/continuous RSM/PCM information reference."""

import numpy as np
import pandas as pd
import pytest
from scipy.special import roots_hermitenorm

import streamlit_app as app
from validation.mml_free_sd_stationarity_adapter import prepare_app_free_sd_problem_from_data
from validation.mml_rasch_integral import RaschIntegral
from validation.mml_rasch_quadrature_review import expansion


def problem(model, sign):
    data = pd.DataFrame([(f"P{p}", f"R{r}", f"C{c}", (p+2*r+c) % 4)
        for p in range(4) for r in range(2) for c in range(2)
        if (p, r, c) != (1, 1, 1)], columns=["Person", "Rater", "Criterion", "Score"])
    _, native = prepare_app_free_sd_problem_from_data(data, person_col="Person", facet_cols=["Rater", "Criterion"],
        score_col="Score", rating_min=0, rating_max=3, model=model, noncenter_facet="Criterion",
        step_facet="Criterion" if model == "PCM" else None, positive_facets=("Rater",) if sign == 1 else ())
    return native, RaschIntegral(data, native.config)


@pytest.mark.parametrize("model", ["RSM", "PCM"])
@pytest.mark.parametrize("sign", [-1, 1])
@pytest.mark.parametrize("sigma", [.7, 1.4])
def test_independent_finite_information_reproduces_app(model, sign, sigma):
    native, reference = problem(model, sign)
    point = np.r_[np.linspace(-.4, .5, reference.d), np.log(sigma)]
    z, w = roots_hermitenorm(31); w /= np.sqrt(2*np.pi)
    independent = reference.finite(point, z, w)
    nll, gradient = native.joint_value_gradient(point)
    assert abs(independent["nll"]-nll) < 1e-12
    np.testing.assert_allclose(independent["gradient"], gradient, atol=1e-11)
    hessian = app._compute_mml_observed_information(point, native.idx, native.config, native.sizes,
        app.gauss_hermite_normal(31, sd=sigma), include_log_sigma=True)
    np.testing.assert_allclose(independent["information"], hessian, atol=2e-8)


def test_continuous_information_and_density_score_have_independent_derivatives():
    _, reference = problem("PCM", 1)
    point = np.r_[np.linspace(-.3, .4, reference.d), np.log(1.1)]
    result = reference.continuous(point, information=True)
    tighter = reference.continuous(point, bound=14, rel_tol=1e-11, abs_tol=1e-13, information=True)
    for key in ("nll", "gradient", "information", "eap", "sd"):
        np.testing.assert_allclose(result[key], tighter[key], atol=1e-8)
    assert abs(result["gradient"][-1]-result["density_log_sigma_score"]) < 1e-9
    finite_hessian = []
    for delta in np.eye(len(point))*1e-4:
        plus, minus = reference.continuous(point+delta), reference.continuous(point-delta)
        finite_hessian.append((plus["gradient"]-minus["gradient"])/2e-4)
    np.testing.assert_allclose(result["information"], np.array(finite_hessian).T, atol=2e-7)


def test_binary_symmetry_control_has_no_population_sd_information():
    data = pd.DataFrame({"Person": ["P0", "P1"], "Item": ["I0", "I0"], "Score": [0, 1]})
    cfg = dict(model="RSM", facet_names=["Item"], facet_levels={"Item": ["I0"]},
               facet_signs={"Item": -1}, noncenter_facet="Item", n_cat=2, rating_min=0)
    reference = RaschIntegral(data, cfg)
    result = reference.continuous(np.array([0., np.log(1.7)]), information=True)
    assert result["nll"] == pytest.approx(2*np.log(2), abs=1e-12)
    assert np.max(abs(result["gradient"])) < 1e-12
    assert abs(result["information"][-1, -1]) < 1e-12
    assert result["eap"][0] == pytest.approx(-result["eap"][1], abs=1e-12)


def test_expanded_se_review_includes_dependent_steps_and_sd_chain_rule():
    config = dict(model="RSM", facet_names=["Rater", "Criterion"],
        facet_levels={"Rater": ["R0", "R1"], "Criterion": ["C0", "C1"]},
        noncenter_facet="Criterion", n_cat=4)
    point = np.array([.2, -.3, .7, -.8, .1, np.log(1.4)])
    labels, values, jacobian = expansion(config, point)
    np.testing.assert_allclose(values, [.2, -.2, -.3, .7, -.8, .1, .7, 1.4])
    covariance = np.diag([.01, .04, .09, .16, .25, .36])
    covariance[3, 4] = covariance[4, 3] = .02
    variance = np.diag(jacobian @ covariance @ jacobian.T)
    assert labels[-2:] == ["step:shared:3", "population_sd"]
    assert variance[-2] == pytest.approx(.16+.25+2*.02)
    assert variance[-1] == pytest.approx(1.4**2*.36)
