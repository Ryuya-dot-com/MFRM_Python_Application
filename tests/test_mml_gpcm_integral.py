"""Independent continuous GPCM information: density score and nonlinear slopes."""

import numpy as np
import pandas as pd
import pytest
from scipy.integrate import quad_vec

from validation.mml_gpcm_integral import GpcmIntegral
from validation.mml_rasch_integral import RaschIntegral


def problem(sign=-1, slope_facet="Criterion", step_facet="Criterion"):
    data = pd.DataFrame([(f"P{p}", f"R{r}", f"C{c}", (p+2*r+c) % 4)
        for p in range(3) for r in range(2) for c in range(3)
        if (p, r, c) != (1, 1, 1)], columns=["Person", "Rater", "Criterion", "Score"])
    cfg = dict(model="GPCM", facet_names=["Rater", "Criterion"],
        facet_levels={"Rater": ["R0", "R1"], "Criterion": ["C0", "C1", "C2"]},
        facet_signs={"Rater": sign, "Criterion": -1}, noncenter_facet="Criterion",
        step_facet=step_facet, slope_facet=slope_facet, n_cat=4, rating_min=0)
    return data, cfg, GpcmIntegral(data, cfg)


@pytest.mark.parametrize("sign", [-1, 1])
@pytest.mark.parametrize("sigma", [.7, 1.4])
def test_unit_slopes_match_rasch_with_distinct_sd_derivation(sign, sigma):
    data, cfg, reference = problem(sign)
    linear = np.linspace(-.3, .4, reference.linear_d)
    point = np.r_[linear, 0., 0., np.log(sigma)]
    actual = reference.continuous(point, information=True)
    pcm = RaschIntegral(data, dict(cfg, model="PCM")).continuous(np.r_[linear, np.log(sigma)], information=True)
    columns = np.r_[np.arange(reference.linear_d), reference.d]
    for key in ["nll", "eap", "sd"]:
        np.testing.assert_allclose(actual[key], pcm[key], atol=1e-11, rtol=0)
    np.testing.assert_allclose(actual["gradient"][columns], pcm["gradient"], atol=1e-11, rtol=0)
    np.testing.assert_allclose(actual["information"][np.ix_(columns, columns)], pcm["information"], atol=1e-10, rtol=0)


def test_continuous_full_gradient_and_information_match_differences():
    _, _, reference = problem(1)
    point = np.r_[np.linspace(-.3, .4, reference.linear_d), np.log([1.8, .8, 1.1])]
    result = reference.continuous(point, information=True)
    tighter = reference.continuous(point, bound=16, rel_tol=1e-12, abs_tol=1e-14, information=True)
    for key in ["nll", "gradient", "information", "eap", "sd"]:
        np.testing.assert_allclose(result[key], tighter[key], atol=1e-10, rtol=0)
    gradients, hessian = [], []
    for delta in np.eye(len(point))*1e-4:
        plus, minus = reference.continuous(point+delta), reference.continuous(point-delta)
        gradients.append((plus["nll"]-minus["nll"])/2e-4)
        hessian.append((plus["gradient"]-minus["gradient"])/2e-4)
    np.testing.assert_allclose(result["gradient"], gradients, atol=2e-7, rtol=0)
    np.testing.assert_allclose(result["information"], np.array(hessian).T, atol=5e-7, rtol=0)
    assert abs(result["information"][-3, -2]) > 1e-3
    assert abs(result["information"][-3, -1]) > 1e-3
    assert tighter["tail_information_abs_bound"] < 1e-30
    assert tighter["numeric_information_error_allowance"] < 1e-7


def test_binary_zero_difficulty_has_no_slope_or_sd_information():
    data = pd.DataFrame(dict(Person=["P0", "P1", "P2", "P3"],
                             Criterion=["C0", "C0", "C1", "C1"], Score=[0, 1, 0, 1]))
    cfg = dict(model="GPCM", facet_names=["Criterion"], facet_levels={"Criterion": ["C0", "C1"]},
        facet_signs={"Criterion": -1}, noncenter_facet="Criterion", step_facet="Criterion",
        slope_facet="Criterion", n_cat=2, rating_min=0)
    result = GpcmIntegral(data, cfg).continuous([0., 0., np.log(2), np.log(1.3)], information=True)
    assert result["nll"] == pytest.approx(4*np.log(2), abs=1e-12)
    np.testing.assert_allclose(result["gradient"], 0., atol=2e-13)
    np.testing.assert_allclose(result["information"][2:, 2:], 0., atol=2e-12)
    np.testing.assert_allclose(result["eap"]**2+result["sd"]**2, 1.3**2, atol=1e-12)


def test_polynomial_tail_enclosures_include_nonlinear_hessian_and_density_score():
    _, _, reference = problem()
    point = np.r_[np.linspace(-.4, .5, reference.linear_d), np.log([2., .8, 1.1])]
    bound = 2.
    score_bound, second_bound, _ = reference._tail_numerators(point, bound)
    def integrand(z):
        blocks = []
        for sign in [-1, 1]:
            ll, score, _, h = reference.conditional(point, sign*z, information=True)
            mass = np.exp(ll-z*z/2)/np.sqrt(2*np.pi)
            second = h-np.einsum("pi,pj->pij", score, score)
            blocks.append(mass[:, None]*np.c_[abs(score), abs(second).reshape(reference.n, -1)])
        return sum(blocks)
    actual, _ = quad_vec(integrand, bound, 12., epsabs=1e-12, epsrel=1e-10)
    assert np.all(actual[:, :reference.d+1] <= score_bound)
    assert np.all(actual[:, reference.d+1:].reshape(second_bound.shape) <= second_bound)


def test_reference_rejects_unsupported_weighting_and_coordinates():
    data, cfg, reference = problem()
    with pytest.raises(ValueError, match="unit-weight"):
        GpcmIntegral(data.assign(Weight=2.), cfg)
    with pytest.raises(ValueError, match="coordinates"):
        reference.continuous(np.zeros(reference.d))
    with pytest.raises(ValueError, match="GPCM"):
        GpcmIntegral(data, dict(cfg, model="PCM"))


@pytest.mark.parametrize("slope_facet", ["Criterion", "Rater"])
def test_continuous_full_information_matches_app_at_accurate_finite_quadrature(slope_facet):
    import streamlit_app as app
    from scipy.special import roots_hermitenorm
    # The app's bounded GPCM contract requires a common step/slope owner.
    data, cfg, reference = problem(1, slope_facet, step_facet=slope_facet)
    cfg.update(method="MML", n_person=reference.n, estimate_population_sd=True,
        facet_specs={f: app.build_facet_constraint(levels, centered=f != cfg["noncenter_facet"])
                     for f, levels in cfg["facet_levels"].items()})
    sizes = app.build_param_sizes(cfg)
    idx = dict(person=reference.person, score_k=reference.y, weight=np.ones(len(data)),
        facets={f: np.array([cfg["facet_levels"][f].index(v) for v in data[f]]) for f in cfg["facet_names"]})
    idx["step_idx"] = idx["facets"][cfg["step_facet"]]; idx["slope_idx"] = idx["facets"][slope_facet]
    point = np.r_[np.linspace(-.3, .4, reference.linear_d),
                  np.linspace(np.log(1.8), np.log(.8), reference.d-reference.linear_d), np.log(1.1)]
    result = reference.continuous(point, information=True)
    coarse_quad = app.gauss_hermite_normal(321, sd=1.1)
    _, coarse_gradient = app.mfrm_loglik_mml_value_grad(point[:-1], idx, cfg, sizes, coarse_quad, include_log_sigma=True)
    # Q321 left a 4.71e-9 gradient error in this control. Preserve that
    # diagnostic; refine the independent rule instead of weakening tolerances.
    # SciPy supports high order without numpy.hermgauss overflow. Exactly-zero
    # floating-point tail weights contribute zero and are omitted before log(w).
    z, w = roots_hermitenorm(641); w /= np.sqrt(2*np.pi)
    positive = w > 0
    quad = dict(nodes=1.1*z[positive], weights=w[positive], sd=1.1)
    assert w.sum() == pytest.approx(1., abs=1e-14)
    nll, gradient = app.mfrm_loglik_mml_value_grad(point[:-1], idx, cfg, sizes, quad, include_log_sigma=True)
    information = app._compute_mml_observed_information(point, idx, cfg, sizes, quad, include_log_sigma=True)
    assert nll == pytest.approx(result["nll"], abs=1e-10)
    np.testing.assert_allclose(gradient, result["gradient"], atol=1e-9, rtol=0)
    np.testing.assert_allclose(information, result["information"], atol=5e-8, rtol=0)
    if slope_facet == "Criterion":
        assert np.max(abs(gradient-result["gradient"])) < np.max(abs(coarse_gradient-result["gradient"]))/100


def test_native_gpcm_keeps_the_existing_shared_step_slope_owner_scope():
    import streamlit_app as app
    _, cfg, reference = problem(1, "Rater", step_facet="Criterion")
    cfg.update(method="MML", n_person=reference.n,
        facet_specs={f: app.build_facet_constraint(levels, centered=f != cfg["noncenter_facet"])
                     for f, levels in cfg["facet_levels"].items()})
    with pytest.raises(ValueError, match="slope_facet == step_facet"):
        app.build_param_sizes(cfg)
