"""Independent finite-sum controls prompted by the Zotero literature audit.

See docs/mml_population_sd_uncertainty_audit.md for sources and limits.
The oracle uses explicit probabilities, not the application's likelihood or
derivative helpers. These are algebra checks, not coverage studies.
"""

import numpy as np
import pandas as pd
import pytest

import streamlit_app as app


@pytest.mark.parametrize("model", ["RSM", "PCM", "GPCM"])
@pytest.mark.parametrize("sign", [-1, 1])
def test_marginal_information_subtracts_posterior_score_variance(model, sign):
    config = {
        "model": model, "method": "MML", "n_cat": 3, "n_person": 2,
        "facet_names": ["Criterion"], "facet_levels": {"Criterion": ["C1", "C2"]},
        "facet_specs": {"Criterion": app.build_facet_constraint(["C1", "C2"], centered=True)},
        "facet_signs": {"Criterion": sign}, "step_facet": "Criterion",
        "slope_facet": "Criterion",
    }
    sizes = app.build_param_sizes(config)
    par = np.array([.35, -.4] if model == "RSM" else [.35, -.4, .2])
    if model == "GPCM":
        par = np.append(par, np.log(2.))
    idx = {"person": np.array([0, 0, 1, 1]), "score_k": np.array([0, 2, 1, 2]),
           "facets": {"Criterion": np.array([0, 1, 0, 1])},
           "step_idx": np.array([0, 1, 0, 1]), "slope_idx": np.array([0, 1, 0, 1]),
           "weight": np.array([1., 2., .5, 1.5])}
    # Deliberate three-point mixture: finite-sum differentiation, not a claim
    # about the accuracy of a normal-distribution quadrature rule.
    quad = {"nodes": np.array([-1.1, .2, 1.7]), "weights": np.array([.2, .5, .3])}
    category = np.arange(3)
    slopes = np.array([2., .5]) if model == "GPCM" else np.ones(2)
    cumulative_steps = np.array([[0., -.4, 0.], [0., -.4 if model == "RSM" else .2, 0.]])
    direction = sign * np.array([1., -1.])
    eta = quad["nodes"][:, None] + .35 * direction
    logits = slopes[None, :, None] * (eta[:, :, None] * category - cumulative_steps)
    probs = np.exp(logits - logits.max(axis=2, keepdims=True))
    probs /= probs.sum(axis=2, keepdims=True)
    expected = probs @ category
    variance = probs @ category**2 - expected**2
    posteriors, means, sds = [], [], []
    marginal_ll = mean_score = complete_information = missing_information = 0.
    for person in range(2):
        rows = np.flatnonzero(idx["person"] == person)
        weights, scores = idx["weight"][rows], idx["score_k"][rows]
        likelihood = np.prod(probs[:, np.arange(2), scores] ** weights, axis=1)
        mass = quad["weights"] * likelihood
        posterior = mass / mass.sum()
        score = np.sum(weights * slopes * direction * (scores - expected), axis=1)
        info = np.sum(weights * (slopes * direction)**2 * variance, axis=1)
        posterior_score = posterior @ score
        marginal_ll += np.log(mass.sum())
        mean_score += posterior_score
        complete_information += posterior @ info
        missing_information += posterior @ (score - posterior_score)**2
        mean = posterior @ quad["nodes"]
        means.append(mean)
        sds.append(np.sqrt(posterior @ (quad["nodes"] - mean)**2))
        posteriors.append(posterior)

    params = app.expand_params(par, sizes, config)
    actual_post, actual_ll = app._e_step_posteriors(idx, config, params, quad)
    nll, gradient = app.mfrm_loglik_mml_value_grad(par, idx, config, sizes, quad)
    hessian = app._compute_mml_observed_information(par, idx, config, sizes, quad)
    np.testing.assert_allclose(actual_post, posteriors, atol=1e-13)
    assert actual_ll == pytest.approx(marginal_ll, abs=1e-12)
    assert nll == pytest.approx(-marginal_ll, abs=1e-12)
    assert gradient[0] == pytest.approx(-mean_score, abs=1e-12)
    assert missing_information > 1e-7
    assert hessian[0, 0] == pytest.approx(complete_information - missing_information, abs=1e-8)
    # The M-step holds posterior weights fixed and therefore lacks the
    # subtraction. It cannot replace the observed marginal information.
    delta = np.zeros_like(par)
    delta[0] = 1e-5
    gradients = [app._m_step_expected_ll_value_grad(
        par + d, idx, config, sizes, quad, actual_post)[1][0] for d in (delta, -delta)]
    assert (gradients[0] - gradients[1]) / 2e-5 == pytest.approx(complete_information, abs=1e-8)
    eap = app.compute_person_eap(idx, config, params, quad)
    np.testing.assert_allclose(eap.Estimate, means, atol=1e-13)
    np.testing.assert_allclose(eap.SD, sds, atol=1e-13)


@pytest.mark.parametrize("model", ["RSM", "PCM", "GPCM"])
@pytest.mark.parametrize("population_mean", [0., 1.25])
def test_grid_map_maximizes_density_not_quadrature_mass(model, population_mean):
    # Four successes and one failure, with N(mu, 0.7^2) prior. For GH5 the
    # mass-maximizing node is mu+0.948938..., but the density-maximizing node
    # is mu. Neither is asserted to be the continuous posterior mode.
    config = {"method": "MML", "model": model, "n_person": 1, "n_cat": 2,
              "facet_names": ["Criterion"], "population_prior_sd": .7,
              "population_model": {"enabled": True, "X": np.ones((1, 1))}}
    params = {"facets": {"Criterion": np.array([population_mean])},
              "population": np.array([population_mean]), "steps": np.zeros(1),
              "steps_mat": np.zeros((1, 1)), "slopes": np.ones(1)}
    idx = {"person": np.zeros(5, dtype=int), "score_k": np.array([1, 1, 1, 1, 0]),
           "facets": {"Criterion": np.zeros(5, dtype=int)},
           "step_idx": np.zeros(5, dtype=int), "slope_idx": np.zeros(5, dtype=int)}
    quad = app.gauss_hermite_normal(5, sd=.7)
    nodes = quad["nodes"]
    log_likelihood = -4*np.logaddexp(0, -nodes) - np.logaddexp(0, nodes)
    log_density = log_likelihood - .5*(nodes/.7)**2
    log_mass = log_likelihood + np.log(quad["weights"])
    assert np.argmax(log_mass) != np.argmax(log_density)
    scores = app.compute_person_posterior_outputs(idx, config, params, quad)["scores"].iloc[0]
    assert scores.PosteriorMAPGrid == pytest.approx(population_mean, abs=1e-14)
    mass = np.exp(log_mass - log_mass.max())
    mass /= mass.sum()
    mean = mass @ nodes
    assert scores.PosteriorMean == pytest.approx(population_mean + mean)
    assert scores.PosteriorSD == pytest.approx(np.sqrt(mass @ (nodes - mean)**2))


def test_eap_uncertainty_names_conditioning_and_normal_interval():
    measures = pd.DataFrame({"Facet": ["Person"], "Level": ["P1"], "Estimate": [.2], "SE": [.7]})
    result = {"config": {"method": "MML", "model": "RSM", "facet_names": []}}
    annotated, _ = app.annotate_measure_uncertainty(measures, result)
    row = annotated.iloc[0]
    assert "conditional on fitted" in row.SE_Detail
    assert "Normal approximation" in row.CI_Method
    assert "calibration uncertainty" in row.UncertaintyCaution
