"""Prepared scalar slices must preserve probabilities, sums and failure states."""
import numpy as np
import pytest

import streamlit_app as app
from validation.mml_population_study import json_bytes


@pytest.mark.parametrize('model', ['RSM', 'PCM', 'GPCM'])
@pytest.mark.parametrize('categories', [2, 6])
@pytest.mark.parametrize('case', ['unit', 'unweighted', 'fractional', 'zero_rows',
                                 'zero_all', 'negative', 'nonfinite', 'minimum',
                                 'maximum', 'limited'])
def test_prepared_cell_matches_existing_model_and_fit_exactly(model, categories, case):
    rng = np.random.default_rng(2026100901)
    eta = rng.normal(size=31)
    scores = rng.integers(categories, size=31)
    # Interleaved groups, an unused step level, and a different slope facet.
    step_idx = rng.choice([0, 2], size=31)
    slope_idx = rng.integers(2, size=31)
    slopes = np.array([.65, 1.6])
    cumulative = np.c_[np.zeros(3), np.cumsum(rng.normal(size=(3, categories-1)), axis=1)]
    weight = np.ones(31)
    bound = 10.
    if case == 'unweighted': weight = None
    if case == 'fractional': weight[::3] = .5
    if case == 'zero_rows': weight[::3] = 0.
    if case == 'zero_all': weight[:] = 0.
    if case == 'negative': weight[0] = -1.
    if case == 'nonfinite': weight[0] = np.nan
    if case == 'minimum': scores[:] = 0
    if case == 'maximum': scores[:] = categories-1
    if case == 'limited': bound = 1e-8
    if model == 'RSM':
        nll = lambda b: -app.loglik_rsm(eta+b, scores, cumulative[0], weight)
        probability = lambda b: app.category_prob_rsm(eta+b, cumulative[0])
    elif model == 'PCM':
        nll = lambda b: -app.loglik_pcm(eta+b, scores, cumulative, step_idx, weight)
        probability = lambda b: app.category_prob_pcm(eta+b, cumulative, step_idx)
    else:
        nll = lambda b: -app.loglik_gpcm(eta+b, scores, cumulative, step_idx, slopes, slope_idx, weight)
        probability = lambda b: app.category_prob_gpcm(eta+b, cumulative, step_idx, slopes, slope_idx)
    prepared_nll, prepared_probability = app._prepare_bias_cell_functions(
        model, eta, scores, weight, cumulative[0], cumulative, step_idx, slopes, slope_idx)
    for b in [-10., -.03, 0., .27, 10.]:
        np.testing.assert_array_equal(prepared_probability(b), probability(b), strict=True)
        assert json_bytes(prepared_nll(b)) == json_bytes(nll(b))
    fit_args = (scores, slopes[slope_idx] if model == 'GPCM' else np.ones(31),
                weight if weight is not None else np.ones(31), bound)
    expected = app._fit_conditional_bias_cell(nll, probability, *fit_args)
    actual = app._fit_conditional_bias_cell(prepared_nll, prepared_probability, *fit_args)
    assert json_bytes(actual) == json_bytes(expected)
    if case in ['zero_all', 'negative', 'nonfinite']:
        assert actual['status'] == 'invalid_input'
    elif case in ['minimum', 'maximum']:
        assert actual['status'] == 'no_finite_mle'
    elif case == 'limited':
        assert actual['status'] == 'search_range_limited'
    else:
        assert actual['status'] == 'interior'


def test_prepared_cells_remain_independent_when_calibration_changes():
    scores = np.array([0, 1, 0, 1, 1, 1])
    args = (scores, np.ones(6), None, np.array([[0., .4]]),
            np.zeros(6, dtype=int), None, None)
    first = app._prepare_bias_cell_functions('PCM', np.zeros(6), *args)
    before = first[1](.2)
    second = app._prepare_bias_cell_functions('PCM', np.ones(6), *args)
    assert not np.array_equal(before, second[1](.2))
    np.testing.assert_array_equal(before, first[1](.2))


@pytest.mark.parametrize('step_idx', [np.array([-1, 0, 0]), np.array([-1, -1, -1])])
def test_unmapped_step_rows_still_fail_probability_validation(step_idx):
    eta, scores, weights = np.zeros(3), np.array([0, 1, 1]), np.ones(3)
    steps = np.array([[0., .4]])
    nll, probability = app._prepare_bias_cell_functions(
        'PCM', eta, scores, weights, None, steps, step_idx, None, None)
    np.testing.assert_array_equal(probability(.2), app.category_prob_pcm(eta+.2, steps, step_idx))
    result = app._fit_conditional_bias_cell(nll, probability, scores, weights, weights, 10.)
    assert result['status'] == 'calculation_failed'
    assert result['error'] == 'ValueError: Invalid category probabilities in bias cell'
    assert np.isnan(result['p'])
