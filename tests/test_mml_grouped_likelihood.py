"""Exact sufficient-statistic reduction and full finite-GH information."""
from dataclasses import replace

import numpy as np
import pandas as pd
import pytest
from scipy.special import roots_hermitenorm

import streamlit_app as app
from tests.test_mfrm_recovery_design import design
from validation.mfrm_recovery_design import generate_population_sample
from validation.mml_grouped_likelihood import GroupedMmlLikelihood


def problem(model, assignment, duplicate=False):
    sample = generate_population_sample(replace(design(model, assignment, persons=6), rating_min=0), 473)
    data = sample['data']
    if duplicate:
        data = pd.concat([data, data.iloc[::3]], ignore_index=True)
    facets = ['Rater', 'Task', 'Criterion']
    config = dict(model=model, method='MML', n_person=6, n_cat=4, rating_min=0, rating_max=3,
        facet_names=facets, facet_levels={f: sample['declared_levels'][f] for f in facets},
        facet_signs={f: -1 for f in facets}, noncenter_facet='Criterion',
        step_facet=None if model == 'RSM' else 'Criterion', slope_facet='Criterion' if model == 'GPCM' else None)
    config['facet_specs'] = {f: app.build_facet_constraint(config['facet_levels'][f], centered=f != 'Criterion') for f in facets}
    sizes = app.build_param_sizes(config)
    idx = dict(person=np.array([sample['declared_levels']['Person'].index(p) for p in data['Person']]),
        score_k=data['Score'].to_numpy(), weight=np.ones(len(data)),
        facets={f: np.array([config['facet_levels'][f].index(v) for v in data[f]]) for f in facets})
    if model != 'RSM': idx['step_idx'] = idx['facets']['Criterion']
    if model == 'GPCM': idx['slope_idx'] = idx['facets']['Criterion']
    return data, config, sizes, idx, GroupedMmlLikelihood(data, config)


@pytest.mark.parametrize('model', ['RSM', 'PCM', 'GPCM'])
@pytest.mark.parametrize('assignment', ['regular', 'sparse', 'unbalanced'])
@pytest.mark.parametrize('duplicate', [False, True])
def test_grouped_objective_gradient_posterior_and_full_information_match_app(model, assignment, duplicate):
    _, config, sizes, idx, grouped = problem(model, assignment, duplicate)
    point = np.r_[np.linspace(-.4, .5, sum(sizes.values())), np.log(1.2)]
    quad = app.gauss_hermite_normal(31, sd=1.2)
    result = grouped.evaluate(point, quad['nodes']/1.2, quad['weights'], information=True)
    nll, gradient = app.mfrm_loglik_mml_value_grad(point[:-1], idx, config, sizes, quad, include_log_sigma=True)
    information = app._compute_mml_observed_information(point, idx, config, sizes, quad, include_log_sigma=True)
    posterior, _ = app._e_step_posteriors(idx, config, app.expand_params(point[:-1], sizes, config), quad)
    assert result['nll'] == pytest.approx(nll, abs=1e-11)
    np.testing.assert_allclose(result['gradient'], gradient, atol=1e-10, rtol=0)
    np.testing.assert_allclose(result['information'], information, atol=5e-7, rtol=0)
    np.testing.assert_allclose(result['posterior'], posterior, atol=2e-14, rtol=0)
    assert result['cells'] <= 18 < result['observations']


def test_analytic_gpcm_information_matches_finite_differences_of_grouped_gradient():
    _, _, sizes, _, grouped = problem('GPCM', 'unbalanced', duplicate=True)
    point = np.r_[np.linspace(-.5, .4, sum(sizes.values())), np.log(.8)]
    z, w = roots_hermitenorm(41); w /= np.sqrt(2*np.pi)
    result = grouped.evaluate(point, z, w, information=True)
    columns = []
    for delta in np.eye(len(point))*1e-5:
        plus = grouped.evaluate(point+delta, z, w); minus = grouped.evaluate(point-delta, z, w)
        columns.append((plus['gradient']-minus['gradient'])/2e-5)
    np.testing.assert_allclose(result['information'], np.array(columns).T, atol=2e-7, rtol=0)
    assert np.max(abs(result['information'][:-1, -1])) > 1e-3


def test_zero_quadrature_weights_have_exactly_zero_mass_without_changing_the_objective():
    _, _, sizes, _, grouped = problem('PCM', 'sparse')
    point = np.r_[np.zeros(sum(sizes.values())), 0.]
    z, w = roots_hermitenorm(641); w /= np.sqrt(2*np.pi)
    assert np.any(w == 0)
    full = grouped.evaluate(point, z, w, information=True)
    positive = grouped.evaluate(point, z[w > 0], w[w > 0], information=True)
    for key in ['nll', 'gradient', 'information', 'eap', 'sd']:
        np.testing.assert_allclose(full[key], positive[key], atol=1e-12, rtol=0)
    np.testing.assert_array_equal(full['posterior'][:, w == 0], 0.)


def test_grouping_does_not_merge_different_persons_before_missing_information():
    data, cfg, _, _, grouped = problem('RSM', 'regular')
    merged = data.copy(); merged['Person'] = 'one-person'
    other = GroupedMmlLikelihood(merged, dict(cfg, n_person=1))
    point = np.r_[np.linspace(-.4, .5, grouped.dimension-1), np.log(1.2)]
    z, w = roots_hermitenorm(31); w /= np.sqrt(2*np.pi)
    a, b = grouped.evaluate(point, z, w, information=True), other.evaluate(point, z, w, information=True)
    assert abs(a['nll']-b['nll']) > .01
    assert np.max(abs(a['information']-b['information'])) > .01


def test_grouped_reference_rejects_power_weights_and_invalid_quadrature():
    data, cfg, _, _, grouped = problem('RSM', 'regular')
    with pytest.raises(ValueError, match='unit-weight'):
        GroupedMmlLikelihood(data.assign(Weight=2.), cfg)
    with pytest.raises(ValueError, match='quadrature'):
        grouped.evaluate(np.zeros(grouped.dimension), [-1., 1.], [.5, -.5])


@pytest.mark.parametrize('model', ['RSM', 'PCM', 'GPCM'])
def test_optimizer_can_omit_moments_without_changing_joint_derivatives(model):
    _, _, _, _, grouped = problem(model, 'unbalanced', duplicate=True)
    point = np.r_[np.linspace(-.5,.5,grouped.dimension-1), np.log(1.7)]
    z,w = roots_hermitenorm(41);w /= np.sqrt(2*np.pi)
    full = grouped.evaluate(point,z,w,information=True)
    fast = grouped.evaluate(point,z,w,information=True,moments=False)
    assert 'eap' not in fast and 'sd' not in fast
    assert set(fast) == set(full)-{'eap','sd'}
    for key in fast: np.testing.assert_array_equal(fast[key],full[key])


def test_quadrature_reuse_preserves_nodes_zero_weights_and_input_guards():
    from validation.mml_grouped_fit import normal_rule
    z,w = normal_rule(641)
    expected_z,expected_w = roots_hermitenorm(641);expected_w /= np.sqrt(2*np.pi)
    np.testing.assert_array_equal(z,expected_z);np.testing.assert_array_equal(w,expected_w)
    assert normal_rule(np.int64(641))[0] is z
    for values in [z,w]:
        with pytest.raises(ValueError): values[0] = 42.
    for invalid in [True,False,2,641.,'641']:
        with pytest.raises(ValueError,match='Quadrature'): normal_rule(invalid)
