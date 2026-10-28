"""Product-path numerical evidence; these tests grant no inference eligibility."""
from copy import deepcopy
import pickle

import numpy as np
import pandas as pd
import pytest

from mfrm_app import engine
from mfrm_app import mml_fixed_sd_review as review


@pytest.fixture(scope='module', params=['RSM', 'PCM'])
def fit(request):
    data = pd.DataFrame([dict(Person=f'{p:03d}', Rater=f'R{r}', Task=f'T{t}',
                              Score=1+(p+2*r+t)%3)
                         for p in range(18) for r in range(3) for t in range(2)])
    # Nonzero score origin, positive facet sign, and a missing cell are deliberate.
    return engine.mfrm_estimate(data=data.iloc[1:], person_col='Person',
        facet_cols=['Rater', 'Task'], score_col='Score', model=request.param,
        method='MML', noncenter_facet='Task', step_facet='Task', positive_facets=['Rater'],
        population_prior_sd=.8, quad_points=15, maxit=60)


def test_actual_result_matches_separate_derivatives_and_full_target_covariance(fit, monkeypatch):
    original = pickle.dumps(fit)
    monkeypatch.setattr(engine, 'minimize', lambda *a, **k: pytest.fail('Review must not refit'))
    prepared = review.prepare_fixed_sd_check(fit)
    record = review.evaluate_fixed_sd_quadrature_check(prepared)
    assert pickle.dumps(fit) == original
    assert record['comparison']['accuracy_pass'] is None and not record['inference_ready']
    assert not record['comparison']['refitted']
    assert not record['structural_boundary']['no_finite_mle_certified']
    assert record['person_order'][0] == '000'
    assert [r['points'] for r in record['comparison']['rows']] == [15, 31, 61, 121]
    j = np.asarray(record['target_jacobian'])
    # Rater's third centered level uses minus both free effects, including covariance.
    np.testing.assert_array_equal(j[2, :2], [-1, -1])
    for row in record['comparison']['rows']:
        quad = engine.make_mml_quadrature(prepared['config'], row['points'])
        h = engine._compute_mml_observed_information(prepared['coordinates'], prepared['idx'],
                prepared['config'], prepared['sizes'], quad, rel_step=1e-4)
        np.testing.assert_allclose(row['information'], h, atol=2e-7, rtol=1e-7)
        v = np.linalg.inv(h)
        actual = row['information_review']
        assert actual['positive_definite'] and not actual['regularized']
        np.testing.assert_allclose(actual['target_se'], np.sqrt(np.diag(j @ v @ j.T)), atol=1e-8)
        assert actual['target_se'][2] == pytest.approx(np.sqrt(v[0,0]+v[1,1]+2*v[0,1]), abs=1e-8)
        assert abs(actual['target_se'][2] - np.sqrt(v[0,0]+v[1,1])) > .01
        # Sigma is FIXED. Inverting a matrix that estimates sigma would be a different SE.
        joint_h = prepared['problem'].evaluate(np.r_[prepared['coordinates'], np.log(.8)],
            quad['nodes']/.8, quad['weights'], information=True)['information']
        assert np.max(np.abs(v - np.linalg.inv(joint_h)[:-1,:-1])) > 1e-4
        g = np.asarray(row['gradient'])
        assert actual['curvature_scaled_gradient'] == pytest.approx(np.sqrt(g @ v @ g), abs=1e-8)
        expected = np.max(np.abs(j @ v @ g)/np.sqrt(np.diag(j @ v @ j.T)))
        assert actual['max_target_newton_over_se'] == pytest.approx(expected, abs=1e-8)
        person = engine.compute_person_posterior_outputs(prepared['idx'], prepared['config'], fit['params'], quad)['scores']
        np.testing.assert_allclose(row['eap'], person.PosteriorMean, atol=1e-12)
        np.testing.assert_allclose(row['posterior_sd'], person.PosteriorSD, atol=1e-12)
    assert record['comparison']['rows'][0]['max_target_relative_se_change'] == 0
    # Full JSON is finite and lossless; numerical status does not bypass public holds.
    import json
    assert json.loads(json.dumps(record, allow_nan=False)) == record
    from mfrm_app.output_qualification import mml_inference_hold_reason
    assert mml_inference_hold_reason(fit['config'])


def test_se_scaling_and_failure_do_not_repair_information():
    j = np.array([[1., 0.], [0., 1.], [-1., -1.]])
    result = review._information_review(np.diag([1e-8, 1e8]), np.array([1e-5, 1e-5]), j)
    assert result['curvature_scaled_gradient'] == pytest.approx(.1)
    assert result['max_target_newton_over_se'] == pytest.approx(.1)
    for h in (np.diag([1., 0.]), np.diag([1., -1e-10]), np.diag([1., np.nan])):
        failed = review._information_review(h, np.ones(2), j)
        assert not failed['positive_definite'] and not failed['regularized']
        assert failed['covariance'] is None and failed['target_se'] is None
        assert failed['max_target_newton_over_se'] is None


def test_synchronous_size_guard_runs_before_grouped_allocation(fit, monkeypatch):
    changed = deepcopy(fit)
    sizes = engine.build_param_sizes(changed['config'])
    changed['config']['facet_specs']['Task']['n_params'] += 65-sum(sizes.values())
    changed['opt'].x = np.zeros(65)
    monkeypatch.setattr(review.mml_joint, 'GroupedMmlLikelihood',
        lambda *args: pytest.fail('Oversized request reached grouped-array allocation'))
    with pytest.raises(ValueError, match='1–64'):
        review.prepare_fixed_sd_check(changed)


@pytest.mark.parametrize('change', ['sd', 'weights', 'scores', 'levels', 'params', 'penalty', 'anchors', 'q', 'regression'])
def test_mismatched_or_unsupported_results_are_rejected_before_evaluation(fit, change, monkeypatch):
    changed = deepcopy(fit)
    if change == 'sd': changed['config']['estimate_population_sd'] = True
    elif change == 'weights': changed['prep']['data'].iloc[0, changed['prep']['data'].columns.get_loc('Weight')] = .5
    elif change == 'scores': changed['prep']['data']['score_k'] = .5
    elif change == 'levels': changed['config']['facet_levels']['Rater'] = ['R1', 'R0', 'R2']
    elif change == 'params': changed['params']['facets']['Rater'][0] += .1
    elif change == 'penalty': changed['config']['facet_regularization'] = {'enabled': True}
    elif change == 'anchors': changed['config']['facet_specs']['Rater']['anchors'][0] = 0.
    elif change == 'q': changed['config']['quad_points'] = True
    elif change == 'regression': changed['config']['population_model'] = {'enabled': True}
    monkeypatch.setattr(engine, 'mfrm_loglik_mml_value_grad', lambda *a, **k: pytest.fail('Unsupported input evaluated'))
    with pytest.raises(ValueError):
        review.prepare_fixed_sd_check(changed)


def test_record_is_tied_to_source_runtime_data_and_returned_point(fit, monkeypatch):
    original = review.prepare_fixed_sd_check(fit)
    changed = deepcopy(fit)
    changed['opt'].message += ' changed'
    assert review.prepare_fixed_sd_check(changed)['input_id'] != original['input_id']
    monkeypatch.setattr(review.importlib.metadata, 'version', lambda _: 'different-runtime')
    assert review.prepare_fixed_sd_check(fit)['input_id'] != original['input_id']
    monkeypatch.undo()
    changed = deepcopy(fit)
    changed['opt'].fun += .1
    with pytest.raises(ValueError, match='Reported NLL'):
        review.evaluate_fixed_sd_quadrature_check(review.prepare_fixed_sd_check(changed))
    assert {'engine.py', 'mml_joint.py', 'mml_fixed_sd_review.py'} == set(original['implementation']['source_sha256'])


def test_separation_is_recorded_even_if_the_local_information_is_positive(fit):
    changed = deepcopy(fit)
    changed['prep']['data']['Score'] = 3
    changed['prep']['data']['score_k'] = 2
    prepared = review.prepare_fixed_sd_check(changed)
    prepared['reported_nll'] = engine.mfrm_loglik_mml_value_grad(prepared['coordinates'], prepared['idx'],
        prepared['config'], prepared['sizes'], engine.make_mml_quadrature(prepared['config']))[0]
    record = review.evaluate_fixed_sd_quadrature_check(prepared)
    assert record['structural_boundary']['no_finite_mle_certified']
    assert not record['inference_ready']


def test_confounded_design_does_not_receive_a_numerical_inverse(fit):
    changed = deepcopy(fit)
    data = changed['prep']['data']
    changed['prep']['data'] = data[(data.Rater == 'R0') == (data.Task == 'T0')].copy()
    prepared = review.prepare_fixed_sd_check(changed)
    prepared['reported_nll'] = engine.mfrm_loglik_mml_value_grad(prepared['coordinates'], prepared['idx'],
        prepared['config'], prepared['sizes'], engine.make_mml_quadrature(prepared['config']))[0]
    record = review.evaluate_fixed_sd_quadrature_check(prepared)
    assert record['structural_boundary']['structural_design_rank'] < len(prepared['coordinates'])
    for row in record['comparison']['rows']:
        assert row['information_review']['status'] == 'rank_deficient_design'
        assert row['information_review']['covariance'] is None
        assert row['information_review']['target_se'] is None
        assert row['max_target_relative_se_change'] is None
