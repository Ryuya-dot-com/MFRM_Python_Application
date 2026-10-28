"""Native-scale uncertainty and replicate-cluster population MCSE controls."""
from dataclasses import replace

import numpy as np
import pytest

import streamlit_app as app
from tests.test_mfrm_recovery_design import design
from validation.mfrm_recovery_design import generate_population_sample
from validation.mfrm_recovery_metrics import (calibration_coordinates, candidate_wald_intervals,
    population_truth_quartile, replicate_ratio, summarize_replicate_errors)
from validation.mml_free_sd_stationarity_adapter import prepare_app_free_sd_problem_from_data
from validation.mml_gpcm_joint_check import prepare_problem


@pytest.mark.parametrize('model', ['RSM', 'PCM', 'GPCM'])
@pytest.mark.parametrize('categories', [2, 4])
def test_exact_native_expansion_and_jacobian_match_app(model, categories):
    spec = replace(design(model, persons=12), rating_min=0)
    if categories == 2: spec = replace(spec, steps=((0.,),)*(1 if model == 'RSM' else 3))
    data = generate_population_sample(spec, 741)['data']
    if model == 'GPCM': _, problem = prepare_problem(data, 31, rating_max=categories-1)
    else:
        _, problem = prepare_app_free_sd_problem_from_data(data, person_col='Person',
            facet_cols=['Rater', 'Task', 'Criterion'], score_col='Score', rating_min=0, rating_max=categories-1,
            model=model, noncenter_facet='Criterion', step_facet=None if model == 'RSM' else 'Criterion')
    point = np.r_[np.linspace(-.3, .4, len(problem.structural_start)), np.log(1.2)]
    rows = calibration_coordinates(problem.config, point)
    native = app.expand_params(point[:-1], problem.sizes, problem.config)
    expected = list(np.concatenate([native['facets'][f] for f in problem.config['facet_names']]))
    expected.extend(native['steps'] if model == 'RSM' else native['steps_mat'].ravel())
    if model == 'GPCM': expected.extend(np.r_[native['log_slopes'], native['slopes']])
    expected.extend([np.log(1.2), 1.2])
    np.testing.assert_allclose([r['value'] for r in rows], expected, atol=1e-15)
    jacobian = np.array([r['jacobian'] for r in rows]); numeric = []
    for delta in np.eye(len(point))*1e-5:
        plus = calibration_coordinates(problem.config, point+delta)
        minus = calibration_coordinates(problem.config, point-delta)
        numeric.append([(a['value']-b['value'])/2e-5 for a,b in zip(plus,minus)])
    np.testing.assert_allclose(jacobian, np.array(numeric).T, atol=1e-10, rtol=0)
    # Criterion's common location must remain visible, unlike mean alignment.
    assert abs(np.mean([r['value'] for r in rows if r['kind']=='Criterion'])) > .01
    covariance = .02*np.eye(len(point))+.01*np.ones((len(point), len(point)))
    intervals = candidate_wald_intervals(rows, covariance)
    for row in intervals:
        if row['recovery_comparable']:
            assert row['se']**2 == pytest.approx(row['jacobian']@covariance@row['jacobian'])
        else: assert row['se'] is row['lower'] is row['upper'] is None
        if 'log_value' in row:
            assert row['lower'] > 0
            assert np.log(row['lower'])+np.log(row['upper']) == pytest.approx(2*row['log_value'])
    if categories == 2:
        assert all(not r['recovery_comparable'] for r in intervals if r['kind']=='Step')
    with pytest.raises(ValueError, match='full joint covariance'):
        candidate_wald_intervals(rows, covariance[:-1, :-1])


def test_ratio_mcse_uses_replicates_with_random_bin_sizes_and_not_independent_rows():
    a = replicate_ratio([1., 9.], [1., 3.])
    assert a['estimate'] == 2.5  # Equal-replicate mean would be 2, a different estimand.
    assert a['mcse'] == pytest.approx(.75)
    b = replicate_ratio([20., 180.], [20., 60.])
    assert b['estimate'] == a['estimate']
    assert b['mcse'] == pytest.approx(a['mcse'])
    c = replicate_ratio([1., 9., 0.], [1., 3., 0.])
    assert c['planned_replicates'] == 3 and c['empty_replicates'] == 1
    assert c['mcse'] == pytest.approx(3*np.sqrt(3)/8)


def test_failed_fits_missing_intervals_and_unknown_truth_keep_distinct_denominators():
    result = summarize_replicate_errors([[1.], [3.,3.,3.], [np.nan, np.nan], None],
        [[1.], [0.,0.,0.], [np.nan,np.nan], None])
    assert result['bias']['estimate'] == 2.5
    assert result['rmse'] == pytest.approx(np.sqrt(7))
    assert result['coverage_given_interval']['estimate'] == .25
    assert result['coverage_per_known_target']['estimate'] == pytest.approx(1/6)
    assert result['interval_availability']['estimate'] == pytest.approx(4/6)
    assert result['truth_unknown_replicates'] == 1
    assert not result['all_planned_truth_known']
    assert result['bias']['planned_replicates'] == 4
    assert not result['scientific_inference_ready']


def test_no_intervals_or_only_one_informative_replicate_do_not_produce_fake_precision():
    a = summarize_replicate_errors([[1.], [np.nan]], [[np.nan], [np.nan]])
    assert a['bias']['mcse'] is None
    assert a['coverage_given_interval']['estimate'] is None
    assert a['coverage_per_known_target']['estimate'] == 0.
    b = summarize_replicate_errors([None, None], [None, None])
    assert b['bias']['estimate'] is b['rmse'] is b['mcse_rmse'] is None
    assert b['truth_unknown_replicates'] == 2


def test_population_quartiles_use_fixed_cutpoints_and_preserve_ties():
    theta = np.array([-3., -.6744897501960817, 0., 0., .8, 3.])
    np.testing.assert_array_equal(population_truth_quartile(theta, 1.), [1,2,3,3,4,4])
    np.testing.assert_array_equal(population_truth_quartile(theta*1.3, 1.3), population_truth_quartile(theta, 1.))
    for value in theta:
        assert population_truth_quartile([value], 1.)[0] == population_truth_quartile(theta, 1.)[np.flatnonzero(theta==value)[0]]


@pytest.mark.parametrize('sums,counts', [([1], [0]), ([1], [-1]), ([np.nan], [1]), ([1], [1,2])])
def test_invalid_ratio_records_fail_closed(sums, counts):
    with pytest.raises(ValueError): replicate_ratio(sums, counts)


def test_a_covered_interval_without_an_estimate_is_invalid():
    with pytest.raises(ValueError): summarize_replicate_errors([[np.nan]], [[1.]])


def test_unrepresentable_log_interval_is_not_returned_as_a_valid_interval():
    row = dict(value=1., jacobian=np.ones(1), log_value=0., log_jacobian=np.ones(1))
    with pytest.raises(ValueError, match='finite representable'):
        candidate_wald_intervals([row], np.array([[1e10]]))
