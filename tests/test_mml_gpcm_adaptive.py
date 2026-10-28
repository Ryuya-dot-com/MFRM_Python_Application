"""Finite adaptive derivatives, independently differentiated and integrated."""

import numpy as np
import pandas as pd
import pytest

from tests.test_mml_gpcm_integral import problem
from validation.mml_gpcm_adaptive import AdaptiveGpcmLikelihood, normal_log_rule
from validation.mml_grouped_fit import normal_rule


@pytest.mark.parametrize('sign', [-1, 1])
@pytest.mark.parametrize('order', [3, 31])
def test_finite_adaptive_gradient_includes_mode_scale_and_sd_derivatives(sign, order):
    data, config, reference = problem(sign)
    # Preserve unequal row counts and repeated cells inside a Person.
    data = pd.concat([data, data.iloc[::4]], ignore_index=True)
    evaluator = AdaptiveGpcmLikelihood(data, config)
    point = np.r_[np.linspace(-.5, .6, reference.linear_d), .8, -.3, np.log(1.7)]
    rule = normal_rule(order)
    result = evaluator.evaluate(point, *rule)
    numeric = []
    for delta in np.eye(len(point))*2e-5:
        plus = evaluator.evaluate(point+delta, *rule)
        minus = evaluator.evaluate(point-delta, *rule)
        numeric.append((plus['nll']-minus['nll'])/4e-5)
    np.testing.assert_allclose(result['gradient'], numeric, atol=3e-7, rtol=0)
    assert np.max(abs(result['mode_score'])) <= 1e-10
    if order == 3:
        # Differentiating while freezing the moving rule is detectably wrong.
        assert np.max(abs(result['gradient']-result['fixed_z_gradient'])) > 1e-3


def test_adaptive_gradient_and_information_converge_to_independent_continuous_integral():
    data, config, reference = problem(1)
    evaluator = AdaptiveGpcmLikelihood(data, config)
    point = np.r_[np.linspace(-.3, .4, reference.linear_d), .5, -.2, np.log(1.2)]
    rule = normal_rule(121)
    actual = evaluator.evaluate(point, *rule)
    target = reference.continuous(point, information=True)
    for key in ['nll', 'gradient', 'eap', 'sd']:
        np.testing.assert_allclose(actual[key], target[key], atol=1e-8, rtol=0)
    matrices = []
    for step in [1e-4, 5e-5]:
        columns = [(evaluator.evaluate(point+d, *rule)['gradient']-
                    evaluator.evaluate(point-d, *rule)['gradient'])/(2*step)
                   for d in np.eye(len(point))*step]
        matrices.append(np.asarray(columns).T)
    np.testing.assert_allclose(matrices[0], matrices[1], atol=2e-7, rtol=0)
    np.testing.assert_allclose(matrices[1], matrices[1].T, atol=1e-7, rtol=0)
    np.testing.assert_allclose(matrices[1], target['information'], atol=2e-7, rtol=0)


def test_zero_tail_weights_and_omitting_moments_preserve_the_finite_objective():
    data, config, reference = problem()
    evaluator = AdaptiveGpcmLikelihood(data, config)
    point = np.r_[np.zeros(reference.linear_d), .7, -.2, np.log(.8)]
    z, w = normal_rule(481)
    assert np.any(w == 0)
    with np.errstate(divide='raise', invalid='raise'):
        full = evaluator.evaluate(point, z, w)
        positive = evaluator.evaluate(point, z[w > 0], w[w > 0])
        fast = evaluator.evaluate(point, z, w, moments=False)
    assert 'eap' not in fast and 'sd' not in fast
    for key in ['nll', 'gradient', 'eap', 'sd']:
        np.testing.assert_allclose(full[key], positive[key], atol=2e-12, rtol=0)
    for key in fast:
        np.testing.assert_array_equal(full[key], fast[key])


def test_preserves_scope_and_rejects_nonrepresentable_or_invalid_inputs():
    data, config, reference = problem()
    evaluator = AdaptiveGpcmLikelihood(data, config)
    point = np.zeros(reference.d+1)
    for invalid in [point[:-1], np.full_like(point, np.nan)]:
        with pytest.raises(ValueError, match='coordinates'):
            evaluator.evaluate(invalid, *normal_rule(15))
    for logs in [1000., -1000.]:
        bad = point.copy(); bad[-1] = logs
        with pytest.raises(FloatingPointError):
            evaluator.evaluate(bad, *normal_rule(15))
    with pytest.raises(ValueError, match='quadrature'):
        evaluator.evaluate(point, [-1., 1.], [-.5, 1.5])
    with pytest.raises(ValueError, match='unit-weight'):
        AdaptiveGpcmLikelihood(data.assign(Weight=2.), config)
    with pytest.raises(ValueError, match='common step/slope'):
        AdaptiveGpcmLikelihood(data, dict(config, slope_facet='Rater'))
    with pytest.raises(ValueError, match='GPCM'):
        AdaptiveGpcmLikelihood(data, dict(config, model='PCM'))


@pytest.mark.parametrize('order', [3, 31, 481, 1921])
def test_logarithmic_normal_rule_matches_representable_weights_and_moments(order):
    z, logw = normal_log_rule(order)
    expected_z, expected_w = normal_rule(order)
    np.testing.assert_array_equal(z, expected_z)
    # Subnormal weights have few significant bits; check relative precision
    # where the original weights are normal, and moments of the full rule.
    normal = expected_w > 1e-280
    np.testing.assert_allclose(logw[normal], np.log(expected_w[normal]), atol=5e-11, rtol=0)
    # An n-point Gaussian rule is exact only through degree 2*n-1.
    degrees = [k for k in [0,2,4,6] if k < 2*order]
    np.testing.assert_allclose(np.exp(logw) @ np.stack([z**k for k in degrees], axis=1),
                               [1.,1.,3.,15.][:len(degrees)], atol=2e-11, rtol=0)
    assert np.isfinite(logw).all()
    assert not z.flags.writeable and not logw.flags.writeable


def test_underflowed_log_weight_matches_independent_high_precision_polynomial():
    from decimal import Decimal, localcontext
    from math import factorial
    order = 1921
    z, logw = normal_log_rule(order)
    assert normal_rule(order)[1][0] == 0
    with localcontext() as context:
        context.prec = 80
        x = Decimal.from_float(float(z[0]))
        previous, current = Decimal(1), x
        # Monic He recurrence, without the floating-point rescaling used above.
        for degree in range(2, order):
            previous, current = current, x*current-(degree-1)*previous
        target = Decimal(factorial(order-1)).ln()-Decimal(order).ln()-2*abs(current).ln()
    assert abs(logw[0]-float(target)) < 5e-11


def test_log_weight_entry_preserves_low_order_derivatives_and_input_guards():
    data, config, reference = problem()
    evaluator = AdaptiveGpcmLikelihood(data, config)
    point = np.r_[np.linspace(-.3, .4, reference.linear_d), .5, -.2, np.log(1.2)]
    ordinary = evaluator.evaluate(point, *normal_rule(31))
    logged = evaluator.evaluate_log_weights(point, *normal_log_rule(31))
    for key in ordinary:
        np.testing.assert_allclose(logged[key], ordinary[key], atol=1e-11, rtol=0)
    for bad in [True, 2, 3.5, '31']:
        with pytest.raises(ValueError, match='Quadrature'):
            normal_log_rule(bad)
    for weights in [[0.,0.], [-np.inf,-np.inf], [0.,np.nan], [0.,np.inf]]:
        with pytest.raises(ValueError, match='quadrature'):
            evaluator.evaluate_log_weights(point, [-1.,1.], weights)


def test_tiny_base_weights_are_not_discarded_after_adaptive_density_correction():
    from validation.mml_gpcm_integral import GpcmIntegral
    data = pd.DataFrame(dict(Person=['P0','P1'], Criterion=['C0','C1'], Score=[1,0]))
    config = dict(model='GPCM', facet_names=['Criterion'],
        facet_levels={'Criterion':['C0','C1']}, facet_signs={'Criterion':-1},
        noncenter_facet='Criterion', step_facet='Criterion', slope_facet='Criterion', n_cat=2, rating_min=0)
    # This finite point is inside the existing log-slope and population-SD
    # box. A sharp mode does not make the whole posterior narrowly Gaussian.
    point = np.array([19.,0.,2.9,np.log(9.5)])
    reference = GpcmIntegral(data,config).continuous(point,bound=16.,rel_tol=1e-12,abs_tol=1e-14)
    evaluator = AdaptiveGpcmLikelihood(data,config)
    rounded = evaluator.evaluate(point,*normal_rule(1921))
    logged = evaluator.evaluate_log_weights(point,*normal_log_rule(1921))
    assert np.max(abs(rounded['eap']-reference['eap'])) > .01
    assert np.max(abs(rounded['sd']-reference['sd'])) > .03
    for key in ['nll','gradient','eap','sd']:
        np.testing.assert_allclose(logged[key],reference[key],atol=2e-7,rtol=0)
