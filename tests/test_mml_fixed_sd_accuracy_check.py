"""Independent controls for the offline fixed-SD integration decomposition."""
import json

import numpy as np
import pandas as pd
import pytest
from scipy.special import ndtr, roots_hermitenorm

import streamlit_app as app
from validation import mml_fixed_sd_accuracy_check as check
from validation import mml_fixed_sd_reference_error_audit as error_audit


@pytest.mark.parametrize('facets,positive', [
    (['Rater', 'Task', 'Criterion'], ['Rater']),
    (['Criterion', 'Task', 'Rater'], ['Criterion', 'Task']),
])
def test_signed_coordinate_map_moments_gradient_and_restore(facets, positive, monkeypatch, tmp_path):
    data = pd.DataFrame([dict(Person=f'{p:03d}', Rater=f'R{r}', Task=f'T{t}',
        Criterion=f'C{c}', Score=(p+2*r+t+c)%4)
        for p in range(8) for r in range(2) for t in range(2) for c in range(2)]).iloc[1:]
    fit = app.mfrm_estimate(data=data, person_col='Person', facet_cols=facets,
        score_col='Score', model='PCM', method='MML', step_facet='Criterion', noncenter_facet='Criterion',
        positive_facets=positive, population_prior_sd=.6, quad_points=121, maxit=1)
    # Deliberately away from a solution: nonzero signs cannot hide in a symmetric fit.
    point = np.array([.2, -.1, .3, -.25, .35, -.15, .12, -.23])
    cfg = fit['config']; sizes = app.build_param_sizes(cfg)
    idx = app.build_indices(fit['prep'], step_facet='Criterion'); quad = app.make_mml_quadrature(cfg)
    fit['opt'].x = point; fit['params'] = app.expand_params(point, sizes, cfg)
    value, gradient = app.mfrm_loglik_mml_value_grad(point, idx, cfg, sizes, quad)
    fit['opt'].fun = value
    assets = app._analysis_input_assets(fit)
    prepared = check.restore(assets['data.csv'], assets['analysis_snapshot.json'])
    transform = check.reference_map(prepared); job = prepared['reference_job']
    np.testing.assert_allclose(transform@point, np.asarray(job['xsi'])[check.PERM], atol=0, rtol=0)
    np.testing.assert_array_equal(transform.T@transform, np.eye(8))
    reference = check.Reference(job)
    z, w = roots_hermitenorm(121); w /= np.sqrt(2*np.pi)
    finite = reference.finite(transform@point, z, w)
    scores = app.compute_person_eap(idx, cfg, fit['params'], quad)
    np.testing.assert_allclose(finite['nll'], value, atol=1e-10, rtol=0)
    np.testing.assert_allclose(transform.T@finite['gradient'], gradient, atol=1e-10, rtol=0)
    np.testing.assert_allclose(finite['eap'], scores.Estimate, atol=1e-10, rtol=0)
    np.testing.assert_allclose(finite['sd'], scores.SD, atol=1e-10, rtol=0)
    continuous = reference.continuous(transform@point, **check.PROTOCOL['continuous'][-1])
    assert max(check.differences(finite, continuous).values()) < 1e-9
    def vg(p):
        v = reference.continuous(transform@p, **check.PROTOCOL['continuous'][-1])
        return v['nll'], transform.T@v['gradient']
    audit = check.audit_joint_gradient(lambda p: vg(p)[0], vg, point)
    assert audit.maximum_absolute_difference < 1e-7
    assert job['person_labels'][0] == '000'
    assert sum(v is None for row in job['responses'] for v in row) == 1
    invalid = json.loads(assets['analysis_snapshot.json']); invalid['config']['estimate_population_sd'] = True
    with pytest.raises(ValueError, match='fixed population SD'):
        check.restore(assets['data.csv'], json.dumps(invalid))
    if facets[0] == 'Rater':
        # Failed optimization must leave all six attempted conditions visible.
        monkeypatch.setattr(check, 'polish_fixed_sd', lambda *a, **kw: dict(completed=False,
            candidate_coordinates=None, stages=[], failure={'stage': 'raw', 'message': 'injected'}))
        report = check.check(prepared, lambda name, obj: check.dump(tmp_path/name, obj))
        assert report['expected_trajectories'] == len(report['trajectories']) == 6
        assert report['evaluated_trajectories'] == 0 and report['cross_start_checks'] == {}
        assert all(r['run']['failure']['message'] == 'injected' for r in report['trajectories'])
        assert len(list(tmp_path.glob('*_optimizer.json'))) == 6
        assert len(list(tmp_path.glob('*_same_point.json'))) == 3
        assert not report['inference_ready'] and not report['original_result_replaced']
        np.testing.assert_array_equal(prepared['coordinates'], point)
        # The follow-up audit rejects a mismatched saved value and retains the
        # six unavailable endpoints instead of dropping them from its denominator.
        inputs = tmp_path/'inputs'; inputs.mkdir()
        outputs = tmp_path/'audit'; outputs.mkdir()
        for name, content in assets.items():
            (inputs/name).write_bytes(content.encode() if isinstance(content, str) else content)
        report['original_continuous']['values'][0]['nll'] += 1.
        check.dump(inputs/'summary.json', report)
        check.dump(inputs/'runtime.json', prepared['implementation'])
        audited = error_audit.audit_case(inputs, lambda name, obj: check.dump(outputs/name, obj))
        assert audited['expected_points'] == len(audited['records']) == 7
        assert 'changed from the retained' in audited['records'][0]['failure']['message']
        assert all(not r['completed'] and not r['localized_diagnostics_meet_thresholds'] for r in audited['records'])
        assert not audited['inference_ready'] and not audited['prior_flags_replaced'] and not audited['refitted']


def test_tail_moment_allowance_and_signed_cancellation():
    # No responses: the posterior is the prior. Exact truncated-normal moments
    # expose a dropped sigma factor or a tail mass mistaken for a moment bound.
    bound, sigma = 4., 1.7
    tail = 2*ndtr(-bound); phi = np.exp(-bound*bound/2)/np.sqrt(2*np.pi)
    truncated_sd = sigma*np.sqrt(1-2*bound*phi/(1-tail))
    value = dict(eap=[0.], sd=[truncated_sd], numeric_relative_mass_error_sum=0.,
                 tail_relative_mass_bound_sum=tail/(1-tail))
    allowance = check.moment_allowances(value, bound, sigma)
    assert 0 < sigma-truncated_sd <= allowance['sd'] < .01
    # Coordinatewise changes can cancel, so the maxima must never be added.
    baseline = dict(nll=10., eap=np.array([0., 0.]), sd=np.array([1., 1.]))
    integral = dict(nll=12., eap=np.array([1., -1.]), sd=np.array([1.2, .8]))
    finite = dict(nll=11., eap=np.array([0., .5]), sd=np.array([.9, 1.1]))
    for result in check.decompose(finite, integral, baseline).values():
        assert result['identity_residual'] < 1e-15
        assert result['maximum_absolute']['total'] < sum(result['maximum_absolute'][k]
                                                       for k in ('integration', 'recalibration'))


def test_reference_capture_preserves_calculation_and_separates_person_and_total_errors():
    # Repeated identical prior-only posteriors isolate aggregation from the
    # integrator. This is an analytic control, not a supported fitted data set.
    settings = dict(bound=4., rel_tol=1e-10, abs_tol=1e-12)
    sigma = 1.7; tail = 2*ndtr(-settings['bound'])
    phi = np.exp(-settings['bound']**2/2)/np.sqrt(2*np.pi)
    expected_sd = sigma*np.sqrt(1-2*settings['bound']*phi/(1-tail))
    results = []
    for people in (1, 32):
        reference = check.Reference(dict(sigma=sigma, responses=[[None]*8]*people))
        direct = reference.continuous(np.zeros(8), **settings)
        inspected = error_audit.inspect_reference(reference, np.zeros(8), settings)
        assert max(check.differences(direct, inspected['value']).values()) == 0.
        assert direct['quadrature_neval'] == inspected['value']['quadrature_neval']
        np.testing.assert_allclose(inspected['scaled_integrals'][:, 0], np.sqrt(2*np.pi)*(1-tail), atol=1e-12, rtol=0)
        np.testing.assert_allclose(inspected['value']['sd'], expected_sd, atol=1e-12, rtol=0)
        assert inspected['relative_error_estimate'].sum() == direct['numeric_relative_mass_error_sum']
        assert inspected['total_nll_allowance'] == inspected['pooled']['nll']
        for a in inspected['per_person']:
            assert 0 < sigma-expected_sd <= a['sd'] <= inspected['pooled']['sd']
        results.append(inspected)
    # Same integrand per person: the shared absolute max-norm error stays fixed,
    # while the sum-based relative mass error (for total NLL) scales by 32.
    assert results[0]['shared_absolute_error_estimate'] == results[1]['shared_absolute_error_estimate']
    np.testing.assert_allclose(results[1]['relative_error_estimate'], results[0]['relative_error_estimate'][0], atol=0, rtol=0)
    assert results[1]['value']['numeric_relative_mass_error_sum'] == 32*results[0]['value']['numeric_relative_mass_error_sum']
    assert results[1]['total_nll_allowance'] > 31*results[0]['total_nll_allowance']
