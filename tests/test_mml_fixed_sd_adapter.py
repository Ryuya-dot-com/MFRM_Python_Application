"""Fixed-SD candidates must preserve the fitted result and the exact finite rule."""
from copy import deepcopy
import pickle

import numpy as np
import pandas as pd
import pytest
from scipy.special import logsumexp, roots_hermitenorm

import streamlit_app as app
from mfrm_app import mml_stationarity as kernel
from mfrm_app.legacy_compat import pcm_native_runner as native
from validation import mml_free_sd_stationarity_adapter as adapter


def fit_example(sigma=1.1, missing=False):
    data = pd.DataFrame([dict(Person=f'{p:03d}', Rater=f'R{r}', Task=f'T{t}',
        Criterion=f'C{c}', Score=(p+2*r+t+c)%4)
        for p in range(16) for r in range(2) for t in range(2) for c in range(2)])
    if missing:
        data = data.iloc[1:]
    return app.mfrm_estimate(data=data, person_col='Person', facet_cols=['Rater','Task','Criterion'],
        score_col='Score', model='PCM', method='MML', step_facet='Criterion', noncenter_facet='Criterion',
        positive_facets=['Rater'], population_prior_sd=sigma, quad_points=121, maxit=60)


@pytest.fixture(scope='module')
def fitted():
    return fit_example()


def input_bytes(result):
    # The optimizer's local class is not picklable; preserve its public values separately.
    return pickle.dumps({k:v for k,v in result.items() if k!='opt'}), pickle.dumps(vars(result['opt']))


def test_disabled_is_a_noop_and_unsupported_inputs_stop_before_optimization(fitted, monkeypatch):
    before = input_bytes(fitted)
    with monkeypatch.context() as m:
        m.setattr(app, 'build_cross_engine_validation_bundle', lambda *a,**k: pytest.fail('disabled exported'))
        assert adapter.run_app_fixed_sd_polish(fitted) is None
        assert adapter.run_app_fixed_sd_polish(None, enabled=False) is None
    assert input_bytes(fitted) == before
    monkeypatch.setattr(adapter, 'polish_fixed_sd', lambda *a,**k: pytest.fail('unsupported input reached optimizer'))
    for key,value in [('model','RSM'),('method','JMLE'),('estimate_population_sd',True),
                      ('facet_regularization_enabled',True),('dummy_facets',['Task']),
                      ('quad_points',15),('quad_points',31.5),('quad_points',True)]:
        changed = deepcopy(fitted); changed['config'][key] = value
        with pytest.raises(ValueError): adapter.run_app_fixed_sd_polish(changed, enabled=True)
    for change in ['anchors','weights','scores','reported_nll','level_order','nested_penalty']:
        changed = deepcopy(fitted)
        if change=='anchors': changed['config']['facet_specs']['Task']['anchors'][0] = .3
        elif change=='weights': changed['prep']['data'].iloc[0, changed['prep']['data'].columns.get_loc('Weight')] = .5
        elif change=='scores': changed['prep']['data'].iloc[0, changed['prep']['data'].columns.get_loc('score_k')] = 3
        elif change=='reported_nll': changed['opt'].fun += .1
        elif change=='nested_penalty': changed['config']['facet_regularization']['enabled'] = True
        else: changed['prep']['data']['Rater'] = changed['prep']['data']['Rater'].cat.reorder_categories(['R1','R0'])
        with pytest.raises(ValueError): adapter.run_app_fixed_sd_polish(changed, enabled=True)
    with pytest.raises(ValueError): adapter.run_app_fixed_sd_polish(fitted, enabled='yes')
    assert input_bytes(fitted) == before


@pytest.mark.parametrize('sigma,missing', [(1.1,False),(.7,True),(1.6,False)])
def test_candidate_preserves_fixed_rule_labels_and_original_and_matches_independent_gh(sigma, missing):
    fit = fit_example(sigma, missing)
    before = input_bytes(fit)
    run = adapter.run_app_fixed_sd_polish(fit, enabled=True)
    assert input_bytes(fit) == before
    assert run['completed'] and run['gradient_tolerance_met'] and run['failure'] is None
    assert run['inference_ready'] is False and run['quadrature_sensitivity_pass'] is None
    assert run['candidate_is_fitted_result'] is False
    quad = app.make_mml_quadrature(fit['config'])
    np.testing.assert_array_equal(run['fixed_quadrature']['nodes'], quad['nodes'])
    np.testing.assert_array_equal(run['fixed_quadrature']['weights'], quad['weights'])
    assert run['fixed_quadrature']['sigma'] == sigma
    assert len(run['candidate_coordinates']) == 8
    assert run['constraint_residual'] < 1e-12
    assert all(i['positive_definite'] for i in run['information'])
    assert run['gradient_audit']['maximum_absolute_difference'] < 1e-6
    assert run['stages'][-1]['nll'] <= run['original_optimizer']['reevaluated_nll']+1e-9
    # Independent category design and SciPy GH, including positive rater signs and missing cells.
    candidate = deepcopy(fit)
    candidate['opt'].x = np.array(run['candidate_coordinates'])
    candidate['params'] = run['candidate_parameters']
    bundle = app.build_cross_engine_validation_bundle(candidate)
    job = native.prepare_inputs(bundle['data.csv'], bundle['analysis_snapshot.json'])
    z, w = roots_hermitenorm(121); theta = sigma*z; w /= np.sqrt(2*np.pi)
    logits = (native.native_design() @ np.asarray(job['xsi']))[:,None,:]+theta[None,:,None]*np.arange(4)
    logp = logits-logsumexp(logits,axis=-1,keepdims=True)
    y = np.asarray(job['responses'],dtype=float); ll = np.zeros((len(y),len(theta)))
    for j in range(8):
        ok = np.isfinite(y[:,j]); ll[ok] += logp[j,:,y[ok,j].astype(int)]
    lm = logsumexp(ll+np.log(w),axis=1); posterior = np.exp(ll+np.log(w)-lm[:,None])
    mean = posterior@theta
    sd = np.sqrt(np.sum(posterior*(theta[None,:]-mean[:,None])**2,axis=1))
    assert abs(-lm.sum()-run['stages'][-1]['nll']) < 1e-9
    np.testing.assert_allclose(mean,run['candidate_person_scores']['Estimate'],atol=1e-10,rtol=0)
    np.testing.assert_allclose(sd,run['candidate_person_scores']['SD'],atol=1e-10,rtol=0)
    job.update(nodes=801, bound=30.)
    reference = native.evaluate(job, job['xsi'])
    if sigma == 1.6:
        # Retained counterexample: finite-GH stationarity does not confer integral accuracy.
        assert abs(reference['raw_nll']-run['stages'][-1]['nll']) > 1e-3
        assert run['quadrature_sensitivity_pass'] is None and not run['inference_ready']
    assert job['person_labels'][0]=='000'
    assert run['candidate_person_scores']['Person'].tolist() == job['person_labels']
    assert sum(v is None for row in job['responses'] for v in row) == int(missing)


def test_kernel_handles_objective_offset_and_preserves_failure_and_large_gradient(monkeypatch):
    target = np.array([1.3, -.7])
    vg = lambda p:(1e8+.5*np.sum((p-target)**2), p-target)
    factory = lambda a:lambda p:float((a-target)@(p-a)+.5*np.sum((p-a)**2))
    run = kernel.polish_fixed_sd(np.zeros(2), vg, factory)
    assert run['completed'] and run['gradient_tolerance_met']
    np.testing.assert_allclose(run['candidate_coordinates'], target, atol=1e-9, rtol=0)
    bad = kernel.polish_fixed_sd(np.zeros(2), vg, lambda a:lambda p:1.)
    assert not bad['completed'] and bad['candidate_coordinates'] is None
    assert bad['stages'][0]['returned_coordinates'] is not None
    assert bad['failure']['stage']=='difference'
    from scipy.optimize import OptimizeResult
    def false_success(fun, x, **kwargs):
        value, gradient = fun(x)
        return OptimizeResult(x=x,fun=value,jac=gradient,success=True,status=0,message='stopped',nit=0,nfev=1)
    monkeypatch.setattr(kernel, 'minimize', false_success)
    stopped = kernel.polish_fixed_sd(np.zeros(2), vg, factory)
    assert stopped['completed'] and not stopped['gradient_tolerance_met']
    assert stopped['inference_ready'] is False
    broken = kernel.polish_fixed_sd(np.zeros(2), lambda p:(np.nan,np.ones(2)), factory)
    assert not broken['completed'] and broken['failure']['stage']=='raw'
