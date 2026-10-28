"""Normal-moment identities and actual MML consumers at high GH orders."""
import math

import numpy as np
import pytest
from numpy.polynomial.hermite import hermgauss

import streamlit_app as app
from tests.test_mml_grouped_likelihood import problem


def test_normal_rules_integrate_known_moments_and_preserve_supported_low_orders():
    for q in [1,2,3,15,31,81,161,321,641,1281]:
        for sd in [.2,1.,2.5]:
            quad=app.gauss_hermite_normal(q,sd)
            z,w=quad['nodes'],quad['weights']
            assert np.isfinite(z).all() and np.isfinite(w).all() and (w>=0).all()
            for power in range(0,min(2*q,13),2):
                exact=math.prod(range(1,power,2))*sd**power
                assert w @ z**power == pytest.approx(exact,rel=5e-13,abs=1e-14)
            assert abs(w@z)<1e-14
            if q<=321:
                old_z,old_w=hermgauss(q)
                np.testing.assert_allclose(z,np.sqrt(2.)*old_z*sd,atol=1e-12,rtol=1e-13)
                np.testing.assert_allclose(w,old_w/np.sqrt(np.pi),atol=2e-15,rtol=1e-12)
    assert (app.gauss_hermite_normal(1281)['weights']==0).any()


def test_reused_standard_rule_keeps_independent_outputs_and_validates_every_order():
    app._standard_normal_quadrature.cache_clear()
    before=app.gauss_hermite_normal(31,1.2)
    expected_z=before['nodes'].copy();expected_w=before['weights'].copy()
    before['nodes'][:]=99.;before['weights'][:]=0.
    again=app.gauss_hermite_normal(np.int64(31),1.2)
    np.testing.assert_array_equal(again['nodes'],expected_z)
    np.testing.assert_array_equal(again['weights'],expected_w)
    assert app._standard_normal_quadrature.cache_info().hits==1
    for q in [0,-1,True,np.bool_(True),31.,2.5,'31',np.nan]:
        with pytest.raises(ValueError):app.gauss_hermite_normal(q)
        with pytest.raises(ValueError):app.make_mml_quadrature({'quad_points':q})
        with pytest.raises(ValueError):app.make_mml_quadrature({},q)
    for sd in [0,-1,np.nan,np.inf,1e308]:
        with pytest.raises(ValueError):app.gauss_hermite_normal(31,sd)
    assert len(app.make_mml_quadrature({})['nodes'])==15
    assert len(app.make_mml_quadrature({'quad_points':31},None)['nodes'])==31


@pytest.mark.parametrize('model',['RSM','PCM','GPCM'])
@pytest.mark.parametrize('q',[641,1281])
def test_high_order_likelihood_gradient_and_person_outputs_match_grouped(model,q):
    _,config,sizes,idx,grouped=problem(model,'unbalanced',duplicate=True)
    sigma=1.2;config['population_prior_sd']=sigma
    point=np.r_[np.linspace(-.4,.5,sum(sizes.values())),np.log(sigma)]
    with np.errstate(divide='raise',invalid='raise',over='raise',under='ignore'):
        quad=app.gauss_hermite_normal(q,sigma)
        expected=grouped.evaluate(point,quad['nodes']/sigma,quad['weights'],information=True)
        params=app.expand_params(point[:-1],sizes,config)
        nll=app.mfrm_loglik_mml(point[:-1],idx,config,sizes,quad)
        grad_nll,gradient=app.mfrm_loglik_mml_value_grad(point[:-1],idx,config,sizes,quad,include_log_sigma=True)
        posterior,ll=app._e_step_posteriors(idx,config,params,quad)
        eap=app.compute_person_eap(idx,config,params,quad)
        scores=app.compute_person_posterior_outputs(idx,config,params,quad,n_plausible_values=3,seed=7241)
        positive=quad['weights']>0
        reduced=dict(quad,nodes=quad['nodes'][positive],weights=quad['weights'][positive])
        reduced_post,reduced_ll=app._e_step_posteriors(idx,config,params,reduced)
        if model=='GPCM' and q==641:
            info=app._compute_mml_observed_information(point,idx,config,sizes,quad,include_log_sigma=True)
            np.testing.assert_allclose(info,expected['information'],atol=5e-7,rtol=0)
    for value in [nll,grad_nll,-ll,-reduced_ll]:
        assert value==pytest.approx(expected['nll'],abs=2e-10)
    np.testing.assert_allclose(gradient,expected['gradient'],atol=2e-10,rtol=0)
    np.testing.assert_allclose(posterior,expected['posterior'],atol=2e-14,rtol=0)
    np.testing.assert_array_equal(posterior[:,~positive],0.)
    np.testing.assert_allclose(posterior[:,positive],reduced_post,atol=2e-14,rtol=0)
    for frame,columns in [(eap,['Estimate','SD']),(scores['scores'],['PosteriorMean','PosteriorSD'])]:
        np.testing.assert_allclose(frame[columns[0]],expected['eap'],atol=1e-12,rtol=0)
        np.testing.assert_allclose(frame[columns[1]],expected['sd'],atol=1e-12,rtol=0)
    assert scores['available']
    assert np.isin(scores['plausible_values'].iloc[:,1:].to_numpy(),quad['nodes'][positive]).all()


def test_invalid_prior_mass_is_rejected_before_all_posterior_consumers():
    _,config,sizes,idx,_=problem('RSM','regular')
    point=np.zeros(sum(sizes.values()));params=app.expand_params(point,sizes,config)
    bad_rules=[dict(nodes=[-1.,1.],weights=w) for w in [[.4,.4],[1.1,-.1],[0.,0.],[np.nan,.5],[1.]]]
    bad_rules.append(dict(nodes=[-1.,np.inf],weights=[.5,.5]))
    for quad in bad_rules:
        calls=[lambda:app.mfrm_loglik_mml(point,idx,config,sizes,quad),
               lambda:app.mfrm_loglik_mml_value_grad(point,idx,config,sizes,quad),
               lambda:app.compute_person_eap(idx,config,params,quad),
               lambda:app.compute_person_posterior_outputs(idx,config,params,quad)]
        for call in calls:
            with pytest.raises(ValueError,match='quadrature'):call()


def test_final_mml_person_table_reuses_posterior_scores(monkeypatch):
    data, _, _, _, _ = problem('RSM', 'regular')

    def duplicate_scoring(*args, **kwargs):
        raise AssertionError('Final scoring must not repeat the EAP likelihood pass')

    monkeypatch.setattr(app, 'compute_person_eap', duplicate_scoring)
    result = app.mfrm_estimate(data, person_col='Person',
        facet_cols=['Rater', 'Task', 'Criterion'], score_col='Score',
        rating_min=0, rating_max=3, keep_original=True, method='MML',
        noncenter_facet='Criterion', quad_points=15, maxit=30,
        compute_plausible_values=True, n_plausible_values=3, plausible_seed=726)
    people = result['facets']['person']
    scores = result['posterior']['scores']
    assert list(people.columns) == ['Person', 'Estimate', 'SD']
    np.testing.assert_array_equal(people.Person, scores.Person)
    np.testing.assert_array_equal(people.Estimate, scores.PosteriorMean)
    np.testing.assert_array_equal(people.SD, scores.PosteriorSD)
    assert len(result['posterior']['plausible_values']) == len(people)
