"""Independent finite-GH information, including structural/log-SD cross terms."""

import numpy as np
import pytest

import streamlit_app as app


@pytest.mark.parametrize("model", ["RSM", "PCM", "GPCM"])
@pytest.mark.parametrize("sign", [-1, 1])
@pytest.mark.parametrize("sigma", [.7, 1.4])
@pytest.mark.parametrize("weighted", [False, True])
def test_joint_information_matches_direct_finite_sum(model, sign, sigma, weighted):
    config = {"model": model, "method": "MML", "n_cat": 3, "n_person": 3,
              "facet_names": ["Criterion"], "facet_levels": {"Criterion": ["C1", "C2"]},
              "facet_specs": {"Criterion": app.build_facet_constraint(["C1", "C2"], centered=True)},
              "facet_signs": {"Criterion": sign}, "step_facet": "Criterion",
              "slope_facet": "Criterion",
              "estimate_population_sd": True, "population_prior_sd": sigma}
    sizes = app.build_param_sizes(config)
    structural = np.array([.35, -.4] if model == "RSM" else [.35, -.4, .2])
    if model == "GPCM":
        structural = np.r_[structural, np.log(2.)]
    point = np.r_[structural, np.log(sigma)]
    idx = {"person": np.repeat(np.arange(3), 2), "score_k": np.array([0, 0, 1, 2, 2, 1]),
           "facets": {"Criterion": np.tile([0, 1], 3)},
           "step_idx": np.tile([0, 1], 3), "slope_idx": np.tile([0, 1], 3),
           "weight": np.array([.5, 2., 1., 1.5, 2., .75]) if weighted else np.ones(6)}
    quad = app.gauss_hermite_normal(19, sd=sigma)
    theta, prior = quad["nodes"], quad["weights"]
    category = np.arange(3)
    direction = sign * np.array([1., -1.])
    steps = np.array([[0., -.4, 0.], [0., -.4 if model == "RSM" else .2, 0.]])
    slopes = np.array([2., .5]) if model == "GPCM" else np.ones(2)
    logits = slopes[None, :, None]*((theta[:, None] + .35*direction)[:, :, None]*category - steps)
    probs = np.exp(logits - logits.max(axis=2, keepdims=True))
    probs /= probs.sum(axis=2, keepdims=True)
    # Derivatives of the explicit category logits, independent of app gradients.
    derivative = np.zeros((*probs.shape, len(point)))
    derivative[:, :, :, 0] = (slopes*direction)[:, None]*category
    for item in range(2):
        derivative[:, item, :, 1 if model == "RSM" else 1+item] = -slopes[item]*(category == 1)
    derivative[:, :, :, -1] = theta[:, None, None]*slopes[None, :, None]*category
    second = np.zeros((*probs.shape, len(point), len(point)))
    second[..., -1, -1] = derivative[..., -1]
    if model == "GPCM":
        slope_sign = np.array([1., -1.])[None, :, None]
        derivative[..., -2] = slope_sign*logits
        second[..., -2, -2] = logits  # second derivative of exp(+/-lambda)
        for j in [0, 1, 2, len(point)-1]:
            second[..., j, -2] = slope_sign*derivative[..., j]
            second[..., -2, j] = second[..., j, -2]
    mean_derivative = np.einsum("qik,qikp->qip", probs, derivative)
    expected = probs @ category
    joint_information = np.zeros((len(point), len(point)))
    joint_gradient = np.zeros(len(point))
    for person in range(3):
        y = idx["score_k"][2*person:2*person+2]
        w = idx["weight"][2*person:2*person+2]
        mass = prior*np.prod(probs[:, np.arange(2), y]**w, axis=1)
        posterior = mass/mass.sum()
        score = np.sum(w[None, :, None]*(mean_derivative - derivative[:, np.arange(2), y]), axis=1)
        centered = derivative - mean_derivative[:, :, None, :]
        conditional_information = np.einsum("i,qik,qikp,qikr->qpr", w, probs, centered, centered)
        # Nonlinear category logits contribute E[d² logit]-d² logit(y),
        # including slope/location, slope/step, slope/log-SD cross terms.
        conditional_information += np.einsum("i,qik,qikpr->qpr", w, probs, second)
        conditional_information -= np.einsum("i,qipr->qpr", w, second[:, np.arange(2), y])
        mean_score = posterior @ score
        joint_gradient += mean_score
        score_residual = score - mean_score
        joint_information += np.einsum("q,qpr->pr", posterior, conditional_information)
        joint_information -= np.einsum("q,qp,qr->pr", posterior, score_residual, score_residual)

    _, gradient = app.mfrm_loglik_mml_value_grad(
        structural, idx, config, sizes, quad, include_log_sigma=True)
    actual = app._compute_mml_observed_information(
        point, idx, config, sizes, quad, include_log_sigma=True)
    np.testing.assert_allclose(gradient, joint_gradient, atol=2e-13)
    np.testing.assert_allclose(actual, joint_information, atol=2e-8)
    assert np.max(np.abs(joint_information[:-1, -1])) > 1e-3
    fixed = app._compute_mml_observed_information(structural, idx, config, sizes, quad)
    np.testing.assert_allclose(actual[:-1, :-1], fixed, atol=2e-9)
    # A different FD step must preserve the full matrix, not only its diagonal.
    refined = app._compute_mml_observed_information(
        point, idx, config, sizes, quad, rel_step=5e-6, include_log_sigma=True)
    np.testing.assert_allclose(refined, joint_information, atol=2e-8)


def test_gpcm_last_slope_generates_cross_terms_between_both_free_log_slopes():
    # Three slopes: exp(l1), exp(l2), exp(-l1-l2). The last discrimination
    # contributes a nonzero mixed second derivative in l1/l2.
    config={'model':'GPCM','method':'MML','n_cat':4,'n_person':2,
        'facet_names':['Criterion'],'facet_levels':{'Criterion':['C1','C2','C3']},
        'facet_specs':{'Criterion':app.build_facet_constraint(['C1','C2','C3'],centered=False)},
        'facet_signs':{'Criterion':-1},'step_facet':'Criterion','slope_facet':'Criterion'}
    sizes=app.build_param_sizes(config)
    point=np.r_[[.2,-.1,.4,-.8,-.1,-1.,.2,-.5,-.2],np.log([2.,.8,1.1])]
    idx={'person':np.repeat([0,1],3),'score_k':np.array([0,2,1,3,1,2]),
        'facets':{'Criterion':np.tile(np.arange(3),2)},'step_idx':np.tile(np.arange(3),2),
        'slope_idx':np.tile(np.arange(3),2),'weight':np.array([1.,2.,.5,.5,1.5,1.])}
    quad=app.gauss_hermite_normal(17,sd=1.1);k=np.arange(4)
    design=np.zeros((3,4,9))
    for item in range(3):
        design[item,:,item]=-k
        for step in range(2):design[item,:,3+2*item+step]=-((k>step)&(k<3)).astype(float)
    transform=np.array([[1.,0.],[0.,1.],[-1.,-1.]])
    slopes=np.exp(transform@point[9:11]); scores=[];conditional=[];ll=[]
    for theta in quad['nodes']:
        logits=slopes[:,None]*(design@point[:9]+theta*k)
        probs=np.exp(logits-logits.max(axis=1,keepdims=True));probs/=probs.sum(axis=1,keepdims=True)
        d=np.zeros((3,4,12));d[:,:,:9]=slopes[:,None,None]*design
        d[:,:,9:11]=logits[:,:,None]*transform[:,None,:]
        d[:,:,-1]=slopes[:,None]*theta*k
        second=np.zeros((3,4,12,12));second[:,:,-1,-1]=d[:,:,-1]
        for j in range(2):
            column=9+j
            for linear in list(range(9))+[11]:
                second[:,:,linear,column]=transform[:,j,None]*d[:,:,linear]
                second[:,:,column,linear]=second[:,:,linear,column]
            for other in range(2):second[:,:,column,9+other]=logits*(transform[:,j]*transform[:,other])[:,None]
        mean=np.einsum('ik,ikd->id',probs,d);center=d-mean[:,None,:]
        qs=[];qh=[];ql=[]
        for person in range(2):
            y=idx['score_k'][3*person:3*person+3];w=idx['weight'][3*person:3*person+3]
            qs.append(np.sum(w[:,None]*(mean-d[np.arange(3),y]),axis=0))
            qh.append(np.einsum('i,ik,ikd,ike->de',w,probs,center,center)
                +np.einsum('i,ik,ikde->de',w,probs,second)-np.einsum('i,ide->de',w,second[np.arange(3),y]))
            ql.append(np.sum(w*np.log(probs[np.arange(3),y])))
        scores.append(qs);conditional.append(qh);ll.append(ql)
    mass=np.exp(np.array(ll).T)*quad['weights'];posterior=mass/mass.sum(axis=1,keepdims=True)
    scores=np.array(scores).transpose(1,0,2);mean=np.einsum('pq,pqd->pd',posterior,scores)
    centered=scores-mean[:,None,:]
    information=np.einsum('pq,qpde->de',posterior,np.array(conditional))-np.einsum('pq,pqd,pqe->de',posterior,centered,centered)
    _,gradient=app.mfrm_loglik_mml_value_grad(point[:-1],idx,config,sizes,quad,include_log_sigma=True)
    actual=app._compute_mml_observed_information(point,idx,config,sizes,quad,include_log_sigma=True)
    np.testing.assert_allclose(gradient,mean.sum(axis=0),atol=2e-12)
    np.testing.assert_allclose(actual,information,atol=5e-8)
    assert abs(information[9,10]) > 1e-3
    assert np.max(abs(information[9:11,-1])) > 1e-3
