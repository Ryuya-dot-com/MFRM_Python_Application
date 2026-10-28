"""Expanded slope constraints must survive joint-SD development construction."""
import numpy as np
import pytest

import streamlit_app as app
from validation.mml_gpcm_joint_check import prepare_problem, joint_constraints, expanded, expansion_jacobian, generate
from validation.mml_free_sd_stationarity_adapter import prepare_app_free_sd_problem_from_data


def test_gpcm_development_problem_keeps_the_dependent_slope_bounded():
    data=generate(2026100439,'sparse',12)
    _,problem=prepare_problem(data,9)
    point=np.r_[problem.structural_start,0.]
    slope=app._build_param_slices(problem.sizes)['log_slopes']
    point[slope]=2.5  # Every free slope is inside; last expanded slope is -5.
    constraint=joint_constraints(problem)[0]
    value=constraint.A@point
    assert np.any(value>constraint.ub)
    assert np.isinf(problem.constraint_residual(point[:-1]))
    point[-1]=np.log(2.)
    np.testing.assert_array_equal(constraint.A@point,value)  # SD is not a slope coordinate.
    point[slope]=[.3,-.2]
    assert problem.constraint_residual(point[:-1]) < 1e-12
    assert np.all(constraint.A@point <= constraint.ub)
    assert np.all(constraint.A@point >= constraint.lb)


def test_gpcm_expansion_jacobian_includes_last_slope_and_exponential_maps():
    _,problem=prepare_problem(generate(2026100438,'unbalanced',12),9)
    point=np.r_[np.linspace(-.2,.3,sum(problem.sizes.values())),np.log(1.3)]
    labels,jac=expansion_jacobian(problem.config,point)
    columns=[]
    for j in range(len(point)):
        d=np.zeros(len(point));d[j]=1e-6
        columns.append((expanded(problem,point+d)-expanded(problem,point-d))/2e-6)
    np.testing.assert_allclose(jac,np.column_stack(columns),atol=3e-10,rtol=1e-9)
    # Full covariance is needed: transformed raw-slope uncertainty uses exp.
    rng=np.random.default_rng(381);a=rng.normal(size=(len(point),len(point)));cov=a@a.T
    variance=np.einsum('ip,pq,iq->i',jac,cov,jac)
    for j in range(3):
        lam=labels.index(f'log_slope:{j}');raw=labels.index(f'slope:{j}')
        assert variance[raw] == pytest.approx(expanded(problem,point)[raw]**2*variance[lam])
    assert variance[-1] == pytest.approx(1.3**2*cov[-1,-1])


def test_rsm_pcm_v2_input_contract_is_not_silently_broadened():
    data=generate(2026100437,'regular',8)
    with pytest.raises(ValueError,match='RSM and PCM only'):
        prepare_app_free_sd_problem_from_data(data,person_col='Person',facet_cols=['Rater','Task','Criterion'],
            score_col='Score',rating_min=0,rating_max=3,model='GPCM',noncenter_facet='Criterion',step_facet='Criterion')
