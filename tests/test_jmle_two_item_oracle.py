import numpy as np
import pytest
from scipy.special import expit
from scipy.optimize import minimize

from validation.jmle_two_item_oracle import estimates,exact_distribution,joint_nll,limit_values


@pytest.mark.parametrize('delta',[-2.,-.5,0.,.5,2.])
def test_conditional_pattern_probability_does_not_depend_on_person(delta):
    theta=np.array([-8.,-2.,0.,1.5,8.])
    p1=expit(theta+delta/2);p2=expit(theta-delta/2)
    p10=p1*(1-p2);p01=(1-p1)*p2
    np.testing.assert_allclose(p10/(p10+p01),expit(delta),rtol=0,atol=2e-13)


@pytest.mark.parametrize('counts',[(1,1),(3,1),(1,3),(40,24)])
def test_full_row_likelihood_has_the_joint_not_conditional_mle(counts):
    answer=estimates(*counts);point=[0.,0.,answer['JMLE']['estimate']]
    fit=minimize(joint_nll,[.3,-.2,0.],args=counts,method='BFGS',options={'gtol':1e-8})
    np.testing.assert_allclose(fit.x,point,rtol=0,atol=2e-6)
    assert abs(joint_nll(point,*counts)-fit.fun)<1e-10
    assert answer['JMLE']['estimate']==2*answer['CML']['estimate']


def test_separation_is_unavailable_and_kept_in_exact_denominator():
    assert estimates(0,8)['status']=='no_finite_mle'
    assert estimates(8,0)['status']=='no_finite_mle'
    result=exact_distribution(8,0.)
    assert result['no_finite_mle_probability']==pytest.approx(2/256)
    for method in result['methods'].values():
        assert method['delivered_coverage']==pytest.approx(method['coverage_given_interval']*254/256)
    assert result['Monte_Carlo_replicates']==0 and not result['scientific_inference_ready']


def test_incidental_bias_and_se_miscalibration_remain_with_many_persons():
    result=exact_distribution(16384,.5);jmle=result['methods']['JMLE'];cml=result['methods']['CML']
    assert abs(jmle['mean_given_finite']-1.)<.0002
    assert abs(cml['mean_given_finite']-.5)<.0001
    assert jmle['coverage_given_interval']<1e-20
    assert .94<cml['coverage_given_interval']<.96
    assert abs(jmle['rms_se_over_sampling_sd']-1/np.sqrt(2))<.001


def test_zero_contrast_does_not_rescue_the_nominal_joint_wald_interval():
    result=exact_distribution(16384,0.)['methods']['JMLE']
    assert abs(result['bias_given_finite'])<1e-12
    assert abs(result['coverage_given_interval']-limit_values(0.)['JMLE_Wald_coverage_at_zero'])<.005


@pytest.mark.parametrize('counts',[(-1,4),(0,0),(1.5,2),(True,4)])
def test_invalid_pattern_counts(counts):
    with pytest.raises(ValueError):estimates(*counts)
