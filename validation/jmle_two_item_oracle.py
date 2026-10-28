"""Exact two-item Rasch illustration of incidental-parameter bias.

P(Y_j=1 | theta)=logistic(theta-b_j), b=(-delta/2,+delta/2).
Condition on M persons with total score 1. K counts pattern 10 and follows
Binomial(M, logistic(delta)) regardless of the Person distribution. Unbounded
joint ML profiles each Person parameter; conditional ML conditions it out.
This is an analytic model, not a qualification run of the bounded app engine.
"""
import numpy as np
from scipy.special import expit,ndtr,ndtri
from scipy.stats import binom


def estimates(n10,n01):
    if any(type(n) is not int or n<0 for n in [n10,n01]) or n10+n01==0:
        raise ValueError('Nonnegative integer pattern counts with positive total required')
    if min(n10,n01)==0:
        return dict(status='no_finite_mle',JMLE=None,CML=None)
    contrast=float(np.log(n10)-np.log(n01));variance=1/n10+1/n01
    return dict(status='interior',
        JMLE=dict(estimate=2*contrast,profile_se=float(np.sqrt(2*variance))),
        CML=dict(estimate=contrast,profile_se=float(np.sqrt(variance))))


def joint_nll(point,n10,n01):
    """Row likelihood with separately profiled abilities for both patterns.

    Persons with the same response pattern share the same unique ability
    maximizer for fixed delta; grouping does not remove that nuisance fit.
    """
    theta=np.asarray(point[:2]);delta=point[2]
    eta=theta[:,None]+np.array([delta/2,-delta/2])
    y=np.array([[1.,0.],[0.,1.]])
    return float(np.dot([n10,n01],np.sum(np.logaddexp(0,eta)-y*eta,axis=1)))


def exact_distribution(informative_persons,delta):
    m=informative_persons
    if type(m) is not int or m<2 or not np.isfinite(delta):raise ValueError('M>=2 and finite delta required')
    p=float(expit(delta))
    if not 0<p<1:raise ValueError('Conditional probability not representable in the open unit interval')
    k=np.arange(m+1);mass=binom.pmf(k,m,p)
    if abs(float(mass.sum())-1)>1e-11:raise ArithmeticError('Enumeration probability mass check failed')
    eligible=(k>0)&(k<m);prob=mass[eligible];available=float(prob.sum());z=ndtri(.975)
    counts=k[eligible];log_odds=np.log(counts)-np.log(m-counts)
    base_variance=1/counts+1/(m-counts);results={}
    for method,multiplier in [('CML',1.),('JMLE',2.)]:
        value=multiplier*log_odds;se=np.sqrt(multiplier*base_variance)
        mean=float(np.dot(prob,value)/available) if available else None
        variance=float(np.dot(prob,(value-mean)**2)/available) if available else None
        rms_se=float(np.sqrt(np.dot(prob,se**2)/available)) if available else None
        covered=float(np.dot(prob,abs(value-delta)<=z*se))
        results[method]=dict(mean_given_finite=mean,bias_given_finite=mean-delta if available else None,
            rmse_given_finite=float(np.sqrt(np.dot(prob,(value-delta)**2)/available)) if available else None,
            sampling_sd_given_finite=float(np.sqrt(variance)) if available else None,
            rms_profile_se_given_finite=rms_se,
            rms_se_over_sampling_sd=rms_se/np.sqrt(variance) if variance is not None and variance>0 else None,
            coverage_given_interval=covered/available if available else None,
            delivered_coverage=covered,
            possible_coverage_completion_bounds=[covered,min(1.,covered+float(mass[~eligible].sum()))])
    return dict(informative_persons=m,true_contrast=float(delta),conditional_pattern10_probability=p,
        probability_mass=float(mass.sum()),finite_estimate_probability=available,
        no_finite_mle_probability=float(mass[~eligible].sum()),methods=results,
        log_no_finite_mle_probability=float(np.logaddexp(m*np.log(p),m*np.log1p(-p))),
        classification='EXACT_BINOMIAL_ENUMERATION_CONDITIONAL_ON_INFORMATIVE_PERSON_COUNT',
        Monte_Carlo_replicates=0,scientific_inference_ready=False,
        interpretation='Complete binomial enumeration evaluated in floating point, with mass check and log failure probability retained. No simulation draws. Bias/SD condition on a finite estimate; delivered coverage retains no-finite-MLE outcomes. Unbounded two-item reference only, not the application optimizer or general MFRM performance.')


def limit_values(delta):
    return dict(CML_contrast=float(delta),JMLE_contrast=2*float(delta),
        JMLE_profile_SE_over_sampling_SD=float(1/np.sqrt(2)),
        JMLE_Wald_coverage_at_zero=float(2*ndtr(ndtri(.975)/np.sqrt(2))-1))
