"""Counts expose selection by unavailable comparisons without fabricating fits."""
from copy import deepcopy

import pytest

from validation.mml_population_analysis import binomial_interval,equivalence_status
from validation.mml_population_missingness_review import review_comparison


def evidence(observed,rejected,count=2000):
    scenario=dict(id='DEVELOPMENT_ONLY_MISSINGNESS_COUNTS',replicates=count,
        paired_rsm_pcm_comparison=True,design=dict(steps=[[0.,0.,0.]]))
    comparison=dict(planned=count,available=observed,unavailable=count-observed,rejected=rejected,
        conditional_rate=rejected/observed if observed else None,all_planned_rate=rejected/count,
        all_planned_mc_interval=binomial_interval(rejected,count),
        null_acceptance=equivalence_status(binomial_interval(rejected,observed),[.035,.065]))
    return dict(status='complete',scenario_id=scenario['id'],planned=count,
        numerically_eligible=count,rsm_pcm_comparison=comparison),scenario


def test_conditional_pass_does_not_hide_missing_comparisons_with_successful_primary_fits():
    summary,scenario=evidence(1600,80);before=deepcopy(summary)
    result=review_comparison(summary,scenario);rate=result['conditional_and_completion_evidence']
    assert result['registered_conditional_acceptance']=='pass'
    assert rate['conditional_rate']==.05
    assert rate['all_planned_identification_interval']==[.04,.24]
    assert rate['all_planned_confidence_envelope'][1]>.24
    assert result['delivered_rejections_per_planned']['estimate']==.04
    assert result['comparison_unavailability']['rate']==.2
    assert not result['all_planned_false_positive_or_power_point_identified']
    assert not result['scientific_inference_ready'] and summary==before


def test_all_missing_is_unknown_not_a_zero_false_positive_rate():
    summary,scenario=evidence(0,0);result=review_comparison(summary,scenario)
    rates=result['conditional_and_completion_evidence']
    assert result['registered_conditional_acceptance']=='unavailable'
    assert rates['conditional_rate'] is None and rates['conditional_mc_interval'] is None
    assert rates['all_planned_identification_interval']==[0.,1.]
    assert rates['all_planned_confidence_envelope']==[0.,1.]
    assert result['delivered_rejections_per_planned']['estimate']==0.


def test_complete_decisions_recover_the_same_rate_and_exact_interval():
    summary,scenario=evidence(2000,100);result=review_comparison(summary,scenario)
    rates=result['conditional_and_completion_evidence']
    assert result['all_planned_false_positive_or_power_point_identified']
    assert rates['all_planned_identification_interval']==[.05,.05]
    assert rates['all_planned_confidence_envelope']==rates['conditional_mc_interval']


@pytest.mark.parametrize('change',[
    dict(unavailable=0),dict(rejected=1601),dict(planned=1999),dict(available=True),
    dict(conditional_rate=.04),dict(all_planned_rate=.05),dict(all_planned_rate=float('nan')),
    dict(all_planned_mc_interval=[0.,1.]),
])
def test_count_or_denominator_mismatch_is_rejected(change):
    summary,scenario=evidence(1600,80);summary['rsm_pcm_comparison'].update(change)
    with pytest.raises(ValueError):review_comparison(summary,scenario)


def test_no_test_rate_cannot_be_filled_with_zero():
    summary,scenario=evidence(0,0);summary['rsm_pcm_comparison']['conditional_rate']=0.
    with pytest.raises(ValueError,match='No observed tests'):review_comparison(summary,scenario)


def test_only_registered_complete_nonequivalent_model_comparisons_are_reviewed():
    summary,scenario=evidence(2000,100)
    incomplete=dict(summary,status='incomplete')
    with pytest.raises(ValueError,match='complete'):review_comparison(incomplete,scenario)
    with pytest.raises(ValueError,match='complete'):review_comparison(summary,dict(scenario,id='other'))
    with pytest.raises(ValueError,match='Binary'):review_comparison(summary,dict(scenario,design=dict(steps=[[0.]])))
    unpaired=dict(scenario,paired_rsm_pcm_comparison=False)
    with pytest.raises(ValueError,match='Unregistered'):review_comparison(summary,unpaired)
    assert review_comparison(dict(summary,rsm_pcm_comparison=None),unpaired)['status']=='not_applicable'
