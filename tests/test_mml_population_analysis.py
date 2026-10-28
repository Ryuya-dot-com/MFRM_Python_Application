"""Synthetic ledgers verify endpoint targets and incomplete/failure handling."""
import numpy as np
import pytest

from tests.test_mml_population_study import fixture_plan
from validation.mml_population_analysis import (binomial_interval, equivalence_status,
    scenario_calibration_truth, summarize_scenario)
from validation.mml_population_study import AttemptLedger


def test_binomial_boundaries_do_not_turn_zero_mcse_into_certainty():
    assert binomial_interval(0,10)[0]==0 and binomial_interval(0,10)[1]>.3
    assert binomial_interval(10,10)[1]==1 and binomial_interval(10,10)[0]<.7
    assert binomial_interval(0,0) is None
    with pytest.raises(ValueError):binomial_interval(True,10)
    assert equivalence_status([.94,.96],[.925,.975])=='pass'
    assert equivalence_status([.90,.94],[.925,.975])=='inconclusive'
    assert equivalence_status([.89,.91],[.925,.975])=='fail'
    assert equivalence_status(None,[.925,.975])=='unavailable'


def test_incomplete_scenario_does_not_release_interim_endpoints(tmp_path):
    plan,scenario,_=fixture_plan();ledger=AttemptLedger(tmp_path/'a.sqlite',plan_sha256='a'*64)
    result=summarize_scenario(ledger,scenario,plan)
    assert result['status']=='incomplete'
    assert 'calibration' not in result and result['terminal']==0
    ledger.close()


def test_native_criterion_location_and_failure_denominators_survive_analysis(tmp_path):
    plan,scenario,_=fixture_plan();scenario['replicates']=3
    ledger=AttemptLedger(tmp_path/'a.sqlite',plan_sha256='a'*64)
    truths=scenario_calibration_truth(scenario);theta=np.linspace(-1.,1.,scenario['design']['persons'])
    for replicate in range(1,4):
        run_id=f"{scenario['id']}__{replicate:04d}";eligible=replicate<3
        ledger.append(run_id,'input_retained',dict(truth=dict(Person=theta)))
        rows=[]
        for truth in truths:
            comparable=bool(np.any(truth['jacobian']!=0));error=.1*replicate if truth['kind']=='Criterion' else 0.
            rows.append(dict(kind=truth['kind'],facet=truth['facet'],level=truth['level'],truth=truth['value'],
                estimate=truth['value']+error,eligible=eligible and comparable,se=.2 if eligible and comparable else None,
                covered=True if eligible and comparable else None))
        endpoints=dict(calibration=rows,person=dict(truth=theta,eligible=eligible,fitted_eap=theta+.1,
            fitted_posterior_sd=np.full(len(theta),.2),known_calibration_numerically_eligible=True,
            known_calibration_eap=.8*theta,known_calibration_posterior_sd=np.full(len(theta),.3)))
        ledger.append(run_id,'PCM_assessment',dict(status='eligible_candidate' if eligible else 'ineligible_candidate',endpoints=endpoints))
        ledger.append(run_id,'terminal',dict(status='completed',comparison=dict(available=eligible,p_value=.01 if eligible else None,
            information_criteria=dict(RSM=dict(AIC=101.,BIC=110.),PCM=dict(AIC=100.,BIC=111.)))))
    result=summarize_scenario(ledger,scenario,plan)
    assert result['planned']==3 and result['failed_or_ineligible']==1
    for row in result['calibration']:
        if row['kind']=='Criterion':
            assert row['metrics']['bias']['estimate']==pytest.approx(.15)
            assert row['metrics']['coverage_given_interval']['estimate']==1.
            assert row['metrics']['coverage_per_known_target']['estimate']==pytest.approx(2/3)
    fitted=next(r for r in result['person'] if r['calibration_basis']=='fitted' and r['truth_quartile']=='all')
    known=next(r for r in result['person'] if r['calibration_basis']=='known' and r['truth_quartile']=='all')
    assert fitted['metrics']['interval_availability']['estimate']==pytest.approx(2/3)
    assert known['metrics']['interval_availability']['estimate']==1.
    assert result['rsm_pcm_comparison']['available']==2 and result['rsm_pcm_comparison']['unavailable']==1
    assert result['rsm_pcm_comparison']['information_criteria']['AIC']['PCM_lower']==2
    assert result['rsm_pcm_comparison']['information_criteria']['BIC']['PCM_lower']==0
    assert result['rsm_pcm_comparison']['all_planned_mcse']==pytest.approx(np.sqrt((2/3)*(1/3)/3))
    assert not result['scientific_inference_ready']
    ledger.close()
    readonly=AttemptLedger.read_only(tmp_path/'a.sqlite',plan_sha256='a'*64)
    assert summarize_scenario(readonly,scenario,plan)['planned']==3
    readonly.close()


def test_generation_failure_leaves_fixed_calibration_truth_known_and_person_truth_unknown(tmp_path):
    plan,scenario,run=fixture_plan();ledger=AttemptLedger(tmp_path/'a.sqlite',plan_sha256='a'*64)
    ledger.append(run['run_id'],'terminal',dict(status='generation_failed',truth_available=False))
    result=summarize_scenario(ledger,scenario,plan)
    assert result['truth_unknown_replicates']==1
    for row in result['calibration']:
        if row['comparable']:
            assert row['metrics']['all_planned_truth_known']
            assert row['metrics']['coverage_per_known_target']['estimate']==0.
            assert row['metrics']['coverage_given_interval']['estimate'] is None
    assert all(not r['metrics']['all_planned_truth_known'] for r in result['person'])
    ledger.close()
