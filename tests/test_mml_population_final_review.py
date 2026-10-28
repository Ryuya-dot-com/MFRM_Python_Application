"""Invented endpoint/reference records test the review, never qualify a model."""
from copy import deepcopy
import hashlib
import json

import numpy as np
import pytest

from validation import mml_population_final_review as review
from validation.mml_population_analysis import scenario_calibration_truth, summarize_scenario
from validation.mml_population_protocol import make_plan
from validation.mml_population_study import AttemptLedger, json_bytes


@pytest.fixture(autouse=True)
def synthetic_reference_arithmetic(monkeypatch):
    # This file tests inventories/assembly using invented passing flags. Actual
    # retained R artifacts and numerical rechecks have their own regression file.
    monkeypatch.setattr(review, '_recheck_r', lambda *args: True)


def plan_fixture():
    plan = make_plan()
    selected = [next(s for s in plan['scenarios'] if s['primary_model']==model and s['cohort']=='population_recovery')
                for model in ['RSM','PCM','GPCM']]
    selected += [next(s for s in plan['scenarios'] if s['cohort']=='negative_control')]
    plan['scenarios'] = deepcopy(selected); plan['negative_controls'] = 1
    return plan


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(json_bytes(value)+b'\n')


def reference_fixture(path, plan, scenario, *, unavailable_tam=False, tam_returncode=0):
    run, inventory = review.reference_inventory(scenario,1,plan)
    cases = [dict(id=key,lane=lane,available=not (unavailable_tam and lane=='TAM')) for key,lane in inventory]
    selection = dict(run_id=run,scenario_id=scenario['id'],study_sha256='DEVELOPMENT_ONLY',
                     required_cases=cases,missing_required_cases=sum(not c['available'] for c in cases))
    write(path/'reference_plan.json',selection)
    write(path/'manifest.json',{'reference_plan.json':hashlib.sha256((path/'reference_plan.json').read_bytes()).hexdigest()})
    output=path/'execution';output.mkdir()
    ledger=AttemptLedger(output/'events.sqlite',plan_sha256='DEVELOPMENT_ONLY')
    rows=[]
    for case in cases:
        key,lane=case['id'],case['lane']
        if not case['available']:
            rows.append(dict(id=key,lane=lane,status='required_input_or_point_unavailable',all_checks_pass=False));continue
        if lane=='TAM':
            files={}
            for name in [key+'_TAM_stdout.txt',key+'_TAM_stderr.txt',f'{key}_TAM/{key}_design.json'] + [f'{key}_TAM/{key}_{g}.json' for _,_,g in review.TAM_GRIDS]:
                write(output/name,dict(DEVELOPMENT_ONLY=True))
                files[name]=hashlib.sha256((output/name).read_bytes()).hexdigest()
            ledger.append(key,'TAM_fit_finished',dict(status='completed',value=dict(returncode=tam_returncode,files=files)))
        suffixes=['_'+g for _,_,g in review.TAM_GRIDS] if lane=='TAM' else ['']
        for suffix in suffixes:
            keys=(review.FINITE_CHECKS if lane=='finite_R' else review.TAM_CHECKS if lane=='TAM' else
                  review.CONTINUOUS_CHECKS | (set() if key.endswith('_truth') else review.INFORMATION_CHECKS | {'stationary','expanded_se'}))
            ledger.append(key+suffix,lane+'_review_finished',dict(status='completed',
                          value=dict(checks={k:True for k in keys},all_checks_pass=True)))
            rows.append(dict(id=key+suffix,lane=lane,status='completed',all_checks_pass=True))
    summary=dict(run_id=run,required_cases=len(cases),comparisons=rows,
                 all_checks_pass=all(r['all_checks_pass'] for r in rows),ledger=ledger.verify())
    ledger.close();write(output/'summary.json',summary)
    return summary


@pytest.mark.parametrize('model',['RSM','PCM','GPCM'])
def test_preselected_reference_inventory_and_read_only_success(tmp_path,model):
    plan=plan_fixture();scenario=next(s for s in plan['scenarios'] if s['primary_model']==model)
    reference_fixture(tmp_path,plan,scenario)
    before={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in tmp_path.rglob('*') if p.is_file()}
    result=review.review_reference(tmp_path,scenario,1,plan,'DEVELOPMENT_ONLY')
    expected=11 if model!='GPCM' else 4
    assert result['required_slots']==result['passed_slots']==expected
    assert result['all_checks_pass'] and not result['scientific_inference_ready']
    assert before=={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in tmp_path.rglob('*') if p.is_file()}


def test_missing_tam_input_preserves_both_grid_slots(tmp_path):
    plan=plan_fixture();scenario=plan['scenarios'][0]
    reference_fixture(tmp_path,plan,scenario,unavailable_tam=True)
    result=review.review_reference(tmp_path,scenario,1,plan,'DEVELOPMENT_ONLY')
    assert result['required_slots']==11 and result['passed_slots']==7
    assert len([s for s in result['slots'] if s['lane']=='TAM' and not s['passed']])==4
    assert not result['all_checks_pass']


def test_nonzero_tam_exit_cannot_be_hidden_by_completed_review(tmp_path):
    plan=plan_fixture();scenario=plan['scenarios'][0]
    reference_fixture(tmp_path,plan,scenario,tam_returncode=1)
    result=review.review_reference(tmp_path,scenario,1,plan,'DEVELOPMENT_ONLY')
    assert result['passed_slots']==7 and not result['all_checks_pass']
    assert all(s['execution_success'] is False for s in result['slots'] if s['lane']=='TAM')


@pytest.mark.parametrize('mutation',['missing_grid','duplicate','overall_flag','head','input','tam_artifact'])
def test_changed_or_missing_reference_evidence_is_rejected(tmp_path,mutation):
    plan=plan_fixture();scenario=plan['scenarios'][0]
    summary=reference_fixture(tmp_path,plan,scenario)
    if mutation=='missing_grid':summary['comparisons'].pop()
    elif mutation=='duplicate':summary['comparisons'].append(deepcopy(summary['comparisons'][0]))
    elif mutation=='overall_flag':summary['all_checks_pass']=False
    elif mutation=='head':summary['ledger']['head_sha256']='0'*64
    elif mutation=='input':(tmp_path/'reference_plan.json').write_text('{}')
    else:next((tmp_path/'execution').rglob('*_design.json')).write_text('{}')
    write(tmp_path/'execution/summary.json',summary)
    with pytest.raises(ValueError):review.review_reference(tmp_path,scenario,1,plan,'DEVELOPMENT_ONLY')


def test_wholly_missing_reference_retains_required_slots(tmp_path):
    plan=plan_fixture();scenario=plan['scenarios'][0]
    result=review.review_reference(tmp_path,scenario,1,plan,'DEVELOPMENT_ONLY')
    assert result['required_slots']==11 and result['passed_slots']==0
    assert result['missing_artifacts'] and not result['all_checks_pass']


def evidence_fixture():
    plan=plan_fixture();summaries={};refs=[]
    for scenario in plan['scenarios']:
        count=scenario['replicates'];key=scenario['id']
        if scenario['cohort']=='negative_control':
            summaries[key]=dict(status='complete',planned=count,expected=scenario['expected_status'],
                                observed=[scenario['expected_status']]*len(plan['estimator']['quadrature']),expected_rejection_pass=True)
            continue
        rows=[]
        for truth in scenario_calibration_truth(scenario):
            variable=bool(np.any(truth['jacobian']!=0))
            row=dict(kind=truth['kind'],facet=truth['facet'],level=truth['level'],truth=truth['value'],comparable=variable)
            if not variable:row['status']='fixed_by_identification'
            else:row.update(monte_carlo_intervals=dict(bias=[-.001,.001],se_ratio=[.99,1.01],coverage_given_interval=[.94,.96]),
                            pointwise_acceptance=dict(bias='pass',se='pass',coverage='pass'))
            rows.append(row)
        from tests.test_mml_population_missingness_review import evidence
        comparison=evidence(2000,100)[0]['rsm_pcm_comparison'] if scenario['paired_rsm_pcm_comparison'] else None
        if comparison is not None:
            comparison['conditional_mc_interval']=comparison['all_planned_mc_interval']
            if scenario['primary_model']!='RSM':comparison['null_acceptance']='not_a_null_calibration_test'
        summaries[key]=dict(status='complete',scenario_id=key,planned=count,numerically_eligible=count,
            failure_mc_interval=[0.,.002],failure_acceptance='pass',calibration=rows,
            rsm_pcm_comparison=comparison,fitted_objective_worse_than_truth=[],
            person=[dict(calibration_basis=b,truth_quartile=q,interpretation='DEVELOPMENT_ONLY')
                    for b in ['fitted','known'] for q in ['all',1,2,3,4]])
        for replicate in sorted(set(scenario['independent_R_replicates']+scenario['TAM_replicates'])):
            run,cases=review.reference_inventory(scenario,replicate,plan)
            slots=sum(2 if lane=='TAM' else 1 for _,lane in cases)
            refs.append(dict(run_id=run,required_slots=slots,passed_slots=slots,all_checks_pass=True))
    return plan,summaries,refs


def test_pointwise_conditional_evidence_never_promotes_scientific_inference():
    plan,summaries,refs=evidence_fixture()
    packet=review.assemble_review(plan,summaries,refs)
    assert all(d['registered_conditional_criteria_met'] for d in packet['decisions'])
    assert not packet['scientific_inference_ready']
    assert all(d['review_required'] and not d['scientific_inference_ready'] for d in packet['decisions'])
    assert packet['summaries']==summaries


@pytest.mark.parametrize('mutation',['one_failure','reference','negative_control','inconclusive'])
def test_no_averaging_away_failed_conditions_or_missing_references(mutation):
    plan,summaries,refs=evidence_fixture();first=summaries[plan['scenarios'][0]['id']]
    if mutation=='one_failure':
        first['calibration'][0]['monte_carlo_intervals']['bias']=[.3,.4]
        first['calibration'][0]['pointwise_acceptance']['bias']='fail'
    elif mutation=='inconclusive':
        first['calibration'][0]['monte_carlo_intervals']['bias']=[.09,.11]
        first['calibration'][0]['pointwise_acceptance']['bias']='inconclusive'
    elif mutation=='reference':refs[0].update(passed_slots=0,all_checks_pass=False)
    else:
        control=summaries[plan['scenarios'][-1]['id']]
        control['observed'][0]='wrong_rejection';control['expected_rejection_pass']=False
    result=review.assemble_review(plan,summaries,refs)
    assert not result['decisions'][0]['registered_conditional_criteria_met']
    if mutation!='negative_control':assert result['decisions'][1]['registered_conditional_criteria_met']


@pytest.mark.parametrize('mutation',['missing_scenario','missing_reference','wrong_denominator','missing_person','changed_limit_flag','missing_control_slot'])
def test_invalid_evidence_inventory_or_decision_is_rejected(mutation):
    plan,summaries,refs=evidence_fixture();first=summaries[plan['scenarios'][0]['id']]
    if mutation=='missing_scenario':summaries.pop(plan['scenarios'][1]['id'])
    elif mutation=='missing_reference':refs.pop()
    elif mutation=='wrong_denominator':refs[0].update(required_slots=10,passed_slots=10)
    elif mutation=='missing_person':first['person'].pop()
    elif mutation=='changed_limit_flag':first['calibration'][0]['pointwise_acceptance']['bias']='fail'
    else:summaries[plan['scenarios'][-1]['id']]['observed'].pop()
    with pytest.raises(ValueError):review.assemble_review(plan,summaries,refs)


def test_incomplete_study_reads_no_endpoint_payloads(tmp_path,monkeypatch):
    plan=plan_fixture();scenario=plan['scenarios'][0]
    path=tmp_path/'attempts'/(scenario['id']+'.sqlite');path.parent.mkdir()
    ledger=AttemptLedger(path,plan_sha256='DEVELOPMENT_ONLY')
    ledger.append(scenario['id']+'__0001','terminal',dict(status='generation_failed'));ledger.close()
    monkeypatch.setattr(AttemptLedger,'get',lambda *a,**k:pytest.fail('Endpoint read during incomplete study'))
    result=review.review_study(tmp_path,plan,'DEVELOPMENT_ONLY')
    assert result['status']=='incomplete' and not result['endpoint_payloads_read']
    assert sum(r['terminal'] for r in result['progress'])==1


def test_complete_ledgers_are_reaggregated_and_summary_mutation_rejected(tmp_path):
    plan=plan_fixture()
    for scenario in plan['scenarios']:scenario['replicates']=1
    for scenario in plan['scenarios']:
        path=tmp_path/'attempts'/(scenario['id']+'.sqlite');path.parent.mkdir(exist_ok=True)
        ledger=AttemptLedger(path,plan_sha256='DEVELOPMENT_ONLY')
        ledger.append(scenario['id']+'__0001','terminal',dict(status='generation_failed'))
        summary=summarize_scenario(ledger,scenario,plan);ledger.close()
        write(tmp_path/'summaries'/(scenario['id']+'.json'),summary)
    result=review.review_study(tmp_path,plan,'DEVELOPMENT_ONLY')
    assert result['registered_analysis_recomputed']
    assert not any(d['registered_conditional_criteria_met'] for d in result['decisions'])
    path=tmp_path/'summaries'/(plan['scenarios'][0]['id']+'.json')
    value=json.loads(path.read_bytes());value['failure_acceptance']='pass';write(path,value)
    with pytest.raises(ValueError,match='reaggregation'):review.review_study(tmp_path,plan,'DEVELOPMENT_ONLY')
