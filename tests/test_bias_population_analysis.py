"""Analytically assigned events, not simulated qualification results."""
import copy

import numpy as np
import pytest

from validation.bias_population_analysis import ORACLE, summarize_bias_scenario
from validation.bias_population_protocol import make_bias_plan
from validation.mml_population_study import AttemptLedger


def fixture_plan(*, effect=False):
    plan=make_bias_plan();scenario=copy.deepcopy(plan['scenarios'][int(effect)])
    scenario['replicates']=4;plan['scenarios']=[scenario]
    return plan,scenario


def cells_for(scenario,basis):
    rows=[]
    for r,shifts in enumerate(scenario['eta_shift']):
        for t,truth in enumerate(shifts):
            rows.append(dict(rater=f'R{r}',task=f'T{t}',DGP_shift=truth,status='interior',
                p_available=True,interval_available=True,oracle_target=basis==ORACLE,
                estimate=0.,conditional_se=.2,p=1.,p_holm=1.,oracle_error=-truth if basis==ORACLE else None,
                oracle_covered=True if basis==ORACLE else None,raw_rejection=False,holm_rejection=False,
                operational_statistical_flag=False,review_label='No flag',review_flag=False,observations=8))
    return rows


def reject(cell,*,estimate=.9,strong=True):
    cell.update(estimate=estimate,p=.001,p_holm=.006,raw_rejection=True,holm_rejection=True,
        operational_statistical_flag=True,review_label='Strong review' if strong else 'Statistical review',review_flag=True)
    if cell['oracle_target']:cell['oracle_error']=estimate-cell['DGP_shift']


def unavailable(cell):
    cell.update(status='no_finite_mle',p_available=False,interval_available=False,estimate=None,conditional_se=None,
        p=None,p_holm=None,oracle_error=None,oracle_covered=None,raw_rejection=None,holm_rejection=None,
        operational_statistical_flag=None,review_label='Inference unavailable',review_flag=True)


def append_sample(ledger,scenario,replicate,lanes,*,terminal=True):
    run=f"{scenario['id']}__{replicate:04d}"
    for basis,cells in lanes.items():
        ledger.append(run,basis,dict(basis=basis,cells=cells,declared_cells=len(cells)))
    if terminal:ledger.append(run,'terminal',dict(status='completed' if lanes else 'generation_failed'))


def test_no_endpoints_until_every_planned_attempt_is_terminal(tmp_path):
    plan,scenario=fixture_plan();ledger=AttemptLedger(tmp_path/'events.sqlite',plan_sha256='fixture')
    append_sample(ledger,scenario,1,{basis:cells_for(scenario,basis) for basis in plan['scope']['lanes']})
    result=summarize_bias_scenario(ledger,scenario,plan)
    assert result==dict(status='incomplete',planned=4,terminal=1,scientific_inference_ready=False)
    ledger.close()


def test_partial_family_null_event_and_oracle_only_recovery(tmp_path):
    plan,scenario=fixture_plan();ledger=AttemptLedger(tmp_path/'events.sqlite',plan_sha256='fixture')
    for replicate in range(1,5):
        lanes={basis:cells_for(scenario,basis) for basis in plan['scope']['lanes']}
        if replicate==1:
            for rows in lanes.values():reject(rows[0]);unavailable(rows[1])
        elif replicate==3:lanes={}
        elif replicate==4:
            # UI practical/unavailable review flags must not become LR rejections.
            for rows in lanes.values():
                rows[0].update(estimate=.8,review_label='Practical review',review_flag=True)
                if rows[0]['oracle_target']:rows[0]['oracle_error']=.8
        append_sample(ledger,scenario,replicate,lanes)
    result=summarize_bias_scenario(ledger,scenario,plan);ledger.close()
    for lane in result['lanes']:
        family=lane['families']['holm']['all_cells']
        assert family['any_rejection']['conditional_rate']==pytest.approx(1/3)
        assert family['any_rejection']['all_planned_identification_interval']==[.25,.5]
        assert family['complete_family_event']['conditional_rate']==0
        assert family['complete_family_event']['observed']==2
        assert lane['incomplete_families']==2
        assert lane['families']['raw']['all_cells']['any_rejection']['pointwise_acceptance']=='descriptive'
        assert lane['cell_status_counts']=={'interior':17,'no_finite_mle':1,'lane_unavailable':6}
        oracle=lane['cells'][0]['oracle_recovery']
        assert (oracle is not None)==(lane['basis']==ORACLE)
        if oracle is not None:
            assert oracle['metrics']['bias']['estimate']==pytest.approx(1.7/3)
            assert oracle['coverage']['all_planned_identification_interval']==[.75,1.]
    assert not result['scientific_inference_ready']


def test_alternative_direction_off_target_and_strong_threshold_are_distinct(tmp_path):
    plan,scenario=fixture_plan(effect=True);ledger=AttemptLedger(tmp_path/'events.sqlite',plan_sha256='fixture')
    for replicate in range(1,5):
        lanes={basis:cells_for(scenario,basis) for basis in plan['scope']['lanes']}
        for rows in lanes.values():
            if replicate==1:reject(rows[0],estimate=-.9);reject(rows[1],estimate=-.9)
            if replicate==2:reject(rows[4])
            if replicate==4:reject(rows[0],estimate=.1,strong=False)
        append_sample(ledger,scenario,replicate,{} if replicate==3 else lanes)
    result=summarize_bias_scenario(ledger,scenario,plan);ledger.close()
    for lane in result['lanes']:
        family=lane['families']['strong_review']
        assert family['any_correct_planted']['any_rejection']['successes']==1
        assert family['any_off_target']['any_rejection']['successes']==1
        assert family['any_off_target']['any_rejection']['pointwise_acceptance']=='descriptive'
        assert lane['cells'][0]['rules']['strong_review']['successes']==0
        assert lane['cells'][0]['rules']['statistical_review']['successes']==1
        assert lane['cells'][4]['target']=='off_target_flag'


def test_oracle_empirical_sd_and_rms_se_use_the_same_replicates(tmp_path):
    plan,scenario=fixture_plan();ledger=AttemptLedger(tmp_path/'events.sqlite',plan_sha256='fixture')
    errors=np.array([-.3,-.1,.1,.3]);ses=np.array([.1,.2,.3,.4])
    for replicate,(error,se) in enumerate(zip(errors,ses),1):
        rows=cells_for(scenario,ORACLE)
        rows[0].update(estimate=float(error),oracle_error=float(error),conditional_se=float(se))
        append_sample(ledger,scenario,replicate,{ORACLE:rows})
    result=summarize_bias_scenario(ledger,scenario,plan);ledger.close()
    recovery=result['lanes'][0]['cells'][0]['oracle_recovery']
    assert recovery['metrics']['bias']['estimate']==pytest.approx(0.,abs=1e-15)
    assert recovery['metrics']['bias']['mcse']==pytest.approx(np.std(errors,ddof=1)/2)
    assert recovery['metrics']['rmse']==pytest.approx(np.sqrt(np.mean(errors**2)))
    assert recovery['se_comparison']['ratio']==pytest.approx(np.sqrt(np.mean(ses**2))/np.std(errors,ddof=1))
    assert result['lanes'][1]['cells'][0]['rules']['raw']['all_planned_confidence_envelope']==[0.,1.]


@pytest.mark.parametrize('mutation',['truth','oracle','rule','strong','p','duplicate'])
def test_changed_targets_and_native_decisions_are_rejected(tmp_path,mutation):
    plan,scenario=fixture_plan();ledger=AttemptLedger(tmp_path/'events.sqlite',plan_sha256='fixture')
    rows=cells_for(scenario,ORACLE)
    if mutation=='truth':rows[0]['DGP_shift']=.1
    if mutation=='oracle':rows[0]['oracle_target']=False
    if mutation=='rule':rows[0]['raw_rejection']=True
    if mutation=='strong':rows[0]['review_label']='Strong review'
    if mutation=='p':rows[0]['p_holm']=.9
    if mutation=='duplicate':rows[1]=copy.deepcopy(rows[0])
    append_sample(ledger,scenario,1,{ORACLE:rows})
    for replicate in range(2,5):append_sample(ledger,scenario,replicate,{})
    with pytest.raises(ValueError):summarize_bias_scenario(ledger,scenario,plan)
    ledger.close()


def test_unexpected_run_and_unregistered_scenario_cannot_be_selected_away(tmp_path):
    plan,scenario=fixture_plan();ledger=AttemptLedger(tmp_path/'events.sqlite',plan_sha256='fixture')
    ledger.append('extra_run','terminal',dict(status='generation_failed'))
    with pytest.raises(ValueError,match='Unexpected RunId'):summarize_bias_scenario(ledger,scenario,plan)
    changed=copy.deepcopy(scenario);changed['replicates']=1
    with pytest.raises(ValueError,match='Scenario differs'):summarize_bias_scenario(ledger,changed,plan)
    ledger.close()


def test_independent_reference_requires_all_eligible_values_and_adjustment_slots():
    from validation.bias_population_check import compare_reference
    limits=dict(estimate=1e-8,se=1e-8,lr=1e-8,p=1e-8,ci=1e-7,adjustment=1e-14)
    row={'Cell Status':'interior','Reference Status':'unit_weight_conditional_lr_candidate',
        'Bias Size':0.,'S.E.':.2,'LR ChiSq':0.,'LR Prob.':1.,
        'Profile CI Lower':-.4,'Profile CI Upper':.4,'Profile CI Status':'ok'}
    py=[dict(id='fixture',row=row)]
    case=dict(id='fixture',status='interior',estimate=0.,se=.2,lr=0.,p=1.,ci=[-.4,.4],ci_status='ok')
    reference=dict(cases=[case],adjustments=[dict(holm=[1.,1.],BH=[1.,1.])])
    checks,adjustments=compare_reference(py,[[1.,1.]],reference,limits)
    assert all(r['pass_'] for r in checks+adjustments)
    del case['p']
    checks,_=compare_reference(py,[[1.,1.]],reference,limits)
    assert checks[0]['missing_reference_fields']==['p'] and not checks[0]['pass_']
    reference['adjustments'][0]['holm']=[1.]
    _,adjustments=compare_reference(py,[[1.,1.]],reference,limits)
    assert not adjustments[0]['pass_']
    with pytest.raises(ValueError):compare_reference(py,[],dict(cases=[],adjustments=[]),limits)


def evidence_fixture():
    """Explicitly invented decision records; these contain no study outcomes."""
    from validation.bias_population_analysis import RULES
    plan=make_bias_plan();plan['scenarios']=plan['scenarios'][:2]
    summaries=[];references=[]
    for scenario in plan['scenarios']:
        scenario.update(replicates=4,independent_reference_replicates=[1,3]);lanes=[]
        for basis in plan['scope']['lanes']:
            cells=[]
            for r in range(3):
                for t in range(2):
                    cells.append(dict(rater=f'R{r}',task=f'T{t}',unavailability_acceptance='pass',
                        rules={rule:dict(pointwise_acceptance='pass') for rule in RULES},
                        oracle_recovery=dict(pointwise_acceptance=dict(bias='pass',se_ratio='pass',coverage='pass')) if basis==ORACLE else None))
            families={rule:dict(all_cells=dict(any_rejection=dict(pointwise_acceptance='pass'))) for rule in RULES}
            if not scenario['null']:
                families['strong_review']['any_correct_planted']=dict(any_rejection=dict(pointwise_acceptance='pass'))
            lanes.append(dict(basis=basis,cells=cells,families=families,family_unavailability_acceptance='pass'))
        summaries.append(dict(scenario_id=scenario['id'],status='complete',planned=4,lanes=lanes))
        for replicate in [1,3]:
            references.append(dict(run_id=f"{scenario['id']}__{replicate:04d}",families=[
                dict(basis=basis,status='completed',required_slots=10,passed_slots=10,all_checks_pass=True)
                for basis in plan['scope']['lanes']]))
    controls=dict(required_native=20,passed_native=20,required_reference=47,passed_reference=47,all_checks_pass=True)
    return plan,summaries,references,controls


def test_failed_condition_is_not_averaged_away_and_no_result_is_publicly_promoted():
    from validation.bias_population_analysis import assemble_evidence
    plan,summaries,references,controls=evidence_fixture()
    summaries[0]['lanes'][3]['cells'][0]['rules']['raw']['pointwise_acceptance']='fail'
    packet=assemble_evidence(plan,summaries,references,controls)
    assert packet['decisions'][3]['statistical_state']=='fail'
    assert not packet['decisions'][3]['registered_condition_criteria_met']
    assert packet['decisions'][7]['statistical_state']=='pointwise_criteria_satisfied'
    assert not packet['scientific_inference_ready'] and not any(d['scientific_inference_ready'] for d in packet['decisions'])


def test_missing_reference_prevents_positive_decision_even_when_operating_criteria_pass():
    from validation.bias_population_analysis import assemble_evidence
    plan,summaries,references,controls=evidence_fixture()
    references[0]=dict(run_id=references[0]['run_id'],status='interrupted_reference_retained',all_checks_pass=False)
    packet=assemble_evidence(plan,summaries,references,controls)
    assert not any(d['registered_condition_criteria_met'] for d in packet['decisions'][:4])
    assert all(d['statistical_state']=='pointwise_criteria_satisfied' for d in packet['decisions'][:4])
    assert all(d['registered_condition_criteria_met'] for d in packet['decisions'][4:])
    with pytest.raises(ValueError,match='Reference evidence inventory'):assemble_evidence(plan,summaries,references[:-1],controls)


def test_control_boolean_alone_does_not_replace_planned_denominator():
    from validation.bias_population_analysis import assemble_evidence
    plan,summaries,references,controls=evidence_fixture();controls['passed_reference']=46
    packet=assemble_evidence(plan,summaries,references,controls)
    assert not any(d['registered_condition_criteria_met'] for d in packet['decisions'])
    controls['required_reference']=46
    with pytest.raises(ValueError,match='control denominator'):assemble_evidence(plan,summaries,references,controls)


def test_plugin_lane_cannot_acquire_oracle_coverage_during_final_assembly():
    from validation.bias_population_analysis import assemble_evidence
    plan,summaries,references,controls=evidence_fixture()
    summaries[0]['lanes'][1]['cells'][0]['oracle_recovery']=dict(pointwise_acceptance=dict(coverage='pass'))
    with pytest.raises(ValueError,match='plugin lane'):assemble_evidence(plan,summaries,references,controls)
