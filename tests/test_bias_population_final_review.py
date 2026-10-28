"""Artificial ledger records, not population qualification evidence."""
import copy
import json

import pytest

from validation import bias_population_final_review as review
from validation.bias_population_protocol import make_bias_plan
from validation.bias_population_references import review_family
from validation.mml_population_study import AttemptLedger, json_bytes


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(json_bytes(value)+b'\n')


def file_hashes(root):
    return {str(p.relative_to(root)):review.digest(p) for p in root.rglob('*') if p.is_file()}


def r_output(root, result, returncode=0):
    write(root/'result.json',result)
    for name in ['stdout.txt','stderr.txt']:(root/name).write_text('')
    files=file_hashes(root)
    write(root/'execution.json',dict(returncode=returncode,files=files))
    return dict(result=result,files=file_hashes(root))


def reference_fixture(tmp_path, *, missing=False, failed=False, mismatch=False, returncode=0):
    plan=make_bias_plan();scenario=copy.deepcopy(plan['scenarios'][0])
    scenario.update(replicates=1,independent_reference_replicates=[1]);plan['scenarios']=[scenario]
    run=scenario['id']+'__0001';root=tmp_path/'reference';identity='synthetic'
    bases=plan['scope']['lanes'];ids=[f'R{r}_T{t}' for r in range(3) for t in range(2)]
    row={'Cell Status':'interior','Reference Status':'unit_weight_conditional_lr_candidate','Bias Size':0.,
         'S.E.':.2,'LR ChiSq':0.,'LR Prob.':1.,'Profile CI Lower':-.4,'Profile CI Upper':.4,'Profile CI Status':'ok'}
    native=dict(cases=[dict(id=key,row=copy.deepcopy(row)) for key in ids],families=[[1.]*6],observed_ids=ids)
    reference=dict(cases=[dict(id=key,status='interior',estimate=0.,se=.2,lr=0.,p=1.,ci=[-.4,.4],ci_status='ok') for key in ids],
        adjustments=[dict(holm=[1.]*6,BH=[1.]*6)],family_observed_ids=ids,derived_adjustments=[dict(holm=[1.]*6,BH=[1.]*6)])
    families=[]
    for i,basis in enumerate(bases):
        available=not (missing and i==0)
        families.append(dict(basis=basis,available=available,input_status='available' if available else 'native_screen_unavailable',
            required_scalar_cells=6,required_adjustments=2,required_derived_adjustments=2))
        if available:write(root/basis/'python.json',native)
    write(root/'selection.json',dict(study_sha256=identity,scenario_id=scenario['id'],run_id=run,replicate=1,families=families))
    write(root/'manifest.json',file_hashes(root))
    export=dict(status='completed',value=dict(manifest_sha256=review.digest(root/'manifest.json')))
    output=root/'execution';output.mkdir()
    ledger=AttemptLedger(output/'events.sqlite',plan_sha256=identity)
    for i,family in enumerate(families):
        if not family['available']:continue
        basis=family['basis'];key=run+'_'+basis
        ledger.append(key,'R_started',dict(single_attempt=True))
        if failed and i==0:
            ledger.append(key,'R_finished',dict(status='reference_failed'));continue
        r=copy.deepcopy(reference)
        if mismatch and i==0:r['cases'][0]['estimate']=.1
        value=r_output(output/basis,r,returncode if i==0 else 0)
        ledger.append(key,'R_finished',dict(status='completed',value=value))
        ledger.append(key,'review_started',dict(single_attempt=True))
        ledger.append(key,'review_finished',dict(status='completed',value=review_family(native,r,plan)))
    ledger.close()
    return root,scenario,plan,identity,export


def run_reference(fixture):
    root,scenario,plan,identity,export=fixture
    return review.review_reference(root,scenario,1,plan,identity,export)


@pytest.mark.parametrize('missing,failed,mismatch,expected',[(False,False,False,40),(True,False,False,30),(False,True,False,30),(False,False,True,39)])
def test_all_lanes_and_failed_slots_are_preserved_without_writes(tmp_path,missing,failed,mismatch,expected):
    fixture=reference_fixture(tmp_path,missing=missing,failed=failed,mismatch=mismatch)
    before=file_hashes(tmp_path);result=run_reference(fixture)
    assert file_hashes(tmp_path)==before
    assert result['required_slots']==40 and result['passed_slots']==expected
    assert result['all_checks_pass'] is (expected==40)
    assert not result['scientific_inference_ready']


@pytest.mark.parametrize('mutation',['hash','lane','denominator','availability','replicate'])
def test_changed_reference_inventory_is_rejected(tmp_path,mutation):
    root,scenario,plan,identity,export=fixture=reference_fixture(tmp_path)
    selection=review.read(root/'selection.json')
    if mutation=='hash':
        (root/plan['scope']['lanes'][0]/'python.json').write_text('{}')
    else:
        if mutation=='lane':selection['families'][1]['basis']=selection['families'][0]['basis']
        if mutation=='denominator':selection['families'][0]['required_adjustments']=1
        if mutation=='availability':selection['families'][0]['available']=1
        if mutation=='replicate':selection['replicate']=2
        write(root/'selection.json',selection)
        manifest=review.read(root/'manifest.json');manifest['selection.json']=review.digest(root/'selection.json')
        write(root/'manifest.json',manifest);export['value']['manifest_sha256']=review.digest(root/'manifest.json')
    with pytest.raises(ValueError):run_reference(fixture)


def test_r_numerical_match_cannot_override_nonzero_exit(tmp_path):
    fixture=reference_fixture(tmp_path,returncode=2)
    with pytest.raises(ValueError,match='nonzero'):run_reference(fixture)


def test_changed_retained_r_artifact_is_rejected(tmp_path):
    fixture=reference_fixture(tmp_path);root,_,plan,_,_=fixture
    (root/'execution'/plan['scope']['lanes'][0]/'stderr.txt').write_text('changed')
    with pytest.raises(ValueError,match='artifact changed'):run_reference(fixture)


def test_changed_summary_does_not_override_recomputed_reference(tmp_path):
    fixture=reference_fixture(tmp_path,mismatch=True);result=run_reference(fixture)
    result['all_checks_pass']=True;write(fixture[0]/'execution/summary.json',result)
    with pytest.raises(ValueError,match='Reference summary'):run_reference(fixture)


def test_changed_review_event_is_detected_by_numeric_recomputation(tmp_path):
    fixture=reference_fixture(tmp_path);root,scenario,plan,identity,_=fixture
    path=root/'execution/events.sqlite';source=AttemptLedger.read_only(path,plan_sha256=identity)
    events=[(run,phase,source.get(run,phase)) for run,phase in source.connection.execute('SELECT run_id,phase FROM event ORDER BY sequence')]
    source.close();path.unlink()
    ledger=AttemptLedger(path,plan_sha256=identity)
    for run,phase,value in events:
        if phase=='review_finished':value['value']['cells'][0]['pass_']=False
        ledger.append(run,phase,value)
    ledger.close()
    with pytest.raises(ValueError,match='Reference comparison'):run_reference(fixture)


def test_missing_export_retains_every_required_slot(tmp_path):
    plan=make_bias_plan();scenario=plan['scenarios'][0]
    result=review.review_reference(tmp_path,scenario,1,plan,'synthetic',None)
    assert result['required_slots']==40 and result['passed_slots']==0 and not result['all_checks_pass']


def test_unfinished_scenario_prevents_all_endpoint_access(tmp_path,monkeypatch):
    plan=make_bias_plan();plan['scenarios']=plan['scenarios'][:2]
    for scenario in plan['scenarios']:scenario['replicates']=1
    scenario=plan['scenarios'][0];path=tmp_path/'attempts'/(scenario['id']+'.sqlite');path.parent.mkdir()
    ledger=AttemptLedger(path,plan_sha256='synthetic');ledger.append(scenario['id']+'__0001','terminal',dict(status='generation_failed'));ledger.close()
    def forbidden(*args,**kwargs):raise AssertionError('Endpoint access while incomplete')
    monkeypatch.setattr(AttemptLedger,'get',forbidden)
    monkeypatch.setattr(review,'review_controls',forbidden)
    result=review.review_study(tmp_path,plan,'synthetic')
    assert result['status']=='incomplete' and result['endpoint_payloads_read'] is False
    assert [r['terminal'] for r in result['progress']]==[1,0]


def test_extra_run_cannot_be_selected_away(tmp_path):
    plan=make_bias_plan();plan['scenarios']=plan['scenarios'][:1]
    scenario=plan['scenarios'][0];path=tmp_path/'attempts'/(scenario['id']+'.sqlite');path.parent.mkdir()
    ledger=AttemptLedger(path,plan_sha256='synthetic');ledger.append('extra','terminal',dict(status='generation_failed'));ledger.close()
    with pytest.raises(ValueError,match='Unexpected RunId'):review.review_study(tmp_path,plan,'synthetic')


def complete_study(tmp_path,monkeypatch):
    from tests.test_bias_population_analysis import fixture_plan,cells_for,append_sample
    from validation.bias_population_analysis import summarize_bias_scenario
    plan,scenario=fixture_plan();scenario['independent_reference_replicates']=[1,3]
    path=tmp_path/'attempts'/(scenario['id']+'.sqlite');path.parent.mkdir()
    ledger=AttemptLedger(path,plan_sha256='synthetic')
    for rep in range(1,5):append_sample(ledger,scenario,rep,{b:cells_for(scenario,b) for b in plan['scope']['lanes']})
    summary=summarize_bias_scenario(ledger,scenario,plan);summary.update(study_sha256='synthetic',source_ledger=ledger.verify());ledger.close()
    write(tmp_path/'summaries'/(scenario['id']+'.json'),summary)
    AttemptLedger(tmp_path/'reference_exports.sqlite',plan_sha256='synthetic').close()
    monkeypatch.setattr(review,'review_controls',lambda *a:dict(required_native=20,passed_native=20,required_reference=47,passed_reference=47,all_checks_pass=True))
    return plan,scenario


def test_complete_study_reaggregates_raw_data_but_missing_references_block_qualification(tmp_path,monkeypatch):
    plan,_=complete_study(tmp_path,monkeypatch);before=file_hashes(tmp_path)
    result=review.review_study(tmp_path,plan,'synthetic')
    assert file_hashes(tmp_path)==before and result['registered_analysis_recomputed']
    assert all(not d['registered_condition_criteria_met'] for d in result['decisions'])
    assert not result['scientific_inference_ready']
    assert result['stored_packet_matches'] is None


@pytest.mark.parametrize('mutation',['summary','packet'])
def test_complete_study_detects_changed_saved_decisions(tmp_path,monkeypatch,mutation):
    plan,scenario=complete_study(tmp_path,monkeypatch)
    if mutation=='summary':
        path=tmp_path/'summaries'/(scenario['id']+'.json');value=review.read(path)
        value['lanes'][0]['cells'][0]['rules']['raw']['pointwise_acceptance']='pass';write(path,value)
    else:write(tmp_path/'evidence_packet.json',dict(scientific_inference_ready=True))
    with pytest.raises(ValueError,match='differs from retained'):review.review_study(tmp_path,plan,'synthetic')


def test_hash_verification_rejects_escape(tmp_path):
    outside=tmp_path/'outside.json';outside.write_text('{}');root=tmp_path/'root';root.mkdir()
    with pytest.raises(ValueError,match='artifact changed'):
        review.verify_files(root,{'../outside.json':review.digest(outside)})
