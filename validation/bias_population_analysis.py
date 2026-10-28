"""Complete-scenario G3 summaries; no interim endpoints or automatic release.

The independent unit is a generated dataset, not a cell or a response. Raw
cell p-values do not imply family control. Only the fully known-calibration,
known-Person lane has scalar recovery targets equal to the planted shifts.
"""
from collections import Counter

import numpy as np

from validation.bias_operating_metrics import binary_operating_summary, family_rejection
from validation.mfrm_recovery_metrics import calibration_se_summary, summarize_replicate_errors
from validation.mml_population_analysis import binomial_interval, equivalence_status, mc_interval
from validation.mml_population_study import json_bytes


RULES=('raw','holm','statistical_review','strong_review')
ORACLE='known_calibration_known_person'


def _finite(value):
    return type(value) in (int,float) and np.isfinite(value)


def _decisions(cell, screen):
    """Audit retained native decisions rather than trusting a generic UI flag."""
    if cell is None:return dict.fromkeys(RULES)
    if type(cell.get('p_available')) is not bool:raise ValueError('Missing cell availability')
    if not cell['p_available']:return dict.fromkeys(RULES)
    if (cell['status']!='interior' or any(not _finite(cell.get(k)) for k in ['p','p_holm','estimate'])
        or not 0<=cell['p']<=cell['p_holm']<=1):
        raise ValueError('Available test lacks finite interior evidence')
    raw=cell['p']<screen['alpha'];holm=cell['p_holm']<screen['alpha']
    support=cell.get('observations')
    if type(support) is not int or support<1:raise ValueError('Invalid unit-weight cell support')
    statistical=bool(holm and support>=screen['min_positive_weight_observations'])
    strong=bool(statistical and abs(cell['estimate'])>=screen['practical_logit'])
    native=[cell.get('raw_rejection'),cell.get('holm_rejection'),cell.get('operational_statistical_flag')]
    if any(type(v) is not bool for v in native) or native!=[raw,holm,statistical]:
        raise ValueError('Retained native decision differs from declared rule')
    if (cell.get('review_label')=='Strong review')!=strong:
        raise ValueError('Retained native strong-review label differs from declared rule')
    return dict(raw=raw,holm=holm,statistical_review=statistical,strong_review=strong)


def _binary(values, *, planned, acceptable=None):
    result=binary_operating_summary(values,planned=planned)
    result['pointwise_acceptance']=('descriptive' if acceptable is None else 'unavailable' if result['observed']==0
        else equivalence_status(result['all_planned_confidence_envelope'],acceptable))
    return result


def _oracle_recovery(cells, truth, plan):
    errors=[];ses=[];coverage=[];estimate_available=[]
    for cell in cells:
        eligible=(cell is not None and cell['status']=='interior' and _finite(cell.get('estimate')))
        estimate_available.append(eligible)
        error=cell['estimate']-truth if eligible else np.nan
        if eligible and (not _finite(cell.get('oracle_error')) or abs(error-cell['oracle_error'])>1e-12):
            raise ValueError('Oracle error differs from retained DGP shift')
        ci=bool(cell is not None and cell.get('interval_available'))
        covered=cell.get('oracle_covered') if ci else None
        if ci and (not eligible or type(covered) is not bool):
            raise ValueError('Oracle interval lacks a corresponding interior estimate')
        errors.append(error);coverage.append(covered)
        ses.append(cell['conditional_se'] if eligible and _finite(cell.get('conditional_se')) else np.nan)
    metrics=summarize_replicate_errors([[e] for e in errors],[[np.nan if v is None else float(v)] for v in coverage])
    standard_error=calibration_se_summary(errors,ses);limits=plan['draft_interpretation_limits']
    bias_ci=mc_interval(metrics['bias']['estimate'],metrics['bias']['mcse'])
    ratio_ci=mc_interval(standard_error['ratio'],standard_error['mcse_ratio'])
    interval=_binary(coverage,planned=len(cells),acceptable=limits['oracle_coverage'])
    available=int(sum(estimate_available));intervals=sum(v is not None for v in coverage)
    return dict(metrics=metrics,se_comparison=standard_error,coverage=interval,
        estimate_available=available,interval_available=intervals,
        estimate_unavailability_mc_interval=binomial_interval(len(cells)-available,len(cells)),
        interval_unavailability_mc_interval=binomial_interval(len(cells)-intervals,len(cells)),
        monte_carlo_intervals=dict(bias=bias_ci,se_ratio=ratio_ci),
        pointwise_acceptance=dict(bias=equivalence_status(bias_ci,[-limits['oracle_absolute_bias'],limits['oracle_absolute_bias']]),
            se_ratio=equivalence_status(ratio_ci,limits['oracle_SE_ratio']),coverage=interval['pointwise_acceptance']),
        interpretation='Bias and SE ratio condition on their reported eligible repetitions. Coverage retains missing intervals as unknown. Coverage-per-known-target counts successful interval delivery, not conditional nominal coverage.')


def _family(decisions, indices, rule, *, planned, acceptable=None):
    records=[family_rejection([replicate[j][rule] for j in indices]) for replicate in decisions]
    complete=[r['fully_available'] for r in records]
    return dict(any_rejection=_binary([r['any_rejection'] for r in records],planned=planned,acceptable=acceptable),
        fully_available_family=_binary(complete,planned=planned),
        complete_family_event=_binary([r['complete_family_event'] for r in records],planned=planned),
        interpretation='The union can be known from one observed rejection in a partial family; its conditional rate is distinct from the rate among fully available families.')


def summarize_bias_scenario(ledger, scenario, plan):
    """One ledger per scenario. Every planned attempt must be terminal first."""
    registered=next((s for s in plan['scenarios'] if s['id']==scenario['id']),None)
    if registered is None or json_bytes(registered)!=json_bytes(scenario):
        raise ValueError('Scenario differs from the declared plan')
    ledger.verify();count=scenario['replicates']
    run_ids=[f"{scenario['id']}__{i:04d}" for i in range(1,count+1)]
    existing={r[0] for r in ledger.connection.execute('SELECT DISTINCT run_id FROM event')}
    if existing-set(run_ids):raise ValueError('Unexpected RunId in scenario ledger')
    terminal=[ledger.get(run_id,'terminal') for run_id in run_ids]
    if any(t is None for t in terminal):
        return dict(status='incomplete',planned=count,terminal=sum(t is not None for t in terminal),scientific_inference_ready=False)
    spec=scenario['design'];shift=np.asarray(scenario['eta_shift'])
    keys=[(f'R{r}',f'T{t}') for r in range(len(spec['rater'])) for t in range(len(spec['task']))]
    if (shift.shape!=(len(spec['rater']),len(spec['task'])) or not np.isfinite(shift).all()
        or scenario['null']!=bool(np.all(shift==0))):raise ValueError('Invalid declared effect targets')
    truth=shift.ravel().tolist();all_indices=list(range(len(keys)))
    planted=[j for j,b in enumerate(truth) if b!=0];off_target=[j for j,b in enumerate(truth) if b==0]
    limits=plan['draft_interpretation_limits'];lanes=[]
    for basis in plan['scope']['lanes']:
        decisions=[];directional=[];rows=[];statuses=Counter();lane_statuses=Counter()
        for run_id,term in zip(run_ids,terminal,strict=True):
            result=ledger.get(run_id,basis)
            if result is None or 'cells' not in result:
                lane_statuses[(result or {}).get('status',term['status'])]+=1
                cells=[None]*len(keys)
            else:
                if result.get('basis')!=basis or result.get('declared_cells')!=len(keys):
                    raise ValueError('Retained lane identity or target count changed')
                cells=result['cells']
                if [(c['rater'],c['task']) for c in cells]!=keys:raise ValueError('Declared target cells changed')
                if any(c['DGP_shift']!=b or c['oracle_target']!=(basis==ORACLE) for c,b in zip(cells,truth,strict=True)):
                    raise ValueError('Retained cell target differs from declared DGP/lane')
                lane_statuses['computed']+=1
            rows.append(cells);decision=[_decisions(c,plan['screen']) for c in cells];decisions.append(decision)
            directional.append([{rule:(None if d[rule] is None else bool(d[rule] and c['estimate']*b>0)) for rule in RULES}
                for c,d,b in zip(cells,decision,truth,strict=True)])
            statuses.update(c['status'] if c is not None else 'lane_unavailable' for c in cells)
        cell_summaries=[]
        for j,(rater,task) in enumerate(keys):
            cells=[replicate[j] for replicate in rows];is_effect=truth[j]!=0
            target='global_null_false_positive' if scenario['null'] else 'correct_direction_power' if is_effect else 'off_target_flag'
            selected=directional if is_effect else decisions
            rule_summaries={rule:_binary([d[j][rule] for d in selected],planned=count,
                acceptable=[0.,limits['global_null_cell_error_upper']] if scenario['null'] else None) for rule in RULES}
            missing=rule_summaries['raw']['unavailable'];failure_ci=binomial_interval(missing,count)
            cell_summaries.append(dict(rater=rater,task=task,DGP_shift=truth[j],target=target,rules=rule_summaries,
                status_counts=dict(Counter(c['status'] if c is not None else 'lane_unavailable' for c in cells)),
                unavailability_mc_interval=failure_ci,unavailability_acceptance=equivalence_status(failure_ci,[0.,limits['unavailability_upper']]),
                oracle_recovery=_oracle_recovery(cells,truth[j],plan) if basis==ORACLE else None))
        families={}
        for rule in RULES:
            if scenario['null']:
                acceptable=[0.,limits['global_null_family_error_upper']] if rule in limits['family_error_rules'] else None
                families[rule]=dict(target='global_null_family_error',
                    all_cells=_family(decisions,all_indices,rule,planned=count,acceptable=acceptable))
            else:
                acceptable=[limits['alternative_any_correct_strong_detection_lower'],1.] if rule=='strong_review' else None
                families[rule]=dict(target='alternative_detection',
                    any_correct_planted=_family(directional,planted,rule,planned=count,acceptable=acceptable),
                    any_off_target=_family(decisions,off_target,rule,planned=count) if off_target else None)
        complete_missing=sum(any(d['raw'] is None for d in replicate) for replicate in decisions)
        failure_ci=binomial_interval(complete_missing,count)
        lanes.append(dict(basis=basis,lane_status_counts=dict(lane_statuses),cell_status_counts=dict(statuses),
            incomplete_families=complete_missing,family_unavailability_mc_interval=failure_ci,
            family_unavailability_acceptance=equivalence_status(failure_ci,[0.,limits['unavailability_upper']]),
            cells=cell_summaries,families=families,
            recovery_target='DGP eta shift with known baseline and Person' if basis==ORACLE else 'Residual screening only; no DGP-effect interval coverage'))
    return dict(status='complete',classification='COMPLETE_SCENARIO_EVIDENCE_PENDING_REFERENCE_AND_REVIEW',
        scenario_id=scenario['id'],planned=count,execution_status_counts=dict(Counter(t['status'] for t in terminal)),
        primary_lane=plan['scope']['primary_lane'],lanes=lanes,scientific_inference_ready=False,
        independent_reference_qualification='pending separate complete reference-attempt review',
        interpretation='Pointwise MC evidence; no simultaneous grid claim, no averaging away failing cells, no automatic qualification')


def assemble_evidence(plan,summaries,references,controls):
    """No average-over-design qualification, missing reference or target switch."""
    expected=[s['id'] for s in plan['scenarios']]
    if [s['scenario_id'] for s in summaries]!=expected:raise ValueError('Scenario evidence inventory changed')
    required_refs=[f"{s['id']}__{r:04d}" for s in plan['scenarios'] for r in s['independent_reference_replicates']]
    if [r['run_id'] for r in references]!=required_refs:raise ValueError('Reference evidence inventory changed')
    native_controls=plan['deterministic_controls']['scalar_cases']+len(plan['deterministic_controls']['family_cases'])
    reference_controls=plan['deterministic_controls']['scalar_cases']+10*len(plan['deterministic_controls']['family_cases'])
    if controls.get('required_native')!=native_controls or controls.get('required_reference')!=reference_controls:
        raise ValueError('Deterministic control denominator differs from protocol')
    controls_ok=(controls.get('all_checks_pass') is True and controls.get('passed_native')==native_controls
        and controls.get('passed_reference')==reference_controls)
    reference_by_run={r['run_id']:r for r in references};decisions=[]
    states={'pass','fail','inconclusive','unavailable'}
    for scenario,summary in zip(plan['scenarios'],summaries,strict=True):
        if summary['status']!='complete' or summary['planned']!=scenario['replicates']:
            raise ValueError('Incomplete scenario cannot enter the final evidence packet')
        if [lane['basis'] for lane in summary['lanes']]!=plan['scope']['lanes']:
            raise ValueError('Missing or changed calibration/Person basis')
        keys=[(f'R{r}',f'T{t}') for r in range(len(scenario['design']['rater'])) for t in range(len(scenario['design']['task']))]
        for lane in summary['lanes']:
            if [(c['rater'],c['task']) for c in lane['cells']]!=keys:raise ValueError('Missing declared target cells')
            checks=[('family_availability',lane['family_unavailability_acceptance'])]
            if scenario['null']:
                for rule in plan['draft_interpretation_limits']['family_error_rules']:
                    checks.append(('family_error_'+rule,lane['families'][rule]['all_cells']['any_rejection']['pointwise_acceptance']))
            else:
                checks.append(('any_correct_strong_power',lane['families']['strong_review']['any_correct_planted']['any_rejection']['pointwise_acceptance']))
            for cell in lane['cells']:
                key=cell['rater']+'_'+cell['task'];checks.append((key+'_availability',cell['unavailability_acceptance']))
                if scenario['null']:
                    for rule in RULES:checks.append((key+'_error_'+rule,cell['rules'][rule]['pointwise_acceptance']))
                oracle=cell['oracle_recovery']
                if lane['basis']==ORACLE:
                    if oracle is None:raise ValueError('Missing oracle recovery evidence')
                    checks.extend((key+'_'+k,v) for k,v in oracle['pointwise_acceptance'].items())
                elif oracle is not None:raise ValueError('DGP-effect recovery cannot be assigned to a plugin lane')
            if any(value not in states for _,value in checks):raise ValueError('Invalid pointwise decision state')
            ref_states=[]
            for replicate in scenario['independent_reference_replicates']:
                ref=reference_by_run[f"{scenario['id']}__{replicate:04d}"]
                family=next((f for f in ref.get('families',[]) if f['basis']==lane['basis']),None)
                required=len(keys)+4
                ok=(family is not None and family.get('required_slots')==required and family.get('passed_slots')==required
                    and family.get('all_checks_pass') is True and family.get('status')=='completed')
                ref_states.append(dict(replicate=replicate,pass_=bool(ok),status=family['status'] if family else ref.get('status','missing')))
            values=[value for _,value in checks]
            status=('fail' if 'fail' in values else 'unavailable' if 'unavailable' in values
                else 'inconclusive' if 'inconclusive' in values else 'pointwise_criteria_satisfied')
            reference_ok=all(r['pass_'] for r in ref_states)
            decisions.append(dict(scenario_id=scenario['id'],basis=lane['basis'],primary=lane['basis']==plan['scope']['primary_lane'],
                statistical_state=status,pointwise_counts=dict(Counter(values)),checks=[dict(endpoint=k,status=v) for k,v in checks],
                required_references=ref_states,independent_references_complete_and_matching=reference_ok,
                registered_condition_criteria_met=status=='pointwise_criteria_satisfied' and reference_ok and controls_ok,
                review_required=True,scientific_inference_ready=False))
    return dict(classification='COMPLETE_EXECUTION_PENDING_SCIENTIFIC_REVIEW',scenarios=summaries,references=references,
        controls=controls,decisions=decisions,scientific_inference_ready=False,
        interpretation='Each decision is limited to the registered condition/basis. No pooled pass rate, simultaneous grid guarantee, automatic public promotion, or scope extension.')
