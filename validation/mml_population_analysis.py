"""Frozen, failure-aware endpoint summaries after a complete scenario ledger."""
import numpy as np
from scipy.special import ndtri
from scipy.stats import beta

from validation.mfrm_recovery_metrics import (calibration_coordinates, calibration_se_summary,
    population_truth_quartile, summarize_replicate_errors)
from validation.mml_population_study import truth_joint


def binomial_interval(successes, trials):
    """Exact two-sided 95% Clopper-Pearson Monte Carlo interval."""
    if type(successes) is not int or type(trials) is not int or not 0 <= successes <= trials:
        raise ValueError('Integer binary counts required')
    if trials == 0: return None
    return [0. if successes==0 else float(beta.ppf(.025, successes, trials-successes+1)),
            1. if successes==trials else float(beta.ppf(.975, successes+1, trials-successes))]


def mc_interval(value, se):
    if value is None or se is None or not np.isfinite([value,se]).all(): return None
    return [float(value-ndtri(.975)*se), float(value+ndtri(.975)*se)]


def equivalence_status(interval, acceptable):
    if interval is None: return 'unavailable'
    lo, hi = interval; a, b = acceptable
    if not np.isfinite([lo,hi,a,b]).all() or lo > hi or a > b: raise ValueError('Invalid acceptance interval')
    if lo >= a and hi <= b: return 'pass'
    if hi < a or lo > b: return 'fail'
    return 'inconclusive'


def scenario_calibration_truth(scenario):
    spec = scenario['design']; levels = dict(Rater=[f'R{i}' for i in range(len(spec['rater']))],
        Task=[f'T{i}' for i in range(len(spec['task']))], Criterion=[f'C{i}' for i in range(len(spec['criterion']))])
    config = dict(model=spec['model'], facet_names=['Rater','Task','Criterion'], facet_levels=levels,
        noncenter_facet='Criterion', step_facet=None if spec['model']=='RSM' else 'Criterion',
        slope_facet='Criterion' if spec['model']=='GPCM' else None, n_cat=len(spec['steps'][0])+1)
    truth = dict(Rater=spec['rater'], Task=spec['task'], Criterion=spec['criterion'], Step=spec['steps'],
        LogSlope=spec['log_slopes'], PopulationSD=spec['population_sd'])
    return calibration_coordinates(config, truth_joint(config, truth))


def summarize_scenario(ledger, scenario, plan):
    """No endpoint peeking on partial replicate sets; no hidden denominator."""
    ledger.verify(); count = scenario['replicates']; model = scenario['primary_model']
    run_ids = [f"{scenario['id']}__{i:04d}" for i in range(1,count+1)]
    existing = {r[0] for r in ledger.connection.execute('SELECT DISTINCT run_id FROM event')}
    if existing-set(run_ids): raise ValueError('Unexpected RunId in scenario ledger')
    terminal = [ledger.get(run_id,'terminal') for run_id in run_ids]
    if any(t is None for t in terminal):
        return dict(status='incomplete', planned=count, terminal=sum(t is not None for t in terminal), scientific_inference_ready=False)
    if scenario['cohort']=='negative_control':
        statuses = []
        for run_id in run_ids:
            for q in plan['estimator']['quadrature']:
                fit = ledger.get(run_id,f'{model}_Q{q}_finished')
                statuses.append(fit['status'] if fit else 'missing_or_execution_failed')
        return dict(status='complete', planned=count, expected=scenario['expected_status'], observed=statuses,
            expected_rejection_pass=all(s==scenario['expected_status'] for s in statuses), scientific_inference_ready=False)
    truths = scenario_calibration_truth(scenario)
    keys = [(r['kind'],r['facet'],r['level']) for r in truths]
    error = np.full((count,len(keys)),np.nan); covered=error.copy(); se=error.copy()
    eligible=np.zeros(count,dtype=bool); available_test=np.zeros(count,dtype=bool); rejects=np.zeros(count,dtype=bool)
    persons = {(basis,quartile):dict(error=[],coverage=[]) for basis in ['fitted','known'] for quartile in ['all',1,2,3,4]}
    execution={}; truth_unknown=0; worse_than_truth=[]
    criteria={name:dict(available=0,PCM_lower=0,tied=0) for name in ['AIC','BIC']}
    for i,run_id in enumerate(run_ids):
        term=terminal[i]; execution[term['status']]=execution.get(term['status'],0)+1
        assessment=ledger.get(run_id,model+'_assessment'); endpoints=(assessment or {}).get('endpoints')
        retained=ledger.get(run_id,'input_retained')
        eligible[i]=assessment is not None and assessment['status']=='eligible_candidate'
        if endpoints:
            difference=endpoints.get('fitted_nll_minus_truth_nll')
            if difference is not None and difference > 1e-8:
                worse_than_truth.append(dict(run_id=run_id,nll_difference=difference,operationally_eligible=bool(eligible[i])))
            rows=endpoints['calibration']
            if [(r['kind'],r['facet'],r['level']) for r in rows] != keys: raise ValueError('Calibration rows differ from registered native coordinates')
            for j,row in enumerate(rows):
                if row['truth'] != truths[j]['value']: raise ValueError('Calibration truth changed across replicates')
                if row['eligible']:
                    if row['estimate'] is None or row['covered'] is None or row['se'] is None: raise ValueError('Eligible row lacks its estimate or interval')
                    error[i,j]=row['estimate']-row['truth'];covered[i,j]=float(row['covered']);se[i,j]=row['se']
        comparison=term.get('comparison')
        if comparison and comparison['available']:
            available_test[i]=True;rejects[i]=comparison['p_value'] < .05
            ic=comparison.get('information_criteria',{})
            for name,record in criteria.items():
                values=[ic.get(m,{}).get(name) for m in ['RSM','PCM']]
                if all(v is not None and np.isfinite(v) for v in values):
                    record['available']+=1;record['PCM_lower']+=int(values[1]<values[0]);record['tied']+=int(values[1]==values[0])
        if retained is None:
            truth_unknown+=1
            for record in persons.values(): record['error'].append(None);record['coverage'].append(None)
            continue
        theta=np.array(retained['truth']['Person']);quartiles=population_truth_quartile(theta,scenario['design']['population_sd'])
        person=(endpoints or {}).get('person',{})
        if person and not np.array_equal(theta,np.array(person['truth'])): raise ValueError('Person truth differs from retained generation')
        for basis in ['fitted','known']:
            ok=bool(person.get('eligible',False)) if basis=='fitted' else bool(person.get('known_calibration_numerically_eligible',False))
            prefix='fitted' if basis=='fitted' else 'known_calibration'
            est=np.array(person[prefix+'_eap']) if ok else np.full(len(theta),np.nan)
            sd=np.array(person[prefix+'_posterior_sd']) if ok else np.full(len(theta),np.nan)
            coverage=np.where(np.isfinite(est)&np.isfinite(sd),(abs(est-theta)<=ndtri(.975)*sd).astype(float),np.nan)
            for quartile in ['all',1,2,3,4]:
                select=np.ones(len(theta),dtype=bool) if quartile=='all' else quartiles==quartile
                persons[basis,quartile]['error'].append((est-theta)[select])
                persons[basis,quartile]['coverage'].append(coverage[select])
    limits=plan['interpretation_limits'];calibration=[]
    for j,truth in enumerate(truths):
        comparable=bool(np.any(truth['jacobian']!=0))
        record=dict(kind=truth['kind'],facet=truth['facet'],level=truth['level'],truth=truth['value'],comparable=comparable)
        if not comparable:
            calibration.append(dict(record,status='fixed_by_identification'));continue
        metrics=summarize_replicate_errors([[v] for v in error[:,j]],[[v] for v in covered[:,j]])
        standard_error=calibration_se_summary(error[:,j],se[:,j])
        successes=int(np.nansum(covered[:,j]));intervals=int(np.isfinite(covered[:,j]).sum())
        coverage_ci=binomial_interval(successes,intervals)
        per_target_ci=binomial_interval(successes,count)
        bias_ci=mc_interval(metrics['bias']['estimate'],metrics['bias']['mcse'])
        positive=truth['kind'] in {'Slope','PopulationSD'}
        bias_limit=limits['positive_parameter_relative_bias']*truth['value'] if positive else limits['calibration_absolute_bias']
        ratio_ci=mc_interval(standard_error['ratio'],standard_error['mcse_ratio'])
        calibration.append(dict(record,metrics=metrics,se_comparison=standard_error,
            monte_carlo_intervals=dict(bias=bias_ci,se_ratio=ratio_ci,coverage_given_interval=coverage_ci,coverage_per_known_target=per_target_ci),
            pointwise_acceptance=dict(bias=equivalence_status(bias_ci,[-bias_limit,bias_limit]),
                se=equivalence_status(ratio_ci,limits['rms_model_se_over_empirical_sd']),
                coverage=equivalence_status(coverage_ci,limits['coverage_95']))))
    person_summaries=[]
    for (basis,quartile),record in persons.items():
        metrics=summarize_replicate_errors(record['error'],record['coverage'])
        person_summaries.append(dict(calibration_basis=basis,truth_quartile=quartile,metrics=metrics,
            coverage_mc_interval=mc_interval(metrics['coverage_given_interval']['estimate'],metrics['coverage_given_interval']['mcse']),
            interpretation='Calibration-fixed posterior normal interval; truth-bin coverage and EAP shrinkage bias are descriptive, not a per-bin 95% or unbiasedness requirement'))
    failed=int((~eligible).sum());failure_ci=binomial_interval(failed,count)
    test_n=int(available_test.sum());rejection_n=int(rejects.sum());test_ci=binomial_interval(rejection_n,test_n)
    test_summary=dict(planned=count,available=test_n,unavailable=count-test_n,rejected=rejection_n,
        conditional_rate=rejection_n/test_n if test_n else None,all_planned_rate=rejection_n/count,
        conditional_mcse=float(np.sqrt((rejection_n/test_n)*(1-rejection_n/test_n)/test_n)) if test_n else None,
        all_planned_mcse=float(np.sqrt((rejection_n/count)*(1-rejection_n/count)/count)),
        conditional_mc_interval=test_ci,all_planned_mc_interval=binomial_interval(rejection_n,count),
        target='null rejection' if model=='RSM' else 'power under declared step variation',
        null_acceptance=equivalence_status(test_ci,limits['null_LRT_rejection']) if model=='RSM' and scenario.get('paired_rsm_pcm_comparison') else 'not_a_null_calibration_test')
    for name,record in criteria.items():
        n=record['available'];hits=record['PCM_lower']
        record.update(planned=count,unavailable=count-n,PCM_lower_conditional_rate=hits/n if n else None,
            PCM_lower_mc_interval=binomial_interval(hits,n),interpretation='descriptive same-input candidate criteria; no automatic model recommendation')
    test_summary['information_criteria']=criteria
    return dict(status='complete',classification='PROSPECTIVE_ENDPOINT_EVIDENCE_PENDING_REFERENCE_AND_REVIEW',
        scenario_id=scenario['id'],planned=count,execution_status_counts=execution,
        numerically_eligible=int(eligible.sum()),failed_or_ineligible=failed,truth_unknown_replicates=truth_unknown,
        fitted_objective_worse_than_truth=worse_than_truth,
        failure_rate=failed/count,failure_mcse=float(np.sqrt((failed/count)*(1-failed/count)/count)),
        failure_mc_interval=failure_ci,failure_acceptance=equivalence_status(failure_ci,[0.,limits['eligible_failure_rate_max']]),
        calibration=calibration,person=person_summaries,
        rsm_pcm_comparison=test_summary if scenario.get('paired_rsm_pcm_comparison') else None,
        scientific_inference_ready=False,independent_reference_qualification='pending separate artifact review',
        scope='Pointwise Monte Carlo evidence at registered cells; not simultaneous qualification or public integration')
