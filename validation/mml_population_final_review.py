"""Read-only completion review of the sealed joint-MML population study.

Run with the frozen source root. No generation, optimizer, reference execution,
retry, new acceptance rule, or automatic scientific qualification is available.
"""
import argparse
from collections import Counter
import fcntl
import hashlib
import json
from pathlib import Path
import sys


TAM_GRIDS = ((12.0, .1, 'b12_h10'), (16.0, .05, 'b16_h05'))
REFERENCE_CHECKS = {'reference_nll', 'reference_gradient', 'reference_moments',
                    'relative_mass_error', 'tail_mass', 'tail_gradient'}
INFORMATION_CHECKS = {'reference_information', 'tail_information', 'numeric_information_allowance'}
FINITE_CHECKS = {'nll_per_person', 'gradient_per_person', 'information_per_person',
                 'hessian_step_per_person', 'posterior_moment', 'expansion', 'rowwise_R',
                 'native_rowwise', 'positive_definite', 'joint_covariance', 'expanded_se'}
CONTINUOUS_CHECKS = REFERENCE_CHECKS | {'r_nll', 'r_moments', 'r_sigma_score', 'r_mass',
                                      'r_tail', 'finite_nll', 'finite_gradient', 'finite_moments'}
TAM_CHECKS = REFERENCE_CHECKS | INFORMATION_CHECKS | {'A', 'B', 'response', 'person_order',
    'terminated', 'grid_nll', 'joint_score', 'fresh_moments', 'refit_point', 'refit_nll', 'expanded_se'}


def reference_inventory(scenario, replicate, plan):
    """Reconstruct preselected cases independently of delivered reference files."""
    if replicate not in set(scenario['independent_R_replicates'] + scenario['TAM_replicates']):
        raise ValueError('Unregistered reference replicate')
    if plan['numerical_reference']['TAM']['grids'] != [dict(bound=b, spacing=h) for b, h, _ in TAM_GRIDS]:
        raise ValueError('Unsupported registered TAM grid inventory')
    run = f"{scenario['id']}__{replicate:04d}"
    rows = []
    for model in scenario['fit_models']:
        for q in plan['estimator']['quadrature']:
            key = f'{run}_{model}_Q{q}'
            if replicate in scenario['independent_R_replicates']:
                rows.append((key, 'finite_R'))
                if q == plan['estimator']['primary_q']:
                    rows.append((key, 'continuous'))
            if (q == plan['estimator']['primary_q'] and replicate in scenario['TAM_replicates']
                    and model in plan['numerical_reference']['TAM']['models']):
                rows.append((key, 'TAM'))
        if model == scenario['primary_model'] and replicate in scenario['independent_R_replicates']:
            rows.append((f'{run}_{model}_truth', 'continuous'))
    return run, rows


def _read(path):
    return json.loads(Path(path).read_bytes())


def _verify_files(root, manifest):
    root = root.resolve()
    for name, wanted in manifest.items():
        path = Path(name)
        target = (root / path).resolve()
        if (path.is_absolute() or '..' in path.parts or not target.is_relative_to(root)
                or hashlib.sha256(target.read_bytes()).hexdigest() != wanted):
            raise ValueError('Retained reference artifact changed: ' + name)


def _recheck_r(bundle, selection, plan, ledger, retained_review):
    """Bind R files to inputs/events and repeat the sealed comparison only.

    No R process, fit or continuous integration is run. The truth-point
    continuous comparison repeats its original finite-rule value evaluation.
    """
    from validation.mml_population_reference_review import finite_review, continuous_review
    from validation.mml_population_study import json_bytes

    def same(actual, expected, label):
        if json_bytes(actual) != json_bytes(expected):
            raise ValueError(label + ' differs from retained reference evidence')

    key, lane = selection['id'], selection['lane']
    event = ledger.get(key, lane+'_R_finished')
    if event is None or event['status'] != 'completed':
        return False
    manifest = _read(bundle/'manifest.json')
    input_name = 'finite_input.json' if lane == 'finite_R' else 'continuous_input.json'
    required = {input_name, 'python_finite.json'} | ({'native_rowwise.json'} if lane == 'finite_R' else set())
    if not required.issubset(manifest):
        raise ValueError('Reference input manifest is incomplete')
    matches = [c for c in _read(bundle/input_name)['cases'] if c['id'] == key]
    if len(matches) != 1:
        raise ValueError('Required R input case is missing or duplicated')
    case = matches[0]
    request = dict(classification='DEVELOPMENT_ONLY', cases=[case])
    if lane == 'continuous':
        request['continuous_settings'] = plan['numerical_reference']['continuous_settings']
    output = bundle/'execution'/(key+'_'+lane+'_R')
    execution = _read(output/'execution.json')
    if not {'input.json', 'result.json', 'stdout.txt', 'stderr.txt'}.issubset(execution['files']):
        raise ValueError('R execution artifact inventory is incomplete')
    _verify_files(output, execution['files'])
    same(_read(output/'input.json'), request, 'R execution input')
    same(_read(output/'result.json')['cases'], [event['value']], 'R returned values')
    if event['value']['id'] != key:
        raise ValueError('R result belongs to another case')
    if type(execution['returncode']) is not int or execution['returncode'] != 0:
        return False
    python_rows = [c for c in _read(bundle/'python_finite.json') if c['id'] == key]
    if len(python_rows) != (0 if key.endswith('_truth') else 1):
        raise ValueError('Required Python comparison point is missing or duplicated')
    python = python_rows[0] if python_rows else None
    if lane == 'finite_R':
        native = [c for c in _read(bundle/'native_rowwise.json') if c['id'] == key]
        if len(native) != 1:
            raise ValueError('Required native comparison point is missing or duplicated')
        recomputed = finite_review(case, python, native[0], event['value'], plan)
    else:
        import numpy as np
        values = ledger.get(key, 'continuous_python_finished')
        if values is None or values['status'] != 'completed':
            return False
        # The sealed comparison originally receives ndarray gradients; JSON
        # retains them as lists. Restore the representation, not the values.
        restored = [dict(v, gradient=np.asarray(v['gradient'], dtype=float)) for v in values['value']]
        recomputed = continuous_review(case, restored, event['value']['continuous'], python, plan)
    same(recomputed, retained_review, 'R numerical comparison')
    return True


def review_reference(bundle, scenario, replicate, plan, identity):
    """Check case/grid denominators, event chains, retained checks and execution."""
    from validation.mml_population_study import AttemptLedger
    run, required = reference_inventory(scenario, replicate, plan)
    bundle = Path(bundle)
    expected_slots = [(key + '_' + grid, lane) for key, lane in required if lane == 'TAM'
                      for _, _, grid in TAM_GRIDS] + [(key, lane) for key, lane in required if lane != 'TAM']
    absent = [str(p.relative_to(bundle)) for p in (bundle/'manifest.json', bundle/'reference_plan.json',
              bundle/'execution/summary.json', bundle/'execution/events.sqlite') if not p.is_file()]
    if absent:
        return dict(run_id=run, required_slots=len(expected_slots), passed_slots=0,
                    missing_artifacts=absent, all_checks_pass=False, scientific_inference_ready=False)
    _verify_files(bundle, _read(bundle/'manifest.json'))
    selection = _read(bundle/'reference_plan.json')
    cases = selection['required_cases']
    if (selection['study_sha256'] != identity or selection['run_id'] != run
            or selection['scenario_id'] != scenario['id']
            or [(c['id'], c['lane']) for c in cases] != required
            or any(type(c['available']) is not bool for c in cases)):
        raise ValueError('Reference selection differs from the frozen plan')
    if selection['missing_required_cases'] != sum(not c['available'] for c in cases):
        raise ValueError('Reference missing-input denominator differs')
    summary = _read(bundle/'execution/summary.json')
    delivered = summary['comparisons']
    expected_delivered = [(c['id']+'_'+grid, c['lane']) for c in cases if c['lane']=='TAM' and c['available']
                          for _, _, grid in TAM_GRIDS]
    expected_delivered += [(c['id'], c['lane']) for c in cases if c['lane']!='TAM' or not c['available']]
    if (summary['run_id'] != run or summary['required_cases'] != len(required)
            or Counter((r['id'], r['lane']) for r in delivered) != Counter(expected_delivered)):
        raise ValueError('Reference comparison inventory differs from required cases/grids')
    records = {(r['id'], r['lane']): r for r in delivered}
    output = bundle/'execution'
    ledger = AttemptLedger.read_only(output/'events.sqlite', plan_sha256=identity)
    slots = []
    try:
        head = ledger.verify()
        if head != summary['ledger']:
            raise ValueError('Reference summary does not bind the current event chain')
        for case in cases:
            key, lane = case['id'], case['lane']
            for suffix in (['_'+grid for _, _, grid in TAM_GRIDS] if lane=='TAM' else ['']):
                slot = key+suffix
                row = records[(slot if case['available'] else key, lane)]
                if not case['available']:
                    if row['status']!='required_input_or_point_unavailable' or row['all_checks_pass'] is not False:
                        raise ValueError('Unavailable reference input was credited')
                    slots.append(dict(id=slot, lane=lane, status=row['status'], passed=False))
                    continue
                review = ledger.get(slot, lane+'_review_finished')
                completed = review is not None and review['status']=='completed'
                passed = False
                if completed:
                    value = review['value']; checks = value['checks']
                    expected = (FINITE_CHECKS if lane=='finite_R' else TAM_CHECKS if lane=='TAM' else
                                CONTINUOUS_CHECKS | (set() if key.endswith('_truth') else
                                                     INFORMATION_CHECKS | {'stationary','expanded_se'}))
                    if set(checks)!=expected or any(type(v) is not bool for v in checks.values()):
                        raise ValueError('Reference check inventory or boolean types differ')
                    passed = all(checks.values())
                    if value['all_checks_pass'] is not passed:
                        raise ValueError('Reference decision differs from its individual checks')
                if row['all_checks_pass'] is not passed or (completed and row['status']!='completed'):
                    raise ValueError('Reference summary differs from its retained review')
                execution_ok = None
                if lane in {'finite_R', 'continuous'} and completed:
                    execution_ok = _recheck_r(bundle, case, plan, ledger, value)
                if lane=='TAM' and completed:
                    fit = ledger.get(key, 'TAM_fit_finished')
                    execution_ok = bool(fit and fit['status']=='completed'
                                        and type(fit['value'].get('returncode')) is int
                                        and fit['value']['returncode']==0)
                    if fit and fit['status']=='completed':
                        required_files = {key+'_TAM_stdout.txt', key+'_TAM_stderr.txt',
                                          f'{key}_TAM/{key}_design.json', f'{key}_TAM/{slot}.json'}
                        if not required_files.issubset(fit['value']['files']):
                            raise ValueError('TAM execution artifact inventory is incomplete')
                        _verify_files(output, fit['value']['files'])
                slots.append(dict(id=slot, lane=lane, status=row['status'], passed=passed and execution_ok is not False,
                                  execution_success=execution_ok))
        native_pass = bool(delivered) and all(r['all_checks_pass'] is True for r in delivered)
        if summary['all_checks_pass'] is not native_pass:
            raise ValueError('Reference summary overall flag differs from its comparisons')
    finally:
        ledger.close()
    return dict(run_id=run, required_slots=len(expected_slots), passed_slots=sum(s['passed'] for s in slots),
                slots=slots, ledger=head, all_checks_pass=all(s['passed'] for s in slots),
                scientific_inference_ready=False)


def assemble_review(plan, summaries, references):
    """Preserve registered conditional criteria and all descriptive endpoints."""
    from validation.mml_population_analysis import equivalence_status, scenario_calibration_truth
    from validation.mml_population_missingness_review import review_comparison
    scenarios = plan['scenarios']
    if set(summaries) != {s['id'] for s in scenarios}:
        raise ValueError('Complete planned scenario inventory required')
    expected_refs = [f"{s['id']}__{r:04d}" for s in scenarios
                     for r in sorted(set(s['independent_R_replicates']+s['TAM_replicates']))]
    if [r['run_id'] for r in references] != expected_refs:
        raise ValueError('Complete planned reference inventory required')
    by_run = {r['run_id']: r for r in references}
    controls = []
    for scenario in scenarios:
        summary = summaries[scenario['id']]
        if summary['status']!='complete' or summary['planned']!=scenario['replicates']:
            raise ValueError('Incomplete scenario cannot enter final review')
        if scenario['cohort']=='negative_control':
            expected = [scenario['expected_status']] * (scenario['replicates'] * len(plan['estimator']['quadrature']))
            if len(summary['observed']) != len(expected):
                raise ValueError('Control observation denominator differs')
            passed = summary['observed']==expected
            if summary['expected']!=scenario['expected_status'] or summary['expected_rejection_pass'] is not passed:
                raise ValueError('Control outcome differs from registered rejection')
            controls.append(dict(scenario_id=scenario['id'], passed=passed, observed=summary['observed']))
    if len(controls)!=plan['negative_controls'] or not controls:
        raise ValueError('Registered negative-control denominator differs')
    decisions = []
    limits = plan['interpretation_limits']
    for scenario in scenarios:
        if scenario['cohort']=='negative_control':
            continue
        summary = summaries[scenario['id']]
        if summary['scenario_id']!=scenario['id']:
            raise ValueError('Scenario summary identity changed')
        truths = scenario_calibration_truth(scenario)
        rows = summary['calibration']
        if [(r['kind'], r['facet'], r['level']) for r in rows] != [(r['kind'], r['facet'], r['level']) for r in truths]:
            raise ValueError('Missing or changed calibration targets')
        checks = [('failure', equivalence_status(summary['failure_mc_interval'], [0., limits['eligible_failure_rate_max']]))]
        if checks[0][1] != summary['failure_acceptance']:
            raise ValueError('Failure acceptance differs from its registered interval')
        for row, truth in zip(rows, truths, strict=True):
            variable = bool(any(truth['jacobian'] != 0))
            if row['truth']!=truth['value'] or row['comparable'] is not variable:
                raise ValueError('Calibration truth or identification changed')
            if not variable:
                if row['status']!='fixed_by_identification':
                    raise ValueError('Fixed coordinates cannot earn recovery credit')
                continue
            width = (limits['positive_parameter_relative_bias']*truth['value']
                     if truth['kind'] in {'Slope','PopulationSD'} else limits['calibration_absolute_bias'])
            intervals = row['monte_carlo_intervals']
            acceptance = dict(bias=equivalence_status(intervals['bias'], [-width,width]),
                              se=equivalence_status(intervals['se_ratio'], limits['rms_model_se_over_empirical_sd']),
                              coverage=equivalence_status(intervals['coverage_given_interval'], limits['coverage_95']))
            if row['pointwise_acceptance']!=acceptance:
                raise ValueError('Calibration acceptance differs from registered intervals')
            label = '/'.join((row['kind'], row['facet'], row['level']))
            checks.extend((label+'/'+name, value) for name,value in acceptance.items())
        comparison = review_comparison(summary, scenario)
        if scenario['paired_rsm_pcm_comparison'] and scenario['primary_model']=='RSM':
            value = summary['rsm_pcm_comparison']
            status = equivalence_status(value['conditional_mc_interval'], limits['null_LRT_rejection'])
            if value['null_acceptance']!=status:
                raise ValueError('Conditional LRT criterion changed')
            checks.append(('conditional_null_LRT',status))
        expected_person = [(basis, q) for basis in ['fitted','known'] for q in ['all',1,2,3,4]]
        if [(p['calibration_basis'],p['truth_quartile']) for p in summary['person']] != expected_person:
            raise ValueError('Missing descriptive Person basis or truth-bin evidence')
        refs = [by_run[f"{scenario['id']}__{r:04d}"] for r in sorted(set(scenario['independent_R_replicates']+scenario['TAM_replicates']))]
        for replicate, record in zip(sorted(set(scenario['independent_R_replicates']+scenario['TAM_replicates'])), refs, strict=True):
            _, inventory = reference_inventory(scenario, replicate, plan)
            required = sum(len(TAM_GRIDS) if lane=='TAM' else 1 for _,lane in inventory)
            if (type(record['required_slots']) is not int or record['required_slots']!=required
                    or type(record['passed_slots']) is not int or not 0<=record['passed_slots']<=required
                    or record['all_checks_pass'] is not (record['passed_slots']==required)):
                raise ValueError('Reference slot denominator or aggregate decision differs')
        states = [state for _,state in checks]
        state = next((s for s in ['fail','unavailable','inconclusive'] if s in states), 'pointwise_criteria_satisfied')
        reference_ok = bool(refs) and all(r['all_checks_pass'] is True and r['passed_slots']==r['required_slots'] for r in refs)
        decisions.append(dict(scenario_id=scenario['id'], statistical_state=state,
            checks=[dict(endpoint=k,status=v) for k,v in checks], comparison_missingness= comparison,
            independent_references_complete_and_matching=reference_ok,
            registered_conditional_criteria_met=state=='pointwise_criteria_satisfied' and reference_ok and all(c['passed'] for c in controls),
            fitted_objective_worse_than_truth=summary['fitted_objective_worse_than_truth'],
            review_required=True,scientific_inference_ready=False))
    return dict(classification='COMPLETE_ENDPOINTS_PENDING_SCIENTIFIC_REVIEW', controls=controls,
        decisions=decisions,summaries=summaries,references=references,scientific_inference_ready=False,
        interpretation='Registered pointwise conditional criteria only. No pooled pass rate, simultaneous claim, posterior-SD sampling interpretation, or automatic release.')


def review_study(study, plan, identity):
    """Check all terminal inventories before reading any endpoint payload."""
    from validation.mml_population_study import AttemptLedger, json_bytes
    from validation.mml_population_analysis import summarize_scenario
    from validation.mml_population_missingness_review import review_calibration_coverage
    study = Path(study); progress = []; heads = {}; summaries = {}
    for scenario in plan['scenarios']:
        expected = {f"{scenario['id']}__{r:04d}" for r in range(1,scenario['replicates']+1)}
        path = study/'attempts'/(scenario['id']+'.sqlite')
        terminal = set()
        if path.is_file():
            ledger = AttemptLedger.read_only(path, plan_sha256=identity)
            try:
                existing = {r[0] for r in ledger.connection.execute('SELECT DISTINCT run_id FROM event')}
                if existing-expected:
                    raise ValueError('Unexpected attempt cannot be selected away')
                terminal = {r[0] for r in ledger.connection.execute("SELECT run_id FROM event WHERE phase='terminal'")}
            finally:
                ledger.close()
        progress.append(dict(scenario_id=scenario['id'],planned=len(expected),terminal=len(terminal)))
    if any(r['planned']!=r['terminal'] for r in progress):
        return dict(status='incomplete',progress=progress,endpoint_payloads_read=False,scientific_inference_ready=False)
    for scenario in plan['scenarios']:
        key = scenario['id']; ledger = AttemptLedger.read_only(study/'attempts'/(key+'.sqlite'),plan_sha256=identity)
        try:
            recomputed = summarize_scenario(ledger,scenario,plan)
            if json_bytes(recomputed)!=json_bytes(_read(study/'summaries'/(key+'.json'))):
                raise ValueError('Stored summary differs from sealed-ledger reaggregation: '+key)
            summaries[key] = recomputed; heads[key] = ledger.verify()
        finally:
            ledger.close()
    references = []
    for scenario in plan['scenarios']:
        for replicate in sorted(set(scenario['independent_R_replicates']+scenario['TAM_replicates'])):
            references.append(review_reference(study/'references'/f"{scenario['id']}__{replicate:04d}",
                                               scenario,replicate,plan,identity))
    return dict(assemble_review(plan,summaries,references),study_sha256=identity,scenario_ledgers=heads,
                registered_analysis_recomputed=True,
                calibration_coverage_interpretation=[review_calibration_coverage(summaries[s['id']], s)
                    for s in plan['scenarios'] if s['cohort']!='negative_control'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('study',type=Path); parser.add_argument('--source-root',required=True,type=Path)
    parser.add_argument('--output',required=True,type=Path); args = parser.parse_args()
    sys.path.insert(0,str(args.source_root.resolve()))
    from validation import mml_population_controller as controller
    from validation.mml_population_study import json_bytes
    if controller.ROOT != args.source_root.resolve():
        raise ValueError('Controller did not load from the requested frozen source')
    # The previously written missingness interpretation is supplemental; it is
    # absent from the original MML seal. Frozen modules retain search priority.
    import validation
    validation.__path__.append(str(Path(__file__).resolve().parent))
    with (args.study/'controller.lock').open('rb') as lock:
        fcntl.flock(lock,fcntl.LOCK_SH|fcntl.LOCK_NB)
        plan,identity = controller.verify_seal(args.study)
        verify = controller.environment_verifier(args.study,plan)
        packet = review_study(args.study,plan,identity)
        verify()
    packet['review_source_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    packet['supplemental_source_sha256'] = {
        name: hashlib.sha256((Path(__file__).parent/name).read_bytes()).hexdigest()
        for name in ['mml_population_missingness_review.py','bias_operating_metrics.py']}
    packet['frozen_seal_and_runtime_verified'] = True
    with args.output.open('xb') as handle:
        handle.write(json_bytes(packet)+b'\n')
    print(json.dumps(dict(output=str(args.output),scientific_inference_ready=False)))


if __name__=='__main__':
    main()
