"""Read-only G3 evidence review using the sealed analysis and reference rules.

No generation, fit, R execution, retry, acceptance change or public promotion.
All planned attempts must be terminal before any population endpoint is read.
"""
import argparse
from collections import Counter
import fcntl
import hashlib
import json
from pathlib import Path
import sys


def read(path):
    return json.loads(Path(path).read_bytes())


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify_files(root, manifest):
    root = Path(root).resolve()
    for name, expected in manifest.items():
        path = (root / name).resolve()
        if not path.is_relative_to(root) or digest(path) != expected:
            raise ValueError('Retained artifact changed: ' + name)


def require_equal(actual, expected, label):
    from validation.mml_population_study import json_bytes
    if json_bytes(actual) != json_bytes(expected):
        raise ValueError(label + ' differs from retained evidence')


def check_r_result(root, value):
    """Confirm retained R execution, without executing a process."""
    files = value['files']
    if not {'stdout.txt', 'stderr.txt', 'execution.json', 'result.json'}.issubset(files):
        raise ValueError('Missing R execution artifacts')
    verify_files(root, files)
    execution = read(root/'execution.json')
    if type(execution['returncode']) is not int or execution['returncode'] != 0:
        raise ValueError('Completed R reference has a nonzero execution status')
    require_equal(execution['files'], {k: v for k, v in files.items() if k != 'execution.json'}, 'R file inventory')
    require_equal(read(root/'result.json'), value['result'], 'R returned values')


def review_controls(root, plan, identity):
    """Recheck all fixed native/R controls; these are not MC repetitions."""
    from validation.bias_negative_controls import scalar_controls
    from validation.bias_population_check import compare_reference
    from validation.bias_population_references import _decode, review_family
    from validation.mml_population_study import AttemptLedger
    root = Path(root)
    definition = plan['deterministic_controls']
    scalar_ids = [c['id'] for c in scalar_controls()]
    if len(scalar_ids) != definition['scalar_cases'] or definition['family_cases'] != ['RSM', 'PCM', 'GPCM']:
        raise ValueError('Unsupported deterministic control inventory')
    ledger = AttemptLedger.read_only(root/'events.sqlite', plan_sha256=identity)
    groups = []
    try:
        head = ledger.verify()
        for key, native_n, reference_n in [('scalar', len(scalar_ids), len(scalar_ids))] + [
                (m, 1, 10) for m in definition['family_cases']]:
            record = ledger.get(key, 'control_finished')
            if record is None:
                raise ValueError('Missing completed control event: ' + key)
            value = record.get('value', {})
            native, refs = value.get('native', []), value.get('references', [])
            verify_files(root, value.get('files', {}))
            if record['status'] == 'completed':
                target = root/key
                rfiles = {str(Path(name).relative_to(key+'/R')): wanted
                          for name, wanted in value['files'].items() if Path(name).is_relative_to(key+'/R')}
                reference = read(target/'R/result.json')
                check_r_result(target/'R', dict(result=reference, files=rfiles))
                retained = _decode(read(target/'native.json'))
                if key == 'scalar':
                    if [r['id'] for r in retained] != scalar_ids:
                        raise ValueError('Native scalar-control inventory changed')
                    for row in retained:
                        if any(type(v) is not bool for v in row['checks'].values()) or row['all_checks_pass'] is not all(row['checks'].values()):
                            raise ValueError('Native scalar-control decision changed')
                    recomputed_native = [dict(id=r['id'], pass_=r['all_checks_pass']) for r in retained]
                    recomputed_refs, _ = compare_reference(retained, [], reference, plan['independent_reference']['limits'])
                else:
                    passed = retained['actual'] == retained['expected']
                    if retained['native_status_pass'] is not passed:
                        raise ValueError('Native family-control decision changed')
                    recomputed_native = [dict(id=key, pass_=passed)]
                    result = review_family(_decode(read(target/'python.json')), reference, plan)
                    recomputed_refs = result['cells'] + result['adjustments'] + result['derived_adjustments']
                require_equal(native, recomputed_native, 'Native control summary')
                require_equal(refs, recomputed_refs, 'R control comparison')
            if any(type(r['pass_']) is not bool for r in native+refs):
                raise ValueError('Control pass flags must be boolean')
            groups.append(dict(id=key, status=record['status'], required_native=native_n, required_reference=reference_n,
                passed_native=sum(r['pass_'] for r in native), passed_reference=sum(r['pass_'] for r in refs),
                all_checks_pass=len(native)==native_n and len(refs)==reference_n and all(r['pass_'] for r in native+refs)))
        summary = dict(classification='DETERMINISTIC_FORMULA_AND_GUARD_CONTROLS', groups=groups,
            required_native=sum(r['required_native'] for r in groups), required_reference=sum(r['required_reference'] for r in groups),
            passed_native=sum(r['passed_native'] for r in groups), passed_reference=sum(r['passed_reference'] for r in groups),
            all_checks_pass=all(r['all_checks_pass'] for r in groups), ledger=head, scientific_inference_ready=False,
            interpretation='Fixed controls, not independent simulation replicates or LR-distribution qualification')
        require_equal(read(root/'summary.json'), summary, 'Control summary')
        return summary
    finally:
        ledger.close()


def review_reference(bundle, scenario, replicate, plan, identity, export):
    from validation.bias_population_references import _decode, review_family
    from validation.mml_population_study import AttemptLedger
    bundle = Path(bundle); run = f"{scenario['id']}__{replicate:04d}"
    if replicate not in scenario['independent_reference_replicates']:
        raise ValueError('Unregistered reference replicate')
    ncells = len(scenario['design']['rater'])*len(scenario['design']['task'])
    required = len(plan['scope']['lanes'])*(ncells+4)
    if export is None or export['status'] != 'completed':
        return dict(run_id=run, status='missing' if export is None else export['status'], required_slots=required,
                    passed_slots=0, all_checks_pass=False, stage='export')
    if digest(bundle/'manifest.json') != export['value']['manifest_sha256']:
        raise ValueError('Reference export manifest changed')
    verify_files(bundle, read(bundle/'manifest.json'))
    selection = read(bundle/'selection.json')
    if (selection['study_sha256'] != identity or selection['scenario_id'] != scenario['id']
            or selection['run_id'] != run or selection['replicate'] != replicate
            or [f['basis'] for f in selection['families']] != plan['scope']['lanes']):
        raise ValueError('Required reference identity or lanes changed')
    output = bundle/'execution'
    ledger = AttemptLedger.read_only(output/'events.sqlite', plan_sha256=identity)
    families = []
    try:
        head = ledger.verify()
        for family in selection['families']:
            if ([family[k] for k in ['required_scalar_cells','required_adjustments','required_derived_adjustments']]
                    != [ncells,2,2] or type(family['available']) is not bool):
                raise ValueError('Reference family denominator or availability changed')
            basis = family['basis']; key = run+'_'+basis
            if not family['available']:
                families.append(dict(basis=basis, status='required_input_unavailable', required_slots=ncells+4,
                    passed_slots=0, all_checks_pass=False, input_status=family['input_status']))
                continue
            ref = ledger.get(key, 'R_finished')
            if ref is None:
                # An interrupted attempt is retained as uncompleted, never rerun.
                ref = dict(status='interrupted_before_result' if ledger.get(key,'R_started') is not None else 'missing')
            checked = dict(status=ref['status'])
            if ref['status'] == 'completed':
                check_r_result(output/basis, ref['value'])
                checked = ledger.get(key, 'review_finished') or dict(status='missing')
                if checked['status'] == 'completed':
                    result = review_family(_decode(read(bundle/basis/'python.json')), ref['value']['result'], plan)
                    require_equal(checked['value'], result, 'Reference comparison')
            value = checked.get('value', {})
            rows = value.get('cells', [])+value.get('adjustments', [])+value.get('derived_adjustments', [])
            if any(type(r['pass_']) is not bool for r in rows):
                raise ValueError('Reference pass flags must be boolean')
            families.append(dict(basis=basis, status=checked['status'], required_slots=ncells+4,
                passed_slots=sum(r['pass_'] for r in rows), all_checks_pass=bool(len(rows)==ncells+4 and value.get('all_checks_pass')),
                reference_status=ref['status']))
        result = dict(run_id=run, families=families, status_counts=dict(Counter(r['status'] for r in families)),
            required_slots=required, passed_slots=sum(r['passed_slots'] for r in families),
            all_checks_pass=bool(families) and all(r['all_checks_pass'] for r in families), ledger=head, scientific_inference_ready=False)
        if (output/'summary.json').is_file():
            require_equal(read(output/'summary.json'), result, 'Reference summary')
        return result
    finally:
        ledger.close()


def review_study(study, plan, identity):
    from validation.bias_population_analysis import assemble_evidence, summarize_bias_scenario
    from validation.mml_population_study import AttemptLedger
    study = Path(study); progress = []; summaries = []; references = []
    for scenario in plan['scenarios']:
        expected = {f"{scenario['id']}__{r:04d}" for r in range(1,scenario['replicates']+1)}
        terminal = set(); path = study/'attempts'/(scenario['id']+'.sqlite')
        if path.is_file():
            ledger = AttemptLedger.read_only(path, plan_sha256=identity)
            try:
                existing = {r[0] for r in ledger.connection.execute('SELECT DISTINCT run_id FROM event')}
                if existing-expected: raise ValueError('Unexpected RunId cannot be omitted')
                terminal = {r[0] for r in ledger.connection.execute("SELECT run_id FROM event WHERE phase='terminal'")}
            finally:
                ledger.close()
        progress.append(dict(scenario_id=scenario['id'], planned=len(expected), terminal=len(terminal)))
    if any(r['planned'] != r['terminal'] for r in progress):
        return dict(status='incomplete', progress=progress, endpoint_payloads_read=False, scientific_inference_ready=False)
    controls = review_controls(study/'controls', plan, identity)
    exports = AttemptLedger.read_only(study/'reference_exports.sqlite', plan_sha256=identity)
    try:
        export_head = exports.verify()
        for scenario in plan['scenarios']:
            key = scenario['id']
            ledger = AttemptLedger.read_only(study/'attempts'/(key+'.sqlite'), plan_sha256=identity)
            try:
                summary = summarize_bias_scenario(ledger, scenario, plan)
                summary.update(study_sha256=identity, source_ledger=ledger.verify())
                require_equal(read(study/'summaries'/(key+'.json')), summary, 'Scenario summary: '+key)
                summaries.append(summary)
            finally:
                ledger.close()
            for replicate in scenario['independent_reference_replicates']:
                run = f'{key}__{replicate:04d}'
                references.append(review_reference(study/'references'/run, scenario, replicate, plan, identity,
                                                   exports.get(run,'export_finished')))
    finally:
        exports.close()
    packet = dict(assemble_evidence(plan, summaries, references, controls), study_sha256=identity)
    stored = study/'evidence_packet.json'
    if stored.is_file():
        require_equal(read(stored), packet, 'Final evidence packet')
    return dict(packet, registered_analysis_recomputed=True, retained_reference_comparisons_rechecked=True,
                reference_export_ledger=export_head, stored_packet_matches=True if stored.is_file() else None)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('study',type=Path); parser.add_argument('--source-root',required=True,type=Path)
    parser.add_argument('--output',required=True,type=Path); args = parser.parse_args()
    sys.path.insert(0,str(args.source_root.resolve()))
    from validation import bias_population_controller as controller
    from validation.mml_population_study import json_bytes
    if controller.ROOT != args.source_root.resolve():
        raise ValueError('Controller did not load from the requested frozen source')
    with (args.study/'controller.lock').open('rb') as lock:
        fcntl.flock(lock,fcntl.LOCK_SH|fcntl.LOCK_NB)
        plan,identity = controller.verify_seal(args.study)
        verify = controller.environment_verifier(args.study,plan)
        packet = review_study(args.study,plan,identity)
        verify()
    packet.update(review_source_sha256=digest(Path(__file__)), frozen_seal_and_runtime_verified=True)
    with args.output.open('xb') as handle:
        handle.write(json_bytes(packet)+b'\n')
    print(json.dumps(dict(output=str(args.output),scientific_inference_ready=False)))


if __name__ == '__main__':
    main()
