"""Single-attempt population-study execution over retained raw input bytes.

The frozen protocol, execution receipt and source/runtime checks are handled
by the command-line controller. This module never generates a plan or silently
retries an interrupted optimizer. SQLite transactions retain compressed inputs
and results without writing large posterior grids that can be reconstructed.
"""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sqlite3
import zlib

import numpy as np
from scipy.stats import chi2

from validation.mfrm_recovery_design import RecoveryDesign, generate_population_sample, replicate_seed
from validation.mfrm_recovery_metrics import calibration_coordinates, candidate_wald_intervals, population_truth_quartile
from validation.mml_free_sd_stationarity_adapter import prepare_app_free_sd_problem_from_data
from validation.mml_free_sd_prefit_transport import _configured_data
from validation.mml_gpcm_joint_check import prepare_problem
from validation.mml_grouped_fit import fit_grouped_joint, normal_rule
from validation.mml_grouped_likelihood import GroupedMmlLikelihood


def json_bytes(value):
    """Nonfinite raw failure evidence is tagged, never changed to a number."""
    def encode(x):
        if isinstance(x, np.ndarray): return encode(x.tolist())
        if isinstance(x, np.generic): return encode(x.item())
        if isinstance(x, float) and not np.isfinite(x): return {'nonfinite': str(x)}
        if isinstance(x, dict): return {k: encode(v) for k, v in x.items()}
        if isinstance(x, (tuple, list)): return [encode(v) for v in x]
        return x
    return json.dumps(encode(value), ensure_ascii=False, sort_keys=True, separators=(',', ':'), allow_nan=False).encode('utf-8')


class AttemptLedger:
    """Append-only events. A unique RunId/phase is durable before its work."""
    def __init__(self, path, *, plan_sha256):
        self.connection = sqlite3.connect(path, timeout=30)
        self.connection.execute('PRAGMA synchronous=FULL')
        self.connection.executescript('''
            CREATE TABLE IF NOT EXISTS binding (plan_sha256 TEXT PRIMARY KEY);
            CREATE TABLE IF NOT EXISTS event (
                sequence INTEGER PRIMARY KEY, run_id TEXT NOT NULL, phase TEXT NOT NULL,
                utc TEXT NOT NULL, payload BLOB NOT NULL, raw_sha256 TEXT NOT NULL,
                parent_sha256 TEXT NOT NULL, event_sha256 TEXT NOT NULL UNIQUE,
                UNIQUE(run_id, phase));
            CREATE TRIGGER IF NOT EXISTS no_event_update BEFORE UPDATE ON event BEGIN SELECT RAISE(ABORT,'append only'); END;
            CREATE TRIGGER IF NOT EXISTS no_event_delete BEFORE DELETE ON event BEGIN SELECT RAISE(ABORT,'append only'); END;
            CREATE TRIGGER IF NOT EXISTS no_binding_update BEFORE UPDATE ON binding BEGIN SELECT RAISE(ABORT,'immutable binding'); END;
            CREATE TRIGGER IF NOT EXISTS no_binding_delete BEFORE DELETE ON binding BEGIN SELECT RAISE(ABORT,'immutable binding'); END;
        ''')
        with self.connection:
            bindings = self.connection.execute('SELECT plan_sha256 FROM binding').fetchall()
            if not bindings: self.connection.execute('INSERT INTO binding VALUES (?)', (plan_sha256,))
            elif bindings != [(plan_sha256,)]: raise ValueError('Attempt ledger belongs to another frozen plan')

    def append(self, run_id, phase, value):
        raw = json_bytes(value); digest = hashlib.sha256(raw).hexdigest()
        utc = datetime.now(timezone.utc).isoformat()
        with self.connection:
            self.connection.execute('BEGIN IMMEDIATE')
            row = self.connection.execute('SELECT event_sha256 FROM event ORDER BY sequence DESC LIMIT 1').fetchone()
            parent = row[0] if row else '0'*64
            event_hash = hashlib.sha256(json_bytes([run_id, phase, utc, digest, parent])).hexdigest()
            self.connection.execute('INSERT INTO event(run_id,phase,utc,payload,raw_sha256,parent_sha256,event_sha256) VALUES(?,?,?,?,?,?,?)',
                (run_id, phase, utc, zlib.compress(raw), digest, parent, event_hash))

    def get(self, run_id, phase):
        row = self.connection.execute('SELECT payload,raw_sha256 FROM event WHERE run_id=? AND phase=?', (run_id, phase)).fetchone()
        if row is None: return None
        raw = zlib.decompress(row[0])
        if hashlib.sha256(raw).hexdigest() != row[1]: raise ValueError('Retained event payload changed')
        return json.loads(raw)

    def started(self, run_id):
        return self.connection.execute('SELECT 1 FROM event WHERE run_id=? LIMIT 1', (run_id,)).fetchone() is not None

    def verify(self):
        parent = '0'*64; count = 0
        for run_id, phase, utc, payload, digest, old_parent, event_hash in self.connection.execute(
                'SELECT run_id,phase,utc,payload,raw_sha256,parent_sha256,event_sha256 FROM event ORDER BY sequence'):
            if old_parent != parent or hashlib.sha256(zlib.decompress(payload)).hexdigest() != digest:
                raise ValueError('Attempt ledger payload or parent changed')
            if hashlib.sha256(json_bytes([run_id, phase, utc, digest, parent])).hexdigest() != event_hash:
                raise ValueError('Attempt ledger event identity changed')
            parent = event_hash; count += 1
        return dict(events=count, head_sha256=parent)

    def close(self): self.connection.close()

    @classmethod
    def read_only(cls, path, *, plan_sha256):
        ledger = cls.__new__(cls)
        ledger.connection = sqlite3.connect(Path(path).resolve().as_uri()+'?mode=ro', uri=True)
        if ledger.connection.execute('SELECT plan_sha256 FROM binding').fetchall() != [(plan_sha256,)]:
            ledger.close(); raise ValueError('Attempt ledger belongs to another frozen plan')
        return ledger


def _error(exc): return dict(type=type(exc).__name__, message=str(exc))


def _mutate(data, scenario):
    data = data.copy(); mutation = scenario.get('mutation')
    if mutation == 'all_minimum': data['Score'] = 0
    elif mutation == 'all_maximum': data['Score'] = len(scenario['design']['steps'][0])
    elif mutation == 'unused_internal_category': data.loc[data.Score == 1, 'Score'] = 0
    elif mutation == 'missing_declared_rater': data = data.loc[data.Rater != 'R2'].reset_index(drop=True)
    elif mutation == 'fractional_weight': data['Weight'] = 1.5
    elif mutation == 'aliased_rater_task':
        data = data.loc[data.Task == data.Rater.map(lambda r: 'T'+str(int(r[1:]) % 2))].reset_index(drop=True)
    elif mutation is not None: raise ValueError('Unregistered negative-control mutation')
    return data


def prepare_retained_problem(retained, model, q):
    raw = retained['csv_utf8'].encode('utf-8')
    if hashlib.sha256(raw).hexdigest() != retained['csv_sha256']: raise ValueError('Raw input digest changed')
    cfg = retained['metadata']['design']; rating_max = cfg['rating_min']+len(cfg['steps'][0])
    data = _configured_data(raw, dict(person_col='Person', facet_cols=['Rater','Task','Criterion'],
        score_col='Score', rating_min=cfg['rating_min'], rating_max=rating_max, model=model,
        noncenter_facet='Criterion', step_facet=None if model=='RSM' else 'Criterion', positive_facets=[],
        weight_col='Weight' if retained.get('mutation')=='fractional_weight' else None, quadrature_points=q, sigma_start=1.))
    for facet, levels in retained['declared_levels'].items():
        if set(data[facet]) != set(levels): return data, None, 'input_missing_declared_level'
    if 'Weight' in data and not (data.Weight == 1).all(): return data, None, 'input_nonunit_weight'
    # Preserve the native G0 guard, which runs before problem construction.
    # In particular, all-minimum/all-maximum data have an improving common
    # Criterion-location ray. Do not misreport this guard as an optimizer error
    # or imply that the later LP recession audit was actually executed.
    if data.Score.nunique() < 2: return data, None, 'input_single_observed_category'
    if model == 'GPCM':
        _, problem = prepare_problem(data, q, rating_min=cfg['rating_min'], rating_max=rating_max)
    else:
        _, problem = prepare_app_free_sd_problem_from_data(data, person_col='Person', facet_cols=['Rater','Task','Criterion'],
            score_col='Score', rating_min=cfg['rating_min'], rating_max=rating_max, model=model,
            noncenter_facet='Criterion', step_facet=None if model=='RSM' else 'Criterion',
            weight_col='Weight' if 'Weight' in data else None, quadrature_points=q)
    return data, problem, None


def compact_fit(fit):
    fit = dict(fit)
    if 'evaluated' in fit:
        fit['evaluated'] = {k:v for k,v in fit['evaluated'].items() if k != 'posterior'}
        fit['posterior_storage'] = 'Reconstruct from retained CSV, point and standardized quadrature; no posterior grid stored'
    return fit


def truth_joint(config, truth):
    blocks = [np.asarray(truth[f])[:len(config['facet_levels'][f])-int(f != config['noncenter_facet'])]
        for f in config['facet_names']]
    blocks.extend(np.asarray(truth['Step'])[:, :-1])
    if config['model']=='GPCM': blocks.append(np.asarray(truth['LogSlope'])[:-1])
    blocks.append([np.log(truth['PopulationSD'])])
    return np.concatenate(blocks)


def compare_orders(data, problem, fits, plan):
    if len(fits) != 2 or not all(f.get('numerical_pass',False) for f in fits):
        return dict(pass_=False, status='one_or_both_finite_rule_fits_ineligible')
    low, high = fits; points = [np.asarray(f['point']) for f in fits]
    cov = [np.asarray(f['covariance']) for f in fits]
    scale = np.sqrt(np.outer(np.diag(cov[1]), np.diag(cov[1])))
    grouped = GroupedMmlLikelihood(data, problem.config)
    cross = [grouped.evaluate(points[1-i], *normal_rule(plan['estimator']['quadrature'][i])) for i in range(2)]
    expanded_se=[]
    for point,matrix in zip(points,cov,strict=True):
        jacobian=np.array([r['jacobian'] for r in calibration_coordinates(problem.config,point)])
        variable=np.any(jacobian!=0,axis=1)
        variance=np.einsum('ij,jk,ik->i',jacobian[variable],matrix,jacobian[variable])
        if np.any(variance<=0): raise FloatingPointError('Nonpositive variance of an estimated expanded coordinate')
        expanded_se.append(np.sqrt(variance))
    errors = dict(coordinate=float(np.max(abs(points[0]-points[1]))),
        cross_score_supnorm=max(float(np.max(abs(v['gradient']))) for v in cross),
        nll_per_person=max(abs(cross[i]['nll']-fits[1-i]['evaluated']['nll']) for i in range(2))/problem.config['n_person'],
        posterior_moment=max(float(np.max(abs(np.asarray(low['evaluated'][key])-high['evaluated'][key]))) for key in ['eap','sd']),
        covariance_scaled=float(np.max(abs(cov[0]-cov[1])/scale)),
        expanded_se_relative=float(np.max(abs(expanded_se[0]/expanded_se[1]-1))))
    limits = plan['estimator']['quadrature_limits']; checks = {k:v <= limits[k] for k,v in errors.items()}
    return dict(pass_=all(checks.values()), status='compared', errors=errors, checks=checks)


def primary_endpoints(retained, problem, fit, eligible, *, q, check_q, limits):
    truth = retained['truth']; config = problem.config
    point = truth_joint(config, truth)
    # Used solely for scoring/recovery; never passed to the optimizer.
    true_rows = calibration_coordinates(config, point)
    fitted_rows = calibration_coordinates(config, fit['point']) if fit.get('point') is not None else None
    intervals = candidate_wald_intervals(fitted_rows, fit['covariance']) if eligible else None
    rows = []
    for j, actual in enumerate(true_rows):
        estimated = fitted_rows[j] if fitted_rows else None; interval = intervals[j] if intervals else None
        comparable = bool(np.any(actual['jacobian'] != 0))
        rows.append(dict(kind=actual['kind'], facet=actual['facet'], level=actual['level'], truth=actual['value'],
            estimate=estimated['value'] if estimated else None, comparable=comparable,
            eligible=bool(eligible and comparable), se=interval['se'] if interval else None,
            lower=interval['lower'] if interval else None, upper=interval['upper'] if interval else None,
            interval_basis=interval['interval_basis'] if interval else 'unavailable',
            covered=bool(interval['lower'] <= actual['value'] <= interval['upper']) if interval and comparable else None))
    data, _, invalid = prepare_retained_problem(retained, config['model'], q)
    if invalid: raise ValueError('Retained input changed before truth scoring')
    grouped = GroupedMmlLikelihood(data, config)
    oracle, oracle_check = [grouped.evaluate(point, *normal_rule(order)) for order in [q, check_q]]
    truth_q_errors = dict(nll_per_person=abs(oracle['nll']-oracle_check['nll'])/config['n_person'],
        cross_score_supnorm=float(np.max(abs(oracle['gradient']-oracle_check['gradient']))),
        posterior_moment=max(float(np.max(abs(oracle[key]-oracle_check[key]))) for key in ['eap','sd']))
    theta = np.asarray(truth['Person']); estimated = fit.get('evaluated', {})
    persons = dict(truth=theta, quartile=population_truth_quartile(theta, truth['PopulationSD']),
        fitted_eap=estimated.get('eap'), fitted_posterior_sd=estimated.get('sd'), eligible=bool(eligible),
        known_calibration_eap=oracle['eap'], known_calibration_posterior_sd=oracle['sd'],
        known_calibration_q_errors=truth_q_errors,
        known_calibration_numerically_eligible=all(v <= limits[k] for k,v in truth_q_errors.items()),
        interval_basis='normal approximation using calibration-fixed posterior SD; not a sampling SE')
    return dict(calibration=rows, person=persons,
        fitted_nll_minus_truth_nll=fit['evaluated']['nll']-oracle['nll'] if fit.get('evaluated') else None,
        scientific_inference_ready=False)


def execute_registered_run(ledger, run, scenario, plan, *, verify_environment):
    """Controller supplies only a source/runtime verifier, never an evaluator."""
    run_id = run['run_id']
    registered = next((s for s in plan['scenarios'] if s['id']==scenario['id']), None)
    replicate = run['replicate']
    if registered is None or json_bytes(registered) != json_bytes(scenario) or type(replicate) is not int or not 1 <= replicate <= scenario['replicates']:
        raise ValueError('Scenario or replicate differs from the frozen plan')
    expected = dict(run_id=f"{scenario['id']}__{replicate:04d}", scenario_id=scenario['id'], replicate=replicate,
        seed=replicate_seed(plan['master_seed'], plan['study_id']+':'+scenario['id'], replicate), cohort=scenario['cohort'])
    if type(run['seed']) is not int or run != expected: raise ValueError('RunId, seed or cohort differs from the frozen plan')
    if ledger.started(run_id): raise ValueError('Registered attempt already started; numerical replay is prohibited')
    implementation = verify_environment()
    ledger.append(run_id, 'generation_started', run)
    if implementation is not None: ledger.append(run_id, 'execution_implementation', implementation)
    try:
        spec = RecoveryDesign(**scenario['design']); sample = generate_population_sample(spec, run['seed'])
        data = _mutate(sample['data'], scenario); csv = data.to_csv(index=False)
        retained = {k:sample[k] for k in ['truth','declared_levels','metadata']}
        retained.update(csv_utf8=csv, csv_sha256=hashlib.sha256(csv.encode('utf-8')).hexdigest(), seed=run['seed'],
            scenario_id=scenario['id'], mutation=scenario.get('mutation'), rows_after_mutation=len(data))
        ledger.append(run_id, 'input_retained', retained)
    except Exception as exc:
        ledger.append(run_id, 'terminal', dict(status='generation_failed', error=_error(exc), truth_available=False))
        return
    # Read the retained bytes back before constructing any fitting problem.
    retained = ledger.get(run_id, 'input_retained')
    if (retained['seed'] != run['seed'] or retained['scenario_id'] != scenario['id']
        or json_bytes(retained['metadata']['design']) != json_bytes(scenario['design'])):
        raise ValueError('Retained generation metadata differs from the frozen run')
    results = {}
    for model in scenario['fit_models']:
        fits = []; problem = None
        for q in plan['estimator']['quadrature']:
            phase = f'{model}_Q{q}'
            try:
                verify_environment()
                data, problem, invalid = prepare_retained_problem(retained, model, q)
                authorization = dict(input_sha256=retained['csv_sha256'], model=model, q=q, input_status=invalid or 'valid',
                    structural_start=problem.structural_start if problem else None, sigma_start=problem.sigma_start if problem else None,
                    seed=run['seed'])
                ledger.append(run_id, phase+'_started', authorization)
                result = dict(status=invalid, numerical_pass=False, point=None, covariance=None) if invalid else fit_grouped_joint(data, problem, q=q)
                result = compact_fit(result)
                ledger.append(run_id, phase+'_finished', result); fits.append(result)
            except Exception as exc:
                failure = dict(status='execution_failed', numerical_pass=False, error=_error(exc))
                ledger.append(run_id, phase+'_failed', failure); fits.append(failure)
        try:
            check = compare_orders(data, problem, fits, plan) if problem else dict(pass_=False, status='input_invalid')
            endpoints = primary_endpoints(retained, problem, fits[0], check['pass_'], q=plan['estimator']['primary_q'],
                check_q=plan['estimator']['check_q'], limits=plan['estimator']['quadrature_limits']) \
                if model==scenario['primary_model'] and problem is not None and scenario['cohort']=='population_recovery' else None
            results[model] = dict(status='eligible_candidate' if check['pass_'] else 'ineligible_candidate',
                q_comparison=check, endpoints=endpoints, primary_fit_status=fits[0]['status'],
                primary_nll=fits[0].get('evaluated',{}).get('nll'), K=len(fits[0]['point']) if fits[0].get('point') is not None else None)
        except Exception as exc:
            results[model] = dict(status='extraction_failed', error=_error(exc))
        ledger.append(run_id, model+'_assessment', results[model])
    comparison = None
    if scenario.get('paired_rsm_pcm_comparison'):
        eligible = all(results[m]['status']=='eligible_candidate' for m in ['RSM','PCM'])
        df = (len(spec.criterion)-1)*(len(spec.steps[0])-1)
        lr = 2*(results['RSM']['primary_nll']-results['PCM']['primary_nll']) if eligible else None
        # Negative LR beyond roundoff is an optimizer/model-consistency failure.
        valid = bool(eligible and lr >= -1e-8 and results['PCM']['K']-results['RSM']['K']==df)
        comparison = dict(available=valid, df=df, raw_lr=lr, p_value=float(chi2.sf(max(lr,0.),df)) if valid else None,
            status='candidate_chi_square_LRT' if valid else 'unavailable', scientific_inference_ready=False)
        # N is the number of independent Person integrals, not response rows.
        # These criteria are retained descriptively; no public recommendation.
        comparison['information_criteria']={}
        for model in ['RSM','PCM']:
            value=results[model];nll=value.get('primary_nll');k=value.get('K')
            comparison['information_criteria'][model]=dict(K=k,n_person=spec.persons,nll=nll,
                AIC=2*nll+2*k if nll is not None and k is not None else None,
                BIC=2*nll+k*np.log(spec.persons) if nll is not None and k is not None else None,
                numerically_eligible=value['status']=='eligible_candidate',interpretation='descriptive candidate only')
    ledger.append(run_id, 'terminal', dict(status='completed', truth_available=True, model_status={m:r['status'] for m,r in results.items()},
        comparison=comparison, negative_control_expected=scenario.get('expected_status'), scientific_inference_ready=False))
