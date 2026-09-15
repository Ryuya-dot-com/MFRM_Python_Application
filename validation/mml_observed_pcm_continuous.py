#!/usr/bin/env python3
"""Continuous-likelihood local optimum for the retained weighted/missing PCM case."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import platform
import subprocess
from types import SimpleNamespace

import numpy as np
import scipy
from scipy.optimize import least_squares

from mml_observed_pcm_refit import ROOT, ObservedPCMReference, discrepancies
from mml_pcm_continuous_refit import dump, sha
from mfrm_app.mml_stationarity import JointPolishOptions, polish_joint_free_sd, information_diagnostics, audit_joint_gradient

PREVIOUS = ROOT/'validation/generated/mml_observed_pcm_refit_20260914'
PREVIOUS_HASH = '9204597c2128374f36793387d3cd9feaab02988a7b89c66fcdda98b39735be99'
PROTOCOL = dict(
    id='observed_pcm_continuous_20260914_v1', classification='OBSERVED_DEVELOPMENT_ONLY',
    scientific_inference_ready=False, qualification_eligible=False,
    question='Do the four retained GH starts reach the same positive-curvature continuous stationary point, and how much do recalibration and integration together change their scores?',
    starts=[31, 61, 121, 181], start_selection='Previous start-SD 0.5 endpoints, all four orders; no best-start selection.',
    data='Unchanged previous 147 observed rows, weight sum 141.5, 24 persons/22 informative; fixed anchors, group mean and regression means unchanged.',
    integration=dict(bound=14., rel_tol=1e-11, abs_tol=1e-13),
    tight_integration=dict(bound=16., rel_tol=1e-12, abs_tol=1e-14),
    preliminary=dict(maxiter=250, gtol=1e-8, ftol=1e-15, maxls=50, log_sigma_relative_step=1e-6),
    score_polish=dict(method='trf', jac='3-point', max_nfev=100, gtol=1e-10, ftol=1e-12, xtol=1e-12),
    design='Two raw L-BFGS-B likelihood minimizations followed by two bounded least-squares solves of the raw continuous NLL gradient. The latter cost is half the squared score norm, not NLL. Assess actual NLL, score and information separately.',
    sigma_bounds=[.05, 10.], information_steps=[1e-4, 3e-5], gradient_steps=[1e-4, 1e-5],
    r_selection='Original GH181 point and final continuous point obtained from GH181, chosen before fitting.',
    r_integration=[dict(bound=12., rel_tol=1e-9, abs_tol=1e-11), dict(bound=16., rel_tol=1e-11, abs_tol=1e-13)],
    limits=dict(gradient=1e-8, standardized_score=1e-7, newton_correction=1e-7,
                cross_start_coordinates=1e-7, restart_coordinates=1e-7, raw_worsening=1e-9,
                information_condition=1e6, information_symmetry=1e-6, gradient_fd=1e-6,
                integration_refinement=1e-8, numeric_relative_mass_error=1e-7, tail_relative_mass=1e-10,
                r_nll=1e-8, r_eap=1e-9, r_sd=1e-9, r_log_sigma_score=1e-8,
                gh_structural=1e-4, gh_log_sigma=1e-4, gh_eap=1e-4, gh_sd=1e-4,
                gh_continuous_nll_loss_per_informative_person=1e-6, replay=1e-9),
    failure_policy='Freeze sources/settings before evaluation; retain all preliminary and score-solver failures and all returned/last-evaluation points; no retries with changed tolerances.',
    scope='One observed development design; a checked local minimum is not proof of a global optimum, recovery, SE/CI/coverage, or new native TAM/ConQuest agreement. No public estimator change.',
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    out = parser.parse_args().output_dir.resolve()
    load = lambda path: json.loads(path.read_text())
    assert sha(PREVIOUS/'summary.json') == PREVIOUS_HASH
    old, old_spec = load(PREVIOUS/'summary.json'), load(PREVIOUS/'protocol.json')
    for path, digest in old['artifact_sha256'].items():
        assert sha(PREVIOUS/path) == digest, path
    for path, digest in old_spec['source_sha256'].items():
        assert sha(ROOT/path) == digest, path
    r_script = Path(__file__).with_suffix('.R')
    sources = [Path(__file__).resolve(), r_script, ROOT/'validation/mml_observed_pcm_refit.py',
               ROOT/'validation/mml_pcm_continuous_refit.py', ROOT/'mfrm_app/mml_stationarity.py']
    hashes = dict(old_spec['source_sha256'], **{str(p.relative_to(ROOT)):sha(p) for p in sources})
    out.mkdir(parents=True, exist_ok=False)
    dump(out/'protocol.json', dict(PROTOCOL, source_sha256=hashes, previous_summary_sha256=PREVIOUS_HASH,
         runtime=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__)))
    data = load(PREVIOUS/'input.json')
    idx = data['indices']
    for key in ('person', 'score_k', 'weight'):
        idx[key] = np.array(idx[key])
    idx['facets'] = {key:np.array(value) for key, value in idx['facets'].items()}
    reference = ObservedPCMReference(SimpleNamespace(idx=idx, config={
        'n_person':data['n_person'], 'population_model':{'X':np.array(data['x'])[:, None]}}))
    assert np.array_equal(reference.design, data['affine_design'])
    assert np.array_equal(reference.offset, data['affine_offset'])
    dump(out/'input.json', dict(data, r_integration=PROTOCOL['r_integration'],
                              classification=PROTOCOL['classification'], scientific_inference_ready=False))
    evaluate = lambda point: reference.continuous(point, **PROTOCOL['integration'])
    tight = lambda point: reference.continuous(point, **PROTOCOL['tight_integration'])
    def joint(point):
        v = evaluate(point)
        return v['nll'], v['gradient']
    raw = lambda parameters, sigma: evaluate(np.r_[parameters, np.log(sigma)])['nll']
    def structural(parameters, sigma):
        v = evaluate(np.r_[parameters, np.log(sigma)])
        return v['nll'], v['gradient'][:-1]
    lower, upper = np.r_[np.full(9, -np.inf), np.log(.05)], np.r_[np.full(9, np.inf), np.log(10.)]
    limits = PROTOCOL['limits']
    results = {}
    for q in PROTOCOL['starts']:
        label = f'from_q{q}'
        old_record = old['results'][f'q{q}_sd0.5']
        start = np.array(old_record['coordinates'])
        initial = evaluate(start)
        initial_replay = discrepancies(initial, old_record['continuous'][-1])
        point = start.copy()
        preliminary = []
        for _ in range(2):
            run = polish_joint_free_sd(point[:-1], np.exp(point[-1]), raw, structural,
                sigma_bounds=tuple(PROTOCOL['sigma_bounds']), joint_value_gradient=joint,
                options=JointPolishOptions(**PROTOCOL['preliminary']))
            preliminary.append(run.to_dict())
            point = np.array(run.joint_coordinates)
        dump(out/(label+'_preliminary.json'), preliminary)
        score_runs = []
        for _ in range(2):
            first = point.copy()
            first_value = evaluate(first)
            last = None
            evaluations = 0
            def score(p):
                nonlocal last, evaluations
                v = evaluate(p)
                last = dict(coordinates=p.copy(), nll=v['nll'], gradient=v['gradient'].copy())
                evaluations += 1
                return v['gradient']
            fitted = least_squares(score, first, bounds=(lower, upper), **PROTOCOL['score_polish'])
            point = fitted.x.copy()
            value = evaluate(point)
            score_runs.append(dict(initial=first, initial_nll=first_value['nll'], returned=point,
                final_nll=value['nll'], final_gradient=value['gradient'], last_residual_evaluation=last,
                success=bool(fitted.success), status=int(fitted.status), message=str(fitted.message),
                optimizer_cost=float(fitted.cost), cost_definition='0.5 * sum(raw NLL gradient squared)',
                optimality=float(fitted.optimality), nfev=int(fitted.nfev), njev=int(fitted.njev),
                actual_residual_evaluations=evaluations))
        dump(out/(label+'_score_runs.json'), score_runs)
        value, tight_value = evaluate(point), tight(point)
        information = [information_diagnostics(joint, point, relative_step=h).to_dict()
                       for h in PROTOCOL['information_steps']]
        audits = [audit_joint_gradient(lambda p:evaluate(p)['nll'], joint, point, relative_step=h).to_dict()
                  for h in PROTOCOL['gradient_steps']]
        refinement = discrepancies(value, tight_value)
        moved = discrepancies(tight_value, initial, keys=('eap', 'sd'))
        total = discrepancies(tight_value, old_record['app'], keys=('eap', 'sd'))
        structural_difference = float(np.max(abs(point[:-1]-start[:-1])))
        log_sigma_difference = abs(float(point[-1]-start[-1]))
        nll_gain = float(initial['nll']-tight_value['nll'])
        replay = []
        for run in preliminary:
            for name, coords, nll in [('initial', run['initial_joint_coordinates'], run['initial_objective']),
                                      ('returned', run['joint_coordinates'], run['final_objective']),
                                      ('last_query', run['last_query_coordinates'], run['last_query_raw_objective'])]:
                replay.append(dict(stage='raw', point=name, error=abs(evaluate(np.array(coords))['nll']-nll)))
        for run in score_runs:
            for name, coords, nll in [('initial', run['initial'], run['initial_nll']),
                                      ('returned', run['returned'], run['final_nll']),
                                      ('last_residual_evaluation', run['last_residual_evaluation']['coordinates'], run['last_residual_evaluation']['nll'])]:
                replay.append(dict(stage='score', point=name, error=abs(evaluate(np.array(coords))['nll']-nll)))
        checks = dict(score_solver_terminated=all(r['success'] for r in score_runs),
            finite_stationary=float(np.max(abs(tight_value['gradient']))) < limits['gradient'],
            sigma_interior=.05 < tight_value['sigma'] < 10.,
            local_information=all(v['positive_definite'] and v['condition_number'] < limits['information_condition']
                and v['relative_symmetry_residual'] < limits['information_symmetry']
                and v['standardized_score_supnorm'] < limits['standardized_score']
                and v['newton_correction_supnorm'] < limits['newton_correction'] for v in information),
            gradient_audit=all(v['maximum_absolute_difference'] < limits['gradient_fd'] for v in audits),
            reference_refinement=max(refinement.values()) < limits['integration_refinement']
                and tight_value['numeric_relative_mass_error_sum'] < limits['numeric_relative_mass_error']
                and tight_value['tail_relative_mass_bound_sum'] < limits['tail_relative_mass'],
            raw_nonworsening=nll_gain >= -limits['raw_worsening']
                and all(r['final_nll']-r['initial_nll'] <= limits['raw_worsening'] for r in score_runs),
            restart=float(np.max(abs(score_runs[-1]['returned']-score_runs[-1]['initial']))) < limits['restart_coordinates'],
            replay=max(v['error'] for v in replay) < limits['replay'] and max(initial_replay.values()) < limits['replay'])
        results[label] = dict(initial_coordinates=start, coordinates=point, initial=initial, final=value, tight=tight_value,
            information=information, gradient_audits=audits, refinement=refinement, replay=replay,
            initial_replay=initial_replay, checks=checks, local_reference_pass=all(checks.values()),
            gh_comparison=dict(structural=structural_difference, log_sigma=log_sigma_difference,
                continuous_nll_improvement=nll_gain, continuous_score_change=moved, total_score_change=total,
                pass_engineering_targets=structural_difference < limits['gh_structural']
                    and log_sigma_difference < limits['gh_log_sigma']
                    and total['eap'] < limits['gh_eap'] and total['sd'] < limits['gh_sd']
                    and nll_gain/22 < limits['gh_continuous_nll_loss_per_informative_person']))
        dump(out/(label+'.json'), results[label])
        print(label, 'sigma', tight_value['sigma'], 'gradient', np.max(abs(tight_value['gradient'])),
              'local_reference', all(checks.values()), 'GH differences', results[label]['gh_comparison'], flush=True)
    r_cases = dict(gh181=dict(coordinates=results['from_q181']['initial_coordinates']),
                   continuous=dict(coordinates=results['from_q181']['coordinates']))
    dump(out/'r_cases.json', r_cases)
    with (out/'r.log').open('x') as log:
        status = subprocess.run(['Rscript', str(r_script), str(out/'input.json'), str(out/'r_cases.json'), str(out)],
                                stdout=log, stderr=subprocess.STDOUT).returncode
    r_comparisons = {}
    if status == 0:
        for name, case in r_cases.items():
            values = load(out/(name+'_r.json'))['continuous']
            python = tight(np.array(case['coordinates']))
            difference = discrepancies(python, values[-1], keys=('nll', 'eap', 'sd'))
            difference['log_sigma_score'] = abs(python['gradient'][-1]-values[-1]['log_sigma_nll_score'])
            r_comparisons[name] = dict(difference=difference,
                refinement=discrepancies(*values, keys=('nll', 'eap', 'sd')),
                passed=all(v < limits['r_'+k] for k,v in difference.items())
                    and max(discrepancies(*values, keys=('nll', 'eap', 'sd')).values()) < limits['integration_refinement'])
    cross_start = float(np.max(np.ptp([v['coordinates'] for v in results.values()], axis=0)))
    checks = dict(local_references=all(r['local_reference_pass'] for r in results.values()),
        independent_r=status == 0 and len(r_comparisons)==2 and all(r['passed'] for r in r_comparisons.values()),
        cross_start=cross_start < limits['cross_start_coordinates'],
        sources_unchanged=all(sha(ROOT/p)==h for p,h in hashes.items()))
    dump(out/'summary.json', dict(classification=PROTOCOL['classification'], scientific_inference_ready=False,
        qualification_eligible=False, checks=checks, reference_checks_pass=all(checks.values()),
        results=results, cross_start_coordinates=cross_start, r_exit_code=status, independent_r=r_comparisons,
        protocol_sha256=sha(out/'protocol.json'),
        artifact_sha256={p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file()}))
    print('Continuous reference checks:', checks, 'cross_start', cross_start, flush=True)
    return 0 if all(checks.values()) else 1


if __name__ == '__main__':
    raise SystemExit(main())
