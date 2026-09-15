#!/usr/bin/env python3
"""Review eight user-run ConQuest grid fits against TAM and independent R."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import subprocess

import numpy as np

from mml_pcm_continuous_refit import ROOT, dump, sha, fd_check
from mml_pcm_independent_conditions import PCMIntegral, difference
from mml_pcm_conquest_mc import rows
from mml_pcm_conquest_history_audit import coordinates
from mml_pcm_external_review import fixed_grid
from mml_pcm_grid_check import nodes, precision


def normalized_grid(problem, par, config):
    theta = nodes(config)
    result = fixed_grid(problem, par, theta)
    sigma = np.exp(par[-1])
    prior = np.exp(-.5*(theta/sigma)**2)
    prior /= prior.sum()
    correction = len(problem.y)*float(np.sum(prior*((theta/sigma)**2-1)))
    gradient = result['gradient'].copy(); gradient[-1] += correction
    return dict(result,raw_nll=result['nll'],raw_gradient=result['gradient'],
        nll=result['nll']+len(problem.y)*np.log(result['raw_prior_mass']),gradient=gradient)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory',type=Path,required=True)
    parser.add_argument('--grid-directory',type=Path,required=True)
    args = parser.parse_args()
    source, grids = args.directory.resolve(),args.grid_directory.resolve()
    execution = json.loads((source/'execution.json').read_text())
    assert len(execution['runs']) == 8 and not execution['remaining_unattempted']
    assert all(r['output_checks_pass'] for r in execution['runs'].values())
    for p,h in execution['artifact_sha256'].items():
        assert sha(source/p) == h,p
    input_spec = json.loads((source/'input.json').read_text())
    spec = input_spec['protocol']
    out = source/'review'; out.mkdir(exist_ok=False)
    dependencies = [Path(__file__).resolve(),ROOT/'validation/mml_pcm_grid_check.R',
        ROOT/'validation/mml_pcm_independent_conditions.py',ROOT/'validation/mml_pcm_external_review.py',
        ROOT/'validation/mml_pcm_conquest_history_audit.py',ROOT/'validation/mml_pcm_grid_check.py']
    source_hashes = {str(p.relative_to(ROOT)):sha(p) for p in dependencies}
    results, inputs = {}, {}
    for dataset in spec['datasets']:
        folder = out/dataset; folder.mkdir()
        original = json.loads((grids/dataset/'input.json').read_text())
        tam_cases_path = grids/dataset/'r_cases.json'
        tam_cases = json.loads(tam_cases_path.read_text())
        inputs[str(tam_cases_path)] = sha(tam_cases_path)
        problem = PCMIntegral(original['data'])
        cases = {'initial':dict(coordinates=np.zeros(24),grids=[],continuous=False)}
        records = {}
        for bound in (12,20):
            name = f'{dataset}_b{bound}'; directory = source/name
            config = dict(id=f'b{bound}_q{bound*20+1}',bound=bound,spacing=.1,q=bound*20+1)
            cases['initial']['grids'].append(config)
            history = rows(directory/'history.csv')
            assert [int(r['Iteration']) for r in history] == list(range(1,len(history)+1))
            last = coordinates(history[-1])
            xsi = np.array([float(r['Estimate']) for r in rows(directory/'parameters.csv')])
            variance = float(rows(directory/'covariance.csv')[0]['Covariance'])
            assert float(rows(directory/'reg_coefficients.csv')[0]['Estimate']) == 0
            returned = np.r_[xsi[5:8],xsi[:5],xsi[8:],.5*np.log(variance)]
            tam = np.array(tam_cases[f'tam_b{bound}_returned']['coordinates'])
            endpoints = {}
            for role,p in dict(returned=returned,last=last).items():
                finite = normalized_grid(problem,p,config)
                continuous = problem.evaluate(p,**original['integration'][1])
                error = difference(dict(finite,nll=finite['raw_nll']),continuous)
                endpoints[role] = dict(coordinates=p,finite=finite,continuous=continuous,
                    raw_grid_integration_error=error,normalized_grid_integration_error=difference(finite,continuous),
                    raw_grid_precision_target_met=precision(error,spec),
                    coordinate_difference_to_tam=float(np.max(abs(p-tam))))
                cases[f'b{bound}_{role}'] = dict(coordinates=p,grids=[config],continuous=True)
            tam_finite = normalized_grid(problem,tam,config)
            initial = normalized_grid(problem,np.zeros(24),config)
            native_scores = rows(directory/'cases.csv')
            assert [r['PID'] for r in native_scores] == [f'P{i:02d}' for i in range(1,81)]
            report = (directory/'review.txt').read_text()
            termination = [line.strip() for line in report.splitlines() if 'Iterations terminated' in line]
            assert len(termination) == 1
            records[str(bound)] = dict(config=config,endpoints=endpoints,native_iterations=len(history),
                termination=termination[0],reported_nll=float(history[-1]['LogLikelihood'])/2,
                reported_nll_difference_to_normalized=abs(float(history[-1]['LogLikelihood'])/2-endpoints['returned']['finite']['nll']),
                reported_nll_difference_to_raw=abs(float(history[-1]['LogLikelihood'])/2-endpoints['returned']['finite']['raw_nll']),
                first_reported_nll_difference_to_zero_start=abs(float(history[0]['LogLikelihood'])/2-initial['nll']),
                returned_last_coordinate_difference=float(np.max(abs(returned-last))),
                deterministic_score_differences_to_tam=difference(endpoints['returned']['finite'],tam_finite),
                native_mc_differences=dict(eap=float(np.max(abs(np.array([float(r['EAP_1']) for r in native_scores])-endpoints['returned']['finite']['eap']))),
                    sd=float(np.max(abs(np.array([float(r['PosteriorSD_1']) for r in native_scores])-endpoints['returned']['finite']['sd'])))))
            if dataset == 'wide_1' and bound == 12:
                records[str(bound)]['normalized_gradient_fd'] = fd_check(lambda p:normalized_grid(problem,p,config),returned)
                assert max(c['max_difference'] for c in records[str(bound)]['normalized_gradient_fd']) < spec['implementation_limits']['gradient_fd']
        dump(folder/'input.json',original)
        dump(folder/'r_cases.json',cases)
        results[dataset] = records
    dump(out/'input.json',dict(protocol=spec,execution_sha256=sha(source/'execution.json'),source_sha256=source_hashes,
        tam_case_sha256=inputs,scientific_inference_ready=False,qualification_eligible=False))
    def run_r(dataset):
        folder = out/dataset
        with (folder/'r_review.log').open('x') as log:
            r = subprocess.run(['Rscript',str(ROOT/'validation/mml_pcm_grid_check.R'),str(folder/'input.json'),str(folder),'evaluate'],stdout=log,stderr=subprocess.STDOUT)
        assert r.returncode == 0,folder/'r_review.log'
        print(dataset,'ConQuest coordinates independently checked',flush=True)
    with ThreadPoolExecutor(max_workers=4) as executor:
        list(executor.map(run_r,spec['datasets']))
    checks = {}
    for dataset,records in results.items():
        folder = out/dataset
        original = json.loads((folder/'input.json').read_text())
        problem = PCMIntegral(original['data'])
        cases = json.loads((folder/'r_cases.json').read_text())
        checks[dataset] = {}
        for key,case in cases.items():
            r = json.loads((folder/f'{key}_r.json').read_text())
            point = np.array(case['coordinates'])
            finite = {}
            for g in case['grids']:
                value = normalized_grid(problem,point,g)
                finite[g['id']] = dict(difference(dict(value,nll=value['raw_nll']),r['finite'][g['id']]),
                    normalized_nll=abs(value['nll']-r['finite'][g['id']]['normalized_prior_nll']))
            continuous = difference(problem.evaluate(point,**original['integration'][1]),r['continuous']) if case['continuous'] else None
            checks[dataset][key] = dict(finite=finite,continuous=continuous)
    passed = all(max(d.values()) < spec['implementation_limits']['same_point'] for v in checks.values()
        for c in v.values() for d in [*c['finite'].values(),*([c['continuous']] if c['continuous'] else [])])
    assert source_hashes == {p:sha(ROOT/p) for p in source_hashes}
    assert inputs == {p:sha(p) for p in inputs}
    dump(out/'summary.json',dict(scientific_inference_ready=False,qualification_eligible=False,results=results,
        independent_r=checks,independent_arithmetic_checks_pass=passed,source_sha256=source_hashes,
        execution_sha256=sha(source/'execution.json'),
        limitations=['Eight observed native fits, not equivalence or inference qualification',
            'Native parameter/history exports are rounded; common-objective gradients retain this uncertainty',
            'Native EAP/SD use 2000 MC draws and are not deterministic quadrature scores'],
        artifact_sha256={str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file()}))
    print('ConQuest independent arithmetic:',passed,'Saved:',out,flush=True)
    return 0 if passed else 1


if __name__ == '__main__':
    raise SystemExit(main())
