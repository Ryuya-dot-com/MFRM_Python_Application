#!/usr/bin/env python3
"""Frozen range/density diagnostics and review of matched-start TAM runs."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import subprocess

import numpy as np

from mml_pcm_independent_conditions import PCMIntegral, fit_pair, difference
from mml_pcm_continuous_refit import ROOT, dump, sha, fd_check
from mml_pcm_external_review import fixed_grid

PROTOCOL = ROOT / 'validation/mml_pcm_grid_protocol.json'
PROTOCOL_HASH = 'd0bc84affcb7ff5eb081935976b9e82436d70d788ebea5bad9ed4d57644cd610'


def cells(spec):
    return [dict(id=f'b{b}_q{round(2*b/h)+1}',bound=b,spacing=h,q=round(2*b/h)+1)
            for b in spec['ranges'] for h in spec['spacings']]


def nodes(config):
    return np.linspace(-config['bound'],config['bound'],config['q'])


def precision(error, spec):
    t = spec['engineering_precision_targets']
    return bool(error['nll']/80 < t['nll_per_person'] and error['eap'] < t['max_eap_logit']
                and error['sd'] < t['max_posterior_sd_logit'])


def verify(spec):
    for p,h in spec['source_sha256'].items():
        assert sha(ROOT/p) == h,p
    for p,h in spec['input_sha256'].items():
        assert sha(Path(spec['directory'])/p) == h,p


def prepare(out):
    assert sha(PROTOCOL) == PROTOCOL_HASH
    spec = json.loads(PROTOCOL.read_text())
    previous = ROOT/spec['source_directory']
    assert sha(previous/'summary.json') == spec['source_summary_sha256']
    old = json.loads((previous/'summary.json').read_text())
    for p,h in old['artifact_sha256'].items():
        assert sha(previous/p) == h,p
    out.mkdir(parents=True,exist_ok=False)
    paths = [Path(__file__).resolve(),Path(__file__).with_suffix('.R'),PROTOCOL,
             ROOT/'validation/mml_pcm_independent_conditions.py',
             ROOT/'validation/mml_pcm_continuous_refit.py',ROOT/'validation/mml_pcm_external_review.py']
    hashes = {str(p.relative_to(ROOT)):sha(p) for p in paths}
    for name in spec['datasets']:
        folder = out/name; folder.mkdir()
        original = json.loads((previous/name/'input.json').read_text())
        dump(folder/'input.json',dict(classification='OBSERVED_DEVELOPMENT_ONLY',
            scientific_inference_ready=False,qualification_eligible=False,data=original['data'],
            integration=original['integration'],reference=old['results'][name]['reference'],
            source_data_sha256=sha(previous/name/'input.json')))
    dump(out/'input.json',dict(spec,directory=str(out),source_sha256=hashes,
        input_sha256={str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file()}))
    print('Prepared range/density and native TAM inputs:',out,flush=True)


def compute(out):
    spec = json.loads((out/'input.json').read_text()); verify(spec)
    for name in spec['datasets']:
        folder = out/name
        original = json.loads((folder/'input.json').read_text())
        problem = PCMIntegral(original['data'])
        reference = original['reference']; par = np.array(reference['coordinates'])
        ref_value = problem.evaluate(par,**original['integration'][1])
        assert max(difference(ref_value,reference['tight']).values()) < spec['implementation_limits']['same_point']
        results = {}
        fd = fd_check(lambda p:fixed_grid(problem,p,np.linspace(-12,12,241)),par)
        assert max(v['max_difference'] for v in fd) < spec['implementation_limits']['gradient_fd']
        for config in cells(spec):
            evaluate = lambda p:fixed_grid(problem,p,nodes(config))
            at_reference = evaluate(par)
            pair = fit_pair(evaluate,dict(zero=np.zeros(24),reference=par),spec)
            continuous = problem.evaluate(pair['coordinates'],**original['integration'][1])
            error = difference(pair['value'],continuous)
            result = dict(config=config,at_reference=at_reference,
                reference_integration_error=difference(at_reference,ref_value),
                reference_precision_target_met=precision(difference(at_reference,ref_value),spec),
                fit=pair,continuous=continuous,integration_error=error,
                precision_target_met=precision(error,spec),
                coordinate_distance_to_reference=float(np.max(abs(pair['coordinates']-par))),
                continuous_nll_excess_to_reference=continuous['nll']-ref_value['nll'])
            results[config['id']] = result
            dump(folder/(config['id']+'.json'),result)
            print(name,config['id'],'SD',np.exp(pair['coordinates'][-1]),'error',error,
                  'target',result['precision_target_met'],flush=True)
        dump(folder/'python.json',dict(results=results,gradient_fd=fd))
    verify(spec)


def native_cases(folder, problem, original, spec):
    records, cases = {}, {}
    for bound in (12,20):
        config = dict(id=f'b{bound}_q{bound*20+1}',bound=bound,spacing=.1,q=bound*20+1)
        native = json.loads((folder/f'tam_b{bound}.json').read_text())
        assert native['beta'] == [[0]] and native['initial']['beta'] == [[0]]
        assert [p['pid'] for p in native['person']] == [f'P{i:02d}' for i in range(1,81)]
        assert [p['score'] for p in native['person']] == problem.total.tolist()
        assert all(i['N'] == 80 for i in native['item'])
        byname = {v['item']:i for i,v in enumerate(native['item'])}
        order = [byname[f'C{c}-raterR{r}'] for r in range(1,5) for c in range(1,6)]
        A = np.array(native['A'])[order]
        assert np.array_equal(np.array(native['B'])[order,:,0],np.tile(np.arange(5),(20,1)))
        transform,_,rank,_ = np.linalg.lstsq(problem.design.reshape(100,23),A.reshape(100,23),rcond=None)
        design_error = float(np.max(abs(problem.design@transform-A)))
        assert rank == 23 and design_error < spec['implementation_limits']['design']
        def point(xsi, variance):
            return np.r_[transform@np.array(xsi),.5*np.log(float(np.asarray(variance).ravel()[0]))]
        returned = point(native['returned_xsi'],native['variance'])
        last = native['after'][-1]
        initial = point(native['initial']['xsi'],native['initial']['variance'])
        assert np.max(abs(initial)) == 0
        endpoints = {}
        for role,p in dict(returned=returned,last=point(last['xsi'],last['variance'])).items():
            finite = fixed_grid(problem,p,nodes(config))
            continuous = problem.evaluate(p,**original['integration'][1])
            error = difference(finite,continuous)
            endpoints[role] = dict(coordinates=p,finite=finite,continuous=continuous,
                integration_error=error,precision_target_met=precision(error,spec),
                coordinate_distance_to_reference=float(np.max(abs(p-original['reference']['coordinates']))))
            cases[f'tam_b{bound}_{role}'] = dict(coordinates=p,grids=[config],continuous=True)
        reported = np.array([r['deviance']/2 for r in native['after']])
        history_points = [point(r['xsi'],r['variance']) for r in native['after']]
        matching = [i+1 for i,p in enumerate(history_points) if np.max(abs(p-returned)) < 1e-12]
        old_nll = [fixed_grid(problem,point(r['xsi'],r['variance']),nodes(config))['nll'] for r in native['before']]
        checks = dict(variance_ever_small=any(r['variance_change'] < spec['tam']['conv'] for r in native['after']),
            beta_ever_small=any(r['beta_change'] < spec['tam']['conv'] for r in native['after']),
            terminal_xsi_small=last['xsi_change'] <= spec['tam']['conv'],
            terminal_deviance_small=last['deviance_change'] <= spec['tam']['convD'])
        records[str(bound)] = dict(config=config,initial=initial,endpoints=endpoints,
            native_iterations=native['iter'],reached_iteration_cap=native['reached_iteration_cap'],
            stopping_checks=checks,returned_matches_iterations=matching,
            design_error=design_error,reported_nll=native['deviance']/2,
            reported_nll_difference_to_returned=abs(native['deviance']/2-endpoints['returned']['finite']['nll']),
            history_nll_difference_to_before=float(np.max(abs(reported-old_nll))),
            stored_score_differences=difference(endpoints['returned']['finite'],
                dict(nll=native['deviance']/2,eap=[p['EAP'] for p in native['person']],sd=[p['SD.EAP'] for p in native['person']])),
            returned_axsi_difference=float(np.max(abs(problem.design@returned[:-1]-np.array(native['AXsi'])[order]))),
            warnings=native['warnings'])
    return records,cases


def review(out):
    spec = json.loads((out/'input.json').read_text()); verify(spec)
    results = {}
    for name in spec['datasets']:
        folder = out/name
        original = json.loads((folder/'input.json').read_text())
        problem = PCMIntegral(original['data'])
        python = json.loads((folder/'python.json').read_text())
        native,extra = native_cases(folder,problem,original,spec)
        cases = dict(reference=dict(coordinates=original['reference']['coordinates'],grids=cells(spec),continuous=True),
                     initial=dict(coordinates=np.zeros(24),grids=[g for g in cells(spec) if g['bound'] in (12,20) and g['spacing']==.1],continuous=False))
        for key,value in python['results'].items():
            cases[key] = dict(coordinates=value['fit']['coordinates'],grids=[value['config']],continuous=True)
        cases.update(extra)
        dump(folder/'r_cases.json',cases)
        results[name] = dict(python=python,native_tam=native)
    def run_r(name):
        folder = out/name
        with (folder/'r_review.log').open('x') as log:
            result = subprocess.run(['Rscript',str(Path(__file__).with_suffix('.R')),
                str(folder/'input.json'),str(folder),'evaluate'],stdout=log,stderr=subprocess.STDOUT)
        assert result.returncode == 0,folder/'r_review.log'
        print(name,'independent R complete',flush=True)
    with ThreadPoolExecutor(max_workers=4) as executor:
        list(executor.map(run_r,spec['datasets']))
    for name,result in results.items():
        folder = out/name
        original = json.loads((folder/'input.json').read_text())
        problem = PCMIntegral(original['data'])
        cases = json.loads((folder/'r_cases.json').read_text())
        checks = {}
        for key,case in cases.items():
            independent = json.loads((folder/f'{key}_r.json').read_text())
            p = np.array(case['coordinates'])
            finite = {g['id']:difference(fixed_grid(problem,p,nodes(g)),independent['finite'][g['id']]) for g in case['grids']}
            continuous = difference(problem.evaluate(p,**original['integration'][1]),independent['continuous']) if case['continuous'] else None
            checks[key] = dict(finite=finite,continuous=continuous)
        arithmetic = all(d['nll'] < spec['implementation_limits']['same_point'] and
            max(d['eap'],d['sd']) < spec['implementation_limits']['independent_moments']
            for c in checks.values() for d in [*c['finite'].values(),*([c['continuous']] if c['continuous'] else [])])
        direct = result['python']['results']
        optimization = all(np.max(abs(np.array(s['value']['gradient']))) < spec['implementation_limits']['stationary_gradient'] and
            all(r['success'] for r in s['runs']) for f in direct.values() for s in f['fit']['starts'].values())
        comparisons = []
        for bound in spec['ranges']:
            a,b = [direct[g['id']] for g in cells(spec) if g['bound']==bound]
            comparisons.append(dict(changed='spacing',bound=bound,coordinate_difference=float(np.max(abs(
                np.array(a['fit']['coordinates'])-b['fit']['coordinates']))),
                same_point_difference=difference(a['at_reference'],b['at_reference']),
                refitted_score_difference=difference(a['fit']['value'],b['fit']['value'])))
        result.update(independent_r=checks,arithmetic_checks_pass=bool(arithmetic),
            direct_optimization_checks_pass=bool(optimization),spacing_comparisons=comparisons)
        dump(folder/'summary.json',result)
        print(name,'arithmetic',arithmetic,'direct optimization',optimization,flush=True)
    verify(spec)
    dump(out/'summary.json',dict(classification=spec['classification'],scientific_inference_ready=False,
        qualification_eligible=False,results=results,protocol_sha256=PROTOCOL_HASH,
        implementation_checks_pass=all(r['arithmetic_checks_pass'] and r['direct_optimization_checks_pass'] for r in results.values()),
        source_sha256=spec['source_sha256'],
        native_stopping_is_not_stationarity=True,
        artifact_sha256={str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file()}))
    print('Review saved:',out/'summary.json',flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['prepare','compute','review'])
    parser.add_argument('--directory',type=Path,required=True)
    args = parser.parse_args()
    dict(prepare=prepare,compute=compute,review=review)[args.action](args.directory.resolve())
