#!/usr/bin/env python3
"""Review checkpoint controls without equating a reported NLL with returned coordinates."""
import argparse
import json
from pathlib import Path
import subprocess

import numpy as np
from scipy.special import roots_hermitenorm

from mml_pcm_conquest_checkpoint_check import SOURCE as DIRECT, SOURCE_HASH as DIRECT_HASH
from mml_pcm_conquest_history_audit import coordinates, finite_gh
from mml_pcm_conquest_mc import SOURCE, SOURCE_HASH, rows
from mml_pcm_continuous_refit import PCMIntegral, PREVIOUS, dump
from run_pcm_conquest_check import ROOT, sha


def selection_summary(nll, selected):
    nll=np.asarray(nll,dtype=float);assert nll.ndim==1 and len(nll)>0 and np.all(np.isfinite(nll))
    assert 1<=selected<=len(nll)
    candidates=np.flatnonzero(np.diff(nll)<0)+1
    best=None if len(candidates)==0 else float(nll[candidates].min())
    return dict(global_minimum_iteration=int(np.argmin(nll))+1,global_minimum_nll=float(nll.min()),
        improving_iterations=(candidates+1).tolist(),minimum_among_improving=best,
        selected_iteration=selected,selected_reported_nll=float(nll[selected-1]),
        selected_is_visible_improving_iteration=bool(selected-1 in candidates),
        selected_matches_candidate_nll=None if best is None else bool(nll[selected-1]==best))


def main():
    c=selection_summary([1,3,2],3)
    assert c['global_minimum_iteration']==1 and c['minimum_among_improving']==2 and c['selected_matches_candidate_nll']
    assert selection_summary([1,2,3],3)['minimum_among_improving'] is None
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory',type=Path,required=True)
    out=parser.parse_args().directory.resolve()
    execution=json.loads((out/'summary.json').read_text());spec=json.loads((out/'input.json').read_text())
    assert len(execution['runs'])==4 and not execution['remaining_unattempted']
    assert all(r['output_checks_pass'] and r['history_matches_previous_prefix'] for r in execution['runs'].values())
    for record,field,base in ((execution,'artifact_sha256',out),(spec,'source_sha256',ROOT),(spec,'input_sha256',out)):
        for name,h in record[field].items():assert sha(base/name)==h,name
    assert sha(DIRECT/'review/summary.json')==DIRECT_HASH and sha(SOURCE/'review_summary.json')==SOURCE_HASH
    review=out/'review';review.mkdir(exist_ok=False)
    original=json.loads((PREVIOUS/'input.json').read_text());problem=PCMIntegral(original['data'])
    r_script=ROOT/'validation/mml_pcm_quadrature_diagnosis.R'
    paths=[Path(__file__).resolve(),r_script]
    hashes=spec['source_sha256']|{str(p.relative_to(ROOT)):sha(p) for p in paths}
    input_spec={k:original[k] for k in ('classification','scientific_inference_ready','qualification_eligible','data','integration')}
    input_spec.update(orders=[61],source_sha256=hashes,execution_sha256=sha(out/'summary.json'))
    dump(review/'evaluation_input.json',input_spec)
    z,w=roots_hermitenorm(61);w/=np.sqrt(2*np.pi)
    results={};cases={}
    for name,run in execution['runs'].items():
        history=rows(out/name/'history.csv');selected=run['selected_iteration'];p=coordinates(history[selected-1])
        cases[name]=dict(coordinates=p,finite=finite_gh(problem,p,z,w),continuous=problem.evaluate(p,**original['integration'][1]))
        results[name]=selection_summary([float(r['LogLikelihood'])/2 for r in history],selected)
        results[name].update(returned_coordinate_nll=cases[name]['finite']['nll'],
                             reported_minus_returned_nll=run['reported_nll']-cases[name]['finite']['nll'])
    legacy={}
    cold=json.loads((SOURCE/'review_summary.json').read_text())
    direct=json.loads((DIRECT/'review/summary.json').read_text())
    for q in (31,61,121,181):
        for name,path,selected in (
            (f'default_q{q}',SOURCE/f'cq_gh{q}',cold['fits'][f'cq_gh{q}']['selected_iteration']),
            (f'direct_q{q}',DIRECT/f'q{q}_iterations2000',direct['results'][str(q)]['selected_iteration'])):
            legacy[name]=selection_summary([float(r['LogLikelihood'])/2 for r in rows(path/'history.csv')],selected)
    dump(review/'python.json',cases)
    subprocess.run(['Rscript',str(r_script),str(review/'evaluation_input.json'),str(review/'python.json'),str(review)],check=True)
    differences={}
    for name,case in cases.items():
        r=json.loads((review/f'{name}_r.json').read_text())
        differences[name]={kind:{k:float(np.max(abs(np.asarray(case[kind][k])-target[k]))) for k in ('nll','eap','sd')}
            for kind,target in (('finite',r['finite']['61']),('continuous',r['continuous'][1]))}
    passed=all(d<1e-8 for r in differences.values() for v in r.values() for d in v.values())
    assert hashes=={p:sha(ROOT/p) for p in hashes}
    dump(review/'summary.json',dict(classification='OBSERVED_DEVELOPMENT_ONLY',scientific_inference_ready=False,
        qualification_eligible=False,independent_arithmetic_checks_pass=passed,results=results,legacy=legacy,
        independent_r_differences=differences,source_sha256=hashes,
        candidate_rule_is_native_contract=False,
        artifact_sha256={p.name:sha(p) for p in review.iterdir() if p.is_file()},
        limitations=['Candidate selection among improving rows is a post-hoc hypothesis, not a documented algorithm',
            'CSV ties can conceal strict improvements; compare values separately from exact selected indices',
            'keep-last controls which point is returned; it does not qualify stationarity or integration precision']))
    print('Independent arithmetic:',passed,'Saved',review)
    return 0 if passed else 1


if __name__=='__main__':raise SystemExit(main())
