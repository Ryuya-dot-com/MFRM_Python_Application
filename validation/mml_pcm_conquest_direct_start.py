#!/usr/bin/env python3
"""Prepare/run ConQuest GH from retained direct finite-GH solutions."""
import argparse
import json
from pathlib import Path
import shutil

import numpy as np
from scipy.special import roots_hermitenorm

from mml_pcm_conquest_history_audit import finite_gh
from mml_pcm_conquest_mc import SOURCE, SOURCE_HASH, rows
from mml_pcm_conquest_variance_replay import variance_updates
from mml_pcm_continuous_refit import PCMIntegral, PREVIOUS, SUMMARY_HASH, dump
from run_pcm_conquest_check import ROOT, launch, sha

ORDERS=(31,61,121,181)


def native_coordinates(par):
    par=np.asarray(par)
    assert par.shape==(24,) and np.all(np.isfinite(par))
    return np.r_[par[3:8],par[:3],par[8:23]]


def check_outputs(directory):
    expected=SOURCE/'cq_gh31'
    assert (directory/'amatrix.csv').read_bytes()==(expected/'amatrix.csv').read_bytes()
    pars=rows(directory/'parameters.csv')
    assert [r['Label'] for r in pars]==[r['Label'] for r in rows(expected/'parameters.csv')]
    assert [int(r['P']) for r in pars]==list(range(1,24))
    assert np.all(np.isfinite([float(r['Estimate']) for r in pars]))
    covariance=rows(directory/'covariance.csv'); beta=rows(directory/'reg_coefficients.csv')
    assert len(covariance)==len(beta)==1 and float(covariance[0]['Covariance'])>0
    assert float(beta[0]['Estimate'])==0
    assert [r['PID'] for r in rows(directory/'cases.csv')]==[f'P{i:02d}' for i in range(1,25)]
    assert rows(directory/'history.csv') and (directory/'review.txt').stat().st_size>0


def prepare(out):
    assert sha(SOURCE/'review_summary.json')==SOURCE_HASH
    assert sha(PREVIOUS/'summary.json')==SUMMARY_HASH
    for base in (SOURCE,PREVIOUS):
        record=json.loads((base/('review_summary.json' if base==SOURCE else 'summary.json')).read_text())
        for name,h in record['artifact_sha256'].items():
            assert sha(base/name)==h,name
    original=json.loads((PREVIOUS/'input.json').read_text())
    fits=json.loads((PREVIOUS/'python.json').read_text())
    problem=PCMIntegral(original['data'])
    out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(SOURCE/'wide.csv',out/'wide.csv')
    starts={};plan=[]
    for q in ORDERS:
        par=np.array(fits[f'q{q}_from_q31']['coordinates'])
        native=native_coordinates(par)
        assert np.array_equal(np.r_[native[5:8],native[:5],native[8:]],par[:-1])
        var=float(np.exp(2*par[-1]))
        (out/f'q{q}_init_parameters.txt').write_text(''.join(f'{i} {v:.17g}\n' for i,v in enumerate(native,1)))
        (out/f'q{q}_init_covariance.txt').write_text(f'1 1 {var:.17g}\n')
        z,w=roots_hermitenorm(q); w/=np.sqrt(2*np.pi)
        initial=finite_gh(problem,par,z,w)
        independent=json.loads((PREVIOUS/f'q{q}_from_q31_r.json').read_text())['finite'][str(q)]
        differences={k:float(np.max(abs(np.asarray(initial[k])-independent[k]))) for k in ('nll','eap','sd')}
        assert max(differences.values())<1e-8 and np.max(abs(initial['gradient']))<1e-4
        starts[str(q)]=dict(coordinates=par,variance=var,finite=initial,independent_r_differences=differences,
                           predicted_variance=variance_updates(initial['eap'],initial['sd']))
        source=(SOURCE/f'cq_gh{q}'/'model.cqc').read_text()
        for iterations in (1,2000):
            name=f'q{q}_iterations{iterations}'
            command=source.replace('title Observed PCM integration check;', 'title Direct finite-GH start check;')
            command=command.replace('p_nodes=2000, exit_on_error=yes;',
                'p_nodes=2000, seed=2, keeplastests=no, progress=no, exit_on_error=yes;')
            model='model criterion + rater + criterion*step;\n'
            assert command.count(model)==1 and 'iterations=2000;' in command
            command=command.replace(model,model+f'import init_parameters << ../q{q}_init_parameters.txt;\n'
                +f'import init_covariance << ../q{q}_init_covariance.txt;\n')
            command=command.replace('iterations=2000;',f'iterations={iterations};')
            (out/name).mkdir(); (out/name/'model.cqc').write_text(command)
            plan.append(dict(id=name,q=q,iterations=iterations))
    helpers=('mml_pcm_conquest_history_audit.py','mml_pcm_conquest_variance_replay.py',
             'mml_pcm_conquest_mc.py','mml_pcm_continuous_refit.py','run_pcm_conquest_check.py')
    paths=[Path(__file__).resolve(),*[ROOT/'validation'/p for p in helpers]]
    dump(out/'input.json',dict(classification='OBSERVED_DEVELOPMENT_ONLY',scientific_inference_ready=False,
        qualification_eligible=False,plan=plan,starts=starts,source_summary_sha256=SOURCE_HASH,
        previous_summary_sha256=SUMMARY_HASH,source_sha256={str(p.relative_to(ROOT)):sha(p) for p in paths},
        binary_sha256=sha(Path('/Applications/ConQuest/ConQuest')),integration=original['integration'][1],
        question='Does a native update leave an independently stationary finite-GH start, and where does continued iteration go?',
        controls=['CASES mean zero, free variance, no explicit mean anchor, identical design and data',
            'Full-precision imported structural coordinates and variance; one step and at most 2000 steps per order',
            'Keep native default best checkpoint and preserve complete history, including the last row'],
        metrics=['First reported NLL versus imported-start NLL',
            'First-step variance versus both prespecified posterior-moment formulas',
            'First/selected/last coordinate movements, finite NLL and gradient, continuous integral',
            'Iteration cap and native stopping reason retained separately from numerical interpretation'],
        input_sha256={str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file()}))


def main():
    assert np.array_equal(native_coordinates(np.arange(24.)),np.r_[np.arange(3.,8),np.arange(3.),np.arange(8.,23)])
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory',type=Path,required=True)
    parser.add_argument('--prepare-only',action='store_true')
    args=parser.parse_args();out=args.directory.resolve()
    if args.prepare_only:
        prepare(out); print('Prepared',out); return 0
    spec=json.loads((out/'input.json').read_text())
    for key,base in (('source_sha256',ROOT),('input_sha256',out)):
        for name,h in spec[key].items():
            assert sha(base/name)==h,name
    assert sha(Path('/Applications/ConQuest/ConQuest'))==spec['binary_sha256']
    (out/'sentinel').mkdir(exist_ok=False)
    sentinel=launch(out/'sentinel','quit;\n',60)
    if sentinel['exit_code']!=0 or not sentinel['end_of_program']:
        print('Startup failed; see',out/'sentinel/console.log');return 1
    results={}
    for condition in spec['plan']:
        name=condition['id'];directory=out/name
        result=launch(directory,(directory/'model.cqc').read_text(),600)
        result['output_checks_pass']=False
        if result['exit_code']==0 and result['end_of_program']:
            try:
                check_outputs(directory); result['output_checks_pass']=True
            except (OSError,KeyError,ValueError,AssertionError) as error:
                result['review_error']=repr(error)
        results[name]=result
        print(name,'Output checks:',result['output_checks_pass'],flush=True)
        if not result['output_checks_pass']: break
    dump(out/'execution.json',dict(scientific_inference_ready=False,qualification_eligible=False,runs=results,
        remaining_unattempted=[c['id'] for c in spec['plan'] if c['id'] not in results],
        artifact_sha256={str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file()}))
    complete=len(results)==8 and all(r['output_checks_pass'] for r in results.values())
    print('Execution complete:',complete,'Saved:',out,'Numerical review pending.')
    return 0 if complete else 1


if __name__=='__main__':
    raise SystemExit(main())
