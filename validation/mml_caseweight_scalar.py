#!/usr/bin/env python3
"""One free criterion separates weighted updates, displayed NLL and returned points."""
import argparse
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import traceback
from types import SimpleNamespace

import numpy as np
from scipy.optimize import brentq
from scipy.special import logsumexp

from mml_observed_pcm_refit import ObservedPCMReference
from mml_pcm_continuous_refit import ROOT, dump, sha
from mml_pcm_conquest_mc import rows
from mml_unit_weight_native import CELLS, native_parameters
from mml_unit_weight_conquest_mc import verify_prior
from mml_weight_semantics import WEIGHTS, verify_sources
from run_pcm_conquest_check import launch

PREVIOUS=ROOT/'validation/generated/mml_weight_semantics_20260914'
PREVIOUS_HASH='7dfacf2e89fbaaec8b2f33fcc07ee8797b8f74489062c7858f16d6a00e517eef'
VARIANTS=dict(one=dict(iterations=1,keep='yes'),last=dict(iterations=100,keep='yes'),best=dict(iterations=100,keep='no'))


def previous_check():
    assert sha(PREVIOUS/'review_summary.json')==PREVIOUS_HASH
    s=json.loads((PREVIOUS/'review_summary.json').read_text())
    for f,h in s['artifact_sha256'].items():assert sha(PREVIOUS/f)==h,f
    verify_sources(PREVIOUS,json.loads((PREVIOUS/'protocol.json').read_text()))
    review=json.loads((PREVIOUS/'review_protocol.json').read_text())
    assert sha(ROOT/'validation/mml_weight_semantics.py')==review['review_source_sha256']
    return verify_prior()[1]


def evaluator(data):
    ref=ObservedPCMReference(SimpleNamespace(idx=data['indices'],config=dict(n_person=32,population_model=dict(X=np.array(data['x'])[:,None]))))
    p=np.array(data['continuous_point']);theta=np.linspace(-12,12,241);mu=ref.x*p[0];sigma=np.exp(p[-1])
    base=(ref.offset+ref.design@p[:-1])[:,None,:]+(theta[None,:,None]-mu[ref.person,None,None])*ref.k
    direction=ref.design[:,:,3][:,None,:]
    lw=-.5*((theta[None,:]-mu[:,None])/sigma)**2-np.log(sigma)-.5*np.log(2*np.pi)+np.log(.1)
    prior_mass=logsumexp(lw,axis=1)
    def evaluate(delta,weights,posterior=None):
        logits=base+delta*direction;lp=logits-logsumexp(logits,axis=2,keepdims=True);prob=np.exp(lp)
        ed=np.sum(prob*direction,axis=2);vd=np.sum(prob*direction**2,axis=2)-ed**2
        ll=np.zeros((32,241));score=ll.copy();info=ll.copy()
        np.add.at(ll,ref.person,lp[np.arange(len(ref.y)),:,ref.y])
        np.add.at(score,ref.person,direction[np.arange(len(ref.y)),0,ref.y][:,None]-ed)
        np.add.at(info,ref.person,vd)
        mass=logsumexp(ll+lw,axis=1);post=np.exp(ll+lw-mass[:,None]) if posterior is None else posterior
        es=np.sum(post*score,axis=1)
        hessian=np.sum(post*info,axis=1)
        if posterior is None:hessian-=np.sum(post*(score-es[:,None])**2,axis=1)
        return dict(nll=float(weights@(-mass+prior_mass)),raw_nll=float(weights@(-mass)),
            unweighted_nll=float(np.sum(-mass+prior_mass)),gradient=float(-weights@es),hessian=float(weights@hessian),
            posterior=post)
    return evaluate


def compact(value):
    return {k:v for k,v in value.items() if k!='posterior'}


def derive(data):
    evaluate=evaluator(data);weights={k:32*w/sum(w) for k,w in WEIGHTS.items()}
    roots={k:brentq(lambda d:evaluate(d,w)['gradient'],-2,2,xtol=1e-13) for k,w in weights.items()}
    start=roots['unit'];post=evaluate(start,weights['unit'])['posterior'];result={}
    for key,w in weights.items():
        em=brentq(lambda d:evaluate(d,w,post)['gradient'],-2,2,xtol=1e-13)
        checks=[]
        for d in (start,em,roots[key]):
            v=evaluate(d,w)
            for h in (1e-4,3e-5):
                plus=evaluate(d+h,w);minus=evaluate(d-h,w)
                checks.append(dict(delta=d,step=h,gradient_error=abs((plus['nll']-minus['nll'])/(2*h)-v['gradient']),
                    hessian_error=abs((plus['gradient']-minus['gradient'])/(2*h)-v['hessian'])))
        assert max(max(c['gradient_error'],c['hessian_error']) for c in checks)<2e-6
        assert abs(evaluate(roots[key],w)['gradient'])<1e-9 and evaluate(roots[key],w)['hessian']>0
        result[key]=dict(root=roots[key],first_em=em,start=compact(evaluate(start,w)),root_value=compact(evaluate(roots[key],w)),derivative_checks=checks)
    return dict(start_delta=start,weights=weights,schemes=result)


def review(out,spec,refs,run_r):
    results={};independent={};controls={}
    for cell in CELLS:
        folder=out/cell;data=json.loads((folder/'input.json').read_text());evaluate=evaluator(data)
        p=np.array(data['continuous_point']);native0=native_parameters(p);start=refs[cell]['start_delta'];results[cell]={};controls[cell]={};points={}
        def point(tag,delta,w):points[tag]=dict(delta=delta,weights=w)
        for scheme,w in refs[cell]['weights'].items():
            w=np.asarray(w);cal=folder/scheme;baseline=PREVIOUS/cell/scheme/'cq';results[cell][scheme]={}
            exact=refs[cell]['schemes'][scheme]
            point(scheme+'_start',start,w);point(scheme+'_root',exact['root'],w);point(scheme+'_em',exact['first_em'],w)
            for variant in VARIANTS:
                cq=cal/variant;execution=json.loads((cq/'launch.json').read_text())
                assert execution['exit_code']==0 and execution['end_of_program'] and not execution['timed_out']
                for filename in ('reg_coefficients.csv','covariance.csv','amatrix.csv','scored.csv','converted.txt'):
                    assert sha(cq/filename)==sha(baseline/filename),cq/filename
                pars=rows(cq/'parameters.csv');oldpars=rows(baseline/'parameters.csv')
                assert len(pars)==8 and pars[1:]==oldpars[1:] and pars[0]['Label']==oldpars[0]['Label']
                returned=float(pars[0]['Estimate'])-native0[0];history=rows(cq/'history.csv');oldrow=rows(baseline/'history.csv')[0]
                assert [int(r['Iteration']) for r in history]==list(range(1,len(history)+1))
                report=(cq/'review.txt').read_text()
                assert int(re.search(r'Total number of estimated parameters:\s*(-?\d+)',report)[1])==1
                assert float(re.search(r'Weighted number of cases in MML/MCMC estimation:\s*([\d.]+)',report)[1])==sum(WEIGHTS[scheme])
                trajectory=[]
                for row in history:
                    assert all(row[k]==oldrow[k] for k in oldrow if k not in ('RowLabels','Iteration','LogLikelihood','xsi 1'))
                    delta=float(row['xsi 1'])-native0[0];v=evaluate(delta,w);unit=evaluate(delta,np.ones(32))
                    nll=float(row['LogLikelihood'])/2;error=abs(nll-v['unweighted_nll'])
                    allowance=spec['limits']['nll_display']+spec['limits']['coordinate_rounding']*abs(unit['gradient'])
                    trajectory.append(dict(iteration=int(row['Iteration']),delta=delta,reported_nll=nll,display_error=error,
                        display_matches_unweighted=bool(error<allowance),weighted_display_error=abs(nll-v['nll']),
                        value=compact(v),unweighted_gradient=unit['gradient'],display_rounding_allowance=allowance))
                v=evaluate(returned,w);selected=[r['iteration'] for r in trajectory if r['delta']==returned]
                # Default return selection can in principle choose the initialization.
                assert selected or abs(returned-start)<spec['limits']['coordinate_rounding']
                final=trajectory[-1];first=trajectory[0]
                results[cell][scheme][variant]=dict(returned_delta=returned,returned=compact(v),selected_iterations=selected,
                    iterations=len(history),termination=[l.strip() for l in report.splitlines() if 'Iterations terminated' in l],
                    reached_cap=len(history)>=VARIANTS[variant]['iterations'],trajectory=trajectory,
                    first_em_error=abs(first['delta']-exact['first_em']),root_error=abs(returned-exact['root']),
                    first_em_pass=abs(first['delta']-exact['first_em'])<spec['limits']['first_em'],
                    scalar_root_pass=abs(returned-exact['root'])<spec['limits']['root'],
                    returned_weighted_nll_change=v['nll']-exact['start']['nll'],
                    returned_unweighted_nll_change=v['unweighted_nll']-exact['start']['unweighted_nll'],
                    last_delta=final['delta'],last_weighted_gradient=final['value']['gradient'])
                point(scheme+'_'+variant+'_returned',returned,w)
            calls=results[cell][scheme]
            assert calls['one']['iterations']==1
            controls[cell][scheme]=dict(first_step_repeat=calls['one']['trajectory'][0]['delta']==calls['last']['trajectory'][0]['delta'],
                return_setting_preserves_history=calls['last']['trajectory']==calls['best']['trajectory'],
                keep_yes_returns_last=calls['last']['returned_delta']==calls['last']['last_delta'])
        for variant in VARIANTS:
            controls[cell]['scale_'+variant]=results[cell]['pattern'][variant]==results[cell]['triple'][variant]
        if run_r:
            dump(folder/'r_points.json',dict(points=points,weights=refs[cell]['weights']))
            with (folder/'r.log').open('x') as log:
                result=subprocess.run(['Rscript',str(Path(__file__).with_suffix('.R')),str(folder/'input.json'),str(folder/'r_points.json'),str(folder/'independent_R.json')],stdout=log,stderr=subprocess.STDOUT)
            assert result.returncode==0,folder/'r.log'
        r=json.loads((folder/'independent_R.json').read_text());independent[cell]={}
        for name,pt in points.items():
            expected=evaluate(pt['delta'],np.asarray(pt['weights']));actual=r['points'][name]
            independent[cell][name]={k:abs(actual[k]-expected[k]) for k in actual}
        root_errors={k:abs(r['roots'][k]['delta']-refs[cell]['schemes'][k]['root']) for k in WEIGHTS}
        independent[cell]['root_errors']=root_errors
        assert max(v for d in independent[cell].values() for v in d.values())<spec['limits']['independent_R']
    return dict(results=results,independent_R=independent,controls=controls)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output-dir',type=Path,required=True);parser.add_argument('--audit',action='store_true')
    args=parser.parse_args();out=args.output_dir.resolve();old_spec=previous_check()
    if args.audit:
        spec=json.loads((out/'protocol.json').read_text());saved=json.loads((out/'summary.json').read_text());refs=json.loads((out/'references.json').read_text())
        for f,h in saved['artifact_sha256'].items():assert sha(out/f)==h,f
        snapshots=json.loads((out/'source_snapshot/manifest.json').read_text()) if (out/'source_snapshot/manifest.json').exists() else {}
        for f,h in spec['source_sha256'].items():
            path=ROOT/f
            if sha(path)!=h and f in snapshots:path=out/snapshots[f]['path']
            assert sha(path)==h,f
        if (out/'audit_implementation.json').exists():
            assert sha(Path(__file__).resolve())==json.loads((out/'audit_implementation.json').read_text())['source_sha256']
        assert sha(out/'protocol.json')==saved['protocol_sha256']
        replay=review(out,spec,refs,False)
        assert json.loads(json.dumps(replay,default=lambda value:value.item()))=={k:saved[k] for k in replay}
        print('Saved audit passed:',len(saved['artifact_sha256']),'artifacts;',len(spec['source_sha256']),'sources; 36 native runs')
        return
    sources=dict(old_spec['source_sha256'])
    for path in (Path(__file__).resolve(),Path(__file__).with_suffix('.R').resolve(),ROOT/'validation/mml_weight_semantics.py',ROOT/'validation/mml_unit_weight_conquest_mc.py'):
        sources[str(path.relative_to(ROOT))]=sha(path)
    out.mkdir(parents=True,exist_ok=False)
    spec=dict(classification='OBSERVED_DEVELOPMENT_ONLY',scientific_inference_ready=False,qualification_eligible=False,
        question='With only native criterion C1 free, do case weights determine the first EM update and stationary point while reported deviance and best-estimate selection follow the unweighted objective?',
        cells=CELLS,case_weights=WEIGHTS,variants=VARIANTS,free_parameter='Native xsi 1; app criterion 1 (coordinate 3, zero-based); all seven other structural values, regression and variance anchored',
        start='Independently solve the unweighted finite-grid scalar score on delta [-2,2], then use the same start for every weight/return setting',
        grid='Physical theta [-12,12], spacing .1, 241 points; fixed regression means and variance, hence prior normalization has zero derivative in the free coordinate',
        controls='4 cells x 3 positive weight vectors x 3 iteration/return settings = 36 calls; preserve all history, returned/final roles, first-step and scale invariance; no replacement',
        reference='Analytic marginal score and Louis curvature; central differences at start/EM prediction/scalar root with steps 1e-4 and 3e-5; scalar EM maximization with old posterior. Independent R formula and uniroot at selected points.',
        limits=dict(nll_display=5.1e-7,coordinate_rounding=5.1e-7,independent_R=1e-8,first_em=2e-6,root=2e-6),
        source_sha256=sources,previous_summary_sha256=PREVIOUS_HASH,binary_sha256=old_spec['binary_sha256'],
        source='https://conquestmanual.acer.org/s4-00.html',
        scope='One free structural coordinate in four paired observed cells, finite range only; no general multi-parameter optimizer, continuous accuracy, response-weight equivalence, native score precision or SE/CI/coverage qualification')
    dump(out/'protocol.json',spec);refs={}
    for cell in CELLS:
        folder=out/cell;folder.mkdir();shutil.copyfile(PREVIOUS/cell/'input.json',folder/'input.json')
        data=json.loads((folder/'input.json').read_text());refs[cell]=derive(data)
        native=native_parameters(np.array(data['continuous_point']));native[0]+=refs[cell]['start_delta']
        for scheme in WEIGHTS:
            cal=folder/scheme;cal.mkdir();old=PREVIOUS/cell/scheme
            for filename in ('data.csv','anchor_beta.txt','anchor_covariance.txt'):shutil.copyfile(old/filename,cal/filename)
            (cal/'anchor_parameters.txt').write_text(''.join((old/'anchor_parameters.txt').read_text().splitlines(keepends=True)[1:]))
            (cal/'initial.txt').write_text(''.join(f'{i+1} {v:.17g}\n' for i,v in enumerate(native)))
            for variant,settings in VARIANTS.items():
                cq=cal/variant;cq.mkdir();command=(old/'cq/model.cqc').read_text()
                assert command.count('import anchor_parameters')==1 and command.count('iterations=3;')==1
                command=command.replace('import anchor_parameters','import init_parameters << ../initial.txt;\nimport anchor_parameters')
                command=command.replace('iterations=3;',f'iterations={settings["iterations"]};').replace('keeplastests=yes',f'keeplastests={settings["keep"]}')
                (cq/'model.cqc').write_text(command)
    dump(out/'references.json',refs)
    dump(out/'prepared.json',dict(artifact_sha256={str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file()}))
    for cell in CELLS:
        for scheme in WEIGHTS:
            for variant in VARIANTS:
                cq=out/cell/scheme/variant;result=launch(cq,(cq/'model.cqc').read_text(),600)
                assert result['exit_code']==0 and result['end_of_program'],cq
                print(cell,scheme,variant,'completed',flush=True)
    record=review(out,spec,refs,True);previous_check()
    for f,h in sources.items():assert sha(ROOT/f)==h,f
    for f,h in json.loads((out/'prepared.json').read_text())['artifact_sha256'].items():assert sha(out/f)==h,f
    record.update(classification=spec['classification'],scientific_inference_ready=False,qualification_eligible=False,
        arithmetic_and_input_checks_pass=True,protocol_sha256=sha(out/'protocol.json'),
        artifact_sha256={str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file()})
    dump(out/'summary.json',record)
    for cell,runs in record['results'].items():
        for scheme,v in runs.items():print(cell,scheme,{k:dict(iterations=z['iterations'],first_error=z['first_em_error'],root_error=z['root_error'],returned_gradient=z['returned']['gradient']) for k,z in v.items()},flush=True)


if __name__=='__main__':
    try:main()
    except Exception as exc:
        if '--output-dir' in sys.argv and '--audit' not in sys.argv:
            out=Path(sys.argv[sys.argv.index('--output-dir')+1])
            if out.is_dir() and not (out/'summary.json').exists():dump(out/'failure.json',dict(error=repr(exc),traceback=traceback.format_exc()))
        raise
