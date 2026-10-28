"""Development controls for conditional scalar bias; not a coverage study."""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import subprocess
import sys
import zipfile

import numpy as np
import pandas as pd
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import streamlit_app as app


def clean(value):
    if isinstance(value, dict):
        return {k: clean(v) for k,v in value.items()}
    if isinstance(value, (list,tuple,np.ndarray)):
        return [clean(v) for v in value]
    if isinstance(value, np.generic):
        return clean(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def save(path, value):
    path.write_text(json.dumps(clean(value),indent=2,allow_nan=False)+'\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    parser.add_argument('--reference-recheck-from', type=Path)
    args = parser.parse_args(); out=args.output
    out.mkdir(parents=True,exist_ok=False)
    limits = dict(probability=1e-12,nll=1e-10,score=1e-10,information=1e-10,
                  gradient_fd=2e-7,information_fd=2e-5,estimate=1e-8,se=1e-8,
                  lr=1e-8,p=1e-8,ci=1e-7,adjustment=1e-14)
    cases=[dict(id=f'{model}_{kind}',model=model,kind=kind,bound=.02 if kind=='range' else 10.)
           for model in ['RSM','PCM','GPCM']
           for kind in ['unit','weighted','extreme_min','extreme_max','range']]
    parent = args.reference_recheck_from
    if parent:
        manifest = json.loads((parent/'manifest.json').read_text())
        for name, digest in manifest.items():
            if hashlib.sha256((parent/name).read_bytes()).hexdigest() != digest:
                raise ValueError(f"Changed parent artifact: {name}")
        original = json.loads((parent/'protocol.json').read_text())
        if limits != original['limits']:
            raise ValueError("Recheck must preserve every original tolerance")
    save(out/'protocol.json',dict(classification='POST_RESULT_REFERENCE_STENCIL_RECHECK' if parent else 'DEVELOPMENT_ONLY',
        parent_manifest_sha256=hashlib.sha256((parent/'manifest.json').read_bytes()).hexdigest() if parent else None,
        scientific_inference_ready=False,
        scope='Scalar likelihood, derivatives, roots, intervals and multiplicity; not sampling coverage',
        limits=limits,cases=cases))
    with zipfile.ZipFile(out/'source.zip','w',zipfile.ZIP_DEFLATED) as archive:
        for p in [Path(__file__),ROOT/'validation/bias_conditional_reference.R',ROOT/'streamlit_app.py',
                  ROOT/'mfrm_app/output_qualification.py']:
            archive.write(p,str(p.relative_to(ROOT)))
    save(out/'runtime.json',dict(python=sys.version,packages={d.metadata['Name']:d.version
                               for d in importlib.metadata.distributions()}))
    results=[]
    for cell in cases:
        n=24; row=np.arange(n); group=row%3
        eta=np.linspace(-1.2,1.4,n)
        steps=np.array([-.7,.1,.6])[None,:]+(group-1)[:,None]*np.array([-.2,.05,.15])
        if cell['model']=='RSM': steps=np.tile([-.7,.1,.6],(n,1))
        slopes=np.array([.5,1.,2.])[group] if cell['model']=='GPCM' else np.ones(n)
        score=(row*3+row//4)%4
        if cell['kind']=='extreme_min': score=np.zeros(n,dtype=int)
        if cell['kind']=='extreme_max': score=np.full(n,3,dtype=int)
        if cell['kind']=='range': score=np.full(n,2,dtype=int)
        weights=.5+(row%4)/3 if cell['kind']=='weighted' else np.ones(n)
        data=pd.DataFrame(dict(eta=eta,score=score,slope=slopes,weight=weights,
                               step1=steps[:,0],step2=steps[:,1],step3=steps[:,2]))
        cell['file']=cell['id']+'.csv'; data.to_csv(out/cell['file'],index=False,float_format='%.17g')
        cum=np.column_stack([np.zeros(n),np.cumsum(steps,axis=1)])
        if cell['model']=='RSM':
            nll=lambda b: -app.loglik_rsm(eta+b,score,cum[0],weight=weights)
            probability=lambda b: app.category_prob_rsm(eta+b,cum[0])
        elif cell['model']=='PCM':
            nll=lambda b: -app.loglik_pcm(eta+b,score,cum,row,weight=weights)
            probability=lambda b: app.category_prob_pcm(eta+b,cum,row)
        else:
            nll=lambda b: -app.loglik_gpcm(eta+b,score,cum,row,slopes,row,weight=weights)
            probability=lambda b: app.category_prob_gpcm(eta+b,cum,row,slopes,row)
        fit=app._fit_conditional_bias_cell(nll,probability,score,slopes,weights,cell['bound'])
        points=[]
        for b in [-.4,0.,.6]:
            prob=probability(b); mean=prob@np.arange(4)
            var=np.sum(prob*(np.arange(4)[None,:]-mean[:,None])**2,axis=1)
            points.append(dict(nll=nll(b),score=np.sum(weights*slopes*(mean-score)),
                               info=np.sum(weights*slopes**2*var),probability=prob))
        results.append(dict(id=cell['id'],fit=fit,points=points))
    families=[[.01,None,.03,.04],[.02,.03,None],[None,None],[.01,.01,1.,0.]]
    save(out/'input.json',dict(cases=cases,families=families));save(out/'python.json',results)
    if parent:
        for name in ['input.json']+[cell['file'] for cell in cases]:
            if (out/name).read_bytes() != (parent/name).read_bytes():
                raise ValueError(f"Recheck inputs changed: {name}")
    run=subprocess.run(['/usr/local/bin/Rscript',str(ROOT/'validation/bias_conditional_reference.R'),
                        str(out/'input.json'),str(out/'R.json')],capture_output=True,text=True)
    (out/'R.log').write_text(run.stdout+run.stderr)
    if run.returncode: raise RuntimeError('Independent R reference failed; see R.log')
    reference=json.loads((out/'R.json').read_text());checks=[]
    for py, r in zip(results,reference['cases'],strict=True):
        error={k:0. for k in ['nll','score','information','probability','gradient_fd','information_fd']}
        for pp,rp in zip(py['points'],r['points'],strict=True):
            for name,source in [('nll','nll'),('score','score'),('information','info'),('probability','probability')]:
                error[name]=max(error[name],float(np.max(abs(np.asarray(pp[source])-np.asarray(rp[source])))))
            error['gradient_fd']=max(error['gradient_fd'],abs(rp['gradient_fd']-pp['score']))
            error['information_fd']=max(error['information_fd'],float(np.max(abs(np.asarray(rp['information_fd'])-pp['info']))))
        same_status=py['fit']['status']==r['status']
        if r['status']=='interior':
            for name in ['estimate','se','lr','p']:
                if name in r: error[name]=abs(py['fit'][name]-r[name])
            if 'ci' in r:
                error['ci']=float(max(abs(np.array([py['fit']['ci'][s] for s in ['lower','upper']])-r['ci'])))
                same_status &= py['fit']['ci']['status']==r['ci_status']
        checks.append(dict(id=py['id'],status_match=same_status,errors=error,
                           passed=same_status and all(v<limits[k] for k,v in error.items())))
    adjustments=[]
    for family,r in zip(families,reference['adjustments'],strict=True):
        p=np.asarray(family,dtype=float)
        for method in ['holm','BH']:
            py=app.adjust_p_values(p,method);rv=np.asarray(r[method],dtype=float)
            passed=np.array_equal(np.isnan(py),np.isnan(rv)) and np.allclose(py,rv,rtol=0,atol=limits['adjustment'],equal_nan=True)
            adjustments.append(dict(method=method,p=p,python=py,R=rv,passed=passed))
    summary=dict(classification='DEVELOPMENT_ONLY',scientific_inference_ready=False,
                 checks=checks,adjustments=adjustments,all_checks_pass=all(c['passed'] for c in checks+adjustments))
    save(out/'summary.json',summary)
    save(out/'manifest.json',{p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                             for p in sorted(out.iterdir()) if p.is_file()})
    print(json.dumps(dict(cases=len(checks),adjustments=len(adjustments),all_checks_pass=summary['all_checks_pass'])))
    return 0 if summary['all_checks_pass'] else 1


if __name__=='__main__': raise SystemExit(main())
