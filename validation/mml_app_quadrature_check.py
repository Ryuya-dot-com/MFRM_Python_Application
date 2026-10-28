"""Bounded app-entrypoint regression, not convergence or sampling validation.

Compare the retained prior app with the new quadrature on six development
fits. Four EM iterations intentionally exercise updates and final outputs;
neither an iteration-limit result nor numerical agreement grants inference.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import streamlit_app as app
from tests.test_mml_grouped_likelihood import problem
from validation.mml_population_study import json_bytes


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('before',type=Path);parser.add_argument('output',type=Path)
    args=parser.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    spec=importlib.util.spec_from_file_location('retained_prior_quadrature_app',args.before.resolve())
    old=importlib.util.module_from_spec(spec);sys.modules[spec.name]=old;spec.loader.exec_module(old)
    def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
    plan=dict(classification='DEVELOPMENT_APP_ENTRYPOINT_REGRESSION',scientific_inference_ready=False,
        cases=[dict(model=m,estimate_population_sd=s) for m in ['RSM','PCM','GPCM'] for s in [False,True]],
        q=15,maxit=4,reltol=1e-8,point_atol=2e-5,nll_per_person_atol=1e-8,
        prior_app_sha256=sha(args.before),candidate_app_sha256=sha(app.__file__),script_sha256=sha(__file__))
    (out/'plan.json').write_bytes(json_bytes(plan));records=[]
    for case in plan['cases']:
        model=case['model'];data,_,_,_,_=problem(model,'unbalanced',duplicate=True)
        name=model+('_free' if case['estimate_population_sd'] else '_fixed')
        (out/(name+'.csv')).write_text(data.to_csv(index=False))
        fits={}
        for label,engine in [('prior',old),('candidate',app)]:
            result=engine.mfrm_estimate(data,person_col='Person',facet_cols=['Rater','Task','Criterion'],
                score_col='Score',rating_min=0,rating_max=3,model=model,method='MML',
                noncenter_facet='Criterion',step_facet=None if model=='RSM' else 'Criterion',
                quad_points=plan['q'],population_prior_sd=1.2,estimate_population_sd=case['estimate_population_sd'],
                maxit=plan['maxit'],reltol=plan['reltol'],compute_plausible_values=True,n_plausible_values=3)
            value=dict(point=result['opt'].x,sigma=result['config']['population_prior_sd'],
                nll=float(result['opt'].raw_fun),optimizer_success=bool(result['opt'].success),
                optimizer_status=getattr(result['opt'],'status',None),
                optimizer_message=str(getattr(result['opt'],'message','')),iterations=int(result['opt'].nit),
                person=result['facets']['person'][['Estimate','SD']].to_numpy(),
                posterior=result['posterior']['scores'].to_dict('records'),
                plausible_values=result['posterior']['plausible_values'].to_dict('records'),
                summary=json.loads(result['summary'].to_json(orient='records')))
            fits[label]=value;(out/(name+'_'+label+'.json')).write_bytes(json_bytes(value))
        a,b=fits['prior'],fits['candidate']
        checks={key:a[key]==b[key] for key in ['optimizer_success','optimizer_status','iterations']}
        errors=dict(point=float(np.max(abs(a['point']-b['point']))),sigma=abs(a['sigma']-b['sigma']),
            person=float(np.max(abs(a['person']-b['person']))),nll_per_person=abs(a['nll']-b['nll'])/data.Person.nunique())
        checks.update({key:np.isfinite(error) and error<=plan['nll_per_person_atol' if key=='nll_per_person' else 'point_atol']
                       for key,error in errors.items()})
        record=dict(case=case,errors=errors,checks=checks,all_checks_pass=all(checks.values()),
                    inference_claim=False,optimizer_success={k:v['optimizer_success'] for k,v in fits.items()})
        records.append(record);print(name,record['all_checks_pass'],flush=True)
    summary=dict(classification=plan['classification'],scientific_inference_ready=False,
        records=records,all_checks_pass=all(r['all_checks_pass'] for r in records))
    (out/'summary.json').write_bytes(json_bytes(summary))
    return 0 if summary['all_checks_pass'] else 1


if __name__=='__main__':raise SystemExit(main())
