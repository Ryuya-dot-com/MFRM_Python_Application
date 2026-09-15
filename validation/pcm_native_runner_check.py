#!/usr/bin/env python3
"""Observed export/runner checks; preserve inputs, every native call and failures."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import streamlit_app as app
from mfrm_app.legacy_compat import pcm_native_runner as native


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    out=parser.parse_args().output.resolve();out.mkdir(parents=True,exist_ok=False)
    prior=ROOT/'validation/generated/mml_independent_observation_designs_20260914'
    assert sha(prior/'summary.json')=='106409fca18611bfdfbfdc7809638f849cbb289f6503e023110131175506f6d8'
    old=json.loads((prior/'summary.json').read_text())
    source=prior/'ordinary_full_input.json'
    assert sha(source)==old['artifact_sha256'][source.name]
    sources=['streamlit_app.py','mfrm_app/legacy_compat/pcm_native_runner.py','mfrm_app/legacy_compat/pcm_native_tam.R',
             'tests/test_pcm_native_runner.py','validation/pcm_native_runner_check.py','locales/en.json','locales/ja.json']
    protocol=dict(classification='OBSERVED_EXPORT_IMPLEMENTATION_CHECK',scientific_inference_ready=False,
        question='Does the downloaded restricted PCM runner preserve observations/model/anchors and gate refits on both fixed-calibration checks?',
        data='Retained 32-person ordinary synthetic responses; complete and exactly the first response removed. No missingness stress or inference study.',
        source_data_sha256=sha(source),source_sha256={name:sha(ROOT/name) for name in sources},
        app_fit=dict(model='PCM',method='MML',mml_engine='EM',population_prior_sd=1.1,quad_points=31,maxit=200,reltol=1e-8,
                     noncenter_facet='Criterion',no_regression=True,unit_weights=True),
        native_runs=[['complete',401,20,True],['one_missing',401,20,True],['complete',241,12,False],['complete',201,20,False]],
        design='Keep 0.1 spacing for the range contrast, and keep +/-20 for the 0.1 versus 0.2 spacing contrast.',
        post_hoc_question='The initial implementation probe showed a gradient at the app returned point; reevaluate that same point at GH31/61/121 to separate integration from stopping.',
        limits='Runner targets are implementation tolerances, not statistical equivalence margins. No universal score, gradient or inference pass.')
    native.write_json(out/'protocol.json',protocol)
    raw=pd.DataFrame(json.loads(source.read_text())['full_responses']);bundles={};gh={}
    for name,data in [('complete',raw),('one_missing',raw.iloc[1:])]:
        fit=app.mfrm_estimate(data=data,person_col='Person',facet_cols=['Rater','Task','Criterion'],score_col='Score',
            rating_min=0,rating_max=3,model='PCM',method='MML',step_facet='Criterion',noncenter_facet='Criterion',
            mml_engine='EM',population_prior_sd=1.1,maxit=200,quad_points=31,reltol=1e-8)
        assets=app.build_cross_engine_validation_bundle(fit)
        assert json.loads(assets['native_runner_status.json'])['available']
        folder=out/f'{name}_bundle';folder.mkdir();bundles[name]=folder
        for filename,value in assets.items(): (folder/filename).write_bytes(value.encode() if isinstance(value,str) else value)
        idx=app.build_indices(fit['prep'],step_facet='Criterion');sizes=app.build_param_sizes(fit['config'])
        gh[name]={}
        for q in [31,61,121]:
            value,gradient=app.mfrm_loglik_mml_value_grad(fit['opt'].x,idx,fit['config'],sizes,app.make_mml_quadrature(fit['config'],q))
            gh[name][str(q)]=dict(nll=float(value),gradient=np.asarray(gradient).tolist(),gradient_supnorm=float(max(abs(gradient))))
        gh[name]['app_message']=str(fit['opt'].message)
        gh[name]['app_success']=bool(fit['opt'].success)
    native.write_json(out/'app_same_point_gh.json',gh)
    results={}
    for name,nodes,bound,refit in protocol['native_runs']:
        key=f'{name}_q{nodes}_b{bound}';folder=out/key
        command=['python3','-B',str(bundles[name]/'run_pcm_native.py'),'--output',str(folder),
                 '--conquest','/Applications/ConQuest/ConQuest','--nodes',str(nodes),'--bound',str(bound)]
        if refit:command.append('--refit')
        with (out/f'{key}.log').open('x') as log:
            process=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT)
        assert process.returncode==0, out/f'{key}.log'
        review=json.loads((folder/'execution_review.json').read_text());assert review['checks_complete']
        results[key]=review['results'];print(key,'complete',flush=True)
    native.write_json(out/'summary.json',dict(protocol_sha256=sha(out/'protocol.json'),results=results,
        scientific_inference_ready=False,qualification_eligible=False,
        artifact_sha256={str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file()}))


if __name__=='__main__':main()
