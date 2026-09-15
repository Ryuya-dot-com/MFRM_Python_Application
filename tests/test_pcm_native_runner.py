"""A native runner must preserve the model or stop before launching an engine."""
from copy import deepcopy
import hashlib
import io
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

import streamlit_app as app
from mfrm_app.legacy_compat import pcm_native_runner as native


@pytest.fixture(scope='module')
def supported_fit():
    data = pd.DataFrame([{'Person':f'{p:03d}', 'Rater':f'R{r}', 'Task':f'T{t}',
        'Criterion':f'C{c}', 'Score':(p+2*r+t+c)%4}
        for p in range(16) for r in range(2) for t in range(2) for c in range(2)])
    return app.mfrm_estimate(data=data,person_col='Person',facet_cols=['Rater','Task','Criterion'],
        score_col='Score',model='PCM',method='MML',step_facet='Criterion',noncenter_facet='Criterion',
        positive_facets=['Rater'],population_prior_sd=1.1,quad_points=13,maxit=60)


def test_supported_export_preserves_logits_gradient_labels_and_hashes(supported_fit, tmp_path):
    bundle = app.build_cross_engine_validation_bundle(supported_fit)
    assert json.loads(bundle['native_runner_status.json'])['available']
    for name, value in bundle.items():
        (tmp_path/name).write_bytes(value.encode() if isinstance(value,str) else value)
    job = native.checked_bundle(tmp_path)
    assert job['person_labels'][0] == '000' and job['pid'][0] == 'P000000'
    job.update(bound=12., nodes=121)
    z = np.array(job['xsi']); value = native.evaluate(job,z)
    fd = [(native.evaluate(job,z+d)['raw_nll']-native.evaluate(job,z-d)['raw_nll'])/2e-5 for d in np.eye(8)*1e-5]
    np.testing.assert_allclose(value['gradient'],fd,atol=2e-8,rtol=0)
    params = supported_fit['params']; signs = supported_fit['config']['facet_signs']
    eta = np.array([.7+signs['Rater']*params['facets']['Rater'][r]+signs['Task']*params['facets']['Task'][t]
                   +signs['Criterion']*params['facets']['Criterion'][c] for r in range(2) for t in range(2) for c in range(2)])
    app_p = app.category_prob_pcm(eta,np.column_stack((np.zeros(2),np.cumsum(params['steps_mat'],axis=1))),np.array([0,1]*4))
    logits = native.native_design() @ z + .7*np.arange(4)
    p = np.exp(logits-logits.max(axis=1,keepdims=True));p /= p.sum(axis=1,keepdims=True)
    np.testing.assert_allclose(app_p,p,atol=1e-14,rtol=0)
    # Tampering cannot keep the original manifest or silently attach its run ID.
    (tmp_path/'data.csv').write_bytes(bundle['data.csv']+b'\n')
    with pytest.raises(ValueError,match='Input hash mismatch'): native.checked_bundle(tmp_path)
    report = json.loads(bundle['comparison_report.json'])
    report['details']['input_sha256']['data.csv'] = hashlib.sha256((tmp_path/'data.csv').read_bytes()).hexdigest()
    (tmp_path/'comparison_report.json').write_text(json.dumps(report))
    with pytest.raises(ValueError,match='manifest digest'): native.checked_bundle(tmp_path)


def test_unsupported_models_constraints_weights_and_duplicate_cells_stop(supported_fit, tmp_path, monkeypatch):
    bundle = app.build_cross_engine_validation_bundle(supported_fit)
    source = json.loads(bundle['analysis_snapshot.json'])
    for key,value in [('model','RSM'),('method','JMLE'),('estimate_population_sd',True),
                      ('noncenter_facet','Person'),('facet_regularization_enabled',True),('dummy_facets',['Task'])]:
        changed = deepcopy(source); changed['config'][key] = value
        with pytest.raises(ValueError): native.prepare_inputs(bundle['data.csv'],json.dumps(changed))
    for spec_change in [dict(anchors=[0.,{'nonfinite':'nan'}]),dict(groups=['A','A'],group_values={'A':0.})]:
        changed = deepcopy(source);changed['config']['facet_specs']['Rater'].update(spec_change)
        with pytest.raises(ValueError,match='Anchors'): native.prepare_inputs(bundle['data.csv'],json.dumps(changed))
    changed = deepcopy(source);changed['params']['steps_mat'][0][0] += .1
    with pytest.raises(ValueError,match='steps disagree'):native.prepare_inputs(bundle['data.csv'],json.dumps(changed))
    data = pd.read_csv(io.BytesIO(bundle['data.csv']),dtype={'Person':str})
    for weight in [0., .5, 2., np.nan]:
        changed = data.copy();changed.loc[0,'Weight'] = weight
        with pytest.raises(ValueError):native.prepare_inputs(app.to_csv_bytes(changed),bundle['analysis_snapshot.json'])
    with pytest.raises(ValueError,match='Repeated'):native.prepare_inputs(app.to_csv_bytes(pd.concat([data,data.iloc[:1]])),bundle['analysis_snapshot.json'])
    sparse = data.loc[~((data.Rater=='R0') & (data.Task=='T0') & (data.Criterion=='C0') & (data.Score==3))]
    with pytest.raises(ValueError,match='All four categories'):native.prepare_inputs(app.to_csv_bytes(sparse),bundle['analysis_snapshot.json'])
    # Missing cells are retained, never turned into zero scores.
    missing = native.prepare_inputs(app.to_csv_bytes(data.iloc[1:]),bundle['analysis_snapshot.json'])
    assert sum(v is None for row in missing['responses'] for v in row)==1
    for name,value in bundle.items():(tmp_path/name).write_bytes(value.encode() if isinstance(value,str) else value)
    (tmp_path/'data.csv').write_bytes(b'changed')
    monkeypatch.setattr(native,'run_process',lambda *a,**k:pytest.fail('Engine must not run'))
    monkeypatch.setattr(sys,'argv',['run_pcm_native.py','--bundle',str(tmp_path),'--output',str(tmp_path/'out'),'--conquest','unused'])
    with pytest.raises(ValueError,match='Input hash mismatch'):native.main()
    assert not (tmp_path/'out').exists()
