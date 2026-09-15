"""Current-fit worksheet identity and its opt-in Streamlit entry."""
from copy import deepcopy
import hashlib
import io
import json
from types import SimpleNamespace
import zipfile

import numpy as np
import pandas as pd
import pytest
from streamlit.testing.v1 import AppTest

import streamlit_app as app


def current_fit():
    return dict(prep=dict(data=pd.DataFrame({'Person':['001','002'], 'Item':['A','A'],
            'Score':[1,2], 'Weight':[1.,1.]}), levels={'Person':['001','002'], 'Item':['A']}),
        config=dict(model='PCM', method='MML', facet_names=['Item'], step_facet='Item',
            quad_points=31, analysis_config_fingerprint='stale-config-id', input_data_fingerprint='stale-input-id',
            population_model={'design':np.zeros((1200,2))},
            anchor_audit={'valid_anchors':pd.DataFrame({'Facet':['Item'],'Level':['A'],'Anchor':[0.]})}),
        params={'steps_mat':np.array([[-1.,1.]])}, opt=SimpleNamespace(x=np.array([.2]),fun=3.,success=True,message='Converged'),
        facets={'person':pd.DataFrame({'Person':['001','002'],'Measure':[.1,.2]})})


def test_worksheet_never_promotes_claims_and_binds_complete_current_content():
    fit = current_fit()
    fit['config']['series'] = pd.Series([.2], index=['A'], name='Named values')
    fit['comparison_report'] = {'scientific_inference_ready':True, 'status':'within_targets'}
    assets = app.build_cross_engine_validation_bundle(fit)
    assert b'MIT License' in assets['LICENSE'] and 'software_citation.md' in assets
    report = json.loads(assets['comparison_report.json'])
    assert json.loads(assets['analysis_snapshot.json'])['config']['series']=={
        'index':['A'], 'data':[.2], 'name':'Named values'}
    assert not report['scientific_inference_ready'] and not report['qualification_eligible']
    assert report['details']['artifact_kind']=='unassessed_template'
    assert not report['details']['external_comparison_run']
    assert {c['status'] for c in report['checks'].values()}=={'not_assessed'}
    assert all(c['metrics']==[] and c['missing_evidence'] for c in report['checks'].values())
    measurements = pd.read_csv(io.BytesIO(assets['comparison_metrics.csv']))
    assert measurements.empty and 'Value' in measurements.columns
    assert b'No qualification of inference' in assets['comparison_report.html']
    for name, digest in report['details']['input_sha256'].items():
        value=assets[name]; value=value.encode() if isinstance(value,str) else value
        assert hashlib.sha256(value).hexdigest()==digest
    with zipfile.ZipFile(io.BytesIO(app._exports.build_mixed_asset_zip(assets))) as z:
        assert set(z.namelist())==set(assets)
        assert z.read('comparison_report.json').decode()==assets['comparison_report.json']
    identity=report['details']['analysis_identity']['analysis_id']
    # Fresh content, even with unchanged dimensions and stale recorded run IDs.
    variants=[deepcopy(fit) for _ in range(7)]
    variants[0]['prep']['data'].loc[0,'Score']=2
    variants[1]['prep']['data'].loc[0,'Weight']=2.
    variants[2]['config']['quad_points']=61
    variants[3]['opt'].x[0]=.3
    variants[4]['config']['population_model']['design'][600,0]=.12345678901234568
    variants[5]['config']['anchor_audit']['valid_anchors'].loc[0,'Anchor']=.1
    variants[6]['params']['steps_mat'][0,0]=-1.1
    for changed in variants:
        bundle=app.build_cross_engine_validation_bundle(changed)
        updated=json.loads(bundle['comparison_report.json'])
        assert updated['details']['analysis_identity']['analysis_id']!=identity
        assert updated['details']['input_bundle_sha256']!=report['details']['input_bundle_sha256']
    design=json.loads(app.build_cross_engine_validation_bundle(variants[4])['analysis_snapshot.json'])['config']['population_model']['design']
    assert len(design)==1200 and design[600][0]==.12345678901234568
    ja=json.loads(app.build_cross_engine_validation_bundle(fit,language='ja')['comparison_report.json'])
    assert ja['details']['analysis_identity']['analysis_id']==identity
    assert ja['title']!=report['title']
    fit['config']['population_model']['design'][600,0]=np.nan
    snapshot=json.loads(app.build_cross_engine_validation_bundle(fit)['analysis_snapshot.json'])
    assert snapshot['config']['population_model']['design'][600][0]=={'nonfinite':'nan'}
    fit['config']['unsupported']=object()
    with pytest.raises(ValueError,match='Unsupported contract value'):
        app.build_cross_engine_validation_bundle(fit)


def test_comparison_entry_is_opt_in_bilingual_and_handles_unexportable_fit():
    def view(lang):
        from unittest.mock import patch
        import pandas as pd
        import streamlit as st
        import streamlit_app as app
        st.session_state['lang']=lang
        fit={'prep':{'data':pd.DataFrame({'Person':['001'], 'Score':[1]})}, 'config':{'method':'JMLE'}}
        if st.session_state.get('broken'):
            fit['config']['unsupported']=object()
        with patch.object(app,'build_publication_gate_summary',return_value=pd.DataFrame()), \
             patch.object(app,'_render_publication_document_section'):
            app._render_guided_report_export_section(fit,{},{},pd.DataFrame(),pd.DataFrame(),None,None,generate_figures=False)
    for lang in ['en','ja']:
        at=AppTest.from_function(view,args=(lang,)).run(timeout=45)
        assert not at.exception and not at.selectbox and not at.get('download_button')
        at.button_group(key='guided_export_task').set_value('review').run()
        assert not at.exception and not at.get('download_button')
        at.selectbox(key='guided_export_resource').select('external_comparison').run()
        assert not at.exception and not at.error and len(at.get('download_button'))==1
        label='未評価' if lang=='ja' else 'Not assessed'
        assert sum(label in m.value for m in at.markdown)==4
        assert not at.expander[-1].proto.expanded
        assert not at.get('file_uploader') and not at.dataframe
        assert not any(x.value.startswith('guided.') for x in at.caption)
        at.session_state['broken']=True
        at.run()
        assert not at.exception and at.error and not at.get('download_button')


def test_supported_native_download_keeps_all_four_questions_unassessed():
    def view(lang):
        from unittest.mock import patch
        import pandas as pd
        import streamlit as st
        import streamlit_app as app
        st.session_state['lang']=lang
        data=pd.DataFrame([{'Person':str(p),'Rater':str(r),'Task':str(t),'Criterion':str(c),
            'Score':(p+r+t+c)%4} for p in range(8) for r in range(2) for t in range(2) for c in range(2)])
        fit=app.mfrm_estimate(data=data,person_col='Person',facet_cols=['Rater','Task','Criterion'],
            score_col='Score',model='PCM',method='MML',step_facet='Criterion',noncenter_facet='Criterion',
            maxit=5,quad_points=13)
        with patch.object(app._pcm_native,'run_process') as external:
            app._render_cross_engine_comparison_section(fit)
            external.assert_not_called()
    for lang in ['en','ja']:
        at=AppTest.from_function(view,args=(lang,)).run(timeout=45)
        assert not at.exception and not at.error and len(at.get('download_button'))==1
        assert sum(('未評価' if lang=='ja' else 'Not assessed') in m.value for m in at.markdown)==4
        assert any(app._load_locale(lang)['guided']['comparison_native_available'] in c.value for c in at.caption)
        assert not at.expander[0].proto.expanded
