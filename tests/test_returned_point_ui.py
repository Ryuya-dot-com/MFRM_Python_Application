"""Returned-point checking is opt-in, current-fit-bound and never a refit."""
import json
import io
from pathlib import Path
import zipfile

import numpy as np
import pytest
from streamlit.testing.v1 import AppTest

import streamlit_app as app


def view(lang):
    from unittest.mock import patch
    import pandas as pd
    import streamlit as st
    import streamlit_app as app
    if 'fit' not in st.session_state:
        data=pd.DataFrame([dict(Person=f'{p:03d}',Rater=f'R{r}',Task=f'T{t}',Criterion=f'C{c}',Score=(p+2*r+t+c)%4)
            for p in range(8) for r in range(2) for t in range(2) for c in range(2)])
        st.session_state['fit']=app.mfrm_estimate(data=data,person_col='Person',facet_cols=['Rater','Task','Criterion'],
            score_col='Score',model='PCM',method='MML',step_facet='Criterion',noncenter_facet='Criterion',
            population_prior_sd=1.1,quad_points=31,maxit=40)
        st.session_state['lang']=lang
        st.session_state['evaluations']=0
    fit=st.session_state['fit']
    original=app._analysis_input_assets(fit)
    actual=app.mfrm_loglik_mml_value_grad
    def evaluate(*args,**kwargs):
        st.session_state['evaluations']+=1
        value,gradient=actual(*args,**kwargs)
        if st.session_state.get('flat'):
            gradient=gradient*0. # UI-only seam; mathematical checks use the actual gradient.
        return value,gradient
    def no_refit(*args,**kwargs):
        raise AssertionError('A returned-point check must not optimize or launch an engine')
    with patch.object(app,'mfrm_loglik_mml_value_grad',evaluate), patch.object(app,'minimize',no_refit), \
         patch.object(app._pcm_native,'run_process',no_refit):
        app.show_convergence_section(fit,key_prefix='test')
        if st.session_state.get('two_contexts'):
            app.show_convergence_section(fit,key_prefix='report')
    assert app._analysis_input_assets(fit)==original


@pytest.mark.parametrize('lang', ['en','ja'])
def test_compact_ui_waits_for_click_and_keeps_integration_unassessed(lang):
    at=AppTest.from_function(view,args=(lang,)).run(timeout=60)
    assert not at.exception and not at.success and at.session_state['evaluations']==0
    assert len(at.expander)==1 and not at.expander[0].proto.expanded
    assert len(at.button)==1 and not at.get('download_button')
    assert len(at.get('popover')) == 1
    # All technical tables are inside the initially closed details panel.
    assert len(at.expander[0].dataframe)==len(at.dataframe)
    at.button(key='test_returned_point_check').click().run()
    assert not at.exception and at.session_state['evaluations']==3 and not at.success
    record=at.session_state['_fixed_sd_returned_point_check']
    assert record['stationarity_pass'] is None and record['quadrature_sensitivity_pass'] is None
    assert record['inference_ready'] is False
    assert record['comparison']['accuracy_pass'] is None
    assert not record['comparison']['refitted'] and not record['comparison']['continuous_reference_checked']
    assert len(at.get('download_button')) == 1
    table = at.dataframe[0].value
    assert table.iloc[:, 0].tolist() == [31, 61, 121]
    np.testing.assert_array_equal(table.iloc[:, 1], [r['nll_change'] for r in record['comparison']['rows']])
    saved = json.loads(app._evidence.canonical_json(record))
    assert saved['comparison'] == record['comparison'] and saved['implementation']['source_sha256']
    text='\n'.join(m.value for m in at.markdown)
    assert app._load_locale(lang)['estimation_subsections']['returned_point']['integration'] in text
    assert ('未確認' if lang=='ja' else 'Not checked') in text
    fit=at.session_state['fit'];prepared=app.prepare_fixed_sd_pcm_check(fit)
    actual=app.evaluate_fixed_sd_pcm_check(prepared)
    np.testing.assert_array_equal(record['gradient'],actual['gradient'])
    bundle = app.build_cross_engine_validation_bundle(fit, language=lang)
    assert app._analysis_input_assets(fit)['analysis_snapshot.json'] == bundle['analysis_snapshot.json']
    expected = (Path(app.__file__).parent/'docs'/f'fixed_sd_quadrature_check.{lang}.md').read_bytes()
    with zipfile.ZipFile(io.BytesIO(app.cached_mixed_asset_zip(bundle, app.bytes_mapping_fingerprint(bundle)))) as archive:
        assert archive.read('README_quadrature_check.md') == expected
        assert archive.read('data.csv') == app._analysis_input_assets(fit)['data.csv']
    guide = app._load_locale(lang)['estimation_subsections']['returned_point']
    assert at.get('popover')[0].proto.popover.label == guide['reading_guide_label']
    at.session_state['lang']='ja' if lang=='en' else 'en'
    at.run()
    assert not at.exception and at.session_state['evaluations']==3
    changed_lang=at.session_state['lang']
    assert at.get('popover')[0].proto.popover.label == app._load_locale(changed_lang)['estimation_subsections']['returned_point']['reading_guide_label']
    at.session_state['flat']=True
    at.button(key='test_returned_point_check').click().run()
    assert not at.exception and not at.success and at.session_state['evaluations']==6
    assert at.session_state['_fixed_sd_returned_point_check']['gradient_supnorm']==0.
    assert at.info and not at.warning
    # Alter the same result object while leaving its recorded run fingerprints untouched.
    fit=at.session_state['fit'];fit['config']['quad_points']=61
    at.run()
    assert not at.exception and at.session_state['evaluations']==6
    assert '_fixed_sd_returned_point_check' not in at.session_state
    at.button(key='test_returned_point_check').click().run()
    assert not at.exception and at.warning
    assert 'error' in at.session_state['_fixed_sd_returned_point_check']
    assert not at.get('download_button')


def test_unsupported_or_inconsistent_inputs_clear_cached_check_and_multiple_contexts_work():
    at=AppTest.from_function(view,args=('ja',)).run(timeout=60)
    at.button(key='test_returned_point_check').click().run()
    calls=at.session_state['evaluations']
    fit=at.session_state['fit']
    fit['prep']['data']['score_k']=3-fit['prep']['data']['score_k']
    at.run()
    assert not at.exception and not at.button and at.session_state['evaluations']==calls
    assert '_fixed_sd_returned_point_check' not in at.session_state
    assert any(app._load_locale('ja')['estimation_subsections']['returned_point']['unavailable'] in c.value for c in at.caption)
    fit['prep']['data']['score_k']=fit['prep']['data']['Score'].to_numpy()
    at.session_state['two_contexts']=True
    at.run()
    assert not at.exception and len(at.button)==2
    assert not at.success


def test_input_identity_includes_current_parameters_and_settings_without_truncated_arrays():
    at=AppTest.from_function(view,args=('en',)).run(timeout=60)
    fit=at.session_state['fit'];original=app.prepare_fixed_sd_pcm_check(fit)
    assets=app._analysis_input_assets(fit)
    snapshot=json.loads(assets['analysis_snapshot.json'])
    np.testing.assert_array_equal(snapshot['returned_free_coordinates'],fit['opt'].x)
    fit['opt'].message+=' changed'
    assert app.prepare_fixed_sd_pcm_check(fit)['input_id']!=original['input_id']
    fit['config']['population_prior_sd']=1.6
    assert app.prepare_fixed_sd_pcm_check(fit)['input_id']!=original['input_id']
    with pytest.raises(ValueError,match='Reported NLL'):
        app.evaluate_fixed_sd_pcm_check(app.prepare_fixed_sd_pcm_check(fit))


def test_normal_workspace_preserves_fit_and_comparison_route_across_navigation(monkeypatch):
    """Use the real app entry, actual fit and actual download data, without UI stubs."""
    import pandas as pd

    captured = {}
    download = app.st.download_button
    def capture(*args, **kwargs):
        key = kwargs.get('key')
        if key in ('dl_cross_engine_comparison', 'results_quadrature_download'):
            captured[key] = kwargs['data']
        return download(*args, **kwargs)
    monkeypatch.setattr(app.st, 'download_button', capture)
    at = AppTest.from_string('import streamlit_app as app\napp.main()', default_timeout=60).run()
    at.button(key='onboarding_skip_guide').click().run()
    at.selectbox(key='data_source_class').set_value('paste').run()
    data = pd.DataFrame([dict(Person=f'{p:03d}', Rater=f'R{r}', Task=f'T{t}', Criterion=f'C{c}', Score=(p+2*r+t+c)%4)
        for p in range(8) for r in range(2) for t in range(2) for c in range(2)])
    csv = data.to_csv(index=False)
    at.text_area(key='paste_data_text').set_value(csv).run()
    at.radio(key='facets_mode_workflow_mode').set_value('Advanced controls').run()
    at.radio(key='facets_mode_model_type').set_value('PCM')
    at.radio(key='facets_mode_estimation_method').set_value('MML').run()
    for prefix in ('facets_mode_step_facet_', 'facets_mode_noncenter_facet_'):
        next(w for w in at.selectbox if w.key and w.key.startswith(prefix)).set_value('Criterion')
    at.number_input(key='facets_mode_quad_points').set_value(31)
    at.number_input(key='facets_mode_population_prior_sd').set_value(1.1)
    at.number_input(key='facets_mode_maxit').set_value(40)
    at.selectbox(key='facets_mode_mml_engine').set_value('Direct')
    at.selectbox(key='facets_mode_analysis_depth').set_value('Fast preview').run()
    next(w for w in at.button if 'Run FACETS' in w.label).click().run()
    assert not at.exception
    out = at.session_state['facets_mode_output']
    fit = out['result']
    original = app._analysis_input_assets(fit)
    scores = out['scorefile'].copy(deep=True)
    def no_refit(*args, **kwargs):
        raise AssertionError('Navigation must not refit or launch an external engine')
    monkeypatch.setitem(app.load_core_namespace(), 'mfrm_estimate', no_refit)
    monkeypatch.setattr(app, 'minimize', no_refit)
    monkeypatch.setattr(app._pcm_native, 'run_process', no_refit)

    at.button_group(key='guided_essential_section').set_value('results').run()
    at.button(key='results_returned_point_check').click().run()
    record = json.loads(captured['results_quadrature_download'])
    assert not record['inference_ready'] and record['comparison']['accuracy_pass'] is None
    assert not record['comparison']['refitted']
    at.button_group(key='guided_essential_section').set_value('report_export').run()
    at.button_group(key='guided_export_task').set_value('review').run()
    at.selectbox(key='guided_export_resource').set_value('external_comparison').run()
    for lang in ('en', 'ja'):
        if lang == 'ja':
            at.radio(key='lang').set_value(lang).run()
            workflow = at.radio(key='facets_mode_workflow_mode')
            assert workflow.proto.set_value
            assert workflow.proto.raw_value == app._load_locale(lang)['sidebar_estimation']['workflow_mode_advanced_display']
        assert not at.exception
        assert at.button_group(key='guided_essential_section').value == 'report_export'
        assert at.button_group(key='guided_export_task').value == 'review'
        assert at.selectbox(key='guided_export_resource').value == 'external_comparison'
        assert at.text_area(key='paste_data_text').value == csv
        with zipfile.ZipFile(io.BytesIO(captured['dl_cross_engine_comparison'])) as archive:
            for name, value in original.items():
                assert archive.read(name) == (value.encode() if isinstance(value, str) else value)
            assert archive.read('README_quadrature_check.md') == (Path(app.__file__).parent/'docs'/f'fixed_sd_quadrature_check.{lang}.md').read_bytes()
            report = json.loads(archive.read('comparison_report.json'))
            assert not report['scientific_inference_ready']
            assert {c['status'] for c in report['checks'].values()} == {'not_assessed'}
        at.button_group(key='guided_essential_section').set_value('results').run()
        assert not at.exception
        assert json.loads(captured['results_quadrature_download']) == record
        at.button_group(key='guided_essential_section').set_value('report_export').run()
        assert at.button_group(key='guided_export_task').value == 'review'
        assert at.selectbox(key='guided_export_resource').value == 'external_comparison'
        assert at.session_state['facets_mode_output'] is out
        assert app._analysis_input_assets(fit) == original
        pd.testing.assert_frame_equal(out['scorefile'], scores)

    at.button(key='mfrm_home_open').click().run()
    at.button(key='mfrm_home_resume').click().run()
    assert not at.exception
    assert at.button_group(key='guided_export_task').value == 'review'
    assert at.selectbox(key='guided_export_resource').value == 'external_comparison'
    assert at.session_state['facets_mode_output'] is out
    assert app._analysis_input_assets(fit) == original

    # A real input edit must still withhold the old result, despite retained navigation.
    at.text_area(key='paste_data_text').set_value(csv.replace('000,R0,T0,C0,0', '000,R0,T0,C0,1')).run()
    assert not at.exception
    assert not at.button_group
    assert not any(w.proto.id.endswith('-dl_cross_engine_comparison') for w in at.get('download_button'))
    assert at.session_state['facets_mode_output'] is out
