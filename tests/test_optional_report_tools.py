"""Optional calculations stay outside reading tables and retain their settings."""
import pytest
from streamlit.testing.v1 import AppTest


@pytest.mark.parametrize('lang', ['en', 'ja'])
def test_only_selected_optional_tool_renders_and_selection_survives_navigation(lang):
    def view(lang):
        from contextlib import ExitStack
        from unittest.mock import patch
        import pandas as pd
        import streamlit as st
        import streamlit_app as app
        st.session_state.setdefault('lang', lang)
        st.radio('Language', ['en', 'ja'], key='lang')
        show = st.checkbox('Show tools', value=True, key='show_tools')
        st.session_state['rendered'] = []
        names = {'model_choice': 'render_model_choice_guidance', 'g_d': 'render_g_d_study_section',
                 'design_network': 'render_design_network_section', 'rater_severity': 'render_rater_severity_section',
                 'rater_halo': 'render_rater_halo_section', 'numerical_check': '_render_fixed_sd_returned_point',
                 'recovery': 'render_parameter_recovery_simulation'}
        def record(key):
            return lambda *a, **kw: st.session_state['rendered'].append(key)
        with ExitStack() as stack:
            for key, name in names.items():
                stack.enter_context(patch.object(app, name, record(key)))
            stack.enter_context(patch.object(app, 'mfrm_estimate', side_effect=AssertionError('No refit on navigation')))
            stack.enter_context(patch.object(app, 'qualified_convergence_table', return_value=pd.DataFrame()))
            if show:
                app._render_additional_report_analyses({'config': {'method': 'MML'}})

    at = AppTest.from_function(view, args=(lang,)).run(timeout=45)
    assert not at.exception and at.session_state['rendered'] == []
    assert at.selectbox(key='report_additional_tool').value == 'none'
    for key in ('model_choice', 'g_d', 'design_network', 'rater_severity', 'rater_halo', 'numerical_check', 'recovery'):
        at.selectbox(key='report_additional_tool').set_value(key).run()
        assert not at.exception and at.session_state['rendered'] == [key]
        assert not any('report_tools.' in c.value for c in at.caption)
    at.radio(key='lang').set_value('ja' if lang == 'en' else 'en').run()
    assert not at.exception and at.session_state['rendered'] == ['recovery']
    at.checkbox(key='show_tools').uncheck().run()
    assert not at.exception and at.session_state['rendered'] == []
    at.checkbox(key='show_tools').check().run()
    assert not at.exception and at.session_state['rendered'] == ['recovery']
    assert at.selectbox(key='report_additional_tool').value == 'recovery'


def test_reading_tables_does_not_open_optional_calculations():
    def view():
        from contextlib import ExitStack
        from unittest.mock import patch
        import streamlit as st
        import streamlit_app as app
        if 'result' not in st.session_state:
            result = app.mfrm_estimate(app._make_self_test_gradient_data(), 'Person', ['Rater', 'Task'], 'Score')
            st.session_state['result'] = result
            st.session_state['diagnostics'] = app.mfrm_diagnostics(result, compute_pca=False)
        with ExitStack() as stack:
            for name in ('render_model_choice_guidance', 'render_g_d_study_section', 'render_design_network_section',
                         'render_rater_severity_section', 'render_rater_halo_section', 'render_parameter_recovery_simulation',
                         'mfrm_estimate', 'evaluate_fixed_sd_quadrature_check'):
                stack.enter_context(patch.object(app, name, side_effect=AssertionError(f'Tables invoked {name}')))
            app._render_report_tables(st.session_state['result'], st.session_state['diagnostics'])
    at = AppTest.from_function(view).run(timeout=60)
    assert not at.exception and at.dataframe
    assert not at.button
    assert at.selectbox(key='apa_presets_table_picker')
    assert any('Estimate' in frame.value.columns for frame in at.dataframe)
    at.run()
    assert not at.exception


@pytest.mark.parametrize('lang', ['en', 'ja'])
def test_fixed_sd_tables_keep_existing_evidence_without_a_run_button(lang):
    def view(lang):
        from types import SimpleNamespace
        from unittest.mock import patch
        import pandas as pd
        import streamlit as st
        import streamlit_app as app
        st.session_state['lang'] = lang
        allow = st.checkbox('Allow check', key='allow')
        st.session_state['_fixed_sd_returned_point_check'] = {
            'input_id': 'same-input', 'gradient_supnorm': 0.4, 'gradient_reference': 0.01,
            'raw_nll': 32.7, 'quadrature': {'points': 31, 'sigma': 1.0},
        }
        with patch.object(app, 'prepare_fixed_sd_check', return_value={'input_id': 'same-input'}), \
             patch.object(app, 'evaluate_fixed_sd_quadrature_check', side_effect=AssertionError('No automatic calculation')):
            app._render_fixed_sd_returned_point(
                {'config': {'method': 'MML'}, 'opt': SimpleNamespace(success=True)},
                pd.DataFrame([{'Converged': True}]), key_prefix='report', allow_check=allow)
        st.session_state['expected_warning'] = app.t('estimation_subsections.returned_point.remaining_notice')
    at = AppTest.from_function(view, args=(lang,)).run(timeout=45)
    assert not at.exception and not at.button
    assert at.session_state['expected_warning'] in [w.value for w in at.warning]
    assert any('4.00e-01' in item.value for item in at.markdown)
    assert any('Additional analyses' in c.value or '追加分析' in c.value for c in at.caption)
    at.checkbox(key='allow').check().run()
    assert not at.exception and len(at.button) == 1
    assert at.session_state['expected_warning'] in [w.value for w in at.warning]


def test_simulation_settings_and_cached_result_survive_hidden_panel_and_language_change():
    def view():
        from unittest.mock import patch
        import streamlit as st
        import streamlit_app as app
        st.radio('Language', ['en', 'ja'], key='lang')
        show = st.checkbox('Show tools', value=True, key='show_tools')
        st.session_state.setdefault('simulation_calls', 0)
        def run(**settings):
            st.session_state['simulation_calls'] += 1
            st.session_state['submitted_settings'] = settings
            return {'available': False, 'reason': 'retained-test-result'}
        with patch.object(app, 'evaluate_parameter_recovery', run), \
             patch.object(app, 'mfrm_estimate', side_effect=AssertionError('No current-data refit')):
            if show:
                app._render_additional_report_analyses({'config': {'method': 'JMLE'}})
    at = AppTest.from_function(view).run(timeout=45)
    at.selectbox(key='report_additional_tool').set_value('recovery').run()
    assert not at.exception and at.session_state['simulation_calls'] == 0
    at.selectbox(key='param_recovery_model').set_value('PCM').run()
    settings = {'param_recovery_reps': 4, 'param_recovery_seed': 17, 'param_recovery_n_person': 40,
                'param_recovery_n_rater': 4, 'param_recovery_n_criterion': 2, 'param_recovery_n_cat': 5}
    for key, value in settings.items():
        at.number_input(key=key).set_value(value)
    at.run()
    assert not at.exception
    before = dict(at.session_state['_parameter_recovery_settings'])
    next(b for b in at.button if b.key.startswith('param_recovery_run::')).click().run()
    assert not at.exception and at.session_state['simulation_calls'] == 1
    at.selectbox(key='report_additional_tool').set_value('none').run()
    assert not at.exception and not at.number_input
    at.radio(key='lang').set_value('ja').run()
    at.selectbox(key='report_additional_tool').set_value('recovery').run()
    assert not at.exception
    assert dict(at.session_state['_parameter_recovery_settings']) == before
    assert any('retained-test-result' in i.value for i in at.info)
    at.checkbox(key='show_tools').uncheck().run()
    at.checkbox(key='show_tools').check().run()
    assert not at.exception and at.session_state['simulation_calls'] == 1
    assert dict(at.session_state['_parameter_recovery_settings']) == before
    assert at.selectbox(key='param_recovery_model').value == 'PCM'
    for key, value in settings.items():
        assert at.number_input(key=key).value == value
    assert any('retained-test-result' in i.value for i in at.info)


@pytest.mark.parametrize('lang', ['en', 'ja'])
def test_guided_optional_entry_keeps_global_hold_and_skips_unrelated_review(lang):
    def view(lang):
        from unittest.mock import patch
        import pandas as pd
        import streamlit as st
        import streamlit_app as app
        st.session_state['lang'] = lang
        st.session_state['guided_export_task'] = 'review'
        st.session_state['guided_export_resource'] = 'additional_analyses'
        st.session_state['opened'] = False
        def opened(result):
            st.session_state['opened'] = True
        with patch.object(app, 'build_publication_gate_summary', return_value=pd.DataFrame([
                {'GateArea': 'Overall manuscript gate', 'GateStatus': 'Not ready'}])), \
             patch.object(app, '_render_additional_report_analyses', opened), \
             patch.object(app, '_render_guided_report_checks', side_effect=AssertionError('Unrelated review')):
            app._render_guided_report_export_section({'config': {'method': 'JMLE'}}, {}, {},
                pd.DataFrame(), pd.DataFrame(), None, None, generate_figures=False)
        st.session_state['expected_hold'] = app.t('guided.export_state_hold')
    at = AppTest.from_function(view, args=(lang,)).run(timeout=45)
    assert not at.exception and at.session_state['opened']
    assert at.session_state['expected_hold'] in [w.value for w in at.warning]


@pytest.mark.parametrize('lang', ['en', 'ja'])
def test_recovery_setup_explains_inputs_and_previews_the_actual_generated_design(lang):
    def view(lang):
        from unittest.mock import patch
        import numpy as np
        import streamlit as st
        import streamlit_app as app
        st.session_state['lang'] = lang
        st.session_state.setdefault('calls', 0)
        def generate_only(**settings):
            st.session_state['calls'] += 1
            st.session_state['submitted'] = settings
            data = app._simulate_recovery_responses(
                theta=np.zeros(settings['n_person']), rater=np.zeros(settings['n_rater']),
                criterion=np.zeros(settings['n_criterion']),
                steps=np.zeros((settings['n_criterion'], settings['n_cat'] - 1)),
                slopes=None, model=settings['model'], rating_min=0,
                rng=np.random.default_rng(settings['seed']))
            st.session_state['generated'] = {'rows': len(data), 'persons': data.Person.nunique(),
                'raters': data.Rater.nunique(), 'criteria': data.Criterion.nunique(),
                'scores': sorted(data.Score.unique().tolist())}
            return {'available': False, 'reason': 'generation-only-test; no fitted simulation'}
        with patch.object(app, 'evaluate_parameter_recovery', generate_only), \
             patch.object(app, 'mfrm_estimate', side_effect=AssertionError('No estimator in this display check')):
            app.render_parameter_recovery_simulation()
    at = AppTest.from_function(view, args=(lang,)).run(timeout=45)
    assert not at.exception and at.session_state['calls'] == 0
    for key in ('param_recovery_model', 'param_recovery_method'):
        assert at.selectbox(key=key).proto.help
    for key in ('param_recovery_reps', 'param_recovery_seed'):
        assert at.number_input(key=key).proto.help
    assert any('270' in c.value for c in at.caption)
    if lang == 'ja':
        labels = {node.label for node in [*at.selectbox, *at.number_input]}
        assert labels == {'生成・推定に使うモデル', '模擬データの推定方法', '受験者数', '評価者数',
                          '評価観点数', '得点カテゴリー数', 'データ生成・推定の反復数', '乱数の種（再現用）'}
        assert '1回の推定内の反復上限とは異なります' in at.number_input(key='param_recovery_reps').proto.help
        assert '十分な精度を保証するものではありません' in at.number_input(key='param_recovery_reps').proto.help
    for key, value in {'param_recovery_n_person': 10, 'param_recovery_n_rater': 2,
                       'param_recovery_n_criterion': 4, 'param_recovery_n_cat': 5,
                       'param_recovery_reps': 3, 'param_recovery_seed': 17}.items():
        at.number_input(key=key).set_value(value)
    at.run()
    assert not at.exception and at.session_state['calls'] == 0
    preview = next(c.value for c in at.caption if '80' in c.value)
    assert ('0〜4' if lang == 'ja' else '0 to 4') in preview
    assert ('3回' if lang == 'ja' else '3 generate-and-fit') in preview
    at.button[0].click().run()
    assert not at.exception and at.session_state['calls'] == 1
    assert dict(at.session_state['generated']) == {'rows': 80, 'persons': 10, 'raters': 2, 'criteria': 4, 'scores': [0, 1, 2, 3, 4]}
    assert at.session_state['submitted']['reps'] == 3
    assert at.session_state['submitted']['seed'] == 17
