"""Beginner disclosure must preserve inputs, evidence and interpretation holds."""
import pytest
from streamlit.testing.v1 import AppTest


@pytest.mark.parametrize("lang", ["en", "ja"])
def test_standard_settings_remain_editable_across_help_and_language(lang):
    def workspace(lang):
        import streamlit as st
        import streamlit_app as app
        st.session_state.setdefault("lang", lang)
        app.main()

    at = AppTest.from_function(workspace, args=(lang,)).run(timeout=45)
    at.button(key="onboarding_skip_guide").click().run(timeout=45)
    settings_label = "分析設定" if lang == "ja" else "Analysis settings"
    settings = next(e for e in at.sidebar.expander if e.label == settings_label)
    assert not settings.proto.expanded
    assert settings.radio(key="facets_mode_model_type").value == "RSM"
    assert settings.radio(key="facets_mode_estimation_method").value == "JMLE"
    at.radio(key="facets_mode_model_type").set_value("PCM").run(timeout=45)
    assert not at.exception
    assert any("PCM · JMLE" in c.value for c in at.main.caption)
    assert at.button(key="facets_mode_run_primary")
    at.button(key="mfrm_help_open_global").click().run(timeout=45)
    assert not at.exception
    other = "ja" if lang == "en" else "en"
    at.radio(key="lang").set_value(other).run(timeout=45)
    assert not at.exception
    assert at.radio(key="facets_mode_model_type").value == "PCM"
    assert "facets_mode_output" not in at.session_state


@pytest.mark.parametrize("lang", ["en", "ja"])
def test_section_navigation_shows_evidence_once_and_keeps_holds(lang):
    def workspace(lang):
        import pandas as pd
        import streamlit as st
        import streamlit_app as app
        from unittest.mock import patch

        st.session_state["lang"] = lang
        plan = pd.DataFrame([{
            "Priority": "P1", "SourceOrder": 1, "Area": "Estimation",
            "Check": "Convergence", "Status": "Do not interpret yet",
            "WhatItSays": "Recorded failure reason", "NextAction": "Inspect the fit.",
            "DetailLocation": "Report",
        }])
        st.session_state["hold_label"] = app.t("guided.interpret_status_pause")
        st.session_state["detail_label"] = app.t("guided.goal_supporting_detail_expander")
        with patch.object(app, "_render_guided_results_section", lambda *a: st.write("Result table")):
            app._render_guided_essential_tabs({}, pd.DataFrame(), {}, {}, {}, pd.DataFrame(),
                pd.DataFrame(), [], "Person", "Score", None, None,
                result_compute_pca=False, result_render_plots=False,
                result_generate_figures=False, action_plan=plan)
        st.session_state["plan"] = plan

    at = AppTest.from_function(workspace, args=(lang,)).run(timeout=45)
    assert not at.exception
    assert len(at.error) == 1 and "Recorded failure reason" in at.error[0].value
    assert at.session_state["hold_label"] in at.error[0].value
    detail = next(e for e in at.expander if e.label == at.session_state["detail_label"])
    assert not detail.proto.expanded
    # Full plan + its legend, with no duplicate plan outside the disclosure.
    assert len(at.dataframe) == len(detail.dataframe) == 2
    assert len(at.get("download_button")) == 1

    original = at.session_state["plan"].copy(deep=True)
    at.button_group(key="guided_essential_section").set_value("results").run(timeout=45)
    assert not at.exception
    assert len(at.error) == 1 and at.session_state["hold_label"] in at.error[0].value
    assert "Recorded failure reason" not in at.error[0].value
    assert any(item.value == "Result table" for item in at.markdown)
    assert not at.get("download_button")
    assert not any(b.key.startswith("guided_action_hub") for b in at.button)
    assert at.session_state["plan"].equals(original)
    at.button_group(key="guided_essential_section").set_value("first_read").run(timeout=45)
    assert not at.exception
    assert "Recorded failure reason" in at.error[0].value
    assert len(at.get("download_button")) == 1


@pytest.mark.parametrize("lang", ["en", "ja"])
def test_person_boundary_remains_visible_when_calculation_details_are_closed(lang):
    def workspace(lang):
        import streamlit as st
        import streamlit_app as app
        from unittest.mock import patch
        st.session_state["lang"] = lang
        st.session_state["expected"] = app.t("guided.person_boundary_hold")
        with patch.object(app, "show_convergence_section"), \
             patch.object(app, "show_person_fit_section"), \
             patch.object(app, "show_posterior_scoring_section"), \
             patch.object(app, "render_guided_measures_reading_guide_tables"):
            app._render_guided_results_section(
                {"person_boundary": {"status": "extreme_persons_present"}}, {}, {})

    at = AppTest.from_function(workspace, args=(lang,)).run(timeout=45)
    assert not at.exception
    assert any(getattr(child, "value", None) == at.session_state["expected"]
               for child in at.main.children.values())


def _boundary_result_and_measures():
    import pandas as pd
    from types import SimpleNamespace
    from tests.test_jmle_person_boundary import _audit
    audit = _audit(estimates=(.25, .4))
    frame = audit['persons'].rename(columns={'Person': 'Level'}).assign(
        Facet='Person', SE=.2, Infit=1., Outfit=1., N_Fit=2,
        CI_Lower=[-.142, .008], CI_Upper=[.642, .792], CI_Level=.95,
        SE_Status='conditional_approximation', CI_Status='conditional_approximation',
        SE_Method='conditional information approximation: other parameters held fixed')
    # Non-person ReportableEstimate is absent by design; it must not blank the facet estimate.
    frame = pd.concat([frame, pd.DataFrame([dict(Facet='Rater', Level='R1', Estimate=.1,
        Infit=1.1, Outfit=1.2, N_Fit=4, SE=.15, CI_Lower=-.194, CI_Upper=.394, CI_Level=.95,
        SE_Status='conditional_approximation', CI_Status='conditional_approximation',
        SE_Method='conditional information approximation: other parameters held fixed')])], ignore_index=True)
    result = {'config': {'method': 'JMLE'}, 'opt': SimpleNamespace(success=True), 'person_boundary': audit}
    return result, frame


@pytest.mark.parametrize('lang', ['en', 'ja'])
def test_boundary_rows_and_real_rating_counts_are_visible_in_basic_and_precision(lang):
    def workspace(lang):
        import streamlit as st
        import streamlit_app as app
        from contextlib import ExitStack
        from unittest.mock import patch
        from tests.test_beginner_ui_simplification import _boundary_result_and_measures
        st.session_state['lang'] = lang
        result, frame = _boundary_result_and_measures()
        with ExitStack() as stack:
            for name in ('show_convergence_section', 'show_person_fit_section', 'show_posterior_scoring_section',
                         'render_eb_shrinkage_section', 'render_guided_measures_reading_guide_tables'):
                stack.enter_context(patch.object(app, name))
            app._render_guided_results_section(result, {'measures': frame}, {})
        st.session_state['source_after'] = frame
    at = AppTest.from_function(workspace, args=(lang,)).run(timeout=45)
    assert not at.exception
    frame = at.dataframe[0].value
    assert frame.Estimate.isna().tolist() == [True, False, False]
    assert frame.EstimateStatus.iloc[0] == ('全て最低点' if lang == 'ja' else 'All minimum scores')
    assert frame.Estimate.iloc[1:].tolist() == [.4, .1]
    assert frame.N_Fit.tolist() == [2, 2, 4]
    assert not {'FiniteJMLEEstimate', 'ReportableEstimate', 'SE_Status'} & set(frame)
    assert at.dataframe[1].value.Estimate.tolist() == [.25, .4, .1]
    at.button_group(key='guided_measures_view').set_value('precision').run(timeout=45)
    assert not at.exception
    precision = at.dataframe[0].value
    assert precision.loc[0, ['Estimate', 'SE', 'CI_Lower', 'CI_Upper']].isna().all()
    assert precision.SE.iloc[1:].tolist() == [.2, .15]
    assert precision.CI_Level.iloc[1:].tolist() == [.95, .95]
    import json
    assert json.loads(at.dataframe[0].proto.columns)['CI_Level']['type_config']['format'] == 'percent'
    assert precision.EstimateStatus.iloc[0] == frame.EstimateStatus.iloc[0]
    assert not {'SE_Status', 'CI_Status', 'EstimateRole'} & set(precision)
    assert at.session_state['source_after'].Estimate.tolist() == [.25, .4, .1]


@pytest.mark.parametrize('hold', ['unfinished', 'unknown', 'fixed_mml', 'free_mml', 'identification'])
def test_current_hold_never_promotes_saved_numbers_to_readable_estimates(hold):
    import pandas as pd
    import streamlit_app as app
    result, frame = _boundary_result_and_measures()
    before = frame.copy(deep=True)
    if hold in ('unfinished', 'unknown'):
        result['opt'].success = False if hold == 'unfinished' else None
    elif hold in ('fixed_mml', 'free_mml'):
        result['config'] = {'method': 'MML', 'estimate_population_sd': hold == 'free_mml'}
    else:
        result['identifiability'] = {'status': 'eta_rank_deficient'}
    for precision in (False, True):
        display = app.guided_measure_display_frame(frame, result, precision=precision)
        assert display.Estimate.isna().all()
        assert display.EstimateStatus.str.len().gt(0).all()
        if precision:
            assert display[['SE', 'CI_Lower', 'CI_Upper']].isna().all().all()
    pd.testing.assert_frame_equal(frame, before)


def test_precision_uses_recorded_basis_and_never_fabricates_missing_intervals():
    import pandas as pd
    import streamlit_app as app
    result, frame = _boundary_result_and_measures()
    frame.loc[1, 'SE_Method'] = 'unknown'
    frame.loc[2, 'CI_Lower'] = float('nan')
    display = app.guided_measure_display_frame(frame, result, precision=True)
    assert pd.isna(display.loc[1, 'SE']) and pd.isna(display.loc[1, 'CI_Lower'])
    assert display.loc[2, 'SE'] == .15
    assert display.loc[2, ['CI_Lower', 'CI_Upper', 'CI_Level']].isna().all()


def test_saved_table_without_person_flags_uses_existing_boundary_audit():
    import streamlit_app as app
    result, frame = _boundary_result_and_measures()
    frame = frame.drop(columns=['FiniteJMLEEstimate', 'PersonInferenceReady', 'ReportableEstimate', 'ExtremeScoreDirection'])
    display = app.guided_measure_display_frame(frame, result)
    assert display.Estimate.isna().tolist() == [True, False, False]


@pytest.mark.parametrize('lang', ['en', 'ja'])
@pytest.mark.parametrize('success', [False, None])
def test_unfinished_next_action_opens_calculation_details_without_document_generation(lang, success):
    def workspace(lang, success):
        import pandas as pd
        import streamlit as st
        import streamlit_app as app
        from types import SimpleNamespace
        from unittest.mock import patch
        st.session_state['lang'] = lang
        result = {'config': {'method': 'JMLE', 'model': 'RSM'},
                  'opt': SimpleNamespace(success=success, message='TEST: iteration limit reached'),
                  'convergence': pd.DataFrame([{'Method': 'JMLE', 'Converged': success,
                                               'Message': 'TEST: iteration limit reached'}])}
        plan = app.build_guided_action_plan(result, {}, {}, localize=True)
        with patch.object(app, '_render_publication_document_section', side_effect=AssertionError('Wrong destination')), \
             patch.object(app, 'show_person_fit_section'), patch.object(app, 'show_posterior_scoring_section'), \
             patch.object(app, 'render_guided_measures_reading_guide_tables'):
            app._render_guided_essential_tabs({}, pd.DataFrame(), result, {}, {}, pd.DataFrame(),
                pd.DataFrame(), [], 'Person', 'Score', None, None, result_compute_pca=False,
                result_render_plots=False, result_generate_figures=False, action_plan=plan)
    at = AppTest.from_function(workspace, args=(lang, success)).run(timeout=45)
    assert not at.exception
    label = '計算の詳細' if lang == 'ja' else 'Calculation details'
    button = at.button(key='guided_action_hub_first_read_primary')
    assert label in button.label
    button.click().run(timeout=45)
    assert not at.exception
    assert at.session_state['guided_essential_section'] == 'results'
    assert at.expander[0].label == label and at.expander[0].proto.expanded
    assert not any(child.type in {'error', 'warning'} for child in at.main.children.values())
    if success is False:
        assert at.expander[0].error
    else:
        assert at.expander[0].warning
    assert 'guided_export_task' not in at.session_state


def test_all_maximum_person_is_identified_without_hiding_interior_or_facet_rows():
    import streamlit_app as app
    from tests.test_jmle_person_boundary import _audit, _boundary_frame
    result, frame = _boundary_result_and_measures()
    scores = _boundary_frame()
    scores.loc[scores.Person.eq('P1'), 'Score'] = 1
    result['person_boundary'] = _audit(estimates=(.25, .4), frame=scores)
    display = app.guided_measure_display_frame(frame, result)
    assert display.Estimate.isna().tolist() == [True, False, False]
    assert display.loc[0, 'EstimateStatus'] == app.t('guided.measure_all_maximum')


def test_weighted_fit_count_keeps_fraction_and_identifies_weight_sum():
    def workspace():
        from contextlib import ExitStack
        from unittest.mock import patch
        import streamlit as st
        import streamlit_app as app
        from tests.test_beginner_ui_simplification import _boundary_result_and_measures
        st.session_state['lang'] = 'ja'
        result, frame = _boundary_result_and_measures()
        result['config']['weight_col'] = 'Weight'
        frame['N_Fit'] = [1.25, 1.5, 2.75]
        with ExitStack() as stack:
            for name in ('show_convergence_section', 'show_person_fit_section', 'show_posterior_scoring_section',
                         'render_eb_shrinkage_section', 'render_guided_measures_reading_guide_tables'):
                stack.enter_context(patch.object(app, name))
            app._render_guided_results_section(result, {'measures': frame}, {})
    import json
    at = AppTest.from_function(workspace).run(timeout=45)
    assert not at.exception
    assert at.dataframe[0].value.N_Fit.tolist() == [1.25, 1.5, 2.75]
    column = json.loads(at.dataframe[0].proto.columns)['N_Fit']
    assert column['label'] == '重み合計'
    assert column['type_config'].get('format') != '%d'
    assert at.dataframe[0].proto.placeholder == '—'


@pytest.mark.parametrize('lang', ['en', 'ja'])
@pytest.mark.parametrize('success', [False, None])
def test_calculation_recovery_leads_with_status_and_keeps_data_without_refitting(lang, success):
    def workspace(lang, success):
        import streamlit as st
        import streamlit_app as app
        import pandas as pd
        from types import SimpleNamespace
        from unittest.mock import patch
        st.session_state['lang'] = lang
        st.session_state.setdefault('guided_essential_section', 'results')
        st.session_state.setdefault('guided_open_calculation_details', True)
        st.session_state.setdefault('facets_mode_maxit', 123)
        st.session_state.setdefault('saved_result', {'config': {'method': 'JMLE'},
            'opt': SimpleNamespace(success=success), 'convergence': pd.DataFrame([{'Method':'JMLE','Converged':success}])})
        with patch.object(app, 'mfrm_estimate', side_effect=AssertionError('Unexpected refit')):
            app._render_guided_calculation_details(st.session_state['saved_result'])
        st.session_state['expected'] = app.t('guided.calculation_unfinished' if success is False else 'guided.calculation_unknown')
        st.session_state['records_label'] = app.t('guided.calculation_records')
    at = AppTest.from_function(workspace, args=(lang,success)).run(timeout=45)
    assert not at.exception
    top = at.main.children[0]
    assert top.type == ('error' if success is False else 'warning')
    assert top.value == at.session_state['expected']
    records=next(e for e in at.expander if e.label==at.session_state['records_label'])
    assert not records.proto.expanded
    assert len(records.dataframe)==len(at.dataframe)
    at.button(key='guided_calculation_review_input').click().run(timeout=45)
    assert not at.exception
    assert at.session_state['guided_essential_section']=='start'
    assert at.session_state['guided_open_calculation_details'] is False
    assert at.session_state['facets_mode_maxit']==123
    assert at.session_state['saved_result']['opt'].success is success


@pytest.mark.parametrize('config', [
    {'method':'MML', 'estimate_population_sd':False},
    {'method':'MML', 'estimate_population_sd':True},
    {'method':'JMLE', 'jmle_facet_boundary_review':True},
])
def test_calculation_summary_does_not_promote_stopping_to_inference(config):
    def workspace(config):
        import streamlit as st
        import streamlit_app as app
        from types import SimpleNamespace
        from unittest.mock import patch
        st.session_state['lang']='ja'
        with patch.object(app, 'show_convergence_section'):
            app._render_guided_calculation_details({'config':config, 'opt':SimpleNamespace(success=True)})
    at=AppTest.from_function(workspace,args=(config,)).run(timeout=45)
    assert not at.exception and not at.success
    assert at.main.children[0].type=='warning'
    assert not at.button


@pytest.mark.parametrize('visible', [False, True])
def test_calculation_summary_deduplicates_only_covered_top_level_alerts(visible):
    def workspace(visible):
        import streamlit as st
        import streamlit_app as app
        import pandas as pd
        from types import SimpleNamespace
        st.session_state['lang']='ja'
        result={'opt':SimpleNamespace(success=False), 'config':{'method':'JMLE',
            'jmle_facet_boundary_review':True, 'anchor_audit':{'overall_status':'error'}}}
        app._show_top_level_warnings(result,{},pd.DataFrame({'Score':[0,1]}),'Person',[],'Score',None,
                                    calculation_summary_visible=visible)
        st.session_state['termination']=app.t('run_warnings.not_converged')
        st.session_state['anchor']=app.t('run_warnings.anchors',count=0)
        st.session_state['identification']=app.t('estimation_subsections.jmle_inference_hold')
    at=AppTest.from_function(workspace,args=(visible,)).run(timeout=45)
    assert not at.exception
    errors='\n'.join(item.value for item in at.error)
    assert (at.session_state['termination'] in errors) is (not visible)
    assert at.session_state['anchor'] in errors
    assert at.session_state['identification'] in '\n'.join(item.value for item in at.warning)
