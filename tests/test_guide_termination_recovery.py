"""Guide progress follows recorded termination; failure review never authorizes inference."""
from copy import deepcopy

import numpy as np
import pandas as pd
import pytest
from streamlit.testing.v1 import AppTest

import streamlit_app as app
from mfrm_app import engine, guidance, output_qualification


@pytest.mark.parametrize('method,maxit,free', [
    ('JMLE', 400, False), ('JMLE', 1, False), ('MML', 100, False), ('MML', 100, True),
])
def test_new_fit_summary_is_qualified_with_the_complete_result(monkeypatch, method, maxit, free):
    contexts = []
    qualify = app.qualified_estimation_summary
    def inspect(result):
        contexts.append(result)
        assert 'opt' in result and 'identifiability' in result and 'person_boundary' in result
        return qualify(result)
    monkeypatch.setattr(engine, 'qualified_estimation_summary', inspect)
    result = app.mfrm_estimate(app._make_self_test_gradient_data(), 'Person', ['Rater', 'Task'],
                               'Score', method=method, maxit=maxit, estimate_population_sd=free)
    assert contexts
    stopped = output_qualification.optimizer_termination_established(result)
    assert stopped == (maxit != 1)
    assert bool(result['summary'].Converged.iloc[0]) == stopped
    assert bool(result['convergence'].Converged.iloc[0]) == stopped
    assert 'SavedConverged' not in result['summary']
    expected = (app.FREE_SD_MML_INFERENCE_HOLD_REASON if free else app.FIXED_SD_MML_INFERENCE_HOLD_REASON
                ) if method == 'MML' else ('' if stopped else output_qualification.OPTIMIZER_TERMINATION_HOLD_REASON)
    assert output_qualification.calibration_inference_hold_reason(result) == expected
    if expected:
        assert not result['summary'].InferenceReady.iloc[0]
        assert result['summary'].InferenceReadinessReason.iloc[0] == expected
    else:
        assert result['summary'].InferenceReady.iloc[0]
        assert 'InferenceReadinessReason' not in result['summary']


@pytest.mark.parametrize('success,allowed', [
    (True, True), (np.bool_(True), True), (False, False), (None, False),
    ('True', False), ('False', False), (1, False),
])
def test_guide_uses_explicit_optimizer_evidence_instead_of_truthy_or_stale_summary(monkeypatch, success, allowed):
    from types import SimpleNamespace
    state = {}
    monkeypatch.setattr(app.st, 'session_state', state)
    app._dispatch_guidance_event(guidance.GuidanceEventType.START_SAMPLE)
    app._dispatch_guidance_event(guidance.GuidanceEventType.NODE_REVIEWED,
                                 data_fingerprint='sample-data', learning_evidence_id='learning.data_roles_reviewed')
    result = {'opt': SimpleNamespace(success=success), 'config': {'method': 'JMLE'},
              'summary': pd.DataFrame([{'Converged': not allowed, 'InferenceReady': False}])}
    original = result['summary'].copy(deep=True)
    assert app._advance_sample_guide_after_fit(result) is allowed
    after = app._get_guidance_state()
    assert after.active_node_id is (guidance.GuideNode.EVIDENCE_REVIEW if allowed else guidance.GuideNode.ESTIMATE)
    assert not after.learning_complete and not after.workflow_complete
    pd.testing.assert_frame_equal(result['summary'], original)


def _run(at):
    at.run(timeout=60)
    assert not at.exception, [e.value for e in at.exception]


def _start():
    at = AppTest.from_string('import streamlit_app as app\napp.main()')
    _run(at)
    at.button(key='onboarding_quickstart').click(); _run(at)
    at.button(key='sample_guide_data_check_continue').click(); _run(at)
    return at


def _inject_first_attempt(monkeypatch, kind):
    app.load_core_namespace.clear()
    baseline = dict(app.load_core_namespace())
    estimate = baseline['mfrm_estimate']
    calls = []
    def attempt(*a, **kw):
        calls.append(dict(kw))
        if len(calls) == 1:
            if kind == 'exception':
                raise RuntimeError('PRIVATE test detail must not be retained by the guide')
            kw = dict(kw, maxit=1)
        return estimate(*a, **kw)
    monkeypatch.setattr(app, 'load_core_namespace', lambda: dict(baseline, mfrm_estimate=attempt))
    return calls


@pytest.mark.parametrize('lang', ['en', 'ja'])
def test_failed_fit_can_be_reviewed_and_explicitly_retried_without_automatic_refits(monkeypatch, lang):
    calls = _inject_first_attempt(monkeypatch, 'termination')
    at = _start()
    at.radio(key='lang').set_value(lang); _run(at)
    at.button(key='sample_guide_estimate_run').click(); _run(at)
    assert len(calls) == 1
    failed = deepcopy(at.session_state['facets_mode_output']['result'])
    assert not failed['opt'].success and not failed['summary'].Converged.iloc[0]
    assert not failed['summary'].InferenceReady.iloc[0]
    assert at.session_state[app._GUIDANCE_STATE_KEY].active_node_id is guidance.GuideNode.ESTIMATE
    assert at.button(key='sample_guide_review_attempt')
    at.radio(key='lang').set_value('ja' if lang == 'en' else 'en'); _run(at)
    assert len(calls) == 1
    at.button(key='sample_guide_review_attempt').click(); _run(at)
    assert at.radio(key='facets_mode_workflow_mode')
    assert at.session_state[app._GUIDANCE_STATE_KEY].active_node_id is guidance.GuideNode.ESTIMATE
    at.button(key='mfrm_help_open_global').click(); _run(at)
    at.button(key='mfrm_help_return_to_source').click(); _run(at)
    assert len(calls) == 1
    pd.testing.assert_frame_equal(at.session_state['facets_mode_output']['result']['summary'], failed['summary'])
    at.radio(key='facets_mode_workflow_mode').set_value('Advanced controls'); _run(at)
    at.number_input(key='facets_mode_maxit').set_value(400); _run(at)
    assert len(calls) == 1
    # The existing settings view exposes its sidebar Run control.
    label = app._load_locale(at.session_state['lang'])['sidebar_perf']['run_button']
    next(b for b in at.sidebar.button if b.label == label).click(); _run(at)
    assert len(calls) == 2
    assert at.session_state['facets_mode_output']['result']['opt'].success
    assert at.session_state[app._GUIDANCE_STATE_KEY].active_node_id is guidance.GuideNode.EVIDENCE_REVIEW
    assert '_facets_mode_guide_attempt_status' not in at.session_state
    assert '_facets_mode_guide_review_attempt' not in at.session_state


def test_exception_keeps_bounded_problem_description_across_locale_and_allows_retry(monkeypatch):
    calls = _inject_first_attempt(monkeypatch, 'exception')
    at = _start()
    at.button(key='sample_guide_estimate_run').click(); _run(at)
    assert len(calls) == 1 and 'facets_mode_output' not in at.session_state
    assert at.session_state['_facets_mode_guide_attempt_status'] == 'exception'
    notice = at.session_state['_facets_mode_guide_attempt_problem']
    assert notice is not None and 'PRIVATE' not in repr(notice)
    at.radio(key='lang').set_value('ja'); _run(at)
    assert len(calls) == 1 and not any(b.key == 'sample_guide_review_attempt' for b in at.button)
    assert any('計算を完了できませんでした' in w.value for w in at.warning)
    assert any(c.value == notice.support_reference for c in at.code)
    at.button(key='sample_guide_estimate_run').click(); _run(at)
    assert len(calls) == 2
    assert at.session_state[app._GUIDANCE_STATE_KEY].active_node_id is guidance.GuideNode.EVIDENCE_REVIEW
    assert '_facets_mode_guide_attempt_problem' not in at.session_state
