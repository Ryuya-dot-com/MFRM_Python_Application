"""App-level acceptance tests for the persistent Help adapter."""

from __future__ import annotations

import ast
from dataclasses import replace
import inspect
import json
from pathlib import Path
import re
import textwrap

from streamlit.testing.v1 import AppTest

import streamlit_app as app
from mfrm_app.help_navigation import HelpFallbackReason


REPO_ROOT = Path(__file__).resolve().parents[1]


def _locale_value(language: str, dotted_key: str) -> str:
    node: object = json.loads(
        (REPO_ROOT / "locales" / f"{language}.json").read_text(encoding="utf-8")
    )
    for part in dotted_key.split("."):
        assert isinstance(node, dict)
        node = node[part]
    assert isinstance(node, str)
    return node


def _app_test() -> AppTest:
    return AppTest.from_string(
        "import streamlit_app as app\napp.main()",
        default_timeout=40,
    )


def _instrument_analysis_calls(monkeypatch) -> dict[str, int]:
    calls = {"estimate": 0, "refit": 0}
    original_load_core = app.load_core_namespace
    targets = (
        ("mfrm_estimate", "estimate"),
        ("simulate_refit_design", "refit"),
        ("evaluate_mml_prior_sd_sensitivity", "refit"),
    )

    # Prime the cached namespace while the production callables are intact.
    # Otherwise the first guarded AppTest can cache the temporary sentinels and
    # leak them into a later test that performs a real fit.
    original_load_core.clear()
    baseline_namespace = dict(original_load_core())

    for name, counter in targets:
        original = getattr(app, name, None)
        if not callable(original):
            continue

        def blocked_global(*args, _counter=counter, **kwargs):
            calls[_counter] += 1
            raise AssertionError("Help navigation initiated an analysis call")

        monkeypatch.setattr(app, name, blocked_global)

    def counted_load_core():
        namespace = dict(baseline_namespace)
        for name, counter in targets:
            original = namespace.get(name)
            if not callable(original):
                continue

            def blocked_core(*args, _counter=counter, **kwargs):
                calls[_counter] += 1
                raise AssertionError("Help navigation initiated a core analysis call")

            namespace[name] = blocked_core
        return namespace

    monkeypatch.setattr(app, "load_core_namespace", counted_load_core)
    return calls


def test_facets_mode_value_widgets_have_explicit_state_keys():
    tree = ast.parse(textwrap.dedent(inspect.getsource(app.run_facets_mode)))
    value_widgets = {
        "checkbox",
        "data_editor",
        "file_uploader",
        "multiselect",
        "number_input",
        "radio",
        "select_slider",
        "selectbox",
        "slider",
        "text_area",
        "text_input",
        "toggle",
    }
    missing = sorted(
        (node.lineno, node.func.attr)
        for node in ast.walk(tree)
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Attribute)
        and node.func.attr in value_widgets
        and not any(keyword.arg == "key" for keyword in node.keywords)
    )
    assert not missing, f"Help-state-unsafe FACETS widgets: {missing}"


def test_no_data_help_locale_and_return_preserve_source_without_refit(
    monkeypatch,
):
    calls = _instrument_analysis_calls(monkeypatch)

    at = _app_test().run()
    assert not at.exception
    at.button(key="onboarding_dismiss").click()
    at.run()
    at.selectbox(key="paste_data_delimiter").select_index(3)
    at.run()
    assert at.session_state["data_source_flat"] == "paste"
    delimiter_before = at.session_state["paste_data_delimiter"]
    assert delimiter_before == "Semicolon (;)"
    draft_before = "Person;Score;Rater;Task\n"
    at.text_area(key="paste_data_text").set_value(draft_before)
    at.run()
    assert at.session_state["paste_data_text"] == draft_before
    at.button(key="mfrm_help_open_global").click()
    at.run()
    assert not at.exception
    route_before_locale = at.session_state[app._HELP_ROUTE_STATE_KEY]
    assert route_before_locale.is_open
    assert route_before_locale.help_topic_id == "help.get_started.overview"
    assert route_before_locale.return_destination == (
        "target.input.start",
        "focus.input.start",
    )
    assert calls["estimate"] == 0
    assert calls["refit"] == 0

    at.radio(key="lang").set_value("ja")
    at.run()
    assert not at.exception
    route_after_locale = at.session_state[app._HELP_ROUTE_STATE_KEY]
    assert route_after_locale == route_before_locale
    assert any(
        element.value
        == _locale_value("ja", "help_topics.get_started.overview.title")
        for element in at.subheader
    ), (
        at.session_state["lang"],
        [element.value for element in at.subheader],
    )
    assert calls["estimate"] == 0
    assert calls["refit"] == 0

    at.button(key="mfrm_help_return_to_source").click()
    at.run()
    assert not at.exception
    returned_route = at.session_state[app._HELP_ROUTE_STATE_KEY]
    assert not returned_route.is_open
    assert at.session_state["data_source_flat"] == "paste"
    assert at.session_state["paste_data_delimiter"] == delimiter_before
    assert at.session_state["paste_data_text"] == draft_before
    assert any(
        element.value == _locale_value("ja", "help_nav.returned_status")
        for element in at.success
    )
    assert any(
        element.value == _locale_value("ja", "app.no_input_info")
        for element in at.info
    )
    assert calls["estimate"] == 0
    assert calls["refit"] == 0


def test_help_locale_and_return_preserve_nondefault_analysis_settings(
    monkeypatch,
):
    calls = _instrument_analysis_calls(monkeypatch)
    at = _app_test().run()
    assert not at.exception
    at.button(key="onboarding_skip_guide").click()
    at.run()

    at.radio(key="facets_mode_workflow_mode").set_value("Advanced controls")
    at.radio(key="facets_mode_model_type").set_value("GPCM")
    at.radio(key="facets_mode_estimation_method").set_value("MML")
    at.run()
    at.number_input(key="facets_mode_quad_points").set_value(21)
    at.number_input(key="facets_mode_maxit").set_value(650)
    at.selectbox(key="facets_mode_analysis_depth").set_value("Custom")
    at.run()

    expected = {
        "facets_mode_workflow_mode": "Advanced controls",
        "facets_mode_model_type": "GPCM",
        "facets_mode_estimation_method": "MML",
        "facets_mode_quad_points": 21,
        "facets_mode_maxit": 650,
        "facets_mode_analysis_depth": "Custom",
    }
    assert {key: at.session_state[key] for key in expected} == expected

    at.button(key="mfrm_help_open_global").click()
    at.run()
    assert not at.exception
    assert {key: at.session_state[key] for key in expected} == expected

    at.radio(key="lang").set_value("ja")
    at.run()
    assert not at.exception
    assert {key: at.session_state[key] for key in expected} == expected

    at.button(key="mfrm_help_return_to_source").click()
    at.run()
    assert not at.exception
    assert {key: at.session_state[key] for key in expected} == expected
    assert calls == {"estimate": 0, "refit": 0}


def test_help_does_not_consume_pending_analysis_triggers(monkeypatch):
    calls = _instrument_analysis_calls(monkeypatch)
    at = _app_test().run()
    assert not at.exception

    at.button(key="mfrm_help_open_global").click()
    at.session_state["_facets_mode_force_rerun"] = True
    at.session_state["_onboarding_quickstart_fired"] = True
    at.run()

    assert not at.exception
    assert at.session_state[app._HELP_ROUTE_STATE_KEY].is_open
    assert at.session_state["_facets_mode_force_rerun"] is True
    assert at.session_state["_onboarding_quickstart_fired"] is True
    assert calls == {"estimate": 0, "refit": 0}

    at.radio(key="lang").set_value("ja")
    at.run()
    assert not at.exception
    assert at.session_state["_facets_mode_force_rerun"] is True
    assert at.session_state["_onboarding_quickstart_fired"] is True
    assert calls == {"estimate": 0, "refit": 0}


def test_simulation_threshold_failure_uses_bounded_problem_notice(monkeypatch):
    calls = _instrument_analysis_calls(monkeypatch)
    secret = "PRIVATE_THRESHOLD /Users/researcher/ratings.csv participant=P-31"
    at = _app_test().run()
    assert not at.exception
    at.button(key="onboarding_skip_guide").click()
    at.run()

    at.selectbox(key="data_source_class").set_value("simulate")
    at.run()
    assert at.session_state["data_source_flat"] == "simulate"
    at.radio(key="sim_threshold_mode").set_value("custom")
    at.run()
    at.text_input(key="sim_threshold_text_5").set_value(secret)
    at.run()

    assert not at.exception
    visible = "\n".join(
        str(element.value)
        for collection in (at.error, at.warning, at.markdown, at.caption, at.code)
        for element in collection
    )
    assert "PRIVATE_THRESHOLD" not in visible
    assert "/Users/researcher/ratings.csv" not in visible
    assert "P-31" not in visible
    assert any(
        re.fullmatch(
            r"MFRM-[0-9A-F]{8}-[0-9A-F]{8}-[0-9A-F]{8}",
            str(element.value),
        )
        for element in at.code
    )
    assert calls == {"estimate": 0, "refit": 0}


def test_unknown_help_request_renders_fallback_and_recovers_home():
    at = AppTest.from_string(
        """
import streamlit_app as app

app._ensure_language_state()
if app._HELP_ROUTE_STATE_KEY not in app.st.session_state:
    app._open_registered_help_link("link.does.not.exist")
app.render_persistent_help_surface()
""",
        default_timeout=20,
    ).run()

    assert not at.exception
    fallback_route = at.session_state[app._HELP_ROUTE_STATE_KEY]
    assert fallback_route.is_fallback
    assert fallback_route.fallback_reason is HelpFallbackReason.LINK_UNAVAILABLE
    assert any(
        element.value == _locale_value("en", "help_nav.unavailable_title")
        for element in at.error
    )
    assert any(
        element.value
        == _locale_value("en", "help_topics.fallback.unavailable.title")
        for element in at.subheader
    )

    at.button(key="mfrm_help_return_home").click()
    at.run()
    assert not at.exception
    home_route = at.session_state[app._HELP_ROUTE_STATE_KEY]
    assert home_route.is_open
    assert not home_route.is_fallback
    assert home_route.help_topic_id == "help.get_started.overview"
    assert len(at.selectbox(key=app._HELP_TOPIC_PICKER_KEY).options) == 22
    assert _locale_value(
        "en", "help_topics.fallback.unavailable.title"
    ) not in at.selectbox(key=app._HELP_TOPIC_PICKER_KEY).options


def test_help_surface_enforces_lifecycle_applicability(monkeypatch):
    original_registry = app._help_topics.HELP_REGISTRY
    original_topic = original_registry.topic("help.get_started.overview")
    assert original_topic is not None
    restricted_topic = replace(original_topic, lifecycle_states=("fitted",))
    restricted_registry = replace(
        original_registry,
        topics=tuple(
            restricted_topic
            if topic.help_topic_id == restricted_topic.help_topic_id
            else topic
            for topic in original_registry.topics
        ),
    )
    monkeypatch.setattr(app._help_topics, "HELP_REGISTRY", restricted_registry)

    harness = """
import streamlit_app as app

app._ensure_language_state()
if app._HELP_ROUTE_STATE_KEY not in app.st.session_state:
    app._open_global_help()
app.render_persistent_help_surface()
"""
    no_data = AppTest.from_string(harness, default_timeout=20).run()
    assert not no_data.exception
    no_data_route = no_data.session_state[app._HELP_ROUTE_STATE_KEY]
    assert no_data_route.is_fallback
    assert no_data_route.fallback_reason is HelpFallbackReason.TOPIC_UNAVAILABLE

    fitted = AppTest.from_string(harness, default_timeout=20)
    fitted.session_state["facets_mode_output"] = {}
    fitted.run()
    assert not fitted.exception
    fitted_route = fitted.session_state[app._HELP_ROUTE_STATE_KEY]
    assert fitted_route.is_open
    assert not fitted_route.is_fallback
    assert fitted_route.help_topic_id == "help.get_started.overview"


def test_safe_problem_notice_hides_exception_and_opens_registered_help(
    monkeypatch,
):
    # A legacy deployment variable must never make raw exceptions public.
    monkeypatch.setenv("MFRM_SHOW_TECHNICAL_ERRORS", "1")
    secret = "PRIVATE_SENTINEL /Users/researcher/confidential.csv participant=P-17"
    harness = f"""
import streamlit_app as app

app._ensure_language_state()
if app._get_help_route_state().is_open:
    app.render_persistent_help_surface()
else:
    app.render_user_problem(
        RuntimeError({secret!r}),
        phase=app._user_problems.UserProblemPhase.APPLICATION,
    )
"""
    at = AppTest.from_string(harness, default_timeout=20).run()
    assert not at.exception

    visible = "\n".join(
        str(element.value)
        for collection in (at.error, at.markdown, at.caption, at.code, at.button)
        for element in collection
    )
    assert "PRIVATE_SENTINEL" not in visible
    assert "/Users/researcher/confidential.csv" not in visible
    assert "P-17" not in visible
    support_references = [
        str(element.value)
        for element in at.code
        if re.fullmatch(
            r"MFRM-[0-9A-F]{8}-[0-9A-F]{8}-[0-9A-F]{8}",
            str(element.value),
        )
    ]
    assert len(support_references) == 1

    at.button(key="mfrm_problem_action_problem_unexpected_1").click()
    at.run()
    assert not at.exception
    help_route = at.session_state[app._HELP_ROUTE_STATE_KEY]
    assert help_route.is_open
    assert help_route.help_topic_id == "help.problem.unexpected"
    assert help_route.return_destination == (
        "target.help.problem_resolution",
        "focus.help.problem_resolution",
    )
    help_visible = "\n".join(
        str(element.value)
        for collection in (at.error, at.markdown, at.caption, at.code, at.button)
        for element in collection
    )
    assert "PRIVATE_SENTINEL" not in help_visible
    assert "/Users/researcher/confidential.csv" not in help_visible
