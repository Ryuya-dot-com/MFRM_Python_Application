"""R3 acceptance tests for exact contextual Help and return intent."""

from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path

from streamlit.testing.v1 import AppTest
import pandas as pd

import streamlit_app as app
from mfrm_app.help_navigation import HelpContextStatus, HelpFallbackReason


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


def _block_analysis_calls(monkeypatch) -> dict[str, int]:
    calls = {"estimate": 0, "refit": 0}
    targets = (
        ("mfrm_estimate", "estimate"),
        ("simulate_refit_design", "refit"),
        ("evaluate_mml_prior_sd_sensitivity", "refit"),
    )
    original_load_core = app.load_core_namespace

    # Prime the cached namespace before replacing module-level callables. This
    # keeps a Help-navigation sentinel from surviving into another real-fit
    # AppTest through ``st.cache_resource``.
    original_load_core.clear()
    baseline_namespace = dict(original_load_core())

    for name, counter in targets:
        if not callable(getattr(app, name, None)):
            continue

        def blocked_global(*args, _counter=counter, **kwargs):
            calls[_counter] += 1
            raise AssertionError("Help return initiated an analysis call")

        monkeypatch.setattr(app, name, blocked_global)

    def guarded_load_core():
        namespace = dict(baseline_namespace)
        for name, counter in targets:
            if not callable(namespace.get(name)):
                continue

            def blocked_core(*args, _counter=counter, **kwargs):
                calls[_counter] += 1
                raise AssertionError("Help return initiated a core analysis call")

            namespace[name] = blocked_core
        return namespace

    monkeypatch.setattr(app, "load_core_namespace", guarded_load_core)
    return calls


CONTEXTUAL_HARNESS = """
import streamlit_app as app
import pandas as pd

app._ensure_language_state()

state = app._get_help_route_state()
if state.is_open:
    app.render_persistent_help_surface(state)
else:
    app._draw_fit_scatter(
        pd.DataFrame(
            {
                "Facet": ["Item"],
                "Level": ["I1"],
                "Infit": [1.0],
                "Outfit": [1.0],
            }
        )
    )
"""


def _return_focus_controllers(at: AppTest) -> list[object]:
    return [
        element.proto
        for element in at.get("html")
        if "mfrm-focus-results-figure-fit-scatter" in element.proto.body
    ]


def _assert_acceptance_states_not_rendered(at: AppTest) -> None:
    rendered = "\n".join(
        str(element.proto)
        for element_type in (
            "title",
            "header",
            "subheader",
            "markdown",
            "caption",
            "text",
            "button",
            "radio",
            "selectbox",
            "success",
            "info",
            "warning",
            "error",
            "html",
        )
        for element in at.get(element_type)
    )
    for internal_state in (
        "APPTEST_ONLY",
        "BROWSER_ACCEPTED",
        "apptest_only",
        "browser_accepted",
    ):
        assert internal_state not in rendered


def test_fit_scatter_opens_exact_section_and_returns_to_full_panel(monkeypatch):
    context = app._help_contract.HelpContextBinding(
        data_context="sample",
        analysis_id="analysis-fit-42",
        analysis_phase="fitted",
    )
    monkeypatch.setattr(app, "_current_fitted_help_context", lambda: context)
    at = AppTest.from_string(CONTEXTUAL_HARNESS, default_timeout=20)
    at.session_state["lang"] = "en"
    at.session_state["app_view_density"] = "Full"
    at.session_state["main_results_panel"] = "fit_details"
    at.session_state["fit_df_method_method"] = "both"
    at.session_state["fit_df_method_cap"] = 12.5
    at.session_state["facets_mode_model_type"] = "GPCM"
    at.session_state["facets_mode_maxit"] = 650
    at.session_state["_facets_mode_force_rerun"] = True
    at.run()
    assert not at.exception
    _assert_acceptance_states_not_rendered(at)
    initial_headers = [
        element
        for element in at.subheader
        if element.value == _locale_value("en", "fit_details.scatter_subheader")
    ]
    assert len(initial_headers) == 1
    assert initial_headers[0].proto.anchor == "mfrm-focus-results-figure-fit-scatter"

    at.button(key="mfrm_popover_open_full_guide_fit_scatter").click()
    at.run()
    assert not at.exception
    _assert_acceptance_states_not_rendered(at)
    opened = at.session_state[app._HELP_ROUTE_STATE_KEY]
    assert opened.is_open
    assert opened.help_link_id == "link.popover.fit_scatter"
    assert opened.help_topic_id == "help.results.fit"
    assert opened.section_id == "fit-scatter"
    assert opened.return_destination == (
        "target.results.figure.fit_scatter",
        "focus.results.figure.fit_scatter",
    )
    assert opened.context_status is HelpContextStatus.CURRENT
    intent = at.session_state[app._HELP_RETURN_INTENT_KEY]
    assert intent.phase == "source_open"
    assert intent.projection_id == "full.fit_details"
    assert intent.source_view_state == (
        ("fit_df_method_method", "both"),
        ("fit_df_method_cap", 12.5),
    )
    assert any(
        _locale_value(
            "en",
            "help_topics.results.fit.sections.fit_scatter.title",
        )
        in str(element.value)
        for element in at.markdown
    )

    route_before_locale = opened
    at.session_state["lang"] = "ja"
    at.run()
    assert not at.exception
    _assert_acceptance_states_not_rendered(at)
    localized = at.session_state[app._HELP_ROUTE_STATE_KEY]
    assert localized == route_before_locale
    assert at.session_state[app._HELP_RETURN_INTENT_KEY] == intent
    assert any(
        _locale_value(
            "ja",
            "help_topics.results.fit.sections.fit_scatter.title",
        )
        in str(element.value)
        for element in at.markdown
    )

    at.button(key="mfrm_help_return_to_source").click()
    at.run()
    assert not at.exception
    _assert_acceptance_states_not_rendered(at)
    assert not at.session_state[app._HELP_ROUTE_STATE_KEY].is_open
    assert app._HELP_RETURN_INTENT_KEY not in at.session_state
    assert at.session_state["app_view_density"] == "Full"
    assert at.session_state["main_results_panel"] == "fit_details"
    assert at.session_state["fit_df_method_method"] == "both"
    assert at.session_state["fit_df_method_cap"] == 12.5
    assert at.session_state["facets_mode_model_type"] == "GPCM"
    assert at.session_state["facets_mode_maxit"] == 650
    assert at.session_state["_facets_mode_force_rerun"] is True
    assert any(
        element.value == _locale_value("ja", "help_nav.returned_status")
        for element in at.caption
    )
    fit_scatter_headers = [
        element
        for element in at.subheader
        if element.value == _locale_value("ja", "fit_details.scatter_subheader")
    ]
    assert len(fit_scatter_headers) == 1
    assert (
        fit_scatter_headers[0].proto.anchor
        == "mfrm-focus-results-figure-fit-scatter"
    )
    controllers = _return_focus_controllers(at)
    assert len(controllers) == 1
    assert controllers[0].body == app._HELP_RETURN_FOCUS_HTML[
        "focus.results.figure.fit_scatter"
    ]
    assert controllers[0].unsafe_allow_javascript
    assert not any(
        element.value == _locale_value("ja", "help_nav.returned_status")
        for element in at.success
    )

    at.run()
    assert not at.exception
    _assert_acceptance_states_not_rendered(at)
    assert app._HELP_RETURN_INTENT_KEY not in at.session_state
    assert not _return_focus_controllers(at)
    assert not any(
        element.value == _locale_value("ja", "help_nav.returned_status")
        for element in at.caption
    )


def test_return_focus_html_is_static_registered_and_accessibility_bounded():
    focus_id = "focus.results.figure.fit_scatter"
    body = app._help_return_focus_html(focus_id)

    assert set(app._HELP_RETURN_FOCUS_HTML) == set(app._HELP_FOCUS_MARKER_IDS)
    assert body == app._HELP_RETURN_FOCUS_HTML[focus_id]
    assert body.count("<script>") == 1
    assert body.count("<style>") == 1
    assert 'document.getElementById("mfrm-focus-results-figure-fit-scatter")' in body
    assert 'target.tagName !== "H3"' in body
    assert 'target.dataset.mfrmFocusHandled === "true"' in body
    assert 'target.dataset.mfrmFocusHandled = "true"' in body
    assert "delete target.dataset.mfrmFocusHandled" not in body
    assert "activeAtMount instanceof HTMLButtonElement" in body
    assert "!activeAtMount.isConnected" in body
    assert "activeNow === document.documentElement" in body
    assert 'target.setAttribute("tabindex", "-1")' in body
    assert 'target.dataset.mfrmFocusRequest = "return-from-help"' in body
    assert 'target.dataset.mfrmFocusStatus = "attempting"' in body
    assert "document.activeElement === target" in body
    assert "target.focus({preventScroll: true})" in body
    assert 'target.dataset.mfrmFocusStatus = "focus-error"' in body
    assert 'scrollIntoView({behavior: "auto"' in body
    assert 'target.dataset.mfrmFocusStatus = "scroll-error"' in body
    assert "scheduleFrame(() => scheduleFrame(restoreFocus))" in body
    assert 'outline: 3px solid currentColor' in body
    assert "aria-label" not in body
    assert "aria-live" not in body
    assert "role=" not in body
    assert "analysis-fit-42" not in body
    assert _locale_value("en", "help_nav.returned_status") not in body
    assert _locale_value("ja", "help_nav.returned_status") not in body
    assert "innerHTML" not in body
    assert "eval(" not in body
    assert "fetch(" not in body
    assert "window.parent" not in body
    assert "window.top" not in body

    try:
        app._help_return_focus_html("focus.not.registered")
    except ValueError as exc:
        assert str(exc) == "focus_id is not registered for Help return"
    else:  # pragma: no cover - contract failure guard
        raise AssertionError("unregistered focus ID was accepted")


def test_fit_scatter_essential_projection_mutates_presentation_keys_only(
    monkeypatch,
):
    context = app._help_contract.HelpContextBinding(
        data_context="sample",
        analysis_id="analysis-fit-42",
        analysis_phase="fitted",
    )
    monkeypatch.setattr(
        app,
        "_capture_help_return_source_view_state",
        lambda _target_id: (
            ("fit_df_method_method", "facets"),
            ("fit_df_method_cap", 9.0),
        ),
    )
    intent = app._new_help_return_intent(
        help_link_id="link.popover.fit_scatter",
        projection_id="essential.fit_details",
        context=context,
    )
    assert intent is not None
    projection = dict(
        app._HELP_RETURN_PROJECTIONS[(intent.target_id, intent.projection_id)]
    )
    assert projection == {
        "app_view_density": "Essential",
        "guided_essential_section": "diagnostics",
        "guided_diagnostics_panel": "fit_details",
    }
    assert app._HELP_RETURN_PRESENTATION_VALUE_ALLOWLIST == {
        "app_view_density": frozenset({"Essential", "Full"}),
        "main_results_panel": frozenset({"fit_details"}),
        "guided_essential_section": frozenset({"diagnostics"}),
        "guided_diagnostics_panel": frozenset({"fit_details"}),
    }
    assert app._HELP_RETURN_PRESENTATION_KEYS == frozenset(
        {
            "app_view_density",
            "main_results_panel",
            "guided_essential_section",
            "guided_diagnostics_panel",
        }
    )
    assert set(projection).issubset(app._HELP_RETURN_PRESENTATION_KEYS)
    assert not set(projection).intersection(
        {
            "data_source_flat",
            "facets_mode_model_type",
            "facets_mode_estimation_method",
            "facets_mode_analysis_depth",
            "_facets_mode_force_rerun",
            "_onboarding_quickstart_fired",
        }
    )
    assert (
        app._validated_help_return_intent(
            replace(
                intent,
                source_view_state=(
                    ("fit_df_method_method", "both"),
                    ("fit_df_method_cap", 50.5),
                ),
            )
        )
        is None
    )
    assert (
        app._validated_help_return_intent(
            replace(
                intent,
                source_view_state=(
                    ("fit_df_method_method", "both"),
                    ("facets_mode_model_type", "GPCM"),
                ),
            )
        )
        is None
    )


def test_invalid_contextual_return_intent_fails_closed_without_success(monkeypatch):
    context = app._help_contract.HelpContextBinding(
        data_context="sample",
        analysis_id="analysis-fit-42",
        analysis_phase="fitted",
    )
    monkeypatch.setattr(app, "_current_fitted_help_context", lambda: context)
    harness = """
import streamlit_app as app

app._ensure_language_state()
if not app.st.session_state.get("opened_once"):
    app.st.session_state["opened_once"] = True
    app._open_registered_help_link(
        "link.popover.fit_scatter",
        "projection.does_not.exist",
    )
app.render_persistent_help_surface()
"""
    at = AppTest.from_string(harness, default_timeout=20).run()
    assert not at.exception
    assert at.session_state[app._HELP_ROUTE_STATE_KEY].is_open
    assert app._HELP_RETURN_INTENT_KEY not in at.session_state

    at.button(key="mfrm_help_return_to_source").click()
    at.run()
    assert not at.exception
    route = at.session_state[app._HELP_ROUTE_STATE_KEY]
    assert route.is_open
    assert route.is_fallback
    assert route.fallback_reason is HelpFallbackReason.RETURN_TARGET_UNAVAILABLE
    assert route.return_destination is None
    assert not any(
        element.value == _locale_value("en", "help_nav.returned_status")
        for element in at.caption
    )
    assert not _return_focus_controllers(at)


def test_only_explicitly_activated_popover_has_full_guide_edge():
    assert {
        key
        for key in app._help_topics.HELP_POPOVER_LINK_IDS
        if app._registered_popover_help_link_id(key) is not None
    } == {"fit_scatter"}
    assert (
        app._registered_popover_help_link_id("fit_scatter")
        == "link.popover.fit_scatter"
    )
    binding = app._ACTIVE_CONTEXTUAL_HELP_ADAPTER_BINDINGS["fit_scatter"]
    assert binding == app._ContextualHelpAdapterBinding(
        topic_key="fit_scatter",
        help_link_id="link.popover.fit_scatter",
        help_topic_id="help.results.fit",
        section_id="fit-scatter",
        target_id="target.results.figure.fit_scatter",
        focus_id="focus.results.figure.fit_scatter",
        surface_id="surface.results",
        panel_id="panel.results.figure.fit_scatter",
        context_policy=app._help_contract.HelpContextPolicy.REQUIRED_CURRENT,
        acceptance=app._ContextualHelpAcceptance.APPTEST_ONLY,
        browser_evidence_id=None,
    )
    assert app._CONTEXTUAL_HELP_APPTEST_ONLY_PILOT_KEYS == frozenset(
        {"fit_scatter"}
    )
    assert app._ACTIVE_CONTEXTUAL_HELP_POPOVER_KEYS == frozenset({"fit_scatter"})
    assert app._ACTIVE_CONTEXTUAL_HELP_RETURN_DESTINATIONS == frozenset(
        {
            (
                "target.results.figure.fit_scatter",
                "focus.results.figure.fit_scatter",
            )
        }
    )
    assert app._ACTIVE_CONTEXTUAL_HELP_POPOVER_KEYS == (
        app._active_contextual_help_popover_keys()
    )
    assert app._contextual_help_adapter_passes_acceptance_gate(
        "fit_scatter",
        binding,
    )
    link = app._help_topics.HELP_REGISTRY.link(binding.help_link_id)
    target = app._help_topics.HELP_REGISTRY.target(binding.target_id)
    assert link is not None
    assert target is not None
    assert link.source_target_id == binding.target_id
    assert target.surface_id == binding.surface_id
    assert target.presentation_state["panel_id"] == binding.panel_id
    assert app._registered_popover_help_link_id("posterior_trace") is None
    assert app._registered_popover_help_link_id("posterior_rhat_ess") is None
    assert app._registered_popover_help_link_id("not_registered") is None

    posterior = AppTest.from_string(
        "import streamlit_app as app\napp.render_help_popover('posterior_trace')",
        default_timeout=20,
    ).run()
    assert not posterior.exception
    assert not any(
        button.key == "mfrm_popover_open_full_guide_posterior_trace"
        for button in posterior.button
    )


def test_contextual_help_acceptance_gate_requires_browser_evidence_after_pilot(
    monkeypatch,
):
    binding = app._ACTIVE_CONTEXTUAL_HELP_ADAPTER_BINDINGS["fit_scatter"]
    non_pilot = replace(binding, topic_key="scree")

    assert not app._contextual_help_adapter_passes_acceptance_gate(
        "scree",
        binding,
    )
    assert not app._contextual_help_adapter_passes_acceptance_gate(
        "fit_scatter",
        object(),
    )
    assert not app._contextual_help_adapter_passes_acceptance_gate(
        "scree",
        non_pilot,
    )
    assert not app._contextual_help_adapter_passes_acceptance_gate(
        "fit_scatter",
        replace(binding, browser_evidence_id="browser.help.fit_scatter.unverified"),
    )
    with monkeypatch.context() as pilot_patch:
        pilot_patch.setattr(
            app,
            "_CONTEXTUAL_HELP_APPTEST_ONLY_PILOT_KEYS",
            frozenset({"fit_scatter", "scree"}),
        )
        assert not app._contextual_help_adapter_passes_acceptance_gate(
            "fit_scatter",
            binding,
        )

    for invalid_evidence in (
        None,
        "",
        " browser.help.scree.chrome_macos.20260724",
        "browser evidence from a manual check",
        "browser.a",
        "browser.help.scree.",
        "browser.help.fit_scatter.r3.20260724",
    ):
        candidate = replace(
            non_pilot,
            acceptance=app._ContextualHelpAcceptance.BROWSER_ACCEPTED,
            browser_evidence_id=invalid_evidence,
        )
        assert not app._contextual_help_adapter_passes_acceptance_gate(
            "scree",
            candidate,
        )

    accepted = replace(
        non_pilot,
        acceptance=app._ContextualHelpAcceptance.BROWSER_ACCEPTED,
        browser_evidence_id="browser.help.scree.r3.20260724",
    )
    assert app._contextual_help_adapter_passes_acceptance_gate(
        "scree",
        accepted,
    )
    assert not app._contextual_help_adapter_passes_acceptance_gate(
        "scree",
        replace(accepted, acceptance="browser_accepted"),
    )

    with monkeypatch.context() as binding_patch:
        binding_patch.setitem(
            app._ACTIVE_CONTEXTUAL_HELP_ADAPTER_BINDINGS,
            "scree",
            accepted,
        )
        assert app._active_contextual_help_popover_keys() == frozenset(
            {"fit_scatter"}
        )

    with monkeypatch.context() as projection_patch:
        projection_patch.setitem(
            app._CONTEXTUAL_HELP_RETURN_PROJECTION_IDS,
            "fit_scatter",
            frozenset({"projection.not.registered"}),
        )
        assert app._registered_popover_help_link_id("fit_scatter") is None
        assert not app._active_contextual_help_popover_keys()

    monkeypatch.setitem(
        app._ACTIVE_CONTEXTUAL_HELP_ADAPTER_BINDINGS,
        "fit_scatter",
        replace(
            binding,
            acceptance=app._ContextualHelpAcceptance.BROWSER_ACCEPTED,
            browser_evidence_id=None,
        ),
    )
    assert app._registered_popover_help_link_id("fit_scatter") is None


def test_browser_runbook_records_pilot_as_not_run_and_locales_keep_states_internal():
    runbook = (REPO_ROOT / "docs" / "help_browser_acceptance.md").read_text(
        encoding="utf-8"
    )
    binding = app._ACTIVE_CONTEXTUAL_HELP_ADAPTER_BINDINGS["fit_scatter"]

    assert binding.acceptance is app._ContextualHelpAcceptance.APPTEST_ONLY
    assert binding.browser_evidence_id is None
    assert "- Current pilot: `fit_scatter`" in runbook
    assert "- Current stage: `APPTEST_ONLY`" in runbook
    assert "- Browser evidence ID: none" in runbook
    assert "- Browser result: **NOT RUN**" in runbook
    for case_id in (
        "CORE-01",
        "CORE-02",
        "CORE-03",
        "CORE-04",
        "REPEAT-01",
        "ABORT-01",
        "STALE-01",
        "CSP-01",
        "MOTION-01",
        "REFLOW-01/02",
    ):
        assert case_id in runbook

    for language in ("en", "ja"):
        locale_text = (REPO_ROOT / "locales" / f"{language}.json").read_text(
            encoding="utf-8"
        )
        assert "APPTEST_ONLY" not in locale_text
        assert "BROWSER_ACCEPTED" not in locale_text
        assert "apptest_only" not in locale_text
        assert "browser_accepted" not in locale_text


def test_legacy_output_never_infers_provenance_from_current_sidebar_source():
    legacy_output = {"result": {"config": {}}}

    assert app._help_data_context_for_output(legacy_output) is None
    assert (
        app._help_data_context_for_output(
            {**legacy_output, "_help_data_context": "not-a-data-context"}
        )
        is None
    )


def test_topic_change_with_missing_intent_fails_closed(monkeypatch):
    context = app._help_contract.HelpContextBinding(
        data_context="sample",
        analysis_id="analysis-fit-42",
        analysis_phase="fitted",
    )
    monkeypatch.setattr(app, "_current_fitted_help_context", lambda: context)
    at = AppTest.from_string(CONTEXTUAL_HARNESS, default_timeout=20)
    at.session_state["lang"] = "en"
    at.session_state["app_view_density"] = "Full"
    at.session_state["main_results_panel"] = "data"
    at.session_state["fit_df_method_method"] = "facets"
    at.session_state["fit_df_method_cap"] = 9.0
    at.run()
    assert not at.exception

    at.button(key="mfrm_popover_open_full_guide_fit_scatter").click()
    at.run()
    assert not at.exception
    assert app._HELP_RETURN_INTENT_KEY in at.session_state

    at.selectbox(key=app._HELP_TOPIC_PICKER_KEY).set_value(
        "help.get_started.overview"
    )
    at.run()
    assert not at.exception
    changed = at.session_state[app._HELP_ROUTE_STATE_KEY]
    assert changed.help_link_id is None
    assert changed.return_destination == (
        "target.results.figure.fit_scatter",
        "focus.results.figure.fit_scatter",
    )

    del at.session_state[app._HELP_RETURN_INTENT_KEY]
    at.button(key="mfrm_help_return_to_source").click()
    at.run()
    assert not at.exception
    route = at.session_state[app._HELP_ROUTE_STATE_KEY]
    assert route.is_open
    assert route.is_fallback
    assert route.fallback_reason is HelpFallbackReason.RETURN_TARGET_UNAVAILABLE
    assert route.return_destination is None
    assert at.session_state["main_results_panel"] == "data"
    assert not any(
        element.value == _locale_value("en", "help_nav.returned_status")
        for element in at.caption
    )
    assert not _return_focus_controllers(at)


def test_topic_change_with_intact_intent_returns_to_exact_source(monkeypatch):
    context = app._help_contract.HelpContextBinding(
        data_context="sample",
        analysis_id="analysis-fit-42",
        analysis_phase="fitted",
    )
    monkeypatch.setattr(app, "_current_fitted_help_context", lambda: context)
    at = AppTest.from_string(CONTEXTUAL_HARNESS, default_timeout=20)
    at.session_state["lang"] = "en"
    at.session_state["app_view_density"] = "Full"
    at.session_state["main_results_panel"] = "data"
    at.session_state["fit_df_method_method"] = "both"
    at.session_state["fit_df_method_cap"] = 12.5
    at.run()
    assert not at.exception

    at.button(key="mfrm_popover_open_full_guide_fit_scatter").click()
    at.run()
    assert not at.exception
    intent = at.session_state[app._HELP_RETURN_INTENT_KEY]
    at.selectbox(key=app._HELP_TOPIC_PICKER_KEY).set_value(
        "help.get_started.overview"
    )
    at.run()
    assert not at.exception
    assert at.session_state[app._HELP_RETURN_INTENT_KEY] == intent
    assert (
        at.session_state[app._HELP_ROUTE_STATE_KEY].context_status
        is HelpContextStatus.STATIC_ONLY
    )

    at.button(key="mfrm_help_return_to_source").click()
    at.run()
    assert not at.exception
    assert not at.session_state[app._HELP_ROUTE_STATE_KEY].is_open
    assert app._HELP_RETURN_INTENT_KEY not in at.session_state
    assert at.session_state["main_results_panel"] == "fit_details"
    assert at.session_state["fit_df_method_method"] == "both"
    assert at.session_state["fit_df_method_cap"] == 12.5
    assert any(
        element.value == _locale_value("en", "help_nav.returned_status")
        for element in at.caption
    )


def test_context_change_in_help_fails_before_restoring_source_state(monkeypatch):
    def current_context():
        return app._help_contract.HelpContextBinding(
            data_context=app.st.session_state.get("test_data_context", "sample"),
            analysis_id=app.st.session_state.get("test_analysis_id", "analysis-A"),
            analysis_phase="fitted",
        )

    monkeypatch.setattr(app, "_current_fitted_help_context", current_context)
    at = AppTest.from_string(CONTEXTUAL_HARNESS, default_timeout=20)
    at.session_state["lang"] = "en"
    at.session_state["test_analysis_id"] = "analysis-A"
    at.session_state["test_data_context"] = "sample"
    at.session_state["app_view_density"] = "Full"
    at.session_state["main_results_panel"] = "fit_details"
    at.session_state["fit_df_method_method"] = "both"
    at.session_state["fit_df_method_cap"] = 12.5
    at.run()
    assert not at.exception

    at.button(key="mfrm_popover_open_full_guide_fit_scatter").click()
    at.run()
    assert not at.exception
    assert at.session_state[app._HELP_RETURN_INTENT_KEY].analysis_id == "analysis-A"

    b_draft_key = app._fit_details_scoped_state_key(
        "fit_df_convention",
        "analysis-B",
    )
    at.session_state["test_analysis_id"] = "analysis-B"
    at.session_state["test_data_context"] = "real"
    at.session_state["main_results_panel"] = "data"
    at.session_state["fit_df_method_method"] = "engine"
    at.session_state["fit_df_method_cap"] = 7.0
    at.session_state[b_draft_key] = {"method": "engine", "cap": 7.0}
    at.run()
    assert not at.exception
    assert (
        at.session_state[app._HELP_ROUTE_STATE_KEY].context_status
        is HelpContextStatus.STALE
    )

    at.button(key="mfrm_help_return_to_source").click()
    at.run()
    assert not at.exception
    route = at.session_state[app._HELP_ROUTE_STATE_KEY]
    assert route.is_open
    assert route.is_fallback
    assert route.fallback_reason is HelpFallbackReason.RETURN_TARGET_UNAVAILABLE
    assert app._HELP_RETURN_INTENT_KEY not in at.session_state
    assert app._HELP_SUPPRESS_ANALYSIS_TRIGGERS_ONCE_KEY not in at.session_state
    assert at.session_state["main_results_panel"] == "data"
    assert at.session_state["fit_df_method_method"] == "engine"
    assert at.session_state["fit_df_method_cap"] == 7.0
    assert at.session_state[b_draft_key] == {"method": "engine", "cap": 7.0}
    assert not any(
        element.value == _locale_value("en", "help_nav.returned_status")
        for element in at.caption
    )
    assert not _return_focus_controllers(at)


def test_new_global_help_discards_unrendered_return_journey(monkeypatch):
    context = app._help_contract.HelpContextBinding(
        data_context="sample",
        analysis_id="analysis-fit-42",
        analysis_phase="fitted",
    )
    monkeypatch.setattr(app, "_current_fitted_help_context", lambda: context)
    harness = """
import streamlit_app as app

app._ensure_language_state()
state = app.render_persistent_help_launcher()
if state.is_open:
    app.render_persistent_help_surface(state)
else:
    if app.st.session_state.get("render_return_target", True):
        app.render_help_return_marker(
            "target.results.figure.fit_scatter",
            "focus.results.figure.fit_scatter",
        )
    app.render_help_popover("fit_scatter")
"""
    at = AppTest.from_string(harness, default_timeout=20)
    at.session_state["lang"] = "en"
    at.session_state["app_view_density"] = "Full"
    at.session_state["fit_df_method_method"] = "facets"
    at.session_state["fit_df_method_cap"] = 9.0
    at.run()
    assert not at.exception

    at.button(key="mfrm_popover_open_full_guide_fit_scatter").click()
    at.run()
    assert not at.exception
    at.session_state["render_return_target"] = False
    at.button(key="mfrm_help_return_to_source").click()
    at.run()
    assert not at.exception
    pending = at.session_state[app._HELP_RETURN_INTENT_KEY]
    assert pending.phase == "awaiting_render"
    assert at.session_state[app._HELP_SUPPRESS_ANALYSIS_TRIGGERS_ONCE_KEY] is True

    at.button(key="mfrm_help_open_global").click()
    at.run()
    assert not at.exception
    global_route = at.session_state[app._HELP_ROUTE_STATE_KEY]
    assert global_route.is_open
    assert not global_route.is_fallback
    assert global_route.help_link_id == app._GLOBAL_HELP_LINK_ID
    assert app._HELP_RETURN_INTENT_KEY not in at.session_state
    assert app._HELP_SUPPRESS_ANALYSIS_TRIGGERS_ONCE_KEY not in at.session_state
    assert app._HELP_RETURNED_STATUS_KEY not in at.session_state

    at.button(key="mfrm_help_return_to_source").click()
    at.run()
    assert not at.exception
    closed = at.session_state[app._HELP_ROUTE_STATE_KEY]
    assert not closed.is_open
    assert not closed.is_fallback


def test_global_help_preserves_fit_detail_view_drafts_across_hidden_rerun():
    harness = """
import pandas as pd
import streamlit_app as app

app._ensure_language_state()
state = app.render_persistent_help_launcher()
if state.is_open:
    app.render_persistent_help_surface(state)
else:
    view_token = app.st.session_state.get(
        "test_fit_view_token",
        "analysis-global-help-test",
    )
    diagnostics = {
        "fit_df_method": "facets",
        "facets_zstd_cap": 9.0,
        "fit": pd.DataFrame(),
        "overall_fit": pd.DataFrame(),
    }
    app.render_fit_df_method_controls(
        diagnostics,
        view_state_token=view_token,
    )
    app._draw_misfit_ranking(
        pd.DataFrame({
            "Facet": ["Rater"] * 8,
            "Level": [f"R{i}" for i in range(8)],
            "InfitZSTD": [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5],
            "OutfitZSTD": [0.1, 0.6, 1.1, 1.6, 2.1, 2.6, 3.1, 3.6],
        }),
        view_state_token=view_token,
    )
"""
    at = AppTest.from_string(harness, default_timeout=30)
    at.session_state["lang"] = "en"
    at.run()
    assert not at.exception
    at.radio(key="fit_df_method_method").set_value("both")
    at.run()
    at.number_input(key="fit_df_method_cap").set_value(12.5)
    at.run()
    at.slider(key="misfit_top_n").set_value(6)
    at.run()
    at.slider(key="misfit_threshold").set_value(3.0)
    at.run()
    assert not at.exception

    at.button(key="mfrm_help_open_global").click()
    at.run()
    assert not at.exception
    assert "fit_df_method_method" not in at.session_state
    assert "misfit_threshold" not in at.session_state
    at.session_state["lang"] = "ja"
    at.run()
    assert not at.exception
    assert "fit_df_method_method" not in at.session_state
    assert "misfit_threshold" not in at.session_state

    at.button(key="mfrm_help_return_to_source").click()
    at.run()
    assert not at.exception
    assert at.session_state["fit_df_method_method"] == "both"
    assert at.session_state["fit_df_method_cap"] == 12.5
    assert at.session_state["misfit_top_n"] == 6
    assert at.session_state["misfit_threshold"] == 3.0

    at.session_state["test_fit_view_token"] = "analysis-B"
    at.run()
    assert not at.exception
    assert at.session_state["fit_df_method_method"] == "facets"
    assert at.session_state["fit_df_method_cap"] == 9.0
    assert at.session_state["misfit_top_n"] == 8
    assert at.session_state["misfit_threshold"] == 2.0


def test_mml_sensitivity_draft_survives_contextual_help_widget_cleanup(
    monkeypatch,
):
    context = app._help_contract.HelpContextBinding(
        data_context="sample",
        analysis_id="analysis-fit-42",
        analysis_phase="fitted",
    )
    monkeypatch.setattr(app, "_current_fitted_help_context", lambda: context)
    monkeypatch.setattr(
        app,
        "build_mml_prior_sensitivity_plan",
        lambda _result: pd.DataFrame(),
    )
    harness = """
import pandas as pd
import streamlit_app as app

app._ensure_language_state()
if "test_mml_result" not in app.st.session_state:
    app.st.session_state["test_mml_result"] = {
        "config": {
            "method": "MML",
            "estimate_population_sd": False,
            "maxit": 200,
            "population_prior_sd": 1.0,
        }
    }
result = app.st.session_state["test_mml_result"]
state = app._get_help_route_state()
if state.is_open:
    app.render_persistent_help_surface(state)
else:
    app._draw_fit_scatter(pd.DataFrame({
        "Facet": ["Rater", "Rater"],
        "Level": ["R1", "R2"],
        "Infit": [0.9, 1.1],
        "Outfit": [0.8, 1.2],
    }))
    app.show_mml_prior_sd_sensitivity_section(result)
"""
    at = AppTest.from_string(harness, default_timeout=30)
    at.session_state["lang"] = "en"
    at.session_state["app_view_density"] = "Full"
    at.session_state["fit_df_method_method"] = "both"
    at.session_state["fit_df_method_cap"] = 12.5
    at.run()
    assert not at.exception
    result = at.session_state["test_mml_result"]
    suffix = str(id(result))

    at.multiselect(
        key=f"mml_prior_sensitivity_multipliers::{suffix}"
    ).set_value([0.5, 1.0, 2.0])
    at.run()
    at.number_input(
        key=f"mml_prior_sensitivity_maxit::{suffix}"
    ).set_value(35)
    at.run()
    at.selectbox(
        key=f"mml_prior_sensitivity_reltol::{suffix}"
    ).set_value(1e-5)
    at.run()
    assert not at.exception
    draft_key = app._mml_prior_sensitivity_draft_key(result)
    assert at.session_state[draft_key] == {
        "multipliers": (0.5, 1.0, 2.0),
        "maxit": 35,
        "reltol": 1e-5,
    }

    at.button(key="mfrm_popover_open_full_guide_fit_scatter").click()
    at.run()
    assert not at.exception
    assert f"mml_prior_sensitivity_multipliers::{suffix}" not in at.session_state
    assert f"mml_prior_sensitivity_maxit::{suffix}" not in at.session_state
    assert f"mml_prior_sensitivity_reltol::{suffix}" not in at.session_state
    assert at.session_state[draft_key]["maxit"] == 35
    at.session_state["lang"] = "ja"
    at.run()
    assert not at.exception
    assert f"mml_prior_sensitivity_multipliers::{suffix}" not in at.session_state
    assert at.session_state[draft_key] == {
        "multipliers": (0.5, 1.0, 2.0),
        "maxit": 35,
        "reltol": 1e-5,
    }

    at.button(key="mfrm_help_return_to_source").click()
    at.run()
    assert not at.exception
    assert at.session_state[f"mml_prior_sensitivity_multipliers::{suffix}"] == [
        0.5,
        1.0,
        2.0,
    ]
    assert at.session_state[f"mml_prior_sensitivity_maxit::{suffix}"] == 35
    assert at.session_state[f"mml_prior_sensitivity_reltol::{suffix}"] == 1e-5


def test_mml_sensitivity_draft_validation_fails_closed():
    config = {"maxit": 200}
    valid = {
        "multipliers": (0.5, 1.0, 2.0),
        "maxit": 35,
        "reltol": 1e-5,
    }
    assert app._validated_mml_prior_sensitivity_draft(valid, config=config) == valid
    assert (
        app._validated_mml_prior_sensitivity_draft(
            {**valid, "multipliers": (0.5, 0.5)},
            config=config,
        )
        is None
    )
    assert app._validated_fit_df_view_draft(
        {"method": "both", "cap": 12.5}
    ) == {"method": "both", "cap": 12.5}
    assert (
        app._validated_fit_df_view_draft(
            {"method": "both", "cap": float("inf")}
        )
        is None
    )
    assert (
        app._validated_mml_prior_sensitivity_draft(
            {**valid, "maxit": 201},
            config=config,
        )
        is None
    )
    assert (
        app._validated_mml_prior_sensitivity_draft(
            {**valid, "reltol": float("nan")},
            config=config,
        )
        is None
    )


def test_real_fit_scatter_help_return_preserves_fit_and_pending_triggers(
    monkeypatch,
):
    at = AppTest.from_string(
        "import streamlit_app as app\napp.main()",
        default_timeout=120,
    ).run()
    assert not at.exception
    at.button(key="onboarding_skip_guide").click()
    at.run()
    at.radio(key="app_view_density").set_value("Full")
    at.run()
    at.button(key="facets_mode_run_primary").click()
    at.run(timeout=120)
    assert not at.exception
    output_before = at.session_state["facets_mode_output"]
    identity_before = app.build_result_analysis_identity(
        output_before["result"]
    ).analysis_id

    at.selectbox(key="main_results_panel").set_value("fit_details")
    at.run()
    assert not at.exception
    at.radio(key="fit_df_method_method").set_value("both")
    at.run()
    assert not at.exception
    at.number_input(key="fit_df_method_cap").set_value(12.5)
    at.run()
    assert not at.exception
    at.slider(key="misfit_top_n").set_value(6)
    at.run()
    assert not at.exception
    at.slider(key="misfit_threshold").set_value(3.0)
    at.run()
    assert not at.exception
    assert any(
        button.key == "mfrm_popover_open_full_guide_fit_scatter"
        for button in at.button
    )

    calls = _block_analysis_calls(monkeypatch)
    at.session_state["_facets_mode_force_rerun"] = True
    at.session_state["_onboarding_quickstart_fired"] = True
    at.button(key="mfrm_popover_open_full_guide_fit_scatter").click()
    at.run()
    assert not at.exception
    assert calls == {"estimate": 0, "refit": 0}
    assert at.session_state[app._HELP_ROUTE_STATE_KEY].context_status is HelpContextStatus.CURRENT
    pending_intent = at.session_state[app._HELP_RETURN_INTENT_KEY]
    assert pending_intent.source_view_state == (
        ("fit_df_method_method", "both"),
        ("fit_df_method_cap", 12.5),
    )
    assert "fit_df_method_method" not in at.session_state
    assert "fit_df_method_cap" not in at.session_state
    assert "misfit_top_n" not in at.session_state
    assert "misfit_threshold" not in at.session_state
    view_state_token = app._fit_details_view_state_token(
        at.session_state["facets_mode_output"]["result"]
    )
    assert view_state_token == identity_before
    assert at.session_state[
        app._fit_details_scoped_state_key(
            app._MISFIT_TOP_N_VIEW_STATE_KEY,
            view_state_token,
        )
    ] == 6
    assert at.session_state[
        app._fit_details_scoped_state_key(
            app._MISFIT_THRESHOLD_VIEW_STATE_KEY,
            view_state_token,
        )
    ] == 3.0

    at.radio(key="lang").set_value("ja")
    at.run()
    assert not at.exception
    assert calls == {"estimate": 0, "refit": 0}

    at.button(key="mfrm_help_return_to_source").click()
    at.run()
    assert not at.exception
    assert calls == {"estimate": 0, "refit": 0}
    assert at.session_state["_facets_mode_force_rerun"] is True
    assert at.session_state["_onboarding_quickstart_fired"] is True
    assert at.session_state["main_results_panel"] == "fit_details"
    assert at.session_state["fit_df_method_method"] == "both"
    assert at.session_state["fit_df_method_cap"] == 12.5
    assert at.session_state["misfit_top_n"] == 6
    assert at.session_state["misfit_threshold"] == 3.0
    assert app._HELP_RETURN_INTENT_KEY not in at.session_state
    assert app._HELP_SUPPRESS_ANALYSIS_TRIGGERS_ONCE_KEY not in at.session_state
    identity_after = app.build_result_analysis_identity(
        at.session_state["facets_mode_output"]["result"]
    ).analysis_id
    assert identity_after == identity_before
    assert any(
        element.value == _locale_value("ja", "help_nav.returned_status")
        for element in at.caption
    )
