"""Regression and responsiveness contract tests for Streamlit P3-A optimizations."""
from __future__ import annotations

import inspect
from copy import deepcopy

import numpy as np
import pandas as pd
import pytest
from streamlit.testing.v1 import AppTest

import streamlit_app as app
from tests.test_help_contextual_app_integration import _block_analysis_calls, _locale_value


@pytest.mark.parametrize(
    "left,right",
    [
        (pd.DataFrame({"x": [None]}, dtype=object), pd.DataFrame({"x": ["<NA>"]}, dtype=object)),
        (pd.DataFrame({"x": [1]}, dtype=object), pd.DataFrame({"x": ["1"]}, dtype=object)),
        (pd.DataFrame({"x": [None]}, dtype=object), pd.DataFrame({"x": [np.nan]}, dtype=object)),
        (pd.DataFrame({"x": [None]}, dtype=object), pd.DataFrame({"x": [pd.NA]}, dtype=object)),
        (pd.DataFrame({1: [1]}), pd.DataFrame({"1": [1]})),
        (pd.DataFrame({"x": [1]}, index=[1]), pd.DataFrame({"x": [1]}, index=["1"])),
        (pd.DataFrame({"x": [1]}), pd.DataFrame({"x": [1]}).rename_axis("row")),
        (
            pd.DataFrame({"x": pd.Categorical(["a"], categories=["a", "b"])}),
            pd.DataFrame({"x": pd.Categorical(["a"], categories=["a", "c"])}),
        ),
        (
            pd.DataFrame({"x": pd.Categorical(["a"], ordered=False)}),
            pd.DataFrame({"x": pd.Categorical(["a"], ordered=True)}),
        ),
    ],
)
def test_fast_dataframe_cache_token_preserves_type_and_category_distinctions(left, right):
    app._cached_canonical_dataframe_fingerprint.clear()
    tok_left = app._fast_dataframe_cache_token(left)
    tok_right = app._fast_dataframe_cache_token(right)
    assert tok_left == app._fast_dataframe_cache_token(left.copy())
    assert tok_left != tok_right
    assert app._cached_canonical_dataframe_fingerprint(left, tok_left, length=64) == app.dataframe_fingerprint(
        left, length=64
    )
    assert app._cached_canonical_dataframe_fingerprint(right, tok_right, length=64) == app.dataframe_fingerprint(
        right, length=64
    )
    assert app._cached_canonical_dataframe_fingerprint(left, tok_left, length=64) != (
        app._cached_canonical_dataframe_fingerprint(right, tok_right, length=64)
    )


def test_cached_rerun_scans_avoid_duplicate_work_and_respect_locale(monkeypatch):
    app._cached_canonical_dataframe_fingerprint.clear()
    app._cached_auto_detect_column_roles.clear()
    app._cached_build_readiness_report.clear()

    df = app._make_self_test_gradient_data()
    tok = app._fast_dataframe_cache_token(df)
    counts = {"fp": 0, "roles": 0, "readiness": 0}

    orig_fp = app.dataframe_fingerprint
    orig_roles = app.auto_detect_column_roles
    orig_readiness = app.build_readiness_report

    def counted_fp(*args, **kwargs):
        counts["fp"] += 1
        return orig_fp(*args, **kwargs)

    def counted_roles(*args, **kwargs):
        counts["roles"] += 1
        return orig_roles(*args, **kwargs)

    def counted_readiness(*args, **kwargs):
        counts["readiness"] += 1
        return orig_readiness(*args, **kwargs)

    monkeypatch.setattr(app, "dataframe_fingerprint", counted_fp)
    monkeypatch.setattr(app, "auto_detect_column_roles", counted_roles)
    monkeypatch.setattr(app, "build_readiness_report", counted_readiness)

    fp_first = app._cached_canonical_dataframe_fingerprint(df, tok)
    fp_second = app._cached_canonical_dataframe_fingerprint(df, tok)
    roles_first = app._cached_auto_detect_column_roles(df, tok)
    roles_second = app._cached_auto_detect_column_roles(df, tok)
    assert fp_first == fp_second == orig_fp(df)
    assert roles_first == roles_second
    assert counts["fp"] == 1
    assert counts["roles"] == 1

    monkeypatch.setitem(app.st.session_state, "lang", "en")
    rep_en_1 = app._cached_build_readiness_report(
        df, tok, person_col="Person", score_col="Score", facet_cols_key=("Rater", "Task"), weight_col=None, lang="en"
    )
    rep_en_2 = app._cached_build_readiness_report(
        df, tok, person_col="Person", score_col="Score", facet_cols_key=("Rater", "Task"), weight_col=None, lang="en"
    )
    assert rep_en_1 == rep_en_2
    assert counts["readiness"] == 1

    monkeypatch.setitem(app.st.session_state, "lang", "ja")
    rep_ja = app._cached_build_readiness_report(
        df, tok, person_col="Person", score_col="Score", facet_cols_key=("Rater", "Task"), weight_col=None, lang="ja"
    )
    assert counts["readiness"] == 2
    assert rep_ja["checks"][0]["headline"] != rep_en_1["checks"][0]["headline"]


def test_presentation_scaling_refreshes_in_place_without_reestimation_while_maxit_stales(monkeypatch):
    at = AppTest.from_string("import streamlit_app as app\napp.main()", default_timeout=120)
    at.run(timeout=120)
    assert not at.exception
    at.button(key="onboarding_dismiss").click()
    at.run(timeout=120)
    data = app._make_self_test_gradient_data()
    at.text_area(key="paste_data_text").set_value(data.to_csv(index=False))
    at.run(timeout=120)
    at.radio(key="facets_mode_workflow_mode").set_value("Advanced controls")
    at.run(timeout=120)
    at.number_input(key="facets_mode_maxit").set_value(100)
    at.radio(key="lang").set_value("en")
    at.run(timeout=120)

    run_label = _locale_value("en", "sidebar_perf.run_button")
    next(b for b in at.sidebar.button if b.label == run_label).click()
    at.run(timeout=120)
    assert not at.exception

    output = at.session_state["facets_mode_output"]
    result = output["result"]
    analysis_id = app.build_result_analysis_identity(result).analysis_id
    base_person_measure = float(output["report_tables"]["Person"]["Measure"].iloc[0])

    calls = _block_analysis_calls(monkeypatch)
    monkeypatch.setattr(
        app,
        "mfrm_diagnostics",
        lambda *a, **k: (_ for _ in ()).throw(AssertionError("Unexpected mfrm_diagnostics rerun")),
    )

    # Presentation-only report scaling updates report_tables in place with zero re-estimation.
    at.number_input(key="facets_mode_umean").set_value(50.0)
    at.number_input(key="facets_mode_uscale").set_value(10.0)
    at.number_input(key="facets_mode_udecimals").set_value(3)
    at.run(timeout=120)
    assert not at.exception
    assert not any(b.key == "stale_rerun" for b in at.button)
    assert _locale_value("en", "app.stale_results_warning") not in [w.value for w in at.warning]
    assert at.session_state["facets_mode_output"] is output
    assert app.build_result_analysis_identity(output["result"]).analysis_id == analysis_id
    scaled_person_measure = float(output["report_tables"]["Person"]["Measure"].iloc[0])
    assert abs(scaled_person_measure - (50.0 + 10.0 * base_person_measure)) < 0.02
    assert output["_umean"] == 50.0
    assert output["_uscale"] == 10.0
    assert output["_udecimals"] == 3
    assert calls == {"estimate": 0, "refit": 0}

    # Estimator parameter changes still require explicit re-estimation and block stale results.
    at.number_input(key="facets_mode_maxit").set_value(120)
    at.run(timeout=120)
    assert not at.exception
    assert any(b.key == "stale_rerun" for b in at.button)
    assert _locale_value("en", "app.stale_results_warning") in [w.value for w in at.warning]
    assert calls == {"estimate": 0, "refit": 0}


def test_compact_history_snapshot_omits_bulky_row_files_and_restores_identically(monkeypatch):
    data = app._make_self_test_gradient_data()
    result = app.mfrm_estimate(data, "Person", ["Rater", "Task"], "Score", method="JMLE", maxit=80)
    diagnostics = app.mfrm_diagnostics(result, compute_pca=False)
    scorefile = app.compute_scorefile(result)
    residuals = app.compute_residual_file(result)
    full_output = {
        "result": result,
        "diagnostics": diagnostics,
        "scorefile": scorefile,
        "residuals": residuals,
        "model_type": "RSM",
        "_est_method": "JMLE",
        "facet_cols": ["Rater", "Task"],
    }

    state = {app._RUN_HISTORY_STORE_SNAPSHOTS_KEY: True}
    monkeypatch.setattr(app.st, "session_state", state)
    entry = app.record_run_in_history(output=deepcopy(full_output), elapsed_sec=0.4)
    snap = entry["output_snapshot"]
    assert entry["snapshot_available"] is True
    assert "scorefile" not in snap
    assert "residuals" not in snap
    assert snap.get("_history_omitted_scorefile") is True
    assert snap.get("_history_omitted_residuals") is True

    calls = _block_analysis_calls(monkeypatch)
    assert app.restore_run_from_history(entry["run_id"]) is True
    restored = state["facets_mode_output"]
    assert "_history_omitted_scorefile" not in restored
    assert "_history_omitted_residuals" not in restored
    pd.testing.assert_frame_equal(restored["scorefile"], scorefile, check_exact=True)
    pd.testing.assert_frame_equal(restored["residuals"], residuals, check_exact=True)
    assert calls == {"estimate": 0, "refit": 0}


def test_interactive_panel_fragments_preserve_wrapped_source_inspection():
    decorated = (
        app._render_guided_figures_section,
        app.show_wright_map_section,
        app.show_visuals_section,
        app.show_bias_section,
        app.show_categories_section,
        app.show_agreement_section,
        app.render_g_d_study_bundle,
    )
    for fn in decorated:
        assert hasattr(fn, "__wrapped__")
        assert f"def {fn.__name__}(" in inspect.getsource(fn)
