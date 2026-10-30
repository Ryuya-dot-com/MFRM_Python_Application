"""Mathematical appropriateness and UX contract tests for Desktop G/D-Study & terminology."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from mfrm_app import generalizability
from mfrm_app.desktop_analysis import (
    _gd_settings,
    _run_gd,
    prepare_gd_input_data,
    preview_gd_design,
    run_analysis,
)
from mfrm_app.desktop_guidance import gd_two_facet_example, writing_example
from mfrm_app.desktop_i18n import EN, GUIDES_EN, TEMPLATES, text
from mfrm_app.desktop_outputs import analysis_figure, current_tables, figure_choices, report_html
from mfrm_app.project import content_fingerprint, load_project, save_project


def test_desktop_i18n_bijection_and_terminology_unification():
    # EN must be a strict 1-to-1 bijection so text(msg, 'ja') reverses without collision.
    en_values = list(EN.values())
    duplicates = [v for v in set(en_values) if en_values.count(v) > 1]
    assert not duplicates, f"Duplicate English values break reverse translation: {duplicates}"

    for ja_msg, en_msg in EN.items():
        assert text(ja_msg, "en") == en_msg
        assert text(en_msg, "ja") == ja_msg
        assert "ファセット" not in ja_msg
        for legacy in ("G/D研究", "G／D研究", "G研究", "D研究", "G/D Study"):
            assert legacy not in ja_msg

    sample_args = {
        "version": "0.2.20-beta",
        "label": "demo.csv",
        "rows": "240",
        "columns": "4",
        "index": "1",
        "facet": "Rater",
        "count": "0",
        "x": "Rater",
        "series": "Task",
        "pair": "Rater × Task",
        "bound": "10",
        "path": "out.csv",
        "target": "Rater",
        "scope": "",
        "level": "95",
        "row": "3",
        "minimum": "0",
        "maximum": "5",
        "excluded": "2",
    }
    for ja_tpl, en_tpl in TEMPLATES:
        assert "ファセット" not in ja_tpl
        ja_filled = ja_tpl.format(**{k: v for k, v in sample_args.items() if "{" + k + "}" in ja_tpl})
        en_filled = text(ja_filled, "en")
        assert en_filled == en_tpl.format(**{k: v for k, v in sample_args.items() if "{" + k + "}" in en_tpl})
        assert text(en_filled, "ja") == ja_filled

    for guide_text in GUIDES_EN.values():
        assert "G/D-Study" in guide_text
        assert "proficiency" not in guide_text.lower()
        assert "strictness" not in guide_text.lower()


def test_gd_two_facet_example_is_balanced_and_immediately_analyzable():
    df = gd_two_facet_example()
    assert list(df.columns) == ["Person", "Rater", "Task", "Score"]
    assert len(df) == 30 * 4 * 2

    preview = preview_gd_design(
        df,
        object_facet="Person",
        score_col="Score",
        random_facets=["Rater", "Task"],
        fixed_condition=None,
    )
    assert preview["ok"] is True
    assert preview["status"] == "ready"
    assert preview["reason_code"] == ""

    analysis = run_analysis(
        df,
        workflow="gd",
        gd_settings={
            "object_facet": "Person",
            "score_col": "Score",
            "random_facets": ["Rater", "Task"],
            "n_max_per_facet": 8,
            "target_population": "Synthetic writing benchmark",
            "selected_plan": None,
        },
    )
    res = analysis["gd_analysis"]["result"]
    assert res["available"] is True
    assert res["design"]["n_observations"] == 240
    assert 0.0 <= res["selected_plan"]["G"] <= 1.0
    assert 0.0 <= res["selected_plan"]["Phi"] <= 1.0


def test_fixed_facet_composite_sum_vs_mean_brennan_scaling_invariance():
    """Brennan (2001, Ch. 2 & §3.4): across K=4 balanced fixed levels,
    SS(sum) = K^2 SS(mean), MS(sum) = K^2 MS(mean), Var(sum) = K^2 Var(mean),
    ErrVar(sum) = K^2 ErrVar(mean), SEM(sum) = K SEM(mean), and G(sum) == G(mean), Phi(sum) == Phi(mean)."""
    df = writing_example()
    k_levels = 4

    sum_entry = _run_gd(
        df,
        _gd_settings(
            {
                "object_facet": "Person",
                "score_col": "Score",
                "random_facets": ["Rater", "Task"],
                "n_max_per_facet": 6,
                "target_population": "",
                "selected_plan": {"Rater": 3, "Task": 2},
                "fixed_condition": {"facet": "Criterion", "mode": "sum", "level": None},
            }
        ),
    )
    mean_entry = _run_gd(
        df,
        _gd_settings(
            {
                "object_facet": "Person",
                "score_col": "Score",
                "random_facets": ["Rater", "Task"],
                "n_max_per_facet": 6,
                "target_population": "",
                "selected_plan": {"Rater": 3, "Task": 2},
                "fixed_condition": {"facet": "Criterion", "mode": "mean", "level": None},
            }
        ),
    )
    res_sum = sum_entry["result"]
    res_mean = mean_entry["result"]
    assert res_sum["available"] is True and res_mean["available"] is True
    assert res_sum["design"]["n_observations"] == 240
    assert res_mean["design"]["n_observations"] == 240

    vc_sum = res_sum["variance_components"].set_index("Source")
    vc_mean = res_mean["variance_components"].set_index("Source")
    assert list(vc_sum.index) == list(vc_mean.index)
    for col in ("SumSquares", "MeanSquares", "RawVariance", "Variance"):
        np.testing.assert_allclose(
            vc_sum[col].to_numpy(dtype=float),
            (k_levels**2) * vc_mean[col].to_numpy(dtype=float),
            rtol=1e-12,
            atol=1e-12,
        )
    np.testing.assert_allclose(
        vc_sum["ProportionVariance"].to_numpy(dtype=float),
        vc_mean["ProportionVariance"].to_numpy(dtype=float),
        rtol=1e-12,
        atol=1e-12,
    )

    sp_sum = res_sum["selected_plan"]
    sp_mean = res_mean["selected_plan"]
    np.testing.assert_allclose(
        sp_sum["RelativeErrorVariance"], (k_levels**2) * sp_mean["RelativeErrorVariance"], rtol=1e-12, atol=1e-12
    )
    np.testing.assert_allclose(
        sp_sum["AbsoluteErrorVariance"], (k_levels**2) * sp_mean["AbsoluteErrorVariance"], rtol=1e-12, atol=1e-12
    )
    np.testing.assert_allclose(
        np.sqrt(sp_sum["RelativeErrorVariance"]), k_levels * np.sqrt(sp_mean["RelativeErrorVariance"]), rtol=1e-12
    )
    np.testing.assert_allclose(
        np.sqrt(sp_sum["AbsoluteErrorVariance"]), k_levels * np.sqrt(sp_mean["AbsoluteErrorVariance"]), rtol=1e-12
    )

    ds_sum = res_sum["d_study"].sort_values(["Rater", "Task"]).reset_index(drop=True)
    ds_mean = res_mean["d_study"].sort_values(["Rater", "Task"]).reset_index(drop=True)
    np.testing.assert_allclose(ds_sum["G"], ds_mean["G"], rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(ds_sum["Phi"], ds_mean["Phi"], rtol=1e-12, atol=1e-12)


def test_fixed_facet_conditioning_rejects_incomplete_or_duplicate_elementary_cells():
    """Composite sum/mean must NEVER silently average unequal fixed-level counts."""
    df = writing_example()

    def _cfg(mode):
        return _gd_settings(
            {
                "object_facet": "Person",
                "score_col": "Score",
                "random_facets": ["Rater", "Task"],
                "n_max_per_facet": 5,
                "target_population": "",
                "selected_plan": None,
                "fixed_condition": {"facet": "Criterion", "mode": mode, "level": None},
            }
        )

    # 1. Incomplete elementary cell: drop a single rating row (959 instead of 960).
    incomplete = df.iloc[1:].reset_index(drop=True)
    for mode in ("sum", "mean"):
        cfg = _cfg(mode)
        _, prep_reason, _ = prepare_gd_input_data(incomplete, cfg)
        assert prep_reason == "missing_cells"

        entry = _run_gd(incomplete, cfg)
        assert entry["result"]["available"] is False
        assert entry["result"]["reason_code"] == "missing_cells"

    # 2. Duplicate elementary cell: append a duplicate row.
    duplicated = pd.concat([df, df.iloc[[0]]], ignore_index=True)
    _, prep_dup, _ = prepare_gd_input_data(duplicated, _cfg("mean"))
    assert prep_dup == "duplicate_cells"

    # 3. Non-unit Weight column must be rejected before aggregation.
    weighted = df.copy()
    weighted["Weight"] = 1.0
    weighted.loc[0, "Weight"] = 2.0
    _, prep_w, _ = prepare_gd_input_data(weighted, _cfg("sum"))
    assert prep_w == "unsupported_weights"

    # 4. Single-level fixed facet cannot form a multi-level composite.
    one_level = df[df["Criterion"].astype(str) == "Content"].reset_index(drop=True)
    _, prep_one, _ = prepare_gd_input_data(one_level, _cfg("sum"))
    assert prep_one == "insufficient_levels"


def test_single_level_fixed_condition_matches_manual_slice_and_preflight():
    df = writing_example()

    # Without fixed_condition, selecting Rater & Task on 3-facet data warns about Criterion (4 levels).
    warn_info = preview_gd_design(
        df,
        object_facet="Person",
        score_col="Score",
        random_facets=["Rater", "Task"],
        fixed_condition=None,
    )
    assert warn_info["ok"] is False
    assert warn_info["reason_code"] == "duplicate_cells"
    assert ("Criterion", 4) in warn_info["extra_facets"]
    assert "Criterion" in warn_info["summary_ja"]
    assert "Criterion" in warn_info["summary_en"]

    ok_info = preview_gd_design(
        df,
        object_facet="Person",
        score_col="Score",
        random_facets=["Rater", "Task"],
        fixed_condition={"facet": "Criterion", "mode": "level", "level": "Content"},
    )
    assert ok_info["ok"] is True

    cond_entry = _run_gd(
        df,
        _gd_settings(
            {
                "object_facet": "Person",
                "score_col": "Score",
                "random_facets": ["Rater", "Task"],
                "n_max_per_facet": 6,
                "target_population": "",
                "selected_plan": None,
                "fixed_condition": {"facet": "Criterion", "mode": "level", "level": "Content"},
            }
        ),
    )
    manual_slice = df[df["Criterion"].astype(str) == "Content"][["Person", "Rater", "Task", "Score"]].reset_index(drop=True)
    manual_entry = _run_gd(
        manual_slice,
        _gd_settings(
            {
                "object_facet": "Person",
                "score_col": "Score",
                "random_facets": ["Rater", "Task"],
                "n_max_per_facet": 6,
                "target_population": "",
                "selected_plan": None,
            }
        ),
    )
    pd.testing.assert_frame_equal(
        cond_entry["result"]["variance_components"],
        manual_entry["result"]["variance_components"],
    )
    pd.testing.assert_frame_equal(
        cond_entry["result"]["d_study"],
        manual_entry["result"]["d_study"],
    )


def test_gd_swap_figure_and_mfrm_project_reopen_without_recomputation(tmp_path, monkeypatch):
    df = writing_example()
    analysis = run_analysis(
        df,
        workflow="gd",
        gd_settings={
            "object_facet": "Person",
            "score_col": "Score",
            "random_facets": ["Rater", "Task"],
            "n_max_per_facet": 6,
            "target_population": "University writing placement",
            "selected_plan": {"Rater": 2, "Task": 3},
            "fixed_condition": {"facet": "Criterion", "mode": "mean", "level": None},
        },
    )
    choices = figure_choices(analysis)
    assert choices == [("D-Study · G / Phi", {"kind": "gd"})]

    fig_default = analysis_figure(analysis, kind="gd", gd_swap=False)
    fig_swapped = analysis_figure(analysis, kind="gd", gd_swap=True)
    assert fig_default.axes[0].get_xlabel() == "Rater count"
    assert fig_swapped.axes[0].get_xlabel() == "Task count"

    saved = dict(
        analysis,
        view_settings={
            "figure": "D-Study · G / Phi",
            "facet": "すべて",
            "statistic": "Infit MnSq",
            "presentation": "explained",
            "language": "en",
            "gd_swap": "1",
        },
    )
    saved["content_fingerprint"] = content_fingerprint(saved)
    project_path = tmp_path / "gd_mean_swapped.mfrm"
    save_project(project_path, saved)

    monkeypatch.setattr(
        generalizability,
        "compute_generalizability_from_data",
        lambda *a, **kw: (_ for _ in ()).throw(AssertionError("Must not recompute G/D on reopen")),
    )
    reloaded = load_project(project_path)
    assert reloaded["gd_analysis"]["settings"]["fixed_condition"] == {
        "facet": "Criterion",
        "mode": "mean",
        "level": None,
    }
    assert reloaded["view_settings"]["gd_swap"] == "1"
    tables = current_tables(reloaded)
    assert len(tables) == 5
    html_doc = report_html(reloaded)
    assert "G/D-Study" in html_doc
    assert "Task count" in html_doc


def test_desktop_app_gd_ux_controls_and_language_roundtrip():
    tk = pytest.importorskip("tkinter")
    from desktop_app import DesktopApp

    try:
        root = tk.Tk()
    except tk.TclError:
        pytest.skip("Tk display unavailable")
    root.withdraw()
    try:
        app = DesktopApp(root)
        # 1. One-click G/D 2-facet sample switches to G/D-Study, excludes Person/Score, and passes preflight.
        app.use_gd_demo()
        assert app.workflow.get() == "G/D-Study"
        facet_items = [app.gd_facets.get(i) for i in range(app.gd_facets.size())]
        assert "Person" not in facet_items
        assert "Score" not in facet_items
        assert [app.gd_facets.get(i) for i in app.gd_facets.curselection()] == ["Rater", "Task"]
        assert "ready" in app.gd_preflight.get().lower() or "ok" in app.gd_preflight.get().lower()
        req = app.settings()
        assert req["workflow"] == "gd"
        assert req["gd_settings"]["random_facets"] == ["Rater", "Task"]

        # 2. Swap random facet order updates settings order immediately.
        app.swap_gd_random_facets()
        assert app.settings()["gd_settings"]["random_facets"] == ["Task", "Rater"]

        # 3. Loading 3-facet sample warns before Analyze when fixed_condition is None, and resolves when Mean selected.
        app.use_demo()
        app.workflow.set("G/D-Study")
        app.change_workflow()
        assert "Criterion" in app.gd_preflight.get()
        app._set_gd_fixed_mode_code("mean")
        app.update_gd_facet_choices()
        assert "ready" in app.gd_preflight.get().lower() or "ok" in app.gd_preflight.get().lower()

        # 4. Language switch en -> ja -> en preserves fixed_condition and valid settings.
        app.language_choice.set("日本語")
        app.change_language()
        assert "設計確認OK" in app.gd_preflight.get()
        assert app.settings()["gd_settings"]["fixed_condition"]["mode"] == "mean"
        app.language_choice.set("English")
        app.change_language()
        assert "ready" in app.gd_preflight.get().lower()
        assert app.settings()["gd_settings"]["fixed_condition"]["mode"] == "mean"
        app.worker.close()
    finally:
        root.destroy()
