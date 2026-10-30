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
        "reason": "欠測セル",
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


@pytest.mark.parametrize(
    "selected_plan,n_max",
    [
        ({"Rater": 7, "Task": 5}, 10),   # Previously failed: 7 absent from discrete grid & Task capped at observed*2=4
        ({"Rater": 12, "Task": 9}, 15),  # Previously failed: grid capped at 10
        ({"Rater": 14, "Task": 6}, 8),   # Custom plan exceeding n_max_per_facet still evaluates & appends cleanly
    ],
)
def test_phase_c1_custom_d_study_plans_evaluate_exact_closed_form(selected_plan, n_max, tmp_path):
    """C-1: Any positive integer custom D-Study plan evaluates Brennan Eq. 3.18 without ValueError."""
    df = gd_two_facet_example()
    analysis = run_analysis(
        df,
        workflow="gd",
        gd_settings={
            "object_facet": "Person",
            "score_col": "Score",
            "random_facets": ["Rater", "Task"],
            "n_max_per_facet": n_max,
            "target_population": "Placement cohort",
            "selected_plan": selected_plan,
        },
    )
    res = analysis["gd_analysis"]["result"]
    assert res["available"] is True
    sel = res["selected_plan"]
    assert sel["Rater"] == selected_plan["Rater"]
    assert sel["Task"] == selected_plan["Task"]

    vc = res["variance_components"].set_index("Source")["Variance"].to_dict()
    n_r, n_t = selected_plan["Rater"], selected_plan["Task"]
    exp_rel = vc["Person:Rater"] / n_r + vc["Person:Task"] / n_t + vc["Residual"] / (n_r * n_t)
    exp_abs = exp_rel + vc["Rater"] / n_r + vc["Task"] / n_t + vc["Rater:Task"] / (n_r * n_t)
    exp_g = vc["Person"] / (vc["Person"] + exp_rel)
    exp_phi = vc["Person"] / (vc["Person"] + exp_abs)

    assert sel["G"] == pytest.approx(exp_g, abs=1e-12)
    assert sel["Phi"] == pytest.approx(exp_phi, abs=1e-12)
    assert sel["RelativeErrorVariance"] == pytest.approx(exp_rel, abs=1e-12)
    assert sel["AbsoluteErrorVariance"] == pytest.approx(exp_abs, abs=1e-12)
    assert sel["RelativeSEM"] == pytest.approx(np.sqrt(exp_rel), abs=1e-12)
    assert sel["AbsoluteSEM"] == pytest.approx(np.sqrt(exp_abs), abs=1e-12)

    hit = res["d_study"][
        res["d_study"]["Rater"].eq(n_r) & res["d_study"]["Task"].eq(n_t)
    ]
    assert len(hit) == 1
    assert float(hit["G"].iloc[0]) == pytest.approx(exp_g, abs=1e-12)

    fig = analysis_figure(analysis, kind="gd")
    assert len(fig.axes[0].get_legend().get_texts()) <= 13  # Up to 6 series x (G, Phi) + Selected plan

    project_path = tmp_path / f"custom_{n_r}_{n_t}.mfrm"
    save_project(project_path, analysis)
    reloaded = load_project(project_path)
    assert reloaded["gd_analysis"]["result"]["selected_plan"]["Rater"] == n_r


def test_phase_c2_orthogonal_tensor_projection_eliminates_ss_cancellation():
    """C-2: High main-effect-to-residual dynamic range and large constant shifts
    must not suffer IEEE-754 subtraction cancellation in 2-way or 3-way ANOVA."""
    import itertools

    # 1) High main-effect-to-residual ratio without base offset:
    # Y = 1e6*p + 3*r + 2*t + 0.01*p*r*t -> SS_total ~ 8e12 while SS_residual = 0.0008 (16-order ratio!).
    # Old subtraction (ss_total - ss_p - ...) suffered >22% cancellation error here;
    # orthogonal projection recovers exact float64 input residual to < 1e-7 relative error.
    rows_3way = [
        {
            "Person": f"P{p:+d}",
            "Rater": f"R{r:+d}",
            "Task": f"T{t:+d}",
            "Score": 1.0e6 * p + 3.0 * r + 2.0 * t + 0.01 * p * r * t,
        }
        for p, r, t in itertools.product((-1, 1), repeat=3)
    ]
    vc_3 = generalizability._crossed_anova_variance_components(
        pd.DataFrame(rows_3way),
        object_facet="Person",
        random_facets=["Rater", "Task"],
        score_col="Score",
    ).set_index("Source")
    assert vc_3.loc["Residual", "SumSquares"] == pytest.approx(0.0008, rel=1e-7)
    assert vc_3.loc["Residual", "RawVariance"] == pytest.approx(0.0008, rel=1e-7)
    assert (vc_3["SumSquares"] >= 0.0).all()

    # 2) Shift invariance under a +1e9 constant offset (where ss_total ~ 8e18 wipes out 0.0008 completely under naive formulas):
    shifted_3way = pd.DataFrame(rows_3way).assign(Score=lambda d: d["Score"] + 1.0e9)
    vc_3_shift = generalizability._crossed_anova_variance_components(
        shifted_3way,
        object_facet="Person",
        random_facets=["Rater", "Task"],
        score_col="Score",
    ).set_index("Source")
    assert vc_3_shift.loc["Residual", "SumSquares"] == pytest.approx(0.0008, rel=1e-5)

    # 3) Exact orthogonal (-1, +1) 2-way construction:
    # Y = 1e6*p + 5*r + 0.01*p*r -> 4 cells, df_residual = 1, SS_residual = 4 * 1e-4 = 0.0004.
    rows_2way = [
        {
            "Person": f"P{p:+d}",
            "Rater": f"R{r:+d}",
            "Score": 1.0e6 * p + 5.0 * r + 0.01 * p * r,
        }
        for p, r in itertools.product((-1, 1), repeat=2)
    ]
    vc_2 = generalizability._crossed_anova_variance_components(
        pd.DataFrame(rows_2way),
        object_facet="Person",
        random_facets=["Rater"],
        score_col="Score",
    ).set_index("Source")
    assert vc_2.loc["Residual", "SumSquares"] == pytest.approx(0.0004, rel=1e-7)
    assert vc_2.loc["Residual", "RawVariance"] == pytest.approx(0.0004, rel=1e-7)
    assert (vc_2["SumSquares"] >= 0.0).all()


def test_phase_c3_web_and_desktop_share_single_gd_engine_and_fixed_conditioning():
    """C-3: streamlit_app delegates directly to mfrm_app.generalizability and matches Desktop on 3-facet conditioning."""
    import streamlit_app as web

    assert web._crossed_anova_variance_components is generalizability._crossed_anova_variance_components
    assert web._g_phi_from_variance_components is generalizability._g_phi_from_variance_components
    assert web.compute_generalizability_from_data is generalizability.compute_generalizability_from_data
    assert web.g_d_study_export_frame is generalizability.g_d_study_export_frame

    df = writing_example()
    fc = {"facet": "Criterion", "mode": "mean", "level": None}
    desk = run_analysis(
        df,
        workflow="gd",
        gd_settings={
            "object_facet": "Person",
            "score_col": "Score",
            "random_facets": ["Rater", "Task"],
            "n_max_per_facet": 8,
            "target_population": "Writing cohort",
            "selected_plan": None,
            "fixed_condition": fc,
        },
    )["gd_analysis"]["result"]
    web_raw = web.compute_generalizability_from_data(
        df,
        object_facet="Person",
        random_facets=["Rater", "Task"],
        score_col="Score",
        n_max_per_facet=8,
        target_population="Writing cohort",
        fixed_condition=fc,
    )
    pd.testing.assert_frame_equal(desk["variance_components"], web_raw["variance_components"])
    pd.testing.assert_frame_equal(desk["d_study"], web_raw["d_study"])
    assert desk["observed_coefficients"] == web_raw["observed_coefficients"]


def test_phase_c4_gd_table_order_localized_notes_and_contextual_banners():
    """C-4: Practical table order, section-specific localized notes, no orange banner on clean G/D,
    and specific banners when unavailable or when negative RawVariance occurs."""
    df = gd_two_facet_example()
    clean = run_analysis(
        df,
        workflow="gd",
        gd_settings={
            "object_facet": "Person",
            "score_col": "Score",
            "random_facets": ["Rater", "Task"],
            "n_max_per_facet": 6,
            "target_population": "Benchmark",
            "selected_plan": {"Rater": 5, "Task": 3},
        },
    )
    tables = current_tables(clean)
    titles = [t for t, _, _ in tables]
    notes = [n for _, _, n in tables]
    assert titles == [
        "G/D: observed_coefficients",
        "G/D: selected_d_study",
        "G/D: variance_components",
        "G/D: d_study",
        "G/D: metadata",
    ]
    assert len(set(notes)) == 5  # Each section has its own concise note (not 5 identical CI warnings)
    assert "Person-only resampling" not in " ".join(notes)
    assert {"RelativeSEM", "AbsoluteSEM"}.issubset(tables[0][1].columns)
    assert {"RelativeSEM", "AbsoluteSEM"}.issubset(tables[1][1].columns)

    ja_html = report_html(clean, language="ja")
    assert "G/D: 観測設計のG・Phiと誤差" in ja_html
    assert "G/D: 分散成分表" in ja_html
    assert "元の観測得点における観測設計の点推定です" in ja_html
    assert "<p>Person-only resampling" not in ja_html
    assert ja_html.count("Person-only resampling has not qualified") == 1  # Only inside the metadata audit Value cell

    tk = pytest.importorskip("tkinter")
    from desktop_app import DesktopApp

    try:
        root = tk.Tk()
    except tk.TclError:
        pytest.skip("Tk display unavailable")
    root.withdraw()
    try:
        app = DesktopApp(root, language="ja")
        # 1. Clean positive-variance G/D run shows NO orange top warning banner.
        app.show_analysis(clean)
        assert app.qualification.get() == ""

        # 2. Unavailable G/D run (missing cell) shows explicit unavailable reason in banner.
        missing = run_analysis(
            df.iloc[1:].reset_index(drop=True),
            workflow="gd",
            gd_settings={
                "object_facet": "Person",
                "score_col": "Score",
                "random_facets": ["Rater", "Task"],
            },
        )
        app.show_analysis(missing)
        assert "計算できませんでした" in app.qualification.get()

        # 3. Positive vs negative RawVariance attaches the diagnostic directly to G/D: variance_components
        #    without cluttering the top-of-window banner.
        pos_valid_df = pd.DataFrame(
            [
                {"Person": "P1", "Rater": "R1", "Score": 22.0},
                {"Person": "P1", "Rater": "R2", "Score": 14.0},
                {"Person": "P2", "Rater": "R1", "Score": 12.0},
                {"Person": "P2", "Rater": "R2", "Score": 8.0},
            ]
        )
        pos_run = run_analysis(
            pos_valid_df,
            workflow="gd",
            gd_settings={"object_facet": "Person", "score_col": "Score", "random_facets": ["Rater"]},
        )
        app.show_analysis(pos_run)
        assert app.qualification.get() == ""
        app.table_selector.current(2)
        app.show_table()
        assert "符号付きのモーメント法推定値" in app.table_note.get()

        neg_valid_df = pd.DataFrame(
            [
                {"Person": "P1", "Rater": "R1", "Score": 20.0},
                {"Person": "P1", "Rater": "R2", "Score": 16.0},
                {"Person": "P2", "Rater": "R1", "Score": 6.0},
                {"Person": "P2", "Rater": "R2", "Score": 10.0},
            ]
        )
        neg_run = run_analysis(
            neg_valid_df,
            workflow="gd",
            gd_settings={"object_facet": "Person", "score_col": "Score", "random_facets": ["Rater"]},
        )
        assert neg_run["gd_analysis"]["result"]["available"] is True
        app.show_analysis(neg_run)
        assert app.qualification.get() == ""
        app.table_selector.current(2)
        app.show_table()
        assert "負のRawVariance推定値が含まれており" in app.table_note.get()
        app.worker.close()
    finally:
        root.destroy()
