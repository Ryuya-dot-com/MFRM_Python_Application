"""Tests for Desktop manuscript diagnostics (agreement/alpha, optional residual PCA,
linear logit rescaling), Web fixed-N(0,1) joint interaction UI, research-option UI
separation, and bidirectional Web <-> Desktop `.mfrm` archive portability."""
from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import subprocess
import sys
from unittest.mock import patch

import matplotlib.figure
import numpy as np
import pandas as pd
import pytest
from streamlit.testing.v1 import AppTest

from desktop_app import DesktopApp
from mfrm_app import (
    basic_analysis as basic,
    desktop_analysis as desktop,
    desktop_examples as examples,
    desktop_i18n as i18n,
    desktop_outputs as outputs,
    engine,
    generalizability as gd,
    joint_interaction as joint,
    project,
    project_portability as portability,
    reporting_facts,
    structural_intervals as intervals,
)
import streamlit_app as web_app

REPO_ROOT = Path(__file__).resolve().parents[1]


def _small_crossed_df(n_person: int = 16, n_rater: int = 3, n_task: int = 2) -> pd.DataFrame:
    rows = []
    for p in range(n_person):
        for r in range(n_rater):
            for t in range(n_task):
                score = 1 + ((2 * p + 3 * r + t) % 4)
                rows.append({
                    "Person": f"P{p + 1:02d}",
                    "Rater": f"R{r + 1}",
                    "Task": f"T{t + 1}",
                    "Score": score,
                })
    return pd.DataFrame(rows)


def _check_project_cli(mfrm_path: Path, out_dir: Path) -> dict:
    subprocess.run(
        [sys.executable, str(REPO_ROOT / "desktop_app.py"), "--check-project", str(mfrm_path), "--output-dir", str(out_dir)],
        cwd=REPO_ROOT,
        check=True,
    )
    return json.loads((out_dir / "reopen.json").read_text(encoding="utf-8"))


def test_interrater_agreement_and_reliability_match_core_and_respect_calibration_holds():
    identical_df = pd.DataFrame([
        {"Person": f"P{p:02d}", "Rater": r, "Task": f"T{t + 1}", "Score": 1 + ((p + t) % 4)}
        for p in range(10)
        for r in ("R1", "R2", "R3")
        for t in range(2)
    ])
    ident_fit = basic.run_basic_analysis(
        identical_df,
        person_col="Person",
        facet_cols=["Rater", "Task"],
        score_col="Score",
        rating_min=1,
        rating_max=4,
        model="RSM",
        method="JMLE",
    )
    ident_tables = dict((name, tbl) for name, tbl, _ in outputs.current_tables(ident_fit))
    agr = ident_tables["Inter-rater Agreement and Krippendorff's Alpha"]
    rater_row = agr.loc[agr["Facet"].eq("Rater")].iloc[0]
    assert float(rater_row["ExactPct"]) == pytest.approx(100.0)
    assert float(rater_row["AdjacentPct"]) == pytest.approx(100.0)
    assert float(rater_row["AlphaOrdinal"]) == pytest.approx(1.0)
    assert float(rater_row["AlphaInterval"]) == pytest.approx(1.0)
    assert float(rater_row["AlphaNominal"]) == pytest.approx(1.0)
    assert np.isfinite(float(rater_row["ExpectedExactPct"]))
    assert "Pairwise Agreement: Rater" in ident_tables

    # Compare JMLE (diagnostic reliability & model-expected exact agreement) vs MML (withheld).
    df = _small_crossed_df()
    jmle_an = desktop.run_analysis(
        df,
        person_col="Person",
        facet_cols=["Rater", "Task"],
        score_col="Score",
        rating_min=1,
        rating_max=4,
        model="RSM",
        method="JMLE",
    )
    mml_an = desktop.run_analysis(
        df,
        person_col="Person",
        facet_cols=["Rater", "Task"],
        score_col="Score",
        rating_min=1,
        rating_max=4,
        model="RSM",
        method="MML",
        quad_points=15,
        population_prior_sd=1.0,
        mml_engine="EM",
    )
    jmle_map = {name: tbl for name, tbl, _ in outputs.current_tables(jmle_an)}
    mml_map = {name: tbl for name, tbl, _ in outputs.current_tables(mml_an)}

    jmle_rel = jmle_map["Separation and Reliability by Facet"]
    mml_rel = mml_map["Separation and Reliability by Facet"]
    assert jmle_rel["Use"].eq("Diagnostic").all()
    assert pd.to_numeric(jmle_rel["Reliability"], errors="coerce").notna().all()
    assert mml_rel["Use"].eq("Withheld").all()
    for col in ("RMSE", "Separation", "Strata", "Reliability"):
        assert mml_rel[col].eq("").all()

    # Exact parity with reporting_facts.publication_reliability_table
    expected_pub_rel, _ = reporting_facts.publication_reliability_table(
        jmle_an["result"], jmle_an["diagnostics"], rel_note=""
    )
    pd.testing.assert_frame_equal(jmle_rel, expected_pub_rel)

    # Under MML calibration hold, observed exact/adjacent/alpha remain available while ExpectedExactPct is masked.
    mml_agr = mml_map["Inter-rater Agreement and Krippendorff's Alpha"]
    assert pd.to_numeric(mml_agr["ExactPct"], errors="coerce").notna().all()
    assert pd.to_numeric(mml_agr["AlphaOrdinal"], errors="coerce").notna().all()
    assert mml_agr["ExpectedExactPct"].eq("").all()
    mml_pairs = mml_map["Pairwise Agreement: Rater"]
    assert mml_pairs["ExpectedExactPct"].eq("").all()


def test_optional_residual_pca_math_scree_plot_and_zero_refit_reopen(tmp_path):
    df = examples.generate_clinical_4facet_df()
    an = desktop.run_analysis(
        df,
        person_col="Examinee",
        facet_cols=["Rater", "Station", "Competency", "Session"],
        score_col="Rating",
        rating_min=1,
        rating_max=5,
        model="RSM",
        method="JMLE",
        compute_pca=True,
    )
    assert an["fit_settings"]["compute_pca"] is True
    assert an["scope"]["pca"] is True
    assert an["diagnostics"]["pca_enabled"] is True

    tbls = {name: tbl for name, tbl, _ in outputs.current_tables(an)}
    assert "Residual PCA: Eigenvalues and Stability" in tbls
    assert "Residual PCA: First-Contrast Loadings" in tbls
    ev_tbl = tbls["Residual PCA: Eigenvalues and Stability"]
    ev_vec = np.asarray(an["diagnostics"]["pca"]["eigenvalues"], dtype=float)
    vp_vec = np.asarray(an["diagnostics"]["pca"]["variance_pct"], dtype=float)
    assert float(ev_vec.sum()) == pytest.approx(float(len(ev_vec)), rel=1e-6, abs=1e-6)
    assert float(vp_vec.sum()) == pytest.approx(100.0, abs=1e-4)
    assert int(ev_tbl.iloc[0]["ResidualColumns"]) == len(ev_vec)
    assert float(ev_tbl.iloc[0]["Eigenvalue"]) == pytest.approx(float(ev_vec[0]), abs=1e-3)

    choices = dict(outputs.figure_choices(an))
    assert "残差PCA · 固有値スクリー" in choices
    fig = outputs.analysis_figure(an, kind="pca")
    assert isinstance(fig, matplotlib.figure.Figure)

    mfrm_path = tmp_path / "pca_run.mfrm"
    project.save_project(mfrm_path, an)
    reopen_info = _check_project_cli(mfrm_path, tmp_path / "reopen_pca")
    assert reopen_info["estimation_calls"] == 0
    restored = project.load_project(mfrm_path)
    assert restored["fit_settings"]["compute_pca"] is True
    assert dict(outputs.figure_choices(restored)).get("残差PCA · 固有値スクリー") == {"kind": "pca"}


def test_linear_logit_rescaling_affine_invariance_hold_protection_and_tamper_guard(tmp_path):
    df = _small_crossed_df()
    base = desktop.run_analysis(
        df,
        person_col="Person",
        facet_cols=["Rater", "Task"],
        score_col="Score",
        rating_min=1,
        rating_max=4,
        model="RSM",
        method="JMLE",
    )
    scaled = desktop.run_analysis(
        df,
        person_col="Person",
        facet_cols=["Rater", "Task"],
        score_col="Score",
        rating_min=1,
        rating_max=4,
        model="RSM",
        method="JMLE",
        umean=500.0,
        uscale=100.0,
    )
    assert "umean" not in base["fit_settings"]
    assert scaled["fit_settings"]["umean"] == 500.0
    assert scaled["fit_settings"]["uscale"] == 100.0

    # Canonical logits in Measures and Fit remain identical; Separation & Reliability are affine-invariant.
    base_tbls = {name: tbl for name, tbl, _ in outputs.current_tables(base)}
    scaled_tbls = {name: tbl for name, tbl, _ in outputs.current_tables(scaled)}
    pd.testing.assert_frame_equal(base_tbls["Measures and Fit"], scaled_tbls["Measures and Fit"])
    pd.testing.assert_frame_equal(
        base_tbls["Separation and Reliability by Facet"],
        scaled_tbls["Separation and Reliability by Facet"],
    )

    rescaled_tbl = scaled_tbls["Rescaled Measures (Umean / Uscale)"]
    assert np.allclose(
        rescaled_tbl["ScaledMeasure"].to_numpy(dtype=float),
        500.0 + 100.0 * rescaled_tbl["LogitEstimate"].to_numpy(dtype=float),
    )
    base_rater_se = pd.to_numeric(base["report_tables"]["Rater"]["ModelSE"], errors="coerce").to_numpy(dtype=float)
    scaled_rater_se = pd.to_numeric(scaled["report_tables"]["Rater"]["ModelSE"], errors="coerce").to_numpy(dtype=float)
    assert np.allclose(scaled_rater_se, 100.0 * base_rater_se)

    # Under MML, ScaledSE / ScaledLower / ScaledUpper must remain blank (`""`), never leaking withheld SEs!
    mml_scaled = desktop.run_analysis(
        df,
        person_col="Person",
        facet_cols=["Rater", "Task"],
        score_col="Score",
        rating_min=1,
        rating_max=4,
        model="RSM",
        method="MML",
        quad_points=15,
        population_prior_sd=1.0,
        mml_engine="EM",
        umean=50.0,
        uscale=10.0,
    )
    mml_tables = dict((n, t) for n, t, _ in outputs.current_tables(mml_scaled))
    mml_rescaled_tbl = mml_tables["Rescaled Measures (Umean / Uscale)"]
    assert mml_rescaled_tbl["ScaledSE"].eq("").all()
    assert mml_rescaled_tbl["ScaledLower"].eq("").all()
    assert mml_rescaled_tbl["ScaledUpper"].eq("").all()
    assert list(mml_rescaled_tbl["Basis"]) == list(mml_tables["Measures and Fit"]["Basis"])

    # Invalid umean / uscale rejected pre-fit.
    for bad_umean, bad_uscale in [(0.0, 0.0), (0.0, -5.0), (float("nan"), 1.0), (0.0, float("inf")), (True, 1.0)]:
        with pytest.raises(ValueError):
            basic.run_basic_analysis(
                df,
                person_col="Person",
                facet_cols=["Rater", "Task"],
                score_col="Score",
                rating_min=1,
                rating_max=4,
                umean=bad_umean,
                uscale=bad_uscale,
            )

    # Tamper guard in project.validate_analysis when Measure disagrees with Estimate * uscale + umean.
    tampered = deepcopy(scaled)
    tampered["report_tables"]["Rater"] = tampered["report_tables"]["Rater"].copy()
    tampered["report_tables"]["Rater"]["Measure"] = tampered["report_tables"]["Rater"]["Measure"] + 1.0
    tampered["content_fingerprint"] = project.content_fingerprint(tampered)
    with pytest.raises(project.ProjectError, match="disagrees with umean/uscale"):
        project.save_project(tmp_path / "tampered.mfrm", tampered)


def test_desktop_gui_rescaling_presets_pca_toggle_and_bilingual_roundtrip():
    tk = pytest.importorskip("tkinter")
    try:
        root = tk.Tk()
    except tk.TclError as exc:
        pytest.skip(f"Tk display unavailable: {exc}")
    root.withdraw()
    try:
        gui = DesktopApp(root)
        gui.set_data(_small_crossed_df(), "sample.csv")
        gui.minimum.set("1")
        gui.maximum.set("4")
        gui.umean_var.set("500")
        gui.uscale_var.set("100")
        gui.compute_pca_var.set(True)
        s = gui.settings()
        assert (s["umean"], s["uscale"], s["compute_pca"]) == (500.0, 100.0, True)

        gui.restore_settings({
            "kind": "mfrm-analysis-v2",
            "fit_settings": {
                "person_col": "Person",
                "facet_cols": ["Rater", "Task"],
                "score_col": "Score",
                "rating_min": 1,
                "rating_max": 4,
                "model": "RSM",
                "method": "JMLE",
                "umean": 50.0,
                "uscale": 10.0,
                "compute_pca": True,
            },
        })
        assert float(gui.umean_var.get()) == pytest.approx(50.0)
        assert float(gui.uscale_var.get()) == pytest.approx(10.0)
        assert gui.compute_pca_var.get() is True

        for en_str in (
            "Separation and Reliability by Facet",
            "Inter-rater Agreement and Krippendorff's Alpha",
            "Rescaled Measures (Umean / Uscale)",
            "Residual PCA: Eigenvalues and Stability",
            "Residual PCA: First-Contrast Loadings",
        ):
            ja_str = i18n.text(en_str, "ja")
            assert ja_str != en_str
            assert i18n.text(ja_str, "en") == en_str
    finally:
        root.destroy()


def test_bidirectional_mfrm_portability_web_to_desktop_and_desktop_to_web(tmp_path, monkeypatch):
    # 1. Web MFRM with custom column names, anchors, PCA, and linear rescaling -> .mfrm -> Desktop zero-refit reopen
    raw_df = _small_crossed_df().rename(columns={"Person": "Examinee", "Score": "EssayScore"})
    anchor_df = pd.DataFrame([{"Facet": "Rater", "Level": "R1", "Anchor": 0.25}])
    group_anchor_df = pd.DataFrame([
        {"Facet": "Task", "Level": "T1", "Group": "G1", "GroupValue": 0.0},
        {"Facet": "Task", "Level": "T2", "Group": "G1", "GroupValue": 0.0},
    ])
    web_res = engine.mfrm_estimate(
        raw_df,
        person_col="Examinee",
        facet_cols=["Rater", "Task"],
        score_col="EssayScore",
        model="RSM",
        method="JMLE",
        anchor_df=anchor_df,
        group_anchor_df=group_anchor_df,
    )
    web_diag = engine.mfrm_diagnostics(web_res, compute_pca=True)
    web_tbls = engine.calc_facets_report_tbls(web_res, web_diag, umean=50.0, uscale=10.0)
    web_sfile = engine.compute_scorefile(web_res)
    web_bias = {"Rater x Task": web_app.estimate_bias_interaction(web_res, web_diag, "Rater", "Task")}

    mfrm_bytes = portability.build_mfrm_project_from_web(
        web_res,
        web_diag,
        web_tbls,
        web_sfile,
        person_col="Examinee",
        score_col="EssayScore",
        all_bias_results=web_bias,
        umean=50.0,
        uscale=10.0,
    )
    web_mfrm_path = tmp_path / "from_web_jmle.mfrm"
    web_mfrm_path.write_bytes(mfrm_bytes)
    reopen_jmle = _check_project_cli(web_mfrm_path, tmp_path / "reopen_jmle")
    assert reopen_jmle["estimation_calls"] == 0

    # 2. Web fixed-N(0,1) MML ("Auto (recommended)") + structural intervals + joint interaction -> .mfrm -> Desktop
    df_mml = _small_crossed_df()
    mml_res = engine.mfrm_estimate(
        df_mml,
        person_col="Person",
        facet_cols=["Rater", "Task"],
        score_col="Score",
        model="PCM",
        step_facet="Task",
        method="MML",
        mml_engine="Auto (recommended)",
        quad_points=15,
        population_prior_sd=1.0,
    )
    mml_diag = engine.mfrm_diagnostics(mml_res, compute_pca=False)
    mml_res["_structural_intervals"] = intervals.compute_structural_intervals(
        mml_res, facet="Rater", difference=["R1", "R2"]
    )
    mml_res["_joint_interaction"] = joint.fit_joint_interaction(
        mml_res, facet_a="Rater", facet_b="Task", difference=[["R1", "T1"], ["R2", "T2"]]
    )
    mml_bytes = portability.build_mfrm_project_from_web(mml_res, mml_diag)
    mml_path = tmp_path / "from_web_mml_joint.mfrm"
    mml_path.write_bytes(mml_bytes)
    reopen_mml = _check_project_cli(mml_path, tmp_path / "reopen_mml")
    assert reopen_mml["estimation_calls"] == 0
    restored_mml = project.load_project(mml_path)
    assert "joint_analysis" in restored_mml and len(restored_mml["interval_analyses"]) == 1

    # 3. Web 3-facet conditioned G/D -> .mfrm -> Desktop zero-refit reopen
    gd_df = pd.DataFrame([
        {
            "Person": f"P{p:02d}",
            "Rater": f"R{r + 1}",
            "Task": f"T{t + 1}",
            "Criterion": f"C{c + 1}",
            "Score": 1 + ((2 * p + r + t + c) % 5),
        }
        for p in range(10)
        for r in range(3)
        for t in range(2)
        for c in range(2)
    ])
    gd_bundle = gd.compute_generalizability_from_data(
        gd_df,
        object_facet="Person",
        random_facets=["Rater", "Task"],
        score_col="Score",
        target_population="Examinees across trained raters and tasks",
        fixed_condition={"facet": "Criterion", "mode": "mean", "level": None},
    )
    assert gd_bundle["available"] is True
    gd_bytes = portability.build_gd_project_from_web(gd_bundle, input_data=gd_df)
    gd_path = tmp_path / "from_web_gd.mfrm"
    gd_path.write_bytes(gd_bytes)
    reopen_gd = _check_project_cli(gd_path, tmp_path / "reopen_gd")
    assert reopen_gd["estimation_calls"] == 0

    # 4. Research / Web-only options cleanly rejected by check_web_mfrm_portability
    gpcm_res = engine.mfrm_estimate(
        df_mml,
        person_col="Person",
        facet_cols=["Rater", "Task"],
        score_col="Score",
        model="GPCM",
        step_facet="Task",
        method="JMLE",
        maxit=40,
    )
    gpcm_chk = portability.check_web_mfrm_portability(gpcm_res)
    assert gpcm_chk["eligible"] is False and any("GPCM" in r for r in gpcm_chk["reasons"])

    free_sd_stub = deepcopy(mml_res)
    free_sd_stub["config"] = dict(free_sd_stub["config"], estimate_population_sd=True)
    free_chk = portability.check_web_mfrm_portability(free_sd_stub)
    assert free_chk["eligible"] is False and any("Free population SD" in r for r in free_chk["reasons"])

    # 4b. Web UI export buttons (`dl_mfrm_project_archive` and `standalone_gd_download_mfrm`)
    def view_exports(eligible: bool):
        import pandas as pd
        import streamlit as st
        import streamlit_app as app

        st.session_state["lang"] = "ja"
        st.session_state["downloads_panel"] = "scripts_config"
        res_obj = st.session_state["_test_res"]
        diag_obj = st.session_state["_test_diag"]
        app._render_downloads(
            res_obj,
            diag_obj,
            {},
            pd.DataFrame(),
            pd.DataFrame(),
            generate_figures=False,
        )
        if eligible:
            app.render_g_d_study_bundle(
                st.session_state["_test_gd_bundle"],
                key_prefix="standalone_gd",
                raw_input_data=st.session_state["_test_gd_df"],
            )

    at_exp = AppTest.from_function(view_exports, args=(True,))
    at_exp.session_state["_test_res"] = mml_res
    at_exp.session_state["_test_diag"] = mml_diag
    at_exp.session_state["_test_gd_bundle"] = gd_bundle
    at_exp.session_state["_test_gd_df"] = gd_df
    at_exp.run(timeout=45)
    assert not at_exp.exception, [e.value for e in at_exp.exception]
    dl_keys = {b.key for b in at_exp.get("download_button")}
    assert "dl_mfrm_project_archive" in dl_keys
    assert "standalone_gd_download_mfrm" in dl_keys

    at_inelig = AppTest.from_function(view_exports, args=(False,))
    at_inelig.session_state["_test_res"] = gpcm_res
    at_inelig.session_state["_test_diag"] = engine.mfrm_diagnostics(gpcm_res, compute_pca=False)
    at_inelig.run(timeout=45)
    assert not at_inelig.exception, [e.value for e in at_inelig.exception]
    assert "dl_mfrm_project_archive" not in {b.key for b in at_inelig.get("download_button")}
    assert at_inelig.info

    # 5. Desktop .mfrm -> Web zero-estimation viewer (AppTest in EN & JA)
    def view_project(payload_bytes: bytes, lang: str):
        import streamlit as st
        import streamlit_app as app
        st.session_state["lang"] = lang
        st.session_state["_mfrm_project_loaded_bytes"] = payload_bytes
        st.session_state["_mfrm_project_loaded_name"] = "desktop_verified.mfrm"
        app.render_mfrm_project_workspace()

    def _forbid(*_a, **_kw):
        raise AssertionError("Reopening .mfrm in Web must never run any estimator")

    for name in ("mfrm_estimate", "mfrm_diagnostics", "compute_generalizability_from_data"):
        monkeypatch.setattr(web_app, name, _forbid)
        if hasattr(engine, name):
            monkeypatch.setattr(engine, name, _forbid)

    for payload in (mml_bytes, gd_bytes):
        for lang in ("en", "ja"):
            at = AppTest.from_function(view_project, args=(payload, lang)).run(timeout=60)
            assert not at.exception, [e.value for e in at.exception]
            assert at.success and at.dataframe
            assert any(b.key == "mfrm_project_download_html" for b in at.get("download_button"))

    # Tampered bytes show user-facing error without uncaught exception.
    at_bad = AppTest.from_function(view_project, args=(b"not-a-valid-mfrm-archive", "ja")).run(timeout=30)
    assert not at_bad.exception and at_bad.error


@pytest.mark.parametrize("lang", ["en", "ja"])
def test_web_joint_interaction_ui_matches_desktop_publication_table_and_restores_without_refit(lang):
    def view_joint(language: str):
        import pandas as pd
        import streamlit as st
        import streamlit_app as app
        from mfrm_app import engine as eng

        if "fit" not in st.session_state:
            df = pd.DataFrame([
                {"Person": f"P{p:02d}", "Rater": f"R{r}", "Task": f"T{t}", "Score": (p + 2 * r + t) % 3}
                for p in range(24)
                for r in range(3)
                for t in range(2)
            ])
            st.session_state["fit"] = eng.mfrm_estimate(
                df,
                "Person",
                ["Rater", "Task"],
                "Score",
                method="MML",
                model="RSM",
                quad_points=31,
                population_prior_sd=1.0,
                maxit=120,
            )
        st.session_state["lang"] = language
        app._render_joint_interaction(st.session_state["fit"], key_prefix="test", expanded=True)

    at = AppTest.from_function(view_joint, args=(lang,)).run(timeout=60)
    assert not at.exception, [e.value for e in at.exception]
    prefix = "test_joint_interaction"
    at.checkbox(key=prefix + "_difference").check().run(timeout=45)
    at.button(key=prefix + "_calculate").click().run(timeout=60)
    assert not at.exception, [e.value for e in at.exception]

    fit = at.session_state["fit"]
    record = fit["_joint_interaction"]
    expected_pub, _ = joint.publication_table(fit, record)
    assert len(at.dataframe) >= 1
    assert len(at.get("download_button")) == 2
    assert record["settings"]["difference"] is not None

    # Language switch and re-render must NOT call fit_joint_interaction again.
    with patch.object(joint, "fit_joint_interaction", side_effect=AssertionError("Must not refit on rerun")):
        at.session_state["lang"] = "ja" if lang == "en" else "en"
        at.run(timeout=45)
        assert not at.exception
        assert len(at.get("download_button")) == 2

    # Web publication report tables include compact qualified Joint interaction table.
    pub_tables = web_app._publication_report_tables(fit, {"measures": pd.DataFrame()}, {}, include_details=False)
    ji_title, ji_tbl, _ = next(item for item in pub_tables if item[0].startswith("Joint interaction:"))
    assert "Rater × Task" in ji_title
    assert list(ji_tbl["Target"]) == list(expected_pub["Target"])

    # Score mutation invalidates saved joint interaction without automatic refit.
    fit["prep"]["data"].iloc[0, fit["prep"]["data"].columns.get_loc("Score")] = 2
    with patch.object(joint, "fit_joint_interaction", side_effect=AssertionError("Must not refit on invalidation")):
        at.run(timeout=45)
    assert not at.exception and at.warning
    assert not at.get("download_button")


def test_web_research_options_are_visually_separated_in_collapsed_sidebar_expander():
    from tests.test_analysis_journeys import _prepare

    at = _prepare("free_sd")
    lang = at.session_state["lang"]
    expected_label = web_app._load_locale(lang)["sidebar_estimation"]["research_options_expander"]
    all_expanders = list(at.expander) + list(at.sidebar.expander)
    research_exp = next((e for e in all_expanders if e.label == expected_label), None)
    assert research_exp is not None
    assert not research_exp.proto.expanded
    assert at.checkbox(key="facets_mode_estimate_population_sd").value is True
