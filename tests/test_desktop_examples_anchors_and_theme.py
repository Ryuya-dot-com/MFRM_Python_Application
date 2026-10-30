"""Desktop suite tests: 8-scenario gallery, Wide->Long conversion, live connectivity &
identifiability diagnosis, fixed/group anchor calibration & .mfrm validation, and Desktop UI theme."""
from __future__ import annotations

from copy import deepcopy
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from desktop_app import DesktopApp, PALETTE
from mfrm_app import basic_analysis as basic
from mfrm_app import desktop_analysis as desktop
from mfrm_app import desktop_examples as examples
from mfrm_app import desktop_outputs as outputs
from mfrm_app import engine, io_tables, project
from mfrm_app.desktop_i18n import text


def test_eight_curated_scenarios_match_bundled_csvs_and_translations():
    keys = [item["key"] for item in examples.SCENARIO_CATALOG]
    assert keys == [
        "writing_3facet",
        "clinical_4facet",
        "anchor_equating",
        "missing_ratings",
        "sparse_spiral",
        "sparse_disconnected",
        "wide_spreadsheet",
        "gd_2facet",
    ]

    generators = {
        "clinical_4facet": examples.generate_clinical_4facet_df,
        "anchor_equating": examples.generate_anchor_equating_df,
        "missing_ratings": examples.generate_missing_ratings_df,
        "sparse_spiral": examples.generate_sparse_spiral_df,
        "sparse_disconnected": examples.generate_disconnected_shifts_df,
        "wide_spreadsheet": examples.generate_wide_spreadsheet_df,
    }
    for item in examples.SCENARIO_CATALOG:
        csv_path = examples.EXAMPLES_DIR / item["filename"]
        assert csv_path.is_file(), f"Missing bundled CSV: {csv_path}"
        loaded = examples.load_scenario_dataframe(item["key"])
        assert len(loaded) > 0
        assert "ファセット" not in item["badge_ja"]
        assert "ファセット" not in item["title_ja"]
        assert "ファセット" not in item["summary_ja"]
        # Title roundtrips cleanly through desktop_i18n
        assert text(item["title_ja"], "en") == item["title_en"]
        assert text(item["title_en"], "ja") == item["title_ja"]
        if item["key"] in generators:
            expected = generators[item["key"]]()
            pd.testing.assert_frame_equal(loaded, expected)


def test_wide_spreadsheet_detection_and_conversion_preserves_ids_and_omits_blanks():
    wide_df = examples.load_scenario_dataframe("wide_spreadsheet")
    assert wide_df.shape == (40, 6)
    det = io_tables.detect_wide_format_columns(wide_df)
    assert det["looks_wide"] is True
    assert det["probable_id_cols"] == ["Person", "Task"]
    assert det["probable_score_cols"] == ["R1", "R2", "R3", "R4"]

    melted, dropped = DesktopApp._convert_wide_dataframe(
        wide_df,
        id_cols=["Person", "Task"],
        score_cols=["R1", "R2", "R3", "R4"],
        new_facet_name="Rater",
        score_col_name="Score",
    )
    assert dropped == 40
    assert melted.shape == (120, 4)
    assert list(melted.columns) == ["Person", "Task", "Rater", "Score"]
    assert all(s.isdigit() for s in melted["Score"])

    # Leading-zero ID preservation check
    custom_wide = pd.DataFrame(
        [("001", "3", "", "4"), ("002", "", "2", "5")],
        columns=["Examinee", "Item1", "Item2", "Item3"],
    )
    m2, d2 = DesktopApp._convert_wide_dataframe(
        custom_wide,
        id_cols=["Examinee"],
        score_cols=["Item1", "Item2", "Item3"],
        new_facet_name="Item",
        score_col_name="Score",
    )
    assert d2 == 2
    assert set(m2["Examinee"]) == {"001", "002"}
    assert m2["Score"].tolist() == ["3", "2", "4", "5"]


def test_connectivity_and_structural_nullity_across_dense_sparse_and_disconnected_designs():
    # 1. 4-facet clinical OSCE (100% density)
    osce = examples.load_scenario_dataframe("clinical_4facet")
    p_osce = desktop.preview_mfrm_connectivity(
        osce, person_col="Examinee", score_col="Rating",
        facet_cols=["Rater", "Station", "Competency", "Session"],
    )
    assert p_osce["ok"] is True
    assert p_osce["status"] == "ok"
    assert p_osce["nullity"] == 0
    assert p_osce["density_pct"] == pytest.approx(100.0)

    # 2. Missing ratings (900 rows, 12 blank Score cells -> 888 valid, 12 excluded)
    miss = examples.load_scenario_dataframe("missing_ratings")
    p_miss = desktop.preview_mfrm_connectivity(
        miss, person_col="Person", score_col="Score",
        facet_cols=["Rater", "Task", "Criterion"],
    )
    assert p_miss["ok"] is True
    assert p_miss["included"] == 888
    assert p_miss["excluded"] == 12
    assert p_miss["nullity"] == 0

    # 3. Sparse spiral design (25.0% cell density, connected via overlapping pairs)
    spiral = examples.load_scenario_dataframe("sparse_spiral")
    p_spiral = desktop.preview_mfrm_connectivity(
        spiral, person_col="Person", score_col="Score",
        facet_cols=["Rater", "Task", "Criterion"],
    )
    assert p_spiral["ok"] is True
    assert p_spiral["status"] == "ok"
    assert p_spiral["nullity"] == 0
    assert p_spiral["density_pct"] == pytest.approx(25.0)
    assert "スパース" in p_spiral["summary_ja"]
    assert "sparse" in p_spiral["summary_en"]

    # 4. Disconnected AM/PM rater shifts sharing Task1 & Task2 (Person x Rater has 2 components -> nullity = 1)
    disc = examples.load_scenario_dataframe("sparse_disconnected")
    p_disc = desktop.preview_mfrm_connectivity(
        disc, person_col="Person", score_col="Score",
        facet_cols=["Rater", "Task"],
    )
    assert p_disc["ok"] is False
    assert p_disc["status"] == "disconnected"
    assert p_disc["nullity"] == 1
    assert p_disc["component_counts"]["Rater"] == 2
    assert p_disc["component_counts"]["Task"] == 1

    # Bridging AM and PM shifts by anchoring one rater in each shift resolves nullity to 0!
    bridge_anchors = [
        {"Facet": "Rater", "Level": "R1", "Anchor": -0.35},
        {"Facet": "Rater", "Level": "R3", "Anchor": -0.20},
    ]
    p_bridged = desktop.preview_mfrm_connectivity(
        disc, person_col="Person", score_col="Score",
        facet_cols=["Rater", "Task"],
        anchors=bridge_anchors,
    )
    assert p_bridged["ok"] is True
    assert p_bridged["status"] == "anchored_ok"
    assert p_bridged["nullity"] == 0


@pytest.mark.parametrize("method", ["JMLE", "MML"])
def test_fixed_and_group_anchor_estimation_archival_and_tamper_rejection(tmp_path, monkeypatch, method):
    df = examples.load_scenario_dataframe("anchor_equating")
    presets = examples.ANCHOR_EQUATING_PRESETS
    extra = dict(quad_points=15, population_prior_sd=1.0) if method == "MML" else {}

    analysis = desktop.run_analysis(
        df,
        person_col="Person",
        score_col="Score",
        facet_cols=["Rater", "Task", "Criterion"],
        rating_min=0,
        rating_max=4,
        model="RSM",
        method=method,
        anchors=presets["anchors"],
        group_anchors=presets["group_anchors"],
        **extra,
    )
    assert analysis["kind"] == "mfrm-analysis-v2"
    measures = analysis["diagnostics"]["measures"]
    facet_specs = analysis["result"]["config"]["facet_specs"]
    facet_levels = analysis["result"]["config"]["facet_levels"]

    # Verify fixed anchors (Status A) hold exact bank values
    for rec in presets["anchors"]:
        row = measures[(measures["Facet"] == rec["Facet"]) & (measures["Level"] == rec["Level"])].iloc[0]
        idx = list(facet_levels[rec["Facet"]]).index(rec["Level"])
        assert float(facet_specs[rec["Facet"]]["anchors"][idx]) == pytest.approx(rec["Anchor"], abs=1e-12)
        assert float(row["Estimate"]) == pytest.approx(rec["Anchor"], abs=1e-10)

    # Verify group anchors (Status G) hold exact group mean 0.0 while allowing free relative spread
    crit_rows = measures[measures["Facet"] == "Criterion"]
    assert list(facet_specs["Criterion"]["groups"]) == ["Rubric", "Rubric", "Rubric"]
    assert float(facet_specs["Criterion"]["group_values"]["Rubric"]) == pytest.approx(0.0, abs=1e-12)
    assert float(crit_rows["Estimate"].mean()) == pytest.approx(0.0, abs=1e-10)
    assert float(crit_rows["Estimate"].std()) > 0.02

    # Connectivity & Scale Linking table is retained and rendered in both JA and EN
    tables_dict = {t[0]: t[1] for t in outputs.current_tables(analysis)}
    assert "Connectivity and Scale Linking" in tables_dict
    link_df = tables_dict["Connectivity and Scale Linking"]
    assert set(link_df["Category"]) == {"Connectivity", "FixedAnchor", "GroupAnchor"}
    html_en = outputs.report_html(analysis, language="en")
    html_ja = outputs.report_html(analysis, language="ja")
    assert "Connectivity and Scale Linking" in html_en
    assert "連結性・アンカー等化サマリ" in html_ja

    # Save and reopen without any estimation calls
    mfrm_path = tmp_path / f"anchored_{method}.mfrm"
    project.save_project(mfrm_path, analysis)

    monkeypatch.setattr(engine, "mfrm_estimate", lambda *a, **k: pytest.fail("Must not refit on reopen"))
    monkeypatch.setattr(engine, "mfrm_diagnostics", lambda *a, **k: pytest.fail("Must not recompute diagnostics"))
    reopened = project.load_project(mfrm_path)
    assert reopened["analysis_id"] == analysis["analysis_id"]
    assert outputs.report_html(reopened, language="en") == html_en

    # Tampering with stored fit_settings["anchors"] must be rejected by validate_analysis
    tampered = deepcopy(reopened)
    tampered["fit_settings"]["anchors"][0]["Anchor"] = 0.99
    tampered["content_fingerprint"] = project.content_fingerprint(tampered)
    with pytest.raises(project.ProjectError, match="anchor"):
        project.validate_analysis(tampered)


def test_desktop_app_ui_gallery_wide_converter_anchor_tab_and_color_banners(tmp_path):
    tk = pytest.importorskip("tkinter")
    try:
        root = tk.Tk()
    except tk.TclError:
        pytest.skip("Tk display unavailable in this environment")
    root.withdraw()
    try:
        app = DesktopApp(root, language="en")
        assert app.step_buttons["data"].cget("style") == "Current.TButton"

        # 1. Sample gallery opens all 8 scenarios
        gallery = app.open_sample_gallery()
        assert len(gallery._scenario_buttons) == 8
        gallery.destroy()

        # 2. Load 4-facet clinical OSCE scenario -> emerald 'ok' preflight banner
        app.load_example_scenario("clinical_4facet")
        assert app.person.get() == "Examinee"
        assert app.score.get() == "Rating"
        assert [app.facets.get(i) for i in app.facets.curselection()] == [
            "Rater", "Station", "Competency", "Session"
        ]
        assert app.mfrm_preflight_banner._stripe.cget("bg") == PALETTE["banner"]["ok"][0]

        # 3. Load disconnected shifts -> amber 'warn' banner -> import anchor CSV -> blue 'info' banner
        app.load_example_scenario("sparse_disconnected")
        assert app.mfrm_preflight_banner._stripe.cget("bg") == PALETTE["banner"]["warn"][0]
        anchor_csv = tmp_path / "bridge_anchors.csv"
        anchor_csv.write_text("Facet,Level,Anchor,Status\nRater,R1,-0.35,A\nRater,R3,-0.20,A\nRater,R2,0.25,F\n", encoding="utf-8")
        app.import_anchor_csv(anchor_csv)
        assert len(app.anchor_records) == 2  # Status F row skipped, two Status A rows imported
        assert app.mfrm_preflight_banner._stripe.cget("bg") == PALETTE["banner"]["info"][0]

        # 4. Load wide_spreadsheet -> Wide->Long converter opens and converts 40x6 -> 120x4
        app.load_example_scenario("wide_spreadsheet")
        # Find the open Wide->Long Toplevel dialog
        dialogs = [w for w in root.winfo_children() if isinstance(w, tk.Toplevel)]
        assert len(dialogs) == 1
        wide_win = dialogs[0]
        wide_win._convert_button.invoke()
        assert app.data.shape == (120, 4)
        assert app.minimum.get() == "0" and app.maximum.get() == "5"

        # 5. Load anchor_equating, switch EN -> JA -> EN, and check settings fingerprint parity
        app.load_example_scenario("anchor_equating")
        assert len(app.anchor_records) == 4
        assert len(app.group_anchor_records) == 3
        app.language_choice.set("日本語")
        app.change_language()
        assert "アンカー・連結確認OK" in app.mfrm_preflight.get()
        app.language_choice.set("English")
        app.change_language()
        assert "Anchor & connectivity OK" in app.mfrm_preflight.get()

        req = app.settings()
        analysis = desktop.run_analysis(app.data, **{k: v for k, v in req.items() if k != "workflow"}, workflow="mfrm")
        app.restore_settings(analysis)
        app.show_analysis(analysis)
        assert app.draft_notice.get() == ""
    finally:
        root.destroy()


def test_palette_meets_wcag_aa_contrast_ratios():
    def _rel_lum(hex_col: str) -> float:
        h = hex_col.lstrip("#")
        rgb = [int(h[i:i + 2], 16) / 255.0 for i in (0, 2, 4)]
        lin = [c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4 for c in rgb]
        return 0.2126 * lin[0] + 0.7152 * lin[1] + 0.0722 * lin[2]

    def _contrast(c1: str, c2: str) -> float:
        l1, l2 = _rel_lum(c1), _rel_lum(c2)
        return (max(l1, l2) + 0.05) / (min(l1, l2) + 0.05)

    for name, fg, bg in [
        ("header_fg/header_bg", PALETTE["header_fg"], PALETTE["header_bg"]),
        ("header_sub/header_bg", PALETTE["header_sub"], PALETTE["header_bg"]),
        ("primary_btn", "#FFFFFF", PALETTE["primary"]),
        ("accent_btn", "#FFFFFF", PALETTE["accent"]),
        ("danger_btn", PALETTE["danger_fg"], PALETTE["danger_bg"]),
        ("step_ready", PALETTE["step_ready_fg"], PALETTE["step_ready_bg"]),
        ("text/card_bg", PALETTE["text"], PALETTE["card_bg"]),
        ("muted/card_bg", PALETTE["muted"], PALETTE["card_bg"]),
        ("table_head", PALETTE["table_head_fg"], PALETTE["table_head_bg"]),
        ("banner_ok", PALETTE["banner"]["ok"][2], PALETTE["banner"]["ok"][1]),
        ("banner_info", PALETTE["banner"]["info"][2], PALETTE["banner"]["info"][1]),
        ("banner_warn", PALETTE["banner"]["warn"][2], PALETTE["banner"]["warn"][1]),
        ("banner_error", PALETTE["banner"]["error"][2], PALETTE["banner"]["error"][1]),
    ]:
        assert _contrast(fg, bg) >= 4.5, f"{name} fell below WCAG AA 4.5:1"

