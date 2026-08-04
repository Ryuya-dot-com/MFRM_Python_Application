"""Custom synthetic-data generator contracts."""

from __future__ import annotations

import inspect

import numpy as np
import pandas as pd
import pytest

import streamlit_app as app


def test_custom_simulation_shape_scores_and_determinism():
    kwargs = dict(
        n_person=12,
        facet_names=("Rater", "Task", "Criterion", "Occasion"),
        facet_level_counts=(4, 2, 3, 2),
        facet_sds=(0.35, 0.25, 0.20, 0.10),
        first_facet_levels_per_person=2,
        n_categories=5,
        thresholds=(-1.5, -0.5, 0.5, 1.5),
        seed=123,
    )
    df1 = app.generate_custom_mfrm_simulation_data(**kwargs)
    df2 = app.generate_custom_mfrm_simulation_data(**kwargs)

    assert isinstance(df1, pd.DataFrame)
    assert df1.equals(df2)
    assert list(df1.columns) == ["Person", "Rater", "Task", "Criterion", "Occasion", "Score"]
    assert df1.shape == (12 * 2 * 2 * 3 * 2, 6)
    assert df1["Score"].between(0, 4).all()
    assert (df1.groupby("Person")["Rater"].nunique() == 2).all()


def test_custom_simulation_missing_rate_is_reproducible():
    full_rows = 40 * 3 * 2 * 2
    df = app.generate_custom_mfrm_simulation_data(
        n_person=40,
        n_rater=4,
        n_task=2,
        n_criterion=2,
        raters_per_person=3,
        n_categories=4,
        missing_rate=0.20,
        seed=99,
    )
    df_again = app.generate_custom_mfrm_simulation_data(
        n_person=40,
        n_rater=4,
        n_task=2,
        n_criterion=2,
        raters_per_person=3,
        n_categories=4,
        missing_rate=0.20,
        seed=99,
    )

    assert df.equals(df_again)
    assert 0 < len(df) < full_rows
    assert abs(len(df) - int(full_rows * 0.80)) <= 40


def test_custom_simulation_keeps_legacy_three_facet_signature():
    df = app.generate_custom_mfrm_simulation_data(
        n_person=10,
        n_rater=3,
        n_task=2,
        n_criterion=2,
        raters_per_person=2,
        n_categories=4,
        seed=77,
    )
    assert list(df.columns) == ["Person", "Rater", "Task", "Criterion", "Score"]
    assert df.shape == (10 * 2 * 2 * 2, 5)


def test_custom_simulation_can_force_zero_count_category():
    df = app.generate_custom_mfrm_simulation_data(
        n_person=30,
        facet_names=("Rater", "Task", "Criterion"),
        facet_level_counts=(4, 2, 3),
        first_facet_levels_per_person=2,
        n_categories=5,
        zero_count_score=3,
        seed=888,
    )
    assert df["Score"].between(0, 4).all()
    assert 3 not in set(df["Score"].unique())


@pytest.mark.parametrize("zero_score", [0, 4])
def test_custom_simulation_can_force_edge_zero_count_category(zero_score):
    df = app.generate_custom_mfrm_simulation_data(
        n_person=30,
        facet_names=("Rater", "Task", "Criterion"),
        facet_level_counts=(4, 2, 3),
        first_facet_levels_per_person=2,
        n_categories=5,
        zero_count_score=zero_score,
        seed=889 + zero_score,
    )
    assert df["Score"].between(0, 4).all()
    assert zero_score not in set(df["Score"].unique())


def test_custom_simulation_bundle_includes_preview_tables():
    bundle = app.generate_custom_mfrm_simulation_bundle(
        n_person=12,
        facet_names=("Judge", "Prompt"),
        facet_level_counts=(3, 2),
        first_facet_levels_per_person=2,
        n_categories=5,
        seed=99,
    )
    assert set(bundle) >= {
        "data",
        "person_truth",
        "facet_truth",
        "threshold_truth",
        "category_counts",
        "sparse_design_audit",
        "sparse_pair_cells",
        "meta",
    }
    assert bundle["meta"]["facet_names"] == ["Judge", "Prompt"]
    assert bundle["category_counts"]["Count"].sum() == len(bundle["data"])
    assert not bundle["sparse_design_audit"].empty
    assert not bundle["sparse_pair_cells"].empty
    assert "thresholds" in bundle["meta"]
    assert "seed" in bundle["meta"]


def test_custom_simulation_sparse_design_audit_flags_sparse_cases():
    bundle = app.generate_custom_mfrm_simulation_bundle(
        n_person=18,
        facet_names=("Rater", "Task", "Criterion"),
        facet_level_counts=(5, 3, 3),
        first_facet_levels_per_person=1,
        n_categories=5,
        zero_count_score=4,
        missing_rate=0.35,
        seed=20260607,
    )

    audit = bundle["sparse_design_audit"]
    assert {
        "Priority",
        "Check",
        "Status",
        "Value",
        "Threshold",
        "Evidence",
        "Action",
        "DownloadFile",
    } <= set(audit.columns)
    assert not audit.empty
    assert set(audit["Check"]) >= {
        "Score category support",
        "Generated missingness",
        "Person row support",
        "Facet-pair cell density",
    }
    assert audit["Status"].isin({"Review", "Hold before reporting"}).any()

    pair_cells = bundle["sparse_pair_cells"]
    assert {
        "FacetPair",
        "Status",
        "PossibleCells",
        "ObservedCells",
        "ZeroCells",
        "LowCountCells",
        "MinCellCount",
        "CellDensity",
        "Action",
    } <= set(pair_cells.columns)
    assert not pair_cells.empty
    assert pair_cells["LowCountCells"].ge(0).all()
    assert pair_cells["ZeroCells"].ge(0).all()


def test_custom_simulation_sparse_audit_is_reproducible():
    kwargs = dict(
        n_person=24,
        facet_names=("Rater", "Task", "Criterion"),
        facet_level_counts=(4, 2, 3),
        first_facet_levels_per_person=2,
        n_categories=5,
        missing_rate=0.15,
        seed=4242,
    )
    bundle_1 = app.generate_custom_mfrm_simulation_bundle(**kwargs)
    bundle_2 = app.generate_custom_mfrm_simulation_bundle(**kwargs)

    pd.testing.assert_frame_equal(
        bundle_1["sparse_design_audit"],
        bundle_2["sparse_design_audit"],
    )
    pd.testing.assert_frame_equal(
        bundle_1["sparse_pair_cells"],
        bundle_2["sparse_pair_cells"],
    )


def test_custom_simulation_preview_exposes_sparse_design_downloads():
    source = inspect.getsource(app.render_custom_simulation_preview_panel)
    assert "sim_preview_tab_sparse" in source
    assert "custom_simulation_sparse_design_audit.csv" in source
    assert "custom_simulation_sparse_pair_cells.csv" in source


def test_custom_simulation_design_presets_cover_common_stress_cases():
    balanced = app.get_custom_simulation_design_preset_settings("balanced")
    assert balanced["n_person"] == 60
    assert balanced["facet_names"] == ["Rater", "Task", "Criterion"]
    assert balanced["facet_level_counts"] == [2, 3, 3]
    assert balanced["first_facet_levels_per_person"] == 2
    assert balanced["zero_count_score"] is None

    sparse = app.get_custom_simulation_design_preset_settings("sparse_coverage")
    assert sparse["facet_level_counts"][0] > balanced["facet_level_counts"][0]
    assert sparse["first_facet_levels_per_person"] == 1
    assert sparse["missing_rate"] > 0

    missing = app.get_custom_simulation_design_preset_settings("planned_missingness")
    assert missing["missing_rate"] >= 0.20

    zero_category = app.get_custom_simulation_design_preset_settings("zero_category")
    assert zero_category["zero_count_score"] == zero_category["n_categories"] - 1

    fallback = app.get_custom_simulation_design_preset_settings("unknown")
    assert fallback == balanced


def test_custom_simulation_sparse_reporting_context_feeds_apa_and_binder():
    bundle = app.generate_custom_mfrm_simulation_bundle(
        n_person=18,
        facet_names=("Rater", "Task", "Criterion"),
        facet_level_counts=(5, 3, 3),
        first_facet_levels_per_person=1,
        n_categories=5,
        zero_count_score=4,
        missing_rate=0.35,
        seed=20260607,
    )
    context = app.build_custom_simulation_sparse_reporting_context(bundle)
    settings = app.build_custom_simulation_settings_table(
        bundle,
        ui_meta={"design_preset": "sparse_coverage", "threshold_mode": "even"},
    )

    assert {
        "SimulationCheck",
        "AuditStatus",
        "ReportStatus",
        "ActionBeforeReporting",
        "CaveatToCarry",
        "EvidenceFiles",
        "DownloadFile",
    } <= set(context.columns)
    assert "Do not report yet" in set(context["ReportStatus"])
    assert context["EvidenceFiles"].astype(str).str.contains("custom_simulation_sparse_design_audit.csv").any()
    assert {
        "Section",
        "Setting",
        "Value",
        "Source",
        "ReportNote",
        "EvidenceFile",
    } <= set(settings.columns)
    assert settings["Setting"].astype(str).str.contains("Seed").any()
    assert settings["Setting"].astype(str).str.contains("Adjacent thresholds").any()
    assert settings["Setting"].astype(str).str.contains("Requested missing rate").any()
    assert settings["Value"].astype(str).str.contains("sparse_coverage").any()

    apa_draft = app.generate_report_ready_apa_results_draft(
        {},
        {},
        simulation_sparse_context=context,
        simulation_settings=settings,
    )
    assert "## Simulation Design Use Conditions" in apa_draft
    assert "## Simulation Settings To Archive" in apa_draft
    assert "custom_simulation_settings.csv" in apa_draft
    assert "custom_simulation_sparse_reporting_context.csv" in apa_draft
    assert "Do not report yet" in apa_draft

    binder_assets = app.build_manuscript_binder_assets(
        {
            "custom_simulation_settings": settings,
            "custom_simulation_sparse_reporting_context": context,
            "custom_simulation_sparse_design_audit": bundle["sparse_design_audit"],
            "custom_simulation_sparse_pair_cells": bundle["sparse_pair_cells"],
        },
        {"apa_results_paragraph_draft.md": apa_draft},
        public_export_mode=True,
    )
    assert "custom_simulation_settings.csv" in binder_assets
    assert "custom_simulation_sparse_reporting_context.csv" in binder_assets
    assert "custom_simulation_sparse_design_audit.csv" in binder_assets
    assert "custom_simulation_sparse_pair_cells.csv" in binder_assets
    assert "custom_simulation_settings.csv" in binder_assets["README_first.md"]
    assert "custom_simulation_sparse_reporting_context.csv" in binder_assets["README_first.md"]

    reanalysis_checklist = app.build_report_ready_reanalysis_checklist(
        {},
        {},
        simulation_sparse_context=context,
        simulation_settings=settings,
    )
    assert {
        "Simulation settings archive",
        "Simulation design reportability",
    }.issubset(set(reanalysis_checklist["Phase"].astype(str)))
    joined_reanalysis = " ".join(reanalysis_checklist.astype(str).to_numpy().ravel().tolist())
    assert "custom_simulation_settings.csv" in joined_reanalysis
    assert "custom_simulation_sparse_reporting_context.csv" in joined_reanalysis


def test_custom_simulation_source_exposes_design_presets():
    source = inspect.getsource(app.render_custom_simulation_source)
    assert "sim_design_preset_label" in source
    assert "sim_design_preset_apply_button" in source
    assert "CUSTOM_SIMULATION_DESIGN_PRESET_KEYS" in source


def test_custom_simulation_runs_existing_analysis_functions():
    bundle = app.generate_custom_mfrm_simulation_bundle(
        n_person=14,
        facet_names=("Rater", "Task", "Criterion"),
        facet_level_counts=(3, 2, 2),
        facet_sds=(0.55, 0.35, 0.45),
        first_facet_levels_per_person=2,
        n_categories=5,
        thresholds=(-1.2, -0.4, 0.4, 1.2),
        zero_count_score=3,
        missing_rate=0.02,
        seed=20260515,
    )
    df = bundle["data"]
    facet_cols = ["Rater", "Task", "Criterion"]

    readiness = app.build_readiness_report(
        data=df,
        person_col="Person",
        score_col="Score",
        facet_cols=facet_cols,
    )
    assert readiness["overall"] in {"ok", "warning"}
    assert readiness["n_issues"] == 0

    result = app.mfrm_estimate(
        df,
        "Person",
        facet_cols,
        "Score",
        rating_min=0,
        rating_max=4,
        keep_original=True,
        model="RSM",
        method="JMLE",
        maxit=60,
        reltol=1e-4,
    )
    assert bool(result["summary"]["Converged"].iloc[0])
    assert int(result["summary"]["Categories"].iloc[0]) == 5

    diagnostics = app.mfrm_diagnostics(
        result,
        compute_pca=True,
        compute_marginal=False,
    )
    assert len(diagnostics["obs"]) == len(df)
    assert not diagnostics["measures"].empty
    assert not diagnostics["reliability"].empty


def test_wright_preview_uses_adjacent_threshold_locations():
    thresholds = (-1.5, -0.5, 0.5, 1.5)
    bundle = app.generate_custom_mfrm_simulation_bundle(
        n_person=12,
        facet_names=("Judge", "Prompt"),
        facet_level_counts=(3, 2),
        first_facet_levels_per_person=2,
        n_categories=5,
        thresholds=thresholds,
        seed=100,
    )
    fig = app.build_custom_simulation_wright_preview_figure(bundle)
    x_lines = [float(shape.x0) for shape in fig.layout.shapes]
    np.testing.assert_allclose(x_lines, np.array(thresholds, dtype=float))


def test_custom_threshold_parser_accepts_common_delimiters():
    parsed = app.parse_custom_simulation_thresholds("-1.2, 0; 1.2", 4)
    np.testing.assert_allclose(parsed, np.array([-1.2, 0.0, 1.2]))


def test_custom_threshold_parser_rejects_wrong_length():
    with pytest.raises(ValueError, match="Expected 3 thresholds"):
        app.parse_custom_simulation_thresholds("-1, 0", 4)


def test_custom_threshold_parser_rejects_non_numeric_and_non_finite():
    with pytest.raises(ValueError, match="numeric"):
        app.parse_custom_simulation_thresholds("-1, nope, 1", 4)
    with pytest.raises(ValueError, match="finite"):
        app.parse_custom_simulation_thresholds("-1, inf, 1", 4)


def test_default_custom_thresholds_binary_are_zero():
    np.testing.assert_allclose(
        app.default_custom_simulation_thresholds(2, step_span=3.0),
        np.array([0.0]),
    )


def test_custom_simulation_rejects_runaway_designs():
    with pytest.raises(ValueError, match="100,000 rows"):
        app.generate_custom_mfrm_simulation_data(
            n_person=300,
            n_rater=20,
            n_task=10,
            n_criterion=10,
            raters_per_person=20,
            n_categories=5,
        )


def test_data_source_options_include_custom_simulation_once():
    labels = [opt["label"] for opt in app.build_data_source_options()]
    assert labels.count("Generate synthetic data") == 1
    assert labels[-3:] == [
        "Generate synthetic data",
        "Paste CSV/TSV text",
        "Upload your own file",
    ]
