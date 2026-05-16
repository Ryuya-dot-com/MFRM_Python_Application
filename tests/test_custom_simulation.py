"""Custom synthetic-data generator contracts."""

from __future__ import annotations

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
        "meta",
    }
    assert bundle["meta"]["facet_names"] == ["Judge", "Prompt"]
    assert bundle["category_counts"]["Count"].sum() == len(bundle["data"])


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
