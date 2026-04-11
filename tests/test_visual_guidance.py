import streamlit_app as app


def test_visual_interpretation_checklist_contract():
    checklist = app.visual_interpretation_checklist()
    assert not checklist.empty
    assert checklist["Priority"].is_monotonic_increasing

    required = {
        "Visualization",
        "Where",
        "PrimaryQuestion",
        "ReadFirst",
        "ReviewTrigger",
        "BeginnerAction",
        "Caveat",
    }
    assert required.issubset(checklist.columns)

    joined = " ".join(checklist["Visualization"].astype(str))
    for expected in [
        "Wright map",
        "Category probability curves",
        "Fit scatter",
        "Residual PCA",
        "Bias heatmap",
        "Strict marginal",
    ]:
        assert expected in joined


def test_visual_method_evidence_table_contract():
    evidence = app.visual_method_evidence_table()
    assert not evidence.empty
    assert {
        "Visualization",
        "Purpose",
        "MethodBasis",
        "AppReadabilityRule",
        "PrimaryReference",
    }.issubset(evidence.columns)
    joined = " ".join(evidence["Visualization"].astype(str))
    for expected in [
        "Wright map",
        "Category probability",
        "Bias heatmap",
        "Strict marginal",
    ]:
        assert expected in joined
    assert evidence["AppReadabilityRule"].str.contains("hover|compact|Limit", case=False, regex=True).any()


def test_category_probability_curve_data_supports_gpcm_level_views():
    result = {
        "params": {
            "steps_mat": app.np.array([[-1.2, 0.3], [-0.4, 1.1]], dtype=float),
            "slopes": app.np.array([0.8, 1.25], dtype=float),
        },
        "config": {"model": "GPCM", "n_cat": 3, "step_facet": "Task", "slope_facet": "Task"},
        "prep": {"rating_min": 1, "levels": {"Task": ["T1", "T2"]}},
    }

    options = app.category_probability_curve_options(result)
    assert options["Label"].tolist() == [
        "Average across Task levels",
        "Task = T1",
        "Task = T2",
    ]

    average = app.build_category_probability_curve_data(result)
    level = app.build_category_probability_curve_data(result, step_level_index=1)
    assert average["available"]
    assert level["available"]
    assert level["metadata"]["Level"] == "T2"
    assert level["metadata"]["Slope"] == 1.25
    fig, fig_expected = app._make_category_probability_curve_figures(level)
    assert fig is not None
    assert fig_expected is not None

    export = app.category_probability_curve_export_table(result)
    assert not export.empty
    assert {"CurveLabel", "Scope", "ExpectedScore", "Slope"}.issubset(export.columns)
    assert set(export["CurveLabel"]) == {"Average across Task levels", "Task = T1", "Task = T2"}
    assert not app.np.allclose(
        average["probability_wide"].drop(columns=["Theta"]).to_numpy(),
        level["probability_wide"].drop(columns=["Theta"]).to_numpy(),
    )
