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
        "SuggestedAction",
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


def test_visual_claim_guardrail_table_contract():
    guardrails = app.visual_claim_guardrail_table()
    assert not guardrails.empty
    assert {
        "GuardrailID",
        "Visualization",
        "WhatCanBeRead",
        "SafeReportWording",
        "DoNotWrite",
        "RequiredEvidence",
        "ReviewTrigger",
        "NextAction",
        "ArchiveFile",
    }.issubset(guardrails.columns)

    joined = " ".join(guardrails["Visualization"].astype(str))
    for expected in [
        "Wright Map",
        "FACETS-style yardstick",
        "Threshold / step lines",
        "Category characteristic curves",
    ]:
        assert expected in joined

    yardstick = guardrails.loc[guardrails["GuardrailID"] == "facets_yardstick_text_map"].iloc[0]
    assert "text map" in yardstick["SafeReportWording"]
    assert "density" in yardstick["DoNotWrite"]
    assert "visual_claim_guardrails.csv" in set(guardrails["ArchiveFile"])


def test_visual_evidence_map_inherits_claim_guardrails():
    figure_manifest = app.pd.DataFrame({
        "FigureName": [
            "wright_map",
            "facets_yardstick_text_map",
            "category_probability_curve_Average",
        ],
        "FormatsAvailable": ["HTML; PNG", "HTML", "HTML; PNG"],
        "Theme": ["Manuscript white", "Manuscript white", "Manuscript white"],
        "BaseFontSize": [13, 13, 13],
        "LabelPolicy": ["Auto", "Auto", "Auto"],
        "CaptionDetail": ["Detailed", "Detailed", "Detailed"],
    })
    visual_map = app.build_visual_evidence_map(figure_manifest, app.pd.DataFrame())

    assert {
        "VisualGuardrailID",
        "SafeReportWording",
        "DoNotWrite",
        "RequiredEvidence",
        "VisualGuardrailArchive",
    }.issubset(visual_map.columns)
    assert "wright_map_targeting" in visual_map["VisualGuardrailID"].tolist()
    assert "facets_yardstick_text_map" in visual_map["VisualGuardrailID"].tolist()
    assert "category_characteristic_curves" in visual_map["VisualGuardrailID"].tolist()
    assert visual_map["DoNotWrite"].str.contains("alone|solely", case=False, regex=True).any()

    captions = app.generate_visual_caption_drafts(visual_map)
    assert "Safe wording" in captions
    assert "Do not write" in captions


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
    assert level["metadata"]["RatingMin"] == 1
    assert level["thresholds"]["DisplayLabel"].tolist() == ["1|2", "2|3"]
    diagnostic = app.category_probability_curve_diagnostic_table(level)
    assert not diagnostic.empty
    assert {
        "PeakTheta",
        "PeakProbability",
        "IsEverMostProbable",
        "ReviewFlag",
    }.issubset(diagnostic.columns)
    fig, fig_expected = app._make_category_probability_curve_figures(level)
    assert fig is not None
    assert fig_expected is not None
    assert any(getattr(trace, "name", "") == "Step threshold lines" for trace in fig.data)
    assert any(str(getattr(trace, "name", "")).endswith("peak") for trace in fig.data)
    annotation_text = " ".join(str(item.text) for item in fig.layout.annotations)
    assert "1|2" in annotation_text and "2|3" in annotation_text

    export = app.category_probability_curve_export_table(result)
    assert not export.empty
    assert {
        "CurveLabel",
        "Scope",
        "ExpectedScore",
        "Slope",
        "PeakTheta",
        "PeakProbability",
        "IsEverMostProbable",
    }.issubset(export.columns)
    assert set(export["CurveLabel"]) == {"Average across Task levels", "Task = T1", "Task = T2"}
    assert not app.np.allclose(
        average["probability_wide"].drop(columns=["Theta"]).to_numpy(),
        level["probability_wide"].drop(columns=["Theta"]).to_numpy(),
    )


def test_rating_scale_dashboard_links_counts_thresholds_and_category_curves():
    result = {
        "params": {
            "steps_mat": app.np.array([[-1.2, 0.3], [-0.4, 1.1]], dtype=float),
            "slopes": app.np.array([0.8, 1.25], dtype=float),
        },
        "config": {"model": "GPCM", "n_cat": 3, "step_facet": "Task", "slope_facet": "Task"},
        "prep": {"rating_min": 1, "levels": {"Task": ["T1", "T2"]}},
        "steps": app.pd.DataFrame({
            "StepFacet": ["T1", "T1", "T2", "T2"],
            "Step": ["Step_1", "Step_2", "Step_1", "Step_2"],
            "Estimate": [-1.2, 0.3, -0.4, 1.1],
        }),
        "slopes": app.pd.DataFrame({
            "StepFacet": ["Task", "Task"],
            "Level": ["T1", "T2"],
            "Slope": [0.8, 1.25],
        }),
    }
    observed = app.np.array([1, 1, 1, 1, 2, 2, 3, 3], dtype=float)
    diagnostics = {
        "obs": app.pd.DataFrame({
            "Observed": observed,
            "Score": observed,
            "Var": app.np.ones(len(observed)),
            "StdSq": app.np.ones(len(observed)),
            "PersonMeasure": [-1.8, -1.5, -1.1, -0.8, 0.1, 0.4, 1.0, 1.3],
            "Expected": observed,
            "Residual": app.np.zeros(len(observed)),
        })
    }

    scope_diag = app.category_curve_diagnostic_scope_table(result)
    assert set(scope_diag["CurveLabel"]) == {"Average across Task levels", "Task = T1", "Task = T2"}

    dashboard = app.rating_scale_functioning_dashboard(result, diagnostics)
    assert not dashboard.empty
    assert {
        "Check",
        "Status",
        "Evidence",
        "InterpretationBoundary",
        "NextAction",
        "EvidenceSource",
    }.issubset(dashboard.columns)
    assert "Overall rating-scale decision" in dashboard["Check"].tolist()
    assert "Category characteristic curves" in dashboard["Check"].tolist()
    assert "Observed category support" in dashboard["Check"].tolist()
    assert (dashboard["Status"] == "Review").any()

    category_evidence = app.rating_scale_category_evidence_table(result, diagnostics)
    assert not category_evidence.empty
    assert {
        "CategoryEvidenceStatus",
        "EvidenceSummary",
        "CurvePeakProbability",
        "CurveEverMostProbable",
    }.issubset(category_evidence.columns)
    assert (category_evidence["CategoryEvidenceStatus"] == "Review").any()

    decision_support = app.rating_scale_decision_support_table(result, diagnostics)
    assert not decision_support.empty
    assert {
        "DecisionArea",
        "DecisionStatus",
        "AffectedCategories",
        "PrimaryEvidence",
        "RecommendedAction",
        "RerunComparison",
        "SafeReportWording",
        "DoNotClaim",
    }.issubset(decision_support.columns)
    assert "Overall reportability" in decision_support["DecisionArea"].tolist()
    assert "Potential recoding scope" in decision_support["DecisionArea"].tolist()
    assert (decision_support["DecisionStatus"].isin(["Do not claim yet", "Review before reporting"])).any()
    assert decision_support["RerunComparison"].str.contains("original scoring", case=False, regex=False).any()
    assert decision_support["DoNotClaim"].str.contains("collapse categories solely", case=False, regex=False).any()

    recode_candidates = app.rating_scale_recode_candidate_table(result, diagnostics)
    assert not recode_candidates.empty
    assert {
        "CandidateID",
        "AnalysisRole",
        "AdjacentBoundary",
        "CollapsedOriginalCategories",
        "OriginalToCandidateScoreMap",
        "Priority",
        "TriggerEvidence",
        "RequiredSubstantiveJustification",
        "RerunComparisonPlan",
        "ReportingBoundary",
    }.issubset(recode_candidates.columns)
    assert "BASE" in recode_candidates["CandidateID"].tolist()
    assert (recode_candidates["AnalysisRole"] == "Sensitivity comparison candidate").any()
    assert recode_candidates["OriginalToCandidateScoreMap"].str.contains("2 -> 1", regex=False).any()
    assert recode_candidates["RerunComparisonPlan"].str.contains("refit the model", case=False, regex=False).any()

    recode_map_long = app.rating_scale_recode_map_long_table(result, diagnostics)
    assert not recode_map_long.empty
    assert {
        "CandidateID",
        "CandidateScoreColumn",
        "OriginalScore",
        "CandidateScore",
        "CandidateFingerprint",
    }.issubset(recode_map_long.columns)
    assert "Score_recode_C01" in recode_map_long["CandidateScoreColumn"].tolist()
    assert (
        recode_map_long
        .loc[recode_map_long["CandidateScoreColumn"] == "Score_recode_C01", "CandidateScore"]
        .astype(float)
        .min()
        == 1.0
    )

    recode_assets = app.rating_scale_recode_script_assets(result, diagnostics)
    assert {
        "README_rating_scale_recode_scripts.md",
        "mfrm_rating_scale_recode_candidates.csv",
        "mfrm_rating_scale_recode_map_long.csv",
        "apply_rating_scale_recodes.py",
        "apply_rating_scale_recodes.R",
        "apply_rating_scale_recodes.jl",
    }.issubset(recode_assets)
    assert "MFRM_RECODE_MAP_CSV" in recode_assets["apply_rating_scale_recodes.py"]
    assert "Score_original" in recode_assets["README_rating_scale_recode_scripts.md"]


def test_rating_scale_help_table_keeps_curve_evidence_in_workflow():
    help_table = app.rating_scale_evidence_help_table()
    assert not help_table.empty
    assert {
        "EvidenceSurface",
        "Where",
        "WhatItAnswers",
        "ReviewTrigger",
        "NextAction",
    }.issubset(help_table.columns)
    joined = " ".join(help_table["EvidenceSurface"].astype(str))
    assert "Category characteristic curves" in joined
    assert "Threshold / step ordering" in joined
