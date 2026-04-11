import io

import streamlit_app as app


def test_cross_package_validation_plan_contract():
    plan = app.cross_package_validation_plan()
    assert not plan.empty
    assert {"TAM", "mirt", "sirt"}.issubset(set(plan["ExternalReference"].str.split().str[0]))
    assert "mml_rsm_latent_regression" in " ".join(plan["PythonScenario"].astype(str))
    assert "mml_gpcm_task_steps" in " ".join(plan["PythonScenario"].astype(str))

    notes = app.cross_package_parameterization_notes()
    assert "Latent variance" in notes["Topic"].tolist()
    assert "GPCM slope identification" in notes["Topic"].tolist()

    tolerance = app.cross_package_tolerance_policy()
    assert "Log-likelihood" in tolerance["EvidenceType"].tolist()
    assert "Plausible-value draws" in tolerance["EvidenceType"].tolist()

    docs = app.external_reference_documentation_table()
    assert {"TAM tam.mml.mfr", "mirt mirt", "sirt rm.facets", "mfrmr 0.1.5 local source"}.issubset(set(docs["Reference"]))

    artifacts = app.external_validation_artifact_checklist()
    assert "External package versions" in artifacts["Artifact"].tolist()

    template = app.external_validation_report_template()
    assert "mfrmr 0.1.5 migration" in template["ClaimArea"].tolist()

    coverage = app.mfrmr_015_migration_coverage_table()
    assert "Bounded GPCM" in coverage["mfrmr015Area"].tolist()
    assert "Latent regression / population_formula" in coverage["mfrmr015Area"].tolist()


def test_final_readiness_uses_five_percent_residual_benchmark():
    class Opt:
        success = True
        message = "ok"

    result = {
        "opt": Opt(),
        "config": {},
        "prep": {},
        "facets": {"person": app.pd.DataFrame(), "others": app.pd.DataFrame()},
    }
    diagnostics_ready = {
        "obs": app.pd.DataFrame({"StdResidual": [2.1] * 5 + [0.0] * 95}),
        "reliability": app.pd.DataFrame(),
        "pca_enabled": False,
    }
    diagnostics_review = {
        "obs": app.pd.DataFrame({"StdResidual": [2.1] * 6 + [0.0] * 94}),
        "reliability": app.pd.DataFrame(),
        "pca_enabled": False,
    }

    ready = app.build_final_report_readiness(result, diagnostics_ready, all_bias_results={})
    review = app.build_final_report_readiness(result, diagnostics_review, all_bias_results={})
    ready_status = ready.loc[ready["Check"] == "Global residual fit", "Status"].iloc[0]
    review_status = review.loc[review["Check"] == "Global residual fit", "Status"].iloc[0]

    assert app.FINAL_RESIDUAL_PCT_GE2_READY == 5.0
    assert ready_status == "Ready"
    assert review_status == "Review"


def test_manuscript_claim_guide_contract():
    class Opt:
        success = True
        message = "ok"

    result = {
        "opt": Opt(),
        "config": {
            "app_version": app.APP_VERSION,
            "model": "RSM",
            "method": "JMLE",
            "facet_names": ["Rater", "Task"],
            "anchor_audit": {"overall_status": "ok", "message": "No anchor issues."},
        },
        "prep": {"n_obs": 100, "n_person": 20, "rating_min": 1, "rating_max": 5},
        "facets": {
            "person": app.pd.DataFrame({"Estimate": [-0.5, 0.5]}),
            "others": app.pd.DataFrame({"Estimate": [-0.2, 0.2]}),
        },
    }
    diagnostics = {
        "obs": app.pd.DataFrame({
            "Observed": [1, 2, 3, 4, 5] * 20,
            "StdResidual": [2.1] * 5 + [0.0] * 95,
        }),
        "reliability": app.pd.DataFrame({"Facet": ["Person"], "Reliability": [0.85]}),
        "pca_enabled": False,
    }

    guide = app.build_manuscript_claim_guide(result, diagnostics, all_bias_results={})
    assert not guide.empty
    assert {
        "ManuscriptArea",
        "ClaimStatus",
        "SafeManuscriptWording",
        "EvidenceToReport",
        "DoNotClaim",
        "NextAction",
    }.issubset(guide.columns)
    assert "External package comparison" in guide["ManuscriptArea"].tolist()
    boundary = guide.loc[guide["ManuscriptArea"] == "External package comparison"].iloc[0]
    assert boundary["ClaimStatus"] == "Boundary"
    assert "Do not force" in boundary["DoNotClaim"]
    bias = guide.loc[guide["ManuscriptArea"] == "Bias / local interaction"].iloc[0]
    assert bias["ClaimStatus"] == "Do not claim"
    fit = guide.loc[guide["ManuscriptArea"] == "Fit and dimensionality"].iloc[0]
    assert fit["ClaimStatus"] == "Report with caveat"

    template = app.generate_manuscript_reporting_template(result, diagnostics, all_bias_results={})
    assert "## Methods Draft" in template
    assert "## Results Draft" in template
    assert "standalone Python" in template
    assert "Claims Requiring Caution" in template
    assert "Do not include an R-vs-Python comparison table" in template

    gate = app.build_publication_gate_summary(result, diagnostics, all_bias_results={})
    assert not gate.empty
    assert {"GateArea", "GateStatus", "Evidence", "ManuscriptAction"}.issubset(gate.columns)
    overall = gate.loc[gate["GateArea"] == "Overall manuscript gate"].iloc[0]
    assert overall["GateStatus"] == "Report with caveat"
    assert "Fit and dimensionality" in gate["GateArea"].astype(str).tolist()
    assert "Bias / local interaction" in gate["GateArea"].astype(str).tolist()


def test_uploaded_file_fingerprint_uses_content_not_only_name_and_size():
    first = io.BytesIO(b"abc")
    second = io.BytesIO(b"abd")
    first.name = "ratings.csv"
    second.name = "ratings.csv"
    first.size = 3
    second.size = 3

    assert app.uploaded_file_fingerprint(first) != app.uploaded_file_fingerprint(second)


def test_population_covariate_type_summary_flags_integer_codes():
    person_data = app.pd.DataFrame({
        "Person": ["P1", "P2", "P3", "P4"],
        "GradeCode": [1, 2, 1, 2],
        "SES": [0.1, -0.4, 0.6, 0.2],
        "Group": ["A", "A", "B", "B"],
    })
    summary = app.summarize_population_covariate_types(
        person_data,
        "Person",
        "~ GradeCode + SES + Group",
        person_levels=["P1", "P2", "P3", "P4"],
    )
    grade_code = summary.loc[summary["Term"] == "GradeCode"].iloc[0]
    ses = summary.loc[summary["Term"] == "SES"].iloc[0]

    assert grade_code["InferredType"] == "numeric"
    assert bool(grade_code["ReviewFlag"])
    assert ses["InferredType"] == "numeric"
    assert not bool(ses["ReviewFlag"])

    forced = app.summarize_population_covariate_types(
        person_data,
        "Person",
        "~ GradeCode + SES + Group",
        categorical_terms=["GradeCode"],
        person_levels=["P1", "P2", "P3", "P4"],
    )
    forced_grade_code = forced.loc[forced["Term"] == "GradeCode"].iloc[0]
    assert forced_grade_code["InferredType"] == "categorical"
    assert forced_grade_code["Override"] == "categorical"
    assert not bool(forced_grade_code["ReviewFlag"])


def test_public_beta_release_contract_tables():
    limitations = app.public_beta_limitations_table()
    readiness = app.public_release_readiness_table()

    assert not limitations.empty
    assert {"Area", "PublicBetaStatus", "SupportedNow", "Boundary", "UserAction"}.issubset(limitations.columns)
    joined = " ".join(limitations["Area"].astype(str).tolist() + limitations["Boundary"].astype(str).tolist())
    assert "GPCM" in joined
    assert "Cross-package" in joined
    assert "privacy" in joined

    assert not readiness.empty
    assert {"Check", "Status", "Evidence", "Action"}.issubset(readiness.columns)
    assert "License" in readiness["Check"].astype(str).tolist()
    assert not (readiness["Status"].astype(str) == "Blocker").any()
