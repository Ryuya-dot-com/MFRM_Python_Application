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
