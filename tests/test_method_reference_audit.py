"""Reference-audit contract for manuscript-facing method claims."""

from __future__ import annotations

import streamlit_app as app


def test_method_reference_audit_is_ready_and_zotero_aligned():
    audit = app.build_method_reference_audit()
    assert not audit.empty
    required = {
        "MethodArea",
        "AppSurface",
        "PrimaryReferenceKeys",
        "CitationTokens",
        "ReferenceCoverageStatus",
        "BibExportStatus",
        "ZoteroAlignment",
        "ManuscriptUse",
        "ClaimBoundary",
    }
    assert required.issubset(audit.columns)
    assert audit["ReferenceCoverageStatus"].eq("Ready").all()
    assert audit["BibExportStatus"].eq("Ready").all()
    assert audit["ZoteroAlignment"].str.contains("Zotero|BibTeX", case=False, regex=True).all()


def test_method_reference_audit_covers_statistical_risk_surfaces():
    audit = app.build_method_reference_audit()
    areas = " ".join(audit["MethodArea"].astype(str))
    for expected in [
        "MML",
        "Fit",
        "Local dependence",
        "Residual PCA",
        "Simulation",
        "External R ecosystem",
        "Recent MFRM extensions",
    ]:
        assert expected in areas
    joined_keys = " ".join(audit["PrimaryReferenceKeys"].astype(str))
    for key in [
        "Rasch_1980",
        "Andersen_1973",
        "Warm_1989",
        "Yen_1984",
        "Christensen_Makransky_Horton_2017",
        "Mair_Hatzinger_2007",
        "Rizopoulos_2006",
        "Buerkner_2021",
        "Uto_Ueno_2020",
        "Uto_2021",
        "Uto_2023",
    ]:
        assert key in joined_keys


def test_result_bundle_includes_method_reference_audit():
    result = {
        "summary": app.pd.DataFrame([{"Model": "RSM"}]),
        "config": {"model": "RSM", "method": "JMLE", "facet_names": ["Person", "Rater"], "n_cat": 5},
        "prep": {"n_obs": 4, "n_person": 2, "rating_min": 0, "rating_max": 4},
    }
    diagnostics = {
        "measures": app.pd.DataFrame([
            {"Facet": "Rater", "Level": "R1", "Infit": 1.0, "SE": 0.2},
        ]),
        "obs": app.pd.DataFrame({"StdResidual": [0.0, 1.0]}),
    }
    frames = app.build_result_bundle_frames(result, diagnostics)
    assert "method_reference_audit" in frames
    assert not frames["method_reference_audit"].empty
    assert "apa_report_sentence_audit" in frames
    assert not frames["apa_report_sentence_audit"].empty


def test_current_run_method_reference_audit_filters_to_active_surfaces():
    result = {
        "config": {
            "model": "RSM",
            "method": "MML",
            "compute_plausible_values": True,
        },
        "prep": {"rating_min": 0, "rating_max": 4},
        "posterior": {
            "available": True,
            "plausible_values": app.pd.DataFrame({"Person": ["P1"], "Draw": [1], "Value": [0.2]}),
        },
    }
    diagnostics = {
        "obs": app.pd.DataFrame({"StdResidual": [0.0, 2.1]}),
        "pca_enabled": True,
        "pca_eigenvalues": app.pd.DataFrame({"Component": [1], "Eigenvalue": [1.4]}),
    }
    audit = app.build_current_run_method_reference_audit(
        result,
        diagnostics,
        all_bias_results={"Rater x Task": {"dff": app.pd.DataFrame({"Cell": ["R1:T1"]})}},
    )
    areas = " ".join(audit["MethodArea"].astype(str).tolist())
    assert "MML, EM, EAP" in areas
    assert "Local dependence" in areas
    assert "Residual PCA" in areas
    assert "Simulation" in areas
    assert {"CurrentRunUse", "TriggerEvidence"}.issubset(audit.columns)
    assert audit["TriggerEvidence"].astype(str).str.len().gt(0).all()


def test_claim_matrix_and_template_surface_reference_boundaries():
    result = {
        "config": {
            "model": "RSM",
            "method": "MML",
            "facet_names": ["Person", "Rater"],
            "n_cat": 5,
            "compute_plausible_values": False,
        },
        "prep": {"n_obs": 20, "n_person": 5, "rating_min": 0, "rating_max": 4},
    }
    diagnostics = {
        "obs": app.pd.DataFrame({"StdResidual": [0.0, 0.5, 2.2]}),
        "measures": app.pd.DataFrame({
            "Facet": ["Person", "Rater"],
            "Level": ["P1", "R1"],
            "Estimate": [0.1, -0.1],
            "SE": [0.2, 0.2],
            "Infit": [1.0, 1.1],
        }),
        "reliability": app.pd.DataFrame({"Facet": ["Person"], "Reliability": [0.8], "Separation": [2.0]}),
        "pca_enabled": False,
    }
    claim_evidence = app.build_claim_to_evidence_matrix(result, diagnostics, all_bias_results={})
    assert {
        "ReferenceAuditAreas",
        "SuggestedCitations",
        "CitationBoundary",
        "CitationEvidenceFile",
    }.issubset(claim_evidence.columns)
    model_row = claim_evidence.loc[
        claim_evidence["ManuscriptArea"].astype(str) == "Model and software scope"
    ].iloc[0]
    assert "Rasch" in model_row["SuggestedCitations"]
    assert model_row["CitationEvidenceFile"] == "method_reference_audit.csv"

    template = app.generate_manuscript_reporting_template(result, diagnostics, all_bias_results={})
    assert "## Method References And Claim Boundaries" in template
    assert "method_reference_audit.csv" in template
    assert "claim_to_evidence_matrix.csv" in template


def test_apa_report_sentence_audit_maps_sentences_to_evidence_and_boundaries():
    result = {
        "config": {
            "model": "RSM",
            "method": "MML",
            "facet_names": ["Person", "Rater"],
            "n_cat": 5,
            "population_prior_sd": 1.0,
        },
        "prep": {"n_obs": 30, "n_person": 10, "rating_min": 0, "rating_max": 4},
        "steps": app.pd.DataFrame({"Step": ["Step1", "Step2", "Step3"], "Estimate": [-1.0, 0.0, 1.0]}),
    }
    diagnostics = {
        "obs": app.pd.DataFrame({"StdResidual": [0.0, 0.5, 2.2]}),
        "measures": app.pd.DataFrame({
            "Facet": ["Person", "Rater"],
            "Level": ["P1", "R1"],
            "Estimate": [0.1, -0.1],
            "SE": [0.2, 0.2],
            "Infit": [1.0, 1.1],
        }),
        "reliability": app.pd.DataFrame({"Facet": ["Person"], "Reliability": [0.8], "Separation": [2.0]}),
        "pca_enabled": True,
        "pca_eigenvalues": app.pd.DataFrame({"Component": [1], "Eigenvalue": [1.5]}),
    }
    audit = app.build_apa_report_sentence_audit(result, diagnostics, all_bias_results={})
    assert not audit.empty
    required = {
        "SentenceOrder",
        "APASection",
        "SentenceRole",
        "LinkedManuscriptArea",
        "DraftSentence",
        "EvidenceFiles",
        "SuggestedCitations",
        "CitationBoundary",
        "DoNotClaim",
        "BeforeSubmissionCheck",
        "CopyDecision",
        "CitationEvidenceFile",
    }
    assert required.issubset(audit.columns)
    assert audit["SentenceOrder"].tolist() == sorted(audit["SentenceOrder"].tolist())
    assert audit["DraftSentence"].astype(str).str.len().gt(0).all()
    assert audit["SuggestedCitations"].astype(str).str.contains("Rasch|Bock|Wright|Linacre", regex=True).any()
    assert audit["CopyDecision"].astype(str).str.len().gt(0).all()
