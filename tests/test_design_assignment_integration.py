from __future__ import annotations

import pandas as pd

import streamlit_app as app


def _result(method: str = "JMLE") -> dict:
    data = pd.DataFrame(
        {
            "Person": ["p1", "p1", "p2", "p2"],
            "Rater": ["r1", "r2", "r1", "r2"],
            "Task": ["t1", "t1", "t1", "t1"],
            "Score": [0, 1, 1, 0],
        }
    )
    return {
        "config": {
            "method": method,
            "model": "RSM",
            "facet_names": ["Rater", "Task"],
            "n_cat": 2,
        },
        "prep": {"data": data, "facet_names": ["Rater", "Task"]},
        "summary": pd.DataFrame([{"Method": method, "Converged": True}]),
        "facets": {},
    }


def test_result_wrapper_uses_fitted_likelihood_rows_and_confirmed_rater_mapping():
    bundle = app.build_design_assignment_audit_for_result(_result())
    row = bundle["summary"].iloc[0]

    assert bundle["available"] is True
    assert row["LikelihoodRowsOnly"]
    assert row["RaterFacet"] == "Rater"
    assert row["RaterMappingConfirmed"]
    assert row["UniquePersonRaterAssignments"] == 4


def test_quick_result_bundle_exports_estimand_assignment_and_sensitivity_contracts():
    frames = app.build_result_bundle_frames(_result("MML"), {})

    required = {
        "estimator_estimand_contract",
        "assignment_design_summary",
        "assignment_rater_exposure",
        "assignment_rater_overlap",
        "fixed_density_sensitivity_plan",
        "informative_assignment_validation_evidence",
        "fixed_density_assignment_runner_gates",
        "fixed_density_assignment_block_profiles",
    }
    assert required.issubset(frames)
    estimand = frames["estimator_estimand_contract"]
    assert estimand.loc[estimand["CurrentRun"], "Method"].tolist() == ["MML"]
    assert frames["assignment_design_summary"].iloc[0]["OutcomeBlindAudit"]
    assert not frames["fixed_density_assignment_runner_gates"].empty


def test_design_assignment_locale_keys_have_english_japanese_parity():
    en = app._load_locale("en")["design_assignment"]
    ja = app._load_locale("ja")["design_assignment"]

    assert set(en) == set(ja)
    assert len(en) >= 20
