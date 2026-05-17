"""Tests for the parameter-recovery (ADEMP) simulation pipeline.

The data-generating mechanism samples a many-facet ground truth with
mean-zero location blocks. ``evaluate_parameter_recovery`` refits the
requested model on each replicate, mean-aligns the location-block
estimates against the truth, and reports per-row error plus aggregate
bias / RMSE / coverage95.

The math contract pinned here covers:

* Refusal on unsupported model / fit_method / n_cat inputs.
* Long-table shape (one row per parameter level per rep) and the
  three location-block ParameterType values (Person, Rater,
  Criterion) when ``include_person=True``.
* The per-rep mean-alignment identity: at every rep,
  ``sum(EstimateAligned) ≈ 0`` and ``sum(ErrorAligned) ≈ 0`` over
  each (rep, ParameterType) subset.
* The summary's Bias is the mean of ErrorAligned, RMSE is the
  root-mean-square, MAE is the mean absolute, and the row count
  matches ``reps * n_levels`` for each location block.
* Coverage95 is the fraction of finite-SE rows whose truth lies in
  ``estimate +/- 1.96 * SE``; the rate must equal the closed-form
  fraction from the recovery rows.
* SE/CI method/status metadata are preserved in the long table, and
  the explicit SE/CI coverage report carries nominal coverage,
  coverage error, interpretation, and recommended action fields.
* Reproducibility: same seed produces the same recovery row for the
  same (rep, ParameterType, Level) triple.
* ADEMP bundle: Aims / DataGenerating / Estimands / Methods /
  PerformanceMeasures fields populated as non-empty strings.

Fits use small fixtures (10-30 persons, 2-3 raters, 3 criteria, 3-4
score levels, 5-10 reps) to keep the test suite cheap; the
correctness contract is what the tests pin, not the absolute RMSE.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy.stats import norm

import streamlit_app as app


# -----------------------------------------------------------------------------
# Unit: refusal on bad inputs
# -----------------------------------------------------------------------------


def test_evaluate_parameter_recovery_refuses_unknown_model():
    out = app.evaluate_parameter_recovery(model="XYZ", reps=1)
    assert out["available"] is False
    assert "Unsupported model" in out["reason"]


def test_evaluate_parameter_recovery_refuses_unknown_fit_method():
    out = app.evaluate_parameter_recovery(model="RSM", fit_method="Bayes", reps=1)
    assert out["available"] is False
    assert "Unsupported fit_method" in out["reason"]


def test_evaluate_parameter_recovery_refuses_one_category():
    out = app.evaluate_parameter_recovery(model="RSM", n_cat=1, reps=1)
    assert out["available"] is False
    assert "n_cat" in out["reason"]


# -----------------------------------------------------------------------------
# End-to-end small fixture
# -----------------------------------------------------------------------------


@pytest.fixture(scope="module")
def small_rsm_recovery_bundle():
    """Reused across the long-table contract tests so the simulation
    cost is paid once."""
    return app.evaluate_parameter_recovery(
        n_person=20, n_rater=2, n_criterion=3, n_cat=3,
        reps=5, model="RSM", fit_method="JMLE",
        seed=20260605, maxit=20, reltol=1e-3,
    )


def test_recovery_bundle_is_available_and_carries_three_parameter_types(
    small_rsm_recovery_bundle,
):
    assert small_rsm_recovery_bundle["available"] is True
    recovery = small_rsm_recovery_bundle["recovery"]
    assert {"Person", "Rater", "Criterion"}.issubset(
        set(recovery["ParameterType"].astype(str).unique())
    )


def test_recovery_rep_overview_has_one_row_per_rep(small_rsm_recovery_bundle):
    rep_overview = small_rsm_recovery_bundle["rep_overview"]
    assert len(rep_overview) == 5
    assert (rep_overview["RunOK"]).all()
    assert (rep_overview["RecoveryRows"] > 0).all()


def test_recovery_long_table_row_counts_match_reps_times_levels(
    small_rsm_recovery_bundle,
):
    """Each (ParameterType, rep) subset must have exactly n_levels rows."""
    recovery = small_rsm_recovery_bundle["recovery"]
    counts = recovery.groupby(["ParameterType", "rep"]).size().unstack(fill_value=0)
    # Person: 20 levels; Rater: 2; Criterion: 3.
    assert (counts.loc["Person"] == 20).all()
    assert (counts.loc["Rater"] == 2).all()
    assert (counts.loc["Criterion"] == 3).all()


# -----------------------------------------------------------------------------
# Math contract: mean-alignment identity within (rep, ParameterType)
# -----------------------------------------------------------------------------


def test_recovery_mean_alignment_identity_holds_per_rep(
    small_rsm_recovery_bundle,
):
    """Mean-aligned estimates and aligned truths must each sum to zero
    over the per-rep, per-parameter-type subset (Rasch identification
    convention; cf. Morris et al. 2019 alignment recipe)."""
    recovery = small_rsm_recovery_bundle["recovery"]
    grouped = recovery.groupby(["rep", "ParameterType"])
    for (rep, ptype), sub in grouped:
        if sub.empty:
            continue
        est_aligned = pd.to_numeric(sub["EstimateAligned"], errors="coerce").to_numpy()
        truth_aligned = pd.to_numeric(sub["Truth_Aligned"], errors="coerce").to_numpy()
        finite_est = est_aligned[np.isfinite(est_aligned)]
        finite_truth = truth_aligned[np.isfinite(truth_aligned)]
        if finite_est.size:
            assert abs(float(finite_est.mean())) < 1e-10, (rep, ptype, finite_est.mean())
        if finite_truth.size:
            assert abs(float(finite_truth.mean())) < 1e-10, (rep, ptype)


# -----------------------------------------------------------------------------
# Summary: Bias / RMSE / MAE match the underlying long-table values
# -----------------------------------------------------------------------------


def test_recovery_summary_bias_matches_mean_error_aligned(
    small_rsm_recovery_bundle,
):
    """The summary Bias is the mean of ErrorAligned across all rows in
    the (ParameterType, Facet, ComparisonScale) group."""
    recovery = small_rsm_recovery_bundle["recovery"]
    summary = small_rsm_recovery_bundle["recovery_summary"]
    for _, row in summary.iterrows():
        sub = recovery[
            (recovery["ParameterType"] == row["ParameterType"])
            & (recovery["Facet"] == row["Facet"])
            & (recovery["ComparisonScale"] == row["ComparisonScale"])
        ]
        err = pd.to_numeric(sub["ErrorAligned"], errors="coerce").dropna()
        assert float(row["Bias"]) == pytest.approx(float(err.mean()), abs=1e-12)
        assert float(row["RMSE"]) == pytest.approx(
            float(np.sqrt((err ** 2).mean())), abs=1e-12
        )
        assert float(row["MAE"]) == pytest.approx(float(err.abs().mean()), abs=1e-12)


def test_recovery_summary_coverage95_matches_closed_form(
    small_rsm_recovery_bundle,
):
    """Coverage95 must equal the fraction of finite-SE rows with
    truth in estimate +/- 1.96 * SE."""
    recovery = small_rsm_recovery_bundle["recovery"]
    summary = small_rsm_recovery_bundle["recovery_summary"]
    z95 = float(norm.ppf(0.975))
    for _, srow in summary.iterrows():
        sub = recovery[
            (recovery["ParameterType"] == srow["ParameterType"])
            & (recovery["Facet"] == srow["Facet"])
            & (recovery["ComparisonScale"] == srow["ComparisonScale"])
        ]
        se = pd.to_numeric(sub["SE"], errors="coerce")
        err = pd.to_numeric(sub["ErrorAligned"], errors="coerce")
        ok = np.isfinite(se) & (se > 0) & np.isfinite(err)
        if not ok.any():
            continue
        expected_coverage = float(
            (np.abs(err[ok]) <= z95 * se[ok]).mean()
        )
        # The summary value may be NaN if the helper found no rows;
        # we already filtered to ok.any() so it must be finite here.
        assert float(srow["Coverage95"]) == pytest.approx(expected_coverage, abs=1e-12)


def test_recovery_long_table_carries_se_ci_metadata(small_rsm_recovery_bundle):
    recovery = small_rsm_recovery_bundle["recovery"]
    expected_cols = {
        "NominalCILevel",
        "CoverageMethod",
        "SE_Method",
        "SE_Status",
        "CI_Method",
        "CI_Status",
        "UncertaintyCaution",
    }
    assert expected_cols.issubset(recovery.columns)
    assert set(recovery["NominalCILevel"].dropna().unique()) == {0.95}
    assert recovery["CoverageMethod"].astype(str).str.contains("mean-aligned Wald").all()
    assert recovery["SE_Status"].astype(str).str.len().gt(0).all()
    assert recovery["CI_Status"].astype(str).str.len().gt(0).all()


def test_recovery_summary_carries_explicit_coverage_diagnostic(
    small_rsm_recovery_bundle,
):
    summary = small_rsm_recovery_bundle["recovery_summary"]
    expected_cols = {
        "CoverageN",
        "NominalCoverage",
        "CoverageErrorFromNominal",
        "PrimarySEStatus",
        "SEStatusSummary",
        "CIStatusSummary",
        "SEBasisRisk",
        "CoverageClaimStatus",
        "CoverageInterpretation",
        "RecommendedAction",
    }
    assert expected_cols.issubset(summary.columns)
    for _, row in summary.iterrows():
        if np.isfinite(float(row["Coverage95"])):
            assert float(row["CoverageErrorFromNominal"]) == pytest.approx(
                float(row["Coverage95"]) - 0.95, abs=1e-12
            )
        assert isinstance(row["CoverageInterpretation"], str)
        assert isinstance(row["RecommendedAction"], str)
        assert isinstance(row["SEBasisRisk"], str)
        assert row["CoverageClaimStatus"] in {
            "Ready",
            "Report with caveat",
            "Screening only",
            "Do not claim",
        }
        assert len(row["RecommendedAction"]) > 0


def test_build_se_ci_coverage_report_matches_recovery_summary(
    small_rsm_recovery_bundle,
):
    report = app.build_se_ci_coverage_report(small_rsm_recovery_bundle)
    summary = small_rsm_recovery_bundle["recovery_summary"]
    assert not report.empty
    assert {
        "CoverageN",
        "NominalCoverage",
        "Coverage95",
        "CoverageError",
        "CoverageAbsError",
        "PrimarySEStatus",
        "SEBasisRisk",
        "CoverageClaimStatus",
        "Interpretation",
        "RecommendedAction",
    }.issubset(report.columns)
    for _, row in report.iterrows():
        match = summary[
            (summary["ParameterType"] == row["ParameterType"])
            & (summary["Facet"] == row["Facet"])
            & (summary["ComparisonScale"] == row["ComparisonScale"])
        ]
        assert len(match) == 1
        srow = match.iloc[0]
        assert int(row["CoverageN"]) == int(srow["CoverageN"])
        assert float(row["Coverage95"]) == pytest.approx(float(srow["Coverage95"]), abs=1e-12)
        assert float(row["CoverageError"]) == pytest.approx(
            float(srow["Coverage95"]) - 0.95, abs=1e-12
        )


def test_evaluate_se_ci_coverage_wrapper_attaches_report():
    out = app.evaluate_se_ci_coverage(
        n_person=12, n_rater=2, n_criterion=2, n_cat=3,
        reps=2, model="RSM", fit_method="JMLE",
        seed=20260607, maxit=15, reltol=1e-3,
    )
    assert out["available"] is True
    assert "coverage_report" in out
    assert isinstance(out["coverage_report"], pd.DataFrame)
    assert not out["coverage_report"].empty


# -----------------------------------------------------------------------------
# Determinism: same seed = same recovery row content
# -----------------------------------------------------------------------------


def test_recovery_is_deterministic_under_fixed_seed():
    out_a = app.evaluate_parameter_recovery(
        n_person=15, n_rater=2, n_criterion=3, n_cat=3,
        reps=3, model="RSM", fit_method="JMLE",
        seed=2026, maxit=15, reltol=1e-3,
    )
    out_b = app.evaluate_parameter_recovery(
        n_person=15, n_rater=2, n_criterion=3, n_cat=3,
        reps=3, model="RSM", fit_method="JMLE",
        seed=2026, maxit=15, reltol=1e-3,
    )
    pd.testing.assert_frame_equal(
        out_a["recovery"].reset_index(drop=True),
        out_b["recovery"].reset_index(drop=True),
    )


# -----------------------------------------------------------------------------
# ADEMP bundle has Morris / White / Crowther fields
# -----------------------------------------------------------------------------


def test_recovery_ademp_bundle_carries_morris_2019_fields(
    small_rsm_recovery_bundle,
):
    ademp = small_rsm_recovery_bundle["ademp"]
    for key in ("Aims", "DataGenerating", "Estimands", "Methods",
                "PerformanceMeasures", "Reference"):
        assert key in ademp
        assert isinstance(ademp[key], str) and len(ademp[key]) > 0
    assert "Morris" in ademp["Reference"]


def test_recovery_settings_round_trip_input_arguments(
    small_rsm_recovery_bundle,
):
    settings = small_rsm_recovery_bundle["settings"]
    assert settings["model"] == "RSM"
    assert settings["fit_method"] == "JMLE"
    assert settings["n_person"] == 20
    assert settings["n_rater"] == 2
    assert settings["n_criterion"] == 3
    assert settings["n_cat"] == 3
    assert settings["reps"] == 5


# -----------------------------------------------------------------------------
# Performance bounds: bias is near zero on a clean RSM fit
# -----------------------------------------------------------------------------


def test_recovery_bias_is_near_zero_on_rsm_data(small_rsm_recovery_bundle):
    """Mean alignment forces Bias to zero by construction at every
    parameter type. Pin to machine precision."""
    summary = small_rsm_recovery_bundle["recovery_summary"]
    for _, row in summary.iterrows():
        assert abs(float(row["Bias"])) < 1e-10


def test_recovery_correlation_is_positive_for_well_identified_blocks():
    """With enough persons / reps, the (truth, aligned estimate)
    correlation must be clearly positive for the Person block under
    a clean RSM fit."""
    out = app.evaluate_parameter_recovery(
        n_person=40, n_rater=3, n_criterion=3, n_cat=4,
        reps=10, model="RSM", fit_method="JMLE",
        seed=20260606, maxit=25, reltol=1e-4,
    )
    summary = out["recovery_summary"]
    person_corr = float(
        summary[summary["ParameterType"] == "Person"]["Correlation"].iloc[0]
    )
    assert person_corr > 0.5  # well-identified Person block
