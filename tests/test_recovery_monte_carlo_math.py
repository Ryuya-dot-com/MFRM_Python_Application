"""Independent repetitions, not parameter rows, define Monte Carlo uncertainty."""

import numpy as np
import pandas as pd
import pytest

import streamlit_app as app


def rows(errors, *, aligned=False):
    records = []
    for rep, values in enumerate(errors, 1):
        for level, error in enumerate(values):
            records.append(dict(rep=rep, ParameterType="Rater", Facet="Rater", Level=str(level),
                ComparisonScale="logit_mean_aligned" if aligned else "identified_logit",
                ErrorAligned=error, ErrorRaw=error, Truth=level, EstimateAligned=level+error,
                TruthQuantile=f"Q{int(4*(level+.5)/len(values))+1}",
                SE=1., Covered95=abs(error) < 1.96, SE_Status="ok", CI_Status="ok"))
    return pd.DataFrame(records)


def test_mcse_does_not_shrink_when_correlated_levels_are_duplicated():
    for levels in (1, 20):
        summary = app._recovery_summarize(rows([[1.] * levels, [3.] * levels])).iloc[0]
        assert summary["Bias"] == 2.
        assert summary["MCSEBias"] == pytest.approx(1.)
        assert summary["RMSE"] == pytest.approx(np.sqrt(5.))
        assert summary["MCSERMSE"] == pytest.approx(2./np.sqrt(5.))
        assert summary["MCSECoverage95"] == pytest.approx(.5)
        assert summary["MCSEUnit"] == "independent simulation replicate"


def test_unequal_rows_do_not_reweight_independent_repetitions():
    summary = app._recovery_summarize(rows([[1.], [3., 3., 3.]])).iloc[0]
    assert summary["Bias"] == 2.
    assert summary["RMSE"] == pytest.approx(np.sqrt(5.))
    assert summary["Coverage95"] == .5
    assert summary["CoverageReps"] == 2


def test_many_levels_do_not_turn_two_repetitions_into_calibration_evidence():
    data = rows([[0.] * 19 + [3.]] * 2)
    summary = app._recovery_summarize(data)
    result = summary.iloc[0]
    assert result["CoverageN"] == 40
    assert result["CoverageReps"] == 2
    assert result["CoverageClaimStatus"] == "Screening only"
    report = app.build_se_ci_coverage_report({"recovery_summary": summary})
    assert "increase replicates" in report.iloc[0]["Interpretation"]


def test_centering_identity_is_separate_from_level_specific_bias():
    data = rows([[-1., 1.], [-2., 2.]], aligned=True)
    block = app._recovery_summarize(data).iloc[0]
    assert np.isnan(block["Bias"])
    assert np.isnan(block["MCSEBias"])
    assert block["MeanAlignedError"] == 0.
    assert block["BiasStatus"] == "not_evaluated_for_centered_block; inspect_level_bias"
    levels = app._recovery_summarize(data, by_level=True).set_index("Level")
    assert levels.loc["0", "Bias"] == -1.5
    assert levels.loc["1", "Bias"] == 1.5
    assert levels.loc["0", "MCSEBias"] == .5


def test_single_repetition_cannot_estimate_mcse():
    result = app._recovery_summarize(rows([[1., 2., 3.]])).iloc[0]
    assert np.isnan(result["MCSEBias"])
    assert np.isnan(result["MCSERMSE"])
    assert np.isnan(result["MCSECoverage95"])


@pytest.mark.parametrize("free_sd", [False, True])
def test_coverage_transforms_covariance_alongside_centering(free_sd):
    config = {"model": "RSM", "method": "MML", "n_cat": 2,
        "facet_names": ["Criterion"], "facet_levels": {"Criterion": ["C1", "C2"]},
        "facet_specs": {"Criterion": app.build_facet_constraint(["C1", "C2"], centered=False)},
        "estimate_population_sd": free_sd}
    res = {"config": config, "facets": {"others": pd.DataFrame({
        "Facet": ["Criterion"] * 2, "Level": ["C1", "C2"], "Estimate": [.2, -.2]})}}
    truth = {"Criterion": {"C1": 0., "C2": 0.}}
    measures = pd.DataFrame({"Facet": ["Criterion"] * 2, "Level": ["C1", "C2"], "SE": [.2, .3]})
    covariance = {"status": "ok", "cov": np.array([[.04, .01], [.01, .09]])}
    output = app._recovery_rows_from_fit(res, truth, rep=1, include_person=False,
                                         measures=measures, covariance=covariance)
    np.testing.assert_allclose(output["SE"], [.2, .3])  # Original diagnostic retained.
    if free_sd:
        assert output["CoverageSE"].isna().all()  # Fixed-SD block cannot stand in for joint covariance.
        assert output["Covered95"].isna().all()
    else:
        np.testing.assert_allclose(output["CoverageSE"], np.sqrt(.0275))
        assert output["Covered95"].all()
    without = app._recovery_rows_from_fit(res, truth, rep=1, include_person=False, measures=measures)
    assert without["CoverageSE"].isna().all()
    assert without["Covered95"].isna().all()


def test_legacy_row_count_does_not_preserve_a_ready_coverage_claim():
    legacy = pd.DataFrame([dict(CoverageN=1000, Coverage95=.95, CoverageClaimStatus="Ready")])
    report = app.build_se_ci_coverage_report({"recovery_summary": legacy})
    assert report.iloc[0]["CoverageClaimStatus"] == "Do not claim"
    with pytest.raises(ValueError, match="recomputing"):
        app.build_se_ci_coverage_report({"recovery_summary": legacy}, nominal=.90)


def test_coverage_point_estimate_alone_does_not_qualify_inference():
    assert app._coverage_claim_status(.95, 1000, "ok") == "Screening only"


@pytest.mark.parametrize("first_converged", [False, True])
def test_recovery_retains_failures_and_never_changes_summary_population(monkeypatch, first_converged):
    calls = 0

    def fit(**kwargs):
        nonlocal calls
        calls += 1
        if calls == 3:
            raise RuntimeError("deliberate fit failure")
        return {"summary": pd.DataFrame([{"Converged": first_converged and calls == 1}])}

    def recovery_rows(res, truth, *, rep, **kwargs):
        data = rows([[float(rep)]])
        data["rep"] = rep
        return data

    monkeypatch.setattr(app, "mfrm_estimate", fit)
    monkeypatch.setattr(app, "mfrm_diagnostics", lambda *args, **kwargs: {})
    monkeypatch.setattr(app, "_recovery_rows_from_fit", recovery_rows)
    bundle = app.evaluate_parameter_recovery(n_person=2, n_rater=2, n_criterion=2, reps=3)
    counts = bundle["execution_summary"].iloc[0]
    assert counts["PlannedReps"] == counts["AttemptedReps"] == 3
    assert counts["CompletedReps"] == 2
    assert counts["FailedReps"] == 1
    assert counts["ConvergedReps"] == int(first_converged)
    assert counts["NonconvergedCompletedReps"] == 2-int(first_converged)
    assert counts["FailureRate"] == pytest.approx(1/3)
    assert counts["MCSEFailureRate"] == pytest.approx(np.sqrt(2/27))
    assert len(bundle["recovery"]) == 2  # Both returned results remain inspectable.
    assert "deliberate fit failure" in bundle["rep_overview"].iloc[2]["Error"]
    assert bundle["summary_basis"] == "converged_only_no_fallback"
    if first_converged:
        summary = bundle["recovery_summary"].iloc[0]
        assert summary["Bias"] == 1.
        assert summary["PlannedReps"] == 3
        assert summary["ExcludedReps"] == 2
        assert summary["ValidReps"] == 1
    else:
        assert bundle["recovery_summary"].empty
        assert bundle["coverage_report"].empty
        # Reconstructing a report must not revive excluded nonconverged rows.
        assert app.build_se_ci_coverage_report(bundle).empty
        assert app._recovery_summarize(bundle["recovery"]).empty


@pytest.mark.parametrize("model,n_cat", [("RSM", 2), ("RSM", 4), ("PCM", 2), ("PCM", 4), ("GPCM", 4)])
def test_recovery_truth_mapping_preserves_every_adjacent_log_odds(monkeypatch, model, n_cat):
    captured = {}
    simulate = app._simulate_recovery_responses

    def capture(**kwargs):
        captured.update(kwargs)
        return simulate(**kwargs)

    monkeypatch.setattr(app, "_simulate_recovery_responses", capture)
    monkeypatch.setattr(app, "mfrm_estimate", lambda **kwargs: {})
    monkeypatch.setattr(app, "mfrm_diagnostics", lambda *args, **kwargs: {})
    bundle = app.evaluate_parameter_recovery(model=model, n_cat=n_cat, reps=1,
        n_person=3, n_rater=2, n_criterion=3, slopes=np.array([.6, 1.7, 2.3]))
    truth = bundle["truth"]
    for i in range(3):
        for r in range(2):
            for c in range(3):
                a = captured["slopes"][c] if model == "GPCM" else 1.
                original = a*(captured["theta"][i]-captured["rater"][r]-captured["criterion"][c]-captured["steps"][c])
                step_truth = truth["Step"]["Shared" if model == "RSM" else f"C{c+1}"]
                mapped = (truth["Person"][f"P{i+1:02d}"]-truth["Rater"][f"R{r+1}"]
                    - truth["Criterion"][f"C{c+1}"]-np.array(list(step_truth.values())))
                mapped *= truth["Slope"][f"C{c+1}"] if model == "GPCM" else 1.
                np.testing.assert_allclose(original, mapped, rtol=1e-14, atol=1e-14)
                assert sum(step_truth.values()) == pytest.approx(0., abs=1e-14)
    if model == "GPCM":
        assert np.prod(list(truth["Slope"].values())) == pytest.approx(1.)
        assert sum(truth["LogSlope"].values()) == pytest.approx(0., abs=1e-14)
    assert bundle["settings"]["person_replication"] == "one_centered_draw_held_fixed"


@pytest.mark.parametrize("slopes", [[1., 0., 1.], [1., -1., 1.], [1., np.nan, 1.], [1., 2.]])
def test_recovery_refuses_invalid_gpcm_slopes(slopes):
    with pytest.raises(ValueError, match="finite positive"):
        app.evaluate_parameter_recovery(model="GPCM", n_criterion=3, reps=1, slopes=np.array(slopes))


def test_generation_failures_remain_in_all_attempts_denominator(monkeypatch):
    def fail(**kwargs):
        raise RuntimeError("deliberate generation failure")

    monkeypatch.setattr(app, "_simulate_recovery_responses", fail)
    bundle = app.evaluate_parameter_recovery(reps=3)
    assert not bundle["available"]
    assert len(bundle["rep_overview"]) == 3
    assert bundle["rep_overview"]["Error"].eq("deliberate generation failure").all()
    counts = bundle["execution_summary"].iloc[0]
    assert counts["PlannedReps"] == counts["AttemptedReps"] == counts["FailedReps"] == 3
    assert counts["CompletedReps"] == counts["ConvergedReps"] == 0


def test_centered_summary_compares_truth_and_estimates_on_the_same_scale():
    result = {"facets":{"person":pd.DataFrame({"Person":["P1","P2"], "Estimate":[2.,4.]})}}
    data = app._recovery_rows_from_fit(result, {"Person":{"P1":3.,"P2":5.}}, rep=1, include_person=True)
    summary = app._recovery_summarize(data).iloc[0]
    assert summary.MeanTruth == summary.MeanEstimate == 0.
    assert summary.MeanTruthRaw == 4.
    assert summary.RawBias == -1.
    assert summary.RMSE == 0.  # Centered target deliberately excludes the common offset.
    assert np.isnan(summary.Bias)
