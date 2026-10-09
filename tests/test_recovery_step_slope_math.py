"""Recovery targets, covariance expansion and positive slope intervals."""

import numpy as np
import pandas as pd
import pytest

import streamlit_app as app


def fixture(model="GPCM", n_cat=4):
    levels = {"Rater": ["R1", "R2"], "Criterion": ["C1", "C2"]}
    config = dict(model=model, method="MML", n_cat=n_cat, facet_names=list(levels), facet_levels=levels,
        step_facet="Criterion", slope_facet="Criterion",
        facet_specs={f:app.build_facet_constraint(v, centered=f == "Rater") for f,v in levels.items()})
    sizes = app.build_param_sizes(config)
    slices = app._build_param_slices(sizes)
    covariance = np.eye(sum(sizes.values()))*.04
    groups = ["Shared"] if model == "RSM" else levels["Criterion"]
    steps, truth_steps = [], {}
    for group in groups:
        values = np.linspace(-.8, .8, n_cat-1) if n_cat > 2 else [0.]
        truth_steps[group] = {f"Step_{j+1}": 0. for j in range(n_cat-1)}
        for j,value in enumerate(values):
            steps.append(dict(Step=f"Step_{j+1}", Estimate=value,
                              **({"StepFacet":group} if model != "RSM" else {})))
    if n_cat > 3:
        start = slices["steps"].start
        covariance[start, start+1] = covariance[start+1, start] = .01
    slopes = np.array([2., .5])
    if model == "GPCM":
        covariance[slices["log_slopes"], slices["log_slopes"]] = .09
    result = dict(config=config, steps=pd.DataFrame(steps).iloc[::-1],
        slopes=pd.DataFrame(dict(Level=levels["Criterion"], Estimate=slopes, LogEstimate=np.log(slopes))).iloc[::-1],
        params={"log_slopes":np.log(slopes)})
    truth = dict(Step=truth_steps, LogSlope={"C1":0., "C2":0.}, Slope={"C1":1., "C2":1.})
    return result, truth, {"status":"ok", "cov":covariance}


@pytest.mark.parametrize("model", ["RSM", "PCM", "GPCM"])
def test_dependent_step_se_uses_off_diagonal_covariance(model):
    result, truth, covariance = fixture(model)
    data = app._recovery_rows_from_fit(result, truth, rep=1, include_person=False, covariance=covariance)
    group = "Shared" if model == "RSM" else "C1"
    steps = data[data.ParameterType.eq("Step") & data.Facet.eq(group)].set_index("Level")
    assert steps.loc["Step_3", "SE"] == pytest.approx(np.sqrt(.04+.04+2*.01))
    assert steps.loc["Step_3", "CoverageSE"] == pytest.approx(np.sqrt(.10))
    assert steps.loc["Step_3", "EstimateAligned"] == .8
    summary = app._recovery_summarize(data)
    constrained = summary.ParameterType.isin(["Step", "LogSlope"])
    assert summary.loc[constrained, "Bias"].isna().all()
    levels = app._recovery_summarize(data, by_level=True)
    assert levels.loc[levels.ParameterType.eq("Step") & levels.Level.eq("Step_3"), "Bias"].eq(.8).all()


def test_slope_delta_se_and_log_wald_interval_are_consistent_under_exponentiation():
    result, truth, covariance = fixture()
    data = app._recovery_rows_from_fit(result, truth, rep=1, include_person=False, covariance=covariance)
    log = data[data.ParameterType.eq("LogSlope")].set_index("Level")
    slope = data[data.ParameterType.eq("Slope")].set_index("Level")
    np.testing.assert_allclose(log.SE, [.3, .3])
    np.testing.assert_allclose(slope.SE, [.6, .15])
    np.testing.assert_allclose(slope.CoverageCILower, np.exp(log.CoverageCILower))
    np.testing.assert_allclose(slope.CoverageCIUpper, np.exp(log.CoverageCIUpper))
    assert (slope.Covered95 == log.Covered95).all()
    assert (slope.CoverageCILower > 0).all()
    # Arithmetic centering would corrupt the identified positive slopes.
    np.testing.assert_allclose(slope.EstimateAligned, [2., .5])
    assert slope.EstimateAligned.prod() == pytest.approx(1.)


@pytest.mark.parametrize("issue", ["free_sd", "regularized", "singular", "boundary", "jmle"])
def test_step_slope_coverage_refuses_incomplete_or_nonregular_covariance(issue):
    result, truth, covariance = fixture()
    if issue == "free_sd":
        result["config"]["estimate_population_sd"] = True
    elif issue == "regularized":
        covariance["regularized"] = True
    elif issue == "singular":
        covariance["cov"][-1, -1] = 0.
    elif issue == "boundary":
        result["params"]["log_slopes"] = np.array([3., -3.])
        result["slopes"]["LogEstimate"] = result["slopes"].Level.map({"C1":3., "C2":-3.})
        result["slopes"]["Estimate"] = np.exp(result["slopes"].LogEstimate)
    elif issue == "jmle":
        result["config"]["method"] = "JMLE"
    data = app._recovery_rows_from_fit(result, truth, rep=1, include_person=False, covariance=covariance)
    assert set(data.ParameterType) == {"Step", "LogSlope", "Slope"}
    assert data.RecoveryComparable.all()
    assert data.CoverageSE.isna().all()
    assert data.Covered95.isna().all()


def test_binary_fixed_step_is_not_counted_as_successful_parameter_recovery():
    result, truth, covariance = fixture("RSM", 2)
    data = app._recovery_rows_from_fit(result, truth, rep=1, include_person=False, covariance=covariance)
    assert len(data) == 1
    assert data.iloc[0].SE == 0.
    assert data.iloc[0].CoverageStatus == "fixed_by_identification"
    assert not data.iloc[0].RecoveryComparable
    assert app._recovery_summarize(data).empty


def test_gpcm_pipeline_retains_all_step_and_slope_targets_on_native_scale():
    bundle = app.evaluate_parameter_recovery(model="GPCM", fit_method="JMLE", reps=2,
        n_person=12, n_rater=2, n_criterion=3, n_cat=3, maxit=5, seed=2026100420,
        slopes=np.array([.6, 1.7, 2.3]))
    data = bundle["recovery"]
    assert len(data[data.ParameterType.eq("Step")]) == 12
    assert len(data[data.ParameterType.eq("LogSlope")]) == 6
    assert len(data[data.ParameterType.eq("Slope")]) == 6
    for _, rows in data[data.ParameterType.eq("Slope")].groupby("rep"):
        assert rows.Estimate.prod() == pytest.approx(1.)
        assert rows.Truth.prod() == pytest.approx(1.)
    assert data.Covered95.isna().all()  # No invented joint JMLE covariance.


@pytest.mark.parametrize("tied_truth", [False, True])
def test_truth_quartiles_use_truth_and_do_not_reintroduce_centering_bias(tied_truth):
    persons = ["P1", "P2", "P3", "P4"]
    values = np.zeros(4) if tied_truth else np.array([-3., -1., 1., 3.])
    result = {"facets":{"person":pd.DataFrame({"Person":persons, "Estimate":[3., 1., -1., -3.]})}}
    data = app._recovery_rows_from_fit(result, {"Person":dict(zip(persons, values))}, rep=1, include_person=True)
    summary = app._recovery_summarize(data, by_quantile=True)
    if tied_truth:
        assert data.TruthQuantile.eq("Q3").all()
        assert summary.Bias.isna().all()  # The complete centered block is still an identity.
    else:
        assert data.TruthQuantile.tolist() == ["Q1", "Q2", "Q3", "Q4"]
        np.testing.assert_allclose(summary.Bias, [6., 2., -2., -6.])


def test_pcm_generator_varies_step_shape_and_binary_span_cannot_shift_difficulty(monkeypatch):
    captured = []
    simulate = app._simulate_recovery_responses

    def capture(**kwargs):
        captured.append(kwargs)
        return simulate(**kwargs)

    monkeypatch.setattr(app, "_simulate_recovery_responses", capture)
    monkeypatch.setattr(app, "mfrm_estimate", lambda **kwargs: {})
    monkeypatch.setattr(app, "mfrm_diagnostics", lambda *args, **kwargs: {})
    for n_cat, span in ((4, 1.4), (2, 1.4), (2, 3.)):
        bundle = app.evaluate_parameter_recovery(model="PCM", n_cat=n_cat, step_span=span,
            reps=1, n_person=3, n_rater=2, n_criterion=3, seed=2026100421)
        assert bundle["settings"]["generator_version"] == "recovery_v4_centered_step_shape"
    np.testing.assert_allclose(captured[0]["steps"].sum(axis=1), 0., atol=1e-14)
    assert not np.allclose(captured[0]["steps"][0], captured[0]["steps"][-1])
    np.testing.assert_array_equal(captured[1]["steps"], np.zeros((3, 1)))
    np.testing.assert_array_equal(captured[1]["steps"], captured[2]["steps"])
    np.testing.assert_array_equal(captured[1]["criterion"], captured[2]["criterion"])


@pytest.mark.parametrize("model", ["RSM", "GPCM"])
@pytest.mark.parametrize("free_sd", [False, True])
def test_mml_recovery_uses_the_declared_identified_fixed_sd_without_changing_free_sd_start(monkeypatch, model, free_sd):
    captured = {}

    def fit(**kwargs):
        captured.update(kwargs)
        return {}

    monkeypatch.setattr(app, "mfrm_estimate", fit)
    monkeypatch.setattr(app, "mfrm_diagnostics", lambda *args, **kwargs: {})
    bundle = app.evaluate_parameter_recovery(model=model, fit_method="MML", reps=1,
        n_person=3, n_rater=2, n_criterion=3, theta_sd=1.4, slopes=np.array([.6,1.7,2.3]),
        estimate_population_sd=free_sd)
    if free_sd:
        assert captured["estimate_population_sd"] is True
        assert "population_prior_sd" not in captured
        assert bundle["settings"]["fitted_fixed_population_sd"] is None
    else:
        expected = 1.4*((.6*1.7*2.3)**(1/3) if model == "GPCM" else 1.)
        assert captured["population_prior_sd"] == pytest.approx(expected)
        assert bundle["settings"]["fitted_fixed_population_sd"] == pytest.approx(expected)
