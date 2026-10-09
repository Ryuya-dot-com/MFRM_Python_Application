"""Closed-form controls for conditional information and delta-method validity."""

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

import streamlit_app as app


def test_conditional_se_uses_gpcm_information_and_reduces_to_rasch():
    # At eta=0, zero steps give uniform probabilities on {0,1,2}:
    # Var(X)=2/3. Slopes 2 and 1/2 satisfy geometric-mean-one identification.
    obs = pd.DataFrame({"Person": ["P1", "P1"], "Criterion": ["C1", "C2"],
                        "Var": [2/3, 2/3], "Weight": [1., 2.],
                        "ScoreInformation": [8/3, 1/6]})
    table = app.calc_facet_se(obs, ["Person", "Criterion"]).set_index(["Facet", "Level"])
    assert table.loc[("Person", "P1"), "SE"] == pytest.approx(1/np.sqrt(3))
    assert table.loc[("Criterion", "C1"), "SE"] == pytest.approx(np.sqrt(3/8))
    assert table.loc[("Criterion", "C2"), "SE"] == pytest.approx(np.sqrt(3))
    rasch = obs.assign(ScoreInformation=obs.Var)
    pd.testing.assert_frame_equal(app.calc_facet_se(rasch, ["Person"]),
                                  app.calc_facet_se(rasch.drop(columns="ScoreInformation"), ["Person"]))


@pytest.mark.parametrize("derivative", ["unavailable", "finite", "zero"])
def test_all_delta_method_outputs_require_every_derivative(monkeypatch, derivative):
    result = {"config": {"method": "MML", "model": "GPCM", "n_cat": 3,
                         "facet_names": ["Criterion"], "step_facet": "Criterion"},
              "prep": {"levels": {"Criterion": ["C1"]}, "rating_min": 0, "rating_max": 2},
              "opt": SimpleNamespace(x=np.zeros(2))}
    covariance = {"status": "ok", "cov": np.array([[1., .25], [.25, 4.]])}

    def value(_result, p, *_args, **_kwargs):
        if derivative == "unavailable" and p[1] != 0:
            return np.nan
        if derivative == "zero":
            return .7
        return .7 + p[0] + 2*p[1]

    for name in ("_compute_measure_value_from_par", "_compute_step_value_from_par",
                 "_compute_fair_average_value_from_par"):
        monkeypatch.setattr(app, name, value)
    measures = app.compute_mml_structural_measure_se_table(result, covariance)
    steps = app.compute_mml_step_se_table(result, covariance)
    fair = app.add_gpcm_fair_average_delta_se(
        {"Criterion": pd.DataFrame({"Level": ["C1"], "FairM": [.7], "FairZ": [.7]})},
        result, covariance)["Criterion"]
    for frame, column, status in ((measures, "StructuralSE", "StructuralSE_Status"),
                                  (steps, "SE", "SE_Status"),
                                  (fair, "FairMSE", "FairM_SE_Status"),
                                  (fair, "FairZSE", "FairZ_SE_Status")):
        if derivative == "unavailable":
            assert frame[column].isna().all()
            assert not frame[status].isin(["ok", "regularized"]).any()
        else:
            # (1,2) @ Cov @ (1,2).T = 18, including the off-diagonal terms.
            expected = 0. if derivative == "zero" else np.sqrt(18)
            np.testing.assert_allclose(frame[column], expected, rtol=1e-10)


def test_gpcm_corrected_se_reaches_measure_ci_and_facets_report():
    data = pd.DataFrame([dict(Person=f"P{p}", Rater=f"R{r}", Criterion=f"C{c}",
                              Score=(p+r+c)%3)
                         for p in range(6) for r in range(2) for c in range(2)])
    result = app.mfrm_estimate(data, "Person", ["Rater", "Criterion"], "Score",
                               model="GPCM", method="JMLE", step_facet="Criterion", maxit=1)
    # A deterministic diagnostic point, not a claimed MLE or recovery study.
    point = np.zeros_like(result["opt"].x)
    sizes = app.build_param_sizes(result["config"])
    point[app._build_param_slices(sizes)["log_slopes"]] = np.log(2.)
    result["opt"].x = point
    result["params"] = app.expand_params(point, sizes, result["config"])
    result["facets"]["person"]["Estimate"] = 0.
    result["facets"]["others"]["Estimate"] = 0.
    diagnostics = app.mfrm_diagnostics(result, compute_pca=False, interaction_pairs=[])
    measures = diagnostics["measures"].set_index(["Facet", "Level"])
    # Each person has two observations at each slope: I=17/3.
    persons = measures.loc["Person"]
    np.testing.assert_allclose(persons.SE, np.sqrt(3/17))
    np.testing.assert_allclose(persons.CI_Upper-persons.Estimate,
                               app._norm.ppf(.975)*np.sqrt(3/17))
    assert persons.SE_Method.str.contains("a^2", regex=False).all()
    facets = app.calc_facets_report_tbls(result, diagnostics, uscale=2.)
    np.testing.assert_allclose(facets["Person"].ModelSE, 2*np.sqrt(3/17))


@pytest.mark.parametrize("scope", [
    {"estimate_population_sd": True},
    {"facet_regularization_enabled": True},
    {"facet_regularization": {"enabled": True}},
])
def test_incomplete_covariance_scope_cannot_supply_structural_intervals(scope):
    result = {"config": {"method": "MML", "model": "GPCM", "n_cat": 3,
                         "facet_names": ["Criterion"], "step_facet": "Criterion", **scope},
              "prep": {"levels": {"Criterion": ["C1"]}},
              "opt": SimpleNamespace(x=np.zeros(2))}
    # An old full-rank covariance cannot override the model's current scope.
    legacy = {"status": "ok", "cov": np.eye(2), "hessian": np.eye(2), "rank": 2}
    for covariance in (app.compute_mml_parameter_covariance(result),
                       app._mml_covariance_for_diagnostics(result, max_params=0)):
        assert covariance["status"] == "unsupported_covariance_scope"
        assert covariance["cov"] is None
    audit = app.build_mml_covariance_audit(legacy, result)
    assert audit.iloc[0].QualificationStatus == "WITHHELD"
    old_audit = pd.DataFrame([{"Status": "ok", "Rank": 2, "ParamCount": 2}])
    outputs = app.build_output_qualification_table(
        result, {"uncertainty": {"covariance_audit": old_audit}}).set_index("OutputID")
    assert outputs.loc["mml.structural_covariance", "QualificationStatus"] == "WITHHELD"
    measures = app.compute_mml_structural_measure_se_table(result, legacy)
    steps = app.compute_mml_step_se_table(result, legacy)
    fair = app.add_gpcm_fair_average_delta_se(
        {"Criterion": pd.DataFrame({"Level": ["C1"], "FairM": [1.]})},
        result, legacy)["Criterion"]
    for frame, column, status in ((measures, "StructuralSE", "StructuralSE_Status"),
                                  (steps, "SE", "SE_Status"),
                                  (fair, "FairMSE", "FairM_SE_Status")):
        assert frame[column].isna().all()
        assert frame[status].eq("unsupported_covariance_scope").all()
    assert legacy["status"] == "ok"  # Saved optimizer evidence is unchanged.


def test_missing_hessian_is_not_evidence_of_rank_deficiency():
    audit = app.build_mml_covariance_audit(
        {"status": "skipped_large_parameter_count", "cov": None},
        {"config": {"method": "MML"}, "opt": SimpleNamespace(x=np.zeros(2))})
    assert audit.iloc[0].ReasonCode == "stat.stat_003.covariance_unavailable"
    assert audit.iloc[0].NextGate == "G2"


@pytest.mark.parametrize("options", [
    {"estimate_population_sd": True},
    {"facet_regularization": {"facets": ["Rater"], "sd": 1.}},
])
def test_scope_hold_reaches_fitted_diagnostics_and_step_table(options):
    result = app.mfrm_estimate(
        app._make_self_test_gradient_data(), "Person", ["Rater", "Task"], "Score",
        model="GPCM", method="MML", mml_engine="EM", quad_points=5, maxit=2, **options)
    diagnostics = app.mfrm_diagnostics(result, compute_pca=False, interaction_pairs=[])
    measures = diagnostics["measures"]
    assert measures.SE_CovarianceStatus.eq("unsupported_covariance_scope").all()
    assert measures.loc[measures.Facet != "Person", "SE_Status"].eq("conditional_approximation").all()
    assert measures.loc[measures.Facet == "Person", "SE_Status"].eq("posterior_eap_sd").all()
    assert result["steps"].SE.isna().all()
    assert result["steps"].SE_Status.eq("unsupported_covariance_scope").all()
