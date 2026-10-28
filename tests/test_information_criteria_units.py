"""Independent-cluster BIC bookkeeping, legacy views, and applicability guards."""
import numpy as np
import pandas as pd
import pytest

import streamlit_app as app


def legacy_result(method="MML", persons=12, n=72, k=5, ll=-123.):
    summary = pd.DataFrame([dict(Model="RSM", Method=method, Persons=persons,
        N=n, KParams=k, LogLik=ll, AIC=-2*ll+2*k, BIC=-2*ll+k*np.log(n))])
    return dict(config=dict(method=method, model="RSM", n_cat=3), summary=summary,
                prep=dict(data=pd.DataFrame(dict(Person=["P1", "P2"], Score=[0, 1]))))


@pytest.mark.parametrize("responses", [24, 72, 120])
def test_mml_bic_penalty_counts_independent_persons_not_rows(responses):
    result = legacy_result(n=responses)
    before = result["summary"].copy(deep=True)
    report = app.qualified_estimation_summary(result)
    row = report.iloc[0]
    assert row["BIC"] == pytest.approx(246+5*np.log(12))
    assert row["BICSampleSize"] == 12
    assert row["N"] == responses
    assert row["SavedBIC"] == pytest.approx(246+5*np.log(responses))
    assert row["BICPerObs"] == pytest.approx(row["BIC"]/responses)
    assert row["BICSampleSizeBasis"] == "independent_person_marginals"
    pd.testing.assert_frame_equal(result["summary"], before)
    pd.testing.assert_frame_equal(app.qualified_estimation_summary(dict(result, summary=report)), report)


def test_mml_missing_person_count_does_not_silently_use_response_count():
    result = legacy_result()
    result["summary"] = result["summary"].drop(columns="Persons")
    report = app.qualified_estimation_summary(result)
    assert report["BIC"].isna().all() and report["BICSampleSize"].isna().all()
    ic = app.build_likelihood_information_criteria(result).set_index("Metric")
    assert np.isnan(ic.loc["BIC", "Value"])
    assert np.isfinite(report.iloc[0]["SavedBIC"])


@pytest.mark.parametrize("persons", [0, -1, 12.5, np.nan])
def test_mml_invalid_person_count_is_unavailable(persons):
    row = app.qualified_estimation_summary(legacy_result(persons=persons)).iloc[0]
    assert np.isnan(row["BIC"]) and np.isnan(row["BICSampleSize"])
    assert row["BICSampleSizeBasis"] == "unavailable_person_count"


def test_jmle_row_count_is_explicitly_unqualified_for_incidental_parameters():
    result = legacy_result(method="JMLE")
    report = app.qualified_estimation_summary(result).iloc[0]
    assert report["BIC"] == pytest.approx(246+5*np.log(72))
    assert report["BICSampleSize"] == 72
    assert "incidental person parameters" in report["BICScope"]


def test_likelihood_table_and_history_view_use_same_mml_unit():
    result = legacy_result()
    ic = app.build_likelihood_information_criteria(result).set_index("Metric")
    history = app.extract_run_likelihood_row(dict(result=result))
    assert ic.loc["BIC", "Value"] == history["BIC"] == pytest.approx(246+5*np.log(12))
    assert ic.loc["BICSampleSize", "Value"] == history["BICSampleSize"] == 12
    assert ic.loc["N", "Value"] == history["N"] == 72
    assert "log(BICSampleSize)" in ic.loc["BIC", "Formula"]


def test_compact_history_without_person_metadata_retains_only_saved_mml_bic():
    entry = dict(method="MML", model="RSM", n_obs=72, k_params=5, loglik=-123., bic=267.3)
    report = app.build_run_history_likelihood_table([entry]).iloc[0]
    assert np.isnan(report["BIC"]) and np.isnan(report["BICPerObs"])
    assert report["SavedBIC"] == 267.3 and entry["bic"] == 267.3
    assert not report["ComparableIC"]


def test_comparison_rebuilds_legacy_mml_bic_for_every_candidate(monkeypatch):
    from types import SimpleNamespace
    original = legacy_result()
    original['config']['facet_names'] = ['Criterion']
    original['prep']['data']['Criterion'] = ['C1', 'C2']
    original["opt"] = SimpleNamespace(success=True)
    def refit(res, model, **kwargs):
        result = legacy_result(k={"PCM": 7, "GPCM": 9}[model])
        result["config"]["model"] = model
        result["opt"] = SimpleNamespace(success=True)
        return result
    monkeypatch.setattr(app, "_refit_candidate_model", refit)
    bundle = app.compute_model_choice_comparison(original)
    comp = bundle["comparison"]
    np.testing.assert_allclose(comp["BIC"], 246+comp["KParams"]*np.log(12))
    np.testing.assert_allclose(comp["DeltaBIC"], [0., 2*np.log(12), 4*np.log(12)])
    assert (comp["BICSampleSize"] == 12).all()
    assert bundle["recommendation"]["tier"] == "withheld"
    assert not bundle["qualification"]["PublicConclusionAllowed"].any()


def test_exports_rebuild_stale_cached_likelihood_table_without_changing_it():
    result = legacy_result()
    result["config"]["facet_names"] = []
    stale = pd.DataFrame([dict(Metric="BIC", Value=999.)])
    before = stale.copy(deep=True)
    result["likelihood_information"] = stale
    frames, _ = app.collect_download_frames(result, {}, {}, pd.DataFrame(), pd.DataFrame(), public_export_mode=True)
    tables = [app.build_result_bundle_frames(result, {})["likelihood_information_criteria"],
              frames["likelihood_information_criteria"],
              app._collect_apa_exportable_tables(result, {})["Likelihood / information criteria"]]
    for table in tables:
        assert table.set_index("Metric").loc["BIC", "Value"] == pytest.approx(246+5*np.log(12))
        assert not table["PublicConclusionAllowed"].any()
    pd.testing.assert_frame_equal(result["likelihood_information"], before)


def test_model_choice_cache_distinguishes_population_models_and_current_ic_policy():
    result = legacy_result()
    fixed = app._model_choice_cache_key(result)
    result["config"]["population_prior_sd"] = 1.7
    other_fixed = app._model_choice_cache_key(result)
    result["config"]["estimate_population_sd"] = True
    free = app._model_choice_cache_key(result)
    assert len({fixed, other_fixed, free}) == 3
    assert fixed.startswith("person_bic_no_aicc_v1|")


def test_fresh_mml_fit_uses_person_bic_and_keeps_response_normalization():
    data = pd.DataFrame([dict(Person=f"P{p}", Rater=f"R{r}", Criterion=f"C{c}", Score=(p+r+c)%3)
                         for p in range(6) for r in range(2) for c in range(2)])
    fit = app.mfrm_estimate(data=data, person_col="Person", facet_cols=["Rater", "Criterion"],
                           score_col="Score", rating_min=0, rating_max=2, model="RSM",
                           method="MML", quad_points=11, maxit=3)
    row = fit["summary"].iloc[0]
    assert row["N"] == 24 and row["Persons"] == row["BICSampleSize"] == 6
    assert row["BIC"] == pytest.approx(-2*row["LogLik"]+row["KParams"]*np.log(6))
    assert row["BICPerObs"] == pytest.approx(row["BIC"]/24)
    assert "SavedBIC" not in fit["summary"]
