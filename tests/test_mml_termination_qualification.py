"""Regression controls for the TAM review; no coverage or general accuracy claim."""
import numpy as np
import pandas as pd
import pytest
from scipy.optimize import OptimizeResult

import streamlit_app as app
from mfrm_app import engine


def small_problem(model="PCM"):
    config = dict(model=model, method="MML", n_cat=3, n_person=2,
                  facet_names=["Criterion"], facet_levels={"Criterion": ["C1", "C2"]},
                  facet_specs={"Criterion": app.build_facet_constraint(["C1", "C2"], centered=True)},
                  facet_signs={"Criterion": -1}, step_facet="Criterion", slope_facet="Criterion")
    sizes = app.build_param_sizes(config)
    point = np.array([.35, -.4, .2] + ([0.] if model == "GPCM" else []))
    idx = dict(person=np.array([0, 0, 1, 1]), score_k=np.array([0, 2, 1, 2]),
               facets={"Criterion": np.array([0, 1, 0, 1])}, step_idx=np.array([0, 1, 0, 1]),
               weight=np.ones(4))
    return point, idx, config, sizes, app.gauss_hermite_normal(31, sd=1.)


@pytest.mark.parametrize("case,reason", [
    ("failed", "optimizer_failed"), ("nan", "nonfinite"),
    ("worse", "objective_increased"), ("infeasible", "infeasible"),
    ("derived_infeasible", "infeasible"),
])
@pytest.mark.parametrize("free_sd", [False, True])
def test_em_rejects_failed_or_invalid_m_step_without_false_convergence(monkeypatch, case, reason, free_sd):
    point, idx, config, sizes, quad = small_problem("GPCM" if "infeasible" in case else "PCM")
    config["estimate_population_sd"] = free_sd
    if case == "derived_infeasible":
        levels = ["C1", "C2", "C3"]
        config["facet_levels"]["Criterion"] = levels
        config["facet_specs"]["Criterion"] = app.build_facet_constraint(levels, centered=True)
        sizes = app.build_param_sizes(config)
        point = np.zeros(sum(sizes.values()))
    calls = []

    def injected(fun, x, *, args, **kwargs):
        value, gradient = fun(x, *args)
        candidate = x.copy()
        if case == "nan":
            candidate[0] = np.nan
        elif case == "worse":
            candidate += 100 * gradient  # climb the negative-Q objective
        elif case == "infeasible":
            candidate[-1] = 4.
        elif case == "derived_infeasible":
            candidate[-2:] = 2.  # each free slope is allowed; the derived slope is -4
        calls.append(candidate)
        return OptimizeResult(x=candidate, fun=value, jac=gradient, success=case != "failed",
                              status=2 if case == "failed" else 0,
                              message="Injected M-step result", nit=0, nfev=1)

    monkeypatch.setattr(engine, "minimize", injected)
    fit = app.mfrm_em_mml(point, idx, config, sizes, quad, maxit=5)
    assert len(calls) == 1
    assert not fit.success and fit.status == 2
    np.testing.assert_array_equal(fit.x, point)
    assert fit.m_step_failure["reason"] == "stat.mml_em.m_step_" + reason
    assert fit.m_step_failure["status"] == (2 if case == "failed" else 0)
    assert "Injected M-step result" in fit.message
    assert np.isfinite(fit.fun)
    if free_sd:
        assert fit.sigma_trace == [1.]
    table = app.build_convergence_summary(fit, config, -fit.raw_fun)
    assert table.MStepFailureReason.iloc[0] == fit.m_step_failure["reason"]
    assert not table.Converged.iloc[0]


def test_later_failure_keeps_last_accepted_point_and_hybrid_records_recovery(monkeypatch):
    point, idx, config, sizes, quad = small_problem()
    minimize = engine.minimize
    accepted = []

    def later_failure(fun, x, **kwargs):
        if not accepted:
            step = minimize(fun, x, **kwargs)
            assert step.success
            accepted.append(step.x.copy())
            return step
        return OptimizeResult(x=x + 100, success=False, status=2,
                              message="later failure", nfev=1)

    monkeypatch.setattr(engine, "minimize", later_failure)
    fit = app.mfrm_em_mml(point, idx, config, sizes, quad, maxit=5, reltol=0.)
    assert not fit.success and fit.m_step_failure["iteration"] == 2
    np.testing.assert_array_equal(fit.x, accepted[0])

    def warmup_failure(fun, x, **kwargs):
        if fun is app._m_step_expected_ll_value_grad:
            return OptimizeResult(x=x.copy(), success=False, status=2,
                                  message="warmup failure", nfev=1)
        return minimize(fun, x, **kwargs)

    monkeypatch.setattr(engine, "minimize", warmup_failure)
    hybrid = app.mfrm_hybrid_mml(point, idx, config, sizes, quad, maxit=200)
    assert hybrid.success
    assert hybrid.em_warmup_failure["message"] == "warmup failure"
    assert "warmup failure" in app.build_convergence_summary(hybrid, config, -hybrid.raw_fun).EMWarmupFailure.iloc[0]


def test_successful_em_keeps_the_marginal_likelihood_and_inference_states_separate():
    point, idx, config, sizes, quad = small_problem()
    fit = app.mfrm_em_mml(point, idx, config, sizes, quad, maxit=200)
    assert fit.success and fit.m_step_failure is None
    assert fit.fun < app.mfrm_loglik_mml(point, idx, config, sizes, quad)
    assert np.min(np.diff(fit.ll_trace)) >= -1e-9
    table = app.build_convergence_summary(fit, config, -fit.raw_fun)
    assert table.Converged.iloc[0]
    assert not table.QuadratureAccuracyQualified.iloc[0]
    assert not table.ScientificInferenceReady.iloc[0]


def test_repeated_qualification_replaces_saved_prose_without_mutating_it():
    original = pd.DataFrame(dict(Estimate=[.2, .4], SE=[.1, .2],
                                 UncertaintyCaution=["Calibration uncertainty excluded.", None]))
    result = dict(config=dict(method="MML"))
    once = app.qualify_measure_uncertainty(original, result)
    # Saved prose can authorize inference; the current hold supplies bounded wording.
    assert once.UncertaintyCaution.str.startswith("Technical terminal-point calculations only").all()
    assert original.UncertaintyCaution.iloc[0] == "Calibration uncertainty excluded."
    assert once.UncertaintyCaution.str.contains("inference is withheld").all()
    pd.testing.assert_frame_equal(once, app.qualify_measure_uncertainty(once, result))
    assert pd.isna(original.UncertaintyCaution.iloc[1])


@pytest.mark.parametrize("missing", [False, True])
def test_reviewed_pcm_cases_retain_estimates_but_withhold_inference(missing):
    data = pd.DataFrame([dict(Person=f"{p:03d}", Rater=f"R{r}", Task=f"T{t}",
                              Criterion=f"C{c}", Score=(p + 2*r + t + c) % 4)
                         for p in range(16) for r in range(2)
                         for t in range(2) for c in range(2)])
    if missing:
        data = data.iloc[1:]
    result = app.mfrm_estimate(data, "Person", ["Rater", "Task", "Criterion"], "Score",
        model="PCM", method="MML", mml_engine="EM", step_facet="Criterion",
        noncenter_facet="Criterion", positive_facets=["Rater"], population_prior_sd=1.1,
        quad_points=31, maxit=100, reltol=1e-6)
    assert result["opt"].success
    assert result["summary"].LogLik.iloc[0] == pytest.approx(
        -195.7026189715229 if missing else -197.3486605636869, abs=1e-10, rel=0)
    assert not result["summary"].InferenceReady.iloc[0]
    assert result["summary"].InferenceReadinessReason.iloc[0] == app.FIXED_SD_MML_INFERENCE_HOLD_REASON
    assert app.build_first_read_guide_rows(result, {})[0]["Status"] == "Do not interpret yet"
    gate = app.build_assignment_sensitivity_preflight(result)["gates"].set_index("Gate")
    assert not gate.loc["Source fit inference readiness", "Passed"]


@pytest.mark.parametrize("free", [False, True])
def test_legacy_success_is_held_in_all_exports_without_mutating_saved_evidence(free):
    result = dict(config=dict(method="MML", model="PCM", estimate_population_sd=free),
                  summary=pd.DataFrame([dict(Method="MML", Converged=True, InferenceReady=True)]))
    expected = app.FREE_SD_MML_INFERENCE_HOLD_REASON if free else app.FIXED_SD_MML_INFERENCE_HOLD_REASON
    frames, _ = app.collect_download_frames(result, {}, {}, pd.DataFrame(), pd.DataFrame(),
                                           public_export_mode=True)
    for summary in (app.qualified_estimation_summary(result), frames["summary"],
                    app.build_result_bundle_frames(result, {})["summary"],
                    app._collect_apa_exportable_tables(result, {})["Estimation summary"]):
        # A legacy table flag without the optimizer record does not establish termination.
        assert not summary.Converged.iloc[0] and summary.SavedConverged.iloc[0]
        assert not summary.InferenceReady.iloc[0]
        assert summary.InferenceReadinessReason.iloc[0] == expected
    assert result["summary"].InferenceReady.iloc[0]
    # A legacy Method column also identifies the estimator if config lacks it.
    assert not app.qualified_estimation_summary(dict(summary=result["summary"])).InferenceReady.iloc[0]
