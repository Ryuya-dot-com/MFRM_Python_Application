"""Contracts for JMLE Person score-boundary reporting."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

import streamlit_app as app
from validation import jmle_extreme_score_integration_check as integration


def _boundary_frame() -> pd.DataFrame:
    return pd.DataFrame({
        "Person": ["P1", "P1", "P2", "P2"],
        "Rater": ["R1", "R2", "R1", "R2"],
        "Task": ["T1", "T1", "T1", "T1"],
        "Criterion": ["C1", "C1", "C1", "C1"],
        "Score": [0, 0, 0, 1],
    })


def _audit(*, method: str = "JMLE", estimates=(0.25, 1000.0), fit_ready=True, frame=None):
    prep = app.prepare_mfrm_data(
        _boundary_frame() if frame is None else frame,
        person_col="Person",
        facet_cols=["Rater", "Task", "Criterion"],
        score_col="Score",
        rating_min=0,
        rating_max=1,
        keep_original=True,
    )
    config = {"method": method, "n_person": 2, "n_cat": 2}
    return app.build_jmle_person_boundary_audit(
        prep,
        config,
        estimates=np.asarray(estimates, dtype=float) if method == "JMLE" else None,
        fit_inference_ready=fit_ready,
    )


def test_boundary_uses_exact_scores_not_terminal_estimate_magnitude():
    audit = _audit()
    persons = audit["persons"].set_index("Person")
    summary = audit["summary"].iloc[0]

    assert audit["scope"] == "jmle_person_score_boundary_only"
    assert audit["status"] == "extreme_persons_present"
    assert int(summary["ExtremePersons"]) == 1
    assert persons.loc["P1", "ExtremeScoreDirection"] == "all_minimum"
    assert not bool(persons.loc["P1", "FiniteJMLEEstimate"])
    assert pd.isna(persons.loc["P1", "ReportableEstimate"])
    assert "not_finite_person_mle" in persons.loc["P1", "EstimateRole"]
    assert bool(persons.loc["P2", "FiniteJMLEEstimate"])
    assert float(persons.loc["P2", "Estimate"]) == 1000.0
    assert float(persons.loc["P2", "ReportableEstimate"]) == 1000.0


def test_interior_person_still_inherits_fit_level_readiness():
    audit = _audit(fit_ready=False)
    persons = audit["persons"].set_index("Person")

    assert bool(persons.loc["P2", "FiniteJMLEEstimate"])
    assert not bool(persons.loc["P2", "PersonInferenceReady"])
    assert pd.isna(persons.loc["P2", "ReportableEstimate"])
    assert persons.loc["P2", "EstimateRole"] == "finite_boundary_status_but_fit_not_inference_ready"


def test_all_maximum_pattern_is_classified_without_estimate_cutoff():
    frame = _boundary_frame()
    frame.loc[frame["Person"].eq("P1"), "Score"] = 1
    audit = _audit(frame=frame, estimates=(-1000.0, 0.25))
    persons = audit["persons"].set_index("Person")

    assert persons.loc["P1", "ExtremeScoreDirection"] == "all_maximum"
    assert not bool(persons.loc["P1", "FiniteJMLEEstimate"])
    assert pd.isna(persons.loc["P1", "ReportableEstimate"])


def test_mml_is_outside_unbounded_jmle_boundary_contract():
    audit = _audit(method="MML")

    assert audit["status"] == "not_applicable_mml"
    assert audit["persons"].empty
    assert pd.isna(audit["summary"].iloc[0]["PersonMeasureInferenceReady"])


def test_registered_integration_plan_forbids_optimizer_or_estimate_changes(tmp_path):
    plan_path = Path(
        "validation/jmle_extreme_score_integration_plan_20260809.json"
    ).resolve()
    plan = integration.load_plan(plan_path)

    assert not plan["authorized_changes"]["change_optimizer_controls"]
    assert not plan["authorized_changes"]["change_parameter_estimates"]
    assert not plan["authorized_changes"]["authorize_precision_polish_in_application_core"]

    invalid = dict(plan)
    invalid["authorized_changes"] = dict(plan["authorized_changes"])
    invalid["authorized_changes"]["change_parameter_estimates"] = True
    invalid_path = tmp_path / "invalid_plan.json"
    invalid_path.write_text(json.dumps(invalid), encoding="utf-8")

    try:
        integration.load_plan(invalid_path)
    except ValueError as exc:
        assert "change_parameter_estimates" in str(exc)
    else:
        raise AssertionError("unsafe extreme-score plan was accepted")


def test_fitted_extreme_person_is_flagged_and_exported_without_replacing_estimate():
    result = app.mfrm_estimate(
        _boundary_frame(),
        person_col="Person",
        facet_cols=["Rater", "Task", "Criterion"],
        score_col="Score",
        rating_min=0,
        rating_max=1,
        model="RSM",
        method="JMLE",
        noncenter_facet="Person",
        maxit=80,
        reltol=1e-6,
        keep_original=True,
    )
    summary = result["summary"].iloc[0]
    persons = result["facets"]["person"].set_index("Person")
    optimizer_theta = result["params"]["theta"]
    quick_frames = app.build_result_bundle_frames(result, {})
    diagnostics = app.mfrm_diagnostics(result, compute_pca=False)
    measure_persons = diagnostics["measures"].loc[
        diagnostics["measures"]["Facet"].eq("Person")
    ].set_index("Level")

    assert int(summary["ExtremeJMLEPersons"]) == 1
    assert int(summary["FiniteJMLEPersonEstimates"]) == 1
    assert not bool(summary["PersonMeasureInferenceReady"])
    assert np.array_equal(
        persons.loc[["P1", "P2"], "Estimate"].to_numpy(),
        optimizer_theta,
    )
    assert pd.isna(persons.loc["P1", "ReportableEstimate"])
    assert not bool(measure_persons.loc["P1", "FiniteJMLEEstimate"])
    assert pd.isna(measure_persons.loc["P1", "ReportableEstimate"])
    assert {
        "jmle_person_boundary_summary",
        "jmle_person_boundary_persons",
    }.issubset(quick_frames)


def test_japanese_convergence_view_explains_extreme_person_boundary():
    audit = _audit()

    class Recorder:
        def __init__(self):
            self.session_state = {"lang": "ja"}
            self.messages = []

        def __getattr__(self, name):
            if name in {"warning", "error", "info", "success", "caption", "markdown", "subheader"}:
                return lambda message="", *args, **kwargs: self.messages.append(str(message))
            if name == "dataframe":
                return lambda *args, **kwargs: None
            if name == "expander":
                def expander(*args, **kwargs):
                    class Context:
                        def __enter__(self):
                            return None

                        def __exit__(self, *exc):
                            return False

                    return Context()

                return expander
            return getattr(app.st, name)

    recorder = Recorder()
    result = {
        "config": {"method": "JMLE"},
        "convergence": pd.DataFrame([{
            "Converged": True,
            "RequestedMmlEngine": "",
            "ResolvedMmlEngine": "",
            "GradientNorm": 0.0,
            "ElapsedSeconds": 0.01,
        }]),
        "person_boundary": audit,
    }
    saved = app.st
    app.st = recorder
    try:
        app.show_convergence_section(result)
    finally:
        app.st = saved

    assert any("有限なJMLE Person最尤推定値は存在しません" in message for message in recorder.messages)
    assert any("表示丸めの閾値は使いません" in message for message in recorder.messages)
