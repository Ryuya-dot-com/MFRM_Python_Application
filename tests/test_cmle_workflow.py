"""Structured pre-optimization CMLE workflow contracts."""

from __future__ import annotations

import pandas as pd

import mfrm_app.cmle_workflow as workflow


def _interior_frame() -> pd.DataFrame:
    scores = {
        "P1": [0, 1, 1, 2],
        "P2": [1, 0, 2, 1],
        "P3": [2, 1, 1, 0],
        "P4": [1, 2, 0, 1],
        "P5": [0, 2, 1, 1],
        "P6": [2, 0, 1, 1],
    }
    units = [("R1", "C1"), ("R1", "C2"), ("R2", "C1"), ("R2", "C2")]
    rows = []
    for person, values in scores.items():
        rows.extend(
            (person, rater, criterion, value)
            for (rater, criterion), value in zip(units, values, strict=True)
        )
    return pd.DataFrame(rows, columns=["Person", "Rater", "Criterion", "Score"])


def _boundary_frame(persons: int = 20) -> pd.DataFrame:
    rows = []
    for index in range(persons):
        rows.extend(
            [(f"P{index:03d}", "R1", 1), (f"P{index:03d}", "R2", 0)]
        )
    return pd.DataFrame(rows, columns=["Person", "Rater", "Score"])


def _structural_frame() -> pd.DataFrame:
    rows = []
    patterns = ([0, 1, 2], [1, 2, 0], [2, 0, 1], [0, 2, 1])
    for person_index in range(12):
        rater = "R1" if person_index < 6 else "R2"
        for criterion, score in zip(
            ["C1", "C2", "C3"], patterns[person_index % 4], strict=True
        ):
            rows.append((f"P{person_index}", rater, criterion, score))
    return pd.DataFrame(rows, columns=["Person", "Rater", "Criterion", "Score"])


def _run(frame: pd.DataFrame, *, binary: bool = False):
    return workflow.run_cmle_calibration_workflow(
        frame,
        person_col="Person",
        facet_cols=["Rater"] if binary else ["Rater", "Criterion"],
        score_col="Score",
        rating_min=0,
        rating_max=1 if binary else 2,
        model="RSM",
        gtol=1e-8,
        maxiter=800,
    )


def _assert_stable_schema(result):
    assert set(result) == {
        "schema_version",
        "summary",
        "stages",
        "availability",
        "design_audit",
        "existence_audit",
        "fit",
        "error",
    }
    assert result["schema_version"] == workflow.CMLE_WORKFLOW_SCHEMA_VERSION
    assert len(result["summary"]) == 1
    assert result["stages"]["Stage"].tolist() == [
        "input_and_design",
        "finite_mle_existence",
        "structural_optimization",
        "downstream_person_scoring",
    ]
    assert result["stages"]["Order"].tolist() == [1, 2, 3, 4]
    assert result["stages"]["HeadlineEn"].astype(str).str.len().gt(0).all()
    assert result["stages"]["HeadlineJa"].astype(str).str.len().gt(0).all()
    summary = result["summary"].iloc[0]
    for column in (
        "HeadlineEn",
        "HeadlineJa",
        "DetailEn",
        "DetailJa",
        "NextActionEn",
        "NextActionJa",
    ):
        assert len(str(summary[column])) > 0
    assert result["availability"]["Artifact"].tolist() == [
        "workflow_diagnostics",
        "boundary_direction",
        "structural_calibration",
        "downstream_person_scoring",
    ]


def test_interior_workflow_reaches_ready_calibration():
    result = _run(_interior_frame())
    _assert_stable_schema(result)
    summary = result["summary"].iloc[0]
    assert summary["WorkflowStatus"] == "calibration_ready"
    assert bool(summary["Ready"])
    assert bool(summary["OptimizationAttempted"])
    assert bool(summary["PersonScoringAvailable"])
    assert not bool(summary["PersonScoringPerformed"])
    assert result["stages"]["Status"].tolist() == [
        "pass",
        "pass",
        "pass",
        "available",
    ]
    assert result["fit"] is not None
    assert bool(result["fit"]["summary"].iloc[0]["InferenceReady"])


def test_boundary_workflow_never_calls_optimizer(monkeypatch):
    def forbidden(*_args, **_kwargs):
        raise AssertionError("fit_cmle must not be called on a boundary")

    monkeypatch.setattr(workflow, "fit_cmle", forbidden)
    result = _run(_boundary_frame(), binary=True)
    _assert_stable_schema(result)
    summary = result["summary"].iloc[0]
    assert summary["WorkflowStatus"] == "finite_mle_boundary"
    assert bool(summary["StoppedBeforeOptimization"])
    assert not bool(summary["OptimizationAttempted"])
    assert result["fit"] is None
    assert result["existence_audit"]["direction"] is not None
    assert result["stages"]["Status"].tolist() == [
        "pass",
        "block",
        "not_run",
        "not_run",
    ]


def test_structural_block_never_calls_oracle_or_optimizer(monkeypatch):
    def forbidden(*_args, **_kwargs):
        raise AssertionError("later stage must not be called")

    monkeypatch.setattr(workflow, "audit_cmle_finite_mle_oracle", forbidden)
    monkeypatch.setattr(workflow, "fit_cmle", forbidden)
    result = _run(_structural_frame())
    _assert_stable_schema(result)
    summary = result["summary"].iloc[0]
    assert summary["WorkflowStatus"] == "design_not_identified"
    assert bool(summary["StoppedBeforeOptimization"])
    assert result["existence_audit"] is None
    assert result["fit"] is None
    assert result["stages"]["Status"].tolist() == [
        "block",
        "not_run",
        "not_run",
        "not_run",
    ]


def test_unavailable_existence_never_calls_optimizer(monkeypatch):
    def unavailable(_design):
        return {
            "summary": pd.DataFrame(
                [
                    {
                        "Status": "tolerance_unstable",
                        "Reason": "injected_test_instability",
                        "ExistenceQualified": False,
                    }
                ]
            ),
            "direction": None,
        }

    def forbidden(*_args, **_kwargs):
        raise AssertionError("fit_cmle must not be called when existence is unavailable")

    monkeypatch.setattr(workflow, "audit_cmle_finite_mle_oracle", unavailable)
    monkeypatch.setattr(workflow, "fit_cmle", forbidden)
    result = _run(_interior_frame())
    _assert_stable_schema(result)
    assert result["summary"].iloc[0]["WorkflowStatus"] == "finite_mle_unavailable"
    assert result["stages"]["Status"].tolist() == [
        "pass",
        "review",
        "not_run",
        "not_run",
    ]


def test_input_error_returns_structured_result():
    result = _run(_interior_frame().drop(columns=["Score"]))
    _assert_stable_schema(result)
    summary = result["summary"].iloc[0]
    assert summary["WorkflowStatus"] == "input_invalid"
    assert "ValueError" in summary["Error"]
    assert result["design_audit"] is None
    assert result["stages"]["Status"].tolist() == [
        "block",
        "not_run",
        "not_run",
        "not_run",
    ]


def test_fit_exception_returns_optimization_not_ready(monkeypatch):
    def failed(*_args, **_kwargs):
        raise FloatingPointError("injected fit failure")

    monkeypatch.setattr(workflow, "fit_cmle", failed)
    result = _run(_interior_frame())
    _assert_stable_schema(result)
    summary = result["summary"].iloc[0]
    assert summary["WorkflowStatus"] == "optimization_not_ready"
    assert bool(summary["OptimizationAttempted"])
    assert not bool(summary["StoppedBeforeOptimization"])
    assert not bool(summary["FitReturned"])
    assert "injected fit failure" in summary["Error"]
    assert result["stages"]["Status"].tolist() == [
        "pass",
        "pass",
        "block",
        "not_run",
    ]


def test_fit_nonready_retains_technical_fit(monkeypatch):
    technical = {
        "summary": pd.DataFrame(
            [
                {
                    "InferenceReady": False,
                    "ReadinessReasons": "injected_nonready",
                }
            ]
        )
    }
    monkeypatch.setattr(workflow, "fit_cmle", lambda *_args, **_kwargs: technical)
    result = _run(_interior_frame())
    _assert_stable_schema(result)
    summary = result["summary"].iloc[0]
    assert summary["WorkflowStatus"] == "optimization_not_ready"
    assert bool(summary["OptimizationAttempted"])
    assert not bool(summary["FitInferenceReady"])
    assert result["fit"] is technical
