"""Contracts for prospective strict-JMLE numerical qualification lanes."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from validation import operating_characteristics_strict_jmle as stage_a
from validation import operating_characteristics_strict_jmle_a2 as stage_a2


REPO = Path(__file__).resolve().parents[1]
STAGE_A_PLAN = REPO / "validation" / "operating_characteristics_strict_jmle_plan_20260809.json"
STAGE_A2_PLAN = REPO / "validation" / "operating_characteristics_strict_jmle_a2_plan_20260809.json"


def test_stage_a_never_authorizes_stage_b_when_terminal_gradient_gate_fails():
    plan = stage_a.load_plan(STAGE_A_PLAN)
    input_checks = pd.DataFrame({"Passed": [True] * 7})
    runs = pd.DataFrame({
        "RequestedAnchorRows": [2] * 8 + [0] * 8,
        "FitReturned": [True] * 16,
        "OptimizerSuccess": [True] * 16,
        "TerminalGradientSupNorm": [2e-4] * 16,
        "GradientGatePassed": [False] * 16,
        "OptimizerJacGatePassed": [True] * 16,
        "FiniteDifferenceGatePassed": [True] * 16,
        "AnchorGatePassed": [True] * 16,
        "MatchedLogLikGatePassed": [True] * 16,
        "NumericalQualificationPassed": [False] * 16,
    })

    gates = stage_a.stage_gates(runs, input_checks, plan)

    assert not gates.loc[gates["Gate"].eq("terminal_gradient_sup_norm"), "GatePassed"].iloc[0]
    assert not gates.loc[gates["Gate"].eq("stage_b_authorization"), "GatePassed"].iloc[0]


def test_stage_a2_plan_rejects_post_hoc_gradient_threshold_change(tmp_path):
    plan = json.loads(STAGE_A2_PLAN.read_text(encoding="utf-8"))
    plan["qualification"]["terminal_gradient_sup_norm_threshold"] = 1e-3
    altered = tmp_path / "altered.json"
    altered.write_text(json.dumps(plan), encoding="utf-8")

    with pytest.raises(ValueError, match="threshold must remain"):
        stage_a2.load_plan(altered)


def _a2_gate_rows(lbfgsb_passed: bool, bfgs_passed: bool) -> pd.DataFrame:
    rows = []
    for priority, candidate_id, passed in (
        (1, "lbfgsb_precision", lbfgsb_passed),
        (2, "bfgs_gradient", bfgs_passed),
    ):
        rows.append({
            "CandidateId": candidate_id,
            "CandidatePriority": priority,
            "Gate": "candidate_qualified",
            "Required": 11,
            "Passed": 11 if passed else 10,
            "GatePassed": passed,
        })
    return pd.DataFrame(rows)


def test_stage_a2_selection_prefers_registered_limited_memory_candidate():
    plan = stage_a2.load_plan(STAGE_A2_PLAN)

    selection = stage_a2.selection_table(_a2_gate_rows(True, True), plan)

    original = selection.loc[selection["Decision"].eq("OriginalStageB")].iloc[0]
    revised = selection.loc[selection["Decision"].eq("RevisedStageB2")].iloc[0]
    integration = selection.loc[selection["Decision"].eq("ApplicationCoreChange")].iloc[0]
    assert not bool(original["Authorized"])
    assert bool(revised["Authorized"])
    assert revised["SelectedCandidate"] == "lbfgsb_precision"
    assert not bool(integration["Authorized"])
    assert not selection["PublicSurfaceEnabled"].any()


def test_stage_a2_selection_falls_back_only_to_fully_qualified_bfgs():
    plan = stage_a2.load_plan(STAGE_A2_PLAN)

    selection = stage_a2.selection_table(_a2_gate_rows(False, True), plan)

    revised = selection.loc[selection["Decision"].eq("RevisedStageB2")].iloc[0]
    assert bool(revised["Authorized"])
    assert revised["SelectedCandidate"] == "bfgs_gradient"


def test_stage_a2_selection_blocks_when_no_candidate_qualifies():
    plan = stage_a2.load_plan(STAGE_A2_PLAN)

    selection = stage_a2.selection_table(_a2_gate_rows(False, False), plan)

    revised = selection.loc[selection["Decision"].eq("RevisedStageB2")].iloc[0]
    assert not bool(revised["Authorized"])
    assert revised["SelectedCandidate"] == ""
