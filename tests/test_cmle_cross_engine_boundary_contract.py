from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from validation.cmle_cross_engine_boundary import (
    EXPECTED_TO_WORKFLOW,
    json_safe,
    r_ready_frame,
    selected_fixture_cases,
)


ROOT = Path(__file__).resolve().parents[1]
PLAN = ROOT / "validation/cmle_cross_engine_boundary_plan_20260810.json"


def _plan() -> dict[str, object]:
    return json.loads(PLAN.read_text(encoding="utf-8"))


def test_registered_case_accounting_and_order_are_stable() -> None:
    plan = _plan()
    cases = selected_fixture_cases(plan)
    assert [case["CaseId"] for case in cases] == plan["fixture_contract"]["case_ids"]
    assert len(cases) == 13
    assert sum(case["ExpectedStatus"] == "interior_finite_cmle_supported" for case in cases) == 5
    assert sum(case["ExpectedStatus"] == "boundary_no_finite_cmle" for case in cases) == 7
    assert sum(case["ExpectedStatus"] == "structural_nonidentification" for case in cases) == 1


def test_r_mapping_preserves_scores_and_uses_event_only_as_response_unit() -> None:
    cases = selected_fixture_cases(_plan())
    for case in cases:
        mapped = r_ready_frame(case)
        assert len(mapped) == len(case["frame"])
        assert mapped["Score"].tolist() == case["frame"]["Score"].astype(int).tolist()
        assert not mapped.duplicated(["CaseId", "Person", "Rater", "Unit"]).any()
        if "Event" in case["frame"].columns:
            assert mapped["Unit"].tolist() == case["frame"]["Event"].astype(str).tolist()
        else:
            assert set(mapped["Unit"]) == {"Response"}


def test_plan_forbids_false_estimator_equivalence_and_rounded_decisions() -> None:
    plan = _plan()
    roles = plan["engine_roles"]
    assert roles["immer"]["role"].startswith("primary same-estimand")
    for engine in ("mfrmr", "TAM", "sirt"):
        assert "different-estimand" in roles[engine]["role"]
        assert "do not score equality" in roles[engine]["interpretation"]
    observation = plan["registered_observations"]["binary64_sensitivity"]
    assert "Never base the primary result on printed or rounded" in observation
    assert EXPECTED_TO_WORKFLOW["boundary_no_finite_cmle"] == "finite_mle_boundary"


def test_decision_json_uses_null_for_nonfinite_descriptive_cells() -> None:
    assert json_safe({"median": np.nan, "maximum": np.float64(np.inf)}) == {
        "median": None,
        "maximum": None,
    }
