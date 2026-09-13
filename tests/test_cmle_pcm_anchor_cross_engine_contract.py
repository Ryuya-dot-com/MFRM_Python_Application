from __future__ import annotations

import json
from pathlib import Path

from validation.cmle_pcm_anchor_cross_engine import registered_cases


ROOT = Path(__file__).resolve().parents[1]
PLAN = ROOT / "validation/cmle_pcm_anchor_cross_engine_plan_20260810.json"
REMEDIATION = (
    ROOT
    / "validation/cmle_pcm_anchor_cross_engine_identity_remediation_plan_20260810.json"
)
REMEDIATION_AMENDMENT = (
    ROOT
    / "validation/cmle_pcm_anchor_cross_engine_identity_remediation_amendment_20260810.json"
)


def _plan() -> dict[str, object]:
    return json.loads(PLAN.read_text(encoding="utf-8"))


def test_registered_pcm_anchor_matrix_is_stable() -> None:
    plan = _plan()
    cases = registered_cases(plan)
    assert len(cases) == 12
    assert [case["CaseId"] for case in cases] == [
        item["case_id"] for item in plan["fixture_contract"]["cases"]
    ]
    assert sum(
        case["ExpectedStatus"] == "interior_finite_cmle_supported"
        for case in cases
    ) == 6
    assert sum(
        case["ExpectedStatus"] == "boundary_no_finite_cmle" for case in cases
    ) == 6


def test_pcm_anchor_inputs_preserve_response_rows_and_scope() -> None:
    for case in registered_cases(_plan()):
        frame = case["frame"]
        args = case["prepare_args"]
        assert list(frame.columns) == ["Person", "Rater", "Criterion", "Score"]
        assert args["model"] == "PCM"
        assert args["step_facet"] == "Criterion"
        assert args["rating_min"] == 0 and args["rating_max"] == 2
        anchors = args["hard_anchors"]
        if anchors is not None:
            assert set(anchors["ParameterType"]) == {"Facet"}
            assert set(anchors["Facet"]) <= {"Rater", "Criterion"}


def test_plan_separates_exact_parity_from_secondary_estimators() -> None:
    plan = _plan()
    matched = plan["matched_python_immer_contract"]
    secondary = plan["secondary_engine_contract"]
    assert "b_const" in matched["fixed_offsets"]
    assert "Never score equality to CMLE" in secondary["mfrmr"]
    assert "typed unsupported" in secondary["TAM"]
    assert "typed unsupported" in secondary["sirt"]
    assert plan["success_gates"]["public_ui"] == "withheld"


def test_identity_remediation_is_result_transparent_and_scope_preserving() -> None:
    remediation = json.loads(REMEDIATION.read_text(encoding="utf-8"))
    assert remediation["failed_trial_identity"]["contract_passed"] is False
    assert remediation["failed_trial_identity"]["only_failed_gate"] == "identity_passed"
    assert remediation["known_results_not_prospective"]["immer_numeric_passed"] is True
    assert "utils::removeSource" in remediation["stable_function_identity_contract"]["algorithm"]
    assert "Do not modify a fixture" in remediation["unchanged_contract"]["forbidden"]


def test_identity_amendment_preserves_failed_replay_and_corrects_only_aggregation() -> None:
    remediation = json.loads(REMEDIATION.read_text(encoding="utf-8"))
    amendment = json.loads(REMEDIATION_AMENDMENT.read_text(encoding="utf-8"))
    failed = amendment["failed_corrected_algorithm_trial"]
    correction = amendment["bookkeeping_correction"]
    frozen = amendment["frozen_inputs"]
    assert failed["contract_passed"] is False
    assert failed["only_failed_gate"] == "identity_passed"
    assert correction["incorrect_per_function_sha256_retained_in_parent"] == (
        remediation["stable_function_identity_contract"]
        ["mfrmr_fit_mfrm_source_stripped_sha256"]
    )
    assert correction["correct_final_named_vector_sha256"] == (
        failed["observed_mfrmr_final_aggregate_sha256"]
    )
    assert "including when that vector has length one" in correction["algorithm"]
    assert not any(
        frozen[key]
        for key in (
            "fixture_changes",
            "tolerance_changes",
            "maxit_changes",
            "estimator_role_changes",
            "numerical_output_changes",
        )
    )
