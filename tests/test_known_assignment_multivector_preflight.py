from __future__ import annotations

import json

import pandas as pd

from validation import known_assignment_multivector_preflight as preflight


def test_frozen_multivector_plan_separates_preflight_from_confirmation():
    plan = json.loads(preflight.PLAN_PATH.read_text(encoding="utf-8"))

    assert plan["status"] == "frozen_before_new_person_or_response_generation"
    assert plan["evidence_status"]["person_vectors"] == 4
    assert plan["evidence_status"]["datasets"] == 12
    assert plan["evidence_status"]["attempt_units"] == 48
    assert plan["evidence_status"]["confirmatory_claims_allowed"] is False
    assert plan["evidence_status"]["screening_observations_pooled"] is False
    assert plan["evidence_status"]["preflight_observations_may_enter_confirmation"] is False
    assert plan["method_comparison_contract"]["estimator_ranking_prohibited"] is True
    assert "two-decimal" in plan["precision_contract"]["prohibition"]


def test_registered_preflight_diagnostics_follow_symmetric_contracts():
    loss_rows = []
    error_rows = []
    run_rows = []
    severity = {"R01": -0.45, "R02": -0.15, "R03": 0.15, "R04": 0.45}
    for vector in range(1, 5):
        for mode in (preflight.MML_FIXED_MODE, preflight.MML_FREE_MODE):
            for gamma, loss in ((-0.8, 0.20), (0.0, 0.10), (0.8, 0.20)):
                loss_rows.append(
                    {
                        "PersonVector": vector,
                        "EstimatorMode": mode,
                        "Gamma": gamma,
                        "RaterRMSE": loss,
                    }
                )
        for gamma, slope in ((-0.8, 0.4), (0.0, 0.0), (0.8, -0.4)):
            for level, truth in severity.items():
                error_rows.append(
                    {
                        "PersonVector": vector,
                        "EstimatorMode": preflight.MML_FREE_MODE,
                        "Gamma": gamma,
                        "Level": level,
                        "ErrorAligned": slope * truth,
                    }
                )
            run_rows.append(
                {
                    "PersonVector": vector,
                    "EstimatorMode": preflight.MML_FREE_MODE,
                    "Gamma": gamma,
                    "EstimatedPopulationSD": 0.8 if gamma == 0 else 0.7,
                }
            )

    diagnostics, summary, loss_wide, slope_wide = (
        preflight.summarize_preflight_diagnostics(
            rater_loss=pd.DataFrame(loss_rows),
            rater_errors=pd.DataFrame(error_rows),
            runs=pd.DataFrame(run_rows),
        )
    )

    assert len(diagnostics) == 16
    assert diagnostics.groupby("DiagnosticId")["PersonVector"].nunique().eq(4).all()
    assert summary["AdvancementSignal"].all()
    assert not summary["PValueComputed"].any()
    assert not summary["ConfirmatoryClaimAllowed"].any()
    assert set(loss_wide["SymmetricStressRaterRMSE"].round(12)) == {0.1}
    assert set(
        slope_wide.loc[
            slope_wide["EstimatorMode"].eq(preflight.MML_FREE_MODE),
            "DirectionAlignedSlopeHalfDifference",
        ].round(12)
    ) == {0.4}


def test_edge_signature_parser_is_lossless_and_score_free():
    signature = "P001::R01|P001::R02|P002::R02"
    assert preflight._parse_edges(signature) == {
        ("P001", "R01"),
        ("P001", "R02"),
        ("P002", "R02"),
    }

