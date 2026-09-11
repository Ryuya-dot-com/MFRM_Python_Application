from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pandas.testing as pdt

from mfrm_app.design_assignment import (
    ASSIGNMENT_EVIDENCE_VERSION,
    build_assignment_design_audit,
    build_estimand_contract,
    build_informative_assignment_evidence_register,
    infer_rater_facet,
)


def _balanced_rows() -> pd.DataFrame:
    rows = []
    for person in ("p1", "p2", "p3"):
        for rater in ("r1", "r2"):
            for task in ("t1", "t2"):
                rows.append(
                    {
                        "Person": person,
                        "Rater": rater,
                        "Task": task,
                        "Score": int(person[-1]) % 2,
                        "Estimate": 100.0,
                    }
                )
    return pd.DataFrame(rows)


def test_assignment_audit_counts_unique_person_rater_pairs_not_repeated_tasks():
    bundle = build_assignment_design_audit(
        _balanced_rows(),
        facet_names=["Rater", "Task"],
    )
    row = bundle["summary"].iloc[0]

    assert bundle["available"] is True
    assert row["AuditStatus"] == "Descriptive audit ready"
    assert row["UniquePersonRaterAssignments"] == 6
    assert row["LikelihoodRows"] == 12
    assert row["AssignmentCellDensity"] == 1.0
    assert row["MedianRatersPerPerson"] == 2.0
    assert row["CoRatedPersonShare"] == 1.0
    assert row["RaterOverlapComponents"] == 1
    assert row["FixedDensityStressTestStatus"] == "Eligible but not run for this dataset"
    assert bundle["sensitivity_plan"]["Status"].eq("Eligible; not run").all()
    assert bundle["rater_overlap"].iloc[0]["SharedPersons"] == 3


def test_assignment_audit_is_outcome_blind_and_fit_blind():
    original = _balanced_rows()
    perturbed = original.copy()
    perturbed["Score"] = [999 if idx % 2 else -999 for idx in range(len(perturbed))]
    perturbed["Estimate"] = range(len(perturbed))
    perturbed["StdResidual"] = 1e9

    first = build_assignment_design_audit(original, facet_names=["Rater", "Task"])
    second = build_assignment_design_audit(perturbed, facet_names=["Rater", "Task"])

    pdt.assert_frame_equal(first["summary"], second["summary"])
    pdt.assert_frame_equal(first["rater_exposure"], second["rater_exposure"])
    pdt.assert_frame_equal(first["rater_overlap"], second["rater_overlap"])
    assert first["summary"].iloc[0]["OutcomeBlindAudit"]


def test_assignment_audit_reports_disconnected_overlap_without_calling_it_mnar():
    data = pd.DataFrame(
        {
            "Person": ["p1", "p2", "p3"],
            "Rater": ["r1", "r2", "r3"],
            "Score": [0, 1, 2],
        }
    )
    bundle = build_assignment_design_audit(data, facet_names=["Rater"])
    row = bundle["summary"].iloc[0]

    assert row["AuditStatus"] == "Review direct rater overlap"
    assert row["RaterOverlapComponents"] == 3
    assert row["OverlappingRaterPairShare"] == 0.0
    assert row["AssignmentMechanismStatus"] == "Not identified from observed assignments"
    assert "do not establish" in row["ClaimBoundary"]
    assert row["FixedDensityStressTestStatus"].startswith("Not ready")
    assert bundle["sensitivity_plan"]["Status"].eq("Not ready").all()


def test_unconfirmed_first_facet_is_a_candidate_not_a_certified_rater_role():
    mapping = infer_rater_facet(["Prompt", "Criterion"])
    assert mapping == {
        "rater_facet": "Prompt",
        "mapping_basis": "first_facet_unconfirmed",
        "mapping_confirmed": False,
        "candidate_only": True,
    }

    data = pd.DataFrame(
        {
            "Person": ["p1", "p1", "p2", "p2"],
            "Prompt": ["a", "b", "a", "b"],
        }
    )
    row = build_assignment_design_audit(data, facet_names=["Prompt"])['summary'].iloc[0]
    assert row["AuditStatus"] == "Review rater mapping"
    assert not row["RaterMappingConfirmed"]
    assert row["FixedDensityStressTestStatus"].startswith("Not ready")


def test_configured_rater_role_takes_precedence_over_column_name_heuristic():
    mapping = infer_rater_facet(
        ["ScoringAgent", "RaterLikeButNotTheRole"],
        configured_rater_facet="ScoringAgent",
    )
    assert mapping["rater_facet"] == "ScoringAgent"
    assert mapping["mapping_basis"] == "configured_role"
    assert mapping["mapping_confirmed"] is True


def test_estimand_contract_marks_one_current_row_and_forbids_cross_basis_ranking():
    contract = build_estimand_contract(
        {"config": {"method": "MML", "estimate_population_sd": True}}
    )

    assert contract["Method"].tolist() == ["JMLE", "MML", "EXACT_CMLE"]
    assert int(contract["CurrentRun"].sum()) == 1
    current = contract.loc[contract["CurrentRun"]].iloc[0]
    assert current["Method"] == "MML"
    assert current["PopulationScale"] == "estimated Gaussian population SD"
    assert contract["ClaimBoundary"].str.contains("cross-basis", case=False).all()
    assert "two-decimal" in contract.loc[contract["Method"].eq("JMLE"), "FitPrecisionBoundary"].iloc[0]
    assert "FACETS 4.5.0" in contract.loc[contract["Method"].eq("JMLE"), "FACETSRelation"].iloc[0]
    assert "not a gold standard" in contract.loc[contract["Method"].eq("MML"), "FACETSRelation"].iloc[0]
    assert "Repository-only" in contract.loc[contract["Method"].eq("EXACT_CMLE"), "Availability"].iloc[0]


def test_validation_evidence_is_a_rationale_not_a_dataset_correction():
    evidence = build_informative_assignment_evidence_register()

    assert len(evidence) == 5
    assert evidence["EvidenceVersion"].eq(ASSIGNMENT_EVIDENCE_VERSION).all()
    assert evidence["EvidenceSource"].str.endswith("assessment_20260811.json").all()
    assert evidence["Replicates"].eq(100).all()
    assert evidence.iloc[0]["Estimate"] == 0.09863414655507183
    assert evidence["UseInCurrentDataset"].str.contains("do not transport", case=False).all()


def test_product_evidence_register_matches_retained_confirmatory_assessment():
    retained = json.loads(
        Path("validation/informative_assignment_confirmatory100_assessment_20260811.json")
        .read_text(encoding="utf-8")
    )
    expected = [retained["primary_result"]["mean_aligned_minus_planned"]]
    expected.extend(
        row.get("mean_aligned_minus_planned", row.get("mean_slope"))
        for row in retained["secondary_holm_family"]["results"]
    )
    evidence = build_informative_assignment_evidence_register()

    assert evidence["Estimate"].tolist() == expected
