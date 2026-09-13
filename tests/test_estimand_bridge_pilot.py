from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from validation import estimand_bridge_pilot as bridge


PARENT = Path("validation/facets_pcm_boundary_pilot_20260811")


@pytest.mark.retained_evidence
def test_parent_identity_and_rank_full_selection_are_locked():
    hashes = bridge.validate_parent_identity(PARENT)
    _, manifest = bridge.selected_manifest(PARENT)

    assert "generated_ratings.csv" in hashes
    assert len(manifest) == 8
    assert set(manifest["Design"]) == {"complete", "planned_connected"}
    assert manifest["ExpectedStructuralNullity"].eq(0).all()


@pytest.mark.retained_evidence
def test_parent_jmle_import_preserves_estimand_labels():
    _, manifest = bridge.selected_manifest(PARENT)
    ledger, recovery, thresholds = bridge.import_parent_jmle(
        PARENT, set(manifest["RunId"])
    )

    assert len(ledger) == 16
    assert len(recovery) == 288
    assert len(thresholds) == 144
    assert set(recovery["EstimatorMode"]) == {
        bridge.PARENT_FACETS_MODE,
        bridge.PARENT_PYTHON_MODE,
    }
    assert recovery["EstimandClass"].eq("fixed_person_joint_likelihood").all()
    assert thresholds["ThresholdCondition"].notna().all()


def test_mml_constraint_audit_excludes_location_bearing_criterion():
    facets = pd.DataFrame({
        "Facet": ["Rater", "Rater", "Task", "Task", "Criterion", "Criterion"],
        "Level": ["R1", "R2", "T1", "T2", "C1", "C2"],
        "Estimate": [-0.2, 0.2, -0.1, 0.1, 0.5, 0.7],
    })
    steps = pd.DataFrame({
        "StepFacetLevel": ["C1", "C1", "C2", "C2"],
        "Estimate": [-0.4, 0.4, -0.3, 0.3],
    })

    audit = bridge._constraint_audit(
        facets,
        steps,
        estimator_mode="PYTHON_MML_FREE_SD_Q31",
        fit_model="PCM",
    )

    assert audit["ConstraintPass"]
    assert "Criterion" not in audit["ConstrainedFacetSums"]


def test_cmle_constraint_audit_includes_all_structural_facets():
    facets = pd.DataFrame({
        "Facet": ["Rater", "Rater", "Task", "Task", "Criterion", "Criterion"],
        "Level": ["R1", "R2", "T1", "T2", "C1", "C2"],
        "Estimate": [-0.2, 0.2, -0.1, 0.1, -0.3, 0.3],
    })
    steps = pd.DataFrame({
        "StepFacetLevel": ["__COMMON__", "__COMMON__"],
        "Estimate": [-0.4, 0.4],
    })

    audit = bridge._constraint_audit(
        facets, steps, estimator_mode=bridge.CMLE_MODE, fit_model="RSM"
    )

    assert audit["ConstraintPass"]
    assert set(audit["ConstrainedFacetSums"]) == {"Rater", "Task", "Criterion"}


def test_model_direction_checks_keep_likelihood_bases_separate():
    parent = pd.DataFrame(columns=[
        "RunId", "Design", "ThresholdCondition", "Replicate", "Rows", "FitModel",
        "ComparisonEligible", "PythonLogLik", "PythonLogLikPerObs", "PythonAIC", "PythonBIC",
    ])
    rows = []
    for condition, gain in (("shared", 0.01), ("heterogeneous", 0.20)):
        for model, loglik in (("RSM", -10.0), ("PCM", -10.0 + gain)):
            rows.append({
                "RunId": f"run-{condition}",
                "Design": "complete",
                "ThresholdCondition": condition,
                "Replicate": 1,
                "Rows": 100,
                "FitModel": model,
                "EstimatorMode": bridge.CMLE_MODE,
                "LikelihoodBasis": "exact_person_total_conditional",
                "IncludedInBridge": True,
                "LogLik": loglik,
                "LogLikPerObs": loglik / 100,
                "AIC": -2 * loglik,
                "BIC": float("nan"),
            })
    pairs, checks = bridge.build_model_direction_checks(parent, pd.DataFrame(rows))

    assert len(pairs) == 2
    assert len(checks) == 1
    assert bool(checks.iloc[0]["DirectionPass"])
    assert checks.iloc[0]["LikelihoodBasis"] == "exact_person_total_conditional"
