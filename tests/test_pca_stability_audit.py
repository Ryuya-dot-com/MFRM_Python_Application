from __future__ import annotations

import numpy as np
import pandas as pd

import streamlit_app as app


def _stable_residual_matrix(n_persons: int = 80) -> pd.DataFrame:
    rng = np.random.default_rng(20260513)
    common = rng.normal(0.0, 1.0, size=n_persons)
    secondary = rng.normal(0.0, 0.15, size=n_persons)
    return pd.DataFrame(
        {
            "Task_R1": common + rng.normal(0.0, 0.08, size=n_persons),
            "Task_R2": common + secondary + rng.normal(0.0, 0.08, size=n_persons),
            "Task_R3": common - secondary + rng.normal(0.0, 0.08, size=n_persons),
            "Task_R4": common + rng.normal(0.0, 0.08, size=n_persons),
        },
        index=[f"P{i:03d}" for i in range(n_persons)],
    )


def test_pca_stability_audit_includes_loo_and_bootstrap_diagnostics():
    bundle = app.compute_pca_bundle(_stable_residual_matrix())

    assert bundle is not None
    stability = bundle["stability_table"]
    row = stability.iloc[0]

    expected_cols = {
        "FullEV1",
        "FullEV2",
        "FullPC1VariancePct",
        "LOOColumnsEvaluated",
        "LOOEV1MaxAbsDelta",
        "LOOEV1ShareMaxAbsDelta",
        "LOOPC1MedianAbsLoadingCorrelation",
        "BootstrapReplicates",
        "BootstrapEV1CV",
        "BootstrapEV1P05",
        "BootstrapEV1P95",
        "BootstrapPC1MedianAbsLoadingCorrelation",
        "BootstrapEV1ThresholdCrossing",
        "SensitivityFlagSummary",
    }
    assert expected_cols.issubset(stability.columns)
    assert row["PCAStabilityStatus"] == "Stable screen"
    assert int(row["LOOColumnsEvaluated"]) == 4
    assert int(row["BootstrapReplicates"]) == app.PCA_STABILITY_BOOTSTRAP_REPS
    assert float(row["LOOPC1MedianAbsLoadingCorrelation"]) >= 0.95
    assert float(row["BootstrapPC1MedianAbsLoadingCorrelation"]) >= 0.95
    assert row["SensitivityFlagSummary"] == "none"


def test_pca_stability_audit_flags_bootstrap_threshold_crossing():
    rng = np.random.default_rng(20260514)
    weak = pd.DataFrame(
        rng.normal(size=(36, 4)),
        columns=["A", "B", "C", "D"],
        index=[f"P{i}" for i in range(36)],
    )

    bundle = app.compute_pca_bundle(weak)

    assert bundle is not None
    row = bundle["stability_table"].iloc[0]
    assert row["PCAStabilityStatus"] == "Review"
    assert "bootstrap" in str(row["Caution"]) or "loading correlation" in str(row["Caution"])


def test_collect_pca_stability_tables_includes_overall_and_facet_scopes():
    overall = app.compute_pca_bundle(_stable_residual_matrix())
    facet = app.compute_pca_bundle(_stable_residual_matrix().iloc[:, :3])
    diagnostics = {
        "pca": overall,
        "pca_by_facet": {"Rater": facet},
    }

    combined = app.collect_pca_stability_tables(diagnostics)

    assert not combined.empty
    assert set(combined["Scope"]) == {"overall", "facet"}
    assert "Rater" in set(combined["Facet"])
    assert "BootstrapEV1CV" in combined.columns
