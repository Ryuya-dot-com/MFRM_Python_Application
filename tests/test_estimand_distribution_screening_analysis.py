import pandas as pd
import pytest

from validation.estimand_distribution_screening_analysis import (
    build_design_pairs,
    build_shape_pairs,
    summarize_contrasts,
)


def _run_metrics():
    rows = []
    for replicate in (1, 2):
        for design, design_shift in (("complete", 0.0), ("planned_connected", 0.2)):
            for distribution, shape_shift in (("normal", 0.0), ("right_skew", 0.1)):
                rows.append({
                    "RunId": f"{replicate}-{design}-{distribution}",
                    "Replicate": replicate,
                    "EstimatorMode": "mode",
                    "PersonDistribution": distribution,
                    "Design": design,
                    "Facet": "Rater",
                    "RMSE": 0.5 + design_shift + shape_shift,
                })
    return pd.DataFrame(rows)


def test_shape_pairs_are_replicate_matched_non_normal_minus_normal():
    pairs, summary = build_shape_pairs(
        _run_metrics(),
        identity_columns=["Replicate", "EstimatorMode", "Design", "Facet"],
        metric_columns=["RMSE"],
    )
    assert len(pairs) == 4
    assert pairs["RMSEContrast"].tolist() == pytest.approx([0.1] * 4)
    assert summary["PairedReplicates"].eq(2).all()
    assert summary["MeanContrast"].tolist() == pytest.approx([0.1, 0.1])


def test_design_pairs_are_replicate_matched_planned_minus_complete():
    pairs, summary = build_design_pairs(
        _run_metrics(),
        identity_columns=["Replicate", "EstimatorMode", "PersonDistribution", "Facet"],
        metric_columns=["RMSE"],
    )
    assert len(pairs) == 4
    assert pairs["RMSEContrast"].tolist() == pytest.approx([0.2] * 4)
    assert summary["PairedReplicates"].eq(2).all()


def test_contrast_summary_reports_mc_se_and_screening_interval():
    pairs = pd.DataFrame({
        "Group": ["A"] * 4,
        "RMSEContrast": [0.1, 0.2, 0.3, 0.4],
    })
    summary = summarize_contrasts(
        pairs, group_columns=["Group"], metric_columns=["RMSE"]
    ).iloc[0]
    assert summary["PairedReplicates"] == 4
    assert summary["MeanContrast"] == pytest.approx(0.25)
    assert summary["MonteCarloSE"] == pytest.approx(pd.Series([0.1, 0.2, 0.3, 0.4]).std() / 2)
    assert summary["EvidenceTier"] == "20-replicate screening; not confirmatory"
