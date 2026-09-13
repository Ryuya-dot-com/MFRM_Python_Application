from __future__ import annotations

import numpy as np
import pandas as pd

from mfrm_app import decision_stability as ds


def test_fit_endpoint_semantics_are_explicit():
    assert ds.classify_fit_mnsq(np.nextafter(0.5, -np.inf)) == "overfit"
    assert ds.classify_fit_mnsq(0.5) == "acceptable"
    assert ds.classify_fit_mnsq(1.5) == "acceptable"
    assert ds.classify_fit_mnsq(np.nextafter(1.5, np.inf)) == "noisy"
    assert ds.classify_fit_mnsq(2.0) == "noisy"
    assert ds.classify_fit_mnsq(np.nextafter(2.0, np.inf)) == "distorting"


def test_exact_thresholds_are_classified_but_marked_as_numerical_boundaries():
    for value in (0.5, 1.5, 2.0):
        evidence = ds.evaluate_fit_mnsq(value)
        assert evidence["BoundaryStatus"] == "numerical_boundary"
        assert evidence["DecisionStable"] is False


def test_raw_decision_never_uses_rounded_display_value():
    evidence = ds.evaluate_fit_mnsq(1.5004, display_decimals=3)

    assert evidence["RawDecision"] == "noisy"
    assert evidence["DisplayValue"] == 1.5
    assert evidence["DisplayDecision"] == "acceptable"
    assert evidence["DisplayDecisionConsistent"] is False
    assert evidence["BoundaryStatus"] == "display_rounding_boundary"


def test_far_from_threshold_is_stable():
    evidence = ds.evaluate_fit_mnsq(1.2344, display_decimals=3)

    assert evidence["RawDecision"] == "acceptable"
    assert evidence["DisplayDecisionConsistent"] is True
    assert evidence["BoundaryStatus"] == "stable"
    assert evidence["DecisionStable"] is True


def test_nonfinite_value_is_typed_as_unavailable():
    evidence = ds.evaluate_fit_mnsq(np.nan)

    assert evidence["BoundaryStatus"] == "unavailable"
    assert evidence["RawDecision"] == "unavailable"
    assert evidence["DecisionStable"] is False


def test_fit_audit_preserves_identity_and_uses_absolute_zstd():
    frame = pd.DataFrame({
        "Facet": ["Rater", "Rater"],
        "Level": ["R1", "R2"],
        "Infit": [1.5004, 1.2],
        "Outfit": [2.0004, 0.4996],
        "InfitZSTD": [-2.0, 0.3],
    })
    audit = ds.audit_fit_decision_stability(frame)

    assert set(audit["Statistic"]) == {"Infit", "Outfit", "InfitZSTD"}
    assert set(audit["Facet"]) == {"Rater"}
    assert set(audit["Level"]) == {"R1", "R2"}
    zstd = audit.loc[(audit["Level"] == "R1") & (audit["Statistic"] == "InfitZSTD")].iloc[0]
    assert zstd["DecisionValue"] == 2.0
    assert zstd["RawDecision"] == "review"
    assert zstd["BoundaryStatus"] == "numerical_boundary"


def test_fit_summary_counts_rounding_mismatches():
    audit = ds.audit_fit_decision_stability(pd.DataFrame({"Infit": [1.5004, 1.2]}))
    summary = ds.summarize_fit_decision_stability(audit).iloc[0]

    assert summary["Status"] == "Review"
    assert summary["DisplayBoundaryStatistics"] == 1
    assert summary["DisplayDecisionMismatches"] == 1


def test_fit_summary_does_not_treat_unavailable_statistics_as_ready():
    all_missing = ds.audit_fit_decision_stability(
        pd.DataFrame({"Infit": [float("nan")], "Outfit": [float("nan")]})
    )
    missing_summary = ds.summarize_fit_decision_stability(all_missing).iloc[0]

    assert missing_summary["Status"] == "Missing"
    assert missing_summary["AvailableStatistics"] == 0
    assert missing_summary["UnavailableStatistics"] == 2

    partial = ds.audit_fit_decision_stability(
        pd.DataFrame({"Infit": [1.0], "Outfit": [float("nan")]})
    )
    partial_summary = ds.summarize_fit_decision_stability(partial).iloc[0]

    assert partial_summary["Status"] == "Review"
    assert partial_summary["AvailableStatistics"] == 1
    assert partial_summary["UnavailableStatistics"] == 1


def test_bias_threshold_audit_separates_raw_and_display_decisions():
    frame = pd.DataFrame({
        "FacetPair": ["Rater x Task"],
        "p_holm": [0.04996],
        "p_bh": [0.04996],
        "AbsBias": [0.49996],
    })
    audit = ds.audit_bias_decision_stability(frame, display_decimals=3)

    assert len(audit) == 3
    assert set(audit["BoundaryStatus"]) == {"display_rounding_boundary"}
    assert (~audit["DisplayDecisionConsistent"]).all()
