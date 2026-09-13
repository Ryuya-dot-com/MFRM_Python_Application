"""Contracts for the post-result known-truth pilot decomposition."""

from __future__ import annotations

from itertools import product

import numpy as np
import pandas as pd
import pytest

from mfrm_app.fit_threshold_operating_characteristics import (
    PERSON_MEASURE_TRUE,
    PERSON_MEASURE_WLE,
)
from mfrm_app.fit_threshold_pilot_addendum import (
    canonical_measure_source_pairs,
    replicate_wle_recovery,
    summarize_measure_source_pairs,
    summarize_wle_recovery,
    threshold_dimension_audit,
    wle_person_recovery,
)


def _person_rows() -> pd.DataFrame:
    rows = []
    for person, truth, wle, wle_infit, true_infit in (
        ("P1", 0.0, 0.2, 1.0, 1.6),
        ("P2", 1.0, 0.6, 1.7, 1.7),
    ):
        common = {
            "ConditionId": "C",
            "Replicate": 1,
            "Person": person,
            "TruthGroup": "clean",
            "TrueTheta": truth,
            "Outfit": 1.0,
            "WLEExtremeScorePattern": False,
        }
        rows.append(
            {
                **common,
                "PersonMeasureSource": PERSON_MEASURE_WLE,
                "PersonMeasure": wle,
                "Infit": wle_infit,
            }
        )
        rows.append(
            {
                **common,
                "PersonMeasureSource": PERSON_MEASURE_TRUE,
                "PersonMeasure": truth,
                "Infit": true_infit,
            }
        )
    return pd.DataFrame(rows)


def test_wle_recovery_and_replicate_summary_use_unrecentered_error():
    recovery = wle_person_recovery(_person_rows())
    assert recovery["WLEError"].tolist() == pytest.approx([0.2, -0.4])
    replicate = replicate_wle_recovery(recovery)
    assert float(replicate.iloc[0]["MeanBias"]) == pytest.approx(-0.1)
    assert float(replicate.iloc[0]["MeanAbsoluteError"]) == pytest.approx(0.3)
    assert float(replicate.iloc[0]["RootMeanSquaredError"]) == pytest.approx(
        np.sqrt(0.1)
    )
    summary = summarize_wle_recovery(replicate)
    assert set(summary["Metric"]) == {
        "MeanBias",
        "MeanAbsoluteError",
        "RootMeanSquaredError",
        "MedianError",
    }
    assert not summary["ConfirmatoryRecoveryClaim"].any()


def test_measure_source_pairing_retains_all_four_direction_labels():
    pairs = canonical_measure_source_pairs(_person_rows())
    comparison = pairs.set_index("Person")["FlagComparison"].to_dict()
    assert comparison == {"P1": "generating_theta_only", "P2": "both_flagged"}
    summary = summarize_measure_source_pairs(pairs).iloc[0]
    assert int(summary["GeneratingThetaOnly"]) == 1
    assert int(summary["BothFlagged"]) == 1
    assert int(summary["FlagDisagreements"]) == 1
    assert int(summary["IdentityUnavailable"]) == 0


def _threshold_rates() -> pd.DataFrame:
    rows = []
    for lower, acceptable, noisy in product(
        [0.45, 0.5], [1.4, 1.5], [1.9, 2.0]
    ):
        for rule in (
            "either_upper",
            "either_distorting",
            "either_nonacceptable",
            "either_overfit",
        ):
            relevant_count = {
                "either_upper": int(round(acceptable * 10)),
                "either_distorting": int(round(noisy * 10)),
                "either_nonacceptable": int(round(lower * 100 + acceptable * 10)),
                "either_overfit": int(round(lower * 100)),
            }[rule]
            rows.append(
                {
                    "ConditionId": "C",
                    "Replicate": 1,
                    "PersonMeasureSource": PERSON_MEASURE_WLE,
                    "TruthGroup": "clean",
                    "OverfitUpper": lower,
                    "AcceptableUpper": acceptable,
                    "NoisyUpper": noisy,
                    "Rule": rule,
                    "PersonsEligible": 100,
                    "Flagged": relevant_count,
                }
            )
    return pd.DataFrame(rows)


def test_threshold_dimension_audit_counts_effective_inputs_and_detects_mutation():
    rates = _threshold_rates()
    audit = threshold_dimension_audit(rates).set_index("Rule")
    assert audit["Passed"].all()
    assert int(audit.loc["either_upper", "NominalThresholdTriplets"]) == 8
    assert int(audit.loc["either_upper", "EffectiveInputConfigurations"]) == 2
    assert int(audit.loc["either_nonacceptable", "EffectiveInputConfigurations"]) == 4

    changed = rates.copy()
    row = changed.index[
        changed["Rule"].eq("either_upper")
        & changed["OverfitUpper"].eq(0.45)
        & changed["AcceptableUpper"].eq(1.4)
        & changed["NoisyUpper"].eq(1.9)
    ][0]
    changed.loc[row, "Flagged"] += 1
    failed = threshold_dimension_audit(changed).set_index("Rule")
    assert not bool(failed.loc["either_upper", "Passed"])
    assert int(failed.loc["either_upper", "IrrelevantDimensionMismatchGroups"]) == 1
