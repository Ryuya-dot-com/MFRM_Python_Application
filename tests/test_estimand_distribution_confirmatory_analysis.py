import json

import numpy as np
import pandas as pd
import pytest

from validation.estimand_distribution_confirmatory_analysis import (
    DEFAULT_PLAN_PATH,
    EXPECTED_ENDPOINT_IDS,
    compute_primary_contrasts,
    holm_adjust,
    summarize_primary,
)


def _synthetic_primary_tables(replicates=3):
    run_rows = []
    threshold_rows = []
    levels = [("C01", 1), ("C01", 2), ("C01", 3), ("C02", 1), ("C02", 2), ("C02", 3)]
    for replicate in range(21, 21 + replicates):
        for distribution, value in (
            ("normal", 0.8),
            ("right_skew", 0.7),
            ("heavy_tail_t3", 0.75),
        ):
            run_rows.append({
                "Replicate": replicate,
                "EstimatorMode": "PYTHON_MML_FREE_SD_Q31",
                "Design": "planned_connected",
                "PersonDistribution": distribution,
                "EstimatedPopulationSD": value + replicate * 1e-4,
                "IncludedInStudy": True,
            })
        for design, error in (("complete", 0.1), ("planned_connected", 0.2)):
            for criterion, category in levels:
                threshold_rows.append({
                    "Replicate": replicate,
                    "EstimatorMode": "PYTHON_JMLE",
                    "Design": design,
                    "PersonDistribution": "normal",
                    "StepFacetLevel": criterion,
                    "Category": category,
                    "TruthError": error,
                    "IncludedInStudy": True,
                })
        for mode in (
            "PYTHON_JMLE",
            "PYTHON_MML_FREE_SD_Q31",
            "PYTHON_EXACT_CMLE",
        ):
            if mode != "PYTHON_JMLE":
                threshold_rows.append({
                    "Replicate": replicate,
                    "EstimatorMode": mode,
                    "Design": "planned_connected",
                    "PersonDistribution": "normal",
                    "StepFacetLevel": "C02",
                    "Category": 2,
                    "TruthError": 0.2,
                    "IncludedInStudy": True,
                })
            threshold_rows.append({
                "Replicate": replicate,
                "EstimatorMode": mode,
                "Design": "planned_connected",
                "PersonDistribution": "symmetric_mixture",
                "StepFacetLevel": "C02",
                "Category": 2,
                "TruthError": 0.1,
                "IncludedInStudy": True,
            })
    return pd.DataFrame(run_rows), pd.DataFrame(threshold_rows)


def test_primary_contrast_extraction_is_exact_and_directional():
    runs, thresholds = _synthetic_primary_tables()
    contrasts = compute_primary_contrasts(runs, thresholds)
    assert tuple(dict.fromkeys(contrasts["EndpointId"])) == EXPECTED_ENDPOINT_IDS
    assert contrasts.groupby("EndpointId").size().eq(3).all()
    expected = {
        EXPECTED_ENDPOINT_IDS[0]: -0.1,
        EXPECTED_ENDPOINT_IDS[1]: -0.05,
        EXPECTED_ENDPOINT_IDS[2]: -0.1,
        EXPECTED_ENDPOINT_IDS[3]: -0.1,
        EXPECTED_ENDPOINT_IDS[4]: -0.1,
        EXPECTED_ENDPOINT_IDS[5]: 0.1,
    }
    observed = contrasts.groupby("EndpointId")["Contrast"].mean().to_dict()
    for endpoint_id, value in expected.items():
        assert observed[endpoint_id] == pytest.approx(value)


def test_threshold_rmse_requires_all_six_registered_elements():
    runs, thresholds = _synthetic_primary_tables()
    drop_index = thresholds[
        thresholds["EstimatorMode"].eq("PYTHON_JMLE")
        & thresholds["Design"].eq("complete")
        & thresholds["PersonDistribution"].eq("normal")
    ].index[0]
    with pytest.raises(ValueError, match="exactly six"):
        compute_primary_contrasts(runs, thresholds.drop(index=drop_index))


def test_holm_adjustment_matches_step_down_definition():
    adjusted = holm_adjust([0.01, 0.04, 0.03])
    assert adjusted.tolist() == pytest.approx([0.03, 0.06, 0.06])


def test_summary_enforces_100_pair_gate_and_separates_precision():
    plan = json.loads(DEFAULT_PLAN_PATH.read_text(encoding="utf-8"))
    rows = []
    rng = np.random.default_rng(8142)
    for index, endpoint_id in enumerate(EXPECTED_ENDPOINT_IDS):
        direction = 1.0 if endpoint_id == EXPECTED_ENDPOINT_IDS[-1] else -1.0
        for replicate, value in enumerate(
            direction * (0.3 + rng.normal(scale=0.01, size=100)), start=21
        ):
            rows.append({
                "EndpointId": endpoint_id,
                "Replicate": replicate,
                "Contrast": value,
            })
    contrasts = pd.DataFrame(rows)
    result = summarize_primary(contrasts, plan)
    assert result["FullPairGate"].all()
    assert result["DirectionConfirmed"].all()
    assert result["PrecisionPass"].all()
    incomplete = summarize_primary(
        contrasts[~(
            contrasts["EndpointId"].eq(EXPECTED_ENDPOINT_IDS[0])
            & contrasts["Replicate"].eq(21)
        )],
        plan,
    )
    first = incomplete.set_index("EndpointId").loc[EXPECTED_ENDPOINT_IDS[0]]
    assert not bool(first["FullPairGate"])
    assert not bool(first["DirectionConfirmed"])
