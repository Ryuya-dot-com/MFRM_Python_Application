from __future__ import annotations

import pandas as pd

from validation.operating_characteristics_facets import (
    facets_python_table8_threshold_agreement,
)
from validation.operating_characteristics_python_table8 import (
    build_jmle_kwargs,
    normalize_steps,
)


def _manifest() -> pd.Series:
    return pd.Series({
        "RunId": "run-1",
        "ConditionId": "condition-1",
        "Design": "balanced",
        "TruthBias": 0.0,
        "TruthPositive": False,
        "BiasCell": "R01 x T01",
        "Replicate": 1,
        "Seed": 123,
        "Categories": 4,
    })


def test_jmle_step_controls_reproduce_registered_pilot():
    kwargs = build_jmle_kwargs(4)
    assert kwargs["method"] == "JMLE"
    assert kwargs["model"] == "RSM"
    assert kwargs["noncenter_facet"] == "Person"
    assert kwargs["maxit"] == 160
    assert kwargs["reltol"] == 1e-6
    assert kwargs["rating_max"] == 3


def test_normalize_steps_maps_step_number_to_facets_category():
    result = {
        "steps": pd.DataFrame({
            "Step": ["Step_1", "Step_2", "Step_3"],
            "Estimate": [-1.2, 0.0, 1.2],
        })
    }
    output = normalize_steps(result, _manifest())
    assert list(output["Category"]) == [1, 2, 3]
    assert list(output["ThresholdEstimate"]) == [-1.2, 0.0, 1.2]
    assert set(output["Mode"]) == {"PYTHON_JMLE_REGISTERED_RSM_STEPS"}


def test_facets_python_threshold_agreement_uses_primary_u6_interval(tmp_path):
    categories = pd.DataFrame({
        "RunId": ["run-1"],
        "ConditionId": ["condition-1"],
        "Design": ["balanced"],
        "TruthBias": [0.0],
        "Replicate": [1],
        "Seed": [123],
        "Category": [1],
        "PrimaryU6ThresholdMeasureToken": ["-1.2345"],
        "PrimaryU6ThresholdMeasureDisplayed": [-1.2345],
        "PrimaryU6ThresholdMeasureDisplayDecimals": [4],
        "PrimaryU6ThresholdMeasureRawLowerBound": [-1.23455],
        "PrimaryU6ThresholdMeasureRawUpperBound": [-1.23445],
        "DirectComparisonEligible": [True],
        "ComparisonClass": ["matched_jmle"],
    })
    pd.DataFrame({
        "RunId": ["run-1"],
        "Category": [1],
        "ThresholdEstimate": [-1.23449],
        "IncludedInComparison": [True],
    }).to_csv(tmp_path / "python_jmle_steps.csv", index=False)

    pairs, summary = facets_python_table8_threshold_agreement(categories, tmp_path)

    assert len(pairs) == 1
    assert pairs.loc[0, "PythonWithinFACETSDisplayInterval"]
    assert pairs.loc[0, "IncludedInDirectComparison"]
    assert summary.loc[0, "WithinFACETSDisplayIntervalRate"] == 1.0
