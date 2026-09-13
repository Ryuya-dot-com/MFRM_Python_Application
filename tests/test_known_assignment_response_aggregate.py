from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from validation import known_assignment_response_aggregate as aggregate
from validation import known_assignment_response_run_v2 as wrapper


def test_screening_interval_helper_has_registered_nonconfirmatory_contract():
    frame = pd.DataFrame(
        {
            "EstimatorMode": ["mode"] * 10,
            "Gamma": [0.8] * 10,
            "Value": np.arange(10, dtype=float),
        }
    )
    result = aggregate._summary_interval(
        frame,
        group_columns=["EstimatorMode", "Gamma"],
        value_column="Value",
        value_name="Value",
    ).iloc[0]

    assert result["N"] == 10
    assert np.isclose(result["MeanValue"], 4.5)
    assert result["ScreeningCI95Low"] < result["MeanValue"] < result["ScreeningCI95High"]
    assert not result["ConfirmatoryClaimAllowed"]


@pytest.mark.retained_evidence
def test_full_registered_denominator_is_complete_before_analysis_registration():
    markers = sorted((wrapper.STUDY_DIR / "attempts").rglob("completion.json"))
    values = [json.loads(path.read_text(encoding="utf-8")) for path in markers]

    assert len(values) == 120
    assert all(value["execution_completed"] for value in values)
    assert all(value["statistical_evidence_ready"] for value in values)
    resilient = [
        value
        for value in values
        if value["attempt_type"] == "RESILIENT_FACETS_PYTHON_JMLE_PCM"
    ]
    assert len(resilient) == 30
    assert all(value["facets_calibration_ready"] for value in resilient)


def test_registered_context_replaces_blank_optional_adapter_metadata():
    frame = pd.DataFrame(
        {
            "RunId": ["run-1", "run-1"],
            "EstimatorMode": ["mml", "mml"],
            "Replicate": [np.nan, np.nan],
            "Gamma": [np.nan, np.nan],
        }
    )
    manifest = pd.DataFrame(
        {
            "RunId": ["run-1"],
            "ConditionId": ["condition"],
            "Design": ["gamma_pos_0p8"],
            "Gamma": [0.8],
            "PersonDistribution": ["fixed_normal_vector"],
            "Replicate": [7],
            "AssignmentCorrelation": [0.4],
            "AssignmentStatistic": [64.0],
        }
    )

    attached = aggregate._attach_registered_context(frame, manifest)

    assert attached["Replicate"].eq(7).all()
    assert attached["Gamma"].eq(0.8).all()
    assert attached["EstimatorMode"].eq("mml").all()
