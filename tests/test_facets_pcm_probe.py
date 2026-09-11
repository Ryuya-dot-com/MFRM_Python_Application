from __future__ import annotations

import pandas as pd
import pytest

from validation import facets_pcm_probe as probe


def _manifest_row() -> pd.Series:
    return pd.Series({
        "RunId": "run-1",
        "ConditionId": "condition-1",
        "Design": "balanced",
        "TruthBias": 0.0,
        "Replicate": 1,
        "Seed": 1,
        "Categories": 4,
    })


def _ratings() -> pd.DataFrame:
    return pd.DataFrame({
        "Person": ["P1", "P1", "P2", "P2"],
        "Rater": ["R1", "R1", "R1", "R1"],
        "Task": ["T1", "T1", "T1", "T1"],
        "Criterion": ["C1", "C2", "C1", "C2"],
        "Score": [0, 1, 2, 3],
    })


def test_pcm_spec_uses_one_hash_facet_and_removes_bias_model(tmp_path):
    spec, maps = probe.build_pcm_spec(
        _manifest_row(),
        _ratings(),
        pd.DataFrame(columns=["Facet", "Level", "Anchor"]),
        score_base=tmp_path / "scores.txt",
    )

    assert "Models=\n?,?,?,#,R3\n*" in spec
    assert "?B" not in spec
    assert maps["Criterion"] == {"C1": 1, "C2": 2}


def test_map_pcm_scale_tables_is_one_to_one():
    categories = pd.DataFrame({
        "TableNumber": ["8.1", "8.1", "8.2", "8.2"],
        "Model": ["?,?,?,1,R3", "?,?,?,1,R3", "?,?,?,2,R3", "?,?,?,2,R3"],
        "Category": [0, 1, 0, 1],
    })

    mapped = probe.map_pcm_scale_tables(
        categories,
        criterion_map={"C1": 1, "C2": 2},
    )

    assert set(mapped.loc[mapped["TableNumber"].eq("8.1"), "StepFacetLevel"]) == {"C1"}
    assert set(mapped.loc[mapped["TableNumber"].eq("8.2"), "StepFacetLevel"]) == {"C2"}


def test_map_pcm_scale_tables_fails_on_missing_scale():
    categories = pd.DataFrame({
        "TableNumber": ["8.1", "8.1"],
        "Model": ["?,?,?,1,R3", "?,?,?,1,R3"],
        "Category": [0, 1],
    })

    with pytest.raises(ValueError, match="not one-to-one"):
        probe.map_pcm_scale_tables(
            categories,
            criterion_map={"C1": 1, "C2": 2},
        )


def test_normalize_python_pcm_steps_extracts_transition_category():
    result = {
        "steps": pd.DataFrame({
            "StepFacet": ["C1", "C1", "C2", "C2"],
            "Step": ["Step_1", "Step_2", "Step_1", "Step_2"],
            "Estimate": [-0.5, 0.5, -0.2, 0.2],
        })
    }

    output = probe.normalize_python_pcm_steps(result)

    assert list(output["Category"]) == [1, 2, 1, 2]
    assert output.groupby("StepFacetLevel")["PythonThreshold"].sum().eq(0).all()
