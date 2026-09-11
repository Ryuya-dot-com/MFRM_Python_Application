from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from validation.operating_characteristics_python_mml import (
    DEFAULT_MODES,
    build_mml_kwargs,
    normalize_mml_recovery,
    validate_registered_bundle_hashes,
)
from validation.operating_characteristics_facets import sha256_file


def _manifest_row() -> pd.Series:
    return pd.Series({
        "RunId": "run-1",
        "ConditionId": "condition-1",
        "Design": "balanced",
        "TruthBias": 0.0,
        "Replicate": 1,
        "Seed": 123,
        "Categories": 4,
    })


def test_registered_default_modes_isolate_population_sd_choice() -> None:
    assert DEFAULT_MODES == ("PYTHON_MML_FIXED_SD1_Q31", "PYTHON_MML_FREE_SD_Q31")
    fixed = build_mml_kwargs(DEFAULT_MODES[0], 4)
    free = build_mml_kwargs(DEFAULT_MODES[1], 4)
    assert fixed["quad_points"] == free["quad_points"] == 31
    assert fixed["population_prior_sd"] == free["population_prior_sd"] == 1.0
    assert fixed["estimate_population_sd"] is False
    assert free["estimate_population_sd"] is True
    assert fixed["noncenter_facet"] == free["noncenter_facet"] == "Criterion"


def test_unknown_mode_fails_closed() -> None:
    with pytest.raises(ValueError, match="Unknown native MML mode"):
        build_mml_kwargs("invented", 4)


def test_unanchored_recovery_is_mean_aligned() -> None:
    measures = pd.DataFrame({
        "Facet": ["Rater", "Rater", "Task", "Criterion"],
        "Level": ["R1", "R2", "T1", "C1"],
        "Estimate": [0.4, -0.2, 0.3, -0.1],
        "SE": [0.1, 0.1, 0.1, 0.1],
    })
    truth = pd.DataFrame({
        "Facet": ["Rater", "Rater", "Task", "Criterion"],
        "Level": ["R1", "R2", "T1", "C1"],
        "Truth": [0.3, -0.3, 0.0, 0.0],
    })
    recovery = normalize_mml_recovery(
        measures, truth, pd.DataFrame(columns=["Facet", "Level", "Anchor"]),
        _manifest_row(), mode=DEFAULT_MODES[0], eligible=True,
        structurally_identified=True,
    )
    rater = recovery[recovery["Facet"].eq("Rater")]
    assert np.isclose(rater["ErrorAligned"].mean(), 0.0)
    assert set(recovery["Engine"]) == {"PythonApp"}
    assert set(recovery["Estimator"]) == {"MML"}
    assert set(recovery["ComparisonClass"]) == {"native_marginal_sensitivity"}


def test_anchored_facet_stays_on_absolute_scale() -> None:
    measures = pd.DataFrame({
        "Facet": ["Rater", "Rater", "Task", "Criterion"],
        "Level": ["R1", "R2", "T1", "C1"],
        "Estimate": [0.3, 0.2, 0.0, 0.0],
        "SE": [np.nan, 0.1, 0.1, 0.1],
    })
    truth = pd.DataFrame({
        "Facet": ["Rater", "Rater", "Task", "Criterion"],
        "Level": ["R1", "R2", "T1", "C1"],
        "Truth": [0.3, -0.3, 0.0, 0.0],
    })
    anchors = pd.DataFrame({"Facet": ["Rater"], "Level": ["R1"], "Anchor": [0.3]})
    recovery = normalize_mml_recovery(
        measures, truth, anchors, _manifest_row(), mode=DEFAULT_MODES[1],
        eligible=True, structurally_identified=True,
    )
    rater = recovery[recovery["Facet"].eq("Rater")]
    assert set(rater["ComparisonScale"]) == {"anchor_identified_absolute"}
    assert np.allclose(rater["EstimateAligned"], rater["Estimate"])


def test_structural_nonidentification_is_annotation_not_mml_exclusion() -> None:
    measures = pd.DataFrame({
        "Facet": ["Rater", "Task", "Criterion"],
        "Level": ["R1", "T1", "C1"],
        "Estimate": [0.0, 0.0, 0.0],
        "SE": [0.1, 0.1, 0.1],
    })
    truth = measures.rename(columns={"Estimate": "Truth"}).drop(columns="SE")
    recovery = normalize_mml_recovery(
        measures, truth, pd.DataFrame(columns=["Facet", "Level", "Anchor"]),
        _manifest_row(), mode=DEFAULT_MODES[0], eligible=True,
        structurally_identified=False,
    )
    assert recovery["IncludedInSummary"].all()
    assert not recovery["ObservedDesignStructurallyIdentified"].any()
    assert set(recovery["MMLIdentificationBasis"]) == {"normal_population_distribution"}


def test_registered_bundle_hash_mismatch_fails_closed(tmp_path: Path) -> None:
    payload = tmp_path / "generated_ratings.csv"
    payload.write_text("a\n1\n", encoding="utf-8")
    pd.DataFrame({
        "File": [payload.name],
        "SHA256": [sha256_file(payload)],
    }).to_csv(tmp_path / "generated_bundle_files.csv", index=False)
    assert validate_registered_bundle_hashes(tmp_path)[payload.name] == sha256_file(payload)
    payload.write_text("a\n2\n", encoding="utf-8")
    with pytest.raises(ValueError, match="SHA256 mismatch"):
        validate_registered_bundle_hashes(tmp_path)
