from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from validation import informative_assignment_confirmatory as confirm
from validation import informative_assignment_screening as screening


def test_registered_range_matches_frozen_confirmation() -> None:
    plan = json.loads(confirm.PLAN_PATH.read_text(encoding="utf-8"))
    assert confirm.registered_range(plan) == (241, 340, 100)
    assert plan["operational_contract"]["attempt_units"] == 800
    assert tuple(plan["observation_designs"]) == confirm.CONFIRMATORY_DESIGNS
    assert plan["primary_endpoint"]["id"] == confirm.PRIMARY_ID
    assert tuple(
        endpoint["id"] for endpoint in plan["secondary_family"]["endpoints"]
    ) == confirm.SECONDARY_IDS


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("replicate_range", [242, 341]),
        ("replicates", 99),
        ("pool_replicates_1_through_240", True),
        ("fixed_sample_size", False),
        ("optional_stopping", True),
        ("interim_scientific_looks", 1),
    ],
)
def test_registered_range_rejects_plan_drift(field: str, value: object) -> None:
    plan = json.loads(confirm.PLAN_PATH.read_text(encoding="utf-8"))
    plan["evidence_status"][field] = value
    with pytest.raises(ValueError):
        confirm.registered_range(plan)


def test_screening_engine_context_restores_shared_module() -> None:
    original_designs = screening.DESIGNS
    original_range = screening._registered_range
    original_validate = screening.validate_study_identity
    with confirm._screening_engine_context():
        assert screening.DESIGNS == confirm.CONFIRMATORY_DESIGNS
        assert screening._registered_range is confirm.registered_range
        assert screening.validate_study_identity is confirm.validate_study_identity
    assert screening.DESIGNS is original_designs
    assert screening._registered_range is original_range
    assert screening.validate_study_identity is original_validate


def test_holm_adjust_is_step_down_and_order_stable() -> None:
    result = confirm.holm_adjust([0.04, 0.001, 0.02, 0.02])
    assert np.allclose(result, [0.06, 0.004, 0.06, 0.06])


def _synthetic_endpoint_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    replicates = np.arange(241, 341)
    paired_rows = []
    for replicate in replicates:
        offset = (replicate - replicates.mean()) / 10000.0
        for mode, contrast in (
            (confirm.MML_FREE_MODE, 0.075 + offset),
            (confirm.MML_FIXED_MODE, 0.040 + offset),
        ):
            paired_rows.append({
                "Replicate": replicate,
                "EstimatorMode": mode,
                "RecoveryDomain": "Facet:Rater",
                "Metric": "RMSE",
                "ContrastAlignedMinusPlanned": contrast,
            })
    rater_rows = []
    for replicate in replicates:
        jitter = (replicate - replicates.mean()) / 100000.0
        for mode, slope in (
            (confirm.MML_FREE_MODE, -0.52 + jitter),
            (confirm.MML_FIXED_MODE, -0.40 + jitter),
        ):
            for level, truth in confirm.RATER_TRUTH.items():
                rater_rows.append({
                    "Replicate": replicate,
                    "EstimatorMode": mode,
                    "Level": level,
                    "ErrorContrastAlignedMinusPlanned": slope * truth,
                })
    free_sd = pd.DataFrame({
        "Replicate": replicates,
        "ContrastAlignedMinusPlanned": -0.125
        + (replicates - replicates.mean()) / 10000.0,
    })
    return pd.DataFrame(paired_rows), pd.DataFrame(rater_rows), free_sd


def test_endpoint_construction_has_five_fresh_hundred_pair_series() -> None:
    paired, rater, free_sd = _synthetic_endpoint_inputs()
    output = confirm.build_endpoint_contrasts(paired, rater, free_sd)
    assert tuple(output["EndpointId"].drop_duplicates()) == (
        confirm.PRIMARY_ID,
        *confirm.SECONDARY_IDS,
    )
    assert output.groupby("EndpointId").size().eq(100).all()
    means = output.groupby("EndpointId")["Contrast"].mean()
    assert np.isclose(means[confirm.PRIMARY_ID], 0.075)
    assert np.isclose(means[confirm.SECONDARY_IDS[0]], 0.040)
    assert np.isclose(means[confirm.SECONDARY_IDS[1]], -0.125)
    assert np.isclose(means[confirm.SECONDARY_IDS[2]], -0.52)
    assert np.isclose(means[confirm.SECONDARY_IDS[3]], -0.40)


def test_endpoint_summary_applies_primary_and_holm_gates() -> None:
    paired, rater, free_sd = _synthetic_endpoint_inputs()
    contrasts = confirm.build_endpoint_contrasts(paired, rater, free_sd)
    primary, secondary = confirm.summarize_endpoints(contrasts)
    assert primary.loc[0, "FinitePairedReplicates"] == 100
    assert bool(primary.loc[0, "DirectionConfirmed"])
    assert bool(primary.loc[0, "PrecisionPass"])
    assert secondary["FinitePairedReplicates"].eq(100).all()
    assert secondary["PrimaryGatePass"].all()
    assert secondary["DirectionConfirmed"].all()
    assert (secondary["HolmAdjustedP"] >= secondary["RawOneSidedP"]).all()


def test_endpoint_summary_requires_all_hundred_primary_pairs() -> None:
    paired, rater, free_sd = _synthetic_endpoint_inputs()
    contrasts = confirm.build_endpoint_contrasts(paired, rater, free_sd)
    contrasts = contrasts[
        ~(
            contrasts["EndpointId"].eq(confirm.PRIMARY_ID)
            & contrasts["Replicate"].eq(340)
        )
    ]
    primary, secondary = confirm.summarize_endpoints(contrasts)
    assert primary.loc[0, "FinitePairedReplicates"] == 99
    assert not bool(primary.loc[0, "DirectionConfirmed"])
    assert not secondary["PrimaryGatePass"].any()
    assert not secondary["DirectionConfirmed"].any()


def test_runner_does_not_exist_in_registered_dependency_before_registration(
    tmp_path: Path,
) -> None:
    missing_registration = tmp_path / "missing.json"
    with pytest.raises(FileNotFoundError):
        confirm.build_dependency_manifest(registration_path=missing_registration)


def test_retained_confirmation_passes_registered_and_independent_gates() -> None:
    study = (
        Path(__file__).resolve().parents[1]
        / "validation"
        / "informative_assignment_confirmatory100_20260811"
    )
    metrics = json.loads(
        (study / "aggregate" / "confirmatory_metrics.json").read_text(
            encoding="utf-8"
        )
    )
    primary = pd.read_csv(study / "aggregate" / "primary_result.csv").iloc[0]
    secondary = pd.read_csv(study / "aggregate" / "secondary_results.csv")
    r_verification = pd.read_csv(study / "aggregate" / "r_verification.csv").iloc[0]

    assert metrics["qualification_pass"] is True
    assert all(metrics["gates"].values())
    assert metrics["facets_calibration"]["ready"] == 200
    assert metrics["facets_calibration"]["facets_total_retries"] == 0
    assert primary["EndpointId"] == confirm.PRIMARY_ID
    assert int(primary["FinitePairedReplicates"]) == 100
    assert bool(primary["DirectionConfirmed"])
    assert bool(primary["PrecisionPass"])
    assert np.isclose(primary["MeanContrast"], 0.09863414655507183)
    assert tuple(secondary["EndpointId"]) == confirm.SECONDARY_IDS
    assert secondary["DirectionConfirmed"].all()
    assert bool(r_verification["VerificationPass"])
    assert r_verification["MaximumAbsolutePythonDifference"] <= 1e-12
