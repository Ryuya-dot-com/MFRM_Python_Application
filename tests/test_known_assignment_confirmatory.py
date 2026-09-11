from __future__ import annotations

import copy
import json

import numpy as np
import pandas as pd
import pytest

from validation import known_assignment_confirmatory as confirm
from validation import known_assignment_mml_crossfit as crossfit


def test_plan_freezes_triplet_unit_precision_and_two_layer_statuses() -> None:
    plan = confirm.validate_plan()
    assert plan["evidence_status"]["independent_person_vector_triplets"] == 200
    assert plan["evidence_status"]["datasets"] == 600
    assert plan["evidence_status"]["optional_stopping"] is False
    assert plan["primary_endpoint"]["id"] == confirm.PRIMARY_ID
    assert plan["primary_endpoint"]["sesoi"].startswith("not defined")
    assert plan["sample_size_contract"]["target_half_width"] == 0.015
    assert tuple(
        plan["execution_lanes"]["facets_external_calibration"]["triplet_subset"]
    ) == confirm.CALIBRATION_TRIPLETS
    assert "ScientificInferenceValidated" in plan["two_layer_terminal_status"]
    assert "FACETSComplementQualification" in plan["two_layer_terminal_status"]
    assert "two-decimal" in plan["facets_contract"]["fit_precision"]


@pytest.mark.parametrize(
    ("path", "replacement"),
    [
        (("evidence_status", "independent_person_vector_triplets"), 201),
        (("evidence_status", "optional_stopping"), True),
        (("evidence_status", "replacement_triplets"), True),
        (("primary_endpoint", "alpha"), 0.025),
        (("sample_size_contract", "target_half_width"), 0.02),
    ],
)
def test_plan_validation_rejects_scientific_drift(
    path: tuple[str, str], replacement: object
) -> None:
    plan = json.loads(confirm.PLAN_PATH.read_text(encoding="utf-8"))
    plan[path[0]][path[1]] = replacement
    with pytest.raises(ValueError):
        confirm.validate_plan(plan)


def test_runtime_seed_expansion_is_complete_disjoint_and_formula_bound() -> None:
    runtime = confirm.runtime_plan(confirm.validate_plan())
    dgm = runtime["data_generating_process"]
    assignment = runtime["assignment"]
    assert dgm["person_vector_ids"] == list(range(1, 201))
    assert dgm["person_vector_seeds"] == list(range(26100001, 26100201))
    assert dgm["response_uniform_seeds"] == list(range(26101001, 26101201))
    assert assignment["assignment_seeds_by_person_vector"]["1"] == {
        "-0.8": 26102001,
        "0.0": 26102002,
        "0.8": 26102003,
    }
    assert assignment["assignment_seeds_by_person_vector"]["200"]["0.8"] == 26102600


def test_attempt_filter_keeps_all_mml_and_only_frozen_external_subset() -> None:
    rows = []
    ordinal = 0
    for vector in range(1, 201):
        for gamma in (-0.8, 0.0, 0.8):
            run_id = f"v{vector:03d}g{gamma:+.1f}"
            for attempt_type in confirm.engine.ATTEMPT_TYPES:
                rows.append(
                    {
                        "AttemptOrdinal": ordinal,
                        "AttemptId": f"{run_id}::{attempt_type}",
                        "RunId": run_id,
                        "AttemptType": attempt_type,
                        "PersonVector": vector,
                        "Gamma": gamma,
                    }
                )
                ordinal += 1
    output = confirm.filter_attempt_manifest(pd.DataFrame(rows))
    assert len(output) == 1272
    counts = output["AttemptType"].value_counts()
    assert counts[confirm.MML_ATTEMPT_TYPES[0]] == 600
    assert counts[confirm.MML_ATTEMPT_TYPES[1]] == 600
    assert counts[confirm.PAIR_ATTEMPT_TYPE] == 36
    assert counts[confirm.CMLE_ATTEMPT_TYPE] == 36
    external = output[output["AttemptType"].isin((confirm.PAIR_ATTEMPT_TYPE, confirm.CMLE_ATTEMPT_TYPE))]
    assert tuple(sorted(external["PersonVector"].unique())) == confirm.CALIBRATION_TRIPLETS


def _synthetic_inputs() -> tuple[pd.DataFrame, pd.DataFrame]:
    truth = {"R01": -0.45, "R02": -0.15, "R03": 0.15, "R04": 0.45}
    errors = []
    runs = []
    for vector in range(1, 201):
        jitter = (vector - 100.5) / 100000.0
        for mode in (confirm.MML_FREE_MODE, confirm.MML_FIXED_MODE):
            for gamma in (-0.8, 0.0, 0.8):
                if gamma == 0.0:
                    base = 0.10 + jitter
                    slope = 0.0
                else:
                    base = (0.20 if mode == confirm.MML_FREE_MODE else 0.18) + jitter
                    slope = 0.4 if gamma < 0 else -0.4
                for level, severity in truth.items():
                    errors.append(
                        {
                            "PersonVector": vector,
                            "Gamma": gamma,
                            "EstimatorMode": mode,
                            "Facet": "Rater",
                            "Level": level,
                            "ErrorAligned": base + slope * severity,
                            "IncludedInStudy": True,
                        }
                    )
                runs.append(
                    {
                        "PersonVector": vector,
                        "Gamma": gamma,
                        "EstimatorMode": mode,
                        "EstimatedPopulationSD": (
                            0.8 if gamma == 0.0 else 0.7
                        )
                        if mode == confirm.MML_FREE_MODE
                        else np.nan,
                    }
                )
    return pd.DataFrame(errors), pd.DataFrame(runs)


def test_registered_endpoints_use_two_hundred_triplets_and_holm_gatekeeping() -> None:
    errors, runs = _synthetic_inputs()
    contrasts = confirm.build_endpoint_contrasts(rater_errors=errors, runs=runs)
    assert tuple(contrasts["EndpointId"].drop_duplicates()) == (
        confirm.PRIMARY_ID,
        *confirm.SECONDARY_IDS,
    )
    assert contrasts.groupby("EndpointId").size().eq(200).all()
    primary, secondary = confirm.summarize_endpoints(contrasts)
    assert int(primary.loc[0, "FiniteTriplets"]) == 200
    assert bool(primary.loc[0, "DirectionConfirmed"])
    assert bool(primary.loc[0, "PrecisionQualified"])
    assert secondary["PrimaryGatePass"].all()
    assert secondary["DirectionConfirmed"].all()
    assert (secondary["HolmAdjustedP"] >= secondary["RawOneSidedP"]).all()


def test_one_missing_gamma_makes_primary_inconclusive_without_replacement() -> None:
    errors, runs = _synthetic_inputs()
    missing = ~(
        errors["PersonVector"].eq(200)
        & errors["Gamma"].eq(0.8)
        & errors["EstimatorMode"].eq(confirm.MML_FREE_MODE)
    )
    contrasts = confirm.build_endpoint_contrasts(rater_errors=errors.loc[missing], runs=runs)
    primary, secondary = confirm.summarize_endpoints(contrasts)
    assert int(primary.loc[0, "FiniteTriplets"]) == 199
    assert not bool(primary.loc[0, "FullTripletGate"])
    assert not bool(primary.loc[0, "DirectionConfirmed"])
    assert not secondary["PrimaryGatePass"].any()
    assert not secondary["DirectionConfirmed"].any()


def test_bias_mse_decomposition_is_distinct_and_satisfies_identity() -> None:
    errors, _ = _synthetic_inputs()
    output = confirm.bias_mse_decomposition(errors)
    assert len(output) == 24
    assert output["FullTripletGate"].all()
    assert output["InferenceStatus"].eq("descriptive_registered_decomposition").all()
    assert output["MSEIdentityResidual"].abs().max() <= 1e-14


def test_facets_path_budget_is_checked_on_short_registered_root(tmp_path) -> None:
    lengths = confirm._worst_case_facets_paths(tmp_path / "kac")
    assert set(lengths) == {"analysis.txt", "report_u6.txt", "scores_u6.txt"}
    assert max(lengths.values()) <= 220


def test_engine_context_restores_frozen_preflight_module(tmp_path) -> None:
    originals = {
        "SCHEMA_VERSION": confirm.engine.SCHEMA_VERSION,
        "PLAN_PATH": confirm.engine.PLAN_PATH,
        "STUDY_DIR": confirm.engine.STUDY_DIR,
        "validate_registration": confirm.engine.validate_registration,
    }
    with confirm._engine_context(tmp_path / "study"):
        assert confirm.engine.SCHEMA_VERSION == confirm.SCHEMA_VERSION
        assert confirm.engine.PLAN_PATH == confirm.PLAN_PATH
        assert confirm.engine.STUDY_DIR == (tmp_path / "study").resolve()
        assert confirm.engine.validate_registration is confirm.validate_registration
    for name, value in originals.items():
        assert getattr(confirm.engine, name) is value


def test_independent_python_protocol_evaluator_reproduces_retained_mml_loglik() -> None:
    study = (
        confirm.ROOT
        / "validation"
        / "known_assignment_multivector_preflight4_20260811"
    )
    run_id = "gamma_neg_0p8__vector-01"
    run = pd.read_csv(study / "work" / "00002" / "run_ledger.csv").iloc[0]
    recovery = pd.read_csv(study / "work" / "00002" / "recovery.csv")
    thresholds = pd.read_csv(study / "work" / "00002" / "thresholds.csv")
    parameter_rows = []
    for block in ("Rater", "Task", "Criterion"):
        for row in recovery.loc[recovery["Facet"].eq(block)].itertuples(index=False):
            parameter_rows.append(
                {"Block": block, "Level": str(row.Level), "PythonEstimate": row.Estimate}
            )
    for row in thresholds.itertuples(index=False):
        parameter_rows.append(
            {
                "Block": "Step",
                "Level": f"{row.StepFacetLevel}::{int(row.Category)}",
                "PythonEstimate": row.Estimate,
            }
        )
    surface = crossfit.parameter_surface(pd.DataFrame(parameter_rows), "PythonEstimate")
    ratings = pd.read_csv(study / "retained_input" / "generated_ratings.csv")
    ratings = ratings.loc[ratings["RunId"].eq(run_id)]
    evaluated = crossfit.marginal_loglik(
        ratings,
        surface,
        sigma=float(run["EstimatedPopulationSD"]),
        points=31,
    )
    assert abs(evaluated - float(run["LogLik"])) <= 1e-9


def test_base_r_crossfit_qualification_passes_twelve_nonconfirmatory_fixtures() -> None:
    assessment_path = (
        confirm.ROOT
        / "validation"
        / "known_assignment_mml_crossfit_preflight12_v2_20260811"
        / "assessment.json"
    )
    assessment = json.loads(assessment_path.read_text(encoding="utf-8"))
    assert assessment["pass"] is True
    assert assessment["datasets"] == 12
    assert assessment["r_version"] == "4.5.1"
    assert all(assessment["gates"].values())
    assert assessment["metrics"]["maximum_cross_language_loglik_difference"] <= 1e-9
    assert assessment["metrics"]["maximum_abs_r_q31_minus_python_parameter"] <= 0.01
    assert assessment["metrics"]["maximum_abs_r_q61_minus_q31_parameter"] <= 0.01
