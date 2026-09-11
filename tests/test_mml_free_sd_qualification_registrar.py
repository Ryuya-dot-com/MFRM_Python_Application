from __future__ import annotations

from copy import deepcopy
from dataclasses import asdict
import hashlib
import json
from pathlib import Path

import pytest

from mfrm_app.mml_engine_v2 import StationarityContract
from mfrm_app.mml_quadrature_sensitivity import QuadratureSensitivityContract
from mfrm_app.mml_qualification_batch import contract_sha256
from mfrm_app.mml_stationarity import JointPolishOptions
from validation import mml_free_sd_qualification_registrar as registrar


def _digest(label: str) -> str:
    return hashlib.sha256(label.encode("utf-8")).hexdigest()


def _stationarity_contract() -> StationarityContract:
    return StationarityContract(
        max_projected_gradient_supnorm=1e-6,
        max_standardized_score_supnorm=1e-6,
        max_newton_correction_supnorm=1e-6,
        max_restart_improvement_per_observation=1e-10,
        max_restart_displacement=1e-5,
        max_gradient_fd_disagreement=1e-7,
        max_objective_value_disagreement=1e-10,
        max_information_relative_symmetry_residual=1e-6,
        max_information_condition_number=1e6,
        max_constraint_residual=1e-10,
        objective_worsening_tolerance_total=1e-12,
        objective_worsening_tolerance_per_observation=1e-14,
        sigma_boundary_tolerance=1e-6,
    )


def _sensitivity_contract() -> QuadratureSensitivityContract:
    return QuadratureSensitivityContract(
        primary_quadrature_points=31,
        sensitivity_quadrature_points=61,
        max_structural_parameter_difference=0.002,
        max_log_sigma_difference=0.002,
        max_sensitivity_optimization_gain_per_observation=1e-5,
        max_negative_sensitivity_gain_per_observation=1e-10,
        max_abs_primary_back_evaluation_change_per_observation=1e-5,
        max_objective_reconstruction_disagreement=1e-10,
    )


def _plan() -> dict[str, object]:
    stationarity = _stationarity_contract()
    sensitivity = _sensitivity_contract()
    optimizer = JointPolishOptions()
    runs = [
        {
            "run_id": "run-001",
            "seed": 1101,
            "dgp_cell": "neutral",
            "replicate": 1,
            "cohort": "POSITIVE_QUALIFICATION",
            "dgp_spec_sha256": _digest("dgp-neutral"),
            "generation_recipe_sha256": _digest("recipe-1"),
            "expected_outcome": "POSITIVE_CASE",
            "expected_gate": "ALL_REQUIRED_GATES_PASS",
            "q31_required": True,
            "q61_required": True,
            "independent_r_required": True,
            "facets_comparison_required": True,
        },
        {
            "run_id": "run-002",
            "seed": 1102,
            "dgp_cell": "negative-stationarity",
            "replicate": 1,
            "cohort": "NEGATIVE_CONTROL",
            "dgp_spec_sha256": _digest("dgp-negative"),
            "generation_recipe_sha256": _digest("recipe-2"),
            "expected_outcome": "NEGATIVE_CONTROL",
            "expected_gate": "STATIONARITY_PASS",
            "q31_required": True,
            "q61_required": True,
            "independent_r_required": False,
            "facets_comparison_required": False,
        },
        {
            "run_id": "run-003",
            "seed": 1103,
            "dgp_cell": "input-invalid",
            "replicate": 1,
            "cohort": "POSITIVE_QUALIFICATION",
            "dgp_spec_sha256": _digest("dgp-invalid"),
            "generation_recipe_sha256": _digest("recipe-3"),
            "expected_outcome": "STRESS_EXPECTED_PASS",
            "expected_gate": "ALL_REQUIRED_GATES_PASS",
            "q31_required": True,
            "q61_required": False,
            "independent_r_required": False,
            "facets_comparison_required": True,
        },
        {
            "run_id": "run-004",
            "seed": 1104,
            "dgp_cell": "generator-failure",
            "replicate": 1,
            "cohort": "POSITIVE_QUALIFICATION",
            "dgp_spec_sha256": _digest("dgp-failure"),
            "generation_recipe_sha256": _digest("recipe-4"),
            "expected_outcome": "STRESS_EXPECTED_PASS",
            "expected_gate": "ALL_REQUIRED_GATES_PASS",
            "q31_required": True,
            "q61_required": False,
            "independent_r_required": True,
            "facets_comparison_required": False,
        },
    ]
    source_manifest = {
        name: _digest(f"source::{name}") for name in registrar.REQUIRED_SOURCE_FILES
    }
    runtime_manifest = {
        "python": "3.14.0",
        "python_implementation": "CPython",
        "os": "Windows-11",
        "architecture": "AMD64",
        "byteorder": "little",
        "numpy": "2.3.0",
        "scipy": "1.16.0",
        "blas": "OpenBLAS",
        "package_lock_sha256": _digest("package-lock"),
        "r_version": "4.5.1",
        "rscript_sha256": _digest("rscript"),
        "facets_version": "4.5.0",
        "facets_executable_sha256": _digest("facets"),
    }
    numerical = {
        "stationarity": asdict(stationarity),
        "stationarity_sha256": contract_sha256(stationarity),
        "sensitivity": asdict(sensitivity),
        "sensitivity_sha256": contract_sha256(sensitivity),
        "optimizer": asdict(optimizer),
        "optimizer_sha256": contract_sha256(optimizer),
        "objective_evaluator_implementation_sha256": _digest("objective-evaluator"),
    }
    plan = {
        "schema_version": registrar.PLAN_SCHEMA,
        "registration_id": "qualification-dev-001",
        "classification": "NUMERICAL_QUALIFICATION_NOT_SCIENTIFIC_INFERENCE",
        "status": "SEALED_BEFORE_QUALIFICATION_DATA_GENERATION",
        "planned_runs": runs,
        "planned_runs_sha256": registrar._records_digest(
            "MFRM_QUALIFICATION_PLANNED_RUNS", "v1", runs
        ),
        "rng_contract": {
            "algorithm": "PCG64DXSM",
            "implementation": "numpy.random.Generator",
            "version": "2.3.0",
            "stream_partition_rule": "one unique scalar seed per RunId",
        },
        "input_canonicalization_sha256": _digest("input-canonical-v1"),
        "input_validity_rules_sha256": _digest("input-validity-v1"),
        "known_fixture_exclusion_registry_sha256": _digest("exclusions-v1"),
        "negative_control_gate_registry_sha256": (
            registrar.negative_control_gate_registry_sha256()
        ),
        "claim_scope": "FINITE_REGISTERED_CASES_ONLY",
        "fresh_confirmatory_data_required_after_qualification": True,
        "source_manifest": source_manifest,
        "source_manifest_sha256": registrar._domain_sha256(
            "MFRM_QUALIFICATION_SOURCE_MANIFEST", "v1", source_manifest
        ),
        "runtime_manifest": runtime_manifest,
        "runtime_manifest_sha256": registrar._domain_sha256(
            "MFRM_QUALIFICATION_RUNTIME_MANIFEST", "v1", runtime_manifest
        ),
        "numerical_contracts": numerical,
        "numerical_contracts_sha256": registrar._domain_sha256(
            "MFRM_QUALIFICATION_NUMERICAL_CONTRACTS", "v1", numerical
        ),
        "failure_policy": deepcopy(registrar.FAILURE_POLICY),
        "lane_separation": deepcopy(registrar.LANE_SEPARATION),
        "independent_r_plan_sha256": _digest("r-plan"),
        "facets_operational_plan_sha256": _digest("facets-plan"),
        "qualification_data_generated_before_registration": False,
        "estimator_attempts_started_before_registration": False,
        "scientific_endpoints_computed_before_registration": False,
        "scientific_inference_enabled": False,
    }
    return plan


def _inputs(plan: dict[str, object]) -> dict[str, object]:
    plan_identity = registrar.pre_generation_plan_identity_sha256(plan)
    records: list[dict[str, object]] = []
    for run in plan["planned_runs"]:
        run_id = run["run_id"]
        seed = run["seed"]
        if run_id == "run-004":
            records.append(
                {
                    "run_id": run_id,
                    "seed": seed,
                    "status": "GENERATION_FAILED",
                    "raw_input_artifact_sha256": None,
                    "canonical_response_sha256": None,
                    "input_instance_sha256": None,
                    "rows": None,
                    "failure_code": "GENERATOR_EXCEPTION",
                    "generation_error_identity_sha256": _digest("generation-error-4"),
                }
            )
            continue
        raw_digest = _digest(f"raw::{run_id}")
        response_digest = _digest(f"response::{run_id}")
        records.append(
            {
                "run_id": run_id,
                "seed": seed,
                "status": (
                    "INPUT_INVALID_UNDER_PREREGISTERED_RULE"
                    if run_id == "run-003"
                    else "GENERATED"
                ),
                "raw_input_artifact_sha256": raw_digest,
                "canonical_response_sha256": response_digest,
                "input_instance_sha256": registrar.input_instance_sha256(
                    plan_identity_sha256=plan_identity,
                    run_id=run_id,
                    seed=seed,
                    raw_input_artifact_sha256=raw_digest,
                    canonical_response_sha256=response_digest,
                ),
                "rows": 960,
                "failure_code": "INPUT_CATEGORY_INVALID" if run_id == "run-003" else None,
                "generation_error_identity_sha256": (
                    _digest("input-invalid-error-3") if run_id == "run-003" else None
                ),
            }
        )
    return {
        "schema_version": registrar.INPUT_SCHEMA,
        "status": "SEALED_AFTER_GENERATION_BEFORE_FIT",
        "plan_identity_sha256": plan_identity,
        "planned_runs_sha256": plan["planned_runs_sha256"],
        "records": records,
        "records_sha256": registrar._records_digest(
            "MFRM_QUALIFICATION_REALIZED_INPUT_RECORDS", "v1", records
        ),
        "estimator_attempts_started_before_manifest": False,
        "scientific_endpoints_computed_before_manifest": False,
        "scientific_inference_enabled": False,
    }


def _authorization(
    plan: dict[str, object], inputs: dict[str, object]
) -> dict[str, object]:
    records: list[dict[str, object]] = []
    input_by_id = {record["run_id"]: record for record in inputs["records"]}
    for run in plan["planned_runs"]:
        run_id = run["run_id"]
        realized = input_by_id[run_id]
        if realized["status"] == "GENERATED":
            records.append(
                {
                    "run_id": run_id,
                    "seed": run["seed"],
                    "input_instance_sha256": realized["input_instance_sha256"],
                    "status": "AUTHORIZED_FOR_SINGLE_ATTEMPT_DEVELOPMENT_ONLY",
                    # Deliberately equal: instance identity, not likelihood content,
                    # distinguishes independently generated coincident datasets.
                    "likelihood_problem_digest": _digest("same-likelihood-problem"),
                    "estimator_input_config_sha256": _digest(f"config::{run_id}"),
                }
            )
        else:
            records.append(
                {
                    "run_id": run_id,
                    "seed": run["seed"],
                    "input_instance_sha256": None,
                    "status": (
                        "NOT_AUTHORIZED_GENERATION_FAILED"
                        if realized["status"] == "GENERATION_FAILED"
                        else "NOT_AUTHORIZED_INPUT_INVALID"
                    ),
                    "likelihood_problem_digest": None,
                    "estimator_input_config_sha256": None,
                }
            )
    return {
        "schema_version": registrar.FIT_SCHEMA,
        "status": "DEVELOPMENT_ONLY_PREFIT_FACTORY_NOT_QUALIFIED",
        "plan_identity_sha256": registrar.pre_generation_plan_identity_sha256(plan),
        "realized_inputs_identity_sha256": registrar.realized_inputs_identity_sha256(
            plan, inputs
        ),
        "planned_runs_sha256": plan["planned_runs_sha256"],
        "problem_digest_factory_implementation_sha256": _digest("prefit-factory-dev"),
        "records": records,
        "records_sha256": registrar._records_digest(
            "MFRM_QUALIFICATION_FIT_AUTHORIZATION_RECORDS", "v1", records
        ),
        "estimator_attempts_started_before_authorization": False,
        "pure_prefit_problem_factory_qualified": False,
        "scientific_inference_enabled": False,
    }


def _ledger(
    plan: dict[str, object],
    inputs: dict[str, object],
    authorization: dict[str, object],
) -> dict[str, object]:
    input_by_id = {record["run_id"]: record for record in inputs["records"]}
    fit_by_id = {record["run_id"]: record for record in authorization["records"]}
    records: list[dict[str, object]] = []
    for run in plan["planned_runs"]:
        run_id = run["run_id"]
        realized = input_by_id[run_id]
        fit = fit_by_id[run_id]
        base = {
            "run_id": run_id,
            "seed": run["seed"],
            "expected_outcome": run["expected_outcome"],
        }
        if realized["status"] != "GENERATED":
            q61_required = bool(run["q61_required"])
            r_required = bool(run["independent_r_required"])
            facets_required = bool(run["facets_comparison_required"])
            base.update(
                {
                    "attempt_id": f"{run_id}::no-fit",
                    "input_instance_sha256": realized["input_instance_sha256"],
                    "attempt_count": 0,
                    "status": (
                        "NOT_ATTEMPTED_GENERATION_FAILED"
                        if realized["status"] == "GENERATION_FAILED"
                        else "NOT_ATTEMPTED_INPUT_INVALID"
                    ),
                    "observed_outcome": "NOT_OBSERVED_NO_FIT",
                    "observed_gate": None,
                    "python_q31_status": "NOT_OBSERVED",
                    "python_q31_evidence_sha256": None,
                    "python_q31_failed_gates": None,
                    "python_q61_status": "NOT_OBSERVED" if q61_required else "NOT_REQUIRED",
                    "python_q61_evidence_sha256": None,
                    "independent_r_status": "NOT_OBSERVED" if r_required else "NOT_REQUIRED",
                    "independent_r_evidence_sha256": None,
                    "facets_status": "NOT_OBSERVED" if facets_required else "NOT_REQUIRED",
                    "facets_evidence_sha256": None,
                    "mml_numerical_status": "NOT_OBSERVED",
                    "facets_operational_status": (
                        "NOT_OBSERVED" if facets_required else "NOT_REQUIRED"
                    ),
                    "mml_scientific_status": "NOT_SCIENTIFIC_INFERENCE",
                    "warm_start_artifact_sha256": None,
                    "qualification_record_sha256": None,
                    "likelihood_problem_digest": None,
                    "failure_stage": (
                        "GENERATION"
                        if realized["status"] == "GENERATION_FAILED"
                        else "INPUT_VALIDITY"
                    ),
                    "failure_code": realized["failure_code"],
                    "error_identity_sha256": realized["generation_error_identity_sha256"],
                    "stdout_sha256": None,
                    "stderr_sha256": None,
                    "started_utc": None,
                    "finished_utc": None,
                }
            )
        else:
            is_negative = run["expected_outcome"] == "NEGATIVE_CONTROL"
            q31_status = "FAIL" if is_negative else "PASS"
            q61_required = bool(run["q61_required"])
            r_required = bool(run["independent_r_required"])
            facets_required = bool(run["facets_comparison_required"])
            base.update(
                {
                    "attempt_id": f"{run_id}::attempt-01",
                    "input_instance_sha256": realized["input_instance_sha256"],
                    "attempt_count": 1,
                    "status": "ASSESSMENT_RETURNED",
                    "observed_outcome": (
                        "EXPECTED_GATE_FAILED"
                        if is_negative
                        else "MML_REQUIRED_GATES_PASS"
                    ),
                    "observed_gate": run["expected_gate"] if is_negative else None,
                    "python_q31_status": q31_status,
                    "python_q31_evidence_sha256": _digest(f"q31::{run_id}"),
                    "python_q31_failed_gates": (
                        ["STATIONARITY_PASS"] if is_negative else []
                    ),
                    "python_q61_status": "PASS" if q61_required else "NOT_REQUIRED",
                    "python_q61_evidence_sha256": (
                        _digest(f"q61::{run_id}") if q61_required else None
                    ),
                    "independent_r_status": "PASS" if r_required else "NOT_REQUIRED",
                    "independent_r_evidence_sha256": (
                        _digest(f"r::{run_id}") if r_required else None
                    ),
                    "facets_status": "PASS" if facets_required else "NOT_REQUIRED",
                    "facets_evidence_sha256": (
                        _digest(f"facets::{run_id}") if facets_required else None
                    ),
                    "mml_numerical_status": "FAIL" if is_negative else "PASS",
                    "facets_operational_status": (
                        "PASS" if facets_required else "NOT_REQUIRED"
                    ),
                    "mml_scientific_status": "NOT_SCIENTIFIC_INFERENCE",
                    "warm_start_artifact_sha256": _digest(f"warm::{run_id}"),
                    "qualification_record_sha256": _digest(f"record::{run_id}"),
                    "likelihood_problem_digest": fit["likelihood_problem_digest"],
                    "failure_stage": None,
                    "failure_code": None,
                    "error_identity_sha256": None,
                    "stdout_sha256": _digest(f"stdout::{run_id}"),
                    "stderr_sha256": _digest(f"stderr::{run_id}"),
                    "started_utc": "2026-08-12T00:00:00Z",
                    "finished_utc": "2026-08-12T00:00:01.25Z",
                }
            )
        records.append(base)
    return {
        "schema_version": registrar.LEDGER_SCHEMA,
        "status": "COMPLETE_DEVELOPMENT_LEDGER_NOT_A_QUALIFICATION",
        "plan_identity_sha256": registrar.pre_generation_plan_identity_sha256(plan),
        "realized_inputs_identity_sha256": registrar.realized_inputs_identity_sha256(
            plan, inputs
        ),
        "fit_authorization_identity_sha256": registrar.fit_authorization_identity_sha256(
            plan, inputs, authorization
        ),
        "planned_runs_sha256": plan["planned_runs_sha256"],
        "records": records,
        "records_sha256": registrar._records_digest(
            "MFRM_QUALIFICATION_ATTEMPT_LEDGER_RECORDS", "v1", records
        ),
        "replacements_performed": False,
        "optional_extensions_performed": False,
        "scientific_inference_enabled": False,
    }


def _artifacts():
    plan = _plan()
    inputs = _inputs(plan)
    authorization = _authorization(plan, inputs)
    ledger = _ledger(plan, inputs, authorization)
    return plan, inputs, authorization, ledger


def _rehash_plan(plan: dict[str, object]) -> None:
    plan["planned_runs_sha256"] = registrar._records_digest(
        "MFRM_QUALIFICATION_PLANNED_RUNS", "v1", plan["planned_runs"]
    )


def _rehash_inputs(inputs: dict[str, object]) -> None:
    inputs["records_sha256"] = registrar._records_digest(
        "MFRM_QUALIFICATION_REALIZED_INPUT_RECORDS", "v1", inputs["records"]
    )


def _rehash_fit(authorization: dict[str, object]) -> None:
    authorization["records_sha256"] = registrar._records_digest(
        "MFRM_QUALIFICATION_FIT_AUTHORIZATION_RECORDS",
        "v1",
        authorization["records"],
    )


def _rehash_ledger(ledger: dict[str, object]) -> None:
    ledger["records_sha256"] = registrar._records_digest(
        "MFRM_QUALIFICATION_ATTEMPT_LEDGER_RECORDS", "v1", ledger["records"]
    )


def test_valid_four_stage_chain_is_development_only_and_keeps_denominator() -> None:
    plan, inputs, authorization, ledger = _artifacts()
    registrar.validate_pre_generation_plan(plan)
    registrar.validate_realized_inputs(plan, inputs)
    registrar.validate_fit_authorization(plan, inputs, authorization)
    registrar.validate_attempt_ledger(plan, inputs, authorization, ledger)

    assert len(inputs["records"]) == len(plan["planned_runs"]) == 4
    assert len(ledger["records"]) == 4
    generated_problem_digests = [
        item["likelihood_problem_digest"]
        for item in authorization["records"]
        if item["likelihood_problem_digest"] is not None
    ]
    assert len(generated_problem_digests) == 2
    assert len(set(generated_problem_digests)) == 1
    assert authorization["pure_prefit_problem_factory_qualified"] is False
    assert ledger["scientific_inference_enabled"] is False


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (lambda plan: plan["planned_runs"][1].update(seed=1101), "seeds must be unique"),
        (
            lambda plan: plan["planned_runs"][1].update(
                cohort="POSITIVE_QUALIFICATION"
            ),
            "negative-control outcome",
        ),
        (
            lambda plan: [
                run.update(q61_required=False) for run in plan["planned_runs"]
            ],
            "subsets must be nonempty",
        ),
        (
            lambda plan: plan.update(
                qualification_data_generated_before_registration=True
            ),
            "pre-generation gate must be false",
        ),
    ],
)
def test_plan_faults_are_rejected(mutation, message: str) -> None:
    plan = _plan()
    mutation(plan)
    _rehash_plan(plan)
    with pytest.raises(ValueError, match=message):
        registrar.validate_pre_generation_plan(plan)


def test_plan_rejects_tampered_numerical_contract_digest() -> None:
    plan = _plan()
    plan["numerical_contracts"]["stationarity_sha256"] = _digest("forged")
    plan["numerical_contracts_sha256"] = registrar._domain_sha256(
        "MFRM_QUALIFICATION_NUMERICAL_CONTRACTS",
        "v1",
        plan["numerical_contracts"],
    )
    with pytest.raises(ValueError, match="contract digest does not reconstruct"):
        registrar.validate_pre_generation_plan(plan)


def test_plan_rejects_free_text_negative_control_gate() -> None:
    plan = _plan()
    plan["planned_runs"][1]["expected_gate"] = "TOTALLY_FAKE_POSTHOC_GATE"
    _rehash_plan(plan)
    with pytest.raises(ValueError, match="registered gate"):
        registrar.validate_pre_generation_plan(plan)

    plan = _plan()
    plan["planned_runs"][1]["expected_gate"] = "INDEPENDENT_R_PASS"
    _rehash_plan(plan)
    with pytest.raises(ValueError, match="requires its evidence lane"):
        registrar.validate_pre_generation_plan(plan)


def test_realized_inputs_reject_copied_instance_and_untyped_failure() -> None:
    plan = _plan()
    inputs = _inputs(plan)
    inputs["records"][1]["input_instance_sha256"] = inputs["records"][0][
        "input_instance_sha256"
    ]
    _rehash_inputs(inputs)
    with pytest.raises(ValueError, match="does not reconstruct"):
        registrar.validate_realized_inputs(plan, inputs)

    inputs = _inputs(plan)
    inputs["records"][2]["failure_code"] = "free-form posthoc excuse"
    _rehash_inputs(inputs)
    with pytest.raises(ValueError, match="unsupported input-validity failure code"):
        registrar.validate_realized_inputs(plan, inputs)

    inputs = _inputs(plan)
    inputs["records"][0]["seed"] = True
    _rehash_inputs(inputs)
    with pytest.raises(ValueError, match="must be an integer"):
        registrar.validate_realized_inputs(plan, inputs)


def test_fit_authorization_cannot_authorize_invalid_or_claim_prefit_qualification() -> None:
    plan = _plan()
    inputs = _inputs(plan)
    authorization = _authorization(plan, inputs)
    authorization["records"][2].update(
        status="AUTHORIZED_FOR_SINGLE_ATTEMPT_DEVELOPMENT_ONLY",
        input_instance_sha256=inputs["records"][2]["input_instance_sha256"],
        likelihood_problem_digest=_digest("invalid-problem"),
        estimator_input_config_sha256=_digest("invalid-config"),
    )
    _rehash_fit(authorization)
    with pytest.raises(ValueError, match="incorrectly authorized"):
        registrar.validate_fit_authorization(plan, inputs, authorization)

    authorization = _authorization(plan, inputs)
    authorization["pure_prefit_problem_factory_qualified"] = True
    with pytest.raises(ValueError, match="not yet qualified"):
        registrar.validate_fit_authorization(plan, inputs, authorization)


def test_attempt_ledger_rejects_retry_relabel_and_wrong_negative_gate() -> None:
    plan, inputs, authorization, ledger = _artifacts()
    ledger["records"][0]["attempt_count"] = 2
    _rehash_ledger(ledger)
    with pytest.raises(ValueError, match="exactly one fit attempt"):
        registrar.validate_attempt_ledger(plan, inputs, authorization, ledger)

    plan, inputs, authorization, ledger = _artifacts()
    ledger["records"][0]["run_id"] = "run-002"
    _rehash_ledger(ledger)
    with pytest.raises(ValueError, match="exactly follow planned RunIds"):
        registrar.validate_attempt_ledger(plan, inputs, authorization, ledger)

    plan, inputs, authorization, ledger = _artifacts()
    ledger["records"][1]["observed_gate"] = "QUADRATURE_SENSITIVITY_PASS"
    _rehash_ledger(ledger)
    with pytest.raises(ValueError, match="outcome does not reconstruct"):
        registrar.validate_attempt_ledger(plan, inputs, authorization, ledger)

    plan, inputs, authorization, ledger = _artifacts()
    ledger["records"][0]["attempt_count"] = True
    _rehash_ledger(ledger)
    with pytest.raises(ValueError, match="must be an integer"):
        registrar.validate_attempt_ledger(plan, inputs, authorization, ledger)


def test_attempt_ledger_rejects_replacement_and_instance_drift() -> None:
    plan, inputs, authorization, ledger = _artifacts()
    ledger["replacements_performed"] = True
    with pytest.raises(ValueError, match="gate must be false"):
        registrar.validate_attempt_ledger(plan, inputs, authorization, ledger)

    plan, inputs, authorization, ledger = _artifacts()
    ledger["records"][0]["input_instance_sha256"] = _digest("replacement-instance")
    _rehash_ledger(ledger)
    with pytest.raises(ValueError, match="input instance differs"):
        registrar.validate_attempt_ledger(plan, inputs, authorization, ledger)


def test_attempt_ledger_requires_r_evidence_but_separates_facets_status() -> None:
    plan, inputs, authorization, ledger = _artifacts()
    ledger["records"][0]["independent_r_status"] = "NOT_OBSERVED"
    ledger["records"][0]["independent_r_evidence_sha256"] = None
    _rehash_ledger(ledger)
    with pytest.raises(ValueError, match="required independent_r status"):
        registrar.validate_attempt_ledger(plan, inputs, authorization, ledger)

    plan, inputs, authorization, ledger = _artifacts()
    ledger["records"][0]["facets_status"] = "FAIL"
    ledger["records"][0]["facets_operational_status"] = "FAIL"
    _rehash_ledger(ledger)
    registrar.validate_attempt_ledger(plan, inputs, authorization, ledger)
    assert ledger["records"][0]["mml_numerical_status"] == "PASS"

    ledger["records"][0]["observed_outcome"] = "NOT_QUALIFIED"
    ledger["records"][0]["observed_gate"] = "FACETS_OPERATIONAL_PASS"
    _rehash_ledger(ledger)
    with pytest.raises(ValueError, match="passing lane evidence"):
        registrar.validate_attempt_ledger(plan, inputs, authorization, ledger)


def test_negative_control_q31_subgate_must_match_registered_expectation() -> None:
    plan, inputs, authorization, ledger = _artifacts()
    ledger["records"][1]["python_q31_failed_gates"] = ["CONSTRAINT_PASS"]
    _rehash_ledger(ledger)
    with pytest.raises(ValueError, match="outcome does not reconstruct"):
        registrar.validate_attempt_ledger(plan, inputs, authorization, ledger)

    plan, inputs, authorization, ledger = _artifacts()
    ledger["records"][1]["python_q31_failed_gates"] = [
        "CONSTRAINT_PASS",
        "STATIONARITY_PASS",
    ]
    ledger["records"][1]["observed_outcome"] = (
        "EXPECTED_AND_OTHER_GATES_FAILED"
    )
    ledger["records"][1]["observed_gate"] = "CONSTRAINT_PASS"
    _rehash_ledger(ledger)
    registrar.validate_attempt_ledger(plan, inputs, authorization, ledger)

    ledger["records"][1]["observed_outcome"] = "EXPECTED_GATE_FAILED"
    ledger["records"][1]["observed_gate"] = "STATIONARITY_PASS"
    _rehash_ledger(ledger)
    with pytest.raises(ValueError, match="outcome does not reconstruct"):
        registrar.validate_attempt_ledger(plan, inputs, authorization, ledger)


def test_negative_control_cannot_claim_unexpected_pass_when_a_gate_failed() -> None:
    plan, inputs, authorization, ledger = _artifacts()
    ledger["records"][1]["observed_outcome"] = "UNEXPECTED_PASS"
    ledger["records"][1]["observed_gate"] = None
    _rehash_ledger(ledger)
    with pytest.raises(ValueError, match="outcome does not reconstruct"):
        registrar.validate_attempt_ledger(plan, inputs, authorization, ledger)


def test_precheck_failure_cannot_retain_warm_start() -> None:
    plan, inputs, authorization, ledger = _artifacts()
    record = ledger["records"][0]
    record.update(
        status="EXECUTION_FAILED",
        observed_outcome="NOT_OBSERVED_EXECUTION_FAILED",
        observed_gate=None,
        python_q31_status="NOT_OBSERVED",
        python_q31_evidence_sha256=None,
        python_q31_failed_gates=None,
        python_q61_status="NOT_OBSERVED",
        python_q61_evidence_sha256=None,
        independent_r_status="NOT_OBSERVED",
        independent_r_evidence_sha256=None,
        facets_status="NOT_OBSERVED",
        facets_evidence_sha256=None,
        mml_numerical_status="NOT_OBSERVED",
        facets_operational_status="NOT_OBSERVED",
        qualification_record_sha256=None,
        failure_stage="PRECHECK",
        failure_code="PRECHECK_IDENTITY_MISMATCH",
        error_identity_sha256=_digest("precheck-error"),
    )
    _rehash_ledger(ledger)
    with pytest.raises(ValueError, match="retained a warm start"):
        registrar.validate_attempt_ledger(plan, inputs, authorization, ledger)


def test_publish_is_identity_last_and_readback_is_strict(tmp_path: Path) -> None:
    plan = _plan()
    output = tmp_path / "t0"
    published = registrar.publish_sealed_artifact(
        output_dir=output,
        document=plan,
        validator=registrar.validate_pre_generation_plan,
        semantic_identity_function=registrar.pre_generation_plan_identity_sha256,
    )
    assert set(path.name for path in output.iterdir()) == {
        "artifact.json",
        "identity.json",
    }
    assert published["identity"]["scientific_inference_ready"] is False
    assert published["identity"]["artifact_bytes"] == len(
        (output / "artifact.json").read_bytes()
    )

    identity_text = (output / "identity.json").read_text(encoding="utf-8")
    duplicate = identity_text.replace(
        '{"artifact_bytes":',
        '{"artifact_bytes":1,"artifact_bytes":',
        1,
    )
    (output / "identity.json").write_text(duplicate, encoding="utf-8")
    with pytest.raises(ValueError, match="duplicate key"):
        registrar.validate_published_artifact(
            output,
            validator=registrar.validate_pre_generation_plan,
            semantic_identity_function=registrar.pre_generation_plan_identity_sha256,
        )


def test_published_identity_rejects_forged_parent_and_publisher(tmp_path: Path) -> None:
    for field in ("parent", "publisher"):
        output = tmp_path / field
        registrar.publish_sealed_artifact(
            output_dir=output,
            document=_plan(),
            validator=registrar.validate_pre_generation_plan,
            semantic_identity_function=registrar.pre_generation_plan_identity_sha256,
        )
        identity = json.loads((output / "identity.json").read_text(encoding="utf-8"))
        if field == "parent":
            identity["parent_identity_sha256"] = {"invented": _digest("invented")}
            message = "parent identities do not reconstruct"
        else:
            identity["publisher_sha256"] = _digest("forged-publisher")
            message = "publisher source does not reconstruct"
        (output / "identity.json").write_text(
            json.dumps(identity, separators=(",", ":")) + "\n",
            encoding="utf-8",
        )
        with pytest.raises(ValueError, match=message):
            registrar.validate_published_artifact(
                output,
                validator=registrar.validate_pre_generation_plan,
                semantic_identity_function=registrar.pre_generation_plan_identity_sha256,
            )


def test_published_validation_detects_concurrent_artifact_mutation(
    tmp_path: Path,
) -> None:
    output = tmp_path / "toctou"
    registrar.publish_sealed_artifact(
        output_dir=output,
        document=_plan(),
        validator=registrar.validate_pre_generation_plan,
        semantic_identity_function=registrar.pre_generation_plan_identity_sha256,
    )

    def mutating_identity(artifact):
        identity = registrar.pre_generation_plan_identity_sha256(artifact)
        (output / "artifact.json").write_text("{}\n", encoding="utf-8")
        return identity

    with pytest.raises(ValueError, match="changed during validation"):
        registrar.validate_published_artifact(
            output,
            validator=registrar.validate_pre_generation_plan,
            semantic_identity_function=mutating_identity,
        )


def test_publish_failure_cleans_owned_partial_directory(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    output = tmp_path / "faulted"
    real_replace = registrar.os.replace
    calls = 0

    def fail_second_replace(source, destination):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise OSError("injected second replace failure")
        return real_replace(source, destination)

    monkeypatch.setattr(registrar.os, "replace", fail_second_replace)
    with pytest.raises(OSError, match="injected second replace failure"):
        registrar.publish_sealed_artifact(
            output_dir=output,
            document=_plan(),
            validator=registrar.validate_pre_generation_plan,
            semantic_identity_function=registrar.pre_generation_plan_identity_sha256,
        )
    assert not output.exists()


def test_published_artifact_rejects_nonfinite_json_token(tmp_path: Path) -> None:
    output = tmp_path / "nonfinite"
    registrar.publish_sealed_artifact(
        output_dir=output,
        document=_plan(),
        validator=registrar.validate_pre_generation_plan,
        semantic_identity_function=registrar.pre_generation_plan_identity_sha256,
    )
    artifact_text = (output / "artifact.json").read_text(encoding="utf-8")
    artifact = json.loads(artifact_text)
    artifact["qualification_data_generated_before_registration"] = float("nan")
    (output / "artifact.json").write_text(
        json.dumps(artifact, allow_nan=True), encoding="utf-8"
    )
    with pytest.raises(ValueError, match="non-finite token"):
        registrar.validate_published_artifact(
            output,
            validator=registrar.validate_pre_generation_plan,
            semantic_identity_function=registrar.pre_generation_plan_identity_sha256,
        )

    with pytest.raises(ValueError, match="non-finite number"):
        registrar._load_strict_json_object(b'{"x":1e999}', "overflow probe")
