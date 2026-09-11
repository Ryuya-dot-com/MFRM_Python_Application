"""Fail-closed schemas for a future prospective free-SD MML qualification.

The registrar is deliberately four-stage:

T0 pre-generation plan -> T1 realized inputs -> T1.5 fit authorization ->
T2 attempt ledger.

The split prevents a fitted result from defining its own planned denominator or
problem identity.  This module validates and publishes development artifacts;
it never enables a registered scope or scientific inference.  T1.5 can become
prospective only after a pure pre-fit likelihood-problem factory is separately
qualified and hash-frozen.
"""

from __future__ import annotations

from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import time
from typing import Any, Callable, Mapping, Sequence

from mfrm_app.mml_engine_v2 import StationarityContract
from mfrm_app.mml_quadrature_sensitivity import QuadratureSensitivityContract
from mfrm_app.mml_qualification_batch import contract_sha256
from mfrm_app.mml_stationarity import JointPolishOptions


PLAN_SCHEMA = "mml-free-sd-qualification-pre-generation-plan-v1"
INPUT_SCHEMA = "mml-free-sd-qualification-realized-inputs-v1"
FIT_SCHEMA = "mml-free-sd-qualification-fit-authorization-v1"
LEDGER_SCHEMA = "mml-free-sd-qualification-attempt-ledger-v1"
IDENTITY_SCHEMA = "mml-free-sd-qualification-artifact-identity-v1"
SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")
UTC_PATTERN = re.compile(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(?:\.\d+)?Z$")

PLAN_KEYS = {
    "schema_version",
    "registration_id",
    "classification",
    "status",
    "planned_runs",
    "planned_runs_sha256",
    "rng_contract",
    "input_canonicalization_sha256",
    "input_validity_rules_sha256",
    "known_fixture_exclusion_registry_sha256",
    "negative_control_gate_registry_sha256",
    "claim_scope",
    "fresh_confirmatory_data_required_after_qualification",
    "source_manifest",
    "source_manifest_sha256",
    "runtime_manifest",
    "runtime_manifest_sha256",
    "numerical_contracts",
    "numerical_contracts_sha256",
    "failure_policy",
    "lane_separation",
    "independent_r_plan_sha256",
    "facets_operational_plan_sha256",
    "qualification_data_generated_before_registration",
    "estimator_attempts_started_before_registration",
    "scientific_endpoints_computed_before_registration",
    "scientific_inference_enabled",
}
PLANNED_RUN_KEYS = {
    "run_id",
    "seed",
    "dgp_cell",
    "replicate",
    "cohort",
    "dgp_spec_sha256",
    "generation_recipe_sha256",
    "expected_outcome",
    "expected_gate",
    "q31_required",
    "q61_required",
    "independent_r_required",
    "facets_comparison_required",
}
RUNTIME_KEYS = {
    "python",
    "python_implementation",
    "os",
    "architecture",
    "byteorder",
    "numpy",
    "scipy",
    "blas",
    "package_lock_sha256",
    "r_version",
    "rscript_sha256",
    "facets_version",
    "facets_executable_sha256",
}
RNG_KEYS = {"algorithm", "implementation", "version", "stream_partition_rule"}
NUMERICAL_CONTRACT_KEYS = {
    "stationarity",
    "stationarity_sha256",
    "sensitivity",
    "sensitivity_sha256",
    "optimizer",
    "optimizer_sha256",
    "objective_evaluator_implementation_sha256",
}
FAILURE_POLICY = {
    "all_planned_required": True,
    "replacement_forbidden": True,
    "optional_extension_forbidden": True,
    "numerical_refit_after_failure_forbidden": True,
    "maximum_fit_attempts_per_generated_run": 1,
    "transport_retry_same_input_only": True,
    "infrastructure_transport_retry_limit": 3,
    "retry_adoption_rule": "FIRST_COMPLETE_IDENTITY_VALID_TRANSPORT",
    "failed_run_disposition": "RETAIN_IN_DENOMINATOR_NOT_READY",
}
LANE_SEPARATION = {
    "facets_is_jmle_only": True,
    "facets_failure_does_not_reduce_mml_denominator": True,
    "facets_display_values_not_raw_inputs": True,
    "independent_r_mml_required": True,
    "scientific_status_separate": True,
}

INPUT_KEYS = {
    "schema_version",
    "status",
    "plan_identity_sha256",
    "planned_runs_sha256",
    "records",
    "records_sha256",
    "estimator_attempts_started_before_manifest",
    "scientific_endpoints_computed_before_manifest",
    "scientific_inference_enabled",
}
INPUT_RECORD_KEYS = {
    "run_id",
    "seed",
    "status",
    "raw_input_artifact_sha256",
    "canonical_response_sha256",
    "input_instance_sha256",
    "rows",
    "failure_code",
    "generation_error_identity_sha256",
}

FIT_KEYS = {
    "schema_version",
    "status",
    "plan_identity_sha256",
    "realized_inputs_identity_sha256",
    "planned_runs_sha256",
    "problem_digest_factory_implementation_sha256",
    "records",
    "records_sha256",
    "estimator_attempts_started_before_authorization",
    "pure_prefit_problem_factory_qualified",
    "scientific_inference_enabled",
}
FIT_RECORD_KEYS = {
    "run_id",
    "seed",
    "input_instance_sha256",
    "status",
    "likelihood_problem_digest",
    "estimator_input_config_sha256",
}

LEDGER_KEYS = {
    "schema_version",
    "status",
    "plan_identity_sha256",
    "realized_inputs_identity_sha256",
    "fit_authorization_identity_sha256",
    "planned_runs_sha256",
    "records",
    "records_sha256",
    "replacements_performed",
    "optional_extensions_performed",
    "scientific_inference_enabled",
}
LEDGER_RECORD_KEYS = {
    "run_id",
    "attempt_id",
    "seed",
    "input_instance_sha256",
    "attempt_count",
    "status",
    "expected_outcome",
    "observed_outcome",
    "observed_gate",
    "python_q31_status",
    "python_q31_evidence_sha256",
    "python_q31_failed_gates",
    "python_q61_status",
    "python_q61_evidence_sha256",
    "independent_r_status",
    "independent_r_evidence_sha256",
    "facets_status",
    "facets_evidence_sha256",
    "mml_numerical_status",
    "facets_operational_status",
    "mml_scientific_status",
    "warm_start_artifact_sha256",
    "qualification_record_sha256",
    "likelihood_problem_digest",
    "failure_stage",
    "failure_code",
    "error_identity_sha256",
    "stdout_sha256",
    "stderr_sha256",
    "started_utc",
    "finished_utc",
}

INPUT_INVALID_FAILURE_CODES = {
    "INPUT_SCHEMA_INVALID",
    "INPUT_CATEGORY_INVALID",
    "INPUT_ASSIGNMENT_INVALID",
    "INPUT_NONFINITE",
}
GENERATION_FAILURE_CODES = {
    "GENERATOR_EXCEPTION",
    "GENERATOR_PROTOCOL_ERROR",
    "GENERATOR_SOURCE_RUNTIME_MISMATCH",
}
ATTEMPT_FAILURE_CODES = {
    "PRECHECK": {
        "PRECHECK_IDENTITY_MISMATCH",
        "PRECHECK_SOURCE_RUNTIME_MISMATCH",
        "PRECHECK_INPUT_MISMATCH",
        "PRECHECK_PROBLEM_MISMATCH",
        "PROCESS_LAUNCH_FAILED",
        "PROCESS_TIMEOUT",
        "PROCESS_NONZERO_EXIT",
        "PROCESS_PROTOCOL_ERROR",
    },
    "WARM_START": {"WARM_START_FAILED", "WARM_START_NONFINITE"},
    "PRIMARY_Q": {"PRIMARY_Q_FAILED", "PRIMARY_Q_NONFINITE"},
    "SENSITIVITY_Q": {"SENSITIVITY_Q_FAILED", "SENSITIVITY_Q_NONFINITE"},
    "CROSS_EVALUATION": {
        "CROSS_EVALUATION_FAILED",
        "CROSS_EVALUATION_NONFINITE",
    },
    "ASSESSMENT": {"ASSESSMENT_FAILED", "ASSESSMENT_SCHEMA_INVALID"},
    "PUBLICATION": {"PUBLICATION_TRANSPORT_FAILED", "PUBLICATION_READBACK_FAILED"},
}
NEGATIVE_CONTROL_GATE_LANES = {
    "ALGORITHM_TERMINATED": "python_q31",
    "STATIONARITY_PASS": "python_q31",
    "CONSTRAINT_PASS": "python_q31",
    "OBJECTIVE_CONSISTENCY_PASS": "python_q31",
    "QUADRATURE_SENSITIVITY_PASS": "python_q61",
    "INDEPENDENT_R_PASS": "independent_r",
    "FACETS_OPERATIONAL_PASS": "facets",
}
NEGATIVE_CONTROL_GATES = frozenset(NEGATIVE_CONTROL_GATE_LANES)
NEGATIVE_CONTROL_GATE_ORDER = (
    "ALGORITHM_TERMINATED",
    "STATIONARITY_PASS",
    "CONSTRAINT_PASS",
    "OBJECTIVE_CONSISTENCY_PASS",
    "QUADRATURE_SENSITIVITY_PASS",
    "INDEPENDENT_R_PASS",
    "FACETS_OPERATIONAL_PASS",
)
PYTHON_Q31_GATES = frozenset(
    gate
    for gate, lane in NEGATIVE_CONTROL_GATE_LANES.items()
    if lane == "python_q31"
)

REQUIRED_SOURCE_FILES = {
    "mfrm_app/mml_stationarity.py",
    "mfrm_app/mml_engine_v2.py",
    "mfrm_app/mml_quadrature_sensitivity.py",
    "mfrm_app/mml_qualification_batch.py",
    "validation/mml_free_sd_stationarity_adapter.py",
    "validation/mml_free_sd_quadrature_adapter.py",
    "validation/mml_free_sd_qualification_registrar.py",
}


def _json_bytes(value: object) -> bytes:
    return (
        json.dumps(
            value,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        )
        + "\n"
    ).encode("utf-8")


def _load_strict_json_object(payload: bytes, label: str) -> Mapping[str, Any]:
    """Decode one JSON object while rejecting duplicate keys and non-finite tokens."""

    def pairs_hook(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"{label} contains a duplicate key: {key}")
            result[key] = value
        return result

    def reject_constant(token: str) -> None:
        raise ValueError(f"{label} contains a non-finite token: {token}")

    def finite_float(token: str) -> float:
        value = float(token)
        if not (-float("inf") < value < float("inf")):
            raise ValueError(f"{label} contains a non-finite number: {token}")
        return value

    try:
        decoded = json.loads(
            payload.decode("utf-8"),
            object_pairs_hook=pairs_hook,
            parse_constant=reject_constant,
            parse_float=finite_float,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"{label} is not strict UTF-8 JSON") from exc
    if not isinstance(decoded, Mapping):
        raise ValueError(f"{label} must contain one JSON object")
    return decoded


def _domain_sha256(domain: str, schema: str, payload: object) -> str:
    return hashlib.sha256(
        _json_bytes(
            {
                "domain": domain,
                "schema_version": schema,
                "payload": payload,
            }
        )
    ).hexdigest()


def _exact_keys(value: object, expected: set[str], label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{label} must be an object")
    actual = set(map(str, value.keys()))
    if actual != expected:
        raise ValueError(
            f"{label} keys differ: missing={sorted(expected - actual)}, "
            f"extra={sorted(actual - expected)}"
        )
    return value


def _sha256(value: object, label: str, *, nullable: bool = False) -> str | None:
    if value is None and nullable:
        return None
    if not isinstance(value, str) or not SHA256_PATTERN.fullmatch(value):
        raise ValueError(f"{label} must be a lowercase SHA256 digest")
    return value


def _text(value: object, label: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{label} must be nonempty text")
    return value


def _bool(value: object, label: str) -> bool:
    if not isinstance(value, bool):
        raise ValueError(f"{label} must be boolean")
    return value


def _integer(value: object, label: str, *, minimum: int = 0) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f"{label} must be an integer >= {minimum}")
    return value


def _utc(value: object, label: str, *, nullable: bool = False) -> str | None:
    if value is None and nullable:
        return None
    if not isinstance(value, str) or not UTC_PATTERN.fullmatch(value):
        raise ValueError(f"{label} must be an ISO-8601 UTC timestamp")
    return value


def _utc_datetime(value: object, label: str) -> datetime:
    text = _utc(value, label)
    assert text is not None
    try:
        parsed = datetime.fromisoformat(text[:-1] + "+00:00")
    except ValueError as exc:
        raise ValueError(f"{label} is not a valid UTC timestamp") from exc
    if parsed.tzinfo != timezone.utc:
        raise ValueError(f"{label} must use UTC")
    return parsed


def _records_digest(domain: str, schema: str, records: Sequence[object]) -> str:
    return _domain_sha256(domain, schema, list(records))


def negative_control_gate_registry_sha256() -> str:
    return _domain_sha256(
        "MFRM_QUALIFICATION_NEGATIVE_CONTROL_GATE_REGISTRY",
        "v1",
        {
            "gate_lanes": NEGATIVE_CONTROL_GATE_LANES,
            "wrong_gate_reporting_priority": NEGATIVE_CONTROL_GATE_ORDER,
        },
    )


def _validate_lane_evidence(
    item: Mapping[str, Any],
    *,
    run_id: str,
    lane: str,
    required: bool,
    attempted: bool,
) -> str:
    status_name = f"{lane}_status"
    evidence_name = f"{lane}_evidence_sha256"
    status = item[status_name]
    if not attempted:
        expected = "NOT_OBSERVED" if required else "NOT_REQUIRED"
        if status != expected or item[evidence_name] is not None:
            raise ValueError(f"unattempted {lane} evidence differs for {run_id}")
        return str(status)
    if required:
        if status not in {"PASS", "FAIL"}:
            raise ValueError(f"required {lane} status differs for {run_id}")
        _sha256(item[evidence_name], f"{run_id} {lane} evidence")
    elif status != "NOT_REQUIRED" or item[evidence_name] is not None:
        raise ValueError(f"unrequired {lane} evidence differs for {run_id}")
    return str(status)


def _gate_failed(item: Mapping[str, Any], gate: str) -> bool:
    try:
        lane = NEGATIVE_CONTROL_GATE_LANES[gate]
    except KeyError as exc:
        raise ValueError(f"unregistered observed gate: {gate}") from exc
    if lane == "python_q31":
        return gate in item["python_q31_failed_gates"]
    return item[f"{lane}_status"] == "FAIL"


def _failed_registered_gates(item: Mapping[str, Any]) -> tuple[str, ...]:
    return tuple(gate for gate in NEGATIVE_CONTROL_GATE_ORDER if _gate_failed(item, gate))


def validate_pre_generation_plan(plan: Mapping[str, Any]) -> Mapping[str, Any]:
    value = _exact_keys(plan, PLAN_KEYS, "pre-generation plan")
    if value["schema_version"] != PLAN_SCHEMA:
        raise ValueError("pre-generation plan schema differs")
    _text(value["registration_id"], "registration_id")
    if value["classification"] != "NUMERICAL_QUALIFICATION_NOT_SCIENTIFIC_INFERENCE":
        raise ValueError("pre-generation plan classification differs")
    if value["status"] != "SEALED_BEFORE_QUALIFICATION_DATA_GENERATION":
        raise ValueError("pre-generation plan status differs")

    records = value["planned_runs"]
    if not isinstance(records, list) or not records:
        raise ValueError("planned_runs must be a nonempty list")
    run_ids: list[str] = []
    seeds: list[int] = []
    cells_and_replicates: list[tuple[str, int]] = []
    lane_counts = {"q61": 0, "r": 0, "facets": 0}
    cohort_counts = {"POSITIVE_QUALIFICATION": 0, "NEGATIVE_CONTROL": 0}
    for index, record in enumerate(records):
        item = _exact_keys(record, PLANNED_RUN_KEYS, f"planned_runs[{index}]")
        run_ids.append(_text(item["run_id"], f"planned_runs[{index}].run_id"))
        seeds.append(_integer(item["seed"], f"planned_runs[{index}].seed"))
        dgp_cell = _text(item["dgp_cell"], f"planned_runs[{index}].dgp_cell")
        replicate = _integer(
            item["replicate"],
            f"planned_runs[{index}].replicate",
            minimum=1,
        )
        cells_and_replicates.append((dgp_cell, replicate))
        _sha256(item["dgp_spec_sha256"], f"planned_runs[{index}].dgp_spec_sha256")
        _sha256(
            item["generation_recipe_sha256"],
            f"planned_runs[{index}].generation_recipe_sha256",
        )
        expected = item["expected_outcome"]
        if expected not in {"POSITIVE_CASE", "STRESS_EXPECTED_PASS", "NEGATIVE_CONTROL"}:
            raise ValueError(f"planned expected outcome is unsupported at index {index}")
        expected_gate = _text(
            item["expected_gate"],
            f"planned_runs[{index}].expected_gate",
        )
        cohort = item["cohort"]
        if cohort not in cohort_counts:
            raise ValueError(f"planned cohort is unsupported at index {index}")
        cohort_counts[str(cohort)] += 1
        if expected == "NEGATIVE_CONTROL":
            if cohort != "NEGATIVE_CONTROL":
                raise ValueError("negative-control outcome must use the negative cohort")
            if expected_gate == "ALL_REQUIRED_GATES_PASS":
                raise ValueError("negative control must name its expected failing gate")
            if expected_gate not in NEGATIVE_CONTROL_GATES:
                raise ValueError("negative control must use a registered gate")
        else:
            if cohort != "POSITIVE_QUALIFICATION":
                raise ValueError("positive/stress outcome must use the positive cohort")
            if expected_gate != "ALL_REQUIRED_GATES_PASS":
                raise ValueError("positive/stress cases must require all gates to pass")
        if not _bool(item["q31_required"], f"planned_runs[{index}].q31_required"):
            raise ValueError("every planned case requires the Q31 lane")
        for field, counter in (
            ("q61_required", "q61"),
            ("independent_r_required", "r"),
            ("facets_comparison_required", "facets"),
        ):
            lane_counts[counter] += int(
                _bool(item[field], f"planned_runs[{index}].{field}")
            )
        if expected == "NEGATIVE_CONTROL":
            required_field = {
                "python_q61": "q61_required",
                "independent_r": "independent_r_required",
                "facets": "facets_comparison_required",
            }.get(NEGATIVE_CONTROL_GATE_LANES[expected_gate])
            if required_field is not None and not item[required_field]:
                raise ValueError("negative-control gate requires its evidence lane")
    if run_ids != sorted(run_ids) or len(set(run_ids)) != len(run_ids):
        raise ValueError("planned RunIds must be unique and canonically sorted")
    if len(set(seeds)) != len(seeds):
        raise ValueError("planned generation seeds must be unique")
    if len(set(cells_and_replicates)) != len(cells_and_replicates):
        raise ValueError("DGP cell/replicate pairs must be unique")
    if any(count == 0 for count in lane_counts.values()):
        raise ValueError("Q61, independent-R, and FACETS subsets must be nonempty")
    if any(count == 0 for count in cohort_counts.values()):
        raise ValueError("positive and negative-control cohorts must both be nonempty")
    expected_runs = _records_digest(
        "MFRM_QUALIFICATION_PLANNED_RUNS",
        "v1",
        records,
    )
    if value["planned_runs_sha256"] != expected_runs:
        raise ValueError("planned_runs_sha256 does not reconstruct")

    rng = _exact_keys(value["rng_contract"], RNG_KEYS, "rng_contract")
    for name in RNG_KEYS:
        _text(rng[name], f"rng_contract.{name}")
    for name in (
        "input_canonicalization_sha256",
        "input_validity_rules_sha256",
        "known_fixture_exclusion_registry_sha256",
    ):
        _sha256(value[name], name)
    if (
        value["negative_control_gate_registry_sha256"]
        != negative_control_gate_registry_sha256()
    ):
        raise ValueError("negative-control gate registry digest differs")
    if value["claim_scope"] not in {
        "FINITE_REGISTERED_CASES_ONLY",
        "OPERATING_CHARACTERISTIC_FAILURE_RATE",
    }:
        raise ValueError("qualification claim_scope is unsupported")
    if not _bool(
        value["fresh_confirmatory_data_required_after_qualification"],
        "fresh_confirmatory_data_required_after_qualification",
    ):
        raise ValueError("fresh confirmatory data must be required after qualification")

    source_manifest = value["source_manifest"]
    if not isinstance(source_manifest, Mapping):
        raise ValueError("source_manifest must be an object")
    if set(map(str, source_manifest.keys())) != REQUIRED_SOURCE_FILES:
        raise ValueError("source_manifest does not contain the exact required files")
    for name, digest in source_manifest.items():
        _sha256(digest, f"source_manifest[{name}]")
    expected_source = _domain_sha256(
        "MFRM_QUALIFICATION_SOURCE_MANIFEST",
        "v1",
        dict(source_manifest),
    )
    if value["source_manifest_sha256"] != expected_source:
        raise ValueError("source_manifest_sha256 does not reconstruct")

    runtime = _exact_keys(value["runtime_manifest"], RUNTIME_KEYS, "runtime_manifest")
    for name in RUNTIME_KEYS - {
        "package_lock_sha256",
        "rscript_sha256",
        "facets_executable_sha256",
    }:
        _text(runtime[name], f"runtime_manifest.{name}")
    for name in (
        "package_lock_sha256",
        "rscript_sha256",
        "facets_executable_sha256",
    ):
        _sha256(runtime[name], f"runtime_manifest.{name}")
    expected_runtime = _domain_sha256(
        "MFRM_QUALIFICATION_RUNTIME_MANIFEST",
        "v1",
        dict(runtime),
    )
    if value["runtime_manifest_sha256"] != expected_runtime:
        raise ValueError("runtime_manifest_sha256 does not reconstruct")

    numerical = _exact_keys(
        value["numerical_contracts"],
        NUMERICAL_CONTRACT_KEYS,
        "numerical_contracts",
    )
    try:
        stationarity = StationarityContract(**numerical["stationarity"])
        sensitivity = QuadratureSensitivityContract(**numerical["sensitivity"])
        optimizer = JointPolishOptions(**numerical["optimizer"])
        stationarity.validate()
        sensitivity.validate()
        optimizer.validate()
    except (TypeError, ValueError) as exc:
        raise ValueError("numerical contracts are malformed") from exc
    expected_contract_hashes = {
        "stationarity_sha256": contract_sha256(stationarity),
        "sensitivity_sha256": contract_sha256(sensitivity),
        "optimizer_sha256": contract_sha256(optimizer),
    }
    for name, expected in expected_contract_hashes.items():
        if numerical[name] != expected:
            raise ValueError(f"numerical contract digest does not reconstruct: {name}")
    _sha256(
        numerical["objective_evaluator_implementation_sha256"],
        "objective evaluator implementation",
    )
    expected_numerical = _domain_sha256(
        "MFRM_QUALIFICATION_NUMERICAL_CONTRACTS",
        "v1",
        dict(numerical),
    )
    if value["numerical_contracts_sha256"] != expected_numerical:
        raise ValueError("numerical_contracts_sha256 does not reconstruct")
    if value["failure_policy"] != FAILURE_POLICY:
        raise ValueError("failure policy differs from the no-replacement contract")
    if value["lane_separation"] != LANE_SEPARATION:
        raise ValueError("FACETS/R/scientific lane separation differs")
    _sha256(value["independent_r_plan_sha256"], "independent R plan")
    _sha256(value["facets_operational_plan_sha256"], "FACETS operational plan")
    for name in (
        "qualification_data_generated_before_registration",
        "estimator_attempts_started_before_registration",
        "scientific_endpoints_computed_before_registration",
        "scientific_inference_enabled",
    ):
        if _bool(value[name], name):
            raise ValueError(f"pre-generation gate must be false: {name}")
    return value


def pre_generation_plan_identity_sha256(plan: Mapping[str, Any]) -> str:
    validate_pre_generation_plan(plan)
    return _domain_sha256("MFRM_QUALIFICATION_T0_PLAN", PLAN_SCHEMA, plan)


def _planned_by_id(plan: Mapping[str, Any]) -> dict[str, Mapping[str, Any]]:
    validate_pre_generation_plan(plan)
    return {str(item["run_id"]): item for item in plan["planned_runs"]}


def input_instance_sha256(
    *,
    plan_identity_sha256: str,
    run_id: str,
    seed: int,
    raw_input_artifact_sha256: str,
    canonical_response_sha256: str,
) -> str:
    _sha256(plan_identity_sha256, "plan identity")
    _text(run_id, "run_id")
    _integer(seed, "seed")
    _sha256(raw_input_artifact_sha256, "raw input artifact")
    _sha256(canonical_response_sha256, "canonical response")
    return _domain_sha256(
        "MFRM_QUALIFICATION_INPUT_INSTANCE",
        "v1",
        {
            "plan_identity_sha256": plan_identity_sha256,
            "run_id": run_id,
            "seed": seed,
            "raw_input_artifact_sha256": raw_input_artifact_sha256,
            "canonical_response_sha256": canonical_response_sha256,
        },
    )


def validate_realized_inputs(
    plan: Mapping[str, Any],
    manifest: Mapping[str, Any],
) -> Mapping[str, Any]:
    planned = _planned_by_id(plan)
    value = _exact_keys(manifest, INPUT_KEYS, "realized-input manifest")
    if value["schema_version"] != INPUT_SCHEMA:
        raise ValueError("realized-input schema differs")
    if value["status"] != "SEALED_AFTER_GENERATION_BEFORE_FIT":
        raise ValueError("realized-input status differs")
    if value["plan_identity_sha256"] != pre_generation_plan_identity_sha256(plan):
        raise ValueError("realized inputs cite a different T0 plan")
    if value["planned_runs_sha256"] != plan["planned_runs_sha256"]:
        raise ValueError("realized-input denominator differs from T0")
    records = value["records"]
    if not isinstance(records, list):
        raise ValueError("realized-input records must be a list")
    if [item.get("run_id") for item in records if isinstance(item, Mapping)] != list(planned):
        raise ValueError("realized-input records must exactly follow planned RunIds")
    input_instances: list[str] = []
    for index, record in enumerate(records):
        item = _exact_keys(record, INPUT_RECORD_KEYS, f"input records[{index}]")
        run_id = _text(item["run_id"], f"input records[{index}].run_id")
        seed = _integer(item["seed"], f"{run_id} seed")
        if seed != planned[run_id]["seed"]:
            raise ValueError(f"realized-input seed differs for {run_id}")
        status = item["status"]
        if status in {"GENERATED", "INPUT_INVALID_UNDER_PREREGISTERED_RULE"}:
            raw_digest = _sha256(item["raw_input_artifact_sha256"], f"{run_id} raw input")
            response_digest = _sha256(
                item["canonical_response_sha256"],
                f"{run_id} canonical response",
            )
            expected_instance = input_instance_sha256(
                plan_identity_sha256=value["plan_identity_sha256"],
                run_id=run_id,
                seed=seed,
                raw_input_artifact_sha256=str(raw_digest),
                canonical_response_sha256=str(response_digest),
            )
            if item["input_instance_sha256"] != expected_instance:
                raise ValueError(f"input instance digest does not reconstruct for {run_id}")
            input_instances.append(str(expected_instance))
            _integer(item["rows"], f"{run_id} rows", minimum=1)
            if status == "GENERATED":
                if (
                    item["failure_code"] is not None
                    or item["generation_error_identity_sha256"] is not None
                ):
                    raise ValueError(f"generated input has failure evidence for {run_id}")
            else:
                if item["failure_code"] not in INPUT_INVALID_FAILURE_CODES:
                    raise ValueError(f"unsupported input-validity failure code for {run_id}")
                _sha256(
                    item["generation_error_identity_sha256"],
                    f"{run_id} input-validity evidence",
                )
        elif status == "GENERATION_FAILED":
            for name in (
                "raw_input_artifact_sha256",
                "canonical_response_sha256",
                "input_instance_sha256",
                "rows",
            ):
                if item[name] is not None:
                    raise ValueError(f"generation failure retained {name} for {run_id}")
            if item["failure_code"] not in GENERATION_FAILURE_CODES:
                raise ValueError(f"unsupported generation failure code for {run_id}")
            _sha256(
                item["generation_error_identity_sha256"],
                f"{run_id} generation error",
            )
        else:
            raise ValueError(f"unsupported generation status for {run_id}")
    if len(set(input_instances)) != len(input_instances):
        raise ValueError("realized input instance digests must be unique")
    expected_records = _records_digest(
        "MFRM_QUALIFICATION_REALIZED_INPUT_RECORDS",
        "v1",
        records,
    )
    if value["records_sha256"] != expected_records:
        raise ValueError("realized-input records_sha256 does not reconstruct")
    for name in (
        "estimator_attempts_started_before_manifest",
        "scientific_endpoints_computed_before_manifest",
        "scientific_inference_enabled",
    ):
        if _bool(value[name], name):
            raise ValueError(f"realized-input pre-fit gate must be false: {name}")
    return value


def realized_inputs_identity_sha256(
    plan: Mapping[str, Any],
    manifest: Mapping[str, Any],
) -> str:
    validate_realized_inputs(plan, manifest)
    return _domain_sha256("MFRM_QUALIFICATION_T1_INPUTS", INPUT_SCHEMA, manifest)


def validate_fit_authorization(
    plan: Mapping[str, Any],
    inputs: Mapping[str, Any],
    authorization: Mapping[str, Any],
) -> Mapping[str, Any]:
    planned = _planned_by_id(plan)
    validate_realized_inputs(plan, inputs)
    input_by_id = {str(item["run_id"]): item for item in inputs["records"]}
    value = _exact_keys(authorization, FIT_KEYS, "fit authorization")
    if value["schema_version"] != FIT_SCHEMA:
        raise ValueError("fit-authorization schema differs")
    if value["status"] != "DEVELOPMENT_ONLY_PREFIT_FACTORY_NOT_QUALIFIED":
        raise ValueError("fit authorization must remain development-only")
    if value["plan_identity_sha256"] != pre_generation_plan_identity_sha256(plan):
        raise ValueError("fit authorization cites a different T0 plan")
    if value["realized_inputs_identity_sha256"] != realized_inputs_identity_sha256(
        plan, inputs
    ):
        raise ValueError("fit authorization cites different T1 inputs")
    if value["planned_runs_sha256"] != plan["planned_runs_sha256"]:
        raise ValueError("fit-authorization denominator differs from T0")
    _sha256(
        value["problem_digest_factory_implementation_sha256"],
        "problem-digest factory implementation",
    )
    records = value["records"]
    if not isinstance(records, list) or [
        item.get("run_id") for item in records if isinstance(item, Mapping)
    ] != list(planned):
        raise ValueError("fit-authorization records must exactly follow planned RunIds")
    for index, record in enumerate(records):
        item = _exact_keys(record, FIT_RECORD_KEYS, f"fit records[{index}]")
        run_id = _text(item["run_id"], f"fit records[{index}].run_id")
        seed = _integer(item["seed"], f"{run_id} seed")
        if seed != planned[run_id]["seed"]:
            raise ValueError(f"fit-authorization seed differs for {run_id}")
        realized = input_by_id[run_id]
        if realized["status"] == "GENERATED":
            if item["status"] != "AUTHORIZED_FOR_SINGLE_ATTEMPT_DEVELOPMENT_ONLY":
                raise ValueError(f"generated run is not singly authorized: {run_id}")
            if item["input_instance_sha256"] != realized["input_instance_sha256"]:
                raise ValueError(f"fit authorization input instance differs for {run_id}")
            problem_digest = _sha256(
                item["likelihood_problem_digest"],
                f"{run_id} likelihood problem",
            )
            _sha256(
                item["estimator_input_config_sha256"],
                f"{run_id} estimator input config",
            )
        else:
            expected_status = (
                "NOT_AUTHORIZED_GENERATION_FAILED"
                if realized["status"] == "GENERATION_FAILED"
                else "NOT_AUTHORIZED_INPUT_INVALID"
            )
            if item["status"] != expected_status:
                raise ValueError(f"invalid input was incorrectly authorized: {run_id}")
            for name in (
                "input_instance_sha256",
                "likelihood_problem_digest",
                "estimator_input_config_sha256",
            ):
                if item[name] is not None:
                    raise ValueError(f"unauthorized run retained {name}: {run_id}")
    expected_records = _records_digest(
        "MFRM_QUALIFICATION_FIT_AUTHORIZATION_RECORDS",
        "v1",
        records,
    )
    if value["records_sha256"] != expected_records:
        raise ValueError("fit-authorization records_sha256 does not reconstruct")
    if _bool(
        value["estimator_attempts_started_before_authorization"],
        "estimator_attempts_started_before_authorization",
    ):
        raise ValueError("fit attempts started before authorization")
    if _bool(
        value["pure_prefit_problem_factory_qualified"],
        "pure_prefit_problem_factory_qualified",
    ):
        raise ValueError("pure pre-fit problem factory is not yet qualified")
    if _bool(value["scientific_inference_enabled"], "scientific_inference_enabled"):
        raise ValueError("fit authorization cannot enable scientific inference")
    return value


def fit_authorization_identity_sha256(
    plan: Mapping[str, Any],
    inputs: Mapping[str, Any],
    authorization: Mapping[str, Any],
) -> str:
    validate_fit_authorization(plan, inputs, authorization)
    return _domain_sha256(
        "MFRM_QUALIFICATION_T1_5_FIT_AUTHORIZATION",
        FIT_SCHEMA,
        authorization,
    )


def validate_attempt_ledger(
    plan: Mapping[str, Any],
    inputs: Mapping[str, Any],
    authorization: Mapping[str, Any],
    ledger: Mapping[str, Any],
) -> Mapping[str, Any]:
    planned = _planned_by_id(plan)
    validate_fit_authorization(plan, inputs, authorization)
    input_by_id = {str(item["run_id"]): item for item in inputs["records"]}
    fit_by_id = {str(item["run_id"]): item for item in authorization["records"]}
    value = _exact_keys(ledger, LEDGER_KEYS, "attempt ledger")
    if value["schema_version"] != LEDGER_SCHEMA:
        raise ValueError("attempt-ledger schema differs")
    if value["status"] != "COMPLETE_DEVELOPMENT_LEDGER_NOT_A_QUALIFICATION":
        raise ValueError("attempt-ledger status differs")
    expected_parents = {
        "plan_identity_sha256": pre_generation_plan_identity_sha256(plan),
        "realized_inputs_identity_sha256": realized_inputs_identity_sha256(plan, inputs),
        "fit_authorization_identity_sha256": fit_authorization_identity_sha256(
            plan, inputs, authorization
        ),
        "planned_runs_sha256": plan["planned_runs_sha256"],
    }
    for name, expected in expected_parents.items():
        if value[name] != expected:
            raise ValueError(f"attempt ledger parent differs: {name}")
    records = value["records"]
    if not isinstance(records, list) or [
        item.get("run_id") for item in records if isinstance(item, Mapping)
    ] != list(planned):
        raise ValueError("attempt-ledger records must exactly follow planned RunIds")
    for index, record in enumerate(records):
        item = _exact_keys(record, LEDGER_RECORD_KEYS, f"ledger records[{index}]")
        run_id = _text(item["run_id"], f"ledger records[{index}].run_id")
        seed = _integer(item["seed"], f"{run_id} seed")
        attempt_count = _integer(item["attempt_count"], f"{run_id} attempt_count")
        if seed != planned[run_id]["seed"]:
            raise ValueError(f"attempt-ledger seed differs for {run_id}")
        if item["expected_outcome"] != planned[run_id]["expected_outcome"]:
            raise ValueError(f"attempt-ledger expected outcome differs for {run_id}")
        realized = input_by_id[run_id]
        fit = fit_by_id[run_id]
        status = item["status"]
        assessment_returned = status == "ASSESSMENT_RETURNED"
        lane_requirements = {
            "python_q31": True,
            "python_q61": _bool(planned[run_id]["q61_required"], "q61_required"),
            "independent_r": _bool(
                planned[run_id]["independent_r_required"],
                "independent_r_required",
            ),
            "facets": _bool(
                planned[run_id]["facets_comparison_required"],
                "facets_comparison_required",
            ),
        }
        lane_statuses = {
            lane: _validate_lane_evidence(
                item,
                run_id=run_id,
                lane=lane,
                required=required,
                attempted=assessment_returned,
            )
            for lane, required in lane_requirements.items()
        }
        q31_failed_gates = item["python_q31_failed_gates"]
        if assessment_returned:
            if not isinstance(q31_failed_gates, list):
                raise ValueError(f"Q31 failed-gate list is malformed for {run_id}")
            if (
                q31_failed_gates != sorted(q31_failed_gates)
                or len(set(q31_failed_gates)) != len(q31_failed_gates)
                or any(gate not in PYTHON_Q31_GATES for gate in q31_failed_gates)
            ):
                raise ValueError(f"Q31 failed-gate list is not canonical for {run_id}")
            if (lane_statuses["python_q31"] == "PASS") != (
                len(q31_failed_gates) == 0
            ):
                raise ValueError(f"Q31 status and failed gates differ for {run_id}")
        elif q31_failed_gates is not None:
            raise ValueError(f"unattempted Q31 retained failed gates for {run_id}")
        if item["mml_scientific_status"] != "NOT_SCIENTIFIC_INFERENCE":
            raise ValueError(f"MML scientific status differs for {run_id}")
        expected_facets_status = (
            lane_statuses["facets"]
            if lane_requirements["facets"]
            else "NOT_REQUIRED"
        )
        if item["facets_operational_status"] != expected_facets_status:
            raise ValueError(f"FACETS operational status does not reconstruct for {run_id}")
        if assessment_returned:
            mml_required_lanes = ["python_q31"]
            if lane_requirements["python_q61"]:
                mml_required_lanes.append("python_q61")
            if lane_requirements["independent_r"]:
                mml_required_lanes.append("independent_r")
            expected_mml_status = (
                "PASS"
                if all(lane_statuses[lane] == "PASS" for lane in mml_required_lanes)
                else "FAIL"
            )
        else:
            expected_mml_status = "NOT_OBSERVED"
        if item["mml_numerical_status"] != expected_mml_status:
            raise ValueError(f"MML numerical status does not reconstruct for {run_id}")
        if realized["status"] != "GENERATED":
            expected_status = (
                "NOT_ATTEMPTED_GENERATION_FAILED"
                if realized["status"] == "GENERATION_FAILED"
                else "NOT_ATTEMPTED_INPUT_INVALID"
            )
            if status != expected_status or attempt_count != 0:
                raise ValueError(f"non-generated attempt state differs for {run_id}")
            if item["error_identity_sha256"] != realized["generation_error_identity_sha256"]:
                raise ValueError(f"non-generated failure identity differs for {run_id}")
            expected_instance = (
                None
                if realized["status"] == "GENERATION_FAILED"
                else realized["input_instance_sha256"]
            )
            if item["input_instance_sha256"] != expected_instance:
                raise ValueError(f"non-generated input instance differs for {run_id}")
            if item["attempt_id"] != f"{run_id}::no-fit":
                raise ValueError(f"non-generated attempt id differs for {run_id}")
            if item["observed_outcome"] != "NOT_OBSERVED_NO_FIT":
                raise ValueError(f"non-generated observed outcome differs for {run_id}")
            if item["observed_gate"] is not None:
                raise ValueError(f"non-generated run retained an observed gate: {run_id}")
            if item["failure_code"] != realized["failure_code"]:
                raise ValueError(f"non-generated failure code differs for {run_id}")
            for name in (
                "qualification_record_sha256",
                "likelihood_problem_digest",
                "warm_start_artifact_sha256",
                "stdout_sha256",
                "stderr_sha256",
                "started_utc",
                "finished_utc",
            ):
                if item[name] is not None:
                    raise ValueError(f"non-generated run retained {name}: {run_id}")
            expected_stage = (
                "GENERATION"
                if realized["status"] == "GENERATION_FAILED"
                else "INPUT_VALIDITY"
            )
            if item["failure_stage"] != expected_stage:
                raise ValueError(f"non-generated failure stage differs for {run_id}")
            continue
        if item["input_instance_sha256"] != realized["input_instance_sha256"]:
            raise ValueError(f"attempt input instance differs for {run_id}")
        if item["likelihood_problem_digest"] != fit["likelihood_problem_digest"]:
            raise ValueError(f"attempt problem digest differs for {run_id}")
        if attempt_count != 1:
            raise ValueError(f"generated run must have exactly one fit attempt: {run_id}")
        if item["attempt_id"] != f"{run_id}::attempt-01":
            raise ValueError(f"fit attempt id differs for {run_id}")
        _sha256(item["stdout_sha256"], f"{run_id} stdout")
        _sha256(item["stderr_sha256"], f"{run_id} stderr")
        started = _utc_datetime(item["started_utc"], f"{run_id} started_utc")
        finished = _utc_datetime(item["finished_utc"], f"{run_id} finished_utc")
        if finished < started:
            raise ValueError(f"fit attempt timestamps are reversed for {run_id}")
        if status == "ASSESSMENT_RETURNED":
            _sha256(item["qualification_record_sha256"], f"{run_id} qualification record")
            _sha256(item["warm_start_artifact_sha256"], f"{run_id} warm start")
            expected_outcome = planned[run_id]["expected_outcome"]
            allowed_observed = (
                {
                    "EXPECTED_GATE_FAILED",
                    "EXPECTED_AND_OTHER_GATES_FAILED",
                    "UNEXPECTED_PASS",
                    "WRONG_GATE_FAILED",
                }
                if expected_outcome == "NEGATIVE_CONTROL"
                else {"MML_REQUIRED_GATES_PASS", "NOT_QUALIFIED"}
            )
            if item["observed_outcome"] not in allowed_observed:
                raise ValueError(f"assessment observed outcome differs for {run_id}")
            if expected_outcome == "NEGATIVE_CONTROL":
                failed_gates = _failed_registered_gates(item)
                expected_gate = str(planned[run_id]["expected_gate"])
                if failed_gates == (expected_gate,):
                    derived_outcome = "EXPECTED_GATE_FAILED"
                    derived_gate: str | None = expected_gate
                elif expected_gate in failed_gates:
                    derived_outcome = "EXPECTED_AND_OTHER_GATES_FAILED"
                    derived_gate = next(
                        gate for gate in failed_gates if gate != expected_gate
                    )
                elif failed_gates:
                    derived_outcome = "WRONG_GATE_FAILED"
                    derived_gate = failed_gates[0]
                else:
                    derived_outcome = "UNEXPECTED_PASS"
                    derived_gate = None
                if (
                    item["observed_outcome"] != derived_outcome
                    or item["observed_gate"] != derived_gate
                ):
                    raise ValueError(
                        f"negative-control outcome does not reconstruct for {run_id}"
                    )
            elif expected_mml_status == "PASS":
                if item["observed_outcome"] != "MML_REQUIRED_GATES_PASS":
                    raise ValueError(
                        f"positive MML outcome contradicts passing lane evidence: {run_id}"
                    )
                if item["observed_gate"] is not None:
                    raise ValueError(f"passing positive case retained a failed gate: {run_id}")
            else:
                if item["observed_outcome"] != "NOT_QUALIFIED":
                    raise ValueError(
                        f"positive MML outcome contradicts failing lane evidence: {run_id}"
                    )
                _text(item["observed_gate"], f"{run_id} failed positive gate")
                if (
                    item["observed_gate"] not in NEGATIVE_CONTROL_GATES
                    or item["observed_gate"] == "FACETS_OPERATIONAL_PASS"
                    or not _gate_failed(item, str(item["observed_gate"]))
                ):
                    raise ValueError(f"positive failed gate is not in lane evidence: {run_id}")
            if (
                item["failure_stage"] is not None
                or item["failure_code"] is not None
                or item["error_identity_sha256"] is not None
            ):
                raise ValueError(f"returned assessment retained failure fields: {run_id}")
        elif status == "EXECUTION_FAILED":
            if item["qualification_record_sha256"] is not None:
                raise ValueError(f"failed attempt retained a qualification record: {run_id}")
            if item["failure_stage"] not in ATTEMPT_FAILURE_CODES:
                raise ValueError(f"failed attempt stage is unsupported: {run_id}")
            if item["failure_stage"] != "PRECHECK":
                _sha256(item["warm_start_artifact_sha256"], f"{run_id} warm start")
            elif item["warm_start_artifact_sha256"] is not None:
                raise ValueError(f"PRECHECK failure retained a warm start for {run_id}")
            if item["failure_code"] not in ATTEMPT_FAILURE_CODES[item["failure_stage"]]:
                raise ValueError(f"attempt failure code/stage differ for {run_id}")
            _sha256(item["error_identity_sha256"], f"{run_id} attempt error")
            if item["observed_outcome"] != "NOT_OBSERVED_EXECUTION_FAILED":
                raise ValueError(f"failed attempt observed outcome differs for {run_id}")
            if item["observed_gate"] is not None:
                raise ValueError(f"failed execution retained an observed gate: {run_id}")
        else:
            raise ValueError(f"attempt status is unsupported for {run_id}")
    expected_records = _records_digest(
        "MFRM_QUALIFICATION_ATTEMPT_LEDGER_RECORDS",
        "v1",
        records,
    )
    if value["records_sha256"] != expected_records:
        raise ValueError("attempt-ledger records_sha256 does not reconstruct")
    for name in (
        "replacements_performed",
        "optional_extensions_performed",
        "scientific_inference_enabled",
    ):
        if _bool(value[name], name):
            raise ValueError(f"attempt-ledger gate must be false: {name}")
    return value


def attempt_ledger_identity_sha256(
    plan: Mapping[str, Any],
    inputs: Mapping[str, Any],
    authorization: Mapping[str, Any],
    ledger: Mapping[str, Any],
) -> str:
    validate_attempt_ledger(plan, inputs, authorization, ledger)
    return _domain_sha256(
        "MFRM_QUALIFICATION_T2_ATTEMPT_LEDGER",
        LEDGER_SCHEMA,
        ledger,
    )


def _durable_create(path: Path, payload: bytes) -> None:
    with path.open("xb") as handle:
        handle.write(payload)
        handle.flush()
        os.fsync(handle.fileno())


def _is_link_or_junction(path: Path) -> bool:
    return path.is_symlink() or bool(
        hasattr(path, "is_junction") and path.is_junction()
    )


def _expected_parent_identities(artifact: Mapping[str, Any]) -> dict[str, str]:
    schema = artifact.get("schema_version")
    if schema == PLAN_SCHEMA:
        return {}
    if schema == INPUT_SCHEMA:
        return {"plan": str(artifact["plan_identity_sha256"])}
    if schema == FIT_SCHEMA:
        return {
            "plan": str(artifact["plan_identity_sha256"]),
            "realized_inputs": str(artifact["realized_inputs_identity_sha256"]),
        }
    if schema == LEDGER_SCHEMA:
        return {
            "plan": str(artifact["plan_identity_sha256"]),
            "realized_inputs": str(artifact["realized_inputs_identity_sha256"]),
            "fit_authorization": str(artifact["fit_authorization_identity_sha256"]),
        }
    raise ValueError("published registrar artifact schema is unsupported")


def validate_published_artifact(
    output_dir: Path,
    *,
    validator: Callable[[Mapping[str, Any]], Mapping[str, Any]],
    semantic_identity_function: Callable[[Mapping[str, Any]], str],
) -> Mapping[str, Any]:
    supplied_output = Path(output_dir)
    if _is_link_or_junction(supplied_output):
        raise ValueError("published artifact directory is a link or junction")
    output = supplied_output.resolve()
    if not output.is_dir() or _is_link_or_junction(output):
        raise ValueError("published artifact directory is invalid")
    if {path.name for path in output.iterdir()} != {"artifact.json", "identity.json"}:
        raise ValueError("published artifact file set differs")
    for name in ("artifact.json", "identity.json"):
        path = output / name
        if not path.is_file() or path.is_symlink():
            raise ValueError(f"published artifact member is invalid: {name}")
    artifact_bytes = (output / "artifact.json").read_bytes()
    identity_bytes = (output / "identity.json").read_bytes()
    artifact = _load_strict_json_object(artifact_bytes, "published artifact")
    validator(artifact)
    identity = _load_strict_json_object(
        identity_bytes,
        "published identity",
    )
    expected_identity_keys = {
        "schema_version",
        "artifact_schema_version",
        "artifact_sha256",
        "artifact_bytes",
        "semantic_identity_sha256",
        "parent_identity_sha256",
        "publisher_sha256",
        "scientific_inference_ready",
    }
    _exact_keys(identity, expected_identity_keys, "artifact identity")
    if identity["schema_version"] != IDENTITY_SCHEMA:
        raise ValueError("artifact identity schema differs")
    if identity["artifact_schema_version"] != artifact.get("schema_version"):
        raise ValueError("artifact identity cites a different schema")
    if identity["artifact_sha256"] != hashlib.sha256(artifact_bytes).hexdigest():
        raise ValueError("artifact byte hash differs")
    if _integer(identity["artifact_bytes"], "artifact_bytes", minimum=1) != len(
        artifact_bytes
    ):
        raise ValueError("artifact byte count differs")
    expected_semantic_identity = semantic_identity_function(artifact)
    _sha256(expected_semantic_identity, "reconstructed semantic identity")
    if identity["semantic_identity_sha256"] != expected_semantic_identity:
        raise ValueError("artifact semantic identity does not reconstruct")
    parents = identity["parent_identity_sha256"]
    if not isinstance(parents, Mapping):
        raise ValueError("parent identities must be an object")
    for name, digest in parents.items():
        _text(name, "parent identity name")
        _sha256(digest, f"parent identity {name}")
    if dict(parents) != _expected_parent_identities(artifact):
        raise ValueError("published parent identities do not reconstruct")
    _sha256(identity["publisher_sha256"], "publisher source")
    publisher_path = Path(__file__).resolve()
    publisher_digest = hashlib.sha256(publisher_path.read_bytes()).hexdigest()
    if identity["publisher_sha256"] != publisher_digest:
        raise ValueError("published publisher source does not reconstruct")
    if _bool(identity["scientific_inference_ready"], "scientific_inference_ready"):
        raise ValueError("published registrar artifact cannot enable scientific inference")
    final_members = list(output.iterdir())
    if (
        {path.name for path in final_members} != {"artifact.json", "identity.json"}
        or any(not path.is_file() or _is_link_or_junction(path) for path in final_members)
        or (output / "artifact.json").read_bytes() != artifact_bytes
        or (output / "identity.json").read_bytes() != identity_bytes
        or hashlib.sha256(publisher_path.read_bytes()).hexdigest() != publisher_digest
    ):
        raise ValueError("published registrar artifact changed during validation")
    return {"artifact": artifact, "identity": identity}


def publish_sealed_artifact(
    *,
    output_dir: Path,
    document: Mapping[str, Any],
    validator: Callable[[Mapping[str, Any]], Mapping[str, Any]],
    semantic_identity_function: Callable[[Mapping[str, Any]], str],
    parent_identity_sha256: Mapping[str, str] | None = None,
) -> Mapping[str, Any]:
    """Publish one new exact two-file directory with identity last."""

    output = Path(output_dir).resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite registrar artifact: {output}")
    validator(document)
    semantic_identity_sha256 = semantic_identity_function(document)
    _sha256(semantic_identity_sha256, "semantic identity")
    parents = _expected_parent_identities(document)
    if parent_identity_sha256 is not None and dict(parent_identity_sha256) != parents:
        raise ValueError("supplied parent identities differ from artifact parents")
    for name, digest in parents.items():
        _text(name, "parent identity name")
        _sha256(digest, f"parent identity {name}")
    artifact_bytes = _json_bytes(document)
    decoded = _load_strict_json_object(artifact_bytes, "staged artifact")
    validator(decoded)
    publisher_path = Path(__file__).resolve()
    publisher_before = hashlib.sha256(publisher_path.read_bytes()).hexdigest()
    identity = {
        "schema_version": IDENTITY_SCHEMA,
        "artifact_schema_version": document["schema_version"],
        "artifact_sha256": hashlib.sha256(artifact_bytes).hexdigest(),
        "artifact_bytes": len(artifact_bytes),
        "semantic_identity_sha256": semantic_identity_sha256,
        "parent_identity_sha256": parents,
        "publisher_sha256": publisher_before,
        "scientific_inference_ready": False,
    }
    identity_bytes = _json_bytes(identity)
    output.mkdir(parents=True, exist_ok=False)
    nonce = f"{os.getpid()}.{time.time_ns()}"
    artifact_stage = output / f".artifact.json.stage.{nonce}"
    identity_stage = output / f".identity.json.stage.{nonce}"
    try:
        _durable_create(artifact_stage, artifact_bytes)
        _durable_create(identity_stage, identity_bytes)
        if hashlib.sha256(publisher_path.read_bytes()).hexdigest() != publisher_before:
            raise ValueError("registrar publisher changed during publication")
        os.replace(artifact_stage, output / "artifact.json")
        os.replace(identity_stage, output / "identity.json")
        published = validate_published_artifact(
            output,
            validator=validator,
            semantic_identity_function=semantic_identity_function,
        )
        if published["artifact"] != decoded or published["identity"] != identity:
            raise ValueError("published registrar artifact differs from staged values")
        return published
    except BaseException:
        for path in (
            output / "identity.json",
            output / "artifact.json",
            identity_stage,
            artifact_stage,
        ):
            try:
                path.unlink(missing_ok=True)
            except OSError:
                pass
        try:
            output.rmdir()
        except OSError:
            pass
        raise


__all__ = [
    "FAILURE_POLICY",
    "FIT_SCHEMA",
    "IDENTITY_SCHEMA",
    "INPUT_SCHEMA",
    "LANE_SEPARATION",
    "LEDGER_SCHEMA",
    "PLAN_SCHEMA",
    "REQUIRED_SOURCE_FILES",
    "attempt_ledger_identity_sha256",
    "fit_authorization_identity_sha256",
    "input_instance_sha256",
    "negative_control_gate_registry_sha256",
    "pre_generation_plan_identity_sha256",
    "publish_sealed_artifact",
    "realized_inputs_identity_sha256",
    "validate_attempt_ledger",
    "validate_fit_authorization",
    "validate_pre_generation_plan",
    "validate_published_artifact",
    "validate_realized_inputs",
]
