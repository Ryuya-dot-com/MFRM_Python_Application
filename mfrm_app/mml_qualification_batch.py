"""Exact-denominator batching for Q-fit/Q-check numerical evidence.

This development-only layer binds every planned RunId to one unique
likelihood-problem digest,
requires one frozen stationarity contract and one frozen quadrature-sensitivity
contract for the whole batch, and reconstructs all pass/fail summaries from raw
runs and independently replayed cross-evaluations.  The replay resolver is an
auditor-side trust boundary: a future registrar must construct it from frozen
raw inputs and an exact source/runtime manifest.  Records cannot supply it.
This module cannot authenticate a registration or create scientific readiness.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, is_dataclass
import hashlib
import json
import re
from typing import Callable, Iterable, Mapping

import numpy as np

from mfrm_app.mml_engine_v2 import (
    StationarityAssessment,
    StationarityContract,
    assess_free_sd_stationarity,
)
from mfrm_app.mml_quadrature_sensitivity import (
    QuadratureSensitivityAssessment,
    QuadratureSensitivityContract,
)
from mfrm_app.mml_stationarity import JointPolishOptions, ValueFunction


SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")


def _json_ready(value: object) -> object:
    if is_dataclass(value) and not isinstance(value, type):
        return _json_ready(asdict(value))
    if isinstance(value, np.ndarray):
        return _json_ready(value.tolist())
    if isinstance(value, np.generic):
        return _json_ready(value.item())
    if isinstance(value, Mapping):
        return {str(key): _json_ready(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_ready(item) for item in value]
    return value


def _canonical_sha256(value: object) -> str:
    encoded = json.dumps(
        _json_ready(value),
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _domain_sha256(domain: str, schema_version: str, payload: object) -> str:
    return _canonical_sha256(
        {
            "domain": domain,
            "schema_version": schema_version,
            "payload": payload,
        }
    )


def contract_sha256(contract: object, *, domain: str | None = None) -> str:
    """Return the canonical digest of a typed or mapping numerical contract."""

    if is_dataclass(contract) and not isinstance(contract, type):
        payload = asdict(contract)
    elif isinstance(contract, Mapping):
        payload = dict(contract)
    else:
        raise TypeError("contract must be a dataclass instance or mapping")
    if domain is None:
        if isinstance(contract, StationarityContract):
            domain = "MFRM_FREE_SD_STATIONARITY_CONTRACT"
        elif isinstance(contract, QuadratureSensitivityContract):
            domain = "MFRM_FREE_SD_QUADRATURE_SENSITIVITY_CONTRACT"
        elif isinstance(contract, JointPolishOptions):
            domain = "MFRM_FREE_SD_JOINT_POLISH_OPTIONS"
        else:
            raise ValueError("mapping contract hashes require an explicit domain")
    return _domain_sha256(domain, "v1", payload)


def batch_contract_sha256(contract: "NumericalQualificationBatchContract") -> str:
    if not isinstance(contract, NumericalQualificationBatchContract):
        raise TypeError("contract must be a NumericalQualificationBatchContract")
    return _domain_sha256(
        "MFRM_FREE_SD_NUMERICAL_QUALIFICATION_BATCH_CONTRACT",
        "v2",
        asdict(contract),
    )


def planned_run_ids_sha256(run_ids: Iterable[str]) -> str:
    values = [str(item) for item in run_ids]
    return _domain_sha256("MFRM_PLANNED_RUN_IDS", "v1", values)


def planned_problem_digests_sha256(
    values: Iterable[tuple[str, str]],
) -> str:
    return _domain_sha256(
        "MFRM_PLANNED_RUN_ID_TO_PROBLEM_DIGEST_MAP",
        "v1",
        [(str(run_id), str(digest)) for run_id, digest in values],
    )


def _is_sha256(value: object) -> bool:
    return bool(SHA256_PATTERN.fullmatch(str(value)))


def _strict_bool(value: object, label: str) -> bool:
    if not isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{label} must be boolean")
    return bool(value)


def _consistent(left: object, right: object) -> bool:
    try:
        left_value = float(left)
        right_value = float(right)
    except (TypeError, ValueError):
        return False
    return bool(
        np.isfinite(left_value)
        and np.isfinite(right_value)
        and abs(left_value - right_value)
        <= 1e-12 * max(1.0, abs(left_value), abs(right_value))
    )


@dataclass(frozen=True)
class NumericalQualificationBatchContract:
    batch_id: str
    scope: str
    registration_identity_sha256: str | None
    planned_run_ids: tuple[str, ...]
    planned_run_ids_sha256: str
    planned_problem_digests: tuple[tuple[str, str], ...]
    planned_problem_digests_sha256: str
    stationarity_contract_sha256: str
    sensitivity_contract_sha256: str
    optimizer_options_sha256: str
    objective_evaluator_implementation_sha256: str
    primary_quadrature_points: int
    sensitivity_quadrature_points: int
    all_planned_required: bool
    replacement_forbidden: bool

    def validate(self) -> None:
        if not str(self.batch_id).strip():
            raise ValueError("batch_id must be nonempty")
        if self.scope != "DEVELOPMENT_ONLY":
            raise ValueError("only DEVELOPMENT_ONLY scope is implemented")
        if self.registration_identity_sha256 is not None:
            raise ValueError("development batches cannot cite a registration identity")
        ids = tuple(str(item) for item in self.planned_run_ids)
        if not ids or any(not item.strip() for item in ids) or len(set(ids)) != len(ids):
            raise ValueError("planned RunIds must be nonempty and unique")
        if ids != tuple(sorted(ids)):
            raise ValueError("planned RunIds must be stored in canonical sorted order")
        expected = planned_run_ids_sha256(ids)
        if self.planned_run_ids_sha256 != expected:
            raise ValueError("planned RunId digest does not reconstruct")

        problem_map = tuple(
            (str(run_id), str(digest))
            for run_id, digest in self.planned_problem_digests
        )
        if tuple(run_id for run_id, _digest in problem_map) != ids:
            raise ValueError("planned problem map must exactly follow planned RunIds")
        problem_digests = tuple(digest for _run_id, digest in problem_map)
        if any(not _is_sha256(digest) for digest in problem_digests):
            raise ValueError("planned problem digests must be lowercase SHA256 values")
        if len(set(problem_digests)) != len(problem_digests):
            raise ValueError("planned problem digests must be unique")
        if (
            self.planned_problem_digests_sha256
            != planned_problem_digests_sha256(problem_map)
        ):
            raise ValueError("planned problem-map digest does not reconstruct")
        for label, digest in (
            ("stationarity contract", self.stationarity_contract_sha256),
            ("sensitivity contract", self.sensitivity_contract_sha256),
            ("optimizer options", self.optimizer_options_sha256),
            (
                "objective evaluator implementation",
                self.objective_evaluator_implementation_sha256,
            ),
        ):
            if not _is_sha256(digest):
                raise ValueError(f"{label} digest must be a lowercase SHA256 value")
        if not isinstance(self.all_planned_required, bool) or not self.all_planned_required:
            raise ValueError("all planned datasets must be required")
        if not isinstance(self.replacement_forbidden, bool) or not self.replacement_forbidden:
            raise ValueError("replacement must be forbidden")
        if (
            isinstance(self.primary_quadrature_points, bool)
            or isinstance(self.sensitivity_quadrature_points, bool)
            or not isinstance(self.primary_quadrature_points, int)
            or not isinstance(self.sensitivity_quadrature_points, int)
            or self.primary_quadrature_points < 3
            or self.sensitivity_quadrature_points <= self.primary_quadrature_points
        ):
            raise ValueError("batch quadrature points are invalid")


@dataclass(frozen=True)
class NumericalQualificationRecord:
    batch_id: str
    scope: str
    registration_identity_sha256: str | None
    run_id: str
    problem_digest: str
    batch_contract_sha256: str
    assessment_sha256: str
    record_sha256: str
    assessment: QuadratureSensitivityAssessment


@dataclass(frozen=True)
class QuadratureObjectiveEvaluator:
    """Auditor-side replay handle; it must never be supplied by a record."""

    problem_digest: str
    primary_quadrature_points: int
    sensitivity_quadrature_points: int
    primary_value_function: ValueFunction
    sensitivity_value_function: ValueFunction

    def validate(self) -> None:
        if not _is_sha256(self.problem_digest):
            raise ValueError("objective evaluator problem digest is malformed")
        if (
            isinstance(self.primary_quadrature_points, bool)
            or isinstance(self.sensitivity_quadrature_points, bool)
            or not isinstance(self.primary_quadrature_points, int)
            or not isinstance(self.sensitivity_quadrature_points, int)
            or self.primary_quadrature_points < 3
            or self.sensitivity_quadrature_points
            <= self.primary_quadrature_points
        ):
            raise ValueError("objective evaluator quadrature bases are invalid")
        if not callable(self.primary_value_function) or not callable(
            self.sensitivity_value_function
        ):
            raise ValueError("objective evaluator functions must be callable")


ObjectiveEvaluatorResolver = Callable[
    [str, str, str],
    QuadratureObjectiveEvaluator,
]


@dataclass(frozen=True)
class NumericalQualificationBatchAssessment:
    batch_id: str
    scope: str
    registration_identity_sha256: str | None
    batch_contract_sha256: str
    planned_run_ids_sha256: str
    planned_problem_digests_sha256: str
    stationarity_contract_sha256: str
    sensitivity_contract_sha256: str
    optimizer_options_sha256: str
    objective_evaluator_implementation_sha256: str
    planned_datasets: int
    returned_datasets: int
    finite_cross_evaluation_datasets: int
    primary_stationarity_pass_datasets: int
    sensitivity_stationarity_pass_datasets: int
    numerical_sensitivity_pass_datasets: int
    missing_run_ids: tuple[str, ...]
    unexpected_run_ids: tuple[str, ...]
    failed_run_ids: tuple[str, ...]
    exact_denominator_pass: bool
    all_numerical_sensitivity_pass: bool
    numerical_batch_pass: bool
    scientific_inference_ready: bool
    status: str
    record_sha256: dict[str, str]
    contract: dict[str, object]

    def to_dict(self) -> dict[str, object]:
        result = asdict(self)
        result.update(
            {
                "NumericalBatchPass": self.numerical_batch_pass,
                "ScientificInferenceReady": self.scientific_inference_ready,
            }
        )
        return result


def make_numerical_qualification_record(
    batch_contract: NumericalQualificationBatchContract,
    run_id: str,
    assessment: QuadratureSensitivityAssessment,
) -> NumericalQualificationRecord:
    batch_contract.validate()
    identifier = str(run_id).strip()
    if not identifier:
        raise ValueError("run_id must be nonempty")
    if not _is_sha256(assessment.problem_digest):
        raise ValueError("assessment problem_digest must be a lowercase SHA256 value")
    contract_digest = batch_contract_sha256(batch_contract)
    assessment_digest = _domain_sha256(
        "MFRM_FREE_SD_QUADRATURE_SENSITIVITY_ASSESSMENT",
        "v1",
        assessment.to_dict(),
    )
    record_payload = {
        "schema_version": "mml-numerical-qualification-record-v2",
        "domain": "MFRM_FREE_SD_QFIT_QCHECK_NUMERICAL_RECORD",
        "batch_id": batch_contract.batch_id,
        "scope": batch_contract.scope,
        "registration_identity_sha256": (
            batch_contract.registration_identity_sha256
        ),
        "run_id": identifier,
        "problem_digest": assessment.problem_digest,
        "batch_contract_sha256": contract_digest,
        "assessment_sha256": assessment_digest,
        "objective_evaluator_implementation_sha256": (
            batch_contract.objective_evaluator_implementation_sha256
        ),
    }
    return NumericalQualificationRecord(
        batch_id=batch_contract.batch_id,
        scope=batch_contract.scope,
        registration_identity_sha256=batch_contract.registration_identity_sha256,
        run_id=identifier,
        problem_digest=assessment.problem_digest,
        batch_contract_sha256=contract_digest,
        assessment_sha256=assessment_digest,
        record_sha256=_canonical_sha256(record_payload),
        assessment=assessment,
    )


def _stationarity_contract(
    assessment: StationarityAssessment,
    label: str,
) -> StationarityContract:
    try:
        contract = StationarityContract(**assessment.contract)
        contract.validate()
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} stationarity contract is malformed") from exc
    return contract


def _same_stationarity_assessment(
    stored: StationarityAssessment,
    reconstructed: StationarityAssessment,
) -> bool:
    return _domain_sha256(
        "MFRM_FREE_SD_STATIONARITY_ASSESSMENT",
        "v1",
        stored.to_dict(),
    ) == _domain_sha256(
        "MFRM_FREE_SD_STATIONARITY_ASSESSMENT",
        "v1",
        reconstructed.to_dict(),
    )


def _validate_quadrature_assessment_summary(
    assessment: QuadratureSensitivityAssessment,
    objective_evaluator: QuadratureObjectiveEvaluator,
    *,
    expected_stationarity_contract_sha256: str,
    expected_sensitivity_contract_sha256: str,
    expected_optimizer_options_sha256: str,
) -> None:
    if not _is_sha256(assessment.problem_digest):
        raise ValueError("quadrature assessment problem digest is malformed")
    objective_evaluator.validate()
    if objective_evaluator.problem_digest != assessment.problem_digest:
        raise ValueError("objective evaluator problem digest differs from assessment")
    if (
        objective_evaluator.primary_quadrature_points
        != assessment.primary_quadrature_points
        or objective_evaluator.sensitivity_quadrature_points
        != assessment.sensitivity_quadrature_points
    ):
        raise ValueError("objective evaluator quadrature bases differ from assessment")
    try:
        sensitivity_contract = QuadratureSensitivityContract(**assessment.contract)
        sensitivity_contract.validate()
    except (TypeError, ValueError) as exc:
        raise ValueError("quadrature assessment contract is malformed") from exc
    if contract_sha256(sensitivity_contract) != expected_sensitivity_contract_sha256:
        raise ValueError("quadrature assessment sensitivity contract differs from batch")
    polish_options = (
        assessment.primary_run.primary_polish.options,
        assessment.primary_run.restart_polish.options,
        assessment.sensitivity_run.primary_polish.options,
        assessment.sensitivity_run.restart_polish.options,
    )
    if any(
        contract_sha256(
            options,
            domain="MFRM_FREE_SD_JOINT_POLISH_OPTIONS",
        )
        != expected_optimizer_options_sha256
        for options in polish_options
    ):
        raise ValueError("quadrature assessment optimizer options differ from batch")

    primary_contract = _stationarity_contract(
        assessment.primary_stationarity,
        "primary",
    )
    sensitivity_stationarity_contract = _stationarity_contract(
        assessment.sensitivity_stationarity,
        "sensitivity",
    )
    if (
        contract_sha256(primary_contract) != expected_stationarity_contract_sha256
        or contract_sha256(sensitivity_stationarity_contract)
        != expected_stationarity_contract_sha256
    ):
        raise ValueError("quadrature assessment stationarity contract differs from batch")
    reconstructed_primary = assess_free_sd_stationarity(
        assessment.primary_run,
        primary_contract,
    )
    reconstructed_sensitivity = assess_free_sd_stationarity(
        assessment.sensitivity_run,
        sensitivity_stationarity_contract,
    )
    if (
        not _same_stationarity_assessment(
            assessment.primary_stationarity,
            reconstructed_primary,
        )
        or not _same_stationarity_assessment(
            assessment.sensitivity_stationarity,
            reconstructed_sensitivity,
        )
    ):
        raise ValueError("nested stationarity assessment does not reconstruct")

    observations = float(assessment.primary_run.observations)
    if (
        not np.isfinite(observations)
        or observations <= 0
        or not _consistent(observations, assessment.sensitivity_run.observations)
        or not _consistent(observations, assessment.observations)
    ):
        raise ValueError("quadrature assessment observation denominator is inconsistent")
    primary_parameters = np.asarray(
        assessment.primary_run.restart_polish.structural_parameters,
        dtype=float,
    )
    sensitivity_parameters = np.asarray(
        assessment.sensitivity_run.restart_polish.structural_parameters,
        dtype=float,
    )
    if (
        primary_parameters.ndim != 1
        or primary_parameters.shape != sensitivity_parameters.shape
        or primary_parameters.size == 0
        or not np.isfinite(primary_parameters).all()
        or not np.isfinite(sensitivity_parameters).all()
    ):
        raise ValueError("quadrature assessment structural evidence is malformed")
    primary_sigma = float(assessment.primary_run.restart_polish.sigma)
    sensitivity_sigma = float(assessment.sensitivity_run.restart_polish.sigma)
    if (
        not np.isfinite(primary_sigma)
        or not np.isfinite(sensitivity_sigma)
        or primary_sigma <= 0
        or sensitivity_sigma <= 0
    ):
        raise ValueError("quadrature assessment sigma evidence is malformed")

    try:
        primary_at_primary = float(
            objective_evaluator.primary_value_function(
                primary_parameters,
                primary_sigma,
            )
        )
        primary_at_sensitivity = float(
            objective_evaluator.primary_value_function(
                sensitivity_parameters,
                sensitivity_sigma,
            )
        )
        sensitivity_at_primary = float(
            objective_evaluator.sensitivity_value_function(
                primary_parameters,
                primary_sigma,
            )
        )
        sensitivity_at_sensitivity = float(
            objective_evaluator.sensitivity_value_function(
                sensitivity_parameters,
                sensitivity_sigma,
            )
        )
    except Exception as exc:
        raise ValueError("objective cross-evaluation replay failed") from exc
    objective_values = np.array(
        [
            primary_at_primary,
            primary_at_sensitivity,
            sensitivity_at_primary,
            sensitivity_at_sensitivity,
        ],
        dtype=float,
    )
    finite = bool(np.isfinite(objective_values).all())
    primary_reconstruction = abs(
        primary_at_primary - float(assessment.primary_run.final_objective)
    )
    sensitivity_reconstruction = abs(
        sensitivity_at_sensitivity - float(assessment.sensitivity_run.final_objective)
    )
    structural_difference = float(
        np.max(np.abs(primary_parameters - sensitivity_parameters))
    )
    log_sigma_difference = abs(float(np.log(primary_sigma) - np.log(sensitivity_sigma)))
    sensitivity_gain = float(
        (sensitivity_at_primary - sensitivity_at_sensitivity) / observations
    )
    primary_back_change = float(
        (primary_at_sensitivity - primary_at_primary) / observations
    )
    structural_pass = bool(
        np.isfinite(structural_difference)
        and structural_difference
        <= sensitivity_contract.max_structural_parameter_difference
    )
    sigma_pass = bool(
        np.isfinite(log_sigma_difference)
        and log_sigma_difference <= sensitivity_contract.max_log_sigma_difference
    )
    gain_pass = bool(
        np.isfinite(sensitivity_gain)
        and sensitivity_gain
        >= -sensitivity_contract.max_negative_sensitivity_gain_per_observation
        and sensitivity_gain
        <= sensitivity_contract.max_sensitivity_optimization_gain_per_observation
    )
    back_pass = bool(
        np.isfinite(primary_back_change)
        and abs(primary_back_change)
        <= sensitivity_contract.max_abs_primary_back_evaluation_change_per_observation
    )
    reconstruction_pass = bool(
        finite
        and primary_reconstruction
        <= sensitivity_contract.max_objective_reconstruction_disagreement
        and sensitivity_reconstruction
        <= sensitivity_contract.max_objective_reconstruction_disagreement
    )
    numerical_pass = all(
        (
            finite,
            reconstructed_primary.stationarity_pass,
            reconstructed_sensitivity.stationarity_pass,
            structural_pass,
            sigma_pass,
            gain_pass,
            back_pass,
            reconstruction_pass,
        )
    )
    expected_status = (
        "NUMERICAL_SENSITIVITY_PASS_REGISTRATION_REQUIRED"
        if numerical_pass
        else "NUMERICAL_SENSITIVITY_NOT_QUALIFIED"
    )
    numeric_checks = (
        _consistent(
            assessment.primary_objective_at_primary_solution,
            primary_at_primary,
        ),
        _consistent(
            assessment.primary_objective_at_sensitivity_solution,
            primary_at_sensitivity,
        ),
        _consistent(
            assessment.sensitivity_objective_at_primary_solution,
            sensitivity_at_primary,
        ),
        _consistent(
            assessment.sensitivity_objective_at_sensitivity_solution,
            sensitivity_at_sensitivity,
        ),
        _consistent(
            assessment.primary_objective_reconstruction_difference,
            primary_reconstruction,
        ),
        _consistent(
            assessment.sensitivity_objective_reconstruction_difference,
            sensitivity_reconstruction,
        ),
        _consistent(
            assessment.maximum_structural_parameter_difference,
            structural_difference,
        ),
        _consistent(assessment.absolute_log_sigma_difference, log_sigma_difference),
        _consistent(
            assessment.sensitivity_optimization_gain_per_observation,
            sensitivity_gain,
        ),
        _consistent(
            assessment.primary_back_evaluation_change_per_observation,
            primary_back_change,
        ),
    )
    boolean_checks = (
        _strict_bool(assessment.finite_cross_evaluation, "finite_cross_evaluation")
        == finite,
        _strict_bool(assessment.primary_stationarity_pass, "primary_stationarity_pass")
        == reconstructed_primary.stationarity_pass,
        _strict_bool(
            assessment.sensitivity_stationarity_pass,
            "sensitivity_stationarity_pass",
        )
        == reconstructed_sensitivity.stationarity_pass,
        _strict_bool(assessment.structural_difference_pass, "structural_difference_pass")
        == structural_pass,
        _strict_bool(assessment.log_sigma_difference_pass, "log_sigma_difference_pass")
        == sigma_pass,
        _strict_bool(assessment.sensitivity_gain_pass, "sensitivity_gain_pass")
        == gain_pass,
        _strict_bool(
            assessment.primary_back_evaluation_pass,
            "primary_back_evaluation_pass",
        )
        == back_pass,
        _strict_bool(
            assessment.objective_reconstruction_pass,
            "objective_reconstruction_pass",
        )
        == reconstruction_pass,
        _strict_bool(assessment.numerical_sensitivity_pass, "numerical_sensitivity_pass")
        == numerical_pass,
    )
    checks = (
        assessment.primary_quadrature_points
        == sensitivity_contract.primary_quadrature_points,
        assessment.sensitivity_quadrature_points
        == sensitivity_contract.sensitivity_quadrature_points,
        *numeric_checks,
        *boolean_checks,
        assessment.scientific_inference_ready is False,
        reconstructed_primary.inference_ready is False,
        reconstructed_sensitivity.inference_ready is False,
        assessment.status == expected_status,
    )
    if not all(checks):
        raise ValueError("quadrature assessment summary does not reconstruct")


def assess_numerical_qualification_batch(
    contract: NumericalQualificationBatchContract,
    records: Iterable[NumericalQualificationRecord],
    *,
    evaluator_resolver: ObjectiveEvaluatorResolver,
) -> NumericalQualificationBatchAssessment:
    contract.validate()
    values = list(records)
    identifiers = [str(record.run_id) for record in values]
    if len(set(identifiers)) != len(identifiers):
        raise ValueError("duplicate returned RunId evidence")
    problem_digests = [str(record.problem_digest) for record in values]
    if len(set(problem_digests)) != len(problem_digests):
        raise ValueError("duplicate returned problem evidence")
    planned_problem_map = dict(contract.planned_problem_digests)
    contract_digest = batch_contract_sha256(contract)
    digest_map: dict[str, str] = {}
    for record in values:
        if (
            record.batch_id != contract.batch_id
            or record.scope != contract.scope
            or record.registration_identity_sha256
            != contract.registration_identity_sha256
        ):
            raise ValueError(f"record batch lineage mismatch for {record.run_id}")
        if record.batch_contract_sha256 != contract_digest:
            raise ValueError(f"record batch-contract digest mismatch for {record.run_id}")
        if not _is_sha256(record.problem_digest):
            raise ValueError(f"invalid problem digest for {record.run_id}")
        if record.problem_digest != record.assessment.problem_digest:
            raise ValueError(f"record problem digest mismatch for {record.run_id}")
        if record.run_id in planned_problem_map and (
            record.problem_digest != planned_problem_map[record.run_id]
        ):
            raise ValueError(f"planned problem digest mismatch for {record.run_id}")
        try:
            objective_evaluator = evaluator_resolver(
                record.run_id,
                record.problem_digest,
                contract.objective_evaluator_implementation_sha256,
            )
        except Exception as exc:
            raise ValueError(
                f"auditor-side objective evaluator resolution failed for {record.run_id}"
            ) from exc
        if not isinstance(objective_evaluator, QuadratureObjectiveEvaluator):
            raise ValueError(
                f"auditor-side evaluator has an invalid type for {record.run_id}"
            )
        _validate_quadrature_assessment_summary(
            record.assessment,
            objective_evaluator,
            expected_stationarity_contract_sha256=(
                contract.stationarity_contract_sha256
            ),
            expected_sensitivity_contract_sha256=(
                contract.sensitivity_contract_sha256
            ),
            expected_optimizer_options_sha256=contract.optimizer_options_sha256,
        )
        if not _is_sha256(record.assessment_sha256):
            raise ValueError(f"invalid assessment digest for {record.run_id}")
        reconstructed_assessment = _domain_sha256(
            "MFRM_FREE_SD_QUADRATURE_SENSITIVITY_ASSESSMENT",
            "v1",
            record.assessment.to_dict(),
        )
        if reconstructed_assessment != record.assessment_sha256:
            raise ValueError(f"assessment digest mismatch for {record.run_id}")
        reconstructed_record = _canonical_sha256(
            {
                "schema_version": "mml-numerical-qualification-record-v2",
                "domain": "MFRM_FREE_SD_QFIT_QCHECK_NUMERICAL_RECORD",
                "batch_id": contract.batch_id,
                "scope": contract.scope,
                "registration_identity_sha256": contract.registration_identity_sha256,
                "run_id": record.run_id,
                "problem_digest": record.problem_digest,
                "batch_contract_sha256": contract_digest,
                "assessment_sha256": reconstructed_assessment,
                "objective_evaluator_implementation_sha256": (
                    contract.objective_evaluator_implementation_sha256
                ),
            }
        )
        if not _is_sha256(record.record_sha256) or (
            reconstructed_record != record.record_sha256
        ):
            raise ValueError(f"record lineage digest mismatch for {record.run_id}")
        if (
            record.assessment.primary_quadrature_points
            != contract.primary_quadrature_points
            or record.assessment.sensitivity_quadrature_points
            != contract.sensitivity_quadrature_points
        ):
            raise ValueError(f"quadrature basis mismatch for {record.run_id}")
        digest_map[record.run_id] = reconstructed_record

    planned = set(contract.planned_run_ids)
    returned = set(identifiers)
    missing = tuple(sorted(planned - returned))
    unexpected = tuple(sorted(returned - planned))
    in_plan = [record for record in values if record.run_id in planned]
    failed = tuple(
        sorted(
            record.run_id
            for record in in_plan
            if not record.assessment.numerical_sensitivity_pass
        )
    )
    exact_denominator = bool(not missing and not unexpected and len(in_plan) == len(planned))
    all_pass = bool(
        exact_denominator
        and all(record.assessment.numerical_sensitivity_pass for record in in_plan)
    )
    numerical_batch_pass = bool(exact_denominator and all_pass)
    return NumericalQualificationBatchAssessment(
        batch_id=contract.batch_id,
        scope=contract.scope,
        registration_identity_sha256=contract.registration_identity_sha256,
        batch_contract_sha256=contract_digest,
        planned_run_ids_sha256=contract.planned_run_ids_sha256,
        planned_problem_digests_sha256=contract.planned_problem_digests_sha256,
        stationarity_contract_sha256=contract.stationarity_contract_sha256,
        sensitivity_contract_sha256=contract.sensitivity_contract_sha256,
        optimizer_options_sha256=contract.optimizer_options_sha256,
        objective_evaluator_implementation_sha256=(
            contract.objective_evaluator_implementation_sha256
        ),
        planned_datasets=len(planned),
        returned_datasets=len(values),
        finite_cross_evaluation_datasets=sum(
            bool(record.assessment.finite_cross_evaluation) for record in in_plan
        ),
        primary_stationarity_pass_datasets=sum(
            bool(record.assessment.primary_stationarity_pass) for record in in_plan
        ),
        sensitivity_stationarity_pass_datasets=sum(
            bool(record.assessment.sensitivity_stationarity_pass) for record in in_plan
        ),
        numerical_sensitivity_pass_datasets=sum(
            bool(record.assessment.numerical_sensitivity_pass) for record in in_plan
        ),
        missing_run_ids=missing,
        unexpected_run_ids=unexpected,
        failed_run_ids=failed,
        exact_denominator_pass=exact_denominator,
        all_numerical_sensitivity_pass=all_pass,
        numerical_batch_pass=numerical_batch_pass,
        scientific_inference_ready=False,
        status=(
            (
                "NUMERICAL_BATCH_PASS_DEVELOPMENT_ONLY"
                if contract.scope == "DEVELOPMENT_ONLY"
                else "NUMERICAL_BATCH_PASS_PROSPECTIVE_REGISTRATION_REQUIRED"
            )
            if numerical_batch_pass
            else "NUMERICAL_BATCH_NOT_QUALIFIED"
        ),
        record_sha256=dict(sorted(digest_map.items())),
        contract=asdict(contract),
    )


__all__ = [
    "NumericalQualificationBatchAssessment",
    "NumericalQualificationBatchContract",
    "NumericalQualificationRecord",
    "ObjectiveEvaluatorResolver",
    "QuadratureObjectiveEvaluator",
    "assess_numerical_qualification_batch",
    "batch_contract_sha256",
    "contract_sha256",
    "make_numerical_qualification_record",
    "planned_problem_digests_sha256",
    "planned_run_ids_sha256",
]
