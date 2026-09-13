from __future__ import annotations

from dataclasses import replace
import hashlib

import numpy as np
import pytest

from mfrm_app.mml_engine_v2 import run_free_sd_stationarity_v2
from mfrm_app.mml_quadrature_sensitivity import assess_quadrature_sensitivity
from mfrm_app.mml_stationarity import JointPolishOptions
from mfrm_app.mml_qualification_batch import (
    NumericalQualificationBatchContract,
    QuadratureObjectiveEvaluator,
    assess_numerical_qualification_batch,
    contract_sha256,
    make_numerical_qualification_record,
    planned_problem_digests_sha256,
    planned_run_ids_sha256,
)
from tests.test_mml_quadrature_sensitivity import (
    objective,
    sensitivity_contract,
    stationarity_contract,
)


RUN_IDS = ("gamma-0__vector-01", "gamma-0__vector-02", "gamma-0__vector-03")
PROBLEM_DIGESTS = tuple(
    hashlib.sha256(f"qualification-problem::{run_id}".encode("utf-8")).hexdigest()
    for run_id in RUN_IDS
)
PROBLEM_MAP = tuple(zip(RUN_IDS, PROBLEM_DIGESTS, strict=True))
EVALUATOR_IMPLEMENTATION_SHA256 = hashlib.sha256(
    b"qualification-test-objective-evaluator-v1"
).hexdigest()


def contract() -> NumericalQualificationBatchContract:
    stationary = stationarity_contract()
    sensitivity = sensitivity_contract()
    return NumericalQualificationBatchContract(
        batch_id="DEVELOPMENT_BATCH_NOT_REGISTERED",
        scope="DEVELOPMENT_ONLY",
        registration_identity_sha256=None,
        planned_run_ids=RUN_IDS,
        planned_run_ids_sha256=planned_run_ids_sha256(RUN_IDS),
        planned_problem_digests=PROBLEM_MAP,
        planned_problem_digests_sha256=planned_problem_digests_sha256(PROBLEM_MAP),
        stationarity_contract_sha256=contract_sha256(stationary),
        sensitivity_contract_sha256=contract_sha256(sensitivity),
        optimizer_options_sha256=contract_sha256(JointPolishOptions()),
        objective_evaluator_implementation_sha256=(
            EVALUATOR_IMPLEMENTATION_SHA256
        ),
        primary_quadrature_points=31,
        sensitivity_quadrature_points=61,
        all_planned_required=True,
        replacement_forbidden=True,
    )


def _objectives(index: int, structural_shift: float):
    primary_target = np.array([0.5 + 0.02 * index, -0.25 - 0.01 * index])
    sensitivity_target = primary_target + np.array([structural_shift, -1e-4])
    primary_log_sigma = float(np.log(1.2 + 0.01 * index))
    sensitivity_log_sigma = primary_log_sigma + 8e-5
    primary_value, primary_gradient = objective(primary_target, primary_log_sigma)
    sensitivity_value, sensitivity_gradient = objective(
        sensitivity_target,
        sensitivity_log_sigma,
    )
    return primary_value, primary_gradient, sensitivity_value, sensitivity_gradient


def _evaluator(
    index: int,
    *,
    structural_shift: float = 2e-4,
    problem_digest: str | None = None,
) -> QuadratureObjectiveEvaluator:
    primary_value, _primary_gradient, sensitivity_value, _sensitivity_gradient = (
        _objectives(index, structural_shift)
    )
    return QuadratureObjectiveEvaluator(
        problem_digest=problem_digest or PROBLEM_DIGESTS[index],
        primary_quadrature_points=31,
        sensitivity_quadrature_points=61,
        primary_value_function=primary_value,
        sensitivity_value_function=sensitivity_value,
    )


def _assessment(
    index: int,
    *,
    structural_shift: float = 2e-4,
    stationarity=None,
    sensitivity=None,
    problem_digest: str | None = None,
):
    (
        primary_value,
        primary_gradient,
        sensitivity_value,
        sensitivity_gradient,
    ) = _objectives(
        index,
        structural_shift,
    )
    primary_run = run_free_sd_stationarity_v2(
        np.array([-1.0, 1.0]),
        0.8,
        primary_value,
        primary_gradient,
        observations=100,
        constraint_residual_function=lambda _parameters: 0.0,
    )
    sensitivity_run = run_free_sd_stationarity_v2(
        primary_run.restart_polish.structural_parameters,
        primary_run.restart_polish.sigma,
        sensitivity_value,
        sensitivity_gradient,
        observations=100,
        constraint_residual_function=lambda _parameters: 0.0,
    )
    return assess_quadrature_sensitivity(
        problem_digest=problem_digest or PROBLEM_DIGESTS[index],
        primary_quadrature_points=31,
        sensitivity_quadrature_points=61,
        primary_run=primary_run,
        sensitivity_run=sensitivity_run,
        primary_value_function=primary_value,
        sensitivity_value_function=sensitivity_value,
        stationarity_contract=stationarity or stationarity_contract(),
        sensitivity_contract=sensitivity or sensitivity_contract(),
    )


@pytest.fixture(scope="module")
def passing_assessments():
    return tuple(_assessment(index) for index in range(len(RUN_IDS)))


def records(assessments):
    batch_contract = contract()
    return [
        make_numerical_qualification_record(
            batch_contract,
            run_id,
            assessment,
        )
        for run_id, assessment in zip(RUN_IDS, assessments, strict=True)
    ]


def assess_batch(batch_contract, values, *, evaluator_overrides=None):
    evaluator_map = {
        run_id: _evaluator(index)
        for index, run_id in enumerate(RUN_IDS)
    }
    evaluator_map.update(evaluator_overrides or {})

    def resolver(run_id, problem_digest, implementation_sha256):
        if implementation_sha256 != EVALUATOR_IMPLEMENTATION_SHA256:
            raise ValueError("unexpected evaluator implementation")
        evaluator = evaluator_map[run_id]
        if evaluator.problem_digest != problem_digest:
            raise ValueError("unexpected evaluator problem")
        return evaluator

    return assess_numerical_qualification_batch(
        batch_contract,
        values,
        evaluator_resolver=resolver,
    )


def test_all_planned_numerical_pass_still_cannot_emit_scientific_readiness(
    passing_assessments,
) -> None:
    result = assess_batch(
        contract(),
        records(passing_assessments),
    )
    assert result.planned_datasets == 3
    assert result.returned_datasets == 3
    assert result.numerical_sensitivity_pass_datasets == 3
    assert result.exact_denominator_pass is True
    assert result.numerical_batch_pass is True
    assert result.scientific_inference_ready is False
    assert result.status == "NUMERICAL_BATCH_PASS_DEVELOPMENT_ONLY"


def test_missing_case_remains_in_denominator(passing_assessments) -> None:
    result = assess_batch(
        contract(),
        records(passing_assessments)[:-1],
    )
    assert result.planned_datasets == 3
    assert result.returned_datasets == 2
    assert result.missing_run_ids == (RUN_IDS[-1],)
    assert result.exact_denominator_pass is False
    assert result.numerical_batch_pass is False


def test_failed_case_is_not_replaced_by_an_unexpected_success(
    passing_assessments,
) -> None:
    failed = _assessment(1, structural_shift=0.01)
    assert failed.numerical_sensitivity_pass is False
    values = records(passing_assessments)
    values[1] = make_numerical_qualification_record(
        contract(),
        RUN_IDS[1],
        failed,
    )
    unexpected_digest = hashlib.sha256(b"unexpected-problem").hexdigest()
    unexpected = _assessment(20, problem_digest=unexpected_digest)
    values.append(
        make_numerical_qualification_record(
            contract(),
            "replacement-success",
            unexpected,
        )
    )
    result = assess_batch(
        contract(),
        values,
        evaluator_overrides={
            RUN_IDS[1]: _evaluator(1, structural_shift=0.01),
            "replacement-success": _evaluator(
                20,
                problem_digest=unexpected_digest,
            ),
        },
    )
    assert result.failed_run_ids == (RUN_IDS[1],)
    assert result.unexpected_run_ids == ("replacement-success",)
    assert result.exact_denominator_pass is False
    assert result.numerical_batch_pass is False


def test_duplicate_and_tampered_evidence_are_rejected(passing_assessments) -> None:
    value = records(passing_assessments)[0]
    with pytest.raises(ValueError, match="duplicate returned RunId"):
        assess_batch(contract(), [value, value])

    tampered = replace(value, assessment_sha256="0" * 64)
    with pytest.raises(ValueError, match="assessment digest mismatch"):
        assess_batch(contract(), [tampered])

    forged = replace(
        passing_assessments[0],
        maximum_structural_parameter_difference=0.0,
    )
    forged_record = make_numerical_qualification_record(
        contract(),
        RUN_IDS[0],
        forged,
    )
    with pytest.raises(ValueError, match="summary does not reconstruct"):
        assess_batch(contract(), [forged_record])


def test_run_id_is_bound_to_its_planned_problem(passing_assessments) -> None:
    copied = make_numerical_qualification_record(
        contract(),
        RUN_IDS[1],
        passing_assessments[0],
    )
    with pytest.raises(ValueError, match="planned problem digest mismatch"):
        assess_batch(contract(), [copied])

    first = records(passing_assessments)[0]
    relabelled = replace(first, run_id=RUN_IDS[1])
    with pytest.raises(ValueError, match="planned problem digest mismatch|lineage digest"):
        assess_batch(contract(), [relabelled])


def test_batch_rejects_mixed_numerical_contracts(passing_assessments) -> None:
    looser_sensitivity = sensitivity_contract(
        max_structural_parameter_difference=0.5,
    )
    mixed_sensitivity = _assessment(0, sensitivity=looser_sensitivity)
    value = make_numerical_qualification_record(
        contract(),
        RUN_IDS[0],
        mixed_sensitivity,
    )
    with pytest.raises(ValueError, match="sensitivity contract differs"):
        assess_batch(contract(), [value])

    looser_stationarity = replace(
        stationarity_contract(),
        max_projected_gradient_supnorm=1e-4,
    )
    mixed_stationarity = _assessment(0, stationarity=looser_stationarity)
    value = make_numerical_qualification_record(
        contract(),
        RUN_IDS[0],
        mixed_stationarity,
    )
    with pytest.raises(ValueError, match="stationarity contract differs"):
        assess_batch(contract(), [value])

    assessment = passing_assessments[0]
    changed_options = dict(assessment.primary_run.primary_polish.options)
    changed_options["maxiter"] = int(changed_options["maxiter"]) + 1
    changed_run = replace(
        assessment.primary_run,
        primary_polish=replace(
            assessment.primary_run.primary_polish,
            options=changed_options,
        ),
    )
    mixed_options = replace(assessment, primary_run=changed_run)
    value = make_numerical_qualification_record(
        contract(),
        RUN_IDS[0],
        mixed_options,
    )
    with pytest.raises(ValueError, match="optimizer options differ"):
        assess_batch(contract(), [value])

    wrong_basis_evaluator = replace(
        _evaluator(0),
        sensitivity_quadrature_points=81,
    )
    with pytest.raises(ValueError, match="quadrature bases differ"):
        assess_batch(
            contract(),
            [records(passing_assessments)[0]],
            evaluator_overrides={RUN_IDS[0]: wrong_basis_evaluator},
        )

    nonfinite_evaluator = replace(
        _evaluator(0),
        primary_value_function=lambda _parameters, _sigma: float("nan"),
    )
    with pytest.raises(ValueError, match="summary does not reconstruct"):
        assess_batch(
            contract(),
            [records(passing_assessments)[0]],
            evaluator_overrides={RUN_IDS[0]: nonfinite_evaluator},
        )


def test_nested_and_cross_evaluation_summaries_are_reconstructed(
    passing_assessments,
) -> None:
    assessment = passing_assessments[0]
    forged_nested = replace(
        assessment,
        primary_stationarity=replace(
            assessment.primary_stationarity,
            stationarity_pass=False,
        ),
        primary_stationarity_pass=False,
        numerical_sensitivity_pass=False,
        status="NUMERICAL_SENSITIVITY_NOT_QUALIFIED",
    )
    value = make_numerical_qualification_record(
        contract(),
        RUN_IDS[0],
        forged_nested,
    )
    with pytest.raises(ValueError, match="nested stationarity assessment"):
        assess_batch(contract(), [value])

    forged_cross = replace(
        assessment,
        primary_objective_at_sensitivity_solution=(
            assessment.primary_objective_at_sensitivity_solution + 1.0
        ),
    )
    value = make_numerical_qualification_record(
        contract(),
        RUN_IDS[0],
        forged_cross,
    )
    with pytest.raises(ValueError, match="summary does not reconstruct"):
        assess_batch(contract(), [value])

    coordinated_cross_forgery = replace(
        assessment,
        primary_objective_at_sensitivity_solution=(
            assessment.primary_objective_at_primary_solution
        ),
        sensitivity_objective_at_primary_solution=(
            assessment.sensitivity_objective_at_sensitivity_solution
        ),
        sensitivity_optimization_gain_per_observation=0.0,
        primary_back_evaluation_change_per_observation=0.0,
        sensitivity_gain_pass=True,
        primary_back_evaluation_pass=True,
        numerical_sensitivity_pass=True,
        status="NUMERICAL_SENSITIVITY_PASS_REGISTRATION_REQUIRED",
    )
    value = make_numerical_qualification_record(
        contract(),
        RUN_IDS[0],
        coordinated_cross_forgery,
    )
    with pytest.raises(ValueError, match="summary does not reconstruct"):
        assess_batch(contract(), [value])


def test_plan_order_problem_map_and_digests_are_exact(passing_assessments) -> None:
    unsorted_ids = tuple(reversed(RUN_IDS))
    unsorted_map = tuple(reversed(PROBLEM_MAP))
    unsorted = replace(
        contract(),
        planned_run_ids=unsorted_ids,
        planned_run_ids_sha256=planned_run_ids_sha256(unsorted_ids),
        planned_problem_digests=unsorted_map,
        planned_problem_digests_sha256=planned_problem_digests_sha256(unsorted_map),
    )
    with pytest.raises(ValueError, match="canonical sorted order"):
        assess_batch(unsorted, records(passing_assessments))

    bad_digest = replace(contract(), planned_run_ids_sha256="f" * 64)
    with pytest.raises(ValueError, match="RunId digest does not reconstruct"):
        assess_batch(bad_digest, records(passing_assessments))

    duplicate_problem_map = (
        PROBLEM_MAP[0],
        (PROBLEM_MAP[1][0], PROBLEM_MAP[0][1]),
        PROBLEM_MAP[2],
    )
    duplicate_problem = replace(
        contract(),
        planned_problem_digests=duplicate_problem_map,
        planned_problem_digests_sha256=planned_problem_digests_sha256(
            duplicate_problem_map
        ),
    )
    with pytest.raises(ValueError, match="problem digests must be unique"):
        assess_batch(
            duplicate_problem,
            records(passing_assessments),
        )

    fake_registered = replace(
        contract(),
        scope="REGISTERED_Q31_Q61_SUBSET",
        registration_identity_sha256=None,
    )
    with pytest.raises(ValueError, match="only DEVELOPMENT_ONLY"):
        assess_batch(
            fake_registered,
            records(passing_assessments),
        )
