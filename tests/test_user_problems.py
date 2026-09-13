from __future__ import annotations

from dataclasses import FrozenInstanceError
import json
from pathlib import Path
import re

import pytest

from mfrm_app import user_problems as problems


EXPECTED_CODES = {
    "input.parse_failed",
    "input.mapping_invalid",
    "input.no_usable_rows",
    "estimation.identification_failed",
    "estimation.nonconvergence",
    "estimation.rating_scale_invalid",
    "estimation.anchor_invalid",
    "resource.hosted_limit",
    "diagnostic.residual_pca_unavailable",
    "diagnostic.residual_pca_failed",
    "problem.unexpected",
}


def test_registry_has_every_initial_code_and_only_stable_targets() -> None:
    problems.validate_user_problem_registry()

    assert set(problems.USER_PROBLEM_CODES) == EXPECTED_CODES
    assert set(problems.USER_PROBLEM_SPECS) == EXPECTED_CODES
    for code, spec in problems.USER_PROBLEM_SPECS.items():
        assert spec.problem_code == code
        assert spec.title_key == f"problems.{code}.title"
        assert spec.body_key == f"problems.{code}.body"
        assert spec.action_keys
        assert len(spec.action_keys) == len(spec.action_target_ids)
        assert all(key.startswith("problems.actions.") for key in spec.action_keys)
        assert all(target.startswith("target.") for target in spec.action_target_ids)
        assert spec.help_topic_id.startswith("help.")


def test_specs_and_notices_are_frozen_and_notices_store_only_code_and_reference(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    spec = problems.get_user_problem_spec("input.parse_failed")
    monkeypatch.setattr(
        problems.secrets,
        "token_bytes",
        lambda size: b"\x00\x11\x22\x33\x44\x55\x66\x77\x88\x99\xaa\xbb",
    )
    notice = problems.build_user_problem_notice(spec.problem_code)

    with pytest.raises(FrozenInstanceError):
        spec.body_key = "problems.problem.unexpected.body"
    with pytest.raises((FrozenInstanceError, AttributeError)):
        notice.problem_code = "problem.unexpected"
    with pytest.raises((FrozenInstanceError, AttributeError, TypeError)):
        notice.exception = RuntimeError("must not attach")

    assert notice.problem_code == "input.parse_failed"
    assert notice.title_key == "problems.input.parse_failed.title"
    assert notice.__slots__ == (
        "problem_code",
        "occurrence_phase",
        "support_reference",
    )

    with pytest.raises(problems.UserProblemContractError, match="build_user_problem_notice"):
        problems.UserProblemNotice()
    with pytest.raises(TypeError):
        problems.UserProblemNotice(  # type: ignore[call-arg]
            problem_code="input.parse_failed",
            support_reference="MFRM-00112233-44556677-8899AABB",
        )


@pytest.mark.parametrize(
    ("message", "expected"),
    [
        ("could not decode the uploaded bytes", "input.parse_failed"),
        ("required column is missing from column mapping", "input.mapping_invalid"),
        ("empty dataframe after filtering; no usable rows", "input.no_usable_rows"),
    ],
)
def test_parse_classifier_returns_stable_codes(message: str, expected: str) -> None:
    assert problems.classify_parse_exception(ValueError(message)).problem_code == expected


@pytest.mark.parametrize(
    ("exception", "expected"),
    [
        (RuntimeError("information matrix is singular"), "estimation.identification_failed"),
        (RuntimeError("maximum iterations reached; did not converge"), "estimation.nonconvergence"),
        (ValueError("rating scale has non-contiguous score categories"), "estimation.rating_scale_invalid"),
        (ValueError("invalid anchor row in anchor table"), "estimation.anchor_invalid"),
        (MemoryError("out of memory"), "resource.hosted_limit"),
    ],
)
def test_estimation_classifier_returns_specific_registered_codes(
    exception: BaseException,
    expected: str,
) -> None:
    assert problems.classify_estimation_exception(exception).problem_code == expected


def test_residual_pca_and_resource_classifiers_keep_unavailable_failure_separate() -> None:
    unavailable = problems.classify_residual_pca_exception(
        RuntimeError("too few comparable residual profiles")
    )
    failed = problems.classify_residual_pca_exception(
        ArithmeticError("eigendecomposition failed")
    )
    hosted = problems.classify_resource_exception(MemoryError("capacity exceeded"))
    unrelated = problems.classify_resource_exception(RuntimeError("unclassified event"))
    numerical_failure = problems.classify_residual_pca_exception(
        ArithmeticError("insufficient workspace for eigendecomposition")
    )

    assert unavailable.problem_code == "diagnostic.residual_pca_unavailable"
    assert unavailable.severity is problems.UserProblemSeverity.INFORMATION
    assert failed.problem_code == "diagnostic.residual_pca_failed"
    assert failed.severity is problems.UserProblemSeverity.CAUTION
    assert hosted.problem_code == "resource.hosted_limit"
    assert unrelated.problem_code == "problem.unexpected"
    assert numerical_failure.problem_code == "diagnostic.residual_pca_failed"
    assert numerical_failure.severity is problems.UserProblemSeverity.CAUTION


def test_unknown_estimation_and_application_exceptions_use_generic_fallback() -> None:
    exception = RuntimeError("a condition with no registered match")

    assert (
        problems.classify_estimation_exception(exception).problem_code
        == "problem.unexpected"
    )
    assert (
        problems.classify_user_problem(
            exception,
            phase=problems.UserProblemPhase.APPLICATION,
        ).problem_code
        == "problem.unexpected"
    )

    # A bare lookup failure during fitting can be an implementation defect;
    # do not misdirect the user to column mapping without a column marker.
    assert (
        problems.classify_estimation_exception(KeyError("unclassified lookup"))
        .problem_code
        == "problem.unexpected"
    )
    assert (
        problems.classify_estimation_exception(TimeoutError("request timed out"))
        .problem_code
        == "problem.unexpected"
    )
    assert (
        problems.classify_parse_exception(KeyError("person"))
        .problem_code
        == "input.mapping_invalid"
    )


def test_support_reference_is_random_format_without_public_value_injection(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    expected = "MFRM-00112233-44556677-8899AABB"
    requested_sizes: list[int] = []

    def fixed_token(size: int) -> bytes:
        requested_sizes.append(size)
        return b"\x00\x11\x22\x33\x44\x55\x66\x77\x88\x99\xaa\xbb"

    monkeypatch.setattr(problems.secrets, "token_bytes", fixed_token)
    generated = problems.generate_support_reference()

    assert generated == expected
    assert requested_sizes == [12]
    assert problems.validate_support_reference(generated) == generated

    with pytest.raises(TypeError):
        problems.generate_support_reference(  # type: ignore[call-arg]
            test_token=b"\x00\x11\x22\x33\x44\x55\x66\x77\x88\x99\xaa\xbb"
        )
    with pytest.raises(TypeError):
        problems.build_user_problem_notice(  # type: ignore[call-arg]
            "input.parse_failed",
            support_reference=expected,
        )
    with pytest.raises(problems.UserProblemContractError):
        problems.validate_support_reference("MFRM-/Users/private/table.csv")


def test_support_references_are_independent_between_notices(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    tokens = iter(
        (
            b"\x00\x11\x22\x33\x44\x55\x66\x77\x88\x99\xaa\xbb",
            b"\xff\xee\xdd\xcc\xbb\xaa\x99\x88\x77\x66\x55\x44",
        )
    )
    monkeypatch.setattr(problems.secrets, "token_bytes", lambda size: next(tokens))

    first = problems.generate_support_reference()
    second = problems.generate_support_reference()

    assert re.fullmatch(r"MFRM-[0-9A-F]{8}-[0-9A-F]{8}-[0-9A-F]{8}", first)
    assert first != second


def test_notice_payload_log_context_and_repr_never_retain_exception_details(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    class ConfidentialPipelineFailure(RuntimeError):
        pass

    sentinels = (
        "Aiko Tanaka",
        "student-secret@example.org",
        "/Users/aiko/private/ratings.csv",
        "ratings.csv",
        "LearnerEmail",
        "raw-value-9917",
        "MFRM_SHOW_TECHNICAL_ERRORS",
        "ana_secret_AnalysisID_8877",
        "ConfidentialPipelineFailure",
    )
    exception = ConfidentialPipelineFailure(
        "Aiko Tanaka student-secret@example.org at "
        "/Users/aiko/private/ratings.csv filename=ratings.csv "
        "column=LearnerEmail value=raw-value-9917 "
        "env=MFRM_SHOW_TECHNICAL_ERRORS AnalysisID=ana_secret_AnalysisID_8877"
    )
    monkeypatch.setattr(
        problems.secrets,
        "token_bytes",
        lambda size: bytes.fromhex("abcdef0123456789abcdef01"),
    )
    notice = problems.notice_from_exception(
        exception,
        phase=problems.UserProblemPhase.ESTIMATION,
    )

    payload = notice.to_payload()
    serialized_surfaces = (
        repr(notice),
        repr(payload),
        notice.to_json(),
        repr(dict(problems.safe_log_context(notice))),
    )
    joined = "\n".join(serialized_surfaces)

    assert notice.problem_code == "problem.unexpected"
    assert notice.phase is problems.UserProblemPhase.ESTIMATION
    assert set(payload) == {
        "schema_version",
        "problem_code",
        "severity",
        "phase",
        "title_key",
        "body_key",
        "action_keys",
        "action_target_ids",
        "help_topic_id",
        "support_reference",
    }
    assert json.loads(notice.to_json()) == payload
    for sentinel in sentinels:
        assert sentinel not in joined


@pytest.mark.parametrize(
    "phase",
    (
        problems.UserProblemPhase.ESTIMATION,
        problems.UserProblemPhase.RESOURCE,
        problems.UserProblemPhase.APPLICATION,
    ),
)
def test_unexpected_notice_preserves_only_the_allowlisted_occurrence_phase(
    phase: problems.UserProblemPhase,
) -> None:
    notice = problems.notice_from_exception(
        RuntimeError("unclassified confidential failure"),
        phase=phase,
    )

    assert notice.problem_code == "problem.unexpected"
    assert notice.phase is phase
    assert notice.to_payload()["phase"] == phase.value
    assert problems.safe_log_context(notice)["phase"] == phase.value


def test_known_problem_rejects_a_mismatched_occurrence_phase() -> None:
    with pytest.raises(problems.UserProblemContractError, match="must match"):
        problems.build_user_problem_notice(
            "estimation.nonconvergence",
            occurrence_phase="parse",
        )


@pytest.mark.parametrize(
    ("exception", "expected_code"),
    (
        (MemoryError("out of memory"), "resource.hosted_limit"),
        (ValueError("required column is missing"), "input.mapping_invalid"),
        (ValueError("no usable rows after filtering"), "input.no_usable_rows"),
    ),
)
def test_estimation_adapter_preserves_occurrence_for_cross_phase_classification(
    exception: BaseException,
    expected_code: str,
) -> None:
    notice = problems.notice_from_exception(exception, phase="estimation")

    assert notice.problem_code == expected_code
    assert notice.phase is problems.UserProblemPhase.ESTIMATION
    assert notice.to_payload()["phase"] == "estimation"


def test_support_reference_does_not_depend_on_exception_or_data_identity(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    tokens = iter((bytes.fromhex("00112233445566778899aabb"), bytes.fromhex("ffeeddccbbaa998877665544")))
    requested_sizes: list[int] = []

    def next_token(size: int) -> bytes:
        requested_sizes.append(size)
        return next(tokens)

    monkeypatch.setattr(problems.secrets, "token_bytes", next_token)
    first = problems.notice_from_exception(
        RuntimeError("secret-message-one data_hash=111 AnalysisID=ana_one"),
        phase="application",
    )
    second = problems.notice_from_exception(
        RuntimeError("secret-message-two data_hash=222 AnalysisID=ana_two"),
        phase="application",
    )

    assert first.problem_code == second.problem_code == "problem.unexpected"
    assert first.support_reference != second.support_reference
    assert first.support_reference == "MFRM-00112233-44556677-8899AABB"
    assert second.support_reference == "MFRM-FFEEDDCC-BBAA9988-77665544"
    assert requested_sizes == [12, 12]


def test_safe_log_context_is_allowlisted_and_read_only(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        problems.secrets,
        "token_bytes",
        lambda size: bytes.fromhex("0123456789abcdef01234567"),
    )
    notice = problems.build_user_problem_notice("estimation.nonconvergence")
    context = problems.safe_log_context(notice)

    assert context == {
        "event": "mfrm.user_problem",
        "schema_version": "mfrm_user_problem_v1",
        "problem_code": "estimation.nonconvergence",
        "severity": "blocked",
        "phase": "estimation",
        "support_reference": "MFRM-01234567-89ABCDEF-01234567",
    }
    with pytest.raises(TypeError):
        context["exception"] = "must not be added"


def test_module_remains_streamlit_and_pandas_independent() -> None:
    source = Path(problems.__file__).read_text(encoding="utf-8")

    assert "import streamlit" not in source
    assert "import pandas" not in source
    assert "os.environ" not in source
