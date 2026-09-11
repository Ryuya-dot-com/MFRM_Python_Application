"""Versioned keyed randomization for prospective MFRM simulation."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
import hashlib
import math
import re

from ._condition_shared import (
    MAX_MONTE_CARLO_REPLICATES,
    SimulationConditionValidationError,
    _fingerprint,
    _json_safe,
    _require_int,
    _require_literal,
    _strict_mapping,
)


RANDOMIZATION_SPEC_VERSION = "mfrm_randomization_spec_v1"

RANDOMIZATION_ALGORITHM_V1 = "sha256_u53_box_muller_v1"
RANDOMIZATION_SHARED_RECORD_SCOPE_TWO_ARM = "two_arm_shared_score_records"
RANDOM_STREAM_PERSON_THETA = "person_theta"
RANDOM_STREAM_RATER_SEVERITY = "rater_severity_raw"
RANDOM_STREAM_RATER_CENTRAL_TENDENCY = "rater_central_tendency"
RANDOM_STREAM_RATER_PERSON_LOCAL_BIAS = "rater_person_local_bias"
RANDOM_STREAM_CRITERION_DIFFICULTY = "criterion_difficulty_raw"
RANDOM_STREAM_RESPONSE = "response_uniform"
RANDOM_STREAM_CHOICES = (
    RANDOM_STREAM_PERSON_THETA,
    RANDOM_STREAM_RATER_SEVERITY,
    RANDOM_STREAM_RATER_CENTRAL_TENDENCY,
    RANDOM_STREAM_RATER_PERSON_LOCAL_BIAS,
    RANDOM_STREAM_CRITERION_DIFFICULTY,
    RANDOM_STREAM_RESPONSE,
)
RANDOM_STREAM_KEY_SCHEMAS = (
    (RANDOM_STREAM_PERSON_THETA, ("person_role", "person_id")),
    (RANDOM_STREAM_RATER_SEVERITY, ("rater_id",)),
    (RANDOM_STREAM_RATER_CENTRAL_TENDENCY, ("rater_id",)),
    (
        RANDOM_STREAM_RATER_PERSON_LOCAL_BIAS,
        ("person_role", "person_id", "rater_id"),
    ),
    (RANDOM_STREAM_CRITERION_DIFFICULTY, ("criterion_id",)),
    (
        RANDOM_STREAM_RESPONSE,
        (
            "person_role",
            "person_id",
            "artifact_id",
            "rater_id",
            "criterion_id",
        ),
    ),
)


@dataclass(frozen=True, slots=True, kw_only=True)
class RandomizationSpecV1:
    """Order-invariant two-arm draw sharing, not inferential comparability."""

    master_seed: int = 20260722
    algorithm: str = RANDOMIZATION_ALGORITHM_V1
    shared_record_scope: str = RANDOMIZATION_SHARED_RECORD_SCOPE_TWO_ARM
    shared_truth: bool = True
    shared_existing_response_draws: bool = True
    schema_version: str = RANDOMIZATION_SPEC_VERSION

    def __post_init__(self) -> None:
        if self.schema_version != RANDOMIZATION_SPEC_VERSION:
            raise SimulationConditionValidationError(
                f"Unsupported randomization version: {self.schema_version!r}"
            )
        _require_int("master_seed", self.master_seed, minimum=0, maximum=2**63 - 1)
        _require_literal("algorithm", self.algorithm, RANDOMIZATION_ALGORITHM_V1)
        _require_literal(
            "shared_record_scope",
            self.shared_record_scope,
            RANDOMIZATION_SHARED_RECORD_SCOPE_TWO_ARM,
        )
        if self.shared_truth is not True:
            raise SimulationConditionValidationError("shared_truth must remain true")
        if self.shared_existing_response_draws is not True:
            raise SimulationConditionValidationError(
                "shared_existing_response_draws must remain true"
            )

    def to_dict(self) -> dict:
        return _json_safe({
            "master_seed": self.master_seed,
            "algorithm": self.algorithm,
            "shared_record_scope": self.shared_record_scope,
            "shared_truth": self.shared_truth,
            "shared_existing_response_draws": self.shared_existing_response_draws,
            "schema_version": self.schema_version,
        })


def normalize_randomization_spec(
    value: RandomizationSpecV1 | Mapping,
) -> RandomizationSpecV1:
    if type(value) is RandomizationSpecV1:
        return value
    return RandomizationSpecV1(
        **_strict_mapping(value, RandomizationSpecV1, label="randomization spec")
    )


def randomization_spec_fingerprint(
    value: RandomizationSpecV1 | Mapping,
    *,
    length: int = 16,
) -> str:
    return _fingerprint(normalize_randomization_spec(value).to_dict(), length=length)


_RANDOM_STREAM_KEY_SCHEMA_MAP = dict(RANDOM_STREAM_KEY_SCHEMAS)
_PERSON_ID_PATTERN = re.compile(r"^(?:P|CA)(?P<index>[0-9]{6})$")
_RATER_ID_PATTERN = re.compile(r"^R(?P<index>[0-9]{6})$")
_CRITERION_ID_PATTERN = re.compile(r"^C(?P<index>[0-9]{3})$")


def _is_positive_canonical_id(value: str, pattern: re.Pattern[str]) -> bool:
    match = pattern.fullmatch(value)
    return match is not None and any(
        character != "0" for character in match.group("index")
    )


def _validate_stream_key_parts(stream: str, key_parts: tuple[str, ...]) -> None:
    field_names = _RANDOM_STREAM_KEY_SCHEMA_MAP[stream]
    if len(key_parts) != len(field_names):
        raise SimulationConditionValidationError(
            f"{stream} requires exactly {len(field_names)} key parts"
        )
    if any(not isinstance(part, str) or not part for part in key_parts):
        raise SimulationConditionValidationError(
            "key_parts must contain non-empty strings"
        )
    try:
        for part in key_parts:
            part.encode("utf-8")
    except UnicodeEncodeError as exc:
        raise SimulationConditionValidationError(
            "key_parts must contain valid UTF-8 text"
        ) from exc

    values = dict(zip(field_names, key_parts))
    role = values.get("person_role")
    person_id = values.get("person_id")
    if role is not None and role not in ("study", "common_anchor"):
        raise SimulationConditionValidationError(
            "person_role must be study or common_anchor"
        )
    if person_id is not None:
        if not _is_positive_canonical_id(person_id, _PERSON_ID_PATTERN):
            raise SimulationConditionValidationError(
                "person_id must use canonical P/CA numbering"
            )
        if role == "study" and not person_id.startswith("P"):
            raise SimulationConditionValidationError(
                "study person_id must use the P prefix"
            )
        if role == "common_anchor" and not person_id.startswith("CA"):
            raise SimulationConditionValidationError(
                "common-anchor person_id must use the CA prefix"
            )
    artifact_id = values.get("artifact_id")
    if artifact_id is not None:
        prefix = f"{person_id}-A"
        if not artifact_id.startswith(prefix):
            raise SimulationConditionValidationError(
                "artifact_id must be nested under person_id"
            )
        suffix = artifact_id[len(prefix):]
        if (
            len(suffix) != 3
            or not suffix.isascii()
            or not suffix.isdigit()
            or not any(character != "0" for character in suffix)
        ):
            raise SimulationConditionValidationError(
                "artifact_id must use canonical A numbering"
            )
        if role == "common_anchor" and suffix != "001":
            raise SimulationConditionValidationError(
                "a common-anchor Person must use its A001 artifact"
            )
    rater_id = values.get("rater_id")
    if rater_id is not None and not _is_positive_canonical_id(
        rater_id, _RATER_ID_PATTERN
    ):
        raise SimulationConditionValidationError(
            "rater_id must use canonical R numbering"
        )
    criterion_id = values.get("criterion_id")
    if (
        criterion_id is not None
        and not _is_positive_canonical_id(criterion_id, _CRITERION_ID_PATTERN)
    ):
        raise SimulationConditionValidationError(
            "criterion_id must use canonical C numbering"
        )


def _length_prefixed(parts: tuple[bytes, ...]) -> bytes:
    encoded = bytearray()
    for part in parts:
        encoded.extend(len(part).to_bytes(8, byteorder="big", signed=False))
        encoded.extend(part)
    return bytes(encoded)


def keyed_uniform01(
    randomization: RandomizationSpecV1 | Mapping,
    *,
    replicate_index: int,
    stream: str,
    key_parts: tuple[str, ...],
    lane: int = 0,
) -> float:
    """Return a deterministic open-interval U(0, 1) variate.

    Inputs are encoded as 8-byte big-endian length-prefixed UTF-8 fields, then
    hashed with SHA-256.  The first 53 digest bits map to
    ``(x + 1) / (2**53 + 1)``.  Both endpoints therefore remain excluded even
    after binary64 rounding.
    No row order, batch size, worker count, or requested replicate total enters
    the key.
    """
    spec = normalize_randomization_spec(randomization)
    rep = _require_int(
        "replicate_index",
        replicate_index,
        minimum=1,
        maximum=MAX_MONTE_CARLO_REPLICATES,
    )
    lane_value = _require_int("lane", lane, minimum=0, maximum=2**31 - 1)
    if stream not in RANDOM_STREAM_CHOICES:
        raise SimulationConditionValidationError(
            f"stream must be one of {RANDOM_STREAM_CHOICES}"
        )
    if stream == RANDOM_STREAM_RESPONSE and lane_value != 0:
        raise SimulationConditionValidationError(
            "response_uniform uses lane 0 only in the v1 score-record contract"
        )
    if not isinstance(key_parts, tuple):
        raise SimulationConditionValidationError("key_parts must be a tuple")
    _validate_stream_key_parts(stream, key_parts)
    parts = (
        spec.algorithm.encode("utf-8"),
        str(spec.master_seed).encode("ascii"),
        str(rep).encode("ascii"),
        stream.encode("utf-8"),
        str(lane_value).encode("ascii"),
        *(part.encode("utf-8") for part in key_parts),
    )
    digest = hashlib.sha256(_length_prefixed(parts)).digest()
    word53 = int.from_bytes(digest[:8], byteorder="big", signed=False) >> 11
    return (word53 + 1) / (2**53 + 1)


def keyed_standard_normal(
    randomization: RandomizationSpecV1 | Mapping,
    *,
    replicate_index: int,
    stream: str,
    key_parts: tuple[str, ...],
) -> float:
    """Return one deterministic standard-normal draw via Box-Muller."""
    if stream == RANDOM_STREAM_RESPONSE:
        raise SimulationConditionValidationError(
            "response_uniform is a uniform-only stream"
        )
    first = keyed_uniform01(
        randomization,
        replicate_index=replicate_index,
        stream=stream,
        key_parts=key_parts,
        lane=0,
    )
    second = keyed_uniform01(
        randomization,
        replicate_index=replicate_index,
        stream=stream,
        key_parts=key_parts,
        lane=1,
    )
    return math.sqrt(-2.0 * math.log(first)) * math.cos(2.0 * math.pi * second)


__all__ = [
    "RANDOMIZATION_ALGORITHM_V1",
    "RANDOMIZATION_SHARED_RECORD_SCOPE_TWO_ARM",
    "RANDOMIZATION_SPEC_VERSION",
    "RANDOM_STREAM_CHOICES",
    "RANDOM_STREAM_KEY_SCHEMAS",
    "RANDOM_STREAM_CRITERION_DIFFICULTY",
    "RANDOM_STREAM_PERSON_THETA",
    "RANDOM_STREAM_RATER_CENTRAL_TENDENCY",
    "RANDOM_STREAM_RATER_PERSON_LOCAL_BIAS",
    "RANDOM_STREAM_RATER_SEVERITY",
    "RANDOM_STREAM_RESPONSE",
    "RandomizationSpecV1",
    "keyed_standard_normal",
    "keyed_uniform01",
    "normalize_randomization_spec",
    "randomization_spec_fingerprint",
]
