"""Deterministic rating-session schedules for prospective MFRM designs."""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import asdict, dataclass, fields
import hashlib
import json
import re

from .cost import DesignWorkloadV1, compute_design_workload
from .schema import (
    ASSIGNMENT_ALGORITHM_V1,
    DISCONNECTED_GROUPS,
    FULLY_CROSSED,
    NESTED_BRIDGE,
    SPIRAL,
    DesignSpecV1,
    design_spec_fingerprint,
    design_spec_to_dict,
    normalize_design_spec,
)


RATING_ASSIGNMENT_VERSION = "mfrm_rating_assignment_v1"
ASSIGNMENT_BUNDLE_VERSION = "mfrm_assignment_bundle_v1"
DEFAULT_MAX_RATING_SESSIONS = 100_000
SEED_NAMESPACE_PERSON = "person"
SEED_NAMESPACE_RATER = "rater"
SEED_NAMESPACE_BRIDGE_ARTIFACT = "bridge_artifact"
SEED_NAMESPACE_CHOICES = (
    SEED_NAMESPACE_PERSON,
    SEED_NAMESPACE_RATER,
    SEED_NAMESPACE_BRIDGE_ARTIFACT,
)

SESSION_PRIMARY = "primary"
SESSION_BRIDGE = "bridge"
SESSION_COMMON_ANCHOR = "common_anchor"
SESSION_TYPE_CHOICES = (
    SESSION_PRIMARY,
    SESSION_BRIDGE,
    SESSION_COMMON_ANCHOR,
)

PERSON_ROLE_STUDY = "study"
PERSON_ROLE_COMMON_ANCHOR = "common_anchor"
PERSON_ROLE_CHOICES = (PERSON_ROLE_STUDY, PERSON_ROLE_COMMON_ANCHOR)

_SESSION_SORT_ORDER = {
    SESSION_PRIMARY: 0,
    SESSION_BRIDGE: 1,
    SESSION_COMMON_ANCHOR: 2,
}
_STUDY_PERSON_PATTERN = re.compile(r"P[0-9]{6,}")
_ANCHOR_PERSON_PATTERN = re.compile(r"CA[0-9]{6,}")
_RATER_PATTERN = re.compile(r"R[0-9]{6,}")
_GROUP_PATTERN = re.compile(r"G[0-9]{3,}")


class AssignmentValidationError(ValueError):
    """Raised when an assignment schedule is invalid or internally inconsistent."""


class AssignmentLimitError(AssignmentValidationError):
    """Raised before materialization when a schedule exceeds its explicit cap."""

    def __init__(self, required_rating_sessions: int, max_rating_sessions: int):
        self.required_rating_sessions = required_rating_sessions
        self.max_rating_sessions = max_rating_sessions
        super().__init__(
            "design requires "
            f"{required_rating_sessions:,} rating sessions, exceeding the "
            f"explicit materialization limit of {max_rating_sessions:,}"
        )


def _has_positive_suffix(value: object, prefix: str) -> bool:
    """Return whether a canonical numeric identifier is one-based."""
    if not isinstance(value, str) or not value.startswith(prefix):
        return False
    suffix = value[len(prefix) :]
    return bool(suffix) and any(character != "0" for character in suffix)


def _canonical_ids(prefix: str, count: int, *, minimum_width: int) -> tuple[str, ...]:
    width = max(minimum_width, len(str(count)))
    return tuple(f"{prefix}{index:0{width}d}" for index in range(1, count + 1))


def _study_person_ids(spec: DesignSpecV1) -> tuple[str, ...]:
    return _canonical_ids("P", spec.n_persons, minimum_width=6)


def _rater_ids(spec: DesignSpecV1) -> tuple[str, ...]:
    return _canonical_ids("R", spec.n_raters, minimum_width=6)


def _group_ids(spec: DesignSpecV1) -> tuple[str, ...]:
    return _canonical_ids("G", spec.n_groups, minimum_width=3)


def _anchor_person_ids(spec: DesignSpecV1) -> tuple[str, ...]:
    return _canonical_ids(
        "CA",
        spec.n_common_anchor_artifacts,
        minimum_width=6,
    )


def _artifact_ids(person_id: str, count: int) -> tuple[str, ...]:
    width = max(3, len(str(count)))
    return tuple(
        f"{person_id}-A{index:0{width}d}" for index in range(1, count + 1)
    )


def _primary_artifacts(spec: DesignSpecV1) -> tuple[tuple[str, str], ...]:
    return tuple(
        (person_id, artifact_id)
        for person_id in _study_person_ids(spec)
        for artifact_id in _artifact_ids(person_id, spec.n_artifacts_per_person)
    )


def _seeded_order(
    identifiers: Sequence[str],
    *,
    seed: int,
    namespace: str,
) -> tuple[str, ...]:
    def key(identifier: str) -> tuple[bytes, str]:
        payload = (
            f"{ASSIGNMENT_ALGORITHM_V1}|{seed}|{namespace}|{identifier}"
        ).encode("utf-8")
        return hashlib.sha256(payload).digest(), identifier

    return tuple(sorted(identifiers, key=key))


def _balanced_partitions(
    values: Sequence[str],
    n_groups: int,
) -> tuple[tuple[str, ...], ...]:
    quotient, remainder = divmod(len(values), n_groups)
    groups: list[tuple[str, ...]] = []
    start = 0
    for group_index in range(n_groups):
        size = quotient + (1 if group_index < remainder else 0)
        groups.append(tuple(values[start : start + size]))
        start += size
    return tuple(groups)


@dataclass(frozen=True, slots=True, kw_only=True)
class RatingAssignmentV1:
    """One artifact-by-rater scoring session; never one long-format score row."""

    person_id: str
    artifact_id: str
    rater_id: str
    group_id: str
    person_role: str
    session_type: str

    def __post_init__(self) -> None:
        if self.person_role not in PERSON_ROLE_CHOICES:
            raise AssignmentValidationError(
                f"Unknown person_role: {self.person_role!r}"
            )
        if self.session_type not in SESSION_TYPE_CHOICES:
            raise AssignmentValidationError(
                f"Unknown session_type: {self.session_type!r}"
            )
        if (
            not isinstance(self.rater_id, str)
            or _RATER_PATTERN.fullmatch(self.rater_id) is None
            or not _has_positive_suffix(self.rater_id, "R")
        ):
            raise AssignmentValidationError(f"Invalid rater_id: {self.rater_id!r}")
        if (
            not isinstance(self.group_id, str)
            or _GROUP_PATTERN.fullmatch(self.group_id) is None
            or not _has_positive_suffix(self.group_id, "G")
        ):
            raise AssignmentValidationError(f"Invalid group_id: {self.group_id!r}")

        if self.person_role == PERSON_ROLE_STUDY:
            if (
                not isinstance(self.person_id, str)
                or _STUDY_PERSON_PATTERN.fullmatch(self.person_id) is None
                or not _has_positive_suffix(self.person_id, "P")
            ):
                raise AssignmentValidationError(
                    f"Invalid study person_id: {self.person_id!r}"
                )
            artifact_pattern = re.compile(
                rf"{re.escape(self.person_id)}-A[0-9]{{3,}}"
            )
            if (
                not isinstance(self.artifact_id, str)
                or artifact_pattern.fullmatch(self.artifact_id) is None
                or not _has_positive_suffix(
                    self.artifact_id,
                    f"{self.person_id}-A",
                )
            ):
                raise AssignmentValidationError(
                    "study artifact_id must be nested under its person_id"
                )
            if self.session_type == SESSION_COMMON_ANCHOR:
                raise AssignmentValidationError(
                    "study persons cannot have common_anchor sessions"
                )
        else:
            if (
                not isinstance(self.person_id, str)
                or _ANCHOR_PERSON_PATTERN.fullmatch(self.person_id) is None
                or not _has_positive_suffix(self.person_id, "CA")
            ):
                raise AssignmentValidationError(
                    f"Invalid common-anchor person_id: {self.person_id!r}"
                )
            if self.artifact_id != f"{self.person_id}-A001":
                raise AssignmentValidationError(
                    "a v1 common-anchor person must own exactly its A001 artifact"
                )
            if self.session_type != SESSION_COMMON_ANCHOR:
                raise AssignmentValidationError(
                    "common-anchor persons require common_anchor sessions"
                )

    @property
    def session_key(self) -> tuple[str, str, str]:
        return self.person_id, self.artifact_id, self.rater_id

    def to_dict(self) -> dict:
        payload = asdict(self)
        json.dumps(payload, allow_nan=False)
        return payload


def _assignment_sort_key(row: RatingAssignmentV1) -> tuple:
    return (
        row.person_id,
        row.artifact_id,
        _SESSION_SORT_ORDER[row.session_type],
        row.rater_id,
        row.group_id,
    )


def assignment_fingerprint(
    assignments: Iterable[RatingAssignmentV1],
    *,
    length: int = 16,
) -> str:
    """Hash a canonical assignment multiset independently of input row order."""
    if (
        isinstance(length, bool)
        or not isinstance(length, int)
        or not 8 <= length <= 64
    ):
        raise ValueError("fingerprint length must be an integer between 8 and 64")
    rows = tuple(assignments)
    if any(not isinstance(row, RatingAssignmentV1) for row in rows):
        raise TypeError("assignments must contain RatingAssignmentV1 rows")
    digest = hashlib.sha256()
    digest.update(f"{RATING_ASSIGNMENT_VERSION}\n".encode("utf-8"))
    for row in sorted(rows, key=_assignment_sort_key):
        encoded = json.dumps(
            row.to_dict(),
            allow_nan=False,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")
        digest.update(encoded)
        digest.update(b"\n")
    return digest.hexdigest()[:length]


@dataclass(frozen=True, slots=True, kw_only=True)
class AssignmentBundleV1:
    """A validated design, workload, and canonical rating-session schedule."""

    schema_version: str
    design_spec: DesignSpecV1
    design_spec_fingerprint: str
    assignment_algorithm: str
    assignment_fingerprint: str
    workload: DesignWorkloadV1
    assignments: tuple[RatingAssignmentV1, ...]

    def __post_init__(self) -> None:
        if self.schema_version != ASSIGNMENT_BUNDLE_VERSION:
            raise AssignmentValidationError(
                f"Unsupported assignment bundle version: {self.schema_version!r}"
            )
        if type(self.design_spec) is not DesignSpecV1:
            raise AssignmentValidationError("design_spec must be a DesignSpecV1")
        expected_design_fingerprint = design_spec_fingerprint(self.design_spec)
        if self.design_spec_fingerprint != expected_design_fingerprint:
            raise AssignmentValidationError(
                "design_spec_fingerprint does not match design_spec"
            )
        if self.assignment_algorithm != self.design_spec.assignment_algorithm:
            raise AssignmentValidationError(
                "assignment_algorithm does not match design_spec"
            )
        if type(self.workload) is not DesignWorkloadV1:
            raise AssignmentValidationError("workload must be a DesignWorkloadV1")
        if self.workload != compute_design_workload(self.design_spec):
            raise AssignmentValidationError("workload does not match design_spec")
        if not isinstance(self.assignments, tuple) or any(
            type(row) is not RatingAssignmentV1 for row in self.assignments
        ):
            raise AssignmentValidationError(
                "assignments must be a tuple of RatingAssignmentV1 rows"
            )
        if tuple(sorted(self.assignments, key=_assignment_sort_key)) != self.assignments:
            raise AssignmentValidationError("assignments must use canonical row order")
        keys = tuple(row.session_key for row in self.assignments)
        if len(keys) != len(set(keys)):
            raise AssignmentValidationError(
                "assignments contain duplicate (Person, Artifact, Rater) sessions"
            )
        if len(self.assignments) != self.workload.total_rating_sessions:
            raise AssignmentValidationError(
                "assignment row count does not match total_rating_sessions"
            )
        expected_types = {
            SESSION_PRIMARY: self.workload.primary_rating_sessions,
            SESSION_BRIDGE: self.workload.bridge_rating_sessions,
            SESSION_COMMON_ANCHOR: self.workload.common_anchor_rating_sessions,
        }
        actual_types = Counter(row.session_type for row in self.assignments)
        if any(actual_types[name] != count for name, count in expected_types.items()):
            raise AssignmentValidationError(
                "assignment session-type counts do not match workload"
            )
        expected_raters = set(_rater_ids(self.design_spec))
        actual_raters = {row.rater_id for row in self.assignments}
        if actual_raters != expected_raters:
            raise AssignmentValidationError(
                "every declared active rater must have at least one session"
            )
        expected_people = set(_study_person_ids(self.design_spec))
        actual_people = {
            row.person_id
            for row in self.assignments
            if row.person_role == PERSON_ROLE_STUDY
        }
        if actual_people != expected_people:
            raise AssignmentValidationError(
                "assignment study-person IDs do not match design_spec"
            )
        expected_artifacts = set(_primary_artifacts(self.design_spec))
        primary_counts = Counter(
            (row.person_id, row.artifact_id)
            for row in self.assignments
            if row.session_type == SESSION_PRIMARY
        )
        if set(primary_counts) != expected_artifacts or any(
            count != self.design_spec.raters_per_artifact
            for count in primary_counts.values()
        ):
            raise AssignmentValidationError(
                "primary assignment counts do not match declared study artifacts"
            )
        expected_anchors = set(_anchor_person_ids(self.design_spec))
        actual_anchors = {
            row.person_id
            for row in self.assignments
            if row.person_role == PERSON_ROLE_COMMON_ANCHOR
        }
        if actual_anchors != expected_anchors:
            raise AssignmentValidationError(
                "assignment common-anchor IDs do not match design_spec"
            )
        expected_groups = set(_group_ids(self.design_spec))
        if {row.group_id for row in self.assignments} != expected_groups:
            raise AssignmentValidationError(
                "assignment group IDs do not match design_spec"
            )
        if self.design_spec.plan == NESTED_BRIDGE:
            primary_raters_by_person: dict[str, set[str]] = {
                person_id: set() for person_id in expected_people
            }
            for row in self.assignments:
                if row.session_type == SESSION_PRIMARY:
                    primary_raters_by_person[row.person_id].add(row.rater_id)
            if any(len(raters) != 1 for raters in primary_raters_by_person.values()):
                raise AssignmentValidationError(
                    "nested_bridge must keep each person's primary artifacts nested"
                )
            bridge_counts = Counter(
                (row.person_id, row.artifact_id)
                for row in self.assignments
                if row.session_type == SESSION_BRIDGE
            )
            if len(bridge_counts) != self.design_spec.n_bridge_artifacts or any(
                count != self.design_spec.bridge_extra_raters_per_artifact
                for count in bridge_counts.values()
            ):
                raise AssignmentValidationError(
                    "bridge target/count structure does not match design_spec"
                )
            if not set(bridge_counts).issubset(expected_artifacts):
                raise AssignmentValidationError(
                    "bridge targets must be declared primary artifacts"
                )
        if self.design_spec.plan == DISCONNECTED_GROUPS:
            ordered_people = _seeded_order(
                _study_person_ids(self.design_spec),
                seed=self.design_spec.assignment_seed,
                namespace=SEED_NAMESPACE_PERSON,
            )
            ordered_raters = _seeded_order(
                _rater_ids(self.design_spec),
                seed=self.design_spec.assignment_seed,
                namespace=SEED_NAMESPACE_RATER,
            )
            person_groups = _balanced_partitions(
                ordered_people,
                self.design_spec.n_groups,
            )
            rater_groups = _balanced_partitions(
                ordered_raters,
                self.design_spec.n_groups,
            )
            person_group = {
                person_id: group_id
                for group_id, people in zip(
                    _group_ids(self.design_spec),
                    person_groups,
                    strict=True,
                )
                for person_id in people
            }
            rater_group = {
                rater_id: group_id
                for group_id, raters in zip(
                    _group_ids(self.design_spec),
                    rater_groups,
                    strict=True,
                )
                for rater_id in raters
            }
            primary_raters_by_person: dict[str, set[str]] = {
                person_id: set() for person_id in expected_people
            }
            for row in self.assignments:
                if row.session_type == SESSION_PRIMARY:
                    if (
                        person_group[row.person_id] != row.group_id
                        or rater_group[row.rater_id] != row.group_id
                    ):
                        raise AssignmentValidationError(
                            "disconnected primary sessions cannot cross groups"
                        )
                    primary_raters_by_person[row.person_id].add(row.rater_id)
                elif row.session_type == SESSION_COMMON_ANCHOR:
                    if rater_group[row.rater_id] != row.group_id:
                        raise AssignmentValidationError(
                            "common-anchor rater does not belong to row group"
                        )
            if any(len(raters) != 1 for raters in primary_raters_by_person.values()):
                raise AssignmentValidationError(
                    "disconnected_groups must nest each person's primary artifacts"
                )
            expected_anchor_groups = {
                (anchor_id, group_id)
                for anchor_id in expected_anchors
                for group_id in expected_groups
            }
            anchor_group_counts = Counter(
                (row.person_id, row.group_id)
                for row in self.assignments
                if row.session_type == SESSION_COMMON_ANCHOR
            )
            if set(anchor_group_counts) != expected_anchor_groups or any(
                count != self.design_spec.common_anchor_raters_per_group
                for count in anchor_group_counts.values()
            ):
                raise AssignmentValidationError(
                    "common-anchor group/count structure does not match design_spec"
                )
        expected_assignment_fingerprint = assignment_fingerprint(self.assignments)
        if self.assignment_fingerprint != expected_assignment_fingerprint:
            raise AssignmentValidationError(
                "assignment_fingerprint does not match assignment rows"
            )
        if self.assignments != _compile_assignment_rows(self.design_spec):
            raise AssignmentValidationError(
                "assignment rows do not match the declared assignment algorithm"
            )

    @property
    def cache_identity(self) -> tuple[str, str, str]:
        """Identity for caches that depend on both design meaning and exact rows."""
        return (
            self.schema_version,
            self.design_spec_fingerprint,
            self.assignment_fingerprint,
        )

    def to_dict(self) -> dict:
        payload = {
            "schema_version": self.schema_version,
            "design_spec": design_spec_to_dict(self.design_spec),
            "design_spec_fingerprint": self.design_spec_fingerprint,
            "assignment_algorithm": self.assignment_algorithm,
            "assignment_fingerprint": self.assignment_fingerprint,
            "workload": self.workload.to_dict(),
            "assignments": [row.to_dict() for row in self.assignments],
        }
        json.dumps(payload, allow_nan=False)
        return payload


def _choose_extra_raters(
    *,
    ordered_raters: Sequence[str],
    active_raters: set[str],
    already_on_artifact: set[str],
    count: int,
    cursor: int,
    inactive_cursor: int,
) -> tuple[tuple[str, ...], int, int]:
    chosen: list[str] = []
    n_raters = len(ordered_raters)
    for _ in range(count):
        selected_index: int | None = None
        while (
            inactive_cursor < n_raters
            and ordered_raters[inactive_cursor] in active_raters
        ):
            inactive_cursor += 1
        if inactive_cursor < n_raters:
            selected_index = inactive_cursor
            inactive_cursor += 1
        else:
            for offset in range(n_raters):
                candidate_index = (cursor + offset) % n_raters
                candidate = ordered_raters[candidate_index]
                if candidate not in already_on_artifact:
                    selected_index = candidate_index
                    break
        if selected_index is None:
            raise AssignmentValidationError(
                "assignment algorithm could not find an eligible extra rater"
            )
        selected = ordered_raters[selected_index]
        chosen.append(selected)
        already_on_artifact.add(selected)
        active_raters.add(selected)
        cursor = (selected_index + 1) % n_raters
    return tuple(chosen), cursor, inactive_cursor


def _compile_fully_crossed(spec: DesignSpecV1) -> list[RatingAssignmentV1]:
    group_id = _group_ids(spec)[0]
    return [
        RatingAssignmentV1(
            person_id=person_id,
            artifact_id=artifact_id,
            rater_id=rater_id,
            group_id=group_id,
            person_role=PERSON_ROLE_STUDY,
            session_type=SESSION_PRIMARY,
        )
        for person_id, artifact_id in _primary_artifacts(spec)
        for rater_id in _rater_ids(spec)
    ]


def _compile_spiral(spec: DesignSpecV1) -> list[RatingAssignmentV1]:
    group_id = _group_ids(spec)[0]
    ordered_raters = _seeded_order(
        _rater_ids(spec),
        seed=spec.assignment_seed,
        namespace=SEED_NAMESPACE_RATER,
    )
    rows: list[RatingAssignmentV1] = []
    for artifact_index, (person_id, artifact_id) in enumerate(
        _primary_artifacts(spec)
    ):
        for slot in range(spec.raters_per_artifact):
            rater_index = (
                artifact_index * spec.raters_per_artifact + slot
            ) % spec.n_raters
            rows.append(
                RatingAssignmentV1(
                    person_id=person_id,
                    artifact_id=artifact_id,
                    rater_id=ordered_raters[rater_index],
                    group_id=group_id,
                    person_role=PERSON_ROLE_STUDY,
                    session_type=SESSION_PRIMARY,
                )
            )
    return rows


def _compile_nested_bridge(spec: DesignSpecV1) -> list[RatingAssignmentV1]:
    group_id = _group_ids(spec)[0]
    people = _study_person_ids(spec)
    ordered_people = _seeded_order(
        people,
        seed=spec.assignment_seed,
        namespace=SEED_NAMESPACE_PERSON,
    )
    ordered_raters = _seeded_order(
        _rater_ids(spec),
        seed=spec.assignment_seed,
        namespace=SEED_NAMESPACE_RATER,
    )
    primary_rater = {
        person_id: ordered_raters[index % spec.n_raters]
        for index, person_id in enumerate(ordered_people)
    }
    rows = [
        RatingAssignmentV1(
            person_id=person_id,
            artifact_id=artifact_id,
            rater_id=primary_rater[person_id],
            group_id=group_id,
            person_role=PERSON_ROLE_STUDY,
            session_type=SESSION_PRIMARY,
        )
        for person_id, artifact_id in _primary_artifacts(spec)
    ]
    active_raters = set(primary_rater.values())
    bridge_candidates = _seeded_order(
        tuple(artifact_id for _, artifact_id in _primary_artifacts(spec)),
        seed=spec.assignment_seed,
        namespace=SEED_NAMESPACE_BRIDGE_ARTIFACT,
    )[: spec.n_bridge_artifacts]
    artifact_to_person = {
        artifact_id: person_id
        for person_id, artifact_id in _primary_artifacts(spec)
    }
    cursor = 0
    inactive_cursor = 0
    for artifact_id in bridge_candidates:
        person_id = artifact_to_person[artifact_id]
        already_on_artifact = {primary_rater[person_id]}
        selected, cursor, inactive_cursor = _choose_extra_raters(
            ordered_raters=ordered_raters,
            active_raters=active_raters,
            already_on_artifact=already_on_artifact,
            count=spec.bridge_extra_raters_per_artifact,
            cursor=cursor,
            inactive_cursor=inactive_cursor,
        )
        rows.extend(
            RatingAssignmentV1(
                person_id=person_id,
                artifact_id=artifact_id,
                rater_id=rater_id,
                group_id=group_id,
                person_role=PERSON_ROLE_STUDY,
                session_type=SESSION_BRIDGE,
            )
            for rater_id in selected
        )
    return rows


def _compile_disconnected_groups(spec: DesignSpecV1) -> list[RatingAssignmentV1]:
    group_ids = _group_ids(spec)
    ordered_people = _seeded_order(
        _study_person_ids(spec),
        seed=spec.assignment_seed,
        namespace=SEED_NAMESPACE_PERSON,
    )
    ordered_raters = _seeded_order(
        _rater_ids(spec),
        seed=spec.assignment_seed,
        namespace=SEED_NAMESPACE_RATER,
    )
    person_groups = _balanced_partitions(ordered_people, spec.n_groups)
    rater_groups = _balanced_partitions(ordered_raters, spec.n_groups)
    rows: list[RatingAssignmentV1] = []
    active_by_group: dict[str, set[str]] = {}
    cursor_by_group: dict[str, int] = {}
    inactive_cursor_by_group: dict[str, int] = {}

    for group_id, people, raters in zip(
        group_ids,
        person_groups,
        rater_groups,
        strict=True,
    ):
        person_to_rater = {
            person_id: raters[index % len(raters)]
            for index, person_id in enumerate(people)
        }
        active_by_group[group_id] = set(person_to_rater.values())
        cursor_by_group[group_id] = 0
        inactive_cursor_by_group[group_id] = 0
        for person_id in people:
            rows.extend(
                RatingAssignmentV1(
                    person_id=person_id,
                    artifact_id=artifact_id,
                    rater_id=person_to_rater[person_id],
                    group_id=group_id,
                    person_role=PERSON_ROLE_STUDY,
                    session_type=SESSION_PRIMARY,
                )
                for artifact_id in _artifact_ids(
                    person_id,
                    spec.n_artifacts_per_person,
                )
            )

    for anchor_person_id in _anchor_person_ids(spec):
        artifact_id = f"{anchor_person_id}-A001"
        for group_id, raters in zip(group_ids, rater_groups, strict=True):
            (
                selected,
                cursor_by_group[group_id],
                inactive_cursor_by_group[group_id],
            ) = _choose_extra_raters(
                ordered_raters=raters,
                active_raters=active_by_group[group_id],
                already_on_artifact=set(),
                count=spec.common_anchor_raters_per_group,
                cursor=cursor_by_group[group_id],
                inactive_cursor=inactive_cursor_by_group[group_id],
            )
            rows.extend(
                RatingAssignmentV1(
                    person_id=anchor_person_id,
                    artifact_id=artifact_id,
                    rater_id=rater_id,
                    group_id=group_id,
                    person_role=PERSON_ROLE_COMMON_ANCHOR,
                    session_type=SESSION_COMMON_ANCHOR,
                )
                for rater_id in selected
            )
    return rows


def _compile_assignment_rows(
    spec: DesignSpecV1,
) -> tuple[RatingAssignmentV1, ...]:
    """Materialize the one canonical row set for the declared algorithm."""
    if spec.plan == FULLY_CROSSED:
        rows = _compile_fully_crossed(spec)
    elif spec.plan == SPIRAL:
        rows = _compile_spiral(spec)
    elif spec.plan == NESTED_BRIDGE:
        rows = _compile_nested_bridge(spec)
    elif spec.plan == DISCONNECTED_GROUPS:
        rows = _compile_disconnected_groups(spec)
    else:  # pragma: no cover - DesignSpecV1 already rejects unknown plans.
        raise AssignmentValidationError(f"Unsupported design plan: {spec.plan!r}")
    return tuple(sorted(rows, key=_assignment_sort_key))


def _validate_materialization_limit(value: object) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError("max_rating_sessions must be a positive integer")
    return value


def compile_rating_design(
    value: DesignSpecV1 | Mapping,
    *,
    max_rating_sessions: int = DEFAULT_MAX_RATING_SESSIONS,
) -> AssignmentBundleV1:
    """Compile one deterministic schedule after a formula-only size preflight."""
    limit = _validate_materialization_limit(max_rating_sessions)
    spec = normalize_design_spec(value)
    required_sessions = (
        spec.n_persons
        * spec.n_artifacts_per_person
        * spec.raters_per_artifact
        + spec.n_bridge_artifacts * spec.bridge_extra_raters_per_artifact
        + spec.n_common_anchor_artifacts
        * spec.n_groups
        * spec.common_anchor_raters_per_group
    )
    if required_sessions > limit:
        raise AssignmentLimitError(required_sessions, limit)
    workload = compute_design_workload(spec)
    if workload.total_rating_sessions != required_sessions:
        raise AssignmentValidationError(
            "formula-only preflight disagrees with the versioned workload"
        )
    assignments = _compile_assignment_rows(spec)
    return AssignmentBundleV1(
        schema_version=ASSIGNMENT_BUNDLE_VERSION,
        design_spec=spec,
        design_spec_fingerprint=design_spec_fingerprint(spec),
        assignment_algorithm=spec.assignment_algorithm,
        assignment_fingerprint=assignment_fingerprint(assignments),
        workload=workload,
        assignments=assignments,
    )


_ROW_FIELDS = tuple(field.name for field in fields(RatingAssignmentV1))
_ROW_FIELD_SET = frozenset(_ROW_FIELDS)
_WORKLOAD_FIELDS = tuple(field.name for field in fields(DesignWorkloadV1))
_WORKLOAD_FIELD_SET = frozenset(_WORKLOAD_FIELDS)
_BUNDLE_FIELDS = tuple(field.name for field in fields(AssignmentBundleV1))
_BUNDLE_FIELD_SET = frozenset(_BUNDLE_FIELDS)


def _require_exact_mapping(
    value: object,
    *,
    field_set: frozenset[str],
    label: str,
) -> Mapping:
    if not isinstance(value, Mapping):
        raise AssignmentValidationError(f"{label} must be a mapping")
    supplied = set(value)
    if any(not isinstance(name, str) for name in supplied):
        raise AssignmentValidationError(f"{label} field names must be strings")
    missing = field_set - supplied
    unknown = supplied - field_set
    if missing:
        raise AssignmentValidationError(
            f"{label} is missing fields: {sorted(missing)}"
        )
    if unknown:
        raise AssignmentValidationError(
            f"{label} contains unknown fields: {sorted(unknown)}"
        )
    return value


def normalize_assignment_bundle(
    value: AssignmentBundleV1 | Mapping,
    *,
    max_rating_sessions: int = DEFAULT_MAX_RATING_SESSIONS,
) -> AssignmentBundleV1:
    """Strictly parse and revalidate a saved assignment bundle."""
    limit = _validate_materialization_limit(max_rating_sessions)
    if isinstance(value, AssignmentBundleV1):
        if type(value) is not AssignmentBundleV1:
            raise AssignmentValidationError(
                "AssignmentBundleV1 subclasses are not accepted as validated bundles"
            )
        if len(value.assignments) > limit:
            raise AssignmentLimitError(len(value.assignments), limit)
        return value
    bundle = _require_exact_mapping(
        value,
        field_set=_BUNDLE_FIELD_SET,
        label="Saved assignment bundle",
    )
    workload_mapping = _require_exact_mapping(
        bundle["workload"],
        field_set=_WORKLOAD_FIELD_SET,
        label="Saved assignment workload",
    )
    assignments_value = bundle["assignments"]
    if not isinstance(assignments_value, (list, tuple)):
        raise AssignmentValidationError(
            "Saved assignment rows must be a list or tuple"
        )
    if len(assignments_value) > limit:
        raise AssignmentLimitError(len(assignments_value), limit)
    rows = []
    for index, value_row in enumerate(assignments_value):
        row = _require_exact_mapping(
            value_row,
            field_set=_ROW_FIELD_SET,
            label=f"Saved assignment row {index}",
        )
        rows.append(RatingAssignmentV1(**{name: row[name] for name in _ROW_FIELDS}))
    try:
        workload = DesignWorkloadV1(
            **{name: workload_mapping[name] for name in _WORKLOAD_FIELDS}
        )
        spec = normalize_design_spec(bundle["design_spec"])
    except (TypeError, ValueError) as exc:
        raise AssignmentValidationError(str(exc)) from exc
    return AssignmentBundleV1(
        schema_version=bundle["schema_version"],
        design_spec=spec,
        design_spec_fingerprint=bundle["design_spec_fingerprint"],
        assignment_algorithm=bundle["assignment_algorithm"],
        assignment_fingerprint=bundle["assignment_fingerprint"],
        workload=workload,
        assignments=tuple(rows),
    )


def assignment_bundle_to_dict(
    value: AssignmentBundleV1 | Mapping,
    *,
    max_rating_sessions: int = DEFAULT_MAX_RATING_SESSIONS,
) -> dict:
    """Return a validated JSON-primitive assignment bundle."""
    return normalize_assignment_bundle(
        value,
        max_rating_sessions=max_rating_sessions,
    ).to_dict()


__all__ = [
    "ASSIGNMENT_BUNDLE_VERSION",
    "DEFAULT_MAX_RATING_SESSIONS",
    "PERSON_ROLE_CHOICES",
    "PERSON_ROLE_COMMON_ANCHOR",
    "PERSON_ROLE_STUDY",
    "RATING_ASSIGNMENT_VERSION",
    "SEED_NAMESPACE_BRIDGE_ARTIFACT",
    "SEED_NAMESPACE_CHOICES",
    "SEED_NAMESPACE_PERSON",
    "SEED_NAMESPACE_RATER",
    "SESSION_BRIDGE",
    "SESSION_COMMON_ANCHOR",
    "SESSION_PRIMARY",
    "SESSION_TYPE_CHOICES",
    "AssignmentBundleV1",
    "AssignmentLimitError",
    "AssignmentValidationError",
    "RatingAssignmentV1",
    "assignment_bundle_to_dict",
    "assignment_fingerprint",
    "compile_rating_design",
    "normalize_assignment_bundle",
]
