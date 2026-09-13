"""Deterministic Person-Rater topology and structural repair contracts.

This module audits a *compiled* prospective rating schedule.  It deliberately
stays narrower than the fitted-data network diagnostics in ``streamlit_app``:
only Person-Rater incidence is represented, and no score, estimator, standard
error, recovery, or classification claim is made.

Two graphs are always reported. ``primary_study`` contains study-person
primary and bridge sessions while retaining every declared rater (including
raters that are active only through common-anchor work) as a possible isolate.
``anchor_augmented`` then adds common-anchor persons and sessions.  Keeping the
two scopes separate prevents an anchor count or plan label from silently being
treated as evidence of an actually connected schedule.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
import hashlib
import json
import re

from .assignments import (
    ASSIGNMENT_BUNDLE_VERSION,
    DEFAULT_MAX_RATING_SESSIONS,
    PERSON_ROLE_COMMON_ANCHOR,
    PERSON_ROLE_STUDY,
    SESSION_PRIMARY,
    AssignmentBundleV1,
    RatingAssignmentV1,
    normalize_assignment_bundle,
)


TOPOLOGY_AUDIT_VERSION = "mfrm_rating_topology_audit_v1"
TOPOLOGY_SCOPE_VERSION = "mfrm_rating_topology_scope_v1"
TOPOLOGY_COMPONENT_VERSION = "mfrm_rating_topology_component_v1"
RATER_LOAD_VERSION = "mfrm_rating_topology_rater_load_v1"
TOPOLOGY_EDGE_VERSION = "mfrm_rating_topology_edge_v1"
STRUCTURAL_REPAIR_VERSION = "mfrm_structural_bridge_repair_v1"
STRUCTURAL_REPAIR_SESSION_VERSION = "mfrm_structural_bridge_session_v1"

PRIMARY_STUDY_SCOPE = "primary_study"
ANCHOR_AUGMENTED_SCOPE = "anchor_augmented"
TOPOLOGY_SCOPE_CHOICES = (PRIMARY_STUDY_SCOPE, ANCHOR_AUGMENTED_SCOPE)

REPAIR_NOT_NEEDED = "not_needed"
REPAIR_PROPOSED = "proposed"
REPAIR_STATUS_CHOICES = (REPAIR_NOT_NEEDED, REPAIR_PROPOSED)
PRECISION_NOT_EVALUATED = "not_evaluated"

_COMPONENT_ID_PATTERN = re.compile(r"C[0-9]{3,}")


class TopologyValidationError(ValueError):
    """Raised when a topology or repair payload contradicts its base design."""


def _require_nonnegative_int(name: str, value: object) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise TopologyValidationError(f"{name} must be a nonnegative integer")
    return value


def _require_positive_int(name: str, value: object) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise TopologyValidationError(f"{name} must be a positive integer")
    return value


def _require_lowercase_fingerprint(name: str, value: object) -> str:
    if (
        not isinstance(value, str)
        or len(value) != 16
        or value != value.lower()
        or any(character not in "0123456789abcdef" for character in value)
    ):
        raise TopologyValidationError(
            f"{name} must be a 16-character lowercase hexadecimal fingerprint"
        )
    return value


def _require_sorted_unique_strings(name: str, value: object) -> tuple[str, ...]:
    if not isinstance(value, tuple) or any(
        not isinstance(item, str) or not item for item in value
    ):
        raise TopologyValidationError(f"{name} must be a tuple of non-empty strings")
    if value != tuple(sorted(set(value))):
        raise TopologyValidationError(f"{name} must be sorted and unique")
    return value


def _json_payload(value: dict) -> dict:
    json.dumps(value, allow_nan=False, ensure_ascii=False)
    return value


@dataclass(frozen=True, slots=True, kw_only=True)
class RaterTopologyLoadV1:
    """Auditable rater workload inside one topology scope."""

    schema_version: str
    rater_id: str
    rating_sessions: int
    unique_persons: int
    study_persons: int
    common_anchor_persons: int

    def __post_init__(self) -> None:
        if self.schema_version != RATER_LOAD_VERSION:
            raise TopologyValidationError(
                f"Unsupported rater-load version: {self.schema_version!r}"
            )
        if not isinstance(self.rater_id, str) or not self.rater_id:
            raise TopologyValidationError("rater_id must be a non-empty string")
        for name in (
            "rating_sessions",
            "unique_persons",
            "study_persons",
            "common_anchor_persons",
        ):
            _require_nonnegative_int(name, getattr(self, name))
        if self.unique_persons != self.study_persons + self.common_anchor_persons:
            raise TopologyValidationError(
                "unique_persons must equal study_persons plus common_anchor_persons"
            )
        if self.rating_sessions < self.unique_persons:
            raise TopologyValidationError(
                "rating_sessions cannot be smaller than unique_persons"
            )

    def to_dict(self) -> dict:
        return _json_payload({
            "schema_version": self.schema_version,
            "rater_id": self.rater_id,
            "rating_sessions": self.rating_sessions,
            "unique_persons": self.unique_persons,
            "study_persons": self.study_persons,
            "common_anchor_persons": self.common_anchor_persons,
        })


@dataclass(frozen=True, slots=True, kw_only=True)
class TopologyIncidenceEdgeV1:
    """One unique Person-Rater edge and its artifact-session multiplicity."""

    schema_version: str
    component_id: str
    person_role: str
    person_id: str
    rater_id: str
    rating_sessions: int

    def __post_init__(self) -> None:
        if self.schema_version != TOPOLOGY_EDGE_VERSION:
            raise TopologyValidationError(
                f"Unsupported topology-edge version: {self.schema_version!r}"
            )
        if (
            not isinstance(self.component_id, str)
            or _COMPONENT_ID_PATTERN.fullmatch(self.component_id) is None
        ):
            raise TopologyValidationError(
                "component_id must use canonical C001-style numbering"
            )
        if self.person_role not in (
            PERSON_ROLE_STUDY,
            PERSON_ROLE_COMMON_ANCHOR,
        ):
            raise TopologyValidationError(
                f"Unknown edge person_role: {self.person_role!r}"
            )
        for name in ("person_id", "rater_id"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value:
                raise TopologyValidationError(f"{name} must be a non-empty string")
        _require_positive_int("rating_sessions", self.rating_sessions)

    @property
    def edge_key(self) -> tuple[str, str, str]:
        return self.person_role, self.person_id, self.rater_id

    def to_dict(self) -> dict:
        return _json_payload({
            "schema_version": self.schema_version,
            "component_id": self.component_id,
            "person_role": self.person_role,
            "person_id": self.person_id,
            "rater_id": self.rater_id,
            "rating_sessions": self.rating_sessions,
        })


@dataclass(frozen=True, slots=True, kw_only=True)
class TopologyComponentV1:
    """One deterministic connected component of a bipartite scope."""

    schema_version: str
    component_id: str
    study_person_ids: tuple[str, ...]
    common_anchor_person_ids: tuple[str, ...]
    rater_ids: tuple[str, ...]
    unique_incidence_edges: int
    rating_sessions: int

    def __post_init__(self) -> None:
        if self.schema_version != TOPOLOGY_COMPONENT_VERSION:
            raise TopologyValidationError(
                f"Unsupported topology-component version: {self.schema_version!r}"
            )
        if (
            not isinstance(self.component_id, str)
            or _COMPONENT_ID_PATTERN.fullmatch(self.component_id) is None
        ):
            raise TopologyValidationError(
                "component_id must use canonical C001-style numbering"
            )
        _require_sorted_unique_strings("study_person_ids", self.study_person_ids)
        _require_sorted_unique_strings(
            "common_anchor_person_ids", self.common_anchor_person_ids
        )
        _require_sorted_unique_strings("rater_ids", self.rater_ids)
        if not (
            self.study_person_ids
            or self.common_anchor_person_ids
            or self.rater_ids
        ):
            raise TopologyValidationError("a topology component cannot be empty")
        _require_nonnegative_int(
            "unique_incidence_edges", self.unique_incidence_edges
        )
        _require_nonnegative_int("rating_sessions", self.rating_sessions)
        if self.rating_sessions < self.unique_incidence_edges:
            raise TopologyValidationError(
                "rating_sessions cannot be smaller than unique_incidence_edges"
            )
        if not self.rater_ids:
            raise TopologyValidationError(
                "every Person-Rater topology component must contain a rater"
            )
        if self.unique_incidence_edges == 0 and (
            self.study_person_ids
            or self.common_anchor_person_ids
            or len(self.rater_ids) != 1
            or self.rating_sessions != 0
        ):
            raise TopologyValidationError(
                "a zero-edge component must be one isolated declared rater"
            )
        if self.unique_incidence_edges > 0 and not (
            self.study_person_ids or self.common_anchor_person_ids
        ):
            raise TopologyValidationError(
                "a nonempty incidence component must contain a person node"
            )

    @property
    def n_nodes(self) -> int:
        return (
            len(self.study_person_ids)
            + len(self.common_anchor_person_ids)
            + len(self.rater_ids)
        )

    def to_dict(self) -> dict:
        return _json_payload({
            "schema_version": self.schema_version,
            "component_id": self.component_id,
            "study_person_ids": list(self.study_person_ids),
            "common_anchor_person_ids": list(self.common_anchor_person_ids),
            "rater_ids": list(self.rater_ids),
            "unique_incidence_edges": self.unique_incidence_edges,
            "rating_sessions": self.rating_sessions,
        })


@dataclass(frozen=True, slots=True, kw_only=True)
class TopologyScopeAuditV1:
    """Summary, components, and rater loads for one graph scope."""

    schema_version: str
    scope: str
    declared_study_persons: int
    declared_common_anchor_persons: int
    declared_raters: int
    nodes: int
    unique_incidence_edges: int
    rating_sessions: int
    repeated_sessions_collapsed: int
    n_components: int
    largest_component_nodes: int
    largest_component_study_persons: int
    largest_component_raters: int
    isolated_raters: int
    connected: bool
    components: tuple[TopologyComponentV1, ...]
    incidence_edges: tuple[TopologyIncidenceEdgeV1, ...]
    rater_loads: tuple[RaterTopologyLoadV1, ...]

    def __post_init__(self) -> None:
        if self.schema_version != TOPOLOGY_SCOPE_VERSION:
            raise TopologyValidationError(
                f"Unsupported topology-scope version: {self.schema_version!r}"
            )
        if self.scope not in TOPOLOGY_SCOPE_CHOICES:
            raise TopologyValidationError(f"Unknown topology scope: {self.scope!r}")
        count_fields = (
            "declared_study_persons",
            "declared_common_anchor_persons",
            "declared_raters",
            "nodes",
            "unique_incidence_edges",
            "rating_sessions",
            "repeated_sessions_collapsed",
            "n_components",
            "largest_component_nodes",
            "largest_component_study_persons",
            "largest_component_raters",
            "isolated_raters",
        )
        for name in count_fields:
            _require_nonnegative_int(name, getattr(self, name))
        for name in (
            "declared_study_persons",
            "declared_raters",
            "nodes",
            "unique_incidence_edges",
            "rating_sessions",
            "n_components",
            "largest_component_nodes",
            "largest_component_raters",
        ):
            _require_positive_int(name, getattr(self, name))
        if not isinstance(self.connected, bool):
            raise TopologyValidationError("connected must be boolean")
        if not isinstance(self.components, tuple) or any(
            type(item) is not TopologyComponentV1 for item in self.components
        ):
            raise TopologyValidationError(
                "components must be a tuple of TopologyComponentV1 values"
            )
        if not isinstance(self.rater_loads, tuple) or any(
            type(item) is not RaterTopologyLoadV1 for item in self.rater_loads
        ):
            raise TopologyValidationError(
                "rater_loads must be a tuple of RaterTopologyLoadV1 values"
            )
        if not isinstance(self.incidence_edges, tuple) or any(
            type(item) is not TopologyIncidenceEdgeV1
            for item in self.incidence_edges
        ):
            raise TopologyValidationError(
                "incidence_edges must be a tuple of TopologyIncidenceEdgeV1 values"
            )
        expected_ids = tuple(
            f"C{index:03d}" for index in range(1, len(self.components) + 1)
        )
        if tuple(component.component_id for component in self.components) != expected_ids:
            raise TopologyValidationError(
                "components must use contiguous canonical IDs in deterministic order"
            )
        if self.n_components != len(self.components):
            raise TopologyValidationError("n_components does not match components")
        if self.connected != (self.n_components == 1):
            raise TopologyValidationError("connected does not match n_components")
        if self.declared_raters != len(self.rater_loads):
            raise TopologyValidationError(
                "declared_raters does not match the rater-load table"
            )
        if tuple(load.rater_id for load in self.rater_loads) != tuple(
            sorted(load.rater_id for load in self.rater_loads)
        ):
            raise TopologyValidationError("rater_loads must be sorted by rater_id")
        if len({load.rater_id for load in self.rater_loads}) != len(self.rater_loads):
            raise TopologyValidationError("rater_loads contain duplicate raters")

        component_study_people_list = [
            person_id
            for component in self.components
            for person_id in component.study_person_ids
        ]
        component_anchor_people_list = [
            person_id
            for component in self.components
            for person_id in component.common_anchor_person_ids
        ]
        component_raters_list = [
            rater_id
            for component in self.components
            for rater_id in component.rater_ids
        ]
        if len(component_study_people_list) != len(set(component_study_people_list)):
            raise TopologyValidationError(
                "study persons cannot appear in more than one component"
            )
        if len(component_anchor_people_list) != len(set(component_anchor_people_list)):
            raise TopologyValidationError(
                "common-anchor persons cannot appear in more than one component"
            )
        if len(component_raters_list) != len(set(component_raters_list)):
            raise TopologyValidationError(
                "raters cannot appear in more than one component"
            )
        component_study_people = set(component_study_people_list)
        component_anchor_people = set(component_anchor_people_list)
        component_raters = set(component_raters_list)
        if self.declared_study_persons != len(component_study_people):
            raise TopologyValidationError(
                "declared_study_persons does not match component membership"
            )
        if self.declared_common_anchor_persons != len(component_anchor_people):
            raise TopologyValidationError(
                "declared_common_anchor_persons does not match component membership"
            )
        if self.declared_raters != len(component_raters):
            raise TopologyValidationError(
                "declared_raters does not match component membership"
            )
        if component_raters != {load.rater_id for load in self.rater_loads}:
            raise TopologyValidationError(
                "component rater membership does not match rater_loads"
            )
        if self.nodes != sum(component.n_nodes for component in self.components):
            raise TopologyValidationError("nodes does not match component membership")
        if self.unique_incidence_edges != sum(
            component.unique_incidence_edges for component in self.components
        ):
            raise TopologyValidationError(
                "unique_incidence_edges does not match component totals"
            )
        if self.unique_incidence_edges != len(self.incidence_edges):
            raise TopologyValidationError(
                "unique_incidence_edges does not match incidence_edges"
            )
        edge_keys = tuple(edge.edge_key for edge in self.incidence_edges)
        if len(edge_keys) != len(set(edge_keys)):
            raise TopologyValidationError("incidence_edges contain duplicate keys")
        expected_edge_order = tuple(sorted(
            edge_keys,
            key=lambda key: (
                0 if key[0] == PERSON_ROLE_STUDY else 1,
                key[1],
                key[2],
            ),
        ))
        if edge_keys != expected_edge_order:
            raise TopologyValidationError(
                "incidence_edges must use deterministic person/rater order"
            )
        component_ids = {component.component_id for component in self.components}
        if any(edge.component_id not in component_ids for edge in self.incidence_edges):
            raise TopologyValidationError(
                "incidence edge refers to an unknown component"
            )
        component_by_person = {
            (PERSON_ROLE_STUDY, person_id): component.component_id
            for component in self.components
            for person_id in component.study_person_ids
        }
        component_by_person.update({
            (PERSON_ROLE_COMMON_ANCHOR, person_id): component.component_id
            for component in self.components
            for person_id in component.common_anchor_person_ids
        })
        component_by_rater = {
            rater_id: component.component_id
            for component in self.components
            for rater_id in component.rater_ids
        }
        edges_by_component: dict[str, list[TopologyIncidenceEdgeV1]] = defaultdict(list)
        for edge in self.incidence_edges:
            if component_by_person.get((edge.person_role, edge.person_id)) != edge.component_id:
                raise TopologyValidationError(
                    "incidence edge person does not belong to its component_id"
                )
            if component_by_rater.get(edge.rater_id) != edge.component_id:
                raise TopologyValidationError(
                    "incidence edge rater does not belong to its component_id"
                )
            edges_by_component[edge.component_id].append(edge)
        role_order = {
            PERSON_ROLE_STUDY: 0,
            PERSON_ROLE_COMMON_ANCHOR: 1,
            "rater": 2,
        }
        component_signatures = []
        for component in self.components:
            component_nodes = {
                *((PERSON_ROLE_STUDY, value) for value in component.study_person_ids),
                *((PERSON_ROLE_COMMON_ANCHOR, value) for value in component.common_anchor_person_ids),
                *(("rater", value) for value in component.rater_ids),
            }
            adjacency = {node: set() for node in component_nodes}
            for edge in edges_by_component[component.component_id]:
                person_node = (edge.person_role, edge.person_id)
                rater_node = ("rater", edge.rater_id)
                adjacency[person_node].add(rater_node)
                adjacency[rater_node].add(person_node)
            start = min(
                component_nodes,
                key=lambda node: (role_order[node[0]], node[1]),
            )
            visited = set()
            pending = [start]
            while pending:
                node = pending.pop()
                if node in visited:
                    continue
                visited.add(node)
                pending.extend(adjacency[node] - visited)
            if visited != component_nodes:
                raise TopologyValidationError(
                    "component membership is not connected by its incidence edges"
                )
            component_signatures.append(tuple(sorted(
                component_nodes,
                key=lambda node: (role_order[node[0]], node[1]),
            )))
        expected_component_order = sorted(
            component_signatures,
            key=lambda members: tuple(
                (role_order[node[0]], node[1]) for node in members
            ),
        )
        if component_signatures != expected_component_order:
            raise TopologyValidationError(
                "component membership does not use deterministic canonical order"
            )
        for component in self.components:
            component_edge_rows = edges_by_component[component.component_id]
            if len(component_edge_rows) != component.unique_incidence_edges:
                raise TopologyValidationError(
                    "component unique-edge count does not match incidence_edges"
                )
            if sum(edge.rating_sessions for edge in component_edge_rows) != component.rating_sessions:
                raise TopologyValidationError(
                    "component session count does not match edge multiplicities"
                )
        if self.rating_sessions != sum(
            component.rating_sessions for component in self.components
        ):
            raise TopologyValidationError(
                "rating_sessions does not match component totals"
            )
        if self.rating_sessions != sum(
            edge.rating_sessions for edge in self.incidence_edges
        ):
            raise TopologyValidationError(
                "rating_sessions does not match incidence-edge multiplicities"
            )
        if self.repeated_sessions_collapsed != (
            self.rating_sessions - self.unique_incidence_edges
        ):
            raise TopologyValidationError(
                "repeated_sessions_collapsed must equal sessions minus unique edges"
            )
        if self.rating_sessions != sum(
            load.rating_sessions for load in self.rater_loads
        ):
            raise TopologyValidationError(
                "rating_sessions does not match rater-load totals"
            )
        edges_by_rater: dict[str, list[TopologyIncidenceEdgeV1]] = defaultdict(list)
        for edge in self.incidence_edges:
            edges_by_rater[edge.rater_id].append(edge)
        for load in self.rater_loads:
            rater_edges = edges_by_rater[load.rater_id]
            if load.rating_sessions != sum(edge.rating_sessions for edge in rater_edges):
                raise TopologyValidationError(
                    "rater session load does not match incidence-edge multiplicities"
                )
            if load.unique_persons != len(rater_edges):
                raise TopologyValidationError(
                    "rater unique-person load does not match incidence edges"
                )
            if load.study_persons != sum(
                edge.person_role == PERSON_ROLE_STUDY for edge in rater_edges
            ):
                raise TopologyValidationError(
                    "rater study-person load does not match incidence edges"
                )
            if load.common_anchor_persons != sum(
                edge.person_role == PERSON_ROLE_COMMON_ANCHOR
                for edge in rater_edges
            ):
                raise TopologyValidationError(
                    "rater common-anchor load does not match incidence edges"
                )
        if self.isolated_raters != sum(
            load.rating_sessions == 0 for load in self.rater_loads
        ):
            raise TopologyValidationError(
                "isolated_raters does not match zero-session rater loads"
            )
        if self.components:
            largest = max(
                self.components,
                key=lambda item: (
                    item.n_nodes,
                    len(item.study_person_ids),
                    len(item.rater_ids),
                    item.component_id,
                ),
            )
            if self.largest_component_nodes != largest.n_nodes:
                raise TopologyValidationError(
                    "largest_component_nodes does not match components"
                )
            if self.largest_component_study_persons != len(largest.study_person_ids):
                raise TopologyValidationError(
                    "largest_component_study_persons does not match the deterministic largest component"
                )
            if self.largest_component_raters != len(largest.rater_ids):
                raise TopologyValidationError(
                    "largest_component_raters does not match the deterministic largest component"
                )

    def to_dict(self) -> dict:
        return _json_payload({
            "schema_version": self.schema_version,
            "scope": self.scope,
            "declared_study_persons": self.declared_study_persons,
            "declared_common_anchor_persons": self.declared_common_anchor_persons,
            "declared_raters": self.declared_raters,
            "nodes": self.nodes,
            "unique_incidence_edges": self.unique_incidence_edges,
            "rating_sessions": self.rating_sessions,
            "repeated_sessions_collapsed": self.repeated_sessions_collapsed,
            "n_components": self.n_components,
            "largest_component_nodes": self.largest_component_nodes,
            "largest_component_study_persons": self.largest_component_study_persons,
            "largest_component_raters": self.largest_component_raters,
            "isolated_raters": self.isolated_raters,
            "connected": self.connected,
            "components": [component.to_dict() for component in self.components],
            "incidence_edges": [edge.to_dict() for edge in self.incidence_edges],
            "rater_loads": [load.to_dict() for load in self.rater_loads],
        })


@dataclass(frozen=True, slots=True, kw_only=True)
class TopologyAuditV1:
    """Versioned topology evidence bound to one exact assignment bundle."""

    schema_version: str
    assignment_bundle_schema_version: str
    design_spec_fingerprint: str
    assignment_fingerprint: str
    graph_semantics: str
    primary_study: TopologyScopeAuditV1
    anchor_augmented: TopologyScopeAuditV1

    def __post_init__(self) -> None:
        if self.schema_version != TOPOLOGY_AUDIT_VERSION:
            raise TopologyValidationError(
                f"Unsupported topology-audit version: {self.schema_version!r}"
            )
        if self.assignment_bundle_schema_version != ASSIGNMENT_BUNDLE_VERSION:
            raise TopologyValidationError(
                "topology audit requires the current assignment-bundle version"
            )
        _require_lowercase_fingerprint(
            "design_spec_fingerprint", self.design_spec_fingerprint
        )
        _require_lowercase_fingerprint(
            "assignment_fingerprint", self.assignment_fingerprint
        )
        if self.graph_semantics != "person_rater_bipartite_v1":
            raise TopologyValidationError(
                f"Unsupported graph semantics: {self.graph_semantics!r}"
            )
        if type(self.primary_study) is not TopologyScopeAuditV1 or (
            self.primary_study.scope != PRIMARY_STUDY_SCOPE
        ):
            raise TopologyValidationError(
                "primary_study must contain the primary-study scope"
            )
        if type(self.anchor_augmented) is not TopologyScopeAuditV1 or (
            self.anchor_augmented.scope != ANCHOR_AUGMENTED_SCOPE
        ):
            raise TopologyValidationError(
                "anchor_augmented must contain the anchor-augmented scope"
            )
        if (
            self.primary_study.declared_study_persons
            != self.anchor_augmented.declared_study_persons
            or self.primary_study.declared_raters
            != self.anchor_augmented.declared_raters
        ):
            raise TopologyValidationError(
                "both topology scopes must retain the same study persons and raters"
            )
        if self.primary_study.declared_common_anchor_persons != 0:
            raise TopologyValidationError(
                "the primary-study scope cannot contain common-anchor persons"
            )
        if (
            self.anchor_augmented.rating_sessions
            < self.primary_study.rating_sessions
            or self.anchor_augmented.unique_incidence_edges
            < self.primary_study.unique_incidence_edges
        ):
            raise TopologyValidationError(
                "adding common anchors cannot remove sessions or incidence edges"
            )
        primary_study_people = {
            person_id
            for component in self.primary_study.components
            for person_id in component.study_person_ids
        }
        augmented_study_people = {
            person_id
            for component in self.anchor_augmented.components
            for person_id in component.study_person_ids
        }
        primary_raters = {
            rater_id
            for component in self.primary_study.components
            for rater_id in component.rater_ids
        }
        augmented_raters = {
            rater_id
            for component in self.anchor_augmented.components
            for rater_id in component.rater_ids
        }
        if primary_study_people != augmented_study_people:
            raise TopologyValidationError(
                "both topology scopes must retain the same study-person IDs"
            )
        if primary_raters != augmented_raters:
            raise TopologyValidationError(
                "both topology scopes must retain the same rater IDs"
            )
        primary_edge_multiplicity = {
            edge.edge_key: edge.rating_sessions
            for edge in self.primary_study.incidence_edges
        }
        augmented_study_edge_multiplicity = {
            edge.edge_key: edge.rating_sessions
            for edge in self.anchor_augmented.incidence_edges
            if edge.person_role == PERSON_ROLE_STUDY
        }
        if primary_edge_multiplicity != augmented_study_edge_multiplicity:
            raise TopologyValidationError(
                "common-anchor augmentation must preserve every study edge multiplicity"
            )

    @property
    def cache_identity(self) -> tuple[str, str, str, str]:
        return (
            self.schema_version,
            self.assignment_bundle_schema_version,
            self.design_spec_fingerprint,
            self.assignment_fingerprint,
        )

    def to_dict(self) -> dict:
        return _json_payload({
            "schema_version": self.schema_version,
            "assignment_bundle_schema_version": self.assignment_bundle_schema_version,
            "design_spec_fingerprint": self.design_spec_fingerprint,
            "assignment_fingerprint": self.assignment_fingerprint,
            "graph_semantics": self.graph_semantics,
            "primary_study": self.primary_study.to_dict(),
            "anchor_augmented": self.anchor_augmented.to_dict(),
        })


@dataclass(frozen=True, slots=True, kw_only=True)
class StructuralBridgeSessionV1:
    """One proposed extra rating session; never part of the base bundle."""

    schema_version: str
    repair_index: int
    person_id: str
    artifact_id: str
    additional_rater_id: str
    source_component_id: str
    target_component_id: str
    target_group_id: str
    added_score_records: int

    def __post_init__(self) -> None:
        if self.schema_version != STRUCTURAL_REPAIR_SESSION_VERSION:
            raise TopologyValidationError(
                f"Unsupported repair-session version: {self.schema_version!r}"
            )
        _require_positive_int("repair_index", self.repair_index)
        for name in (
            "person_id",
            "artifact_id",
            "additional_rater_id",
            "target_group_id",
        ):
            value = getattr(self, name)
            if not isinstance(value, str) or not value:
                raise TopologyValidationError(f"{name} must be a non-empty string")
        for name in ("source_component_id", "target_component_id"):
            value = getattr(self, name)
            if (
                not isinstance(value, str)
                or _COMPONENT_ID_PATTERN.fullmatch(value) is None
            ):
                raise TopologyValidationError(
                    f"{name} must use canonical component numbering"
                )
        if self.source_component_id == self.target_component_id:
            raise TopologyValidationError(
                "a structural bridge must join two different base components"
            )
        _require_positive_int("added_score_records", self.added_score_records)

    def to_dict(self) -> dict:
        return _json_payload({
            "schema_version": self.schema_version,
            "repair_index": self.repair_index,
            "person_id": self.person_id,
            "artifact_id": self.artifact_id,
            "additional_rater_id": self.additional_rater_id,
            "source_component_id": self.source_component_id,
            "target_component_id": self.target_component_id,
            "target_group_id": self.target_group_id,
            "added_score_records": self.added_score_records,
        })


@dataclass(frozen=True, slots=True, kw_only=True)
class StructuralBridgeRepairV1:
    """A K-1 structural overlay bound to, but separate from, a base schedule."""

    schema_version: str
    design_spec_fingerprint: str
    assignment_fingerprint: str
    topology_audit_schema_version: str
    topology_scope: str
    status: str
    target_person_id: str
    target_artifact_id: str
    target_existing_rater_ids: tuple[str, ...]
    components_before: int
    base_component_ids: tuple[str, ...]
    projected_components_after: int
    structural_minimum_rating_sessions: int
    added_rating_sessions: int
    added_score_records: int
    added_calibration_events: int
    overlay_sessions: tuple[StructuralBridgeSessionV1, ...]
    repair_fingerprint: str
    structural_connectivity_only: bool
    precision_evaluation_status: str
    expected_se_improvement_percent: None
    classification_agreement_improvement_percent: None

    def __post_init__(self) -> None:
        if self.schema_version != STRUCTURAL_REPAIR_VERSION:
            raise TopologyValidationError(
                f"Unsupported structural-repair version: {self.schema_version!r}"
            )
        _require_lowercase_fingerprint(
            "design_spec_fingerprint", self.design_spec_fingerprint
        )
        _require_lowercase_fingerprint(
            "assignment_fingerprint", self.assignment_fingerprint
        )
        if self.topology_audit_schema_version != TOPOLOGY_AUDIT_VERSION:
            raise TopologyValidationError(
                "repair proposal requires the current topology-audit version"
            )
        if self.topology_scope != ANCHOR_AUGMENTED_SCOPE:
            raise TopologyValidationError(
                "v1 repairs must target the anchor-augmented topology"
            )
        if self.status not in REPAIR_STATUS_CHOICES:
            raise TopologyValidationError(f"Unknown repair status: {self.status!r}")
        for name in ("target_person_id", "target_artifact_id"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value:
                raise TopologyValidationError(f"{name} must be a non-empty string")
        _require_sorted_unique_strings(
            "target_existing_rater_ids",
            self.target_existing_rater_ids,
        )
        if not self.target_existing_rater_ids:
            raise TopologyValidationError(
                "a repair target must already have at least one rater"
            )
        for name in (
            "components_before",
            "projected_components_after",
            "structural_minimum_rating_sessions",
            "added_rating_sessions",
            "added_score_records",
            "added_calibration_events",
        ):
            _require_nonnegative_int(name, getattr(self, name))
        if self.components_before < 1:
            raise TopologyValidationError("components_before must be positive")
        if not isinstance(self.base_component_ids, tuple) or any(
            not isinstance(value, str) or not value
            for value in self.base_component_ids
        ):
            raise TopologyValidationError(
                "base_component_ids must be a tuple of non-empty strings"
            )
        if len(self.base_component_ids) != len(set(self.base_component_ids)):
            raise TopologyValidationError("base_component_ids must be unique")
        expected_component_ids = tuple(
            f"C{index:03d}" for index in range(1, self.components_before + 1)
        )
        if self.base_component_ids != expected_component_ids:
            raise TopologyValidationError(
                "base_component_ids must exactly enumerate the audited base components"
            )
        if self.projected_components_after != 1:
            raise TopologyValidationError(
                "a complete v1 repair must project one connected component"
            )
        expected_minimum = max(0, self.components_before - 1)
        if self.structural_minimum_rating_sessions != expected_minimum:
            raise TopologyValidationError(
                "structural_minimum_rating_sessions must equal components_before - 1"
            )
        if not isinstance(self.overlay_sessions, tuple) or any(
            type(item) is not StructuralBridgeSessionV1
            for item in self.overlay_sessions
        ):
            raise TopologyValidationError(
                "overlay_sessions must be a tuple of StructuralBridgeSessionV1 values"
            )
        if tuple(item.repair_index for item in self.overlay_sessions) != tuple(
            range(1, len(self.overlay_sessions) + 1)
        ):
            raise TopologyValidationError(
                "repair sessions must use contiguous one-based indices"
            )
        if self.added_rating_sessions != len(self.overlay_sessions):
            raise TopologyValidationError(
                "added_rating_sessions does not match overlay_sessions"
            )
        if self.added_rating_sessions != self.structural_minimum_rating_sessions:
            raise TopologyValidationError(
                "v1 repair must use the K-1 structural minimum"
            )
        if self.added_score_records != sum(
            item.added_score_records for item in self.overlay_sessions
        ):
            raise TopologyValidationError(
                "added_score_records does not match overlay session totals"
            )
        if self.added_calibration_events != 0:
            raise TopologyValidationError(
                "structural bridge sessions cannot create calibration events"
            )
        expected_status = REPAIR_NOT_NEEDED if self.components_before == 1 else REPAIR_PROPOSED
        if self.status != expected_status:
            raise TopologyValidationError("repair status does not match components_before")
        source_ids = {item.source_component_id for item in self.overlay_sessions}
        target_ids = [item.target_component_id for item in self.overlay_sessions]
        if self.overlay_sessions and len(source_ids) != 1:
            raise TopologyValidationError(
                "all v1 repair sessions must use one deterministic hub component"
            )
        if len(target_ids) != len(set(target_ids)):
            raise TopologyValidationError(
                "each non-hub base component may appear only once in a v1 repair"
            )
        for session in self.overlay_sessions:
            if (
                session.person_id != self.target_person_id
                or session.artifact_id != self.target_artifact_id
            ):
                raise TopologyValidationError(
                    "every v1 repair session must use the declared common target artifact"
                )
            if session.additional_rater_id in self.target_existing_rater_ids:
                raise TopologyValidationError(
                    "an additional rater cannot already rate the target artifact"
                )
        if self.overlay_sessions:
            source_id = next(iter(source_ids))
            if source_id not in self.base_component_ids:
                raise TopologyValidationError(
                    "repair hub is not present in base_component_ids"
                )
            if set(target_ids) != set(self.base_component_ids) - {source_id}:
                raise TopologyValidationError(
                    "repair targets must cover every non-hub base component exactly once"
                )
        elif self.base_component_ids != ("C001",):
            raise TopologyValidationError(
                "an empty repair overlay is valid only for one base component"
            )
        _require_lowercase_fingerprint("repair_fingerprint", self.repair_fingerprint)
        expected_fingerprint = _repair_fingerprint(
            design_spec_fingerprint=self.design_spec_fingerprint,
            assignment_fingerprint=self.assignment_fingerprint,
            base_component_ids=self.base_component_ids,
            target_person_id=self.target_person_id,
            target_artifact_id=self.target_artifact_id,
            target_existing_rater_ids=self.target_existing_rater_ids,
            sessions=self.overlay_sessions,
        )
        if self.repair_fingerprint != expected_fingerprint:
            raise TopologyValidationError(
                "repair_fingerprint does not match the base assignment and overlay"
            )
        if self.structural_connectivity_only is not True:
            raise TopologyValidationError(
                "structural_connectivity_only must remain true in v1"
            )
        if self.precision_evaluation_status != PRECISION_NOT_EVALUATED:
            raise TopologyValidationError(
                "v1 structural repair cannot claim evaluated precision"
            )
        if self.expected_se_improvement_percent is not None:
            raise TopologyValidationError(
                "expected SE improvement is not evaluated by a structural repair"
            )
        if self.classification_agreement_improvement_percent is not None:
            raise TopologyValidationError(
                "classification agreement is not evaluated by a structural repair"
            )

    def to_dict(self) -> dict:
        return _json_payload({
            "schema_version": self.schema_version,
            "design_spec_fingerprint": self.design_spec_fingerprint,
            "assignment_fingerprint": self.assignment_fingerprint,
            "topology_audit_schema_version": self.topology_audit_schema_version,
            "topology_scope": self.topology_scope,
            "status": self.status,
            "target_person_id": self.target_person_id,
            "target_artifact_id": self.target_artifact_id,
            "target_existing_rater_ids": list(self.target_existing_rater_ids),
            "components_before": self.components_before,
            "base_component_ids": list(self.base_component_ids),
            "projected_components_after": self.projected_components_after,
            "structural_minimum_rating_sessions": self.structural_minimum_rating_sessions,
            "added_rating_sessions": self.added_rating_sessions,
            "added_score_records": self.added_score_records,
            "added_calibration_events": self.added_calibration_events,
            "overlay_sessions": [item.to_dict() for item in self.overlay_sessions],
            "repair_fingerprint": self.repair_fingerprint,
            "structural_connectivity_only": self.structural_connectivity_only,
            "precision_evaluation_status": self.precision_evaluation_status,
            "expected_se_improvement_percent": self.expected_se_improvement_percent,
            "classification_agreement_improvement_percent": (
                self.classification_agreement_improvement_percent
            ),
        })


class _UnionFind:
    def __init__(self, nodes: Sequence[tuple[str, str]]) -> None:
        self.parent = {node: node for node in nodes}
        self.rank = {node: 0 for node in nodes}

    def find(self, node: tuple[str, str]) -> tuple[str, str]:
        parent = self.parent[node]
        if parent != node:
            self.parent[node] = self.find(parent)
        return self.parent[node]

    def union(self, left: tuple[str, str], right: tuple[str, str]) -> None:
        root_left = self.find(left)
        root_right = self.find(right)
        if root_left == root_right:
            return
        if self.rank[root_left] < self.rank[root_right]:
            root_left, root_right = root_right, root_left
        self.parent[root_right] = root_left
        if self.rank[root_left] == self.rank[root_right]:
            self.rank[root_left] += 1


def _node_sort_key(node: tuple[str, str]) -> tuple[int, str]:
    role_order = {
        PERSON_ROLE_STUDY: 0,
        PERSON_ROLE_COMMON_ANCHOR: 1,
        "rater": 2,
    }
    return role_order[node[0]], node[1]


def _scope_rows(
    bundle: AssignmentBundleV1,
    scope: str,
) -> tuple[RatingAssignmentV1, ...]:
    if scope == PRIMARY_STUDY_SCOPE:
        return tuple(
            row
            for row in bundle.assignments
            if row.person_role == PERSON_ROLE_STUDY
        )
    if scope == ANCHOR_AUGMENTED_SCOPE:
        return bundle.assignments
    raise TopologyValidationError(f"Unknown topology scope: {scope!r}")


def _build_scope_audit(
    bundle: AssignmentBundleV1,
    scope: str,
) -> TopologyScopeAuditV1:
    rows = _scope_rows(bundle, scope)
    all_study_people = tuple(sorted({
        row.person_id
        for row in bundle.assignments
        if row.person_role == PERSON_ROLE_STUDY
    }))
    all_raters = tuple(sorted({row.rater_id for row in bundle.assignments}))
    anchor_people = tuple(sorted({
        row.person_id
        for row in rows
        if row.person_role == PERSON_ROLE_COMMON_ANCHOR
    }))
    nodes = tuple(
        [(PERSON_ROLE_STUDY, person_id) for person_id in all_study_people]
        + [
            (PERSON_ROLE_COMMON_ANCHOR, person_id)
            for person_id in anchor_people
        ]
        + [("rater", rater_id) for rater_id in all_raters]
    )
    union_find = _UnionFind(nodes)
    edge_multiplicity: Counter[
        tuple[tuple[str, str], tuple[str, str]]
    ] = Counter()
    for row in rows:
        person_node = (row.person_role, row.person_id)
        rater_node = ("rater", row.rater_id)
        union_find.union(person_node, rater_node)
        edge_multiplicity[(person_node, rater_node)] += 1

    members_by_root: dict[tuple[str, str], list[tuple[str, str]]] = defaultdict(list)
    for node in nodes:
        members_by_root[union_find.find(node)].append(node)
    ordered_members = sorted(
        (tuple(sorted(members, key=_node_sort_key)) for members in members_by_root.values()),
        key=lambda members: tuple(_node_sort_key(node) for node in members),
    )
    component_id_by_node: dict[tuple[str, str], str] = {}
    for index, members in enumerate(ordered_members, start=1):
        component_id = f"C{index:03d}"
        for node in members:
            component_id_by_node[node] = component_id

    edge_count_by_component = Counter(
        component_id_by_node[person_node]
        for person_node, _ in edge_multiplicity
    )
    session_count_by_component = Counter(
        component_id_by_node[(row.person_role, row.person_id)] for row in rows
    )
    components: list[TopologyComponentV1] = []
    for index, members in enumerate(ordered_members, start=1):
        component_id = f"C{index:03d}"
        components.append(TopologyComponentV1(
            schema_version=TOPOLOGY_COMPONENT_VERSION,
            component_id=component_id,
            study_person_ids=tuple(
                node_id for role, node_id in members if role == PERSON_ROLE_STUDY
            ),
            common_anchor_person_ids=tuple(
                node_id
                for role, node_id in members
                if role == PERSON_ROLE_COMMON_ANCHOR
            ),
            rater_ids=tuple(
                node_id for role, node_id in members if role == "rater"
            ),
            unique_incidence_edges=edge_count_by_component[component_id],
            rating_sessions=session_count_by_component[component_id],
        ))

    incidence_edges = tuple(
        TopologyIncidenceEdgeV1(
            schema_version=TOPOLOGY_EDGE_VERSION,
            component_id=component_id_by_node[person_node],
            person_role=person_node[0],
            person_id=person_node[1],
            rater_id=rater_node[1],
            rating_sessions=multiplicity,
        )
        for (person_node, rater_node), multiplicity in sorted(
            edge_multiplicity.items(),
            key=lambda item: (
                _node_sort_key(item[0][0]),
                item[0][1][1],
            ),
        )
    )

    people_by_rater: dict[str, set[tuple[str, str]]] = defaultdict(set)
    session_loads = Counter(row.rater_id for row in rows)
    for row in rows:
        people_by_rater[row.rater_id].add((row.person_role, row.person_id))
    rater_loads = tuple(
        RaterTopologyLoadV1(
            schema_version=RATER_LOAD_VERSION,
            rater_id=rater_id,
            rating_sessions=session_loads[rater_id],
            unique_persons=len(people_by_rater[rater_id]),
            study_persons=sum(
                role == PERSON_ROLE_STUDY
                for role, _ in people_by_rater[rater_id]
            ),
            common_anchor_persons=sum(
                role == PERSON_ROLE_COMMON_ANCHOR
                for role, _ in people_by_rater[rater_id]
            ),
        )
        for rater_id in all_raters
    )
    largest = max(
        components,
        key=lambda item: (
            item.n_nodes,
            len(item.study_person_ids),
            len(item.rater_ids),
            item.component_id,
        ),
    )
    return TopologyScopeAuditV1(
        schema_version=TOPOLOGY_SCOPE_VERSION,
        scope=scope,
        declared_study_persons=len(all_study_people),
        declared_common_anchor_persons=len(anchor_people),
        declared_raters=len(all_raters),
        nodes=len(nodes),
        unique_incidence_edges=len(edge_multiplicity),
        rating_sessions=len(rows),
        repeated_sessions_collapsed=len(rows) - len(edge_multiplicity),
        n_components=len(components),
        largest_component_nodes=largest.n_nodes,
        largest_component_study_persons=len(largest.study_person_ids),
        largest_component_raters=len(largest.rater_ids),
        isolated_raters=sum(load.rating_sessions == 0 for load in rater_loads),
        connected=len(components) == 1,
        components=tuple(components),
        incidence_edges=incidence_edges,
        rater_loads=rater_loads,
    )


def audit_rating_design_topology(
    value: AssignmentBundleV1 | Mapping,
    *,
    max_rating_sessions: int = DEFAULT_MAX_RATING_SESSIONS,
) -> TopologyAuditV1:
    """Audit primary and anchor-augmented Person-Rater incidence graphs."""
    bundle = normalize_assignment_bundle(
        value,
        max_rating_sessions=max_rating_sessions,
    )
    return TopologyAuditV1(
        schema_version=TOPOLOGY_AUDIT_VERSION,
        assignment_bundle_schema_version=bundle.schema_version,
        design_spec_fingerprint=bundle.design_spec_fingerprint,
        assignment_fingerprint=bundle.assignment_fingerprint,
        graph_semantics="person_rater_bipartite_v1",
        primary_study=_build_scope_audit(bundle, PRIMARY_STUDY_SCOPE),
        anchor_augmented=_build_scope_audit(bundle, ANCHOR_AUGMENTED_SCOPE),
    )


def _repair_fingerprint(
    *,
    design_spec_fingerprint: str,
    assignment_fingerprint: str,
    base_component_ids: Sequence[str],
    target_person_id: str,
    target_artifact_id: str,
    target_existing_rater_ids: Sequence[str],
    sessions: Sequence[StructuralBridgeSessionV1],
) -> str:
    payload = {
        "schema_version": STRUCTURAL_REPAIR_VERSION,
        "design_spec_fingerprint": design_spec_fingerprint,
        "assignment_fingerprint": assignment_fingerprint,
        "topology_scope": ANCHOR_AUGMENTED_SCOPE,
        "base_component_ids": list(base_component_ids),
        "target_person_id": target_person_id,
        "target_artifact_id": target_artifact_id,
        "target_existing_rater_ids": list(target_existing_rater_ids),
        "overlay_sessions": [session.to_dict() for session in sessions],
    }
    encoded = json.dumps(
        payload,
        allow_nan=False,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()[:16]


def propose_structural_bridge_repair(
    value: AssignmentBundleV1 | Mapping,
    *,
    topology_audit: TopologyAuditV1 | None = None,
    max_rating_sessions: int = DEFAULT_MAX_RATING_SESSIONS,
) -> StructuralBridgeRepairV1:
    """Propose the deterministic K-1 overlay for the augmented graph.

    Every overlay row sends one rater from a non-hub base component to an
    already-rated study artifact in the hub component.  This is a graph repair
    only.  The returned null improvement fields are intentional: expected SE,
    recovery, and classification effects require paired simulation under the
    same data-generating truths and random draws.
    """
    bundle = normalize_assignment_bundle(
        value,
        max_rating_sessions=max_rating_sessions,
    )
    expected_audit = audit_rating_design_topology(
        bundle,
        max_rating_sessions=max_rating_sessions,
    )
    audit = topology_audit or expected_audit
    if type(audit) is not TopologyAuditV1:
        raise TypeError("topology_audit must be a TopologyAuditV1 or None")
    if (
        audit.design_spec_fingerprint != bundle.design_spec_fingerprint
        or audit.assignment_fingerprint != bundle.assignment_fingerprint
        or audit.assignment_bundle_schema_version != bundle.schema_version
    ):
        raise TopologyValidationError(
            "topology_audit does not belong to the supplied assignment bundle"
        )
    if audit != expected_audit:
        raise TopologyValidationError(
            "topology_audit does not exactly reproduce the supplied assignment bundle"
        )

    scope = audit.anchor_augmented
    components = scope.components
    hub = min(
        components,
        key=lambda component: (
            -len(component.study_person_ids),
            -len(component.rater_ids),
            component.component_id,
        ),
    )
    if not hub.study_person_ids:
        raise TopologyValidationError(
            "repair hub must contain at least one study person"
        )
    target_person_id = hub.study_person_ids[0]
    primary_rows = sorted(
        (
            row
            for row in bundle.assignments
            if row.person_id == target_person_id
            and row.person_role == PERSON_ROLE_STUDY
            and row.session_type == SESSION_PRIMARY
        ),
        key=lambda row: (row.artifact_id, row.rater_id, row.group_id),
    )
    if not primary_rows:
        raise TopologyValidationError(
            "repair hub study person has no primary artifact"
        )
    target_artifact_id = primary_rows[0].artifact_id
    existing_rater_ids = tuple(sorted({
        row.rater_id
        for row in bundle.assignments
        if row.person_id == target_person_id
        and row.artifact_id == target_artifact_id
        and row.person_role == PERSON_ROLE_STUDY
    }))
    group_ids_by_rater: dict[str, set[str]] = defaultdict(set)
    for row in bundle.assignments:
        group_ids_by_rater[row.rater_id].add(row.group_id)

    overlay_sessions: list[StructuralBridgeSessionV1] = []
    for component in components:
        if component.component_id == hub.component_id:
            continue
        if not component.rater_ids:
            raise TopologyValidationError(
                "every repair target component must contain a declared rater"
            )
        additional_rater_id = component.rater_ids[0]
        target_groups = sorted(group_ids_by_rater[additional_rater_id])
        if not target_groups:
            raise TopologyValidationError(
                "repair target rater has no group identity in the base schedule"
            )
        overlay_sessions.append(StructuralBridgeSessionV1(
            schema_version=STRUCTURAL_REPAIR_SESSION_VERSION,
            repair_index=len(overlay_sessions) + 1,
            person_id=target_person_id,
            artifact_id=target_artifact_id,
            additional_rater_id=additional_rater_id,
            source_component_id=hub.component_id,
            target_component_id=component.component_id,
            target_group_id=target_groups[0],
            added_score_records=bundle.design_spec.score_records_per_session,
        ))

    sessions = tuple(overlay_sessions)
    fingerprint = _repair_fingerprint(
        design_spec_fingerprint=bundle.design_spec_fingerprint,
        assignment_fingerprint=bundle.assignment_fingerprint,
        base_component_ids=tuple(component.component_id for component in components),
        target_person_id=target_person_id,
        target_artifact_id=target_artifact_id,
        target_existing_rater_ids=existing_rater_ids,
        sessions=sessions,
    )
    return StructuralBridgeRepairV1(
        schema_version=STRUCTURAL_REPAIR_VERSION,
        design_spec_fingerprint=bundle.design_spec_fingerprint,
        assignment_fingerprint=bundle.assignment_fingerprint,
        topology_audit_schema_version=audit.schema_version,
        topology_scope=ANCHOR_AUGMENTED_SCOPE,
        status=REPAIR_NOT_NEEDED if scope.connected else REPAIR_PROPOSED,
        target_person_id=target_person_id,
        target_artifact_id=target_artifact_id,
        target_existing_rater_ids=existing_rater_ids,
        components_before=scope.n_components,
        base_component_ids=tuple(component.component_id for component in components),
        projected_components_after=1,
        structural_minimum_rating_sessions=max(0, scope.n_components - 1),
        added_rating_sessions=len(sessions),
        added_score_records=sum(item.added_score_records for item in sessions),
        added_calibration_events=0,
        overlay_sessions=sessions,
        repair_fingerprint=fingerprint,
        structural_connectivity_only=True,
        precision_evaluation_status=PRECISION_NOT_EVALUATED,
        expected_se_improvement_percent=None,
        classification_agreement_improvement_percent=None,
    )


def topology_audit_to_dict(value: TopologyAuditV1) -> dict:
    if type(value) is not TopologyAuditV1:
        raise TypeError("value must be a TopologyAuditV1")
    return value.to_dict()


def structural_repair_to_dict(value: StructuralBridgeRepairV1) -> dict:
    if type(value) is not StructuralBridgeRepairV1:
        raise TypeError("value must be a StructuralBridgeRepairV1")
    return value.to_dict()


__all__ = [
    "ANCHOR_AUGMENTED_SCOPE",
    "PRECISION_NOT_EVALUATED",
    "PRIMARY_STUDY_SCOPE",
    "RATER_LOAD_VERSION",
    "REPAIR_NOT_NEEDED",
    "REPAIR_PROPOSED",
    "REPAIR_STATUS_CHOICES",
    "STRUCTURAL_REPAIR_SESSION_VERSION",
    "STRUCTURAL_REPAIR_VERSION",
    "TOPOLOGY_AUDIT_VERSION",
    "TOPOLOGY_COMPONENT_VERSION",
    "TOPOLOGY_EDGE_VERSION",
    "TOPOLOGY_SCOPE_CHOICES",
    "TOPOLOGY_SCOPE_VERSION",
    "RaterTopologyLoadV1",
    "StructuralBridgeRepairV1",
    "StructuralBridgeSessionV1",
    "TopologyAuditV1",
    "TopologyComponentV1",
    "TopologyIncidenceEdgeV1",
    "TopologyScopeAuditV1",
    "TopologyValidationError",
    "audit_rating_design_topology",
    "propose_structural_bridge_repair",
    "structural_repair_to_dict",
    "topology_audit_to_dict",
]
