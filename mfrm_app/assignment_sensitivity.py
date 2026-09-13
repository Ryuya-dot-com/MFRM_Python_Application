"""Fail-closed fixed-density assignment perturbation primitives.

This module changes *design rows*, never observed outcomes.  It supports a
strict initial scope in which every observed Person-Rater assignment carries
the same non-rater context signature.  Under that condition a bipartite
degree-preserving 2-switch also preserves response-row exposure exactly.

Counterfactual scores must be simulated after perturbation by a caller that
owns the fitted response model.  Returning the observed Score on a changed
Rater row is intentionally prohibited.
"""

from __future__ import annotations

from itertools import combinations
from typing import Iterable, Mapping

import numpy as np
import pandas as pd


SENSITIVITY_SCHEMA_VERSION = "fixed_density_assignment_sensitivity_v1"
MAX_EXACT_SWITCH_EDGES = 500
MAX_EXACT_SWITCH_ROWS = 100_000


def _clean_string_frame(data: pd.DataFrame, columns: Iterable[str]) -> pd.DataFrame:
    out = data[list(columns)].copy()
    for column in columns:
        out[column] = out[column].astype(str)
    return out


def _rater_overlap_component_count(edges: set[tuple[str, str]]) -> int:
    raters = sorted({rater for _, rater in edges})
    if not raters:
        return 0
    people_by_rater = {
        rater: {person for person, candidate in edges if candidate == rater}
        for rater in raters
    }
    adjacency = {rater: set() for rater in raters}
    for left, right in combinations(raters, 2):
        if people_by_rater[left] & people_by_rater[right]:
            adjacency[left].add(right)
            adjacency[right].add(left)
    remaining = set(raters)
    components = 0
    while remaining:
        components += 1
        stack = [remaining.pop()]
        while stack:
            node = stack.pop()
            unseen = adjacency[node] & remaining
            remaining.difference_update(unseen)
            stack.extend(unseen)
    return components


def _rank_coordinate(values: Mapping[str, float], levels: list[str]) -> dict[str, float]:
    series = pd.Series({level: values.get(level, np.nan) for level in levels}, dtype=float)
    if series.isna().any() or not np.isfinite(series.to_numpy(dtype=float)).all():
        missing = series.index[~np.isfinite(series.to_numpy(dtype=float))].tolist()
        raise ValueError("Non-finite ordering coordinate for: " + ", ".join(map(str, missing[:10])))
    if series.nunique() < 2:
        raise ValueError("Ordering coordinates must contain at least two distinct values.")
    ranks = series.rank(method="average", pct=True)
    return {str(level): float(ranks.loc[level] - 0.5) for level in levels}


def _assignment_objective(
    edges: set[tuple[str, str]],
    person_coordinate: Mapping[str, float],
    rater_coordinate: Mapping[str, float],
) -> float:
    if not edges:
        return float("nan")
    return float(
        np.mean(
            [person_coordinate[person] * rater_coordinate[rater] for person, rater in edges]
        )
    )


def _assignment_correlation(
    edges: set[tuple[str, str]],
    person_coordinate: Mapping[str, float],
    rater_coordinate: Mapping[str, float],
) -> float:
    if len(edges) < 2:
        return float("nan")
    left = np.asarray([person_coordinate[person] for person, _ in sorted(edges)], dtype=float)
    right = np.asarray([rater_coordinate[rater] for _, rater in sorted(edges)], dtype=float)
    if np.std(left) <= 0 or np.std(right) <= 0:
        return float("nan")
    return float(np.corrcoef(left, right)[0, 1])


def _context_signature_table(
    work: pd.DataFrame,
    *,
    person_col: str,
    rater_col: str,
    context_cols: list[str],
) -> tuple[pd.DataFrame, bool]:
    signature_rows: list[dict[str, object]] = []
    duplicates_found = False
    for (person, rater), block in work.groupby([person_col, rater_col], observed=True, sort=True):
        if context_cols:
            contexts = [tuple(row) for row in block[context_cols].astype(str).itertuples(index=False, name=None)]
        else:
            contexts = [("__single_context__",)] * len(block)
        duplicates = len(contexts) - len(set(contexts))
        duplicates_found = duplicates_found or duplicates > 0
        signature_rows.append(
            {
                "Person": str(person),
                "Rater": str(rater),
                "Rows": int(len(block)),
                "UniqueContexts": int(len(set(contexts))),
                "DuplicateContexts": int(duplicates),
                "ContextSignature": repr(tuple(sorted(set(contexts)))),
            }
        )
    return pd.DataFrame(signature_rows), duplicates_found


def _has_connected_improving_switch(
    edges: set[tuple[str, str]],
    person_coordinate: Mapping[str, float] | None = None,
    rater_coordinate: Mapping[str, float] | None = None,
    *,
    direction: str = "aligned",
) -> bool:
    sign = 1.0 if direction == "aligned" else -1.0
    ordered_edges = sorted(edges)
    for (p1, r1), (p2, r2) in combinations(ordered_edges, 2):
        if p1 == p2 or r1 == r2 or (p1, r2) in edges or (p2, r1) in edges:
            continue
        if person_coordinate is not None and rater_coordinate is not None:
            delta = sign * (
                person_coordinate[p1] * rater_coordinate[r2]
                + person_coordinate[p2] * rater_coordinate[r1]
                - person_coordinate[p1] * rater_coordinate[r1]
                - person_coordinate[p2] * rater_coordinate[r2]
            )
            if delta <= 1e-15:
                continue
        candidate = set(edges)
        candidate.remove((p1, r1))
        candidate.remove((p2, r2))
        candidate.add((p1, r2))
        candidate.add((p2, r1))
        if _rater_overlap_component_count(candidate) == 1:
            return True
    return False


def evaluate_fixed_density_perturbation_feasibility(
    data: pd.DataFrame,
    *,
    person_col: str = "Person",
    rater_col: str = "Rater",
    context_cols: Iterable[str] = (),
    weight_col: str | None = "Weight",
    rater_mapping_confirmed: bool = False,
) -> dict[str, object]:
    """Gate the exact exposure-preserving v1 perturbation scope."""

    context_cols = [str(column) for column in context_cols]
    gate_rows: list[dict[str, object]] = []

    def add_gate(gate: str, passed: bool, evidence: str, action: str) -> None:
        gate_rows.append(
            {
                "Gate": gate,
                "Passed": bool(passed),
                "Evidence": evidence,
                "ActionIfFailed": action,
            }
        )

    frame_ready = isinstance(data, pd.DataFrame) and not data.empty
    add_gate(
        "Likelihood-row data",
        frame_ready,
        f"rows={len(data) if isinstance(data, pd.DataFrame) else 0}",
        "Fit the model and retain likelihood-included rows.",
    )
    if not frame_ready:
        return {
            "available": False,
            "reason": "Likelihood-row data are unavailable.",
            "gates": pd.DataFrame(gate_rows),
            "block_profiles": pd.DataFrame(),
        }

    required = [person_col, rater_col, *context_cols]
    columns_ready = all(column in data.columns for column in required)
    add_gate(
        "Required design columns",
        columns_ready,
        "required=" + ",".join(required),
        "Map Person, the confirmed Rater facet, and all non-rater context facets.",
    )
    if not columns_ready:
        missing = [column for column in required if column not in data.columns]
        return {
            "available": False,
            "reason": "Missing design columns: " + ", ".join(missing),
            "gates": pd.DataFrame(gate_rows),
            "block_profiles": pd.DataFrame(),
        }

    add_gate(
        "Confirmed Rater role",
        bool(rater_mapping_confirmed),
        f"rater_col={rater_col}; confirmed={bool(rater_mapping_confirmed)}",
        "Confirm which facet represents raters/severity-bearing agents.",
    )

    work = data[required + ([weight_col] if weight_col and weight_col in data.columns else [])].copy()
    nonblank = work[required].notna().all(axis=1)
    add_gate(
        "Nonblank design identifiers",
        bool(nonblank.all()),
        f"nonblank={int(nonblank.sum())}/{len(work)}",
        "Resolve blank Person/facet values before sensitivity analysis.",
    )
    if not nonblank.all():
        work = work.loc[nonblank].copy()
    work = _clean_string_frame(work, required).join(
        work.drop(columns=required), how="left"
    )

    weights_ready = True
    weight_evidence = "No frequency-weight column."
    if weight_col and weight_col in work.columns:
        weights = pd.to_numeric(work[weight_col], errors="coerce")
        weights_ready = bool(weights.notna().all() and np.allclose(weights.to_numpy(dtype=float), 1.0))
        weight_evidence = (
            f"unit_weights={int(np.isclose(weights.fillna(np.nan), 1.0).sum())}/{len(weights)}"
        )
    add_gate(
        "Uncompressed row topology",
        weights_ready,
        weight_evidence,
        "Expand frequency-weighted rows before constructing assignment counterfactuals.",
    )

    persons = sorted(work[person_col].unique().tolist())
    raters = sorted(work[rater_col].unique().tolist())
    level_ready = len(persons) >= 2 and len(raters) >= 2
    add_gate(
        "Multiple Persons and Raters",
        level_ready,
        f"persons={len(persons)}; raters={len(raters)}",
        "At least two Persons and two confirmed Raters are required.",
    )

    profiles, duplicate_contexts = _context_signature_table(
        work,
        person_col=person_col,
        rater_col=rater_col,
        context_cols=context_cols,
    )
    add_gate(
        "Unique row per assignment context",
        not duplicate_contexts,
        f"duplicate_context_rows={int(profiles['DuplicateContexts'].sum()) if not profiles.empty else 0}",
        "Add the omitted context/occasion facet or aggregate the intended analysis unit explicitly.",
    )
    signatures = profiles["ContextSignature"].nunique() if not profiles.empty else 0
    signature_ready = signatures == 1
    add_gate(
        "Common context signature per Person-Rater block",
        signature_ready,
        f"distinct_block_signatures={int(signatures)}",
        "Use the future generalized flow runner; v1 only moves exchangeable equal-context blocks.",
    )

    edges = set(zip(work[person_col], work[rater_col]))
    edge_limit_ready = len(edges) <= MAX_EXACT_SWITCH_EDGES and len(work) <= MAX_EXACT_SWITCH_ROWS
    add_gate(
        "Exact-switch computation envelope",
        edge_limit_ready,
        f"edges={len(edges)}/{MAX_EXACT_SWITCH_EDGES}; rows={len(work)}/{MAX_EXACT_SWITCH_ROWS}",
        "Use a sampled/flow optimization runner for a larger design.",
    )
    component_count = _rater_overlap_component_count(edges)
    connected = component_count == 1
    add_gate(
        "Connected direct Rater overlap",
        connected,
        f"components={component_count}",
        "Confirm links or redesign before preserving a connected counterfactual.",
    )
    switch_exists = _has_connected_improving_switch(edges) if connected else False
    add_gate(
        "At least one connected degree-preserving 2-switch",
        switch_exists,
        f"switch_exists={switch_exists}",
        "A complete or structurally fixed design has no within-degree assignment perturbation in v1.",
    )

    gates = pd.DataFrame(gate_rows)
    available = bool(gates["Passed"].all())
    failed = gates.loc[~gates["Passed"], "Gate"].astype(str).tolist()
    return {
        "available": available,
        "reason": (
            "Exact row-exposure-preserving assignment perturbation is ready."
            if available
            else "Blocked by: " + "; ".join(failed)
        ),
        "gates": gates,
        "block_profiles": profiles,
        "edges": edges,
        "persons": persons,
        "raters": raters,
        "context_cols": context_cols,
        "schema_version": SENSITIVITY_SCHEMA_VERSION,
    }


def _invariant_table(
    original: pd.DataFrame,
    counterfactual: pd.DataFrame,
    *,
    person_col: str,
    rater_col: str,
    context_cols: list[str],
) -> pd.DataFrame:
    original_edges = original[[person_col, rater_col]].drop_duplicates()
    counter_edges = counterfactual[[person_col, rater_col]].drop_duplicates()

    def series_equal(left: pd.Series, right: pd.Series) -> bool:
        return left.sort_index().astype(int).equals(right.sort_index().astype(int))

    person_degree_before = original_edges.groupby(person_col, observed=True)[rater_col].nunique()
    person_degree_after = counter_edges.groupby(person_col, observed=True)[rater_col].nunique()
    rater_person_before = original_edges.groupby(rater_col, observed=True)[person_col].nunique()
    rater_person_after = counter_edges.groupby(rater_col, observed=True)[person_col].nunique()
    rater_rows_before = original.groupby(rater_col, observed=True).size()
    rater_rows_after = counterfactual.groupby(rater_col, observed=True).size()
    cell_cols = [person_col, rater_col, *context_cols]
    duplicates_after = int(counterfactual.duplicated(cell_cols).sum())
    after_edges = set(map(tuple, counter_edges[[person_col, rater_col]].astype(str).itertuples(index=False, name=None)))
    changed_edges = len(
        set(map(tuple, original_edges[[person_col, rater_col]].astype(str).itertuples(index=False, name=None)))
        ^ after_edges
    ) // 2
    rows = [
        ("Row count", len(original) == len(counterfactual), f"{len(original)} -> {len(counterfactual)}"),
        ("Per-Person Rater degree", series_equal(person_degree_before, person_degree_after), "exact equality"),
        ("Per-Rater Person exposure", series_equal(rater_person_before, rater_person_after), "exact equality"),
        ("Per-Rater response-row exposure", series_equal(rater_rows_before, rater_rows_after), "exact equality"),
        ("Unique Person-Rater-context rows", duplicates_after == 0, f"duplicates={duplicates_after}"),
        ("Connected direct Rater overlap", _rater_overlap_component_count(after_edges) == 1, f"components={_rater_overlap_component_count(after_edges)}"),
        ("Assignment changed", changed_edges > 0, f"changed_edges={changed_edges}"),
    ]
    return pd.DataFrame(rows, columns=["Invariant", "Passed", "Evidence"])


def build_degree_preserving_assignment_perturbation(
    data: pd.DataFrame,
    *,
    person_scores: Mapping[str, float],
    rater_scores: Mapping[str, float],
    person_col: str = "Person",
    rater_col: str = "Rater",
    facet_cols: Iterable[str] = (),
    context_cols: Iterable[str] = (),
    weight_col: str | None = "Weight",
    rater_mapping_confirmed: bool = False,
    direction: str = "aligned",
    max_switches: int | None = None,
) -> dict[str, object]:
    """Create a deterministic connected 2-switch path and score-free design."""

    direction = str(direction).lower()
    if direction not in {"aligned", "anti_aligned"}:
        raise ValueError("direction must be 'aligned' or 'anti_aligned'.")
    context_cols = [str(column) for column in context_cols]
    facet_cols = [str(column) for column in facet_cols]
    preflight = evaluate_fixed_density_perturbation_feasibility(
        data,
        person_col=person_col,
        rater_col=rater_col,
        context_cols=context_cols,
        weight_col=weight_col,
        rater_mapping_confirmed=rater_mapping_confirmed,
    )
    if not preflight.get("available"):
        return {
            **preflight,
            "design": pd.DataFrame(),
            "assignment_map": pd.DataFrame(),
            "switch_ledger": pd.DataFrame(),
            "trajectory": pd.DataFrame(),
            "invariants": pd.DataFrame(),
        }

    work_cols = [person_col, rater_col, *facet_cols, *context_cols]
    if weight_col and weight_col in data.columns:
        work_cols.append(weight_col)
    work_cols = list(dict.fromkeys(work_cols))
    source = data[work_cols].copy().reset_index(drop=True)
    source["_SourceRow"] = np.arange(len(source), dtype=int)
    for column in dict.fromkeys([person_col, rater_col, *facet_cols, *context_cols]):
        source[column] = source[column].astype(str)

    persons = list(preflight["persons"])
    raters = list(preflight["raters"])
    person_coordinate = _rank_coordinate(
        {str(key): float(value) for key, value in person_scores.items()}, persons
    )
    rater_coordinate = _rank_coordinate(
        {str(key): float(value) for key, value in rater_scores.items()}, raters
    )
    sign = 1.0 if direction == "aligned" else -1.0

    edges: set[tuple[str, str]] = set(preflight["edges"])
    block_at = {(person, rater): rater for person, rater in edges}
    objective = _assignment_objective(edges, person_coordinate, rater_coordinate)
    trajectory_rows = [
        {
            "Switch": 0,
            "Objective": objective,
            "DirectionAdjustedObjective": sign * objective,
            "AssignmentRankCorrelation": _assignment_correlation(
                edges, person_coordinate, rater_coordinate
            ),
            "ChangedBlocks": 0,
        }
    ]
    ledger_rows: list[dict[str, object]] = []
    switch_limit = min(
        max(1, int(max_switches if max_switches is not None else len(edges) * 4)),
        2_000,
    )
    for switch_index in range(1, switch_limit + 1):
        best: tuple[float, str, str, str, str, set[tuple[str, str]]] | None = None
        for (p1, r1), (p2, r2) in combinations(sorted(edges), 2):
            if p1 == p2 or r1 == r2 or (p1, r2) in edges or (p2, r1) in edges:
                continue
            raw_delta = (
                person_coordinate[p1] * rater_coordinate[r2]
                + person_coordinate[p2] * rater_coordinate[r1]
                - person_coordinate[p1] * rater_coordinate[r1]
                - person_coordinate[p2] * rater_coordinate[r2]
            ) / len(edges)
            adjusted_delta = sign * raw_delta
            if adjusted_delta <= 1e-15:
                continue
            candidate = set(edges)
            candidate.remove((p1, r1))
            candidate.remove((p2, r2))
            candidate.add((p1, r2))
            candidate.add((p2, r1))
            if _rater_overlap_component_count(candidate) != 1:
                continue
            key = (adjusted_delta, p1, r1, p2, r2, candidate)
            if best is None or key[:5] > best[:5]:
                best = key
        if best is None:
            break
        adjusted_delta, p1, r1, p2, r2, candidate = best
        source_rater_1 = block_at.pop((p1, r1))
        source_rater_2 = block_at.pop((p2, r2))
        block_at[(p1, r2)] = source_rater_1
        block_at[(p2, r1)] = source_rater_2
        before = objective
        edges = candidate
        objective = _assignment_objective(edges, person_coordinate, rater_coordinate)
        ledger_rows.append(
            {
                "Switch": switch_index,
                "PersonA": p1,
                "PersonB": p2,
                "RaterAFrom": r1,
                "RaterATo": r2,
                "RaterBFrom": r2,
                "RaterBTo": r1,
                "ObjectiveBefore": before,
                "ObjectiveAfter": objective,
                "DirectionAdjustedGain": adjusted_delta,
                "OverlapComponentsAfter": _rater_overlap_component_count(edges),
            }
        )
        mapping_now = {
            (person, source_rater): current_rater
            for (person, current_rater), source_rater in block_at.items()
        }
        changed = sum(current != original for (_, original), current in mapping_now.items())
        trajectory_rows.append(
            {
                "Switch": switch_index,
                "Objective": objective,
                "DirectionAdjustedObjective": sign * objective,
                "AssignmentRankCorrelation": _assignment_correlation(
                    edges, person_coordinate, rater_coordinate
                ),
                "ChangedBlocks": int(changed),
            }
        )

    mapping = {
        (person, source_rater): current_rater
        for (person, current_rater), source_rater in block_at.items()
    }
    assignment_map = pd.DataFrame(
        [
            {
                "Person": person,
                "SourceRater": source_rater,
                "CounterfactualRater": current_rater,
                "Changed": bool(source_rater != current_rater),
            }
            for (person, source_rater), current_rater in sorted(mapping.items())
        ]
    )
    counterfactual = source.copy()
    counterfactual[rater_col] = [
        mapping[(str(person), str(rater))]
        for person, rater in zip(source[person_col], source[rater_col])
    ]
    invariants = _invariant_table(
        source,
        counterfactual,
        person_col=person_col,
        rater_col=rater_col,
        context_cols=context_cols,
    )
    available = bool(
        ledger_rows
        and not invariants.empty
        and invariants["Passed"].all()
    )
    return {
        "available": available,
        "reason": (
            "A score-free fixed-density counterfactual design was constructed."
            if available
            else "No improving invariant-preserving switch path was constructed."
        ),
        "schema_version": SENSITIVITY_SCHEMA_VERSION,
        "direction": direction,
        "design": counterfactual,
        "assignment_map": assignment_map,
        "switch_ledger": pd.DataFrame(ledger_rows),
        "trajectory": pd.DataFrame(trajectory_rows),
        "invariants": invariants,
        "gates": preflight["gates"],
        "block_profiles": preflight["block_profiles"],
        "person_coordinate": pd.DataFrame(
            {"Person": list(person_coordinate), "RankCoordinate": list(person_coordinate.values())}
        ),
        "rater_coordinate": pd.DataFrame(
            {"Rater": list(rater_coordinate), "SeverityRankCoordinate": list(rater_coordinate.values())}
        ),
        "claim_boundary": (
            "This is a deterministic assignment counterfactual. It contains no observed Score and does not "
            "identify the actual assignment mechanism. Scores must be simulated from an explicit fitted model."
        ),
    }


def build_perturbation_path_snapshots(
    data: pd.DataFrame,
    perturbation: Mapping[str, object],
    *,
    requested_doses: Iterable[float] = (0.0, 1.0),
    person_col: str = "Person",
    rater_col: str = "Rater",
    facet_cols: Iterable[str] = (),
    context_cols: Iterable[str] = (),
    weight_col: str | None = "Weight",
) -> dict[str, object]:
    """Reconstruct score-free snapshots along a validated greedy switch path.

    Dose is normalized progress in the direction-adjusted assignment objective,
    not the fraction of switches and not a claim about the true assignment
    mechanism. Requested doses that map to the same discrete switch are fitted
    only once and remain visible in ``dose_table``.
    """

    if not isinstance(perturbation, Mapping) or not perturbation.get("available"):
        return {
            "available": False,
            "reason": "An available assignment perturbation is required.",
            "designs": {},
            "dose_table": pd.DataFrame(),
            "path_invariants": pd.DataFrame(),
            "path_assignment_map": pd.DataFrame(),
        }
    trajectory = perturbation.get("trajectory", pd.DataFrame())
    ledger = perturbation.get("switch_ledger", pd.DataFrame())
    if not isinstance(trajectory, pd.DataFrame) or trajectory.empty:
        raise ValueError("Perturbation trajectory is unavailable.")
    if not isinstance(ledger, pd.DataFrame) or ledger.empty:
        raise ValueError("Perturbation switch ledger is unavailable.")

    doses: list[float] = []
    for value in requested_doses:
        dose = float(value)
        if not np.isfinite(dose) or dose < 0.0 or dose > 1.0:
            raise ValueError("requested_doses must contain finite values in [0, 1].")
        if dose not in doses:
            doses.append(dose)
    doses = sorted(set([0.0, 1.0, *doses]))

    trajectory = trajectory.copy().sort_values("Switch").reset_index(drop=True)
    adjusted = pd.to_numeric(
        trajectory["DirectionAdjustedObjective"], errors="coerce"
    ).to_numpy(dtype=float)
    if not np.isfinite(adjusted).all() or len(adjusted) < 2:
        raise ValueError("Perturbation trajectory objective is incomplete.")
    objective_gain = float(adjusted[-1] - adjusted[0])
    if objective_gain <= 0:
        raise ValueError("Perturbation trajectory has no positive alignment gain.")
    achieved = np.clip((adjusted - adjusted[0]) / objective_gain, 0.0, 1.0)
    trajectory["AchievedAlignmentDose"] = achieved

    selected_by_switch: dict[int, list[float]] = {}
    for dose in doses:
        candidates = np.flatnonzero(achieved >= dose - 1e-12)
        row_index = int(candidates[0]) if candidates.size else len(trajectory) - 1
        switch = int(trajectory.iloc[row_index]["Switch"])
        selected_by_switch.setdefault(switch, []).append(float(dose))

    facet_cols = [str(column) for column in facet_cols]
    context_cols = [str(column) for column in context_cols]
    work_cols = list(dict.fromkeys([person_col, rater_col, *facet_cols, *context_cols]))
    if weight_col and weight_col in data.columns:
        work_cols.append(weight_col)
    missing = [column for column in work_cols if column not in data.columns]
    if missing:
        raise ValueError("Snapshot source is missing columns: " + ", ".join(missing))
    source = data[work_cols].copy().reset_index(drop=True)
    source[person_col] = source[person_col].astype(str)
    source[rater_col] = source[rater_col].astype(str)
    for column in dict.fromkeys([*facet_cols, *context_cols]):
        source[column] = source[column].astype(str)
    source["_SourceRow"] = np.arange(len(source), dtype=int)
    current = source.copy()

    designs: dict[str, pd.DataFrame] = {}
    dose_rows: list[dict[str, object]] = []
    invariant_parts: list[pd.DataFrame] = []
    map_parts: list[pd.DataFrame] = []
    target_switches = sorted(selected_by_switch)
    ledger_by_switch = ledger.set_index("Switch", drop=False)
    last_applied = 0
    total_switches = int(trajectory["Switch"].max())

    for switch in target_switches:
        for switch_index in range(last_applied + 1, switch + 1):
            if switch_index not in ledger_by_switch.index:
                raise ValueError(f"Switch {switch_index} is missing from the perturbation ledger.")
            step = ledger_by_switch.loc[switch_index]
            if isinstance(step, pd.DataFrame):
                raise ValueError(f"Switch {switch_index} is duplicated in the perturbation ledger.")
            p1 = str(step["PersonA"])
            p2 = str(step["PersonB"])
            r1 = str(step["RaterAFrom"])
            r2 = str(step["RaterBFrom"])
            mask1 = current[person_col].eq(p1) & current[rater_col].eq(r1)
            mask2 = current[person_col].eq(p2) & current[rater_col].eq(r2)
            if not mask1.any() or not mask2.any() or (mask1 & mask2).any():
                raise RuntimeError(
                    f"Switch {switch_index} cannot be replayed on the current assignment."
                )
            current.loc[mask1, rater_col] = r2
            current.loc[mask2, rater_col] = r1
        last_applied = switch

        trajectory_row = trajectory.loc[trajectory["Switch"].eq(switch)].iloc[0]
        dose = float(trajectory_row["AchievedAlignmentDose"])
        scenario = (
            "observed_assignment"
            if switch == 0
            else (
                "aligned_counterfactual"
                if switch == total_switches
                else f"alignment_switch_{switch:04d}"
            )
        )
        snapshot = current.copy()
        designs[scenario] = snapshot
        requested = selected_by_switch[switch]
        dose_rows.append({
            "Scenario": scenario,
            "RequestedAlignmentDoses": "|".join(f"{value:.6f}" for value in requested),
            "RequestedDoseCount": len(requested),
            "SwitchesApplied": switch,
            "TotalSwitches": total_switches,
            "AchievedAlignmentDose": dose,
            "AssignmentRankCorrelation": float(trajectory_row["AssignmentRankCorrelation"]),
            "Objective": float(trajectory_row["Objective"]),
            "ChangedBlocks": int(trajectory_row["ChangedBlocks"]),
            "IsObserved": switch == 0,
            "IsEndpoint": switch == total_switches,
        })

        invariants = _invariant_table(
            source,
            snapshot,
            person_col=person_col,
            rater_col=rater_col,
            context_cols=context_cols,
        )
        invariants.insert(0, "AchievedAlignmentDose", dose)
        invariants.insert(0, "Scenario", scenario)
        invariants["RequiredAtDose"] = ~invariants["Invariant"].eq("Assignment changed") | (switch > 0)
        invariants["PassedForDose"] = invariants["Passed"] | ~invariants["RequiredAtDose"]
        invariant_parts.append(invariants)

        source_raters = source[["_SourceRow", person_col, rater_col]].rename(
            columns={rater_col: "SourceRater"}
        )
        current_raters = snapshot[["_SourceRow", rater_col]].rename(
            columns={rater_col: "CounterfactualRater"}
        )
        assignment_map = source_raters.merge(
            current_raters, on="_SourceRow", how="left", validate="one_to_one"
        ).drop(columns=["_SourceRow"])
        assignment_map = assignment_map.drop_duplicates(
            [person_col, "SourceRater", "CounterfactualRater"]
        ).reset_index(drop=True)
        assignment_map.insert(0, "AchievedAlignmentDose", dose)
        assignment_map.insert(0, "Scenario", scenario)
        assignment_map["Changed"] = assignment_map["SourceRater"].ne(
            assignment_map["CounterfactualRater"]
        )
        map_parts.append(assignment_map)

    path_invariants = pd.concat(invariant_parts, ignore_index=True)
    available = bool(
        designs
        and path_invariants.loc[path_invariants["RequiredAtDose"], "PassedForDose"].all()
    )
    return {
        "available": available,
        "reason": (
            "Score-free invariant-preserving assignment-dose snapshots were constructed."
            if available
            else "At least one required path invariant failed."
        ),
        "schema_version": SENSITIVITY_SCHEMA_VERSION,
        "designs": designs,
        "dose_table": pd.DataFrame(dose_rows),
        "path_invariants": path_invariants,
        "path_assignment_map": pd.concat(map_parts, ignore_index=True),
        "trajectory": trajectory,
        "claim_boundary": (
            "Alignment dose is normalized progress along one deterministic greedy switch path. "
            "It is neither an estimated assignment propensity nor a global worst-case scale."
        ),
    }
