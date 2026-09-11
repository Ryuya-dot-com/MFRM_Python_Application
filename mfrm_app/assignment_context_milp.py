"""Score-free assignment perturbations with exact Rater-by-context margins.

The strict fixed-density runner in :mod:`mfrm_app.assignment_sensitivity`
uses degree-preserving 2-switches and therefore requires exchangeable
Person-Rater blocks.  This module covers the harder case in which blocks have
different non-rater context profiles.

Each original Person-Rater block is assigned as an indivisible unit by a
binary MILP.  The model preserves Person degree, per-Rater Person exposure,
and every per-Rater context-cell row count.  Direct Rater-overlap connectivity
is guaranteed by locking one observed witness Person for every edge of a
deterministic spanning tree.  Observed outcomes are never copied into the
returned design.
"""

from __future__ import annotations

from collections import Counter
from itertools import combinations
from typing import Iterable, Mapping

import numpy as np
import pandas as pd
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix


CONTEXT_MILP_SCHEMA_VERSION = "context_margin_assignment_milp_v1"
MAX_MILP_VARIABLES = 10_000
MAX_MILP_CONSTRAINTS = 50_000
MAX_MILP_NONZEROS = 1_000_000


def _rank_coordinate(values: Mapping[str, float], levels: list[str]) -> dict[str, float]:
    series = pd.Series({level: values.get(level, np.nan) for level in levels}, dtype=float)
    finite = np.isfinite(series.to_numpy(dtype=float))
    if not finite.all():
        missing = series.index[~finite].astype(str).tolist()
        raise ValueError("Non-finite ordering coordinate for: " + ", ".join(missing[:10]))
    if series.nunique() < 2:
        raise ValueError("Ordering coordinates must contain at least two distinct values.")
    ranks = series.rank(method="average", pct=True)
    return {level: float(ranks.loc[level] - 0.5) for level in levels}


def _overlap_components(edges: set[tuple[str, str]]) -> int:
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


def _spanning_tree_witnesses(edges: set[tuple[str, str]]) -> pd.DataFrame:
    raters = sorted({rater for _, rater in edges})
    people_by_rater = {
        rater: {person for person, candidate in edges if candidate == rater}
        for rater in raters
    }
    parent = {rater: rater for rater in raters}

    def find(value: str) -> str:
        while parent[value] != value:
            parent[value] = parent[parent[value]]
            value = parent[value]
        return value

    def union(left: str, right: str) -> bool:
        root_left, root_right = find(left), find(right)
        if root_left == root_right:
            return False
        parent[root_right] = root_left
        return True

    rows: list[dict[str, object]] = []
    for left, right in combinations(raters, 2):
        shared = sorted(people_by_rater[left] & people_by_rater[right])
        if shared and union(left, right):
            rows.append(
                {
                    "TreeEdge": len(rows) + 1,
                    "RaterA": left,
                    "RaterB": right,
                    "WitnessPerson": shared[0],
                    "ObservedSharedPersons": len(shared),
                }
            )
            if len(rows) == max(0, len(raters) - 1):
                break
    return pd.DataFrame(rows)


def _context_key_rows(
    frame: pd.DataFrame, context_cols: list[str]
) -> list[tuple[str, ...]]:
    if not context_cols:
        return [("__single_context__",)] * len(frame)
    return [
        tuple(row)
        for row in frame[context_cols].astype(str).itertuples(index=False, name=None)
    ]


def evaluate_context_margin_milp_feasibility(
    data: pd.DataFrame,
    *,
    person_col: str = "Person",
    rater_col: str = "Rater",
    context_cols: Iterable[str] = (),
    weight_col: str | None = "Weight",
    rater_mapping_confirmed: bool = False,
) -> dict[str, object]:
    """Audit whether an exact context-margin MILP can be constructed."""

    context_cols = [str(column) for column in context_cols]
    gate_rows: list[dict[str, object]] = []

    def gate(name: str, passed: bool, evidence: str, action: str) -> None:
        gate_rows.append(
            {
                "Gate": name,
                "Passed": bool(passed),
                "Evidence": evidence,
                "ActionIfFailed": action,
            }
        )

    ready = isinstance(data, pd.DataFrame) and not data.empty
    gate(
        "Likelihood-row data",
        ready,
        f"rows={len(data) if isinstance(data, pd.DataFrame) else 0}",
        "Fit the model and retain likelihood-included rows.",
    )
    if not ready:
        return {
            "available": False,
            "reason": "Likelihood-row data are unavailable.",
            "gates": pd.DataFrame(gate_rows),
            "block_profiles": pd.DataFrame(),
            "witnesses": pd.DataFrame(),
        }

    required = list(dict.fromkeys([person_col, rater_col, *context_cols]))
    columns_ready = all(column in data.columns for column in required)
    gate(
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
            "witnesses": pd.DataFrame(),
        }

    gate(
        "Confirmed Rater role",
        bool(rater_mapping_confirmed),
        f"rater_col={rater_col}; confirmed={bool(rater_mapping_confirmed)}",
        "Confirm which facet represents raters/severity-bearing agents.",
    )
    work_cols = required + (
        [weight_col] if weight_col and weight_col in data.columns else []
    )
    work = data[work_cols].copy()
    nonblank = work[required].notna().all(axis=1)
    gate(
        "Nonblank design identifiers",
        bool(nonblank.all()),
        f"nonblank={int(nonblank.sum())}/{len(work)}",
        "Resolve blank Person/facet values before sensitivity analysis.",
    )
    weights_ready = True
    weight_evidence = "No frequency-weight column."
    if weight_col and weight_col in work.columns:
        weights = pd.to_numeric(work[weight_col], errors="coerce")
        weights_ready = bool(
            weights.notna().all()
            and np.allclose(weights.to_numpy(dtype=float), 1.0)
        )
        weight_evidence = (
            f"unit_weights={int(np.isclose(weights.fillna(np.nan), 1.0).sum())}/{len(weights)}"
        )
    gate(
        "Uncompressed row topology",
        weights_ready,
        weight_evidence,
        "Expand frequency-weighted rows before constructing assignment counterfactuals.",
    )
    if not nonblank.all():
        gates = pd.DataFrame(gate_rows)
        return {
            "available": False,
            "reason": "Blocked by: Nonblank design identifiers",
            "gates": gates,
            "block_profiles": pd.DataFrame(),
            "witnesses": pd.DataFrame(),
        }

    for column in required:
        work[column] = work[column].astype(str)
    persons = sorted(work[person_col].unique().tolist())
    raters = sorted(work[rater_col].unique().tolist())
    level_ready = len(persons) >= 2 and len(raters) >= 2
    gate(
        "Multiple Persons and Raters",
        level_ready,
        f"persons={len(persons)}; raters={len(raters)}",
        "At least two Persons and two confirmed Raters are required.",
    )

    profile_rows: list[dict[str, object]] = []
    for (person, rater), block in work.groupby(
        [person_col, rater_col], observed=True, sort=True
    ):
        context_counts = Counter(_context_key_rows(block, context_cols))
        profile_rows.append(
            {
                "Person": str(person),
                "Rater": str(rater),
                "Rows": int(len(block)),
                "UniqueContexts": int(len(context_counts)),
                "DuplicateContextRows": int(len(block) - len(context_counts)),
                "ContextCountSignature": repr(tuple(sorted(context_counts.items()))),
            }
        )
    profiles = pd.DataFrame(profile_rows)
    block_count = len(profiles)
    variable_count = block_count * len(raters)
    context_count = len(set(_context_key_rows(work, context_cols)))
    witness_upper = max(0, len(raters) - 1) * 2
    constraint_estimate = (
        block_count
        + len(persons) * len(raters)
        + len(raters)
        + len(raters) * context_count
        + witness_upper
    )
    # Every variable enters the block, Person-Rater, Rater-count, and at least
    # one context constraint. Blocks with multiple contexts enter more rows.
    context_entries = int(profiles["UniqueContexts"].sum()) * len(raters)
    nonzero_estimate = variable_count * 3 + context_entries + witness_upper * block_count
    envelope_ready = bool(
        variable_count <= MAX_MILP_VARIABLES
        and constraint_estimate <= MAX_MILP_CONSTRAINTS
        and nonzero_estimate <= MAX_MILP_NONZEROS
    )
    gate(
        "MILP computation envelope",
        envelope_ready,
        (
            f"variables={variable_count}/{MAX_MILP_VARIABLES}; "
            f"constraints_estimate={constraint_estimate}/{MAX_MILP_CONSTRAINTS}; "
            f"nonzeros_estimate={nonzero_estimate}/{MAX_MILP_NONZEROS}"
        ),
        "Use an offline/sharded optimization study or reduce the fitted design.",
    )

    edges = set(
        map(
            tuple,
            work[[person_col, rater_col]].drop_duplicates().itertuples(
                index=False, name=None
            ),
        )
    )
    components = _overlap_components(edges)
    connected = components == 1
    gate(
        "Connected direct Rater overlap",
        connected,
        f"components={components}",
        "Confirm links or redesign before preserving a connected counterfactual.",
    )
    witnesses = _spanning_tree_witnesses(edges) if connected else pd.DataFrame()
    witnesses_ready = bool(
        connected and len(witnesses) == max(0, len(raters) - 1)
    )
    gate(
        "Connectivity witness spanning tree",
        witnesses_ready,
        f"witness_edges={len(witnesses)}/{max(0, len(raters) - 1)}",
        "A direct-overlap spanning tree could not be certified.",
    )

    gates = pd.DataFrame(gate_rows)
    available = bool(not gates.empty and gates["Passed"].all())
    failed = gates.loc[~gates["Passed"], "Gate"].astype(str).tolist()
    return {
        "available": available,
        "reason": (
            "Exact Rater-by-context margin optimization is constructible."
            if available
            else "Blocked by: " + "; ".join(failed)
        ),
        "gates": gates,
        "block_profiles": profiles,
        "witnesses": witnesses,
        "persons": persons,
        "raters": raters,
        "edges": edges,
        "context_cols": context_cols,
        "schema_version": CONTEXT_MILP_SCHEMA_VERSION,
    }


def _context_margin_audit(
    original: pd.DataFrame,
    counterfactual: pd.DataFrame,
    *,
    rater_col: str,
    context_cols: list[str],
) -> pd.DataFrame:
    group_cols = [rater_col, *context_cols]
    if context_cols:
        before = original.groupby(group_cols, observed=True).size().rename("RowsBefore")
        after = counterfactual.groupby(group_cols, observed=True).size().rename("RowsAfter")
    else:
        before = original.groupby([rater_col], observed=True).size().rename("RowsBefore")
        after = counterfactual.groupby([rater_col], observed=True).size().rename("RowsAfter")
    audit = pd.concat([before, after], axis=1).fillna(0).reset_index()
    audit["RowsBefore"] = audit["RowsBefore"].astype(int)
    audit["RowsAfter"] = audit["RowsAfter"].astype(int)
    audit["Difference"] = audit["RowsAfter"] - audit["RowsBefore"]
    audit["ExactMatch"] = audit["Difference"].eq(0)
    return audit


def _series_exact(left: pd.Series, right: pd.Series) -> bool:
    return left.sort_index().astype(int).equals(right.sort_index().astype(int))


def build_context_margin_assignment_perturbation(
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
    time_limit_seconds: float = 30.0,
) -> dict[str, object]:
    """Optimize a score-free endpoint while preserving exact design margins."""

    direction = str(direction).lower()
    if direction not in {"aligned", "anti_aligned"}:
        raise ValueError("direction must be 'aligned' or 'anti_aligned'.")
    context_cols = [str(column) for column in context_cols]
    facet_cols = [str(column) for column in facet_cols]
    preflight = evaluate_context_margin_milp_feasibility(
        data,
        person_col=person_col,
        rater_col=rater_col,
        context_cols=context_cols,
        weight_col=weight_col,
        rater_mapping_confirmed=rater_mapping_confirmed,
    )
    empty = {
        "design": pd.DataFrame(),
        "assignment_map": pd.DataFrame(),
        "trajectory": pd.DataFrame(),
        "invariants": pd.DataFrame(),
        "context_margin_audit": pd.DataFrame(),
        "solver_audit": pd.DataFrame(),
    }
    if not preflight.get("available"):
        return {**preflight, **empty}

    work_cols = list(dict.fromkeys([person_col, rater_col, *facet_cols, *context_cols]))
    if weight_col and weight_col in data.columns:
        work_cols.append(weight_col)
    source = data[work_cols].copy().reset_index(drop=True)
    string_cols = list(dict.fromkeys([person_col, rater_col, *facet_cols, *context_cols]))
    for column in string_cols:
        source[column] = source[column].astype(str)
    source["_SourceRow"] = np.arange(len(source), dtype=int)

    persons = list(preflight["persons"])
    raters = list(preflight["raters"])
    person_coordinate = _rank_coordinate(
        {str(key): float(value) for key, value in person_scores.items()}, persons
    )
    rater_coordinate = _rank_coordinate(
        {str(key): float(value) for key, value in rater_scores.items()}, raters
    )
    sign = 1.0 if direction == "aligned" else -1.0

    blocks: list[dict[str, object]] = []
    for block_id, ((person, source_rater), block) in enumerate(
        source.groupby([person_col, rater_col], observed=True, sort=True)
    ):
        counts = Counter(_context_key_rows(block, context_cols))
        blocks.append(
            {
                "Block": int(block_id),
                "Person": str(person),
                "SourceRater": str(source_rater),
                "Rows": int(len(block)),
                "ContextCounts": counts,
                "ContextCountSignature": repr(tuple(sorted(counts.items()))),
            }
        )
    block_count = len(blocks)
    rater_index = {rater: index for index, rater in enumerate(raters)}

    def variable(block_index: int, target_rater: str) -> int:
        return block_index * len(raters) + rater_index[target_rater]

    row_indices: list[int] = []
    col_indices: list[int] = []
    coefficients: list[float] = []
    lower: list[float] = []
    upper: list[float] = []
    constraint_names: list[str] = []

    def add_constraint(
        entries: Iterable[tuple[int, float]],
        lb: float,
        ub: float,
        name: str,
    ) -> None:
        row = len(lower)
        for column, value in entries:
            if value:
                row_indices.append(row)
                col_indices.append(int(column))
                coefficients.append(float(value))
        lower.append(float(lb))
        upper.append(float(ub))
        constraint_names.append(name)

    # Every source block is moved exactly once.
    for block_index in range(block_count):
        add_constraint(
            ((variable(block_index, rater), 1.0) for rater in raters),
            1.0,
            1.0,
            f"source_block_once:{block_index}",
        )

    blocks_by_person = {
        person: [index for index, block in enumerate(blocks) if block["Person"] == person]
        for person in persons
    }
    # A Person cannot contribute two source blocks to the same target Rater.
    for person in persons:
        for rater in raters:
            add_constraint(
                ((variable(index, rater), 1.0) for index in blocks_by_person[person]),
                0.0,
                1.0,
                f"unique_person_target:{person}:{rater}",
            )

    source_blocks_per_rater = Counter(block["SourceRater"] for block in blocks)
    for rater in raters:
        target = float(source_blocks_per_rater[rater])
        add_constraint(
            ((variable(index, rater), 1.0) for index in range(block_count)),
            target,
            target,
            f"rater_person_exposure:{rater}",
        )

    all_contexts = sorted(
        {context for block in blocks for context in block["ContextCounts"]}
    )
    original_context_targets: Counter[tuple[str, tuple[str, ...]]] = Counter()
    for block in blocks:
        for context, count in block["ContextCounts"].items():
            original_context_targets[(str(block["SourceRater"]), context)] += int(count)
    for rater in raters:
        for context in all_contexts:
            target = float(original_context_targets[(rater, context)])
            add_constraint(
                (
                    (variable(index, rater), float(block["ContextCounts"].get(context, 0)))
                    for index, block in enumerate(blocks)
                ),
                target,
                target,
                f"rater_context_rows:{rater}:{context!r}",
            )

    witnesses = preflight["witnesses"].copy()
    for witness in witnesses.itertuples(index=False):
        for rater in (str(witness.RaterA), str(witness.RaterB)):
            person = str(witness.WitnessPerson)
            add_constraint(
                ((variable(index, rater), 1.0) for index in blocks_by_person[person]),
                1.0,
                1.0,
                f"connectivity_witness:{person}:{rater}",
            )

    variable_count = block_count * len(raters)
    constraint_count = len(lower)
    nonzero_count = len(coefficients)
    if (
        variable_count > MAX_MILP_VARIABLES
        or constraint_count > MAX_MILP_CONSTRAINTS
        or nonzero_count > MAX_MILP_NONZEROS
    ):
        reason = (
            "Constructed MILP exceeded its fail-closed envelope: "
            f"variables={variable_count}, constraints={constraint_count}, nonzeros={nonzero_count}."
        )
        return {**preflight, **empty, "available": False, "reason": reason}

    matrix = coo_matrix(
        (coefficients, (row_indices, col_indices)),
        shape=(constraint_count, variable_count),
        dtype=float,
    ).tocsr()
    lower_array = np.asarray(lower, dtype=float)
    upper_array = np.asarray(upper, dtype=float)
    observed = np.zeros(variable_count, dtype=float)
    for index, block in enumerate(blocks):
        observed[variable(index, str(block["SourceRater"]))] = 1.0
    observed_lhs = np.asarray(matrix @ observed, dtype=float)
    observed_violation = np.maximum(
        np.maximum(lower_array - observed_lhs, 0.0),
        np.maximum(observed_lhs - upper_array, 0.0),
    )
    observed_feasible = bool(np.max(observed_violation, initial=0.0) <= 1e-9)
    if not observed_feasible:
        worst = int(np.argmax(observed_violation))
        reason = (
            "The observed assignment did not reproduce the constructed constraints; "
            f"worst={constraint_names[worst]}, violation={observed_violation[worst]:.3g}."
        )
        solver_audit = pd.DataFrame(
            [{
                "ObservedBaselineFeasible": False,
                "ObservedMaxConstraintViolation": float(np.max(observed_violation)),
                "SolverCalled": False,
                "Variables": variable_count,
                "Constraints": constraint_count,
                "Nonzeros": nonzero_count,
            }]
        )
        return {
            **preflight,
            **empty,
            "available": False,
            "reason": reason,
            "solver_audit": solver_audit,
        }

    objective = np.empty(variable_count, dtype=float)
    for index, block in enumerate(blocks):
        person = str(block["Person"])
        for rater in raters:
            # scipy minimizes. Scaling keeps the reported function on the same
            # mean-product scale as the strict runner.
            objective[variable(index, rater)] = (
                -sign * person_coordinate[person] * rater_coordinate[rater] / block_count
            )
    solution = milp(
        c=objective,
        integrality=np.ones(variable_count, dtype=np.int8),
        bounds=Bounds(np.zeros(variable_count), np.ones(variable_count)),
        constraints=LinearConstraint(matrix, lower_array, upper_array),
        options={
            "disp": False,
            "presolve": True,
            "time_limit": max(1.0, float(time_limit_seconds)),
            "mip_rel_gap": 0.0,
        },
    )
    optimal = bool(solution.status == 0 and solution.success and solution.x is not None)
    solver_audit = pd.DataFrame(
        [{
            "ObservedBaselineFeasible": True,
            "ObservedMaxConstraintViolation": float(np.max(observed_violation, initial=0.0)),
            "SolverCalled": True,
            "Solver": "scipy.optimize.milp (HiGHS)",
            "SolverStatus": int(solution.status),
            "SolverSuccess": bool(solution.success),
            "SolverMessage": str(solution.message),
            "GlobalOptimumCertified": optimal,
            "MIPGap": float(getattr(solution, "mip_gap", np.nan)),
            "MIPNodeCount": float(getattr(solution, "mip_node_count", np.nan)),
            "Variables": variable_count,
            "Constraints": constraint_count,
            "Nonzeros": nonzero_count,
            "TimeLimitSeconds": max(1.0, float(time_limit_seconds)),
        }]
    )
    if not optimal:
        return {
            **preflight,
            **empty,
            "available": False,
            "reason": "MILP did not certify a global optimum: " + str(solution.message),
            "solver_audit": solver_audit,
        }

    rounded = np.rint(np.asarray(solution.x, dtype=float))
    solution_lhs = np.asarray(matrix @ rounded, dtype=float)
    solution_violation = np.maximum(
        np.maximum(lower_array - solution_lhs, 0.0),
        np.maximum(solution_lhs - upper_array, 0.0),
    )
    maximum_solution_violation = float(np.max(solution_violation, initial=0.0))
    integral_residual = float(np.max(np.abs(np.asarray(solution.x) - rounded), initial=0.0))
    solver_audit["SolutionMaxConstraintViolation"] = maximum_solution_violation
    solver_audit["SolutionMaxIntegralityResidual"] = integral_residual
    numerical_ready = bool(maximum_solution_violation <= 1e-7 and integral_residual <= 1e-7)
    if not numerical_ready:
        return {
            **preflight,
            **empty,
            "available": False,
            "reason": "MILP optimum failed the post-solve numerical audit.",
            "solver_audit": solver_audit,
        }

    mapping: dict[tuple[str, str], str] = {}
    assignment_rows: list[dict[str, object]] = []
    for index, block in enumerate(blocks):
        selected = [rater for rater in raters if rounded[variable(index, rater)] > 0.5]
        if len(selected) != 1:
            raise RuntimeError(f"Block {index} did not have exactly one target Rater.")
        target = selected[0]
        person = str(block["Person"])
        source_rater = str(block["SourceRater"])
        mapping[(person, source_rater)] = target
        assignment_rows.append(
            {
                "Block": index,
                "Person": person,
                "SourceRater": source_rater,
                "CounterfactualRater": target,
                "Rows": int(block["Rows"]),
                "ContextCountSignature": block["ContextCountSignature"],
                "Changed": bool(source_rater != target),
            }
        )
    assignment_map = pd.DataFrame(assignment_rows)
    counterfactual = source.copy()
    counterfactual[rater_col] = [
        mapping[(str(person), str(rater))]
        for person, rater in zip(source[person_col], source[rater_col])
    ]

    original_edges = set(
        map(
            tuple,
            source[[person_col, rater_col]].drop_duplicates().itertuples(
                index=False, name=None
            ),
        )
    )
    counter_edges = set(
        map(
            tuple,
            counterfactual[[person_col, rater_col]].drop_duplicates().itertuples(
                index=False, name=None
            ),
        )
    )
    before_objective = float(
        np.mean(
            [person_coordinate[person] * rater_coordinate[rater] for person, rater in original_edges]
        )
    )
    after_objective = float(
        np.mean(
            [person_coordinate[person] * rater_coordinate[rater] for person, rater in counter_edges]
        )
    )

    def correlation(edges: set[tuple[str, str]]) -> float:
        ordered = sorted(edges)
        left = np.asarray([person_coordinate[person] for person, _ in ordered], dtype=float)
        right = np.asarray([rater_coordinate[rater] for _, rater in ordered], dtype=float)
        if len(ordered) < 2 or np.std(left) <= 0 or np.std(right) <= 0:
            return float("nan")
        return float(np.corrcoef(left, right)[0, 1])

    moved_source_blocks = int(assignment_map["Changed"].sum())
    changed_edges = len(original_edges.symmetric_difference(counter_edges)) // 2
    improvement = float(sign * (after_objective - before_objective))
    context_audit = _context_margin_audit(
        source,
        counterfactual,
        rater_col=rater_col,
        context_cols=context_cols,
    )
    person_before = source[[person_col, rater_col]].drop_duplicates().groupby(person_col)[rater_col].nunique()
    person_after = counterfactual[[person_col, rater_col]].drop_duplicates().groupby(person_col)[rater_col].nunique()
    rater_people_before = source[[person_col, rater_col]].drop_duplicates().groupby(rater_col)[person_col].nunique()
    rater_people_after = counterfactual[[person_col, rater_col]].drop_duplicates().groupby(rater_col)[person_col].nunique()
    rater_rows_before = source.groupby(rater_col, observed=True).size()
    rater_rows_after = counterfactual.groupby(rater_col, observed=True).size()
    witness_ok = True
    present_a: list[bool] = []
    present_b: list[bool] = []
    for witness in witnesses.itertuples(index=False):
        person = str(witness.WitnessPerson)
        has_a = (person, str(witness.RaterA)) in counter_edges
        has_b = (person, str(witness.RaterB)) in counter_edges
        present_a.append(has_a)
        present_b.append(has_b)
        witness_ok = witness_ok and has_a and has_b
    if not witnesses.empty:
        witnesses["PresentAtRaterAAfter"] = present_a
        witnesses["PresentAtRaterBAfter"] = present_b
        witnesses["LockSatisfied"] = witnesses[
            ["PresentAtRaterAAfter", "PresentAtRaterBAfter"]
        ].all(axis=1)
    invariant_rows = [
        ("Row count", len(source) == len(counterfactual), f"{len(source)} -> {len(counterfactual)}"),
        ("Each source block assigned once", len(mapping) == block_count, f"blocks={len(mapping)}/{block_count}"),
        ("Unique Person-target Rater block", len(counter_edges) == block_count, f"edges={len(counter_edges)}/{block_count}"),
        ("Per-Person Rater degree", _series_exact(person_before, person_after), "exact equality"),
        ("Per-Rater Person exposure", _series_exact(rater_people_before, rater_people_after), "exact equality"),
        ("Per-Rater response-row exposure", _series_exact(rater_rows_before, rater_rows_after), "exact equality"),
        ("Per-Rater context-cell row margins", bool(context_audit["ExactMatch"].all()), f"cells={int(context_audit['ExactMatch'].sum())}/{len(context_audit)}"),
        ("Connectivity witness locks", bool(witness_ok), f"tree_edges={len(witnesses)}"),
        ("Connected direct Rater overlap", _overlap_components(counter_edges) == 1, f"components={_overlap_components(counter_edges)}"),
        (
            "Assignment changed",
            changed_edges > 0,
            f"changed_edges={changed_edges}; moved_source_blocks={moved_source_blocks}",
        ),
        ("Direction-adjusted objective improved", improvement > 1e-12, f"gain={improvement:.12g}"),
        ("Score-free design", "Score" not in counterfactual.columns, "observed outcome omitted"),
    ]
    invariants = pd.DataFrame(invariant_rows, columns=["Invariant", "Passed", "Evidence"])
    available = bool(invariants["Passed"].all())
    trajectory = pd.DataFrame(
        [
            {
                "Scenario": "observed_assignment",
                "Objective": before_objective,
                "DirectionAdjustedObjective": sign * before_objective,
                "AssignmentRankCorrelation": correlation(original_edges),
                "ChangedBlocks": 0,
                "ChangedEdges": 0,
            },
            {
                "Scenario": "context_margin_counterfactual",
                "Objective": after_objective,
                "DirectionAdjustedObjective": sign * after_objective,
                "AssignmentRankCorrelation": correlation(counter_edges),
                "ChangedBlocks": moved_source_blocks,
                "ChangedEdges": changed_edges,
            },
        ]
    )
    return {
        **preflight,
        "available": available,
        "reason": (
            "A globally optimized score-free endpoint with exact Rater-by-context margins was constructed."
            if available
            else "The optimal MILP endpoint did not provide a changed, improving invariant-preserving contrast."
        ),
        "schema_version": CONTEXT_MILP_SCHEMA_VERSION,
        "direction": direction,
        "design": counterfactual if available else pd.DataFrame(),
        "assignment_map": assignment_map,
        "trajectory": trajectory,
        "invariants": invariants,
        "context_margin_audit": context_audit,
        "solver_audit": solver_audit,
        "witnesses": witnesses,
        "person_coordinate": pd.DataFrame(
            {"Person": list(person_coordinate), "RankCoordinate": list(person_coordinate.values())}
        ),
        "rater_coordinate": pd.DataFrame(
            {"Rater": list(rater_coordinate), "SeverityRankCoordinate": list(rater_coordinate.values())}
        ),
        "claim_boundary": (
            "This endpoint is globally optimal only within the declared binary assignment model, exact margin "
            "constraints, and locked connectivity witnesses. It contains no observed Score, does not identify "
            "the actual assignment mechanism, is not an estimated propensity, and is not a causal-bias estimate."
        ),
    }


def build_context_margin_endpoint_pair(
    data: pd.DataFrame,
    perturbation: Mapping[str, object],
    *,
    person_col: str = "Person",
    rater_col: str = "Rater",
    facet_cols: Iterable[str] = (),
    context_cols: Iterable[str] = (),
    weight_col: str | None = "Weight",
) -> dict[str, object]:
    """Package the observed design and one MILP endpoint for paired refits.

    Unlike the strict 2-switch runner, a meaningful nested MILP dose path has
    not been defined.  This function therefore exposes only 0 and 1 and labels
    the latter as an optimization endpoint, not as a propensity or switch dose.
    """

    if not isinstance(perturbation, Mapping) or not perturbation.get("available"):
        return {
            "available": False,
            "reason": "An available context-margin MILP perturbation is required.",
            "designs": {},
            "dose_table": pd.DataFrame(),
            "path_invariants": pd.DataFrame(),
            "path_assignment_map": pd.DataFrame(),
            "trajectory": pd.DataFrame(),
        }
    facet_cols = [str(column) for column in facet_cols]
    context_cols = [str(column) for column in context_cols]
    work_cols = list(dict.fromkeys([person_col, rater_col, *facet_cols, *context_cols]))
    if weight_col and weight_col in data.columns:
        work_cols.append(weight_col)
    missing = [column for column in work_cols if column not in data.columns]
    if missing:
        raise ValueError("Endpoint source is missing columns: " + ", ".join(missing))
    observed = data[work_cols].copy().reset_index(drop=True)
    for column in dict.fromkeys([person_col, rater_col, *facet_cols, *context_cols]):
        observed[column] = observed[column].astype(str)
    observed["_SourceRow"] = np.arange(len(observed), dtype=int)
    endpoint = perturbation.get("design", pd.DataFrame())
    if not isinstance(endpoint, pd.DataFrame) or endpoint.empty:
        raise ValueError("MILP endpoint design is unavailable.")
    endpoint = endpoint.copy().sort_values("_SourceRow").reset_index(drop=True)
    if not observed["_SourceRow"].equals(endpoint["_SourceRow"]):
        raise ValueError("Observed and endpoint source-row identities differ.")

    trajectory = perturbation.get("trajectory", pd.DataFrame()).copy()
    if not isinstance(trajectory, pd.DataFrame) or len(trajectory) != 2:
        raise ValueError("MILP endpoint trajectory must contain exactly two rows.")
    trajectory = trajectory.reset_index(drop=True)
    trajectory.loc[0, "Scenario"] = "observed_assignment"
    trajectory.loc[1, "Scenario"] = "aligned_counterfactual"
    trajectory["AchievedAlignmentDose"] = [0.0, 1.0]
    dose_table = pd.DataFrame(
        [
            {
                "Scenario": "observed_assignment",
                "RequestedAlignmentDoses": "0.000000",
                "RequestedDoseCount": 1,
                "SwitchesApplied": 0,
                "TotalSwitches": 0,
                "AchievedAlignmentDose": 0.0,
                "AssignmentRankCorrelation": float(trajectory.iloc[0]["AssignmentRankCorrelation"]),
                "Objective": float(trajectory.iloc[0]["Objective"]),
                "ChangedBlocks": 0,
                "IsObserved": True,
                "IsEndpoint": False,
                "OptimizationEndpoint": False,
            },
            {
                "Scenario": "aligned_counterfactual",
                "RequestedAlignmentDoses": "1.000000",
                "RequestedDoseCount": 1,
                "SwitchesApplied": 0,
                "TotalSwitches": 0,
                "AchievedAlignmentDose": 1.0,
                "AssignmentRankCorrelation": float(trajectory.iloc[1]["AssignmentRankCorrelation"]),
                "Objective": float(trajectory.iloc[1]["Objective"]),
                "ChangedBlocks": int(trajectory.iloc[1]["ChangedBlocks"]),
                "IsObserved": False,
                "IsEndpoint": True,
                "OptimizationEndpoint": True,
            },
        ]
    )

    endpoint_invariants = perturbation.get("invariants", pd.DataFrame()).copy()
    invariant_parts: list[pd.DataFrame] = []
    for scenario, dose in (("observed_assignment", 0.0), ("aligned_counterfactual", 1.0)):
        current = endpoint_invariants.copy()
        current.insert(0, "AchievedAlignmentDose", dose)
        current.insert(0, "Scenario", scenario)
        if scenario == "observed_assignment":
            identity_exceptions = current["Invariant"].isin(
                ["Assignment changed", "Direction-adjusted objective improved"]
            )
            current.loc[~identity_exceptions, "Passed"] = True
            current.loc[~identity_exceptions, "Evidence"] = "observed baseline identity"
            current.loc[identity_exceptions, "Passed"] = False
            current.loc[identity_exceptions, "Evidence"] = "not required at observed baseline"
            current["RequiredAtDose"] = ~identity_exceptions
        else:
            current["RequiredAtDose"] = True
        current["PassedForDose"] = current["Passed"] | ~current["RequiredAtDose"]
        invariant_parts.append(current)
    path_invariants = pd.concat(invariant_parts, ignore_index=True)

    endpoint_map = perturbation.get("assignment_map", pd.DataFrame()).copy()
    observed_map = endpoint_map.copy()
    observed_map["CounterfactualRater"] = observed_map["SourceRater"]
    observed_map["Changed"] = False
    observed_map.insert(0, "AchievedAlignmentDose", 0.0)
    observed_map.insert(0, "Scenario", "observed_assignment")
    endpoint_map.insert(0, "AchievedAlignmentDose", 1.0)
    endpoint_map.insert(0, "Scenario", "aligned_counterfactual")
    path_map = pd.concat([observed_map, endpoint_map], ignore_index=True)

    available = bool(
        path_invariants.loc[path_invariants["RequiredAtDose"], "PassedForDose"].all()
    )
    return {
        "available": available,
        "reason": (
            "Score-free observed and context-margin MILP endpoint designs were packaged."
            if available
            else "At least one required endpoint invariant failed."
        ),
        "schema_version": CONTEXT_MILP_SCHEMA_VERSION,
        "designs": {
            "observed_assignment": observed,
            "aligned_counterfactual": endpoint,
        },
        "dose_table": dose_table,
        "path_invariants": path_invariants,
        "path_assignment_map": path_map,
        "trajectory": trajectory,
        "claim_boundary": (
            "Only observed versus one constrained MILP endpoint is defined. Intermediate dose values are "
            "not available because nested feasible context-margin assignments have not been qualified."
        ),
    }
