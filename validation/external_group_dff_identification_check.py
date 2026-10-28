#!/usr/bin/env python3
"""Bounded design checks for the proposed two-school RSM interaction model.

No scores, estimation, MML integration or operating-characteristic study.
Run with: python -m validation.external_group_dff_identification_check OUTPUT
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

import streamlit_app as app


CASES = ("crossed", "nested_references", "missing_group_cell")
MODELS = ("reference_relative", "symmetric")
RATERS = ("A_R1", "A_R2", "B_R1", "B_R2")


def teaching_design(case, *, full_lesson=False):
    """Keep the lesson's IDs; use two persons per school and two references."""
    data = app.sample_mfrm_data_by_key("two_school_linking")
    groups = app._two_school_sample_groups(data).set_index("Person")["School"]
    if not full_lesson:
        data = data.loc[data.Person.isin(["A01", "A02", "B01", "B02", "C01", "C02"])]
    data = data[["Person", "Rater", "Criterion"]].drop_duplicates()
    groups = groups.reindex(sorted(data.Person.unique())).replace("", pd.NA)
    if case in {"crossed", "missing_group_cell"}:
        data = pd.MultiIndex.from_product(
            [groups.index, RATERS, sorted(data.Criterion.unique())], names=data.columns,
        ).to_frame(index=False)
        if case == "missing_group_cell":
            data = data.loc[~(data.Person.map(groups).eq("B") & data.Rater.eq("A_R1"))]
    elif case == "one_cross_rating":
        extra = data.loc[data.Person.eq("A01") & data.Rater.eq("A_R1")].copy()
        extra["Rater"] = "B_R1"
        data = pd.concat([data, extra])
    elif case != "nested_references":
        raise ValueError(f"Unknown bounded design: {case}")
    return data.reset_index(drop=True), groups


def model_matrices(data, groups, interaction="reference_relative"):
    """Full adjacent-logit design and the app's sum-zero free coordinates.

    Four RSM transitions are included. Person locations are free; rater,
    criterion and step locations sum to zero. Each group's gamma sums to zero
    over raters; references have no interaction term. The optional symmetric
    model additionally imposes gamma[r,A] = -gamma[r,B], a substantive
    restriction once reference observations define the rater baseline.
    """
    levels = {name: sorted(data[name].unique()) for name in data.columns}
    levels.update(Step=list(range(4)), Gamma=[(r, g) for r in RATERS for g in ("A", "B")])
    columns, blocks, expansions = {}, {}, []
    offset = 0
    for name, values in levels.items():
        blocks[name] = slice(offset, offset + len(values))
        columns.update({(name, value): offset + i for i, value in enumerate(values)})
        offset += len(values)
        if name == "Gamma":
            if interaction not in MODELS:
                raise ValueError(f"Unknown interaction specification: {interaction}")
            group_expansion = np.eye(2) if interaction == "reference_relative" else np.array([[1], [-1]])
            expansion = np.kron(expansions[1], group_expansion)
        else:
            spec = app.build_facet_constraint(values, centered=name != "Person")
            expansion = app._constraint_expansion_matrix_for_identifiability(spec)
        expansions.append(expansion)
    expansion = np.zeros((offset, sum(e.shape[1] for e in expansions)))
    start = 0
    for block, part in zip(blocks.values(), expansions):
        expansion[block, start:start + part.shape[1]] = part
        start += part.shape[1]

    raw = np.zeros((len(data) * 4, offset), dtype=int)
    for obs, row in enumerate(data.itertuples(index=False)):
        selected = slice(4 * obs, 4 * obs + 4)
        for name, sign in (("Person", 1), ("Rater", -1), ("Criterion", -1)):
            raw[selected, columns[name, getattr(row, name)]] = sign
        raw[np.arange(4 * obs, 4 * obs + 4), np.arange(blocks["Step"].start, blocks["Step"].stop)] = -1
        group = groups.loc[row.Person]
        if pd.notna(group):
            raw[selected, columns["Gamma", (row.Rater, group)]] = 1

    constraints = []
    for name in ("Rater", "Criterion", "Step"):
        constraint = np.zeros(offset, dtype=int)
        constraint[blocks[name]] = 1
        constraints.append(constraint)
    interaction_constraints = [[(r, g) for r in RATERS] for g in ("A", "B")]
    if interaction == "symmetric":
        interaction_constraints += [[(r, g) for g in ("A", "B")] for r in RATERS]
    for pairs in interaction_constraints:
        constraint = np.zeros(offset, dtype=int)
        for pair in pairs:
            constraint[columns["Gamma", pair]] = 1
        constraints.append(constraint)
    return raw, np.asarray(constraints), expansion, columns


def exact_row_basis(matrix):
    """Rational row elimination: independent of floating point SVD/coding."""
    basis = {}
    for values in matrix:
        row = [Fraction(int(value)) for value in values]
        for pivot, known in sorted(basis.items()):
            if row[pivot]:
                factor = row[pivot]
                row = [a - factor * b for a, b in zip(row, known)]
        pivot = next((i for i, value in enumerate(row) if value), None)
        if pivot is not None:
            scale = row[pivot]
            basis[pivot] = [value / scale for value in row]
    return basis


def exact_estimable(contrast, basis):
    row = [Fraction(int(value)) for value in contrast]
    for pivot, known in sorted(basis.items()):
        factor = row[pivot]
        row = [a - factor * b for a, b in zip(row, known)]
    return not any(row)


def planned_contrasts(columns):
    differences = {}
    for rater in RATERS:
        row = np.zeros(len(columns), dtype=int)
        row[columns["Gamma", (rater, "A")]] = 1
        row[columns["Gamma", (rater, "B")]] = -1
        differences[rater] = row
    contrasts = [(f"{r}: A minus B", row, [r], ("A", "B")) for r, row in differences.items()]
    for left, right in (("A_R1", "A_R2"), ("B_R1", "B_R2"), ("A_R1", "B_R1"), ("A_R2", "B_R1")):
        contrasts.append((f"({left}: A-B) minus ({right}: A-B)",
                          differences[left] - differences[right], [left, right], ("A", "B")))
    for group, left, right in (("A", "A_R1", "A_R2"), ("B", "B_R1", "B_R2")):
        row = np.zeros(len(columns), dtype=int)
        row[columns["Gamma", (left, group)]] = 1
        row[columns["Gamma", (right, group)]] = -1
        contrasts.append((f"Within {group}: {left} minus {right}", row, [left, right], (group,)))
    return contrasts


def audit_case(case, interaction="reference_relative"):
    data, groups = teaching_design(case)
    raw, constraints, expansion, columns = model_matrices(data, groups, interaction)
    design = raw @ expansion
    _, singular_values, vh = np.linalg.svd(design, full_matrices=False)
    tolerance = max(design.shape) * np.finfo(float).eps * singular_values[0]
    rank = int(np.sum(singular_values > tolerance))
    null_basis = vh[rank:]
    # The full-coordinate check does not use the app's constraint expansion.
    row_basis = exact_row_basis(np.vstack([raw, constraints]))
    constraint_rank = len(exact_row_basis(constraints))
    exact_rank = len(row_basis) - constraint_rank
    assert expansion.shape[1] == raw.shape[1] - constraint_rank
    assert np.linalg.matrix_rank(expansion) == expansion.shape[1]
    assert np.array_equal(constraints @ expansion, np.zeros((len(constraints), expansion.shape[1])))
    assert rank == exact_rank
    observed = {(r.Rater, groups.loc[r.Person]) for r in data.itertuples()
                if pd.notna(groups.loc[r.Person])}
    contrasts = []
    for name, contrast, raters, compared_groups in planned_contrasts(columns):
        null_projection = float(np.linalg.norm(null_basis @ (contrast @ expansion)))
        estimable = null_projection < 1e-10
        exact = exact_estimable(contrast, row_basis)
        assert estimable == exact
        contrasts.append(dict(
            Case=case, Model=interaction, Contrast=name, Estimable=estimable, ExactEstimable=exact,
            ComparisonScope="between_groups" if len(compared_groups) == 2 else "within_group",
            NullProjection=null_projection,
            AllComparedCellsObserved=all((r, g) in observed for r in raters for g in compared_groups),
        ))
    summary = dict(
        Case=case, Model=interaction, Persons=len(groups), Ratings=len(data), AdjacentLogits=len(design),
        Estimator="JML", RatingModel="RSM",
        ObservedStudentCells=len(observed), PlannedStudentCells=8,
        FreeParameters=design.shape[1], Rank=rank, ExactRank=exact_rank,
        Nullity=design.shape[1] - rank, SVDTolerance=float(tolerance),
    )
    return summary, contrasts, data.join(groups.rename("School"), on="Person")


def build_output_plan(summary, contrasts, *, targets=None):
    """Exercise a design-only output specification on these bounded fixtures.

    This is not a public-data audit or an estimator qualification function.
    Unknown or inconsistent evidence raises instead of certifying a design.
    """
    if (summary.get("Case") not in CASES or summary.get("Model") not in MODELS
            or summary.get("Estimator") != "JML" or summary.get("RatingModel") != "RSM"):
        raise ValueError("The output specification covers only the declared design fixtures.")
    rank, size, nullity = (summary.get(k) for k in ("Rank", "FreeParameters", "Nullity"))
    if (any(type(v) is not int for v in (rank, size, nullity))
            or not 0 <= rank <= size or size == 0 or nullity != size - rank
            or summary.get("ExactRank") != rank):
        raise ValueError("A consistent constrained-rank audit is required.")
    by_name = {}
    for row in contrasts:
        if (row.get("Case") != summary["Case"] or row.get("Model") != summary["Model"]
                or row.get("ComparisonScope") not in {"between_groups", "within_group"}
                or any(type(row.get(k)) is not bool for k in
                       ("Estimable", "ExactEstimable", "AllComparedCellsObserved"))
                or row["Estimable"] != row["ExactEstimable"]
                or (nullity == 0 and not row["Estimable"])):
            raise ValueError("Contrast evidence is incomplete or inconsistent with the design audit.")
        name = row.get("Contrast")
        if not isinstance(name, str) or not name.strip() or name in by_name:
            raise ValueError("Contrasts must have unique, explicit names.")
        by_name[name] = row
    if targets is None:
        targets = [name for name, row in by_name.items() if row["ComparisonScope"] == "between_groups"]
    if (not isinstance(targets, (list, tuple)) or not targets
            or any(not isinstance(name, str) for name in targets)
            or len(set(targets)) != len(targets)
            or any(name not in by_name or by_name[name]["ComparisonScope"] != "between_groups"
                   for name in targets)):
        raise ValueError("Select explicit between-group contrasts; within-group contrasts are a separate target.")
    selected = [by_name[name] for name in targets]
    n_estimable = sum(row["Estimable"] for row in selected)
    state = ("target_unidentified" if n_estimable == 0 else
             "model_identified" if nullity == 0 else "partially_identified")
    return dict(
        Case=summary["Case"], Model=summary["Model"], State=state,
        Estimator="JML", RatingModel="RSM", Stage="design_only",
        TargetCount=len(selected), EstimableTargetCount=n_estimable, ModelNullity=nullity,
        PointEstimatesReleased=False, IntervalsReleased=False, PValuesReleased=False,
        GroupCurvesReleased=False,
        Comparisons=[dict(
            Case=row["Case"], Model=row["Model"], Contrast=row["Contrast"],
            StructurallyEstimable=row["Estimable"], AllComparedCellsObserved=row["AllComparedCellsObserved"],
            OutputScope="Design review only; no effect estimate",
        ) for row in selected],
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    summaries, contrasts, output_plans = [], [], []
    for case in CASES:
        for interaction in MODELS:
            summary, rows, data = audit_case(case, interaction)
            summaries.append(summary)
            contrasts.extend(rows)
            if interaction == "reference_relative":
                output_plans.append(build_output_plan(summary, rows))
        data.to_csv(args.output / f"{case}_assignments.csv", index=False)
    pd.DataFrame(summaries).to_csv(args.output / "design_ranks.csv", index=False)
    pd.DataFrame(contrasts).to_csv(args.output / "contrast_estimability.csv", index=False)
    (args.output / "output_contract_cases.json").write_text(json.dumps(output_plans, indent=2) + "\n")
    pd.DataFrame([row for plan in output_plans for row in plan["Comparisons"]]).to_csv(
        args.output / "planned_group_comparisons.csv", index=False,
    )
    repo = Path(__file__).resolve().parents[1]
    sources = [Path(__file__).resolve(), repo / "streamlit_app.py",
               repo / "tests/test_external_group_dff_identification.py"]
    report = dict(
        scope="Three deterministic JML/RSM design checks; no estimation or MML qualification",
        assumptions=["Free person locations; rater, criterion and RSM steps sum to zero",
                     "Group-specific sum-zero rater departures, compared with an additional A/B symmetry restriction",
                     "Reference recordings outside A/B; interaction absent on references"],
        cases=summaries, contrasts=contrasts, formal_studies="remain paused",
        sources={str(p.relative_to(repo)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
    )
    (args.output / "identification_review.json").write_text(json.dumps(report, indent=2) + "\n")
    print(pd.DataFrame(summaries).to_string(index=False))


if __name__ == "__main__":
    main()
