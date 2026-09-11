#!/usr/bin/env python3
"""Run the prospectively registered finite-CMLE existence stress."""

from __future__ import annotations

import argparse
import hashlib
from itertools import product
import json
from pathlib import Path
import shutil
import sys
import time

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from mfrm_app.cmle import fit_cmle, prepare_cmle_design  # noqa: E402
from mfrm_app.cmle_existence import audit_cmle_finite_mle  # noqa: E402
from validation.cmle_native_anchor_connectivity_stress import (  # noqa: E402
    anchor_levels,
    hard_anchors,
)


DEFAULT_PLAN = ROOT / "validation/cmle_finite_mle_existence_plan_20260810.json"
DEFAULT_OUTPUT = ROOT / "validation/cmle_finite_mle_existence_20260810"
CONNECTIVITY_PLAN = ROOT / "validation/cmle_native_anchor_connectivity_plan_20260810.json"
CONNECTIVITY_RESPONSES = ROOT / "validation/cmle_native_anchor_connectivity_20260810/connectivity_responses.csv"
REMEDIATION_DIR = ROOT / "validation/cmle_finite_domain_optimizer_remediation_20260810"
PRIMARY_REPLAY_MAX_CONFIGURATIONS = 50_000
PRIMARY_REPLAY_MAX_CONSTRAINTS = 50_000
SENTINEL_MAX_CONFIGURATIONS = 300_000
SENTINEL_MAX_CONSTRAINTS = 150_000


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, index=False, lineterminator="\n", float_format="%.17g")


def validate_plan(path: Path) -> tuple[dict[str, object], dict[str, object]]:
    plan = json.loads(path.read_text(encoding="utf-8"))
    if plan.get("study_id") != "cmle_finite_mle_existence_v1":
        raise ValueError("Unexpected finite-MLE existence plan.")
    mismatches = [
        relative
        for relative, digest in plan["parent_identity"].items()
        if not (ROOT / relative).exists() or sha256_file(ROOT / relative) != digest
    ]
    if mismatches:
        raise ValueError(f"Finite-MLE parent identity failed: {mismatches}")
    remediation = json.loads(
        (REMEDIATION_DIR / "remediation_decision.json").read_text(encoding="utf-8")
    )
    output_mismatches = [
        name
        for name, digest in remediation["output_sha256"].items()
        if not (REMEDIATION_DIR / name).exists()
        or sha256_file(REMEDIATION_DIR / name) != digest
    ]
    if output_mismatches:
        raise ValueError(f"Remediation evidence identity failed: {output_mismatches}")
    connectivity_plan = json.loads(CONNECTIVITY_PLAN.read_text(encoding="utf-8"))
    return plan, connectivity_plan


def _anchor_table(levels: list[str], values: list[float] | None = None) -> pd.DataFrame | None:
    if not levels:
        return None
    supplied = values if values is not None else [0.0] * len(levels)
    return pd.DataFrame(
        [
            {
                "ParameterType": "Facet",
                "Facet": "Rater",
                "Level": level,
                "Value": float(value),
            }
            for level, value in zip(levels, supplied, strict=True)
        ]
    )


def _binary_two_rater_frame(patterns: list[tuple[int, int]]) -> pd.DataFrame:
    rows = []
    for index, (first, second) in enumerate(patterns):
        rows.extend(
            [
                (f"P{index:04d}", "R1", first),
                (f"P{index:04d}", "R2", second),
            ]
        )
    return pd.DataFrame(rows, columns=["Person", "Rater", "Score"])


def _three_rater_frame(mechanism: str) -> pd.DataFrame:
    rows: list[tuple[str, str, str, int]] = []
    person_index = 0

    def add_edge(left: str, right: str, patterns: list[tuple[int, int]]) -> None:
        nonlocal person_index
        for first, second in patterns:
            person = f"P{person_index:04d}"
            rows.extend(
                [(person, left, "E1", first), (person, right, "E2", second)]
            )
            person_index += 1

    mixed = [(1, 0)] * 5 + [(0, 1)] * 5
    if mechanism == "finite":
        add_edge("R1", "R2", mixed)
        add_edge("R2", "R3", mixed)
    elif mechanism == "one_separated_edge":
        add_edge("R1", "R2", [(1, 0)] * 10)
        add_edge("R2", "R3", mixed)
    elif mechanism == "opposing_edges":
        add_edge("R1", "R2", [(1, 0)] * 10)
        add_edge("R2", "R3", [(0, 1)] * 10)
    elif mechanism == "monotone_edges":
        add_edge("R1", "R2", [(1, 0)] * 10)
        add_edge("R2", "R3", [(1, 0)] * 10)
    elif mechanism == "disconnected":
        add_edge("R1", "R2", mixed)
        for _ in range(10):
            person = f"P{person_index:04d}"
            rows.extend(
                [(person, "R3", "E1", 1), (person, "R3", "E2", 0)]
            )
            person_index += 1
    else:  # pragma: no cover - frozen manifest controls this
        raise ValueError(f"Unknown three-Rater mechanism: {mechanism}")
    return pd.DataFrame(rows, columns=["Person", "Rater", "Event", "Score"])


def _rsm_frame(mechanism: str) -> pd.DataFrame:
    if mechanism == "finite":
        patterns = [
            pair for pair in product(range(4), repeat=2) if 0 < sum(pair) < 6
        ]
    elif mechanism == "separated":
        choices = [(1, 0), (2, 0), (3, 0), (3, 1), (3, 2)]
        patterns = [choices[index % len(choices)] for index in range(50)]
    elif mechanism == "unused_declared_top_category":
        patterns = [
            pair for pair in product(range(3), repeat=2) if 0 < sum(pair) < 6
        ]
    else:  # pragma: no cover
        raise ValueError(f"Unknown RSM mechanism: {mechanism}")
    return _binary_two_rater_frame(patterns)


def _pcm_frame(mechanism: str) -> pd.DataFrame:
    units = [("R1", "C1"), ("R1", "C2"), ("R2", "C1"), ("R2", "C2")]
    if mechanism == "finite":
        configurations = [
            values for values in product(range(3), repeat=4) if sum(values) == 4
        ]
    elif mechanism == "separated":
        configurations = [(2, 2, 0, 0)] * 20
    elif mechanism == "unused_declared_top_category":
        configurations = [
            values
            for values in product(range(2), repeat=4)
            if 0 < sum(values) < 8
        ]
    else:  # pragma: no cover
        raise ValueError(f"Unknown PCM mechanism: {mechanism}")
    rows = []
    for person_index, values in enumerate(configurations):
        for (rater, criterion), score in zip(units, values, strict=True):
            rows.append((f"P{person_index:04d}", rater, criterion, score))
    return pd.DataFrame(
        rows, columns=["Person", "Rater", "Criterion", "Score"]
    )


def fixture_cases() -> list[dict[str, object]]:
    cases: list[dict[str, object]] = []

    def add(
        case_id: str,
        family: str,
        mechanism: str,
        frame: pd.DataFrame,
        prepare_args: dict[str, object],
        expected: str,
        anchor_label: str,
    ) -> None:
        cases.append(
            {
                "CaseId": case_id,
                "Family": family,
                "Mechanism": mechanism,
                "AnchorLabel": anchor_label,
                "ExpectedStatus": expected,
                "frame": frame,
                "prepare_args": prepare_args,
            }
        )

    for persons in (2, 5, 10, 20, 100, 300):
        balanced = [(1, 0)] * max(1, persons // 2) + [(0, 1)] * (
            persons - max(1, persons // 2)
        )
        mechanisms = {
            "separated_forward": ([(1, 0)] * persons, "boundary_no_finite_cmle"),
            "separated_reverse": ([(0, 1)] * persons, "boundary_no_finite_cmle"),
            "balanced_finite": (balanced, "interior_finite_cmle_supported"),
            "near_forward": ([(1, 0)] * (persons - 1) + [(0, 1)], "interior_finite_cmle_supported"),
            "near_reverse": ([(0, 1)] * (persons - 1) + [(1, 0)], "interior_finite_cmle_supported"),
        }
        for mechanism, (patterns, expected) in mechanisms.items():
            for anchor_label, anchors in (
                ("unanchored", None),
                ("R1_correct", _anchor_table(["R1"])),
            ):
                add(
                    f"binary2_n{persons:03d}_{mechanism}_{anchor_label}",
                    "binary_two_rater",
                    mechanism,
                    _binary_two_rater_frame(patterns),
                    {
                        "person_col": "Person",
                        "facet_cols": ["Rater"],
                        "score_col": "Score",
                        "rating_min": 0,
                        "rating_max": 1,
                        "model": "RSM",
                        "hard_anchors": anchors,
                    },
                    expected,
                    anchor_label,
                )

    three_expected = {
        "finite": {
            "unanchored": "interior_finite_cmle_supported",
            "R1_correct": "interior_finite_cmle_supported",
            "R1_R3_correct": "interior_finite_cmle_supported",
        },
        "one_separated_edge": {
            "unanchored": "boundary_no_finite_cmle",
            "R1_correct": "boundary_no_finite_cmle",
            "R1_R3_correct": "interior_finite_cmle_supported",
        },
        "opposing_edges": {
            "unanchored": "boundary_no_finite_cmle",
            "R1_correct": "boundary_no_finite_cmle",
            "R1_R3_correct": "boundary_no_finite_cmle",
        },
        "monotone_edges": {
            "unanchored": "boundary_no_finite_cmle",
            "R1_correct": "boundary_no_finite_cmle",
            "R1_R3_correct": "interior_finite_cmle_supported",
        },
        "disconnected": {
            "unanchored": "structural_nonidentification",
            "R1_correct": "structural_nonidentification",
            "R1_R3_correct": "interior_finite_cmle_supported",
        },
    }
    for mechanism, expected_by_anchor in three_expected.items():
        for anchor_label, anchors in (
            ("unanchored", None),
            ("R1_correct", _anchor_table(["R1"])),
            ("R1_R3_correct", _anchor_table(["R1", "R3"])),
        ):
            add(
                f"binary3_{mechanism}_{anchor_label}",
                "binary_three_rater_graph",
                mechanism,
                _three_rater_frame(mechanism),
                {
                    "person_col": "Person",
                    "facet_cols": ["Rater"],
                    "score_col": "Score",
                    "rating_min": 0,
                    "rating_max": 1,
                    "model": "RSM",
                    "response_unit_col": "Event",
                    "hard_anchors": anchors,
                },
                expected_by_anchor[anchor_label],
                anchor_label,
            )

    for family, maker, extra in (
        ("polytomous_rsm", _rsm_frame, {}),
        ("polytomous_pcm", _pcm_frame, {"step_facet": "Criterion"}),
    ):
        for mechanism, expected in (
            ("finite", "interior_finite_cmle_supported"),
            ("separated", "boundary_no_finite_cmle"),
            ("unused_declared_top_category", "boundary_no_finite_cmle"),
        ):
            for anchor_label, anchors in (
                ("unanchored", None),
                ("R1_correct", _anchor_table(["R1"], [0.0])),
                ("R1_contaminated", _anchor_table(["R1"], [0.25])),
            ):
                facet_cols = ["Rater"] if family == "polytomous_rsm" else ["Rater", "Criterion"]
                add(
                    f"{family}_{mechanism}_{anchor_label}",
                    family,
                    mechanism,
                    maker(mechanism),
                    {
                        "person_col": "Person",
                        "facet_cols": facet_cols,
                        "score_col": "Score",
                        "rating_min": 0,
                        "rating_max": 3 if family == "polytomous_rsm" else 2,
                        "model": "RSM" if family == "polytomous_rsm" else "PCM",
                        "hard_anchors": anchors,
                        **extra,
                    },
                    expected,
                    anchor_label,
                )
    return cases


def _append_table(
    frames: list[pd.DataFrame], table: pd.DataFrame, identity: dict[str, object]
) -> None:
    if table.empty:
        return
    value = table.copy()
    for column, item in reversed(list(identity.items())):
        value.insert(0, column, item)
    frames.append(value)


def run_fixture_case(case: dict[str, object]) -> tuple[
    dict[str, object], pd.DataFrame, pd.DataFrame, pd.DataFrame
]:
    started = time.perf_counter()
    design = prepare_cmle_design(case["frame"], **case["prepare_args"])
    audit = audit_cmle_finite_mle(design)
    summary = audit["summary"].iloc[0].to_dict()
    fit_returned = False
    fit_ready = False
    fit_failure = ""
    max_estimate = np.nan
    max_se = np.nan
    if bool(design.audit["eligible"]):
        try:
            fit = fit_cmle(case["frame"], **case["prepare_args"], gtol=1e-8, maxiter=800)
            fit_returned = True
            fit_summary = fit["summary"].iloc[0]
            fit_ready = bool(fit_summary["InferenceReady"])
            structural = pd.concat(
                [fit["facets"]["others"], fit["steps"]], ignore_index=True
            )
            max_estimate = float(structural["Estimate"].abs().max())
            max_se = float(structural["SE"].abs().max())
        except Exception as exc:  # pragma: no cover - retained as evidence
            fit_failure = f"{type(exc).__name__}: {str(exc)[:500]}"
    row = {
        "CaseId": case["CaseId"],
        "Family": case["Family"],
        "Mechanism": case["Mechanism"],
        "AnchorLabel": case["AnchorLabel"],
        "ExpectedStatus": case["ExpectedStatus"],
        "ObservedStatus": summary["Status"],
        "ExpectedStatusMatch": summary["Status"] == case["ExpectedStatus"],
        "PrefitEligible": bool(design.audit["eligible"]),
        "PrefitRank": int(design.audit["conditional_rank"]),
        "PrefitNullity": int(design.audit["conditional_nullity"]),
        "KParams": int(design.n_parameters),
        "EnumerationComplete": bool(summary.get("EnumerationComplete", False)),
        "LPComplete": bool(summary.get("LPComplete", False)),
        "BoundaryDetected": summary.get("BoundaryDetected", pd.NA),
        "ExistenceQualified": bool(summary.get("ExistenceQualified", False)),
        "Configurations": summary.get("Configurations", np.nan),
        "UniqueNormalizedConstraints": summary.get("UniqueNormalizedConstraints", np.nan),
        "DirectionalObjectiveNonincreasing": summary.get("DirectionalObjectiveNonincreasing", pd.NA),
        "CurrentFitReturned": fit_returned,
        "CurrentInferenceReady": fit_ready,
        "ExistenceQualifiedInferenceReady": bool(fit_ready and summary.get("ExistenceQualified", False)),
        "CurrentFitFailure": fit_failure,
        "MaxAbsoluteEstimate": max_estimate,
        "MaxSE": max_se,
        "ElapsedSeconds": float(time.perf_counter() - started),
    }
    identity = {"CaseId": case["CaseId"], "Family": case["Family"]}
    lp = audit["lp_tolerances"].copy()
    trace = audit["directional_trace"].copy()
    strata = audit["strata"].copy()
    for table in (lp, trace, strata):
        for column, item in reversed(list(identity.items())):
            if not table.empty:
                table.insert(0, column, item)
    return row, lp, trace, strata


def run_retained_case(
    run_frame: pd.DataFrame,
    scenario: dict[str, object],
    connectivity_plan: dict[str, object],
    current_fit: pd.DataFrame,
    *,
    max_configurations: int,
    max_constraints: int,
) -> tuple[dict[str, object], pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    replicate = int(run_frame["Replicate"].iloc[0])
    run_id = str(run_frame["RunId"].iloc[0])
    scenario_name = str(scenario["scenario"])
    levels = anchor_levels(connectivity_plan, replicate, scenario_name)
    anchors = hard_anchors(connectivity_plan, levels)
    analysis = run_frame[["Person", "Rater", "Criterion", "ObservedCategory"]].copy()
    started = time.perf_counter()
    design = prepare_cmle_design(
        analysis,
        person_col="Person",
        facet_cols=["Rater", "Criterion"],
        score_col="ObservedCategory",
        rating_min=0,
        rating_max=3,
        model="RSM",
        hard_anchors=anchors,
    )
    audit = audit_cmle_finite_mle(
        design,
        max_configurations=max_configurations,
        max_constraints=max_constraints,
    )
    summary = audit["summary"].iloc[0].to_dict()
    fitted = current_fit.loc[
        current_fit["RunId"].eq(run_id)
        & current_fit["Scenario"].eq(scenario_name)
    ]
    current_returned = bool(fitted.iloc[0]["Returned"]) if len(fitted) else False
    current_ready = bool(fitted.iloc[0]["InferenceReady"]) if len(fitted) else False
    identity = {
        "Topology": str(run_frame["Topology"].iloc[0]),
        "Replicate": replicate,
        "RunId": run_id,
        "Scenario": scenario_name,
    }
    row = {
        **identity,
        "AnchorCount": len(levels),
        "AnchorLevels": ";".join(levels),
        "PrefitEligible": bool(design.audit["eligible"]),
        "PrefitRank": int(design.audit["conditional_rank"]),
        "PrefitNullity": int(design.audit["conditional_nullity"]),
        "KParams": int(design.n_parameters),
        "ExistenceStatus": summary["Status"],
        "ExistenceReason": summary["Reason"],
        "EnumerationComplete": bool(summary.get("EnumerationComplete", False)),
        "LPComplete": bool(summary.get("LPComplete", False)),
        "BoundaryDetected": summary.get("BoundaryDetected", pd.NA),
        "ExistenceQualified": bool(summary.get("ExistenceQualified", False)),
        "Configurations": summary.get("Configurations", np.nan),
        "UniqueNormalizedConstraints": summary.get("UniqueNormalizedConstraints", np.nan),
        "DirectionalObjectiveNonincreasing": summary.get("DirectionalObjectiveNonincreasing", pd.NA),
        "CurrentFitReturned": current_returned,
        "CurrentInferenceReady": current_ready,
        "ExistenceQualifiedInferenceReady": bool(
            current_ready and summary.get("ExistenceQualified", False)
        ),
        "ElapsedSeconds": float(time.perf_counter() - started),
    }
    for table in (audit["lp_tolerances"], audit["directional_trace"], audit["strata"]):
        for column, item in reversed(list(identity.items())):
            if not table.empty:
                table.insert(0, column, item)
    return row, audit["lp_tolerances"], audit["directional_trace"], audit["strata"]


def report_text(decision: dict[str, object], retained_summary: pd.DataFrame) -> str:
    summary_lines = [
        "| Existence status | Cases | Current ready | Existence-qualified ready |",
        "| --- | ---: | ---: | ---: |",
    ]
    for row in retained_summary.itertuples(index=False):
        summary_lines.append(
            f"| {row.ExistenceStatus} | {int(row.Cases)} | {int(row.CurrentReady)} | {int(row.ExistenceQualifiedReady)} |"
        )
    table = "\n".join(summary_lines)
    return f"""# Native exact-CMLE finite-MLE existence stress

> Prospectively registered research audit. It qualifies conditional convex-support geometry under explicit enumeration caps; it does not establish model fit, estimator performance, anchor validity, or public UI readiness.

## Contract

- Contract passed: `{decision['contract_passed']}`
- Registered fixtures: `{decision['fixture_cases']}`; expected-status matches: `{decision['fixture_expected_status_matches']}`
- Retained cases ledgered: `{decision['retained_cases']}/350`
- Retained prefit eligible: `{decision['retained_prefit_eligible']}/287`
- Boundary / interior / structural / unavailable: `{decision['retained_boundary']}` / `{decision['retained_interior']}` / `{decision['retained_structural_nonidentification']}` / `{decision['retained_unavailable']}`
- Current inference-ready / existence-qualified ready: `{decision['retained_current_ready']}` / `{decision['retained_existence_qualified_ready']}`
- Directional verification failures: `{decision['directional_verification_failures']}`

## Retained same-byte transition

{table}

## Boundary

Structural rank deficiency is reported before separation. A boundary direction with monotone exact likelihood trace withholds finite-CMLE existence; an enumeration/LP cap yields unavailable, not finite. The initial exact enumerator is deliberately capped and is not yet a scalable public algorithm. Current optimizer outputs and Person fit values are unchanged; this evidence projects an additional fail-closed gate only.
"""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--plan", type=Path, default=DEFAULT_PLAN)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    plan_path = args.plan.resolve()
    output = args.output.resolve()
    plan, connectivity_plan = validate_plan(plan_path)
    if output.exists() and any(output.iterdir()):
        if not args.overwrite:
            raise FileExistsError(f"Output is non-empty: {output}; pass --overwrite explicitly.")
        shutil.rmtree(output)
    output.mkdir(parents=True, exist_ok=True)
    shutil.copy2(plan_path, output / plan_path.name)

    fixture_rows = []
    fixture_lp: list[pd.DataFrame] = []
    fixture_trace: list[pd.DataFrame] = []
    fixture_strata: list[pd.DataFrame] = []
    response_frames = []
    cases = fixture_cases()
    for case in cases:
        response = case["frame"].copy()
        response.insert(0, "CaseId", case["CaseId"])
        response_frames.append(response)
        row, lp, trace, strata = run_fixture_case(case)
        fixture_rows.append(row)
        _append_table(fixture_lp, lp, {})
        _append_table(fixture_trace, trace, {})
        _append_table(fixture_strata, strata, {})
    fixtures = pd.DataFrame(fixture_rows)
    fixture_responses = pd.concat(response_frames, ignore_index=True, sort=False)

    retained_responses = pd.read_csv(CONNECTIVITY_RESPONSES, float_precision="round_trip")
    current_fit = pd.read_csv(REMEDIATION_DIR / "candidate_fit_ledger.csv")
    retained_rows = []
    retained_lp: list[pd.DataFrame] = []
    retained_trace: list[pd.DataFrame] = []
    retained_strata: list[pd.DataFrame] = []
    for _, run_frame in retained_responses.groupby("RunId", sort=False):
        for scenario in connectivity_plan["anchor_scenarios"]:
            row, lp, trace, strata = run_retained_case(
                run_frame,
                scenario,
                connectivity_plan,
                current_fit,
                max_configurations=PRIMARY_REPLAY_MAX_CONFIGURATIONS,
                max_constraints=PRIMARY_REPLAY_MAX_CONSTRAINTS,
            )
            retained_rows.append(row)
            _append_table(retained_lp, lp, {})
            _append_table(retained_trace, trace, {})
            _append_table(retained_strata, strata, {})
    retained = pd.DataFrame(retained_rows)

    sentinel_keys = {
        ("chain_strong", 1, "a0_unanchored"),
        ("chain_strong", 1, "a3_random_correct"),
        ("chain_strong", 1, "a6_all_correct"),
        ("hub_strong", 1, "a0_unanchored"),
        ("hub_strong", 1, "a3_random_correct"),
        ("hub_strong", 1, "a6_all_correct"),
    }
    sentinel_rows = []
    for topology, replicate, scenario_name in sorted(sentinel_keys):
        run_frame = retained_responses.loc[
            retained_responses["Topology"].eq(topology)
            & retained_responses["Replicate"].eq(replicate)
        ]
        scenario = next(
            row
            for row in connectivity_plan["anchor_scenarios"]
            if row["scenario"] == scenario_name
        )
        row, _, _, _ = run_retained_case(
            run_frame,
            scenario,
            connectivity_plan,
            current_fit,
            max_configurations=SENTINEL_MAX_CONFIGURATIONS,
            max_constraints=SENTINEL_MAX_CONSTRAINTS,
        )
        sentinel_rows.append(row)
    sentinels = pd.DataFrame(sentinel_rows)

    concat = lambda frames: pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    outputs = {
        "existence_fixture_responses.csv": fixture_responses,
        "existence_fixture_ledger.csv": fixtures,
        "existence_fixture_lp.csv": concat(fixture_lp),
        "existence_fixture_directional_trace.csv": concat(fixture_trace),
        "existence_fixture_strata.csv": concat(fixture_strata),
        "existence_retained_ledger.csv": retained,
        "existence_retained_lp.csv": concat(retained_lp),
        "existence_retained_directional_trace.csv": concat(retained_trace),
        "existence_retained_strata.csv": concat(retained_strata),
        "existence_high_cap_sentinels.csv": sentinels,
    }
    for name, frame in outputs.items():
        write_csv(frame, output / name)

    retained_summary = (
        retained.groupby("ExistenceStatus", dropna=False)
        .agg(
            Cases=("RunId", "size"),
            CurrentReady=("CurrentInferenceReady", "sum"),
            ExistenceQualifiedReady=("ExistenceQualifiedInferenceReady", "sum"),
        )
        .reset_index()
    )
    write_csv(retained_summary, output / "existence_retained_summary.csv")
    outputs["existence_retained_summary.csv"] = retained_summary

    boundary_traces = pd.concat(
        [concat(fixture_trace), concat(retained_trace)], ignore_index=True
    )
    directional_failures = 0
    if not boundary_traces.empty:
        group_columns = [column for column in ["CaseId", "RunId", "Scenario"] if column in boundary_traces.columns]
        for _, group in boundary_traces.groupby(group_columns, dropna=False):
            finite = group.loc[group["Finite"]]
            if len(finite) < 2:
                directional_failures += 1
                continue
            objective = finite.sort_values("Radius")["Objective"].to_numpy(dtype=float)
            tolerance = 1e-10 * max(1.0, float(np.max(np.abs(objective))))
            if np.any(np.diff(objective) > tolerance):
                directional_failures += 1

    fixture_matches = int(fixtures["ExpectedStatusMatch"].sum())
    fixture_contract = fixture_matches == len(fixtures)
    rank_controls = fixtures.loc[
        fixtures["ExpectedStatus"].eq("structural_nonidentification")
    ]
    separation_controls = fixtures.loc[
        fixtures["ExpectedStatus"].eq("boundary_no_finite_cmle")
    ]
    interior_controls = fixtures.loc[
        fixtures["ExpectedStatus"].eq("interior_finite_cmle_supported")
    ]
    fixture_contract = bool(
        fixture_contract
        and rank_controls["ObservedStatus"].eq("structural_nonidentification").all()
        and separation_controls["ObservedStatus"].eq("boundary_no_finite_cmle").all()
        and interior_controls["ObservedStatus"].eq("interior_finite_cmle_supported").all()
    )
    retained_accounting = bool(
        len(retained) == 350
        and int(retained["PrefitEligible"].sum()) == 287
        and int(retained["CurrentInferenceReady"].sum()) == 274
    )
    sentinel_passed = bool(
        len(sentinels) == len(sentinel_keys)
        and sentinels["ExistenceStatus"].ne("unavailable").all()
    )
    contract_passed = bool(
        fixture_contract
        and retained_accounting
        and directional_failures == 0
        and sentinel_passed
    )
    status_counts = retained["ExistenceStatus"].value_counts()
    decision = {
        "schema_version": "mfrm-cmle-finite-mle-existence-decision-v1",
        "study_id": plan["study_id"],
        "contract_passed": contract_passed,
        "plan_sha256": sha256_file(plan_path),
        "core_sha256": sha256_file(ROOT / "mfrm_app/cmle.py"),
        "existence_module_sha256": sha256_file(ROOT / "mfrm_app/cmle_existence.py"),
        "existence_test_sha256": sha256_file(ROOT / "tests/test_cmle_existence.py"),
        "same_retained_response_bytes": sha256_file(CONNECTIVITY_RESPONSES)
        == plan["parent_identity"]["validation/cmle_native_anchor_connectivity_20260810/connectivity_responses.csv"],
        "fixture_cases": len(fixtures),
        "fixture_expected_status_matches": fixture_matches,
        "retained_cases": len(retained),
        "retained_prefit_eligible": int(retained["PrefitEligible"].sum()),
        "retained_boundary": int(status_counts.get("boundary_no_finite_cmle", 0)),
        "retained_interior": int(status_counts.get("interior_finite_cmle_supported", 0)),
        "retained_structural_nonidentification": int(status_counts.get("structural_nonidentification", 0)),
        "retained_unavailable": int(
            len(retained)
            - status_counts.get("boundary_no_finite_cmle", 0)
            - status_counts.get("interior_finite_cmle_supported", 0)
            - status_counts.get("structural_nonidentification", 0)
        ),
        "retained_current_ready": int(retained["CurrentInferenceReady"].sum()),
        "retained_existence_qualified_ready": int(
            retained["ExistenceQualifiedInferenceReady"].sum()
        ),
        "current_ready_withheld_by_boundary": int(
            (
                retained["CurrentInferenceReady"]
                & retained["ExistenceStatus"].eq("boundary_no_finite_cmle")
            ).sum()
        ),
        "current_ready_withheld_by_unavailable": int(
            (
                retained["CurrentInferenceReady"]
                & ~retained["ExistenceStatus"].isin(
                    ["boundary_no_finite_cmle", "interior_finite_cmle_supported", "structural_nonidentification"]
                )
            ).sum()
        ),
        "directional_verification_failures": directional_failures,
        "high_cap_sentinels": len(sentinels),
        "high_cap_sentinels_complete": sentinel_passed,
        "primary_replay_max_configurations": PRIMARY_REPLAY_MAX_CONFIGURATIONS,
        "primary_replay_max_constraints": PRIMARY_REPLAY_MAX_CONSTRAINTS,
        "performance_value_success_gate": None,
        "public_ui_ready": False,
        "output_sha256": {},
    }
    hash_names = list(outputs)
    decision["output_sha256"] = {
        name: sha256_file(output / name) for name in hash_names
    }
    (output / "existence_decision.json").write_text(
        json.dumps(decision, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (output / "CMLE_FINITE_MLE_EXISTENCE_STRESS.md").write_text(
        report_text(decision, retained_summary), encoding="utf-8"
    )
    print(json.dumps(decision, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
