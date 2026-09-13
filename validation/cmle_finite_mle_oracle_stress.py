#!/usr/bin/env python3
"""Validate the prospectively registered scalable finite-CMLE support oracle."""

from __future__ import annotations

import argparse
import hashlib
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

from mfrm_app.cmle import prepare_cmle_design  # noqa: E402
from mfrm_app.cmle_existence import (  # noqa: E402
    _score_configurations,
    audit_cmle_finite_mle_oracle,
    cmle_score_stratum_support_maximum,
)
from validation.cmle_finite_mle_existence_stress import (  # noqa: E402
    CONNECTIVITY_PLAN,
    CONNECTIVITY_RESPONSES,
    REMEDIATION_DIR,
    fixture_cases,
)
from validation.cmle_native_anchor_connectivity_stress import (  # noqa: E402
    anchor_levels,
    hard_anchors,
)


DEFAULT_PLAN = ROOT / "validation/cmle_finite_mle_oracle_plan_20260810.json"
DEFAULT_OUTPUT = ROOT / "validation/cmle_finite_mle_oracle_20260810"
EXHAUSTIVE_DIR = ROOT / "validation/cmle_finite_mle_existence_20260810"
VALUE_TOLERANCE = 1e-10
RANDOM_SEED = 20_260_810


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
    if plan.get("study_id") != "cmle_finite_mle_support_oracle_v1":
        raise ValueError("Unexpected support-oracle plan.")
    direct_parents = [
        "validation/cmle_finite_mle_existence_stress.py",
        "validation/cmle_finite_mle_existence_20260810/existence_decision.json",
        "validation/cmle_finite_mle_existence_20260810/existence_retained_ledger.csv",
    ]
    mismatches = [
        relative
        for relative in direct_parents
        if not (ROOT / relative).exists()
        or sha256_file(ROOT / relative) != plan["parent_identity"][relative]
    ]
    prior_decision = json.loads(
        (EXHAUSTIVE_DIR / "existence_decision.json").read_text(encoding="utf-8")
    )
    if (
        prior_decision["existence_module_sha256"]
        != plan["parent_identity"]["mfrm_app/cmle_existence.py"]
    ):
        mismatches.append("mfrm_app/cmle_existence.py@registered_parent")
    if (
        prior_decision["existence_test_sha256"]
        != plan["parent_identity"]["tests/test_cmle_existence.py"]
    ):
        mismatches.append("tests/test_cmle_existence.py@registered_parent")
    if mismatches:
        raise ValueError(f"Support-oracle parent identity failed: {mismatches}")
    connectivity_plan = json.loads(CONNECTIVITY_PLAN.read_text(encoding="utf-8"))
    return plan, connectivity_plan


def concat(frames: list[pd.DataFrame]) -> pd.DataFrame:
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def identify(table: pd.DataFrame, identity: dict[str, object]) -> pd.DataFrame:
    if table.empty:
        return table
    result = table.copy()
    for column, value in reversed(list(identity.items())):
        result.insert(0, column, value)
    return result


def direction_suite(width: int, case_id: str, pattern_index: int, raw_score: int):
    directions: list[tuple[str, np.ndarray]] = []
    for coordinate in range(width):
        basis = np.zeros(width, dtype=float)
        basis[coordinate] = 1.0
        directions.append((f"coordinate_{coordinate}_positive", basis))
        directions.append((f"coordinate_{coordinate}_negative", -basis))
    alternating = np.where(np.arange(width) % 2 == 0, 1.0, -1.0)
    directions.append(("alternating", alternating))
    identity = f"{case_id}|{pattern_index}|{raw_score}|{RANDOM_SEED}".encode()
    seed = int.from_bytes(hashlib.sha256(identity).digest()[:8], "big")
    rng = np.random.default_rng(seed)
    for index in range(3):
        directions.append((f"seeded_random_{index + 1}", rng.normal(size=width)))
    return directions


def fixture_value_comparisons(
    case: dict[str, object], design
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for pattern_index, pattern in enumerate(design.patterns):
        for raw_score_value in np.flatnonzero(pattern.score_frequencies > 0):
            raw_score = int(raw_score_value)
            configurations = _score_configurations(
                pattern.design.shape[0], design.n_categories, raw_score
            )
            statistics = np.zeros(
                (len(configurations), design.n_parameters), dtype=float
            )
            for unit in range(pattern.design.shape[0]):
                statistics += pattern.design[
                    unit, configurations[:, unit], :
                ]
            for label, direction in direction_suite(
                design.n_parameters,
                str(case["CaseId"]),
                pattern_index,
                raw_score,
            ):
                observed = cmle_score_stratum_support_maximum(
                    design,
                    pattern_index=pattern_index,
                    raw_score=raw_score,
                    direction=direction,
                )
                exhaustive = float(np.max(statistics @ direction))
                oracle = float(observed["Maximum"])
                tolerance = VALUE_TOLERANCE * max(1.0, abs(exhaustive))
                statistic_in_support = bool(
                    np.any(
                        np.all(
                            statistics
                            == np.asarray(observed["Statistic"], dtype=float),
                            axis=1,
                        )
                    )
                )
                rows.append(
                    {
                        "CaseId": case["CaseId"],
                        "Family": case["Family"],
                        "PatternIndex": pattern_index,
                        "RawScore": raw_score,
                        "DirectionLabel": label,
                        "Configurations": len(configurations),
                        "ExhaustiveMaximum": exhaustive,
                        "OracleMaximum": oracle,
                        "AbsoluteDifference": abs(exhaustive - oracle),
                        "Tolerance": tolerance,
                        "ValueMatch": abs(exhaustive - oracle) <= tolerance,
                        "ReturnedStatisticInSupport": statistic_in_support,
                        "ReturnedCategoryScore": int(
                            np.sum(observed["Categories"])
                        ),
                        "ReturnedCategoryScoreMatch": int(
                            np.sum(observed["Categories"])
                        )
                        == raw_score,
                    }
                )
    return rows


def run_fixture(
    case: dict[str, object], exhaustive_status: str
) -> tuple[dict[str, object], dict[str, pd.DataFrame], list[dict[str, object]]]:
    started = time.perf_counter()
    design = prepare_cmle_design(case["frame"], **case["prepare_args"])
    value_rows = (
        fixture_value_comparisons(case, design)
        if bool(design.audit["eligible"])
        else []
    )
    audit = audit_cmle_finite_mle_oracle(design)
    elapsed = time.perf_counter() - started
    summary = audit["summary"].iloc[0]
    row = {
        "CaseId": case["CaseId"],
        "Family": case["Family"],
        "Mechanism": case["Mechanism"],
        "AnchorLabel": case["AnchorLabel"],
        "ExhaustiveStatus": exhaustive_status,
        "OracleStatus": summary["Status"],
        "StatusMatch": summary["Status"] == exhaustive_status,
        "PrefitEligible": bool(design.audit["eligible"]),
        "PrefitRank": int(design.audit["conditional_rank"]),
        "PrefitNullity": int(design.audit["conditional_nullity"]),
        "KParams": int(design.n_parameters),
        "TheoreticalConfigurations": summary.get(
            "TheoreticalConfigurations", np.nan
        ),
        "OracleStateCells": summary.get("OracleStateCells", np.nan),
        "GeneratedConstraintsMax": summary.get(
            "GeneratedConstraintsMax", np.nan
        ),
        "CuttingPlaneRoundsTotal": summary.get(
            "CuttingPlaneRoundsTotal", np.nan
        ),
        "OracleCallsTotal": summary.get("OracleCallsTotal", np.nan),
        "ExistenceQualified": bool(summary.get("ExistenceQualified", False)),
        "ElapsedSeconds": elapsed,
    }
    identity = {"CaseId": case["CaseId"], "Family": case["Family"]}
    tables = {
        "lp": identify(audit["lp_tolerances"], identity),
        "cuts": identify(audit["cut_history"], identity),
        "strata": identify(audit["strata"], identity),
        "trace": identify(audit["directional_trace"], identity),
    }
    return row, tables, value_rows


def run_retained(
    run_frame: pd.DataFrame,
    scenario: dict[str, object],
    connectivity_plan: dict[str, object],
    current_fit: pd.DataFrame,
) -> tuple[dict[str, object], dict[str, pd.DataFrame]]:
    topology = str(run_frame["Topology"].iloc[0])
    replicate = int(run_frame["Replicate"].iloc[0])
    run_id = str(run_frame["RunId"].iloc[0])
    scenario_name = str(scenario["scenario"])
    levels = anchor_levels(connectivity_plan, replicate, scenario_name)
    anchors = hard_anchors(connectivity_plan, levels)
    analysis = run_frame[
        ["Person", "Rater", "Criterion", "ObservedCategory"]
    ].copy()
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
    audit = audit_cmle_finite_mle_oracle(design)
    elapsed = time.perf_counter() - started
    summary = audit["summary"].iloc[0]
    fitted = current_fit.loc[
        current_fit["RunId"].eq(run_id)
        & current_fit["Scenario"].eq(scenario_name)
    ]
    current_returned = bool(fitted.iloc[0]["Returned"]) if len(fitted) else False
    current_ready = bool(fitted.iloc[0]["InferenceReady"]) if len(fitted) else False
    identity = {
        "Topology": topology,
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
        "OracleStatus": summary["Status"],
        "OracleReason": summary["Reason"],
        "OracleComplete": bool(summary.get("OracleComplete", False)),
        "BoundaryDetected": summary.get("BoundaryDetected", pd.NA),
        "ExistenceQualified": bool(summary.get("ExistenceQualified", False)),
        "TheoreticalConfigurations": summary.get(
            "TheoreticalConfigurations", np.nan
        ),
        "OracleStateCells": summary.get("OracleStateCells", np.nan),
        "GeneratedConstraintsMax": summary.get(
            "GeneratedConstraintsMax", np.nan
        ),
        "CuttingPlaneRoundsTotal": summary.get(
            "CuttingPlaneRoundsTotal", np.nan
        ),
        "OracleCallsTotal": summary.get("OracleCallsTotal", np.nan),
        "DirectionalObjectiveNonincreasing": summary.get(
            "DirectionalObjectiveNonincreasing", pd.NA
        ),
        "CurrentFitReturned": current_returned,
        "CurrentInferenceReady": current_ready,
        "ExistenceQualifiedInferenceReady": bool(
            current_ready and summary.get("ExistenceQualified", False)
        ),
        "ElapsedSeconds": elapsed,
    }
    tables = {
        "lp": identify(audit["lp_tolerances"], identity),
        "cuts": identify(audit["cut_history"], identity),
        "strata": identify(audit["strata"], identity),
        "trace": identify(audit["directional_trace"], identity),
    }
    return row, tables


def report_text(decision: dict[str, object], transition: pd.DataFrame) -> str:
    lines = [
        "| Oracle status | Cases | Current ready | Existence-qualified ready |",
        "| --- | ---: | ---: | ---: |",
    ]
    for row in transition.itertuples(index=False):
        lines.append(
            f"| {row.OracleStatus} | {int(row.Cases)} | {int(row.CurrentReady)} | {int(row.ExistenceQualifiedReady)} |"
        )
    return f"""# Scalable finite-CMLE support-oracle stress

> Prospectively registered repository-only research evidence. Public CMLE routing remains withheld.

## Result

- Contract passed: `{decision['contract_passed']}`
- Oracle/exhaustive value comparisons: `{decision['oracle_value_matches']}/{decision['oracle_value_comparisons']}`
- Fixture status matches: `{decision['fixture_status_matches']}/{decision['fixture_cases']}`
- Retained exact-status matches: `{decision['retained_exact_status_matches']}/{decision['retained_exact_status_comparisons']}`
- High-cap sentinel matches: `{decision['high_cap_sentinel_matches']}/{decision['high_cap_sentinels']}`
- Retained structural / boundary / interior / unavailable: `{decision['retained_structural']}` / `{decision['retained_boundary']}` / `{decision['retained_interior']}` / `{decision['retained_unavailable']}`
- Current inference-ready / existence-qualified ready: `{decision['retained_current_ready']}` / `{decision['retained_existence_qualified_ready']}`
- Eligible audit seconds P95 / max; all-350 total: `{decision['eligible_case_p95_seconds']:.6f}` / `{decision['eligible_case_max_seconds']:.6f}`; `{decision['all_350_total_seconds']:.6f}`
- Directional verification failures: `{decision['directional_verification_failures']}`

## Same-byte transition

{chr(10).join(lines)}

## Interpretation

The support oracle changes the computational route, not the likelihood or finite-existence criterion. It removes exhaustive configuration materialization by adding only oracle-certified support cuts. An interior result remains conditional on the implemented model and anchors; it is not evidence of fit, fairness, anchor validity, or cross-engine parity.
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
            raise FileExistsError(
                f"Output is non-empty: {output}; pass --overwrite explicitly."
            )
        shutil.rmtree(output)
    output.mkdir(parents=True, exist_ok=True)
    shutil.copy2(plan_path, output / plan_path.name)

    exhaustive_fixtures = pd.read_csv(
        EXHAUSTIVE_DIR / "existence_fixture_ledger.csv"
    ).set_index("CaseId")
    fixture_rows: list[dict[str, object]] = []
    value_rows: list[dict[str, object]] = []
    fixture_tables = {name: [] for name in ("lp", "cuts", "strata", "trace")}
    for case in fixture_cases():
        exhaustive_status = str(
            exhaustive_fixtures.loc[case["CaseId"], "ObservedStatus"]
        )
        row, tables, comparisons = run_fixture(case, exhaustive_status)
        fixture_rows.append(row)
        value_rows.extend(comparisons)
        for name, table in tables.items():
            if not table.empty:
                fixture_tables[name].append(table)
    fixtures = pd.DataFrame(fixture_rows)
    values = pd.DataFrame(value_rows)

    retained_responses = pd.read_csv(
        CONNECTIVITY_RESPONSES, float_precision="round_trip"
    )
    current_fit = pd.read_csv(REMEDIATION_DIR / "candidate_fit_ledger.csv")
    retained_rows: list[dict[str, object]] = []
    retained_tables = {name: [] for name in ("lp", "cuts", "strata", "trace")}
    retained_started = time.perf_counter()
    for _, run_frame in retained_responses.groupby("RunId", sort=False):
        for scenario in connectivity_plan["anchor_scenarios"]:
            row, tables = run_retained(
                run_frame, scenario, connectivity_plan, current_fit
            )
            retained_rows.append(row)
            for name, table in tables.items():
                if not table.empty:
                    retained_tables[name].append(table)
    retained_total_seconds = time.perf_counter() - retained_started
    retained = pd.DataFrame(retained_rows)

    primary_exact = pd.read_csv(
        EXHAUSTIVE_DIR / "existence_retained_ledger.csv"
    )[
        ["RunId", "Scenario", "ExistenceStatus"]
    ].rename(columns={"ExistenceStatus": "PrimaryExhaustiveStatus"})
    retained = retained.merge(primary_exact, on=["RunId", "Scenario"], how="left")
    retained["PrimaryExhaustiveComparable"] = retained[
        "PrimaryExhaustiveStatus"
    ].ne("unavailable")
    retained["PrimaryExhaustiveStatusMatch"] = (
        retained["OracleStatus"].eq(retained["PrimaryExhaustiveStatus"])
        & retained["PrimaryExhaustiveComparable"]
    )
    high_cap = pd.read_csv(EXHAUSTIVE_DIR / "existence_high_cap_sentinels.csv")[
        ["Topology", "Replicate", "Scenario", "ExistenceStatus"]
    ].rename(columns={"ExistenceStatus": "HighCapExhaustiveStatus"})
    retained = retained.merge(
        high_cap, on=["Topology", "Replicate", "Scenario"], how="left"
    )
    retained["HighCapSentinel"] = retained["HighCapExhaustiveStatus"].notna()
    retained["HighCapStatusMatch"] = (
        retained["OracleStatus"].eq(retained["HighCapExhaustiveStatus"])
        & retained["HighCapSentinel"]
    )

    transition = (
        retained.groupby("OracleStatus", dropna=False)
        .agg(
            Cases=("RunId", "size"),
            CurrentReady=("CurrentInferenceReady", "sum"),
            ExistenceQualifiedReady=("ExistenceQualifiedInferenceReady", "sum"),
        )
        .reset_index()
    )
    status_counts = retained["OracleStatus"].value_counts()
    eligible_times = retained.loc[retained["PrefitEligible"], "ElapsedSeconds"]
    directional_failures = int(
        (
            retained["OracleStatus"].eq("boundary_no_finite_cmle")
            & ~retained["DirectionalObjectiveNonincreasing"].eq(True)
        ).sum()
    )
    value_matches = int(
        (
            values["ValueMatch"]
            & values["ReturnedStatisticInSupport"]
            & values["ReturnedCategoryScoreMatch"]
        ).sum()
    )
    fixture_matches = int(fixtures["StatusMatch"].sum())
    primary_comparisons = int(retained["PrimaryExhaustiveComparable"].sum())
    primary_matches = int(retained["PrimaryExhaustiveStatusMatch"].sum())
    sentinel_count = int(retained["HighCapSentinel"].sum())
    sentinel_matches = int(retained["HighCapStatusMatch"].sum())
    p95_seconds = float(eligible_times.quantile(0.95))
    max_seconds = float(eligible_times.max())
    tolerance_unstable = int(
        retained["OracleStatus"].eq("tolerance_unstable").sum()
    )
    unavailable = int(
        (
            ~retained["OracleStatus"].isin(
                [
                    "structural_nonidentification",
                    "boundary_no_finite_cmle",
                    "interior_finite_cmle_supported",
                ]
            )
        ).sum()
    )
    performance_passed = bool(
        p95_seconds <= 2.0
        and max_seconds <= 10.0
        and retained_total_seconds <= 120.0
    )
    contract_passed = bool(
        value_matches == len(values)
        and fixture_matches == len(fixtures)
        and primary_matches == primary_comparisons
        and sentinel_matches == sentinel_count == 6
        and len(retained) == 350
        and int(retained["PrefitEligible"].sum()) == 287
        and directional_failures == 0
        and tolerance_unstable == 0
        and unavailable == 0
        and performance_passed
    )

    outputs = {
        "oracle_value_equivalence.csv": values,
        "oracle_fixture_ledger.csv": fixtures,
        "oracle_fixture_lp.csv": concat(fixture_tables["lp"]),
        "oracle_fixture_cut_history.csv": concat(fixture_tables["cuts"]),
        "oracle_fixture_strata.csv": concat(fixture_tables["strata"]),
        "oracle_fixture_directional_trace.csv": concat(fixture_tables["trace"]),
        "oracle_retained_ledger.csv": retained,
        "oracle_retained_lp.csv": concat(retained_tables["lp"]),
        "oracle_retained_cut_history.csv": concat(retained_tables["cuts"]),
        "oracle_retained_strata.csv": concat(retained_tables["strata"]),
        "oracle_retained_directional_trace.csv": concat(retained_tables["trace"]),
        "oracle_retained_transition.csv": transition,
    }
    for name, frame in outputs.items():
        write_csv(frame, output / name)

    decision = {
        "schema_version": "mfrm-cmle-finite-mle-support-oracle-decision-v1",
        "study_id": plan["study_id"],
        "contract_passed": contract_passed,
        "plan_sha256": sha256_file(plan_path),
        "core_sha256": sha256_file(ROOT / "mfrm_app/cmle.py"),
        "oracle_module_sha256": sha256_file(ROOT / "mfrm_app/cmle_existence.py"),
        "oracle_test_sha256": sha256_file(ROOT / "tests/test_cmle_existence.py"),
        "runner_sha256": sha256_file(Path(__file__).resolve()),
        "same_retained_response_bytes": sha256_file(CONNECTIVITY_RESPONSES)
        == json.loads(
            (
                ROOT
                / "validation/cmle_finite_mle_existence_plan_20260810.json"
            ).read_text(encoding="utf-8")
        )["parent_identity"][
            "validation/cmle_native_anchor_connectivity_20260810/connectivity_responses.csv"
        ],
        "oracle_value_comparisons": len(values),
        "oracle_value_matches": value_matches,
        "fixture_cases": len(fixtures),
        "fixture_status_matches": fixture_matches,
        "retained_cases": len(retained),
        "retained_prefit_eligible": int(retained["PrefitEligible"].sum()),
        "retained_exact_status_comparisons": primary_comparisons,
        "retained_exact_status_matches": primary_matches,
        "high_cap_sentinels": sentinel_count,
        "high_cap_sentinel_matches": sentinel_matches,
        "retained_structural": int(
            status_counts.get("structural_nonidentification", 0)
        ),
        "retained_boundary": int(
            status_counts.get("boundary_no_finite_cmle", 0)
        ),
        "retained_interior": int(
            status_counts.get("interior_finite_cmle_supported", 0)
        ),
        "retained_unavailable": unavailable,
        "retained_current_ready": int(retained["CurrentInferenceReady"].sum()),
        "retained_existence_qualified_ready": int(
            retained["ExistenceQualifiedInferenceReady"].sum()
        ),
        "directional_verification_failures": directional_failures,
        "tolerance_unstable": tolerance_unstable,
        "eligible_case_p95_seconds": p95_seconds,
        "eligible_case_max_seconds": max_seconds,
        "all_350_total_seconds": retained_total_seconds,
        "performance_contract_passed": performance_passed,
        "public_ui_ready": False,
        "output_sha256": {
            name: sha256_file(output / name) for name in outputs
        },
    }
    (output / "oracle_decision.json").write_text(
        json.dumps(decision, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (output / "CMLE_FINITE_MLE_SUPPORT_ORACLE_STRESS.md").write_text(
        report_text(decision, transition), encoding="utf-8"
    )
    print(json.dumps(decision, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
