#!/usr/bin/env python3
"""Validate the prospectively registered structured CMLE workflow."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from mfrm_app.cmle_workflow import run_cmle_calibration_workflow  # noqa: E402
from validation.cmle_finite_domain_optimizer_remediation import (  # noqa: E402
    boolean_mismatches,
    numeric_pair_evidence,
)
from validation.cmle_finite_mle_existence_stress import (  # noqa: E402
    CONNECTIVITY_PLAN,
    CONNECTIVITY_RESPONSES,
    fixture_cases,
)
from validation.cmle_native_anchor_connectivity_stress import (  # noqa: E402
    anchor_levels,
    hard_anchors,
    write_csv,
)


DEFAULT_PLAN = ROOT / "validation/cmle_structured_workflow_plan_20260810.json"
DEFAULT_OUTPUT = ROOT / "validation/cmle_structured_workflow_20260810"
INTEGRATION_DIR = ROOT / "validation/cmle_finite_mle_readiness_integration_20260810"
ORACLE_DIR = ROOT / "validation/cmle_finite_mle_oracle_20260810"
FIT_KEY = ["Topology", "Replicate", "RunId", "Scenario"]
STRUCTURAL_KEY = FIT_KEY + ["ParameterType", "Facet", "Level", "StepKey"]
EXPECTED_WORKFLOW = {
    "structural_nonidentification": "design_not_identified",
    "boundary_no_finite_cmle": "finite_mle_boundary",
    "interior_finite_cmle_supported": "calibration_ready",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_plan(path: Path) -> tuple[dict[str, object], dict[str, object]]:
    plan = json.loads(path.read_text(encoding="utf-8"))
    if plan.get("study_id") != "cmle_structured_preoptimization_workflow_v1":
        raise ValueError("Unexpected structured-workflow plan.")
    mismatches = [
        relative
        for relative, digest in plan["parent_identity"].items()
        if not (ROOT / relative).exists()
        or sha256_file(ROOT / relative) != digest
    ]
    integration_decision = json.loads(
        (INTEGRATION_DIR / "integration_decision.json").read_text(encoding="utf-8")
    )
    for name, digest in integration_decision["output_sha256"].items():
        if not (INTEGRATION_DIR / name).exists() or sha256_file(
            INTEGRATION_DIR / name
        ) != digest:
            mismatches.append(f"integration_output:{name}")
    if mismatches:
        raise ValueError(f"Structured-workflow parent identity failed: {mismatches}")
    return plan, json.loads(CONNECTIVITY_PLAN.read_text(encoding="utf-8"))


def _identified(table: pd.DataFrame, identity: dict[str, object]) -> pd.DataFrame:
    if table.empty:
        return table
    value = table.copy()
    for column, item in reversed(list(identity.items())):
        value.insert(0, column, item)
    return value


def _fixture_run(case: dict[str, object]):
    started = time.perf_counter()
    result = run_cmle_calibration_workflow(
        case["frame"],
        **case["prepare_args"],
        gtol=1e-8,
        maxiter=800,
    )
    summary = result["summary"].iloc[0]
    expected = EXPECTED_WORKFLOW[str(case["ExpectedStatus"])]
    identity = {"CaseId": case["CaseId"], "Family": case["Family"]}
    row = {
        **identity,
        "Mechanism": case["Mechanism"],
        "AnchorLabel": case["AnchorLabel"],
        "ExpectedWorkflowStatus": expected,
        "WorkflowStatus": summary["WorkflowStatus"],
        "TerminalMatch": summary["WorkflowStatus"] == expected,
        "OptimizationAttempted": bool(summary["OptimizationAttempted"]),
        "FitReturned": bool(summary["FitReturned"]),
        "Ready": bool(summary["Ready"]),
        "ElapsedSeconds": float(time.perf_counter() - started),
    }
    return row, _identified(result["stages"], identity)


def _retained_run(
    run_frame: pd.DataFrame,
    scenario: dict[str, object],
    connectivity_plan: dict[str, object],
    expected_status: str,
):
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
    result = run_cmle_calibration_workflow(
        analysis,
        person_col="Person",
        facet_cols=["Rater", "Criterion"],
        score_col="ObservedCategory",
        rating_min=0,
        rating_max=3,
        model="RSM",
        hard_anchors=anchors,
        gtol=float(connectivity_plan["fit_contract"]["gtol"]),
        maxiter=int(connectivity_plan["fit_contract"]["maxiter"]),
        newton_polish_maxiter=int(
            connectivity_plan["fit_contract"]["newton_polish_maxiter"]
        ),
    )
    elapsed = time.perf_counter() - started
    summary = result["summary"].iloc[0]
    identity = {
        "Topology": topology,
        "Replicate": replicate,
        "RunId": run_id,
        "Scenario": scenario_name,
    }
    expected_workflow = EXPECTED_WORKFLOW[expected_status]
    fit = result["fit"]
    fit_summary = fit["summary"].iloc[0] if fit is not None else None
    anchor_exact = False
    structural = pd.DataFrame()
    if fit is not None:
        fitted_anchors = fit["facets"]["others"].loc[
            fit["facets"]["others"]["Anchored"]
        ]
        truth = connectivity_plan["simulation"]["rater_truth"]
        anchor_exact = bool(
            len(fitted_anchors) == len(levels)
            and all(
                float(row.Estimate) == float(truth[str(row.Level)])
                and float(row.SE) == 0.0
                for row in fitted_anchors.itertuples(index=False)
            )
        )
        structural = pd.concat(
            [fit["facets"]["others"], fit["steps"]], ignore_index=True
        )
        for column, item in reversed(list(identity.items())):
            structural.insert(0, column, item)
    row = {
        **identity,
        "AnchorCount": len(levels),
        "AnchorLevels": ";".join(levels),
        "ExpectedExistenceStatus": expected_status,
        "ExpectedWorkflowStatus": expected_workflow,
        "WorkflowStatus": summary["WorkflowStatus"],
        "TerminalMatch": summary["WorkflowStatus"] == expected_workflow,
        "PrefitEligible": bool(summary["PrefitEligible"]),
        "PrefitRank": summary["PrefitRank"],
        "PrefitNullity": summary["PrefitNullity"],
        "ExistenceStatus": summary["ExistenceStatus"],
        "OptimizationAttempted": bool(summary["OptimizationAttempted"]),
        "StoppedBeforeOptimization": bool(summary["StoppedBeforeOptimization"]),
        "FitReturned": bool(summary["FitReturned"]),
        "Ready": bool(summary["Ready"]),
        "AnchorExact": anchor_exact,
        "ConditionalLogLik": (
            float(fit_summary["ConditionalLogLik"])
            if fit_summary is not None
            else np.nan
        ),
        "ReadinessReasons": (
            str(fit_summary["ReadinessReasons"])
            if fit_summary is not None
            else ""
        ),
        "ElapsedSeconds": elapsed,
    }
    tables = {
        "stages": _identified(result["stages"], identity),
        "availability": _identified(result["availability"], identity),
    }
    return row, tables, structural


def report_text(decision: dict[str, object]) -> str:
    return f"""# Structured pre-optimization CMLE workflow stress

> Prospectively registered repository-only UX orchestration evidence. No public CMLE selector is enabled.

## Result

- Contract passed: `{decision['contract_passed']}`
- Fixture terminal matches: `{decision['fixture_terminal_matches']}/{decision['fixture_cases']}`
- Retained terminal matches: `{decision['retained_terminal_matches']}/{decision['retained_cases']}`
- Design blocks / boundary stops / calibration ready: `{decision['design_blocks']}` / `{decision['boundary_stops']}` / `{decision['calibration_ready']}`
- Optimization attempts: `{decision['optimization_attempts']}`; forbidden early-stop attempts: `{decision['forbidden_optimizer_attempts']}`
- Structural estimate / SE maximum differences: `{decision['structural_estimate_max_abs_difference']:.3g}` / `{decision['structural_se_max_abs_difference']:.3g}`
- Conditional-loglikelihood maximum difference: `{decision['conditional_loglik_max_abs_difference']:.3g}`
- P95 seconds (design block / boundary / ready): `{decision['design_block_p95_seconds']:.6f}` / `{decision['boundary_stop_p95_seconds']:.6f}` / `{decision['interior_workflow_p95_seconds']:.6f}`
- Selected tests passed: `{decision['selected_tests_passed']}`

## Interpretation

The workflow stops before optimization for structural nonidentification and a certified finite-CMLE boundary, while preserving the direct `fit_cmle` technical-audit route separately. Its bilingual messages and artifact availability are machine-readable but have not yet undergone a user-comprehension study.
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

    fixture_rows = []
    fixture_stages = []
    for case in fixture_cases():
        row, stages = _fixture_run(case)
        fixture_rows.append(row)
        fixture_stages.append(stages)
    fixtures = pd.DataFrame(fixture_rows)
    fixture_stage_table = pd.concat(fixture_stages, ignore_index=True)

    responses = pd.read_csv(CONNECTIVITY_RESPONSES, float_precision="round_trip")
    oracle = pd.read_csv(ORACLE_DIR / "oracle_retained_ledger.csv")
    expected = {
        (str(row.RunId), str(row.Scenario)): str(row.OracleStatus)
        for row in oracle.itertuples(index=False)
    }
    retained_rows = []
    retained_stages = []
    retained_availability = []
    structural_frames = []
    replay_started = time.perf_counter()
    for _, run_frame in responses.groupby("RunId", sort=False):
        run_id = str(run_frame["RunId"].iloc[0])
        for scenario in connectivity_plan["anchor_scenarios"]:
            scenario_name = str(scenario["scenario"])
            row, tables, structural = _retained_run(
                run_frame,
                scenario,
                connectivity_plan,
                expected[(run_id, scenario_name)],
            )
            retained_rows.append(row)
            retained_stages.append(tables["stages"])
            retained_availability.append(tables["availability"])
            if not structural.empty:
                structural_frames.append(structural)
    replay_seconds = time.perf_counter() - replay_started
    retained = pd.DataFrame(retained_rows)
    stages = pd.concat(retained_stages, ignore_index=True)
    availability = pd.concat(retained_availability, ignore_index=True)
    structural = pd.concat(structural_frames, ignore_index=True)

    baseline_fit = pd.read_csv(INTEGRATION_DIR / "integration_fit_ledger.csv")
    baseline_ready = baseline_fit.loc[baseline_fit["InferenceReady"]].copy()
    candidate_fit = retained.loc[retained["FitReturned"]].copy()
    fit_pairs, fit_evidence = numeric_pair_evidence(
        baseline_ready,
        candidate_fit,
        key=FIT_KEY,
        fields=["ConditionalLogLik"],
    )
    fit_boolean_mismatches = boolean_mismatches(
        fit_pairs, ["AnchorExact"]
    )
    baseline_structural = pd.read_csv(
        INTEGRATION_DIR / "integration_structural_results.csv"
    ).merge(
        baseline_ready[FIT_KEY], on=FIT_KEY, how="inner", validate="many_to_one"
    )
    for frame in (baseline_structural, structural):
        frame["StepKey"] = (
            pd.to_numeric(frame["Step"], errors="coerce").fillna(-1).astype(int)
        )
    structural_pairs, structural_evidence = numeric_pair_evidence(
        baseline_structural.drop(columns=["Step"], errors="ignore"),
        structural.drop(columns=["Step"], errors="ignore"),
        key=STRUCTURAL_KEY,
        fields=["Estimate", "SE"],
    )

    test_files = [
        "tests/test_cmle_workflow.py",
        "tests/test_cmle.py",
        "tests/test_cmle_existence.py",
        "tests/test_cmle_hard_anchors.py",
    ]
    completed = subprocess.run(
        [sys.executable, "-m", "pytest", "-q", *test_files],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    (output / "selected_pytest.txt").write_text(
        completed.stdout + completed.stderr, encoding="utf-8"
    )
    tests_passed = completed.returncode == 0

    design_block = retained["WorkflowStatus"].eq("design_not_identified")
    boundary_stop = retained["WorkflowStatus"].eq("finite_mle_boundary")
    ready = retained["WorkflowStatus"].eq("calibration_ready")
    forbidden_attempts = int(
        ((design_block | boundary_stop) & retained["OptimizationAttempted"]).sum()
    )
    stage_schema_passed = bool(
        len(stages) == 4 * len(retained)
        and stages.groupby(FIT_KEY)["Stage"].nunique().eq(4).all()
        and stages["HeadlineEn"].astype(str).str.len().gt(0).all()
        and stages["HeadlineJa"].astype(str).str.len().gt(0).all()
    )
    p95 = lambda mask: float(retained.loc[mask, "ElapsedSeconds"].quantile(0.95))
    design_p95 = p95(design_block)
    boundary_p95 = p95(boundary_stop)
    ready_p95 = p95(ready)
    tolerance = plan["success_gates"]["numeric_invariance"]
    numeric_passed = bool(
        structural_evidence["Estimate_max_abs_difference"]
        <= float(tolerance["structural_estimate_absolute"])
        and structural_evidence["SE_max_abs_difference"]
        <= float(tolerance["structural_se_absolute"])
        and fit_evidence["ConditionalLogLik_max_abs_difference"]
        <= float(tolerance["conditional_loglik_absolute"])
        and all(value == 0 for value in fit_boolean_mismatches.values())
    )
    performance = plan["success_gates"]["performance"]
    performance_passed = bool(
        design_p95 <= float(performance["structural_block_p95_seconds"])
        and boundary_p95 <= float(performance["boundary_stop_p95_seconds"])
        and ready_p95 <= float(performance["interior_workflow_p95_seconds"])
    )
    same_bytes = sha256_file(CONNECTIVITY_RESPONSES) == json.loads(
        (
            ROOT / "validation/cmle_finite_mle_existence_plan_20260810.json"
        ).read_text(encoding="utf-8")
    )["parent_identity"][
        "validation/cmle_native_anchor_connectivity_20260810/connectivity_responses.csv"
    ]
    contract_passed = bool(
        len(fixtures) == 93
        and fixtures["TerminalMatch"].all()
        and len(retained) == 350
        and retained["TerminalMatch"].all()
        and int(design_block.sum()) == 63
        and int(boundary_stop.sum()) == 13
        and int(ready.sum()) == 274
        and int(retained["OptimizationAttempted"].sum()) == 274
        and forbidden_attempts == 0
        and stage_schema_passed
        and numeric_passed
        and performance_passed
        and tests_passed
        and same_bytes
    )

    outputs = {
        "workflow_fixture_ledger.csv": fixtures,
        "workflow_fixture_stages.csv": fixture_stage_table,
        "workflow_retained_ledger.csv": retained,
        "workflow_retained_stages.csv": stages,
        "workflow_retained_availability.csv": availability,
        "workflow_retained_structural.csv": structural.drop(columns=["StepKey"]),
        "workflow_fit_pair_comparison.csv": fit_pairs,
        "workflow_structural_pair_comparison.csv": structural_pairs,
    }
    for name, frame in outputs.items():
        write_csv(frame, output / name)
    comparison = {
        "fit": fit_evidence,
        "fit_boolean_mismatches": fit_boolean_mismatches,
        "structural": structural_evidence,
    }
    (output / "numeric_comparison.json").write_text(
        json.dumps(comparison, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    decision = {
        "schema_version": "mfrm-cmle-structured-workflow-decision-v1",
        "study_id": plan["study_id"],
        "contract_passed": contract_passed,
        "plan_sha256": sha256_file(plan_path),
        "workflow_module_sha256": sha256_file(ROOT / "mfrm_app/cmle_workflow.py"),
        "workflow_test_sha256": sha256_file(ROOT / "tests/test_cmle_workflow.py"),
        "runner_sha256": sha256_file(Path(__file__).resolve()),
        "same_retained_response_bytes": same_bytes,
        "fixture_cases": len(fixtures),
        "fixture_terminal_matches": int(fixtures["TerminalMatch"].sum()),
        "retained_cases": len(retained),
        "retained_terminal_matches": int(retained["TerminalMatch"].sum()),
        "design_blocks": int(design_block.sum()),
        "boundary_stops": int(boundary_stop.sum()),
        "calibration_ready": int(ready.sum()),
        "optimization_attempts": int(retained["OptimizationAttempted"].sum()),
        "forbidden_optimizer_attempts": forbidden_attempts,
        "stage_schema_passed": stage_schema_passed,
        "structural_estimate_max_abs_difference": float(
            structural_evidence["Estimate_max_abs_difference"]
        ),
        "structural_se_max_abs_difference": float(
            structural_evidence["SE_max_abs_difference"]
        ),
        "conditional_loglik_max_abs_difference": float(
            fit_evidence["ConditionalLogLik_max_abs_difference"]
        ),
        "anchor_exact_mismatches": int(sum(fit_boolean_mismatches.values())),
        "design_block_p95_seconds": design_p95,
        "boundary_stop_p95_seconds": boundary_p95,
        "interior_workflow_p95_seconds": ready_p95,
        "performance_contract_passed": performance_passed,
        "replay_seconds": replay_seconds,
        "selected_tests_passed": tests_passed,
        "public_ui_ready": False,
        "output_sha256": {},
    }
    hash_names = [*outputs, "numeric_comparison.json", "selected_pytest.txt"]
    decision["output_sha256"] = {
        name: sha256_file(output / name) for name in hash_names
    }
    (output / "workflow_decision.json").write_text(
        json.dumps(decision, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (output / "CMLE_STRUCTURED_WORKFLOW_STRESS.md").write_text(
        report_text(decision), encoding="utf-8"
    )
    print(json.dumps(decision, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
