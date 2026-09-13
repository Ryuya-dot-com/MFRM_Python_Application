#!/usr/bin/env python3
"""Verify the JMLE Person boundary guard against the frozen 160-run evidence."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import platform
import sys

import numpy as np
import pandas as pd


REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

import streamlit_app as app  # noqa: E402
from validation import operating_characteristics_strict_jmle as stage_a  # noqa: E402


PLAN_SCHEMA = "mfrm-jmle-extreme-score-integration-plan-v1"
DEFAULT_INPUT = REPO / "validation" / "operating_characteristics_pilot20_20260809"
DEFAULT_B2 = REPO / "validation" / "operating_characteristics_strict_jmle_b2_20260809"
DEFAULT_MOVEMENT = REPO / "validation" / "operating_characteristics_jmle_movement_20260809"
DEFAULT_IDENTIFIABILITY = REPO / "validation" / "operating_characteristics_identifiability_20260809"
DEFAULT_PLAN = REPO / "validation" / "jmle_extreme_score_integration_plan_20260809.json"
DEFAULT_OUTPUT = REPO / "validation" / "jmle_extreme_score_integration_20260809"


def as_bool(value) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    return str(value).strip().lower() in {"true", "1", "yes"}


def markdown_table(frame: pd.DataFrame) -> str:
    columns = [str(column) for column in frame.columns]
    rows = [columns] + [
        [str(value).replace("|", "\\|").replace("\n", " ") for value in values]
        for values in frame.itertuples(index=False, name=None)
    ]
    widths = [max(len(row[index]) for row in rows) for index in range(len(columns))]
    rendered = [
        "| " + " | ".join(value.ljust(widths[index]) for index, value in enumerate(rows[0])) + " |",
        "| " + " | ".join("-" * width for width in widths) + " |",
    ]
    rendered.extend(
        "| " + " | ".join(value.ljust(widths[index]) for index, value in enumerate(row)) + " |"
        for row in rows[1:]
    )
    return "\n".join(rendered)


def load_plan(path: Path) -> dict:
    plan = json.loads(path.read_text(encoding="utf-8"))
    if plan.get("schema_version") != PLAN_SCHEMA:
        raise ValueError(f"integration plan schema must be {PLAN_SCHEMA}")
    changes = plan.get("authorized_changes", {})
    if not bool(changes.get("add_person_boundary_audit_before_result_construction", False)):
        raise ValueError("plan must authorize the Person boundary audit")
    for forbidden in (
        "change_optimizer_controls",
        "change_parameter_estimates",
        "automatically_switch_estimator",
        "substitute_finite_extreme_score_correction",
        "authorize_precision_polish_in_application_core",
    ):
        if bool(changes.get(forbidden, True)):
            raise ValueError(f"plan must keep {forbidden}=false")
    if plan.get("application_contract", {}).get("scope_label") != "jmle_person_score_boundary_only":
        raise ValueError("Person boundary audit scope changed")
    return plan


def validate_evidence(
    b2_dir: Path,
    movement_dir: Path,
    ident_dir: Path,
    plan: dict,
) -> pd.DataFrame:
    expected = plan["evidence_identity"]
    paths = {
        "stage_b2_plan_sha256": REPO / "validation" / "operating_characteristics_strict_jmle_b2_plan_20260809.json",
        "stage_b2_runs_sha256": b2_dir / "strict_jmle_b2_runs.csv",
        "stage_b2_identity_sha256": b2_dir / "strict_jmle_b2_identity.json",
        "movement_plan_sha256": REPO / "validation" / "operating_characteristics_jmle_movement_plan_20260809.json",
        "movement_adapter_sha256": REPO / "validation" / "operating_characteristics_jmle_movement_audit.py",
        "movement_run_summary_sha256": movement_dir / "jmle_movement_run_summary.csv",
        "movement_expanded_parameters_sha256": movement_dir / "jmle_movement_expanded_parameters.csv",
        "movement_identity_sha256": movement_dir / "jmle_movement_identity.json",
        "identifiability_plan_sha256": REPO / "validation" / "operating_characteristics_identifiability_plan_20260809.json",
        "identifiability_adapter_sha256": REPO / "validation" / "operating_characteristics_identifiability_audit.py",
        "identifiability_runs_sha256": ident_dir / "jmle_identifiability_runs.csv",
        "identifiability_identity_sha256": ident_dir / "jmle_identifiability_identity.json",
        "pre_change_identifiability_integration_identity_sha256": (
            REPO / "validation" /
            "jmle_extreme_score_prechange_identifiability_identity_20260809.json"
        ),
    }
    rows = []
    for key, path in paths.items():
        actual = stage_a.sha256_file(path)
        target = str(expected[key])
        rows.append({
            "Check": key,
            "Passed": actual == target,
            "Evidence": f"actual={actual}; expected={target}",
        })
    prior_identity_path = paths["pre_change_identifiability_integration_identity_sha256"]
    prior_identity = json.loads(prior_identity_path.read_text(encoding="utf-8"))
    prior_source = str(prior_identity.get("post_change_streamlit_app_sha256", ""))
    rows.append({
        "Check": "pre_change_streamlit_identity_chain",
        "Passed": prior_source == str(expected["pre_change_streamlit_app_sha256"]),
        "Evidence": (
            f"prior integration source={prior_source}; "
            f"expected={expected['pre_change_streamlit_app_sha256']}"
        ),
    })
    checks = pd.DataFrame(rows)
    if not checks["Passed"].all():
        failed = checks.loc[~checks["Passed"], "Check"].astype(str).tolist()
        raise ValueError(f"extreme-score integration evidence identity failed: {failed}")
    return checks


def build_gates(
    comparison: pd.DataFrame,
    run_summary: pd.DataFrame,
    checks: pd.DataFrame,
) -> pd.DataFrame:
    structurally_identified = run_summary["StructurallyIdentified"].astype(bool)
    definitions = [
        ("evidence_identity", len(checks), int(checks["Passed"].sum())),
        ("runids_checked", 160, int(comparison["RunId"].nunique())),
        ("person_rows_checked", len(comparison), int(comparison["BoundaryClassMatched"].sum())),
        ("existing_estimate_retained", len(comparison), int(comparison["EstimateRetained"].sum())),
        ("reportability_contract", len(comparison), int(comparison["ReportabilityContractPassed"].sum())),
        ("person_readiness_contract", len(comparison), int(comparison["PersonReadinessMatched"].sum())),
        ("structurally_identified_runs", 120, int(structurally_identified.sum())),
        (
            "identified_large_non_theta_movement_absent",
            0,
            int(run_summary.loc[structurally_identified, "LargeNonThetaCoordinates"].sum()),
        ),
        (
            "identified_large_nonextreme_theta_movement_absent",
            0,
            int(run_summary.loc[structurally_identified, "LargeNonExtremeThetaCoordinates"].sum()),
        ),
    ]
    rows = [
        {
            "Gate": gate,
            "Required": required,
            "Observed": observed,
            "GatePassed": required == observed,
        }
        for gate, required, observed in definitions
    ]
    rows.append({
        "Gate": "person_boundary_guard_verified",
        "Required": len(rows),
        "Observed": sum(bool(row["GatePassed"]) for row in rows),
        "GatePassed": all(bool(row["GatePassed"]) for row in rows),
    })
    return pd.DataFrame(rows)


def run(
    input_dir: Path,
    b2_dir: Path,
    movement_dir: Path,
    ident_dir: Path,
    output: Path,
    plan_path: Path,
) -> None:
    input_dir = input_dir.resolve()
    b2_dir = b2_dir.resolve()
    movement_dir = movement_dir.resolve()
    ident_dir = ident_dir.resolve()
    output = output.resolve()
    plan_path = plan_path.resolve()
    plan = load_plan(plan_path)
    checks = validate_evidence(b2_dir, movement_dir, ident_dir, plan)

    manifest = pd.read_csv(input_dir / "manifest.csv")
    ratings = pd.read_csv(input_dir / "generated_ratings.csv")
    truth = pd.read_csv(input_dir / "generated_facet_truth.csv")
    anchors = pd.read_csv(input_dir / "generated_anchors.csv")
    b2 = pd.read_csv(b2_dir / "strict_jmle_b2_runs.csv").set_index("RunId", drop=False)
    movement_runs = pd.read_csv(movement_dir / "jmle_movement_run_summary.csv")
    expanded = pd.read_csv(movement_dir / "jmle_movement_expanded_parameters.csv")
    expanded_theta = expanded.loc[expanded["Block"].eq("theta")].copy()
    ident = pd.read_csv(ident_dir / "jmle_identifiability_runs.csv").set_index("RunId", drop=False)

    rows = []
    for _, manifest_row in manifest.iterrows():
        run_id = str(manifest_row["RunId"])
        generated = stage_a.generated_design_for_run(manifest_row, ratings, truth, anchors)
        prep = app.prepare_mfrm_data(
            generated.data,
            person_col="Person",
            facet_cols=["Rater", "Task", "Criterion"],
            score_col="Score",
            rating_min=0,
            rating_max=int(manifest_row["Categories"]) - 1,
            keep_original=True,
        )
        reference = expanded_theta.loc[expanded_theta["RunId"].eq(run_id)].copy()
        reference["Level"] = reference["Level"].astype(str)
        reference = reference.set_index("Level", drop=False)
        person_levels = [str(value) for value in prep["levels"]["Person"]]
        if set(reference.index) != set(person_levels):
            raise ValueError(f"Person levels changed for {run_id}")
        estimates = reference.loc[person_levels, "PolishedValue"].to_numpy(dtype=float)
        structural_ready = as_bool(ident.loc[run_id, "StructurallyIdentified"])
        numerical_ready = as_bool(b2.loc[run_id, "NumericalQualificationPassed"])
        fit_ready = bool(structural_ready and numerical_ready)
        audit = app.build_jmle_person_boundary_audit(
            prep,
            {"method": "JMLE", "n_person": len(person_levels), "n_cat": int(manifest_row["Categories"])},
            estimates=estimates,
            fit_inference_ready=fit_ready,
        )
        current = audit["persons"].set_index("Person", drop=False)
        for person in person_levels:
            cur = current.loc[person]
            ref = reference.loc[person]
            reference_extreme = as_bool(ref["ExtremeScorePattern"])
            reference_minimum = as_bool(ref["ExtremeAllMinimum"])
            reference_maximum = as_bool(ref["ExtremeAllMaximum"])
            expected_direction = (
                "all_minimum" if reference_minimum else
                "all_maximum" if reference_maximum else
                "interior"
            )
            boundary_match = bool(
                bool(cur["ExtremeScorePattern"]) == reference_extreme
                and str(cur["ExtremeScoreDirection"]) == expected_direction
            )
            estimate_retained = bool(float(cur["Estimate"]) == float(ref["PolishedValue"]))
            expected_ready = bool(not reference_extreme and fit_ready and np.isfinite(ref["PolishedValue"]))
            actual_reportable = pd.to_numeric(pd.Series([cur["ReportableEstimate"]]), errors="coerce").iloc[0]
            reportability_pass = bool(
                (expected_ready and actual_reportable == float(ref["PolishedValue"]))
                or (not expected_ready and pd.isna(actual_reportable))
            )
            rows.append({
                "RunId": run_id,
                "ConditionId": str(manifest_row["ConditionId"]),
                "Design": str(manifest_row["Design"]),
                "TruthBias": float(manifest_row["TruthBias"]),
                "Replicate": int(manifest_row["Replicate"]),
                "Person": person,
                "StructurallyIdentified": structural_ready,
                "NumericallyQualified": numerical_ready,
                "FitReadyForPersonBoundary": fit_ready,
                "ReferenceExtremeScorePattern": reference_extreme,
                "ReferenceExtremeScoreDirection": expected_direction,
                "CurrentExtremeScorePattern": bool(cur["ExtremeScorePattern"]),
                "CurrentExtremeScoreDirection": str(cur["ExtremeScoreDirection"]),
                "BoundaryClassMatched": boundary_match,
                "PolishedValue": float(ref["PolishedValue"]),
                "CurrentEstimate": float(cur["Estimate"]),
                "EstimateRetained": estimate_retained,
                "ExpectedPersonInferenceReady": expected_ready,
                "CurrentPersonInferenceReady": bool(cur["PersonInferenceReady"]),
                "PersonReadinessMatched": bool(cur["PersonInferenceReady"]) == expected_ready,
                "CurrentReportableEstimate": actual_reportable,
                "ReportabilityContractPassed": reportability_pass,
                "EstimateRole": str(cur["EstimateRole"]),
                "ObservedCount": int(cur["ObservedCount"]),
                "WeightedScoreTotal": float(cur["WeightedScoreTotal"]),
                "MinimumPossibleWeightedTotal": float(cur["MinimumPossibleWeightedTotal"]),
                "MaximumPossibleWeightedTotal": float(cur["MaximumPossibleWeightedTotal"]),
            })

    comparison = pd.DataFrame(rows)
    ident_flat = ident.reset_index(drop=True)[["RunId", "StructurallyIdentified", "StructuralNullity"]]
    run_summary = movement_runs.merge(ident_flat, on="RunId", validate="one_to_one")
    run_summary["StructurallyIdentified"] = run_summary["StructurallyIdentified"].map(as_bool)
    per_run_boundary = (
        comparison.groupby(["RunId", "ConditionId", "Design", "TruthBias", "Replicate"], sort=False)
        .agg(
            Persons=("Person", "size"),
            ExtremePersons=("CurrentExtremeScorePattern", "sum"),
            BoundaryMatches=("BoundaryClassMatched", "sum"),
            EstimateRetained=("EstimateRetained", "sum"),
            ReportabilityContracts=("ReportabilityContractPassed", "sum"),
            StructurallyIdentified=("StructurallyIdentified", "first"),
            NumericallyQualified=("NumericallyQualified", "first"),
        )
        .reset_index()
    )
    condition_summary = (
        per_run_boundary.groupby(["ConditionId", "Design", "TruthBias"], sort=False)
        .agg(
            Runs=("RunId", "size"),
            Persons=("Persons", "sum"),
            ExtremePersons=("ExtremePersons", "sum"),
            StructurallyIdentifiedRuns=("StructurallyIdentified", "sum"),
            NumericallyQualifiedRuns=("NumericallyQualified", "sum"),
            CompleteBoundaryMatches=("BoundaryMatches", "sum"),
            CompleteReportabilityContracts=("ReportabilityContracts", "sum"),
        )
        .reset_index()
    )
    gate_table = build_gates(comparison, run_summary, checks)
    verified = bool(
        gate_table.loc[gate_table["Gate"].eq("person_boundary_guard_verified"), "GatePassed"].iloc[0]
    )
    identified_runs = run_summary.loc[run_summary["StructurallyIdentified"]]
    evidence_summary = pd.DataFrame([{
        "RetainedRunIds": int(comparison["RunId"].nunique()),
        "RetainedPersonRows": int(len(comparison)),
        "StructurallyIdentifiedRuns": int(run_summary["StructurallyIdentified"].sum()),
        "ExtremePersonRowsInIdentifiedRuns": int(
            comparison.loc[comparison["StructurallyIdentified"], "CurrentExtremeScorePattern"].sum()
        ),
        "IdentifiedRunsWithLargeMovement": int(identified_runs["LargeMovementCoordinates"].gt(0).sum()),
        "IdentifiedLargeNonThetaCoordinates": int(identified_runs["LargeNonThetaCoordinates"].sum()),
        "IdentifiedLargeNonExtremeThetaCoordinates": int(identified_runs["LargeNonExtremeThetaCoordinates"].sum()),
        "IdentifiedMaxAbsNonThetaMovement": float(identified_runs["MaxAbsNonThetaMovement"].max()),
        "IdentifiedMaxAbsNonExtremeThetaMovement": float(identified_runs["MaxAbsNonExtremeThetaMovement"].max()),
        "IdentifiedMaxAbsExtremeThetaMovement": float(identified_runs["MaxAbsExtremeThetaMovement"].max()),
    }])
    decision = pd.DataFrame([{
        "Decision": "JMLEPersonBoundaryGuard",
        "Status": "Verified" if verified else "Blocked",
        "Evidence": (
            f"{int(comparison['BoundaryClassMatched'].sum())}/{len(comparison)} Person rows matched; "
            f"{int(comparison.loc[comparison['StructurallyIdentified'], 'CurrentExtremeScorePattern'].sum())} "
            "extreme Person rows occurred in structurally identified runs."
        ),
        "OptimizerChanged": False,
        "EstimatesChanged": False,
        "AutomaticEstimatorSwitch": False,
        "FiniteExtremeCorrectionSelected": False,
        "PrecisionPolishCoreIntegrationAuthorized": False,
        "PublicPerformanceClaimAuthorized": False,
    }])

    output.mkdir(parents=True, exist_ok=True)
    (output / "registered_jmle_extreme_score_integration_plan.json").write_bytes(plan_path.read_bytes())
    comparison.to_csv(output / "jmle_extreme_score_person_comparison.csv", index=False)
    per_run_boundary.to_csv(output / "jmle_extreme_score_run_summary.csv", index=False)
    condition_summary.to_csv(output / "jmle_extreme_score_condition_summary.csv", index=False)
    evidence_summary.to_csv(output / "jmle_extreme_score_evidence_summary.csv", index=False)
    gate_table.to_csv(output / "jmle_extreme_score_gates.csv", index=False)
    checks.to_csv(output / "jmle_extreme_score_input_checks.csv", index=False)
    decision.to_csv(output / "jmle_extreme_score_decision.csv", index=False)
    report = f"""# JMLE Person score-boundary integration check

## Decision

Status: `{'Verified' if verified else 'Blocked'}`. The application boundary
classification and reportability contract matched
{int(comparison['BoundaryClassMatched'].sum())}/{len(comparison)} retained Person rows.
No optimizer controls or estimates changed, no finite correction or automatic
estimator switch was selected, and precision-polish core integration remains
unauthorized.

## Evidence summary

{markdown_table(evidence_summary)}

## Condition summary

{markdown_table(condition_summary)}

## Interpretation boundary

All-minimum/all-maximum status is defined from retained integer score rows, not
from a terminal-estimate magnitude or rounded display. `Estimate` remains the
technical optimizer/constraint value for reproducibility; `ReportableEstimate`
is withheld for boundary Persons and whenever the fit-level structural and
optimizer gates do not pass. This does not select a finite extreme-score
correction or qualify the precision polish for application-core use.
"""
    (output / "JMLE_EXTREME_SCORE_INTEGRATION.md").write_text(report, encoding="utf-8")
    artifacts = [
        "jmle_extreme_score_person_comparison.csv",
        "jmle_extreme_score_run_summary.csv",
        "jmle_extreme_score_condition_summary.csv",
        "jmle_extreme_score_evidence_summary.csv",
        "jmle_extreme_score_gates.csv",
        "jmle_extreme_score_input_checks.csv",
        "jmle_extreme_score_decision.csv",
        "JMLE_EXTREME_SCORE_INTEGRATION.md",
    ]
    identity = {
        "schema_version": PLAN_SCHEMA,
        "platform": platform.platform(),
        "python": platform.python_version(),
        "verified": verified,
        "optimizer_changed": False,
        "estimates_changed": False,
        "automatic_estimator_switch": False,
        "finite_extreme_correction_selected": False,
        "precision_polish_core_integration_authorized": False,
        "public_performance_claim_authorized": False,
        "integration_plan_sha256": stage_a.sha256_file(plan_path),
        "integration_adapter_sha256": stage_a.sha256_file(Path(__file__).resolve()),
        "post_change_streamlit_app_sha256": stage_a.sha256_file(REPO / "streamlit_app.py"),
        "artifacts": {name: stage_a.sha256_file(output / name) for name in artifacts},
    }
    (output / "jmle_extreme_score_identity.json").write_text(
        json.dumps(identity, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--b2", type=Path, default=DEFAULT_B2)
    parser.add_argument("--movement", type=Path, default=DEFAULT_MOVEMENT)
    parser.add_argument("--identifiability", type=Path, default=DEFAULT_IDENTIFIABILITY)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--plan", type=Path, default=DEFAULT_PLAN)
    args = parser.parse_args()
    run(args.input, args.b2, args.movement, args.identifiability, args.output, args.plan)


if __name__ == "__main__":
    main()
