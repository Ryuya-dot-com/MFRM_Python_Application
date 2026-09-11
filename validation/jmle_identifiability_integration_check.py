#!/usr/bin/env python3
"""Verify the integrated JMLE eta-rank guard against all 160 retained audits."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import platform
import sys

import pandas as pd


REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

import streamlit_app as app  # noqa: E402
from validation import operating_characteristics_strict_jmle as stage_a  # noqa: E402


PLAN_SCHEMA = "mfrm-jmle-identifiability-integration-plan-v1"
DEFAULT_INPUT = REPO / "validation" / "operating_characteristics_pilot20_20260809"
DEFAULT_REFERENCE = REPO / "validation" / "operating_characteristics_identifiability_20260809"
DEFAULT_PLAN = REPO / "validation" / "jmle_identifiability_integration_plan_20260809.json"
DEFAULT_OUTPUT = REPO / "validation" / "jmle_identifiability_integration_20260809"


def markdown_table(frame: pd.DataFrame) -> str:
    """Render a small evidence table without pandas' optional tabulate dependency."""
    columns = [str(column) for column in frame.columns]
    rows = [columns]
    for values in frame.itertuples(index=False, name=None):
        rows.append([
            str(value).replace("|", "\\|").replace("\n", " ")
            for value in values
        ])
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
    if not bool(changes.get("add_pre_optimizer_eta_design_rank_audit", False)):
        raise ValueError("integration plan must authorize the eta-rank audit")
    if bool(changes.get("change_optimizer_controls", True)):
        raise ValueError("integration plan must not change optimizer controls")
    if bool(changes.get("automatically_switch_jmle_to_mml", True)):
        raise ValueError("integration plan must not auto-switch estimators")
    if bool(changes.get("change_parameter_estimates", True)):
        raise ValueError("integration plan must not authorize estimate changes")
    if plan.get("contracts", {}).get("scope_label") != "eta_person_facet_blocks_only":
        raise ValueError("integration audit scope changed")
    return plan


def validate_evidence(reference_dir: Path, plan: dict) -> pd.DataFrame:
    expected = plan["evidence_identity"]
    reference_identity = json.loads(
        (reference_dir / "jmle_identifiability_identity.json").read_text(encoding="utf-8")
    )
    paths = {
        "identifiability_plan_sha256": REPO / "validation" / "operating_characteristics_identifiability_plan_20260809.json",
        "identifiability_adapter_sha256": REPO / "validation" / "operating_characteristics_identifiability_audit.py",
        "identifiability_runs_sha256": reference_dir / "jmle_identifiability_runs.csv",
        "identifiability_identity_sha256": reference_dir / "jmle_identifiability_identity.json",
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
    rows.append({
        "Check": "pre_change_streamlit_identity",
        "Passed": reference_identity.get("streamlit_app_sha256") == expected["pre_change_streamlit_app_sha256"],
        "Evidence": (
            f"reference={reference_identity.get('streamlit_app_sha256')}; "
            f"expected={expected['pre_change_streamlit_app_sha256']}"
        ),
    })
    checks = pd.DataFrame(rows)
    if not checks["Passed"].all():
        failed = checks.loc[~checks["Passed"], "Check"].astype(str).tolist()
        raise ValueError(f"integration evidence identity failed: {failed}")
    return checks


def integrated_audit(manifest_row: pd.Series, generated) -> dict[str, object]:
    prep = app.prepare_mfrm_data(
        generated.data,
        person_col="Person",
        facet_cols=["Rater", "Task", "Criterion"],
        score_col="Score",
        rating_min=0,
        rating_max=int(manifest_row["Categories"]) - 1,
        keep_original=True,
    )
    specs = app.prepare_constraint_specs(
        prep,
        anchor_df=generated.anchors if not generated.anchors.empty else None,
        noncenter_facet="Person",
    )
    config = {
        "model": "RSM",
        "method": "JMLE",
        "n_cat": int(manifest_row["Categories"]),
        "n_person": len(prep["levels"]["Person"]),
        "facet_names": prep["facet_names"],
        "facet_levels": {facet: prep["levels"][facet] for facet in prep["facet_names"]},
        "facet_signs": {"Rater": -1, "Task": -1, "Criterion": -1},
        "theta_spec": specs["theta_spec"],
        "facet_specs": specs["facet_specs"],
        "population_model": {"enabled": False, "n_params": 0},
    }
    audit = app.build_jmle_eta_identifiability_audit(prep, config)
    summary = audit["summary"].iloc[0]
    rater = audit["connectivity"].loc[
        audit["connectivity"]["Facet"].eq("Rater")
    ].iloc[0]
    return {
        "Status": str(audit["status"]),
        "Scope": str(audit["scope"]),
        "EtaDesignRank": int(summary["EtaDesignRank"]),
        "EtaDesignColumns": int(summary["EtaDesignColumns"]),
        "EtaStructuralNullity": int(summary["EtaStructuralNullity"]),
        "EtaStructurallyIdentified": bool(summary["EtaStructurallyIdentified"]),
        "InferenceReady": bool(summary["InferenceReady"]),
        "RaterComponents": int(rater["Components"]),
        "MinimumRatersPerPerson": int(rater["MinimumRatersOrLevelsPerPerson"]),
    }


def gates(comparison: pd.DataFrame, checks: pd.DataFrame) -> pd.DataFrame:
    expected = 160
    definitions = [
        ("evidence_identity", len(checks), int(checks["Passed"].sum())),
        ("runids_checked", expected, len(comparison)),
        ("rank_match", expected, int(comparison["RankMatched"].sum())),
        ("nullity_match", expected, int(comparison["NullityMatched"].sum())),
        ("connectivity_match", expected, int(comparison["ConnectivityMatched"].sum())),
        ("scope_match", expected, int(comparison["ScopeMatched"].sum())),
        ("readiness_match", expected, int(comparison["ReadinessMatched"].sum())),
        ("complete_contract_match", expected, int(comparison["CompleteContractMatched"].sum())),
    ]
    rows = [
        {"Gate": gate, "Required": required, "Passed": passed, "GatePassed": required == passed}
        for gate, required, passed in definitions
    ]
    rows.append({
        "Gate": "guard_integration_verified",
        "Required": len(rows),
        "Passed": sum(bool(row["GatePassed"]) for row in rows),
        "GatePassed": all(bool(row["GatePassed"]) for row in rows),
    })
    return pd.DataFrame(rows)


def run(input_dir: Path, reference_dir: Path, output: Path, plan_path: Path) -> None:
    input_dir = input_dir.resolve()
    reference_dir = reference_dir.resolve()
    output = output.resolve()
    plan_path = plan_path.resolve()
    plan = load_plan(plan_path)
    checks = validate_evidence(reference_dir, plan)
    manifest = pd.read_csv(input_dir / "manifest.csv")
    ratings = pd.read_csv(input_dir / "generated_ratings.csv")
    truth = pd.read_csv(input_dir / "generated_facet_truth.csv")
    anchors = pd.read_csv(input_dir / "generated_anchors.csv")
    reference = pd.read_csv(reference_dir / "jmle_identifiability_runs.csv").set_index("RunId", drop=False)
    rows = []
    for _, manifest_row in manifest.iterrows():
        run_id = str(manifest_row["RunId"])
        generated = stage_a.generated_design_for_run(manifest_row, ratings, truth, anchors)
        current = integrated_audit(manifest_row, generated)
        retained = reference.loc[run_id]
        rank_match = current["EtaDesignRank"] == int(retained["EtaDesignRank"])
        nullity_match = current["EtaStructuralNullity"] == int(retained["StructuralNullity"])
        connectivity_match = current["RaterComponents"] == int(retained["PersonRaterComponents"])
        scope_match = current["Scope"] == str(plan["contracts"]["scope_label"])
        readiness_match = current["InferenceReady"] == bool(retained["StructurallyIdentified"])
        rows.append({
            "RunId": run_id,
            "ConditionId": str(manifest_row["ConditionId"]),
            "Design": str(manifest_row["Design"]),
            "TruthBias": float(manifest_row["TruthBias"]),
            **current,
            "ReferenceEtaDesignRank": int(retained["EtaDesignRank"]),
            "ReferenceStructuralNullity": int(retained["StructuralNullity"]),
            "ReferenceRaterComponents": int(retained["PersonRaterComponents"]),
            "RankMatched": rank_match,
            "NullityMatched": nullity_match,
            "ConnectivityMatched": connectivity_match,
            "ScopeMatched": scope_match,
            "ReadinessMatched": readiness_match,
            "CompleteContractMatched": bool(
                rank_match and nullity_match and connectivity_match and scope_match and readiness_match
            ),
        })
    comparison = pd.DataFrame(rows)
    gate_table = gates(comparison, checks)
    verified = bool(
        gate_table.loc[gate_table["Gate"].eq("guard_integration_verified"), "GatePassed"].iloc[0]
    )
    condition_summary = (
        comparison.groupby(["ConditionId", "Design", "TruthBias"], sort=False)
        .agg(
            Runs=("RunId", "size"),
            ContractMatches=("CompleteContractMatched", "sum"),
            EtaStructuralNullity=("EtaStructuralNullity", "first"),
            RaterComponents=("RaterComponents", "first"),
            InferenceReady=("InferenceReady", "all"),
        )
        .reset_index()
    )
    decision = pd.DataFrame([{
        "Decision": "IdentifiabilityGuardIntegration",
        "Status": "Verified" if verified else "Blocked",
        "Evidence": f"{int(comparison['CompleteContractMatched'].sum())}/{len(comparison)} retained RunIds matched.",
        "OptimizerChanged": False,
        "AutomaticEstimatorSwitch": False,
        "ApplicationCoreChangeAuthorized": verified,
        "PublicPerformanceClaimAuthorized": False,
    }])
    output.mkdir(parents=True, exist_ok=True)
    (output / "registered_jmle_identifiability_integration_plan.json").write_bytes(plan_path.read_bytes())
    comparison.to_csv(output / "jmle_identifiability_integration_comparison.csv", index=False)
    condition_summary.to_csv(output / "jmle_identifiability_integration_conditions.csv", index=False)
    gate_table.to_csv(output / "jmle_identifiability_integration_gates.csv", index=False)
    checks.to_csv(output / "jmle_identifiability_integration_input_checks.csv", index=False)
    decision.to_csv(output / "jmle_identifiability_integration_decision.csv", index=False)
    report = f"""# JMLE identifiability-guard integration check

## Decision

Status: `{'Verified' if verified else 'Blocked'}`. The integrated app audit
matched {int(comparison['CompleteContractMatched'].sum())}/{len(comparison)}
retained RunIds for exact eta rank, nullity, Person-Rater connectivity, scope,
and readiness. Optimizer controls and estimates were not changed, and no
automatic JMLE-to-MML switch or public performance claim is authorized.

## Condition summary

{markdown_table(condition_summary)}

## Retained artifacts

- `jmle_identifiability_integration_comparison.csv`
- `jmle_identifiability_integration_conditions.csv`
- `jmle_identifiability_integration_gates.csv`
- `jmle_identifiability_integration_input_checks.csv`
- `jmle_identifiability_integration_decision.csv`
- `jmle_identifiability_integration_identity.json`
"""
    (output / "JMLE_IDENTIFIABILITY_INTEGRATION.md").write_text(report, encoding="utf-8")
    artifact_names = [
        "jmle_identifiability_integration_comparison.csv",
        "jmle_identifiability_integration_conditions.csv",
        "jmle_identifiability_integration_gates.csv",
        "jmle_identifiability_integration_input_checks.csv",
        "jmle_identifiability_integration_decision.csv",
        "JMLE_IDENTIFIABILITY_INTEGRATION.md",
    ]
    identity = {
        "schema_version": PLAN_SCHEMA,
        "platform": platform.platform(),
        "python": platform.python_version(),
        "verified": verified,
        "optimizer_changed": False,
        "automatic_estimator_switch": False,
        "public_performance_claim_authorized": False,
        "integration_plan_sha256": stage_a.sha256_file(plan_path),
        "integration_adapter_sha256": stage_a.sha256_file(Path(__file__).resolve()),
        "post_change_streamlit_app_sha256": stage_a.sha256_file(REPO / "streamlit_app.py"),
        "artifacts": {name: stage_a.sha256_file(output / name) for name in artifact_names},
    }
    (output / "jmle_identifiability_integration_identity.json").write_text(
        json.dumps(identity, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--plan", type=Path, default=DEFAULT_PLAN)
    args = parser.parse_args()
    run(args.input, args.reference, args.output, args.plan)


if __name__ == "__main__":
    main()
