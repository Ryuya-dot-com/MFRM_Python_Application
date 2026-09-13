#!/usr/bin/env python3
"""Replay retained CMLE fits after the finite-MLE readiness integration."""

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

from validation.cmle_finite_domain_optimizer_remediation import (  # noqa: E402
    boolean_mismatches,
    numeric_pair_evidence,
)
from validation.cmle_native_anchor_connectivity_stress import (  # noqa: E402
    analyze_one,
    realized_graph,
    write_csv,
)


DEFAULT_PLAN = ROOT / "validation/cmle_finite_mle_readiness_integration_plan_20260810.json"
DEFAULT_OUTPUT = ROOT / "validation/cmle_finite_mle_readiness_integration_20260810"
REMEDIATION_DIR = ROOT / "validation/cmle_finite_domain_optimizer_remediation_20260810"
ORACLE_DIR = ROOT / "validation/cmle_finite_mle_oracle_20260810"
CONNECTIVITY_DIR = ROOT / "validation/cmle_native_anchor_connectivity_20260810"
CONNECTIVITY_PLAN = ROOT / "validation/cmle_native_anchor_connectivity_plan_20260810.json"
FIT_KEY = ["Topology", "Replicate", "RunId", "Scenario"]
STRUCTURAL_KEY = FIT_KEY + ["ParameterType", "Facet", "Level", "StepKey"]
PERSON_KEY = FIT_KEY + ["Person"]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_plan(path: Path) -> tuple[dict[str, object], dict[str, object]]:
    plan = json.loads(path.read_text(encoding="utf-8"))
    if plan.get("study_id") != "cmle_finite_mle_readiness_integration_v1":
        raise ValueError("Unexpected finite-MLE readiness integration plan.")
    direct = [
        "validation/cmle_finite_mle_oracle_20260810/oracle_decision.json",
        "validation/cmle_finite_mle_oracle_20260810/oracle_retained_ledger.csv",
    ]
    mismatches = [
        relative
        for relative in direct
        if not (ROOT / relative).exists()
        or sha256_file(ROOT / relative) != plan["parent_identity"][relative]
    ]
    oracle_decision = json.loads(
        (ORACLE_DIR / "oracle_decision.json").read_text(encoding="utf-8")
    )
    remediation_decision = json.loads(
        (REMEDIATION_DIR / "remediation_decision.json").read_text(
            encoding="utf-8"
        )
    )
    recorded = {
        "mfrm_app/cmle.py": oracle_decision["core_sha256"],
        "mfrm_app/cmle_existence.py": oracle_decision["oracle_module_sha256"],
        "tests/test_cmle.py": remediation_decision["candidate_test_sha256"],
        "tests/test_cmle_existence.py": oracle_decision["oracle_test_sha256"],
    }
    mismatches.extend(
        relative
        for relative, digest in recorded.items()
        if digest != plan["parent_identity"][relative]
    )
    for name, digest in remediation_decision["output_sha256"].items():
        if not (REMEDIATION_DIR / name).exists() or sha256_file(
            REMEDIATION_DIR / name
        ) != digest:
            mismatches.append(f"remediation_output:{name}")
    if mismatches:
        raise ValueError(f"Readiness-integration parent identity failed: {mismatches}")
    return plan, json.loads(CONNECTIVITY_PLAN.read_text(encoding="utf-8"))


def report_text(decision: dict[str, object]) -> str:
    return f"""# Finite-CMLE readiness integration stress

> Prospectively registered repository-only integration evidence. The public estimator selector remains unchanged.

## Result

- Contract passed: `{decision['contract_passed']}`
- Audits / eligible fit attempts / returned: `{decision['audits']}` / `{decision['eligible_attempts']}` / `{decision['returned_fits']}`
- Interior / boundary / other existence statuses: `{decision['interior_fits']}` / `{decision['boundary_fits']}` / `{decision['other_existence_statuses']}`
- Inference-ready: `{decision['inference_ready_fits']}`
- Boundary fits withheld with typed reason: `{decision['boundary_typed_reason_fits']}/{decision['boundary_fits']}`
- Oracle ledger status matches: `{decision['oracle_status_matches']}/{decision['eligible_attempts']}`
- Person raw/rounded flag mismatches: `{decision['person_flag_mismatches_total']}`
- Selected regression tests passed: `{decision['selected_tests_passed']}`

## Numeric preservation

- Structural estimate / SE maximum absolute difference: `{decision['structural_estimate_max_abs_difference']:.3g}` / `{decision['structural_se_max_abs_difference']:.3g}`
- Conditional log-likelihood maximum absolute difference: `{decision['conditional_loglik_max_abs_difference']:.3g}`
- WLE estimate / SE maximum absolute difference: `{decision['wle_estimate_max_abs_difference']:.3g}` / `{decision['wle_se_max_abs_difference']:.3g}`
- Infit / Outfit maximum absolute difference: `{decision['infit_max_abs_difference']:.3g}` / `{decision['outfit_max_abs_difference']:.3g}`

## Interpretation

The default native exact-CMLE fit now requires an oracle-supported interior before inference readiness can be true. Technical optimizer values remain available for reproducibility on a boundary, but are explicitly non-ready. This integration does not validate fit thresholds or anchors and does not add a penalty, prior, clipping rule, or estimator switch.
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

    responses_path = CONNECTIVITY_DIR / "connectivity_responses.csv"
    responses = pd.read_csv(responses_path, float_precision="round_trip")
    topologies = {
        str(row["topology"]): row for row in connectivity_plan["topologies"]
    }
    audits: list[dict[str, object]] = []
    fits: list[dict[str, object]] = []
    structural_frames: list[pd.DataFrame] = []
    person_frames: list[pd.DataFrame] = []
    replay_started = time.perf_counter()
    for _, run_frame in responses.groupby("RunId", sort=False):
        topology = topologies[str(run_frame["Topology"].iloc[0])]
        graph = realized_graph(run_frame)
        for scenario in connectivity_plan["anchor_scenarios"]:
            audit, fit, structural, persons = analyze_one(
                run_frame, topology, scenario, connectivity_plan, graph
            )
            audits.append(audit)
            if fit is not None:
                fits.append(fit)
            if not structural.empty:
                structural_frames.append(structural)
            if not persons.empty:
                person_frames.append(persons)
    replay_seconds = time.perf_counter() - replay_started
    candidate_audit = pd.DataFrame(audits).sort_values(
        ["TopologyOrder", "Replicate", "ScenarioOrder"]
    ).reset_index(drop=True)
    candidate_fit = pd.DataFrame(fits).sort_values(
        ["TopologyOrder", "Replicate", "ScenarioOrder"]
    ).reset_index(drop=True)
    candidate_structural = pd.concat(structural_frames, ignore_index=True)
    candidate_person = pd.concat(person_frames, ignore_index=True)

    baseline_fit = pd.read_csv(REMEDIATION_DIR / "candidate_fit_ledger.csv")
    baseline_structural = pd.read_csv(
        REMEDIATION_DIR / "candidate_structural_results.csv"
    )
    baseline_person = pd.read_csv(REMEDIATION_DIR / "candidate_person_results.csv")
    for frame in (candidate_structural, baseline_structural):
        frame["StepKey"] = (
            pd.to_numeric(frame["Step"], errors="coerce").fillna(-1).astype(int)
        )
    fit_pairs, fit_evidence = numeric_pair_evidence(
        baseline_fit,
        candidate_fit,
        key=FIT_KEY,
        fields=[
            "ConditionalLogLik",
            "GradientSupNorm",
            "InformationRank",
            "InformationNullity",
        ],
    )
    fit_boolean_mismatches = boolean_mismatches(
        fit_pairs, ["Returned", "Converged", "InferenceReady", "AnchorExact"]
    )
    structural_pairs, structural_evidence = numeric_pair_evidence(
        baseline_structural.drop(columns=["Step"], errors="ignore"),
        candidate_structural.drop(columns=["Step"], errors="ignore"),
        key=STRUCTURAL_KEY,
        fields=["Estimate", "SE"],
    )
    person_pairs, person_evidence = numeric_pair_evidence(
        baseline_person,
        candidate_person,
        key=PERSON_KEY,
        fields=["WLEEstimate", "WLEStandardError", "Infit", "Outfit"],
    )
    flag_fields = [
        "FitEligibleRaw",
        "EitherUpperRaw",
        "EitherUpperRounded3",
        "RawRounded3Disagreement",
        "EitherUpperRounded6",
        "RawRounded6Disagreement",
    ]
    person_flag_mismatches = boolean_mismatches(person_pairs, flag_fields)

    oracle = pd.read_csv(ORACLE_DIR / "oracle_retained_ledger.csv")[
        ["RunId", "Scenario", "OracleStatus", "ExistenceQualified"]
    ]
    candidate_fit = candidate_fit.merge(
        oracle, on=["RunId", "Scenario"], how="left", validate="one_to_one"
    )
    candidate_fit = candidate_fit.merge(
        baseline_fit[FIT_KEY + ["InferenceReady"]].rename(
            columns={"InferenceReady": "BaselineInferenceReady"}
        ),
        on=FIT_KEY,
        how="left",
        validate="one_to_one",
    )
    candidate_fit["OracleStatusMatch"] = candidate_fit["FiniteMLEStatus"].eq(
        candidate_fit["OracleStatus"]
    )
    candidate_fit["ExpectedIntegratedReadiness"] = (
        candidate_fit["OracleStatus"].eq("interior_finite_cmle_supported")
        & candidate_fit["BaselineInferenceReady"].astype(bool)
    )
    candidate_fit["IntegratedReadinessMatch"] = candidate_fit[
        "InferenceReady"
    ].eq(candidate_fit["ExpectedIntegratedReadiness"])

    test_files = [
        "tests/test_cmle.py",
        "tests/test_cmle_existence.py",
        "tests/test_cmle_hard_anchors.py",
        "tests/test_cmle_person_scoring.py",
        "tests/test_cmle_wle_fit.py",
        "tests/test_cmle_wle_uncertainty.py",
        "tests/test_person_scoring_wle.py",
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

    tolerance = plan["success_gates"]["numeric_invariance"]
    numeric_passed = bool(
        structural_evidence["Estimate_max_abs_difference"]
        <= float(tolerance["structural_estimate_absolute"])
        and structural_evidence["SE_max_abs_difference"]
        <= float(tolerance["structural_se_absolute"])
        and fit_evidence["ConditionalLogLik_max_abs_difference"]
        <= float(tolerance["conditional_loglik_absolute"])
        and person_evidence["WLEEstimate_max_abs_difference"]
        <= float(tolerance["person_wle_absolute"])
        and person_evidence["WLEStandardError_max_abs_difference"]
        <= float(tolerance["person_wle_absolute"])
        and person_evidence["Infit_max_abs_difference"]
        <= float(tolerance["person_fit_absolute"])
        and person_evidence["Outfit_max_abs_difference"]
        <= float(tolerance["person_fit_absolute"])
        and all(value == 0 for value in person_flag_mismatches.values())
    )
    boundary = candidate_fit["FiniteMLEStatus"].eq("boundary_no_finite_cmle")
    interior = candidate_fit["FiniteMLEStatus"].eq(
        "interior_finite_cmle_supported"
    )
    boundary_reason = candidate_fit["ReadinessReasons"].astype(str).str.contains(
        "finite_mle_boundary_no_finite_cmle", regex=False
    )
    other_statuses = int((~(boundary | interior)).sum())
    contract_passed = bool(
        len(candidate_audit) == 350
        and int(candidate_audit["ExactEligible"].sum()) == 287
        and len(candidate_fit) == 287
        and int(candidate_fit["Returned"].sum()) == 287
        and int(interior.sum()) == 274
        and int(boundary.sum()) == 13
        and other_statuses == 0
        and int(candidate_fit["InferenceReady"].sum()) == 274
        and candidate_fit["FiniteMLEGateEnabled"].all()
        and candidate_fit["OracleStatusMatch"].all()
        and candidate_fit["IntegratedReadinessMatch"].all()
        and int((boundary & ~candidate_fit["InferenceReady"] & boundary_reason).sum())
        == 13
        and all(value == 0 for value in fit_boolean_mismatches.values())
        and numeric_passed
        and tests_passed
        and sha256_file(responses_path)
        == json.loads(
            (
                ROOT
                / "validation/cmle_finite_mle_existence_plan_20260810.json"
            ).read_text(encoding="utf-8")
        )["parent_identity"][
            "validation/cmle_native_anchor_connectivity_20260810/connectivity_responses.csv"
        ]
    )

    outputs = {
        "integration_audit_ledger.csv": candidate_audit,
        "integration_fit_ledger.csv": candidate_fit,
        "integration_structural_results.csv": candidate_structural.drop(
            columns=["StepKey"]
        ),
        "integration_person_results.csv": candidate_person,
        "fit_pair_comparison.csv": fit_pairs,
        "structural_pair_comparison.csv": structural_pairs,
        "person_pair_comparison.csv": person_pairs,
    }
    for name, frame in outputs.items():
        write_csv(frame, output / name)
    comparison = {
        "fit": fit_evidence,
        "fit_boolean_mismatches": fit_boolean_mismatches,
        "structural": structural_evidence,
        "person": person_evidence,
        "person_flag_mismatches": person_flag_mismatches,
    }
    (output / "numeric_comparison.json").write_text(
        json.dumps(comparison, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    decision = {
        "schema_version": "mfrm-cmle-finite-mle-readiness-integration-decision-v1",
        "study_id": plan["study_id"],
        "contract_passed": contract_passed,
        "plan_sha256": sha256_file(plan_path),
        "core_sha256": sha256_file(ROOT / "mfrm_app/cmle.py"),
        "existence_module_sha256": sha256_file(
            ROOT / "mfrm_app/cmle_existence.py"
        ),
        "test_sha256": sha256_file(ROOT / "tests/test_cmle.py"),
        "runner_sha256": sha256_file(Path(__file__).resolve()),
        "same_retained_response_bytes": True,
        "audits": len(candidate_audit),
        "eligible_attempts": len(candidate_fit),
        "returned_fits": int(candidate_fit["Returned"].sum()),
        "interior_fits": int(interior.sum()),
        "boundary_fits": int(boundary.sum()),
        "other_existence_statuses": other_statuses,
        "inference_ready_fits": int(candidate_fit["InferenceReady"].sum()),
        "boundary_typed_reason_fits": int(
            (boundary & ~candidate_fit["InferenceReady"] & boundary_reason).sum()
        ),
        "oracle_status_matches": int(candidate_fit["OracleStatusMatch"].sum()),
        "integrated_readiness_matches": int(
            candidate_fit["IntegratedReadinessMatch"].sum()
        ),
        "fit_boolean_mismatches_total": int(
            sum(fit_boolean_mismatches.values())
        ),
        "person_rows": len(candidate_person),
        "person_flag_mismatches_total": int(
            sum(person_flag_mismatches.values())
        ),
        "structural_estimate_max_abs_difference": float(
            structural_evidence["Estimate_max_abs_difference"]
        ),
        "structural_se_max_abs_difference": float(
            structural_evidence["SE_max_abs_difference"]
        ),
        "conditional_loglik_max_abs_difference": float(
            fit_evidence["ConditionalLogLik_max_abs_difference"]
        ),
        "wle_estimate_max_abs_difference": float(
            person_evidence["WLEEstimate_max_abs_difference"]
        ),
        "wle_se_max_abs_difference": float(
            person_evidence["WLEStandardError_max_abs_difference"]
        ),
        "infit_max_abs_difference": float(
            person_evidence["Infit_max_abs_difference"]
        ),
        "outfit_max_abs_difference": float(
            person_evidence["Outfit_max_abs_difference"]
        ),
        "replay_seconds": replay_seconds,
        "selected_tests_passed": tests_passed,
        "public_ui_ready": False,
        "output_sha256": {},
    }
    hash_names = [*outputs, "numeric_comparison.json", "selected_pytest.txt"]
    decision["output_sha256"] = {
        name: sha256_file(output / name) for name in hash_names
    }
    (output / "integration_decision.json").write_text(
        json.dumps(decision, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (output / "CMLE_FINITE_MLE_READINESS_INTEGRATION.md").write_text(
        report_text(decision), encoding="utf-8"
    )
    print(json.dumps(decision, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
