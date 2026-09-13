#!/usr/bin/env python3
"""Replay retained connectivity bytes after the finite-domain optimizer guard."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from validation.cmle_native_anchor_connectivity_stress import (  # noqa: E402
    analyze_one,
    realized_graph,
    write_csv,
)


DEFAULT_PLAN = ROOT / "validation/cmle_finite_domain_optimizer_remediation_plan_20260810.json"
DEFAULT_OUTPUT = ROOT / "validation/cmle_finite_domain_optimizer_remediation_20260810"
BASELINE_DIR = ROOT / "validation/cmle_native_anchor_connectivity_20260810"
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


def load_and_validate(plan_path: Path) -> tuple[dict[str, object], dict[str, object]]:
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    if plan.get("study_id") != "cmle_finite_domain_optimizer_remediation_v1":
        raise ValueError("Unexpected finite-domain remediation plan.")
    immutable = {
        relative: digest
        for relative, digest in plan["baseline_identity"].items()
        if relative not in {
            "mfrm_app/cmle.py",
            "tests/test_cmle.py",
            "validation/cmle_native_anchor_connectivity_stress.py",
        }
    }
    mismatches = [
        relative
        for relative, digest in immutable.items()
        if not (ROOT / relative).exists() or sha256_file(ROOT / relative) != digest
    ]
    if mismatches:
        raise ValueError(f"Immutable remediation baseline identity failed: {mismatches}")
    baseline_decision = json.loads(
        (BASELINE_DIR / "connectivity_decision.json").read_text(encoding="utf-8")
    )
    output_mismatches = [
        name
        for name, digest in baseline_decision["output_sha256"].items()
        if not (BASELINE_DIR / name).exists()
        or sha256_file(BASELINE_DIR / name) != digest
    ]
    if output_mismatches:
        raise ValueError(f"Baseline evidence output identity failed: {output_mismatches}")
    connectivity_plan = json.loads(CONNECTIVITY_PLAN.read_text(encoding="utf-8"))
    return plan, connectivity_plan


def numeric_pair_evidence(
    baseline: pd.DataFrame,
    candidate: pd.DataFrame,
    *,
    key: list[str],
    fields: list[str],
) -> tuple[pd.DataFrame, dict[str, object]]:
    paired = baseline.merge(
        candidate,
        on=key,
        how="outer",
        suffixes=("Baseline", "Candidate"),
        indicator=True,
        validate="one_to_one",
    )
    evidence: dict[str, object] = {
        "baseline_rows": len(baseline),
        "candidate_rows": len(candidate),
        "paired_rows": int(paired["_merge"].eq("both").sum()),
        "baseline_only_rows": int(paired["_merge"].eq("left_only").sum()),
        "candidate_only_rows": int(paired["_merge"].eq("right_only").sum()),
    }
    for field in fields:
        left = pd.to_numeric(paired[f"{field}Baseline"], errors="coerce")
        right = pd.to_numeric(paired[f"{field}Candidate"], errors="coerce")
        finite = np.isfinite(left) & np.isfinite(right)
        evidence[f"{field}_finite_pairs"] = int(finite.sum())
        evidence[f"{field}_nan_pattern_mismatches"] = int(
            (left.isna() != right.isna()).sum()
        )
        evidence[f"{field}_max_abs_difference"] = (
            float(np.max(np.abs(left[finite] - right[finite]))) if finite.any() else 0.0
        )
    return paired, evidence


def boolean_mismatches(
    paired: pd.DataFrame, fields: list[str]
) -> dict[str, int]:
    return {
        field: int(
            (
                paired[f"{field}Baseline"].astype("boolean")
                != paired[f"{field}Candidate"].astype("boolean")
            )
            .fillna(True)
            .sum()
        )
        for field in fields
    }


def report_text(decision: dict[str, object]) -> str:
    return f"""# CMLE finite-domain optimizer remediation

> Prospectively registered same-byte numerical remediation. Returnability is not convergence, finite-MLE existence, inference readiness, estimator performance, or public UI readiness.

## Contract

- Contract passed: `{decision['contract_passed']}`
- Same retained response SHA-256: `{decision['same_retained_response_bytes']}`
- Audits / prefit eligible attempts: `{decision['candidate_audits']}/350` and `{decision['candidate_fit_attempts']}/287`
- Candidate returned / inference-ready: `{decision['candidate_returned_fits']}/287` and `{decision['candidate_inference_ready_fits']}/287`
- Baseline returned / inference-ready: `286/287` and `274/287`
- Invalid optimizer trials: `{decision['optimizer_invalid_evaluations_total']}` across `{decision['fits_with_invalid_optimizer_trials']}` fits
- Baseline-returned readiness mismatches: `{decision['baseline_readiness_mismatches']}`
- Shared Person raw/rounded flag mismatches: `{decision['person_flag_mismatches_total']}`
- Selected regression tests passed: `{decision['selected_tests_passed']}`

## Target failure

`chain_minimal::rep-00006 / a0_unanchored` now returned `{decision['target_returned']}` and remained inference-ready `{decision['target_inference_ready']}`. It recorded `{decision['target_invalid_optimizer_evaluations']}` invalid trial evaluations, final rank/nullity `{decision['target_information_rank']}/{decision['target_information_nullity']}`, and gradient sup norm `{decision['target_gradient_sup_norm']:.6g}`.

## Numeric preservation

- Structural estimate maximum absolute difference: `{decision['structural_estimate_max_abs_difference']:.3g}`
- Structural SE maximum absolute difference: `{decision['structural_se_max_abs_difference']:.3g}`
- Conditional log-likelihood maximum absolute difference: `{decision['conditional_loglik_max_abs_difference']:.3g}`
- WLE estimate / SE maximum absolute differences: `{decision['wle_estimate_max_abs_difference']:.3g}` / `{decision['wle_se_max_abs_difference']:.3g}`
- Infit / Outfit maximum absolute differences: `{decision['infit_max_abs_difference']:.3g}` / `{decision['outfit_max_abs_difference']:.3g}`

## Boundary

The optimizer-only adapter rejects non-finite trial points and retains the best finite point; the public exact objective still fails closed on non-finite direct evaluations. No stationarity, rank, covariance, anchor, or fit-decision gate was relaxed. The target changed from an exception to a typed returned-but-nonready result, not to a successful estimate. Public UI, separation detection, and performance claims remain withheld.
"""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--plan", type=Path, default=DEFAULT_PLAN)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    plan_path = args.plan.resolve()
    output = args.output.resolve()
    plan, connectivity_plan = load_and_validate(plan_path)
    if output.exists() and any(output.iterdir()):
        if not args.overwrite:
            raise FileExistsError(f"Output is non-empty: {output}; pass --overwrite explicitly.")
        shutil.rmtree(output)
    output.mkdir(parents=True, exist_ok=True)
    shutil.copy2(plan_path, output / plan_path.name)

    response_path = BASELINE_DIR / "connectivity_responses.csv"
    responses = pd.read_csv(response_path, float_precision="round_trip")
    baseline_fit = pd.read_csv(BASELINE_DIR / "connectivity_fit_ledger.csv")
    baseline_structural = pd.read_csv(
        BASELINE_DIR / "connectivity_structural_results.csv"
    )
    baseline_person = pd.read_csv(BASELINE_DIR / "connectivity_person_results.csv")
    topologies = {
        str(row["topology"]): row for row in connectivity_plan["topologies"]
    }

    audits: list[dict[str, object]] = []
    fits: list[dict[str, object]] = []
    structural_frames: list[pd.DataFrame] = []
    person_frames: list[pd.DataFrame] = []
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

    candidate_audit = pd.DataFrame(audits).sort_values(
        ["TopologyOrder", "Replicate", "ScenarioOrder"]
    ).reset_index(drop=True)
    candidate_fit = pd.DataFrame(fits).sort_values(
        ["TopologyOrder", "Replicate", "ScenarioOrder"]
    ).reset_index(drop=True)
    candidate_structural = pd.concat(structural_frames, ignore_index=True)
    candidate_person = pd.concat(person_frames, ignore_index=True)
    candidate_structural["StepKey"] = (
        pd.to_numeric(candidate_structural["Step"], errors="coerce").fillna(-1).astype(int)
    )
    baseline_structural["StepKey"] = (
        pd.to_numeric(baseline_structural["Step"], errors="coerce").fillna(-1).astype(int)
    )

    baseline_returned = baseline_fit.loc[baseline_fit["Returned"]].copy()
    candidate_for_baseline = candidate_fit.merge(
        baseline_returned[FIT_KEY], on=FIT_KEY, how="inner", validate="one_to_one"
    )
    fit_pairs, fit_evidence = numeric_pair_evidence(
        baseline_returned,
        candidate_for_baseline,
        key=FIT_KEY,
        fields=["ConditionalLogLik", "GradientSupNorm", "InformationRank", "InformationNullity"],
    )
    readiness_mismatches = boolean_mismatches(
        fit_pairs, ["Returned", "Converged", "InferenceReady", "AnchorExact"]
    )
    baseline_structural_compare = baseline_structural.drop(columns=["Step"], errors="ignore")
    candidate_structural_compare = (
        candidate_structural.merge(
            baseline_returned[FIT_KEY],
            on=FIT_KEY,
            how="inner",
            validate="many_to_one",
        )
        .drop(columns=["Step"], errors="ignore")
    )
    structural_pairs, structural_evidence = numeric_pair_evidence(
        baseline_structural_compare,
        candidate_structural_compare,
        key=STRUCTURAL_KEY,
        fields=["Estimate", "SE"],
    )
    person_pairs, person_evidence = numeric_pair_evidence(
        baseline_person,
        candidate_person,
        key=PERSON_KEY,
        fields=["WLEEstimate", "WLEStandardError", "Infit", "Outfit"],
    )
    person_flag_fields = [
        "FitEligibleRaw",
        "EitherUpperRaw",
        "EitherUpperRounded3",
        "RawRounded3Disagreement",
        "EitherUpperRounded6",
        "RawRounded6Disagreement",
    ]
    person_flag_mismatches = boolean_mismatches(person_pairs, person_flag_fields)

    target_contract = plan["baseline_contract"]["target_failure"]
    target = candidate_fit.loc[
        candidate_fit["RunId"].eq(target_contract["run_id"])
        & candidate_fit["Scenario"].eq(target_contract["scenario"])
    ].iloc[0]
    test_files = [
        "tests/test_cmle.py",
        "tests/test_cmle_hard_anchors.py",
        "tests/test_cmle_person_scoring.py",
        "tests/test_cmle_wle_bootstrap.py",
        "tests/test_cmle_wle_fit.py",
        "tests/test_cmle_wle_uncertainty.py",
        "tests/test_decision_stability.py",
        "tests/test_person_scoring_wle.py",
    ]
    completed = subprocess.run(
        [sys.executable, "-m", "pytest", "-q", *test_files],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    test_log = completed.stdout + completed.stderr
    (output / "selected_pytest.txt").write_text(test_log, encoding="utf-8")
    selected_tests_passed = completed.returncode == 0

    tolerance = plan["full_same_byte_replay"]["primary_numeric_tolerance"]
    numeric_passed = bool(
        structural_evidence["Estimate_max_abs_difference"]
        <= float(tolerance["structural_estimate_max_abs"])
        and structural_evidence["SE_max_abs_difference"]
        <= float(tolerance["structural_se_max_abs"])
        and fit_evidence["ConditionalLogLik_max_abs_difference"]
        <= float(tolerance["conditional_loglik_max_abs"])
        and person_evidence["WLEEstimate_max_abs_difference"]
        <= float(tolerance["wle_estimate_max_abs"])
        and person_evidence["WLEStandardError_max_abs_difference"]
        <= float(tolerance["wle_se_max_abs"])
        and person_evidence["Infit_max_abs_difference"]
        <= float(tolerance["infit_outfit_max_abs"])
        and person_evidence["Outfit_max_abs_difference"]
        <= float(tolerance["infit_outfit_max_abs"])
        and all(
            structural_evidence[f"{field}_nan_pattern_mismatches"] == 0
            for field in ("Estimate", "SE")
        )
        and all(
            person_evidence[f"{field}_nan_pattern_mismatches"] == 0
            for field in ("WLEEstimate", "WLEStandardError", "Infit", "Outfit")
        )
    )
    readiness_identical = all(value == 0 for value in readiness_mismatches.values())
    flags_identical = all(value == 0 for value in person_flag_mismatches.values())
    target_passed = bool(
        target["Returned"]
        and not target["InferenceReady"]
        and int(target["OptimizerInvalidEvaluations"]) > 0
    )
    ledger_passed = bool(
        len(candidate_audit) == 350
        and int(candidate_audit["ExactEligible"].sum()) == 287
        and len(candidate_fit) == 287
        and int(candidate_fit["Returned"].sum()) == 287
        and int(candidate_fit["InferenceReady"].sum()) == 274
    )
    same_bytes = sha256_file(response_path) == plan["baseline_identity"][
        "validation/cmle_native_anchor_connectivity_20260810/connectivity_responses.csv"
    ]
    contract_passed = bool(
        ledger_passed
        and same_bytes
        and numeric_passed
        and readiness_identical
        and flags_identical
        and target_passed
        and selected_tests_passed
    )

    output_frames = {
        "candidate_audit_ledger.csv": candidate_audit,
        "candidate_fit_ledger.csv": candidate_fit,
        "candidate_structural_results.csv": candidate_structural.drop(columns="StepKey"),
        "candidate_person_results.csv": candidate_person,
        "fit_pair_comparison.csv": fit_pairs,
        "structural_pair_comparison.csv": structural_pairs,
        "person_pair_comparison.csv": person_pairs,
    }
    for name, frame in output_frames.items():
        write_csv(frame, output / name)
    comparison = {
        "fit": fit_evidence,
        "fit_boolean_mismatches": readiness_mismatches,
        "structural": structural_evidence,
        "person": person_evidence,
        "person_flag_mismatches": person_flag_mismatches,
    }
    (output / "numeric_comparison.json").write_text(
        json.dumps(comparison, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    decision = {
        "schema_version": "mfrm-cmle-finite-domain-optimizer-remediation-decision-v1",
        "study_id": plan["study_id"],
        "contract_passed": contract_passed,
        "plan_sha256": sha256_file(plan_path),
        "candidate_core_sha256": sha256_file(ROOT / "mfrm_app/cmle.py"),
        "candidate_test_sha256": sha256_file(ROOT / "tests/test_cmle.py"),
        "candidate_connectivity_runner_sha256": sha256_file(
            ROOT / "validation/cmle_native_anchor_connectivity_stress.py"
        ),
        "same_retained_response_bytes": same_bytes,
        "candidate_audits": len(candidate_audit),
        "candidate_prefit_eligible": int(candidate_audit["ExactEligible"].sum()),
        "candidate_fit_attempts": len(candidate_fit),
        "candidate_returned_fits": int(candidate_fit["Returned"].sum()),
        "candidate_inference_ready_fits": int(candidate_fit["InferenceReady"].sum()),
        "optimizer_invalid_evaluations_total": int(
            candidate_fit["OptimizerInvalidEvaluations"].sum()
        ),
        "fits_with_invalid_optimizer_trials": int(
            candidate_fit["OptimizerInvalidEvaluations"].gt(0).sum()
        ),
        "best_finite_fallbacks": int(
            candidate_fit["OptimizerBestFiniteFallbackUsed"].sum()
        ),
        "baseline_readiness_mismatches": int(sum(readiness_mismatches.values())),
        "person_flag_mismatches_total": int(sum(person_flag_mismatches.values())),
        "selected_tests_passed": selected_tests_passed,
        "target_returned": bool(target["Returned"]),
        "target_inference_ready": bool(target["InferenceReady"]),
        "target_invalid_optimizer_evaluations": int(target["OptimizerInvalidEvaluations"]),
        "target_information_rank": int(target["InformationRank"]),
        "target_information_nullity": int(target["InformationNullity"]),
        "target_gradient_sup_norm": float(target["GradientSupNorm"]),
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
        "infit_max_abs_difference": float(person_evidence["Infit_max_abs_difference"]),
        "outfit_max_abs_difference": float(person_evidence["Outfit_max_abs_difference"]),
        "performance_value_success_gate": None,
        "public_ui_ready": False,
        "output_sha256": {},
    }
    hash_names = [*output_frames, "numeric_comparison.json", "selected_pytest.txt"]
    decision["output_sha256"] = {
        name: sha256_file(output / name) for name in hash_names
    }
    decision_path = output / "remediation_decision.json"
    decision_path.write_text(
        json.dumps(decision, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (output / "CMLE_FINITE_DOMAIN_OPTIMIZER_REMEDIATION.md").write_text(
        report_text(decision), encoding="utf-8"
    )
    print(json.dumps(decision, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
