#!/usr/bin/env python3
"""Run the registered PCM hard-anchor cross-engine study."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from mfrm_app.cmle import (  # noqa: E402
    cmle_objective_value_grad,
    prepare_cmle_design,
)
from mfrm_app.cmle_workflow import run_cmle_calibration_workflow  # noqa: E402
from validation.cmle_cross_engine_boundary import (  # noqa: E402
    install_mfrmr,
    json_safe,
    mfrmr_source_identity_without_git,
    sha256_file,
    write_csv,
)
from validation.cmle_finite_mle_existence_stress import _pcm_frame  # noqa: E402


PLAN = ROOT / "validation/cmle_pcm_anchor_cross_engine_plan_20260810.json"
REMEDIATION_PLAN = (
    ROOT
    / "validation/cmle_pcm_anchor_cross_engine_identity_remediation_plan_20260810.json"
)
REMEDIATION_AMENDMENT = (
    ROOT
    / "validation/cmle_pcm_anchor_cross_engine_identity_remediation_amendment_20260810.json"
)
OUTPUT = ROOT / "validation/cmle_pcm_anchor_cross_engine_20260810"
R_ADAPTER = ROOT / "validation/cmle_pcm_anchor_cross_engine.R"
SIRT_WORKER = ROOT / "validation/cmle_cross_engine_sirt_worker.R"
EXPECTED_WORKFLOW = {
    "interior_finite_cmle_supported": "calibration_ready",
    "boundary_no_finite_cmle": "finite_mle_boundary",
}


def validate_plan(path: Path) -> dict[str, object]:
    plan = json.loads(path.read_text(encoding="utf-8"))
    if plan.get("study_id") != "cmle_pcm_hard_anchor_cross_engine_v1":
        raise ValueError("Unexpected PCM-anchor plan.")
    mismatches = [
        relative
        for relative, digest in plan["parent_identity"].items()
        if not (ROOT / relative).exists()
        or sha256_file(ROOT / relative) != digest
    ]
    if mismatches:
        raise ValueError(f"PCM-anchor parent identity failed: {mismatches}")
    remediation = json.loads(REMEDIATION_PLAN.read_text(encoding="utf-8"))
    if (
        remediation.get("study_id")
        != "cmle_pcm_hard_anchor_cross_engine_identity_remediation_v1"
        or remediation["parent_plan"]["sha256"] != sha256_file(path)
    ):
        raise ValueError("PCM-anchor identity-remediation plan is invalid.")
    amendment = json.loads(REMEDIATION_AMENDMENT.read_text(encoding="utf-8"))
    correction = amendment["bookkeeping_correction"]
    if (
        amendment.get("study_id")
        != "cmle_pcm_hard_anchor_cross_engine_identity_remediation_amendment_v1"
        or amendment["parent_remediation_plan"]["sha256"]
        != sha256_file(REMEDIATION_PLAN)
        or correction["incorrect_per_function_sha256_retained_in_parent"]
        != remediation["stable_function_identity_contract"][
            "mfrmr_fit_mfrm_source_stripped_sha256"
        ]
    ):
        raise ValueError("PCM-anchor identity-remediation amendment is invalid.")
    snapshot = Path(plan["r_identity"]["mfrmr_frozen_snapshot_path"])
    if (
        not snapshot.is_dir()
        or mfrmr_source_identity_without_git(snapshot)
        != plan["r_identity"]["mfrmr_frozen_snapshot_relative_sha256"]
    ):
        raise ValueError("Frozen mfrmr snapshot is missing or changed.")
    plan["_identity_remediation"] = remediation
    plan["_identity_remediation_amendment"] = amendment
    return plan


def hard_anchor_frame(rows: list[dict[str, object]]) -> pd.DataFrame | None:
    if not rows:
        return None
    return pd.DataFrame(
        [
            {
                "ParameterType": "Facet",
                "Facet": row["Facet"],
                "Level": row["Level"],
                "Value": float(row["Value"]),
            }
            for row in rows
        ]
    )


def registered_cases(plan: dict[str, object]) -> list[dict[str, object]]:
    cases: list[dict[str, object]] = []
    for item in plan["fixture_contract"]["cases"]:
        anchors = hard_anchor_frame(item["anchors"])
        cases.append(
            {
                "CaseId": item["case_id"],
                "Mechanism": item["mechanism"],
                "ExpectedStatus": item["expected_status"],
                "AnchorLabel": item["case_id"].removeprefix(
                    f"polytomous_pcm_{item['mechanism']}_"
                ),
                "frame": _pcm_frame(item["mechanism"]),
                "prepare_args": {
                    "person_col": "Person",
                    "facet_cols": ["Rater", "Criterion"],
                    "score_col": "Score",
                    "rating_min": 0,
                    "rating_max": 2,
                    "model": "PCM",
                    "step_facet": "Criterion",
                    "hard_anchors": anchors,
                },
            }
        )
    expected = plan["fixture_contract"]["expected_accounting"]
    counts = pd.Series([case["ExpectedStatus"] for case in cases]).value_counts()
    if len(cases) != expected["total"] or any(
        int(counts.get(status, 0)) != count
        for status, count in expected.items()
        if status != "total"
    ):
        raise ValueError("PCM-anchor fixture accounting changed.")
    return cases


def build_bundle(
    cases: list[dict[str, object]], output: Path
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    manifests: list[dict[str, object]] = []
    ratings: list[pd.DataFrame] = []
    anchor_rows: list[dict[str, object]] = []
    runs: list[dict[str, object]] = []
    prefit_rows: list[dict[str, object]] = []
    coefficients: list[pd.DataFrame] = []
    for case in cases:
        anchors = case["prepare_args"]["hard_anchors"]
        anchor_count = 0 if anchors is None else len(anchors)
        manifests.append(
            {
                "CaseId": case["CaseId"],
                "Family": "polytomous_pcm",
                "Mechanism": case["Mechanism"],
                "ExpectedStatus": case["ExpectedStatus"],
                "ExpectedWorkflowStatus": EXPECTED_WORKFLOW[case["ExpectedStatus"]],
                "AnchorLabel": case["AnchorLabel"],
                "AnchorCount": anchor_count,
                "RatingMax": 2,
                "Rows": len(case["frame"]),
                "Persons": case["frame"]["Person"].nunique(),
            }
        )
        frame = case["frame"].copy()
        frame.insert(0, "CaseId", case["CaseId"])
        frame["Unit"] = frame["Criterion"]
        ratings.append(frame[["CaseId", "Person", "Rater", "Criterion", "Unit", "Score"]])
        if anchors is not None:
            for row in anchors.itertuples(index=False):
                anchor_rows.append(
                    {
                        "CaseId": case["CaseId"],
                        "Facet": row.Facet,
                        "Level": row.Level,
                        "Value": float(row.Value),
                    }
                )
        design = prepare_cmle_design(case["frame"], **case["prepare_args"])
        zero = np.zeros(design.n_parameters, dtype=float)
        nll_zero, gradient_zero = cmle_objective_value_grad(zero, design)
        for parameter, gradient in zip(
            design.parameter_names, gradient_zero, strict=True
        ):
            prefit_rows.append(
                {
                    "CaseId": case["CaseId"],
                    "Parameter": parameter,
                    "PythonGradientZero": float(gradient),
                    "PythonNLLZero": float(nll_zero),
                }
            )
        started = time.perf_counter()
        result = run_cmle_calibration_workflow(
            case["frame"], **case["prepare_args"], gtol=1e-8, maxiter=800
        )
        elapsed = time.perf_counter() - started
        summary = result["summary"].iloc[0]
        fit = result["fit"]
        fit_summary = fit["summary"].iloc[0] if fit is not None else None
        anchor_exact = True
        if fit is not None and anchors is not None:
            facets = fit["facets"]["others"]
            merged = anchors.merge(
                facets[["Facet", "Level", "Estimate", "SE"]],
                on=["Facet", "Level"],
                how="left",
                validate="one_to_one",
            )
            anchor_exact = bool(
                np.array_equal(
                    merged["Estimate"].to_numpy(float),
                    merged["Value"].to_numpy(float),
                )
                and np.array_equal(merged["SE"].to_numpy(float), np.zeros(len(merged)))
            )
        runs.append(
            {
                "CaseId": case["CaseId"],
                "ExpectedStatus": case["ExpectedStatus"],
                "WorkflowStatus": summary["WorkflowStatus"],
                "TerminalMatch": summary["WorkflowStatus"]
                == EXPECTED_WORKFLOW[case["ExpectedStatus"]],
                "AnchorCount": anchor_count,
                "ParameterNames": "|".join(design.parameter_names),
                "KParams": design.n_parameters,
                "PrefitRank": int(design.audit["conditional_rank"]),
                "PrefitNullity": int(design.audit["conditional_nullity"]),
                "PythonNLLZero": float(nll_zero),
                "OptimizationAttempted": bool(summary["OptimizationAttempted"]),
                "FitReturned": bool(summary["FitReturned"]),
                "InferenceReady": bool(summary["FitInferenceReady"]),
                "AnchorExact": anchor_exact,
                "ConditionalLogLik": (
                    float(fit_summary["ConditionalLogLik"])
                    if fit_summary is not None
                    else np.nan
                ),
                "ElapsedSeconds": elapsed,
            }
        )
        if fit is not None:
            table = fit["coefficients"].copy()
            table.insert(0, "CaseId", case["CaseId"])
            coefficients.append(table)
    manifest = pd.DataFrame(manifests)
    rating_table = pd.concat(ratings, ignore_index=True)
    anchors_table = pd.DataFrame(
        anchor_rows, columns=["CaseId", "Facet", "Level", "Value"]
    )
    python_runs = pd.DataFrame(runs)
    python_prefit = pd.DataFrame(prefit_rows)
    python_coefficients = pd.concat(coefficients, ignore_index=True)
    tables = {
        "manifest.csv": manifest,
        "ratings.csv": rating_table,
        "anchors.csv": anchors_table,
        "python_runs.csv": python_runs,
        "python_prefit.csv": python_prefit,
        "python_coefficients.csv": python_coefficients,
    }
    for name, table in tables.items():
        write_csv(table, output / name)
    inventory = pd.DataFrame(
        [{"File": name, "SHA256": sha256_file(output / name)} for name in tables]
    )
    write_csv(inventory, output / "bundle_inventory.csv")
    return manifest, python_runs, python_prefit, python_coefficients


def r_environment(library: Path, current_paths: list[str]) -> dict[str, str]:
    return {
        **os.environ,
        "R_LIBS": os.pathsep.join([str(library), *current_paths]),
    }


def run_r_adapter(
    output: Path,
    library: Path,
    current_paths: list[str],
    plan: dict[str, object],
) -> None:
    result = subprocess.run(
        [
            "Rscript",
            str(R_ADAPTER),
            "--input",
            str(output),
            "--output",
            str(output),
            "--plan-sha256",
            sha256_file(PLAN),
            "--mfrmr-source-sha256",
            plan["r_identity"]["mfrmr_frozen_snapshot_relative_sha256"],
        ],
        env=r_environment(library, current_paths),
        capture_output=True,
        text=True,
    )
    (output / "r_adapter_stdout.txt").write_text(result.stdout, encoding="utf-8")
    (output / "r_adapter_stderr.txt").write_text(result.stderr, encoding="utf-8")
    if result.returncode:
        raise RuntimeError("PCM-anchor R adapter failed; inspect retained logs.")


def run_sirt_lane(
    output: Path, library: Path, current_paths: list[str]
) -> None:
    manifest = pd.read_csv(output / "manifest.csv")
    python_runs = pd.read_csv(output / "python_runs.csv")
    worker_dir = output / "sirt_workers"
    worker_dir.mkdir()
    rows: list[pd.DataFrame] = []
    for item in manifest.itertuples(index=False):
        if int(item.AnchorCount) > 0:
            rows.append(
                pd.DataFrame(
                    [
                        {
                            "CaseId": item.CaseId,
                            "Mechanism": item.Mechanism,
                            "ExpectedStatus": item.ExpectedStatus,
                            "AnchorLabel": item.AnchorLabel,
                            "AnchorCount": item.AnchorCount,
                            "PythonWorkflowStatus": python_runs.loc[
                                python_runs["CaseId"].eq(item.CaseId),
                                "WorkflowStatus",
                            ].iloc[0],
                            "Engine": "sirt",
                            "Estimator": "MML",
                            "Mode": "rm.facets_PCM_isolated",
                            "FitAttempted": False,
                            "FitReturned": False,
                            "Converged": False,
                            "MaxAbsEstimate": np.nan,
                            "MaxAbsSE": np.nan,
                            "Warnings": "",
                            "Error": "anchor_scope_unsupported: equivalent PCM fixed-offset mapping not registered",
                        }
                    ]
                )
            )
            continue
        path = worker_dir / f"{item.CaseId}.csv"
        result = subprocess.run(
            [
                "Rscript",
                str(SIRT_WORKER),
                "--input",
                str(output),
                "--output",
                str(path),
                "--case-id",
                item.CaseId,
            ],
            env=r_environment(library, current_paths),
            capture_output=True,
            text=True,
        )
        (worker_dir / f"{item.CaseId}.stdout.txt").write_text(
            result.stdout, encoding="utf-8"
        )
        (worker_dir / f"{item.CaseId}.stderr.txt").write_text(
            result.stderr, encoding="utf-8"
        )
        if result.returncode or not path.exists():
            rows.append(
                pd.DataFrame(
                    [
                        {
                            "CaseId": item.CaseId,
                            "Mechanism": item.Mechanism,
                            "ExpectedStatus": item.ExpectedStatus,
                            "AnchorLabel": item.AnchorLabel,
                            "AnchorCount": 0,
                            "Engine": "sirt",
                            "Estimator": "MML",
                            "Mode": "rm.facets_PCM_isolated",
                            "FitAttempted": True,
                            "FitReturned": False,
                            "Converged": False,
                            "Error": f"isolated_process_exit_{result.returncode}",
                        }
                    ]
                )
            )
        else:
            worker = pd.read_csv(path)
            worker["AnchorLabel"] = item.AnchorLabel
            worker["AnchorCount"] = 0
            rows.append(worker)
    sirt = pd.concat(rows, ignore_index=True, sort=False)
    secondary = pd.read_csv(output / "secondary_engine_runs.csv")
    write_csv(
        pd.concat([secondary, sirt], ignore_index=True, sort=False),
        output / "secondary_engine_runs.csv",
    )


def analyze(
    plan: dict[str, object],
    manifest: pd.DataFrame,
    python_runs: pd.DataFrame,
    python_prefit: pd.DataFrame,
    python_coefficients: pd.DataFrame,
    output: Path,
) -> dict[str, object]:
    immer = pd.read_csv(output / "immer_runs.csv")
    r_prefit = pd.read_csv(output / "immer_prefit_gradient.csv")
    r_coefficients = pd.read_csv(output / "immer_coefficients.csv")
    secondary = pd.read_csv(output / "secondary_engine_runs.csv")
    identity = pd.read_csv(output / "r_engine_identity.csv").iloc[0]
    prefit = python_prefit.merge(
        r_prefit,
        on=["CaseId", "Parameter"],
        how="outer",
        indicator=True,
    )
    prefit["GradientAbsoluteDifference"] = (
        prefit["PythonGradientZero"] - prefit["RGradientZero"]
    ).abs()
    prefit["NLLAbsoluteDifference"] = (
        prefit["PythonNLLZero"] - prefit["RNLLZero"]
    ).abs()
    write_csv(prefit, output / "immer_prefit_comparison.csv")
    primary = immer.loc[immer["MaxIt"].eq(2000)].copy()
    rank = python_runs[
        ["CaseId", "ParameterNames", "KParams", "PrefitRank", "PrefitNullity"]
    ].merge(
        primary[
            [
                "CaseId",
                "ParameterNames",
                "KParams",
                "RPrefitRank",
                "RPrefitNullity",
            ]
        ],
        on="CaseId",
        suffixes=("Python", "Immer"),
        validate="one_to_one",
    )
    rank["ParameterNamesMatch"] = rank["ParameterNamesPython"].eq(
        rank["ParameterNamesImmer"]
    )
    rank["KMatch"] = rank["KParamsPython"].eq(rank["KParamsImmer"])
    rank["RankMatch"] = rank["PrefitRank"].eq(rank["RPrefitRank"])
    rank["NullityMatch"] = rank["PrefitNullity"].eq(rank["RPrefitNullity"])
    write_csv(rank, output / "immer_prefit_rank_and_names.csv")
    interior_ids = set(
        manifest.loc[
            manifest["ExpectedStatus"].eq("interior_finite_cmle_supported"),
            "CaseId",
        ]
    )
    boundary_ids = set(manifest["CaseId"]) - interior_ids
    py_coeff = python_coefficients.loc[
        python_coefficients["CaseId"].isin(interior_ids)
    ]
    r_coeff = r_coefficients.loc[
        r_coefficients["CaseId"].isin(interior_ids)
        & r_coefficients["MaxIt"].eq(2000)
    ]
    coefficient_comparison = py_coeff.merge(
        r_coeff,
        on=["CaseId", "Parameter"],
        how="outer",
        suffixes=("Python", "Immer"),
        indicator=True,
    )
    coefficient_comparison["EstimateAbsoluteDifference"] = (
        coefficient_comparison["EstimatePython"]
        - coefficient_comparison["EstimateImmer"]
    ).abs()
    coefficient_comparison["SEAbsoluteDifference"] = (
        coefficient_comparison["SEPython"] - coefficient_comparison["SEImmer"]
    ).abs()
    write_csv(
        coefficient_comparison, output / "immer_interior_coefficient_comparison.csv"
    )
    loglik = python_runs.loc[
        python_runs["CaseId"].isin(interior_ids),
        ["CaseId", "ConditionalLogLik"],
    ].merge(
        primary.loc[
            primary["CaseId"].isin(interior_ids),
            ["CaseId", "ConditionalLogLik"],
        ],
        on="CaseId",
        suffixes=("Python", "Immer"),
        validate="one_to_one",
    )
    loglik["AbsoluteDifference"] = (
        loglik["ConditionalLogLikPython"] - loglik["ConditionalLogLikImmer"]
    ).abs()
    write_csv(loglik, output / "immer_interior_loglik_comparison.csv")
    boundary = immer.loc[immer["CaseId"].isin(boundary_ids)].copy()
    boundary["FinalReadinessAfterOracle"] = False
    boundary["FinalDecision"] = "blocked_no_finite_cmle"
    write_csv(boundary, output / "immer_boundary_maxit_sensitivity.csv")
    secondary_summary = (
        secondary.assign(HasError=secondary["Error"].fillna("").ne(""))
        .groupby(["Engine", "ExpectedStatus"], as_index=False)
        .agg(
            Cases=("CaseId", "nunique"),
            FitAttempted=("FitAttempted", "sum"),
            FitReturned=("FitReturned", "sum"),
            Converged=("Converged", "sum"),
            Errors=("HasError", "sum"),
            MedianMaxAbsEstimate=("MaxAbsEstimate", "median"),
            MaximumMaxAbsEstimate=("MaxAbsEstimate", "max"),
        )
    )
    write_csv(secondary_summary, output / "secondary_engine_by_status.csv")
    interior_primary = primary.loc[primary["CaseId"].isin(interior_ids)]
    expected = plan["fixture_contract"]["expected_accounting"]
    stable_identity = plan["_identity_remediation"][
        "stable_function_identity_contract"
    ]
    corrected_mfrmr_identity = plan["_identity_remediation_amendment"][
        "bookkeeping_correction"
    ]["correct_final_named_vector_sha256"]
    results = {
        "identity_passed": bool(
            identity["PlanSHA256"] == sha256_file(PLAN)
            and str(identity["ImmerVersion"]) == plan["r_identity"]["immer_version"]
            and identity["ImmerFunctionSHA256"]
            == stable_identity["immer_source_stripped_sha256"]
            and str(identity["MfrmrVersion"]) == plan["r_identity"]["mfrmr_version"]
            and identity["MfrmrFunctionSHA256"]
            == corrected_mfrmr_identity
            and identity["TAMFunctionSHA256"]
            == stable_identity["TAM_source_stripped_sha256"]
            and identity["SirtFunctionSHA256"]
            == stable_identity["sirt_source_stripped_sha256"]
            and identity["MfrmrSourceSHA256"]
            == plan["r_identity"]["mfrmr_frozen_snapshot_relative_sha256"]
        ),
        "python_passed": bool(
            len(python_runs) == expected["total"]
            and python_runs["TerminalMatch"].all()
            and int(python_runs["OptimizationAttempted"].sum())
            == expected["interior_finite_cmle_supported"]
            and python_runs.loc[python_runs["InferenceReady"], "AnchorExact"].all()
        ),
        "immer_prefit_passed": bool(
            prefit["_merge"].eq("both").all()
            and rank[
                ["ParameterNamesMatch", "KMatch", "RankMatch", "NullityMatch"]
            ].all().all()
            and prefit["GradientAbsoluteDifference"].max() <= 1e-10
            and prefit["NLLAbsoluteDifference"].max() <= 1e-10
        ),
        "immer_interior_passed": bool(
            len(interior_primary) == expected["interior_finite_cmle_supported"]
            and interior_primary["FitReturned"].all()
            and interior_primary["ConvergenceCode"].eq(0).all()
            and interior_primary["AllCoefficientsFinite"].all()
            and interior_primary["AllSEFinite"].all()
            and interior_primary["InformationNullity"].eq(0).all()
            and interior_primary["GradientSupNorm"].le(1e-6).all()
        ),
        "immer_numeric_passed": bool(
            coefficient_comparison["_merge"].eq("both").all()
            and coefficient_comparison["EstimateAbsoluteDifference"].max() <= 1e-7
            and loglik["AbsoluteDifference"].max() <= 1e-8
        ),
        "boundary_complete": bool(
            len(boundary) == expected["boundary_no_finite_cmle"] * 3
            and not boundary["FinalReadinessAfterOracle"].any()
        ),
        "secondary_complete": bool(
            len(secondary) == expected["total"] * 3
            and secondary.groupby("Engine")["CaseId"]
            .nunique()
            .eq(expected["total"])
            .all()
        ),
        "prefit_zero_nll_max_absolute_difference": float(
            prefit["NLLAbsoluteDifference"].max()
        ),
        "prefit_zero_gradient_max_absolute_difference": float(
            prefit["GradientAbsoluteDifference"].max()
        ),
        "estimate_max_absolute_difference": float(
            coefficient_comparison["EstimateAbsoluteDifference"].max()
        ),
        "se_max_absolute_difference_descriptive": float(
            coefficient_comparison["SEAbsoluteDifference"].max()
        ),
        "conditional_loglik_max_absolute_difference": float(
            loglik["AbsoluteDifference"].max()
        ),
        "secondary_by_status": secondary_summary.to_dict(orient="records"),
    }
    return results


def run_tests(output: Path) -> dict[str, object]:
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "tests/test_cmle_pcm_anchor_cross_engine_contract.py",
            "tests/test_cmle.py",
            "tests/test_cmle_hard_anchors.py",
            "tests/test_cmle_existence.py",
            "tests/test_cmle_workflow.py",
        ],
        cwd=ROOT,
        env={**os.environ, "MPLCONFIGDIR": "/private/tmp/matplotlib-cmle-pcm-anchor"},
        capture_output=True,
        text=True,
    )
    (output / "selected_tests_stdout.txt").write_text(result.stdout, encoding="utf-8")
    (output / "selected_tests_stderr.txt").write_text(result.stderr, encoding="utf-8")
    return {"passed": result.returncode == 0, "returncode": result.returncode}


def critical_review(output: Path, decision: dict[str, object]) -> None:
    results = decision["results"]
    text = f"""# PCM Hard-Anchor Cross-Engine Critical Review

## Decision

The prospectively registered contract **{'passed' if decision['contract_passed'] else 'did not pass'}**. Direct numerical parity is limited to native Python exact PCM CMLE and a matched `immer::immer_cml` W/`b_const` parameterization.

## Matched CMLE evidence

Across six oracle-interior cases, including Rater, Criterion, contaminated, and mixed differential anchors, the maximum Python–immer free-coordinate difference was `{results['estimate_max_absolute_difference']:.17g}` and the conditional-log-likelihood difference was `{results['conditional_loglik_max_absolute_difference']:.17g}`. Before optimization, the maximum zero-point NLL and gradient differences were `{results['prefit_zero_nll_max_absolute_difference']:.17g}` and `{results['prefit_zero_gradient_max_absolute_difference']:.17g}`. This prefit check is essential: it tests the fixed-offset sign rather than trusting a final optimizer coincidence.

## Boundaries and anchors

All six registered separation/unused-support cases remained blocked by the exact support oracle at all three retained immer iteration caps. A hard anchor changes the free-coordinate cone and can change an existence result, but exact implementation does not validate anchor content. Criterion main-effect anchors also do not fix Criterion-specific PCM step transitions.

## Different estimators

The mfrmr JML, TAM MML, and sirt MML ledgers are descriptive sensitivity evidence only: {results['secondary_by_status']}. TAM and sirt anchor conditions are explicitly unsupported rather than silently approximated. Their values must not be compared as expected PCM-CMLE equality.

## Remaining risks

- The mfrmr 0.2.3 snapshot is dirty, unreleased, temporary, and identified only by content hash; a clean release rerun is mandatory.
- The study has two Raters, two Criteria, and category support 0..2. It is not a high-dimensional sparse-PCM stress.
- Step anchors, group anchors, slopes/GPCM, new-Person scoring, and calibration uncertainty are outside scope.
- Binary64 tolerances remain in optimization, finite differences, LP feasibility, and rank checks; raw values, not displayed rounding, drive decisions.
- Public UI and bilingual comprehension remain withheld.

## Product implication

The app can represent hard anchors in matched exact CMLE through a fixed-offset contract, but the one-click result must show anchored main effects, free effects, PCM steps, and finite-existence status separately. A successful anchored fit cannot be rendered as evidence that the anchor was unbiased.
"""
    (output / "CMLE_PCM_ANCHOR_CROSS_ENGINE_CRITICAL_REVIEW.md").write_text(
        text, encoding="utf-8"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--plan", type=Path, default=PLAN)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    plan = validate_plan(args.plan)
    if args.output.exists():
        raise FileExistsError(f"Evidence output exists: {args.output}")
    args.output.mkdir(parents=True)
    cases = registered_cases(plan)
    manifest, python_runs, python_prefit, python_coefficients = build_bundle(
        cases, args.output
    )
    snapshot = Path(plan["r_identity"]["mfrmr_frozen_snapshot_path"])
    source_identity = {
        "git_head": plan["r_identity"]["mfrmr_git_head_at_snapshot"],
        "git_status_sha256": "retained_in_parent_boundary_study",
        "build_content_sha256": plan["r_identity"][
            "mfrmr_frozen_snapshot_relative_sha256"
        ],
        "build_relative_content_sha256": plan["r_identity"][
            "mfrmr_frozen_snapshot_relative_sha256"
        ],
        "build_file_count": "546",
    }
    library, current_paths = install_mfrmr(snapshot, source_identity, args.output)
    run_r_adapter(args.output, library, current_paths, plan)
    run_sirt_lane(args.output, library, current_paths)
    results = analyze(
        plan,
        manifest,
        python_runs,
        python_prefit,
        python_coefficients,
        args.output,
    )
    tests = run_tests(args.output)
    contract_passed = bool(
        all(
            results[key]
            for key in (
                "identity_passed",
                "python_passed",
                "immer_prefit_passed",
                "immer_interior_passed",
                "immer_numeric_passed",
                "boundary_complete",
                "secondary_complete",
            )
        )
        and tests["passed"]
    )
    decision: dict[str, object] = {
        "study_id": plan["study_id"],
        "plan_sha256": sha256_file(args.plan),
        "identity_remediation_plan_sha256": sha256_file(REMEDIATION_PLAN),
        "identity_remediation_amendment_sha256": sha256_file(
            REMEDIATION_AMENDMENT
        ),
        "contract_passed": contract_passed,
        "implementation_sha256": {
            "validation/cmle_pcm_anchor_cross_engine.py": sha256_file(Path(__file__)),
            "validation/cmle_pcm_anchor_cross_engine.R": sha256_file(R_ADAPTER),
            "validation/cmle_cross_engine_sirt_worker.R": sha256_file(SIRT_WORKER),
            "tests/test_cmle_pcm_anchor_cross_engine_contract.py": sha256_file(
                ROOT / "tests/test_cmle_pcm_anchor_cross_engine_contract.py"
            ),
        },
        "results": results,
        "selected_tests": tests,
        "interpretation": {
            "direct_parity": "Python exact PCM CMLE versus matched immer W+b_const only",
            "secondary": "mfrmr JML and TAM/sirt MML are sensitivity lanes",
            "public_ui": "withheld",
        },
    }
    critical_review(args.output, decision)
    names = sorted(
        path.relative_to(args.output).as_posix()
        for path in args.output.rglob("*")
        if path.is_file() and path.name != "decision.json"
    )
    decision["output_sha256"] = {
        name: sha256_file(args.output / name) for name in names
    }
    decision = json_safe(decision)
    (args.output / "decision.json").write_text(
        json.dumps(decision, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(decision, indent=2, ensure_ascii=False, allow_nan=False))
    if not contract_passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
