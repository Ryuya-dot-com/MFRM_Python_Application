#!/usr/bin/env python3
"""Run the prospectively registered cross-engine CMLE boundary study."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from mfrm_app.cmle import prepare_cmle_design  # noqa: E402
from mfrm_app.cmle_workflow import run_cmle_calibration_workflow  # noqa: E402
from validation.cmle_finite_mle_existence_stress import fixture_cases  # noqa: E402


DEFAULT_PLAN = ROOT / "validation/cmle_cross_engine_boundary_plan_20260810.json"
DEFAULT_AMENDMENT = (
    ROOT / "validation/cmle_cross_engine_boundary_plan_amendment_20260810.json"
)
DEFAULT_AMENDMENT2 = (
    ROOT / "validation/cmle_cross_engine_boundary_plan_amendment2_20260810.json"
)
DEFAULT_OUTPUT = ROOT / "validation/cmle_cross_engine_boundary_20260810"
R_ADAPTER = ROOT / "validation/cmle_cross_engine_boundary.R"
SIRT_WORKER = ROOT / "validation/cmle_cross_engine_sirt_worker.R"
MFRMR_SOURCE = Path("/Users/tohokusla/Dropbox/MFRM_Application/mfrmr/development")
EXPECTED_TO_WORKFLOW = {
    "interior_finite_cmle_supported": "calibration_ready",
    "boundary_no_finite_cmle": "finite_mle_boundary",
    "structural_nonidentification": "design_not_identified",
}
BUILD_ROOTS = ("DESCRIPTION", "NAMESPACE", "R", "src", "inst", "data", "man")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def write_csv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, index=False, lineterminator="\n", float_format="%.17g")


def json_safe(value: object) -> object:
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def validate_plan(path: Path) -> dict[str, object]:
    plan = json.loads(path.read_text(encoding="utf-8"))
    if plan.get("study_id") != "cmle_cross_engine_boundary_semantics_v1":
        raise ValueError("Unexpected cross-engine boundary plan.")
    mismatches = [
        relative
        for relative, digest in plan["parent_identity"].items()
        if not (ROOT / relative).exists()
        or sha256_file(ROOT / relative) != digest
    ]
    if mismatches:
        raise ValueError(f"Cross-engine parent identity failed: {mismatches}")
    amendment = json.loads(DEFAULT_AMENDMENT.read_text(encoding="utf-8"))
    if (
        amendment.get("study_id")
        != "cmle_cross_engine_boundary_semantics_v1_identity_amendment_1"
        or amendment["parent_plan"]["sha256"] != sha256_file(path)
    ):
        raise ValueError("Cross-engine identity amendment is invalid or stale.")
    amendment2 = json.loads(DEFAULT_AMENDMENT2.read_text(encoding="utf-8"))
    if (
        amendment2.get("study_id")
        != "cmle_cross_engine_boundary_semantics_v1_identity_amendment_2"
        or amendment2["parent_amendment"]["sha256"]
        != sha256_file(DEFAULT_AMENDMENT)
    ):
        raise ValueError("Cross-engine frozen-snapshot amendment is invalid or stale.")
    plan["_effective_identity_amendment"] = amendment2
    return plan


def selected_fixture_cases(plan: dict[str, object]) -> list[dict[str, object]]:
    requested = list(plan["fixture_contract"]["case_ids"])
    available = {str(case["CaseId"]): case for case in fixture_cases()}
    missing = [case_id for case_id in requested if case_id not in available]
    if missing:
        raise ValueError(f"Registered fixture IDs are unavailable: {missing}")
    cases = [available[case_id] for case_id in requested]
    counts = pd.Series([case["ExpectedStatus"] for case in cases]).value_counts()
    expected = plan["fixture_contract"]["expected_accounting"]
    if len(cases) != int(expected["total"]) or any(
        int(counts.get(status, 0)) != int(count)
        for status, count in expected.items()
        if status != "total"
    ):
        raise ValueError("Registered fixture accounting changed.")
    return cases


def r_ready_frame(case: dict[str, object]) -> pd.DataFrame:
    frame = case["frame"].copy()
    unit = (
        frame["Event"].astype(str)
        if "Event" in frame.columns
        else pd.Series("Response", index=frame.index, dtype="object")
    )
    output = pd.DataFrame(
        {
            "CaseId": str(case["CaseId"]),
            "Person": frame["Person"].astype(str),
            "Rater": frame["Rater"].astype(str),
            "Unit": unit,
            "Score": pd.to_numeric(frame["Score"], errors="raise").astype(int),
        }
    )
    if output.duplicated(["CaseId", "Person", "Rater", "Unit"]).any():
        raise ValueError(f"Duplicate R virtual unit in {case['CaseId']}")
    return output


def mfrmr_source_identity(source: Path) -> dict[str, str]:
    head = subprocess.run(
        ["git", "-C", str(source), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
    ).stdout.decode().strip()
    status_bytes = subprocess.run(
        ["git", "-C", str(source), "status", "--porcelain=v1"],
        check=True,
        capture_output=True,
    ).stdout
    files: list[Path] = []
    for relative in BUILD_ROOTS:
        value = source / relative
        if value.is_file():
            files.append(value)
        elif value.is_dir():
            files.extend(path for path in value.rglob("*") if path.is_file())
    inventory = bytearray()
    relative_inventory = bytearray()
    for path in sorted(files, key=lambda value: os.fsencode(str(value))):
        digest = sha256_file(path)
        inventory.extend(f"{digest}  {path}\n".encode("utf-8"))
        relative_inventory.extend(
            f"{digest}  {path.relative_to(source)}\n".encode("utf-8")
        )
    return {
        "git_head": head,
        "git_status_sha256": sha256_bytes(status_bytes),
        "build_content_sha256": sha256_bytes(bytes(inventory)),
        "build_relative_content_sha256": sha256_bytes(
            bytes(relative_inventory)
        ),
        "build_file_count": str(len(files)),
    }


def assert_mfrmr_identity(
    observed: dict[str, str], plan: dict[str, object]
) -> None:
    registered = plan["_effective_identity_amendment"][
        "effective_mfrmr_source_identity"
    ]
    pairs = {
        "git_head": "git_head_at_snapshot",
        "git_status_sha256": "git_status_sha256_at_snapshot",
        "build_content_sha256": "snapshot_relative_content_sha256",
        "build_relative_content_sha256": "snapshot_relative_content_sha256",
    }
    mismatches = [
        label
        for label, registered_label in pairs.items()
        if observed[label] != registered[registered_label]
    ]
    if mismatches:
        raise ValueError(
            "mfrmr 0.2.3 development source changed after registration: "
            f"{mismatches}"
        )


def build_python_bundle(
    cases: list[dict[str, object]], output: Path
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    manifest_rows: list[dict[str, object]] = []
    rating_frames: list[pd.DataFrame] = []
    run_rows: list[dict[str, object]] = []
    coefficient_frames: list[pd.DataFrame] = []
    for case in cases:
        design = prepare_cmle_design(case["frame"], **case["prepare_args"])
        started = time.perf_counter()
        result = run_cmle_calibration_workflow(
            case["frame"],
            **case["prepare_args"],
            gtol=1e-8,
            maxiter=800,
        )
        elapsed = time.perf_counter() - started
        summary = result["summary"].iloc[0]
        expected_status = str(case["ExpectedStatus"])
        expected_workflow = EXPECTED_TO_WORKFLOW[expected_status]
        fit = result["fit"]
        fit_summary = fit["summary"].iloc[0] if fit is not None else None
        manifest_rows.append(
            {
                "CaseId": case["CaseId"],
                "Family": case["Family"],
                "Mechanism": case["Mechanism"],
                "ExpectedStatus": expected_status,
                "ExpectedWorkflowStatus": expected_workflow,
                "Model": case["prepare_args"]["model"],
                "RatingMin": case["prepare_args"]["rating_min"],
                "RatingMax": case["prepare_args"]["rating_max"],
                "Rows": len(case["frame"]),
                "Persons": case["frame"]["Person"].nunique(),
                "KParams": design.n_parameters,
            }
        )
        run_rows.append(
            {
                "CaseId": case["CaseId"],
                "ExpectedStatus": expected_status,
                "ExpectedWorkflowStatus": expected_workflow,
                "WorkflowStatus": summary["WorkflowStatus"],
                "TerminalMatch": summary["WorkflowStatus"] == expected_workflow,
                "PrefitEligible": bool(design.audit["eligible"]),
                "PrefitRank": int(design.audit["conditional_rank"]),
                "PrefitNullity": int(design.audit["conditional_nullity"]),
                "KParams": design.n_parameters,
                "OptimizationAttempted": bool(summary["OptimizationAttempted"]),
                "FitReturned": bool(summary["FitReturned"]),
                "InferenceReady": bool(summary["FitInferenceReady"]),
                "ConditionalLogLik": (
                    float(fit_summary["ConditionalLogLik"])
                    if fit_summary is not None
                    else np.nan
                ),
                "GradientSupNorm": (
                    float(fit_summary["GradientSupNorm"])
                    if fit_summary is not None
                    else np.nan
                ),
                "ElapsedSeconds": elapsed,
            }
        )
        if fit is not None:
            coefficients = fit["coefficients"].copy()
            coefficients.insert(0, "CaseId", case["CaseId"])
            coefficient_frames.append(coefficients)
        rating_frames.append(r_ready_frame(case))
    manifest = pd.DataFrame(manifest_rows)
    ratings = pd.concat(rating_frames, ignore_index=True)
    python_runs = pd.DataFrame(run_rows)
    python_coefficients = pd.concat(coefficient_frames, ignore_index=True)
    write_csv(manifest, output / "manifest.csv")
    write_csv(ratings, output / "ratings.csv")
    write_csv(python_runs, output / "python_runs.csv")
    write_csv(python_coefficients, output / "python_coefficients.csv")
    inventory = pd.DataFrame(
        [
            {"File": name, "SHA256": sha256_file(output / name)}
            for name in (
                "manifest.csv",
                "ratings.csv",
                "python_runs.csv",
                "python_coefficients.csv",
            )
        ]
    )
    write_csv(inventory, output / "bundle_inventory.csv")
    return manifest, ratings, python_runs, python_coefficients


def install_mfrmr(
    source: Path, source_identity: dict[str, str], output: Path
) -> tuple[Path, list[str]]:
    temp_root = Path(
        tempfile.mkdtemp(prefix="mfrmr023-cross-engine-", dir="/private/tmp")
    )
    library = temp_root / "library"
    library.mkdir()
    snapshot = temp_root / "mfrmr-0.2.3-snapshot"
    shutil.copytree(
        source,
        snapshot,
        ignore=shutil.ignore_patterns(".git", ".Rproj.user"),
    )
    snapshot_relative_hash = mfrmr_source_identity_without_git(snapshot)
    if (
        snapshot_relative_hash
        != source_identity["build_relative_content_sha256"]
    ):
        raise ValueError("mfrmr temporary snapshot content differs from source.")
    current_paths = subprocess.run(
        ["Rscript", "-e", 'cat(.libPaths(), sep="\\n")'],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.splitlines()
    result = subprocess.run(
        [
            "R",
            "CMD",
            "INSTALL",
            "--no-multiarch",
            "--with-keep.source",
            "-l",
            str(library),
            str(snapshot),
        ],
        capture_output=True,
        text=True,
    )
    (output / "mfrmr_install_stdout.txt").write_text(
        result.stdout, encoding="utf-8"
    )
    (output / "mfrmr_install_stderr.txt").write_text(
        result.stderr, encoding="utf-8"
    )
    if result.returncode != 0:
        raise RuntimeError(
            "Current mfrmr 0.2.3 snapshot installation failed; see retained logs."
        )
    return library, current_paths


def mfrmr_source_identity_without_git(source: Path) -> str:
    files: list[Path] = []
    for relative in BUILD_ROOTS:
        value = source / relative
        if value.is_file():
            files.append(value)
        elif value.is_dir():
            files.extend(path for path in value.rglob("*") if path.is_file())
    inventory = bytearray()
    for path in sorted(
        files, key=lambda value: os.fsencode(str(value.relative_to(source)))
    ):
        inventory.extend(
            f"{sha256_file(path)}  {path.relative_to(source)}\n".encode("utf-8")
        )
    return sha256_bytes(bytes(inventory))


def run_r_adapter(
    output: Path,
    plan_sha256: str,
    source_identity: dict[str, str],
    r_library: Path,
    current_paths: list[str],
) -> None:
    environment = os.environ.copy()
    environment["R_LIBS"] = os.pathsep.join(
        [str(r_library), *current_paths]
    )
    result = subprocess.run(
        [
            "Rscript",
            str(R_ADAPTER),
            "--input",
            str(output),
            "--output",
            str(output),
            "--plan-sha256",
            plan_sha256,
            "--mfrmr-git-head",
            source_identity["git_head"],
            "--mfrmr-status-sha256",
            source_identity["git_status_sha256"],
            "--mfrmr-content-sha256",
            source_identity["build_content_sha256"],
        ],
        env=environment,
        capture_output=True,
        text=True,
    )
    (output / "r_adapter_stdout.txt").write_text(result.stdout, encoding="utf-8")
    (output / "r_adapter_stderr.txt").write_text(result.stderr, encoding="utf-8")
    if result.returncode != 0:
        raise RuntimeError("R adapter failed; see retained stdout/stderr.")


def run_sirt_workers(
    output: Path,
    r_library: Path,
    current_paths: list[str],
) -> None:
    environment = os.environ.copy()
    environment["R_LIBS"] = os.pathsep.join(
        [str(r_library), *current_paths]
    )
    manifest = pd.read_csv(output / "manifest.csv")
    python_runs = pd.read_csv(output / "python_runs.csv")
    worker_dir = output / "sirt_workers"
    worker_dir.mkdir()
    rows: list[pd.DataFrame] = []
    for manifest_row in manifest.itertuples(index=False):
        case_id = str(manifest_row.CaseId)
        result_path = worker_dir / f"{case_id}.csv"
        result = subprocess.run(
            [
                "Rscript",
                str(SIRT_WORKER),
                "--input",
                str(output),
                "--output",
                str(result_path),
                "--case-id",
                case_id,
            ],
            env=environment,
            capture_output=True,
            text=True,
        )
        (worker_dir / f"{case_id}.stdout.txt").write_text(
            result.stdout, encoding="utf-8"
        )
        (worker_dir / f"{case_id}.stderr.txt").write_text(
            result.stderr, encoding="utf-8"
        )
        if result.returncode == 0 and result_path.exists():
            rows.append(pd.read_csv(result_path))
            continue
        python_row = python_runs.loc[
            python_runs["CaseId"].eq(case_id)
        ].iloc[0]
        rows.append(
            pd.DataFrame(
                [
                    {
                        "CaseId": case_id,
                        "Family": manifest_row.Family,
                        "Mechanism": manifest_row.Mechanism,
                        "ExpectedStatus": manifest_row.ExpectedStatus,
                        "PythonWorkflowStatus": python_row["WorkflowStatus"],
                        "Engine": "sirt",
                        "Estimator": "MML",
                        "Mode": "rm.facets_isolated_worker",
                        "FitAttempted": True,
                        "FitReturned": False,
                        "Converged": False,
                        "Iterations": np.nan,
                        "LogLik": np.nan,
                        "MaxAbsEstimate": np.nan,
                        "MaxAbsSE": np.nan,
                        "NonfiniteEstimateCount": np.nan,
                        "NonfiniteSECount": np.nan,
                        "ObservedSupportMax": np.nan,
                        "DeclaredRatingMax": manifest_row.RatingMax,
                        "SupportMismatch": False,
                        "Warnings": "",
                        "Error": (
                            f"isolated_process_exit_{result.returncode}: "
                            + result.stderr.strip().replace("\n", " | ")[:2000]
                        ),
                        "ElapsedSeconds": np.nan,
                    }
                ]
            )
        )
    sirt = pd.concat(rows, ignore_index=True)
    secondary = pd.read_csv(output / "secondary_engine_runs.csv")
    combined = pd.concat([secondary, sirt], ignore_index=True, sort=False)
    write_csv(combined, output / "secondary_engine_runs.csv")


def analyze_results(
    plan: dict[str, object],
    manifest: pd.DataFrame,
    python_runs: pd.DataFrame,
    python_coefficients: pd.DataFrame,
    output: Path,
) -> tuple[dict[str, object], pd.DataFrame, pd.DataFrame]:
    immer = pd.read_csv(output / "immer_runs.csv")
    immer_coefficients = pd.read_csv(output / "immer_coefficients.csv")
    secondary = pd.read_csv(output / "secondary_engine_runs.csv")
    identity = pd.read_csv(output / "r_engine_identity.csv")
    primary = immer.loc[immer["MaxIt"].eq(2000)].copy()
    interior_ids = set(
        manifest.loc[
            manifest["ExpectedStatus"].eq("interior_finite_cmle_supported"),
            "CaseId",
        ]
    )
    boundary_ids = set(
        manifest.loc[
            manifest["ExpectedStatus"].eq("boundary_no_finite_cmle"),
            "CaseId",
        ]
    )
    structural_ids = set(
        manifest.loc[
            manifest["ExpectedStatus"].eq("structural_nonidentification"),
            "CaseId",
        ]
    )
    interior_primary = primary.loc[primary["CaseId"].isin(interior_ids)].copy()
    boundary = immer.loc[immer["CaseId"].isin(boundary_ids)].copy()
    structural = immer.loc[immer["CaseId"].isin(structural_ids)].copy()
    false_reassurance = (
        boundary.assign(
            CodeZero=boundary["ConvergenceCode"].eq(0),
            FiniteVector=boundary["AllCoefficientsFinite"].fillna(False),
            FiniteSE=boundary["AllSEFinite"].fillna(False),
            SmallGradient1e6=boundary["GradientReadyAt1e6"].fillna(False),
            FullInformationRank=boundary["InformationNullity"].eq(0),
            NominalComposite=(
                boundary["ConvergenceCode"].eq(0)
                & boundary["AllCoefficientsFinite"].fillna(False)
                & boundary["AllSEFinite"].fillna(False)
                & boundary["GradientReadyAt1e6"].fillna(False)
                & boundary["InformationNullity"].eq(0)
            ),
        )
        .groupby("MaxIt", as_index=False)
        .agg(
            BoundaryCases=("CaseId", "size"),
            CodeZero=("CodeZero", "sum"),
            FiniteVector=("FiniteVector", "sum"),
            FiniteSE=("FiniteSE", "sum"),
            SmallGradient1e6=("SmallGradient1e6", "sum"),
            FullInformationRank=("FullInformationRank", "sum"),
            NominalComposite=("NominalComposite", "sum"),
            MedianMaxAbsEstimate=("MaxAbsEstimate", "median"),
            MaximumMaxAbsEstimate=("MaxAbsEstimate", "max"),
        )
    )
    write_csv(false_reassurance, output / "boundary_false_reassurance.csv")
    boundary_detail = boundary[
        [
            "CaseId",
            "MaxIt",
            "PythonWorkflowStatus",
            "FitAttempted",
            "FitReturned",
            "ConvergenceCode",
            "AllCoefficientsFinite",
            "AllSEFinite",
            "GradientSupNorm",
            "GradientReadyAt1e4",
            "GradientReadyAt1e5",
            "GradientReadyAt1e6",
            "GradientReadyAt1e8",
            "InformationRank",
            "InformationNullity",
            "InformationMinEigenvalue",
            "MaxAbsEstimate",
            "MaxAbsSE",
            "Warnings",
            "Error",
        ]
    ].copy()
    boundary_detail["NominalCompositeAt1e6"] = (
        boundary_detail["ConvergenceCode"].eq(0)
        & boundary_detail["AllCoefficientsFinite"].fillna(False)
        & boundary_detail["AllSEFinite"].fillna(False)
        & boundary_detail["GradientReadyAt1e6"].fillna(False)
        & boundary_detail["InformationNullity"].eq(0)
    )
    boundary_detail["FinalReadinessAfterOracle"] = False
    boundary_detail["FinalDecision"] = "blocked_no_finite_cmle"
    write_csv(boundary_detail, output / "boundary_maxit_sensitivity.csv")

    primary_coefficients = immer_coefficients.loc[
        immer_coefficients["MaxIt"].eq(2000)
        & immer_coefficients["CaseId"].isin(interior_ids)
    ].copy()
    comparison = python_coefficients.merge(
        primary_coefficients,
        on=["CaseId", "Parameter"],
        how="outer",
        suffixes=("Python", "Immer"),
        indicator=True,
    )
    comparison["EstimateAbsoluteDifference"] = (
        comparison["EstimatePython"] - comparison["EstimateImmer"]
    ).abs()
    comparison["SEAbsoluteDifference"] = (
        comparison["SEPython"] - comparison["SEImmer"]
    ).abs()
    write_csv(comparison, output / "immer_interior_coefficient_comparison.csv")
    loglik = python_runs.loc[
        python_runs["CaseId"].isin(interior_ids),
        ["CaseId", "ConditionalLogLik"],
    ].merge(
        interior_primary[["CaseId", "ConditionalLogLik"]],
        on="CaseId",
        suffixes=("Python", "Immer"),
        validate="one_to_one",
    )
    loglik["AbsoluteDifference"] = (
        loglik["ConditionalLogLikPython"] - loglik["ConditionalLogLikImmer"]
    ).abs()
    write_csv(loglik, output / "immer_interior_loglik_comparison.csv")

    expected = plan["fixture_contract"]["expected_accounting"]
    identity_row = identity.iloc[0]
    effective_identity = plan["_effective_identity_amendment"][
        "effective_mfrmr_source_identity"
    ]
    identity_passed = bool(
        identity_row["PlanSHA256"] == sha256_file(DEFAULT_PLAN)
        and str(identity_row["MfrmrVersion"]) == "0.2.3"
        and identity_row["MfrmrSourceGitHead"]
        == effective_identity["git_head_at_snapshot"]
        and identity_row["MfrmrSourceStatusSHA256"]
        == effective_identity["git_status_sha256_at_snapshot"]
        and identity_row["MfrmrSourceContentSHA256"]
        == effective_identity["snapshot_relative_content_sha256"]
    )
    python_passed = bool(
        len(python_runs) == int(expected["total"])
        and python_runs["TerminalMatch"].all()
        and int(python_runs["OptimizationAttempted"].sum())
        == int(expected["interior_finite_cmle_supported"])
    )
    immer_prefit_passed = bool(
        len(immer) == int(expected["total"]) * 3
        and immer["PrefitRankAgreement"].all()
        and not structural["FitAttempted"].any()
    )
    interior_diagnostics_passed = bool(
        len(interior_primary) == int(expected["interior_finite_cmle_supported"])
        and interior_primary["FitReturned"].all()
        and interior_primary["ConvergenceCode"].eq(0).all()
        and interior_primary["AllCoefficientsFinite"].all()
        and interior_primary["AllSEFinite"].all()
        and interior_primary["InformationNullity"].eq(0).all()
        and interior_primary["GradientSupNorm"].le(1e-6).all()
    )
    coefficient_max = float(comparison["EstimateAbsoluteDifference"].max())
    se_max = float(comparison["SEAbsoluteDifference"].max())
    loglik_max = float(loglik["AbsoluteDifference"].max())
    numeric_passed = bool(
        comparison["_merge"].eq("both").all()
        and coefficient_max <= 1e-7
        and loglik_max <= 1e-8
    )
    secondary_complete = bool(
        len(secondary) == int(expected["total"]) * 3
        and secondary.groupby("Engine")["CaseId"].nunique().eq(int(expected["total"])).all()
    )
    secondary_by_status = (
        secondary.assign(
            HasError=secondary["Error"].fillna("").ne(""),
            HasWarning=secondary["Warnings"].fillna("").ne(""),
        )
        .groupby(["Engine", "ExpectedStatus"], as_index=False)
        .agg(
            Cases=("CaseId", "nunique"),
            FitAttempted=("FitAttempted", "sum"),
            FitReturned=("FitReturned", "sum"),
            Converged=("Converged", "sum"),
            Errors=("HasError", "sum"),
            Warnings=("HasWarning", "sum"),
            MedianMaxAbsEstimate=("MaxAbsEstimate", "median"),
            MaximumMaxAbsEstimate=("MaxAbsEstimate", "max"),
        )
    )
    write_csv(secondary_by_status, output / "secondary_engine_by_status.csv")
    boundary_complete = len(boundary) == int(expected["boundary_no_finite_cmle"]) * 3
    return (
        {
            "identity_passed": identity_passed,
            "python_passed": python_passed,
            "immer_prefit_passed": immer_prefit_passed,
            "immer_interior_diagnostics_passed": interior_diagnostics_passed,
            "immer_numeric_passed": numeric_passed,
            "boundary_complete": bool(boundary_complete),
            "secondary_complete": secondary_complete,
            "estimate_max_absolute_difference": coefficient_max,
            "se_max_absolute_difference_descriptive": se_max,
            "conditional_loglik_max_absolute_difference": loglik_max,
            "boundary_false_reassurance": false_reassurance.to_dict(orient="records"),
            "secondary_summary": (
                secondary.groupby("Engine", as_index=False)
                .agg(
                    Cases=("CaseId", "nunique"),
                    FitReturned=("FitReturned", "sum"),
                    Converged=("Converged", "sum"),
                    Errors=("Error", lambda value: int(value.fillna("").ne("").sum())),
                )
                .to_dict(orient="records")
            ),
            "secondary_by_python_status": secondary_by_status.to_dict(
                orient="records"
            ),
        },
        comparison,
        loglik,
    )


def run_tests(output: Path) -> dict[str, object]:
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "tests/test_cmle_cross_engine_boundary_contract.py",
            "tests/test_cmle_workflow.py",
            "tests/test_cmle_existence.py",
            "tests/test_cmle.py",
        ],
        cwd=ROOT,
        capture_output=True,
        text=True,
        env={**os.environ, "MPLCONFIGDIR": "/private/tmp/matplotlib-cmle-boundary"},
    )
    (output / "selected_tests_stdout.txt").write_text(result.stdout, encoding="utf-8")
    (output / "selected_tests_stderr.txt").write_text(result.stderr, encoding="utf-8")
    return {"passed": result.returncode == 0, "returncode": result.returncode}


def write_critical_review(output: Path, decision: dict[str, object]) -> None:
    false_rows = decision["results"]["boundary_false_reassurance"]
    primary = next(row for row in false_rows if int(row["MaxIt"]) == 2000)
    secondary = decision["results"]["secondary_summary"]
    secondary_by_status = decision["results"]["secondary_by_python_status"]
    mfrmr_boundary = next(
        row
        for row in secondary_by_status
        if row["Engine"] == "mfrmr"
        and row["ExpectedStatus"] == "boundary_no_finite_cmle"
    )
    text = f"""# CMLE Cross-Engine Boundary Semantics: Critical Review

## Decision

The registered contract **{'passed' if decision['contract_passed'] else 'did not pass'}**. The primary comparison is deliberately narrow: native Python exact CMLE versus a matched `immer::immer_cml` design on five oracle-qualified interior cases. TAM, sirt, and mfrmr are not CMLE parity engines in this experiment.

## What the boundary controls show

At the primary `immer` cap (`maxit=2000`), {int(primary['CodeZero'])}/7 boundary cases returned optimizer code zero, {int(primary['FiniteVector'])}/7 returned finite coefficient vectors, and {int(primary['FiniteSE'])}/7 returned finite standard errors. Even the combined code-zero/finite-vector/finite-SE/gradient-at-1e-6/full-information rule was satisfied in {int(primary['NominalComposite'])}/7. These observations do not overturn the exact support-oracle result that no finite CMLE exists. The key product implication is that a green optimizer message cannot be the public readiness gate.

Increasing `maxit` is a diagnostic, not a proof. Coefficient growth, information deterioration, or gradient changes support a boundary interpretation, but the finite-existence decision comes from the conditional-support geometry before optimization.

## Matched interior comparison

Across matched free coordinates, the maximum absolute Python–immer estimate difference was {decision['results']['estimate_max_absolute_difference']:.17g}; the maximum conditional-log-likelihood difference was {decision['results']['conditional_loglik_max_absolute_difference']:.17g}. Standard-error differences are retained descriptively because covariance construction and floating-point Hessian details can differ even when the estimand and optimum match.

## Why the other engines are not an equivalence test

The secondary ledger is {secondary}. In particular, the frozen `mfrmr::fit_mfrm(method='JML')` snapshot returned and reported convergence for {int(mfrmr_boundary['Converged'])}/{int(mfrmr_boundary['Cases'])} Python-CMLE boundary cases; the largest retained absolute estimate in that group was {float(mfrmr_boundary['MaximumMaxAbsEstimate']):.8g}. This is not a contradiction: mfrmr JML estimates Person parameters jointly, while TAM and sirt integrate over a Person distribution. Their convergence, warnings, or finite outputs answer how alternative workflows behave on the same response patterns, not whether they reproduce exact CMLE. Treating their numbers as expected CMLE equality would be a category error.

## Floating-point and decision boundaries

All raw binary64 values are retained. `immer` gradient classifications are shown at 1e-4, 1e-5, 1e-6, and 1e-8. No displayed rounding is used for the primary decision, and MnSq is not used to infer finite-MLE existence. Later UI work must display a review band near any fit threshold and expose raw-versus-rounded decision sensitivity separately.

## Remaining threats

- The suite is intentionally small and additive-RSM only; PCM needs a separately registered matched threshold design.
- No hard-anchor behavior is compared here. Anchors change the free-coordinate cone and therefore the existence question.
- The tested mfrmr 0.2.3 tree is dirty and unreleased; conclusions attach only to the retained content hash.
- A mathematically correct gate can still be poorly understood. Bilingual strings have not undergone a user-comprehension study.
- Sparse high-dimensional topology, biased anchors, DIF/bias procedures, visualization, and one-click Person scoring remain separate validation layers.

## Product recommendation

Keep the pre-optimization support-oracle gate mandatory and fail closed. In a future one-click screen, show one terminal status first, then a concise reason and next action; place optimizer diagnostics and cross-engine sensitivity under an expandable technical-evidence section. Do not display boundary-case coefficients as ordinary estimates.
"""
    (output / "CMLE_CROSS_ENGINE_BOUNDARY_CRITICAL_REVIEW.md").write_text(
        text, encoding="utf-8"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--plan", type=Path, default=DEFAULT_PLAN)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    plan = validate_plan(args.plan)
    if args.output.exists():
        raise FileExistsError(f"Evidence output already exists: {args.output}")
    args.output.mkdir(parents=True)
    try:
        cases = selected_fixture_cases(plan)
        effective_identity = plan["_effective_identity_amendment"][
            "effective_mfrmr_source_identity"
        ]
        frozen_source = Path(effective_identity["temporary_snapshot_path"])
        frozen_hash = mfrmr_source_identity_without_git(frozen_source)
        source_before = {
            "git_head": effective_identity["git_head_at_snapshot"],
            "git_status_sha256": effective_identity[
                "git_status_sha256_at_snapshot"
            ],
            "build_content_sha256": frozen_hash,
            "build_relative_content_sha256": frozen_hash,
            "build_file_count": str(effective_identity["build_relevant_file_count"]),
        }
        assert_mfrmr_identity(source_before, plan)
        manifest, _, python_runs, python_coefficients = build_python_bundle(
            cases, args.output
        )
        library, current_paths = install_mfrmr(
            frozen_source, source_before, args.output
        )
        if mfrmr_source_identity_without_git(frozen_source) != frozen_hash:
            raise ValueError("Frozen mfrmr snapshot changed during installation.")
        (args.output / "mfrmr_source_identity.json").write_text(
            json.dumps(source_before, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        run_r_adapter(
            args.output,
            sha256_file(args.plan),
            source_before,
            library,
            current_paths,
        )
        run_sirt_workers(args.output, library, current_paths)
        results, _, _ = analyze_results(
            plan, manifest, python_runs, python_coefficients, args.output
        )
        tests = run_tests(args.output)
        contract_passed = bool(
            all(
                results[key]
                for key in (
                    "identity_passed",
                    "python_passed",
                    "immer_prefit_passed",
                    "immer_interior_diagnostics_passed",
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
            "identity_amendment_sha256": sha256_file(DEFAULT_AMENDMENT),
            "frozen_snapshot_amendment_sha256": sha256_file(
                DEFAULT_AMENDMENT2
            ),
            "implementation_sha256": {
                "validation/cmle_cross_engine_boundary.py": sha256_file(
                    Path(__file__).resolve()
                ),
                "validation/cmle_cross_engine_boundary.R": sha256_file(
                    R_ADAPTER
                ),
                "validation/cmle_cross_engine_sirt_worker.R": sha256_file(
                    SIRT_WORKER
                ),
                "tests/test_cmle_cross_engine_boundary_contract.py": sha256_file(
                    ROOT / "tests/test_cmle_cross_engine_boundary_contract.py"
                ),
            },
            "contract_passed": contract_passed,
            "results": results,
            "selected_tests": tests,
            "interpretation": {
                "matched_parity_scope": "Python native exact CMLE versus immer CMLE, interior cases only",
                "boundary_rule": "support-oracle status cannot be overwritten by finite optimizer output or nominal convergence",
                "secondary_scope": "mfrmr JML, TAM MML, and sirt MML are descriptive sensitivity analyses only",
                "public_ui": "withheld",
            },
        }
        write_critical_review(args.output, decision)
        output_names = sorted(
            path.relative_to(args.output).as_posix()
            for path in args.output.rglob("*")
            if path.is_file() and path.name != "decision.json"
        )
        decision["output_sha256"] = {
            name: sha256_file(args.output / name) for name in output_names
        }
        decision = json_safe(decision)
        (args.output / "decision.json").write_text(
            json.dumps(decision, indent=2, ensure_ascii=False, allow_nan=False)
            + "\n",
            encoding="utf-8",
        )
        print(json.dumps(decision, indent=2, ensure_ascii=False, allow_nan=False))
        if not contract_passed:
            raise SystemExit(1)
    except Exception:
        # Keep typed installation/adapter logs for diagnosis but never replace a
        # failed run with a partial decision artifact.
        raise


if __name__ == "__main__":
    main()
