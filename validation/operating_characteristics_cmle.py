#!/usr/bin/env python3
"""Fit native exact CMLE to the byte-validated OC rating bundle.

The adapter deliberately estimates an unanchored additive RSM structural
calibration.  Runs whose generating condition requested anchors are still
useful for byte-matched Python/immer numerical parity, but are marked
ineligible for the *requested anchored condition* because CMLE Phase 0 does
not implement anchors.  Local Rater x Task bias is also outside this additive
CMLE estimand and is never manufactured from a main-effect fit.
"""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
import sys
import time

import numpy as np
import pandas as pd


REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from mfrm_app.cmle import (  # noqa: E402
    CMLE_SCHEMA_VERSION,
    audit_cmle_eligibility,
    fit_cmle,
)


MODE = "PYTHON_EXACT_CMLE_UNANCHORED"
FACETS = ["Rater", "Task", "Criterion"]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_bundle(input_dir: Path) -> tuple[pd.DataFrame, str]:
    inventory_path = input_dir / "generated_bundle_files.csv"
    inventory = pd.read_csv(inventory_path)
    for row in inventory.itertuples(index=False):
        path = input_dir / str(row.File)
        if not path.is_file() or sha256_file(path).lower() != str(row.SHA256).lower():
            raise RuntimeError(f"Generated bundle byte identity failed: {row.File}")
        if len(pd.read_csv(path)) != int(row.Rows):
            raise RuntimeError(f"Generated bundle row count failed: {row.File}")
    inventory_hash = sha256_file(inventory_path)
    bridge = pd.read_csv(input_dir / "bridge_validation_r.csv")
    if not bridge["Passed"].astype(bool).all():
        raise RuntimeError("R bridge validation contains a failed check.")
    if not bridge["BundleInventorySHA256"].astype(str).eq(inventory_hash).all():
        raise RuntimeError("R bridge validation is stale for the generated bundle.")
    return inventory, inventory_hash


def base_run(manifest_row: pd.Series, rows: int, anchor_rows: int, bundle_hash: str) -> dict[str, object]:
    return {
        "SchemaVersion": str(manifest_row["SchemaVersion"]),
        "CMLESchemaVersion": CMLE_SCHEMA_VERSION,
        "RunId": str(manifest_row["RunId"]),
        "ConditionId": str(manifest_row["ConditionId"]),
        "Design": str(manifest_row["Design"]),
        "TruthBias": float(manifest_row["TruthBias"]),
        "TruthPositive": bool(manifest_row["TruthPositive"]),
        "Replicate": int(manifest_row["Replicate"]),
        "Seed": int(manifest_row["Seed"]),
        "Engine": "Python",
        "Estimator": "CMLE",
        "Mode": MODE,
        "Model": "RSM",
        "Estimand": "unanchored additive structural parameters after conditioning out Person",
        "Rows": int(rows),
        "RequestedAnchorRows": int(anchor_rows),
        "RequestedAnchorSupported": anchor_rows == 0,
        "BiasEstimandSupported": False,
        "ConditionalDesignEligible": False,
        "FitReturned": False,
        "Converged": False,
        "InferenceReady": False,
        "ParityEligible": False,
        "RequestedConditionEligible": False,
        "FailureStage": "prefit",
        "FailureReason": "",
        "ConditionalRank": np.nan,
        "ConditionalNullity": np.nan,
        "PersonsTotal": np.nan,
        "PersonsInformative": np.nan,
        "PersonsExtreme": np.nan,
        "Patterns": np.nan,
        "KParams": np.nan,
        "ConditionalLogLik": np.nan,
        "ConditionalDeviance": np.nan,
        "Iterations": np.nan,
        "FunctionEvaluations": np.nan,
        "OptimizerSuccess": False,
        "GradientSupNorm": np.nan,
        "StationarityTolerance": np.nan,
        "InformationConditionNumber": np.nan,
        "InformationMinEigenvalue": np.nan,
        "ReadinessReasons": "",
        "AuditReasonCodes": "",
        "ElapsedSeconds": np.nan,
        "BundleInventorySHA256": bundle_hash,
        "CMLESourceSHA256": sha256_file(REPO / "mfrm_app" / "cmle.py"),
    }


def recovery_rows(
    manifest_row: pd.Series,
    fit: dict[str, object],
    truth: pd.DataFrame,
    parity_eligible: bool,
) -> pd.DataFrame:
    estimates = fit["facets"]["others"][["Facet", "Level", "Estimate", "SE"]].copy()
    expected = truth.loc[
        truth["RunId"].astype(str).eq(str(manifest_row["RunId"])),
        ["Facet", "Level", "Truth"],
    ]
    merged = expected.merge(estimates, on=["Facet", "Level"], how="left", validate="one_to_one")
    merged["EstimateAligned"] = np.nan
    merged["TruthAligned"] = np.nan
    for _, indices in merged.groupby("Facet", sort=False).groups.items():
        idx = list(indices)
        merged.loc[idx, "EstimateAligned"] = (
            merged.loc[idx, "Estimate"] - merged.loc[idx, "Estimate"].mean()
        )
        merged.loc[idx, "TruthAligned"] = (
            merged.loc[idx, "Truth"] - merged.loc[idx, "Truth"].mean()
        )
    merged["ErrorAligned"] = merged["EstimateAligned"] - merged["TruthAligned"]
    merged["ComparisonScale"] = "mean_aligned_unanchored_location"
    merged["IncludedInSummary"] = bool(parity_eligible) & merged["ErrorAligned"].notna()
    merged.insert(0, "RunId", str(manifest_row["RunId"]))
    merged.insert(1, "ConditionId", str(manifest_row["ConditionId"]))
    merged.insert(2, "Engine", "Python")
    merged.insert(3, "Estimator", "CMLE")
    merged.insert(4, "Mode", MODE)
    merged["Design"] = str(manifest_row["Design"])
    merged["TruthBias"] = float(manifest_row["TruthBias"])
    merged["Replicate"] = int(manifest_row["Replicate"])
    merged["Seed"] = int(manifest_row["Seed"])
    return merged


def run_adapter(input_dir: Path, output_dir: Path) -> None:
    _, bundle_hash = validate_bundle(input_dir)
    manifest = pd.read_csv(input_dir / "manifest.csv")
    ratings = pd.read_csv(input_dir / "generated_ratings.csv")
    truth = pd.read_csv(input_dir / "generated_facet_truth.csv")
    anchors = pd.read_csv(input_dir / "generated_anchors.csv")
    if manifest["RunId"].duplicated().any():
        raise RuntimeError("Manifest RunId values must be unique.")

    run_rows: list[dict[str, object]] = []
    coefficient_frames: list[pd.DataFrame] = []
    recovery_frames: list[pd.DataFrame] = []
    for _, manifest_row in manifest.iterrows():
        run_id = str(manifest_row["RunId"])
        data = ratings.loc[
            ratings["RunId"].astype(str).eq(run_id),
            ["Person", *FACETS, "Score"],
        ].copy()
        run_anchor_rows = int(anchors["RunId"].astype(str).eq(run_id).sum())
        run = base_run(manifest_row, len(data), run_anchor_rows, bundle_hash)
        started = time.perf_counter()
        audit = audit_cmle_eligibility(
            data,
            person_col="Person",
            facet_cols=FACETS,
            score_col="Score",
            rating_min=0,
            rating_max=int(manifest_row["Categories"]) - 1,
            model="RSM",
        )
        run["ConditionalDesignEligible"] = bool(audit.get("eligible", False))
        run["ConditionalRank"] = audit.get("conditional_rank", np.nan)
        run["ConditionalNullity"] = audit.get("conditional_nullity", np.nan)
        run["PersonsTotal"] = audit.get("persons_total", np.nan)
        run["PersonsInformative"] = audit.get("persons_informative", np.nan)
        run["PersonsExtreme"] = audit.get("persons_extreme", np.nan)
        run["Patterns"] = audit.get("unique_missingness_design_patterns", np.nan)
        run["KParams"] = audit.get("parameters", np.nan)
        issues = audit.get("issues_table", pd.DataFrame())
        reason_codes = issues.get("Code", pd.Series(dtype=str)).astype(str).tolist()
        run["AuditReasonCodes"] = ";".join(reason_codes)
        if not bool(audit.get("eligible", False)):
            blocking = issues.loc[issues.get("Severity", pd.Series(dtype=str)).eq("Block"), "Message"]
            run["FailureReason"] = " | ".join(blocking.astype(str).tolist()) or "CMLE prefit audit failed"
            run["ElapsedSeconds"] = time.perf_counter() - started
            run_rows.append(run)
            continue

        try:
            fit = fit_cmle(
                data,
                person_col="Person",
                facet_cols=FACETS,
                score_col="Score",
                rating_min=0,
                rating_max=int(manifest_row["Categories"]) - 1,
                model="RSM",
                maxiter=2000,
                gtol=1e-9,
                newton_polish_maxiter=12,
            )
        except Exception as exc:  # fail closed and retain the attempted run
            run["FailureStage"] = "fit"
            run["FailureReason"] = f"{type(exc).__name__}: {exc}"
            run["ElapsedSeconds"] = time.perf_counter() - started
            run_rows.append(run)
            continue

        summary = fit["summary"].iloc[0]
        run.update(
            {
                "FitReturned": True,
                "Converged": bool(summary["Converged"]),
                "InferenceReady": bool(summary["InferenceReady"]),
                "ConditionalLogLik": float(summary["ConditionalLogLik"]),
                "ConditionalDeviance": float(summary["ConditionalDeviance"]),
                "Iterations": int(summary["Iterations"]),
                "FunctionEvaluations": int(summary["FunctionEvaluations"]),
                "OptimizerSuccess": bool(summary["OptimizerSuccess"]),
                "GradientSupNorm": float(summary["GradientSupNorm"]),
                "StationarityTolerance": float(summary["StationarityTolerance"]),
                "InformationConditionNumber": float(summary["InformationConditionNumber"]),
                "InformationMinEigenvalue": float(summary["InformationMinEigenvalue"]),
                "ReadinessReasons": str(summary["ReadinessReasons"]),
            }
        )
        parity_eligible = bool(summary["Converged"] and summary["InferenceReady"])
        run["ParityEligible"] = parity_eligible
        run["RequestedConditionEligible"] = parity_eligible and run_anchor_rows == 0
        if not parity_eligible:
            run["FailureStage"] = "readiness"
            run["FailureReason"] = str(summary["ReadinessReasons"]) or "CMLE inference readiness withheld"
        elif run_anchor_rows:
            run["FailureStage"] = "requested_scope"
            run["FailureReason"] = "anchors requested by condition but unsupported by CMLE Phase 0; unanchored parity fit retained"
        else:
            run["FailureStage"] = ""
            run["FailureReason"] = ""

        coefficients = fit["coefficients"].copy()
        coefficients.insert(0, "RunId", run_id)
        coefficients.insert(1, "ConditionId", str(manifest_row["ConditionId"]))
        coefficients.insert(2, "Engine", "Python")
        coefficients.insert(3, "Estimator", "CMLE")
        coefficients.insert(4, "Mode", MODE)
        coefficients["IncludedInParity"] = parity_eligible
        coefficient_frames.append(coefficients)
        recovery_frames.append(recovery_rows(manifest_row, fit, truth, parity_eligible))
        run["ElapsedSeconds"] = time.perf_counter() - started
        run_rows.append(run)

    output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(run_rows).to_csv(output_dir / "python_cmle_runs.csv", index=False)
    pd.concat(coefficient_frames, ignore_index=True).to_csv(
        output_dir / "python_cmle_coefficients.csv", index=False
    )
    pd.concat(recovery_frames, ignore_index=True).to_csv(
        output_dir / "python_cmle_parameter_recovery.csv", index=False
    )
    identity = pd.DataFrame(
        [
            {
                "BundleInventorySHA256": bundle_hash,
                "CMLESchemaVersion": CMLE_SCHEMA_VERSION,
                "CMLESourceSHA256": sha256_file(REPO / "mfrm_app" / "cmle.py"),
                "AdapterScriptSHA256": sha256_file(Path(__file__).resolve()),
                "Mode": MODE,
                "FacetOrder": "Rater|Task|Criterion",
                "RatingSupport": "0..Categories-1",
                "AnchorScope": "unsupported; unanchored parity fit only",
                "BiasScope": "unsupported by additive CMLE estimand",
            }
        ]
    )
    identity.to_csv(output_dir / "python_cmle_adapter_identity.csv", index=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run_adapter(args.input.resolve(), args.output.resolve())


if __name__ == "__main__":
    main()
