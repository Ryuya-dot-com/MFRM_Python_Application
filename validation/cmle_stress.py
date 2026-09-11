#!/usr/bin/env python3
"""Run the native exact CMLE core on the frozen 2026-08-09 stress datasets.

This runner reuses the deterministic RSM/PCM datasets and generating-truth
surfaces created by ``cross_engine_stress.py``.  It is a repository-only
feasibility and failure-behavior pilot.  One or two replicates per cell do not
support estimator rankings, coverage claims, or sample-size recommendations.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import platform
import subprocess
import sys
import time

import numpy as np
import pandas as pd
import scipy

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from mfrm_app.cmle import (
    CMLEEligibilityError,
    DEFAULT_RANK_AUDIT_MAX_BYTES,
    fit_cmle,
    prepare_cmle_design,
)


DEFAULT_INPUT = Path(__file__).resolve().parent / "cross_engine_stress_20260809"
DEFAULT_OUTPUT = Path(__file__).resolve().parent / "cmle_phase0_20260809"


def _git_value(*args: str) -> str:
    completed = subprocess.run(
        ["git", *args], capture_output=True, text=True, check=False
    )
    return completed.stdout.strip() if completed.returncode == 0 else "unavailable"


def _surface_metrics(fit: dict, truth: pd.DataFrame) -> tuple[float, float, int]:
    estimated = fit["surfaces"][
        ["Rater", "Criterion", "Category", "CumulativeDifficulty"]
    ].copy()
    joined = truth.merge(
        estimated,
        on=["Rater", "Criterion", "Category"],
        how="left",
        validate="one_to_one",
    )
    if joined["CumulativeDifficulty"].isna().any():
        raise RuntimeError("CMLE stress surface merge contains missing estimates.")
    error = joined["CumulativeDifficulty"] - joined["Truth"]
    return (
        float(np.sqrt(np.mean(np.square(error)))),
        float(np.max(np.abs(error))),
        int(len(joined)),
    )


def run(
    input_dir: Path,
    output_dir: Path,
    rank_audit_max_work: int,
    rank_audit_max_bytes: int,
) -> None:
    manifest_path = input_dir / "manifest.csv"
    ratings_path = input_dir / "simulated_ratings.csv"
    truth_path = input_dir / "truth_surface.csv"
    for path in (manifest_path, ratings_path, truth_path):
        if not path.exists():
            raise FileNotFoundError(f"Required stress input is missing: {path}")

    manifest = pd.read_csv(manifest_path)
    ratings = pd.read_csv(ratings_path)
    truth = pd.read_csv(truth_path)
    expected_ids = manifest["DatasetId"].astype(str).tolist()
    if set(ratings["DatasetId"].astype(str)) != set(expected_ids):
        raise RuntimeError("Ratings DatasetId coverage does not match manifest.")
    if set(truth["DatasetId"].astype(str)) != set(expected_ids):
        raise RuntimeError("Truth DatasetId coverage does not match manifest.")

    output_dir.mkdir(parents=True, exist_ok=True)
    fit_rows: list[dict[str, object]] = []
    for manifest_row in manifest.itertuples(index=False):
        dataset_id = str(manifest_row.DatasetId)
        model = str(manifest_row.Model)
        frame = ratings.loc[ratings["DatasetId"] == dataset_id].copy()
        truth_frame = truth.loc[truth["DatasetId"] == dataset_id].copy()
        started = time.perf_counter()
        common = dict(
            person_col="Person",
            facet_cols=["Rater", "Criterion"],
            score_col="Score",
            rating_min=0,
            rating_max=int(manifest_row.Categories) - 1,
            model=model,
            step_facet="Criterion" if model == "PCM" else None,
            rank_audit_max_work=int(rank_audit_max_work),
            rank_audit_max_bytes=int(rank_audit_max_bytes),
        )
        row: dict[str, object] = {
            "DatasetId": dataset_id,
            "Model": model,
            "Scenario": str(manifest_row.Scenario),
            "Replicate": int(manifest_row.Replicate),
            "ObservedRatings": int(len(frame)),
            "ExpectedConnected": bool(manifest_row.ExpectedConnected),
            "RankAuditMaxWork": int(rank_audit_max_work),
            "RankAuditMaxBytes": int(rank_audit_max_bytes),
            "InformationMethod": "exact conditional sufficient-statistic covariance",
        }
        try:
            fit = fit_cmle(
                frame,
                maxiter=1_000,
                gtol=1e-8,
                **common,
            )
            summary = fit["summary"].iloc[0]
            rmse, maximum, coordinates = _surface_metrics(fit, truth_frame)
            row.update(
                {
                    "Status": "ok",
                    "Converged": bool(summary["Converged"]),
                    "OptimizerSuccess": bool(summary["OptimizerSuccess"]),
                    "OptimizerGradientSupNorm": float(
                        summary["OptimizerGradientSupNorm"]
                    ),
                    "NewtonPolishSteps": int(summary["NewtonPolishSteps"]),
                    "NewtonPolishAcceptedSteps": int(
                        summary["NewtonPolishAcceptedSteps"]
                    ),
                    "NewtonPolishReason": str(summary["NewtonPolishReason"]),
                    "InferenceReady": bool(summary["InferenceReady"]),
                    "PersonsInformative": int(summary["PersonsInformative"]),
                    "PersonsExtreme": int(summary["PersonsExtreme"]),
                    "Patterns": int(summary["Patterns"]),
                    "KParams": int(summary["KParams"]),
                    "ConditionalLogLik": float(summary["ConditionalLogLik"]),
                    "GradientSupNorm": float(summary["GradientSupNorm"]),
                    "InformationRank": int(summary["InformationRank"]),
                    "InformationNullity": int(summary["InformationNullity"]),
                    "InformationMethod": str(summary["InformationMethod"]),
                    "InformationConditionNumber": float(
                        summary["InformationConditionNumber"]
                    ),
                    "InformationMinEigenvalue": float(
                        summary["InformationMinEigenvalue"]
                    ),
                    "RankAuditWorkProxy": int(
                        fit["audit"]["rank_audit_work_proxy"]
                    ),
                    "RankAuditPeakMiBProxy": float(
                        fit["audit"]["rank_audit_peak_bytes_proxy"] / 1024**2
                    ),
                    "SurfaceCoordinates": coordinates,
                    "SurfaceRMSE": rmse,
                    "SurfaceMaxAbsError": maximum,
                    "BlockingCodes": "",
                    "Message": str(summary["OptimizerMessage"]),
                }
            )
        except CMLEEligibilityError as exc:
            design = prepare_cmle_design(frame, **common)
            blocking = [
                issue["Code"]
                for issue in design.audit["issues"]
                if issue["Severity"] == "Block"
            ]
            row.update(
                {
                    "Status": "blocked",
                    "Converged": False,
                    "OptimizerSuccess": False,
                    "InferenceReady": False,
                    "PersonsInformative": int(design.audit["persons_informative"]),
                    "PersonsExtreme": int(design.audit["persons_extreme"]),
                    "Patterns": int(design.audit["unique_missingness_design_patterns"]),
                    "KParams": int(design.audit["parameters"]),
                    "ConditionalLogLik": np.nan,
                    "GradientSupNorm": np.nan,
                    "InformationRank": int(design.audit["conditional_rank"]),
                    "InformationNullity": int(design.audit["conditional_nullity"]),
                    "InformationConditionNumber": design.audit[
                        "conditional_condition_number"
                    ],
                    "InformationMinEigenvalue": design.audit[
                        "conditional_min_eigenvalue"
                    ],
                    "RankAuditWorkProxy": int(
                        design.audit["rank_audit_work_proxy"]
                    ),
                    "RankAuditPeakMiBProxy": float(
                        design.audit["rank_audit_peak_bytes_proxy"] / 1024**2
                    ),
                    "SurfaceCoordinates": 0,
                    "SurfaceRMSE": np.nan,
                    "SurfaceMaxAbsError": np.nan,
                    "BlockingCodes": ";".join(blocking),
                    "Message": str(exc),
                }
            )
        except Exception as exc:  # retained as a failed denominator
            row.update(
                {
                    "Status": "failed",
                    "Converged": False,
                    "OptimizerSuccess": False,
                    "InferenceReady": False,
                    "BlockingCodes": "unexpected_failure",
                    "Message": repr(exc),
                }
            )
        row["ElapsedSeconds"] = float(time.perf_counter() - started)
        fit_rows.append(row)

    fit_runs = pd.DataFrame(fit_rows)
    if fit_runs["DatasetId"].duplicated().any():
        raise RuntimeError("CMLE stress output has duplicate DatasetId rows.")
    if set(fit_runs["DatasetId"]) != set(expected_ids):
        raise RuntimeError("CMLE stress output lost expected DatasetId rows.")
    fit_runs.sort_values(["Model", "Scenario", "Replicate", "DatasetId"]).to_csv(
        output_dir / "stress_fit_runs.csv", index=False
    )

    scenario = (
        fit_runs.groupby(["Model", "Scenario"], dropna=False)
        .agg(
            ExpectedFits=("DatasetId", "size"),
            ReturnedFits=("Status", lambda values: int((values == "ok").sum())),
            BlockedFits=("Status", lambda values: int((values == "blocked").sum())),
            FailedFits=("Status", lambda values: int((values == "failed").sum())),
            InferenceReadyFits=("InferenceReady", lambda values: int(values.fillna(False).sum())),
            MedianInformativePersons=("PersonsInformative", "median"),
            MedianPatterns=("Patterns", "median"),
            MedianSurfaceRMSE=("SurfaceRMSE", "median"),
            MaximumSurfaceRMSE=("SurfaceRMSE", "max"),
            MedianElapsedSeconds=("ElapsedSeconds", "median"),
        )
        .reset_index()
    )
    scenario.to_csv(output_dir / "stress_scenario_summary.csv", index=False)

    identity = {
        "schema_version": "native-exact-cmle-stress-v2",
        "run_date": "2026-08-09",
        "git_head": _git_value("rev-parse", "HEAD"),
        "git_status_short": _git_value("status", "--short"),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "input_dir": str(input_dir.resolve()),
        "manifest_rows": int(len(manifest)),
        "output_rows": int(len(fit_runs)),
        "rank_audit_max_work": int(rank_audit_max_work),
        "rank_audit_max_bytes": int(rank_audit_max_bytes),
        "information_method": "exact_conditional_second_moment_dp",
    }
    (output_dir / "stress_runtime_identity.json").write_text(
        json.dumps(identity, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--rank-audit-max-work", type=int, default=500_000_000)
    parser.add_argument(
        "--rank-audit-max-bytes",
        type=int,
        default=DEFAULT_RANK_AUDIT_MAX_BYTES,
    )
    args = parser.parse_args()
    run(
        args.input_dir,
        args.output_dir,
        args.rank_audit_max_work,
        args.rank_audit_max_bytes,
    )


if __name__ == "__main__":
    main()
