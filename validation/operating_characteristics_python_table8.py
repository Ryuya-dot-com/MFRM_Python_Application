#!/usr/bin/env python3
"""Extract native Python JMLE RSM steps for FACETS Table 8 comparison.

This adapter consumes an existing byte-retained operating-characteristics
bundle.  It does not generate or modify ratings.  Its fit controls reproduce
the registered Python JMLE pilot so that RSM thresholds can be joined to the
same RunId/category rows parsed from FACETS Table 8.
"""

from __future__ import annotations

import argparse
import json
import platform
import sys
import time
import warnings
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from validation.operating_characteristics_facets import sha256_file, validate_bundle  # noqa: E402
from validation.operating_characteristics_python_mml import (  # noqa: E402
    validate_registered_bundle_hashes,
)


SCHEMA_VERSION = "mfrm-python-jmle-table8-adapter-v1"
OUTPUT_FILES = (
    "python_jmle_step_runs.csv",
    "python_jmle_steps.csv",
    "python_jmle_step_adapter_identity.json",
)


def build_jmle_kwargs(categories: int) -> dict[str, Any]:
    return {
        "person_col": "Person",
        "facet_cols": ["Rater", "Task", "Criterion"],
        "score_col": "Score",
        "rating_min": 0,
        "rating_max": int(categories) - 1,
        "model": "RSM",
        "method": "JMLE",
        "noncenter_facet": "Person",
        "anchor_policy": "warn",
        "min_common_anchors": 2,
        "maxit": 160,
        "reltol": 1e-6,
        "keep_original": True,
    }


def normalize_steps(result: dict[str, Any], manifest_row: pd.Series) -> pd.DataFrame:
    source = result.get("steps", pd.DataFrame()).copy()
    if source.empty or not {"Step", "Estimate"}.issubset(source.columns):
        raise ValueError("Python JMLE result does not contain RSM Step/Estimate rows")
    output = pd.DataFrame({
        "RunId": manifest_row["RunId"],
        "ConditionId": manifest_row["ConditionId"],
        "Design": manifest_row["Design"],
        "TruthBias": manifest_row["TruthBias"],
        "TruthPositive": manifest_row.get("TruthPositive", np.nan),
        "BiasCell": manifest_row.get("BiasCell", ""),
        "Replicate": manifest_row["Replicate"],
        "Seed": manifest_row["Seed"],
        "Engine": "PythonApp",
        "Estimator": "JMLE",
        "Mode": "PYTHON_JMLE_REGISTERED_RSM_STEPS",
        "Step": source["Step"].astype(str),
        "Category": source["Step"].astype(str).str.extract(r"(\d+)$", expand=False).astype(int),
        "ThresholdEstimate": pd.to_numeric(source["Estimate"], errors="coerce"),
    })
    expected = set(range(1, int(manifest_row["Categories"])))
    observed = set(output["Category"])
    if observed != expected:
        raise ValueError(f"Python RSM step categories differ from expected: {observed} vs {expected}")
    if not np.isfinite(output["ThresholdEstimate"]).all():
        raise ValueError("Python RSM step estimates include non-finite values")
    return output


def run_adapter(args: argparse.Namespace) -> None:
    input_dir = args.input_dir.resolve()
    output_dir = args.output_dir.resolve()
    tables = validate_bundle(input_dir)
    registered_hashes = validate_registered_bundle_hashes(input_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    existing = [output_dir / filename for filename in OUTPUT_FILES if (output_dir / filename).exists()]
    if existing and not args.overwrite:
        raise FileExistsError(
            "Python JMLE step output exists; pass --overwrite to replace only: "
            + ", ".join(path.name for path in existing)
        )

    manifest = tables["manifest.csv"].copy()
    if args.run_id:
        requested = set(args.run_id)
        manifest = manifest[manifest["RunId"].astype(str).isin(requested)]
        missing = requested.difference(manifest["RunId"].astype(str))
        if missing:
            raise ValueError(f"Requested RunIds are absent from manifest: {sorted(missing)}")
    if args.limit is not None:
        manifest = manifest.head(args.limit)
    if manifest.empty:
        raise ValueError("No manifest rows selected")

    import streamlit_app as app  # pylint: disable=import-outside-toplevel

    ratings_all = tables["generated_ratings.csv"]
    anchors_all = tables["generated_anchors.csv"]
    run_rows: list[dict[str, Any]] = []
    step_parts: list[pd.DataFrame] = []
    for _, manifest_row in manifest.iterrows():
        started = time.perf_counter()
        run_id = str(manifest_row["RunId"])
        ratings = ratings_all[ratings_all["RunId"].astype(str).eq(run_id)].copy()
        anchors = anchors_all[anchors_all["RunId"].astype(str).eq(run_id)].copy()
        run: dict[str, Any] = {
            **{column: manifest_row[column] for column in manifest_row.index},
            "Engine": "PythonApp",
            "Estimator": "JMLE",
            "Mode": "PYTHON_JMLE_REGISTERED_RSM_STEPS",
            "FitReturned": False,
            "Converged": False,
            "AnalysisEligible": False,
            "FailureReason": "",
            "Warnings": "",
            "ElapsedSeconds": np.nan,
        }
        caught_messages: list[str] = []
        try:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                result = app.mfrm_estimate(
                    data=ratings[["Person", "Rater", "Task", "Criterion", "Score"]].copy(),
                    anchor_df=(
                        anchors[["Facet", "Level", "Anchor"]].copy()
                        if not anchors.empty else None
                    ),
                    **build_jmle_kwargs(int(manifest_row["Categories"])),
                )
                caught_messages = [str(item.message) for item in caught]
            summary = result["summary"].iloc[0]
            steps = normalize_steps(result, manifest_row)
            converged = bool(summary.get("Converged", False))
            run.update({
                "FitReturned": True,
                "Converged": converged,
                "AnalysisEligible": bool(converged and len(steps) == int(manifest_row["Categories"]) - 1),
                "Iterations": pd.to_numeric(summary.get("Iterations"), errors="coerce"),
                "GradientNorm": pd.to_numeric(summary.get("GradientNorm"), errors="coerce"),
                "LogLik": pd.to_numeric(summary.get("LogLik"), errors="coerce"),
            })
            steps["IncludedInComparison"] = bool(run["AnalysisEligible"])
            step_parts.append(steps)
        except Exception as exc:  # complete attempt ledger
            run["FailureReason"] = f"{type(exc).__name__}: {exc}"
        finally:
            run["Warnings"] = " | ".join(caught_messages)
            run["ElapsedSeconds"] = time.perf_counter() - started
            run_rows.append(run)
        print(
            f"[{len(run_rows):03d}/{len(manifest):03d}] {run_id}: "
            f"returned={run['FitReturned']} converged={run['Converged']} "
            f"eligible={run['AnalysisEligible']}",
            flush=True,
        )

    pd.DataFrame(run_rows).to_csv(output_dir / "python_jmle_step_runs.csv", index=False)
    steps = pd.concat(step_parts, ignore_index=True) if step_parts else pd.DataFrame()
    steps.to_csv(output_dir / "python_jmle_steps.csv", index=False)
    identity = {
        "schema_version": SCHEMA_VERSION,
        "adapter_sha256": sha256_file(Path(__file__).resolve()),
        "app_sha256": sha256_file(REPO_ROOT / "streamlit_app.py"),
        "python_version": sys.version,
        "platform": platform.platform(),
        "selected_runs": int(len(manifest)),
        "fit_controls": build_jmle_kwargs(4),
        "registered_generated_file_sha256": registered_hashes,
        "interpretation_contract": (
            "RSM step thresholds are directly compared only for structurally identified, "
            "matched JMLE runs; FACETS Table 8 display intervals remain explicit."
        ),
    }
    (output_dir / "python_jmle_step_adapter_identity.json").write_text(
        json.dumps(identity, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-id", action="append")
    parser.add_argument("--limit", type=int)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args(argv)


def main() -> None:
    run_adapter(parse_args())


if __name__ == "__main__":
    main()
