#!/usr/bin/env python3
"""Run a fail-closed FACETS 4.5 versus Python JMLE PCM syntax probe.

The probe consumes one retained additive-RSM pilot run but fits it as a PCM
with Criterion as the step facet.  It qualifies FACETS ``#`` syntax, Table 8
multi-scale parsing, parameter signs/constraints, and exact count mapping.  It
is not a known-truth PCM or operating-characteristics study.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import platform
import re
import sys
import time
from typing import Any, Iterable

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from validation.operating_characteristics_facets import (  # noqa: E402
    FACET_COLUMNS,
    RECOVERY_FACETS,
    align_recovery,
    build_facets_spec,
    invoke_facets,
    parse_iteration_report,
    parse_score_file,
    parse_table8_categories,
    sha256_file,
    validate_bundle,
)
from validation.operating_characteristics_python_table8 import (  # noqa: E402
    build_jmle_kwargs,
)


SCHEMA_VERSION = "mfrm-facets-pcm-syntax-probe-v1"
DEFAULT_RUN_ID = "balanced_small__bias_0p0::rep-00001"
PRIMARY_DECIMALS = 6
AUXILIARY_DECIMALS = 2


def build_pcm_spec(
    manifest_row: pd.Series,
    ratings: pd.DataFrame,
    anchors: pd.DataFrame,
    *,
    score_base: Path,
) -> tuple[str, dict[str, dict[str, int]]]:
    """Convert the registered four-facet RSM specification to Criterion PCM."""

    base, level_maps = build_facets_spec(
        manifest_row,
        ratings,
        anchors,
        score_base=score_base,
    )
    source = "Models=\n?,?,?,?,R3\n?,?B,?B,?,R3\n*"
    replacement = "Models=\n?,?,?,#,R3\n*"
    if base.count(source) != 1:
        raise ValueError("Registered RSM Models block changed; PCM conversion fails closed")
    output = base.replace(source, replacement)
    output = output.replace("Title=FACETS replay ", "Title=FACETS PCM syntax probe ", 1)
    return output, level_maps


def map_pcm_scale_tables(
    categories: pd.DataFrame,
    *,
    criterion_map: dict[str, int],
) -> pd.DataFrame:
    """Map every FACETS PCM Table 8 scale to one Criterion level."""

    if categories.empty:
        raise ValueError("No FACETS PCM Table 8 category rows")
    reverse = {int(number): str(level) for level, number in criterion_map.items()}
    if len(reverse) != len(criterion_map):
        raise ValueError("Criterion element numbers are not one-to-one")
    output = categories.copy()

    def element_number(model: object) -> int:
        fields = [field.strip() for field in str(model).split(",")]
        if len(fields) < 5:
            raise ValueError(f"Unexpected FACETS PCM model heading: {model!r}")
        token = fields[3]
        match = re.fullmatch(r"(?P<number>\d+)#?", token)
        if not match:
            raise ValueError(f"Could not map FACETS PCM scale token: {token!r}")
        return int(match.group("number"))

    table_models = output[["TableNumber", "Model"]].drop_duplicates()
    if table_models["TableNumber"].duplicated().any():
        raise ValueError("One FACETS Table 8 number has multiple model headings")
    table_models["CriterionElementNumber"] = table_models["Model"].map(element_number)
    table_models["StepFacetLevel"] = table_models["CriterionElementNumber"].map(reverse)
    if table_models["StepFacetLevel"].isna().any():
        raise ValueError("FACETS PCM Table 8 references an unknown Criterion element")
    observed = set(table_models["StepFacetLevel"])
    expected = set(criterion_map)
    if observed != expected or table_models["StepFacetLevel"].duplicated().any():
        raise ValueError(
            f"FACETS PCM scale mapping is not one-to-one: observed={observed}, expected={expected}"
        )
    return output.merge(table_models, on=["TableNumber", "Model"], validate="many_to_one")


def normalize_python_pcm_steps(result: dict[str, Any]) -> pd.DataFrame:
    source = result.get("steps", pd.DataFrame()).copy()
    required = {"StepFacet", "Step", "Estimate"}
    missing = required.difference(source.columns)
    if missing:
        raise ValueError(f"Python PCM steps are missing columns: {sorted(missing)}")
    output = source.rename(columns={"StepFacet": "StepFacetLevel"})
    output["StepFacetLevel"] = output["StepFacetLevel"].astype(str)
    output["Category"] = (
        output["Step"].astype(str).str.extract(r"(\d+)$", expand=False).astype(int)
    )
    output["PythonThreshold"] = pd.to_numeric(output["Estimate"], errors="coerce")
    if not np.isfinite(output["PythonThreshold"]).all():
        raise ValueError("Python PCM thresholds include non-finite values")
    return output[["StepFacetLevel", "Category", "Step", "PythonThreshold"]]


def pcm_threshold_pairs(
    primary_categories: pd.DataFrame,
    python_steps: pd.DataFrame,
) -> pd.DataFrame:
    facets = primary_categories.loc[
        pd.to_numeric(primary_categories["Category"], errors="coerce").gt(0),
        [
            "TableNumber",
            "Model",
            "CriterionElementNumber",
            "StepFacetLevel",
            "Category",
            "ThresholdMeasureToken",
            "ThresholdMeasureDisplayed",
            "ThresholdMeasureDisplayDecimals",
            "ThresholdMeasureRawLowerBound",
            "ThresholdMeasureRawUpperBound",
        ],
    ].copy()
    pairs = facets.merge(
        python_steps,
        on=["StepFacetLevel", "Category"],
        how="outer",
        validate="one_to_one",
        indicator=True,
    )
    if not pairs["_merge"].eq("both").all():
        raise ValueError("FACETS and Python PCM threshold keys do not match")
    pairs = pairs.drop(columns="_merge")
    pairs["ThresholdDifference"] = (
        pd.to_numeric(pairs["ThresholdMeasureDisplayed"], errors="coerce")
        - pd.to_numeric(pairs["PythonThreshold"], errors="coerce")
    )
    pairs["AbsoluteThresholdDifference"] = pairs["ThresholdDifference"].abs()
    pairs["PythonWithinFACETSDisplayInterval"] = (
        pairs["PythonThreshold"].ge(pairs["ThresholdMeasureRawLowerBound"])
        & pairs["PythonThreshold"].le(pairs["ThresholdMeasureRawUpperBound"])
    )
    return pairs


def pcm_category_count_audit(
    categories: pd.DataFrame,
    ratings: pd.DataFrame,
) -> pd.DataFrame:
    expected = (
        ratings.groupby(["Criterion", "Score"], observed=False)
        .size()
        .rename("GeneratedCount")
        .reset_index()
        .rename(columns={"Criterion": "StepFacetLevel", "Score": "Category"})
    )
    audit = categories.merge(
        expected,
        on=["StepFacetLevel", "Category"],
        how="outer",
        validate="one_to_one",
        indicator=True,
    )
    audit["CountMatches"] = (
        audit["_merge"].eq("both")
        & pd.to_numeric(audit["TotalCount"], errors="coerce")
        .eq(pd.to_numeric(audit["GeneratedCount"], errors="coerce"))
    )
    return audit


def python_pcm_fit(
    ratings: pd.DataFrame,
    anchors: pd.DataFrame,
    *,
    categories: int,
) -> dict[str, Any]:
    import streamlit_app as app  # pylint: disable=import-outside-toplevel

    kwargs = build_jmle_kwargs(categories)
    kwargs.update({"model": "PCM", "step_facet": "Criterion", "maxit": 400})
    return app.mfrm_estimate(
        data=ratings[["Person", "Rater", "Task", "Criterion", "Score"]].copy(),
        anchor_df=(anchors[["Facet", "Level", "Anchor"]].copy() if not anchors.empty else None),
        **kwargs,
    )


def run_probe(args: argparse.Namespace) -> None:
    input_dir = args.input_dir.resolve()
    output_dir = args.output_dir.resolve()
    facets_exe = args.facets_exe.resolve()
    if output_dir.exists():
        raise FileExistsError(f"Output directory already exists: {output_dir}")
    if not facets_exe.is_file():
        raise FileNotFoundError(f"FACETS executable is missing: {facets_exe}")
    tables = validate_bundle(input_dir)
    manifest = tables["manifest.csv"]
    selected = manifest[manifest["RunId"].astype(str).eq(args.run_id)]
    if len(selected) != 1:
        raise ValueError(f"Expected exactly one manifest row for {args.run_id!r}, found {len(selected)}")
    manifest_row = selected.iloc[0]
    run_id = str(manifest_row["RunId"])
    ratings = tables["generated_ratings.csv"]
    ratings = ratings[ratings["RunId"].astype(str).eq(run_id)].copy()
    truth = tables["generated_facet_truth.csv"]
    truth = truth[truth["RunId"].astype(str).eq(run_id)].copy()
    anchors = tables["generated_anchors.csv"]
    anchors = anchors[anchors["RunId"].astype(str).eq(run_id)].copy()
    if ratings.empty:
        raise ValueError("Selected PCM probe run has no rating rows")

    output_dir.mkdir(parents=True)
    spec_path = output_dir / "pcm_analysis.txt"
    report_path = output_dir / "pcm_report_u6.txt"
    aux_report_path = output_dir / "pcm_report_u2.txt"
    score_base = output_dir / "pcm_scores.txt"
    aux_score_base = output_dir / "pcm_scores_u2.txt"
    spec, level_maps = build_pcm_spec(
        manifest_row,
        ratings,
        anchors,
        score_base=score_base,
    )
    spec_path.write_text(spec, encoding="utf-8", newline="\n")

    started = time.perf_counter()
    primary = invoke_facets(
        facets_exe,
        spec_path,
        report_path,
        timeout_seconds=args.timeout_seconds,
    )
    (output_dir / "pcm_stdout_u6.txt").write_text(primary.stdout or "", encoding="utf-8")
    (output_dir / "pcm_stderr_u6.txt").write_text(primary.stderr or "", encoding="utf-8")
    if primary.returncode != 0 or not report_path.is_file():
        raise RuntimeError(f"FACETS PCM primary pass failed: exit={primary.returncode}")
    aux = invoke_facets(
        facets_exe,
        spec_path,
        aux_report_path,
        timeout_seconds=args.timeout_seconds,
        # Keep the auxiliary two-decimal score files separate: the primary
        # Umean=6 pass is authoritative for direct parameter comparisons.
        extra_specs=(
            f"Umean=0,1,{AUXILIARY_DECIMALS}",
            f"Scorefile={aux_score_base}",
        ),
    )
    (output_dir / "pcm_stdout_u2.txt").write_text(aux.stdout or "", encoding="utf-8")
    (output_dir / "pcm_stderr_u2.txt").write_text(aux.stderr or "", encoding="utf-8")
    if aux.returncode != 0 or not aux_report_path.is_file():
        raise RuntimeError(f"FACETS PCM auxiliary pass failed: exit={aux.returncode}")

    report_info = parse_iteration_report(report_path)
    score_parts = [
        parse_score_file(
            output_dir / f"pcm_scores.{facet_number}.txt",
            facet_number=facet_number,
            facet_name=facet_name,
        )
        for facet_number, facet_name in FACET_COLUMNS
    ]
    scores = pd.concat(score_parts, ignore_index=True)
    facets_recovery = align_recovery(scores, truth, anchors, manifest_row)

    primary_categories = map_pcm_scale_tables(
        parse_table8_categories(report_path),
        criterion_map=level_maps["Criterion"],
    )
    auxiliary_categories = map_pcm_scale_tables(
        parse_table8_categories(aux_report_path),
        criterion_map=level_maps["Criterion"],
    )
    count_audit = pcm_category_count_audit(auxiliary_categories, ratings)
    if not count_audit["CountMatches"].all():
        raise ValueError("FACETS PCM Table 8 category counts do not match retained data")

    python_result = python_pcm_fit(
        ratings,
        anchors,
        categories=int(manifest_row["Categories"]),
    )
    python_summary = python_result["summary"].iloc[0]
    python_estimates = python_result["facets"]["others"].copy()
    python_estimates["SE"] = np.nan
    python_estimates["Status"] = np.nan
    python_recovery = align_recovery(python_estimates, truth, anchors, manifest_row)
    recovery_keys = ["RunId", "Facet", "Level"]
    recovery_pairs = facets_recovery.merge(
        python_recovery[
            recovery_keys + ["EstimateAligned", "ErrorAligned", "ComparisonScale"]
        ].rename(columns={
            "EstimateAligned": "PythonEstimateAligned",
            "ErrorAligned": "PythonErrorAligned",
            "ComparisonScale": "PythonComparisonScale",
        }),
        on=recovery_keys,
        how="inner",
        validate="one_to_one",
    )
    recovery_pairs["EstimateDifference"] = (
        recovery_pairs["EstimateAligned"] - recovery_pairs["PythonEstimateAligned"]
    )
    recovery_pairs["AbsoluteEstimateDifference"] = recovery_pairs["EstimateDifference"].abs()
    recovery_pairs["ComparisonEligible"] = (
        recovery_pairs["Facet"].isin(RECOVERY_FACETS)
        & bool(report_info.get("Converged", False))
        & bool(python_summary.get("Converged", False))
        & bool(python_summary.get("InferenceReady", False))
    )

    python_steps = normalize_python_pcm_steps(python_result)
    threshold_pairs = pcm_threshold_pairs(primary_categories, python_steps)
    threshold_pairs["ComparisonEligible"] = (
        bool(report_info.get("Converged", False))
        & bool(python_summary.get("Converged", False))
        & bool(python_summary.get("InferenceReady", False))
    )

    primary_categories.to_csv(output_dir / "facets_pcm_table8_u6.csv", index=False)
    auxiliary_categories.to_csv(output_dir / "facets_pcm_table8_u2.csv", index=False)
    count_audit.to_csv(output_dir / "facets_pcm_category_count_audit.csv", index=False)
    facets_recovery.to_csv(output_dir / "facets_pcm_recovery.csv", index=False)
    python_recovery.to_csv(output_dir / "python_pcm_recovery.csv", index=False)
    recovery_pairs.to_csv(output_dir / "facets_python_pcm_parameter_pairs.csv", index=False)
    python_steps.to_csv(output_dir / "python_pcm_steps.csv", index=False)
    threshold_pairs.to_csv(output_dir / "facets_python_pcm_threshold_pairs.csv", index=False)

    direct_parameters = recovery_pairs[recovery_pairs["ComparisonEligible"]]
    direct_thresholds = threshold_pairs[threshold_pairs["ComparisonEligible"]]
    metrics = {
        "schema_version": SCHEMA_VERSION,
        "run_id": run_id,
        "data_generating_model": "RSM_negative_control_fitted_as_PCM",
        "facets_version": report_info.get("FacetsVersion"),
        "facets_returned": True,
        "facets_converged": bool(report_info.get("Converged", False)),
        "python_returned": True,
        "python_converged": bool(python_summary.get("Converged", False)),
        "python_inference_ready": bool(python_summary.get("InferenceReady", False)),
        "table8_scale_tables": int(primary_categories["TableNumber"].nunique()),
        "category_rows": int(len(auxiliary_categories)),
        "category_count_matches": int(count_audit["CountMatches"].sum()),
        "direct_parameters": int(len(direct_parameters)),
        "parameter_mae_logits": (
            float(direct_parameters["AbsoluteEstimateDifference"].mean())
            if len(direct_parameters) else None
        ),
        "parameter_max_abs_difference_logits": (
            float(direct_parameters["AbsoluteEstimateDifference"].max())
            if len(direct_parameters) else None
        ),
        "direct_thresholds": int(len(direct_thresholds)),
        "threshold_mae_logits": (
            float(direct_thresholds["AbsoluteThresholdDifference"].mean())
            if len(direct_thresholds) else None
        ),
        "threshold_max_abs_difference_logits": (
            float(direct_thresholds["AbsoluteThresholdDifference"].max())
            if len(direct_thresholds) else None
        ),
        "threshold_sums_facets": {
            str(level): float(group["ThresholdMeasureDisplayed"].sum())
            for level, group in direct_thresholds.groupby("StepFacetLevel")
        },
        "threshold_sums_python": {
            str(level): float(group["PythonThreshold"].sum())
            for level, group in direct_thresholds.groupby("StepFacetLevel")
        },
        "elapsed_seconds": time.perf_counter() - started,
        "claim_limit": "Syntax/parser and RSM-negative-control PCM implementation probe only.",
    }
    (output_dir / "facets_pcm_probe_metrics.json").write_text(
        json.dumps(metrics, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    identity = {
        "schema_version": SCHEMA_VERSION,
        "adapter_sha256": sha256_file(Path(__file__).resolve()),
        "app_sha256": sha256_file(REPO_ROOT / "streamlit_app.py"),
        "facets_executable": str(facets_exe),
        "facets_executable_sha256": sha256_file(facets_exe),
        "python_version": sys.version,
        "platform": platform.platform(),
        "input_sha256": {
            name: sha256_file(input_dir / name)
            for name in (
                "manifest.csv",
                "generated_ratings.csv",
                "generated_facet_truth.csv",
                "generated_anchors.csv",
            )
        },
        "specification_sha256": sha256_file(spec_path),
        "primary_report_sha256": sha256_file(report_path),
        "auxiliary_report_sha256": sha256_file(aux_report_path),
        "model_contract": "FACETS Models=?,?,?,#,R3; Python PCM step_facet=Criterion",
    }
    (output_dir / "facets_pcm_probe_identity.json").write_text(
        json.dumps(identity, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-id", default=DEFAULT_RUN_ID)
    parser.add_argument("--facets-exe", type=Path, default=Path(r"C:\Facets\Facets.exe"))
    parser.add_argument("--timeout-seconds", type=float, default=120.0)
    return parser.parse_args(argv)


def main() -> None:
    run_probe(parse_args())


if __name__ == "__main__":
    main()
