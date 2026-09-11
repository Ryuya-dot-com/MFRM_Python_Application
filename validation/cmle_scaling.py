#!/usr/bin/env python3
"""Deterministic scaling envelope for the repository-only exact CMLE core.

Each benchmark case runs in a fresh Python process so peak resident memory is
case-local.  The matrix varies Person count, virtual response-unit count,
PCM parameter count, and unique planned-missingness patterns.  Two additional
audit-only cases exercise the default work and memory guards.

This is an implementation benchmark, not a sample-size recommendation or a
claim about inferential operating characteristics.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from itertools import combinations
import json
from pathlib import Path
import platform
import resource
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
    DEFAULT_RANK_AUDIT_MAX_BYTES,
    cmle_objective_value_grad,
    cmle_objective_value_grad_hessian,
    fit_cmle,
    prepare_cmle_design,
)


SCHEMA_VERSION = "native-exact-cmle-scaling-v1"
DEFAULT_OUTPUT = Path(__file__).resolve().parent / "cmle_scaling_20260809"
BENCHMARK_MAX_WORK = 1_000_000_000
BENCHMARK_MAX_BYTES = 2 * 1024**3
DEFAULT_MAX_WORK = 50_000_000


@dataclass(frozen=True)
class ScalingCase:
    case_id: str
    axis: str
    model: str
    persons: int
    raters: int
    criteria: int
    categories: int = 4
    patterns: int = 1
    mode: str = "fit"
    seed: int = 0


def scaling_cases() -> list[ScalingCase]:
    return [
        ScalingCase("persons-0100", "persons", "RSM", 100, 4, 4, seed=101),
        ScalingCase("persons-1000", "persons", "RSM", 1_000, 4, 4, seed=102),
        ScalingCase("persons-5000", "persons", "RSM", 5_000, 4, 4, seed=103),
        ScalingCase("units-0008", "units", "RSM", 400, 2, 4, seed=201),
        ScalingCase("units-0016", "units", "RSM", 400, 4, 4, seed=202),
        ScalingCase("units-0036", "units", "RSM", 400, 6, 6, seed=203),
        ScalingCase("units-0064", "units", "RSM", 400, 8, 8, seed=204),
        ScalingCase("pcm-params-0012", "parameters", "PCM", 400, 2, 4, seed=301),
        ScalingCase("pcm-params-0014", "parameters", "PCM", 400, 4, 4, seed=302),
        ScalingCase("pcm-params-0022", "parameters", "PCM", 400, 6, 6, seed=303),
        ScalingCase("pcm-params-0030", "parameters", "PCM", 400, 8, 8, seed=304),
        ScalingCase(
            "patterns-0010", "patterns", "PCM", 600, 4, 4, patterns=10, seed=401
        ),
        ScalingCase(
            "patterns-0050", "patterns", "PCM", 600, 4, 4, patterns=50, seed=402
        ),
        ScalingCase(
            "patterns-0100", "patterns", "PCM", 600, 4, 4, patterns=100, seed=403
        ),
        ScalingCase(
            "guard-work",
            "guard",
            "PCM",
            30,
            10,
            10,
            categories=5,
            mode="audit_only",
            seed=501,
        ),
        ScalingCase(
            "guard-memory",
            "guard",
            "PCM",
            20,
            4,
            40,
            categories=5,
            mode="audit_only",
            seed=502,
        ),
    ]


def _retained_unit_patterns(n_units: int, requested: int) -> list[tuple[int, ...]]:
    if requested < 1:
        raise ValueError("patterns must be positive")
    all_units = tuple(range(n_units))
    retained = [all_units]
    if requested == 1:
        return retained
    for omitted_count in range(1, n_units - 1):
        for omitted in combinations(range(n_units), omitted_count):
            omitted_set = set(omitted)
            retained.append(tuple(unit for unit in all_units if unit not in omitted_set))
            if len(retained) == requested:
                return retained
    raise ValueError(f"Cannot construct {requested} patterns from {n_units} units")


def generate_frame(case: ScalingCase) -> pd.DataFrame:
    rng = np.random.default_rng(case.seed)
    units = [
        (f"R{r + 1:02d}", f"C{c + 1:02d}")
        for r in range(case.raters)
        for c in range(case.criteria)
    ]
    patterns = _retained_unit_patterns(len(units), case.patterns)
    rows: list[tuple[str, str, str, int]] = []
    width = max(4, len(str(case.persons)))
    for person_index in range(case.persons):
        retained = patterns[person_index % len(patterns)]
        scores = rng.integers(0, case.categories, size=len(retained))
        if np.all(scores == 0) or np.all(scores == case.categories - 1):
            scores[0] = 1
        person = f"P{person_index + 1:0{width}d}"
        for unit_index, score in zip(retained, scores):
            rater, criterion = units[unit_index]
            rows.append((person, rater, criterion, int(score)))
    return pd.DataFrame(rows, columns=["Person", "Rater", "Criterion", "Score"])


def _common(case: ScalingCase, *, benchmark_caps: bool) -> dict[str, object]:
    return {
        "person_col": "Person",
        "facet_cols": ["Rater", "Criterion"],
        "score_col": "Score",
        "rating_min": 0,
        "rating_max": case.categories - 1,
        "model": case.model,
        "step_facet": "Criterion" if case.model == "PCM" else None,
        "rank_audit_max_work": (
            BENCHMARK_MAX_WORK if benchmark_caps else DEFAULT_MAX_WORK
        ),
        "rank_audit_max_bytes": (
            BENCHMARK_MAX_BYTES
            if benchmark_caps
            else DEFAULT_RANK_AUDIT_MAX_BYTES
        ),
    }


def _peak_rss_mib() -> float:
    raw = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    if sys.platform == "darwin":
        return raw / 1024**2
    return raw / 1024.0


def run_case(case: ScalingCase) -> dict[str, object]:
    total_started = time.perf_counter()
    generated_started = time.perf_counter()
    frame = generate_frame(case)
    data_seconds = time.perf_counter() - generated_started

    benchmark_caps = case.mode == "fit"
    prepared_started = time.perf_counter()
    design = prepare_cmle_design(frame, **_common(case, benchmark_caps=benchmark_caps))
    prepare_seconds = time.perf_counter() - prepared_started
    audit = design.audit
    blocking_codes = [
        row["Code"] for row in audit["issues"] if row["Severity"] == "Block"
    ]
    row: dict[str, object] = {
        **asdict(case),
        "rows": int(len(frame)),
        "units": int(case.raters * case.criteria),
        "patterns_actual": int(audit["unique_missingness_design_patterns"]),
        "parameters": int(design.n_parameters),
        "informative_persons": int(audit["persons_informative"]),
        "eligible": bool(audit["eligible"]),
        "blocking_codes": ";".join(blocking_codes),
        "rank_audit_work_proxy": int(audit["rank_audit_work_proxy"]),
        "rank_audit_peak_bytes_proxy": int(audit["rank_audit_peak_bytes_proxy"]),
        "rank_audit_peak_mib_proxy": float(
            audit["rank_audit_peak_bytes_proxy"] / 1024**2
        ),
        "default_work_cap_exceeded": bool(
            audit["rank_audit_work_proxy"] > DEFAULT_MAX_WORK
        ),
        "default_memory_cap_exceeded": bool(
            audit["rank_audit_peak_bytes_proxy"] > DEFAULT_RANK_AUDIT_MAX_BYTES
        ),
        "data_seconds": float(data_seconds),
        "prepare_seconds": float(prepare_seconds),
        "objective_median_seconds": np.nan,
        "information_seconds": np.nan,
        "fit_seconds": np.nan,
        "converged": False,
        "inference_ready": False,
        "iterations": 0,
        "function_evaluations": 0,
        "optimizer_success": False,
        "optimizer_gradient_sup_norm": np.nan,
        "newton_polish_steps": 0,
        "newton_polish_accepted_steps": 0,
        "newton_polish_reason": "not_run",
        "gradient_sup_norm": np.nan,
        "conditional_loglik": np.nan,
    }

    if case.mode == "fit":
        if not audit["eligible"]:
            raise RuntimeError(
                f"Benchmark case {case.case_id} was unexpectedly blocked: {blocking_codes}"
            )
        zero = np.zeros(design.n_parameters, dtype=float)
        objective_times = []
        for _ in range(5):
            started = time.perf_counter()
            cmle_objective_value_grad(zero, design)
            objective_times.append(time.perf_counter() - started)
        information_started = time.perf_counter()
        cmle_objective_value_grad_hessian(zero, design)
        information_seconds = time.perf_counter() - information_started
        fit_started = time.perf_counter()
        fit = fit_cmle(
            frame,
            maxiter=1_000,
            gtol=1e-8,
            **_common(case, benchmark_caps=True),
        )
        fit_seconds = time.perf_counter() - fit_started
        summary = fit["summary"].iloc[0]
        row.update(
            {
                "objective_median_seconds": float(np.median(objective_times)),
                "information_seconds": float(information_seconds),
                "fit_seconds": float(fit_seconds),
                "converged": bool(summary["Converged"]),
                "inference_ready": bool(summary["InferenceReady"]),
                "iterations": int(summary["Iterations"]),
                "function_evaluations": int(summary["FunctionEvaluations"]),
                "optimizer_success": bool(summary["OptimizerSuccess"]),
                "optimizer_gradient_sup_norm": float(
                    summary["OptimizerGradientSupNorm"]
                ),
                "newton_polish_steps": int(summary["NewtonPolishSteps"]),
                "newton_polish_accepted_steps": int(
                    summary["NewtonPolishAcceptedSteps"]
                ),
                "newton_polish_reason": str(summary["NewtonPolishReason"]),
                "gradient_sup_norm": float(summary["GradientSupNorm"]),
                "conditional_loglik": float(summary["ConditionalLogLik"]),
            }
        )

    row["peak_rss_mib"] = _peak_rss_mib()
    row["total_seconds"] = float(time.perf_counter() - total_started)
    return row


def _git_value(*args: str) -> str:
    completed = subprocess.run(
        ["git", *args], capture_output=True, text=True, check=False
    )
    return completed.stdout.strip() if completed.returncode == 0 else "unavailable"


def run_all(output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    script = Path(__file__).resolve()
    rows = []
    for case in scaling_cases():
        completed = subprocess.run(
            [sys.executable, str(script), "--case", case.case_id],
            capture_output=True,
            text=True,
            check=False,
        )
        if completed.returncode != 0:
            raise RuntimeError(
                f"Scaling case {case.case_id} failed:\n{completed.stderr}"
            )
        rows.append(json.loads(completed.stdout))
    results = pd.DataFrame(rows)
    results.to_csv(output_dir / "scaling_runs.csv", index=False)
    fit_results = results.loc[results["mode"] == "fit"].copy()
    axis_summary = (
        fit_results.groupby("axis", dropna=False)
        .agg(
            Cases=("case_id", "size"),
            InferenceReadyCases=("inference_ready", "sum"),
            MaximumRows=("rows", "max"),
            MaximumUnits=("units", "max"),
            MaximumPatterns=("patterns_actual", "max"),
            MaximumParameters=("parameters", "max"),
            MaximumInformationSeconds=("information_seconds", "max"),
            MaximumFitSeconds=("fit_seconds", "max"),
            MaximumPeakMiBProxy=("rank_audit_peak_mib_proxy", "max"),
            MaximumPeakRSSMiB=("peak_rss_mib", "max"),
        )
        .reset_index()
    )
    axis_summary.to_csv(output_dir / "scaling_axis_summary.csv", index=False)

    identity = {
        "schema_version": SCHEMA_VERSION,
        "run_date": "2026-08-09",
        "git_head": _git_value("rev-parse", "HEAD"),
        "git_status_short": _git_value("status", "--short"),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "cases": int(len(results)),
        "fit_cases": int((results["mode"] == "fit").sum()),
        "audit_only_cases": int((results["mode"] == "audit_only").sum()),
        "inference_ready_fit_cases": int(fit_results["inference_ready"].sum()),
        "default_rank_audit_max_work": DEFAULT_MAX_WORK,
        "default_rank_audit_max_bytes": DEFAULT_RANK_AUDIT_MAX_BYTES,
        "benchmark_rank_audit_max_work": BENCHMARK_MAX_WORK,
        "benchmark_rank_audit_max_bytes": BENCHMARK_MAX_BYTES,
        "timing_scope": "single-machine implementation benchmark",
    }
    (output_dir / "runtime_identity.json").write_text(
        json.dumps(identity, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--case", choices=[case.case_id for case in scaling_cases()])
    args = parser.parse_args()
    if args.case:
        selected = next(case for case in scaling_cases() if case.case_id == args.case)
        print(json.dumps(run_case(selected), allow_nan=True))
    else:
        run_all(args.output_dir)


if __name__ == "__main__":
    main()
