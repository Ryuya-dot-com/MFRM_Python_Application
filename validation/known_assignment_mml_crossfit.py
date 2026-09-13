#!/usr/bin/env python3
"""Cross-evaluate the independent base-R free-SD PCM MML refits."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
from typing import Any, Iterable

import numpy as np
from numpy.polynomial.hermite import hermgauss
import pandas as pd
from scipy.special import logsumexp


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from validation.operating_characteristics_facets import sha256_file  # noqa: E402


SCHEMA_VERSION = "known_assignment_mml_crossfit_v1"
R_SCRIPT = ROOT / "validation" / "known_assignment_mml_crossfit.R"
PARAMETER_BLOCK_LEVELS = {
    "Rater": ("R01", "R02", "R03", "R04"),
    "Task": ("T01", "T02", "T03"),
    "Criterion": ("C01", "C02"),
    "Step": ("C01::1", "C01::2", "C01::3", "C02::1", "C02::2", "C02::3"),
}


def _json_dump(path: Path, value: Any) -> None:
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def quadrature(points: int, sigma: float) -> tuple[np.ndarray, np.ndarray]:
    nodes, weights = hermgauss(int(points))
    return np.sqrt(2.0) * nodes * float(sigma), weights / np.sqrt(np.pi)


def parameter_surface(frame: pd.DataFrame, column: str) -> dict[str, Any]:
    output: dict[str, Any] = {}
    for block, levels in PARAMETER_BLOCK_LEVELS.items():
        selected = frame.loc[frame["Block"].astype(str).eq(block)].set_index("Level")
        if set(selected.index.astype(str)) != set(levels):
            raise ValueError(f"Cross-fit parameter levels changed for {block}")
        values = pd.to_numeric(selected[column], errors="raise")
        if block == "Step":
            output[block] = {
                criterion: np.array(
                    [float(values.loc[f"{criterion}::{step}"]) for step in (1, 2, 3)],
                    dtype=float,
                )
                for criterion in ("C01", "C02")
            }
        else:
            output[block] = {
                level: float(values.loc[level]) for level in levels
            }
    return output


def marginal_loglik(
    ratings: pd.DataFrame,
    surface: dict[str, Any],
    *,
    sigma: float,
    points: int,
) -> float:
    """Protocol-level Python evaluator independent of the app optimizer."""

    if not np.isfinite(sigma) or sigma <= 0:
        raise ValueError("Population SD must be positive and finite")
    data = ratings.copy()
    score = pd.to_numeric(data["Score"], errors="raise").astype(int).to_numpy()
    if not np.isin(score, (0, 1, 2, 3)).all():
        raise ValueError("Cross-fit score map must remain 0..3")
    base = -(
        data["Rater"].astype(str).map(surface["Rater"]).to_numpy(dtype=float)
        + data["Task"].astype(str).map(surface["Task"]).to_numpy(dtype=float)
        + data["Criterion"].astype(str).map(surface["Criterion"]).to_numpy(dtype=float)
    )
    if not np.isfinite(base).all():
        raise ValueError("Cross-fit facet map contains missing levels")
    criteria = data["Criterion"].astype(str).to_numpy()
    cumulative = {
        criterion: np.concatenate([[0.0], np.cumsum(surface["Step"][criterion])])
        for criterion in ("C01", "C02")
    }
    thresholds = np.vstack([cumulative[criterion] for criterion in criteria])
    person_codes, people = pd.factorize(data["Person"].astype(str), sort=True)
    nodes, weights = quadrature(points, sigma)
    person_node = np.zeros((len(people), points), dtype=float)
    categories = np.arange(4, dtype=float)
    row_index = np.arange(len(data))
    for node_index, node in enumerate(nodes):
        logits = np.outer(node + base, categories) - thresholds
        row_ll = logits[row_index, score] - logsumexp(logits, axis=1)
        person_node[:, node_index] = np.bincount(
            person_codes, weights=row_ll, minlength=len(people)
        )
    return float(np.sum(logsumexp(np.log(weights)[None, :] + person_node, axis=1)))


def _maximum_absolute(frame: pd.DataFrame, column: str) -> float:
    return float(pd.to_numeric(frame[column], errors="coerce").abs().max())


def assess_crossfit(
    *,
    study_dir: Path,
    output_dir: Path,
    expected_datasets: int,
    tolerance: dict[str, float],
) -> dict[str, Any]:
    study_dir = study_dir.resolve()
    output_dir = output_dir.resolve()
    runs = pd.read_csv(output_dir / "crossfit_runs.csv")
    parameters = pd.read_csv(output_dir / "crossfit_parameters.csv")
    gradients = pd.read_csv(output_dir / "crossfit_gradients.csv")
    runtime = pd.read_csv(output_dir / "runtime_identity.csv")
    ratings_all = pd.read_csv(study_dir / "retained_input" / "generated_ratings.csv")
    if len(runs) != expected_datasets or runs["RunId"].nunique() != expected_datasets:
        raise ValueError("R cross-fit dataset denominator changed")
    if len(parameters) != expected_datasets * 15:
        raise ValueError("R cross-fit expanded-parameter denominator changed")
    if len(gradients) != expected_datasets * 12:
        raise ValueError("R cross-fit gradient denominator changed")

    evaluation_rows: list[dict[str, Any]] = []
    for run in runs.itertuples(index=False):
        run_id = str(run.RunId)
        parameter = parameters.loc[parameters["RunId"].astype(str).eq(run_id)]
        ratings = ratings_all.loc[
            ratings_all["RunId"].astype(str).eq(run_id),
            ["Person", "Rater", "Task", "Criterion", "Score"],
        ]
        python_surface = parameter_surface(parameter, "PythonEstimate")
        r31_surface = parameter_surface(parameter, "REstimateQ31")
        r61_surface = parameter_surface(parameter, "REstimateQ61")
        values = {
            "PythonSolutionPythonLogLikQ31": marginal_loglik(
                ratings, python_surface, sigma=float(run.PythonSigma), points=31
            ),
            "PythonSolutionPythonLogLikQ61": marginal_loglik(
                ratings, python_surface, sigma=float(run.PythonSigma), points=61
            ),
            "RQ31SolutionPythonLogLikQ31": marginal_loglik(
                ratings, r31_surface, sigma=float(run.RSigma), points=31
            ),
            "RQ31SolutionPythonLogLikQ61": marginal_loglik(
                ratings, r31_surface, sigma=float(run.RSigma), points=61
            ),
            "RQ61SolutionPythonLogLikQ61": marginal_loglik(
                ratings, r61_surface, sigma=float(run.RQ61Sigma), points=61
            ),
        }
        evaluation_rows.append({"RunId": run_id, **values})
    evaluations = pd.DataFrame(evaluation_rows).merge(
        runs, on="RunId", how="left", validate="one_to_one"
    )
    evaluations["AbsRAtPythonVsRecordedQ31"] = (
        evaluations["PythonSolutionRLogLikQ31"] - evaluations["PythonRecordedLogLik"]
    ).abs()
    evaluations["AbsPythonAtPythonVsRecordedQ31"] = (
        evaluations["PythonSolutionPythonLogLikQ31"] - evaluations["PythonRecordedLogLik"]
    ).abs()
    evaluations["AbsCrossLanguageAtPythonQ31"] = (
        evaluations["PythonSolutionPythonLogLikQ31"]
        - evaluations["PythonSolutionRLogLikQ31"]
    ).abs()
    evaluations["AbsCrossLanguageAtRQ31"] = (
        evaluations["RQ31SolutionPythonLogLikQ31"] - evaluations["RLogLikQ31"]
    ).abs()
    evaluations["AbsCrossLanguageAtRQ61"] = (
        evaluations["RQ61SolutionPythonLogLikQ61"]
        - evaluations["ROptimizedLogLikQ61"]
    ).abs()
    evaluations["RQ31ImprovementOverPythonQ31"] = (
        evaluations["RLogLikQ31"] - evaluations["PythonRecordedLogLik"]
    )
    evaluations["PythonQ61MinusQ31"] = (
        evaluations["PythonSolutionPythonLogLikQ61"]
        - evaluations["PythonSolutionPythonLogLikQ31"]
    )
    evaluations["RQ61OptimizedMinusRQ31EvaluatedQ61"] = (
        evaluations["ROptimizedLogLikQ61"]
        - evaluations["RQ31SolutionPythonLogLikQ61"]
    )

    parameter_q31 = _maximum_absolute(parameters, "DifferenceRQ31MinusPython")
    parameter_q61 = _maximum_absolute(parameters, "DifferenceRQ61MinusRQ31")
    q31_sigma = float((evaluations["RSigma"] - evaluations["PythonSigma"]).abs().max())
    q61_sigma = float((evaluations["RQ61Sigma"] - evaluations["RSigma"]).abs().max())
    gates = {
        "r_version_4_5_1": str(runtime.loc[0, "RVersion"]) == "4.5.1",
        "datasets_complete": len(evaluations) == expected_datasets,
        "r_q31_convergence_all": bool(evaluations["RConvergenceCode"].eq(0).all()),
        "r_q61_convergence_all": bool(evaluations["RQ61ConvergenceCode"].eq(0).all()),
        "r_at_python_matches_recorded_q31": float(
            evaluations["AbsRAtPythonVsRecordedQ31"].max()
        )
        <= tolerance["cross_language_loglik"],
        "python_at_python_matches_recorded_q31": float(
            evaluations["AbsPythonAtPythonVsRecordedQ31"].max()
        )
        <= tolerance["cross_language_loglik"],
        "cross_language_python_solution_q31": float(
            evaluations["AbsCrossLanguageAtPythonQ31"].max()
        )
        <= tolerance["cross_language_loglik"],
        "cross_language_r_solution_q31": float(
            evaluations["AbsCrossLanguageAtRQ31"].max()
        )
        <= tolerance["cross_language_loglik"],
        "cross_language_r_solution_q61": float(
            evaluations["AbsCrossLanguageAtRQ61"].max()
        )
        <= tolerance["cross_language_loglik"],
        "r_q31_terminal_gradient": float(evaluations["RGradientSupNormQ31"].max())
        <= tolerance["r_terminal_gradient_q31"],
        "r_q61_terminal_gradient": float(evaluations["RGradientSupNormQ61"].max())
        <= tolerance["r_terminal_gradient_q61"],
        "python_sigma_fixed_point": float(
            evaluations["PythonSigmaFixedPointResidual"].abs().max()
        )
        <= tolerance["python_sigma_fixed_point"],
        "r_q31_vs_python_parameter": parameter_q31 <= tolerance["r_q31_parameter"],
        "r_q31_vs_python_sigma": q31_sigma <= tolerance["r_q31_sigma"],
        "q61_vs_q31_parameter": parameter_q61 <= tolerance["q61_parameter"],
        "q61_vs_q31_sigma": q61_sigma <= tolerance["q61_sigma"],
        "r_q31_loglik_improvement_bounded": float(
            evaluations["RQ31ImprovementOverPythonQ31"].max()
        )
        <= tolerance["r_q31_loglik_improvement"],
    }
    metrics = {
        "maximum_abs_r_at_python_vs_recorded_q31": float(
            evaluations["AbsRAtPythonVsRecordedQ31"].max()
        ),
        "maximum_abs_python_at_python_vs_recorded_q31": float(
            evaluations["AbsPythonAtPythonVsRecordedQ31"].max()
        ),
        "maximum_cross_language_loglik_difference": float(
            evaluations[
                [
                    "AbsCrossLanguageAtPythonQ31",
                    "AbsCrossLanguageAtRQ31",
                    "AbsCrossLanguageAtRQ61",
                ]
            ].to_numpy(dtype=float).max()
        ),
        "maximum_r_q31_gradient_sup_norm": float(evaluations["RGradientSupNormQ31"].max()),
        "maximum_r_q61_gradient_sup_norm": float(evaluations["RGradientSupNormQ61"].max()),
        "maximum_abs_python_sigma_fixed_point_residual": float(
            evaluations["PythonSigmaFixedPointResidual"].abs().max()
        ),
        "maximum_abs_r_q31_minus_python_parameter": parameter_q31,
        "maximum_abs_r_q31_minus_python_sigma": q31_sigma,
        "maximum_abs_r_q61_minus_q31_parameter": parameter_q61,
        "maximum_abs_r_q61_minus_q31_sigma": q61_sigma,
        "maximum_r_q31_loglik_improvement_over_python": float(
            evaluations["RQ31ImprovementOverPythonQ31"].max()
        ),
        "maximum_abs_python_q61_minus_q31_loglik": float(
            evaluations["PythonQ61MinusQ31"].abs().max()
        ),
    }
    evaluations.to_csv(output_dir / "crossfit_cross_evaluations.csv", index=False, lineterminator="\n")
    identity = {
        "schema_version": f"{SCHEMA_VERSION}_identity_v1",
        "study_dir": str(study_dir),
        "expected_datasets": expected_datasets,
        "r_script_sha256": sha256_file(R_SCRIPT),
        "python_verifier_sha256": sha256_file(Path(__file__).resolve()),
        "source_sha256": {
            filename: sha256_file(output_dir / filename)
            for filename in (
                "crossfit_runs.csv",
                "crossfit_parameters.csv",
                "crossfit_gradients.csv",
                "runtime_identity.csv",
            )
        },
        "retained_input_sha256": {
            filename: sha256_file(study_dir / "retained_input" / filename)
            for filename in ("manifest.csv", "generated_ratings.csv", "attempt_manifest.csv")
        },
    }
    _json_dump(output_dir / "verification_identity.json", identity)
    assessment = {
        "schema_version": SCHEMA_VERSION,
        "pass": bool(all(gates.values())),
        "gates": {key: bool(value) for key, value in gates.items()},
        "tolerance": tolerance,
        "metrics": metrics,
        "datasets": expected_datasets,
        "r_version": str(runtime.loc[0, "RVersion"]),
        "r_script_sha256": identity["r_script_sha256"],
        "python_verifier_sha256": identity["python_verifier_sha256"],
        "claim_boundary": "Same-estimand numerical cross-fit; not FACETS parity, endpoint confirmation, DGP replication, or estimator ranking.",
    }
    _json_dump(output_dir / "assessment.json", assessment)
    return assessment


def default_tolerance() -> dict[str, float]:
    """Candidate tolerances; freeze them in registration after fixture qualification."""

    return {
        "cross_language_loglik": 1e-9,
        "r_terminal_gradient_q31": 1e-3,
        "r_terminal_gradient_q61": 1e-2,
        "python_sigma_fixed_point": 1e-3,
        "r_q31_parameter": 0.01,
        "r_q31_sigma": 0.01,
        "q61_parameter": 0.01,
        "q61_sigma": 0.01,
        "r_q31_loglik_improvement": 0.01,
    }


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected-datasets", type=int, required=True)
    return parser.parse_args(argv)


def main() -> None:
    args = parse_args()
    result = assess_crossfit(
        study_dir=args.study,
        output_dir=args.output,
        expected_datasets=args.expected_datasets,
        tolerance=default_tolerance(),
    )
    print(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
