"""Reproducible development check of joint MML curvature against independent R.

All generated cases are excluded from prospective qualification. The reference
shares the declared quadrature rule, not Python's probability/derivative code.
No coverage, continuous-integral accuracy, or public inference claim is made.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, replace
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import subprocess
import sys
import zipfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd

import streamlit_app as app
from mfrm_app.mml_engine_v2 import run_free_sd_stationarity_v2
from mfrm_app.mml_stationarity import JointPolishOptions
from validation.mml_free_sd_quadrature_adapter import app_free_sd_problem_digest
from validation.mml_free_sd_stationarity_adapter import prepare_app_free_sd_problem_from_data


def save(path, value):
    def encode(item):
        if isinstance(item, np.ndarray):
            return item.tolist()
        if isinstance(item, np.generic):
            return item.item()
        raise TypeError(type(item).__name__)
    path.write_text(json.dumps(value, default=encode, indent=2, allow_nan=False) + "\n")


def generate(seed, model, design, persons):
    """Explicit DGP, independent of the app probability and simulation helpers."""
    rng = np.random.default_rng(seed)
    theta = rng.normal(0.0, 1.2, persons)
    rater, task, criterion = [-.4, 0., .4], [-.2, .2], [-.4, .1, .7]
    steps = np.array([[-.8, -.1, .9]] * 3 if model == "RSM" else
                     [[-.8, -.1, .9], [-1., .2, .8], [-.5, -.2, .7]])
    rows = []
    for p in range(persons):
        for r in range(3):
            for t in range(2):
                for c in range(3):
                    if design == "sparse" and (r == p % 3 or t != p % 2 or c == (p+1) % 3):
                        continue
                    if design == "unbalanced" and (r, t, c) != (p % 3, p % 2, p % 3):
                        if rng.random() > [.25, .6, .95][r] * [.45, .8, 1.][c]:
                            continue
                    logits = (theta[p] - rater[r] - task[t] - criterion[c]) * np.arange(4)
                    logits -= np.r_[0., np.cumsum(steps[c])]
                    probability = np.exp(logits - logits.max()); probability /= probability.sum()
                    rows.append((f"P{p:03d}", f"R{r}", f"T{t}", f"C{c}", rng.choice(4, p=probability)))
    return pd.DataFrame(rows, columns=["Person", "Rater", "Task", "Criterion", "Score"])


def expanded_coordinates(problem, point):
    expanded = app.expand_params(point[:-1], problem.sizes, problem.config)
    return np.r_[np.concatenate([expanded["facets"][f] for f in problem.config["facet_names"]]),
                 expanded["steps"] if problem.config["model"] == "RSM" else expanded["steps_mat"].ravel(),
                 np.exp(point[-1])]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--rscript", default="/usr/local/bin/Rscript")
    args = parser.parse_args()
    out = args.output.resolve(); out.mkdir(parents=True, exist_ok=False)
    options = JointPolishOptions(maxiter=300, gtol=1e-8, ftol=1e-15)
    contract = {
        "classification": "DEVELOPMENT_ONLY", "scientific_inference_ready": False,
        "persons": 48, "models": ["RSM", "PCM"], "designs": ["regular", "sparse", "unbalanced"],
        "seeds": list(range(2026100410, 2026100416)), "quadrature_points": [31, 61],
        "optimizer": asdict(options), "hessian_relative_steps": [1e-5, 5e-6],
        "reference_derivatives": {"coordinates": "point + max(1, abs(point)) * u; u=0",
            "method": "numDeriv Richardson", "hessian_eps": [.01, .005], "gradient_eps": .001,
            "richardson_levels": 4},
        "tolerances": {"nll": 1e-8, "gradient": 2e-5, "information": 1e-4,
            "information_step": 1e-4, "covariance": 1e-5, "expanded": 1e-10,
            "expanded_jacobian": 1e-8, "expanded_covariance": 1e-5},
        "scope": "same-point finite-GH probability, gradient, full joint curvature and coordinate transformation",
        "inversion": "direct inverse only when positive definite; no eigenvalue clipping",
        "not_qualified": ["sampling covariance", "coverage", "continuous integral", "prospective evidence"],
    }
    save(out / "contract.json", contract)  # Before data generation or fitting.
    (out / "reference.R").write_bytes((ROOT / "validation/mml_joint_information_reference.R").read_bytes())
    paths = sorted(set([Path(__file__), ROOT / "streamlit_app.py",
        * (ROOT / "mfrm_app").glob("*.py"),
        ROOT / "validation/mml_free_sd_stationarity_adapter.py",
        ROOT / "validation/mml_free_sd_quadrature_adapter.py",
        ROOT / "validation/mml_joint_information_reference.R"]))
    with zipfile.ZipFile(out / "sources.zip", "w", zipfile.ZIP_DEFLATED) as archive:
        for path in paths:
            archive.write(path, str(path.relative_to(ROOT)))
    save(out / "source_manifest.json", {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    save(out / "runtime.json", {"python": sys.version, "platform": platform.platform(),
        "packages": {d.metadata["Name"]: d.version for d in importlib.metadata.distributions()}})
    cases, prefit = [], []
    for model in contract["models"]:
        for design in contract["designs"]:
            seed = contract["seeds"][len(cases)]
            data = generate(seed, model, design, contract["persons"])
            case_id = f"{model}_{design}"
            path = out / f"{case_id}.csv"; data.to_csv(path, index=False)
            config = dict(person_col="Person", facet_cols=["Rater", "Task", "Criterion"],
                score_col="Score", rating_min=0, rating_max=3, model=model,
                noncenter_facet="Criterion", step_facet="Criterion" if model == "PCM" else None)
            prep, problem = prepare_app_free_sd_problem_from_data(data, **config)
            prefit.append({"id": case_id, "seed": seed, "input_config": config,
                "input_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                "likelihood_problem_digest": app_free_sd_problem_digest({"prep": prep}, problem),
                "structural_start": problem.structural_start, "sigma_start": problem.sigma_start})
            cases.append((case_id, data, problem))
    save(out / "prefit_manifest.json", prefit)
    python_cases, reference_cases, failures = [], [], []
    for case_id, data, base in cases:
        for q in contract["quadrature_points"]:
            run_id = f"{case_id}_Q{q}"
            try:
                problem = replace(base, quadrature_points=q)
                run = run_free_sd_stationarity_v2(
                    problem.structural_start, problem.sigma_start, problem.value, problem.value_gradient,
                    observations=problem.observations, structural_bounds=problem.structural_bounds,
                    sigma_bounds=problem.sigma_bounds, options=options,
                    constraint_residual_function=problem.constraint_residual,
                    joint_value_gradient=problem.joint_value_gradient)
                point = np.asarray(run.restart_polish.joint_coordinates)
                quad = app.gauss_hermite_normal(q, sd=np.exp(point[-1]))
                information = [app._compute_mml_observed_information(point, problem.idx, problem.config,
                    problem.sizes, quad, rel_step=h, include_log_sigma=True)
                    for h in contract["hessian_relative_steps"]]
                positive = np.linalg.eigvalsh(information[1]).min() > 0
                covariance = np.linalg.inv(information[1]) if positive else None
                expanded = expanded_coordinates(problem, point)
                jacobian = np.zeros((len(expanded), len(point)))
                for j in range(len(point)-1):
                    delta = np.eye(len(point))[j]
                    jacobian[:, j] = (expanded_coordinates(problem, point+delta) - expanded_coordinates(problem, point-delta))/2
                jacobian[-1, -1] = np.exp(point[-1])
                record = {"id": run_id, "run": run.to_dict(), "point": point,
                    "nll": problem.joint_value_gradient(point)[0],
                    "gradient": problem.joint_value_gradient(point)[1], "information": information[1],
                    "information_step_disagreement": np.max(np.abs(information[1]-information[0])),
                    "positive_definite": positive, "covariance": covariance,
                    "expanded": expanded, "expanded_jacobian": jacobian,
                    "expanded_covariance": jacobian @ covariance @ jacobian.T if positive else None}
                save(out / f"{run_id}.json", record)
                python_cases.append(record)
                reference_cases.append({"id": run_id, "config": {k: problem.config[k] for k in
                    ("model", "facet_names", "facet_levels", "noncenter_facet", "facet_signs",
                     "step_facet", "n_cat", "rating_min", "rating_max")},
                    "data": data.to_dict("records"), "point": point,
                    "nodes": app.gauss_hermite_normal(q)["nodes"], "weights": quad["weights"]})
                print(run_id, "Python complete", flush=True)
            except Exception as exc:
                failures.append({"id": run_id, "type": type(exc).__name__, "message": str(exc)})
                save(out / "failures.json", failures)
    save(out / "reference_input.json", {"classification": "DEVELOPMENT_ONLY", "cases": reference_cases})
    with (out / "R_execution.log").open("w") as log:
        result = subprocess.run([args.rscript, str(out / "reference.R"),
            str(out / "reference_input.json"), str(out / "R_results.json")], stdout=log, stderr=subprocess.STDOUT)
    differences = []
    if result.returncode == 0:
        reference = json.loads((out / "R_results.json").read_text())
        by_id = {c["id"]: c for c in reference["cases"]}
        for case in python_cases:
            other = by_id[case["id"]]
            errors = {k: float(np.max(np.abs(np.asarray(case[k])-np.asarray(other[k]))))
                for k in ("nll", "gradient", "information", "expanded", "expanded_jacobian")}
            errors["information_step"] = max(case["information_step_disagreement"], other["hessian_step_disagreement"])
            for k in ("covariance", "expanded_covariance"):
                errors[k] = float(np.max(np.abs(np.asarray(case[k])-np.asarray(other[k])))) if case[k] is not None and other[k] is not None else None
            checks = {k: v is not None and v <= contract["tolerances"][k] for k, v in errors.items()}
            differences.append({"id": case["id"], "errors": errors, "checks": checks,
                "all_checks_pass": all(checks.values()),
                "optimizer_terminated": case["run"]["algorithm_terminated"],
                "returned_score_supnorm": float(np.max(np.abs(case["gradient"])))})
    summary = {"classification": "DEVELOPMENT_ONLY", "scientific_inference_ready": False,
        "planned_cases": 12, "completed_python_cases": len(python_cases), "R_returncode": result.returncode,
        "failures": failures, "comparisons": differences,
        "all_checks_pass": len(differences) == 12 and not failures and all(x["all_checks_pass"] for x in differences)}
    save(out / "summary.json", summary)
    save(out / "manifest.json", {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(out.iterdir()) if p.is_file()})
    print(json.dumps(summary, indent=2), flush=True)
    return 0 if summary["all_checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
