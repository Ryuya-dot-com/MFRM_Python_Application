#!/usr/bin/env python3
"""Cross-check the observed 24-Person RSM/PCM fixture; no inference promotion.

Reuses retained Q31/Q61 points, the app's v2 kernel, and the independent R
likelihood. Additional Q121/Q181 fits and continuous evaluations are development
diagnostics. Requires R/statmod/jsonlite and the existing Python test environment.
"""
from dataclasses import replace
import argparse
import hashlib
import json
from pathlib import Path
import platform
import runpy
import subprocess
import sys

import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
import streamlit_app as app
from mfrm_app.mml_engine_v2 import run_free_sd_stationarity_v2
from mfrm_app.mml_stationarity import JointPolishOptions, make_joint_free_sd_functions
from validation.mml_free_sd_stationarity_adapter import prepare_app_free_sd_problem

CLASSIFICATION = "OBSERVED_DEVELOPMENT_ONLY"
ORDERS = (31, 61, 121, 181)
# Implementation checks fixed before this probe's results; not scientific margins.
LIMITS = {"nll": 1e-9, "gradient": 1e-6, "probability": 1e-12,
          "eap": 1e-10, "sd": 1e-10, "node": 1e-12,
          "relative_weight": 1e-10, "moment": 1e-10,
          "refit_coordinate": 1e-4, "refit_nll": 1e-8, "refit_score": 1e-4}


def dump(path, value):
    with path.open("x") as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write("\n")


def difference(left, right):
    a, b = np.asarray(left, dtype=float), np.asarray(right, dtype=float)
    if a.shape != b.shape or not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError("Nonfinite or misaligned cross-language evidence")
    return float(np.max(np.abs(a - b)))


def evaluate_python(problem, coordinates):
    point = np.asarray(coordinates, dtype=float)
    _, value_gradient = make_joint_free_sd_functions(problem.value, problem.value_gradient)
    value, gradient = value_gradient(point)
    params = app.expand_params(point[:-1], problem.sizes, problem.config)
    quad = app.gauss_hermite_normal(problem.quadrature_points, sd=np.exp(point[-1]))
    scores = app.compute_person_eap(problem.idx, problem.config, params, quad)
    base = app.compute_base_eta(problem.idx, params, problem.config)
    if problem.config["model"] == "RSM":
        cumulative = np.r_[0, np.cumsum(params["steps"])]
        probabilities = [app.category_prob_rsm(base + t, cumulative)[:8] for t in (-3, 0, 3)]
    else:
        cumulative = np.column_stack([np.zeros(2), np.cumsum(params["steps_mat"], axis=1)])
        probabilities = [app.category_prob_pcm(base + t, cumulative, problem.idx["step_idx"])[:8]
                         for t in (-3, 0, 3)]
    return {"nll": float(value), "gradient": gradient.tolist(),
            "eap": scores.Estimate.tolist(), "sd": scores.SD.tolist(),
            "probabilities": np.asarray(probabilities).tolist()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    retained_path = ROOT / "validation/mml_free_sd_development_probe_20260911.json"
    retained = json.loads(retained_path.read_text())
    hashes = dict(retained["source_sha256"])
    for name, expected in hashes.items():
        if hashlib.sha256((ROOT / name).read_bytes()).hexdigest() != expected:
            raise ValueError(f"Retained point source changed: {name}")
    for path in (Path(__file__), Path(__file__).with_suffix(".R"), retained_path,
                 ROOT / "validation/known_assignment_mml_crossfit.R"):
        hashes[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
    test = runpy.run_path(str(ROOT / "tests/test_mml_free_sd_quadrature_adapter.py"))
    source = test["native_result"].__wrapped__()
    data = source["prep"]["data"][["Person", "Rater", "Task", "Criterion", "Score"]].copy()
    assert len(data) == 192 and data.Person.nunique() == 24
    assert data.Person.is_monotonic_increasing
    assert len(data.iloc[:8].drop_duplicates(["Rater", "Task", "Criterion"])) == 8
    problems, python_cases, python_fits, cases, fits = {}, {}, {}, [], []
    for record in retained["records"]:
        model = record["model"]
        native = source if model == "RSM" else app.mfrm_estimate(
            data, person_col="Person", facet_cols=["Rater", "Task", "Criterion"],
            score_col="Score", rating_min=0, rating_max=3, model="PCM", method="MML",
            noncenter_facet="Criterion", step_facet="Criterion", mml_engine="EM",
            quad_points=9, estimate_population_sd=True, population_prior_sd=1.0,
            maxit=300, reltol=1e-6, min_obs_per_element=1, min_obs_per_category=1)
        problem = prepare_app_free_sd_problem(native, quadrature_points=31)
        assert list(problem.sizes.values()) == ([0, 0, 1, 1, 2, 2] if model == "RSM"
                                                else [0, 0, 1, 1, 2, 4])
        assert problem.config["facet_names"] == ["Rater", "Task", "Criterion"]
        assert all(problem.config["facet_signs"][f] == -1 for f in problem.config["facet_names"])
        assert np.all(problem.idx["weight"] == 1)
        for f, levels in (("Rater", ["R1", "R2"]), ("Task", ["T1", "T2"]),
                          ("Criterion", ["C1", "C2"])):
            assert list(problem.config["facet_levels"][f]) == levels
        old_runs = {31: record["raw"]["primary_run"], 61: record["raw"]["sensitivity_run"]}
        legacy = np.asarray(old_runs[31]["primary_polish"]["initial_joint_coordinates"])
        points = {"legacy": legacy}
        points.update({f"q{q}": np.asarray(run["restart_polish"]["joint_coordinates"])
                       for q, run in old_runs.items()})
        points["perturbed"] = points["q61"] + np.linspace(-0.15, 0.15, len(legacy))
        for q in ORDERS:
            current = replace(problem, quadrature_points=q)
            problems[model, q] = current
            for label, point in points.items():
                key = f"{model}:{label}:Q{q}"
                cases.append({"id": key, "model": model, "q": q, "coordinates": point.tolist()})
                python_cases[key] = evaluate_python(current, point)
            if q in old_runs:
                run = old_runs[q]
            else:
                start = points["q61"]
                run = run_free_sd_stationarity_v2(start[:-1], np.exp(start[-1]),
                    current.value, current.value_gradient, observations=current.observations,
                    structural_bounds=current.structural_bounds, sigma_bounds=current.sigma_bounds,
                    options=JointPolishOptions(maxiter=250, gtol=1e-8, ftol=1e-15),
                    constraint_residual_function=current.constraint_residual).to_dict()
            point = run["restart_polish"]["joint_coordinates"]
            python_fits[f"{model}:Q{q}"] = {"run": run, "values": evaluate_python(current, point)}
            fits.append({"model": model, "q": q, "start": legacy.tolist(), "python_coordinates": point})
            print(f"{model} Q{q}: Python point ready", flush=True)
    dump(output / "input.json", {"classification": CLASSIFICATION, "limits": LIMITS, "source_sha256": hashes,
        "data": data.to_dict(orient="records"), "cases": cases, "fits": fits})
    dump(output / "python.json", {"cases": python_cases, "fits": python_fits})
    subprocess.run(["Rscript", "--vanilla", str(Path(__file__).with_suffix(".R")),
        str(output / "input.json"), str(ROOT / "validation/known_assignment_mml_crossfit.R"),
        str(output / "r.json")], check=True, cwd=ROOT)
    r = json.loads((output / "r.json").read_text())
    assert r["classification"] == CLASSIFICATION and r["scientific_inference_ready"] is False
    assert [row["id"] for row in r["cases"]] == [row["id"] for row in cases]
    assert [(row["model"], row["q"]) for row in r["refits"]] == list(problems)
    rule_checks, cross_checks, refit_checks, movements = [], [], [], []
    for q in ORDERS:
        rq = r["rules"][str(q)]
        py = app.gauss_hermite_normal(q)
        w, x = np.asarray(rq["weights"]), np.asarray(rq["nodes"])
        assert np.all(w > 0)
        check = {"q": q, "node": difference(x, py["nodes"]),
                 "relative_weight": float(np.max(np.abs(w / py["weights"] - 1))),
                 "moment": difference([w.sum(), w @ x, w @ x**2, w @ x**4], [1, 0, 1, 3])}
        check["pass"] = all(check[k] <= LIMITS[k] for k in ("node", "relative_weight", "moment"))
        rule_checks.append(check)
    for row in r["cases"]:
        py = python_cases[row["id"]]
        check = {"id": row["id"], "nll": abs(py["nll"] - row["nll"]),
                 "gradient": max(difference(py["gradient"], row["gradient_h"]),
                                 difference(py["gradient"], row["gradient_half"])),
                 "eap": difference(py["eap"], row["eap"]), "sd": difference(py["sd"], row["sd"]),
                 "probability": difference(py["probabilities"], row["probabilities"])
                                if row["probabilities"] is not None else 0.0}
        check["pass"] = all(check[k] <= LIMITS[k] for k in ("nll", "gradient", "eap", "sd", "probability"))
        cross_checks.append(check)
    for row in r["refits"]:
        model, q = row["model"], row["q"]
        py = python_fits[f"{model}:Q{q}"]
        at_r = evaluate_python(problems[model, q], row["coordinates"])
        check = {"model": model, "q": q, "r_convergence": row["convergence"],
                 "refit_coordinate": difference(row["coordinates"], py["run"]["restart_polish"]["joint_coordinates"]),
                 "refit_nll": abs(at_r["nll"] - py["values"]["nll"]),
                 "refit_score": max(abs(np.asarray(at_r["gradient"]))),
                 "r_point_nll_parity": abs(at_r["nll"] - row["at_r_fit"]["nll"]),
                 "r_point_gradient_parity": difference(at_r["gradient"], row["at_r_fit"]["gradient_half"])}
        check["pass"] = bool(row["convergence"] == 0 and all(check[k] <= LIMITS[k] for k in
            ("refit_coordinate", "refit_nll", "refit_score")) and check["r_point_nll_parity"] <= LIMITS["nll"]
            and check["r_point_gradient_parity"] <= LIMITS["gradient"])
        check["numeric_checks_pass"] = bool(all(check[k] <= LIMITS[k] for k in
            ("refit_coordinate", "refit_nll", "refit_score")) and check["r_point_nll_parity"] <= LIMITS["nll"]
            and check["r_point_gradient_parity"] <= LIMITS["gradient"])
        restart = row["diagnostic_restart"]
        restart_py = evaluate_python(problems[model, q], restart["coordinates"])
        check["diagnostic_restart"] = {"convergence": restart["convergence"],
            "displacement": restart["displacement"], "nll_gain": restart["nll_gain"],
            "python_score": float(np.max(np.abs(restart_py["gradient"]))),
            "nll_parity": abs(restart_py["nll"] - restart["at_point"]["nll"]),
            "gradient_parity": difference(restart_py["gradient"], restart["at_point"]["gradient_half"])}
        refit_checks.append(check)
        continuous = row["continuous_at_python_fit"]
        movements.append({"model": model, "q": q,
            "sigma": py["run"]["final_sigma"], "joint_score": py["run"]["final_projected_gradient_supnorm"],
            "nll": py["values"]["nll"], "continuous_nll_at_same_point": continuous["nll"],
            "nll_vs_continuous": abs(py["values"]["nll"] - continuous["nll"]),
            "eap_vs_continuous": difference(py["values"]["eap"], continuous["eap"]),
            "sd_vs_continuous": difference(py["values"]["sd"], continuous["sd"])})
    for name, expected in hashes.items():
        assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == expected, name
    result = {"classification": CLASSIFICATION, "qualification_eligible": False,
        "scientific_inference_ready": False, "limits": LIMITS, "source_sha256": hashes,
        "runtime": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__, **r["runtime"]},
        "fixture": {"generation_model": "RSM", "seed": 20260811, "persons": 24, "ratings": 192},
        "quadrature": rule_checks, "cross_evaluation": cross_checks, "refits": refit_checks,
        "integration_at_python_fits": movements,
        "artifact_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in output.iterdir()}}
    result["implementation_checks_pass"] = all(c["pass"] for c in rule_checks + cross_checks + refit_checks)
    result["same_point_checks_pass"] = all(c["pass"] for c in rule_checks + cross_checks)
    result["refit_numeric_checks_pass"] = all(c["numeric_checks_pass"] for c in refit_checks)
    dump(output / "summary.json", result)
    print(json.dumps({"implementation_checks_pass": result["implementation_checks_pass"],
                      "integration_at_python_fits": movements}, indent=2), flush=True)
    if not result["implementation_checks_pass"]:
        raise SystemExit("Development comparison needs review; all raw results retained.")


if __name__ == "__main__":
    main()
