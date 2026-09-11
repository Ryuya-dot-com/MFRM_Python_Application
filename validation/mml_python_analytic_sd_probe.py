#!/usr/bin/env python3
"""Check the Python analytic log-SD path; development evidence, never qualification."""
import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import platform
import sys

import numpy as np
import pandas as pd
import scipy

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
import streamlit_app as app
from mfrm_app.mml_engine_v2 import run_free_sd_stationarity_v2
from mfrm_app.mml_stationarity import JointPolishOptions, audit_joint_gradient, make_joint_free_sd_functions
from validation.mml_free_sd_cross_language_probe import dump, difference
from validation.mml_free_sd_stationarity_adapter import prepare_app_free_sd_problem
from validation.mml_free_sd_stationarity_known_fixture import _json_value

REFERENCE = "validation/mml_free_sd_analytic_score_20260912.json"
REFERENCE_SHA256 = "ea8d041a63fe7bb4a1018cc70d1c24720559622798da38ce1403a8b36a52300c"
ORDERS = (31, 61, 121)
LIMITS = {"r_gradient": 1e-9, "fd_gradient": 1e-5,
          "coordinate": 1e-4, "nll": 1e-8, "score": 1e-4}
OPTIONS = JointPolishOptions(maxiter=250, gtol=1e-8, ftol=1e-15)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def native_fit(data, model):
    return app.mfrm_estimate(
        data, person_col="Person", facet_cols=[c for c in ("Rater", "Task", "Criterion") if c in data],
        score_col="Score", rating_min=0, rating_max=int(data.Score.max()), model=model, method="MML",
        noncenter_facet="Criterion", step_facet="Criterion", mml_engine="EM", quad_points=9,
        estimate_population_sd=True, population_prior_sd=1.0, maxit=300, reltol=1e-6,
        min_obs_per_element=1, min_obs_per_category=1,
    )


def new_dataset(model):
    # Distinct seeds/designs; PCM extremes are deliberate stress edits, not recovery data.
    seed = 2026091201 if model == "RSM" else 2026091202
    rng = np.random.default_rng(seed)
    n_person, n_rater, n_criterion = (30, 3, 3) if model == "RSM" else (24, 4, 5)
    theta = rng.normal(0, 1.4, n_person)
    rater = np.linspace(-0.4, 0.4, n_rater)
    criterion = np.linspace(-0.3, 0.3, n_criterion)
    steps = np.tile([-1.5, -0.5, 0.5, 1.5], (n_criterion, 1))
    if model == "PCM":
        steps += np.linspace(-1, 1, n_criterion)[:, None] * [-0.3, 0.1, 0.3, -0.1]
    data = app._simulate_recovery_responses(
        theta=theta, rater=rater, criterion=criterion, steps=steps,
        slopes=None, model=model, rating_min=0, rng=rng,
    )
    if model == "PCM":
        data.loc[data.Person == "P01", "Score"] = 0
        data.loc[data.Person == "P02", "Score"] = 4
    return data, {"seed": seed, "model": model, "theta": theta.tolist(),
        "rater": rater.tolist(), "criterion": criterion.tolist(), "adjacent_steps": steps.tolist(),
        "forced_extreme_persons": ["P01", "P02"] if model == "PCM" else [],
        "data": data.to_dict("records")}


def scores(problem, point):
    quad = app.gauss_hermite_normal(problem.quadrature_points, sd=np.exp(point[-1]))
    expanded = app.expand_params(point[:-1], problem.sizes, problem.config)
    result = app.compute_person_eap(problem.idx, problem.config, expanded, quad)
    return {"nll": problem.value(point[:-1], np.exp(point[-1])),
            "eap": result.Estimate.tolist(), "sd": result.SD.tolist()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    out = parser.parse_args().output_dir
    out.mkdir(parents=True, exist_ok=False)
    assert sha256(ROOT / REFERENCE) == REFERENCE_SHA256
    reference = json.loads((ROOT / REFERENCE).read_text())
    old_input_path = ROOT / "validation/mml_free_sd_cross_language_20260912/input.json"
    assert sha256(old_input_path) == reference["prior_artifact_sha256"]["input.json"]
    old = json.loads(old_input_path.read_text())
    paths = set(reference["source_sha256"]) | {
        str(Path(__file__).relative_to(ROOT)),
        "validation/mml_free_sd_stationarity_known_fixture.py",
    }
    hashes = {p: sha256(ROOT / p) for p in sorted(paths)}
    datasets = {m: new_dataset(m) for m in ("RSM", "PCM")}
    dump(out / "input.json", {"classification": "OBSERVED_DEVELOPMENT_ONLY",
        "qualification_eligible": False, "source_sha256": hashes, "limits": LIMITS,
        "orders": ORDERS, "options": OPTIONS.__dict__,
        "reference_sha256": {REFERENCE: REFERENCE_SHA256, str(old_input_path.relative_to(ROOT)): sha256(old_input_path)},
        "datasets": {m: spec for m, (_, spec) in datasets.items()}})

    # Compare changed Python code with untouched, hash-bound independent R scores.
    problems = {m: prepare_app_free_sd_problem(native_fit(pd.DataFrame(old["data"]), m))
                for m in ("RSM", "PCM")}
    r_cases = {c["id"]: c for c in reference["cases"]}
    parity = []
    for case in old["cases"]:
        problem = replace(problems[case["model"]], quadrature_points=case["q"])
        point = np.array(case["coordinates"])
        _, gradient = problem.joint_value_gradient(point)
        error = difference(gradient, r_cases[case["id"]]["gradient"])
        parity.append({"id": case["id"], "gradient": gradient.tolist(), "max_difference": error,
                       "pass": error < LIMITS["r_gradient"]})
    dump(out / "r_parity.json", parity)
    print("Independent R gradient parity:", sum(c["pass"] for c in parity), "/", len(parity), flush=True)

    all_checks, fits = [], []
    for model, (data, _) in datasets.items():
        native = native_fit(data, model)
        assert native["prep"]["n_obs"] == len(data)
        assert not bool(native["summary"].iloc[0]["InferenceReady"])
        base = prepare_app_free_sd_problem(native)
        start = np.r_[base.structural_start, np.log(base.sigma_start)]
        # Off-stationary points, with two additional SD edge checks at Q31.
        for q, sigma in [(q, 1.4) for q in ORDERS] + [(31, 0.08), (31, 6.0)]:
            problem = replace(base, quadrature_points=q)
            point = start + np.r_[np.linspace(-0.15, 0.2, len(start) - 1), 0]
            point[-1] = np.log(sigma)
            objective, analytic = make_joint_free_sd_functions(
                problem.value, problem.value_gradient, joint_value_gradient=problem.joint_value_gradient,
            )
            audits = [audit_joint_gradient(objective, analytic, point, relative_step=h)
                      for h in (1e-5, 5e-6)]
            all_checks.append({"model": model, "q": q, "sigma": sigma, "point": point.tolist(),
                "audits": [a.to_dict() for a in audits],
                "pass": all(a.maximum_absolute_difference < LIMITS["fd_gradient"] for a in audits)})
        for q in ORDERS:
            problem = replace(base, quadrature_points=q)
            runs = {}
            for method in ("central_difference", "analytic"):
                run = run_free_sd_stationarity_v2(
                    base.structural_start, base.sigma_start, problem.value, problem.value_gradient,
                    observations=problem.observations, structural_bounds=problem.structural_bounds,
                    sigma_bounds=problem.sigma_bounds, options=OPTIONS,
                    constraint_residual_function=problem.constraint_residual,
                    joint_value_gradient=problem.joint_value_gradient if method == "analytic" else None,
                )
                runs[method] = run
                dump(out / f"{model}_q{q}_{method}.json", _json_value(run.to_dict()))
                print(model, q, method, "terminated:", run.algorithm_terminated,
                      "score:", run.final_projected_gradient_supnorm, flush=True)
            a, b = runs["analytic"], runs["central_difference"]
            point = np.array(a.restart_polish.joint_coordinates)
            own = scores(problem, point)
            higher = scores(replace(problem, quadrature_points=181), point)
            coordinate = difference(point, b.restart_polish.joint_coordinates)
            nll = abs(a.final_objective - b.final_objective)
            fits.append({"model": model, "q": q, "coordinate_difference": coordinate,
                "nll_difference": nll, "analytic_sigma": a.final_sigma,
                "numeric_pass": coordinate < LIMITS["coordinate"] and nll < LIMITS["nll"] and
                    all(r.final_projected_gradient_supnorm < LIMITS["score"] for r in runs.values()),
                "status_pass": all(r.algorithm_terminated for r in runs.values()),
                "at_analytic_fit": own, "q181_at_same_point": higher,
                "q181_movement": {k: difference(own[k], higher[k]) for k in ("nll", "eap", "sd")}})
    passed = all(c["pass"] for c in parity + all_checks) and all(
        f["numeric_pass"] and f["status_pass"] for f in fits)
    assert hashes == {p: sha256(ROOT / p) for p in hashes}, "sources changed during run"
    dump(out / "summary.json", {"classification": "OBSERVED_DEVELOPMENT_ONLY",
        "scientific_inference_ready": False, "qualification_eligible": False,
        "implementation_checks_pass": passed, "source_sha256": hashes, "limits": LIMITS,
        "r_parity_max_difference": max(c["max_difference"] for c in parity),
        "new_data_gradient_checks": all_checks, "fits": fits,
        "runtime": {"python": platform.python_version(), "numpy": np.__version__,
                    "scipy": scipy.__version__, "platform": platform.platform()},
        "artifact_sha256": {p.name: sha256(p) for p in sorted(out.iterdir())},
        "limitations": ["Two development datasets; no scientific acceptance or recovery study",
            "Q181 comparison is a finite-rule comparison, not a continuous reference",
            "All seeds/designs are now observed and excluded from qualification"]})
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
