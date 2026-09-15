#!/usr/bin/env python3
"""Separate finite-GH optimization and integration error on one observed PCM fixture."""
import argparse
from dataclasses import replace
import json
from pathlib import Path
import platform
import subprocess
import sys

import numpy as np
import pandas as pd
import scipy

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
import streamlit_app as app
from mfrm_app.mml_engine_v2 import run_free_sd_stationarity_v2
from validation.mml_python_analytic_sd_probe import native_fit, scores, sha256, OPTIONS
from validation.mml_free_sd_cross_language_probe import dump, difference
from validation.mml_free_sd_stationarity_adapter import prepare_app_free_sd_problem
from validation.mml_free_sd_stationarity_known_fixture import _json_value

PREVIOUS = ROOT / "validation/mml_python_analytic_sd_20260912"
SUMMARY_SHA256 = "02315c5b9f10c9a9ebd475be25bee81028632dcec1d46a14fb83c0ef8b65e3f2"
START_ORDERS = (31, 61, 121)
ORDERS = (*START_ORDERS, 181)
# Implementation/replay checks, fixed before this run; no scientific margins.
LIMITS = {"replay": 1e-8, "r_nll": 1e-9, "r_score": 1e-10,
          "r_probability": 1e-12, "continuous_refinement": 1e-8,
          "relative_mass_error": 1e-8, "gradient": 1e-4}
INTEGRATION = [
    {"bound": 12, "rel_tol": 1e-10, "abs_tol": 1e-12},
    {"bound": 14, "rel_tol": 1e-12, "abs_tol": 1e-13},
]


def evaluate(base, point):
    expanded = app.expand_params(point[:-1], base.sizes, base.config)
    # The independent R code reconstructs these same coordinates, not these arrays.
    eta = app.compute_base_eta(base.idx, expanded, base.config)
    cumulative = np.column_stack([np.zeros(5), np.cumsum(expanded["steps_mat"], axis=1)])
    probabilities = [app.category_prob_pcm(eta + theta, cumulative, base.idx["step_idx"])[:20]
                     for theta in (-3, 0, 3)]
    return {"coordinates": point.tolist(), "sigma": float(np.exp(point[-1])),
            "finite": {str(q): scores(replace(base, quadrature_points=q), point) for q in ORDERS},
            "probabilities": np.asarray(probabilities).tolist()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    out = parser.parse_args().output_dir.resolve()
    out.mkdir(parents=True, exist_ok=False)
    assert sha256(PREVIOUS / "summary.json") == SUMMARY_SHA256
    previous = json.loads((PREVIOUS / "summary.json").read_text())
    assert not previous["scientific_inference_ready"] and not previous["qualification_eligible"]
    for name, expected in previous["artifact_sha256"].items():
        assert sha256(PREVIOUS / name) == expected, name
    spec = json.loads((PREVIOUS / "input.json").read_text())["datasets"]["PCM"]
    data = pd.DataFrame(spec["data"])
    assert spec["seed"] == 2026091202 and len(data) == 480 and data.Person.nunique() == 24
    paths = set(previous["source_sha256"]) | {
        str(Path(__file__).relative_to(ROOT)),
        str(Path(__file__).with_suffix(".R").relative_to(ROOT)),
        "validation/mml_python_analytic_sd_probe.py",
        "validation/mml_free_sd_cross_language_probe.py",
    }
    hashes = {p: sha256(ROOT / p) for p in sorted(paths)}
    starts = {q: np.asarray(json.loads((PREVIOUS / f"PCM_q{q}_analytic.json").read_text())
                           ["restart_polish"]["joint_coordinates"]) for q in START_ORDERS}
    dump(out / "input.json", {"classification": "OBSERVED_DEVELOPMENT_ONLY",
        "scientific_inference_ready": False, "qualification_eligible": False,
        "source_sha256": hashes, "previous_summary_sha256": SUMMARY_SHA256,
        "previous_artifact_sha256": previous["artifact_sha256"],
        "data": spec["data"], "seed": spec["seed"], "forced_extreme_persons": spec["forced_extreme_persons"],
        "starts": {str(q): p.tolist() for q, p in starts.items()}, "orders": ORDERS,
        "options": OPTIONS.__dict__, "limits": LIMITS, "integration": INTEGRATION,
        "design": "3 retained starts x 4 finite GH orders; one polish and one restart per cell; no retries",
        "continuous_scope": "Same-point independent R log-PCM integrals; no continuous-objective refitting"})

    native = native_fit(data, "PCM")  # Prepare the app's existing identified design.
    base = prepare_app_free_sd_problem(native)
    assert native["prep"]["n_obs"] == 480 and not bool(native["summary"].iloc[0]["InferenceReady"])
    assert base.config["facet_names"] == ["Rater", "Criterion"]
    assert base.config["facet_levels"] == {"Rater": ["R1", "R2", "R3", "R4"],
                                           "Criterion": ["C1", "C2", "C3", "C4", "C5"]}
    assert [(k, v) for k, v in base.sizes.items() if v] == [("Rater", 3), ("Criterion", 5), ("steps", 15)]
    assert np.all(base.idx["weight"] == 1) and base.constraint_residual(starts[31][:-1]) < 1e-12
    assert list(native["prep"]["levels"]["Person"]) == [f"P{i:02d}" for i in range(1, 25)]
    cases, runs, replay = {}, {}, []
    for q, point in starts.items():
        case = evaluate(base, point)
        case.update(kind="retained", q=q)
        cases[f"retained_q{q}"] = case
        old = next(f for f in previous["fits"] if f["model"] == "PCM" and f["q"] == q)
        differences = {k: difference(case["finite"][str(q)][k], old["at_analytic_fit"][k])
                       for k in ("nll", "eap", "sd")}
        replay.append({"q": q, "differences": differences,
                       "pass": max(differences.values()) < LIMITS["replay"]})
    dump(out / "replay.json", replay)
    assert all(r["pass"] for r in replay), "Retained-point replay failed"
    print("Retained-point replay: 3/3", flush=True)

    for q in ORDERS:
        problem = replace(base, quadrature_points=q)
        for source_q, start in starts.items():
            name = f"q{q}_from_q{source_q}"
            run = run_free_sd_stationarity_v2(
                start[:-1], np.exp(start[-1]), problem.value, problem.value_gradient,
                observations=problem.observations, structural_bounds=problem.structural_bounds,
                sigma_bounds=problem.sigma_bounds, options=OPTIONS,
                constraint_residual_function=problem.constraint_residual,
                joint_value_gradient=problem.joint_value_gradient,
            )
            record = _json_value(run.to_dict())
            dump(out / f"{name}.json", record)
            runs[name] = record
            case = evaluate(base, np.asarray(run.restart_polish.joint_coordinates))
            case.update(kind="cross_start", q=q, source_q=source_q)
            cases[name] = case
            print(name, "sigma", run.final_sigma, "nll", run.final_objective,
                  "terminated", run.algorithm_terminated, "score", run.final_projected_gradient_supnorm, flush=True)
    dump(out / "python.json", cases)
    subprocess.run(["Rscript", str(Path(__file__).with_suffix(".R")),
                    str(out / "input.json"), str(out / "python.json"), str(out)], check=True)

    comparisons = []
    for name, case in cases.items():
        r = json.loads((out / f"{name}_r.json").read_text())
        same = {str(q): {k: difference(case["finite"][str(q)][k], r["finite"][str(q)][k])
                        for k in ("nll", "eap", "sd")} for q in ORDERS}
        probability = difference(case["probabilities"], r["probabilities"])
        primary, tight = r["continuous"]
        refinement = {k: difference(primary[k], tight[k]) for k in ("nll", "eap", "sd")}
        integration_ok = all(
            c["numeric_relative_mass_error_sum"] < LIMITS["relative_mass_error"] and
            c["tail_relative_mass_bound_sum"] < LIMITS["relative_mass_error"]
            for c in (primary, tight))
        errors = {str(q): {k: difference(case["finite"][str(q)][k], tight[k])
                           for k in ("nll", "eap", "sd")} for q in ORDERS}
        passed = (probability < LIMITS["r_probability"] and integration_ok and
                  max(refinement.values()) < LIMITS["continuous_refinement"] and
                  all(d["nll"] < LIMITS["r_nll"] and max(d["eap"], d["sd"]) < LIMITS["r_score"]
                      for d in same.values()))
        comparisons.append({"id": name, "kind": case["kind"], "q": case["q"],
            "sigma": case["sigma"], "same_point_r_differences": same,
            "probability_difference": probability, "continuous_refinement": refinement,
            "continuous_nll": tight["nll"], "continuous_log_sigma_nll_score": tight["log_sigma_nll_score"],
            "continuous_errors": errors, "integration_error_checks_pass": integration_ok,
            "pass": passed})
    spread = []
    for q in ORDERS:
        names = [f"q{q}_from_q{s}" for s in START_ORDERS]
        group = [cases[n] for n in names]
        spread.append({"q": q, "ids": names,
            "sigma_range": [min(c["sigma"] for c in group), max(c["sigma"] for c in group)],
            "finite_nll_range": [min(c["finite"][str(q)]["nll"] for c in group),
                                 max(c["finite"][str(q)]["nll"] for c in group)],
            "max_coordinate_distance": max(difference(a["coordinates"], b["coordinates"])
                                           for a in group for b in group)})
    terminated = all(r["algorithm_terminated"] for r in runs.values())
    stationary = all(r["final_projected_gradient_supnorm"] < LIMITS["gradient"] for r in runs.values())
    assert hashes == {p: sha256(ROOT / p) for p in hashes}, "Sources changed during run"
    assert sha256(PREVIOUS / "summary.json") == SUMMARY_SHA256
    assert all(sha256(PREVIOUS / n) == h for n, h in previous["artifact_sha256"].items())
    passed = all(c["pass"] for c in comparisons) and terminated and stationary
    dump(out / "summary.json", {"classification": "OBSERVED_DEVELOPMENT_ONLY",
        "scientific_inference_ready": False, "qualification_eligible": False,
        "implementation_checks_pass": passed, "same_point_checks_pass": all(c["pass"] for c in comparisons),
        "all_refits_terminated": terminated, "all_refits_small_projected_gradient": stationary,
        "comparisons": comparisons, "cross_start_spread": spread,
        "limits": LIMITS, "source_sha256": hashes,
        "runtime": {"python": platform.python_version(), "numpy": np.__version__,
                    "scipy": scipy.__version__, "platform": platform.platform()},
        "artifact_sha256": {p.name: sha256(p) for p in sorted(out.iterdir()) if p.is_file()},
        "limitations": ["One observed stress dataset; not a recovery or qualification study",
            "Adaptive integral error estimates are not rigorous interval certificates",
            "Three observed starts do not prove global optimality",
            "No continuous-objective optimum, SE, CI or coverage qualification",
            "All primary and restart termination codes are retained, including failures"]})
    print("Implementation checks:", passed, "scientific inference ready: False", flush=True)
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
