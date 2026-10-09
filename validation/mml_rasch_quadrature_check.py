"""High-order GH refits and independent continuous-information development checks.

Uses every retained regular/sparse/unbalanced RSM/PCM dataset. Prior failures
are preserved. Sources, settings and tolerances are frozen before new fits.
"""

from dataclasses import asdict, replace
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import subprocess
import sys
import zipfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd

import streamlit_app as app
from mfrm_app.mml_engine_v2 import run_free_sd_two_stage
from mfrm_app.mml_stationarity import JointPolishOptions
from validation.mml_free_sd_stationarity_adapter import prepare_app_free_sd_problem_from_data
from validation.mml_joint_information_check import save, expanded_coordinates
from validation.mml_rasch_integral import RaschIntegral


PROTOCOL = {
    "classification": "OBSERVED_DEVELOPMENT_ONLY", "scientific_inference_ready": False,
    "orders": [121, 181, 241], "primary_start": "retained Q61 endpoint",
    "additional_start": "retained Q31 endpoint at Q241 for every dataset; no best-start selection",
    "continuous_settings": [dict(bound=12., rel_tol=1e-9, abs_tol=1e-11),
                            dict(bound=16., rel_tol=1e-11, abs_tol=1e-13)],
    "preliminary": asdict(JointPolishOptions(maxiter=300, gtol=1e-8, ftol=1e-15)),
    "refinement": asdict(JointPolishOptions(maxiter=300, gtol=1e-10, ftol=0.)),
    "anchor_gradient_limit": 1e-3,
    "limits": dict(finite_nll=1e-8, finite_gradient=1e-8, finite_information=1e-5,
        continuous_refinement_nll=1e-8, continuous_refinement_gradient=1e-8,
        continuous_refinement_information=1e-7, continuous_refinement_moments=1e-9,
        relative_mass_error=1e-7, tail_mass=1e-10, tail_gradient=1e-8, tail_information=1e-7,
        numeric_information_allowance=1e-5, density_score_identity=1e-8,
        r_nll=1e-8, r_moments=1e-9, r_log_sigma_score=1e-8,
        gh_nll_per_person=1e-6, gh_score_supnorm=1e-4, gh_moments=1e-4,
        gh_se_relative=0.005, same_q_cross_start_coordinates=1e-6),
    "failure_policy": "Retain every planned point and solver status; no retuning, replacement or extra start selection",
    "scope": "Six observed development datasets, not recovery/coverage or prospective qualification",
}


def difference(a, b, key):
    return float(np.max(np.abs(np.asarray(a[key])-np.asarray(b[key]))))


def evaluate(problem, reference, point):
    q = app.gauss_hermite_normal(problem.quadrature_points)
    if not np.all(q["weights"] > 0):
        raise FloatingPointError("GH weights are not positive")
    finite = reference.finite(point, q["nodes"], q["weights"])
    nll, gradient = problem.joint_value_gradient(point)
    info = app._compute_mml_observed_information(point, problem.idx, problem.config, problem.sizes,
        app.gauss_hermite_normal(problem.quadrature_points, sd=np.exp(point[-1])), include_log_sigma=True)
    continuous = [reference.continuous(point, **s, information=True) for s in PROTOCOL["continuous_settings"]]
    tight = continuous[-1]; limits = PROTOCOL["limits"]
    refinement = {key: difference(*continuous, key) for key in ("nll", "gradient", "information", "eap", "sd")}
    integration = {key: difference(finite, tight, key) for key in ("nll", "gradient", "information", "eap", "sd")}
    covariance, se_error = {}, None
    for name, matrix in (("finite", info), ("continuous", tight["information"])):
        eigen = np.linalg.eigvalsh(matrix)
        if eigen.min() > 0:
            covariance[name] = np.linalg.inv(matrix)
    if len(covariance) == 2:
        se = {name: np.sqrt(np.diag(cov)) for name, cov in covariance.items()}
        se_error = float(np.max(abs(se["finite"]/se["continuous"]-1)))
    checks = dict(
        finite_nll=abs(nll-finite["nll"]) < limits["finite_nll"],
        finite_gradient=np.max(abs(gradient-finite["gradient"])) < limits["finite_gradient"],
        finite_information=np.max(abs(info-finite["information"])) < limits["finite_information"],
        reference_nll=refinement["nll"] < limits["continuous_refinement_nll"],
        reference_gradient=refinement["gradient"] < limits["continuous_refinement_gradient"],
        reference_information=refinement["information"] < limits["continuous_refinement_information"],
        reference_moments=max(refinement["eap"], refinement["sd"]) < limits["continuous_refinement_moments"],
        reference_mass=tight["numeric_relative_mass_error_sum"] < limits["relative_mass_error"],
        reference_tail=tight["tail_relative_mass_bound_sum"] < limits["tail_mass"] and
            tight["tail_gradient_abs_bound"] < limits["tail_gradient"] and tight["tail_information_abs_bound"] < limits["tail_information"],
        reference_information_allowance=tight["numeric_information_error_allowance"] < limits["numeric_information_allowance"],
        density_score=abs(tight["gradient"][-1]-tight["density_log_sigma_score"]) < limits["density_score_identity"],
        gh_nll=integration["nll"]/reference.n < limits["gh_nll_per_person"],
        gh_score=integration["gradient"] < limits["gh_score_supnorm"],
        gh_moments=max(integration["eap"], integration["sd"]) < limits["gh_moments"],
        gh_se=se_error is not None and se_error < limits["gh_se_relative"],
    )
    return dict(point=point, finite=finite, continuous=continuous, covariance=covariance,
        app_information=info, app_gradient=gradient, refinement=refinement, integration=integration,
        se_relative_difference=se_error, checks=checks)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("parent", type=Path); parser.add_argument("output", type=Path)
    args = parser.parse_args(); parent, out = args.parent.resolve(), args.output.resolve()
    for name, digest in json.loads((parent/"manifest.json").read_text()).items():
        if hashlib.sha256((parent/name).read_bytes()).hexdigest() != digest:
            raise ValueError(f"Retained artifact changed: {name}")
    out.mkdir(parents=True, exist_ok=False)
    save(out/"protocol.json", dict(PROTOCOL, parent_manifest_sha256=hashlib.sha256((parent/"manifest.json").read_bytes()).hexdigest()))
    paths = sorted(set([Path(__file__).resolve(), ROOT/"streamlit_app.py", *(ROOT/"mfrm_app").glob("*.py"),
        *(ROOT/"validation").glob("mml_free_sd_*adapter.py"), ROOT/"validation/mml_joint_information_check.py",
        ROOT/"validation/mml_joint_information_reference.R", ROOT/"validation/mml_rasch_integral.py"]))
    with zipfile.ZipFile(out/"sources.zip", "w", zipfile.ZIP_DEFLATED) as archive:
        for path in paths: archive.write(path, str(path.relative_to(ROOT)))
    save(out/"source_manifest.json", {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    save(out/"runtime.json", dict(python=sys.version, packages={d.metadata["Name"]: d.version for d in importlib.metadata.distributions()}))
    (out/"reference.R").write_bytes((ROOT/"validation/mml_joint_information_reference.R").read_bytes())
    planned = json.loads((parent/"prefit_manifest.json").read_text())
    records, failures, r_cases = {}, [], []
    for case in planned:
        case_id = case["id"]
        data = pd.read_csv(parent/(case_id+".csv"))
        (out/(case_id+".csv")).write_bytes((parent/(case_id+".csv")).read_bytes())
        _, base = prepare_app_free_sd_problem_from_data(data, **case["input_config"])
        reference = RaschIntegral(data, base.config)
        starts = {q:np.array(json.loads((parent/f"{case_id}_Q{q}.json").read_text())["point"]) for q in (31, 61)}
        jobs = [(q, "retained", starts[q]) for q in (31, 61)]
        jobs += [(q, "from61", starts[61]) for q in PROTOCOL["orders"]]
        jobs += [(241, "from31", starts[31])]
        for q, label, start in jobs:
            key = f"{case_id}_Q{q}_{label}"
            try:
                problem = replace(base, quadrature_points=q); point = start
                runs = None
                if label != "retained":
                    runs = run_free_sd_two_stage(start[:-1], np.exp(start[-1]), problem.value, problem.value_gradient,
                        observations=problem.observations, structural_bounds=problem.structural_bounds,
                        sigma_bounds=problem.sigma_bounds, difference_factory=problem.likelihood_difference,
                        joint_value_gradient=problem.joint_value_gradient, constraint_residual_function=problem.constraint_residual,
                        preliminary_options=JointPolishOptions(**PROTOCOL["preliminary"]),
                        refinement_options=JointPolishOptions(**PROTOCOL["refinement"]), anchor_gradient_limit=PROTOCOL["anchor_gradient_limit"])
                    save(out/(key+"_runs.json"), runs.to_dict())
                    final = runs.refinement if runs.refinement is not None else runs.preliminary
                    point = np.array(final.restart_polish.joint_coordinates)
                record = evaluate(problem, reference, point)
                record["solver_terminated"] = None if runs is None else runs.refinement is not None and runs.refinement.algorithm_terminated
                record["expanded"] = expanded_coordinates(problem, point)
                records[key] = record; save(out/(key+".json"), record)
                if q in (31, 61, 241):
                    r_cases.append(dict(id=key, config={k:base.config[k] for k in ("model", "facet_names", "facet_levels",
                        "noncenter_facet", "facet_signs", "step_facet", "n_cat", "rating_min", "rating_max")},
                        data=data.to_dict("records"), point=point, continuous_only=True))
                print(key, "complete; failed checks:", [k for k,v in record["checks"].items() if not v], flush=True)
            except Exception as exc:
                failures.append(dict(id=key, type=type(exc).__name__, message=str(exc)))
                save(out/"failures.json", failures)
                print(key, "FAILED", repr(exc), flush=True)
    save(out/"R_input.json", dict(classification="DEVELOPMENT_ONLY", continuous_settings=PROTOCOL["continuous_settings"], cases=r_cases))
    with (out/"R_execution.log").open("w") as log:
        result = subprocess.run(["/usr/local/bin/Rscript", str(out/"reference.R"), str(out/"R_input.json"), str(out/"R_results.json")], stdout=log, stderr=subprocess.STDOUT)
    r_comparisons = []
    if result.returncode == 0:
        for case in json.loads((out/"R_results.json").read_text())["cases"]:
            python = records[case["id"]]["continuous"][-1]; r = case["continuous"][-1]
            errors = {k:difference(python, r, k) for k in ("nll", "eap", "sd", "density_log_sigma_score")}
            r_comparisons.append(dict(id=case["id"], errors=errors,
                passed=errors["nll"] < PROTOCOL["limits"]["r_nll"] and max(errors["eap"], errors["sd"]) < PROTOCOL["limits"]["r_moments"] and errors["density_log_sigma_score"] < PROTOCOL["limits"]["r_log_sigma_score"]))
    cross_start = {}
    for case in planned:
        prefix = case["id"]+"_Q241_"
        if all(prefix+x in records for x in ("from31", "from61")):
            cross_start[case["id"]] = difference(records[prefix+"from31"], records[prefix+"from61"], "point")
    summary = dict(classification=PROTOCOL["classification"], scientific_inference_ready=False,
        expected_points=36, completed_points=len(records), failures=failures, R_returncode=result.returncode,
        R_comparisons=r_comparisons, Q241_cross_start_coordinates=cross_start,
        cases={k:dict(checks=v["checks"], solver_terminated=v["solver_terminated"], integration=v["integration"],
            se_relative_difference=v["se_relative_difference"], sigma=float(np.exp(v["point"][-1]))) for k,v in records.items()})
    save(out/"summary.json", summary)
    save(out/"manifest.json", {p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(out.iterdir()) if p.is_file()})
    print("Completed", len(records), "/ 36 points; R confirmations:", len(r_comparisons), flush=True)
    return 0 if len(records) == 36 and len(r_comparisons) == 24 and not failures and all(r["passed"] for r in r_comparisons) else 1


if __name__ == "__main__":
    raise SystemExit(main())
