#!/usr/bin/env python3
"""Independent complete-crossing PCM numerical probes; no inference qualification."""
import argparse
from dataclasses import replace
import json
from pathlib import Path
import platform
import subprocess
import sys
import time

import numpy as np
import pandas as pd
import scipy
from scipy.integrate import quad_vec
from scipy.optimize import brentq, minimize
from scipy.special import logsumexp, ndtr, roots_hermitenorm

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from mml_pcm_continuous_refit import dump, sha, fd_check
from mml_pcm_conquest_history_audit import finite_gh
from mml_pcm_external_review import fixed_grid
from validation.mml_python_analytic_sd_probe import native_fit, scores
from validation.mml_free_sd_stationarity_adapter import prepare_app_free_sd_problem

class PCMIntegral:
    # ponytail: complete crossings only; use the app observation design for missing-data validation.
    # Adapt the preserved historical integrator; never mutate its hash-bound source.
    def __init__(self, data):
        n = len({d["Person"] for d in data})
        assert n > 0
        expected = [(f"P{p:02d}", f"R{r}", f"C{c}")
                    for p in range(1, n + 1) for r in range(1, 5) for c in range(1, 6)]
        assert [(d["Person"], d["Rater"], d["Criterion"]) for d in data] == expected
        y = np.asarray([d["Score"] for d in data])
        assert np.all(np.isin(y, range(5)))
        self.y = y.astype(int).reshape(n, 20)
        self.total = self.y.sum(axis=1)
        self.k = np.arange(5.)
        self.design = np.zeros((20, 5, 23))
        step = np.array([[0, 0, 0], [-1, 0, 0], [-1, -1, 0],
                         [-1, -1, -1], [0, 0, 0]])
        for r in range(4):
            for c in range(5):
                i = r * 5 + c
                if r < 3:
                    self.design[i, :, r] = -self.k
                else:
                    self.design[i, :, :3] = self.k[:, None]
                self.design[i, :, 3 + c] = -self.k
                self.design[i, :, 8 + c * 3:11 + c * 3] = step
        self.observed = self.design[np.arange(20)[None, :], self.y].sum(axis=1)

    def evaluate(self, par, bound=12, rel_tol=1e-10, abs_tol=1e-12):
        par = np.asarray(par)
        assert par.shape == (24,) and np.all(np.isfinite(par))
        sigma = np.exp(par[-1])
        offset = self.design @ par[:-1]
        observed_offset = self.observed @ par[:-1]

        def conditional(z):
            logits = offset + sigma * z * self.k
            normalizer = logsumexp(logits, axis=1)
            prob = np.exp(logits - normalizer[:, None])
            ll = sigma * z * self.total + observed_offset - normalizer.sum()
            return ll - z * z / 2 - np.log(2 * np.pi) / 2, prob

        def mode_score(z, total):
            return sigma * (total - (conditional(z)[1] @ self.k).sum()) - z

        modes = np.array([brentq(mode_score, -bound, bound, args=(s,), xtol=1e-13)
                          for s in self.total])
        assert np.max(np.abs(modes)) < bound - .1
        centers = np.array([conditional(z)[0][i] for i, z in enumerate(modes)])

        def integrand(z):
            log_joint, prob = conditional(z)
            mass = np.exp(log_joint - centers)
            structural = self.observed - np.einsum("ik,ikd->d", prob, self.design)
            log_sd = sigma * z * (self.total - (prob @ self.k).sum())
            return mass[:, None] * np.column_stack(
                [np.ones(len(self.y)), structural, log_sd, np.full(len(self.y), z), np.full(len(self.y), z*z)])

        values, error, info = quad_vec(integrand, -bound, bound, epsrel=rel_tol,
            epsabs=abs_tol, norm="max", points=np.unique(np.r_[0., modes]),
            quadrature="gk21", workers=1, full_output=True)
        assert info.success and np.all(np.isfinite(values)) and np.all(values[:, 0] > 0)
        normalized = values / values[:, :1]
        log_mass = centers + np.log(values[:, 0])
        ez, ez2 = normalized[:, -2], normalized[:, -1]
        assert np.all(ez2 > ez**2)
        return dict(nll=float(-log_mass.sum()), gradient=-normalized[:, 1:25].sum(axis=0),
            eap=sigma*ez, sd=sigma*np.sqrt(ez2-ez**2), sigma=float(sigma),
            density_log_sigma_score=float(-np.sum(ez2-1)),
            numeric_relative_mass_error_sum=float(np.sum(error / values[:, 0])),
            tail_relative_mass_bound_sum=float(np.sum(np.exp(np.log(2*ndtr(-bound))-log_mass))),
            quadrature_neval=int(info.neval), quadrature_status=int(info.status))


def difference(a, b, keys=("nll", "eap", "sd")):
    return {k: float(np.max(abs(np.asarray(a[k]) - b[k]))) for k in keys}


def simulate(spec, condition):
    rng = np.random.Generator(np.random.PCG64(condition["seed"]))
    theta = rng.normal(0, condition["sigma"], spec["persons"])
    raters, criteria = np.array(spec["raters"]), np.array(spec["criteria"])
    steps = (np.array(spec["base_steps"]) + np.outer(
        spec["criterion_step_multipliers"], spec["step_adjustment"]))
    assert abs(raters.sum()) < 1e-14 and np.max(abs(steps.sum(axis=1))) < 1e-14
    cumulative = np.column_stack([np.zeros(5), np.cumsum(steps, axis=1)])
    logits = ((theta[:, None, None] - raters[None, :, None] - criteria[None, None, :])
              [:, :, :, None] * np.arange(5) - cumulative)
    prob = np.exp(logits - logsumexp(logits, axis=-1, keepdims=True))
    cdf = prob.cumsum(axis=-1)
    cdf[..., -1] = 1.0
    y = (rng.random((spec["persons"], 4, 5, 1)) > cdf).sum(axis=-1)
    data = [dict(Person=f"P{p+1:02d}", Rater=f"R{r+1}", Criterion=f"C{c+1}",
                 Score=int(y[p, r, c]))
            for p in range(spec["persons"]) for r in range(4) for c in range(5)]
    return data, theta, np.r_[raters[:3], criteria, steps[:, :3].ravel(), np.log(condition["sigma"])]


def fit_pair(evaluate, starts, spec):
    results = {}
    bounds = [(None, None)] * 23 + [tuple(np.log(spec["sigma_bounds"]))]
    for name, start in starts.items():
        point, runs = np.array(start), []
        for _ in range(spec["passes_per_start"]):
            begin = time.monotonic()
            def objective(p):
                v = evaluate(p)
                return v["nll"], v["gradient"]
            opt = minimize(objective, point, jac=True, method="L-BFGS-B",
                           bounds=bounds, options=spec["optimizer"])
            runs.append(dict(initial=point, returned=opt.x.copy(), success=bool(opt.success),
                status=int(opt.status), message=str(opt.message), nit=int(opt.nit), nfev=int(opt.nfev),
                elapsed=time.monotonic()-begin, reported_nll=float(opt.fun)))
            point = opt.x.copy()
        value = evaluate(point)
        results[name] = dict(coordinates=point, runs=runs, value=value)
    selected = min(results, key=lambda name: results[name]["value"]["nll"])
    return dict(starts=results, selected=selected, coordinates=results[selected]["coordinates"],
        value=results[selected]["value"],
        cross_start_coordinates=float(np.max(np.ptp([v["coordinates"] for v in results.values()], axis=0))),
        cross_start_nll=float(np.ptp([v["value"]["nll"] for v in results.values()])))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    out = parser.parse_args().output_dir.resolve()
    protocol = Path(__file__).with_name("mml_pcm_independent_protocol.json")
    assert sha(protocol) == "b732c535860e554e9b6a3028ac28c3fd41c2992d48baa87c06517e5592a76b05"
    spec = json.loads(protocol.read_text())
    limits = spec["implementation_limits"]
    sources = [Path(__file__), Path(__file__).with_suffix(".R"), protocol,
        ROOT / "streamlit_app.py", ROOT / "validation/mml_pcm_continuous_refit.py",
        ROOT / "validation/mml_pcm_conquest_history_audit.py",
        ROOT / "validation/mml_pcm_external_review.py",
        ROOT / "validation/mml_python_analytic_sd_probe.py",
        ROOT / "validation/mml_free_sd_stationarity_adapter.py"]
    hashes = {str(p.relative_to(ROOT)): sha(p) for p in sources}
    out.mkdir(parents=True, exist_ok=False)
    dump(out / "protocol.json", dict(spec, source_sha256=hashes,
        runtime=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__)))
    neutral = np.r_[np.zeros(8), np.tile([-1.5, -.5, .5], 5), 0.]
    all_results = {}
    for condition in spec["datasets"]:
        name = condition["id"]
        folder = out / name
        folder.mkdir()
        data, theta, truth = simulate(spec, condition)
        dump(folder / "input.json", dict(spec, condition=condition, data=data,
            latent_theta=theta, truth=truth, source_sha256=hashes, forced_extreme_persons=[]))
        problem = PCMIntegral(data)
        control = problem.evaluate(np.r_[np.zeros(23), np.log(1e-10)])
        assert abs(control["nll"] - len(data)*np.log(5)) < 1e-8
        assert np.max(abs(control["sd"] - 1e-10)) < 1e-18
        # Same design as the historical check; this constructor accepts actual new data.
        primary = lambda p: problem.evaluate(p, **spec["integration"][0])
        tight = lambda p: problem.evaluate(p, **spec["integration"][1])
        checks = {}
        if name.endswith("_1"):
            checks["continuous_fd"] = fd_check(primary, truth)
        print(name, "continuous fits starting", flush=True)
        reference = fit_pair(primary, {str(s): np.r_[neutral[:-1], np.log(s)]
                                      for s in spec["continuous_starts"]}, spec)
        reference["tight"] = tight(reference["coordinates"])
        reference["refinement"] = difference(reference["value"], reference["tight"],
                                               ("nll", "eap", "sd", "gradient"))
        dump(folder / "reference.json", reference)
        print(name, "reference sigma", np.exp(reference["coordinates"][-1]),
              "gradient", np.max(abs(reference["tight"]["gradient"])),
              "start distance", reference["cross_start_coordinates"], flush=True)

        native = native_fit(pd.DataFrame(data), "PCM")  # Only prepares the real app likelihood.
        base = prepare_app_free_sd_problem(native)
        assert native["prep"]["n_obs"] == len(data)
        assert not bool(native["summary"].iloc[0]["InferenceReady"])
        assert list(native["prep"]["levels"]["Person"]) == [f"P{i:02d}" for i in range(1, spec["persons"]+1)]
        assert [(k, v) for k, v in base.sizes.items() if v] == [("Rater", 3), ("Criterion", 5), ("steps", 15)]
        assert np.all(base.idx["weight"] == 1)
        cases = {"reference": dict(coordinates=reference["coordinates"], continuous=reference["tight"], finite={})}
        fits, grids = {}, {}
        for q in spec["orders"]:
            z, w = roots_hermitenorm(q)
            w /= np.sqrt(2*np.pi)
            assert np.all(w > 0) and abs(w.sum()-1) < 1e-14
            evaluate = lambda p: finite_gh(problem, p, z, w)
            pair = fit_pair(evaluate, dict(neutral=neutral, reference=reference["coordinates"]), spec)
            point, finite = pair["coordinates"], pair["value"]
            continuous = tight(point)
            app_problem = replace(base, quadrature_points=q)
            app_nll, app_grad = app_problem.joint_value_gradient(point)
            checks[f"app_q{q}"] = dict(difference(scores(app_problem, point), finite),
                gradient=float(np.max(abs(app_grad-finite["gradient"]))),
                value=abs(app_nll-finite["nll"]), constraint=base.constraint_residual(point[:-1]))
            if name.endswith("_1") and q == 31:
                checks["gh_fd"] = fd_check(evaluate, point)
            pair.update(same_point_continuous=continuous,
                integration_error=difference(finite, continuous),
                continuous_nll_excess_to_reference=continuous["nll"]-reference["tight"]["nll"],
                coordinate_distance_to_reference=float(np.max(abs(point-reference["coordinates"]))),
                sigma_difference_to_reference=float(np.exp(point[-1])-np.exp(reference["coordinates"][-1])))
            targets = spec["engineering_precision_targets"]
            error = pair["integration_error"]
            pair["engineering_precision_target_met"] = bool(
                error["nll"]/spec["persons"] < targets["nll_per_person"] and
                error["eap"] < targets["max_eap_logit"] and error["sd"] < targets["max_posterior_sd_logit"])
            fits[str(q)] = pair
            cases[f"q{q}"] = dict(coordinates=point, continuous=continuous, finite={str(q): finite})
            cases["reference"]["finite"][str(q)] = evaluate(reference["coordinates"])
            dump(folder / f"q{q}.json", pair)
            print(name, q, "sigma", np.exp(point[-1]), "start distance", pair["cross_start_coordinates"],
                  "integration errors", error, flush=True)
        for grid in spec["fixed_grids"]:
            bound, spacing = grid["bound"], grid["spacing"]
            value = fixed_grid(problem, reference["coordinates"],
                               np.linspace(-bound, bound, round(2*bound/spacing)+1))
            grids[str(bound)] = dict(config=grid, value=value,
                                     integration_error=difference(value, reference["tight"]))
        dump(folder / "python.json", cases)
        subprocess.run(["Rscript", str(Path(__file__).with_suffix(".R")),
                        str(folder / "input.json"), str(folder / "python.json"), str(folder)], check=True)
        r_checks = {}
        for case_name, case in cases.items():
            r = json.loads((folder / f"{case_name}_r.json").read_text())
            r_checks[case_name] = dict(continuous=difference(case["continuous"], r["continuous"][1]),
                finite={q: difference(v, r["finite"][q]) for q, v in case["finite"].items()},
                refinement=difference(*r["continuous"]),
                numeric_error=max(v["numeric_relative_mass_error_sum"] for v in r["continuous"]),
                tail_bound=max(v["tail_relative_mass_bound_sum"] for v in r["continuous"]))
        numerical = all(max(v.values()) < limits["app_arithmetic"] for k, v in checks.items() if k.startswith("app_"))
        numerical &= all(c["max_difference"] < limits["gradient_fd"]
                         for k, v in checks.items() if k.endswith("_fd") for c in v)
        numerical &= all(d["nll"] < limits["r_nll"] and max(d["eap"], d["sd"]) < limits["r_moments"]
                         for c in r_checks.values() for d in [c["continuous"], *c["finite"].values()])
        numerical &= all(max(c["refinement"].values()) < limits["reference_refinement"] and
            max(c["numeric_error"], c["tail_bound"]) < limits["relative_mass_error"] for c in r_checks.values())
        reference_ok = (reference["cross_start_coordinates"] < limits["reference_start_coordinates"] and
            reference["cross_start_nll"] < limits["reference_start_nll"] and
            np.max(abs(reference["tight"]["gradient"])) < limits["stationary_gradient"] and
            max(reference["refinement"].values()) < limits["reference_refinement"])
        all_fits = [reference, *fits.values()]
        optimization_ok = all(np.max(abs(s["value"]["gradient"])) < limits["stationary_gradient"] and
            all(run["success"] for run in s["runs"]) for pair in all_fits for s in pair["starts"].values())
        result = dict(condition=condition, reference=reference, fits=fits, fixed_grids=grids,
            checks=checks, independent_r=r_checks, arithmetic_checks_pass=bool(numerical),
            reference_checks_pass=bool(reference_ok), optimization_checks_pass=bool(optimization_ok),
            observed_extreme_persons=int(np.sum(np.isin(problem.total, [0, 80]))))
        dump(folder / "summary.json", result)
        all_results[name] = result
        print(name, "arithmetic/reference/optimization", numerical, reference_ok, optimization_ok, flush=True)
    assert hashes == {p: sha(ROOT/p) for p in hashes}
    passed = all(r["arithmetic_checks_pass"] and r["reference_checks_pass"] and r["optimization_checks_pass"]
                 for r in all_results.values())
    dump(out / "summary.json", dict(classification=spec["classification"], scientific_inference_ready=False,
        qualification_eligible=False, implementation_checks_pass=passed, results=all_results,
        source_sha256=hashes, protocol_sha256=sha(protocol),
        artifact_sha256={str(p.relative_to(out)): sha(p) for p in sorted(out.rglob("*")) if p.is_file()},
        limitations=spec["limitations"]))
    print("Implementation checks:", passed, "Saved:", out, flush=True)
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
