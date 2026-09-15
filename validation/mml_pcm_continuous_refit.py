#!/usr/bin/env python3
"""Adaptive-integral refits of the retained, fully crossed 24-person PCM fixture.

Development check only. Coordinates: three sum-zero rater contrasts, five
criterion locations, three sum-zero step contrasts per criterion, log(SD).
"""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import time

import numpy as np
import scipy
from scipy.integrate import quad_vec
from scipy.optimize import brentq, minimize
from scipy.special import logsumexp, ndtr

ROOT = Path(__file__).resolve().parents[1]
PREVIOUS = ROOT / "validation/generated/mml_pcm_quadrature_diagnosis_20260913"
SUMMARY_HASH = "14a3aa94f98f611ab07ad91f193fe01b22abb2068be1d1052002f20923da9567"
OPTIONS = dict(maxiter=250, gtol=1e-8, ftol=1e-15, maxls=50)
LIMITS = dict(r_nll=1e-8, r_moments=1e-9, gradient_fd=1e-6,
              gradient=1e-4, refinement=1e-8, relative_mass_error=1e-8)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, obj):
    with Path(path).open("x") as f:
        json.dump(obj, f, indent=2, allow_nan=False,
                  default=lambda x: x.tolist() if isinstance(x, np.ndarray) else x.item())
        f.write("\n")


class PCMIntegral:
    # ponytail: this check requires a complete 24 x 4 x 5 crossing; use the app's
    # general observation design if a later validation needs missing responses.
    def __init__(self, data):
        expected = [(f"P{p:02d}", f"R{r}", f"C{c}")
                    for p in range(1, 25) for r in range(1, 5) for c in range(1, 6)]
        assert [(d["Person"], d["Rater"], d["Criterion"]) for d in data] == expected
        y = np.asarray([d["Score"] for d in data])
        assert np.all(np.isin(y, range(5))) and np.all(y[:20] == 0) and np.all(y[20:40] == 4)
        self.y = y.astype(int).reshape(24, 20)
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
                [np.ones(24), structural, log_sd, np.full(24, z), np.full(24, z*z)])

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


def fd_check(evaluate, par):
    results = []
    analytic = evaluate(par)["gradient"]
    for h in (1e-4, 3e-5):
        fd = np.array([(evaluate(par+np.eye(24)[j]*h)["nll"] -
                        evaluate(par-np.eye(24)[j]*h)["nll"])/(2*h) for j in range(24)])
        results.append(dict(step=h, gradient=fd, max_difference=float(np.max(abs(fd-analytic)))))
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    out = parser.parse_args().output_dir.resolve()
    out.mkdir(parents=True, exist_ok=False)
    assert sha(PREVIOUS / "summary.json") == SUMMARY_HASH
    old = json.loads((PREVIOUS / "summary.json").read_text())
    for name, expected in old["artifact_sha256"].items():
        assert sha(PREVIOUS / name) == expected, name
    original = json.loads((PREVIOUS / "input.json").read_text())
    cases = json.loads((PREVIOUS / "python.json").read_text())
    starts = {f"from_q{q}": np.asarray(cases[f"q{q}_from_q31"]["coordinates"])
              for q in (31, 61, 121, 181)}
    r_script = ROOT / "validation/mml_pcm_quadrature_diagnosis.R"
    hashes = {str(p.relative_to(ROOT)): sha(p) for p in (Path(__file__), r_script)}
    spec = dict(classification="OBSERVED_DEVELOPMENT_ONLY", scientific_inference_ready=False,
        qualification_eligible=False, data=original["data"], orders=original["orders"],
        integration=original["integration"], starts=starts, options=OPTIONS, limits=LIMITS,
        previous_summary_sha256=SUMMARY_HASH, source_sha256=hashes,
        design="Four retained GH starts, one continuous refit and one restart each; no retries",
        sigma_bounds=[.05, 10], mean=0, slope=1,
        runtime=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__))
    dump(out / "input.json", spec)
    problem = PCMIntegral(original["data"])
    checks = {}
    for name, par in starts.items():
        result = problem.evaluate(par)
        r = json.loads((PREVIOUS / f"{name[5:]}_from_q31_r.json").read_text())
        # name 'from_q31' -> retained case 'q31_from_q31'.
        differences = {k: float(np.max(np.abs(np.asarray(result[k])-r["continuous"][1][k])))
                       for k in ("nll", "eap", "sd")}
        checks[name] = dict(result=result, r_differences=differences)
        assert differences["nll"] < LIMITS["r_nll"]
        assert max(differences["eap"], differences["sd"]) < LIMITS["r_moments"]
    checks["gradient_fd"] = fd_check(problem.evaluate, starts["from_q61"])
    assert max(c["max_difference"] for c in checks["gradient_fd"]) < LIMITS["gradient_fd"]
    dump(out / "preflight.json", checks)
    print("Independent R replay and full-coordinate gradient checks passed", flush=True)

    fits = {}
    def objective(par):
        value = problem.evaluate(par)
        return value["nll"], value["gradient"]
    bounds = [(None, None)] * 23 + [(np.log(.05), np.log(10))]
    for name, start in starts.items():
        runs = []
        par = start
        for _ in range(2):
            begin = time.monotonic()
            fit = minimize(objective, par, jac=True, method="L-BFGS-B", bounds=bounds, options=OPTIONS)
            runs.append(dict(success=bool(fit.success), status=int(fit.status), message=str(fit.message),
                nit=int(fit.nit), nfev=int(fit.nfev), coordinates=fit.x, nll=float(fit.fun),
                gradient=fit.jac, elapsed=time.monotonic()-begin))
            par = fit.x
        primary, tight = [problem.evaluate(par, **s) for s in spec["integration"]]
        fits[name] = dict(coordinates=par, runs=runs, primary=primary, tight=tight,
            coordinate_movement=float(np.max(abs(par-start))), sigma_movement=float(np.exp(par[-1])-np.exp(start[-1])),
            continuous_nll_improvement=checks[name]["result"]["nll"]-tight["nll"],
            refinement={k: float(np.max(abs(np.asarray(primary[k])-tight[k])))
                        for k in ("nll", "gradient", "eap", "sd")})
        dump(out / f"{name}.json", fits[name])
        print(name, "sigma", tight["sigma"], "nll", tight["nll"],
              "gradient", np.max(abs(tight["gradient"])), "success", [r["success"] for r in runs], flush=True)
    dump(out / "python.json", fits)
    subprocess.run(["Rscript", str(r_script), str(out / "input.json"),
                    str(out / "python.json"), str(out)], check=True)
    par = fits["from_q181"]["coordinates"]
    tight_eval = lambda p: problem.evaluate(p, **spec["integration"][1])
    final_fd = fd_check(tight_eval, par)
    hessians = []
    for h in (1e-4, 3e-5):
        raw = np.column_stack([(tight_eval(par+np.eye(24)[j]*h)["gradient"]-
                                tight_eval(par-np.eye(24)[j]*h)["gradient"])/(2*h) for j in range(24)])
        hessians.append(dict(step=h, matrix=raw, max_asymmetry=float(np.max(abs(raw-raw.T))),
                             eigenvalues=np.linalg.eigvalsh((raw+raw.T)/2)))
    comparisons = {}
    for name, fit in fits.items():
        r = json.loads((out / f"{name}_r.json").read_text())["continuous"][1]
        comparisons[name] = {k: float(np.max(abs(np.asarray(fit["tight"][k])-r[k])))
                             for k in ("nll", "eap", "sd")}
    passed = (all(r["success"] for f in fits.values() for r in f["runs"]) and
        all(np.max(abs(f["tight"]["gradient"])) < LIMITS["gradient"] and
            max(f["refinement"].values()) < LIMITS["refinement"] and
            max(f["tight"][k] for k in ("numeric_relative_mass_error_sum", "tail_relative_mass_bound_sum"))
                < LIMITS["relative_mass_error"] for f in fits.values()) and
        all(c["nll"] < LIMITS["r_nll"] and max(c["eap"], c["sd"]) < LIMITS["r_moments"]
            for c in comparisons.values()) and
        max(c["max_difference"] for c in final_fd) < LIMITS["gradient_fd"] and
        all(h["eigenvalues"][0] > 0 for h in hessians))
    assert hashes == {p: sha(ROOT / p) for p in hashes}
    assert sha(PREVIOUS / "summary.json") == SUMMARY_HASH
    dump(out / "summary.json", dict(classification=spec["classification"],
        scientific_inference_ready=False, qualification_eligible=False, implementation_checks_pass=passed,
        fits=fits, independent_r_differences=comparisons, final_gradient_fd=final_fd, hessians=hessians,
        cross_start_max_coordinate_distance=float(np.max(np.ptp([f["coordinates"] for f in fits.values()], axis=0))),
        limits=LIMITS, source_sha256=hashes, previous_summary_sha256=SUMMARY_HASH,
        artifact_sha256={p.name: sha(p) for p in sorted(out.iterdir()) if p.is_file()},
        limitations=["One observed stress fixture, two persons forced extreme; no recovery claim",
            "Bounded adaptive integration with estimated numerical error, not exact integration",
            "Four starts and local positive Hessian do not prove a global optimum",
            "No SE, confidence interval, coverage or scientific-inference qualification"]))
    print("Implementation checks:", passed, flush=True)
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
