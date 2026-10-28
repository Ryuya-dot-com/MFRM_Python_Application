"""Bounded, separately frozen fixed-SD writing study; no app qualification.

Run this file from the retained frozen source directory. The protocol is an
input, never rewritten here. Every input, attempt and returned result is kept.
"""
import argparse
from contextlib import contextmanager
import hashlib
import json
from pathlib import Path
import signal
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import softmax
from scipy.stats import norm

import streamlit_app as app
from mfrm_app.mml_joint import GroupedMmlLikelihood, normal_rule, OPTIONS
from validation.mml_free_sd_stationarity_adapter import prepare_app_free_sd_problem_from_data
from validation.mml_rasch_integral import RaschIntegral


PROTOCOL_SHA256 = "93cc5367ab099c4e55330d25a259fe8ba3c2f581c33b10870e31e4eabe178f27"


def plain(value):
    if isinstance(value, dict):
        return {str(k): plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [plain(v) for v in value]
    if isinstance(value, np.generic):
        return plain(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x") as handle:
        json.dump(plain(value), handle, indent=2, allow_nan=False)
        handle.write("\n")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class StudyTimeout(TimeoutError):
    pass


@contextmanager
def deadline(seconds):
    def timeout(signum, frame):
        raise StudyTimeout(f"Prespecified time limit ({seconds:g} seconds) exceeded")
    previous = signal.signal(signal.SIGALRM, timeout)
    signal.setitimer(signal.ITIMER_REAL, seconds)
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous)


def case_id(n, design, delta):
    return f"N{n}_{design}_d{delta:.1f}"


def cases(plan):
    return [(n, design, delta) for n in plan["persons"]
            for design in plan["designs"] for delta in plan["deltas"]]


def assignment(n, design):
    p = np.arange(n)
    if design == "rotating":
        pairs = np.c_[p % 5, 5 + (p // 5 + p) % 5]
    elif design == "weak_bridge":
        h, block = p // 2, p % 2
        a = h % 5
        bridge = h % 10 == (h // 10) % 5
        pairs = np.c_[5 * block + a, 5 * np.where(bridge, 1 - block, block) + (a + 1) % 5]
    else:
        raise ValueError("Unknown assignment design")
    if n not in (200, 500) or not np.all(np.bincount(pairs.ravel(), minlength=10) == 2 * n // 10):
        raise ValueError("Assignment is outside the prespecified exactly balanced sizes")
    if np.any(pairs[:, 0] == pairs[:, 1]):
        raise ValueError("A person was assigned the same rater twice")
    expected = n if design == "rotating" else n // 10
    if np.sum(pairs[:, 0] // 5 != pairs[:, 1] // 5) != expected:
        raise ValueError("Cross-block assignment count differs from protocol")
    return pairs


def generate(n, design, delta, replicate, seed_entropy):
    pairs = assignment(n, design)
    seed = [int(seed_entropy), n, ("rotating", "weak_bridge").index(design),
            (0., .4).index(float(delta)), int(replicate)]
    rng = np.random.default_rng(np.random.SeedSequence(seed))
    theta = rng.normal(size=n)  # The realized sample is deliberately not recentered.
    rater = np.r_[np.full(5, delta / 2), np.full(5, -delta / 2)]
    task = np.array([-.3, 0., .3])
    person = np.repeat(np.arange(n), 6)
    r = np.repeat(pairs, 3, axis=1).ravel()
    t = np.tile(np.arange(3), 2 * n)
    eta = theta[person] + rater[r] - task[t]
    # Independent direct RSM generation: cumulative steps [0,-1,-1,0].
    probability = softmax(eta[:, None] * np.arange(4) - np.array([0., -1., -1., 0.]), axis=1)
    score = (rng.random(len(person))[:, None] > probability.cumsum(axis=1)).sum(axis=1)
    data = pd.DataFrame(dict(Person=[f"P{v:03d}" for v in person],
        Rater=[f"R{v}" for v in r], Task=[f"T{v}" for v in t], Score=score))
    assert len(data) == 6 * n and data.groupby(["Rater", "Task"]).size().eq(2 * n // 10).all()
    truth = dict(seed_components=seed, theta=theta, rater=rater, task=task,
                 steps=[-1., 0., 1.], primary=delta, secondary=0., scientific_inference_ready=False)
    return data, truth


def prepare(data, q):
    prep, problem = prepare_app_free_sd_problem_from_data(data, person_col="Person",
        facet_cols=["Rater", "Task"], score_col="Score", rating_min=0, rating_max=3,
        model="RSM", noncenter_facet="Task", positive_facets=("Rater",), quadrature_points=q)
    problem.config["estimate_population_sd"] = False
    assert problem.sizes == app.build_param_sizes(problem.config)
    assert len(problem.structural_start) == 14
    assert prep["levels"]["Rater"] == [f"R{r}" for r in range(10)]
    assert prep["levels"]["Task"] == [f"T{t}" for t in range(3)]
    return prep, problem


def contrast_rows():
    expansion = np.zeros((10, 14))
    expansion[:9, :9] = np.eye(9)
    expansion[9, :9] = -1.
    weights = {"primary": np.r_[np.full(5, .2), np.full(5, -.2)],
               "secondary": np.r_[1., -1., np.zeros(8)]}
    return expansion, weights


def contrast_results(point, covariance, truth):
    expansion, weights = contrast_rows()
    result = {}
    for name, w in weights.items():
        c = w @ expansion
        row = dict(estimate=float(c @ point), truth=float(truth[name]))
        if covariance is not None:
            rater_cov = expansion @ covariance @ expansion.T
            full = float(np.sqrt(c @ covariance @ c))
            ignored = float(np.sqrt((w * w) @ np.diag(rater_cov)))
            for label, se in (("full", full), ("ignored_covariance", ignored)):
                lower, upper = row["estimate"] + np.array([-1., 1.]) * norm.ppf(.975) * se
                row[label] = dict(se=se, lower=float(lower), upper=float(upper),
                    covered=bool(lower <= truth[name] <= upper), excludes_zero=bool(lower > 0 or upper < 0))
        result[name] = row
    return result


def fit(data, truth, q, plan, engine="grouped"):
    if engine not in ("grouped", "row_native"):
        raise ValueError("Unknown likelihood engine")
    started = time.perf_counter()
    phases = {}
    result = dict(engine=engine, q=q, numerical_pass=False, scientific_inference_ready=False,
                  point=None, covariance=None, runs=[])
    phase = "prepare"
    try:
        before = time.perf_counter()
        prep, problem = prepare(data, q)
        phases[phase] = time.perf_counter() - before
        phase = "construct"
        before = time.perf_counter()
        grouped = GroupedMmlLikelihood(prep["data"], problem.config)
        nodes, weights = normal_rule(q)
        phases[phase] = time.perf_counter() - before
        phase = "recession_audit"
        before = time.perf_counter()
        audit = grouped.structural_recession_audit()
        phases[phase] = time.perf_counter() - before
        result["recession"] = audit
        if audit["status"] == "unavailable" or audit["no_finite_mle_certified"]:
            result["status"] = "prefit_" + audit["status"]
            return result
        if audit["structural_design_rank"] != audit["linear_dimension"]:
            result["status"] = "prefit_structural_rank_deficient"
            return result
        quad = dict(nodes=nodes, weights=weights, sd=1.)
        calls, evaluation_seconds = 0, 0.
        def value_gradient(point):
            nonlocal calls, evaluation_seconds
            before = time.perf_counter()
            try:
                if engine == "row_native":
                    return app.mfrm_loglik_mml_value_grad(point, problem.idx, problem.config, problem.sizes, quad)
                evaluated = grouped.evaluate(np.r_[point, 0.], nodes, weights, moments=False)
                return evaluated["nll"], evaluated["gradient"][:-1]
            finally:
                calls += 1
                evaluation_seconds += time.perf_counter() - before
        point = np.array(problem.structural_start)
        assert list(OPTIONS) == plan["estimator"]["mandatory_stages"]
        phase = "optimization"
        before = time.perf_counter()
        for options in OPTIONS:
            initial = point.copy()
            opt = minimize(value_gradient, initial, jac=True, method="SLSQP", options=options)
            point = np.asarray(opt.x)
            result["point"] = point.copy()
            result["runs"].append(dict(initial=initial, point=point, nll=float(opt.fun),
                gradient=opt.jac, success=bool(opt.success), status=int(opt.status),
                message=str(opt.message), iterations=int(opt.nit), evaluations=int(opt.nfev)))
            if not np.isfinite(point).all():
                raise FloatingPointError("Optimizer returned nonfinite coordinates")
        phases[phase] = time.perf_counter() - before
        phases["value_gradient"] = dict(calls=calls, seconds=evaluation_seconds)
        phase = "information_and_moments"
        before = time.perf_counter()
        evaluated = grouped.evaluate(np.r_[point, 0.], nodes, weights, information=True)
        # Fixed SD: discard its row/column BEFORE inversion; no nuisance-SD covariance.
        raw = evaluated["information"][:-1, :-1]
        symmetry = float(np.max(abs(raw - raw.T)))
        information = (raw + raw.T) / 2
        eigen = np.linalg.eigvalsh(information)
        covariance = np.linalg.inv(information) if eigen[0] > 0 else None
        gradient = evaluated["gradient"][:-1]
        correction = float(np.max(abs(covariance @ gradient))) if covariance is not None else np.inf
        condition = float(eigen[-1] / eigen[0]) if eigen[0] > 0 else np.inf
        phases[phase] = time.perf_counter() - before
        phase = "native_same_point_check"
        before = time.perf_counter()
        native_nll, native_gradient = app.mfrm_loglik_mml_value_grad(point, problem.idx, problem.config, problem.sizes, quad)
        parity = dict(nll_abs=float(abs(native_nll - evaluated["nll"])),
                      gradient_abs=float(np.max(abs(native_gradient - gradient))))
        phases[phase] = time.perf_counter() - before
        limits = plan["finite_numerical_limits"]
        restart = float(np.max(abs(result["runs"][1]["point"] - result["runs"][0]["point"])))
        constraints = float(problem.constraint_residual(point))
        checks = dict(optimizer_terminated=all(run["success"] for run in result["runs"]),
            score=float(np.max(abs(gradient))) < limits["score_supnorm"],
            newton_correction=correction < limits["newton_correction_supnorm"],
            restart=restart < limits["restart_difference_supnorm"], constraints=constraints <= limits["constraints"],
            positive_definite=covariance is not None, symmetry=symmetry < limits["symmetry"],
            condition=condition < limits["condition"],
            objective_nonworsening=result["runs"][1]["nll"] - result["runs"][0]["nll"] <= limits["objective_worsening"])
        eligible = all(checks.values())
        result.update(status="finite_rule_stationary" if eligible else "numerical_checks_failed",
            numerical_pass=eligible, point=point, information=information, information_eigenvalues=eigen,
            covariance=covariance if eligible else None, checks=checks, parity=parity,
            gradient=gradient, gradient_supnorm=float(np.max(abs(gradient))), newton_correction=correction,
            condition_number=condition, symmetry_residual=symmetry, restart_difference=restart,
            constraint_residual=constraints, nll=float(evaluated["nll"]), eap=evaluated["eap"],
            posterior_sd=evaluated["sd"], person_labels=grouped.person_labels,
            contrasts=contrast_results(point, covariance if eligible else None, truth))
    except Exception as exc:
        result.update(status="timeout" if isinstance(exc, StudyTimeout) else "exception",
                      error=dict(type=type(exc).__name__, message=str(exc), phase=phase))
    finally:
        result["seconds"] = time.perf_counter() - started
        result["phases"] = phases
        result["reason_code"] = "stat.medium_writing." + result.get("status", "incomplete")
    return result


def retained_input(folder, n, design, delta, rep, entropy):
    if (folder / "input_manifest.json").exists():
        manifest = json.loads((folder / "input_manifest.json").read_text())
        for name, expected in manifest.items():
            if digest(folder / name) != expected:
                raise ValueError("Retained input changed: " + str(folder / name))
        return pd.read_csv(folder / "data.csv", dtype={"Person": str}), json.loads((folder / "truth.json").read_text())
    data, truth = generate(n, design, delta, rep, entropy)
    folder.mkdir(parents=True, exist_ok=True)
    with (folder / "data.csv").open("x") as handle:
        data.to_csv(handle, index=False)
    save(folder / "truth.json", truth)
    save(folder / "input_manifest.json", {name: digest(folder / name) for name in ("data.csv", "truth.json")})
    return data, truth


def attempt(path, seconds, function):
    if path.exists():
        return json.loads(path.read_text())
    marker = path.with_suffix(".started.json")
    if marker.exists():
        result = dict(status="interrupted_attempt_not_retried", numerical_pass=False,
                      scientific_inference_ready=False, reason_code="stat.medium_writing.interrupted_attempt_not_retried")
        save(path, result)
        return result
    save(marker, dict(started_unix=time.time(), seconds_limit=seconds))
    started = time.perf_counter()
    try:
        with deadline(seconds):
            result = function()
    except Exception as exc:
        result = dict(status="timeout" if isinstance(exc, StudyTimeout) else "exception", numerical_pass=False,
                      error=dict(type=type(exc).__name__, message=str(exc)))
    result.setdefault("scientific_inference_ready", False)
    result.setdefault("reason_code", "stat.medium_writing." + result["status"])
    result.setdefault("seconds", time.perf_counter() - started)
    save(path, result)
    return result


def run_fits(out, plan, *, pilot=False):
    start = time.monotonic()
    budget = plan["resource_limits"]["fit_stage_wall_seconds"]
    selected = [(n, "weak_bridge", .4) for n in plan["persons"]] if pilot else cases(plan)
    for n, design, delta in selected:
        for rep in range(1 if pilot else plan["replicates_per_cell"]):
            if time.monotonic() - start >= budget:
                return
            folder = out / ("pilot" if pilot else "cases") / case_id(n, design, delta) / f"rep_{rep:03d}"
            entropy = plan["development_seed_entropy" if pilot else "seed_entropy"]
            data, truth = retained_input(folder, n, design, delta, rep, entropy)
            for q in plan["quadrature_orders"]:
                remaining = budget - (time.monotonic() - start)
                if remaining <= 0:
                    return
                result = attempt(folder / f"q{q}_grouped.json",
                    min(plan["resource_limits"]["fit_timeout_seconds"], remaining),
                    lambda: fit(data, truth, q, plan))
                print(case_id(n, design, delta), rep, q, result["status"], result.get("seconds"), flush=True)


def continuous_case(folder, n, plan):
    data = pd.read_csv(folder / "data.csv", dtype={"Person": str})
    truth = json.loads((folder / "truth.json").read_text())
    saved = json.loads((folder / "q81_grouped.json").read_text())
    if saved.get("point") is None:
        return dict(status="source_point_unavailable", numerical_pass=False)
    prep, problem = prepare(data, 81)
    reference = RaschIntegral(prep["data"], problem.config)
    point = np.asarray(saved["point"])
    spec = plan["continuous_reference"]
    settings = [dict(bound=spec["bounds"], rel_tol=spec["relative_tolerance"], abs_tol=spec["absolute_tolerance"])]
    if n == 200:
        settings.append(dict(bound=14., rel_tol=1e-12, abs_tol=1e-14))
    values = [reference.continuous(np.r_[point, 0.], information=True, **setting) for setting in settings]
    exact = values[-1]
    information = exact["information"][:-1, :-1]
    eigen = np.linalg.eigvalsh((information + information.T) / 2)
    covariance = np.linalg.inv(information) if eigen[0] > 0 else None
    contrasts = contrast_results(point, covariance, truth)
    errors = dict(nll_per_person=abs(exact["nll"] - saved["nll"]) / n,
                  structural_gradient_abs=float(np.max(abs(exact["gradient"][:-1] - saved["gradient"]))),
                  continuous_score_supnorm=float(np.max(abs(exact["gradient"][:-1]))))
    se_errors = [abs(contrasts[name]["full"]["se"] / saved["contrasts"][name]["full"]["se"] - 1)
                 for name in contrasts if "full" in contrasts[name] and "full" in saved["contrasts"][name]]
    errors["target_se_relative"] = max(se_errors) if len(se_errors) == 2 else np.inf
    if len(values) == 2:
        errors.update(reference_refinement_nll_abs=abs(values[0]["nll"] - exact["nll"]),
            reference_refinement_info_abs=float(np.max(abs(values[0]["information"] - exact["information"]))))
    checks = {key: errors[key] <= limit for key, limit in spec["accuracy_limits"].items() if key in errors}
    checks["continuous_score"] = errors["continuous_score_supnorm"] <= spec["accuracy_limits"]["structural_gradient_abs"]
    checks["positive_definite"] = covariance is not None
    return dict(status="continuous_subset_checked", point=point, source_numerical_pass=saved["numerical_pass"],
                settings=settings, evaluations=values, contrasts=contrasts, errors=errors, checks=checks,
                all_checks_pass=all(checks.values()), refinement_performed=len(values) == 2)


def run_references(out, plan):
    start = time.monotonic()
    for n, design, delta in cases(plan):
        remaining = plan["resource_limits"]["reference_stage_wall_seconds"] - (time.monotonic() - start)
        if remaining <= 0:
            return
        name = case_id(n, design, delta)
        folder = out / "cases" / name / "rep_000"
        value = attempt(out / "references" / (name + ".json"),
            min(remaining, plan["continuous_reference"]["per_case_timeout_seconds"]),
            lambda: continuous_case(folder, n, plan))
        print("reference", name, value["status"], value.get("all_checks_pass"), flush=True)


def matched_comparison(a, b, n, limits):
    checks = dict(both_numerical_pass=bool(a.get("numerical_pass") and b.get("numerical_pass")))
    errors = {}
    if checks["both_numerical_pass"]:
        errors.update(point_abs=float(np.max(abs(np.asarray(a["point"]) - b["point"]))),
            nll_per_person=abs(a["nll"] - b["nll"]) / n,
            eap_abs=float(np.max(abs(np.asarray(a["eap"]) - b["eap"]))),
            posterior_sd_abs=float(np.max(abs(np.asarray(a["posterior_sd"]) - b["posterior_sd"]))))
        ac, bc = np.asarray(a["covariance"]), np.asarray(b["covariance"])
        errors["covariance_scaled"] = float(np.max(abs(ac - bc) / np.sqrt(np.outer(np.diag(ac), np.diag(ac)))))
    checks.update({key: key in errors and errors[key] <= limit for key, limit in limits.items()})
    return dict(errors=errors, checks=checks, all_checks_pass=all(checks.values()))


def run_benchmark(out, plan):
    start = time.monotonic()
    spec = plan["matched_benchmark"]
    for n in plan["persons"]:
        for design in plan["designs"]:
            name = case_id(n, design, .4)
            folder = out / "cases" / name / "rep_000"
            data, truth = retained_input(folder, n, design, .4, 0, plan["seed_entropy"])
            for repeat in range(spec["pairs_per_case"]):
                pair = {}
                order = ("grouped", "row_native") if repeat % 2 == 0 else ("row_native", "grouped")
                for engine in order:
                    remaining = plan["resource_limits"]["fit_stage_wall_seconds"] - (time.monotonic() - start)
                    if remaining <= 0:
                        return
                    path = out / "benchmarks" / f"{name}_pair{repeat}_{engine}.json"
                    pair[engine] = attempt(path, min(remaining, plan["resource_limits"]["fit_timeout_seconds"]),
                                          lambda: fit(data, truth, spec["q"], plan, engine))
                    print("benchmark", name, repeat, engine, pair[engine]["status"], pair[engine].get("seconds"), flush=True)
                path = out / "benchmarks" / f"{name}_pair{repeat}_comparison.json"
                if not path.exists():
                    save(path, matched_comparison(pair["grouped"], pair["row_native"], n, spec["acceptance"]))


def binary_summary(successes, eligible, planned):
    if eligible:
        p = successes / eligible
        z = norm.ppf(.975)
        center = (p + z * z / (2 * eligible)) / (1 + z * z / eligible)
        half = z * np.sqrt(p * (1 - p) / eligible + z * z / (4 * eligible * eligible)) / (1 + z * z / eligible)
        mcse, lower, upper = np.sqrt(p * (1 - p) / eligible), center - half, center + half
    else:
        p = mcse = lower = upper = None
    return dict(successes=successes, eligible=eligible, planned=planned, proportion=p, mcse=mcse,
                wilson_lower=lower, wilson_upper=upper, unresolved=planned - eligible,
                planned_lower=successes / planned, planned_upper=(successes + planned - eligible) / planned)


def analyze(out, plan):
    rows, pairs = [], []
    planned = plan["replicates_per_cell"]
    for n, design, delta in cases(plan):
        records = {}
        for q in plan["quadrature_orders"]:
            records[q] = {}
            statuses = {}
            for rep in range(planned):
                path = out / "cases" / case_id(n, design, delta) / f"rep_{rep:03d}" / f"q{q}_grouped.json"
                record = json.loads(path.read_text()) if path.exists() else dict(status="unresolved", numerical_pass=False)
                records[q][rep] = record
                statuses[record["status"]] = statuses.get(record["status"], 0) + 1
            valid = [v for v in records[q].values() if v.get("numerical_pass")]
            for target in ("primary", "secondary"):
                estimates = np.array([v["contrasts"][target]["estimate"] for v in valid])
                truth = delta if target == "primary" else 0.
                errors = estimates - truth
                row = dict(case=case_id(n, design, delta), persons=n, design=design, delta=delta, q=q,
                    target=target, planned=planned, eligible=len(valid), statuses=statuses,
                    bias=float(errors.mean()) if len(valid) else None,
                    bias_mcse=float(errors.std(ddof=1) / np.sqrt(len(valid))) if len(valid) > 1 else None,
                    rmse=float(np.sqrt(np.mean(errors ** 2))) if len(valid) else None,
                    empirical_sd=float(estimates.std(ddof=1)) if len(valid) > 1 else None)
                for method in ("full", "ignored_covariance"):
                    values = [v["contrasts"][target][method] for v in valid]
                    row[method] = dict(coverage=binary_summary(sum(v["covered"] for v in values), len(values), planned),
                        zero_exclusion=binary_summary(sum(v["excludes_zero"] for v in values), len(values), planned),
                        mean_se=float(np.mean([v["se"] for v in values])) if values else None)
                rows.append(row)
        low, high = plan["quadrature_orders"]
        common = [rep for rep in range(planned) if records[low][rep].get("numerical_pass") and records[high][rep].get("numerical_pass")]
        for target in ("primary", "secondary"):
            differences = []
            for rep in common:
                a, b = [records[q][rep]["contrasts"][target] for q in (low, high)]
                differences.append(dict(replicate=rep, estimate_difference=a["estimate"] - b["estimate"],
                    se_difference=a["full"]["se"] - b["full"]["se"],
                    decision_disagrees=a["full"]["excludes_zero"] != b["full"]["excludes_zero"]))
            pairs.append(dict(case=case_id(n, design, delta), target=target, planned=planned,
                              both_eligible=len(common), differences=differences))
    benchmarks = []
    for n in plan["persons"]:
        for design in plan["designs"]:
            name = case_id(n, design, .4)
            times = {engine: [] for engine in ("grouped", "row_native")}
            checks = []
            for rep in range(plan["matched_benchmark"]["pairs_per_case"]):
                for engine in times:
                    path = out / "benchmarks" / f"{name}_pair{rep}_{engine}.json"
                    if path.exists():
                        value = json.loads(path.read_text())
                        if value.get("seconds") is not None:
                            times[engine].append(value["seconds"])
                path = out / "benchmarks" / f"{name}_pair{rep}_comparison.json"
                checks.append(json.loads(path.read_text())["all_checks_pass"] if path.exists() else False)
            medians = {key: float(np.median(value)) if value else None for key, value in times.items()}
            ratio = medians["grouped"] / medians["row_native"] if all(v is not None for v in medians.values()) else None
            benchmarks.append(dict(case=name, seconds=times, median_seconds=medians, grouped_over_native=ratio,
                matched_accuracy_pass=all(checks), reduction_20_percent=all(checks) and ratio is not None and ratio <= .8))
    summary = dict(classification=plan["classification"], scientific_inference_ready=False,
                   cells=rows, paired_quadrature=pairs, benchmarks=benchmarks)
    save(out / "analysis.json", summary)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    action = parser.add_mutually_exclusive_group(required=True)
    for name in ("pilot", "run", "references", "benchmark", "analyze"):
        action.add_argument("--" + name, action="store_true")
    args = parser.parse_args()
    out = args.output.resolve()
    if digest(out / "protocol.json") != PROTOCOL_SHA256:
        raise ValueError("Protocol differs from the prespecified cohort")
    plan = json.loads((out / "protocol.json").read_text())
    manifest_path = out / "numerical_source_manifest.json"
    if digest(manifest_path) != plan["source_manifest_sha256"]:
        raise ValueError("Frozen numerical manifest changed")
    manifest = json.loads(manifest_path.read_text())
    def verify():
        if Path(app.__file__).resolve() != out / "frozen" / "streamlit_app.py":
            raise ValueError("This study must import the frozen application")
        for name, expected in manifest.items():
            if digest(out / "frozen" / name) != expected:
                raise ValueError("Frozen source changed: " + name)
    verify()
    if args.pilot or args.run:
        run_fits(out, plan, pilot=args.pilot)
    elif args.references:
        run_references(out, plan)
    elif args.benchmark:
        run_benchmark(out, plan)
    else:
        analyze(out, plan)
    verify()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
