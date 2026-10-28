"""Separate bounded robustness study; reuse the retained fixed-SD fit engine."""
import argparse
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
import pandas as pd
from scipy.special import softmax
from scipy.stats import norm

from validation import medium_writing_check as base
from validation.mml_rasch_block_reference import BlockRaschIntegral


def cells(plan):
    return [(design, mechanism) for design in plan["designs"] for mechanism in plan["mechanisms"]]


def case_id(design, mechanism):
    return f"N200_{design}_{mechanism}"


def generate(design, mechanism, replicate, seed_entropy):
    """Conditional balanced groups; identical latent pool and uniforms across cells."""
    n = 200
    rng = np.random.default_rng(np.random.SeedSequence([int(seed_entropy), int(replicate)]))
    z, task_shock = rng.normal(size=n), rng.normal(size=(n, 3))
    permutation = rng.permutation(n)
    uniform = rng.random((n, 2, 3))
    group = np.where(np.arange(n) % 2 == 0, 1., -1.)
    original_pairs = base.assignment(n, design)
    pairs = original_pairs if mechanism == "mixture_selected" else original_pairs[permutation]
    if mechanism == "normal_independent":
        person, interaction = z, np.zeros((n, 3))
    elif mechanism in ("mixture_independent", "mixture_selected"):
        person, interaction = .8 * z + .6 * group, np.zeros((n, 3))
    elif mechanism == "shared_task_total_variance":
        person, interaction = .8 * z, .6 * task_shock
    else:
        raise ValueError("Unknown prespecified generating mechanism")
    eta = person[:, None] + interaction - np.array([-.3, 0., .3])
    probability = softmax(eta[:, :, None] * np.arange(4) - np.array([0., -1., -1., 0.]), axis=2)
    score = (uniform[:, :, :, None] > probability.cumsum(axis=2)[:, None, :, :]).sum(axis=3)
    data = pd.DataFrame(dict(Person=np.repeat([f"P{p:03d}" for p in range(n)], 6),
        Rater=[f"R{r}" for r in np.repeat(pairs, 3, axis=1).ravel()],
        Task=np.tile(["T0", "T1", "T2"], 2*n), Score=score.ravel()))
    workload = np.bincount(pairs.ravel(), minlength=10)
    positive = np.bincount(pairs.ravel(), weights=np.repeat(group > 0, 2), minlength=10).astype(int)
    assert np.all(workload == 40) and len(data) == 1200
    assert data.groupby(["Rater", "Task"]).size().eq(40).all()
    assert np.isin(score, np.arange(4)).all()
    if mechanism == "mixture_selected":
        expected = np.full(10, 20) if design == "rotating" else np.r_[np.full(5, 38), np.full(5, 2)]
        np.testing.assert_array_equal(positive, expected)
    truth = dict(primary=0., secondary=0., task=.6, seed_components=[int(seed_entropy), int(replicate)],
        mechanism=mechanism, group=group, z=z, task_shock=task_shock, permutation=permutation,
        person_effect=person, person_task_effect=interaction, generating_latent_predictor=eta,
        positive_group_persons_per_rater=positive, persons_per_rater=workload,
        fixed_balanced_group_sampling=True, scientific_inference_ready=False)
    return data, truth


def task_contrast(point, covariance):
    """T2 minus T0: task coordinates 9..11 are unconstrained in this model."""
    c = np.zeros(14)
    c[9], c[11] = -1., 1.
    result = dict(estimate=float(c @ point), truth=.6)
    if covariance is not None:
        for method, variance in (("full", c @ covariance @ c),
                                 ("ignored_covariance", (c*c) @ np.diag(covariance))):
            se = float(np.sqrt(variance))
            lower, upper = result["estimate"] + norm.ppf(.975) * se * np.array([-1., 1.])
            result[method] = dict(se=se, lower=float(lower), upper=float(upper),
                                 covered=bool(lower <= .6 <= upper), excludes_zero=bool(lower > 0 or upper < 0))
    return result


def fit(data, truth, plan):
    result = base.fit(data, truth, plan["quadrature_order"], plan)
    point = result.get("point")
    if point is not None and np.isfinite(point).all():
        result.setdefault("contrasts", {})["task"] = task_contrast(np.asarray(point), result.get("covariance"))
    return result


def retained_input(folder, design, mechanism, replicate, entropy):
    manifest_path = folder / "input_manifest.json"
    if manifest_path.exists():
        for name, expected in json.loads(manifest_path.read_text()).items():
            if base.digest(folder / name) != expected:
                raise ValueError("Retained input changed: " + str(folder / name))
        return pd.read_csv(folder / "data.csv", dtype={"Person": str}), json.loads((folder / "truth.json").read_text())
    data, truth = generate(design, mechanism, replicate, entropy)
    folder.mkdir(parents=True, exist_ok=True)
    with (folder / "data.csv").open("x") as handle:
        data.to_csv(handle, index=False)
    base.save(folder / "truth.json", truth)
    base.save(manifest_path, {name: base.digest(folder/name) for name in ("data.csv", "truth.json")})
    return data, truth


def run_fits(out, plan):
    started = time.monotonic()
    for design, mechanism in cells(plan):
        for rep in range(plan["replicates_per_cell"]):
            remaining = plan["resource_limits"]["fit_stage_wall_seconds"] - (time.monotonic() - started)
            if remaining <= 0:
                return
            folder = out / "cases" / case_id(design, mechanism) / f"rep_{rep:03d}"
            data, truth = retained_input(folder, design, mechanism, rep, plan["seed_entropy"])
            result = base.attempt(folder / "q81_grouped.json",
                min(remaining, plan["resource_limits"]["fit_timeout_seconds"]), lambda: fit(data, truth, plan))
            print("fit", case_id(design, mechanism), rep, result["status"], result.get("seconds"), flush=True)


def continuous_case(folder, plan):
    base.RaschIntegral = BlockRaschIntegral
    result = base.continuous_case(folder, 200, plan)
    if result["status"] != "continuous_subset_checked":
        return result
    exact = result["evaluations"][-1]
    information = np.asarray(exact["information"])[:-1, :-1]
    covariance = np.linalg.inv(information) if np.linalg.eigvalsh((information+information.T)/2)[0] > 0 else None
    task = task_contrast(np.asarray(result["point"]), covariance)
    saved = json.loads((folder / "q81_grouped.json").read_text())
    target = saved.get("contrasts", {}).get("task", {}).get("full")
    error = abs(task["full"]["se"] / target["se"] - 1) if "full" in task and target else np.inf
    result["contrasts"]["task"] = task
    result["errors"]["task_se_relative"] = error
    result["errors"]["target_se_relative"] = max(result["errors"]["target_se_relative"], error)
    result["checks"]["target_se_relative"] = result["errors"]["target_se_relative"] <= plan["continuous_reference"]["accuracy_limits"]["target_se_relative"]
    result["all_checks_pass"] = all(result["checks"].values())
    result["interpretation"] = "Numerical integration of the fitted working model only; not generating-model or interval qualification."
    return result


def run_references(out, plan):
    started = time.monotonic()
    for design, mechanism in cells(plan):
        remaining = plan["resource_limits"]["reference_stage_wall_seconds"] - (time.monotonic() - started)
        if remaining <= 0:
            return
        name = case_id(design, mechanism)
        folder = out / "cases" / name / "rep_000"
        result = base.attempt(out / "references" / (name + ".json"),
            min(remaining, plan["continuous_reference"]["per_case_timeout_seconds"]),
            lambda: continuous_case(folder, plan))
        print("reference", name, result["status"], result.get("all_checks_pass"), flush=True)


def analyze(out, plan):
    planned = plan["replicates_per_cell"]
    records, rows, pairs = {}, [], []
    for design, mechanism in cells(plan):
        records[design, mechanism], statuses = {}, {}
        for rep in range(planned):
            path = out / "cases" / case_id(design, mechanism) / f"rep_{rep:03d}" / "q81_grouped.json"
            value = json.loads(path.read_text()) if path.exists() else dict(status="unresolved", numerical_pass=False)
            records[design, mechanism][rep] = value
            statuses[value["status"]] = statuses.get(value["status"], 0) + 1
        valid = [value for value in records[design, mechanism].values() if value.get("numerical_pass")]
        for target, target_spec in plan["targets"].items():
            estimates = np.array([v["contrasts"][target]["estimate"] for v in valid])
            errors = estimates - target_spec["truth"]
            row = dict(design=design, mechanism=mechanism, target=target, truth=target_spec["truth"],
                planned=planned, eligible=len(valid), statuses=statuses,
                bias=float(errors.mean()) if len(valid) else None,
                bias_mcse=float(errors.std(ddof=1)/np.sqrt(len(valid))) if len(valid)>1 else None,
                rmse=float(np.sqrt(np.mean(errors**2))) if len(valid) else None,
                empirical_sd=float(estimates.std(ddof=1)) if len(valid)>1 else None)
            for method in ("full", "ignored_covariance"):
                values = [v["contrasts"][target][method] for v in valid]
                stats = dict(generating_truth_containment=base.binary_summary(sum(v["covered"] for v in values), len(values), planned),
                             mean_se=float(np.mean([v["se"] for v in values])) if values else None)
                if target != "task":
                    stats["generating_null_zero_exclusion"] = base.binary_summary(sum(v["excludes_zero"] for v in values), len(values), planned)
                row[method] = stats
            rows.append(row)
    for design in plan["designs"]:
        for target in plan["targets"]:
            differences, pair_rows = [], []
            for rep in range(planned):
                a, b = records[design, "mixture_selected"][rep], records[design, "mixture_independent"][rep]
                eligible = bool(a.get("numerical_pass") and b.get("numerical_pass"))
                difference = a["contrasts"][target]["estimate"] - b["contrasts"][target]["estimate"] if eligible else None
                pair_rows.append(dict(replicate=rep, both_eligible=eligible, estimate_difference=difference))
                if eligible:
                    differences.append(difference)
            pairs.append(dict(design=design, target=target, planned=planned, both_eligible=len(differences),
                selected_minus_independent_mean=float(np.mean(differences)) if differences else None,
                paired_mcse=float(np.std(differences,ddof=1)/np.sqrt(len(differences))) if len(differences)>1 else None,
                rows=pair_rows))
    result = dict(classification=plan["classification"], scientific_inference_ready=False, cells=rows,
                  paired_selected_minus_independent=pairs, replication_unit="independent seed; scenarios within seed are paired")
    base.save(out / "analysis.json", result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    action = parser.add_mutually_exclusive_group(required=True)
    for name in ("run", "references", "analyze"):
        action.add_argument("--"+name, action="store_true")
    args = parser.parse_args()
    out = args.output.resolve()
    seal = json.loads((out / "study_seal.json").read_text())
    plan = json.loads((out / "protocol.json").read_text())
    manifest = json.loads((out / "numerical_source_manifest.json").read_text())
    def verify():
        if base.digest(out/"protocol.json") != seal["protocol_sha256"]:
            raise ValueError("Sealed protocol changed")
        if base.digest(out/"numerical_source_manifest.json") != seal["numerical_manifest_sha256"]:
            raise ValueError("Sealed source manifest changed")
        if Path(base.app.__file__).resolve() != out/"frozen/streamlit_app.py":
            raise ValueError("Run the separately retained frozen study source")
        for relative, expected in manifest.items():
            if base.digest(out/"frozen"/relative) != expected:
                raise ValueError("Retained source changed: "+relative)
    verify()
    if args.run:
        run_fits(out,plan)
    elif args.references:
        run_references(out,plan)
    else:
        analyze(out,plan)
    verify()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
