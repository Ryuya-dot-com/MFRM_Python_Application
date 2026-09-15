#!/usr/bin/env python3
"""Review retained TAM endpoints; keep ConQuest launch failures separate."""
import argparse
import json
from pathlib import Path
import subprocess

import numpy as np
from scipy.special import logsumexp

from mml_pcm_continuous_refit import PCMIntegral, PREVIOUS, ROOT, dump, sha


def fixed_grid(problem, par, theta):
    sigma = np.exp(par[-1])
    step = np.diff(theta)
    assert np.allclose(step, step[0], rtol=0, atol=1e-12)
    offset = problem.design @ par[:-1]
    logits = offset[:, :, None] + problem.k[None, :, None] * theta
    norm = logsumexp(logits, axis=1)
    prob = np.exp(logits-norm[:, None, :])
    ll = (theta[None, :]*problem.total[:, None] +
          (problem.observed@par[:-1])[:, None] - norm.sum(axis=0))
    # TAM's reported deviance uses the density * spacing, without renormalizing
    # the truncated grid's prior mass. Posterior weights are normalized per case.
    logw = -.5*(theta/sigma)**2-np.log(sigma)-.5*np.log(2*np.pi)+np.log(step[0])
    joint = ll+logw
    lm = logsumexp(joint, axis=1)
    post = np.exp(joint-lm[:, None])
    eap = post@theta
    structural = -(problem.observed - post@np.einsum("ikq,ikd->qd", prob, problem.design)).sum(axis=0)
    return dict(nll=float(-lm.sum()), eap=eap,
        sd=np.sqrt(np.sum(post*(theta-eap[:, None])**2, axis=1)),
        gradient=np.r_[structural, np.sum(1-post@(theta/sigma)**2)],
        raw_prior_mass=float(np.exp(logsumexp(logw))))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    args = parser.parse_args()
    out = args.directory.resolve()
    spec = json.loads((out/"input.json").read_text())
    assert not spec["scientific_inference_ready"] and not spec["qualification_eligible"]
    for name, expected in spec["input_sha256"].items():
        assert sha(out/name) == expected, name
    for name, expected in spec["source_sha256"].items():
        assert sha(ROOT/name) == expected, name
    continuous = ROOT/"validation/generated/mml_pcm_continuous_refit_20260913"
    assert sha(continuous/"summary.json") == "d92f7359e826772f7a21b498fa1f44277c612165924c3766cff7c2dabda8630f"
    ref = json.loads((continuous/"from_q181.json").read_text())
    original = json.loads((PREVIOUS/"input.json").read_text())
    assert spec["data"] == original["data"]
    evaluation_spec = {k: original[k] for k in
        ("classification", "scientific_inference_ready", "qualification_eligible", "data", "orders", "integration")}
    hashes = {str(p.relative_to(ROOT)): sha(p) for p in
        (Path(__file__).resolve(), ROOT/"validation/mml_pcm_continuous_refit.py",
         ROOT/"validation/mml_pcm_quadrature_diagnosis.R", ROOT/"validation/mml_pcm_external_refit.R")}
    evaluation_spec.update(source_sha256=hashes, limits=dict(nll=1e-8, moments=1e-8, design=1e-10),
        scope="Same-point arithmetic checks, not equivalence or convergence acceptance margins")
    dump(out/"evaluation_input.json", evaluation_spec)
    problem = PCMIntegral(spec["data"])
    results = {}
    for config in spec["tam"]:
        name = config["id"]
        native = json.loads((out/f"{name}.json").read_text())
        assert native["config"] == config and native["beta"] == [[0]]
        assert len(native["item"]) == 20 and all(i["N"] == 24 for i in native["item"])
        assert [p["pid"] for p in native["person"]] == [f"P{i:02d}" for i in range(1,25)]
        assert [p["score"] for p in native["person"]] == problem.total.tolist()
        byname = {row["item"]: i for i,row in enumerate(native["item"])}
        order = [byname[f"C{c}-raterR{r}"] for r in range(1,5) for c in range(1,6)]
        offset = np.array(native["AXsi"])[order]
        slopes = np.array(native["B"])[order]
        assert slopes.shape == (20,5,1) and np.array_equal(slopes[:,:,0], np.tile(np.arange(5),(20,1)))
        coords, _, rank, _ = np.linalg.lstsq(problem.design.reshape(100,23), offset.ravel(), rcond=None)
        residual = float(np.max(abs(problem.design@coords-offset)))
        assert rank == 23 and residual < 1e-10
        par = np.r_[coords, .5*np.log(native["variance"][0][0])]
        finite = fixed_grid(problem, par, np.array(native["control"]["nodes"]))
        adaptive = problem.evaluate(par, **original["integration"][1])
        native_values = dict(nll=native["deviance"]/2, eap=[p["EAP"] for p in native["person"]],
                             sd=[p["SD.EAP"] for p in native["person"]])
        delta = lambda a,b: {k: float(np.max(abs(np.asarray(a[k])-b[k]))) for k in ("nll","eap","sd")}
        results[name] = dict(coordinates=par, config=config, rank=int(rank), design_residual=residual,
            native_iterations=native["iter"], native_iteration_cap=native["control"]["maxiter"],
            native_reached_iteration_cap=native["iter"] >= native["control"]["maxiter"],
            native_warnings=native["warnings"], native_fixed_grid=finite,
            native_replay_differences=delta(finite,native_values), continuous=adaptive,
            same_point_integration_differences=delta(finite,adaptive),
            differences_to_continuous_refit=delta(finite,ref["tight"]),
            coordinate_distance_to_continuous_refit=float(np.max(abs(par-ref["coordinates"]))))
    dump(out/"python.json", results)
    subprocess.run(["Rscript", str(ROOT/"validation/mml_pcm_quadrature_diagnosis.R"),
                    str(out/"evaluation_input.json"), str(out/"python.json"), str(out)], check=True)
    independent = {}
    for name,result in results.items():
        r = json.loads((out/f"{name}_r.json").read_text())["continuous"][1]
        independent[name] = {k: float(np.max(abs(np.asarray(result["continuous"][k])-r[k])))
                             for k in ("nll","eap","sd")}
    passed = all(max(r["native_replay_differences"].values()) < 1e-8 for r in results.values()) and \
             all(max(d.values()) < 1e-8 for d in independent.values())
    launch = json.loads((out/"conquest_launch.json").read_text())
    assert hashes == {p: sha(ROOT/p) for p in hashes}
    dump(out/"summary.json", dict(classification=spec["classification"], scientific_inference_ready=False,
        qualification_eligible=False, same_point_arithmetic_checks_pass=passed, tam=results,
        independent_r_differences=independent, conquest=launch, source_sha256=hashes,
        limitations=["One observed stress fixture; every TAM fit reached the 2000-iteration cap",
            "External reported SEs were retained, but not validated or used for inference",
            "No ConQuest model was executed successfully; no current ConQuest agreement claim",
            "Range and density changes are not independent except the two grids with spacing 0.1"],
        artifact_sha256={str(p.relative_to(out)):sha(p) for p in sorted(out.rglob('*')) if p.is_file()}))
    print("Same-point arithmetic checks:",passed,"TAM iteration caps:",sum(r["native_reached_iteration_cap"] for r in results.values()))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
