#!/usr/bin/env python3
"""Separate TAM stopping, grid range and density on the retained PCM fixture."""
import argparse
import json
from pathlib import Path
import subprocess

import numpy as np

from mml_pcm_continuous_refit import PCMIntegral, PREVIOUS, ROOT, dump, sha, fd_check
from mml_pcm_external_review import fixed_grid


def distance(a, b):
    return float(np.max(abs(np.asarray(a)-b)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    out = parser.parse_args().directory.resolve()
    spec = json.loads((out/"input.json").read_text())
    assert not spec["scientific_inference_ready"] and not spec["qualification_eligible"]
    for p, h in spec["source_sha256"].items():
        assert sha(ROOT/p) == h, p
    old = ROOT/"validation/generated/mml_pcm_external_refit_20260913"
    assert sha(old/"summary.json") == spec["previous_summary_sha256"]
    reference = ROOT/"validation/generated/mml_pcm_continuous_refit_20260913"
    assert sha(reference/"summary.json") == "d92f7359e826772f7a21b498fa1f44277c612165924c3766cff7c2dabda8630f"
    ref = json.loads((reference/"from_q181.json").read_text())
    original = json.loads((PREVIOUS/"input.json").read_text())
    assert spec["data"] == original["data"]
    problem = PCMIntegral(spec["data"])
    evaluation = {k: original[k] for k in ("classification", "scientific_inference_ready",
        "qualification_eligible", "data", "orders", "integration")}
    hashes = dict(spec["source_sha256"])
    hashes.update({str(p.relative_to(ROOT)):sha(p) for p in
        (Path(__file__).resolve(), ROOT/"validation/mml_pcm_quadrature_diagnosis.R")})
    evaluation.update(source_sha256=hashes, same_point_limit=spec["limits"]["same_point"])
    dump(out/"evaluation_input.json", evaluation)
    fits, replay = {}, {}
    for config in spec["tam"]:
        name = config["id"]
        native = json.loads((out/f"{name}.json").read_text())
        assert config == native["config"] and native["beta"] == [[0]]
        assert len(native["item"]) == 20 and all(i["N"] == 24 for i in native["item"])
        assert [p["pid"] for p in native["person"]] == [f"P{i:02d}" for i in range(1,25)]
        assert [p["score"] for p in native["person"]] == problem.total.tolist()
        byname = {r["item"]:i for i,r in enumerate(native["item"])}
        order = [byname[f"C{c}-raterR{r}"] for r in range(1,5) for c in range(1,6)]
        offset = np.array(native["AXsi"])[order]
        assert np.array_equal(np.array(native["B"])[order,:,0], np.tile(np.arange(5),(20,1)))
        coord, _, rank, _ = np.linalg.lstsq(problem.design.reshape(100,23),offset.ravel(),rcond=None)
        residual = distance(problem.design@coord, offset)
        assert rank == 23 and residual < 1e-10
        par = np.r_[coord,.5*np.log(native["variance"][0][0])]
        finite = fixed_grid(problem, par, np.array(native["control"]["nodes"]))
        continuous = problem.evaluate(par, **original["integration"][1])
        fits[name] = dict(config=config, coordinates=par, design_residual=residual,
            sigma=float(np.exp(par[-1])), iter=native["iter"],
            loop_checks=native["terminal_loop_checks"], reached_cap=native["reached_iteration_cap"],
            last_updates=native["updates"][-1], finite=finite, continuous=continuous,
            native_replay=dict(nll=distance(finite["nll"],native["deviance"]/2),
                eap=distance(finite["eap"],[p["EAP"] for p in native["person"]]),
                sd=distance(finite["sd"],[p["SD.EAP"] for p in native["person"]])),
            integration_error={k:distance(finite[k],continuous[k]) for k in ("nll","eap","sd")},
            error_to_reference={k:distance(finite[k],ref["tight"][k]) for k in ("nll","eap","sd")},
            coordinate_distance_to_reference=distance(par,ref["coordinates"]))
        if config["role"] == "trace_previous_settings":
            previous = json.loads((old/f"tam_b{config['bound']}_q{config['q']}.json").read_text())
            replay[name] = {k:distance(native[k],previous[k]) for k in ("AXsi","variance","deviance")}
    # Independent check of the fixed-grid score, including the prior SD score.
    point = np.array(ref["coordinates"])
    gradient_fd = fd_check(lambda p:fixed_grid(problem,p,np.linspace(-6,6,61)),point)
    assert max(c["max_difference"] for c in gradient_fd) < spec["limits"]["gradient_fd"]
    dump(out/"python.json",fits)
    subprocess.run(["Rscript",str(ROOT/"validation/mml_pcm_quadrature_diagnosis.R"),
        str(out/"evaluation_input.json"),str(out/"python.json"),str(out)],check=True)
    independent = {}
    for name,fit in fits.items():
        r = json.loads((out/f"{name}_r.json").read_text())["continuous"][1]
        independent[name] = {k:distance(fit["continuous"][k],r[k]) for k in ("nll","eap","sd")}
    contrasts = []
    def contrast(a,b,changed):
        x,y = fits[a],fits[b]
        contrasts.append(dict(first=a,second=b,changed=changed,
            sigma_difference=abs(x["sigma"]-y["sigma"]),
            coordinate_difference=distance(x["coordinates"],y["coordinates"]),
            finite_output_difference={k:distance(x["finite"][k],y["finite"][k]) for k in ("nll","eap","sd")}))
    for b in (6,8,10,12):
        contrast(f"aligned_b{b}_q{10*b+1}",f"aligned_b{b}_q{20*b+1}","density_only")
    for step in (.2,.1):
        for b in (6,8,10):
            contrast(f"aligned_b{b}_q{round(2*b/step)+1}",
                     f"aligned_b{b+2}_q{round(2*(b+2)/step)+1}","range_only")
    contrast("aligned_b10_q201","aligned_b10_q401","density_only")
    aligned = [f for f in fits.values() if f["config"]["role"] != "trace_previous_settings"]
    passed = (all(max(d.values()) < spec["limits"]["replay"] for d in replay.values()) and
        all(max(f["native_replay"].values()) < spec["limits"]["same_point"] for f in fits.values()) and
        all(max(d.values()) < spec["limits"]["same_point"] for d in independent.values()) and
        all(not f["reached_cap"] and all(f["loop_checks"].values()) for f in aligned))
    assert hashes == {p:sha(ROOT/p) for p in hashes}
    dump(out/"summary.json",dict(classification=spec["classification"],scientific_inference_ready=False,
        qualification_eligible=False, implementation_checks_pass=passed, fits=fits,
        previous_settings_replay=replay, independent_r_differences=independent,
        fixed_grid_gradient_fd=gradient_fd, contrasts=contrasts, source_sha256=hashes,
        limitations=["Observed single PCM fixture with two forced extreme persons",
            "New stopping settings are exploratory diagnostics, not retrospective success relabeling",
            "Loop termination is not a gradient, integration-error or inference certificate"],
        artifact_sha256={p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file()}))
    print("Implementation checks:",passed,"Aligned normal terminations:",
          sum(not f["reached_cap"] and all(f["loop_checks"].values()) for f in aligned),"/",len(aligned))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
