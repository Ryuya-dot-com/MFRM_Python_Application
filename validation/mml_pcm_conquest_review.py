#!/usr/bin/env python3
"""Review actual ConQuest exports at the same PCM coordinates and scoring target."""
import argparse
import csv
import json
from pathlib import Path
import re
import subprocess

import numpy as np

from mml_pcm_continuous_refit import PCMIntegral, PREVIOUS, ROOT, dump, sha
from mml_pcm_external_review import fixed_grid


def read_csv(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def distance(a,b):
    return float(np.max(abs(np.asarray(a)-b)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory",type=Path,required=True)
    out = parser.parse_args().directory.resolve()
    execution = json.loads((out/"execution.json").read_text())
    assert not execution["scientific_inference_ready"] and not execution["qualification_eligible"]
    for name,h in execution["artifact_sha256"].items():
        assert sha(out/name) == h,name
    spec = json.loads((out/"input.json").read_text())
    assert not execution["remaining_unattempted"] and len(execution["runs"]) == 7
    assert all(r["execution_complete"] for r in execution["runs"].values())
    original = json.loads((PREVIOUS/"input.json").read_text())
    refdir = ROOT/"validation/generated/mml_pcm_continuous_refit_20260913"
    assert sha(refdir/"summary.json") == "d92f7359e826772f7a21b498fa1f44277c612165924c3766cff7c2dabda8630f"
    ref = json.loads((refdir/"from_q181.json").read_text())
    tam = json.loads((ROOT/"validation/generated/mml_pcm_tam_grid_diagnosis_20260913/python.json").read_text())
    problem = PCMIntegral(original["data"])
    wide = read_csv(out/"wide.csv")
    assert [r["Person"] for r in wide] == [f"P{i:02d}" for i in range(1,25)]
    assert np.array_equal([[int(r[f"Y{i}"]) for i in range(1,21)] for r in wide],problem.y)
    hashes = {str(p.relative_to(ROOT)):sha(p) for p in (Path(__file__).resolve(),
        ROOT/"validation/mml_pcm_continuous_refit.py",ROOT/"validation/mml_pcm_external_review.py",
        ROOT/"validation/mml_pcm_quadrature_diagnosis.R")}
    evaluation = {k:original[k] for k in ("classification","scientific_inference_ready",
        "qualification_eligible","data","orders","integration")}
    evaluation.update(source_sha256=hashes,scope="Observed, rounded exports; no scientific equivalence margins",
        independent_same_point_limit=1e-8,execution_sha256=sha(out/"execution.json"))
    dump(out/"evaluation_input.json",evaluation)
    fits = {}
    for config in spec["plan"]:
        name = config["id"]
        directory = out/name
        native = read_csv(directory/"parameters.csv")
        estimates = np.array([float(r["Estimate"]) for r in native])
        arows = read_csv(directory/"amatrix.csv")
        assert len(native) == 23 and len(arows) == 100
        labels = [r["Label"] for r in native]
        assert list(arows[0])[2:] == labels
        assert [(int(r["GIN"]),int(r["Category"])) for r in arows] == [(i,k) for i in range(1,21) for k in range(1,6)]
        A = np.array([[float(r[l]) for l in labels] for r in arows])
        assert np.array_equal(A,problem.design.reshape(100,23)[:,[3,4,5,6,7,0,1,2,*range(8,23)]])
        variance = read_csv(directory/"covariance.csv")
        beta = read_csv(directory/"reg_coefficients.csv")
        assert len(variance) == len(beta) == 1 and float(beta[0]["Estimate"]) == 0
        var = float(variance[0]["Covariance"])
        par = np.r_[estimates[5:8],estimates[:5],estimates[8:],.5*np.log(var)]
        review = (directory/"review.txt").read_text()
        assert "Location constraint was: CASES" in review and "Slopes are fixed" in review
        assert re.search(r"Cases in MML/MCMC estimation:\s*24\b",review)
        iteration = int(re.search(r"The number of iterations:\s*(\d+)",review)[1])
        termination = re.search(r"Iterations terminated[^\n]+",review)[0]
        history = read_csv(directory/"history.csv")
        selected = next(r for r in history if int(r["Iteration"]) == iteration)
        assert np.array_equal([float(selected[f"xsi {i}"]) for i in range(1,24)],estimates)
        assert float(selected["wvar 1 1"]) == var
        # Despite its column heading, this field contains deviance in this export.
        native_nll = float(selected["LogLikelihood"])/2
        assert abs(2*native_nll-float(re.search(r"Final Deviance:\s*([\d.]+)",review)[1])) <= 5.1e-6
        cases = read_csv(directory/"cases.csv")
        assert [r["PID"] for r in cases] == [f"P{i:02d}" for i in range(1,25)]
        eap = np.array([float(r["EAP_1"]) for r in cases])
        sd = np.array([float(r["PosteriorSD_1"]) for r in cases])
        continuous = problem.evaluate(par,**original["integration"][1])
        fits[name] = dict(config=config,coordinates=par,sigma=float(np.sqrt(var)),
            selected_iteration=iteration,executed_history_rows=len(history),termination=termination,
            native_nll=native_nll,native_eap=eap,native_posterior_sd=sd,
            native_estimate_tokens=[r["Estimate"] for r in native],native_variance_token=variance[0]["Covariance"],
            continuous=continuous,coordinate_distance_to_reference=distance(par,ref["coordinates"]),
            continuous_nll_excess=continuous["nll"]-ref["tight"]["nll"],
            conditional_score_difference_to_reference={k:distance(continuous[k],ref["tight"][k]) for k in ("eap","sd")},
            native_mc_scoring_difference=dict(eap=distance(eap,continuous["eap"]),sd=distance(sd,continuous["sd"])))
        if config["method"] == "quadrature":
            finite = fixed_grid(problem,par,np.linspace(-config["bound"],config["bound"],config["q"]))
            target = tam[f"aligned_b{config['bound']}_q{config['q']}"]
            fits[name].update(fixed_grid=finite,
                normalized_grid_nll=finite["nll"]+24*np.log(finite["raw_prior_mass"]),
                tam_coordinate_difference=distance(par,target["coordinates"]),
                tam_sigma_difference=abs(np.sqrt(var)-target["sigma"]),
                tam_continuous_scoring_difference={k:distance(continuous[k],target["continuous"][k]) for k in ("eap","sd")})
    dump(out/"python.json",fits)
    subprocess.run(["Rscript",str(ROOT/"validation/mml_pcm_quadrature_diagnosis.R"),
        str(out/"evaluation_input.json"),str(out/"python.json"),str(out)],check=True)
    independent = {}
    for name,fit in fits.items():
        r = json.loads((out/f"{name}_r.json").read_text())
        independent[name] = {k:distance(fit["continuous"][k],r["continuous"][1][k]) for k in ("nll","eap","sd")}
        if fit["config"]["method"] == "gauss":
            finite = r["finite"][str(fit["config"]["q"])]
            fit["independent_finite_gh"] = finite
            fit["finite_gh_native_nll_difference"] = finite["nll"]-fit["native_nll"]
            fit["finite_gh_continuous_difference"] = {k:distance(finite[k],fit["continuous"][k]) for k in ("nll","eap","sd")}
    passed = all(max(d.values())<1e-8 for d in independent.values())
    assert hashes == {p:sha(ROOT/p) for p in hashes}
    dump(out/"review_summary.json",dict(classification="OBSERVED_DEVELOPMENT_ONLY",scientific_inference_ready=False,
        qualification_eligible=False,independent_arithmetic_checks_pass=passed,fits=fits,
        independent_r_differences=independent,source_sha256=hashes,
        limitations=["One observed fixture and one native start per ConQuest condition",
            "Native parameter exports rounded to six decimal places; no bitwise equivalence claim",
            "ConQuest EAP and posterior SD use separate Monte Carlo scoring with p_nodes=2000",
            "GH Q61 and Q121 stopped without further deviance improvement; selected and executed iterations differ",
            "No ConQuest SE/CI/coverage qualification or user-facing estimator change"],
        artifact_sha256={str(p.relative_to(out)):sha(p) for p in sorted(out.rglob('*')) if p.is_file()}))
    print("Independent same-point arithmetic checks:",passed,"Native numerical agreement remains condition-specific.")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
