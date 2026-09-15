#!/usr/bin/env python3
"""Summarize ConQuest scoring variation with calibration held fixed."""
import argparse
import json
from pathlib import Path
import re

import numpy as np

from mml_pcm_conquest_mc import COUNTS, SEEDS, SOURCE, SOURCE_HASH, rows, validate_calibration
from run_pcm_conquest_check import ROOT, sha, write


def error_summary(errors):
    errors=np.asarray(errors,dtype=float)
    assert errors.shape==(10,24) and np.all(np.isfinite(errors))
    rmse=np.sqrt(np.mean(errors**2,axis=1))
    maxima=np.max(abs(errors),axis=1)
    return dict(aggregate_rmse=float(np.sqrt(np.mean(errors**2))),
        seed_rmse=rmse.tolist(),seed_max_absolute_error=maxima.tolist(),
        median_seed_max_absolute_error=float(np.median(maxima)),
        worst_seed_max_absolute_error=float(np.max(maxima)),
        per_person_mean_signed_error=errors.mean(axis=0).tolist(),
        per_person_across_seed_sd=errors.std(axis=0,ddof=1).tolist())


def main():
    # Small analytic controls remain runnable with every review.
    assert error_summary(np.zeros((10,24)))["aggregate_rmse"]==0
    control=error_summary(np.tile(np.arange(10.)[:,None],(1,24)))
    assert abs(control["aggregate_rmse"]-np.sqrt(28.5))<1e-12
    assert control["seed_max_absolute_error"]==list(np.arange(10.))
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory",type=Path,required=True)
    args=parser.parse_args()
    out=args.directory.resolve()
    spec=json.loads((out/"input.json").read_text())
    execution=json.loads((out/"execution.json").read_text())
    ref=json.loads((out/"reference.json").read_text())
    assert not spec["scientific_inference_ready"] and not spec["qualification_eligible"]
    assert not execution["remaining_unattempted"] and len(execution["runs"])==31
    assert all(r["calibration_unchanged"] for r in execution["runs"].values())
    assert sha(SOURCE/"review_summary.json")==ref["source_summary_sha256"]==SOURCE_HASH
    assert sha(SOURCE/"cq_b10_q201_r.json")==ref["source_R_sha256"]
    assert ref["coordinates"]==spec["calibration"]
    for name,h in execution["artifact_sha256"].items():
        assert sha(out/name)==h,name
    for name,h in spec["input_sha256"].items():
        assert sha(out/name)==h,name
    for key in ("nll","eap","sd"):
        assert np.max(abs(np.asarray(ref["python"][key])-ref["independent_R"][key])) < 1e-8
    exact={"eap":np.asarray(ref["independent_R"]["eap"]),"sd":np.asarray(ref["independent_R"]["sd"])}
    all_scores={}
    native_parameter_counts={}
    for config in spec["plan"]:
        name=config["id"]
        directory=out/name
        validate_calibration(directory,spec)
        review=(directory/"review.txt").read_text()
        assert int(re.search(r"Number of nodes used when drawing PVs:\s*(\d+)",review)[1])==config["p_nodes"]
        assert float(re.search(r"Random number generation seed:\s*([\d.]+)",review)[1])==config["seed"]
        native_parameter_counts[name]=int(re.search(r"Total number of estimated parameters:\s*(-?\d+)",review)[1])
        data=rows(directory/"cases.csv")
        values={"eap":np.array([float(r["EAP_1"]) for r in data]),
                "sd":np.array([float(r["PosteriorSD_1"]) for r in data])}
        assert all(np.all(np.isfinite(a)) for a in values.values()) and np.all(values["sd"]>0)
        all_scores[name]={k:v.tolist() for k,v in values.items()}
    repeated={k:bool(np.array_equal(all_scores['n2000_seed2'][k],all_scores['n2000_seed2_repeat'][k]))
              for k in exact}
    groups={}
    for n in COUNTS:
        groups[str(n)]={k:error_summary(np.array([all_scores[f"n{n}_seed{s}"][k] for s in SEEDS])-target)
                        for k,target in exact.items()}
    changes={k:{"aggregate_rmse_ratio_20000_to_2000":groups['20000'][k]['aggregate_rmse']/groups['2000'][k]['aggregate_rmse'],
        "aggregate_rmse_ratio_200000_to_2000":groups['200000'][k]['aggregate_rmse']/groups['2000'][k]['aggregate_rmse'],
        "log_log_slope":float(np.polyfit(np.log(COUNTS),
            np.log([groups[str(n)][k]['aggregate_rmse'] for n in COUNTS]),1)[0])} for k in exact}
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(9.2,4.5),layout="constrained")
    for ax,key,label in zip(axes,("eap","sd"),("EAP","Posterior SD")):
        aggregate=np.array([groups[str(n)][key]['aggregate_rmse'] for n in COUNTS])
        seed_rmse=np.array([groups[str(n)][key]['seed_rmse'] for n in COUNTS])
        ax.fill_between(COUNTS,seed_rmse.min(axis=1),seed_rmse.max(axis=1),color="#dce9ef",label="Range across 10 seeds")
        ax.loglog(COUNTS,aggregate,"o-",color="#17617d",linewidth=2,label="Aggregate RMSE")
        ax.loglog(COUNTS,aggregate[0]*np.sqrt(COUNTS[0]/np.array(COUNTS)),"--",color="#797979",label="1 / sqrt(N), anchored at 2,000")
        ax.set_title(label)
        ax.set_xticks(COUNTS,["2,000","20,000","200,000"])
        ax.set_xlabel("ConQuest scoring draws (p_nodes)")
        ax.set_ylabel("Error against adaptive integration (logits)")
        ax.grid(True,which="major",alpha=.2)
    axes[0].legend(fontsize=7,loc="lower left")
    fig.suptitle("Monte Carlo scoring error at one fixed PCM calibration",fontsize=13)
    for suffix in ("png","svg"):
        fig.savefig(out/f"mc_scoring_error.{suffix}",dpi=180)
    plt.close(fig)
    write(out/"review_summary.json",dict(classification="OBSERVED_DEVELOPMENT_ONLY",scientific_inference_ready=False,
        qualification_eligible=False,all_calibrations_unchanged=True,repeat_exact=repeated,
        native_reported_parameter_counts=native_parameter_counts,
        native_parameter_count_expected=0,native_parameter_count_checks_pass=all(n==0 for n in native_parameter_counts.values()),
        native_information_criteria_used=False,
        counts=COUNTS,seeds=SEEDS,groups=groups,scaling=changes,scores=all_scores,
        reference_sha256=sha(out/"reference.json"),execution_sha256=sha(out/"execution.json"),
        source_sha256={str(p.relative_to(ROOT)):sha(p) for p in
            (Path(__file__).resolve(),ROOT/"validation/mml_pcm_conquest_mc.py",ROOT/"validation/run_pcm_conquest_check.py")},
        limitations=["One observed fixture, one fixed rounded calibration and ten seeds per count",
            "Conditional scoring only: calibration uncertainty and recovery are not evaluated",
            "Native parameter-count bookkeeping is retained separately from verified fixed calibration; no native IC or SE validation",
            "Error need not decrease for every seed or person; no retrospective error acceptance threshold",
            "The log-log slope is descriptive, not proof of a universal Monte Carlo rate"],
        artifact_sha256={str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file()}))
    print("All 31 calibrations unchanged; repeat exact:",repeated)
    for n,g in groups.items():
        print(n,{k:dict(rmse=v['aggregate_rmse'],worst=v['worst_seed_max_absolute_error']) for k,v in g.items()})
    return 0 if all(repeated.values()) else 1


if __name__=="__main__":
    raise SystemExit(main())
