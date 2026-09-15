#!/usr/bin/env python3
"""Isolate the explicit mean anchor with two otherwise identical ConQuest runs."""
import argparse
import json
from pathlib import Path
import re
import shutil

from mml_pcm_conquest_mc import rows, validate_calibration
from run_pcm_conquest_check import ROOT, launch, sha, write

SOURCE = ROOT / "validation/generated/mml_pcm_conquest_mc_index1_20260913"
SOURCE_HASH = "279ec456a8e37332d6b8e56b89795b82ac23c1c5de8cc01a439c065bc7e0c0c4"
ANCHOR = "import anchor_reg_coefficients << ../anchor_reg_coefficients.txt;\n"
CONDITIONS = ("explicit_mean_anchor", "cases_mean_only")


def report(text):
    return {key: float(re.search(pattern, text)[1]) for key, pattern in {
        "parameter_count": r"Total number of estimated parameters:\s*(-?\d+)",
        "deviance": r"Final Deviance:\s*([\d.]+)",
        "aic": r"Akaike Information Criterion \(AIC\):\s*([\d.]+)",
    }.items()}


def main():
    assert report("Total number of estimated parameters: -1\nFinal Deviance: 10\n"
                  "Akaike Information Criterion (AIC): 8")==dict(parameter_count=-1., deviance=10., aic=8.)
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--prepare-only", action="store_true")
    args=parser.parse_args()
    out=args.directory.resolve()
    assert sha(SOURCE/"review_summary.json")==SOURCE_HASH
    previous=json.loads((SOURCE/"review_summary.json").read_text())
    for name,h in previous["artifact_sha256"].items():
        assert sha(SOURCE/name)==h, name
    if args.prepare_only:
        spec=json.loads((SOURCE/"input.json").read_text())
        original=(SOURCE/"n2000_seed2/model.cqc").read_text()
        assert original.count(ANCHOR)==1 and "lconstraints=cases" in original
        out.mkdir(parents=True, exist_ok=False)
        for name in ("wide.csv", "anchor_parameters.txt", "anchor_covariance.txt", "anchor_reg_coefficients.txt"):
            shutil.copyfile(SOURCE/name,out/name)
        for name in CONDITIONS:
            (out/name).mkdir()
            command=original if name==CONDITIONS[0] else original.replace(ANCHOR, "")
            (out/name/"model.cqc").write_text(command)
        spec={k:spec[k] for k in ("expected_parameter_tokens", "expected_variance_token", "expected_mean", "binary_sha256")}
        spec.update(classification="OBSERVED_DEVELOPMENT_ONLY", scientific_inference_ready=False,
            qualification_eligible=False, source_summary_sha256=SOURCE_HASH,
            expected_counts=dict(explicit_mean_anchor=-1, cases_mean_only=0),
            contrast="Remove only the explicit mean anchor; CASES, all other anchors, seed and p_nodes unchanged",
            source_sha256={str(p.relative_to(ROOT)):sha(p) for p in
                (Path(__file__).resolve(), ROOT/"validation/mml_pcm_conquest_mc.py", ROOT/"validation/run_pcm_conquest_check.py")},
            input_sha256={str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file()})
        write(out/"input.json",spec)
        print("Prepared",out)
        return 0
    spec=json.loads((out/"input.json").read_text())
    for name,h in spec["source_sha256"].items():
        assert sha(ROOT/name)==h, name
    for name,h in spec["input_sha256"].items():
        assert sha(out/name)==h, name
    assert sha(Path("/Applications/ConQuest/ConQuest"))==spec["binary_sha256"]
    (out/"sentinel").mkdir(exist_ok=False)
    sentinel=launch(out/"sentinel", "quit;\n", 60)
    if sentinel["exit_code"]!=0 or not sentinel["end_of_program"]:
        print("Startup failed; see",out/"sentinel/console.log")
        return 1
    results={}
    prior_scores=rows(SOURCE/"n2000_seed2/cases.csv")
    for name in CONDITIONS:
        directory=out/name
        result=launch(directory,(directory/"model.cqc").read_text(),600)
        result["calibration_unchanged"]=False
        if result["exit_code"]==0 and result["end_of_program"]:
            try:
                validate_calibration(directory,spec)
                result["calibration_unchanged"]=True
                result["native_report"]=report((directory/"review.txt").read_text())
                scores=rows(directory/"cases.csv")
                result["max_score_difference_to_previous"]={k:max(abs(float(a[k])-float(b[k]))
                    for a,b in zip(scores,prior_scores)) for k in ("EAP_1","PosteriorSD_1")}
                result["expected_parameter_count_matches"]=(result["native_report"]["parameter_count"]==spec["expected_counts"][name])
            except (OSError,KeyError,ValueError,TypeError,AssertionError) as error:
                result["review_error"]=repr(error)
        results[name]=result
        print(name,result,flush=True)
        if not result["calibration_unchanged"] or "review_error" in result:
            break
    write(out/"summary.json",dict(classification="OBSERVED_DEVELOPMENT_ONLY", scientific_inference_ready=False,
        qualification_eligible=False, native_information_criteria_used_for_model_selection=False,
        runs=results, remaining_unattempted=[n for n in CONDITIONS if n not in results],
        artifact_sha256={str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file()}))
    complete=len(results)==2 and all(r["calibration_unchanged"] and "review_error" not in r for r in results.values())
    print("Execution complete:",complete,"Saved:",out,"Numerical interpretation pending.")
    return 0 if complete else 1


if __name__=="__main__":
    raise SystemExit(main())
