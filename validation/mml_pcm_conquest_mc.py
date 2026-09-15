#!/usr/bin/env python3
"""Prepare/run fixed-calibration ConQuest posterior Monte Carlo checks."""
import argparse
import csv
import json
from pathlib import Path
import shutil

from run_pcm_conquest_check import ROOT, launch, sha, write

SOURCE = ROOT / "validation/generated/mml_pcm_conquest_terminal_20260913_233939"
SOURCE_HASH = "a1e317deffe4b245e18421b7bbcea9a71aaa7999e224a9c72c18e2907dc7cca5"
CALIBRATION = "cq_b10_q201"
COUNTS = (2000, 20000, 200000)
SEEDS = tuple(range(2, 12))


def rows(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def command(n, seed):
    original = (SOURCE / CALIBRATION / "model.cqc").read_text().splitlines()
    start = original[:original.index("model criterion + rater + criterion*step;")+1]
    start[0] = "title Fixed PCM calibration posterior MC check;"
    start[2] = (f"set lconstraints=cases, sconstraint=none, nodefilter=0, p_nodes={n}, "
                f"seed={seed}, exit_on_error=yes, progress=no;")
    return "\n".join(start + [
        "import anchor_parameters << ../anchor_parameters.txt;",
        "import anchor_reg_coefficients << ../anchor_reg_coefficients.txt;",
        "import anchor_covariance << ../anchor_covariance.txt;",
        # Keep the previous integration method; all calibration parameters are
        # anchored. Verify exported values after every call, not just the syntax.
        "estimate ! method=quadrature, nodes=201, minnode=-10, maxnode=10, distribution=normal, "
        "fit=no, stderr=quick, abilities=eap, matrixout=check, iterations=10;",
        "export parameters ! filetype=csv >> parameters.csv;",
        "export amatrix ! filetype=csv >> amatrix.csv;",
        "export reg_coefficients ! filetype=csv >> reg_coefficients.csv;",
        "export covariance ! filetype=csv >> covariance.csv;",
        "show cases ! estimates=eap, filetype=csv >> cases.csv;",
        "write check_history ! filetype=csv >> history.csv;",
        "show parameters ! tables=1:2:3, estimates=eap >> review.txt;",
        "quit;", ""])


def prepare(out):
    assert sha(SOURCE / "review_summary.json") == SOURCE_HASH
    source = json.loads((SOURCE / "review_summary.json").read_text())
    assert not source["scientific_inference_ready"] and not source["qualification_eligible"]
    for name,h in source["artifact_sha256"].items():
        assert sha(SOURCE/name) == h,name
    calibration = rows(SOURCE/CALIBRATION/"parameters.csv")
    variance = rows(SOURCE/CALIBRATION/"covariance.csv")[0]["Covariance"]
    assert [int(r["P"]) for r in calibration] == list(range(1,24))
    assert float(rows(SOURCE/CALIBRATION/"reg_coefficients.csv")[0]["Estimate"]) == 0
    out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(SOURCE/"wide.csv",out/"wide.csv")
    (out/"anchor_parameters.txt").write_text("".join(f"{r['P']} {r['Estimate']}\n" for r in calibration))
    # This installed 5.47.5 rejects intercept index 0 despite the manual's text.
    # Use its exported index 1 and verify the anchored mean after execution.
    (out/"anchor_reg_coefficients.txt").write_text("1 1 0\n")
    (out/"anchor_covariance.txt").write_text(f"1 1 {variance}\n")
    plan = [dict(id=f"n{n}_seed{s}",p_nodes=n,seed=s,role="seed_count_comparison")
            for n in COUNTS for s in SEEDS]
    plan.append(dict(id="n2000_seed2_repeat",p_nodes=2000,seed=2,role="exact_repeat_control"))
    for config in plan:
        directory=out/config["id"]
        directory.mkdir()
        (directory/"model.cqc").write_text(command(config["p_nodes"],config["seed"]))
    write(out/"input.json",dict(classification="OBSERVED_DEVELOPMENT_ONLY",
        scientific_inference_ready=False,qualification_eligible=False,plan=plan,
        source_summary_sha256=SOURCE_HASH,calibration=source["fits"][CALIBRATION]["coordinates"],
        expected_parameter_tokens=[r["Estimate"] for r in calibration],expected_variance_token=variance,
        expected_mean=0,script_sha256=sha(Path(__file__).resolve()),
        launcher_sha256=sha(ROOT/"validation/run_pcm_conquest_check.py"),
        binary_sha256=sha(Path("/Applications/ConQuest/ConQuest")),
        metrics=["Per-person EAP/posterior-SD signed error against adaptive integration",
            "Per-seed RMSE and maximum absolute error; spread across 10 seeds",
            "Aggregate RMS error scaling with Monte Carlo count; not a per-run monotonicity gate"],
        controls=["All 23 coordinates, population mean and variance anchored in every process",
            "New process and explicit seed per condition; one duplicate seed/count control",
            "Same rounded calibration as the previous ConQuest exports; calibration uncertainty excluded"],
        input_sha256={str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file()}))


def validate_calibration(directory, spec):
    parameters = rows(directory/"parameters.csv")
    variance = rows(directory/"covariance.csv")
    beta = rows(directory/"reg_coefficients.csv")
    assert [int(r["P"]) for r in parameters] == list(range(1,24))
    assert [float(r["Estimate"]) for r in parameters] == [float(v) for v in spec["expected_parameter_tokens"]]
    assert len(variance)==len(beta)==1 and float(variance[0]["Covariance"])==float(spec["expected_variance_token"])
    assert float(beta[0]["Estimate"])==spec["expected_mean"]
    assert (directory/"amatrix.csv").read_bytes()==(SOURCE/CALIBRATION/"amatrix.csv").read_bytes()
    assert [r["PID"] for r in rows(directory/"cases.csv")]==[f"P{i:02d}" for i in range(1,25)]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory",type=Path,required=True)
    parser.add_argument("--prepare-only",action="store_true")
    args=parser.parse_args()
    out=args.directory.resolve()
    if args.prepare_only:
        prepare(out)
        print("Prepared",out)
        return 0
    spec=json.loads((out/"input.json").read_text())
    assert not spec["scientific_inference_ready"] and not spec["qualification_eligible"]
    assert sha(Path(__file__).resolve())==spec["script_sha256"]
    assert sha(ROOT/"validation/run_pcm_conquest_check.py")==spec["launcher_sha256"]
    assert sha(Path("/Applications/ConQuest/ConQuest"))==spec["binary_sha256"]
    for name,h in spec["input_sha256"].items():
        assert sha(out/name)==h,name
    sentinel=out/"sentinel"
    sentinel.mkdir(exist_ok=False)
    result=launch(sentinel,"quit;\n",60)
    if result["exit_code"]!=0 or not result["end_of_program"]:
        print("Startup failed; see",sentinel/"console.log")
        return 1
    results={}
    for config in spec["plan"]:
        name=config["id"]
        directory=out/name
        print("Scoring",name,flush=True)
        result=launch(directory,(directory/"model.cqc").read_text(),600)
        result["calibration_unchanged"]=False
        if result["exit_code"]==0 and result["end_of_program"]:
            try:
                validate_calibration(directory,spec)
                result["calibration_unchanged"]=True
            except (OSError,KeyError,ValueError,AssertionError) as error:
                result["validation_error"]=repr(error)
        results[name]=result
        if not result["calibration_unchanged"]:
            print("Stopped: calibration or output check failed for",name,flush=True)
            break
    write(out/"execution.json",dict(scientific_inference_ready=False,qualification_eligible=False,runs=results,
        remaining_unattempted=[c["id"] for c in spec["plan"] if c["id"] not in results],
        artifact_sha256={str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file()}))
    complete=len(results)==len(spec["plan"]) and all(r["calibration_unchanged"] for r in results.values())
    print("Saved:",out,"All fixed-calibration checks:",complete,flush=True)
    return 0 if complete else 1


if __name__=="__main__":
    raise SystemExit(main())
