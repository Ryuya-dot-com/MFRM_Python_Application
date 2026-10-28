"""Review all retained quadrature checks, including expanded-parameter SEs.

Post-result development audit; never a prospective qualification decision.
Reads frozen artifacts only and does not refit or select a preferred start.
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def read(path):
    return json.loads(path.read_text())


def verify(directory):
    manifest = read(directory / "manifest.json")
    for name, expected in manifest.items():
        if hashlib.sha256((directory / name).read_bytes()).hexdigest() != expected:
            raise ValueError(f"Changed artifact: {directory / name}")


def expansion(config, point):
    """Explicit sum-zero expansion; last coordinate is population SD."""
    labels, blocks = [], []
    for facet in config["facet_names"]:
        levels = config["facet_levels"][facet]
        centered = facet != config["noncenter_facet"]
        matrix = np.eye(len(levels), len(levels)-int(centered))
        if centered:
            matrix[-1] = -1.
        blocks.append(matrix)
        labels.extend(f"{facet}:{level}" for level in levels)
    steps = config["n_cat"]-1
    step_levels = ["shared"] if config["model"] == "RSM" else config["facet_levels"][config["step_facet"]]
    for level in step_levels:
        matrix = np.eye(steps, steps-1); matrix[-1] = -1.
        blocks.append(matrix)
        labels.extend(f"step:{level}:{k+1}" for k in range(steps))
    jacobian = np.zeros((len(labels)+1, len(point)))
    row, column = 0, 0
    for matrix in blocks:
        nr, nc = matrix.shape
        jacobian[row:row+nr, column:column+nc] = matrix
        row += nr; column += nc
    if column != len(point)-1:
        raise ValueError("Expanded coordinates do not match the retained point")
    jacobian[-1, -1] = np.exp(point[-1])
    expanded = jacobian @ point
    expanded[-1] = np.exp(point[-1])
    return labels + ["population_sd"], expanded, jacobian


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("parent", type=Path)
    parser.add_argument("study", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    verify(args.parent); verify(args.study)
    protocol, summary = read(args.study/"protocol.json"), read(args.study/"summary.json")
    if protocol["parent_manifest_sha256"] != hashlib.sha256((args.parent/"manifest.json").read_bytes()).hexdigest():
        raise ValueError("Study belongs to a different parent")
    configs = {case["id"].rsplit("_Q", 1)[0]: case["config"]
               for case in read(args.parent/"reference_input.json")["cases"]}
    r_checks = {case["id"]: case["passed"] for case in summary["R_comparisons"]}
    records = []
    for key, brief in summary["cases"].items():
        case_id, q_label = key.rsplit("_Q", 1)
        q = int(q_label.split("_")[0]); result = read(args.study/(key+".json"))
        point = np.array(result["point"])
        labels, expanded, jacobian = expansion(configs[case_id], point)
        if np.max(abs(expanded-np.array(result["expanded"]))) > 1e-10:
            raise ValueError(f"Coordinate expansion disagrees: {key}")
        ses = {}
        for kind, cov in result["covariance"].items():
            variance = np.einsum("ij,jk,ik->i", jacobian, np.array(cov), jacobian)
            if np.any(variance < 0):
                raise ValueError(f"Negative expanded variance: {key}")
            ses[kind] = np.sqrt(variance)
        relative, worst = None, None
        if set(ses) == {"finite", "continuous"}:
            valid = ses["continuous"] > 0
            if np.any((~valid) & (ses["finite"] != 0)):
                raise ValueError(f"Fixed-coordinate variance changed: {key}")
            errors = abs(ses["finite"][valid]/ses["continuous"][valid]-1.)
            relative = float(errors.max())
            worst = np.array(labels)[valid][errors.argmax()]
        expanded_pass = relative is not None and relative < protocol["limits"]["gh_se_relative"]
        checks = dict(result["checks"], expanded_se=expanded_pass)
        if key in r_checks:
            checks["independent_R"] = r_checks[key]
        if brief["solver_terminated"] is not None:
            checks["solver_terminated"] = brief["solver_terminated"]
        if q == 241:
            checks["cross_start"] = summary["Q241_cross_start_coordinates"].get(case_id, float("inf")) < protocol["limits"]["same_q_cross_start_coordinates"]
        records.append(dict(id=key, q=q, checks=checks, all_development_checks_pass=all(checks.values()),
            expanded_se_relative_difference=relative, worst_expanded_se_coordinate=worst,
            expanded_labels=labels, expanded_se={k:v.tolist() for k,v in ses.items()},
            continuous_score_supnorm=float(np.max(abs(np.array(result["continuous"][-1]["gradient"])))),
            finite_score_supnorm=float(np.max(abs(np.array(result["app_gradient"])))),
            integration=brief["integration"]))
    complete = (summary["completed_points"] == 36 and not summary["failures"]
                and len(r_checks) == 24 and summary["R_returncode"] == 0)
    highest = [r for r in records if r["q"] == 241]
    report = dict(classification="POST_RESULT_DEVELOPMENT_REVIEW", scientific_inference_ready=False,
        scope="Expanded SE agreement at retained points plus all original gates; no recovery or coverage evidence",
        complete=complete, highest_order_all_checks_pass=complete and len(highest) == 12 and all(r["all_development_checks_pass"] for r in highest),
        parent_manifest_sha256=protocol["parent_manifest_sha256"],
        study_manifest_sha256=hashlib.sha256((args.study/"manifest.json").read_bytes()).hexdigest(),
        cases=records)
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output/"review.py").write_bytes(Path(__file__).read_bytes())
    (args.output/"review.json").write_text(json.dumps(report, indent=2)+"\n")
    (args.output/"manifest.json").write_text(json.dumps({p.name:hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(args.output.iterdir()) if p.is_file()}, indent=2)+"\n")
    print(json.dumps({k:report[k] for k in ("complete", "highest_order_all_checks_pass", "scientific_inference_ready")}))
    return 0 if report["highest_order_all_checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
