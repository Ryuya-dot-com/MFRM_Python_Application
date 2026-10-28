"""Independent same-point and refit comparison of frozen native TAM results."""

import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import sys
import zipfile

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from validation.mml_rasch_integral import RaschIntegral
from validation.mml_rasch_quadrature_review import expansion, read, verify


def save(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False,
        default=lambda x: x.tolist() if isinstance(x, np.ndarray) else x.item())+"\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("parent", "tam", "gh", "output"):
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    for folder in (args.parent, args.tam, args.gh):
        verify(folder)
    tam_protocol, gh_protocol = read(args.tam/"protocol.json"), read(args.gh/"protocol.json")
    parent_digest = hashlib.sha256((args.parent/"manifest.json").read_bytes()).hexdigest()
    if any(p["parent_manifest_sha256"] != parent_digest for p in (tam_protocol, gh_protocol)):
        raise ValueError("Unmatched parent datasets")
    args.output.mkdir(parents=True, exist_ok=False)
    save(args.output/"protocol.json", dict(classification="POST_RESULT_DEVELOPMENT_ONLY", scientific_inference_ready=False,
        tam_limits=tam_protocol["limits"], continuous_settings=gh_protocol["continuous_settings"],
        additional_limits=dict(reference_information=1e-7, reference_gradient=1e-8,
            reference_nll=1e-8, reference_moments=1e-9, tam_returned_moments=1e-8),
        se_definition="TAM native SE retained separately; inference comparison uses the full continuous joint information",
        comparison="Both native TAM grids versus Q241_from61 for every case; Q241 cross-start audit retained separately",
        input_manifests={n:hashlib.sha256((getattr(args,n)/"manifest.json").read_bytes()).hexdigest() for n in ("parent","tam","gh")}))
    with zipfile.ZipFile(args.output/"sources.zip", "w", zipfile.ZIP_DEFLATED) as archive:
        for p in (Path(__file__), ROOT/"validation/mml_rasch_integral.py", ROOT/"validation/mml_rasch_quadrature_review.py"):
            archive.write(p, str(p.relative_to(ROOT)))
    save(args.output/"runtime.json", dict(python=sys.version,
        packages={d.metadata["Name"]:d.version for d in importlib.metadata.distributions()}))
    planned = read(args.parent/"prefit_manifest.json")
    records, failures = [], []
    for case in planned:
        case_id = case["id"]
        data = pd.read_csv(args.parent/(case_id+".csv"))
        design = read(args.tam/(case_id+"_design.json")); config = design["config"]
        reference = RaschIntegral(data, config)
        item_map = {tuple(row[f] for f in config["facet_names"]): i for i,row in enumerate(design["items"])}
        item_indices = [item_map[tuple(row[f] for f in config["facet_names"])] for row in data.to_dict("records")]
        np.testing.assert_allclose(np.array(design["A"])[item_indices], reference.design, rtol=0, atol=1e-14)
        np.testing.assert_allclose(np.array(design["B"])[:,:,0], np.broadcast_to(reference.k, (len(design["items"]),len(reference.k))), rtol=0, atol=0)
        response = np.array(design["response"], dtype=float)
        assert np.isfinite(response).sum() == len(data)
        np.testing.assert_array_equal(response[reference.person, item_indices], reference.y)
        assert design["persons"] == reference.person_labels
        app_point = read(args.gh/(case_id+"_Q241_from61.json"))
        for grid in tam_protocol["grids"]:
            key = case_id+"_"+grid["id"]
            try:
                tam = read(args.tam/(key+".json")); point = np.array(tam["point"])
                sigma = np.exp(point[-1]); theta = np.array(tam["theta"]).ravel()
                weights = np.exp(-.5*(theta/sigma)**2)/sigma/np.sqrt(2*np.pi)*grid["spacing"]
                # Values/moments only: a fixed-theta grid has different SD
                # derivatives from a moving-node GH rule.
                fixed_grid = reference.finite(point, theta/sigma, weights/weights.sum())
                continuous = [reference.continuous(point, **s, information=True) for s in gh_protocol["continuous_settings"]]
                tight = continuous[-1]; info = tight["information"]
                if np.linalg.eigvalsh(info).min() <= 0:
                    raise ValueError("Continuous joint information is not positive definite")
                covariance = np.linalg.inv(info)
                labels, expanded, jacobian = expansion(config, point)
                se = np.sqrt(np.einsum("ij,jk,ik->i", jacobian, covariance, jacobian))
                _, _, app_jacobian = expansion(config, np.array(app_point["point"]))
                app_se = np.sqrt(np.einsum("ij,jk,ik->i", app_jacobian, np.array(app_point["covariance"]["continuous"]), app_jacobian))
                errors = dict(grid_nll_per_person=abs(fixed_grid["nll"]-tam["deviance"]/2)/reference.n,
                    coordinates=float(np.max(abs(point-np.array(app_point["point"])))),
                    continuous_nll_per_person=abs(tight["nll"]-app_point["continuous"][-1]["nll"])/reference.n,
                    continuous_score_supnorm=float(np.max(abs(tight["gradient"]))),
                    expanded_joint_se_relative=float(np.max(abs(se/app_se-1))))
                person = pd.DataFrame(tam["person"]).set_index("pid").loc[reference.person_labels]
                errors["returned_eap"] = float(np.max(abs(person["EAP"].to_numpy()-fixed_grid["eap"])))
                errors["returned_sd"] = float(np.max(abs(person["SD.EAP"].to_numpy()-fixed_grid["sd"])))
                refinement = {k:float(np.max(abs(np.array(continuous[0][k])-np.array(tight[k])))) for k in ("nll","gradient","information","eap","sd")}
                limits = tam_protocol["limits"]
                checks = dict(grid_nll=errors["grid_nll_per_person"] < limits["same_point_grid_nll_per_person"],
                    refit_coordinates=errors["coordinates"] < limits["point_difference"],
                    refit_nll=errors["continuous_nll_per_person"] < limits["continuous_nll_per_person"],
                    stationary=errors["continuous_score_supnorm"] < limits["continuous_score_supnorm"],
                    joint_se=errors["expanded_joint_se_relative"] < limits["continuous_joint_se_relative"],
                    returned_moments=max(errors["returned_eap"],errors["returned_sd"]) < 1e-8,
                    reference=refinement["nll"] < 1e-8 and refinement["gradient"] < 1e-8 and refinement["information"] < 1e-7 and max(refinement["eap"],refinement["sd"]) < 1e-9,
                    reference_information_error=tight["numeric_information_error_allowance"] < gh_protocol["limits"]["numeric_information_allowance"],
                    reference_tail=tight["tail_information_abs_bound"] < gh_protocol["limits"]["tail_information"],
                    terminated=tam["terminated_below_cap"])
                native_se = np.array([r["se.xsi"] for r in tam["xsi"]])
                result = dict(id=key, point=point, continuous=continuous, covariance=covariance,
                    expanded=expanded, expanded_labels=labels, expanded_se=se, errors=errors,
                    refinement=refinement, checks=checks, all_checks_pass=all(checks.values()),
                    native_se=native_se, continuous_free_se=np.sqrt(np.diag(covariance)),
                    native_vs_joint_se_relative=float(np.max(abs(native_se/np.sqrt(np.diag(covariance))[:-1]-1))))
                save(args.output/(key+".json"), result)
                records.append({k:result[k] for k in ("id","errors","checks","all_checks_pass","native_vs_joint_se_relative")})
                print(key, "failed checks:", [k for k,v in checks.items() if not v], flush=True)
            except Exception as exc:
                failures.append(dict(id=key, error=repr(exc)))
                save(args.output/"failures.json", failures)
                print(key, "FAILED", repr(exc), flush=True)
    summary = dict(classification="POST_RESULT_DEVELOPMENT_ONLY", scientific_inference_ready=False,
        expected_points=12, completed_points=len(records), cases=records, failures=failures,
        all_checks_pass=len(records)==12 and not failures and all(r["all_checks_pass"] for r in records))
    save(args.output/"summary.json", summary)
    save(args.output/"manifest.json", {p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(args.output.iterdir()) if p.is_file()})
    return 0 if summary["all_checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
