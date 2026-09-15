"""Downloadable, offline PCM comparison. Importing this file never runs an engine.

First supported model: three two-level facets, four categories, unit response
weights, fixed zero-mean normal prior, unanchored centered nuisance facets and
an uncentered step facet. No inference qualification or Streamlit result import.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path
import platform
import re
import subprocess
import time

import numpy as np
import pandas as pd
from scipy.special import logsumexp


SCOPE = "PCM; 3 facets × 2 levels; scores 0:3; unit weights; fixed SD; no regression, anchors or penalty; only the step facet uncentered; all four categories observed in each facet combination."
TAM_VERSION = "4.3.25"
CONQUEST_VERSION = "5.47.5"
TAM_ORDER = [0, 1, 2, 3, 4, 6, 5, 7]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def prepare_inputs(data_csv, snapshot_json):
    """Validate the exported model, then recode labels without aggregating ratings."""
    saved = json.loads(snapshot_json)
    cfg, params, levels = saved["config"], saved["params"], saved["levels"]
    require(cfg.get("model") == "PCM" and cfg.get("method") == "MML", "Requires PCM estimated by MML.")
    require(cfg.get("n_cat") == 4 and cfg.get("rating_min") == 0 and cfg.get("rating_max") == 3,
            "Requires the four score categories 0, 1, 2, 3.")
    require(cfg.get("estimate_population_sd") is False and cfg.get("population_model", {}).get("enabled") is False,
            "Requires fixed population SD and zero prior mean, without latent regression.")
    sigma = cfg["population_prior_sd"]
    require(isinstance(sigma, (int, float)) and not isinstance(sigma, bool) and np.isfinite(sigma) and sigma > 0
            and np.isfinite(sigma*sigma) and sigma*sigma > 0,
            "Population SD must be finite and positive.")
    require(cfg.get("facet_regularization_enabled") is False and not cfg.get("dummy_facets"),
            "Penalties and dummy facets are unsupported.")
    facets = cfg["facet_names"]
    c = cfg.get("step_facet")
    require(len(facets) == len(set(facets)) == 3 and c in facets and cfg.get("noncenter_facet") == c,
            "Requires three distinct facets; only the step facet may be uncentered.")
    r, t = [f for f in facets if f != c]
    require(cfg.get("parameterization", {}).get("step_coordinate") == "last_derived_negative_sum",
            "Unsupported step parameterization.")
    free = np.asarray(saved["returned_free_coordinates"], dtype=float)
    require(free.shape == (8,) and np.isfinite(free).all(), "Requires eight finite returned structural coordinates.")
    position, expanded = 0, {}
    for f in facets:
        lev, spec = levels[f], cfg["facet_specs"][f]
        centered = f != c
        require(len(lev) == len(set(lev)) == 2 and all(isinstance(v, str) and v for v in lev), "Each facet needs two distinct labels.")
        require(lev == cfg["facet_levels"][f] == spec["levels"], "Facet level orders disagree.")
        require(spec.get("anchors") == [{"nonfinite": "nan"}]*2 and spec.get("groups") == [None]*2
                and spec.get("group_values") == {} and spec.get("centered") is centered
                and spec.get("origin_constraint") == ("sum_to_zero" if centered else "free")
                and spec.get("n_params") == (1 if centered else 2), "Anchors or unsupported facet constraints are present.")
        count = 1 if centered else 2
        expected = np.r_[free[position], -free[position]] if centered else free[position:position+2]
        expanded[f] = np.asarray(params["facets"][f], dtype=float)
        require(np.array_equal(expanded[f], expected), "Expanded facet estimates disagree with returned coordinates.")
        position += count
    steps = np.column_stack((free[position:].reshape(2, 2), -free[position:].reshape(2, 2).sum(axis=1)))
    require(np.array_equal(steps, np.asarray(params["steps_mat"], dtype=float)), "Expanded PCM steps disagree with returned coordinates.")
    signs = cfg["facet_signs"]
    require(all(type(signs[f]) is int and signs[f] in (-1, 1) for f in facets), "Unsupported facet signs.")
    z = np.r_[-signs[c]*expanded[c], -signs[r]*expanded[r][0], -signs[t]*expanded[t][0], steps[:, :2].ravel()]
    # All 32 category logits, including unobserved person/cell combinations.
    expected = np.array([np.arange(4)*(signs[r]*expanded[r][ri]+signs[t]*expanded[t][ti]+signs[c]*expanded[c][ci])
                         -np.r_[0, np.cumsum(steps[ci])] for ri in range(2) for ti in range(2) for ci in range(2)])
    require(np.max(abs(native_design() @ z - expected)) < 1e-12, "Native and app category logits disagree.")
    data = pd.read_csv(io.BytesIO(data_csv), dtype=str, keep_default_na=False)
    require({"Person", "Score", *facets}.issubset(data), "Required response columns are missing.")
    require(not data.duplicated(["Person", *facets]).any(), "Repeated person/facet cells are unsupported; no aggregation is performed.")
    score = pd.to_numeric(data["Score"], errors="raise").to_numpy(dtype=float)
    require(np.isfinite(score).all() and np.isin(score, [0,1,2,3]).all(), "Responses must be integers 0:3.")
    if "Weight" in data:
        require(np.all(pd.to_numeric(data["Weight"], errors="raise").to_numpy() == 1), "Only unit response weights are supported.")
    people = levels["Person"]
    require(0 < len(people) <= 999999 and len(set(people)) == len(people) and set(data.Person) == set(people),
            "Person labels must be unique in the level map; empty persons are unsupported.")
    require(len(people) == cfg["n_person"], "Person count disagrees with the fitted configuration.")
    codes = {f: data[f].map({v:i for i,v in enumerate(levels[f])}) for f in ["Person", *facets]}
    require(all(v.notna().all() for v in codes.values()), "Unknown person or facet labels.")
    column = (4*codes[r]+2*codes[t]+codes[c]).to_numpy(dtype=int)
    response = np.full((len(people),8), np.nan)
    response[codes["Person"].to_numpy(dtype=int), column] = score
    require(all(set(response[np.isfinite(response[:,j]),j]) == {0,1,2,3} for j in range(8)),
            "All four categories must be observed in each facet combination; native category inference is otherwise unsupported.")
    return dict(schema="mfrm_pcm_native_input_v1", scope=SCOPE, sigma=float(sigma), xsi=z.tolist(),
        responses=[[None if np.isnan(v) else int(v) for v in row] for row in response],
        pid=[f"P{i:06d}" for i in range(len(people))], person_labels=people,
        facet_mapping=dict(criterion=c, rater=r, task=t), facet_levels={f:levels[f] for f in facets},
        source_app_quad_points=cfg["quad_points"], scientific_inference_ready=False)


def native_design():
    a = np.zeros((8,4,8)); k = np.arange(4)
    for r in range(2):
        for t in range(2):
            for c in range(2):
                j = 4*r+2*t+c
                a[j,:,c] = -k
                a[j,:,2] = k*(2*r-1); a[j,:,3] = k*(2*t-1)
                a[j,:,4+2*c:6+2*c] = [[0,0],[-1,0],[-1,-1],[0,0]]
    return a


def evaluate(job, xsi):
    """Independent finite-grid NLL, structural gradient and deterministic scores."""
    theta = np.linspace(-job["bound"],job["bound"],job["nodes"])
    logw = -.5*(theta/job["sigma"])**2-np.log(job["sigma"])-.5*np.log(2*np.pi)+np.log(theta[1]-theta[0])
    y = np.asarray(job["responses"],dtype=float); a = native_design()
    logits = (a @ np.asarray(xsi))[:,None,:]+theta[None,:,None]*np.arange(4)
    logp = logits-logsumexp(logits,axis=2,keepdims=True)
    ll = np.zeros((len(y),len(theta)))
    for j in range(8):
        ok = np.isfinite(y[:,j]); ll[ok] += logp[j,:,y[ok,j].astype(int)]
    lm = logsumexp(ll+logw,axis=1); post = np.exp(ll+logw-lm[:,None])
    eap = post @ theta
    gradient = np.zeros(8)
    for j in range(8):
        ok = np.isfinite(y[:,j])
        gradient += post[ok].sum(axis=0) @ (np.exp(logp[j]) @ a[j]) - a[j,y[ok,j].astype(int)].sum(axis=0)
    return dict(raw_nll=float(-lm.sum()), normalized_nll=float(-lm.sum()+len(y)*logsumexp(logw)),
        gradient=gradient.tolist(), gradient_supnorm=float(max(abs(gradient))), eap=eap.tolist(),
        sd=np.sqrt(np.sum(post*(theta[None,:]-eap[:,None])**2,axis=1)).tolist())


def checked_bundle(folder):
    record = json.loads((folder/"comparison_report.json").read_text())
    require(record.get('schema_version') == 'mfrm_numerical_comparison_v1'
            and record.get('scientific_inference_ready') is False and record.get('qualification_eligible') is False,
            "Unsupported comparison record or qualification claim.")
    hashes = record["details"]["input_sha256"]
    digest = hashlib.sha256(json.dumps(hashes,sort_keys=True,ensure_ascii=False,separators=(',',':'),allow_nan=False).encode()).hexdigest()
    require(digest == record['details']['input_bundle_sha256'], "Input manifest digest mismatch.")
    for name, expected in hashes.items():
        path = folder/name
        require(Path(name).name == name and not path.is_symlink(), "Unsafe input path.")
        require(re.fullmatch(r"[0-9a-f]{64}", expected) is not None and hashlib.sha256(path.read_bytes()).hexdigest() == expected,
                f"Input hash mismatch: {name}")
    require({"data.csv", "analysis_snapshot.json", "run_pcm_native.py", "run_pcm_tam.R"}.issubset(hashes), "Missing native-runner input identities.")
    require(hashlib.sha256(Path(__file__).read_bytes()).hexdigest() == hashes["run_pcm_native.py"], "Running script differs from the exported script.")
    job = prepare_inputs((folder/"data.csv").read_bytes(), (folder/"analysis_snapshot.json").read_text())
    job["input_sha256"] = hashes
    job["input_bundle_sha256"] = digest
    job["exported_analysis_identity"] = record["details"]["analysis_identity"]  # Preserved label, not independently authenticated.
    return job


def write_json(path, value):
    with path.open("x") as f:
        json.dump(value,f,indent=2,ensure_ascii=False,allow_nan=False);f.write("\n")


def conquest_control(job, mode):
    verb = "anchor" if mode == "fixed" else "init"
    return f'''title Exported PCM comparison ({mode});
export logfile >> internal.log;
set lconstraints=cases, sconstraint=none, nodefilter=0, p_nodes=2000, seed=2, keeplastests=yes, progress=no, exit_on_error=yes;
datafile ../native_data.csv ! filetype=csv, header=yes, columnlabels=no, pid=Person, pidwidth=7, responses=Y1 to Y8, facets=criterion(2) task(2) rater(2) >> converted.txt;
codes 0,1,2,3;
labels 1 C1 ! criterion;
labels 2 C2 ! criterion;
labels 1 R1 ! rater;
labels 2 R2 ! rater;
labels 1 T1 ! task;
labels 2 T2 ! task;
model criterion + rater + task + criterion*step;
import {verb}_parameters << ../parameters.txt;
import anchor_covariance << ../variance.txt;
estimate ! method=quadrature, nodes={job['nodes']}, minnode=-{job['bound']:.17g}, maxnode={job['bound']:.17g}, distribution=normal, fit=no, stderr=quick, abilities=eap, matrixout=check, convergence=0.000000001, deviancechange=0.00000000001, iterations={3 if mode == 'fixed' else 2000};
export parameters ! filetype=csv >> parameters.csv;
export amatrix ! filetype=csv >> amatrix.csv;
export reg_coefficients ! filetype=csv >> reg_coefficients.csv;
export covariance ! filetype=csv >> covariance.csv;
show cases ! estimates=eap, filetype=csv >> cases.csv;
write check_history ! filetype=csv >> history.csv;
show parameters ! tables=1:2:3:4:5, estimates=eap >> review.txt;
export scoreddata ! filetype=csv >> scored.csv;
quit;
'''


def run_process(command, folder, *, input_text=None):
    started = time.monotonic()
    with (folder/"console.log").open("x") as log:
        try:
            run = subprocess.run(command,input=input_text,text=True,cwd=folder,stdout=log,stderr=subprocess.STDOUT,timeout=600)
            status = dict(exit_code=run.returncode,timed_out=False)
        except subprocess.TimeoutExpired:
            status = dict(exit_code=None,timed_out=True)
    status["elapsed_seconds"] = time.monotonic()-started
    write_json(folder/"launch.json",status)
    require(status["exit_code"] == 0 and not status["timed_out"], f"Execution failed; inspect {folder/'console.log'}")


def review_outputs(job, folder, engine, mode):
    """Verify native coding/design/anchors; keep score MC separate from quadrature."""
    y = np.asarray(job["responses"],dtype=float); a = native_design()
    if engine == "tam":
        out = json.loads((folder/"result.json").read_text())
        require(out["runtime"]["TAM"] == TAM_VERSION, "Unreviewed TAM version.")
        require(out["xsi_names"] == ['C1','C2','raterR1','taskT1','C1:step1','C2:step1','C1:step2','C2:step2'], "TAM parameter labels changed.")
        z = np.asarray(out["xsi"])[TAM_ORDER]; variance = float(np.asarray(out["variance"]).ravel()[0])
        beta = np.asarray(out["beta"]).ravel()
        reported = out["deviance"]/2
        require(len(out["item_names"]) == 8 and np.array_equal(np.asarray(out['B']),
                np.broadcast_to(np.arange(4)[None,:,None], (8,4,1))), "TAM slope design changed.")
        order = []
        for j,name in enumerate(out["item_names"]):
            c,r,t = map(int,re.fullmatch(r'C(\d)-raterR(\d)-taskT(\d)',name).groups()); k = 4*(r-1)+2*(t-1)+c-1
            order.append(k)
            require(np.array_equal(np.asarray(out["resp"],dtype=float)[:,j],y[:,k],equal_nan=True), "TAM responses changed.")
            require(np.array_equal(np.asarray(out["A"])[j][:,TAM_ORDER],a[k]), "TAM category design changed.")
        require(sorted(order) == list(range(8)), "TAM omitted or duplicated a facet combination.")
        people = pd.DataFrame(out["person"])
        native_eap, native_sd = people.EAP.to_numpy(), people['SD.EAP'].to_numpy()
        require(people.pid.tolist() == job["pid"], "TAM person order changed.")
        anchor_tolerance = 2e-10; criterion = 'raw_nll'
        history = out["deviance_history"]
        last = None  # Native history lacks retained coordinate rows in this runner.
    else:
        console = (folder/"console.log").read_text()
        require(f"ConQuest version: {CONQUEST_VERSION}" in console and "End of Program" in console, "Unreviewed or incomplete ConQuest execution.")
        pars = pd.read_csv(folder/"parameters.csv")
        require(pars.Label.tolist() == [' criterion C1',' criterion C2',' rater R1',' task T1',
                ' criterion C1 category 1',' criterion C1 category 2',' criterion C2 category 1',' criterion C2 category 2'], "ConQuest parameter labels changed.")
        z = pars.Estimate.to_numpy(); variance = pd.read_csv(folder/"covariance.csv").Covariance.iloc[0]
        beta = pd.read_csv(folder/"reg_coefficients.csv").Estimate.to_numpy()
        scored = pd.read_csv(folder/"scored.csv",keep_default_na=False)
        require(scored.PID.tolist() == job["pid"], "ConQuest person order changed.")
        headers = list(scored)[2:]; order = []
        for name in headers:
            c,r,t = map(int,re.match(r'criterion:(\d).*rater:(\d).*task:(\d)',name).groups()); k = 4*(r-1)+2*(t-1)+c-1;order.append(k)
            actual = pd.to_numeric(scored[name].replace({'.':np.nan,'':np.nan,'NA':np.nan})).to_numpy()
            require(np.array_equal(actual,y[:,k],equal_nan=True), "ConQuest responses changed.")
        require(sorted(order) == list(range(8)), "ConQuest omitted or duplicated a facet combination.")
        design = pd.read_csv(folder/"amatrix.csv")
        require(len(design) == 32 and set(zip(design.GIN,design.Category)) == {(j,k) for j in range(1,9) for k in range(1,5)},
                "Incomplete ConQuest design.")
        for _,row in design.iterrows():
            require(np.array_equal(row.iloc[2:].to_numpy(dtype=float),a[order[int(row.GIN)-1],int(row.Category)-1]), "ConQuest category design changed.")
        h = pd.read_csv(folder/"history.csv")
        last = [float(h.iloc[-1][f'xsi {i}']) for i in range(1,9)]
        require(np.array_equal(z,last), "ConQuest did not return the last saved coordinate row.")
        history = h.to_dict('records'); reported = float(h.LogLikelihood.iloc[-1])/2
        people = pd.read_csv(folder/"cases.csv")
        require(people.PID.tolist() == job["pid"], "ConQuest score order changed.")
        native_eap,native_sd = people.EAP_1.to_numpy(),people.PosteriorSD_1.to_numpy()
        anchor_tolerance = 5.1e-7; criterion = 'normalized_nll'
    require(np.isfinite(z).all() and z.shape == (8,) and np.max(abs(beta)) == 0, "Invalid calibration or nonzero prior mean.")
    require(abs(variance-job["sigma"]**2) <= anchor_tolerance, "Fixed variance changed.")
    point = job["xsi"] if mode == "fixed" else z.tolist()
    value = evaluate(job,point)
    if mode == "fixed":
        require(max(abs(z-job["xsi"])) <= anchor_tolerance, "Fixed calibration changed.")
        tolerance = (1e-7 if engine == 'tam' else 2e-6)+1e-9*len(y)
        require(abs(reported-value[criterion]) < tolerance, "Fixed-calibration likelihood mismatch.")
        if engine == "tam":
            require(max(abs(native_eap-value['eap'])) < 1e-8 and max(abs(native_sd-value['sd'])) < 1e-8, "TAM fixed-calibration quadrature scores differ.")
    return dict(engine=engine,mode=mode,initial=job["xsi"],returned=z.tolist(),last_saved=last,
        exported_variance=float(variance),reference_variance=job['sigma']**2,
        fixed_calibration_checked=mode=='fixed',reported_nll=reported,
        reported_nll_point='fixed' if mode=='fixed' else 'native history timing; not assumed to be returned',
        recalculated=value,native_eap_max_difference=float(max(abs(native_eap-value['eap']))),
        native_sd_max_difference=float(max(abs(native_sd-value['sd']))),
        native_score_method='quadrature' if engine=='tam' else 'Monte Carlo: p_nodes=2000, seed=2; stability unassessed',
        history=history,scientific_inference_ready=False,qualification_eligible=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundle',type=Path,default=Path(__file__).resolve().parent)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--conquest',type=Path,required=True)
    parser.add_argument('--nodes',type=int,default=401)
    parser.add_argument('--bound',type=float,default=20.)
    parser.add_argument('--refit',action='store_true',help='Refit only after both fixed-calibration checks pass.')
    args = parser.parse_args(); bundle = args.bundle.resolve(); out = args.output.resolve()
    require(args.nodes >= 3 and args.nodes % 2 == 1 and np.isfinite(args.bound) and args.bound > 0,
            "Use an odd node count >=3 and a finite positive bound.")
    job = checked_bundle(bundle)
    job.update(nodes=args.nodes,bound=args.bound,rule='equally spaced density-times-spacing grid; not Gauss-Hermite',
               scoring_nodes=2000,scoring_seed=2,TAM_version=TAM_VERSION,ConQuest_version=CONQUEST_VERSION)
    binary = args.conquest.resolve(); require(binary.is_file(), "ConQuest executable is missing.")
    job['conquest_binary_sha256'] = hashlib.sha256(binary.read_bytes()).hexdigest()
    out.mkdir(parents=True,exist_ok=False)
    retained = out/'bundle'; retained.mkdir()
    for name, expected in job['input_sha256'].items():
        data = (bundle/name).read_bytes()
        require(hashlib.sha256(data).hexdigest() == expected, f"Input changed while preparing: {name}")
        (retained/name).write_bytes(data)
    (retained/'comparison_report.json').write_bytes((bundle/'comparison_report.json').read_bytes())
    bundle = retained
    write_json(out/'input.json',job)
    frame = pd.DataFrame(job['responses'],columns=[f'Y{i}' for i in range(1,9)])
    frame.insert(0,'Person',job['pid']);frame.to_csv(out/'native_data.csv',index=False,na_rep='.')
    (out/'parameters.txt').write_text(''.join(f'{i+1} {v:.17g}\n' for i,v in enumerate(job['xsi'])))
    (out/'variance.txt').write_text(f"1 1 {job['sigma']**2:.17g}\n")
    command = [str(binary)]
    if platform.system() == 'Darwin' and platform.machine() == 'arm64':
        command = ['/usr/bin/arch','-x86_64',*command]
    results = {}; failure = None
    try:
        for mode in (['fixed','refit'] if args.refit else ['fixed']):
            for engine in ['tam','conquest']:
                folder = out/f'{engine}_{mode}'; folder.mkdir()
                if engine == 'tam':
                    run_process(['Rscript',str(bundle/'run_pcm_tam.R'),str(out/'input.json'),mode],folder)
                else:
                    control = conquest_control(job,mode); (folder/'model.cqc').write_text(control)
                    run_process(command,folder,input_text=control)
                results[f'{engine}_{mode}'] = review_outputs(job,folder,engine,mode)
                write_json(folder/'review.json',results[f'{engine}_{mode}'])
                print(engine,mode,'output checks complete',flush=True)
    except Exception as exc:
        failure = f'{type(exc).__name__}: {exc}'
    write_json(out/'execution_review.json',dict(results=results,failure=failure,
        checks_complete=failure is None,scientific_inference_ready=False,qualification_eligible=False,
        source_sha256={n:hashlib.sha256((bundle/n).read_bytes()).hexdigest() for n in ['run_pcm_native.py','run_pcm_tam.R']},
        artifact_sha256={str(p.relative_to(out)):hashlib.sha256(p.read_bytes()).hexdigest() for p in out.rglob('*') if p.is_file()}))
    require(failure is None, failure)
    print('Saved raw outputs and numerical checks:',out,'; four-question interpretation remains separate.')


if __name__ == '__main__':
    main()
