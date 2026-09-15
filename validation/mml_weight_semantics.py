#!/usr/bin/env python3
"""Distinguish response powers from native case weights at fixed PCM calibration."""
import argparse
import copy
import csv
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import traceback
from types import SimpleNamespace

import numpy as np
from scipy.special import logsumexp

from mml_observed_pcm_refit import ObservedPCMReference
from mml_pcm_continuous_refit import ROOT, dump, sha
from mml_pcm_conquest_mc import rows
from mml_unit_weight_native import CELLS, grid as prior_grid
from mml_unit_weight_conquest_mc import PRIOR, verify_prior
from run_pcm_conquest_check import launch

MFRMR = Path('/Users/tohokusla/Dropbox/MFRM_Application/mfrmr/development')
WEIGHTS = dict(unit=np.ones(32), pattern=np.tile([.5, 1., 1.5, 2.], 8), triple=np.tile([1.5, 3., 4.5, 6.], 8))


def self_check():
    # Same-person repeated evidence and two weighted copies of a person differ.
    likelihood = np.array([.2, .8]); prior = np.array([.5, .5])
    assert abs(prior @ likelihood**2 - .34) < 1e-15
    assert abs((prior @ likelihood)**2 - .25) < 1e-15
    assert abs((prior @ likelihood**2) - (prior @ likelihood)**2) > .08


def verify_sources(out, spec):
    manifest_path=out/'source_snapshot/manifest.json'
    manifest=json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    for filename, expected in spec['source_sha256'].items():
        path=Path(filename)
        if sha(path)!=expected and filename in manifest:
            path=out/manifest[filename]['path']
        assert sha(path)==expected,filename


def reference(input_data, response_weights):
    idx = copy.deepcopy(input_data['indices']); idx['weight'] = response_weights
    ref = ObservedPCMReference(SimpleNamespace(idx=idx, config=dict(n_person=32, population_model=dict(X=np.array(input_data['x'])[:, None]))))
    p = np.array(input_data['continuous_point']); sigma = np.exp(p[-1])
    theta = np.linspace(-12, 12, 241)
    conditional_point = p.copy(); conditional_point[0] = 0  # Physical theta already includes the regression mean.
    ll = np.array([ref.conditional(conditional_point, v/sigma)[0] for v in theta]).T
    logw = -.5*((theta[None, :]-ref.x[:, None]*p[0])/sigma)**2 - np.log(sigma) - .5*np.log(2*np.pi) + np.log(.1)
    mass = logsumexp(ll+logw, axis=1); post = np.exp(ll+logw-mass[:, None]); eap = post @ theta
    result = dict(person_nll=-mass, log_prior_mass=logsumexp(logw, axis=1), eap=eap,
                  sd=np.sqrt(np.sum(post*(theta[None, :]-eap[:, None])**2, axis=1)))
    if np.all(np.asarray(response_weights)==1):
        old = prior_grid(ref, p, 12)
        assert abs(sum(result['person_nll'])-old['nll']) < 1e-10
        assert max(np.max(abs(result[k]-old[k])) for k in ('eap', 'sd')) < 1e-10
    return result


def compare(out, spec):
    results = {}
    for cell in CELLS:
        folder = out/cell; data = json.loads((folder/'input.json').read_text())
        scalar = json.loads((folder/'scalar_R.json').read_text())
        unit = reference(data, np.ones(len(data['indices']['person'])))
        results[cell] = {}
        for scheme, weights in WEIGHTS.items():
            run = folder/scheme; native = run/'cq'; old = PRIOR/cell/'cq_b12_fixed'
            launched = json.loads((native/'launch.json').read_text())
            assert launched['exit_code']==0 and not launched['timed_out'] and launched['end_of_program']
            for filename in ('parameters.csv', 'reg_coefficients.csv', 'covariance.csv', 'amatrix.csv', 'scored.csv'):
                assert sha(native/filename)==sha(old/filename), native/filename
            text = (native/'converted.txt').read_text().splitlines()
            assert len(text)==32
            old_text = (old/'converted.txt').read_text().splitlines()
            assert all(line[:36]==old_text[i][:36] and abs(float(line[36:60])-weights[i])<1e-15 for i, line in enumerate(text))
            report = (native/'review.txt').read_text()
            assert int(re.search(r'Total number of estimated parameters:\s*(-?\d+)', report)[1])==0
            assert int(re.search(r'Number of nodes used when drawing PVs:\s*(\d+)', report)[1])==2000
            assert float(re.search(r'Random number generation seed:\s*([\d.]+)', report)[1])==2
            history = rows(native/'history.csv'); cq_nll = float(history[-1]['LogLikelihood'])/2
            # All history calibration columns must match; only the objective may change.
            baseline = rows(old/'history.csv')[-1]
            assert len(history)==1
            assert all(history[0][key]==baseline[key] for key in baseline if key!='LogLikelihood')
            native_scores = rows(native/'cases.csv'); previous_scores = rows(old/'cases.csv')
            assert [r['PID'] for r in native_scores]==[f'P{i:03d}' for i in range(32)]
            cq_score_change = {key:float(max(abs(float(a[col])-float(b[col])) for a,b in zip(native_scores, previous_scores)))
                               for key,col in [('eap','EAP_1'), ('sd','PosteriorSD_1')]}
            tam = json.loads((run/'tam.json').read_text()); assert all(tam['checks'].values())
            scaled = 32*weights/sum(weights)
            raw = float(weights @ unit['person_nll'])
            normalized = float(weights @ (unit['person_nll']+unit['log_prior_mass']))
            expected_tam = float(scaled @ unit['person_nll'])
            cq_candidates = dict(raw_case_weights=normalized, mean_one_case_weights=normalized*32/sum(weights))
            cq_errors = {key:abs(cq_nll-v) for key,v in cq_candidates.items()}
            matches = [key for key,v in cq_errors.items() if v<spec['limits']['cq_nll']]
            predeclared_match = len(matches)==(2 if scheme=='unit' else 1)
            # Post-hoc diagnosis after the original weighted-NLL hypotheses failed.
            # Preserve their failed status; unweighted display is a separate finding.
            unweighted_nll=float(sum(unit['person_nll']+unit['log_prior_mass']))
            amatrix=np.zeros((8,4,8))
            for row in rows(native/'amatrix.csv'):
                amatrix[int(row['GIN'])-1,int(row['Category'])-1]=[float(row[k]) for k in list(row)[2:]]
            observed=np.zeros((32,8))
            for i,row in enumerate(rows(native/'scored.csv')):
                for gin,key in enumerate(list(row)[2:]):
                    if row[key] not in ('','NA','.'):
                        observed[i]+=amatrix[gin,int(row[key])]
            section=(native/'internal.log').read_text().split('Original Sufficient Statistics\n')[1].split('Modified Sufficient Statistics')[0]
            parsed=re.findall(r'^\s*(\d+)\s+(-?[\d.]+)\s*/\*',section,re.M)
            assert [int(i) for i,_ in parsed]==list(range(1,9))
            statistics=np.array([float(v) for _,v in parsed])
            statistic_errors=dict(raw=float(max(abs(statistics-weights@observed))),
                mean_one=float(max(abs(statistics-scaled@observed))),unit=float(max(abs(statistics-observed.sum(axis=0)))))
            assert abs(cq_nll-unweighted_nll)<spec['limits']['cq_nll'] and statistic_errors['mean_one']<2e-5
            assert float(re.search(r'Weighted number of cases in MML/MCMC estimation:\s*([\d.]+)',report)[1])==sum(weights)
            assert abs(tam['deviance']/2-expected_tam)<spec['limits']['tam_nll']
            assert max(tam['score_change_to_previous'].values())<spec['limits']['tam_score']
            inner = reference(data, weights[np.array(data['indices']['person'])])
            r_key = 'unit' if scheme=='unit' else scheme
            differences = {key:float(np.max(abs(np.array(scalar[r_key][key])-inner[key]))) for key in inner}
            assert max(differences.values())<spec['limits']['independent_R']
            wrong = dict(nll=float(sum(inner['person_nll'])), normalized_nll=float(sum(inner['person_nll']+inner['log_prior_mass'])),
                         max_eap_change=float(max(abs(inner['eap']-unit['eap']))), max_sd_change=float(max(abs(inner['sd']-unit['sd']))))
            if scheme!='unit':
                assert abs(wrong['normalized_nll']-normalized)>1e-3 and wrong['max_eap_change']>1e-3
            results[cell][scheme] = dict(raw_case_weight_nll=raw, normalized_prior_case_weight_nll=normalized,
                expected_tam_nll=expected_tam, tam_reported_nll=tam['deviance']/2, tam_nll_error=abs(tam['deviance']/2-expected_tam),
                cq_reported_nll=cq_nll, cq_candidates=cq_candidates, cq_candidate_errors=cq_errors, cq_matching_candidates=matches,
                predeclared_cq_hypothesis_pass=predeclared_match, cq_unweighted_nll_error=abs(cq_nll-unweighted_nll),
                cq_observed_sufficient_statistics=statistics, cq_sufficient_statistic_errors=statistic_errors,
                weighted_log_prior_mass=float(weights @ unit['log_prior_mass']), tam_score_change=tam['score_change_to_previous'],
                cq_score_change=cq_score_change, tam_iterations=tam['iter'], tam_reached_cap=tam['reached_cap'],
                cq_termination=[line.strip() for line in report.splitlines() if 'Iterations terminated' in line],
                independent_R_difference=differences, inner_response_weight_reference=inner, outer_unit_reference=unit,
                wrong_weight_substitution=wrong)
    return results


def main():
    self_check(); parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True); parser.add_argument('--audit', action='store_true')
    parser.add_argument('--review-retained',action='store_true',help='Diagnose the retained first hypothesis failure without launching engines')
    args=parser.parse_args(); out=args.output_dir.resolve(); _, old_spec=verify_prior()
    if args.audit:
        spec=json.loads((out/'protocol.json').read_text()); summary=json.loads((out/'review_summary.json').read_text())
        for f,h in summary['artifact_sha256'].items(): assert sha(out/f)==h,f
        verify_sources(out,spec)
        review_spec=json.loads((out/'review_protocol.json').read_text())
        assert sha(Path(__file__).resolve())==review_spec['review_source_sha256']
        assert sha(out/'protocol.json')==summary['protocol_sha256']
        replay=json.loads(json.dumps(compare(out,spec),default=lambda a:a.tolist()))
        assert replay==summary['results']
        print('Saved audit passed:',len(summary['artifact_sha256']),'artifacts;',len(spec['source_sha256']),'sources; 24 native calls')
        return
    if args.review_retained:
        spec=json.loads((out/'protocol.json').read_text());verify_sources(out,spec)
        assert (out/'failure.json').is_file() and (out/'source_snapshot/manifest.json').is_file()
        for f,h in json.loads((out/'prepared.json').read_text())['artifact_sha256'].items(): assert sha(out/f)==h,f
        dump(out/'review_protocol.json',dict(classification='POST_HOC_NUMERICAL_DIAGNOSIS',
            original_failure_sha256=sha(out/'failure.json'),original_protocol_sha256=sha(out/'protocol.json'),
            review_source_sha256=sha(Path(__file__).resolve()),source_snapshot_manifest_sha256=sha(out/'source_snapshot/manifest.json'),
            question='After both predeclared weighted-NLL hypotheses failed, does the fixed ConQuest display match unweighted NLL while its original sufficient statistics match mean-one case weights?',
            limits=dict(unweighted_display_nll=spec['limits']['cq_nll'],sufficient_statistics=2e-5),
            scope='Saved-output reanalysis only, no new native runs; never promote the failed predeclared hypotheses to passed'))
        results=compare(out,spec)
        dump(out/'review_summary.json',dict(classification='POST_HOC_NUMERICAL_DIAGNOSIS',scientific_inference_ready=False,
            qualification_eligible=False,arithmetic_and_native_input_checks_pass=True,
            predeclared_cq_hypotheses_pass=all(v['predeclared_cq_hypothesis_pass'] for r in results.values() for v in r.values()),
            results=results,protocol_sha256=sha(out/'protocol.json'),
            artifact_sha256={str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file()}))
        print('Retained diagnosis complete; predeclared ConQuest weighted-display hypotheses remain failed')
        return
    source_paths=[Path(__file__).resolve(),Path(__file__).with_suffix('.R').resolve(), ROOT/'validation/mml_unit_weight_conquest_mc.py',
        MFRMR/'DESCRIPTION',MFRMR/'R/mfrm_core.R',MFRMR/'R/core-data-prep.R',MFRMR/'src/mml_backend.cpp']
    hashes={str(ROOT/f):h for f,h in old_spec['source_sha256'].items()}
    hashes.update({str(p):sha(p) for p in source_paths})
    out.mkdir(parents=True,exist_ok=False)
    spec=dict(classification='OBSERVED_DEVELOPMENT_ONLY',scientific_inference_ready=False,qualification_eligible=False,
        question='At fixed PCM calibration, are native person weights applied outside the integral and rescaled? Can even constant-within-person response weights be substituted?',
        cells=CELLS,case_weights=WEIGHTS,prior=str(PRIOR),prior_summary_sha256=sha(PRIOR/'summary.json'),
        bound=12,spacing=.1,tam_version='4.3.25',conquest_version='5.47.5',mfrmr_source_version='0.2.4.9000',
        design='Three positive case-weight vectors (sum 32, 40, 120), four retained response datasets, both engines: 24 fixed calls; no refits. The latter two vectors differ by global factor 3. No new responses.',
        hypotheses='Replay the original raw-vs-mean-one weighted ConQuest display hypotheses, retaining any rejection. Also check the unweighted display and mean-one sufficient-statistic explanation identified from the first retained failure. ConQuest finite prior normalized per person, TAM raw density-times-spacing.',
        negative_control='Independently integrate the same values as powers of each conditional response likelihood. Even constant within person, this differs from weighting the marginal log likelihood. Not a native response-weight fit.',
        limits=dict(cq_nll=2e-6,tam_nll=1e-7,tam_score=1e-8,independent_R=1e-8),
        source_sha256=hashes,binary_sha256=old_spec['binary_sha256'],
        zotero=dict(query='ConQuest',key='SKHZWZD5',title='ACER ConQuest Manual',year=None,other_queries={'pseudo-likelihood':[], 'sampling weights':[]},scope='Metadata search only; mathematical/software claims use official manual and inspected implementation.'),
        sources=['https://conquestmanual.acer.org/s4-00.html','https://alexanderrobitzsch.github.io/TAM/reference/tam.mml.html'],
        limitations=['Two response sets paired across missingness, not four independent replications',
          'Positive person weights only; zero case weights and arbitrary varying within-person native response weights are not tested',
          'Finite-range semantics audit, not continuous accuracy, optimizer, native MC scoring, SE/CI/coverage or survey inference qualification',
          'mfrmr source inspection only; no mfrmr fit and no change to its quadrature-zero finding',
          'No automatic conversion of response weights to person weights; preserve all failed attempts without replacement'])
    dump(out/'protocol.json',spec)
    for cell in CELLS:
        folder=out/cell;folder.mkdir();data=json.loads((PRIOR/cell/'input.json').read_text());dump(folder/'input.json',data)
        for scheme,w in WEIGHTS.items():
            run=folder/scheme;run.mkdir();cq=run/'cq';cq.mkdir()
            for name in ('anchor_parameters.txt','anchor_beta.txt','anchor_covariance.txt'):shutil.copyfile(PRIOR/cell/name,run/name)
            table=rows(PRIOR/cell/'data.csv')
            with (run/'data.csv').open('x',newline='') as f:
                writer=csv.DictWriter(f,fieldnames=[*table[0],'W']);writer.writeheader()
                writer.writerows(dict(row,W=format(w[i],'.17g')) for i,row in enumerate(table))
            command=(PRIOR/cell/'cq_b12_fixed/model.cqc').read_text()
            assert command.count('keeps=x, keepswidth=24')==1 and command.count('regression x;')==1
            command=command.replace('keeps=x, keepswidth=24','keeps=x W, keepswidth=24').replace('regression x;','regression x;\ncaseweight W;')
            (cq/'model.cqc').write_text(command)
    dump(out/'prepared.json',dict(artifact_sha256={str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file()}))
    with (out/'tam.log').open('x') as log:
        status=subprocess.run(['Rscript',str(Path(__file__).with_suffix('.R')),str(out)],stdout=log,stderr=subprocess.STDOUT)
    assert status.returncode==0,out/'tam.log'
    for cell in CELLS:
        for scheme in WEIGHTS:
            cq=out/cell/scheme/'cq';result=launch(cq,(cq/'model.cqc').read_text(),600)
            assert result['exit_code']==0 and result['end_of_program'],cq
            print(cell,scheme,'ConQuest completed',flush=True)
    results=compare(out,spec);verify_prior()
    for f,h in hashes.items():assert sha(Path(f))==h,f
    for f,h in json.loads((out/'prepared.json').read_text())['artifact_sha256'].items():assert sha(out/f)==h,f
    dump(out/'summary.json',dict(classification=spec['classification'],scientific_inference_ready=False,qualification_eligible=False,
        arithmetic_and_native_input_checks_pass=True,
        predeclared_cq_hypotheses_pass=all(v['predeclared_cq_hypothesis_pass'] for r in results.values() for v in r.values()),
        results=results,protocol_sha256=sha(out/'protocol.json'),
        artifact_sha256={str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file()}))
    for cell,runs in results.items():
        for scheme,v in runs.items():print(cell,scheme,'CQ matches',v['cq_matching_candidates'],'TAM NLL error',v['tam_nll_error'],'wrong EAP change',v['wrong_weight_substitution']['max_eap_change'],flush=True)


if __name__=='__main__':
    try:main()
    except Exception as exc:
        if '--output-dir' in sys.argv and '--audit' not in sys.argv and '--review-retained' not in sys.argv:
            folder=Path(sys.argv[sys.argv.index('--output-dir')+1])
            if folder.is_dir() and not (folder/'summary.json').exists():dump(folder/'failure.json',dict(error=repr(exc),traceback=traceback.format_exc()))
        raise
