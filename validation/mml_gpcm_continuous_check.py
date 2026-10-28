"""Independent continuous-normal checks at all retained finite-GH GPCM points.

A development accuracy diagnosis, not another fit or a qualification study.
Includes a binary zero-difficulty closed-form marginal-likelihood control.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))

import numpy as np
from validation.mml_joint_information_check import save
from validation.mml_rasch_quadrature_review import verify


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('parent',type=Path);parser.add_argument('output',type=Path)
    args=parser.parse_args();verify(args.parent)
    out=args.output;out.mkdir(parents=True,exist_ok=False)
    original=json.loads((args.parent/'reference_input.json').read_text())
    parent_summary=json.loads((args.parent/'summary.json').read_text())
    settings=[dict(bound=12.,rel_tol=1e-10,abs_tol=1e-12),dict(bound=16.,rel_tol=1e-12,abs_tol=1e-14)]
    limits=dict(reference_nll=1e-8,reference_moments=1e-9,reference_sigma_score=1e-8,
        tail_relative_mass=1e-10,numeric_relative_mass=1e-7,
        gh_nll_per_person=1e-6,gh_eap=1e-4,gh_sd=1e-4,gh_log_sigma_score=1e-4,closed_form=1e-10)
    contract=dict(classification='DEVELOPMENT_ONLY',scientific_inference_ready=False,
        parent_manifest_sha256=hashlib.sha256((args.parent/'manifest.json').read_bytes()).hexdigest(),
        settings=settings,limits=limits,scope='continuous likelihood, posterior moments and density log-SD score; no continuous Hessian',
        planned_cases=parent_summary['planned'],parent_all_checks_pass=parent_summary['all_checks_pass'])
    save(out/'protocol.json',contract)
    (out/'check.py').write_bytes(Path(__file__).read_bytes())
    (out/'reference.R').write_bytes(Path(__file__).with_name('mml_joint_information_reference.R').read_bytes())
    cases=[dict(c,continuous_only=True) for c in original['cases']]
    cases.append(dict(id='binary_zero_difficulty_control',continuous_only=True,nodes=[],weights=[],
        config=dict(model='GPCM',facet_names=['Criterion'],facet_levels={'Criterion':['C1','C2']},
            facet_signs={'Criterion':-1},noncenter_facet='Criterion',step_facet='Criterion',slope_facet='Criterion',
            n_cat=2,rating_min=0,rating_max=1),
        data=[dict(Person='P1',Criterion='C1',Score=0),dict(Person='P2',Criterion='C2',Score=1)],
        point=[0.,0.,float(np.log(2)),float(np.log(1.3))]))
    save(out/'input.json',dict(classification='DEVELOPMENT_ONLY',continuous_settings=settings,cases=cases))
    with (out/'R_execution.log').open('w') as log:
        result=subprocess.run(['/usr/local/bin/Rscript',str(out/'reference.R'),str(out/'input.json'),str(out/'R_results.json')],stdout=log,stderr=subprocess.STDOUT)
    if result.returncode:raise RuntimeError('Continuous R reference failed; artifacts retained')
    reference=json.loads((out/'R_results.json').read_text());records=[];control={}
    for case in reference['cases']:
        loose,tight=case['continuous']
        error=dict(reference_nll=abs(loose['nll']-tight['nll']),
            reference_moments=max(float(np.max(abs(np.array(loose[k])-np.array(tight[k])))) for k in ['eap','sd']),
            reference_sigma_score=abs(loose['density_log_sigma_score']-tight['density_log_sigma_score']),
            tail_relative_mass=tight['tail_relative_mass_bound_sum'],numeric_relative_mass=tight['numeric_relative_mass_error_sum'])
        checks={k:v<limits[k] for k,v in error.items()}
        if case['id']=='binary_zero_difficulty_control':
            control=dict(nll_error=abs(tight['nll']-2*np.log(2)),sigma_score_error=abs(tight['density_log_sigma_score']),
                second_moment_error=float(np.max(abs(np.array(tight['eap'])**2+np.array(tight['sd'])**2-1.3**2))))
            control['passed']=all(v<limits['closed_form'] for v in control.values()) and all(checks.values())
            continue
        point=json.loads((args.parent/(case['id']+'.json')).read_text());persons=len(point['eap'])
        error.update(gh_nll_per_person=abs(point['nll']-tight['nll'])/persons,
            gh_eap=float(np.max(abs(np.array(point['eap'])-np.array(tight['eap'])))),
            gh_sd=float(np.max(abs(np.array(point['sd'])-np.array(tight['sd'])))),
            gh_log_sigma_score=abs(point['gradient'][-1]-tight['density_log_sigma_score']))
        checks.update({k:v<limits[k] for k,v in error.items()})
        records.append(dict(id=case['id'],errors=error,checks=checks,all_checks_pass=all(checks.values())))
    summary=dict(classification='DEVELOPMENT_ONLY',scientific_inference_ready=False,comparisons=records,
        closed_form=control,complete=len(records)==contract['planned_cases'],
        all_checks_pass=len(records)==contract['planned_cases'] and control.get('passed',False) and all(r['all_checks_pass'] for r in records))
    save(out/'summary.json',summary)
    save(out/'manifest.json',{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(out.iterdir()) if p.is_file()})
    print((out/'summary.json').read_text(),flush=True)
    return 0 if summary['all_checks_pass'] else 1


if __name__=='__main__':raise SystemExit(main())
