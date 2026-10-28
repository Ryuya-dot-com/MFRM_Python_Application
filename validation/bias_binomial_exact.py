"""Exact finite-sample controls for the binary conditional bias screen.

Enumerate all binomial counts: no random data or Monte Carlo uncertainty.
The same-data intercept control shows complete absorption of a cell effect,
not a general approximation to the MFRM fitted-calibration distribution.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))

import numpy as np
from scipy.special import expit,logit
from scipy.stats import binom
import streamlit_app as app


def enumerate_cell(n,baseline_probability,shift):
    if type(n) is not int or n<2 or not 0<baseline_probability<1 or not np.isfinite(shift):
        raise ValueError('Integer sample size, interior probability and finite shift required')
    baseline=float(logit(baseline_probability));actual=float(expit(baseline+shift))
    rows=[]
    for successes in range(n+1):
        scores=np.r_[np.ones(successes),np.zeros(n-successes)]
        mass=float(binom.pmf(successes,n,actual));results={}
        for calibration in ['known_baseline','same_data_intercept']:
            if calibration=='same_data_intercept' and successes in [0,n]:
                results[calibration]=dict(status='calibration_has_no_finite_mle',p=None,estimate=None,
                    interval_available=False,coverage_of_DGP_shift=None);continue
            eta=baseline if calibration=='known_baseline' else float(logit(successes/n))
            nll=lambda b:float(n*np.logaddexp(0.,eta+b)-successes*(eta+b))
            probability=lambda b:np.tile([expit(-eta-b),expit(eta+b)],(n,1))
            fit=app._fit_conditional_bias_cell(nll,probability,scores,np.ones(n),np.ones(n),10.)
            interval=fit['ci'];available=interval['status']=='ok'
            value={k:fit[k] for k in ['status','estimate','se','lr','p','reference_status']}
            value={k:None if isinstance(v,float) and not np.isfinite(v) else v for k,v in value.items()}
            value.update(interval_available=available,
                interval_lower=interval['lower'] if np.isfinite(interval['lower']) else None,
                interval_upper=interval['upper'] if np.isfinite(interval['upper']) else None,
                # Only the known-baseline scalar has the DGP shift as target.
                coverage_of_DGP_shift=bool(interval['lower']<=shift<=interval['upper'])
                    if available and calibration=='known_baseline' else None)
            results[calibration]=value
        rows.append(dict(successes=successes,probability=mass,results=results))
    summaries={}
    for calibration in ['known_baseline','same_data_intercept']:
        availability=sum(r['probability'] for r in rows if r['results'][calibration]['p'] is not None)
        rejection=sum(r['probability'] for r in rows if r['results'][calibration]['p'] is not None and r['results'][calibration]['p']<.05)
        ci_available=sum(r['probability'] for r in rows if r['results'][calibration]['interval_available'])
        covered=sum(r['probability'] for r in rows if r['results'][calibration]['coverage_of_DGP_shift'] is True)
        summaries[calibration]=dict(test_availability=availability,test_unavailability=1-availability,
            rejection_per_generated_sample=rejection,rejection_given_available=rejection/availability if availability else None,
            interval_availability=ci_available,
            coverage_per_generated_sample=covered if calibration=='known_baseline' else None,
            coverage_given_interval=covered/ci_available if calibration=='known_baseline' and ci_available else None,
            interpretation='DGP shift is the scalar target' if calibration=='known_baseline' else
                'The fitted intercept absorbs the DGP shift; this residual screen has no separate identifiable effect parameter')
    return dict(n=n,baseline_probability=baseline_probability,DGP_shift=shift,response_probability=actual,
        probability_sum=sum(r['probability'] for r in rows),rows=rows,summaries=summaries,
        calculation='exact binomial enumeration up to floating-point evaluation; no Monte Carlo error',scientific_inference_ready=False)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('output',type=Path);args=parser.parse_args()
    out=args.output;out.mkdir(parents=True,exist_ok=False)
    def write(name,value):
        with (out/name).open('x') as f:json.dump(value,f,indent=2,allow_nan=False);f.write('\n')
    write('protocol.json',dict(classification='ANALYTIC_ENUMERATION_DEVELOPMENT',sample_sizes=[4,8,16,32,64],
        probabilities=[.1,.3,.5,.8],shifts=[0.,.8],alpha=.05,max_abs=10.,
        all_counts='0..n; unavailable extreme-cell fits retained in exact probability denominator',
        qualification_cohort=False,scientific_inference_ready=False))
    write('source_manifest.json',{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
        for p in [Path(__file__),ROOT/'streamlit_app.py']})
    records=[]
    for n in [4,8,16,32,64]:
        for probability in [.1,.3,.5,.8]:
            for shift in [0.,.8]:
                result=enumerate_cell(n,probability,shift);key=f'n{n}_p{round(100*probability)}_b{round(100*shift)}'
                write(key+'.json',result);records.append({k:v for k,v in result.items() if k!='rows'})
    write('summary.json',dict(classification='ANALYTIC_ENUMERATION_DEVELOPMENT',cases=records,scientific_inference_ready=False))
    write('manifest.json',{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir() if p.is_file()})
    print('Exact count enumeration complete:',len(records),'conditions',flush=True)


if __name__=='__main__':main()
