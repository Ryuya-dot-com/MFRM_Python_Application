#!/usr/bin/env python3
"""Compare two variance-update formulas against every saved GH transition."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import roots_hermitenorm

from mml_pcm_conquest_history_audit import coordinates, finite_gh
from mml_pcm_conquest_mc import SOURCE, SOURCE_HASH, rows
from mml_pcm_continuous_refit import PCMIntegral, PREVIOUS, dump
from run_pcm_conquest_check import ROOT, sha

HISTORY = ROOT/'validation/generated/mml_pcm_conquest_history_audit_20260914'
HISTORY_HASH = '3f79c3cd253d047048a833140b8b9eefeebf8f5553ae6d732535324dfa206e6d'


def variance_updates(eap, sd):
    eap,sd=np.asarray(eap),np.asarray(sd)
    second=float(np.mean(sd**2+eap**2))
    return dict(fixed_zero_mean=second, centered_mean=second-float(eap.mean())**2)


def main():
    assert variance_updates([1,3],[1,1])==dict(fixed_zero_mean=6.,centered_mean=2.)
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,required=True)
    out=parser.parse_args().output_dir.resolve()
    assert sha(HISTORY/'summary.json')==HISTORY_HASH
    assert sha(SOURCE/'review_summary.json')==SOURCE_HASH
    for base in (HISTORY,SOURCE):
        record=json.loads((base/('summary.json' if base==HISTORY else 'review_summary.json')).read_text())
        for name,h in record['artifact_sha256'].items():
            assert sha(base/name)==h,name
    out.mkdir(parents=True,exist_ok=False)
    problem=PCMIntegral(json.loads((PREVIOUS/'input.json').read_text())['data'])
    hashes={str(p.relative_to(ROOT)):sha(p) for p in (Path(__file__).resolve(),
        ROOT/'validation/mml_pcm_conquest_history_audit.py',ROOT/'validation/mml_pcm_continuous_refit.py',
        ROOT/'validation/mml_pcm_conquest_mc.py',ROOT/'validation/run_pcm_conquest_check.py')}
    dump(out/'input.json',dict(classification='OBSERVED_DEVELOPMENT_ONLY',scientific_inference_ready=False,
        qualification_eligible=False,history_summary_sha256=HISTORY_HASH,source_summary_sha256=SOURCE_HASH,
        source_sha256=hashes,posthoc=True,
        question='Does the preceding GH posterior reproduce the next native variance with or without mean centering?',
        scope='All 1,539 saved adjacent transitions; one fixture; rounded coordinates; no native execution'))
    results={}
    for q in (31,61,121,181):
        history=rows(SOURCE/f'cq_gh{q}'/'history.csv')
        z,w=roots_hermitenorm(q); w/=np.sqrt(2*np.pi)
        transitions=[]
        for before,after in zip(history,history[1:]):
            value=finite_gh(problem,coordinates(before),z,w)
            predicted=variance_updates(value['eap'],value['sd'])
            actual=float(after['wvar 1 1'])
            transitions.append(dict(from_iteration=int(before['Iteration']),to_iteration=int(after['Iteration']),
                actual_variance=actual,predicted=predicted,residual={k:actual-v for k,v in predicted.items()}))
        dump(out/f'q{q}.json',transitions)
        results[str(q)]={k:dict(max_absolute_residual=float(max(abs(t['residual'][k]) for t in transitions)),
            median_absolute_residual=float(np.median([abs(t['residual'][k]) for t in transitions]))) for k in predicted}
        results[str(q)]['transitions']=len(transitions)
        print(q,results[str(q)],flush=True)
    assert hashes=={p:sha(ROOT/p) for p in hashes}
    dump(out/'summary.json',dict(classification='OBSERVED_DEVELOPMENT_ONLY',scientific_inference_ready=False,
        qualification_eligible=False,results=results,source_sha256=hashes,
        artifact_sha256={p.name:sha(p) for p in out.iterdir() if p.is_file()},
        limitations=['Post-hoc numerical reconstruction, not source-level proof of ConQuest internals',
            'Residuals include six-decimal export rounding; no fitted error acceptance threshold',
            'Changing native initial values or constraints has not been tested here']))


if __name__=='__main__':
    main()
