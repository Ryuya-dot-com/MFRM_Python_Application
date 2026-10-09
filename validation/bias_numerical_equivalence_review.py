"""Supplementary, explicitly post-hoc review of the retained G3 comparison.

The original blanket comparison remains unchanged, including its failures.
Stationarity residuals at different solver endpoints need not agree as raw
numbers: their change is approximately H times the endpoint displacement.
Evaluate both implementations at BOTH retained endpoints and both Q orders.
This is implementation evidence, never population operating characteristics.
"""
import argparse
import json
from pathlib import Path
import statistics

import numpy as np

from validation.mml_numerical_efficiency_check import digest, module
from validation.mml_population_study import AttemptLedger, json_bytes, prepare_retained_problem
from validation.bias_population_design import prepare_bias_problem
from validation.mfrm_recovery_design import RecoveryDesign
from validation.mml_grouped_fit import normal_rule
from validation.mml_grouped_likelihood import GroupedMmlLikelihood


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('parent',type=Path);parser.add_argument('output',type=Path)
    args=parser.parse_args();parent=args.parent.resolve();out=args.output.resolve()
    out.mkdir(parents=True,exist_ok=False)
    root=Path(__file__).resolve().parents[1]
    plan=json.loads((parent/'plan.json').read_bytes())
    summary=json.loads((parent/'summary.json').read_bytes())
    # Reuse the earlier fixed-point test tolerances, not a relaxed fit tolerance.
    review_plan=dict(classification='POST_HOC_IMPLEMENTATION_REVIEW',scientific_inference_ready=False,
        original_all_checks_pass=summary['all_checks_pass'],parent_summary_sha256=digest(parent/'summary.json'),
        source_manifest_sha256=digest(parent/'source_manifest.json'),
        common_point_tolerance=dict(atol=1e-9,rtol=1e-12),
        original_quadrature_limits=plan['estimator']['quadrature_limits'],
        allowed_original_discrepancy='calibration_assessment/errors/cross_score_supnorm only, with both original acceptance checks true; no tolerance or formal decision is changed',
        derivative_identity='g(new)-g(old) = integral_0^1 H(old+t*(new-old)) dt * (new-old); midpoint approximation is diagnostic only',
        script_sha256=digest(Path(__file__)))
    (out/'plan.json').write_bytes(json_bytes(review_plan))
    (out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    hashes=json.loads((parent/'source_manifest.json').read_bytes())
    for name in ['validation/mml_grouped_likelihood.py','validation/mml_grouped_fit.py',
                 'validation/mml_grouped_likelihood_baseline.py','streamlit_app.py']:
        if digest(root/name)!=hashes[name]:raise ValueError('Compared implementation changed: '+name)
    old_engine=module(root/'validation/mml_grouped_likelihood_baseline.py','review_old_engine').GroupedMmlLikelihood
    ledgers=[AttemptLedger.read_only(parent/(label+'.sqlite'),plan_sha256=digest(parent/'plan.json'))
             for label in ['frozen','candidate']]
    records=[]
    try:
        for record,scenario in zip(summary['records'],plan['scenarios'],strict=True):
            key=record['run_id'];retained=ledgers[0].get(key,'input_retained')
            if retained['scenario_id']!=scenario['id']:raise ValueError('Scenario ordering differs')
            qs=plan['estimator']['quadrature'];spec=RecoveryDesign(**scenario['design'])
            data,_,invalid=prepare_retained_problem(retained,spec.model,qs[0])
            if invalid:raise ValueError(invalid)
            _,problem=prepare_bias_problem(data,spec,qs[0])
            engines=[old_engine(data,problem.config),GroupedMmlLikelihood(data,problem.config)]
            common_checks=[];maximum={};diagnostics=[]
            for fitted_q in qs:
                fits=[l.get(key,f'Q{fitted_q}_finished') for l in ledgers]
                points=[np.asarray(f['point']) for f in fits]
                for evaluated_q in qs:
                    nodes,weights=normal_rule(evaluated_q)
                    endpoint_values=[]
                    for point in points:
                        values=[e.evaluate(point,nodes,weights,information=True) for e in engines]
                        endpoint_values.append(values[1])
                        for name in values[0]:
                            a,b=[np.asarray(v[name]) for v in values]
                            common_checks.append(bool(np.allclose(a,b,**review_plan['common_point_tolerance'])))
                            maximum[name]=max(maximum.get(name,0.),float(np.max(abs(a-b))))
                    delta=points[1]-points[0]
                    middle=engines[1].evaluate((points[0]+points[1])/2,nodes,weights,information=True)
                    observed=endpoint_values[1]['gradient']-endpoint_values[0]['gradient']
                    predicted=middle['information']@delta
                    diagnostics.append(dict(fitted_q=fitted_q,evaluated_q=evaluated_q,
                        endpoint_displacement=float(np.max(abs(delta))),
                        score_change=float(np.max(abs(observed))),
                        midpoint_H_delta_residual=float(np.max(abs(observed-predicted)))))
            comparisons=[l.get(key,'calibration_assessment') for l in ledgers]
            differences=record['screen_differences']
            allowed=all(d['path']=='calibration_assessment/errors/cross_score_supnorm'
                        and d['reason']=='numeric_tolerance' for d in differences)
            same_checks=comparisons[0]['checks']==comparisons[1]['checks']
            original_acceptance=all(c['pass_'] and c['checks']['cross_score_supnorm'] for c in comparisons)
            result=dict(run_id=key,original_all_checks_pass=record['all_checks_pass'],
                original_discrepancies=differences,common_point_checks_pass=all(common_checks),
                common_point_max_errors=maximum,endpoint_gradient_diagnostics=diagnostics,
                original_acceptance_unchanged=same_checks and original_acceptance,
                implementation_equivalence_supported=record['input_match'] and len(record['fits'])==len(qs)
                and all(f['all_checks_pass'] and f['both_numerical_pass'] for f in record['fits'])
                and allowed and same_checks and original_acceptance and all(common_checks))
            records.append(result)
            (out/(scenario['id']+'.json')).write_bytes(json_bytes(result))
            print(key,result['implementation_equivalence_supported'],flush=True)
        result=dict(classification=review_plan['classification'],scientific_inference_ready=False,
            original_all_checks_pass=summary['all_checks_pass'],records=records,
            implementation_equivalence_supported=all(r['implementation_equivalence_supported'] for r in records),
            median_full_pipeline_speed_ratio=statistics.median(r['time_ratio'] for r in summary['records']),
            ledger_chains={label:l.verify() for label,l in zip(['frozen','candidate'],ledgers)})
        if result['ledger_chains']!=summary['ledger_chains']:raise ValueError('Original benchmark ledger changed')
        (out/'summary.json').write_bytes(json_bytes(result))
        return 0 if result['implementation_equivalence_supported'] else 1
    finally:
        for ledger in ledgers:ledger.close()


if __name__=='__main__':raise SystemExit(main())
