"""Paired, fixed-point development benchmark; no fits or sampling claims.

Run with the candidate and a retained baseline in the import path. The fixed
matrix is written before generating benchmark-only inputs. Neither input nor
results come from an in-progress population study. Timing on a shared machine
is diagnostic, and tracemalloc peaks are not whole-process RSS measurements.
"""
import argparse
from dataclasses import asdict
import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import statistics
import sys
import time
import tracemalloc

import numpy as np
from scipy.special import roots_hermitenorm

from validation.mfrm_recovery_design import RecoveryDesign, generate_population_sample
from validation.mml_grouped_likelihood import GroupedMmlLikelihood


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('baseline', type=Path); parser.add_argument('output', type=Path)
    args = parser.parse_args(); out = args.output.resolve(); out.mkdir(parents=True, exist_ok=False)
    spec = importlib.util.spec_from_file_location('retained_grouped_baseline', args.baseline)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    matrix = [dict(model=m, assignment=a, persons=n, categories=k, q=q)
        for m in ['RSM', 'PCM', 'GPCM'] for a in ['regular', 'sparse', 'unbalanced']
        for n, k, q in [(48,2,31), (192,4,641), (192,6,1281)]]
    plan = dict(classification='DEVELOPMENT_FIXED_POINT_EFFICIENCY_ONLY', cases=matrix,
        seeds=[2**117+2026100600+i for i in range(len(matrix))], timing_repeats=7,
        comparison_atol=1e-9, comparison_rtol=1e-12, scientific_inference_ready=False,
        timing_scope='Shared Mac while two formal workers run; paired alternating order; no total-fit or total-study speed claim',
        memory_scope='tracemalloc allocations during one evaluation, including returned arrays; inputs constructed beforehand; not process RSS',
        source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
            [args.baseline, Path(sys.modules[GroupedMmlLikelihood.__module__].__file__), Path(__file__)]},
        runtime=dict(python=platform.python_version(), numpy=np.__version__, platform=platform.platform()))
    (out/'plan.json').write_text(json.dumps(plan, indent=2)+'\n')
    records = []
    for i, case in enumerate(matrix):
        m, a, n, k, q = (case[key] for key in ['model','assignment','persons','categories','q'])
        steps = np.linspace(-.9,.9,k-1) if k>2 else np.zeros(1)
        shape = np.linspace(-1.,1.,k-1)**2; shape -= shape.mean()
        design = RecoveryDesign(model=m, assignment=a, persons=n, population_sd=1.2,
            rater=(-.45,0.,.45), task=(-.25,.25), criterion=(-.5,.1,.7),
            steps=tuple(tuple(v) for v in ([steps] if m=='RSM' else [steps+.15*(j-1)*shape for j in range(3)])),
            log_slopes=(-.4,0.,.4) if m=='GPCM' else ())
        sample = generate_population_sample(design, plan['seeds'][i]); data = sample['data']
        facets = ['Rater','Task','Criterion']
        config = dict(model=m, n_cat=k, rating_min=0, facet_names=facets,
            facet_levels={f:sample['declared_levels'][f] for f in facets},
            facet_signs={f:-1 for f in facets}, noncenter_facet='Criterion',
            step_facet=None if m=='RSM' else 'Criterion', slope_facet='Criterion' if m=='GPCM' else None)
        name = f'{m}_{a}_N{n}_K{k}_Q{q}'
        data.to_csv(out/(name+'.csv'), index=False)
        engines = [module.GroupedMmlLikelihood(data,config), GroupedMmlLikelihood(data,config)]
        point = np.r_[np.linspace(-.3,.4,engines[0].dimension-1), np.log(1.2)]
        nodes, weights = roots_hermitenorm(q); weights /= np.sqrt(2*np.pi)
        errors = {}
        # A nonzero point and the all-zero point include nonlinear slopes/SD.
        for label, value in [('nonzero',point), ('zero',np.zeros_like(point))]:
            for information in [False, True]:
                old, new = [e.evaluate(value,nodes,weights,information=information) for e in engines]
                assert set(old)==set(new)
                for key in old:
                    np.testing.assert_allclose(new[key],old[key],atol=plan['comparison_atol'],rtol=plan['comparison_rtol'])
                    errors[f'{label}_{information}_{key}'] = float(np.max(abs(np.asarray(new[key])-np.asarray(old[key]))))
        timings = [[],[]]
        for e in engines: e.evaluate(point,nodes,weights)
        for repeat in range(plan['timing_repeats']):
            for j in ([0,1] if repeat%2==0 else [1,0]):
                start=time.perf_counter(); engines[j].evaluate(point,nodes,weights)
                timings[j].append(time.perf_counter()-start)
        peaks=[]
        for e in engines:
            tracemalloc.start(); result=e.evaluate(point,nodes,weights)
            peaks.append(tracemalloc.get_traced_memory()[1]); tracemalloc.stop(); del result
        record = dict(case=case,id=name,input_sha256=hashlib.sha256((out/(name+'.csv')).read_bytes()).hexdigest(),
            design=asdict(design), config=config, point=point.tolist(), checks_pass=True, errors=errors,
            baseline_seconds=timings[0], candidate_seconds=timings[1],
            baseline_median=statistics.median(timings[0]), candidate_median=statistics.median(timings[1]),
            baseline_traced_peak_bytes=peaks[0], candidate_traced_peak_bytes=peaks[1])
        record['median_time_ratio']=record['baseline_median']/record['candidate_median']
        record['traced_peak_ratio']=peaks[0]/peaks[1]
        records.append(record); (out/(name+'.json')).write_text(json.dumps(record,indent=2)+'\n')
        print(name, 'checks pass; time ratio', round(record['median_time_ratio'],2), flush=True)
    summary=dict(classification=plan['classification'],cases=len(records),all_checks_pass=True,
        scientific_inference_ready=False,records=records)
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    return 0


if __name__=='__main__':
    raise SystemExit(main())
