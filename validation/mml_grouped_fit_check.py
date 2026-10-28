"""Cold-start development refits of nine retained data sets at Q321/Q641."""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import sys
import zipfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
import streamlit_app as app
from validation.mml_free_sd_stationarity_adapter import prepare_app_free_sd_problem_from_data
from validation.mml_gpcm_joint_check import prepare_problem
from validation.mml_grouped_fit import fit_grouped_joint, normal_rule, OPTIONS, LIMITS
from validation.mml_joint_information_check import save
from validation.mml_rasch_quadrature_review import verify


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    base = ROOT/'validation/generated'
    parents = {name: base/f'{name}_20261004' for name in ['g4_joint_information', 'g4_gpcm_joint_information',
        'g4_rasch_quadrature', 'g4_gpcm_joint_refinement']}
    for path in parents.values(): verify(path)
    out = args.output; out.mkdir(parents=True, exist_ok=False)
    limits = dict(native_nll=1e-8, native_gradient=1e-8, prior_coordinate=2e-5,
        prior_nll_per_person=1e-6, prior_covariance_relative=.005,
        q_coordinate=1e-5, q_nll_per_person=1e-6, q_eap=1e-4, q_sd=1e-4)
    save(out/'protocol.json', dict(classification='DEVELOPMENT_ONLY', scientific_inference_ready=False,
        planned=18, data_status='Previously observed; no qualification data generated',
        models=['RSM', 'PCM', 'GPCM'], designs=['regular', 'sparse', 'unbalanced'], q=[321, 641],
        initialization='Native input-only cold start separately for each order, followed by one mandatory polish',
        estimator_options=OPTIONS, estimator_limits=LIMITS, comparison_limits=limits,
        parent_manifests={k: hashlib.sha256((p/'manifest.json').read_bytes()).hexdigest() for k, p in parents.items()},
        limitation='Finite-rule checks and retained-point proximity do not qualify recovery, coverage or all future integral errors'))
    sources = [ROOT/'streamlit_app.py', *sorted((ROOT/'mfrm_app').glob('*.py')),
        *[ROOT/'validation'/n for n in ['mml_grouped_fit_check.py', 'mml_grouped_fit.py', 'mml_grouped_likelihood.py',
            'mml_free_sd_stationarity_adapter.py', 'mml_gpcm_joint_check.py', 'mml_gpcm_integral.py',
            'mml_rasch_integral.py', 'mml_joint_information_check.py', 'mml_rasch_quadrature_review.py']]]
    hashes = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    save(out/'source_manifest.json', hashes)
    with zipfile.ZipFile(out/'sources.zip', 'w', zipfile.ZIP_DEFLATED) as z:
        for p in sources: z.write(p, str(p.relative_to(ROOT)))
    save(out/'runtime.json', dict(python=sys.version, packages={d.metadata['Name']: d.version for d in importlib.metadata.distributions()}))
    references, records, failures = [], [], []
    for model in ['RSM', 'PCM', 'GPCM']:
        source = parents['g4_gpcm_joint_information' if model == 'GPCM' else 'g4_joint_information']
        cases = {c['id']: c for c in json.loads((source/'reference_input.json').read_text())['cases']}
        for design in ['regular', 'sparse', 'unbalanced']:
            identifier = model+'_'+design
            case = cases[identifier+'_Q31']; data = pd.DataFrame(case['data']); cfg = case['config']
            if model == 'GPCM':
                _, problem = prepare_problem(data, 321)
                prior = json.loads((parents['g4_gpcm_joint_refinement']/(identifier+'_Q321_from61.json')).read_text())
                old_nll, old_covariance = prior['nll'], np.array(prior['covariance'])
            else:
                _, problem = prepare_app_free_sd_problem_from_data(data, person_col='Person', facet_cols=cfg['facet_names'],
                    score_col='Score', rating_min=cfg['rating_min'], rating_max=cfg['rating_max'], model=model,
                    noncenter_facet=cfg['noncenter_facet'], step_facet=cfg['step_facet'], quadrature_points=321)
                prior = json.loads((parents['g4_rasch_quadrature']/(identifier+'_Q241_from61.json')).read_text())
                old_nll, old_covariance = prior['finite']['nll'], np.array(prior['covariance']['finite'])
            pair = []
            for q in [321, 641]:
                run_id = identifier+'_Q'+str(q)
                try:
                    fit = fit_grouped_joint(data, problem, q=q)
                    save(out/(run_id+'_fit.json'), fit)
                    if fit['point'] is None:
                        raise ValueError('Retained regular fixture had prefit failure: '+fit['status'])
                    point = np.asarray(fit['point']); nodes, weights = normal_rule(q)
                    positive = weights > 0
                    # Zero tail weights contribute exactly zero; native log(w)
                    # needs only the remaining nodes, without renormalization.
                    quad = dict(nodes=np.exp(point[-1])*nodes[positive], weights=weights[positive], sd=np.exp(point[-1]))
                    native_nll, native_g = app.mfrm_loglik_mml_value_grad(point[:-1], problem.idx, problem.config,
                        problem.sizes, quad, include_log_sigma=True)
                    evaluated = fit['evaluated']
                    errors = dict(native_nll=abs(native_nll-evaluated['nll']),
                        native_gradient=float(np.max(abs(native_g-evaluated['gradient']))),
                        prior_coordinate=float(np.max(abs(point-np.array(prior['point'])))),
                        prior_nll_per_person=abs(evaluated['nll']-old_nll)/cfg.get('n_person', data.Person.nunique()),
                        prior_covariance_relative=None)
                    if fit['covariance'] is not None:
                        scale = np.sqrt(np.outer(np.diag(old_covariance), np.diag(old_covariance)))
                        errors['prior_covariance_relative'] = float(np.max(abs(fit['covariance']-old_covariance)/scale))
                    checks = {k: v is not None and v < limits[k] for k, v in errors.items()}
                    checks['finite_rule_numerical'] = fit['numerical_pass']
                    record = dict(id=run_id, errors=errors, checks=checks, all_checks_pass=all(checks.values()))
                    save(out/(run_id+'_comparison.json'), record); records.append(record); pair.append(fit)
                    references.append(dict(id=run_id, data=case['data'], config=cfg, point=point, nodes=nodes, weights=weights))
                    print(run_id, fit['status'], 'failed', [k for k,v in checks.items() if not v], flush=True)
                except Exception as exc:
                    failures.append(dict(id=run_id, type=type(exc).__name__, message=str(exc)))
                    save(out/'failures.json', failures)
                    print(run_id, 'FAILED', repr(exc), flush=True)
            if len(pair) == 2:
                low, high = pair
                errors = dict(q_coordinate=float(np.max(abs(low['point']-high['point']))),
                    q_nll_per_person=abs(low['evaluated']['nll']-high['evaluated']['nll'])/data.Person.nunique(),
                    **{'q_'+k: float(np.max(abs(low['evaluated'][k]-high['evaluated'][k]))) for k in ['eap', 'sd']})
                checks = {k: v < limits[k] for k, v in errors.items()}
                save(out/(identifier+'_q_comparison.json'), dict(errors=errors, checks=checks, all_checks_pass=all(checks.values())))
    if any(hashlib.sha256((ROOT/n).read_bytes()).hexdigest() != h for n,h in hashes.items()):
        raise ValueError('Source changed during the retained study')
    save(out/'reference_input.json', dict(classification='DEVELOPMENT_ONLY', cases=references))
    q_records = {p.name: json.loads(p.read_text()) for p in out.glob('*_q_comparison.json')}
    summary = dict(classification='DEVELOPMENT_ONLY', scientific_inference_ready=False,
        planned=18, completed=len(records), failures=failures, comparisons=records, q_comparisons=q_records,
        all_checks_pass=len(records)==18 and len(q_records)==9 and not failures
            and all(r['all_checks_pass'] for r in records) and all(r['all_checks_pass'] for r in q_records.values()))
    save(out/'summary.json', summary)
    save(out/'manifest.json', {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir() if p.is_file()})
    return 0 if summary['all_checks_pass'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
