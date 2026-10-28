"""One retained input: independent IC arithmetic, likelihood and numerical sensitivity.

This is a local development check, not a model-selection performance study.
The independent evaluator below does not call the app's probability or expansion code.
"""
from pathlib import Path
import hashlib
import json
import sys
import time

import numpy as np
import pandas as pd
from scipy.integrate import quad
from scipy.optimize import minimize
from scipy.special import logsumexp, roots_hermitenorm
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
SOURCE = ROOT / "validation/generated/mml_fixed_sd_api_20260915/checked_v1/complete_input"
OUT = ROOT / "validation/generated/information_criteria_validation_20261027"
MODELS = {"RSM": 6, "PCM": 8, "GPCM": 9}
PROTOCOL = {
    "scope": "One retained complete input, 32 persons x 8 responses; fixed SD 1.1, unit weights, no anchors or penalties",
    "question": "Are K, person-based BIC, likelihood values and model differences consistent with independent calculations?",
    "default_fit": {"mml_engine": "EM", "quad_points": 15, "maxit": 400, "reltol": 1e-6},
    "reference_fit": {"nodes": 81, "start": "zero", "method": "L-BFGS-B", "maxiter": 500, "ftol": 1e-13, "gtol": 1e-7},
    "checks": {"same_grid_nll": 1e-9, "ic_arithmetic": 1e-9, "reference_integral_refinement": 1e-7,
               "reference_gradient_supnorm": 1e-4, "native_printed_ic": 5e-4},
    "continuous_integral": {"bounds_standard_normal": [10, 12], "epsabs": 1e-13, "epsrel": 1e-11},
    "failure_policy": "Retain all solver outcomes; do not change limits, replace failures or promote inference qualifications.",
    "inference_qualified": False,
}


def write(name, value):
    (OUT / name).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class Reference:
    """Explicit 2-rater / 2-task / 2-criterion, four-category model."""
    def __init__(self, data, model):
        self.model = model
        self.person, levels = pd.factorize(data.Person, sort=True)
        self.nperson = len(levels)
        self.rater = data.Rater.map({"R1": 0, "R2": 1}).to_numpy()
        self.task = data.Task.map({"T1": 0, "T2": 1}).to_numpy()
        self.criterion = data.Criterion.map({"C1": 0, "C2": 1}).to_numpy()
        self.score = data.Score.to_numpy(dtype=int)
        assert self.nperson == 32 and len(data) == 256 and data.Weight.eq(1).all()
        assert np.bincount(self.person).tolist() == [8] * 32
        self.rules = {}

    def log_response(self, point, nodes):
        point = np.asarray(point)
        assert len(point) == MODELS[self.model]
        rater = np.array([point[0], -point[0]])
        task = np.array([point[1], -point[1]])
        location = -rater[self.rater] - task[self.task] - point[2:4][self.criterion]
        free_steps = point[4:6][None, :] if self.model == "RSM" else point[4:8].reshape(2, 2)
        steps = np.column_stack([free_steps, -free_steps.sum(axis=1)])
        cumulative = np.column_stack([np.zeros(len(steps)), np.cumsum(steps, axis=1)])
        row_steps = cumulative[np.zeros(len(location), dtype=int) if self.model == "RSM" else self.criterion]
        slopes = np.exp(np.array([point[-1], -point[-1]])[self.criterion]) if self.model == "GPCM" else np.ones(len(location))
        logits = ((np.asarray(nodes)[:, None, None] + location[None, :, None]) * np.arange(4)
                  - row_steps[None, :, :]) * slopes[None, :, None]
        return logits[:, np.arange(len(location)), self.score] - logsumexp(logits, axis=2)

    def nll(self, point, order=81):
        if order not in self.rules:
            nodes, weights = roots_hermitenorm(order)
            self.rules[order] = (nodes * 1.1, np.log(weights / np.sqrt(2*np.pi)))
        nodes, logweights = self.rules[order]
        rows = self.log_response(point, nodes)
        grouped = np.stack([rows[:, self.person == p].sum(axis=1) for p in range(self.nperson)], axis=1)
        return float(-logsumexp(grouped + logweights[:, None], axis=0).sum())

    def continuous(self, point, bound):
        values, errors = [], []
        # Scale each person's integral by a finite-grid mode to avoid tiny absolute integrals.
        grid = np.linspace(-8, 8, 81)
        row_lp = self.log_response(point, grid*1.1)
        shifts = [float(np.max(row_lp[:, self.person == p].sum(axis=1) + norm.logpdf(grid)))
                  for p in range(self.nperson)]
        for p, shift in enumerate(shifts):
            def integrand(z):
                lp = self.log_response(point, np.array([z*1.1]))[0, self.person == p].sum()
                return float(np.exp(lp + norm.logpdf(z) - shift))
            value, error = quad(integrand, -bound, bound, epsabs=1e-13, epsrel=1e-11, limit=150)
            values.append(np.log(value) + shift)
            errors.append(error/value)
        return {"nll": -float(np.sum(values)), "relative_error_sum_estimate": float(sum(errors))}


def main():
    import streamlit_app as app
    from mfrm_app import engine
    manifest = {str(p.relative_to(ROOT)): digest(p) for p in [SOURCE/'data.csv', SOURCE/'analysis_snapshot.json',
                ROOT/'streamlit_app.py', ROOT/'mfrm_app/engine.py', ROOT/'mfrm_app/output_qualification.py', Path(__file__)]}
    PROTOCOL['source_sha256'] = manifest
    if (OUT/'protocol.json').exists():
        assert json.loads((OUT/'protocol.json').read_text()) == PROTOCOL
    else:
        write('protocol.json', PROTOCOL)
    data = pd.read_csv(SOURCE/'data.csv')
    if not (OUT/'app_fits.json').exists():
        fits = {}
        original_refit = app._refit_candidate_model
        def capture(res, candidate, **kwargs):
            fit = original_refit(res, candidate, **kwargs)
            fits[candidate] = fit
            print('app refit', candidate, getattr(fit.get('opt'), 'success', None) if fit else 'failed', flush=True)
            return fit
        fits['RSM'] = app.mfrm_estimate(data, person_col='Person', facet_cols=['Rater', 'Task', 'Criterion'],
            score_col='Score', weight_col='Weight', rating_min=0, rating_max=3, noncenter_facet='Criterion',
            model='RSM', method='MML', population_prior_sd=1.1, **PROTOCOL['default_fit'])
        app._refit_candidate_model = capture
        try:
            bundle = app.compute_model_choice_comparison(fits['RSM'], step_facet='Criterion', slope_facet='Criterion')
        finally:
            app._refit_candidate_model = original_refit
        bundle['comparison'].to_csv(OUT/'app_comparison.csv', index=False)
        saved = {}
        for model, fit in fits.items():
            assert fit is not None, model
            s = fit['summary'].iloc[0]
            saved[model] = dict(point=fit['opt'].x.tolist(), success=bool(fit['opt'].success),
                message=str(fit['opt'].message), iterations=int(fit['opt'].nit),
                gradient_norm=float(fit['opt'].gradient_norm),
                summary={k: float(s[k]) for k in ['LogLik','AIC','BIC','KParams','Persons','N','BICSampleSize']},
                sizes=dict(engine.build_param_sizes(fit['config'])))
        write('app_fits.json', saved)
    saved = json.loads((OUT/'app_fits.json').read_text())
    checks = []
    for model, k in MODELS.items():
        if (OUT/f'{model}.json').exists():
            checks.append(json.loads((OUT/f'{model}.json').read_text())); continue
        start = time.monotonic()
        r = Reference(data, model); fit = saved[model]; s = fit['summary']; point = np.array(fit['point'])
        assert s['KParams'] == k and s['Persons'] == s['BICSampleSize'] == 32 and s['N'] == 256
        assert abs(s['AIC'] - (-2*s['LogLik']+2*k)) < 1e-9
        assert abs(s['BIC'] - (-2*s['LogLik']+k*np.log(32))) < 1e-9
        replay = r.nll(point, 15)
        assert abs(replay+s['LogLik']) < 1e-9, (model,replay,s['LogLik'])
        bounds = [(None,None)]*k
        if model == 'GPCM': bounds[-1] = (-3,3)
        options = {key:PROTOCOL['reference_fit'][key] for key in ['maxiter','ftol','gtol']}
        refit = minimize(r.nll, np.zeros(k), jac='3-point', method='L-BFGS-B', bounds=bounds, options=options)
        default_integrals = [r.continuous(point, b) for b in [10,12]]
        reference_integrals = [r.continuous(refit.x, b) for b in [10,12]]
        ref_nll = reference_integrals[-1]['nll']
        record = dict(model=model, K=k, same_grid_nll_difference=replay+s['LogLik'],
            default=s, default_terminated=fit['success'], default_integrals=default_integrals,
            reference=dict(point=refit.x.tolist(), success=bool(refit.success), message=str(refit.message),
                iterations=int(refit.nit), gradient_supnorm=float(np.max(np.abs(refit.jac))),
                nll_81=float(refit.fun), nll_161=r.nll(refit.x,161), integrals=reference_integrals,
                AIC=2*ref_nll+2*k, BIC=2*ref_nll+k*np.log(32)),
            integration_ic_shift=2*(default_integrals[-1]['nll']+s['LogLik']),
            optimization_ic_shift=2*(ref_nll-default_integrals[-1]['nll']),
            elapsed_seconds=time.monotonic()-start)
        record['reference_checks_pass'] = bool(refit.success and np.max(abs(refit.jac))<1e-4
            and abs(reference_integrals[0]['nll']-ref_nll)<1e-7
            and abs(r.nll(refit.x,161)-ref_nll)<1e-7)
        write(f'{model}.json',record); checks.append(record)
        print(model, 'K',k,'default IC',s['AIC'],s['BIC'],'reference IC',record['reference']['AIC'],record['reference']['BIC'],
              'reference pass',record['reference_checks_pass'],flush=True)
    rows=[]
    for lane in ['default','reference']:
        for check in checks:
            row=dict(lane=lane, Model=check['model'], K=check['K'], **{ic:check[lane][ic] for ic in ['AIC','BIC']})
            rows.append(row)
    table=pd.DataFrame(rows)
    for metric in ['AIC','BIC']:
        table['Delta'+metric]=table[metric]-table.groupby('lane')[metric].transform('min')
    table.to_csv(OUT/'independent_comparison.csv',index=False)
    assert all(digest(ROOT/name)==h for name,h in manifest.items()), 'Source changed during validation'
    write('numerical_summary.json',dict(checks=checks, all_reference_checks_pass=all(c['reference_checks_pass'] for c in checks),
        inference_qualified=False, scope=PROTOCOL['scope']))


if __name__ == '__main__':
    main()
