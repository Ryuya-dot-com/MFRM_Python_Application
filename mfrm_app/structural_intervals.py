"""Bounded, pointwise fixed-SD MML intervals and replay without optimization.

Numerical availability is separate from finite-sample coverage qualification.
These objects do not change the eligibility of Person scores or other tests.
"""
import hashlib
import html
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

from . import engine, mml_fixed_sd_review as review, output_qualification


POLICY = dict(version='fixed_sd_structural_intervals_v1', scaled_gradient=.01,
              inverse_residual=1e-6, quadrature_score_shift=.01, covariance_change=.01)
METHOD = 'Pointwise asymptotic normal interval from full fixed-SD MML observed information'
CAUTION = ('Conditional on the specified fixed normal population and response model, independent Persons, '
           'ignorable assignment/missingness, an identified interior parameter and information growth. '
           'Numerical checks are sensitivity checks, not continuous-integration error bounds or finite-sample '
           'coverage guarantees. Not simultaneous intervals or equivalence tests. Person calibration '
           'uncertainty and estimated population SD are outside this output.')


def _source_hash():
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def _fingerprint(record):
    return engine.stable_json_fingerprint({k: v for k, v in record.items() if k != 'record_hash'}, 64)


def _assess(record):
    checks = []
    def check(name, passed, value=None, limit=None):
        checks.append(dict(check=name, passed=bool(passed), value=value, limit=limit))
    boundary = record['boundary']
    check('optimizer_termination', record['optimizer_success'] is True)
    check('structural_identification', boundary['structural_design_rank'] == len(record['coordinates']))
    check('structural_boundary', boundary['status'] == 'no_structural_ray_detected'
          and not boundary['no_finite_mle_certified'])
    rows = record['calculations']
    inverses, gradients = [], []
    for row in rows:
        h, g = np.asarray(row['information'], dtype=float), np.asarray(row['gradient'], dtype=float)
        p = len(record['coordinates'])
        valid = h.shape == (p, p) and g.shape == (p,) and np.isfinite(h).all() and np.isfinite(g).all()
        information = review._information_review(h, g, np.eye(p)) if valid else {}
        prefix = f"GH{row['points']}:"
        ok = information.get('positive_definite', False)
        check(prefix+'positive_information', ok)
        if not ok:
            inverses.append(None); gradients.append(g)
            continue
        v = np.asarray(information['covariance'])
        residual = h @ v - np.eye(p)
        residual = float(max(np.linalg.norm(residual, 1), np.linalg.norm(residual, np.inf)))
        scaled = information['curvature_scaled_gradient']
        check(prefix+'inverse_residual', np.isfinite(residual) and residual <= POLICY['inverse_residual'], residual, POLICY['inverse_residual'])
        check(prefix+'scaled_gradient', np.isfinite(scaled) and scaled <= POLICY['scaled_gradient'], scaled, POLICY['scaled_gradient'])
        inverses.append(v); gradients.append(g)
    if all(v is not None for v in inverses):
        v, fine = inverses
        delta = gradients[1]-gradients[0]
        shift = float(np.sqrt(max(0., delta @ v @ delta)))
        root = np.linalg.cholesky(rows[0]['information'])
        change = float(np.linalg.norm(root.T @ (fine-v) @ root, 2))
        check('quadrature_score_shift', shift <= POLICY['quadrature_score_shift'], shift, POLICY['quadrature_score_shift'])
        check('full_covariance_change', change <= POLICY['covariance_change'], change, POLICY['covariance_change'])
    return checks, inverses[0]


def _table(record):
    checks, covariance = _assess(record)
    passed = all(check['passed'] for check in checks)
    reason = '; '.join(check['check'] for check in checks if not check['passed'])
    matrix = np.asarray(record['target_jacobian'])
    point = np.asarray(record['coordinates'])
    estimates = matrix @ point
    target_covariance = matrix @ covariance @ matrix.T if covariance is not None else None
    z = norm.ppf((1+record['level'])/2)
    rows = []
    for i, target in enumerate(record['targets']):
        fixed = not np.any(matrix[i])
        variance = None if target_covariance is None else float(target_covariance[i, i])
        usable = passed and not fixed and variance is not None and np.isfinite(variance) and variance > 0
        se = float(np.sqrt(variance)) if usable else None
        rows.append(dict(**target, Estimate=float(estimates[i]), SE=se,
            Lower=float(estimates[i]-z*se) if usable else None,
            Upper=float(estimates[i]+z*se) if usable else None, CI_Level=record['level'],
            CIEligible=bool(usable), CoverageQualified=False,
            Status='model_approximation' if usable else ('fixed' if fixed else 'withheld'),
            Reason='' if usable else ('constraint_fixed' if fixed else reason or 'target_variance'),
            Method=METHOD, Caution=CAUTION))
    return rows, checks, target_covariance.tolist() if passed and target_covariance is not None else None


def compute_structural_intervals(result, *, facet, difference=None, level=.95):
    """Intervals for one facet's levels and an explicitly selected A-minus-B.

    Uses the returned point unchanged. Insufficient stationarity refuses
    intervals; there is no hidden refit, model change or fallback covariance.
    """
    if isinstance(level, (bool, np.bool_)) or not np.isfinite(level) or not 0 < level < 1:
        raise ValueError('Interval level must be between zero and one')
    prepared = review.prepare_fixed_sd_check(result)
    cfg, sizes = prepared['config'], prepared['sizes']
    if facet not in cfg['facet_names']:
        raise ValueError('Select a fitted structural facet')
    labels, jacobian = review._targets(cfg, sizes)
    positions = [i for i, row in enumerate(labels) if row['kind'] == 'facet' and row['facet'] == facet]
    levels = [labels[i]['level'] for i in positions]
    matrix = jacobian[positions].copy()
    targets = [dict(Target=level_name, Kind='level', Facet=facet, Level=level_name, Reference='') for level_name in levels]
    if difference is not None:
        if (not isinstance(difference, (list, tuple)) or len(difference) != 2
                or any(x not in levels for x in difference) or difference[0] == difference[1]):
            raise ValueError('A difference requires two distinct observed levels of the selected facet')
        a, b = difference
        matrix = np.vstack([matrix, matrix[levels.index(a)]-matrix[levels.index(b)]])
        targets.append(dict(Target=f'{a} − {b}', Kind='difference', Facet=facet, Level=a, Reference=b))
    q = cfg['quad_points']
    evaluated = review.evaluate_fixed_sd_quadrature_check(prepared, comparison_points=2*q-1)
    rows = evaluated['comparison']['rows']
    record = dict(schema=POLICY['version'], policy=dict(POLICY), source_hash=_source_hash(),
        source_input_id=prepared['input_id'], implementation=prepared['implementation'],
        optimizer_success=output_qualification.optimizer_termination_established(result),
        coordinates=prepared['coordinates'].tolist(), targets=targets, target_jacobian=matrix.tolist(),
        level=float(level), boundary=evaluated['structural_boundary'],
        calculations=[{k: row[k] for k in ('points', 'raw_nll', 'gradient', 'information')} for row in rows],
        scope=dict(model=cfg['model'], population_sd=engine.get_population_prior_sd(cfg),
                   noncenter_facet=cfg['noncenter_facet'], ratings=len(prepared['prep']['data']),
                   persons=cfg['n_person'], coverage_qualified=False))
    record['table'], record['checks'], record['target_covariance'] = _table(record)
    record['record_hash'] = _fingerprint(record)
    return record


def read_structural_intervals(result, record):
    """Validate source/content and derive the table again, without a refit."""
    if not isinstance(record, dict) or record.get('schema') != POLICY['version']:
        raise ValueError('Unsupported structural interval record')
    if record.get('record_hash') != _fingerprint(record):
        raise ValueError('Structural interval record changed')
    if record.get('policy') != POLICY or record.get('source_hash') != _source_hash():
        raise ValueError('Structural interval policy or numerical source changed')
    prepared = review.prepare_fixed_sd_check(result)
    if record['source_input_id'] != prepared['input_id']:
        raise ValueError('Structural intervals belong to different inputs, estimates or runtime')
    level = record.get('level')
    if isinstance(level, bool) or not isinstance(level, (int, float)) or not np.isfinite(level) or not 0 < level < 1:
        raise ValueError('Invalid stored interval level')
    q = prepared['config']['quad_points']
    if ([row['points'] for row in record['calculations']] != [q, 2*q-1]
            or not all(np.isfinite(row['raw_nll']) for row in record['calculations'])
            or not np.array_equal(record['coordinates'], prepared['coordinates'])
            or record['optimizer_success'] is not output_qualification.optimizer_termination_established(result)
            or abs(record['calculations'][0]['raw_nll']-prepared['reported_nll']) > 1e-9):
        raise ValueError('Stored numerical review disagrees with the fitted point or rule')
    rows, checks, covariance = _table(record)
    if engine.stable_json_fingerprint([rows, checks, covariance], 64) != engine.stable_json_fingerprint(
            [record['table'], record['checks'], record['target_covariance']], 64):
        raise ValueError('Stored interval table or qualification disagrees with its numerical evidence')
    return pd.DataFrame(rows)


def interval_json(record):
    return json.dumps(record, ensure_ascii=False, allow_nan=False, indent=2)


def interval_figure(table):
    """A common figure for screen and HTML; never draw a withheld interval."""
    import plotly.graph_objects as go
    available = table[table.CIEligible.eq(True)]
    return go.Figure(go.Scatter(x=available.Estimate, y=available.Target, mode='markers',
        error_x=dict(type='data', symmetric=False, array=available.Upper-available.Estimate,
                     arrayminus=available.Estimate-available.Lower))).update_layout(
        xaxis_title='Logits', yaxis_title=str(table.Facet.iloc[0]), template='plotly_white')


def interval_html(result, record):
    table = read_structural_intervals(result, record)
    columns = ['Target', 'Estimate', 'SE', 'Lower', 'Upper', 'CIEligible', 'CoverageQualified', 'Reason']
    return ('<!doctype html><meta charset="utf-8"><title>Structural intervals</title>'
            '<h1>Specified structural intervals</h1><p>'+str(100*record['level'])+'%: '+html.escape(METHOD)+'</p><p>'+html.escape(CAUTION)+'</p>'
            +table[columns].to_html(index=False, na_rep='Unavailable', escape=True)
            +interval_figure(table).to_html(full_html=False, include_plotlyjs=True)
            +'<p>Source: '+html.escape(record['source_input_id'])+'</p>')


def publication_table(result, record):
    table = read_structural_intervals(result, record)
    out = table[['Target', 'Estimate', 'SE']].copy()
    out['Interval'] = [f'[{row.Lower:.3f}, {row.Upper:.3f}]' if row.CIEligible else ''
                       for row in table.itertuples()]
    out['Use'] = table.Status.map({'model_approximation':'Approximate', 'fixed':'Fixed', 'withheld':'Withheld'})
    out[['Estimate', 'SE']] = out[['Estimate', 'SE']].round(3)
    caption = f"{100*record['level']:g}%: {METHOD}. {CAUTION}"
    if table.Reason.ne('').any():
        caption += ' Unavailable: '+', '.join(sorted(set(table.loc[table.Reason.ne(''), 'Reason'])))+'.'
    return out, caption
