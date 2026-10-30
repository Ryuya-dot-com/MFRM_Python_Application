"""A separate, jointly refitted fixed-population RSM/PCM interaction model.

The shared marginal likelihood supplies all nuisance cross-information. The
baseline fit and its fixed-calibration bias screens are never overwritten.
"""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import logsumexp
from scipy.stats import norm

from . import engine, mml_fixed_sd_review as review, mml_joint, structural_intervals
from .bias_analysis import adjust_p_values


VERSION = 'joint_interaction_fixed_normal_v1'
METHOD = 'Pointwise Wald approximation from full joint fixed-N(0,1) marginal information'
CAUTION = ('A separately refitted main-effects, steps and interaction model. Positive interaction '
           'raises expected coded scores. Independent Persons, the normal population and conditional '
           'response model, ignorable assignment/missingness, identification, an interior parameter '
           'and information growth are assumed. Intervals are pointwise; Holm adjusts the recorded '
           'family of requested cell/contrast tests, not the intervals. Numerical eligibility is not '
           'a finite-sample coverage/error-rate guarantee or a fairness verdict.')


def settings(values, facets, method):
    if values is None:
        return None
    if not isinstance(values, dict) or set(values) - {'facet_a', 'facet_b', 'difference', 'level'}:
        raise ValueError('Unknown joint interaction settings')
    a, b = values.get('facet_a'), values.get('facet_b')
    level, difference = values.get('level', .95), values.get('difference')
    if method != 'MML' or a not in facets or b not in facets or a == b:
        raise ValueError('Joint interactions require fixed-SD MML and two distinct structural facets')
    if isinstance(level, (bool, np.bool_)) or not np.isfinite(level) or not 0 < level < 1:
        raise ValueError('Interval level must be between zero and one')
    if difference is not None:
        if (not isinstance(difference, (list, tuple)) or len(difference) != 2
                or any(not isinstance(cell, (list, tuple)) or len(cell) != 2
                       or any(not isinstance(x, str) for x in cell) for cell in difference)
                or list(difference[0]) == list(difference[1])):
            raise ValueError('A joint contrast requires two distinct cells, each with an A and B level')
        difference = [list(cell) for cell in difference]
    return dict(facet_a=a, facet_b=b, difference=difference, level=float(level))


def _design(prepared, request):
    cfg, data = prepared['config'], prepared['prep']['data']
    if engine.get_population_prior_sd(cfg) != 1.:
        raise ValueError('Joint interaction inference currently requires a fixed N(0,1) population')
    a, b = request['facet_a'], request['facet_b']
    levels_a, levels_b = [list(map(str, cfg['facet_levels'][f])) for f in (a, b)]
    na, nb = len(levels_a), len(levels_b)
    if min(na, nb) < 2:
        raise ValueError('Each interaction facet requires at least two observed levels')
    A = np.vstack([np.eye(na-1), -np.ones(na-1)])
    B = np.vstack([np.eye(nb-1), -np.ones(nb-1)])
    # Column-major cells: A varies first, as in mfrmr's zero-margin coding.
    expansion = np.kron(B, A)
    map_a = {lvl: i for i, lvl in enumerate(levels_a)}
    map_b = {lvl: i for i, lvl in enumerate(levels_b)}
    cell_index = np.fromiter(
        (map_a[str(x)] + na * map_b[str(y)] for x, y in zip(data[a], data[b], strict=True)),
        dtype=int,
        count=len(data),
    )
    rows = expansion[cell_index]
    base_d = len(prepared['coordinates'])
    dimension = base_d + rows.shape[1]
    if dimension > 64:
        raise ValueError('Joint information calculation exceeds the 64-coordinate synchronous scope')
    if max(121, 2*cfg['quad_points']-1)*len(data[cfg['facet_names']].drop_duplicates())*cfg['n_cat']*(dimension+1) > 8_000_000:
        raise ValueError('Joint information calculation exceeds the synchronous memory scope')
    # Reuse the existing grouped kernel without changing its default path or
    # invalidating unrelated saved structural-interval source fingerprints.
    problem = prepared['problem']
    native = mml_joint._NativeRows(data, cfg, allow_all_centered=True)
    _, first, group = np.unique(native.design.reshape(len(data), -1), axis=0,
                               return_index=True, return_inverse=True)
    if not np.array_equal(rows, rows[first][group]):
        raise ValueError('Baseline grouping combines different interaction cells')
    problem.design = np.concatenate([problem.design, rows[first, None, :]*problem.k[None, :, None]], axis=2)
    problem.linear_d = dimension
    problem.dimension = dimension+1
    counts = np.bincount(cell_index, minlength=na*nb)
    cells = [[x, y] for y in levels_b for x in levels_a]
    labels = [dict(Target=json.dumps(cell, ensure_ascii=False), Kind='cell', Facet=a+' × '+b,
                   LevelA=cell[0], LevelB=cell[1], Reference='', Ratings=int(counts[i]))
              for i, cell in enumerate(cells)]
    matrix = np.c_[np.zeros((len(cells), base_d)), expansion]
    if request['difference'] is not None:
        if any(cell not in cells for cell in request['difference']):
            raise ValueError('Select observed levels of both interaction facets for the cell contrast')
        i, j = [cells.index(cell) for cell in request['difference']]
        matrix = np.vstack([matrix, matrix[i]-matrix[j]])
        labels.append(dict(Target=labels[i]['Target']+' − '+labels[j]['Target'], Kind='difference',
            Facet=a+' × '+b, LevelA='', LevelB='', Reference=labels[j]['Target'], Ratings=None))
    return dict(problem=problem, row_design=rows, base_dimension=base_d,
                targets=labels, target_jacobian=matrix, cells=cells, expansion=expansion)


def _row_nll(point, prepared, design, q):
    """Direct per-person/category evaluation, separate from grouped counts."""
    cfg, idx = prepared['config'], prepared['idx']
    params = engine.expand_params(point[:design['base_dimension']], prepared['sizes'], cfg)
    eta = engine.compute_base_eta(idx, params, cfg) + design['row_design'] @ point[design['base_dimension']:]
    steps = params['steps'] if cfg['model'] == 'RSM' else params['steps_mat'][idx['step_idx']]
    cumulative = np.r_[0., np.cumsum(steps)] if cfg['model'] == 'RSM' else np.c_[np.zeros(len(eta)), np.cumsum(steps, axis=1)]
    nodes, weights = mml_joint.normal_rule(q)
    active = weights > 0
    logits = (eta[:, None, None]+nodes[None, active, None])*np.arange(cfg['n_cat'])
    logits -= cumulative[:, None, :] if cfg['model'] == 'PCM' else cumulative
    log_p = logits-logsumexp(logits, axis=2, keepdims=True)
    observed = np.take_along_axis(log_p, np.broadcast_to(idx['score_k'][:, None, None], (*log_p.shape[:2], 1)), axis=2)[..., 0]
    person_red = engine._get_person_reduction(idx, len(idx['rows_by_person']))
    person_ll = engine._sum_2d_by_group(observed, person_red)[person_red['nonempty']]
    if person_ll.size == 0:
        return 0.0
    return -float(np.sum(logsumexp(person_ll + np.log(weights[active])[None, :], axis=1)))


def _implementation(prepared):
    return dict(**prepared['implementation'], joint_source=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        interval_source=structural_intervals._source_hash(),
        probability_source=hashlib.sha256(Path(adjust_p_values.__code__.co_filename).read_bytes()).hexdigest())


def _fingerprint(record):
    return engine.stable_json_fingerprint({k: v for k, v in record.items() if k != 'record_hash'}, 64)


def _tables(record):
    rows, checks, covariance = structural_intervals._table(record)
    for row in rows:
        row.update(Method=METHOD, Caution=CAUTION,
            WaldP=float(2*norm.sf(abs(row['Estimate']/row['SE']))) if row['CIEligible'] else None,
            PMethod='Two-sided asymptotic normal Wald; joint marginal covariance',
            FamilySize=len(rows), Family='All cells and the prespecified contrast in this request')
    adjusted = adjust_p_values([row['WaldP'] if row['WaldP'] is not None else np.nan for row in rows], 'holm')
    for row, p in zip(rows, adjusted, strict=True):
        row['HolmP'] = float(p) if np.isfinite(p) else None
    return rows, checks, covariance


def fit_joint_interaction(result, *, facet_a, facet_b, difference=None, level=.95):
    """Explicitly refit one interaction; retain the baseline result unchanged."""
    request = settings(dict(facet_a=facet_a, facet_b=facet_b, difference=difference, level=level),
                       result['config']['facet_names'], result['config']['method'])
    prepared = review.prepare_fixed_sd_check(result)
    design = _design(prepared, request)
    problem = design['problem']
    boundary = problem.structural_recession_audit()
    if boundary['structural_design_rank'] != problem.linear_d:
        raise ValueError('Joint interaction design is not identified; no intervals or p-values')
    if boundary['status'] != 'no_structural_ray_detected':
        raise ValueError('Joint interaction boundary review: '+boundary['status'])
    q = prepared['config']['quad_points']
    nodes, weights = mml_joint.normal_rule(q)
    start = np.r_[prepared['coordinates'], np.zeros(problem.linear_d-design['base_dimension'])]
    def objective(point):
        value = problem.evaluate(np.r_[point, 0.], nodes, weights, moments=False)
        return value['nll'], value['gradient'][:-1]
    opt = minimize(objective, start, jac=True, method='L-BFGS-B',
                   options=dict(maxiter=int(prepared['config']['maxit']), gtol=1e-8, ftol=1e-14, maxls=50))
    if not np.isfinite(opt.x).all() or not np.isfinite(opt.fun):
        raise ValueError('Joint optimizer did not return a finite candidate')
    calculations = []
    for order in (q, 2*q-1):
        nodes, weights = mml_joint.normal_rule(order)
        value = problem.evaluate(np.r_[opt.x, 0.], nodes, weights, information=True, moments=False)
        row_nll = _row_nll(opt.x, prepared, design, order)
        if not np.isfinite(row_nll) or abs(value['nll']-row_nll) > 1e-8:
            raise ValueError('Joint row and grouped marginal likelihoods disagree')
        calculations.append(dict(points=order, raw_nll=value['nll'], row_nll=row_nll,
            gradient=value['gradient'][:-1].tolist(), information=value['information'][:-1, :-1].tolist()))
    record = dict(schema=VERSION, settings=request, policy=dict(structural_intervals.POLICY),
        implementation=_implementation(prepared), source_input_id=prepared['input_id'],
        optimizer_success=bool(opt.success), optimizer_message=str(opt.message), iterations=int(opt.nit),
        coordinates=opt.x.tolist(), base_dimension=design['base_dimension'],
        targets=design['targets'], target_jacobian=design['target_jacobian'].tolist(), level=request['level'],
        boundary=boundary, calculations=calculations, baseline_nll=prepared['reported_nll'],
        model=prepared['config']['model'], estimation='Separate joint refit; fixed N(0,1); L-BFGS-B',
        scope='One prespecified two-facet interaction; unit weights; no anchors, penalties or regression')
    labels, jacobian = review._targets(prepared['config'], prepared['sizes'])
    record['main_effects_and_steps'] = [dict(**label, Estimate=float(value)) for label, value in
        zip(labels, jacobian @ opt.x[:design['base_dimension']], strict=True)]
    record['table'], record['checks'], record['target_covariance'] = _tables(record)
    record['record_hash'] = _fingerprint(record)
    return record


def read_joint_interaction(result, record):
    """Validate and re-derive display without optimization or likelihood calls."""
    if not isinstance(record, dict) or record.get('schema') != VERSION or record.get('record_hash') != _fingerprint(record):
        raise ValueError('Joint interaction record is unsupported or changed')
    prepared = review.prepare_fixed_sd_check(result)
    if record['source_input_id'] != prepared['input_id'] or record['implementation'] != _implementation(prepared):
        raise ValueError('Joint interaction inputs, numerical source or runtime changed')
    request = settings(record['settings'], prepared['config']['facet_names'], prepared['config']['method'])
    design = _design(prepared, request)
    q = prepared['config']['quad_points']
    if (record['policy'] != structural_intervals.POLICY or record['level'] != request['level']
            or record['targets'] != design['targets'] or record['target_jacobian'] != design['target_jacobian'].tolist()
            or record['base_dimension'] != design['base_dimension'] or record['model'] != prepared['config']['model']
            or record['baseline_nll'] != prepared['reported_nll']
            or np.asarray(record['coordinates']).shape != (design['problem'].linear_d,)
            or not np.isfinite(record['coordinates']).all()
            or [row['points'] for row in record['calculations']] != [q, 2*q-1]):
        raise ValueError('Joint interaction targets or numerical evidence disagree with their source')
    rows, checks, covariance = _tables(record)
    if engine.stable_json_fingerprint([rows, checks, covariance], 64) != engine.stable_json_fingerprint(
            [record['table'], record['checks'], record['target_covariance']], 64):
        raise ValueError('Joint interaction table disagrees with its numerical evidence')
    labels, jacobian = review._targets(prepared['config'], prepared['sizes'])
    expected = [dict(**label, Estimate=float(value)) for label, value in zip(labels,
        jacobian @ np.asarray(record['coordinates'])[:design['base_dimension']], strict=True)]
    if record['main_effects_and_steps'] != expected:
        raise ValueError('Joint main effects and steps disagree with the retained coordinates')
    return pd.DataFrame(rows)


def publication_table(result, record):
    table = read_joint_interaction(result, record)
    columns = ['Target', 'Estimate', 'SE', 'Lower', 'Upper', 'CI_Level', 'WaldP', 'HolmP',
               'CIEligible', 'CoverageQualified', 'Status', 'Reason', 'FamilySize']
    return table[columns], METHOD+'. '+CAUTION
