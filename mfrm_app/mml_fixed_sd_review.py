"""Same-point numerical review for fixed-SD RSM/PCM, without UI or refitting.

Observed marginal information includes all estimated structural coordinates.
The population SD is fixed; its row/column is therefore excluded BEFORE
inversion. No ridge, eigenvalue floor, or inference qualification is applied.
"""
from copy import deepcopy
import hashlib
import importlib.metadata
import operator
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import cho_solve

from . import engine, mml_joint


def prepare_fixed_sd_check(result):
    """Validate and identify a fitted result without evaluating its likelihood."""
    cfg, prep = deepcopy(result['config']), deepcopy(result['prep'])
    if cfg.get('method') != 'MML' or cfg.get('model') not in ('RSM', 'PCM'):
        raise ValueError('The same-point review requires MML RSM or PCM')
    if cfg.get('estimate_population_sd'):
        raise ValueError('Population SD must be fixed')
    q = cfg.get('quad_points')
    if isinstance(q, (bool, np.bool_)):
        raise ValueError('Quadrature order must be an integer')
    q = operator.index(q)
    if not 3 <= q <= 121:
        raise ValueError('The synchronous review supports fitted GH orders 3–121')
    sigma = engine.get_population_prior_sd(cfg)
    if not np.isfinite(sigma) or not .05 <= sigma <= 10:
        raise ValueError('The same-point review requires fixed SD in [0.05, 10]')
    if prep['facet_names'] != cfg['facet_names']:
        raise ValueError('Prepared facets disagree with the fitted model')
    data = prep['data']
    if not len(data):
        raise ValueError('Same-point GH comparison requires observed responses')
    for facet in ['Person', *cfg['facet_names']]:
        levels = list(prep['levels'][facet])
        if (list(data[facet].cat.categories) != levels or data[facet].isna().any()
                or set(data[facet].astype(str)) != set(levels)):
            raise ValueError('Prepared category order or observed levels disagree with the fitted level map')
        if facet != 'Person' and (levels != list(cfg['facet_levels'][facet])
                                  or levels != list(cfg['facet_specs'][facet]['levels'])):
            raise ValueError('Configuration and prepared facet levels disagree')
    if cfg['n_person'] != len(prep['levels']['Person']):
        raise ValueError('Person count disagrees with observed levels')
    raw = np.asarray(data['Score'], dtype=float) - cfg['rating_min']
    if not np.isin(raw, np.arange(cfg['n_cat'])).all() or not np.array_equal(raw, data['score_k']):
        raise ValueError('Prepared scores disagree with original responses')
    sizes = engine.build_param_sizes(cfg)
    point = np.asarray(result['opt'].x, dtype=float).copy()
    if point.shape != (sum(sizes.values()),) or not np.isfinite(point).all():
        raise ValueError('Structural coordinates are nonfinite or inconsistent')
    if not 0 < len(point) <= 64:
        raise ValueError('The synchronous review requires 1–64 structural coordinates')
    # Bound arrays BEFORE constructing the grouped counts/design. Distinct
    # observed facet combinations bound the number of grouped cells from above.
    cells = len(data[cfg['facet_names']].drop_duplicates())
    if (max(121, 2*q-1) * cells * cfg['n_cat'] * (len(point)+1) > 8_000_000
            or cfg['n_person'] * cells * cfg['n_cat'] > 2_000_000):
        raise ValueError('Information calculation exceeds the synchronous memory scope')
    # This constructor checks unit weights, simple centering, anchors, dummy
    # facets, regression and penalties. It evaluates no likelihood or optimizer.
    problem = mml_joint.GroupedMmlLikelihood(data, cfg, allow_all_centered=True)
    if problem.dimension != len(point)+1:
        raise ValueError('Grouped design disagrees with structural coordinate dimensions')
    params = engine.expand_params(point, sizes, cfg)
    for key in ('facets', 'steps', 'steps_mat'):
        actual, expected = result['params'].get(key), params[key]
        pairs = [(actual[f], expected[f]) for f in cfg['facet_names']] if key == 'facets' else [(actual, expected)]
        for left, right in pairs:
            if right is None:
                consistent = left is None
            else:
                consistent = np.array_equal(np.asarray(left), right)
            if not consistent:
                raise ValueError('Stored structural parameters disagree with returned coordinates')
    implementation = dict(
        source_sha256={Path(module.__file__).name: hashlib.sha256(Path(module.__file__).read_bytes()).hexdigest()
                       for module in (sys.modules[__name__], engine, mml_joint)},
        python=sys.version, packages={name: importlib.metadata.version(name) for name in ('numpy', 'scipy', 'pandas')})
    inputs = dict(config=cfg, prep=prep, params=result['params'], coordinates=point,
                  reported_nll=float(result['opt'].fun),
                  optimizer_success=getattr(result['opt'], 'success', None),
                  optimizer_message=str(getattr(result['opt'], 'message', '')))
    input_hash = engine.stable_json_fingerprint(inputs, 64)
    return dict(config=cfg, prep=prep, sizes=sizes, coordinates=point, problem=problem,
                idx=engine.build_indices(prep, step_facet=cfg.get('step_facet')),
                reported_nll=inputs['reported_nll'], implementation=implementation,
                input_sha256={'fitted_result': input_hash},
                input_id=engine.stable_json_fingerprint(dict(inputs=input_hash, implementation=implementation), 64))


def _targets(config, sizes):
    """Exact linear maps, including the dependent last facet level and step."""
    def expand(point):
        p = engine.expand_params(point, sizes, config)
        values, labels = [], []
        for facet in config['facet_names']:
            for level, value in zip(config['facet_levels'][facet], p['facets'][facet], strict=True):
                labels.append(dict(kind='facet', facet=facet, level=str(level)))
                values.append(value)
        steps = [p['steps']] if config['model'] == 'RSM' else p['steps_mat']
        levels = ['shared'] if config['model'] == 'RSM' else config['facet_levels'][config['step_facet']]
        for level, row in zip(levels, steps, strict=True):
            for index, value in enumerate(row):
                labels.append(dict(kind='step', facet=config.get('step_facet') or 'shared',
                                   level=str(level), step=index + 1))
                values.append(value)
        return np.array(values), labels
    origin, labels = expand(np.zeros(sum(sizes.values())))
    jacobian = np.column_stack([expand(basis)[0]-origin for basis in np.eye(sum(sizes.values()))])
    return labels, jacobian


def _information_review(information, gradient, targets):
    """Retain failure evidence; invert only positive-definite raw information."""
    h = np.asarray(information, dtype=float)
    g, j = np.asarray(gradient, dtype=float), np.asarray(targets, dtype=float)
    if h.shape != (g.size, g.size) or j.ndim != 2 or j.shape[1] != g.size:
        raise ValueError('Inconsistent information, gradient or target dimensions')
    out = dict(status='nonfinite_information', positive_definite=False, regularized=False,
               covariance=None, target_se=None, newton_correction=None,
               target_newton_correction=None, max_target_newton_over_se=None,
               curvature_scaled_gradient=None, inverse_residual=None)
    if not np.isfinite(h).all() or not np.isfinite(g).all() or not np.isfinite(j).all():
        return out
    out['symmetry_residual'] = float(np.max(np.abs(h-h.T)))
    h = (h+h.T)/2
    try:
        eigenvalues = np.linalg.eigvalsh(h)
        out.update(minimum_eigenvalue=float(eigenvalues[0]), maximum_eigenvalue=float(eigenvalues[-1]))
        if eigenvalues[0] <= 0:
            return dict(out, status='nonpositive_information')
        lower = np.linalg.cholesky(h)
    except np.linalg.LinAlgError:
        return dict(out, status='nonpositive_information')
    covariance = cho_solve((lower, True), np.eye(g.size))
    correction = cho_solve((lower, True), g)
    variance = np.einsum('ij,jk,ik->i', j, covariance, j)
    residual = float(np.max(np.abs(h @ covariance - np.eye(g.size))))
    if not np.isfinite(covariance).all() or not np.isfinite(correction).all() or np.any(variance < 0):
        return dict(out, status='unstable_inverse')
    se = np.sqrt(variance)
    target_correction = j @ correction
    estimated = np.any(j != 0, axis=1)
    if not np.isfinite(se).all() or np.any(se[estimated] <= 0) or not np.isfinite(residual):
        return dict(out, status='unstable_inverse')
    return dict(out, status='positive_definite', positive_definite=True,
                condition_number=float(eigenvalues[-1]/eigenvalues[0]),
                covariance=covariance.tolist(), target_se=se.tolist(),
                newton_correction=correction.tolist(), target_newton_correction=target_correction.tolist(),
                max_target_newton_over_se=float(np.max(np.abs(target_correction[estimated])/se[estimated])),
                curvature_scaled_gradient=float(np.linalg.norm(np.linalg.solve(lower, g))),
                inverse_residual=residual)


def evaluate_fixed_sd_quadrature_check(prepared, *, comparison_points=None):
    """Evaluate all rules at the same point; no pass threshold or refitting."""
    cfg, idx, sizes, point, problem = (prepared[k] for k in ('config', 'idx', 'sizes', 'coordinates', 'problem'))
    primary = cfg['quad_points']
    labels, targets = _targets(cfg, sizes)
    sigma = engine.get_population_prior_sd(cfg)
    joint = np.r_[point, np.log(sigma)]
    person_order = list(prepared['prep']['levels']['Person'])
    person_label_map = {lbl: i for i, lbl in enumerate(problem.person_labels)}
    order = [person_label_map[str(label)] for label in person_order]
    boundary = problem.structural_recession_audit()
    rows = []
    orders = {primary, 31, 61, 121} if comparison_points is None else {primary, comparison_points}
    if any(isinstance(q, (bool, np.bool_)) or operator.index(q) != q or not 3 <= q <= 241 for q in orders):
        raise ValueError('Comparison GH order must be an integer in [3, 241]')
    if max(orders) > max(121, 2*primary-1):
        raise ValueError('Comparison order exceeds the prepared array budget')
    for q in sorted(orders):
        quad = engine.make_mml_quadrature(cfg, q)
        nll, gradient = engine.mfrm_loglik_mml_value_grad(point, idx, cfg, sizes, quad)
        evaluated = problem.evaluate(joint, quad['nodes']/sigma, quad['weights'], information=True)
        if (not np.isfinite(nll) or not np.isfinite(gradient).all()
                or abs(nll-evaluated['nll']) > 1e-9
                or np.max(np.abs(gradient-evaluated['gradient'][:-1])) > 1e-8):
            raise ValueError('Grouped and row likelihood/gradient disagree or are nonfinite')
        if q == primary and (not np.isfinite(prepared['reported_nll']) or abs(nll-prepared['reported_nll']) > 1e-9):
            raise ValueError('Reported NLL does not match the current observations, parameters and fixed rule')
        # Sigma is fixed, not an estimated nuisance parameter. Invert the
        # complete STRUCTURAL matrix, including all cross-facet/step terms.
        hessian = evaluated['information'][:-1, :-1]
        if boundary['structural_design_rank'] < len(point):
            review = dict(status='rank_deficient_design', positive_definite=False, regularized=False,
                          covariance=None, target_se=None, newton_correction=None,
                          target_newton_correction=None, max_target_newton_over_se=None,
                          curvature_scaled_gradient=None, inverse_residual=None)
        else:
            review = _information_review(hessian, gradient, targets)
        eap, sd = evaluated['eap'][order], evaluated['sd'][order]
        if not np.isfinite(eap).all() or not np.isfinite(sd).all():
            raise ValueError('Nonfinite posterior scores')
        rows.append(dict(points=int(q), raw_nll=float(nll), gradient=np.asarray(gradient).tolist(),
                         gradient_supnorm=float(np.max(np.abs(gradient))),
                         nodes=quad['nodes'].tolist(), weights=quad['weights'].tolist(),
                         eap=eap.tolist(), posterior_sd=sd.tolist(),
                         information=hessian.tolist() if np.isfinite(hessian).all() else None,
                         information_review=review))
    baseline = next(row for row in rows if row['points'] == primary)
    base_se = baseline['information_review']['target_se']
    for row in rows:
        se = row['information_review']['target_se']
        relative_se_change = None
        if base_se is not None and se is not None:
            selected = np.asarray(base_se) > 0
            relative_se_change = float(np.max(np.abs(np.asarray(se)[selected]/np.asarray(base_se)[selected]-1)))
        row.update(nll_change=row['raw_nll']-baseline['raw_nll'],
                   max_eap_change=float(np.max(np.abs(np.asarray(row['eap'])-baseline['eap']))),
                   max_posterior_sd_change=float(np.max(np.abs(np.asarray(row['posterior_sd'])-baseline['posterior_sd']))),
                   max_target_relative_se_change=relative_se_change)
    return dict(schema='mfrm_fixed_sd_same_point_gh_v2', input_id=prepared['input_id'],
                input_sha256=prepared['input_sha256'], implementation=prepared['implementation'],
                coordinates=point.tolist(), coordinate_blocks=[dict(name=k, count=int(v)) for k,v in sizes.items()],
                raw_nll=baseline['raw_nll'], reported_nll=prepared['reported_nll'],
                gradient=baseline['gradient'], gradient_supnorm=baseline['gradient_supnorm'], gradient_reference=1e-8,
                quadrature=dict(points=int(primary), sigma=sigma), person_order=person_order,
                targets=labels, target_jacobian=targets.tolist(), target_estimates=(targets @ point).tolist(),
                structural_boundary=boundary,
                comparison=dict(scope='same_point_fixed_sd', primary_points=int(primary), rows=rows,
                                refitted=False, continuous_reference_checked=False, accuracy_pass=None),
                stationarity_pass=None, quadrature_sensitivity_pass=None, inference_ready=False)
