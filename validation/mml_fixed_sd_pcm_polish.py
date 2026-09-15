#!/usr/bin/env python3
"""Separate EM stopping, finite-GH stationarity and integration on two retained fits."""
import argparse
import json
from pathlib import Path
import platform
import sys
from unittest.mock import patch
import zipfile

import numpy as np
import pandas as pd
import scipy
from scipy.integrate import quad_vec
from scipy.optimize import brentq, minimize
from scipy.special import logsumexp, ndtr, roots_hermitenorm

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import streamlit_app as app
from mfrm_app.legacy_compat import pcm_native_runner as native
from mfrm_app.mml_stationarity import audit_joint_gradient, information_diagnostics
from mml_free_sd_stationarity_adapter import _log_average_exp
from mml_pcm_continuous_refit import dump, sha

PRIOR = ROOT / 'validation/generated/pcm_native_runner_20260915/checked_v1'
PRIOR_SHA = '225340194c5a086d475c303a70d4fc8587c0b9f38ad1fa52b0c59fb505305919'
PERM = [2, 3, 0, 1, 4, 5, 6, 7]  # Both directions: app [R,T,C,C,steps] <-> native [C,C,R,T,steps].
OPTIONS = dict(maxiter=500, gtol=1e-9, ftol=0., maxls=50)
SPEC = dict(
    classification='OBSERVED_DEVELOPMENT_ONLY', scientific_inference_ready=False, qualification_eligible=False,
    question='How much of the retained fixed-SD PCM discrepancy is EM stopping, and how much is finite-GH integration?',
    cases=['complete', 'one_missing'], orders=[31, 61, 121], starts=['saved_em', 'zero'],
    model='Unchanged 32-person, three binary facets, PCM 0:3, unit weights, mean 0 and fixed SD 1.1; no anchors/regression/penalty.',
    optimization='For every case/order/start: raw app-NLL L-BFGS-B, independent stable fixed-anchor NLL-difference L-BFGS-B using the same app gradient, then one restart with the same anchor. No SD coordinate is optimized.',
    options=OPTIONS,
    em_resume=dict(maxit=200, reltol=1e-8, record='Every inner M-step result and observed-likelihood trace; no altered inner settings.'),
    continuous=[dict(bound=12., rel_tol=1e-10, abs_tol=1e-12), dict(bound=16., rel_tol=1e-12, abs_tol=1e-14)],
    continuous_selection='Original EM and saved-EM-start polished points at all three orders, both adaptive settings; no continuous refit.',
    independence='Native affine category design, SciPy probabilists GH and adaptive standard-normal integration; no app probability/expansion/gradient helper in the reference.',
    information_steps=[1e-4, 3e-5], gradient_steps=[1e-4, 1e-5],
    engineering_limits=dict(finite_replay=1e-8, gradient_fd=1e-6, stationary_gradient=1e-8,
        reconstruction=1e-9, raw_worsening=1e-9, restart_coordinates=1e-7, cross_start_coordinates=1e-7,
        information_symmetry=1e-6, information_condition=1e6, newton_correction=1e-7,
        continuous_refinement=1e-8, numeric_relative_mass_error=1e-7, tail_relative_mass=1e-10),
    failure_policy='Freeze this specification and sources before evaluation; retain all solver statuses, returned and last queried coordinates, failed checks and artifacts without retuning.',
    limits='Engineering diagnostics, not scientific equivalence margins. Two observed examples cannot establish global optimality, SE/CI/coverage, new native score agreement or public estimator readiness.',
)


class Reference:
    # ponytail: only the preflight-validated eight-coordinate PCM; use a general
    # observation adapter if a later study actually needs other model structures.
    def __init__(self, job):
        self.sigma = job['sigma']
        self.y = np.asarray(job['responses'], dtype=float)
        self.mask = np.isfinite(self.y).astype(float)
        self.score = np.nan_to_num(self.y).astype(int)
        self.a = native.native_design()[:, :, PERM]
        self.observed = (self.a[np.arange(8), self.score] * self.mask[:, :, None]).sum(axis=1)
        self.total = (self.score * self.mask).sum(axis=1)
        self.k = np.arange(4.)

    def conditional(self, point, z):
        logits = self.a @ point + self.sigma*z*self.k
        logp = logits - logsumexp(logits, axis=-1, keepdims=True)
        prob = np.exp(logp)
        ll = (logp[np.arange(8), self.score] * self.mask).sum(axis=1)
        gradient = self.observed - self.mask @ np.einsum('ik,ikd->id', prob, self.a)
        eta_score = self.total - self.mask @ (prob @ self.k)
        return ll, gradient, eta_score

    def finite(self, point, z, w):
        evaluated = [self.conditional(point, t) for t in z]
        joint = np.array([v[0] for v in evaluated]).T + np.log(w)
        marginal = logsumexp(joint, axis=1)
        post = np.exp(joint-marginal[:, None]); mean = post @ z
        return dict(nll=float(-marginal.sum()),
            gradient=-np.einsum('pq,qpd->d', post, np.array([v[1] for v in evaluated])),
            eap=self.sigma*mean, sd=self.sigma*np.sqrt(np.sum(post*(z[None, :]-mean[:, None])**2, axis=1)))

    def difference(self, anchor, z, w):
        logits = (self.a @ anchor)[:, None, :] + self.sigma*z[None, :, None]*self.k
        logp = logits-logsumexp(logits, axis=-1, keepdims=True)
        joint = np.tile(np.log(w), (len(self.y), 1))
        for j in range(8):
            joint += self.mask[:, j, None]*logp[j, :, self.score[:, j]]

        def delta(point):
            change = (self.a @ (point-anchor))[:, None, :]
            normalization = _log_average_exp(logp, np.broadcast_to(change, logp.shape))
            per_person = np.zeros_like(joint)
            for j in range(8):
                per_person += self.mask[:, j, None]*(change[j, 0, self.score[:, j], None]-normalization[j])
            return -float(np.sum(_log_average_exp(joint, per_person)))
        return delta

    def continuous(self, point, bound, rel_tol, abs_tol):
        n = len(self.y)
        modes = np.array([brentq(lambda z: self.sigma*self.conditional(point, z)[2][p]-z,
                                -bound, bound, xtol=1e-13) for p in range(n)])
        assert np.max(abs(modes)) < bound-.1
        centers = np.array([self.conditional(point, z)[0][p]-z*z/2-np.log(2*np.pi)/2
                            for p, z in enumerate(modes)])
        def integrand(z):
            ll, gradient, _ = self.conditional(point, z)
            mass = np.exp(ll-z*z/2-np.log(2*np.pi)/2-centers)
            return mass[:, None]*np.column_stack((np.ones(n), gradient, np.full(n, z), np.full(n, z*z)))
        value, error, info = quad_vec(integrand, -bound, bound, epsrel=rel_tol, epsabs=abs_tol,
            norm='max', points=np.unique(np.r_[0., modes]), quadrature='gk21', workers=1, full_output=True)
        assert info.success and np.isfinite(value).all() and np.all(value[:, 0] > 0)
        normalized = value/value[:, :1]; marginal = centers + np.log(value[:, 0])
        mean, second = normalized[:, -2], normalized[:, -1]
        assert np.all(second > mean*mean)
        return dict(nll=float(-marginal.sum()), gradient=-normalized[:, 1:9].sum(axis=0),
            eap=self.sigma*mean, sd=self.sigma*np.sqrt(second-mean*mean),
            numeric_relative_mass_error_sum=float(np.sum(error/value[:, 0])),
            tail_relative_mass_bound_sum=float(np.sum(np.exp(np.log(2*ndtr(-bound))-marginal))),
            quadrature_neval=int(info.neval), quadrature_status=int(info.status))


def differences(a, b, keys=('nll', 'gradient', 'eap', 'sd')):
    return {k:float(np.max(abs(np.asarray(a[k])-np.asarray(b[k])))) for k in keys}


def restore(folder):
    text = (folder/'analysis_snapshot.json').read_text()
    job = native.prepare_inputs((folder/'data.csv').read_bytes(), text)
    saved = json.loads(text, object_hook=lambda d: float(d['nonfinite']) if set(d)=={'nonfinite'} else d)
    config = saved['config']
    assert config['facet_names'] == ['Rater', 'Task', 'Criterion']
    assert config['positive_facets'] == [] and job['sigma'] == 1.1
    for spec in [config['theta_spec'], *config['facet_specs'].values()]:
        spec['anchors'] = np.asarray(spec['anchors'], dtype=float)
        spec['groups'] = np.asarray(spec['groups'], dtype=object)
    config['population_model']['X'] = np.asarray(config['population_model']['X'], dtype=float)
    data = pd.read_csv(folder/'data.csv', dtype={'Person':str})
    prep = app.prepare_mfrm_data(data, 'Person', config['facet_names'], 'Score', rating_min=0, rating_max=3)
    assert prep['levels'] == saved['levels']
    idx = app.build_indices(prep, step_facet='Criterion'); sizes = app.build_param_sizes(config)
    assert sum(sizes.values()) == 8
    start = np.asarray(saved['returned_free_coordinates'])
    assert np.array_equal(start[PERM], job['xsi'])
    return saved, job, config, idx, sizes, start


def optimize(value_gradient, start, difference=None, anchor_objective=None):
    initial = value_gradient(start)[0]; last = None; reconstruction = 0.
    origin = (initial if anchor_objective is None else anchor_objective) if difference is not None else 0.
    def objective(point):
        nonlocal last, reconstruction
        raw, gradient = value_gradient(point)
        value = raw if difference is None else difference(point)
        reconstruction = max(reconstruction, abs(value+origin-raw))
        last = dict(coordinates=point.copy(), optimizer_objective=float(value), nll=float(raw))
        return value, gradient
    fit = minimize(objective, start, jac=True, method='L-BFGS-B', options=OPTIONS)
    value, gradient = value_gradient(fit.x)
    return dict(initial=start.copy(), returned=fit.x.copy(), initial_nll=float(initial), nll=float(value),
        gradient=gradient, gradient_supnorm=float(max(abs(gradient))), last_query=last,
        success=bool(fit.success), status=int(fit.status), message=str(fit.message), nit=int(fit.nit),
        nfev=int(fit.nfev), optimizer_objective=float(fit.fun), reconstruction_error=float(reconstruction))


def main():
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--output', type=Path, required=True)
    out = parser.parse_args().output.resolve()
    assert sha(PRIOR/'summary.json') == PRIOR_SHA
    old = json.loads((PRIOR/'summary.json').read_text())
    for name, digest in old['artifact_sha256'].items():
        assert sha(PRIOR/name) == digest, name
    sources = ['validation/mml_fixed_sd_pcm_polish.py', 'streamlit_app.py',
        'mfrm_app/legacy_compat/pcm_native_runner.py', 'mfrm_app/mml_stationarity.py',
        'mfrm_app/mml_engine_v2.py', 'validation/mml_free_sd_stationarity_adapter.py',
        'validation/mml_pcm_continuous_refit.py']
    out.mkdir(parents=True, exist_ok=False)
    dump(out/'protocol.json', dict(SPEC, previous_summary_sha256=PRIOR_SHA,
        source_sha256={p:sha(ROOT/p) for p in sources},
        runtime=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__)))
    with zipfile.ZipFile(out/'sources.zip', 'x', compression=zipfile.ZIP_DEFLATED) as archive:
        for p in sources: archive.write(ROOT/p, p)
    limits = SPEC['engineering_limits']; results = {}
    previous_gh = json.loads((PRIOR/'app_same_point_gh.json').read_text())
    for name in SPEC['cases']:
        saved, job, cfg, idx, sizes, start = restore(PRIOR/f'{name}_bundle')
        reference = Reference(job); case = dict(input=job, original_coordinates=start, runs={})
        def app_vg(q):
            quad = app.make_mml_quadrature(cfg, q)
            return lambda p: app.mfrm_loglik_mml_value_grad(p, idx, cfg, sizes, quad)
        vg31 = app_vg(31); inner = []
        def capture(*args, **kwargs):
            fitted = minimize(*args, **kwargs)
            observed = vg31(fitted.x)
            inner.append(dict(initial=np.asarray(args[1]).copy(), returned=fitted.x.copy(),
                success=bool(fitted.success), message=str(fitted.message), options=kwargs['options'],
                expected_nll=float(fitted.fun), expected_gradient=fitted.jac.copy(),
                observed_nll=float(observed[0]), observed_gradient=observed[1], nit=int(fitted.nit)))
            return fitted
        with patch.object(app, 'minimize', capture):
            em = app.mfrm_em_mml(start, idx, cfg, sizes, app.make_mml_quadrature(cfg, 31),
                maxit=SPEC['em_resume']['maxit'], reltol=SPEC['em_resume']['reltol'])
        case['em_resume'] = dict(returned=em.x, nll=em.fun, success=em.success, message=em.message,
            nit=em.nit, ll_trace=em.ll_trace, inner=inner, gradient=vg31(em.x)[1])
        case['continuous_original'] = [reference.continuous(start, **s) for s in SPEC['continuous']]
        continuous = lambda p: reference.continuous(p, **SPEC['continuous'][-1])
        def continuous_vg(p):
            v = continuous(p)
            return v['nll'], v['gradient']
        case['continuous_gradient_audit'] = audit_joint_gradient(lambda p:continuous(p)['nll'],
            continuous_vg, start, relative_step=1e-5).to_dict()
        for q in SPEC['orders']:
            z, w = roots_hermitenorm(q); w /= np.sqrt(2*np.pi)
            assert np.all(w > 0) and abs(w.sum()-1) < 1e-12
            app_quad = app.make_mml_quadrature(cfg, q)
            vg = app_vg(q); original = reference.finite(start, z, w)
            replay = differences(original, previous_gh[name][str(q)], keys=('nll', 'gradient'))
            assert max(replay.values()) < limits['finite_replay']
            for label, first in [('saved_em', start), ('zero', np.zeros(8))]:
                key = f'q{q}_{label}'; raw = optimize(vg, first)
                anchor = raw['returned'].copy()
                # The optimizer uses the app's exact stored nodes/weights; the
                # separate finite reference above uses SciPy's independent rule.
                delta = reference.difference(anchor, app_quad['nodes']/job['sigma'], app_quad['weights'])
                assert delta(anchor) == 0.
                shifted = optimize(vg, anchor, delta)
                restarted = optimize(vg, shifted['returned'], delta, anchor_objective=raw['nll'])
                point = restarted['returned']; finite = reference.finite(point, z, w)
                app_value, app_gradient = vg(point)
                moments = app.compute_person_eap(idx, cfg, app.expand_params(point, sizes, cfg), app.make_mml_quadrature(cfg, q))
                checked = dict(nll=app_value, gradient=app_gradient, eap=moments.Estimate, sd=moments.SD)
                info = [information_diagnostics(vg, point, relative_step=h).to_dict() for h in SPEC['information_steps']]
                audits = [audit_joint_gradient(lambda p:vg(p)[0], vg, first, relative_step=h).to_dict() for h in SPEC['gradient_steps']]
                delta_audit = audit_joint_gradient(delta, lambda p:(delta(p), vg(p)[1]), point, relative_step=1e-5).to_dict()
                record = dict(raw=raw, shifted=shifted, restart=restarted, coordinates=point,
                    finite=finite, finite_app_difference=differences(finite, checked), original_replay=replay,
                    information=info, gradient_audits=audits, difference_gradient_audit=delta_audit,
                    original_to_final_coordinates=float(max(abs(point-start))))
                record['checks'] = dict(finite_replay=max(record['finite_app_difference'].values()) < limits['finite_replay'],
                    stationary=max(abs(app_gradient)) < limits['stationary_gradient'],
                    information=all(i['positive_definite'] and i['condition_number'] < limits['information_condition']
                        and i['relative_symmetry_residual'] < limits['information_symmetry']
                        and i['newton_correction_supnorm'] < limits['newton_correction'] for i in info),
                    gradient_audit=max([a['maximum_absolute_difference'] for a in audits]+[delta_audit['maximum_absolute_difference']]) < limits['gradient_fd'],
                    reconstruction=max(r['reconstruction_error'] for r in [raw, shifted, restarted]) < limits['reconstruction'],
                    nonworsening=all(r['nll']-r['initial_nll'] <= limits['raw_worsening'] for r in [raw, shifted, restarted]),
                    restart=max(abs(restarted['returned']-restarted['initial'])) < limits['restart_coordinates'])
                if label == 'saved_em':
                    record['continuous'] = [reference.continuous(point, **s) for s in SPEC['continuous']]
                    record['continuous_refinement'] = differences(*record['continuous'])
                    record['finite_vs_continuous'] = differences(finite, record['continuous'][-1])
                    record['continuous_reference_checks'] = dict(
                        refinement=max(record['continuous_refinement'].values()) < limits['continuous_refinement'],
                        mass_error=all(v['numeric_relative_mass_error_sum'] < limits['numeric_relative_mass_error'] for v in record['continuous']),
                        tail_mass=all(v['tail_relative_mass_bound_sum'] < limits['tail_relative_mass'] for v in record['continuous']),
                        gradient=case['continuous_gradient_audit']['maximum_absolute_difference'] < limits['gradient_fd'])
                case['runs'][key] = record
                dump(out/f'{name}_{key}.json', record)
                print(name, key, 'gradient', restarted['gradient_supnorm'], 'checks', record['checks'], flush=True)
        case['cross_start'] = {str(q):float(max(abs(case['runs'][f'q{q}_saved_em']['coordinates']-case['runs'][f'q{q}_zero']['coordinates']))) for q in SPEC['orders']}
        case['native_comparison'] = {}
        endpoint = case['runs']['q121_saved_em']['coordinates'][PERM]
        for engine in ['tam', 'conquest']:
            retained = old['results'][f'{name}_q401_b20'][f'{engine}_refit']
            case['native_comparison'][engine] = dict(retained_coordinates=retained['returned'],
                coordinate_difference=float(max(abs(endpoint-np.asarray(retained['returned'])))),
                native_score_method=retained['native_score_method'],
                note='Reuses retained calibration only; no new native execution or Monte Carlo score claim.')
        case['q31_to_q121'] = dict(coordinates=float(max(abs(case['runs']['q31_saved_em']['coordinates']-case['runs']['q121_saved_em']['coordinates']))),
            scores=differences(case['runs']['q31_saved_em']['finite'], case['runs']['q121_saved_em']['finite'], keys=('eap', 'sd')))
        results[name] = case
        dump(out/f'{name}.json', case)
    dump(out/'summary.json', dict(protocol_sha256=sha(out/'protocol.json'), results=results,
        scientific_inference_ready=False, qualification_eligible=False,
        artifact_sha256={str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file()}))


if __name__ == '__main__':
    main()
