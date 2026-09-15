#!/usr/bin/env python3
"""Refit the retained weighted/missing PCM design; separate GH and integral accuracy."""
import argparse
from dataclasses import asdict, replace
from decimal import Decimal, localcontext
import json
from pathlib import Path
import platform
import runpy
import sys

import numpy as np
import scipy
from scipy.integrate import quad_vec
from scipy.optimize import brentq
from scipy.special import logsumexp, ndtr, roots_hermitenorm

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import streamlit_app as app
from mml_pcm_continuous_refit import dump, sha
from mfrm_app.mml_engine_v2 import StationarityContract, assess_free_sd_stationarity, run_free_sd_two_stage
from mfrm_app.mml_quadrature_sensitivity import QuadratureSensitivityContract, assess_quadrature_sensitivity
from mfrm_app.mml_stationarity import JointPolishOptions, audit_joint_gradient

PRIOR = ROOT / 'validation/generated/mml_app_observation_difference_20260914'
API_SPEC = ROOT / 'validation/generated/mml_pcm_two_stage_api_20260914/protocol.json'
TEST_SOURCE = ROOT / 'tests/test_mml_free_sd_stationarity_adapter.py'
PROTOCOL = dict(
    id='observed_pcm_refit_20260914_v1', classification='OBSERVED_DEVELOPMENT_ONLY',
    scientific_inference_ready=False, qualification_eligible=False,
    question='On the retained missing/weighted PCM design, do starts agree at each GH order, and do returned finite-GH solutions approximate continuous likelihood and posterior moments?',
    data='Exactly the previous 147-row, weight-sum 141.5 PCM test design; unchanged responses, deletion and weights; 24 persons including one empty and one wholly zero-weight person.',
    orders=[31, 61, 121, 181], starts=[0.5, 1.0, 3.0],
    structural_start=[0., 0., 0., 0., 0., -1., 0., -1., 0.],
    independent='Explicit affine logits for this two-rater/task/criterion design; SciPy probabilists GH and adaptive standard-normal integration, without app expansion, probability or gradient helpers.',
    continuous=[dict(bound=10., rel_tol=1e-9, abs_tol=1e-11), dict(bound=14., rel_tol=1e-11, abs_tol=1e-13)],
    limits=dict(finite_nll=1e-8, finite_gradient=1e-8, finite_eap=1e-9, finite_sd=1e-9,
                continuous_refinement=1e-8, continuous_fd=1e-6, density_score_identity=1e-8,
                relative_mass_error=1e-7, tail_relative_mass_bound=1e-10,
                cross_start_coordinates=1e-7, cross_start_nll=1e-8,
                integration_nll_per_informative_person=1e-6, integration_eap=1e-4, integration_sd=1e-4,
                continuous_stationarity_gradient=1e-4, no_information_prior=1e-10,
                raw_replay=1e-9, decimal_delta_absolute=1e-18, decimal_delta_relative=1e-11, decimal_refinement=1e-45),
    decimal_selection='Start SD 0.5 at Q31 and Q181: each run own fixed anchor and final endpoint, 60 and 90 digits using the previous independent scalar formula.',
    continuous_selection='Every returned endpoint, two adaptive settings; no continuous refitting or best-start selection.',
    sensitivity='Existing Q121/Q181 raw-NLL assessment at each start; use response-weight sum for its explicitly engineering normalization, not an effective sample size.',
    failure_policy='Freeze sources and settings before preparation/fitting; retain every failed call and failed check without retuning or replacing outputs.',
    scope='Observed development case only. No new native TAM/ConQuest fit, missingness or survey-weight inference, recovery, SE/CI/coverage, public UI or default-estimator change.',
)


class ObservedPCMReference:
    """Independent formula for the specified nine structural coordinates only."""
    def __init__(self, problem):
        idx = problem.idx
        self.person = np.array(idx['person'])
        self.y, self.weight = np.array(idx['score_k']), np.array(idx['weight'])
        self.n = problem.config['n_person']
        self.x = np.array(problem.config['population_model']['X'])[:, 0]
        self.k = np.arange(4.)
        self.design = np.zeros((len(self.y), 4, 9))
        self.offset = np.zeros((len(self.y), 4))
        for i, p in enumerate(self.person):
            r, t, c = (idx['facets'][name][i] for name in ('Rater', 'Task', 'Criterion'))
            assert r in (0, 1) and t in (0, 1) and c in (0, 1)
            self.offset[i] = self.k * (0.35 * (r == 0) - 0.4 * (t == 1))
            self.design[i, :, 0] = self.k * self.x[p]
            self.design[i, :, 1] = self.k * (r == 1)
            self.design[i, :, 2] = self.k * (-1 if t == 0 else 1)
            self.design[i, :, 3+c] = -self.k
            self.design[i, :, 5+2*c:7+2*c] = [[0, 0], [-1, 0], [-1, -1], [0, 0]]
        self.person_weight = np.bincount(self.person, weights=self.weight, minlength=self.n)

    def conditional(self, point, z):
        sigma = np.exp(point[-1])
        logits = self.offset + self.design @ point[:-1] + sigma*z*self.k
        log_prob = logits - logsumexp(logits, axis=1, keepdims=True)
        prob = np.exp(log_prob)
        row_score = self.weight * (self.y - prob @ self.k)
        eta_score = np.bincount(self.person, weights=row_score, minlength=self.n)
        structural = self.weight[:, None] * (self.design[np.arange(len(self.y)), self.y]
                      - np.einsum('ik,ikd->id', prob, self.design))
        gradient = np.zeros((self.n, 10))
        np.add.at(gradient[:, :-1], self.person, structural)
        gradient[:, -1] = sigma*z*eta_score
        ll = np.bincount(self.person, weights=self.weight*log_prob[np.arange(len(self.y)), self.y], minlength=self.n)
        return ll, gradient, eta_score

    def finite(self, point, z, w):
        evaluated = [self.conditional(point, node) for node in z]
        log_joint = np.array([v[0] for v in evaluated]).T + np.log(w)
        marginal = logsumexp(log_joint, axis=1)
        posterior = np.exp(log_joint - marginal[:, None])
        ez = posterior @ z
        variance = np.sum(posterior*(z[None, :]-ez[:, None])**2, axis=1)
        return dict(nll=float(-marginal.sum()),
                    gradient=-np.einsum('pq,qpd->d', posterior, np.array([v[1] for v in evaluated])),
                    eap=self.x*point[0] + np.exp(point[-1])*ez,
                    sd=np.exp(point[-1])*np.sqrt(variance),
                    sigma=float(np.exp(point[-1])))

    def continuous(self, point, bound, rel_tol, abs_tol):
        sigma = np.exp(point[-1])
        modes = np.array([brentq(lambda z: sigma*self.conditional(point, z)[2][p]-z,
                                -bound, bound, xtol=1e-13) for p in range(self.n)])
        assert np.max(abs(modes)) < bound-0.1
        centers = np.array([self.conditional(point, z)[0][p]-z*z/2-np.log(2*np.pi)/2
                            for p, z in enumerate(modes)])

        def integrand(z):
            ll, gradient, _ = self.conditional(point, z)
            mass = np.exp(ll-z*z/2-np.log(2*np.pi)/2-centers)
            return mass[:, None]*np.column_stack((np.ones(self.n), gradient,
                                                  np.full(self.n, z), np.full(self.n, z*z)))

        value, error, info = quad_vec(integrand, -bound, bound, epsrel=rel_tol, epsabs=abs_tol,
                                     norm='max', points=np.unique(np.r_[0., modes]),
                                     quadrature='gk21', workers=1, full_output=True)
        assert info.success and np.isfinite(value).all() and np.all(value[:, 0] > 0)
        normalized = value / value[:, :1]
        marginal = centers + np.log(value[:, 0])
        ez, ez2 = normalized[:, -2], normalized[:, -1]
        assert np.all(ez2 > ez*ez)
        return dict(nll=float(-marginal.sum()), gradient=-normalized[:, 1:11].sum(axis=0),
                    eap=self.x*point[0]+sigma*ez, sd=sigma*np.sqrt(ez2-ez*ez), sigma=float(sigma),
                    density_log_sigma_score=float(-np.sum(ez2-1)),
                    numeric_relative_mass_error_sum=float(np.sum(error/value[:, 0])),
                    tail_relative_mass_bound_sum=float(np.sum(np.exp(np.log(2*ndtr(-bound))-marginal))),
                    quadrature_neval=int(info.neval), quadrature_status=int(info.status))


def discrepancies(left, right, keys=('nll', 'eap', 'sd', 'gradient')):
    return {key: float(np.max(abs(np.asarray(left[key])-np.asarray(right[key])))) for key in keys}


def app_evaluate(problem, point):
    nll, gradient = problem.joint_value_gradient(point)
    quad = app.gauss_hermite_normal(problem.quadrature_points, sd=np.exp(point[-1]))
    expanded = app.expand_params(point[:-1], problem.sizes, problem.config)
    score = app.compute_person_eap(problem.idx, problem.config, expanded, quad)
    return dict(nll=nll, gradient=gradient, eap=score.Estimate.to_numpy(), sd=score.SD.to_numpy())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    out = parser.parse_args().output_dir.resolve()
    load = lambda path: json.loads(path.read_text())
    assert sha(PRIOR/'summary.json') == 'f8b263e206cb6f825c0918896990df71ceabc62c1f109026aeb1dcceac7fbd79'
    assert sha(API_SPEC) == '5d83b66f287a8b75fa23a87b90622885f600cc37b64d80057eaa862788acee5a'
    for path, digest in load(PRIOR/'protocol.json')['sources'].items():
        assert sha(ROOT/path) == digest, path
    api_spec = load(API_SPEC)
    stationarity_spec = dict(api_spec['stationarity'], objective_worsening_tolerance_per_observation=1e-18/141.5)
    stationarity = StationarityContract(**stationarity_spec)
    sensitivity = QuadratureSensitivityContract(**api_spec['sensitivity'])
    paths = [Path(__file__).resolve(), TEST_SOURCE, ROOT/'streamlit_app.py',
             ROOT/'validation/mml_free_sd_stationarity_adapter.py', ROOT/'validation/mml_pcm_continuous_refit.py',
             ROOT/'mfrm_app/mml_stationarity.py', ROOT/'mfrm_app/mml_engine_v2.py',
             ROOT/'mfrm_app/mml_quadrature_sensitivity.py', PRIOR/'summary.json', API_SPEC]
    hashes = {str(p.relative_to(ROOT)): sha(p) for p in paths}
    out.mkdir(parents=True, exist_ok=False)
    dump(out/'protocol.json', dict(PROTOCOL, source_sha256=hashes, stationarity=stationarity_spec,
         sensitivity=asdict(sensitivity), preliminary=api_spec['preliminary'], refinement=api_spec['refinement'],
         anchor_gradient_limit=api_spec['anchor_gradient_limit'],
         runtime=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__)))
    helpers = runpy.run_path(str(TEST_SOURCE))
    native = helpers['native_free_sd_result'].__wrapped__()
    base = helpers['_observed_problem'](native, 'PCM')
    neutral = np.array(PROTOCOL['structural_start'])
    assert sum(base.sizes.values()) == 9 and len(base.idx['score_k']) == 147 and base.observations == 141.5
    base = replace(base, structural_start=tuple(neutral), structural_bounds=((None, None),)*9)
    reference = ObservedPCMReference(base)
    dump(out/'input.json', dict(indices=base.idx, x=reference.x, n_person=reference.n,
         parameter_sizes=dict(base.sizes), person_weight=reference.person_weight,
         affine_design=reference.design, affine_offset=reference.offset, structural_start=neutral))
    rules = {}
    for q in PROTOCOL['orders']:
        z, w = roots_hermitenorm(q)
        w /= np.sqrt(2*np.pi)
        rules[q] = dict(nodes=z, weights=w, app=app.gauss_hermite_normal(q))
    dump(out/'rules.json', rules)
    limits = PROTOCOL['limits']
    control_point = np.r_[np.linspace(.8, -.3, 9), np.log(1.7)]
    continuous = lambda point: reference.continuous(point, **PROTOCOL['continuous'][-1])
    def value_gradient(point):
        value = continuous(point)
        return value['nll'], value['gradient']
    gradient_audit = audit_joint_gradient(lambda point: continuous(point)['nll'], value_gradient,
                                          control_point, relative_step=1e-5)
    decimal_nll = helpers['_decimal_observed_nll'](base, control_point, 60)
    independent_control = reference.finite(control_point, rules[31]['nodes'], rules[31]['weights'])
    controls = dict(continuous_gradient_audit=asdict(gradient_audit),
                    decimal_nll=str(decimal_nll), independent_nll=independent_control['nll'],
                    decimal_nll_error=abs(float(decimal_nll)-independent_control['nll']))
    controls['pass'] = bool(gradient_audit.maximum_absolute_difference < limits['continuous_fd']
                             and controls['decimal_nll_error'] < limits['finite_nll'])
    dump(out/'controls.json', controls)
    results, runs, problems, decimal_checks = {}, {}, {}, []
    for q in PROTOCOL['orders']:
        problem = replace(base, quadrature_points=q)
        problems[q] = problem
        for sd in PROTOCOL['starts']:
            key = f'q{q}_sd{sd:g}'
            run = run_free_sd_two_stage(neutral, sd, problem.value, problem.value_gradient,
                observations=problem.observations, structural_bounds=problem.structural_bounds,
                sigma_bounds=problem.sigma_bounds, difference_factory=problem.likelihood_difference,
                joint_value_gradient=problem.joint_value_gradient,
                constraint_residual_function=problem.constraint_residual,
                anchor_gradient_limit=api_spec['anchor_gradient_limit'],
                preliminary_options=JointPolishOptions(**api_spec['preliminary']),
                refinement_options=JointPolishOptions(**api_spec['refinement']))
            dump(out/(key+'_run.json'), run.to_dict())
            final = run.refinement if run.refinement is not None else run.preliminary
            runs[key] = final
            point = np.array(final.restart_polish.joint_coordinates)
            app_value = app_evaluate(problem, point)
            finite = reference.finite(point, rules[q]['nodes'], rules[q]['weights'])
            integrals = [reference.continuous(point, **settings) for settings in PROTOCOL['continuous']]
            tight = integrals[-1]
            parity = discrepancies(app_value, finite)
            refinement = discrepancies(*integrals)
            integration_error = discrepancies(finite, tight)
            information_free = np.flatnonzero(reference.person_weight == 0)
            prior_errors = [max(float(np.max(abs(v['eap'][information_free]-reference.x[information_free]*point[0]))),
                                float(np.max(abs(v['sd'][information_free]-np.exp(point[-1])))))
                            for v in (app_value, finite, *integrals)]
            assessment = assess_free_sd_stationarity(final, stationarity)
            replay = []
            for stage in (run.preliminary, run.refinement):
                if stage is None:
                    continue
                for polish in (stage.primary_polish, stage.restart_polish):
                    shifted = polish.objective_shift
                    delta = problem.likelihood_difference(np.array(shifted.anchor_coordinates)) if shifted else None
                    for label, coords, expected_raw, expected_delta in (
                        ('initial', polish.initial_joint_coordinates, polish.initial_objective, shifted.initial_difference if shifted else None),
                        ('returned', polish.joint_coordinates, polish.final_objective, shifted.final_difference if shifted else None),
                        ('last_query', polish.last_query_coordinates, polish.last_query_raw_objective, polish.last_query_optimizer_objective if shifted else None)):
                        p = np.array(coords)
                        raw_error = abs(problem.joint_value_gradient(p)[0]-expected_raw)
                        delta_error = abs(delta(p)-expected_delta) if delta else 0.
                        replay.append(dict(stage='refinement' if shifted else 'preliminary', point=label,
                                           raw_error=raw_error, delta_error=delta_error))
            count = int(np.count_nonzero(reference.person_weight))
            checks = dict(finite_arithmetic=all(parity[k] < limits['finite_'+k] for k in parity),
                continuous_reference=all(v < limits['continuous_refinement'] for v in refinement.values())
                    and tight['numeric_relative_mass_error_sum'] < limits['relative_mass_error']
                    and tight['tail_relative_mass_bound_sum'] < limits['tail_relative_mass_bound']
                    and abs(tight['gradient'][-1]-tight['density_log_sigma_score']) < limits['density_score_identity'],
                no_information_prior=max(prior_errors) < limits['no_information_prior'],
                replay=all(max(r['raw_error'], r['delta_error']) < limits['raw_replay'] for r in replay),
                integration_accuracy=integration_error['nll']/count < limits['integration_nll_per_informative_person']
                    and integration_error['eap'] < limits['integration_eap'] and integration_error['sd'] < limits['integration_sd'],
                continuous_stationarity=float(np.max(abs(tight['gradient']))) < limits['continuous_stationarity_gradient'])
            record = dict(coordinates=point, anchor_admitted=run.anchor_admitted, assessment=assessment.to_dict(),
                app=app_value, independent_finite=finite, continuous=integrals, parity=parity,
                continuous_refinement=refinement, integration_error=integration_error,
                informative_persons=count, no_information_prior_errors=prior_errors, replay=replay, checks=checks)
            results[key] = record
            dump(out/(key+'.json'), record)
            print(key, 'sigma', np.exp(point[-1]), 'stationarity', assessment.stationarity_pass,
                  'integration', checks['integration_accuracy'], 'checks', checks, flush=True)
            if sd == .5 and q in (31, 181) and run.refinement is not None:
                anchor = np.array(final.primary_polish.objective_shift.anchor_coordinates)
                precise = {}
                for digits in (60, 90):
                    with localcontext() as context:
                        context.prec = digits
                        a = helpers['_decimal_observed_nll'](problem, anchor, digits)
                        b = helpers['_decimal_observed_nll'](problem, point, digits)
                        precise[str(digits)] = dict(anchor=str(a), final=str(b), delta=str(b-a))
                with localcontext() as context:
                    context.prec = 90
                    target = Decimal(precise['90']['delta'])
                    difference = Decimal.from_float(problem.likelihood_difference(anchor)(point))
                    error = abs(difference-target)
                    tolerance = Decimal(str(limits['decimal_delta_absolute']))+Decimal(str(limits['decimal_delta_relative']))*abs(target)
                    precision_error = max(abs(Decimal(precise['60'][k])-Decimal(precise['90'][k])) for k in ('anchor', 'final', 'delta'))
                check = dict(q=q, precise=precise, error=str(error), precision_error=str(precision_error),
                             passed=bool(error <= tolerance and precision_error < Decimal(str(limits['decimal_refinement']))))
                decimal_checks.append(check)
                dump(out/(key+'_decimal.json'), check)
    cross_start = {}
    for q in PROTOCOL['orders']:
        rows = [results[f'q{q}_sd{sd:g}'] for sd in PROTOCOL['starts']]
        distance = float(np.max(np.ptp([r['coordinates'] for r in rows], axis=0)))
        nll = float(np.ptp([r['app']['nll'] for r in rows]))
        cross_start[str(q)] = dict(coordinates=distance, nll=nll,
            passed=distance < limits['cross_start_coordinates'] and nll < limits['cross_start_nll'])
    sensitivities = {}
    for sd in PROTOCOL['starts']:
        check = assess_quadrature_sensitivity(problem_digest=sha(out/'input.json'),
            primary_quadrature_points=121, sensitivity_quadrature_points=181,
            primary_run=runs[f'q121_sd{sd:g}'], sensitivity_run=runs[f'q181_sd{sd:g}'],
            primary_value_function=problems[121].value, sensitivity_value_function=problems[181].value,
            stationarity_contract=stationarity, sensitivity_contract=sensitivity)
        sensitivities[str(sd)] = check.to_dict()
    checks = dict(all_series_retained=len(results)==12, controls=controls['pass'],
        arithmetic=all(r['checks']['finite_arithmetic'] for r in results.values()),
        references=all(r['checks']['continuous_reference'] for r in results.values()),
        prior_scoring=all(r['checks']['no_information_prior'] for r in results.values()),
        replay=all(r['checks']['replay'] for r in results.values()),
        decimal=len(decimal_checks)==2 and all(c['passed'] for c in decimal_checks),
        raw_reconstruction=all(c['objective_reconstruction_pass'] for c in sensitivities.values()),
        inference_withheld=all(not r['assessment']['InferenceReady'] for r in results.values())
            and all(not c['scientific_inference_ready'] for c in sensitivities.values()),
        sources_unchanged=all(sha(ROOT/p)==h for p, h in hashes.items()))
    dump(out/'summary.json', dict(classification=PROTOCOL['classification'], scientific_inference_ready=False,
        qualification_eligible=False, implementation_checks=checks, implementation_checks_pass=all(checks.values()),
        results=results, cross_start=cross_start, sensitivities=sensitivities, decimal_checks=decimal_checks,
        protocol_sha256=sha(out/'protocol.json'),
        artifact_sha256={str(p.relative_to(out)):sha(p) for p in sorted(out.rglob('*')) if p.is_file()}))
    print('Implementation checks:', checks, 'cross_start:', cross_start, flush=True)
    return 0 if all(checks.values()) else 1


if __name__ == '__main__':
    raise SystemExit(main())
