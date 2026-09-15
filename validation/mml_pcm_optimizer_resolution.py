#!/usr/bin/env python3
"""Diagnose one retained line-search failure; never promote its qualification."""
import argparse
from collections import Counter
from dataclasses import asdict
import json
from pathlib import Path
import platform
from unittest.mock import patch

import mpmath as mp
import numpy as np
import scipy
from scipy.optimize import minimize, _lbfgsb_py

from mml_pcm_continuous_refit import ROOT, dump, sha, fd_check
from mml_pcm_independent_conditions import PCMIntegral, difference
from mml_pcm_external_review import fixed_grid
from mfrm_app.mml_stationarity import information_diagnostics

SOURCE = ROOT / 'validation/generated/mml_pcm_grid_range_density_20260914'
SOURCE_HASH = '53e68b8425a6f40347c2b57a0a35a4e229e51a4ba5757970c722e08a0736656d'
PROTOCOL = dict(
    id='pcm_optimizer_resolution_20260914_v1',
    classification='OBSERVED_DEVELOPMENT_ONLY', scientific_inference_ready=False,
    qualification_eligible=False,
    question='Does objective rounding obscure local improvement at the retained ABNORMAL point?',
    replay='Exactly one unchanged L-BFGS-B restart; record all callback queries and native tasks; no retries',
    optimizer=dict(maxiter=250, gtol=1e-8, ftol=1e-15, maxls=50),
    sigma_bounds=[0.05, 10.0], hessian_steps=[1e-4, 3e-5], mp_digits=[60, 90],
    cases=['failed ordinary_1 b20_q801 reference', 'successful ordinary_1 b20_q801 zero',
           'off-point control: failed point plus 0.01 in coordinate 0',
           'integration control: successful wide_1 b12_q241 zero'],
    proposals='One diagnostic Newton point from the 1e-4 Hessian, and the noninitial replay query with smallest quadratic predicted NLL change; no proposal is adopted as an estimate',
    multiprecision='Independent PCM scalar likelihood at initial and both proposals, at 60 and 90 decimal digits; exact binary input coordinates/nodes and first float64 grid spacing; unnormalized prior',
    checks=dict(mp_refinement=1e-45, float64_absolute_nll=1e-8,
                gradient_fd=1e-5, hessian_summary=1e-9),
    reporting='Raw success/status/message stay immutable; report requested gtol separately from historical 1e-4 diagnostic gradient budget, local curvature, integration targets and deterministic scoring. All inference flags remain false.',
    scope='Diagnostic extension on observed data, not a new stopping policy or a qualification experiment',
)


def mp_nll(problem, point, theta):
    """Independent PCM reconstruction, avoiding the float64 design/evaluator."""
    # ponytail: this probe covers the retained complete 4-rater, 5-criterion PCM only.
    p = [mp.mpf(float(x)) for x in point]
    raters = p[:3] + [-mp.fsum(p[:3])]
    offsets = []
    for r in range(4):
        for c in range(5):
            steps = p[8+3*c:11+3*c]
            steps = steps + [-mp.fsum(steps)]
            offsets.append([-k*(raters[r]+p[3+c])-mp.fsum(steps[:k]) for k in range(5)])
    sigma = mp.exp(p[-1])
    nodes = [mp.mpf(float(x)) for x in theta]
    log_constant = mp.log(mp.mpf(float(theta[1]-theta[0]))) - p[-1] - mp.log(2*mp.pi)/2
    base = [log_constant-t*t/(2*sigma*sigma) - mp.fsum(
        mp.log(mp.fsum(mp.exp(o[k]+k*t) for k in range(5))) for o in offsets) for t in nodes]
    totals = Counter(map(int, problem.y.sum(axis=1)))
    observed = mp.fsum(offsets[i][int(y[i])] for y in problem.y for i in range(20))
    return -observed - mp.fsum(count*mp.log(mp.fsum(
        mp.exp(total*t+b) for t,b in zip(nodes,base))) for total,count in totals.items())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    out = parser.parse_args().output_dir.resolve()
    assert scipy.__version__ == '1.17.1', 'Private task recorder is version-bound'
    assert sha(SOURCE/'summary.json') == SOURCE_HASH
    old = json.loads((SOURCE/'summary.json').read_text())
    assert old['implementation_checks_pass'] is False
    for name,h in old['artifact_sha256'].items():
        assert sha(SOURCE/name) == h, name
    for name,h in old['source_sha256'].items():
        assert sha(ROOT/name) == h, name
    paths = [Path(__file__).resolve(), ROOT/'mfrm_app/mml_stationarity.py',
             ROOT/'validation/mml_pcm_continuous_refit.py',
             ROOT/'validation/mml_pcm_independent_conditions.py',
             ROOT/'validation/mml_pcm_external_review.py', Path(_lbfgsb_py.__file__),
             Path(_lbfgsb_py._lbfgsb.__file__)]
    hashes = {str(p):sha(p) for p in paths}
    out.mkdir(parents=True, exist_ok=False)
    # Freeze protocol, inputs and implementation before any new numerical evaluation.
    dump(out/'input.json', dict(PROTOCOL, source_summary_sha256=SOURCE_HASH,
        source_sha256=hashes, environment=dict(python=platform.python_version(),
            numpy=np.__version__, scipy=scipy.__version__, mpmath=mp.__version__)))
    load = lambda name: json.loads((SOURCE/name).read_text())
    original = load('ordinary_1/input.json')
    problem = PCMIntegral(original['data'])
    cell = load('ordinary_1/b20_q801.json')
    failure = cell['fit']['starts']['reference']['runs'][1]
    x = np.array(failure['returned'])
    assert not failure['success'] and np.array_equal(x, failure['initial'])
    theta = np.linspace(-20, 20, 801)
    evaluate = lambda p: fixed_grid(problem, p, theta)
    trace, native_tasks = [], []
    def objective(p):
        v = evaluate(p)
        trace.append(dict(coordinates=p.copy(), nll=v['nll'], gradient=v['gradient']))
        return v['nll'], v['gradient']
    original_setulb = _lbfgsb_py._lbfgsb.setulb
    def record_setulb(*args):
        result = original_setulb(*args)
        native_tasks.append(dict(task=args[11].copy(), isave=args[13].copy(),
                                 dsave=args[14].copy(), line_task=args[16].copy()))
        return result
    with patch.object(_lbfgsb_py._lbfgsb, 'setulb', record_setulb):
        opt = minimize(objective, x, jac=True, method='L-BFGS-B',
            bounds=[(None,None)]*23+[tuple(np.log(PROTOCOL['sigma_bounds']))],
            options=PROTOCOL['optimizer'])
    replay = dict(success=bool(opt.success), status=int(opt.status), message=str(opt.message),
        nit=int(opt.nit), nfev=int(opt.nfev), initial=x, returned=opt.x,
        reported_nll=float(opt.fun), returned_value=evaluate(opt.x), queries=trace,
        native_tasks=native_tasks)
    dump(out/'replay.json', replay)
    print('Unchanged replay:', replay['message'], 'queries', opt.nfev, flush=True)

    wide = load('wide_1/b12_q241.json')
    wide_input = load('wide_1/input.json')
    wide_problem = PCMIntegral(wide_input['data'])
    displaced = x.copy(); displaced[0] += .01
    cases = dict(failed=(problem,x,theta,original,failure),
        successful=(problem,np.array(cell['fit']['starts']['zero']['coordinates']),theta,
                    original,cell['fit']['starts']['zero']['runs'][1]),
        displaced=(problem,displaced,theta,original,None),
        integration_control=(wide_problem,np.array(wide['fit']['starts']['zero']['coordinates']),
            np.linspace(-12,12,241),wide_input,wide['fit']['starts']['zero']['runs'][1]))
    diagnostics = {}
    for name,(pr,p,t,inp,raw) in cases.items():
        ev = lambda a: fixed_grid(pr,a,t)
        def value_gradient(a):
            v = ev(a)
            return v['nll'],v['gradient']
        value = ev(p)
        continuous = pr.evaluate(p,**inp['integration'][1])
        error = difference(value,continuous)
        diagnostics[name] = dict(coordinates=p, raw_optimizer=raw,
            gradient_supnorm=float(np.max(abs(value['gradient']))),
            requested_gtol_met=bool(np.max(abs(value['gradient'])) <= 1e-8),
            historical_gradient_budget_met=bool(np.max(abs(value['gradient'])) < 1e-4),
            curvature=[asdict(information_diagnostics(value_gradient,p,relative_step=h))
                       for h in PROTOCOL['hessian_steps']],
            gradient_fd=fd_check(ev,p), integration_error=error,
            engineering_integration_target_met=bool(error['nll']/len(pr.y) < 1e-6 and
                                                    error['eap'] < 1e-4 and error['sd'] < 1e-4),
            scoring_method='deterministic posterior on the stated grid',
            scientific_inference_ready=False)
    dump(out/'diagnostics.json', diagnostics)

    # The existing helper retains summaries; reconstruct H solely for this one proposal.
    steps = 1e-4*np.maximum(1,abs(x))
    directions = np.diag(steps)
    jac = np.column_stack([(evaluate(x+d)['gradient']-evaluate(x-d)['gradient'])/(2*h)
                           for d,h in zip(directions,steps)])
    H = (jac+jac.T)/2
    g = evaluate(x)['gradient']
    correction = np.linalg.solve(H,g)
    assert abs(max(abs(correction))-diagnostics['failed']['curvature'][0]['newton_correction_supnorm']) < 1e-9
    noninitial = [i for i,r in enumerate(trace) if not np.array_equal(r['coordinates'],x)]
    def predicted(i):
        d = trace[i]['coordinates']-x
        return float(g@d+.5*d@H@d)
    selected = min(noninitial,key=predicted)
    proposals = dict(initial=x, newton=x-correction, line_query=trace[selected]['coordinates'])
    dump(out/'proposals.json', dict(coordinates=proposals, line_query_index=selected,
        hessian=H, newton_correction=correction, predicted_newton_gain=float(.5*g@correction),
        predicted_line_gain=-predicted(selected), float64_nll_ulp=float(np.spacing(evaluate(x)['nll']))))
    precise = {}
    for digits in PROTOCOL['mp_digits']:
        with mp.workdps(digits):
            values = {name:mp_nll(problem,p,theta) for name,p in proposals.items()}
            precise[str(digits)] = dict(nll={k:mp.nstr(v,digits) for k,v in values.items()},
                gain={k:mp.nstr(values['initial']-v,digits) for k,v in values.items()},
                float64_gain={k:evaluate(x)['nll']-evaluate(p)['nll'] for k,p in proposals.items()})
        dump(out/f'mp_{digits}.json',precise[str(digits)])
        print('Precision',digits,'Newton gain',precise[str(digits)]['gain']['newton'],flush=True)
    with mp.workdps(90):
        refinement = max(abs(mp.mpf(precise['60']['nll'][k])-mp.mpf(precise['90']['nll'][k])) for k in proposals)
        float_error = max(abs(mp.mpf(float(evaluate(p)['nll']))-mp.mpf(precise['90']['nll'][k])) for k,p in proposals.items())
    checks = dict(mp_refinement=bool(refinement < PROTOCOL['checks']['mp_refinement']),
        mp_float64_formula_agreement=bool(float_error < PROTOCOL['checks']['float64_absolute_nll']),
        gradient_fd=all(v['max_difference'] < PROTOCOL['checks']['gradient_fd']
                        for d in diagnostics.values() for v in d['gradient_fd']),
        failure_retained=not failure['success'] and not old['implementation_checks_pass'],
        gradient_negative_control=not diagnostics['displaced']['historical_gradient_budget_met'],
        integration_negative_control=not diagnostics['integration_control']['engineering_integration_target_met'],
        input_unchanged=sha(SOURCE/'summary.json') == SOURCE_HASH,
        source_unchanged=all(sha(Path(p))==h for p,h in hashes.items()))
    dump(out/'summary.json', dict(classification=PROTOCOL['classification'],
        scientific_inference_ready=False, qualification_eligible=False,
        diagnostic_checks=checks, diagnostic_checks_pass=all(checks.values()),
        historical_implementation_checks_pass=False,
        replay_matches_retained_failure=bool(not opt.success and opt.status==failure['status'] and
            str(opt.message)==failure['message'] and opt.nit==failure['nit'] and
            opt.nfev==failure['nfev'] and np.array_equal(opt.x,x)),
        mp_refinement_max=float(refinement), mp_float64_nll_difference_max=float(float_error),
        protocol_sha256=sha(out/'input.json'),
        artifact_sha256={p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file()}))
    print('Diagnostic checks:',checks,'Historical overall pass: False',flush=True)
    return 0 if all(checks.values()) else 1


if __name__ == '__main__':
    raise SystemExit(main())
