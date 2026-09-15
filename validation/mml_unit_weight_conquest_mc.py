#!/usr/bin/env python3
"""Retain ConQuest scoring variation at four fixed regression PCM calibrations."""
import argparse
import json
from pathlib import Path
import re
import shutil
import sys
import traceback

import numpy as np

from mml_pcm_continuous_refit import dump, sha
from mml_pcm_conquest_mc import COUNTS, SEEDS, rows
from mml_unit_weight_native import CELLS, native_parameters
from run_pcm_conquest_check import ROOT, launch

PRIOR = ROOT / 'validation/generated/mml_unit_weight_native_v2_20260914'
PRIOR_SHA = '0324b33650f22b0d2f0c04acf491cf5ad2e083d32efdae22c8a775845847a54f'
BINARY = Path('/Applications/ConQuest/ConQuest')
FIXED_FILES = ('parameters.csv', 'reg_coefficients.csv', 'covariance.csv',
               'amatrix.csv', 'scored.csv', 'converted.txt', 'history.csv')
INPUT_FILES = ('input.json', 'data.csv', 'anchor_parameters.txt',
               'anchor_beta.txt', 'anchor_covariance.txt')


def errors(values, target):
    delta = np.asarray(values, dtype=float) - np.asarray(target, dtype=float)
    assert delta.shape == (10, 32) and np.all(np.isfinite(delta))
    maxima = np.max(abs(delta), axis=1)
    return dict(aggregate_rmse=float(np.sqrt(np.mean(delta**2))),
                seed_rmse=np.sqrt(np.mean(delta**2, axis=1)),
                seed_max_absolute_error=maxima,
                median_seed_max_absolute_error=float(np.median(maxima)),
                worst_seed_max_absolute_error=float(max(maxima)),
                per_person_mean_signed_error=delta.mean(axis=0),
                per_person_across_seed_sd=delta.std(axis=0, ddof=1))


def self_check():
    assert errors(np.zeros((10, 32)), np.zeros(32))['aggregate_rmse'] == 0
    control = errors(np.tile(np.arange(10.)[:, None], (1, 32)), np.zeros(32))
    assert abs(control['aggregate_rmse'] - np.sqrt(28.5)) < 1e-12
    assert np.array_equal(control['seed_max_absolute_error'], np.arange(10.))
    assert np.allclose(control['per_person_across_seed_sd'], np.std(np.arange(10.), ddof=1))


def verify_hashes(base, hashes):
    for name, expected in hashes.items():
        assert sha(base / name) == expected, name


def verify_prior():
    assert sha(PRIOR / 'summary.json') == PRIOR_SHA
    prior = json.loads((PRIOR / 'summary.json').read_text())
    protocol = json.loads((PRIOR / 'protocol.json').read_text())
    assert prior['implementation_checks_pass']
    assert sha(PRIOR / 'protocol.json') == prior['protocol_sha256']
    verify_hashes(PRIOR, prior['artifact_sha256'])
    verify_hashes(ROOT, protocol['source_sha256'])
    assert sha(BINARY) == protocol['binary_sha256']
    return prior, protocol


def check_run(out, config):
    name = config['cell']; folder = out / config['id']
    previous = PRIOR / name / 'cq_b20_fixed'
    executed = json.loads((folder / 'launch.json').read_text())
    assert executed['exit_code'] == 0 and not executed['timed_out'] and executed['end_of_program'], folder
    # All calibration, design, converted responses, regression means and likelihood
    # history must be byte-identical to the already audited fixed-calibration run.
    for filename in FIXED_FILES:
        assert sha(folder / filename) == sha(previous / filename), folder / filename
    text = (folder / 'review.txt').read_text()
    assert int(re.search(r'Number of nodes used when drawing PVs:\s*(\d+)', text)[1]) == config['p_nodes']
    assert float(re.search(r'Random number generation seed:\s*([\d.]+)', text)[1]) == config['seed']
    assert int(re.search(r'Total number of estimated parameters:\s*(-?\d+)', text)[1]) == 0
    scores = rows(folder / 'cases.csv')
    assert [r['PID'] for r in scores] == [f'P{i:03d}' for i in range(32)]
    values = {key: np.array([float(r[column]) for r in scores])
              for key, column in [('eap', 'EAP_1'), ('sd', 'PosteriorSD_1')]}
    assert all(np.all(np.isfinite(a)) for a in values.values()) and np.all(values['sd'] > 0)
    return dict(scores=values, elapsed=executed['elapsed'], calibration_data_design_history_unchanged=True,
                reported_parameter_count=0,
                termination=[line.strip() for line in text.splitlines() if 'Iterations terminated' in line])


def summarize(out, spec):
    reference = json.loads((out / 'reference.json').read_text())
    runs = {c['id']: check_run(out, c) for c in spec['plan']}
    groups = {}; repeats = {}; scaling = {}
    for cell in CELLS:
        a = runs[f'{cell}/n2000_seed2']['scores']; b = runs[f'{cell}/n2000_seed2_repeat']['scores']
        old_scores = rows(PRIOR / cell / 'cq_b20_fixed/cases.csv')
        repeats[cell] = {key: bool(np.array_equal(a[key], b[key])) for key in a}
        repeats[cell]['previous_baseline_exact'] = all(np.array_equal(a[key], [float(r[col]) for r in old_scores])
            for key, col in [('eap', 'EAP_1'), ('sd', 'PosteriorSD_1')])
        groups[cell] = {}
        for n in COUNTS:
            values = {key: [runs[f'{cell}/n{n}_seed{s}']['scores'][key] for s in SEEDS] for key in a}
            groups[cell][str(n)] = {target: {key: errors(values[key], reference[cell][target][key]) for key in a}
                                   for target in ('continuous_R', 'finite_grid')}
        scaling[cell] = {key: dict(
            aggregate_rmse_ratio_20000_to_2000=groups[cell]['20000']['continuous_R'][key]['aggregate_rmse'] / groups[cell]['2000']['continuous_R'][key]['aggregate_rmse'],
            aggregate_rmse_ratio_200000_to_2000=groups[cell]['200000']['continuous_R'][key]['aggregate_rmse'] / groups[cell]['2000']['continuous_R'][key]['aggregate_rmse'],
            log_log_slope=float(np.polyfit(np.log(COUNTS), np.log([groups[cell][str(n)]['continuous_R'][key]['aggregate_rmse'] for n in COUNTS]), 1)[0])) for key in a}
    assert all(all(r.values()) for r in repeats.values())
    return dict(classification=spec['classification'], scientific_inference_ready=False,
                qualification_eligible=False, implementation_checks_pass=True,
                repeat_exact=repeats, groups=groups, scaling=scaling, runs=runs,
                limitations=spec['limitations'])


def main():
    self_check()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--audit', action='store_true', help='Replay checks from saved output without launching ConQuest')
    args = parser.parse_args(); out = args.output_dir.resolve()
    prior, old_spec = verify_prior()
    if args.audit:
        spec = json.loads((out / 'protocol.json').read_text())
        summary = json.loads((out / 'summary.json').read_text())
        verify_hashes(out, summary['artifact_sha256'])
        verify_hashes(ROOT, spec['source_sha256'])
        assert sha(out / 'protocol.json') == summary['protocol_sha256']
        replay = summarize(out, spec)
        # JSON round trip handles only serialization, with no numerical tolerance.
        assert json.loads(json.dumps(replay, default=lambda x: x.tolist())) == {k: summary[k] for k in replay}
        print('Saved audit passed:', len(summary['runs']), 'runs,', len(summary['artifact_sha256']), 'artifacts,', len(spec['source_sha256']), 'sources')
        return
    sources = dict(old_spec['source_sha256'])
    for filename in (Path(__file__).resolve(), ROOT / 'validation/mml_pcm_continuous_refit.py'):
        sources[str(filename.relative_to(ROOT))] = sha(filename)
    plan = [dict(id=f'{cell}/n{n}_seed{s}', cell=cell, p_nodes=n, seed=s)
            for cell in CELLS for n in COUNTS for s in SEEDS]
    plan += [dict(id=f'{cell}/n2000_seed2_repeat', cell=cell, p_nodes=2000, seed=2) for cell in CELLS]
    out.mkdir(parents=True, exist_ok=False)
    spec = dict(classification='OBSERVED_DEVELOPMENT_ONLY', scientific_inference_ready=False, qualification_eligible=False,
        question='With calibration held fixed, how do ConQuest native EAP and posterior SD vary with posterior draw count and seed on the four retained regression/complete-or-missing unit-weight cells?',
        plan=plan, counts=COUNTS, seeds=SEEDS, prior_summary_sha256=PRIOR_SHA, source_sha256=sources,
        binary_sha256=old_spec['binary_sha256'], conquest_version=old_spec['conquest_version'],
        fixed='Exact continuous-reference anchors, cases mean constraint and regression slope anchor 2; unchanged finite MML grid +/-20 with spacing 0.1 (401 points). Fresh process for every run. Only p_nodes and seed vary.',
        reference='Reuse hash-verified independent scalar R continuous integration at the exact anchors; verify agreement with Python and R refinement below 1e-8. Also retain +/-20 finite-grid score discrepancies.',
        scoring='Summarize 10 by 32 errors per cell/count as aggregate RMSE, per-seed RMSE/max, per-person mean signed error and across-seed SD. Exclude repeat controls from aggregates.',
        failure_policy='Retain every attempt and stop at failed execution/calibration/data/design/history controls; no replacements or retrospective scoring-precision acceptance threshold.',
        sources=['https://conquestmanual.acer.org/s4-00.html'],
        limitations=['Four paired numerical cells from two independent response sets, not four independent replications',
          'Conditional scoring at one fixed calibration per cell; calibration uncertainty, recovery, global optimum, SE/CI/coverage and nonunit weights not evaluated',
          'Native score export rounds to six decimal places; finite-grid truncation is reported separately from the continuous reference',
          'Ten seeds and three counts describe these runs, not a universal 1/sqrt(N) rate or guarantee for every seed/person',
          'All calibration exports/history must match audited baselines; no inference qualification or application-default changes'])
    dump(out / 'protocol.json', spec)
    references = {}
    for cell in CELLS:
        folder = out / cell; folder.mkdir()
        for filename in INPUT_FILES:
            shutil.copyfile(PRIOR / cell / filename, folder / filename)
        value = json.loads((folder / 'input.json').read_text()); p = np.array(value['continuous_point'])
        cases = json.loads((PRIOR / cell / 'r_cases.json').read_text())
        assert np.array_equal(p, cases['cq_b20_fixed']['coordinates'])
        assert np.array_equal(np.loadtxt(folder / 'anchor_parameters.txt')[:, 1], native_parameters(p))
        assert np.array_equal(np.loadtxt(folder / 'anchor_beta.txt'), [1, 2, p[0]])
        assert np.array_equal(np.loadtxt(folder / 'anchor_covariance.txt'), [1, 1, np.exp(2*p[-1])])
        endpoint = prior['results'][cell]['b20_fixed']['endpoints']['cq']
        rc = json.loads((PRIOR / cell / 'cq_b20_fixed_r.json').read_text())['continuous']
        differences = {key: float(np.max(abs(np.asarray(rc[-1][key]) - endpoint['continuous'][key]))) for key in ('nll', 'eap', 'sd')}
        refinement = {key: float(np.max(abs(np.asarray(rc[0][key]) - rc[-1][key]))) for key in differences}
        assert max(differences.values()) < 1e-8 and max(refinement.values()) < 1e-8
        references[cell] = dict(coordinates=p, continuous_R=rc[-1], continuous_python=endpoint['continuous'],
            finite_grid=endpoint['grid'], R_python_difference=differences, R_refinement=refinement,
            finite_continuous_difference={key: float(np.max(abs(np.asarray(endpoint['grid'][key])-rc[-1][key]))) for key in ('eap', 'sd')},
            prior_R_sha256=sha(PRIOR / cell / 'cq_b20_fixed_r.json'))
    dump(out / 'reference.json', references)
    for config in plan:
        folder = out / config['id']; folder.mkdir()
        template = (PRIOR / config['cell'] / 'cq_b20_fixed/model.cqc').read_text()
        assert template.count('p_nodes=2000, seed=2,') == 1
        (folder / 'model.cqc').write_text(template.replace('p_nodes=2000, seed=2,', f'p_nodes={config["p_nodes"]}, seed={config["seed"]},'))
    dump(out / 'prepared.json', dict(artifact_sha256={str(p.relative_to(out)): sha(p) for p in out.rglob('*') if p.is_file()}))
    for i, config in enumerate(plan, 1):
        folder = out / config['id']
        launch(folder, (folder / 'model.cqc').read_text(), 600)
        checked = check_run(out, config)
        print(i, '/', len(plan), config['id'], 'fixed controls pass', round(checked['elapsed'], 3), 's', flush=True)
    verify_prior(); verify_hashes(ROOT, sources)
    verify_hashes(out, json.loads((out / 'prepared.json').read_text())['artifact_sha256'])
    summary = summarize(out, spec)
    summary.update(protocol_sha256=sha(out / 'protocol.json'), artifact_sha256={str(p.relative_to(out)): sha(p) for p in out.rglob('*') if p.is_file()})
    dump(out / 'summary.json', summary)
    print('All fixed controls and exact repeats pass.', flush=True)
    for cell, group in summary['groups'].items():
        for n, g in group.items():
            print(cell, n, {k: dict(rmse=v['aggregate_rmse'], worst=v['worst_seed_max_absolute_error']) for k, v in g['continuous_R'].items()}, flush=True)


if __name__ == '__main__':
    try:
        main()
    except Exception as exc:
        if '--output-dir' in sys.argv and '--audit' not in sys.argv:
            out = Path(sys.argv[sys.argv.index('--output-dir') + 1])
            if out.is_dir() and not (out / 'summary.json').exists():
                dump(out / 'failure.json', dict(error=repr(exc), traceback=traceback.format_exc()))
        raise
