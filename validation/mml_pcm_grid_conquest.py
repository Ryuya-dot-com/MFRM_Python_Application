#!/usr/bin/env python3
"""Prepare and run eight matched-start, fixed-grid ConQuest checks."""
import argparse
import csv
import json
import math
from pathlib import Path

from run_pcm_conquest_check import ROOT, launch, sha, write
from mml_pcm_conquest_mc import SOURCE, SOURCE_HASH, rows

PROTOCOL = ROOT / 'validation/mml_pcm_grid_protocol.json'


def prepare(out):
    spec = json.loads(PROTOCOL.read_text())
    source = ROOT / spec['source_directory']
    assert sha(source/'summary.json') == spec['source_summary_sha256']
    previous = json.loads((source/'summary.json').read_text())
    for name, h in previous['artifact_sha256'].items():
        assert sha(source/name) == h, name
    assert sha(SOURCE/'review_summary.json') == SOURCE_HASH
    assert sha(Path('/Applications/ConQuest/ConQuest')) == spec['conquest']['binary_sha256']
    template = (SOURCE/'cq_b8_q161/model.cqc').read_text()
    out.mkdir(parents=True, exist_ok=False)
    (out/'init_parameters.txt').write_text(''.join(f'{i} 0\n' for i in range(1,24)))
    (out/'init_covariance.txt').write_text('1 1 1\n')
    plan = []
    for dataset in spec['datasets']:
        data = json.loads((source/dataset/'input.json').read_text())['data']
        assert len(data) == 1600
        with (out/f'{dataset}.csv').open('x', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['Person', *[f'Y{i}' for i in range(1,21)]])
            for p in range(80):
                block = data[p*20:(p+1)*20]
                assert all(row['Person'] == f'P{p+1:02d}' for row in block)
                assert [(r['Rater'],r['Criterion']) for r in block] == [(f'R{r}',f'C{c}') for r in range(1,5) for c in range(1,6)]
                writer.writerow([block[0]['Person'], *[row['Score'] for row in block]])
        for bound in (12,20):
            name = f'{dataset}_b{bound}'
            folder = out/name
            folder.mkdir()
            command = template.replace('title Observed PCM integration check;', 'title Independent PCM fixed-grid check;')
            command = command.replace('../wide.csv', f'../{dataset}.csv')
            command = command.replace('p_nodes=2000, exit_on_error=yes;',
                'p_nodes=2000, seed=2, keeplastests=yes, progress=no, exit_on_error=yes;')
            command = command.replace('nodes=161, minnode=-8, maxnode=8',
                                      f'nodes={bound*20+1}, minnode=-{bound}, maxnode={bound}')
            model = 'model criterion + rater + criterion*step;\n'
            assert command.count(model) == 1
            command = command.replace(model, model + 'import init_parameters << ../init_parameters.txt;\n'
                + 'import init_covariance << ../init_covariance.txt;\n')
            assert 'anchor_' not in command and 'keeplastests=yes' in command
            (folder/'model.cqc').write_text(command)
            plan.append(dict(id=name, dataset=dataset, bound=bound, spacing=.1, q=bound*20+1))
    write(out/'input.json', dict(protocol=spec, protocol_sha256=sha(PROTOCOL), plan=plan,
        scientific_inference_ready=False, qualification_eligible=False,
        expected_amatrix_sha256=sha(SOURCE/'cq_b8_q161/amatrix.csv'),
        expected_labels=[r['Label'] for r in rows(SOURCE/'cq_b8_q161/parameters.csv')],
        source_sha256={str(p.relative_to(ROOT)):sha(p) for p in (Path(__file__).resolve(),PROTOCOL,
            ROOT/'validation/run_pcm_conquest_check.py',ROOT/'validation/mml_pcm_conquest_mc.py')},
        input_sha256={str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file()}))
    print('Prepared eight runs:',out,flush=True)


def run(out):
    spec = json.loads((out/'input.json').read_text())
    for key, base in (('source_sha256',ROOT),('input_sha256',out)):
        for p,h in spec[key].items():
            assert sha(base/p) == h,p
    assert sha(Path('/Applications/ConQuest/ConQuest')) == spec['protocol']['conquest']['binary_sha256']
    (out/'sentinel').mkdir(exist_ok=False)
    sentinel = launch(out/'sentinel','quit;\n',60)
    if sentinel['exit_code'] != 0 or not sentinel['end_of_program']:
        print('Startup failed; see',out/'sentinel/console.log'); return 1
    results = {}
    for condition in spec['plan']:
        name = condition['id']; directory = out/name
        result = launch(directory,(directory/'model.cqc').read_text(),600)
        result['output_checks_pass'] = False
        if result['exit_code'] == 0 and result['end_of_program']:
            try:
                assert sha(directory/'amatrix.csv') == spec['expected_amatrix_sha256']
                pars = rows(directory/'parameters.csv')
                assert [r['Label'] for r in pars] == spec['expected_labels']
                assert [int(r['P']) for r in pars] == list(range(1,24))
                assert all(math.isfinite(float(r['Estimate'])) for r in pars)
                beta, variance = rows(directory/'reg_coefficients.csv'), rows(directory/'covariance.csv')
                assert len(beta) == len(variance) == 1 and float(beta[0]['Estimate']) == 0
                assert math.isfinite(float(variance[0]['Covariance'])) and float(variance[0]['Covariance']) > 0
                assert [r['PID'] for r in rows(directory/'cases.csv')] == [f'P{i:02d}' for i in range(1,81)]
                assert rows(directory/'history.csv') and (directory/'review.txt').stat().st_size > 0
                result['output_checks_pass'] = True
            except (OSError,ValueError,KeyError,AssertionError) as error:
                result['check_error'] = repr(error)
        results[name] = result
        print(name,result,flush=True)
        if not result['output_checks_pass']: break
    write(out/'execution.json',dict(scientific_inference_ready=False,qualification_eligible=False,runs=results,
        remaining_unattempted=[c['id'] for c in spec['plan'] if c['id'] not in results],
        artifact_sha256={str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file()}))
    complete = len(results) == 8 and all(r['output_checks_pass'] for r in results.values())
    print('Execution complete:',complete,'Saved:',out,'Numerical review pending.',flush=True)
    return 0 if complete else 1


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory',type=Path,required=True)
    parser.add_argument('--prepare-only',action='store_true')
    args = parser.parse_args()
    if args.prepare_only: prepare(args.directory.resolve())
    else: raise SystemExit(run(args.directory.resolve()))
