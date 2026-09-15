#!/usr/bin/env python3
"""Isolate iteration-cap and keep-last effects on one retained GH61 start."""
import argparse
import json
from pathlib import Path
import re
import shutil

import numpy as np

from mml_pcm_conquest_direct_start import check_outputs
from mml_pcm_conquest_history_audit import coordinates
from mml_pcm_conquest_mc import rows
from mml_pcm_continuous_refit import dump
from run_pcm_conquest_check import ROOT, launch, sha

SOURCE=ROOT/'validation/generated/mml_pcm_conquest_direct_start_20260914'
SOURCE_HASH='3bd9c5d2982f18b6183295452bc0d315b140fb4cd0733192cb6921ff53c43338'
PLAN=[dict(id='full_keep_no',iterations=2000,keep='no'),
      dict(id='full_keep_yes',iterations=2000,keep='yes'),
      dict(id='three_keep_no',iterations=3,keep='no'),
      dict(id='five_keep_no',iterations=5,keep='no')]


def command(original, config):
    assert original.count('keeplastests=no')==original.count('iterations=2000;')==1
    return original.replace('keeplastests=no',f"keeplastests={config['keep']}").replace(
        'iterations=2000;',f"iterations={config['iterations']};")


def main():
    assert command('keeplastests=no; iterations=2000;',dict(keep='yes',iterations=3))=='keeplastests=yes; iterations=3;'
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory',type=Path,required=True)
    parser.add_argument('--prepare-only',action='store_true')
    args=parser.parse_args();out=args.directory.resolve()
    assert sha(SOURCE/'review/summary.json')==SOURCE_HASH
    source=json.loads((SOURCE/'input.json').read_text())
    execution=json.loads((SOURCE/'execution.json').read_text())
    for name,h in execution['artifact_sha256'].items(): assert sha(SOURCE/name)==h,name
    if args.prepare_only:
        out.mkdir(parents=True,exist_ok=False)
        for name in ('wide.csv','q61_init_parameters.txt','q61_init_covariance.txt'):
            shutil.copyfile(SOURCE/name,out/name)
        original=(SOURCE/'q61_iterations2000/model.cqc').read_text()
        for config in PLAN:
            (out/config['id']).mkdir()
            (out/config['id']/'model.cqc').write_text(command(original,config))
        dump(out/'input.json',dict(classification='OBSERVED_DEVELOPMENT_ONLY',scientific_inference_ready=False,
            qualification_eligible=False,plan=PLAN,source_summary_sha256=SOURCE_HASH,
            binary_sha256=source['binary_sha256'],
            source_sha256=source['source_sha256']|{str(Path(__file__).resolve().relative_to(ROOT)):sha(Path(__file__).resolve())},
            question='Does keep-last only change exports, and what is returned before/after the first improving iteration?',
            controls='Same GH61 imported start, mean constraint, data, seed and stopping tolerances; change only cap or keep-last',
            input_sha256={str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file()}))
        print('Prepared',out);return 0
    spec=json.loads((out/'input.json').read_text())
    for key,base in (('source_sha256',ROOT),('input_sha256',out)):
        for name,h in spec[key].items():assert sha(base/name)==h,name
    assert sha(Path('/Applications/ConQuest/ConQuest'))==spec['binary_sha256']
    (out/'sentinel').mkdir(exist_ok=False)
    sentinel=launch(out/'sentinel','quit;\n',60)
    if sentinel['exit_code']!=0 or not sentinel['end_of_program']:
        print('Startup failed; see',out/'sentinel/console.log');return 1
    original_history=rows(SOURCE/'q61_iterations2000/history.csv');results={}
    for config in spec['plan']:
        name=config['id'];directory=out/name
        result=launch(directory,(directory/'model.cqc').read_text(),600);result['output_checks_pass']=False
        if result['exit_code']==0 and result['end_of_program']:
            try:
                check_outputs(directory)
                history=rows(directory/'history.csv');text=(directory/'review.txt').read_text()
                selected=int(re.search(r'The number of iterations:\s*(\d+)',text)[1])
                assert 1<=selected<=len(history)
                native=np.array([float(r['Estimate']) for r in rows(directory/'parameters.csv')])
                assert np.array_equal(np.r_[native[5:8],native[:5],native[8:]],coordinates(history[selected-1])[:-1])
                assert float(rows(directory/'covariance.csv')[0]['Covariance'])==float(history[selected-1]['wvar 1 1'])
                result.update(output_checks_pass=True,selected_iteration=selected,executed_rows=len(history),
                    history_matches_previous_prefix=history==original_history[:len(history)],
                    reported_nll=float(history[selected-1]['LogLikelihood'])/2,
                    termination=re.search(r'Iterations terminated[^\n]+',text)[0])
            except (OSError,KeyError,ValueError,TypeError,AssertionError) as error:result['review_error']=repr(error)
        results[name]=result;print(name,result,flush=True)
        if not result['output_checks_pass']:break
    dump(out/'summary.json',dict(classification='OBSERVED_DEVELOPMENT_ONLY',scientific_inference_ready=False,
        qualification_eligible=False,runs=results,remaining_unattempted=[c['id'] for c in spec['plan'] if c['id'] not in results],
        artifact_sha256={str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file()}))
    complete=len(results)==4 and all(r['output_checks_pass'] for r in results.values())
    print('Execution complete:',complete,'Saved:',out,'Numerical interpretation pending.')
    return 0 if complete else 1


if __name__=='__main__':raise SystemExit(main())
