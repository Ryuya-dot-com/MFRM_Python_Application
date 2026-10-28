"""Prepare an implementation-only amendment over two retained study seals.

Uses their original source archives so unrelated workspace edits cannot enter
the formal studies. Preparation never activates a version or sends a request.
"""
import argparse
import ast
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import zipfile

from validation.mml_execution_amendment import digest


ROOT=Path(__file__).resolve().parents[1]
STUDIES={
    'bias_population_operating_20261004_v1':'08544e10ecf5d5b2e6c5af67506bdac26fc17ebebb0722054b4e682f694875d0',
    'joint_mml_population_20261004_v1':'f49dce11d5d4492b3fd790e743de25505f2e3396a8fbb7500576bbc780984eed',
}
CORE=['validation/mml_grouped_likelihood.py','validation/mml_grouped_fit.py']
PATCHES={
 'validation/mml_population_controller.py':[
    ('from validation.mml_population_analysis import summarize_scenario',
     'from validation.mml_population_analysis import summarize_scenario\nfrom validation.mml_execution_amendment import expected_root, verified_sources'),
    ("def environment_verifier(out,plan,*,check_R=True):\n    hashes=json.loads((out/'source_manifest.json').read_bytes())",
     "def environment_verifier(out,plan,*,check_R=True):\n    hashes=json.loads((out/'source_manifest.json').read_bytes())\n    amended=verified_sources(out,ROOT,OPENSSL,TRUSTED_CA_SHA256)\n    binding=None\n    if amended is not None:hashes,binding=amended"),
    ("        if json_bytes(python_runtime())!=json_bytes(expected_python): raise ValueError('Python source/runtime changed')",
     "        if json_bytes(python_runtime())!=json_bytes(expected_python): raise ValueError('Python source/runtime changed')\n        if binding is not None:\n            if digest(out/'implementation_activation.json')!=binding['activation_sha256']:\n                raise ValueError('Execution activation changed during a run')\n            return binding"),
    ("ROOT!=Path('/private/tmp')/('mfrm-population-'+identity[:16])",
     "ROOT!=expected_root(out,identity,Path('/private/tmp')/('mfrm-population-'+identity[:16]))"),
 ],
 'validation/bias_population_controller.py':[
    ('from validation.mml_population_reference_execution import _once',
     'from validation.mml_population_reference_execution import _once\nfrom validation.mml_execution_amendment import expected_root'),
    ('ROOT!=isolated_root(identity)', 'ROOT!=expected_root(out,identity,isolated_root(identity))'),
 ],
 'validation/mml_population_study.py':[
    ("    verify_environment()\n    ledger.append(run_id, 'generation_started', run)",
     "    implementation=verify_environment()\n    ledger.append(run_id, 'generation_started', run)\n    if implementation is not None:ledger.append(run_id,'execution_implementation',implementation)"),
 ],
 'validation/bias_population_study.py':[
    ("    verify_environment();ledger.append(run_id,'generation_started',run)",
     "    implementation=verify_environment();ledger.append(run_id,'generation_started',run)\n    if implementation is not None:ledger.append(run_id,'execution_implementation',implementation)"),
 ],
}


def patched_sources(original):
    candidate=dict(original)
    for name in CORE+['validation/mml_execution_amendment.py']:
        raw=(ROOT/name).read_bytes()
        if len(raw)!=(ROOT/name).stat().st_size:raise OSError('Incomplete workspace read')
        candidate[name]=raw
    for name,replacements in PATCHES.items():
        if name not in original:continue
        text=original[name].decode()
        for before,after in replacements:
            expected=3 if 'ROOT!=Path' in before else 2 if before=='ROOT!=isolated_root(identity)' else 1
            if text.count(before)!=expected:raise ValueError('Unexpected original source: '+name)
            text=text.replace(before,after)
        candidate[name]=text.encode()
    # Cold starts, optimizer stages and all stopping/acceptance constants stay.
    def settings(raw):
        return {s.targets[0].id:ast.dump(s.value) for s in ast.parse(raw).body
                if isinstance(s,ast.Assign) and isinstance(s.targets[0],ast.Name)
                and s.targets[0].id in {'OPTIONS','LIMITS'}}
    if settings(original[CORE[1]])!=settings(candidate[CORE[1]]):raise ValueError('Solver settings changed')
    for name,raw in candidate.items():
        if name.endswith('.py'):ast.parse(raw)
    return candidate


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('output',type=Path)
    args=parser.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    evidence={}
    review=ROOT/'validation/generated/g4_bias_numerical_equivalence_review_20261006'
    reviewed=json.loads((review/'summary.json').read_bytes())
    review_plan=json.loads((review/'plan.json').read_bytes())
    if not reviewed['implementation_equivalence_supported']:
        raise ValueError('Supplementary common-point equivalence review failed')
    for label,name in [('mml','g4_numerical_efficiency_20261006'),('g3','g4_bias_numerical_efficiency_20261006')]:
        parent=ROOT/'validation/generated'/name
        summary=json.loads((parent/'summary.json').read_bytes())
        if label=='mml' and not summary['all_checks_pass']:
            raise ValueError('Implementation comparison failed: '+name)
        if label=='g3' and (review_plan['parent_summary_sha256']!=digest(parent/'summary.json')
                or review_plan['source_manifest_sha256']!=digest(parent/'source_manifest.json')):
            raise ValueError('Supplementary review belongs to different G3 evidence')
        hashes=json.loads((parent/'source_manifest.json').read_bytes())
        if any(hashes[n]!=digest(ROOT/n) for n in CORE):raise ValueError('Candidate differs from checked implementation')
        for file in ['summary.json','plan.json','source_manifest.json']:
            dst=out/'evidence'/label/file;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(parent/file,dst)
        evidence[label]=dict(directory=name,summary_sha256=digest(parent/'summary.json'),
                             source_manifest_sha256=digest(parent/'source_manifest.json'),
                             original_all_checks_pass=summary['all_checks_pass'])
    for file in ['summary.json','plan.json','bias_numerical_equivalence_review.py']:
        dst=out/'evidence/g3_supplementary_review'/file
        dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(review/file,dst)
    evidence['g3_supplementary_review']=dict(directory=review.name,
        summary_sha256=digest(review/'summary.json'),classification=reviewed['classification'],
        original_strict_comparison_remains_failed=not reviewed['original_all_checks_pass'],
        explanation='Two raw cross-score residual differences at distinct solver endpoints; all original decisions and fitted-output comparisons agree. Common-point likelihood/gradient/information/posterior identity is checked for all 72 cases. This post-hoc review does not erase those differences or qualify sampling inference.')
    for file in ['validation/mml_prepare_efficiency_amendment.py','tests/test_mml_execution_amendment.py',
                 'tests/test_mml_amended_execution.py','tests/test_mml_grouped_likelihood.py']:
        dst=out/'evidence/preparation'/file
        dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(ROOT/file,dst)
    for name in ['mfrm-amendment-guards-20261006','mfrm-amended-execution-g3-20261006',
                 'mfrm-amended-execution-mml-20261006','mfrm-efficient-mml-grouped-regression-20261006',
                 'mfrm-numerical-efficiency-final-regression-20261006']:
        shutil.copyfile(Path('/private/tmp')/(name+'.log'),out/'evidence/preparation'/(name+'.log'))
    amendment=dict(schema='mfrm-implementation-amendment-plan-v1',
        change_scope='implementation_only_same_estimand_design_seeds_solver_options_and_criteria',
        scientific_inference_ready=False,evidence=evidence,studies={},
        activation_policy='Stop the original coordinator before activation. Retain every started attempt, including interruption failures. Resume only unstarted RunIds in the same original ledgers; append implementation binding before generating new input. No retry, sample replacement, seed change or criterion change.',
        qualification_policy='Mixed implementation evidence remains subject to all original gates and independent references; no automatic qualification.',
        reference_policy='Original selected RunIds and original reference/analysis code are retained; references evaluate retained fit points regardless of implementation.')
    for name,identity in STUDIES.items():
        study=ROOT/'validation/generated'/name
        if digest(study/'T0.json')!=identity:raise ValueError('Original seal changed')
        seal=json.loads((study/'T0.json').read_bytes())
        for file,sha in seal['files'].items():
            if digest(study/file)!=sha:raise ValueError('Original artifact changed')
        manifest=json.loads((study/'source_manifest.json').read_bytes())
        with zipfile.ZipFile(study/'sources.zip') as z:
            if set(z.namelist())!=set(manifest):raise ValueError('Original source inventory differs')
            original={file:z.read(file) for file in manifest}
        for file,sha in manifest.items():
            if hashlib.sha256(original[file]).hexdigest()!=sha:raise ValueError('Original source changed')
        candidate=patched_sources(original);directory=out/identity;directory.mkdir()
        hashes={file:hashlib.sha256(raw).hexdigest() for file,raw in candidate.items()}
        (directory/'source_manifest.json').write_text(json.dumps(hashes,sort_keys=True,indent=2)+'\n')
        with zipfile.ZipFile(directory/'sources.zip','x',zipfile.ZIP_DEFLATED) as z:
            for file,raw in candidate.items():z.writestr(file,raw)
        changed={file:dict(before=manifest.get(file),after=sha) for file,sha in hashes.items() if manifest.get(file)!=sha}
        (directory/'changes.json').write_text(json.dumps(changed,sort_keys=True,indent=2)+'\n')
        amendment['studies'][identity]=dict(study_id=name,original_source_manifest_sha256=digest(study/'source_manifest.json'),
            candidate_source_manifest=identity+'/source_manifest.json',candidate_source_archive=identity+'/sources.zip',changes=changed)
    (out/'amendment.json').write_text(json.dumps(amendment,sort_keys=True,indent=2)+'\n')
    shutil.copyfile(study/'cacert.pem',out/'cacert.pem')
    if digest(out/'cacert.pem')!=seal['trusted_ca_sha256']:raise ValueError('Trusted CA differs')
    openssl=Path(seal['openssl_executable'])
    if digest(openssl)!=seal['openssl_sha256']:raise ValueError('OpenSSL changed')
    record=dict(schema='mfrm-numerical-implementation-amendment-v1',created_utc=datetime.now(timezone.utc).isoformat(),
        scientific_inference_ready=False,openssl_sha256=seal['openssl_sha256'],
        files={str(p.relative_to(out)):digest(p) for p in out.rglob('*') if p.is_file()})
    (out/'T0.json').write_text(json.dumps(record,sort_keys=True,indent=2)+'\n')
    subprocess.run([str(openssl),'ts','-query','-data',str(out/'T0.json'),'-sha256','-cert','-out',str(out/'query.tsq')],check=True)
    print(json.dumps(dict(directory=str(out),T0_sha256=digest(out/'T0.json'),status='prepared_not_active')))


if __name__=='__main__':main()
