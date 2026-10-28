"""Retain exact count enumeration and independent R likelihood controls."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import zipfile

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from validation.jmle_two_item_oracle import estimates,exact_distribution,joint_nll,limit_values


def write(path,value):
    with path.open('x') as handle:json.dump(value,handle,indent=2,sort_keys=True,allow_nan=False);handle.write('\n')


def main(out):
    import numpy as np
    import scipy
    out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    definition=dict(classification='DETERMINISTIC_ANALYTIC_CONTROL_NOT_QUALIFICATION_DATA',
        informative_persons=[8,16,32,64,128,256,512,1024,4096,16384],delta=[-.5,0.,.5,1.],
        R_counts=[[1,1],[3,1],[1,3],[40,24]],
        R_limits=dict(point=1e-6,se=1e-6,nll=1e-8,gradient=1e-6),
        population_replicates=0,application_estimator_executed=False,
        scope='Unbounded two-item Rasch profile JMLE and exact CML conditional on informative-person count',
        references=[dict(zotero_key='HYGDM2L9',title='The Numerical Solution of a Set of Conditional Estimation Equations',
            author='E. B. Andersen',year=1972,read_scope='Zotero metadata and public publisher abstract',
            url='https://academic.oup.com/jrsssb/article/34/1/42/7027082'),
            dict(title='The Role of Conditional Likelihoods in Latent Variable Modeling',
                authors=['Anders Skrondal','Sophia Rabe-Hesketh'],year=2022,
                read_scope='Public publisher HTML, Section 6; not a Zotero attachment',
                url='https://link.springer.com/article/10.1007/s11336-021-09816-8')])
    write(out/'definition_before_execution.json',definition)
    sources=['validation/jmle_two_item_oracle.py','validation/jmle_two_item_check.py',
        'validation/jmle_two_item_reference.R','tests/test_jmle_two_item_oracle.py','validation/jmle_two_item_plot.R']
    hashes={}
    with zipfile.ZipFile(out/'sources.zip','x',zipfile.ZIP_DEFLATED) as archive:
        for name in sources:
            data=(ROOT/name).read_bytes();hashes[name]=hashlib.sha256(data).hexdigest();archive.writestr(name,data)
    write(out/'source_manifest.json',hashes)
    with (out/'R_stdout.txt').open('x') as stdout,(out/'R_stderr.txt').open('x') as stderr:
        completed=subprocess.run(['/usr/local/bin/Rscript','--vanilla',str(ROOT/sources[2]),str(out/'R_result.json')],stdout=stdout,stderr=stderr)
    write(out/'R_execution.json',dict(returncode=completed.returncode))
    if completed.returncode:raise RuntimeError('Independent R failed; output retained')
    reference=json.loads((out/'R_result.json').read_bytes());checks=[]
    if [[r['n10'],r['n01']] for r in reference['cases']]!=definition['R_counts']:raise ValueError('R case inventory differs')
    limits=definition['R_limits']
    for row in reference['cases']:
        counts=[row['n10'],row['n01']];closed=estimates(*counts)['JMLE'];point=[0.,0.,closed['estimate']]
        errors=dict(point=float(np.max(abs(np.array(row['point'])-point))),
            se=abs(row['joint_se'][2]-closed['profile_se']),nll=abs(row['nll']-joint_nll(point,*counts)),
            gradient=float(np.max(abs(np.array(row['gradient'])))))
        checks.append(dict(counts=counts,errors=errors,pass_=row['convergence']==0 and
            row['information_min_eigenvalue']>0 and all(errors[k]<=limits[k] for k in limits)))
    write(out/'R_comparison.json',dict(cases=checks,all_checks_pass=all(r['pass_'] for r in checks)))
    if not all(r['pass_'] for r in checks):raise RuntimeError('Independent R mismatch; limits unchanged')
    rows=[exact_distribution(m,d) for d in definition['delta'] for m in definition['informative_persons']]
    write(out/'enumeration.json',dict(rows=rows,limits=limit_values(.5),scientific_inference_ready=False))
    write(out/'runtime.json',dict(Python=sys.version,executable=sys.executable,numpy=np.__version__,
        scipy=scipy.__version__,R=reference['R_version']))
    subprocess.run(['/usr/local/bin/Rscript','--vanilla',str(ROOT/'validation/jmle_two_item_plot.R'),str(out)],check=True)
    manifest={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(out.iterdir()) if p.is_file()}
    write(out/'manifest.json',manifest)
    print(json.dumps(dict(output=str(out),R_controls=len(checks),enumerated_conditions=len(rows),all_R_checks_pass=True,
        Monte_Carlo_replicates=0,scientific_inference_ready=False)))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('output',type=Path)
    main(parser.parse_args().output)
