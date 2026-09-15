#!/usr/bin/env python3
"""Common-calibration and native-refit checks on retained unit-weight PCM data."""
import argparse
import csv
import json
from pathlib import Path
import re
import subprocess
import sys
from types import SimpleNamespace
import numpy as np
from scipy.special import logsumexp
from mml_observed_pcm_refit import ROOT, ObservedPCMReference, discrepancies
from mml_pcm_continuous_refit import dump, sha
from run_pcm_conquest_check import launch
from mml_pcm_conquest_mc import rows

PRIOR=ROOT/'validation/generated/mml_independent_observation_designs_20260914'
PROBE=ROOT/'validation/generated/mml_unit_native_probe_20260914'
CELLS=[f'{data}_{design}_unit' for data in ('ordinary','wide') for design in ('complete','missing')]
SPEC=dict(classification='OBSERVED_DEVELOPMENT_ONLY',scientific_inference_ready=False,qualification_eligible=False,
    question='With regression means and identification aligned, do native TAM/ConQuest likelihoods and scores agree at common calibration, and do their refits approximate the continuous reference?',
    cells=CELLS,bounds=[12,20],spacing=.1,tam_version='4.3.25',conquest_version='5.47.5',maxiter=2000,
    start=[0.,0.,0.,0.,0.,-1.,0.,-1.,0.,0.],
    phases='Fixed continuous-reference calibration in all four cells, both ranges and both engines. Admit native refits from shared neutral coordinates only after all fixed-case calibration/data/design/NLL checks pass. No new responses or selected replacements.',
    mapping='Native centered rater1=(app_rater2-.35)/2; task1=app_task1-.2; native criterion_c=app_criterion_c+.2-(.35+app_rater2)/2. Native intercept=0, slope=app beta, variance=app sigma squared. PCM steps unchanged with explicit label reordering.',
    tam_fixed='est.variance=TRUE with explicit variance.fixed and variance.inits; est.variance=FALSE standardizes variance to one. Retain the installed 1e-10 diagonal adjustment and verify returned anchors.',
    conquest_fixed='Cases constraint fixes intercept; anchor only regression index 2 plus all structural parameters and variance, avoiding redundant mean anchor. CSV parameter/covariance output is rounded to six decimals; fixed-point likelihood uses requested anchors after checking exports within 5.1e-7.',
    conquest_scores='Native EAP/SD uses posterior MC; retain with seed 2 and 2000 nodes, report discrepancies without treating these draws as deterministic quadrature scores.',
    objective='TAM reports raw density-times-spacing NLL; ConQuest normalizes the finite prior separately for each regression mean. Record raw NLL plus sum_p log M_p and both score comparisons.',
    limits=dict(design=1e-10,tam_anchor=2e-10,cq_anchor=5.1e-7,tam_fixed_nll=1e-7,cq_fixed_nll=2e-6,
      independent=1e-8,score=1e-4,structural=1e-4,log_sigma=1e-4,continuous_loss_per_person=1e-6),
    r_integration=[dict(bound=12.,rel_tol=1e-9,abs_tol=1e-11),dict(bound=16.,rel_tol=1e-11,abs_tol=1e-13)],
    continuous=dict(bound=16.,rel_tol=1e-12,abs_tol=1e-14),
    mean_export_check='TAM verifies exact input Y internally at 1e-15; JSON export is checked at 1e-12 because jsonlite digits=NA serializes about 15 significant digits. Original overstrict export comparison and source retained.',
    sources=['https://alexanderrobitzsch.github.io/TAM/reference/tam.mml.html','https://conquestmanual.acer.org/s4-00.html'],
    failure_policy='Preserve all outputs and unmet numerical targets. Syntax probes are retained separately; no new-engine default, SE/CI/coverage or weighted inference qualification.',
)


def native_parameters(p):
    return np.r_[p[3:5]+.2-(.35+p[1])/2,(p[1]-.35)/2,p[2]-.2,p[5:9]]


def app_parameters(xsi,beta,variance,tam=False):
    z=np.asarray(xsi,dtype=float)
    if tam:z=z[[0,1,2,3,4,6,5,7]]
    r=.35+2*z[2]
    return np.r_[np.asarray(beta,dtype=float).ravel()[1],r,z[3]+.2,z[:2]-.2+(.35+r)/2,z[4:],.5*np.log(float(np.asarray(variance).ravel()[0]))]


def grid(reference,p,bound):
    theta=np.linspace(-bound,bound,int(round(2*bound/SPEC['spacing']))+1)
    mu=reference.x*p[0];sigma=np.exp(p[-1]);k=reference.k
    base=reference.offset+reference.design@p[:-1]-mu[reference.person,None]*k
    logits=base[:,None,:]+theta[None,:,None]*k
    logp=logits-logsumexp(logits,axis=2,keepdims=True)
    ll=np.zeros((reference.n,len(theta)))
    np.add.at(ll,reference.person,logp[np.arange(len(reference.y)),:,reference.y])
    logw=-.5*((theta[None,:]-mu[:,None])/sigma)**2-np.log(sigma)-.5*np.log(2*np.pi)+np.log(SPEC['spacing'])
    lm=logsumexp(ll+logw,axis=1);prior=logsumexp(logw,axis=1);post=np.exp(ll+logw-lm[:,None]);eap=post@theta
    return dict(nll=float(-lm.sum()),normalized_nll=float(-lm.sum()+prior.sum()),log_prior_mass=prior,eap=eap,
                sd=np.sqrt(np.sum(post*(theta[None,:]-eap[:,None])**2,axis=1)))


def run_r(out,cell,mode):
    with (cell/f'r_{mode}.log').open('x') as log:
        result=subprocess.run(['Rscript',str(Path(__file__).with_suffix('.R')),str(out/'protocol.json'),str(cell),mode],stdout=log,stderr=subprocess.STDOUT)
    assert result.returncode==0,cell/f'r_{mode}.log'


def cq_command(kind,bound):
    text=(PROBE/'cq_fixed/model.cqc').read_text()
    # Preserve the already checked data declaration, facet labels and model.
    text=text.replace('Unit-weight regression PCM syntax probe','Retained unit-weight PCM native comparison')
    text=text.replace('nodes=401, minnode=-20, maxnode=20',f'nodes={bound*20+1}, minnode=-{bound}, maxnode={bound}')
    text=text.replace('iterations=3;',f'convergence=0.000000001, deviancechange=0.00000000001, iterations={3 if kind=="fixed" else SPEC["maxiter"]};')
    if kind=='refit':text=text.replace('anchor_parameters','init_parameters').replace('anchor_reg_coefficients','init_reg_coefficients').replace('anchor_covariance','init_covariance').replace('../anchor_beta.txt','../init_beta.txt')
    return text


def review(cell,kind,bound,reference,input_data):
    key=f'b{bound}_{kind}';p0=np.array(input_data['continuous_point'] if kind=='fixed' else SPEC['start'])
    tam=json.loads((cell/f'tam_{key}.json').read_text());cq=cell/f'cq_{key}'
    names=tam['xsi_names'];assert names==['C1','C2','raterR1','taskT1','C1:step1','C2:step1','C1:step2','C2:step2']
    tp=app_parameters(tam['xsi'],tam['beta'],tam['variance'],True)
    ti=app_parameters(tam['initial']['xsi'],tam['initial']['beta'],tam['initial']['variance'],True)
    tl=app_parameters(tam['after'][-1]['xsi'],tam['after'][-1]['beta'],tam['after'][-1]['variance'],True)
    pars=rows(cq/'parameters.csv');coef=rows(cq/'reg_coefficients.csv');cov=rows(cq/'covariance.csv')
    assert [r['Label'] for r in pars]==[r['Label'] for r in rows(PROBE/'cq_structure/parameters.csv')]
    assert [(r['Dimension'],r['Regressor']) for r in coef]==[('1','1'),('1','2')]
    cp=app_parameters([r['Estimate'] for r in pars],[r['Estimate'] for r in coef],float(cov[0]['Covariance']))
    history=rows(cq/'history.csv');last=history[-1]
    cl=app_parameters([last[f'xsi {i}'] for i in range(1,9)],[last['Dim 1 Var 1'],last['Dim 2 Var 2']],float(last['wvar 1 1']))
    assert all(float(r['Estimate'])==0 for r in coef if r['Regressor']=='1') and float(np.array(tam['beta']).ravel()[0])==0
    converted=(cq/'converted.txt').read_text().splitlines()
    assert len(converted)==32
    assert all(line[:4]==f'P{i:03d}' and abs(float(line[12:36])-reference.x[i])<1e-15 for i,line in enumerate(converted))
    scores=rows(cq/'scored.csv');assert [r['PID'] for r in scores]==[f'P{i:03d}' for i in range(32)]==[r['pid'] for r in tam['person']]
    full=np.array(input_data['responses'],dtype=float).reshape(32,8)
    native_z=native_parameters(p0)
    fixed_anchor_checks=dict(tam_initial=float(np.max(abs(ti-p0))),
      tam_parameters=float(np.max(abs(np.array(tam['xsi'])-native_z[[0,1,2,3,4,6,5,7]]))) if kind=='fixed' else None,
      tam_variance=abs(float(np.asarray(tam['variance']).ravel()[0])-np.exp(2*p0[-1])) if kind=='fixed' else None,
      cq_parameters=float(np.max(abs(np.array([float(r['Estimate']) for r in pars])-native_z))) if kind=='fixed' else None,
      cq_beta=abs(float(coef[1]['Estimate'])-p0[0]) if kind=='fixed' else None,
      cq_variance=abs(float(cov[0]['Covariance'])-np.exp(2*p0[-1])) if kind=='fixed' else None)
    # Match responses by actual exported facet labels, never by native column position.
    for j,name in enumerate(tam['item_names']):
        c,r,t=map(int,re.fullmatch(r'C(\d)-raterR(\d)-taskT(\d)',name).groups()); col=(r-1)*4+(t-1)*2+c-1
        assert np.array_equal(np.array(tam['resp'],dtype=float)[:,j],full[:,col],equal_nan=True)
    for name in list(scores[0])[2:]:
        c,r,t=map(int,re.match(r'criterion:(\d).*rater:(\d).*task:(\d)',name).groups());col=(r-1)*4+(t-1)*2+c-1
        actual=np.array([np.nan if row[name] in ('','.', 'NA') else float(row[name]) for row in scores])
        assert np.array_equal(actual,full[:,col],equal_nan=True)
    assert np.max(abs(np.asarray(tam['Y'])[:,1]-reference.x))<1e-12
    # Native A values must generate the same category logits after the affine mapping.
    def base(r,t,c,p):
        cumulative=np.r_[0,p[5+2*c],p[5+2*c]+p[6+2*c],0]
        return np.arange(4)*(([.35,p[1]][r])-([p[2],.4-p[2]][t])-p[3+c])-cumulative
    tam_design=max(float(np.max(abs(np.asarray(tam['A'])[j]@np.asarray(tam['xsi'])-base(r-1,t-1,c-1,tp))))
        for j,name in enumerate(tam['item_names']) for c,r,t in [tuple(map(int,re.fullmatch(r'C(\d)-raterR(\d)-taskT(\d)',name).groups()))])
    am=rows(cq/'amatrix.csv');headers=list(scores[0])[2:];cq_design=0.
    for row in am:
        c,r,t=map(int,re.match(r'criterion:(\d).*rater:(\d).*task:(\d)',headers[int(row['GIN'])-1]).groups())
        values=np.array([float(row[k]) for k in list(row)[2:]])
        cq_design=max(cq_design,abs(values@np.array([float(z['Estimate']) for z in pars])-base(r-1,t-1,c-1,cp)[int(row['Category'])-1]))
    endpoints={}
    for engine,point,lastpoint,reported in [('tam',tp,tl,tam['deviance']/2),('cq',cp,cl,float(history[-1]['LogLikelihood'])/2)]:
        chosen=p0 if kind=='fixed' else point
        value=grid(reference,chosen,bound);continuous=reference.continuous(chosen,**SPEC['continuous'])
        total=discrepancies(value,input_data['continuous_value'],keys=('eap','sd'))
        delta=np.abs(chosen-np.asarray(input_data['continuous_point']))
        loss=(continuous['nll']-input_data['continuous_value']['nll'])/32
        endpoints[engine]=dict(coordinates=chosen,exported_coordinates=point,last_coordinates=lastpoint,initial_coordinates=ti if engine=='tam' else p0,
            grid=value,continuous=continuous,same_point_error=discrepancies(value,continuous,keys=('nll','eap','sd')),
            total_score_error=total,structural=float(max(delta[:-1])),log_sigma=float(delta[-1]),continuous_nll_loss_per_person=loss,
            accuracy_pass=bool(max(total.values())<1e-4 and max(delta)<1e-4 and -1e-9/32<=loss<1e-6),
            reported_nll=reported,reported_difference_to_raw=abs(reported-value['nll']),reported_difference_to_normalized=abs(reported-value['normalized_nll']),
            native_score_difference=discrepancies(value,dict(eap=[r['EAP'] for r in tam['person']],sd=[r['SD.EAP'] for r in tam['person']]) if engine=='tam' else dict(eap=[float(r['EAP_1']) for r in rows(cq/'cases.csv')],sd=[float(r['PosteriorSD_1']) for r in rows(cq/'cases.csv')]),keys=('eap','sd')))
    anchor_pass=(fixed_anchor_checks['tam_initial']<2e-10 and (kind!='fixed' or
        max(fixed_anchor_checks[k] for k in ('tam_parameters','tam_variance'))<SPEC['limits']['tam_anchor'] and
        max(fixed_anchor_checks[k] for k in ('cq_parameters','cq_beta','cq_variance'))<SPEC['limits']['cq_anchor']))
    checks=dict(data_and_means=True,design=max(tam_design,cq_design)<SPEC['limits']['design'],anchors=anchor_pass,
        fixed_scores=kind!='fixed' or max(endpoints['tam']['native_score_difference'].values())<1e-8,
        fixed_likelihood=kind!='fixed' or endpoints['tam']['reported_difference_to_raw']<SPEC['limits']['tam_fixed_nll'] and endpoints['cq']['reported_difference_to_normalized']<SPEC['limits']['cq_fixed_nll'])
    return dict(endpoints=endpoints,checks=checks,implementation_pass=bool(all(checks.values())),anchors=fixed_anchor_checks,
        design_error=dict(tam=tam_design,cq=cq_design),tam_iterations=tam['iter'],tam_reached_cap=tam['reached_cap'],
        cq_iterations=len(history),cq_termination=[line.strip() for line in (cq/'review.txt').read_text().splitlines() if 'Iterations terminated' in line],
        native_coordinate_difference=float(np.max(abs(tp-cp))),deterministic_score_difference=discrepancies(endpoints['tam']['grid'],endpoints['cq']['grid'],keys=('eap','sd')))


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output-dir',type=Path,required=True);out=parser.parse_args().output_dir.resolve()
    assert sha(PRIOR/'summary.json')=='106409fca18611bfdfbfdc7809638f849cbb289f6503e023110131175506f6d8'
    old=json.loads((PRIOR/'summary.json').read_text());assert old['implementation_checks_pass']
    for f,h in old['artifact_sha256'].items():assert sha(PRIOR/f)==h,f
    sources=[Path(__file__).resolve(),Path(__file__).with_suffix('.R'),ROOT/'validation/mml_observed_pcm_refit.py',ROOT/'validation/mml_independent_observation_designs.R',ROOT/'validation/run_pcm_conquest_check.py',ROOT/'validation/mml_pcm_conquest_mc.py',PROBE/'cq_fixed/model.cqc',PROBE/'cq_structure/parameters.csv']
    hashes=dict(json.loads((PRIOR/'protocol.json').read_text())['source_sha256'],**{str(p.relative_to(ROOT)):sha(p) for p in sources})
    for f,h in hashes.items():assert sha(ROOT/f)==h,f
    binary=Path('/Applications/ConQuest/ConQuest')
    spec=dict(SPEC,source_sha256=hashes,prior_summary_sha256=sha(PRIOR/'summary.json'),binary_sha256=sha(binary))
    out.mkdir(parents=True,exist_ok=False);dump(out/'protocol.json',spec)
    subprocess.run(['Rscript','-e',f'writeLines(c(R.version.string,as.character(packageVersion("TAM"))),"{out}/tam_environment.txt"); writeLines(deparse(get("tam_mml_mstep_regression",asNamespace("TAM"))),"{out}/tam_installed_regression.txt")'],check=True)
    inputs={};refs={};results={}
    for name in CELLS:
        folder=out/name;folder.mkdir();original=json.loads((PRIOR/name/'input.json').read_text());source=old['results'][name]
        data=json.loads((PRIOR/(name.split('_')[0]+'_full_input.json')).read_text());flat=[None]*256
        for pos in original['original_rows']:flat[pos]=data['full_responses'][pos]['Score']
        assert set(original['indices']['weight'])=={1.}
        value=dict(original,responses=flat,continuous_point=source['continuous']['from_q181']['coordinates'],continuous_value=source['continuous']['from_q181']['tight'])
        dump(folder/'input.json',value);inputs[name]=value
        refs[name]=ObservedPCMReference(SimpleNamespace(idx=value['indices'],config=dict(n_person=32,population_model=dict(X=np.array(value['x'])[:,None]))))
        with (folder/'data.csv').open('x',newline='') as f:
            writer=csv.writer(f);writer.writerow(['Person','x',*[f'Y{i}' for i in range(1,9)]])
            for i in range(32):writer.writerow([f'P{i:03d}',format(value['x'][i],'.17g'),*['.' if v is None else v for v in flat[8*i:8*i+8]]])
        for kind,p in [('anchor',np.array(value['continuous_point'])),('init',np.array(SPEC['start']))]:
            (folder/f'{kind}_parameters.txt').write_text(''.join(f'{i+1} {v:.17g}\n' for i,v in enumerate(native_parameters(p))))
            (folder/f'{kind}_beta.txt').write_text(f'1 2 {p[0]:.17g}\n')
            (folder/f'{kind}_covariance.txt').write_text(f'1 1 {np.exp(2*p[-1]):.17g}\n')
        results[name]={}
    dump(out/'prepared.json',dict(artifact_sha256={str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file()}))
    for kind in ('fixed','refit'):
        for name in CELLS:
            folder=out/name;run_r(out,folder,kind)
            for bound in SPEC['bounds']:
                key=f'b{bound}_{kind}';cq=folder/f'cq_{key}';cq.mkdir();command=cq_command(kind,bound);(cq/'model.cqc').write_text(command)
                execution=launch(cq,command,600);assert execution['exit_code']==0 and execution['end_of_program'],cq
                result=review(folder,kind,bound,refs[name],inputs[name]);results[name][key]=result;dump(folder/(key+'_review.json'),result)
                print(name,key,'checks',result['checks'],'native coord difference',result['native_coordinate_difference'],flush=True)
        if kind=='fixed':assert all(r['implementation_pass'] for cell in results.values() for r in cell.values()),'Fixed-calibration admission failed; native refits not attempted'
    rchecks={}
    for name in CELLS:
        folder=out/name;cases={f'{engine}_{key}':dict(coordinates=v['coordinates'],bound=int(key.split('_')[0][1:])) for key,r in results[name].items() for engine,v in r['endpoints'].items()}
        dump(folder/'r_cases.json',cases);run_r(out,folder,'grid')
        with (folder/'r_continuous.log').open('x') as log:
            status=subprocess.run(['Rscript',str(ROOT/'validation/mml_independent_observation_designs.R'),str(folder/'input.json'),str(folder/'r_cases.json'),str(folder)],stdout=log,stderr=subprocess.STDOUT)
        assert status.returncode==0,folder/'r_continuous.log'
        rg=json.loads((folder/'r_grid.json').read_text());rchecks[name]={}
        for key,r in results[name].items():
            for engine,v in r['endpoints'].items():
                tag=f'{engine}_{key}';rc=json.loads((folder/(tag+'_r.json')).read_text())['continuous']
                rchecks[name][tag]=dict(grid=discrepancies(v['grid'],rg[tag],keys=('nll','normalized_nll','eap','sd')),
                  continuous=discrepancies(v['continuous'],rc[-1],keys=('nll','eap','sd')),refinement=discrepancies(*rc,keys=('nll','eap','sd')))
    checks=dict(implementation=all(r['implementation_pass'] for c in results.values() for r in c.values()),
        independent_r=all(max(d.values())<SPEC['limits']['independent'] for c in rchecks.values() for r in c.values() for d in r.values()),
        sources_unchanged=all(sha(ROOT/f)==h for f,h in hashes.items()) and sha(binary)==spec['binary_sha256'])
    dump(out/'summary.json',dict(classification=SPEC['classification'],scientific_inference_ready=False,qualification_eligible=False,
        results=results,independent_r=rchecks,checks=checks,implementation_checks_pass=bool(all(checks.values())),protocol_sha256=sha(out/'protocol.json'),
        artifact_sha256={str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file()}))
    print('Study checks:',checks,flush=True)
    return 0 if all(checks.values()) else 1


if __name__=='__main__':
    try:raise SystemExit(main())
    except Exception as exc:
        import traceback
        folder=Path(sys.argv[sys.argv.index('--output-dir')+1])
        if folder.is_dir():dump(folder/'failure.json',dict(error=repr(exc),traceback=traceback.format_exc()))
        raise
