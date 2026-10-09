"""G3 generated effects, target labels and the actual native screen path."""
import json
from pathlib import Path
import subprocess
import numpy as np
import pandas as pd
import pytest

from validation.bias_population_design import (contrast_matrix,shifted_probabilities,
    generate_shifted_sample,screen_from_point)
from validation.mfrm_recovery_design import generate_population_sample,category_probabilities
from validation.mml_population_protocol import design
from validation.mml_population_study import truth_joint
from validation.mml_gpcm_joint_check import prepare_problem
from validation.mml_free_sd_stationarity_adapter import prepare_app_free_sd_problem_from_data


def problem_for(sample,spec):
    if spec.model=='GPCM':return prepare_problem(sample['data'],31,rating_min=spec.rating_min,rating_max=spec.rating_min+len(spec.steps[0]))
    return prepare_app_free_sd_problem_from_data(sample['data'],person_col='Person',facet_cols=['Rater','Task','Criterion'],
        score_col='Score',rating_min=spec.rating_min,rating_max=spec.rating_min+len(spec.steps[0]),model=spec.model,
        noncenter_facet='Criterion',step_facet=None if spec.model=='RSM' else 'Criterion',quadrature_points=31)


@pytest.mark.parametrize('model',['RSM','PCM','GPCM'])
@pytest.mark.parametrize('assignment',['regular','sparse','unbalanced'])
def test_zero_shift_reproduces_base_response_stream_exactly_and_effect_preserves_assignment(model,assignment):
    spec=design(model,assignment,persons=12);seed=913713
    base=generate_population_sample(spec,seed)
    null=generate_shifted_sample(spec,seed,contrast_matrix(3,2,0.))
    effect=generate_shifted_sample(spec,seed,contrast_matrix(3,2,.8))
    pd.testing.assert_frame_equal(base['data'],null['data'])
    pd.testing.assert_frame_equal(base['data'].drop(columns='Score'),effect['data'].drop(columns='Score'))
    np.testing.assert_array_equal(base['truth']['Person'],effect['truth']['Person'])
    assert effect['missing_levels']==base['missing_levels']


@pytest.mark.parametrize('model',['RSM','PCM','GPCM'])
def test_effect_changes_adjacent_logits_by_slope_times_eta_shift(model):
    spec=design(model,'regular',persons=4,categories=6);theta=np.array([-1.,-.1,.4,1.2]);shift=contrast_matrix(3,2,.6)
    null=category_probabilities(spec,theta);effect=shifted_probabilities(spec,theta,shift)
    difference=np.diff(np.log(effect),axis=-1)-np.diff(np.log(null),axis=-1)
    slopes=np.exp(spec.log_slopes) if model=='GPCM' else np.ones(3)
    expected=shift[None,:,:,None,None]*slopes[None,None,None,:,None]
    np.testing.assert_allclose(difference,np.broadcast_to(expected,difference.shape),rtol=0,atol=5e-15)
    with pytest.raises(ValueError,match='sum-zero'):shifted_probabilities(spec,theta,shift+1.)


@pytest.mark.parametrize('model',['RSM','PCM','GPCM'])
def test_native_screen_recovery_target_is_only_attached_to_fully_known_calibration(model):
    spec=design(model,'regular',persons=24);sample=generate_shifted_sample(spec,90311,contrast_matrix(3,2,.8))
    prep,problem=problem_for(sample,spec);point=truth_joint(problem.config,sample['truth']);theta=sample['truth']['Person']
    oracle=screen_from_point(prep,problem,point,theta,sample,basis='known_calibration_known_person')
    plugin=screen_from_point(prep,problem,point,.8*theta,sample,basis='fitted_calibration_EAP')
    assert oracle['observed_family_size']==oracle['declared_cells']==6 and len(oracle['pairwise'])==3
    assert all(c['oracle_target'] for c in oracle['cells'])
    assert all(not c['oracle_target'] and c.get('oracle_error') is None and c.get('oracle_covered') is None for c in plugin['cells'])
    assert [r['DGP_local_contrast'] for r in oracle['pairwise']]==pytest.approx([-1.6,1.6,0.])
    assert all(not r['PublicConclusionAllowed'] for r in oracle['table'])
    assert not oracle['scientific_inference_ready']
    with pytest.raises(ValueError,match='Known Person'):
        screen_from_point(prep,problem,point,.8*theta,sample,basis='known_calibration_known_person')
    changed=point.copy();changed[0]+=.01
    with pytest.raises(ValueError,match='Known calibration'):
        screen_from_point(prep,problem,changed,theta,sample,basis='known_calibration_EAP')
    altered=dict(sample,data=sample['data'].iloc[::-1].reset_index(drop=True))
    with pytest.raises(ValueError,match='Prepared rows'):
        screen_from_point(prep,problem,point,theta,altered,basis='known_calibration_known_person')
    with pytest.raises(ValueError,match='Unknown calibration'):
        screen_from_point(prep,problem,point,theta,sample,basis='independent_calibration_EAP')


def test_weighted_screen_cannot_acquire_a_p_value_and_unobserved_cells_remain_declared():
    spec=design('PCM','regular',persons=24);sample=generate_shifted_sample(spec,90311,contrast_matrix(3,2,0.))
    sample['data']=sample['data'].loc[~(sample['data'].Rater.eq('R0')&sample['data'].Task.eq('T0'))].reset_index(drop=True)
    prep,problem=problem_for(sample,spec);point=truth_joint(problem.config,sample['truth'])
    result=screen_from_point(prep,problem,point,sample['truth']['Person'],sample,
        basis='known_calibration_known_person',weights=np.full(len(sample['data']),1.5))
    assert result['declared_cells']==6 and result['observed_family_size']==5 and not result['family_complete']
    assert sum(c['status']=='unobserved_declared_cell' for c in result['cells'])==1
    assert not any(c['p_available'] for c in result['cells'])
    # A review request for unavailable evidence is not a statistical rejection.
    assert any(c.get('review_flag') for c in result['cells'])
    assert all(c.get('holm_rejection') is None for c in result['cells'])


@pytest.mark.skipif(not Path('/usr/local/bin/Rscript').exists(),reason='Independent installed R reference unavailable')
@pytest.mark.parametrize('model',['RSM','PCM','GPCM'])
@pytest.mark.parametrize('categories',[2,6])
def test_independent_R_scalar_reference_handles_binary_and_six_categories(tmp_path,model,categories):
    spec=design(model,'regular',persons=24,categories=categories)
    sample=generate_shifted_sample(spec,90311,contrast_matrix(3,2,.4));prep,problem=problem_for(sample,spec)
    point=truth_joint(problem.config,sample['truth'])
    screen=screen_from_point(prep,problem,point,sample['truth']['Person'],sample,basis='known_calibration_known_person')
    data=pd.DataFrame(screen['reference_rows']);data.loc[data.Rater.eq('R0')&data.Task.eq('T0')].to_csv(tmp_path/'cell.csv',index=False,float_format='%.17g')
    request=dict(cases=[dict(id='DEVELOPMENT_TEST_ONLY',file='cell.csv',bound=10.)],families=[])
    (tmp_path/'input.json').write_text(json.dumps(request))
    source=Path(__file__).resolve().parents[1]/'validation/bias_conditional_reference.R'
    subprocess.run(['/usr/local/bin/Rscript','--vanilla',str(source),str(tmp_path/'input.json'),str(tmp_path/'R.json')],check=True,capture_output=True,text=True)
    r=json.loads((tmp_path/'R.json').read_bytes())['cases'][0]
    py=next(row for row in screen['table'] if row['FacetA_Level']=='R0' and row['FacetB_Level']=='T0')
    assert py['Cell Status']==r['status']=='interior'
    for key,column in [('estimate','Bias Size'),('se','S.E.'),('lr','LR ChiSq'),('p','LR Prob.')]:
        assert py[column]==pytest.approx(r[key],abs=1e-8)
    np.testing.assert_allclose([py['Profile CI Lower'],py['Profile CI Upper']],r['ci'],rtol=0,atol=1e-7)
