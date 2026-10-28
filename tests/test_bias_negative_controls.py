"""Scalar target existence differs from whole-model threshold existence."""
from pathlib import Path

import pytest

from validation.bias_negative_controls import scalar_controls,native_scalar,execute_controls
from validation.bias_population_protocol import make_bias_plan

ROOT=Path(__file__).resolve().parents[1]
RSCRIPT=Path('/usr/local/bin/Rscript')


def test_closed_form_controls_include_fixed_threshold_and_range_distinctions():
    cases=scalar_controls();assert len(cases)==17
    results={c['id']:native_scalar(c) for c in cases}
    assert all(r['all_checks_pass'] for r in results.values())
    assert results['fixed_steps_missing_middle_categories']['row']['Cell Status']=='interior'
    assert results['root_outside_range']['row']['Cell Status']=='search_range_limited'
    for side in ['lower','upper','both']:
        row=results[f'interval_{side}_outside_range']['row']
        assert row['Cell Status']=='interior' and 'limited by search range' in row['Profile CI Status']


@pytest.mark.skipif(not RSCRIPT.exists(),reason='Independent R unavailable')
def test_independent_R_guards_and_complete_family_adjustments(tmp_path):
    result=execute_controls(tmp_path/'controls',make_bias_plan(),'fixture',ROOT,RSCRIPT,lambda:None)
    assert result['required_native']==result['passed_native']==20
    assert result['required_reference']==result['passed_reference']==47
    assert result['all_checks_pass'] and not result['scientific_inference_ready']
    again=execute_controls(tmp_path/'controls',make_bias_plan(),'fixture',ROOT,RSCRIPT,lambda:None)
    assert result==again


def test_required_reference_counts_match_every_prespecified_model_design_and_lane():
    plan=make_bias_plan();spec=plan['independent_reference'];scenarios=plan['scenarios']
    families=sum(len(s['independent_reference_replicates']) for s in scenarios)*len(plan['scope']['lanes'])
    scalar=sum(len(s['independent_reference_replicates'])*len(s['design']['rater'])*len(s['design']['task'])
        for s in scenarios)*len(plan['scope']['lanes'])
    assert families==spec['nominal_required_families']==576
    assert scalar==spec['nominal_required_scalar_cells']==3456
    assert 2*families==spec['nominal_required_adjustments']==spec['nominal_required_derived_adjustments']
