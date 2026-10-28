"""No missing-decision-to-negative shortcut in family error/power evidence."""
import pytest
from validation.bias_operating_metrics import binary_operating_summary,family_rejection


def test_partial_family_with_a_rejection_is_known_but_not_a_complete_family():
    a=family_rejection([True,None,False]);b=family_rejection([False,None,False])
    assert a['any_rejection'] is True and a['complete_family_event'] is None
    assert b['any_rejection'] is None and b['complete_family_event'] is None
    assert a['unavailable_tests']==b['unavailable_tests']==1
    assert family_rejection([False,False])['complete_family_event'] is False


def test_all_planned_interval_accounts_for_every_missing_replicate():
    result=binary_operating_summary([True,False,None,None],planned=4)
    assert result['conditional_rate']==.5
    assert result['all_planned_identification_interval']==[.25,.75]
    lo,hi=result['all_planned_confidence_envelope']
    assert lo<.25 and hi>.75
    assert result['observed']==2 and result['unavailable']==2


def test_complete_data_reduces_to_binomial_and_no_outputs_remain_fully_unknown():
    complete=binary_operating_summary([True,False,False,False],planned=4)
    assert complete['all_planned_identification_interval']==[.25,.25]
    assert complete['all_planned_confidence_envelope']==complete['conditional_mc_interval']
    missing=binary_operating_summary([None,None],planned=2)
    assert missing['conditional_rate'] is None and missing['conditional_mcse'] is None
    assert missing['all_planned_identification_interval']==[0.,1.]
    assert missing['all_planned_confidence_envelope']==[0.,1.]


def test_planned_denominator_and_boolean_domain_cannot_be_silently_changed():
    with pytest.raises(ValueError):binary_operating_summary([False],planned=2)
    with pytest.raises(ValueError):binary_operating_summary([0,1],planned=2)
    with pytest.raises(ValueError):family_rejection([])
    with pytest.raises(ValueError):family_rejection([.01,.04])


def test_draft_plan_has_a_distinct_scope_seed_and_does_not_authorize_execution():
    import numpy as np
    from validation.bias_population_protocol import make_bias_plan
    from validation.mml_population_protocol import make_plan
    plan=make_bias_plan();population=make_plan()
    assert len(plan['scenarios'])==72
    assert sum(s['replicates'] for s in plan['scenarios'])==144000
    assert plan['study_id']!=population['study_id'] and plan['master_seed']!=population['master_seed']
    assert not plan['qualification_execution_allowed'] and not plan['scientific_inference_ready']
    for scenario in plan['scenarios']:
        matrix=np.array(scenario['eta_shift'])
        assert np.all(matrix.sum(axis=0)==0) and np.all(matrix.sum(axis=1)==0)
        assert scenario['null']==bool(np.all(matrix==0))
