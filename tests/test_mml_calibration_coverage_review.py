"""Artificial counts verify interpretation, not any MFRM coverage property."""
from copy import deepcopy
from itertools import product

import numpy as np
import pytest

from validation.mfrm_recovery_metrics import summarize_replicate_errors
from validation.mml_population_analysis import binomial_interval, equivalence_status
from validation.mml_population_missingness_review import review_calibration_coverage


def evidence(planned=2000, available=1600, hits=1520):
    errors = [[0.]]*available+[[np.nan]]*(planned-available)
    covered = [[1.]]*hits+[[0.]]*(available-hits)+[[np.nan]]*(planned-available)
    interval = binomial_interval(hits, available)
    row = dict(kind='Step', facet='Shared', level='1', truth=0., comparable=True,
        metrics=summarize_replicate_errors(errors, covered),
        monte_carlo_intervals=dict(coverage_given_interval=interval,
                                  coverage_per_known_target=binomial_interval(hits, planned)),
        pointwise_acceptance=dict(coverage=equivalence_status(interval, [.925, .975])))
    summary = dict(status='complete', scenario_id='DEVELOPMENT_ONLY', planned=planned,
                   calibration=[row, dict(comparable=False, status='fixed_by_identification')])
    scenario = dict(id='DEVELOPMENT_ONLY', replicates=planned)
    return summary, scenario


def test_conditional_coverage_is_distinct_from_delivery_and_completion_bounds():
    summary, scenario = evidence()
    before = deepcopy(summary)
    packet = review_calibration_coverage(summary, scenario)
    assert summary == before
    assert len(packet['targets']) == 1
    target = packet['targets'][0]
    values = target['conditional_and_completion_evidence']
    assert values['planned'] == 2000 and values['observed'] == 1600 and values['unavailable'] == 400
    assert values['conditional_rate'] == .95
    assert target['delivered_covering_intervals_per_planned']['estimate'] == .76
    assert values['all_planned_identification_interval'] == [.76, .96]
    assert target['registered_conditional_acceptance'] == 'pass'
    assert target['registered_acceptance_unchanged']
    assert not packet['scientific_inference_ready'] and not target['scientific_inference_ready']


@pytest.mark.parametrize('available,hits', [(0, 0), (10, 0), (10, 10), (10, 9), (1, 1)])
def test_zero_and_complete_availability_keep_the_full_denominator(available, hits):
    summary, scenario = evidence(10, available, hits)
    target = review_calibration_coverage(summary, scenario)['targets'][0]
    value = target['conditional_and_completion_evidence']
    assert value['conditional_rate'] == (hits/available if available else None)
    assert value['all_planned_identification_interval'] == [hits/10, (hits+10-available)/10]
    assert target['delivered_covering_intervals_per_planned']['estimate'] == hits/10
    if available == 0:
        assert value['all_planned_confidence_envelope'] == [0., 1.]
        assert target['registered_conditional_acceptance'] == 'unavailable'


def test_completion_bounds_match_all_small_binary_completions():
    summary, scenario = evidence(5, 3, 2)
    value = review_calibration_coverage(summary, scenario)['targets'][0]['conditional_and_completion_evidence']
    means = [(2+sum(missing))/5 for missing in product([0, 1], repeat=2)]
    assert value['all_planned_identification_interval'] == [min(means), max(means)]
    assert value['all_planned_confidence_envelope'][0] == binomial_interval(2, 5)[0]
    assert value['all_planned_confidence_envelope'][1] == binomial_interval(4, 5)[1]


@pytest.mark.parametrize('mutation', ['partial', 'scenario', 'planned', 'fractional_count',
    'fractional_hits', 'delivery_rate', 'target_denominator', 'availability', 'unknown_truth',
    'conditional_interval', 'delivered_interval', 'boolean_count'])
def test_inconsistent_or_person_target_denominators_are_rejected(mutation):
    summary, scenario = evidence()
    row = summary['calibration'][0]
    metrics = row['metrics']
    if mutation == 'partial': summary['status'] = 'incomplete'
    elif mutation == 'scenario': summary['scenario_id'] = 'OTHER'
    elif mutation == 'planned': summary['planned'] += 1
    elif mutation == 'fractional_count': metrics['coverage_given_interval']['denominator'] += .5
    elif mutation == 'fractional_hits': metrics['coverage_given_interval']['estimate'] += .0001
    elif mutation == 'delivery_rate': metrics['coverage_per_known_target']['estimate'] = .95
    elif mutation == 'target_denominator': metrics['coverage_per_known_target']['denominator'] = 1600
    elif mutation == 'availability': metrics['interval_availability']['estimate'] = 1.
    elif mutation == 'unknown_truth': metrics['truth_unknown_replicates'] = 1
    elif mutation == 'conditional_interval': row['monte_carlo_intervals']['coverage_given_interval'] = [0., 1.]
    elif mutation == 'delivered_interval': row['monte_carlo_intervals']['coverage_per_known_target'] = [0., 1.]
    else: metrics['coverage_given_interval']['denominator'] = True
    with pytest.raises(ValueError):
        review_calibration_coverage(summary, scenario)
