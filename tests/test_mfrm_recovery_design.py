"""DGP controls, excluded from every future qualification cohort."""
from dataclasses import replace

import numpy as np
import pandas as pd
import pytest

from validation.mfrm_recovery_design import RecoveryDesign, category_probabilities, generate_population_sample, replicate_seed


def design(model='GPCM', assignment='regular', persons=12):
    return RecoveryDesign(model=model, persons=persons, population_sd=1.2,
        rater=(-.4, 0., .4), task=(-.2, .2), criterion=(-.4, .1, .7),
        steps=((-1., .2, .8),) if model == 'RSM' else ((-1., .2, .8), (-.7, -.2, .9), (-.6, -.1, .7)),
        log_slopes=(-.4, .1, .3) if model == 'GPCM' else (), assignment=assignment, rating_min=1)


@pytest.mark.parametrize('model', ['RSM', 'PCM', 'GPCM'])
def test_generator_probabilities_obey_every_adjacent_category_logit(model):
    spec = design(model); theta = np.linspace(-1.5, 1.5, spec.persons)
    p = category_probabilities(spec, theta)
    np.testing.assert_allclose(p.sum(axis=-1), 1., atol=2e-15, rtol=0)
    steps = np.array(spec.steps if model != 'RSM' else spec.steps*3)
    slopes = np.exp(spec.log_slopes) if model == 'GPCM' else np.ones(3)
    for person, rater, task, criterion in np.ndindex(p.shape[:-1]):
        expected = slopes[criterion]*(theta[person]-spec.rater[rater]-spec.task[task]-spec.criterion[criterion]-steps[criterion])
        np.testing.assert_allclose(np.log(p[person,rater,task,criterion,1:]/p[person,rater,task,criterion,:-1]), expected, atol=2e-15)


def test_equal_steps_and_unit_slopes_collapse_to_the_same_dgp():
    rsm = design('RSM'); theta = np.linspace(-1., 1., rsm.persons)
    pcm = replace(rsm, model='PCM', steps=rsm.steps*3)
    gpcm = replace(pcm, model='GPCM', log_slopes=(0., 0., 0.))
    np.testing.assert_array_equal(category_probabilities(rsm, theta), category_probabilities(pcm, theta))
    np.testing.assert_array_equal(category_probabilities(pcm, theta), category_probabilities(gpcm, theta))


@pytest.mark.parametrize('assignment', ['regular', 'sparse', 'unbalanced'])
def test_each_design_draws_population_persons_and_keeps_declared_category_support(assignment):
    spec = design(assignment=assignment)
    a = generate_population_sample(spec, 2**100+42); b = generate_population_sample(spec, 2**100+43)
    assert not np.array_equal(a['theta'], b['theta'])
    assert abs(a['theta'].mean()) > 1e-6  # No per-sample location fixing.
    assert set(a['data']['Person']) == set(a['declared_levels']['Person'])
    assert a['data']['Score'].between(1, 4).all()
    assert a['metadata']['rows_before_missing'] == spec.persons*3*2*3
    assert not a['metadata']['scientific_inference_ready']
    assert a['truth']['PopulationSD'] == 1.2
    np.testing.assert_array_equal(a['truth']['LogSlope'], spec.log_slopes)
    np.testing.assert_allclose(a['truth']['Step'].sum(axis=1), 0., atol=1e-15)


def test_assignment_changes_do_not_change_person_or_response_streams():
    spec = design(); seed = 2**100+7
    regular = generate_population_sample(spec, seed)
    for name in ['sparse', 'unbalanced']:
        sample = generate_population_sample(replace(spec, assignment=name), seed)
        np.testing.assert_array_equal(sample['theta'], regular['theta'])
        common = sample['data'].merge(regular['data'], on=['Person', 'Rater', 'Task', 'Criterion'], validate='one_to_one')
        np.testing.assert_array_equal(common['Score_x'], common['Score_y'])
        assert len(sample['data']) < len(regular['data'])


def test_replicates_are_addressed_independently_of_loop_order_or_requested_count():
    master = 2**120+19
    first = {i: replicate_seed(master, 'GPCM_sparse', i) for i in range(1, 5)}
    extended = {i: replicate_seed(master, 'GPCM_sparse', i) for i in reversed(range(1, 9))}
    assert first == {i: extended[i] for i in first}
    assert len(set(extended.values())) == 8
    assert replicate_seed(master, 'PCM_sparse', 1) != first[1]
    assert replicate_seed(master+2**80, 'GPCM_sparse', 1) != first[1]
    spec = design(assignment='sparse')
    a = generate_population_sample(spec, first[1]); b = generate_population_sample(spec, extended[1])
    pd.testing.assert_frame_equal(a['data'], b['data'])
    np.testing.assert_array_equal(a['theta'], b['theta'])


def test_prefix_of_population_sample_and_design_is_stable_when_person_count_increases():
    small = generate_population_sample(design(assignment='unbalanced', persons=8), 73)
    large = generate_population_sample(design(assignment='unbalanced', persons=16), 73)
    np.testing.assert_array_equal(small['theta'], large['theta'][:8])
    subset = large['data'][large['data']['Person'].isin(small['declared_levels']['Person'])].reset_index(drop=True)
    pd.testing.assert_frame_equal(small['data'], subset)


def test_sparse_missing_facet_levels_are_reported_without_replacement():
    result = generate_population_sample(design(assignment='sparse', persons=1), 3)
    assert len(result['data']) == 8
    assert len(result['missing_levels']['Rater']) == len(result['missing_levels']['Criterion']) == 1
    assert len(result['declared_levels']['Rater']) == len(result['declared_levels']['Criterion']) == 3


@pytest.mark.parametrize('changes', [dict(rater=(-.3, 0., .4)), dict(steps=((-1., .2, .9),)*3),
    dict(log_slopes=(.1, .2, .3)), dict(log_slopes=(-4., 1., 3.)), dict(population_sd=True), dict(persons=True)])
def test_unidentified_or_invalid_truth_is_not_silently_recentered(changes):
    with pytest.raises(ValueError): generate_population_sample(replace(design(), **changes), 42)


def test_binary_steps_are_fixed_zero_and_high_seed_bits_are_not_truncated():
    spec = replace(design(), steps=((0.,),)*3)
    a = generate_population_sample(spec, 42); b = generate_population_sample(spec, 42+2**100)
    assert a['data']['Score'].isin([1, 2]).all()
    np.testing.assert_array_equal(a['truth']['Step'], np.zeros((3, 1)))
    assert not np.array_equal(a['theta'], b['theta'])
    with pytest.raises(ValueError): replicate_seed(True, 'scenario', 1)
