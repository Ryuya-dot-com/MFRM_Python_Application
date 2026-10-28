"""Global monotonicity certificates, distinct from MML optimizer convergence."""
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from scipy.special import roots_hermitenorm

from tests.test_mml_grouped_likelihood import problem
from validation.mml_grouped_likelihood import GroupedMmlLikelihood


def balanced(model, missing=None, only_first_criterion=False):
    data, config, *_ = problem(model, 'regular')
    data = pd.concat([data.assign(Score=k) for k in range(4)], ignore_index=True)
    if missing is not None:
        selected = data.Score == missing
        if only_first_criterion:
            selected &= data.Criterion == 'C0'
        data.loc[selected, 'Score'] = (missing+1) % 4
    return GroupedMmlLikelihood(data, config)


@pytest.mark.parametrize('model', ['RSM', 'PCM', 'GPCM'])
@pytest.mark.parametrize('missing', [0, 1, 2, 3])
def test_unused_category_has_an_exact_structural_ray_and_increasing_marginal_likelihood(model, missing):
    grouped = balanced(model, missing)
    audit = grouped.structural_recession_audit()
    assert audit['no_finite_mle_certified']
    assert audit['minimum_exact_margin'] >= 0 < audit['maximum_exact_margin']
    assert not audit['finite_mle_existence_established']
    ray = np.r_[audit['direction'], np.zeros(grouped.slope_d+1)]
    ray /= np.max(abs(ray))
    point = np.r_[np.linspace(-.4, .5, grouped.dimension-1), np.log(1.2)]
    z, w = roots_hermitenorm(61); w /= np.sqrt(2*np.pi)
    values = [grouped.evaluate(point+t*ray, z, w) for t in [0., 1., 5.]]
    assert all(b['nll'] < a['nll']-1e-5 for a, b in zip(values, values[1:]))
    assert all(v['gradient']@ray < 0 for v in values)


@pytest.mark.parametrize('model', ['PCM', 'GPCM'])
def test_category_absent_within_step_owner_is_detected_even_if_global_support_is_complete(model):
    grouped = balanced(model, missing=1, only_first_criterion=True)
    assert (grouped.counts.sum(axis=(0, 1)) > 0).all()
    assert grouped.structural_recession_audit()['no_finite_mle_certified']


@pytest.mark.parametrize('model', ['RSM', 'PCM', 'GPCM'])
def test_balanced_support_prohibits_a_structural_ray_but_does_not_prove_joint_existence(model):
    audit = balanced(model).structural_recession_audit()
    assert audit['status'] == 'no_structural_ray_detected'
    assert not audit['no_finite_mle_certified']
    assert not audit['finite_mle_existence_established']
    assert audit['structural_design_rank'] == audit['linear_dimension']


def test_rsm_shared_steps_do_not_require_every_category_within_every_criterion():
    assert balanced('RSM', missing=1, only_first_criterion=True).structural_recession_audit()['status'] == 'no_structural_ray_detected'


def test_a_tiny_directional_score_does_not_establish_a_finite_mle():
    grouped = balanced('PCM', missing=1)
    audit = grouped.structural_recession_audit()
    ray = np.r_[audit['direction'], np.zeros(grouped.slope_d+1)]
    ray /= np.max(abs(ray))
    z, w = roots_hermitenorm(61); w /= np.sqrt(2*np.pi)
    point = 100*ray
    result = grouped.evaluate(point, z, w)
    assert abs(result['gradient']@ray) < 1e-8
    assert audit['no_finite_mle_certified']


def test_lp_failure_and_inexact_candidate_cannot_be_promoted(monkeypatch):
    import validation.mml_grouped_likelihood as module
    grouped = balanced('PCM')
    monkeypatch.setattr(module, 'linprog', lambda *a, **k: SimpleNamespace(success=False, status=2, message='test', x=None))
    assert grouped.structural_recession_audit()['status'] == 'unavailable'
    # Pretend a solver returned a feasible nonzero direction in balanced data.
    monkeypatch.setattr(module, 'linprog', lambda *a, **k: SimpleNamespace(success=True, status=0, message='test', x=np.ones(grouped.linear_d)))
    audit = grouped.structural_recession_audit()
    assert audit['reason'] == 'candidate_failed_exact_inequality_check'
    assert not audit['no_finite_mle_certified']
