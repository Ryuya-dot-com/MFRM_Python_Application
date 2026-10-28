"""Identified MFRM population DGP with independently addressable replicates.

No fitting, fitted truth, replacement, or automatic repair of sparse designs.
Every call draws new iid Person effects from the declared population. Facet,
step and slope truths are fixed in the scenario's identified coordinates.
"""
from dataclasses import asdict, dataclass
import hashlib
import json

import numpy as np
import pandas as pd
from scipy.special import softmax


GENERATOR_VERSION = 'mfrm_population_recovery_v1_independent_streams'
STREAMS = dict(person=0, assignment=1, response=2)


def _integer(value, label, *, minimum=0, maximum=None):
    if type(value) is not int or value < minimum or (maximum is not None and value > maximum):
        raise ValueError(label+' must be an integer in the declared range')
    return value


def replicate_seed(master_seed, scenario_id, replicate):
    """A stable 128-bit seed; independent of requested run count and run order."""
    _integer(master_seed, 'master_seed', maximum=2**128-1)
    _integer(replicate, 'replicate', minimum=1)
    if not isinstance(scenario_id, str) or not scenario_id:
        raise ValueError('scenario_id must be a nonempty string')
    payload = json.dumps([GENERATOR_VERSION, master_seed, scenario_id, replicate],
                         ensure_ascii=False, separators=(',', ':')).encode('utf-8')
    return int.from_bytes(hashlib.sha256(payload).digest()[:16], 'big')


@dataclass(frozen=True)
class RecoveryDesign:
    model: str
    persons: int
    population_sd: float
    rater: tuple[float, ...]
    task: tuple[float, ...]
    criterion: tuple[float, ...]
    steps: tuple[tuple[float, ...], ...]
    log_slopes: tuple[float, ...] = ()
    assignment: str = 'regular'
    rating_min: int = 0

    def validate(self):
        if self.model not in {'RSM', 'PCM', 'GPCM'} or self.assignment not in {'regular', 'sparse', 'unbalanced'}:
            raise ValueError('Unsupported response model or assignment')
        _integer(self.persons, 'persons', minimum=1)
        _integer(self.rating_min, 'rating_min', minimum=-(2**31), maximum=2**31-1)
        if isinstance(self.population_sd, (bool, np.bool_)) or not np.isfinite(self.population_sd) or self.population_sd <= 0:
            raise ValueError('population_sd must be positive and finite')
        for name in ['rater', 'task', 'criterion']:
            values = np.asarray(getattr(self, name), dtype=float)
            if values.ndim != 1 or not len(values) or not np.isfinite(values).all():
                raise ValueError(name+' must be a nonempty finite vector')
            if name != 'criterion' and abs(values.sum()) > 1e-12:
                raise ValueError(name+' truth must use native sum-zero coordinates')
        step = np.asarray(self.steps, dtype=float)
        expected_rows = 1 if self.model == 'RSM' else len(self.criterion)
        if step.ndim != 2 or step.shape[0] != expected_rows or step.shape[1] < 1 or not np.isfinite(step).all():
            raise ValueError('steps must have one shared row or one row per criterion')
        if np.max(abs(step.sum(axis=1))) > 1e-12:
            raise ValueError('Step truth must have sum zero in every criterion')
        log_slope = np.asarray(self.log_slopes, dtype=float)
        if self.model == 'GPCM':
            if (log_slope.shape != (len(self.criterion),) or not np.isfinite(log_slope).all()
                or abs(log_slope.sum()) > 1e-12 or np.any(abs(log_slope) > 3)):
                raise ValueError('GPCM log-slopes require sum zero and native bounds [-3, 3]')
        elif log_slope.size:
            raise ValueError('RSM/PCM have fixed unit slopes')


def category_probabilities(design, theta):
    """Explicit adjacent-category logits; no app probability helpers."""
    design.validate()
    theta = np.asarray(theta, dtype=float)
    if theta.shape != (design.persons,) or not np.isfinite(theta).all():
        raise ValueError('Person values differ from the declared population sample')
    steps = np.asarray(design.steps, dtype=float)
    if design.model == 'RSM':
        steps = np.repeat(steps, len(design.criterion), axis=0)
    eta = (theta[:, None, None, None]-np.asarray(design.rater)[None, :, None, None]
           -np.asarray(design.task)[None, None, :, None]-np.asarray(design.criterion)[None, None, None, :])
    cumulative = np.c_[np.zeros(len(steps)), np.cumsum(steps, axis=1)]
    slopes = np.exp(design.log_slopes) if design.model == 'GPCM' else np.ones(len(design.criterion))
    logits = slopes[None, None, None, :, None]*(eta[..., None]*np.arange(steps.shape[1]+1)
              -cumulative[None, None, None, :, :])
    if not np.isfinite(logits).all():
        raise FloatingPointError('Nonfinite category logits in the declared DGP')
    return softmax(logits, axis=-1)


def _assignment_mask(design, rng):
    counts = (len(design.rater), len(design.task), len(design.criterion))
    shape = (design.persons, *counts)
    if design.assignment == 'regular':
        return np.ones(shape, dtype=bool)
    keep = np.zeros(shape, dtype=bool)
    for person in range(design.persons):
        selected = [rng.choice(n, size=min(2, n), replace=False) for n in counts]
        if design.assignment == 'sparse':
            keep[person][np.ix_(*selected)] = True
        else:
            # Two prespecified backbone cells span the available factors for
            # each Person; their selection is independent of theta and scores.
            for j in range(2):
                keep[(person, *(s[j % len(s)] for s in selected))] = True
            probabilities = (np.linspace(.2, .95, counts[0])[:, None, None]
                *np.linspace(.5, 1., counts[1])[None, :, None]
                *np.linspace(.35, 1., counts[2])[None, None, :])
            keep[person] |= rng.random(counts) < probabilities
    return keep


def generate_population_sample(design, seed):
    """Generate one registered replicate, preserving unsupported assignments.

    Missing facet levels are reported and never filled in or resampled. All
    three designs retain every planned Person; sparse/unbalanced generation
    alone does not prove the resulting model is identified.
    """
    design.validate()
    _integer(seed, 'replicate seed', maximum=2**128-1)
    rng = {name: np.random.Generator(np.random.PCG64DXSM(np.random.SeedSequence(seed, spawn_key=(key,))))
           for name, key in STREAMS.items()}
    theta = rng['person'].normal(0., design.population_sd, design.persons)
    probability = category_probabilities(design, theta)
    cdf = probability.cumsum(axis=-1); cdf[..., -1] = 1.
    uniform = rng['response'].random(probability.shape[:-1])
    scores = np.sum(uniform[..., None] >= cdf, axis=-1)+design.rating_min
    keep = _assignment_mask(design, rng['assignment'])
    labels = dict(Person=[f'P{i:08d}' for i in range(design.persons)],
        Rater=[f'R{i}' for i in range(len(design.rater))], Task=[f'T{i}' for i in range(len(design.task))],
        Criterion=[f'C{i}' for i in range(len(design.criterion))])
    coordinates = np.array(np.nonzero(keep)).T
    data = pd.DataFrame({name: np.asarray(levels)[coordinates[:, axis]] for axis, (name, levels) in enumerate(labels.items())})
    data['Score'] = scores[keep]
    missing = {name: sorted(set(levels)-set(data[name])) for name, levels in labels.items()}
    slopes = np.exp(design.log_slopes) if design.model == 'GPCM' else np.ones(len(design.criterion))
    return dict(data=data, theta=theta, declared_levels=labels, missing_levels=missing,
        truth=dict(Person=theta, Rater=np.array(design.rater), Task=np.array(design.task), Criterion=np.array(design.criterion),
            Step=np.array(design.steps), LogSlope=np.array(design.log_slopes) if design.model == 'GPCM' else np.zeros(len(slopes)),
            Slope=slopes, PopulationSD=design.population_sd),
        metadata=dict(generator_version=GENERATOR_VERSION, seed=seed, bit_generator='PCG64DXSM', numpy_version=np.__version__,
            stream_spawn_keys=STREAMS.copy(), design=asdict(design), rows_before_missing=int(keep.size),
            rows_after_missing=int(keep.sum()), person_sampling='fresh_iid_normal_per_replicate_without_sample_centering',
            scientific_inference_ready=False))
