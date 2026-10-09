"""Replicate-level binary operating characteristics with unavailable tests.

None is an unknown decision, never a nonrejection. Missing decisions identify
a probability interval. A conditional rate describes the explicitly observed
subset, which need not be representative of all planned repetitions.
"""
import numpy as np
from validation.mml_population_analysis import binomial_interval


def binary_operating_summary(outcomes,*,planned):
    values=list(outcomes)
    if type(planned) is not int or planned<1 or len(values)!=planned:
        raise ValueError('Exactly one outcome slot per planned independent replicate required')
    if any(v is not None and type(v) is not bool for v in values):
        raise ValueError('Use boolean decisions or None for unavailable outcomes')
    observed=sum(v is not None for v in values);hits=sum(v is True for v in values);missing=planned-observed
    rate=hits/observed if observed else None
    lower,upper=hits/planned,(hits+missing)/planned
    # Each bound is a mean of independent replicate indicators. Combining
    # a 97.5% lower and a 97.5% upper bound covers the unidentified event
    # probability of any binary completion consistent with the observed
    # outcomes with at least 95% confidence by the union bound. This does not
    # define a nominal test for samples where the method returned no test.
    lower_ci=binomial_interval(hits,planned);upper_ci=binomial_interval(hits+missing,planned)
    return dict(planned=planned,observed=observed,unavailable=missing,successes=hits,
        conditional_rate=rate,conditional_mcse=float(np.sqrt(rate*(1-rate)/observed)) if observed else None,
        conditional_mc_interval=binomial_interval(hits,observed),
        all_planned_identification_interval=[lower,upper],
        all_planned_confidence_envelope=[lower_ci[0],upper_ci[1]],
        availability=observed/planned,availability_mc_interval=binomial_interval(observed,planned),
        interpretation='Conditional rate is conditional on outcome observability; all-planned bounds cover possible completions of unavailable decisions, not a nominal FWER claim for an undefined test',
        scientific_inference_ready=False)


def family_rejection(decisions):
    """An observed rejection proves the union event despite other missing tests.

    Also return the event restricted to fully available families. Do not pool
    partially observed positive events into a complete-family conditional rate.
    The empty family has no requested operating-characteristic estimand.
    """
    values=list(decisions)
    if not values:raise ValueError('An explicitly targeted nonempty family required')
    if any(v is not None and type(v) is not bool for v in values):
        raise ValueError('Use boolean decisions or None for unavailable tests')
    complete=all(v is not None for v in values)
    event=True if any(v is True for v in values) else False if complete else None
    return dict(any_rejection=event,fully_available=complete,
        complete_family_event=event if complete else None,
        tests=len(values),unavailable_tests=sum(v is None for v in values))
