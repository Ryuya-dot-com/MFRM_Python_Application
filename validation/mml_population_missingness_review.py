"""Interpret frozen MML comparison/coverage counts without changing endpoints.

The frozen conditional rejection criterion remains conditional. Its
all_planned_rate counts delivered rejections, not unknown test decisions.
This review adds completion bounds; it neither reruns a fit nor supplies a
reference distribution for a sample on which the test was unavailable.
"""
import math

from validation.bias_operating_metrics import binary_operating_summary
from validation.mml_population_analysis import binomial_interval


def review_calibration_coverage(summary, scenario):
    """One fixed calibration target per replicate, not random Person targets.

    A missing interval has no defined coverage event. Delivered coverage and
    hypothetical binary-completion bounds do not assign such an event to it.
    This supplement never changes the registered conditional acceptance rule.
    """
    count = scenario['replicates']
    if (summary.get('status') != 'complete' or summary.get('scenario_id') != scenario['id']
            or type(count) is not int or count < 1 or summary.get('planned') != count):
        raise ValueError('Matching complete calibration evidence required')
    output = []
    for row in summary['calibration']:
        if not row['comparable']:
            continue
        metrics = row['metrics']
        observed = metrics['coverage_given_interval']
        delivered = metrics['coverage_per_known_target']
        availability = metrics['interval_availability']
        denominator = observed['denominator']
        if (type(denominator) not in (int, float) or not math.isfinite(denominator)
                or denominator != int(denominator) or not 0 <= denominator <= count):
            raise ValueError('One binary calibration outcome per replicate required')
        available = int(denominator)
        rate = observed['estimate']
        if available == 0:
            if rate is not None:
                raise ValueError('No available intervals cannot have a conditional coverage rate')
            hits = 0
        else:
            if type(rate) not in (int, float) or not math.isfinite(rate) or not 0 <= rate <= 1:
                raise ValueError('Invalid conditional coverage rate')
            hits = round(rate*available)
            if not math.isclose(rate, hits/available, rel_tol=0., abs_tol=1e-15):
                raise ValueError('Conditional coverage does not represent integer binary counts')
        for record, expected_rate in [(delivered, hits/count), (availability, available/count)]:
            if (record['denominator'] != count or record['planned_replicates'] != count
                    or type(record['estimate']) not in (int, float)
                    or not math.isclose(record['estimate'], expected_rate, rel_tol=0., abs_tol=1e-15)):
                raise ValueError('Calibration coverage/availability denominator differs')
        if (observed['planned_replicates'] != count
                or metrics['truth_unknown_replicates'] != 0 or metrics['all_planned_truth_known'] is not True):
            raise ValueError('Calibration truth is fixed by the scenario, including generation failures')
        intervals = row['monte_carlo_intervals']
        if (intervals['coverage_given_interval'] != binomial_interval(hits, available)
                or intervals['coverage_per_known_target'] != binomial_interval(hits, count)):
            raise ValueError('Stored calibration coverage intervals differ from their binary counts')
        evidence = binary_operating_summary(
            [True]*hits+[False]*(available-hits)+[None]*(count-available), planned=count)
        evidence['interpretation'] = (
            'Conditional coverage describes returned intervals. Completion bounds describe possible binary '
            'completions only; the procedure returned no coverage event for an unavailable interval.')
        output.append(dict(kind=row['kind'], facet=row['facet'], level=row['level'],
            conditional_and_completion_evidence=evidence,
            delivered_covering_intervals_per_planned=dict(estimate=hits/count,
                monte_carlo_interval=intervals['coverage_per_known_target'],
                target='Probability that the procedure returns an interval containing the fixed calibration truth'),
            registered_conditional_acceptance=row['pointwise_acceptance']['coverage'],
            registered_acceptance_unchanged=True, scientific_inference_ready=False))
    return dict(scenario_id=scenario['id'], targets=output, independent_unit='generated dataset',
        interpretation='Fixed calibration targets only; Person truth bins and calibration-fixed posterior intervals have separate descriptive summaries.',
        scientific_inference_ready=False)


def review_comparison(summary, scenario):
    """Review a previously verified complete-scenario summary, not its seal.

    Callers still have to verify the study, raw ledger and numerical reference
    evidence. This function audits count semantics only and cannot qualify an
    estimator or alter the study's registered conditional acceptance decision.
    """
    count=scenario['replicates']
    if (summary.get('status')!='complete' or summary.get('scenario_id')!=scenario['id']
        or type(count) is not int or count<1 or summary.get('planned')!=count):
        raise ValueError('Matching complete scenario evidence required')
    comparison=summary.get('rsm_pcm_comparison')
    paired=scenario.get('paired_rsm_pcm_comparison',False)
    if not paired:
        if comparison is not None:raise ValueError('Unregistered model comparison')
        return dict(status='not_applicable',scientific_inference_ready=False)
    if len(scenario['design']['steps'][0])+1<=2:
        raise ValueError('Binary equal models do not have this LRT')
    if not isinstance(comparison,dict):raise ValueError('Missing registered comparison counts')
    planned=comparison.get('planned');observed=comparison.get('available');missing=comparison.get('unavailable')
    rejected=comparison.get('rejected')
    if (any(type(v) is not int for v in [planned,observed,missing,rejected]) or planned!=count
        or not 0<=rejected<=observed<=count or missing!=count-observed):
        raise ValueError('Inconsistent planned/available/unavailable/rejected counts')
    conditional=rejected/observed if observed else None
    for key,expected in [('conditional_rate',conditional),('all_planned_rate',rejected/count)]:
        actual=comparison.get(key)
        if expected is None:
            if actual is not None:raise ValueError('No observed tests cannot have a conditional rate')
        elif (type(actual) not in (int,float) or not math.isfinite(actual)
              or not math.isclose(actual,expected,rel_tol=0.,abs_tol=1e-15)):
            raise ValueError('Stored rate differs from its count denominator')
    # Only counts are needed for these summaries. This list does not invent
    # replicate identities, response data, or a pairing with other endpoints.
    outcomes=[True]*rejected+[False]*(observed-rejected)+[None]*missing
    reviewed=binary_operating_summary(outcomes,planned=count)
    delivered_interval=binomial_interval(rejected,count)
    if comparison.get('all_planned_mc_interval')!=delivered_interval:
        raise ValueError('Stored delivered-rejection interval differs from its counts')
    return dict(status='comparison_counts_reviewed',scenario_id=scenario['id'],
        registered_conditional_acceptance=comparison['null_acceptance'],
        registered_acceptance_unchanged=True,
        conditional_and_completion_evidence=reviewed,
        delivered_rejections_per_planned=dict(estimate=comparison['all_planned_rate'],
            monte_carlo_interval=delivered_interval,
            target='Probability that this procedure returns an available rejection; missing tests are not assigned a decision'),
        comparison_unavailability=dict(count=missing,rate=missing/count,
            interpretation='Failure of either fitted model can make this comparison unavailable; primary-model fit availability alone does not establish test availability'),
        all_planned_false_positive_or_power_point_identified=missing==0,
        interpretation='A conditional acceptance result describes available tests. Completion bounds describe uncertainty from missing decisions and do not define a test on an unavailable sample. No new acceptance rule or automatic qualification is applied.',
        scientific_inference_ready=False)
