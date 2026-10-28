"""Population-recovery mathematics in native identified coordinates.

No mean alignment, no promotion of finite-rule diagnostics to inference, and
no treatment of parameter rows within one replicate as independent samples.
"""
import numpy as np
from scipy.special import ndtr, ndtri


def calibration_coordinates(config, point):
    """Exact output/Jacobian map, including fixed and dependent coordinates."""
    point = np.asarray(point, dtype=float)
    if (config['model'] not in {'RSM', 'PCM', 'GPCM'} or point.ndim != 1 or not len(point)
        or not np.isfinite(point).all() or config['noncenter_facet'] not in config['facet_names']):
        raise ValueError('Native identified model and finite joint coordinates required')
    width, cursor, rows = len(point), 0, []

    def block(kind, facet, levels, centered):
        nonlocal cursor
        size = len(levels); free = size-int(centered)
        jacobian = np.zeros((size, width)); jacobian[:free, cursor:cursor+free] = np.eye(free)
        if centered: jacobian[-1, cursor:cursor+free] = -1.
        cursor += free
        for level, jac in zip(levels, jacobian, strict=True):
            rows.append(dict(kind=kind, facet=facet, level=str(level), value=float(jac@point), jacobian=jac))
        return jacobian

    for facet in config['facet_names']:
        block(facet, facet, config['facet_levels'][facet], facet != config['noncenter_facet'])
    owners = ['Shared'] if config['model'] == 'RSM' else config['facet_levels'][config['step_facet']]
    for owner in owners:
        block('Step', str(owner), range(1, config['n_cat']), True)
    if config['model'] == 'GPCM':
        levels = config['facet_levels'][config['slope_facet']]
        jacobian = block('LogSlope', config['slope_facet'], levels, True)
        for level, jac in zip(levels, jacobian, strict=True):
            value = float(np.exp(jac@point))
            rows.append(dict(kind='Slope', facet=config['slope_facet'], level=str(level),
                value=value, jacobian=value*jac, log_value=float(jac@point), log_jacobian=jac))
    if cursor != width-1:
        raise ValueError('Joint coordinates differ from the native identified model dimension')
    jac = np.zeros(width); jac[-1] = 1.
    rows.append(dict(kind='LogPopulationSD', facet='Population', level='SD', value=float(point[-1]), jacobian=jac))
    value = float(np.exp(point[-1]))
    rows.append(dict(kind='PopulationSD', facet='Population', level='SD', value=value,
        jacobian=value*jac, log_value=float(point[-1]), log_jacobian=jac))
    if not np.isfinite(point).all() or any(not np.isfinite(row['value']) for row in rows):
        raise ValueError('Nonfinite calibration coordinate')
    return rows


def candidate_wald_intervals(coordinates, covariance):
    """Candidate full-joint Wald intervals; positive quantities use log CIs.

    This is a research calculation, not an inference-readiness decision. The
    caller must retain numerical eligibility, missing intervals and failures.
    Fixed coordinates have no sampling interval or recovery-success credit.
    """
    cov = np.asarray(covariance, dtype=float)
    width = len(coordinates[0]['jacobian'])
    if (cov.shape != (width, width) or not np.isfinite(cov).all()
        or np.max(abs(cov-cov.T)) > 1e-10 or np.linalg.eigvalsh(cov).min() <= 0):
        raise ValueError('A positive-definite full joint covariance is required')
    z = float(ndtri(.975)); output = []
    for row in coordinates:
        jac = row['jacobian']; variable = bool(np.any(jac != 0))
        value = row['value']; se = float(np.sqrt(jac@cov@jac)) if variable else None
        lower = upper = None
        if variable and 'log_value' in row:
            log_jac = row['log_jacobian']; log_se = np.sqrt(log_jac@cov@log_jac)
            with np.errstate(over='ignore', under='ignore'):
                lower, upper = np.exp([row['log_value']-z*log_se, row['log_value']+z*log_se])
        elif variable:
            lower, upper = value-z*se, value+z*se
        if variable and (not np.isfinite([se, lower, upper]).all() or ('log_value' in row and lower <= 0)):
            raise ValueError('Candidate interval is outside finite representable coordinates')
        basis = 'candidate_joint_information_log_wald' if 'log_value' in row else 'candidate_joint_information_wald'
        output.append(dict(row, se=se, lower=lower, upper=upper, recovery_comparable=variable,
            interval_basis=basis if variable else 'fixed_by_identification',
            scientific_inference_ready=False))
    return output


def population_truth_quartile(theta, population_sd):
    """Fixed population-normal cutpoints; no data-dependent ranking or bins."""
    theta = np.asarray(theta, dtype=float)
    if (not np.isfinite(theta).all() or isinstance(population_sd, (bool, np.bool_))
        or not np.isfinite(population_sd) or population_sd <= 0):
        raise ValueError('Finite Person truth and positive population SD required')
    return np.searchsorted([.25, .5, .75], ndtr(theta/population_sd), side='right')+1


def replicate_ratio(numerators, denominators):
    """Ratio-of-totals with a replicate-cluster delta-method MCSE.

    Arrays have one entry per PLANNED independent replicate. A replicate with
    no eligible rows contributes (0,0), rather than disappearing. For a random
    truth-bin size this estimates the per-Person target, distinct from the
    equal-replicate mean already used by the fixed-Person UI simulation.
    """
    sums, counts = np.asarray(numerators, dtype=float), np.asarray(denominators, dtype=float)
    if (sums.ndim != 1 or not len(sums) or counts.shape != sums.shape or not np.isfinite(sums).all()
        or not np.isfinite(counts).all() or np.any(counts < 0) or np.any((counts == 0) & (sums != 0))):
        raise ValueError('Finite per-replicate sums and nonnegative denominators required')
    eligible = int(np.count_nonzero(counts)); total = float(counts.sum())
    estimate = float(sums.sum()/total) if total > 0 else None
    mcse = None
    if len(sums) > 1 and eligible > 1:
        influence = sums-estimate*counts
        mcse = float(np.std(influence, ddof=1)/(np.sqrt(len(sums))*np.mean(counts)))
    return dict(estimate=estimate, mcse=mcse, planned_replicates=len(sums),
        nonempty_replicates=eligible, empty_replicates=len(sums)-eligible, denominator=total,
        mcse_unit='independent simulation replicate', estimand='ratio_of_expected_sums_to_expected_counts')


def summarize_replicate_errors(errors, coverage):
    """Arrays retain every generated target, with NaN for unavailable output.

    There must be one slot per planned replicate. None means truth itself is
    unavailable (e.g. generation failure); an empty array means a known empty
    truth bin. A failed fit has its known target count and all-NaN errors/CIs.
    Coverage uses 0/1/NaN, and failures remain in the known-target denominator.
    """
    if len(errors) != len(coverage) or not len(errors):
        raise ValueError('Every planned replicate requires error and coverage slots')
    sums, squares, valid, hits, intervals, targets = ([] for _ in range(6))
    unknown = 0
    for error, covered in zip(errors, coverage, strict=True):
        if error is None or covered is None:
            if error is not None or covered is not None: raise ValueError('Unknown truth must use two None slots')
            unknown += 1; error = covered = np.array([], dtype=float)
        error, covered = np.asarray(error, dtype=float), np.asarray(covered, dtype=float)
        if (error.ndim != 1 or covered.shape != error.shape or np.isinf(error).any()
            or not np.isin(covered[np.isfinite(covered)], [0., 1.]).all() or np.isinf(covered).any()
            or np.any(np.isfinite(covered) & ~np.isfinite(error))):
            raise ValueError('Coverage must be 0/1/NaN and correspond to a returned estimate')
        ok = np.isfinite(error); sums.append(float(error[ok].sum())); squares.append(float((error[ok]**2).sum()))
        valid.append(int(ok.sum())); hits.append(float(np.nansum(covered)))
        intervals.append(int(np.isfinite(covered).sum())); targets.append(len(error))
    bias = replicate_ratio(sums, valid); mse = replicate_ratio(squares, valid)
    rmse = None if mse['estimate'] is None else float(np.sqrt(mse['estimate']))
    return dict(bias=bias, mse=mse, rmse=rmse,
        mcse_rmse=mse['mcse']/(2*rmse) if rmse is not None and rmse > 0 and mse['mcse'] is not None else None,
        coverage_given_interval=replicate_ratio(hits, intervals),
        coverage_per_known_target=replicate_ratio(hits, targets),
        estimate_availability=replicate_ratio(valid, targets), interval_availability=replicate_ratio(intervals, targets),
        truth_unknown_replicates=unknown, all_planned_truth_known=unknown == 0, scientific_inference_ready=False)


def calibration_se_summary(errors, standard_errors):
    """One fixed calibration target per replicate; joint-moment MCSE.

    Compare sqrt(mean(model SE squared)) with empirical sampling SD on the
    same eligible repetitions. The delta-method MCSE includes their covariance;
    it does not assume independent model variances and estimation errors, or
    a normal distribution of estimates. Availability is reported separately.
    """
    error, se = np.asarray(errors, dtype=float), np.asarray(standard_errors, dtype=float)
    if error.ndim != 1 or se.shape != error.shape or np.isinf(error).any() or np.isinf(se).any():
        raise ValueError('One finite-or-missing estimate/SE per planned replicate required')
    valid = np.isfinite(error) & np.isfinite(se) & (se > 0)
    n = int(valid.sum()); result = dict(planned_replicates=len(error), eligible_replicates=n,
        empirical_sd=None, rms_model_se=None, ratio=None, mcse_ratio=None,
        mcse_empirical_sd=None, mcse_rms_model_se=None, scientific_inference_ready=False)
    if n < 2: return result
    error, se = error[valid], se[valid]
    # Center first to avoid cancellation of E(error^2)-E(error)^2 when bias
    # is large; the empirical-variance influence is translation-invariant.
    centered = error-error.mean(); correction = n/(n-1)
    moments = np.c_[centered**2, se**2]
    empirical_variance = float(moments[:,0].mean()*correction)
    model_variance = float(moments[:,1].mean())
    if empirical_variance <= 0: return result
    empirical = np.sqrt(empirical_variance); model = np.sqrt(model_variance)
    ratio = model/empirical
    covariance = np.cov(moments, rowvar=False, ddof=1)/n
    gradient = np.array([-.5*ratio*correction/empirical_variance, .5*ratio/model_variance])
    variance = float(gradient@covariance@gradient)
    # A theoretically zero delta variance can acquire signed roundoff. Report
    # unavailable if materially negative, and bound only machine-level noise.
    allowance = 32*np.finfo(float).eps*float(abs(gradient)@abs(covariance)@abs(gradient))
    if variance < -allowance: raise FloatingPointError('Invalid joint-moment MC variance')
    result.update(empirical_sd=float(empirical), rms_model_se=float(model), ratio=float(ratio),
        mcse_ratio=float(np.sqrt(max(0.,variance))),
        mcse_empirical_sd=float(correction*np.std(moments[:,0],ddof=1)/(2*np.sqrt(n)*empirical)),
        mcse_rms_model_se=float(np.std(moments[:,1],ddof=1)/(2*np.sqrt(n)*model)))
    return result
