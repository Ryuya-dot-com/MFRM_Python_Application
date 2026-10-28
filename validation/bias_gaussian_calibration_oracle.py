"""Exact Gaussian control for fixed, reused and independent calibration.

This is a separate analytic illustration, not an MFRM performance estimate or
a replicate in either frozen population study. It generates no response data
and runs no optimizer. Rational row coefficients independently check the
closed-form variances and local effect targets.
"""
from fractions import Fraction as F
import json
import math
from statistics import NormalDist


def check(n1, n0, calibration_n):
    """Y_i = location + effect * group_i + N(0, 1), independent rows."""
    n = n1 + n0
    # Order: group 1, group 0, external calibration at the baseline location.
    row_coefficients = {
        'known_location': [F(1, n1)] * n1 + [F(0)] * (n0 + calibration_n),
        'same_data_null_location': [F(1, n1) - F(1, n)] * n1
            + [-F(1, n)] * n0 + [F(0)] * calibration_n,
        'independent_estimated_location': [F(1, n1)] * n1 + [F(0)] * n0
            + [-F(1, calibration_n)] * calibration_n,
        'joint_profile_effect': [F(1, n1)] * n1 + [-F(1, n0)] * n0
            + [F(0)] * calibration_n,
    }
    # These formulas are checked against the independent row representation.
    targets = dict(known_location=F(1), same_data_null_location=F(n0, n),
        independent_estimated_location=F(1), joint_profile_effect=F(1))
    variances = dict(known_location=F(1, n1),
        same_data_null_location=F(1, n1)-F(1, n),
        independent_estimated_location=F(1, n1)+F(1, calibration_n),
        joint_profile_effect=F(1, n1)+F(1, n0))
    z = NormalDist().inv_cdf(.975)
    rows = []
    for basis, coefficients in row_coefficients.items():
        target = sum(coefficients[:n1]); variance = sum(c*c for c in coefficients)
        assert target == targets[basis] and variance == variances[basis]
        # The known-location estimate also subtracts the true location.
        assert sum(coefficients) == (1 if basis == 'known_location' else 0)
        information = F(n1*n0, n) if basis == 'joint_profile_effect' else F(n1)
        scale = variance * information
        null_rejection = math.erfc(z / math.sqrt(2 * float(scale)))
        rows.append(dict(basis=basis, effect_multiplier=str(target),
            sampling_variance=str(variance), likelihood_curvature_variance=str(1/information),
            lr_chi_square_scale=str(scale),
            null_rejection_at_nominal_005=null_rejection,
            null_interval_coverage=1-null_rejection))
    assert rows[1]['null_rejection_at_nominal_005'] < .05
    assert rows[2]['null_rejection_at_nominal_005'] > .05
    for index in (0, 3):
        assert math.isclose(rows[index]['null_rejection_at_nominal_005'], .05, abs_tol=1e-14)
    return dict(n1=n1, n0=n0, calibration_n=calibration_n, rows=rows)


if __name__ == '__main__':
    print(json.dumps(dict(classification='EXACT_GAUSSIAN_ANALYTIC_CONTROL',
        scientific_inference_ready=False, independent_simulation_replicates=0,
        assumptions='Independent Gaussian errors, known unit variance, fixed group sizes; no latent Person/EAP, boundaries, weights or MFRM fit',
        interpretation='The reused-location residual targets an attenuated effect. Its fixed-curvature LR is not the jointly profiled LR. Neither correction is asserted for the application.',
        cases=[check(*sizes) for sizes in [(8,8,8),(8,24,8),(24,8,24),(1,1,1)]]), indent=2))
