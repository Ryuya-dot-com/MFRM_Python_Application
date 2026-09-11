# Known-assignment confirmation: adversarial review record

Date: 2026-08-11  
Status: incorporated before confirmatory Person, response, assignment, or fit generation

## Bottom line

The preceding evidence chain is operationally strong, but its main preflight
statistic must not be called “bias.” It is the expected change in a
within-dataset RMSE across three assignment conditions. The next confirmation
therefore freezes that estimand explicitly and adds a separate Monte Carlo
bias/variance/MSE decomposition.

The independent replicate is one fresh Person-vector triplet containing
`gamma=-0.8, 0.0, +0.8`. The 600 datasets produced by 200 triplets are not 600
independent observations for the primary test.

## Binding corrections

1. **Estimand language.** The primary statistic is
   `0.5 * (RMSE[-0.8] + RMSE[+0.8]) - RMSE[0]`, with each RMSE calculated over
   the same four fixed Rater truth levels. It is not `E(error)`,
   `sqrt(E(error^2))`, or an estimator ranking.
2. **Sample size.** The fixed `R=200` is justified by Monte Carlo precision,
   not the favorable preflight mean. Using the n=10 screening SD and its 95%
   variance upper bound gives an expected 95% half-width of about `0.0143`,
   below the frozen `0.015` target. Achieved precision is a separate gate and
   cannot trigger more runs.
3. **R verification.** The prior R 4.5.1 work independently verified endpoint
   aggregation only. The confirmation requires free-SD MML likelihood and
   cross-fit checks on 36 preselected datasets, including cross-evaluated
   likelihoods, gradients, constraints, and Q31/Q61 sensitivity.
4. **FACETS separation.** FACETS/Python JMLE parity uses a fixed 36-dataset
   subset and a separate denominator. A GUI, report, or path failure does not
   delete valid MML evidence. A parity discrepancy that implicates shared data
   encoding or generation quarantines the affected scientific artifacts.
5. **Failure accounting.** The primary requires 200/200 complete triplets.
   There are no replacement seeds, success-until-fit loops, or complete-case
   confirmatory claims.
6. **FACETS precision.** Umean=6 measures and thresholds may support numerical
   calibration; ordinary fit fields remain two-decimal displays. Fit values
   are never treated as raw inputs. Any display-level comparison is an
   individual quantization-interval check, not an average-error check.

## Deliberate compromise

The n=4 fresh-vector dispersion estimate is too uncertain to guarantee the
`0.015` precision target at `R=200`; its variance upper bound would favor about
400 triplets. The selected 200-triplet design is an explicit compute/precision
tradeoff. It reports `DirectionConfirmed` and `PrecisionQualified` separately
and prohibits extension after seeing either result.

## Generalization boundary

Within-vector standardization means the realized Persons are not an iid sample
from `N(0, 0.8^2)`. The resulting claim is deliberately narrow: one finite
80-Person design, four fixed Raters, one fixed-margin connected assignment
family, three gamma values, and one heterogeneous PCM configuration. Broader
operating-characteristic work belongs in later, separately registered cells.

The complete binding protocol is
`validation/known_assignment_confirmatory_plan_20260811.json`.
