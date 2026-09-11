# CMLE-WLE bootstrap design

## Status and purpose

This is a prospective repository-only design. No bootstrap result, confidence
interval, total standard error, Streamlit control, or public capability claim
is authorized by this document. The purpose is to separate two different
questions that the current asymptotic coefficient-draw sensitivity cannot
answer simultaneously:

1. How does response-pattern variation at each observed Person total move the
   exact-CMLE calibration and the dependent fit-sample WLE scores?
2. How does joint response variation under plug-in CMLE + WLE parameters move
   the full two-stage estimator, including Person totals?

Both lanes preserve the observed administered-row and missingness design.
Neither silently changes estimator, imputes a missing response, adds an unseen
unit, or repairs a rank-deficient replicate.

## Lane A: fixed-score conditional-pattern bootstrap

For Person `p`, let `s_p` be the retained observed total and `q_uk(beta_hat)`
the fitted exact-CMLE category log kernel. Sample a response vector from

```text
P(X_p = x | sum(x) = s_p, beta_hat)
    proportional to exp(sum_u q_u,x_u(beta_hat)).
```

A forward/backward coefficient recurrence must sample categories sequentially
without enumerating all response vectors. At row `u`, the probability of
category `k` is proportional to the current row kernel times the suffix
coefficient for remaining score `s_remaining - k`. Impossible categories have
exact zero probability. The sampled total must equal `s_p` exactly before the
replicate can be fitted.

Each valid replicate refits exact CMLE and then applies the repository Warm WLE
bridge. This lane preserves Person totals, so it isolates pattern-driven
calibration variation and its same-sample coupling to WLE. It is not a total
Person-score sampling distribution. Exact all-minimum and all-maximum totals
have only one admissible response vector and are negative controls, not
evidence of zero total uncertainty.

## Lane B: joint plug-in parametric bootstrap

For each retained fit-sample Person and administered row, form

```text
P(X_pu = k | theta_hat_p, beta_hat)
    = softmax(q_uk(beta_hat) + k theta_hat_p).
```

Draw categories independently conditional on the fitted Person and structural
parameters, refit exact CMLE, and re-score all retained Persons with Warm WLE.
Person totals may change, and a baseline exact-extreme Person need not remain
extreme. This lane carries response-score variation and calibration/WLE
dependence under the fitted plug-in model. Its percentile and standard-
deviation summaries remain candidate sensitivity outputs until repeated-truth
simulation establishes interval coverage and failure behavior. Plug-in WLE
values are estimates, not known generating truths.

## Required outputs and non-equivalences

Every attempted replicate receives a row with its seed, lane, model, CMLE
convergence/readiness, rank/nullity, WLE availability, failure stage, and
reason. Failed or structurally unidentified replicates remain in the
denominator and are never counted as negative decisions.

Person output must retain the baseline estimate and exact-extreme status,
successful/attempted replicate counts, bootstrap mean, SD, quantiles, sign and
fit-zone transitions, raw/display threshold mismatches, and extreme-state
transitions. Model output must report numerical failure and inference-readiness
rates with Monte Carlo uncertainty.

The following quantities are not interchangeable:

- asymptotic coefficient-draw SD: local normal calibration sensitivity with
  same-sample dependence omitted;
- fixed-score bootstrap SD: conditional pattern/calibration coupling at the
  observed Person total;
- joint plug-in bootstrap SD: model-based joint two-stage sampling sensitivity;
- conditional WLE SE: Person information with calibration fixed.

In particular, subtracting one bootstrap variance from another must not be
labelled a variance decomposition without a proved covariance identity. No
bootstrap quantile is a confidence interval merely because it is a percentile.

## Implementation and pilot gates

Before inspecting a 200-replicate implementation pilot for each model/lane:

1. enumerate small RSM and PCM conditional response-pattern probabilities and
   match the dynamic recurrence to at most `1e-12`;
2. require exact Person-total preservation for every fixed-score draw;
3. match joint-lane empirical category frequencies to the analytic kernels
   within a prospectively computed Monte Carlo tolerance;
4. reproduce retained outputs exactly with the same seed and differ with a
   different seed;
5. reject unsupported model, non-inference-ready calibration, parameter-
   identity mismatch, impossible score, material covariance/design mismatch,
   excessive work, or non-finite kernel without fallback;
6. retain every refit failure and distinguish optimizer convergence from
   inference readiness; and
7. make all decisions from unrounded values while exporting display-boundary
   mismatches.

The 200-replicate matrix is an implementation pilot, not precision-qualified
bootstrap evidence. It must not set a confirmatory success threshold after its
results are seen. A later repeated-truth ADEMP study will freeze sample-size,
facet, category, missingness, sparse/connectivity, anchor, extreme-score, and
bias-boundary conditions; it will assess bias, RMSE, interval coverage,
failure/readiness, threshold decisions, and Monte Carlo uncertainty.

## Resource and UX contract

The core must require an explicit seed and apply independent maximum-replicate,
row-replicate-work, coefficient-state, memory, and elapsed-time guards. Partial
progress must be archivable. Cancellation must stop safely without converting
partial output into a passed result.

A later one-click UI may orchestrate the two lanes, but the first-read panel
must show separate cards for CMLE refit readiness, fixed-score sensitivity,
joint plug-in sensitivity, exact-extreme transitions, interval-coverage status,
and public-use status. “Bootstrap completed” cannot by itself create a green
inference card.

## Method anchors

- [`immer_cml` documentation](https://rdrr.io/cran/immer/man/immer_cml.html)
  describes linear PCM conditional calibration, coefficient covariance,
  sufficient statistics, score frequencies, and retained Person scores.
- [Alexandrowicz and Draxler (2016)](https://link.springer.com/article/10.1186/s40488-016-0039-y)
  describes fixed-marginal conditional bootstrap sampling for the binary Rasch
  model and explicitly distinguishes it from plug-in response generation. The
  RSM/PCM recurrence in this design is an extension that must pass the stated
  enumeration gate rather than being presumed equivalent.
- [Warm (1989)](https://doi.org/10.1007/BF02294627) defines weighted likelihood
  Person estimation; it does not turn a conditional calibration bootstrap into
  a total uncertainty result by itself.
