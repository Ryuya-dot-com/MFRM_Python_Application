# Native Exact CMLE Structural-Calibration Design

- Status: repository-only numerical core; not exposed in the Streamlit estimator selector
- Date: 2026-08-09
- Native hard-anchor increment: 2026-08-10
- Implementation: `mfrm_app/cmle.py`
- Validation: `tests/test_cmle.py` and `validation/cmle_phase0_20260809/RESULTS.md`

## Decision

The application will investigate native exact conditional maximum-likelihood
estimation as a differentiated strict-Rasch calibration route. The defensible
claim is not that conditional maximum likelihood is new. The candidate value
is a Python-native, long-format, arbitrary-facet workflow that refuses
conditionally unidentified designs before optimization and retains category,
missingness, boundary, and estimator-scope evidence in the result.

This research core is deliberately below the public product boundary. A numerically
successful fit does not by itself authorize a new UI option, public parity
claim, or reporting recommendation.

## Estimand and likelihood

For response unit `u`, category `k`, Person `p`, and structural free-coordinate
vector `beta`, the supported Rasch-family kernel is

```text
log P(Y_pu = k | theta_p, beta) =
    k theta_p + q_uk(beta) - log sum_h exp(h theta_p + q_uh(beta))

q_uk(beta) = a_uk + X_uk beta
```

`X_uk` contains signed facet main-effect contrasts and cumulative threshold
contrasts. `a_uk` is zero for the legacy sum-to-zero design and contains the
fixed category-kernel contribution of any native hard facet anchor.
Conditioning on the Person raw score `S_p = sum_u Y_pu` removes `theta_p`:

```text
log L_C(beta) =
    sum_p [sum_u q_u,y_pu(beta) - log Z_pattern(p),S_p(beta)]
```

`Z_pattern,S` is the coefficient for total score `S` in the product of the
response-unit category polynomials. A scaled dynamic program carries that
coefficient and the unnormalized first moment of the structural sufficient
statistic for every attainable score. Its gradient is the observed statistic
minus its exact conditional expectation.

For rank and uncertainty, the same recurrence also carries the unnormalized
second moment. For sufficient statistic `T`, this gives

```text
H(beta) = sum_pattern,sum_score n_pattern,score
          * Cov_beta(T | pattern, raw score)
```

which is the exact Hessian of the negative conditional log likelihood for the
linear Rasch-family kernel. Pre-fit rank, fitted rank, positive-definiteness,
and covariance therefore use one common analytical information object; none
is obtained by finite-differencing the gradient. Per-unit log-kernel shifts
and per-stage coefficient scaling protect the recurrence against ordinary
underflow/overflow. The preflight work proxy reflects the quadratic
parameter cost of the second-moment states and remains fail-closed above its
configured cap.

Persons with all-minimum or all-maximum scores have a single response pattern
conditional on their score and therefore contribute zero structural
information. They remain in the typed Person-status table rather than being
silently re-expressed as finite ability estimates.

## Supported structural-calibration contract

| Area | Structural-calibration contract |
| --- | --- |
| Response families | RSM and rectangular PCM, including their binary reductions |
| Discrimination | Fixed unit discrimination only |
| Structural effects | Additive fixed facet main effects and RSM/step-facet PCM thresholds |
| Identification | Exact sum-to-zero coordinates for unanchored facet/threshold blocks; native hard-anchored facets use an affine fixed-offset map and direct free coordinates for unanchored levels |
| Input | Long-format observations with explicit `rating_min` and `rating_max` |
| Missingness | Planned incomplete response-unit sets, grouped by common design pattern |
| Weights | Unit row weights only |
| Person parameters | Conditioned out; no research-core Person estimate |
| Uncertainty | Exact observed conditional information from conditional sufficient-statistic covariance |
| Information criteria | `ConditionalAIC` only, labelled as comparable only across matched conditional likelihoods |
| Resource guards | 200 free parameters, 50 million work-proxy units, and 512 MiB conservative core peak-memory proxy by default |

The explicit rating support is a substantive contract. A top or intermediate
category with zero observations is not dropped merely because it is absent
from one virtual response unit's observed values.

## Fail-closed exclusions

The structural-calibration core blocks rather than approximates:

- bounded or unrestricted GPCM with estimated slopes;
- latent regression and population-distribution terms;
- facet regularization or penalized conditional likelihood;
- non-unit observation weights;
- duplicated Person × response-unit cells without an explicit Event/Occasion
  identity or an explicit research override;
- malformed, duplicated, unknown, non-finite, Person, or Step hard-anchor rows;
- no informative Persons; and
- deficient conditional-information rank.

The conditional-rank check is stricter than an ordinary global facet-network
check. A facet that is constant across every response unit within each Person
can cancel after conditioning even when its levels appear in the dataset. Such
a parameter is blocked with `conditional_information_rank_deficient`.

## Numerical readiness

The result separates the optimizer's raw success flag from numerical
stationarity. SciPy BFGS can report precision loss after reaching a small
analytical gradient. The result therefore records both `OptimizerSuccess` and
`Converged`; inferential readiness additionally requires:

- a finite objective;
- gradient sup norm within the recorded stationarity tolerance;
- full fitted conditional-information rank; and
- a finite inverse information matrix.

`InformationMethod` and the audit field `rank_audit_method` identify this path
as `exact_conditional_second_moment_dp`. BFGS is still the optimizer, but its
internal inverse-Hessian approximation is not used for reported standard
errors or readiness. If BFGS stops above the fixed stationarity tolerance, a
full-rank positive-definite fit may receive up to eight exact-Newton polishing
steps. Every step requires Armijo decrease, or—when the predicted change is
below objective resolution—an objective difference within a documented
floating-point allowance plus material gradient reduction. Raw BFGS success,
its gradient, polishing steps, and the final gradient remain separately
recorded.

## Scaling evidence and limit

The deterministic scaling matrix in
`validation/cmle_scaling_20260809/RESULTS.md` covers 100--5,000 Persons,
8--64 complete response units, 6--30 free parameters, and 1--100 missingness
patterns. All 14 eligible cases were inference-ready; separate 100- and
160-unit PCM cases demonstrated the default work and memory blocks.

The memory proxy is deliberately conservative for the large numerical arrays,
but it is neither a process-RSS estimate nor a guaranteed upper bound on Python,
pandas, or Streamlit overhead. This is enough to retain the exact engine as a
bounded research core, not to
remove the public scaling gate. The next scaling evidence must include
multi-platform and in-Streamlit memory measurement, adversarial extreme
parameters, concurrent deployment assumptions, and cancellation/interruption
behavior. A user-visible override must never be automatic.

This does not establish model fit, local independence, missing-at-random
behavior, or validity of a proposed score use.

## Repository-only CMLE plus WLE bridge

The first explicit post-calibration Person-scoring bridge is now implemented
in `mfrm_app/cmle_person_scoring.py`. It is a two-stage estimator, not a new
form of conditional likelihood:

```text
exact CMLE structural calibration
    -> reconstruct q_uk from the fitted CMLE row design
    -> fixed-calibration Warm WLE Person score
```

For each retained fitted row and category, the bridge computes
`q_uk = row_offset[u,k] + row_design[u,k,:] beta_hat`; the Person coefficient
is the internal category `k` under the supported unit-discrimination RSM/PCM kernel. It then
solves the separately validated Warm adjusted score. The combined label is
`CMLE_exact_structural_plus_WLE_fixed_calibration_person`, and exact-extreme
Persons retain both their CMLE conditional status and their finite WLE status.
They are never relabelled as finite JMLE MLEs.

The prospectively registered deterministic bridge evidence passed for RSM and
PCM: 14 Persons including two exact extremes, zero bridge-versus-generic WLE
theta/SE difference, zero CMLE category-surface/intercept difference, and a
maximum adjusted-score residual of `6.18e-16`. Coefficient row permutation did
not affect scores, while a missing coefficient failed closed. This establishes
the computational handoff only. WLE SE treats `beta_hat` as fixed and therefore
does not include CMLE calibration uncertainty; fit-sample scoring is not
new-Person prediction. See
`validation/cmle_wle_bridge_20260810/CMLE_WLE_BRIDGE_RESULTS.md`.

## Native hard facet anchors

The hard-anchor research increment accepts a fail-closed `hard_anchors` table with
`ParameterType`, `Facet`, `Level`, and `Value`. Version 1 supports finite fixed
values for observed non-Person facet levels only; Step/group/soft anchors and
automatic anchor selection remain unsupported. Anchor values are on the
expanded structural-estimate scale before the facet sign is applied. For an
anchored facet,

```text
expanded facet effect = fixed_offset + J beta
category kernel       = row_offset + row_design beta
```

Anchored rows have a zero Jacobian, the supplied estimate exactly, SE zero,
`Anchored=True`, and `Constraint=hard_anchor`. Unanchored levels in the same
facet are direct free coordinates; the old sum-to-zero constraint is not also
imposed. Facets without anchors and all threshold blocks keep their original
sum-to-zero coordinates. The conditional-rank audit is performed after this
reduction and uses the fixed offsets in the conditional probabilities.

The prospectively specified two-replicate engineering smoke reused 24,000
retained response rows across four dense/sparse clean/random conditions and
four anchor scenarios. All 32 fits returned inference-ready; input identity,
free-parameter counts, exact fixed values/zero SEs, and raw fit recomputation
passed. This is implementation evidence, not an anchor-share or robustness
claim.

A transparently post-result diagnostic confirmed a critical limitation. If
every Rater anchor is contaminated by the same constant `delta`, that common
facet-origin shift cancels from the conditional likelihood at fixed Person
total. With `delta=0.25`, the maximum conditional-log-likelihood difference
was `1.82e-12`; WLE moved by `+0.25` within `3.61e-08`, while the maximum
Infit/Outfit difference was `7.37e-08` and no raw either-upper decision
changed. Thus Person fit cannot validate the absolute anchor origin. External
anchor provenance and scale validation are mandatory even when fit is
unchanged.

The next prospectively registered stress reused 120,000 retained responses in
40 RunIds and completed 320/320 inference-ready fits: four dense/sparse,
clean/random-response conditions crossed with zero, one, two, three, or five
correct anchors and three common/differential/focal contamination scenarios.
This larger debugging run still does not identify a universal anchor share.
Correct nested anchor sets did not monotonically improve Person recovery, and
their share was confounded with which deterministic Rater levels were fixed.
The one-anchor coordinate system was also less well-conditioned than the
unanchored system in these data (mean information condition number `41.45`
versus `9.68`); two, three, and five anchors reduced it to `20.15`, `12.40`,
and `9.35`. All remained full rank, illustrating why full rank alone is not an
anchor-validity certificate.

Relative anchor errors mattered even when their mean effect looked small. The
zero-mean differential pattern moved condition/group mean WLE by only
`0.00497`--`0.01307` logits relative to three correct anchors, but the maximum
absolute Person movement was `0.287525` and raw either-upper decisions changed
for `2.4%`--`15.8%` of paired rows. Wrong fixed values occasionally had higher
in-sample conditional likelihood in sparse samples, so likelihood cannot act
as an automatic anchor oracle either. Across 64,000 Person-scenario rows, a
counterfactual 3-decimal MnSq decision changed 12 raw decisions; 6 decimals
changed none in this dataset, without establishing a universal safe display
precision. All decisions and exports must therefore retain raw values and
boundary distance.

See `validation/cmle_native_hard_anchor_smoke_20260810/` and
`validation/cmle_native_hard_anchor_common_shift_20260810/`, plus the registered
plan, evidence report, and critical review under
`validation/cmle_native_anchor_differential_stress_20260810/`.

The following prospectively registered connectivity stress crossed seven
fully disconnected, two-component, chain, and hub designs with five randomized
correct-anchor counts over ten replicates. A transparent rule based on the
realized graph of non-extreme within-Person Rater contrasts agreed with exact
zero-point conditional eligibility in 350/350 cases. With no bridge, fixing
zero through three Raters left nullity; fixing all six removed the Rater block
rather than learning it. In the strong two-component design, two random
anchors were eligible in exactly the 7/10 replicates where the selected levels
covered both components.

Prefit rank did not guarantee final inference. All 250 connected cases passed
prefit, but the one-Person-per-edge minimal chain produced one non-finite
conditional-normalizer failure and 12 returned fits with deficient final
information rank. The other 200 weak/strong chain/hub fits were inference-
ready. Condition numbers depended strongly on bridge support, topology,
anchor content, and coordinates; a one-anchor direct-free map was often less
well conditioned than the unanchored sum-to-zero map. Because condition number
is not parameterization invariant, this is an operational optimizer warning,
not an intrinsic information ranking. See
`validation/cmle_native_anchor_connectivity_20260810/`.

The next prospectively registered numerical remediation changed only the
primary optimizer adapter. Non-finite BFGS trial points are now counted and
returned to the line search as rejected positive-infinity trials; the public
exact objective remains fail-closed, and no likelihood parameter is clipped.
The adapter retains the best finite point and forces raw optimizer success
false if that fallback is needed. The fallback branch has a deliberate unit
fixture; it was not used in the retained replay.

On identical connectivity bytes, candidate returnability rose from 286/287 to
287/287 while inference readiness remained 274/287. The target recorded 12
invalid trials and returned nonconverged with rank 8/10 and nullity 2. All 286
baseline-returned ready states and all 82,200 shared Person raw/rounded fit
flags were identical; numerical differences remained at binary64 noise. This
does not diagnose finite-MLE existence or conditional separation. A post-
remediation binary perfect-separation check can still stop at large finite
effects with a small gradient and full finite-point rank. Finite-MLE/existence
evidence is therefore a separate required gate. See
`validation/cmle_finite_domain_optimizer_remediation_20260810/`.

That separate research gate is now implemented in
`mfrm_app/cmle_existence.py`. The first prospectively registered route exactly
enumerated fixed-score support points. It matched all 93 registered controls,
but a deliberately low research cap retained 227/350 connectivity cases as
unavailable. A second prospective route uses fixed-score dynamic programming
as an exact support-function oracle and generates only LP constraints violated
by a candidate direction. It matched 1,184/1,184 exhaustive support maxima,
93/93 fixture statuses, 123/123 completed retained statuses, and 6/6 high-cap
sentinels. The same-byte retained replay resolved 63 structural, 13 boundary,
and 274 interior cases with no unavailable or tolerance-unstable result; all
274 current inference-ready cases were interior.

For 287 eligible cases, theoretical response configurations had median 97,040
and maximum 247,668, while generated constraints had median 163 and maximum
304. Eligible-case P95/max time was 0.521/0.568 seconds and the registered
350-case total was 109.216 seconds. The oracle is exact in discrete support
search but the LP remains binary64, so cross-tolerance disagreement must stay
fail-closed. Anchor-induced removal of a separating direction is conditional
on the fixed anchor values and cannot validate those anchors. See
`validation/cmle_finite_mle_oracle_20260810/`.

The next prospective integration now makes that oracle result a default
`fit_cmle` readiness prerequisite. The same-byte replay returned 287/287 fits:
274 interior cases were ready and 13 boundary cases were non-ready with a
dedicated reason. Oracle status and integrated readiness matched 287/287,
102 selected tests passed, and 82,200 Person raw/rounded fit flags were
unchanged. Maximum numerical movement across structural estimates, likelihood,
WLE, and fit was `2.27e-13`. The current research core still completes the
optimizer on a boundary to preserve technical comparison values; public
orchestration should instead return a structured boundary result before
optimization. See
`validation/cmle_finite_mle_readiness_integration_20260810/`.

The structured early-stop layer is now implemented separately in
`mfrm_app/cmle_workflow.py`. Its prospectively registered replay matched all
93 fixture and 350 retained terminal states, stopped 63 structural and 13
boundary cases before optimization, and optimized only 274 interior cases.
No early-stopped case invoked the optimizer. P95 time was 0.068 seconds for a
design block, 0.306 for a boundary stop, and 1.084 for the ready interior path;
ready estimates/SE stayed within `4.44e-16` and conditional likelihood within
`2.27e-13` of the direct integrated path. Workflow v1 records downstream
Person scoring as available/not-run. Bilingual text is machine-checked but has
not undergone comprehension testing. See
`validation/cmle_structured_workflow_20260810/`.

## Calibration-sensitivity boundary

The first prospectively frozen uncertainty increment is implemented in
`mfrm_app/cmle_wle_uncertainty.py`. It draws exact-CMLE free coordinates from
the fitted asymptotic covariance and re-scores the same fit-sample Persons. It
is deliberately labelled a calibration-draw sensitivity diagnostic, not a
total uncertainty estimator:

```text
beta^(b) ~ Normal(beta_hat, scale^2 V_hat)
    -> reconstruct q_uk^(b)
    -> fixed-calibration Warm WLE theta_p^(b)
    -> report dispersion as sensitivity only
```

The frozen RSM/PCM fixture used 2,000 coefficient draws per model. All 32,000
draw-Person scores returned; model-median calibration-draw SD was
`0.0616`--`0.0801` logits, the maximum Person value was `0.518` logits, and the
maximum ratio to the fixed-calibration conditional WLE SE was `0.427`. Exact
extremes had the largest sensitivities. Zero covariance scale returned
numerical zero, same-seed runs were exactly reproducible, and covariance-scale
response was strictly increasing for both models. Material covariance
non-symmetry or non-positive-semidefiniteness fails closed.

The coefficient-draw quantiles are neither confidence nor credible intervals.
`QuadratureSensitivitySE` is a planning quantity only, not an inference-
qualified total SE: calibration and WLE reuse the same responses, but the
coefficient draws do not represent that dependence. The maximum nonlinear
draw-mean displacement was `0.239` logits, reinforcing the need to retain the
full Person-level diagnostic rather than showing only one pooled scalar. See
`validation/cmle_wle_calibration_sensitivity_20260810/CMLE_WLE_CALIBRATION_SENSITIVITY_RESULTS.md`.

The first prospectively registered cross-engine boundary study is now complete
for 13 unanchored additive-RSM controls. Five oracle-interior cases matched
`immer::immer_cml()` within `1.78e-9` for free coordinates and `5.33e-15` for
conditional log likelihood. Seven oracle-boundary cases remained blocked even
though immer returned code zero in 3/7 and finite coefficient/SE vectors in
4/7 at `maxit=2000`. The frozen mfrmr 0.2.3 JML snapshot labelled all 7/7
boundary fits converged; this is different-estimator behavior, not CMLE
parity. TAM and sirt are likewise descriptive MML sensitivity lanes. The UI
contract therefore treats the support oracle as authoritative for finite-CMLE
existence and places cross-engine optimizer evidence only in technical detail.
See `validation/cmle_cross_engine_boundary_20260810/`.

The separately registered 12-case PCM hard-anchor cross-engine study is also
complete. Six oracle-interior cases covered unanchored, Rater-main-effect,
Criterion-main-effect, contaminated, and differential anchors. Native Python
and the matched `immer::immer_cml()` affine `W`/`b_const` mapping agreed within
`2.45e-9` for free coordinates and `3.55e-14` for conditional log likelihood;
the pre-optimization objective/gradient comparison agreed within `4.27e-14`.
All six separation or unused-declared-category cases stayed support-oracle
boundary results under every retained iteration cap. Exact anchor values and
zero anchored SEs verify the Python coordinate implementation only; they do not
validate anchor provenance or bias. mfrmr JML and TAM/sirt MML remain
different-estimand sensitivity lanes, and anchored TAM/sirt cases are typed
unsupported. See `validation/cmle_pcm_anchor_cross_engine_20260810/`.

## One-click UX contract for a later public gate

If this estimator stack clears the remaining evidence gates, the user should
still need only one explicit analysis action. A single button may orchestrate
the stages, but it must not collapse their meanings:

1. preflight declared support, duplicate units, conditional rank, informative
   Persons, missingness patterns, hard-anchor identity/provenance, and
   work/memory limits;
2. certify finite-CMLE existence with the support oracle or stop with a typed
   structural/boundary/computational repair action;
3. run exact CMLE structural calibration or stop with a numerical repair action;
4. run fit-sample WLE Person scoring only when CMLE is inference-ready; and
5. render one first-read panel with separate cards for structural calibration,
   Person scoring, exact extremes, uncertainty scope, and threshold-sensitive
   downstream results.

The result title should say “Exact CMLE calibration + fixed-calibration WLE
Person scores,” never only “CMLE results.” The first view should expose counts
of informative/conditional-extreme Persons, WLE unavailable rows, raw fit-zone
switches, raw/display mismatches, and whether calibration uncertainty is
propagated. When anchors are active it must also show anchored/free levels,
the supplied scale origin, and a persistent warning that unchanged fit cannot
validate a common anchor-origin shift. It must also show anchor observation
support and conditional rank, distinguish common-origin from relative-anchor
sensitivity, and report raw fit transitions plus raw/display disagreements.
Nominal graph connection, realized informative edge support, prefit rank, and
final fitted rank must remain separate states; complete fixed calibration must
not be rendered as data-supported recovery of disconnected free Raters.
Rejected non-finite trials and eventual finite-MLE/separation status must also
remain distinct from convergence and rank.
Expert tables and downloads may follow underneath. Any blocked
stage ends with its reason and next action; the button must not silently fall
back to JMLE, MML, approximate CMLE, or a clamped Person score. The current
public app does not yet expose this button.

The first Streamlit-free implementation of this view-model contract is now
retained in `validation/cmle_one_click_result_contract_20260810/`. Five
registered analytical states returned the same six bilingual cards. Only two
ready calibrations reached fixed-calibration WLE/MnSq, with exact reproduction
of the direct handoff for 14 Persons; the three blocked states made no
downstream scoring attempt. Its constructed `1.5004` Infit probe remains
raw-classified as `noisy` even though three-decimal display is `1.500`. This
qualifies orchestration and raw-value control only. It does not qualify a
Streamlit button or show that users understand the warnings.

The current-schema private archive increment is retained in
`validation/cmle_one_click_archive_identity_20260810/`. Its three deterministic
ready/anchored/boundary ZIPs passed all 21 replay identity checks, and 56 raw
WLE/SE/Infit/Outfit values survived CSV round-trip with identical binary64 hex
representations. Ordered input, semantic input, settings, result assets, core
archive content, and ZIP transport have separate hashes; changed, missing,
extra, or duplicate entries fail before loading. This is private
controlled-access reproduction only. Public construction, de-identification,
encryption, schema migration, and cross-platform bitwise replay remain
unqualified.

The next non-public instrument-readiness gate is retained in
`validation/cmle_one_click_comprehension_readiness_20260810/`. It generates 10
English/Japanese static six-card previews and splits 100 participant task rows
from a 100-row researcher key. Six dangerous misconceptions are
non-compensatory, including the constructed raw `1.5004` versus displayed
`1.500` MnSq decision. Thirty synthetic packets passed the registered scorer
contract, but this is not a human study: participants remain zero,
comprehension and language equivalence are unassessed, and real-browser plus
assistive-technology QA was not completed. The public surface stays withheld.

Private cognitive-interview operations are subsequently frozen in
`validation/cmle_one_click_cognitive_interview_operations_20260810/`. The
schedule exhausts all 10 state pairs within each language and balances case
position and target experience stratum. Twenty answer-free session packets,
a response-lock-first private moderator guide, a 400-row blank record schema,
and typed fail-closed validation are prepared. This repository work does not
constitute ethics approval, consent, recruitment authorization, accessibility
evidence, or a completed cognitive interview. Human participants remain zero.

Instrument content identity and directional confirmatory mechanics are then
frozen in `validation/cmle_one_click_versioned_confirmatory_gate_20260810/`.
The safety numerator is case-eligible and directional rather than identical to
all comprehension errors; each future slot contributes at most once per danger
domain. Six domains are kept separate by language. The strict one-sided 95%
Wilson zero-error boundary is above 0.10 at `n=24` and below it at `n=25`, but
the 12-cell Bonferroni sensitivity is still `0.217782` at `n=25`. This
arithmetic does not select sample size. A future frozen protocol must justify
minimum n, power/precision, invalidity/attrition, clustering, and recruitment
before data.

The subsequent no-human-data sensitivity surface preserves that boundary.
Exact pass probability inherits discrete changes in the maximum allowed error
count, so it can fall between adjacent n values. For a true dangerous-error
probability of `0.05`, the single-cell 90% first crossing is `n=224`, while the
first n sustained through the prospectively registered search maximum of 5000
is `n=260`. A conservative 12-cell union-bound projection gives `n=422` and
`n=456`. Independent retention assurance, heuristic cluster design effects,
and exposure to three distinct blocked mechanisms are reported separately.
The retained evidence therefore supports planning critique, not a selected
sample size, exact correlated-response inference, recruitment, or a public
claim. Human participants remain zero.

The registered dependence/MNAR stress then replaces the single design-effect
number with 64,000 conditional synthetic joint-gate replicates. Participant
latent dependence, fixed-size moderator/site-like clusters, heterogeneous
domains and blocked mechanisms, and outcome-dependent retention are varied
without changing the frozen Wilson rule. At n=500, the danger-under-retained
scenario turns complete-data truth near 0.10 into observed valid-record risk
`0.060315` and false reassurance `0.4035`; the blocked-mechanism hotspot gives
`0.4955`, and the combined adversarial condition `0.9675`. Thus a larger n can
increase confidence in a biased observed estimand or conceal an unacceptable
mechanism. These are constructed Gaussian-threshold stressors, not estimates
of realistic users, correlations, missingness, language effects, or
accessibility conditions. They do not select n or authorize recruitment.

The subsequent private protocol preflight converts these limitations into a
machine-enforced recruitment stop. Five parent decisions are fixed, while 13
scientific and external-authority decisions remain unresolved. A categorical
attempt ledger preserves every frozen slot and rejects PII-like/free-text
columns, denominator drift, unregistered invalidity reasons, and post-outcome
exclusions. A 72-row partial-identification surface shows that, with valid
n=500 and observed 5% error, invalid/valid ratio 5% produces an all-invalid-
dangerous Wilson upper bound of `0.118434`; excluding those records would hide
the lack of worst-case robustness. A deterministic nine-member private ZIP
binds the templates and blocked first read to content hashes. Passing this
software contract means the stop works correctly—not that ethics, recruitment,
sample size, language equivalence, accessibility evidence, or public release
is ready.

The private decision workbench then expands the 13 unresolved rows into 39
non-ranked alternatives. It exposes each option's strength, primary risk,
required evidence, and qualitative impact without issuing a repository
recommendation. A 32-edge acyclic graph places numeric sample-size decision
SCI-02 after seven scientific prerequisites. The blank worksheet is valid but
incomplete; premature dependencies, missing numeric assumptions, outcome-
informed decisions, and repository-generated external evidence fail closed.
Even a synthetically complete worksheet remains substantively unverified and
cannot make recruitment ready. The bundled HTML is bilingual, read-only,
self-contained, and network-free; no public UI or writable decision surface is
introduced.

SCI-01 now has a narrower private comparison layer. It formalizes the three
unranked estimand alternatives and projects the same 72 exact-count
partial-identification scenarios through each one. The resulting 216 rows
include 36 cases where the observed-valid analysis passes but the
all-invalid-dangerous diagnostic proxy fails and the dual view is not robust.
These differences are caused by the estimand and denominator assumptions, not
MnSq rounding or display precision. The exercise also prevents a subtle
overclaim: the existing surface does not contain scheduled-slot totals or a
prospectively adjudicated invalidity-code map, so it cannot directly calculate
the scheduled-slot composite estimand. The comparison remains unselected and
private. Its static bilingual HTML is not evidence that a public button,
assistive-technology path, or user comprehension has passed acceptance.

The following invalidity-adjudication layer keeps the 15 registered causes
mechanistically separate. Only `none` is fixed; every other code remains blank
pending the appropriate scientific, statistical, operational, accessibility,
governance, or ethics owner. A 72-by-5 diagnostic allocation surface has 36
source scenarios whose raw decision changes as more invalid records enter the
composite numerator. This is a warning against post-result classification,
not a data-driven way to choose the composite. Withdrawal cannot bypass ethics
data-use authority, accessibility barriers cannot be silently pooled,
duplicates cannot become dangerous events, and unknown reasons require a
protocol amendment. No resolved code map or public control is introduced.

The retained `cmle_wle_first_read_projection.csv` is the evidence-only view
model for this contract. It has six ordered cards: structural calibration and
Person scoring can be `research_ready`; exact extremes and calibration
sensitivity remain `caution`; interval claims and public UI remain `withheld`.
This prevents a successful numerical run from being rendered as an unsupported
green inferential conclusion.

## Public-integration gate

Before CMLE can enter the Streamlit estimator selector, all of the following
remain required:

1. complete the remaining multi-platform, in-app, extreme-parameter, and
   concurrency scale/memory checks beyond the repository scaling matrix;
2. broaden the completed small-fixture one-click calibration, Person-scoring,
   and fit/rounding card contract to high-dimensional sparse and resource-
   limited states, keeping the direct technical-comparison route distinct from
   reportable estimates and fitted-rank qualification; then extend the initial
   differential-contamination evidence to drift-distribution and prospectively
   powered designs; add group-anchor support only under its own contract;
3. extend the now-passing fit-sample CMLE-WLE bridge to an explicitly gated
   new-Person/unseen-unit design contract and replace the initial asymptotic
   calibration-sensitivity diagnostic with a qualified conditional-bootstrap
   or independent-calibration uncertainty contract;
4. downstream diagnostic qualification for conditionally calibrated
   structures;
5. extend the completed current-schema private AnalysisIdentity/archive
   contract to schema migration, independently governed public disclosure,
   localization, Help, and Streamlit upload/download surfaces;
6. RSM/PCM truth-recovery simulations beyond deterministic implementation
   parity, high-dimensional sparse-PCM and clean released-mfrmr reproduction
   beyond the completed additive-RSM and 12-case PCM/hard-anchor cross-engine
   studies; and
7. a UI preflight that shows informative Persons, unique design patterns,
   conditional rank/nullity, declared category support, and estimated
   computational work before fitting; and
8. bilingual cognitive interviews, an instrument-frozen pilot, and a fresh
   confirmatory study in which every critical item within each language clears
   its registered dangerous-error bound, together with real-browser and
   assistive-technology accessibility QA.

Composite conditional likelihood, if later added, must be a separately named
estimator with its own uncertainty and information-criterion contract. It must
not silently replace exact CMLE when an exact run is too large.
