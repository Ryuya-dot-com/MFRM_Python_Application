# MFRM Streamlit

Development checkout: `integration/unified-app`. The public `main` baseline is
`d8558fc` (`0.2.16-beta`); the inherited version label below is not a release-order
indicator. See the [development integration ledger](docs/development_integration.md)
for preserved branches, completed transfers, and pending release gates.

Standalone Python beta application for Many-Facet Rasch Model estimation in Streamlit.

This app estimates, diagnoses, stress-tests, and reports MFRM analyses in
Python. Other software may be cited as methodological history, but it is not
called, imported, handed a job, or used as a public release gate.

## Status

- Release status: public beta / research preview (**v0.2.15-beta**)
- Runtime engine: standalone Python
- Primary entrypoint: `streamlit_app.py`
- Intended use: exploratory analysis, teaching, reporting support, and research workflow prototyping
- Not intended as: a validated drop-in replacement for FACETS, TAM, sirt, mirt, or `mfrmr`

The result workflow makes the FACETS-complementary boundary explicit.  Before
technical tables it names the current estimator's estimand and likelihood
basis; JMLE, MML, and exact CMLE are never ranked by cross-basis likelihood,
AIC, or BIC.  The Data view also includes an outcome-blind Person-by-Rater
assignment audit.  It describes density, exposure, co-rating, and overlap but
does not diagnose MCAR/MAR/MNAR or informative assignment.  The associated
assignment sensitivity workflow is now an on-demand, fail-closed paired
parametric runner for eligible JMLE/MML RSM/PCM fits.  It preserves Person and
Rater assignment degrees, Rater response-row exposure, and direct-overlap
connectivity; changes no observed Score; generates new paired responses; and
refits only the current estimator.  Exchangeable equal-context blocks use a
connected degree-preserving 2-switch path: requested 0/25/50/75/100% alignment
doses are mapped to achievable points, with one common row-level random draw
shared across every dose in a replicate.  Unequal-context blocks fall back to
a binary MILP endpoint that preserves every Rater-by-context row margin exactly
and locks observed spanning-tree overlap witnesses.  The MILP first proves the
observed assignment is feasible and exposes an endpoint only when HiGHS
certifies a zero-gap optimum and all post-solve invariants pass.  Intermediate
MILP doses are intentionally unavailable.  Dose/path labels are not estimated
assignment propensities or unconstrained global worst-case scales.  It is a local fitted-model
stress test, not an MNAR diagnosis, bias correction, causal estimate, fairness
verdict, or cross-estimator ranking.  FACETS two-decimal display values are not
treated as raw precision for fit recalculation or boundary adjudication.

MML assignment refits now expose two non-interchangeable Person generators.
The default holds source posterior/EAP measures fixed and remains a local
conditional reference analysis.  The optional population mode draws fresh
normal order statistics using the fitted fixed or free MML population SD and
maps them to the strict source-EAP ordering.  Thus the same generated Person
coordinates are shared across assignment scenarios and the rank-to-assignment
relation is preserved, while population-scale realizations vary by replicate.
This is a rank-preserving fitted-assignment stress—not an unconditional new
sample, an empirical assignment model, or evidence that EAP values are truth.
JMLE cannot select this MML population generator.

Repository research note: a native exact RSM/PCM CMLE structural-calibration core
is under deterministic validation in [`mfrm_app/cmle.py`](mfrm_app/cmle.py).
It is not exposed in the Streamlit estimator selector and is not part of the
current public beta capability claim. Its scope and promotion gate are
documented in [`docs/cmle_phase0_design.md`](docs/cmle_phase0_design.md). The
current research core obtains rank and covariance from exact conditional
sufficient-statistic moments rather than a finite-difference Hessian. Its
repository scaling envelope also enforces independent work and conservative
512 MiB core-memory guards; these checks do not yet authorize public UI use.
The repository-only repeated-simulation lane is separately specified in
[`validation/OPERATING_CHARACTERISTICS_PROTOCOL_20260809.md`](validation/OPERATING_CHARACTERISTICS_PROTOCOL_20260809.md);
its retained smoke profile verifies failure accounting and negative controls,
not estimator performance or cross-package parity. Repository-only companion
results now include strict mfrmr JMLE, matched unanchored
Python/`immer_cml` numerical comparisons, and separately labelled TAM/sirt MML
sensitivity adapters. TAM retains the additive RSM surface while integrating a
normal Person population; sirt uses virtual-item PCM thresholds. Neither is
described as JMLE/CMLE parity. All lanes retain optimizer/readiness,
floating-boundary, constraint, and unsupported-scope differences rather than
changing the public capability claim.

The separately registered PCM hard-anchor mapping has also passed its
repository-only contract. Across six oracle-interior unanchored, Rater-anchor,
Criterion-anchor, contaminated-anchor, and differential-anchor cases, native
Python and a matched `immer_cml` `W`/`b_const` parameterization differed by at
most `2.45e-9` in free coordinates and `3.55e-14` in conditional log
likelihood. Six separation or unused-declared-category cases remained blocked
by the exact support oracle at all retained iteration caps. The mfrmr 0.2.3
JML and TAM/sirt MML outputs are retained only as different-estimand
sensitivity evidence; anchored TAM/sirt cases are explicitly unsupported.
Exact anchor enforcement does not establish anchor validity, and this evidence
does not enable the public CMLE UI. See the retained
[`critical review`](validation/cmle_pcm_anchor_cross_engine_20260810/CMLE_PCM_ANCHOR_CROSS_ENGINE_CRITICAL_REVIEW.md).

A Streamlit-free one-click result contract now joins the guarded calibration
workflow to fixed-calibration Warm WLE and untrimmed Person MnSq without an
automatic estimator fallback. Five registered states all returned the same six
bilingual cards. Only the two ready calibrations reached Person scoring; the
input, structural, and finite-existence stops did not call downstream scoring.
The 14 ready Person handoffs reproduced direct WLE/SE/Infit/Outfit values
exactly. A constructed raw Infit of `1.5004` remained `noisy` even though its
three-decimal display was `1.500` (`acceptable` if misused). This is a
presentation-data contract, not task-comprehension evidence, and the Streamlit
surface remains disabled. See the
[`one-click critical review`](validation/cmle_one_click_result_contract_20260810/CMLE_ONE_CLICK_RESULT_CONTRACT_CRITICAL_REVIEW.md).

The same research result can now be saved as a deterministic **private
reproduction** ZIP and reopened or replayed without Streamlit. Three registered
ready/boundary archives passed 21/21 replay identity checks; 56 retained
WLE/SE/Infit/Outfit values survived 17-significant-digit CSV round-trip with
identical binary64 hexadecimal values. Input order, semantic row-order
identity, configuration, results, archive content, and ZIP transport have
separate hashes. Changed, missing, extra, or duplicate entries are rejected
before loading. These ZIPs contain raw rows and identifiers: public-mode
construction fails closed and no anonymization, encryption, or sharing approval
is implied. See the
[`private-archive critical review`](validation/cmle_one_click_archive_identity_20260810/CMLE_ONE_CLICK_PRIVATE_ARCHIVE_CRITICAL_REVIEW.md).

A non-public bilingual comprehension instrument is now prepared around those
same five analytical states. Ten standalone English/Japanese HTML previews
retain six textual-status cards each; participant and separately held
researcher-key packets contain 100 rows apiece. Six dangerous misconceptions
are non-compensatory: false-green calibration, invented Person scores,
display-rounded MnSq reclassification, ignored calibration uncertainty,
anchor-validity overclaim, and treating a private archive as publicly
shareable. Thirty synthetic packets verified fail-closed scoring, including
duplicate, missing, unknown-option, and invalid-time responses; 40 selected
tests passed. This is instrument-readiness evidence only. Human participants
remain zero, comprehension and Japanese/English equivalence are not assessed,
real-browser/assistive-technology QA was not completed in this run, and the
public surface remains disabled. See the
[`comprehension-readiness critical review`](validation/cmle_one_click_comprehension_readiness_20260810/HTML_PREVIEW_INSTRUMENT_CRITICAL_REVIEW.md).

The subsequent private cognitive-interview operations kit is also prepared,
without recruiting anyone. Within each language, 10 planned slots cover every
pair of the five analytical states once; each state appears four times, twice
in each presentation position and twice in each target experience stratum.
Twenty answer-free session packets contain 400 assigned tasks in total. A
separate private moderator guide requires response lock before probing, and an
exact-column validator protects a 400-row blank record template against direct
identifier columns, assignment drift, missing/duplicate items, unknown codes,
unreviewed notes, early probes, and invalid timing/confidence values. Fourteen
synthetic validator cases and 47 selected tests passed. The repository does
not supply ethics approval or recruitment authorization; human participants
remain zero. See the
[`operations critical review`](validation/cmle_one_click_cognitive_interview_operations_20260810/COGNITIVE_INTERVIEW_OPERATIONS_CRITICAL_REVIEW.md).

Instrument version identity and future confirmatory gate mechanics are now
frozen separately, still with no human data. The composite identity covers the
participant task bytes, private key bytes, 10 preview hashes, and directional
dangerous-error rules; any mutation changes the identity and cross-content
pooling is rejected. The safety numerator is narrower than any comprehension
error: `display_acceptable` is dangerous for raw `1.5004` versus displayed
`1.500`, while `cannot_decide` is incorrect but is not relabelled a dangerous
display-driven reclassification. Six domains are evaluated separately in
English and Japanese (12 cells). At zero errors, the one-sided 95% Wilson upper
bound is `0.101310` for `n=24` and `0.097654` for `n=25`; the formative pilot
therefore cannot promote the UI. The `n=25` arithmetic is not a sample-size
recommendation: its 12-cell Bonferroni sensitivity bound is still `0.217782`.
The future minimum valid n remains deliberately unregistered. See the
[`versioned-gate critical review`](validation/cmle_one_click_versioned_confirmatory_gate_20260810/VERSIONED_CONFIRMATORY_GATE_CRITICAL_REVIEW.md).

A prospective, no-human-data sample-size sensitivity surface now makes that
non-selection explicit. Because the maximum Wilson-passing error count changes
discretely, exact pass probability is sawtoothed rather than monotone. Under an
assumed true dangerous-error probability of `0.05`, a single primary cell first
reaches 90% pass probability at `n=224`, but does not remain at or above 90%
through the registered `n=5000` search until `n=260`. The conservative 12-cell
union-bound projection gives `n=422` and `n=456`, respectively. Separate
sensitivities show that 200 valid observations require 231 planned slots under
90% homogeneous retention with 95% assurance, while an effective `n=100`
becomes a heuristic nominal `n=145` for mean cluster size 10 and ICC 0.05.
Seventy-five slots per language expose each of the three blocked danger
mechanisms only 25 times. These are planning diagnostics, not a chosen sample
size, power claim, exact clustered correction, or authorization to recruit.
Human participants and the public surface remain zero/withheld. See the
[`sample-size sensitivity critical review`](validation/cmle_one_click_confirmatory_sample_size_sensitivity_20260810/CONFIRMATORY_SAMPLE_SIZE_SENSITIVITY_CRITICAL_REVIEW.md).

The next registered synthetic stress replaces the cluster heuristic with an
explicit, conditional Gaussian-threshold generator for participant dependence,
fixed-size moderator/site-like clusters, domain and blocked-mechanism
heterogeneity, and outcome-dependent invalidity. Eight scenarios by four
planned-slot values and 2,000 replicates produced 64,000 synthetic joint-gate
replicates. At 300 slots per language with marginal `p=0.05`, all-12 primary
pass rates were `0.5480` when independent, `0.5660` with participant latent
dependence, and `0.5125` with added clustering. More importantly, at `n=500`
an MNAR stress with complete-data truth near `0.10` but lower retention of
dangerous responses reduced the observed valid-record rate to `0.060315` and
produced false reassurance in `0.4035` of replicates (Monte Carlo Wilson 95%
`0.3822`–`0.4252`). A pooled blocked-mechanism hotspot gave `0.4955`, and the
combined adversarial condition `0.9675`. These are constructed stress results,
not estimated user behavior or chosen n; they show that larger samples do not
repair informative invalidity or a misspecified pooled estimand. See the
[`dependence/MNAR critical review`](validation/cmle_one_click_confirmatory_dependence_stress_20260810/CONFIRMATORY_DEPENDENCE_STRESS_CRITICAL_REVIEW.md).

A deterministic private pre-recruitment bundle now turns those warnings into a
fail-closed protocol workflow. Its 18-row decision register retains five fixed
repository decisions and 13 unresolved scientific, ethical, governance,
cluster, mechanism, language, and accessibility decisions. A zero-row exact-
schema denominator ledger rejects direct identifier/free-text columns, frozen-
slot gaps or duplicates, unknown invalidity reasons, post-outcome exclusions,
and incomplete PII/audit attestations. A 72-row partial-identification surface
retains all-invalid-safe and all-invalid-dangerous bounds instead of silently
dropping invalid records. With valid `n=500` and observed dangerous-error rate
5%, the registered grid is worst-case robust at invalid/valid ratio 2.5%, but
not at 5%, where the all-invalid-dangerous Wilson upper bound is `0.118434`.
The nine-member private ZIP is byte-deterministic and tamper checked. The
software contract passed with 91 selected tests precisely because the first
read remains `BLOCKED`: ethics approval, recruitment authorization, minimum n,
human data, confirmatory results, and public UI are all absent. See the
[`protocol-preflight critical review`](validation/cmle_one_click_confirmatory_protocol_preflight_20260810/CONFIRMATORY_PROTOCOL_PREFLIGHT_CRITICAL_REVIEW.md)
and the
[`private preflight ZIP`](validation/cmle_one_click_confirmatory_protocol_preflight_20260810/private_confirmatory_protocol_preflight.zip).

The 13 blockers now have a private read-only decision workbench rather than an
opaque checklist. Thirty-nine options expose strengths, primary risks,
required evidence, and qualitative effects on n, estimand, participant burden,
accessibility, and public claims. None is ranked, recommended, defaulted, or
auto-selectable. Thirty-two prerequisite edges are acyclic; numeric sample-size
decision `SCI-02` cannot resolve until seven scientific prerequisites do.
Unknown or cross-decision options, premature child decisions, missing numeric
assumptions, repository-generated external evidence, and decisions after
outcome inspection fail closed. Even a fully completed synthetic worksheet
remains `SubstantiveEvidenceVerified=false` and `RecruitmentReady=false`. The
self-contained bilingual HTML has no forms, scripts, external resources, or
writable controls; its deterministic seven-member ZIP passed tamper checks and
105 selected tests. See the
[`offline workbench`](validation/cmle_one_click_confirmatory_decision_workbench_20260810/protocol_decision_workbench.html),
[`critical review`](validation/cmle_one_click_confirmatory_decision_workbench_20260810/CONFIRMATORY_DECISION_WORKBENCH_CRITICAL_REVIEW.md),
and
[`private workbench ZIP`](validation/cmle_one_click_confirmatory_decision_workbench_20260810/private_confirmatory_decision_workbench.zip).

SCI-01 now has a separate private estimand-comparison memo. It maps the same
72 exact-count sensitivity scenarios to all three registered alternatives,
creating 216 projections without ranking or selecting an option. In 36
scenarios, the observed-valid gate passes while the all-invalid-dangerous
proxy fails and the dual analysis is explicitly non-robust; raw integer counts
and unrounded Wilson bounds, never displayed rounding, drive those labels. A
critical denominator gap is now machine-visible: the existing surface lacks a
complete scheduled-slot count and a prospective invalidity-code-to-composite
map, so Option B's `(D+I)/(V+I)` output is only a diagnostic proxy, not a
calculated scheduled-slot estimand. The offline memo and deterministic ZIP
remain private, read-only, and `BLOCKED`; no option, sample size, recruitment,
or public button is enabled. See the
[`SCI-01 offline memo`](validation/cmle_one_click_confirmatory_sci01_estimand_memo_20260810/sci01_estimand_comparison.html),
[`critical review`](validation/cmle_one_click_confirmatory_sci01_estimand_memo_20260810/SCI01_ESTIMAND_MEMO_CRITICAL_REVIEW.md),
and
[`same-data figure`](validation/cmle_one_click_confirmatory_sci01_estimand_memo_20260810/sci01_estimand_projection.png).

The shared SCI-01/SCI-05 invalidity bottleneck now has a separate private
adjudication workbook. It preserves all 15 registered reason codes and fixes
only `none`; the other 14 remain unresolved. Withdrawal requires external
ethics/data-use authority, accessibility barriers require an accessibility
owner and separate-stratum sensitivity, duplicate records cannot be dangerous
composite events, and unregistered reasons require amendment. Crossing the 72
count scenarios with five diagnostic fractions of invalid records classified
as composite events gives 360 rows; 36 scenarios cross the raw `< 0.10` gate
within that grid. This is evidence that post-outcome classification could
change the conclusion, not a basis for choosing a fraction. The
[`private invalidity workbench`](validation/cmle_one_click_confirmatory_invalidity_adjudication_20260810/invalidity_adjudication_workbench.html)
remains unselected and recruitment-blocking.

A prospective machine-readable precision plan and a separate 20-replicate
Python pilot are now retained under `validation/`. The pilot completed 160/160
fits in about 82 recorded fit-seconds, but it blocks rather than promotes the
study: only 124/160 focal decisions met the registered app eligibility rule,
the two sparse conditions each produced 4/20 eligible decisions, and 0/160
fits passed the post-pilot `1e-4` terminal-gradient sup-norm sensitivity.
Row-level floating-point auditing also found five fit classifications that
would change if 3-decimal display values drove the decision. See
[`PILOT20_ASSESSMENT.md`](validation/operating_characteristics_pilot20_20260809/PILOT20_ASSESSMENT.md).
The finding requires a separately qualified strict Python JMLE mode and sparse
design review before a fixed confirmatory manifest. The prospective strict
lane subsequently passed its frozen `1e-4` numerical gate in 160/160 runs, but
the accompanying movement and exact eta-design-rank audits found that all 40
sparse runs were structurally nonidentified (nullity 7; eight Person-Rater
components). The app now reports optimizer `Converged` separately from
`InferenceReady`, fails JMLE inference readiness when this pre-fit rank audit
fails, and withholds conditional bias output for the affected fit. A retained
160/160 integration check verifies the guard without changing optimizer
controls, estimates, or automatically switching JMLE to MML. This still does
not enable a public evidence button or performance claim; the audit is limited
to the Person/facet eta block and does not qualify thresholds or GPCM slopes.
An additional exact-score boundary guard now distinguishes an optimizer or
constraint value from a finite JMLE Person measure. In the frozen evidence it
matched 8,320/8,320 Person rows; four boundary Persons occurred among the 120
structurally identified runs. Their technical `Estimate` is retained for
reproducibility, but `ReportableEstimate` is withheld. The rule uses exact
all-minimum/all-maximum integer patterns, never an estimate-magnitude or display
rounding cutoff, and does not silently select a finite correction.

A repository-only fixed-calibration Warm WLE lane now provides the first
explicit correction comparison without relabelling WLE as a finite JMLE MLE.
The generic RSM/PCM/GPCM scorer matched `TAM::tam.mml.wle2` for all 16 frozen
Person fixtures (maximum theta and SE differences `1.68e-11` and `1.18e-11`).
On the 120 structurally identified Stage-B2 runs, 7,600 Person scores were
rescored with facets and steps held fixed; the four exact-extreme patterns
moved by about 24--25 logits, while the median absolute interior shift was
`0.009` logits. Raw Person-fit zones changed in 51/17,360 comparisons, and
three-significant-digit display values disagreed with raw fit zones 54 times
across the two estimator surfaces. The focal Holm and combined strong-bias
decisions did not change, but the practical `|bias| >= 0.50` rule changed once
(`0.500839` to `0.498068`). All 40 structurally rank-deficient sparse runs
remain withheld: finite WLE Person scores cannot identify their facet
calibration. These findings do not enable a Streamlit WLE option, automatic
fallback, or public equivalence claim; see
[`FIXED_CALIBRATION_WLE_DOWNSTREAM_RESULTS.md`](validation/fixed_calibration_wle_downstream_20260809/FIXED_CALIBRATION_WLE_DOWNSTREAM_RESULTS.md).
The first explicit exact-CMLE-to-WLE fit-sample bridge also passes its frozen
RSM/PCM handoff checks with zero bridge/direct-score and category-surface
differences; it remains repository-only and its WLE SE treats CMLE calibration
as fixed. See
[`CMLE_WLE_BRIDGE_RESULTS.md`](validation/cmle_wle_bridge_20260810/CMLE_WLE_BRIDGE_RESULTS.md).
The first frozen calibration-sensitivity gate now perturbs the exact-CMLE free
coordinates with 2,000 asymptotic covariance draws per RSM/PCM model. All
32,000 draw-Person scores returned. Median Person calibration-draw SD was
`0.0616`--`0.0801` logits by model, while the largest value was `0.518` logits
and the largest ratio to conditional WLE SE was `0.427`; the largest effects
occurred at exact-extreme patterns. These draw quantiles are not confidence or
credible intervals, and the quadrature sum is not an inference-qualified total
SE because the fit-sample CMLE and WLE stages reuse responses. A six-card
machine-readable first-read projection therefore marks numerical research
success separately from caution and withheld inference/UI states. See
[`CMLE_WLE_CALIBRATION_SENSITIVITY_RESULTS.md`](validation/cmle_wle_calibration_sensitivity_20260810/CMLE_WLE_CALIBRATION_SENSITIVITY_RESULTS.md).
The next bootstrap stage was designed and frozen before implementation. Its
repository-only core separates a fixed-Person-total conditional-pattern lane
from a joint plug-in response lane, retains refit rank/readiness failures, and
applies seed/work/output guards; neither lane is pre-labelled a total SE or
confidence interval. The frozen 200-replicate RSM/PCM matrix has now completed
all 800 attempted refits with full rank/readiness and WLE availability (the
per-cell Wilson 95% lower bound is `0.981`, not proof of a population rate of
one). Fixed-score totals had zero mismatches. Median Person bootstrap SD was
`0.0602` logits in the fixed-score lane and `0.615` in the joint lane; the
joint lane produced 441 exact-extreme state transitions and the maximum
bootstrap-mean displacement was `0.370` logits. The sampler contract passed,
but these dispersions are not a variance decomposition and percentile columns
are not confidence intervals. See
[`docs/cmle_wle_bootstrap_design.md`](docs/cmle_wle_bootstrap_design.md) and
[`CMLE_WLE_BOOTSTRAP_PILOT_RESULTS.md`](validation/cmle_wle_bootstrap_pilot_20260810/CMLE_WLE_BOOTSTRAP_PILOT_RESULTS.md).

A separately registered downstream gate now validates untrimmed fixed-theta
CMLE-WLE Person Infit/Outfit against `sirt::pcm.fit`. Across RSM and PCM, 32
Person rows and 162 administered observations—including four exact-extreme
Persons and non-rectangular missingness—matched within `4.89e-15` for Infit and
`4.00e-15` for Outfit; expected categories and variances also matched at
binary64 noise. TAM uses the same Infit and pre-trim residual contributions,
but its default JML fit path trims unusually large Outfit contributions, so
default TAM Outfit is not claimed identical. The corrected 800-replicate fit
extension retained 12,800 successful Person-replicates, zero fixed-score total
mismatches, and exact persisted WLE identity with the original pilot. It found
2,372 Infit and 2,454 Outfit raw-class transitions (2,554 Person-replicates had
at least one), plus 24 display-rounding boundary statistics and nine cases in
which a three-decimal display implied a different class. All decisions used
unrounded MnSq. The first extension run remains retained as failed because an
asymmetric CSV-versus-memory comparison introduced a one-ULP identity
difference; a registered amendment kept the exact-zero gate unchanged and the
full corrected rerun passed under symmetric persisted comparison. ZSTD,
p-values, threshold optimality, coverage intervals, total SE, and Streamlit
integration remain withheld. See
[`CMLE_WLE_PERSON_FIT_RESULTS.md`](validation/cmle_wle_person_fit_20260810/CMLE_WLE_PERSON_FIT_RESULTS.md)
and
[`CMLE_WLE_BOOTSTRAP_FIT_EXTENSION_RESULTS.md`](validation/cmle_wle_bootstrap_fit_extension_corrected_20260810/CMLE_WLE_BOOTSTRAP_FIT_EXTENSION_RESULTS.md).
The accompanying
[`critical decomposition`](validation/cmle_wle_bootstrap_fit_extension_corrected_20260810/CMLE_WLE_BOOTSTRAP_FIT_CRITICAL_REVIEW.md)
shows why the aggregate transition count is not a misfit-prevalence estimate.

A post-result sensitivity audit now separates display rounding from the
substantive threshold convention. Recomputing from raw values reproduced every
canonical class and transition. Counterfactually classifying rounded values
produced nine replicate-statistic class mismatches and five either-statistic
transition disagreements at three decimals, but none at six decimals. Across
the frozen 125 threshold triplets, aggregate transition counts ranged from
1,789--3,461 for Infit, 1,877--3,550 for Outfit, and 1,975--3,675 for either
statistic. A registered diagnostic addendum also found that moving only the
noisy/distorting boundary from 1.90 to 2.10 left the binary transition totals
exactly unchanged while reallocating 37 Infit and 51 Outfit replicate labels
from `distorting` to `noisy`. Thus neither rounded display nor a binary
changed/not-changed projection is an adequate audit trail. These are
descriptive post-result surfaces, not threshold validation. See the
[`threshold surface`](validation/cmle_wle_fit_threshold_surface_20260810/CMLE_WLE_FIT_THRESHOLD_SURFACE_RESULTS.md)
and
[`class-decomposition addendum`](validation/cmle_wle_fit_threshold_surface_severity_20260810/CMLE_WLE_FIT_THRESHOLD_SEVERITY_ADDENDUM.md).

The first prospectively registered known-truth design pilot now retains 160
runs, 624,000 response rows, and 32,000 Persons for both fit-sample WLE and
generating-theta sensitivity. Under the raw `Infit > 1.50 or Outfit > 1.50`
rule, clean-Person mean replicate rates were 3.4%--3.7% with 24 observations
but 10.1%--10.8% with six; affected-Person detection ranged from 3.8% for the
frozen threshold-heterogeneity mechanism to 46.6% for dense random response.
WLE and generating-theta fit flags disagreed for 2,333/32,000 Persons, and
non-extreme WLE RMSE ranged from 0.263 to 0.720 logits. The nominal 125-triplet
grid reduces to 5, 5, 25, and 5 effective input configurations for the four
declared rules. These 20-replicate results debug the design; they do not
validate a threshold or authorize a public traffic light. See the
[`known-truth pilot`](validation/cmle_wle_fit_known_truth_pilot_20260810/CMLE_WLE_FIT_KNOWN_TRUTH_PILOT_RESULTS.md),
[`bias/dimension addendum`](validation/cmle_wle_fit_known_truth_pilot_addendum_20260810/CMLE_WLE_FIT_KNOWN_TRUTH_PILOT_ADDENDUM.md),
and
[`critical review`](validation/cmle_wle_fit_known_truth_pilot_addendum_20260810/CMLE_WLE_FIT_KNOWN_TRUTH_CRITICAL_REVIEW.md).

The next repository-only Phase-B increment now implements native hard facet
anchors inside the exact CMLE likelihood rather than as a post-fit recentering.
Its prospectively amended engineering smoke reused 24,000 retained response
rows for 32 same-byte fits; all returned inference-ready, preserved supplied
anchors exactly with SE zero, and passed rank, parameter-count, persistence,
and raw-fit recomputation gates. Two replicates per condition do not establish
an anchor percentage, robustness, or cross-engine performance. A separately
labelled post-result diagnostic found a more fundamental warning: shifting all
Rater anchors by the same +0.25 moved WLE by +0.25 but left conditional
log-likelihood and Person fit invariant to numerical noise (zero raw flag
mismatches). Fit indices therefore cannot validate the absolute anchor origin.
See the [`Phase-B smoke`](validation/cmle_native_hard_anchor_smoke_20260810/CMLE_NATIVE_HARD_ANCHOR_SMOKE.md)
and [`common-shift diagnostic`](validation/cmle_native_hard_anchor_common_shift_20260810/CMLE_HARD_ANCHOR_COMMON_SHIFT_DIAGNOSTIC.md).

The subsequent prospectively frozen differential-anchor stress reused 120,000
retained response rows for 320 exact-CMLE fits across four dense/sparse,
clean/random-response conditions and eight anchor scenarios. All fits were
inference-ready, but the substantive result is cautionary. Correct nested
anchor sets did not improve WLE recovery monotonically; the anchor share was
inseparable from which Rater levels were fixed. A zero-mean differential
anchor error changed raw either-upper Person-fit decisions for 2.4%--15.8% of
paired rows by condition/group, while a common +0.25 origin error changed none.
Across 64,000 Person-scenario rows, using 3-decimal displayed MnSq instead of
raw values would reverse 12 decisions; 6 decimals reproduced the raw decisions
in this dataset. This is not a universal precision guarantee. The evidence
therefore supports provenance, connectivity, raw-boundary, and contamination
sensitivity reporting—not an anchor-percentage recommendation or public CMLE
button. See the
[`registered stress report`](validation/cmle_native_anchor_differential_stress_20260810/CMLE_NATIVE_ANCHOR_DIFFERENTIAL_STRESS.md)
and its
[`critical review`](validation/cmle_native_anchor_differential_stress_20260810/CMLE_NATIVE_ANCHOR_DIFFERENTIAL_CRITICAL_REVIEW.md).

A second prospectively registered stress then varied the Rater bridge topology
and randomized which correct levels were anchored. Its 350/350 comparisons
matched a transparent realized-component prediction to the exact prefit rank:
partial anchors did not identify isolated Raters, while two random anchors
identified a two-component design only in the 7/10 replicates where both
components were covered. Of 287 prefit-eligible cases, 286 returned and 274
were inference-ready. All 13 shortfalls occurred in the one-Person-per-edge
minimal chain: one non-finite conditional-normalizer failure and 12 returned
fits with deficient final information rank. Thus graph connectedness and
prefit rank are necessary but not sufficient when bridge support is extremely
thin. Across 82,200 inference-ready Person rows, 3-decimal display-driven MnSq
changed 13 raw decisions; 6 decimals changed none in this dataset. See the
[`connectivity report`](validation/cmle_native_anchor_connectivity_20260810/CMLE_NATIVE_ANCHOR_CONNECTIVITY_STRESS.md)
and
[`critical review`](validation/cmle_native_anchor_connectivity_20260810/CMLE_NATIVE_ANCHOR_CONNECTIVITY_CRITICAL_REVIEW.md).

The prospectively registered finite-domain remediation now rejects and counts
non-finite BFGS trial points without clipping the exact likelihood or relaxing
readiness. On the same retained bytes, all 287 prefit-eligible fits returned
and the same 274 remained inference-ready. The former exception recorded 12
invalid trials and returned fail-closed at rank 8/10, nullity 2, and gradient
sup norm `0.00296`. The 286 previously returned fits retained identical ready
states; 82,200 shared Person rows retained every raw/rounded fit flag. Maximum
structural-estimate and conditional-loglikelihood differences were
`3.55e-15` and `2.27e-13`. This fixes returnability, not finite-MLE existence:
a separate adversarial binary separation check shows that full finite-point
rank and a small gradient can still mask a likelihood maximum at infinity.
See the
[`remediation report`](validation/cmle_finite_domain_optimizer_remediation_20260810/CMLE_FINITE_DOMAIN_OPTIMIZER_REMEDIATION.md)
and
[`critical review`](validation/cmle_finite_domain_optimizer_remediation_20260810/CMLE_FINITE_DOMAIN_OPTIMIZER_REMEDIATION_CRITICAL_REVIEW.md).

The finite-MLE existence gate has now passed a separately prospective research
stress without entering the public estimator. An exhaustive convex-support
audit first matched all 93 registered binary/RSM/PCM/anchor controls, but its
50,000-configuration research cap left 227/350 retained cases unavailable. A
subsequently registered fixed-score support oracle replaced materialization
with exact dynamic-programming separation and LP constraint generation. It
matched 1,184/1,184 exhaustive support maxima, 93/93 fixture statuses, 123/123
completed retained exhaustive statuses, and all six high-cap sentinels. The
same-byte 350-case replay classified 63 structural, 13 boundary, and 274
interior cases with zero unavailable or tolerance-unstable results. All 274
current inference-ready fits were interior; the 13 boundary cases were
non-ready. Median theoretical configurations were 97,040 versus 163 generated
constraints, and per-eligible-case time was 0.521 seconds at P95 (0.568 maximum)
on the evidence machine. This certifies only conditional finite-existence
geometry under the implemented coordinates and anchors: it does not validate
fit, anchor bias, MnSq thresholds, or R-package equality. See the
[`oracle stress report`](validation/cmle_finite_mle_oracle_20260810/CMLE_FINITE_MLE_SUPPORT_ORACLE_STRESS.md)
and
[`critical review`](validation/cmle_finite_mle_oracle_20260810/CMLE_FINITE_MLE_SUPPORT_ORACLE_CRITICAL_REVIEW.md).
The oracle is now connected to the repository-only `fit_cmle` readiness path
under a third prospective contract. All 287 same-byte eligible fits returned;
274 oracle-interior fits remained ready and all 13 boundary fits remained
non-ready with a dedicated reason. Independent oracle status matched 287/287,
102 selected tests passed, and 82,200 Person rows retained every raw and
rounded-counterfactual fit decision. Estimate, SE, likelihood, WLE, Infit, and
Outfit changes were at most `2.27e-13`. The core still retains technical
optimizer outputs on a boundary for research reproducibility; a future UI
orchestrator should stop before optimization and explain the boundary. See the
[`integration report`](validation/cmle_finite_mle_readiness_integration_20260810/CMLE_FINITE_MLE_READINESS_INTEGRATION.md)
and
[`integration critical review`](validation/cmle_finite_mle_readiness_integration_20260810/CMLE_FINITE_MLE_READINESS_INTEGRATION_CRITICAL_REVIEW.md).
A fourth prospective layer now provides a structured, bilingual pre-
optimization workflow in
[`mfrm_app/cmle_workflow.py`](mfrm_app/cmle_workflow.py). It returns the same
four-stage schema for invalid input, structural nonidentification, a finite-MLE
boundary, an unavailable existence audit, numerical non-readiness, and a ready
calibration. On the same 350 cases it stopped 63 structural and 13 boundary
cases before optimization and optimized only the 274 interior cases; all
terminal states matched the frozen evidence and no forbidden optimizer call
occurred. P95 time was 0.068 seconds for a design block, 0.306 for a boundary
stop, and 1.084 for a ready interior workflow. Japanese/English text is
machine-checked but not comprehension-tested. Person scoring is only marked
available and is not run in workflow v1. See the
[`workflow report`](validation/cmle_structured_workflow_20260810/CMLE_STRUCTURED_WORKFLOW_STRESS.md)
and
[`workflow critical review`](validation/cmle_structured_workflow_20260810/CMLE_STRUCTURED_WORKFLOW_CRITICAL_REVIEW.md).

A separately registered boundary-semantics replay now distinguishes matched
CMLE parity from different-estimator sensitivity. On five oracle-interior
additive-RSM controls, native Python and `immer::immer_cml()` free coordinates
agreed within `1.78e-9` and conditional log likelihood within `5.33e-15`.
On seven oracle-boundary controls, however, `immer` returned optimizer code
zero in 3/7 cases and finite coefficient/SE vectors in 4/7 at `maxit=2000`;
none became reportable after the finite-existence gate. The frozen dirty
`mfrmr 0.2.3` JML snapshot returned and labelled all 7/7 boundary fits
converged, with category-support warnings and a maximum absolute estimate of
`103.75`. TAM and sirt were retained only as MML sensitivity lanes, not
numerical parity targets. Thus neither finite output nor nominal convergence
may override the support oracle, and no rounded value or MnSq was used to
classify existence. See the
[`cross-engine decision`](validation/cmle_cross_engine_boundary_20260810/decision.json)
and
[`critical review`](validation/cmle_cross_engine_boundary_20260810/CMLE_CROSS_ENGINE_BOUNDARY_CRITICAL_REVIEW.md).

Before using results for high-stakes scoring, placement, certification, employment, or institutional decisions, cross-check the analysis with an established workflow and document the model assumptions.

## Active Checkout

The current public-deployment line for the latest UI, help, FACETS-style
yardstick, and reproducibility refinements is `main`. The earlier
`slope-aware-bias-inference` development branch has been merged into `main`.

If you keep more than one local checkout of this repository, verify the one you
are running before starting Streamlit:

```bash
git status --short --branch
git branch --show-current
git rev-parse --show-toplevel
```

The app should be launched from the repository root that contains the updated
`streamlit_app.py`, `locales/en.json`, and `locales/ja.json`. Seeing the old
"Try another scenario" sidebar buttons or the old "Show yardstick labels
directly on plot" control means the process is running an older checkout or the
wrong commit, not the current `main`.

The app header now includes a short `source commit: <sha>` badge. Use it to
confirm that a hosted Streamlit app is actually running the intended GitHub
commit after a push or reboot.

## What's new in the current beta line

The current app label is v0.2.15-beta. This beta line adds the v0.2.14
sample-scenario and quick-results-bundle work, then layers on the
release-readiness and UX/documentation pass documented in `CHANGELOG.md`.
Key current-line improvements include:

- JMLE results now distinguish numerical convergence from structural
  inference readiness. Before optimization, an exact free-coordinate eta
  design audit reports rank, column count, nullity, Person-facet connectivity,
  and null-space concentration. Rank-deficient fits remain visible but cannot
  produce a ready inference label or conditional bias table; downloadable
  evidence preserves the reason. The app does not silently substitute MML.
- JMLE Person tables now distinguish finite interior measures from
  all-minimum/all-maximum score boundaries. Boundary rows carry an explicit
  direction, readiness flag, reportable estimate, and value-role label; the
  Japanese and English convergence views explain why a finite optimizer value
  is not a finite Person MLE. The same audit is included in quick, full, demo,
  and publication-oriented downloads.
- A shared floating-point decision contract for fit and conditional bias
  screens. Classifications use unrounded values; exact MNSQ endpoints are
  explicit (`0.50` and `1.50` acceptable, `2.00` noisy), and numerical or
  3-decimal display boundaries are labelled instead of silently changing the
  decision. Fully unavailable fit statistics are not reported as stable. Both
  the quick bundle and full table download include row-level and summary
  stability audits.
- Descriptive anchor coverage now reports unique anchored levels, unanchored
  levels, and share without treating any percentage as a universal adequacy
  cutoff. Sparse simulation previews, fit scatter, bias heatmaps, and anchor
  coverage charts keep the relevant design caveats next to the visualization.
- FACETS-primary fit standardisation for RSM/PCM: Wright-Masters
  fourth-moment d.f., Wilson-Hilferty ZSTD, and the FACETS +/-9 cap are the
  default. Historical engine d.f./ZSTD remain in `*_ENGINE` sidecar columns.
  GPCM uses the same calculation only as an explicitly labelled
  FACETS-style approximation. See
  [`docs/fit_df_zstd_conventions.md`](docs/fit_df_zstd_conventions.md).
- Modular helper boundaries under `mfrm_app/` for CLI/doctor checks,
  export-frame collection, distribution metadata, privacy text, preflight
  contracts, and table IO.
- Guided Essential-view reading order and localized Japanese help popover
  overlays, so first-time users have clearer local guidance without hiding
  advanced evidence from full-view users.
- README-visible public verification commands: `--doctor`, `--release-check`,
  `--self-test`, and the Python-native demo report export.
- Stronger reproducibility and publication-export contracts for download
  frames, public-release readiness, runtime dependency boundaries, and
  privacy/caching expectations.

Earlier v0.2.0-beta (2026-04-17) shipped four major feature tracks; v0.2.1 and
v0.2.2 were post-release hotfixes landing the findings from a parallel UX audit.
Those foundational features are summarized below; see `CHANGELOG.md` for the
full per-version breakdown.

### New in v0.2.0

- **Publication Document downloads** — one click to produce a
  manuscript-ready Word (.docx), PDF, or HTML file with auto-generated
  abstract, exhaustive Methods, results tables, embedded figures, and an
  APA 7 reference list. Accessible from **Report → Exports**.
- Historical Posterior Viewer and Stan-generator code remains a legacy
  compatibility surface while its deprecation/removal decision is completed.
  It is outside the standalone product boundary and is not included in default
  downloads, demo archives, release evidence, the sidebar, Help, or Report
  routes.
- **Pattern-matched estimation errors** — instead of a generic "common
  causes" checklist, failures now surface a specific diagnosis + action
  (singular matrix → non-centered facet, maxit reached → raise to 1000,
  rating scale error → check Score column, …).
- **Pre-estimation readiness panel** — [OK] / [CAUTION] / [ISSUE] statuses across
  eight data-quality checks before the Run button fires.
- **Run history + two-run comparison** — up to 5 recent runs sit behind a
  Restore button (v0.2.1 lowered the cap from 10 → 5 for Community Cloud
  memory safety); pick any two to compare convergence, element-level
  Pearson r + RMSE, and reliability in a side-by-side panel.
- **Step-by-step progress** — `st.status` accordion reports Estimate →
  Diagnostics → Report tables → Bias interactions instead of a silent
  spinner, with a completion toast.
- **Surface-reason diagnostics** — residual PCA, bias interaction, and
  reliability now explain *why* they were skipped instead of returning
  a blank "not available".

### New in v0.2.1 (first UX-audit hotfix pass)

- **Download sample CSV** — inspect the built-in demo dataset's
  column structure before preparing your own upload.
- **Upload size preflight** — ≥50 MB warns, ≥100 MB errors, before
  pandas parsing consumes memory on Streamlit Community Cloud.
- **Clear-history confirmation** — the Clear history button now
  requires an explicit "Yes, delete all N snapshots" confirmation.
- **Downloads tab signpost** — a clear pointer to Publication Document
  so users can locate the Word / PDF export.
- **Config JSON import** — upload a previously downloaded config to
  inspect its settings as a table (read-only).
- **Accessibility CSS** — `:focus-visible` outline for keyboard users,
  `prefers-reduced-motion` support, narrow-screen (<900 px) layout.
- **Readiness status labels** — all readiness indicators use
  [OK] / [CAUTION] / [ISSUE] labels instead of colour-only signals.

### New in v0.2.2 (second UX-audit hotfix pass)

- **MFRM / Bayesian quick-reference glossary** — 20-term table
  (logit / Infit / Outfit / MnSq / ZSTD / step facet / Rhat / ESS / …)
  surfaced under Help → Glossary.
- **Run breadcrumb** — a persistent one-line banner above the results
  showing model / method / convergence / counts / fingerprint so users
  always know which run they are viewing after a history restore.
- **Measure table reorder** — Facet / Level / Estimate / SE / CI /
  Infit / Outfit columns are surfaced first; Anchor / N / metadata
  move to the tail.
- **Fit-cell colour highlighting** — Infit / Outfit mean-square cells
  colour-coded to the Wright & Linacre (1994) interpretation bands
  (0.5–1.5 acceptable · >1.5–2.0 noisy · >2.0 distorting ·
  < 0.5 over-fit).
- **Scree plot EV = 3 reference line** — previously only EV = 1 and
  EV = 2 were drawn; the "EV > 3 = strong secondary residual signal"
  threshold is now visible on the plot itself.
- **Help tab refresh** — Quick Start documents all v0.2.x features;
  Troubleshooting cites the new diagnostic helpers.

See `CHANGELOG.md` for the per-commit breakdown and
`docs/v0.2.0_beta_plan.md` for the original v0.2.0 release checklist.

## Preview

![MFRM Streamlit sample-data result overview](docs/images/app-data-overview.png)

The screenshot uses the built-in synthetic sample data and shows the
post-estimation workspace after the optional sample guide. On a new session,
the application first offers three routes: learn with the sample, start with
your own data, or continue without the guide. Raw response rows remain behind
a collapsed preview instead of being exposed by default. In Guided defaults,
the primary **Run this analysis / この設定で分析する** action now follows the
column mapping and readiness checks in the main workspace; the sidebar no
longer duplicates a technical `FACETS-mode` Run action. Advanced controls keep
the established sidebar route for expert workflows. The data source is also
grouped into four user tasks—built-in example, generate, paste, or upload.
Choosing an example keeps its scenario selector visible directly below, while
the unchanged internal IDs preserve saved-state and analysis compatibility.
Guided defaults also limits setup to first-run decisions: column roles, model,
estimation method, and analysis coverage. Weighting, anchors, optimizer and
identification controls, regularization, population models, score-scale
overrides, report scaling, and visualization styling remain available after
switching to Advanced controls. Returning to Guided restores documented
standard values instead of retaining an invisible advanced customization.

The current cross-cutting UX baseline and the planned workflow-shell work are
documented in [`docs/ux_adversarial_audit.md`](docs/ux_adversarial_audit.md).
The pre-browser accessibility gate is defined in
[`docs/browser_accessibility_acceptance.md`](docs/browser_accessibility_acceptance.md).
It generates 52 English/Japanese task cases and 482 required evidence rows;
templates, CSS inspection, and AppTest correctly remain `NOT_ACCEPTED` until
all version-matched browser records pass with retained evidence.

## Data Privacy

Rating data can contain person IDs, rater IDs, school or institution identifiers, subgroup labels, and other sensitive information.

For confidential data, run the app locally on a controlled machine:

```bash
streamlit run streamlit_app.py
```

If the app is deployed to Streamlit Community Cloud or another hosted platform, uploaded data may be processed on remote infrastructure. Confirm your institutional data-handling requirements before upload. Prefer removing direct identifiers, limiting uploaded columns to analysis variables, and avoiding hosted deployments for regulated or confidential datasets.

The app is intended to keep uploaded data in memory during the session unless the user explicitly downloads or exports outputs.

## Install

Use a virtual environment:

```bash
git clone https://github.com/Ryuya-dot-com/MFRM_Python_Application.git
cd MFRM_Python_Application
git checkout main
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

For local verification and CI-equivalent smoke tests, install the development dependencies:

```bash
python -m pip install -r requirements-dev.txt
```

## Run

```bash
streamlit run streamlit_app.py
```

For a first local smoke check, keep **Sample data (built-in)** selected (the
default Writing essay scenario is loaded automatically), leave the
guided defaults unchanged, click **Run FACETS-mode estimation**, then open the
**What should I look at first?**, **Data**, **Visuals**, and **Help** tabs.
Other built-in scenarios switchable from the **Data source** radio:
*Large-scale writing* (PCA-ready), *L2 speaking* (analytic-rubric /
PCM), *Clinical OSCE* (station-dominant), *Reading testlet — binary*
(0/1 scoring with item-text-person nesting). Each scenario's sidebar
**About this dataset** expander lists its APA 7 references.

The same **Data source** radio also includes **Generate synthetic data**.
Use it to create public long-format MFRM simulation data by setting the number
of facets, each facet's level count and spread, the score-category support,
zero-count categories, threshold spacing, missingness, and seed values. The
generated view includes reactive score-histogram, Wright-preview, and
pathway-preview plots, plus a generated CSV download.

The generator samples from an RSM-style many-facet adjacent-category model:
`log(P_k/P_{k-1}) = eta - tau_k`, where
`eta = person measure - sum(facet effects)`. This is the Andrich (1978)
rating-scale form extended to many-facet main effects in the sense of Linacre
(1989). The normal draws for effects, rotating first-facet assignment, row
deletion for missingness, and zero-count category recoding are simulation
conveniences for demonstrations and stress tests, not empirical defaults from
those sources. Row deletion is uniform MCAR deletion. Observation noise and
zero-count recoding deliberately depart from the pure RSM-style data-generating
model; the reported truth thresholds describe the pre-recoding sampling model.

For deployment notes, including Streamlit Community Cloud settings and privacy gates for hosted use, see `DEPLOYMENT.md`.

## Verify

```bash
python -m py_compile streamlit_app.py
python streamlit_app.py --doctor
python streamlit_app.py --release-check
python streamlit_app.py --self-test
make apptest
python streamlit_app.py --benchmark-quick --benchmark-csv validation/generated/benchmark_smoke.csv
python streamlit_app.py --export-demo-report validation/generated/demo_report
```

For public distribution, run at least `--doctor`, `--release-check`, and
`--self-test` from a fresh checkout before tagging or deploying. `--doctor`
confirms Python/package floors, bundled assets, the standalone runtime boundary,
and the privacy/cache boundary. `--self-test` exercises the app's built-in
statistical and export contracts. The export command above generates a
sanitized, Python-native demo-report package without uploading private data.
The legacy parity-fixture generator and frozen fixtures are compatibility-only
maintenance code. They have no CLI or Make target, are not called by
`--self-test`, and tests marked `legacy_compat` are excluded from `make
apptest`, `make verify`, and GitHub CI. They are outside the standalone release
gate.

Optional Make shortcuts:

```bash
make verify
make clean
make run
```

## Data Format

Use long-format rating data. The interactive upload path accepts CSV, TSV, TXT,
Excel `.xlsx` / `.xlsm` (first worksheet), Apache Parquet, and JSON / JSON-lines.
Pasted spreadsheet text supports comma, tab, or semicolon delimiters.

```csv
Person,Rater,Task,Criterion,Score
P01,R1,T1,C1,4
P01,R2,T1,C1,5
P02,R1,T2,C2,2
P02,R2,T2,C2,3
```

Required:

- Person column
- Score column with ordered integer categories
- At least two non-person facet columns such as rater, task, criterion, prompt, form, or occasion

Supported score handling includes:

- Intended rating boundaries, such as 1 to 5
- Zero-count intermediate categories, such as a 1-5 scale where category 3 is absent
- Optional recoding of non-consecutive observed categories
- Rejection of fractional categories with a clear error

### Paste From a Spreadsheet

Teachers and workshop users can paste directly from Excel, Google Sheets, or a
gradebook export. Keep the first row as column names, then select the cells and
copy/paste them into **Data source → Paste CSV/TSV text** in the sidebar.

```csv
Student,Rater,Assignment,Criterion,Score
S01,TeacherA,Essay1,Content,4
S01,TeacherA,Essay1,Organization,3
S01,TeacherB,Essay1,Content,5
S02,TeacherA,Essay1,Content,2
```

When mapping columns in the app, choose `Student` as **Person**, `Score` as
**Score**, and use `Rater`, `Assignment`, and `Criterion` as facets. Do not paste
total scores if rubric-level ratings are available; one row should represent one
rating event.

## Model Scope

Implemented in the standalone Python engine:

- RSM
- PCM
- bounded GPCM
- JMLE
- MML via EM, Direct, Hybrid, and Auto engine selection
- latent regression via constrained `population_formula`
- fixed user-set population prior SD for the current MML population model path
- on-demand MML prior-SD sensitivity refits with fit, measure-shift, rank, and
  latent-regression coefficient diagnostics, plus prespecified AnalysisID,
  SensitivityPlan, variant-record, and SensitivityDecision lineage
- MML observed-information covariance for non-person facet SE/CI when the
  fitted parameter vector is small enough for routine diagnostics
- MML covariance audit exports with Hessian rank, condition number,
  eigenvalue regularization, fallback/skip status, claim status, and the S0
  output-qualification reason code
- EAP posterior scoring
- plausible values
- strict marginal diagnostics
- residual PCA diagnostics
- residual PCA stability audit for sparse or weakly overlapped residual matrices,
  leave-one-column-out EV1/loading sensitivity, and row-bootstrap EV1/loading
  stability
- bias / local interaction screening with DFF-style sparse-cell, Holm, BH/FDR,
  and practical-logit review flags
- bias-inference audit exports that label conditionality, multiplicity family,
  sparse cells, profile-CI status, and connected-scale caveats
- versioned `output_qualification.csv` gates that keep unqualified information
  criteria, automatic model recommendations, LR decisions, structural
  covariance claims, and pairwise bias measures out of public conclusions while
  retaining explicitly marked technical values for audit
- exact identified step and GPCM log-slope optimizer coordinates, with a
  `parameterization_audit.csv` that preserves expanded public estimates while
  recording free dimensions, constraint residuals, coordinate rank, and slope
  bounds
- anchor audit and linking review
- anchor drift and equating-chain summaries
- anchor/equating workflow checklist for current-run linking evidence
- prediction for fitted, held-out, and scenario rows
- simulation and design evaluation
- final-report readiness checklist with legacy display columns and a canonical
  EvidenceRecord ledger that separates computation state from conclusion
  stability
- visual interpretation checklist for guided figure reading
- visual evidence binder with figure files, caption drafts, and figure-to-claim mapping
- visualization preferences for theme, label density, font size, dimensions, and caption detail
- claim-to-evidence matrix for manuscript and reviewer-response planning
- method-reference audit that maps model, estimation, fit, dimensionality, bias, simulation, and external-validation claims to APA/Zotero-aligned references
- latent-regression covariate type preview for numeric IDs/codes that may need categorical coding
- downloadable tables, configuration, scripts, method appendix, manuscript handoff, and manuscript binder
- reproducibility/config fingerprints for analysis exports

## Guided Reporting Workflow

To begin an article, open **Report & Export → Start a paper**. Download the
English APA Word template, fill in the study-specific prompts, and consult the
optional Results draft with its evidence checks. See the
[manuscript guide and reusable template](docs/apa_manuscript_template.md).

After fitting a model, inspect results in this order:

1. Output qualification: open `output_qualification.csv`; a
   `PublicConclusionAllowed = False` row overrides any raw numeric value shown
   elsewhere in the bundle.
2. Parameterization: open `parameterization_audit.csv`; require `PASS`, zero
   artificial coordinate null directions, and satisfied expanded GPCM bounds.
3. Convergence: do not interpret final measures until the optimizer converges.
4. SE/CI basis: review `SE_Method`, `SE_Status`, `CI_Method`, and
   `CI_Status`; for MML also review `mml_covariance_audit.csv` before
   using observed-information intervals in prose.
5. Category functioning: check sparse categories, monotonic average measures, and threshold order.
6. Reliability / separation: confirm whether the design supports stable person and facet ordering.
7. Wright map targeting: check whether person locations and facet difficulty/severity ranges overlap.
8. Fit diagnostics: review large standardized residuals and misfitting elements.
9. Bias / local interaction: use `bias_inference_audit.csv` with the DFF tables; treat flags as review prompts, not automatic proof of bias.
10. PCA / dimensionality: check both the residual PCA result and `pca_stability_audit.csv`, including missingness, pairwise overlap, leave-one-column-out, and bootstrap stability.
11. MML prior SD: for MML reports, state the fixed population prior SD and use the Fit Details prior-SD sensitivity screen, or justify why the fixed scale is part of the design, when population-scale claims matter.
12. Anchor / linking review: check connectedness and anchor stability before comparing runs or groups.
13. Strict marginal diagnostics: use for final MML reports when feasible.
14. Publication gate: check whether APA-style conclusions are ready, caveated, or blocked.
15. Submission action plan: fix prioritized blockers, caveats, boundaries, and wording repairs before manuscript use.
16. Case interpretation guidance: review common interpretation traps and safer wording repairs detected in the current run.
17. Final-report readiness: use the generated checklist before writing conclusions.
18. Manuscript claim guide: check what is safe to claim, what requires a caveat, and what should not be claimed yet.
19. Claim-to-evidence matrix: map each manuscript claim to exported tables, figures, diagnostics, caveats, reviewer questions, and archive files.
20. Method-reference audit: check which APA/Zotero-aligned references support each method surface before writing the literature-backed Methods and Limitations.
21. Manuscript template: adapt the generated Methods, Results, limitations, and reviewer preflight scaffold after resolving claim-guide cautions.
22. Visual evidence binder: review figure files, figure-to-claim links, caption drafts, and visual reviewer questions.
23. Manuscript handoff and binder: download the final-result guide, checklist, and curated writing packet for coauthor review or submission prep.

The final-report readiness checklist, publication gate, submission action plan,
first-read guide, manuscript template, and generated report text use the same main thresholds:
about 5% or fewer observation residuals with `|z| >= 2`, person reliability at
least 0.80 when person separation is the goal, and residual PCA first eigenvalue
below 2.0 for a clean residual-structure screen. Residual PCA is not a
standalone proof of unidimensionality and is not a DIMTEST/UNIDIM
implementation; report it with fit and local-dependence evidence.

The complete readiness download adds deterministic `AnalysisID`, `EvidenceID`,
`ComputationState`, `StabilityState`, `ReasonCode`, scope, and source-artifact
fields without changing the original five readiness columns. A displayed
`Review` therefore no longer conflates an executed caution, a failed fit, and
an analysis that was not assessable. Standard ZIPs also include reconstructable
AnalysisIdentity, EvidenceRecord, MML sensitivity decision, and settings JSON
sidecars; public mode uses non-identifying aggregate residual/linking evidence
instead of restoring excluded row-level tables.

For MML latent regression, inspect the covariate type preview before fitting.
Integer-like columns such as `GradeCode = 1, 2, 3` are flagged so you can decide
whether they are continuous predictors or category labels that should be forced
categorical.

The app and demo export also include a manuscript claim guide, claim-to-evidence
matrix, method-reference audit, visual evidence binder, final-result handoff,
manuscript binder, public-beta limitations, and release readiness tables. The submission action
plan combines these sources into a prioritized first-read table so public claims
stay aligned with what the standalone Python engine currently supports.
In the app UI, wide reporting tables show the most important columns first,
wrap short guide tables for reading, and place full-detail tables in expanders;
downloads still contain the complete columns. On desktop screens, long result-tab
bars wrap instead of forcing users to rely on horizontal scrolling.
On long result pages, the authoritative Essential or All-panels selector stays
at the top of the viewport while content scrolls. Narrow screens keep the
Essential selector on one touch-scrollable line so the navigation does not
turn into a tall overlay. The app does not add custom keyboard shortcuts.
For a new session, the optional sample guide uses five focused steps: choose a
route, check the sample's Person/Score/facet roles, run the production
estimator, distinguish a bounded interpretation from an overclaim, and review
a supported/limited/next-action checkpoint. It is skippable, resumable, and
restartable within the session. Guide completion records learning progress
only; it does not validate a data set, design, claim, or intended use.
Default table, manuscript, and demo archives now contain only Python-native
analysis, diagnostics, evidence contracts, figures, and reproduction assets.
Legacy cross-package inventories, R/Julia scripts, Stan/Posterior handoff
packages, and mfrmr migration maps remain compatibility-maintenance surfaces;
they are not included in default downloads.

The Visuals tab also includes a downloadable visual interpretation checklist.
It maps each figure to the first signal to read, the review trigger, and the
recommended next action for guided review.
It also includes a visual method evidence table that links each plot family to
its Rasch/MFRM diagnostic role and explains the app's readability rules.
Figure exports use a manuscript profile: white background, consistent font,
compact margins, 300 DPI PNG when static export is available, and matching
interactive HTML for inspection. Users can also set the figure theme, label
density, base font size, static width, minimum height, and caption-detail level
from the sidebar. These choices are recorded in `visualization_settings.json`
and `visualization_settings.csv` and propagated into `figure_manifest.csv`,
`visual_evidence_map.csv`, caption
drafts, and the visual QA preflight table. The figure bundle includes
`figure_manifest.csv` with the recommended manuscript use and reporting caution
for each figure.
The visual evidence binder adds `visual_evidence_map.csv`,
`visual_caption_drafts.md`, `visual_qa_preflight.csv`, the figure files,
visual interpretation checklist, method evidence table, and claim-to-evidence
matrix in one review packet.
For PCM and bounded GPCM, category probability curves can be read either as an
averaged overview or for each selected step-facet level. The downloadable table
bundle also includes long-form curve data for all available curve scopes.

To generate a synthetic guided report without uploading data:

```bash
python streamlit_app.py --export-demo-report validation/generated/demo_report
```

Open `validation/generated/demo_report/manuscript_handoff.md` first, then read
`manuscript_handoff_checklist.csv`, `claim_to_evidence_matrix.csv`,
`apa_report_sentence_audit.csv`,
`method_reference_audit.csv`,
`visual_evidence_map.csv`, `visual_qa_preflight.csv`,
`visualization_settings.json`, `visualization_settings.csv`,
`MFRM_Demo_Visual_Evidence_Binder.zip`,
`MFRM_Demo_Manuscript_Binder.zip`, `MFRM_Demo_Report.html`,
`export_privacy_manifest.csv`,
`final_report_readiness_analysis_identity.json`,
`final_report_readiness_evidence_contract.json`,
`publication_gate_summary.csv`, `submission_action_plan.csv`,
`case_interpretation_guidance.csv`,
`final_report_readiness.csv`,
`manuscript_claim_guide.csv`,
`manuscript_template.md`,
`visual_interpretation_checklist.csv`, `visual_method_evidence.csv`,
`public_beta_limitations.csv`, `public_release_readiness.csv`,
`mfrm_app_engine_runner.py`, `mfrm_local_batch_workflow.md`,
`mfrm_jmle_self_contained.py`,
`visual_caption_drafts.md`, `figure_manifest.csv`,
`MFRM_Demo_Publication_Figures.zip`, and the interactive diagnostic figures in
`figures_html/`.

## Statistical Caveats

- GPCM is not a strict Rasch model. Its slope parameters change the interpretation of invariance and should be reported explicitly.
- For RSM/PCM, primary Infit/Outfit ZSTD now uses the FACETS/Wright-Masters fourth-moment d.f., Wilson-Hilferty transform, and an absolute cap of 9. The `DF_*_ENGINE` and `*ZSTD_ENGINE` sidecars preserve the prior homogeneous-variance convention for sensitivity checks. For GPCM, FACETS-primary columns are labelled `facets_style_approximation_for_gpcm`; do not claim exact FACETS equivalence.
- `KParams` now counts exact identified optimizer coordinates, and `parameterization_audit.csv` records the expanded/free dimensions, sum-zero residuals, coordinate rank, and bounded GPCM log-slope check. AIC/AICc/BIC and model-choice LR quantities nevertheless remain technical audit values until the G2 comparison-scope and calibration gates clear. Automatic recommendations remain withheld; GPCM-involving LR pairs must not use an ordinary chi-square decision.
- Measure SE/CI columns now carry method/status metadata. Non-person MML facet SEs use observed-information delta-method covariance when available; otherwise conditional information approximations are labelled as such. MML person SEs are EAP posterior SDs, not structural fixed-effect ML SEs. Rank-deficient or regularized MML covariance is `WITHHELD`, and even full numerical rank remains `TECHNICAL_ONLY` until interval coverage is validated. Archive `mml_covariance_audit.csv` with its reason code.
- ADEMP parameter-recovery exports include an explicit SE/CI coverage diagnostic table with `SEBasisRisk` and `CoverageClaimStatus`. Treat it as design-specific Monte Carlo evidence: cite the generator, fit method, replicate count, seed, and SE/CI status before making interval-calibration claims.
- The current latent regression path uses the app's documented fixed population prior SD behavior. Use the versioned MML prior-SD sensitivity decision before treating population-scale or latent-regression coefficients as robust; do not generalize beyond the implemented quadrature, variance treatment, constraints, and tested SD grid.
- Bias and differential functioning outputs are conditional screening tools. Report the exported bias-inference audit with DFF tables; do not make no-bias or confirmatory bias claims unless linking, common-scale evidence, sample size, multiplicity review, and precision support that scope. Pairwise local-measure output is withheld until its direction and contrast-SE defects clear G3.
- Residual PCA is a sparse-matrix diagnostic. Use the exported PCA stability audit before making dimensionality claims from eigenvalues or loadings; review leave-one-column-out and bootstrap sensitivity when the first residual component is near EV = 2 or EV = 3.
- Cross-package equality is not expected by default because FACETS, TAM, sirt, mirt, and this app can use different parameterizations, constraints, latent variance handling, and optimization details.
- Treat external R cross-check outputs as validation evidence, not as a runtime dependency.
- Reproducibility fingerprints help confirm that a report came from the same data/settings, but they are not a privacy guarantee or encryption method.

## Caching

- `st.cache_resource` is used only for the read-only core function namespace.
- `st.cache_data` is used for built-in sample data, bundled anchor templates/guidelines, and short-lived export bytes keyed by reproducibility fingerprints.
- Uploaded or pasted rating data are not cached as a resource. Export-byte caches are bounded and short-lived, but confidential analyses should still be run locally.

## Figure Exports

Interactive HTML figure exports are always attempted. Static PNG exports use
Plotly/Kaleido and require Chrome or Chromium in the runtime when using
`kaleido >= 1.0`; if that browser dependency is unavailable, the app falls back
to interactive HTML figures instead of blocking the analysis.

## Methodological Reference Roles (Not Integrations)

- TAM: faceted MML design, latent regression, EAP, and multifacet reference checks.
- mirt: GPCM, EAP/factor scores, plausible values, and broader IRT diagnostics.
- sirt: rater-facet and hierarchical rater model reference checks.
- mfrmr: historical functional-capability reference.

For the frozen historical comparison matrix, see `validation/README.md`; it is
not a current app workflow or release gate. For the archived optional R smoke status, see
`validation/R_CROSSCHECK_STATUS.md`.
For the archived external numerical-validation inventory, see
`validation/SIMULATION_REFERENCE_STATUS.md`.
These archived compatibility materials are not default downloads, demo
artifacts, release-check inputs, or release evidence for the standalone Python
application.

## Continuous Integration

The app-specific workflow is in:

```text
.github/workflows/python-streamlit.yml
```

It runs:

- dependency installation
- `python -m py_compile streamlit_app.py`
- `python streamlit_app.py --doctor`
- `python streamlit_app.py --release-check`
- `python streamlit_app.py --self-test`
- `python streamlit_app.py --benchmark-quick --benchmark-csv validation/generated/benchmark_smoke.csv`
- Streamlit AppTest smoke check
- demo report export smoke check

If this directory is used as a standalone GitHub repository, the workflow will be discovered normally. If it remains a subdirectory inside a larger repository, copy or mirror the workflow into the repository root `.github/workflows/` directory.

## License

MIT License. See `LICENSE` and `LICENSE_NOTICE.md`.

Commercial use is permitted, including paid teaching, consulting, internal
training, hosted demonstrations, and product evaluation. The app and
documentation are provided as-is, without warranty. Users are responsible for
data privacy, model assumptions, interpretation, external validation, and
institutional approval before relying on outputs in operational or high-stakes
settings.

We intentionally do not use CC BY-NC because the NonCommercial restriction would
conflict with the intended permission for commercial use. If documentation or
teaching excerpts are reused separately, keep attribution and do not imply author
endorsement.

## Publication Note

This repository is the publication target for the standalone Streamlit app.
Keep generated validation outputs, private datasets, and machine-specific paths
out of the public branch unless they have been explicitly sanitized.

## Files

```text
Makefile
streamlit_app.py
requirements.txt
requirements-dev.txt
LICENSE
README.md
DEPLOYMENT.md
CHANGELOG.md
CONTRIBUTING.md
SECURITY.md
MFRM_STREAMLIT_RELEASE_PLAN.md
RELEASE_CHECKLIST.md
anchor_templates_and_guideline/
locales/
mfrm_app/
tests/
.streamlit/config.toml
.github/workflows/python-streamlit.yml
.github/ISSUE_TEMPLATE/
docs/images/app-data-overview.png
```

Generated files such as `__pycache__`, benchmark CSVs, local logs, and validation output directories are ignored.

## Known assignment-mechanism validation

Repository-only validation now includes an exact small-state oracle for a
degree-conditioned informative Person-Rater assignment mechanism. The
Metropolis 2-switch sampler and its score-free row materializer passed every
frozen oracle, margin, connectivity, numerical, and outcome-blind gate for
`gamma=-0.8, 0, +0.8`; see
`validation/KNOWN_ASSIGNMENT_MECHANISM_PILOT_20260811.md`.

This does not expose an inferred missingness model in Streamlit and does not
qualify any estimator under informative assignment. It establishes a clean
simulation layer for the next response-generation study. FACETS 4.5 remains a
same-estimand external check for JMLE only; MML and exact CMLE remain explicitly
different-estimand sensitivities. Rounded FACETS fit display values are never
promoted to raw calculation inputs.

The subsequent fixed-Person, 10-replicate response screening is now complete;
see `validation/KNOWN_ASSIGNMENT_RESPONSE_SCREENING10_20260811.md`. All 30
FACETS/Python JMLE pairs and all 90 MML/CMLE fits qualified. Under the frozen
PCM, only the fixed/free normal-person MML Rater RMSE and MAE screening
intervals excluded zero for both `gamma=-0.8` and `gamma=+0.8`; JMLE, exact
CMLE, Task, Criterion, and threshold intervals did not. This is a candidate
assignment-sensitivity pattern conditional on one fixed Person vector, not a
confirmatory result or estimator ranking, and remains repository-only.
