# Validation Notes

This directory primarily preserves historical, optional compatibility records.
It is not part of the standalone Python product journey or default downloads,
and the Streamlit app does not call an external estimation engine.

## FACETS 4.5.0 Windows launch diagnosis (2026-08-11)

The earlier conclusion that FACETS 4.5.0 was generally blocked on this host is
superseded.  Direct and Python hidden-window routes reproduced `0xC0000005` in
`XojoGUIFramework64.dll`, but interactive FACETS and visible `BATCH=NO`
completed.  The qualified Python adapter now keeps the analysis window visible,
waits for a stable nonempty report and all expected Scorefiles, and posts
`WM_CLOSE`; forced termination is never accepted as success.  A later auxiliary
pass stalled because its generated Scorefile path was exactly 260 characters.
The shared invocation code now fails before launch above a qualified
220-character path budget.  See `FACETS_450_NATIVE_CRASH_DIAGNOSIS_20260811.md`,
`facets_shared_helper_qualification_20260811/result.json`, and
`fx2_20260811/result.json`.  The retained dumps remain useful evidence about
the hidden-window path, but are not evidence that visible FACETS 4.5.0 cannot
run.

## FACETS 4.5 PCM/JMLE qualification (2026-08-11)

The registered four-category PCM/JMLE scope with `Criterion` as the sole step
facet is now externally matched to FACETS 4.5.0.  The 20-replicate paired
pipeline pilot completed 40/40 directly eligible runs: all 320 FACETS Table 8
scale/category counts matched the retained inputs, 360 main-facet parameters
had weighted MAE 0.001230 logits (maximum 0.003059), and 240 criterion-specific
thresholds had weighted MAE 0.003128 logits (maximum 0.006815).  All expanded
PCM designs had nullity zero and all 40 paired heterogeneity-direction checks
passed.

This qualifies only the locked complete, balanced, common-category-support
scope under `Models=?,?,?,#,R3`.  It does not qualify named General/Specific
scales, mixed category supports, recoding, threshold anchors, sparse/missing
designs, more than one `#` facet, or raw fit parity.  FACETS measure/threshold
comparisons use a separate `Umean=0,1,6` primary pass; two-decimal auxiliary
output is never used to invent unreported raw fit precision.  See
`FACETS_PCM_QUALIFICATION_20260811.md`, the three preregistered PCM plan JSON
files, `facets_pcm_syntax_probe_20260811/`, and
`facets_pcm_pipeline_pilot20_20260811/`.

The follow-on missingness/misspecification boundary pilot is retained under
`facets_pcm_boundary_pilot_20260811/`.  All 16 rank-full PCM/RSM fit contracts
remained directly matched, and 144/144 Table 8 counts agreed.  The deliberately
disconnected Person-Rater control was more informative: the independent design
and the application both found nullity one and withheld inference in 8/8 fits,
while FACETS converged and printed `Subset connection O.K.` in all eight.  The
FACETS wording is therefore retained as product output, not promoted to a full
constrained-design identification guarantee.  See
`FACETS_PCM_BOUNDARY_PILOT_20260811.md` and
`facets_pcm_boundary_pilot_plan_20260811.json`.

## Same-row JMLE/MML/CMLE estimand bridge (2026-08-11)

The eight rank-full complete/planned-connected boundary runs now feed an
explicit same-row estimand bridge under `estimand_bridge_pilot_20260811/`.
It reuses 32 hash-bound FACETS/Python JMLE contracts and adds 32 fixed/free-SD
Q31 MML fits plus 16 native exact-CMLE fits.  All 48 new fits were inference
ready, all constraints passed at a maximum residual of `1.11e-16`, and all 16
within-likelihood-basis RSM-versus-PCM direction checks passed.  Recovery is
retained for 720 main-facet and 360 threshold rows.

The bridge does not rank estimators: JMLE uses fixed Persons, MML assumes a
normal Person population, and CMLE conditions on Person totals.  Joint,
marginal, and conditional likelihood/AIC values are kept in separate bases.
At two-replicate depth, the larger estimator shifts and free-SD changes under
planned missingness are hypotheses for a repeated study, not bias or RMSE
ordering claims.  See `ESTIMAND_BRIDGE_PILOT_20260811.md` and
`estimand_bridge_pilot_plan_20260811.json`.

## Person-population shape preflight (2026-08-11)

The next bridge layer fixes every realized Person vector to mean zero and
population SD 0.8 while varying only its shape across normal, right-skewed,
symmetric-mixture, and heavy-tailed conditions.  Each shape is crossed with
complete and ring-connected planned-missing designs under the correct
criterion-specific PCM.  The qualified v3 preflight completed 32/32 attempt
units: eight FACETS/Python JMLE pairs, 16 fixed/free-SD Q31 MML fits, and eight
native exact-CMLE fits.  Every fit and constraint gate passed; the worst
FACETS/Python main-effect difference was 0.010392 logits and the worst PCM
threshold difference was 0.006188 logits.

The worker is immutable-input, modulo-sharded, and attempt-resumable.  A replay
skipped all 32 completed attempts, performed zero refits, and changed zero
completion-marker hashes.  Two result-blind infrastructure amendments retain
the discovered Windows/Dropbox directory-finalization and legacy path-length
boundaries.  A later mechanical aggregation amendment joins condition metadata
from the retained manifest; only `aggregate_corrected/` is qualified, while
the original aggregate remains superseded diagnostic evidence.

Free-SD values of 0.687--0.811 and the single extreme-Person case are one-seed
hypotheses only.  The registered 20-replicate extension remains screening
depth below the repository's 100-eligible-replicate operating-characteristic
threshold.  See `ESTIMAND_DISTRIBUTION_PREFLIGHT_20260811.md`,
`estimand_distribution_screening_plan_20260811.json`, and
`estimand_distribution_preflight_v3_20260811/aggregate_corrected/`.

## Person-population shape screening, 20 replicates (2026-08-11)

The preregistered extension under `estimand_distribution_screening20_20260811/`
completed 160 datasets and 640/640 successful attempt units.  FACETS/Python
JMLE passed direct agreement in 160/160 pairs; 320/320 MML and 160/160 exact-
CMLE fits were inference-ready.  An eight-shard run took about 26.5 minutes,
and a full resume replay skipped all 640 attempts without changing a completion
marker.

At matched realized Person mean/SD, no global non-normal-versus-normal RMSE
screening interval excluded zero for main facets (0/90) or PCM thresholds
(0/30).  Planned missingness was much more consequential: Rater and Task RMSE
increased in 20/20 estimator-by-shape contrasts and threshold RMSE increased in
20/20.  Free-SD MML nevertheless showed shape sensitivity, especially for
right-skewed and heavy-tailed Persons.  Local ledgers also found a symmetric-
mixture C02 threshold pattern shared by all five modes despite stable aggregate
threshold RMSE.  This is why centered facet-average signed errors are not used
as bias evidence; level- and step-specific errors are retained.

The screen does not rank estimators and remains below confirmatory depth.  Any
follow-up must use 100 fresh replicates rather than reuse the 20 replicates from
which contrasts were selected.  See
`ESTIMAND_DISTRIBUTION_SCREENING20_20260811.md` and
`estimand_distribution_screening20_20260811/screening_analysis_v2/`.

## Person-population shape confirmation, 100 fresh replicates (2026-08-11)

The fixed-N confirmatory extension uses only replicates 21--120 and retains 800
datasets plus 3,200 hash-bound attempt units.  The eight-shard run completed all
markers in 135.2 minutes; 3,199 attempts succeeded.  A full resume audit skipped
3,200/3,200 valid markers, performed zero refits, and preserved the combined
marker SHA-256.  MML was operational in 1,600/1,600 fits and exact CMLE in
800/800 fits.

Two registered free-SD MML effects reproduced after Holm adjustment.  In the
planned-connected design, right-skewed minus normal Person shape was `-0.08983`
(95% CI `[-0.10380, -0.07585]`) and heavy-tailed minus normal was `-0.06502`
(`[-0.08532, -0.04472]`).  The latter direction was confirmed but narrowly
missed its separate precision target.  By contrast, the symmetric-mixture C02
step-2 screening signal failed to reproduce in Python JMLE, free-SD MML, and
exact CMLE; all three point estimates changed to small positive values.

The registered JMLE threshold-RMSE design effect remained positive and precise,
but one FACETS exit-0/missing-report event left 99 rather than 100 pairs.  It is
therefore confirmatorily inconclusive, and strict batch/workbench qualification
is false.  Conditional on report production, numerical calibration remained
strong: 799/799 FACETS/Python pairs passed direct agreement and all within-facet
rank correlations were 1.0.  Three isolated post-analysis diagnostic replays
succeeded without replacing the failed evidence row, supporting a transient
report-I/O classification but not proving concurrency as its cause.

Base-R 4.5.1 independently reproduced all paired statistics, intervals, Holm
adjustments, and logical gates with maximum numerical difference `1.78e-15`.
See `ESTIMAND_DISTRIBUTION_CONFIRMATORY100_20260811.md`,
`estimand_distribution_confirmatory_plan_20260811.json`,
`estimand_distribution_confirmatory100_assessment_20260811.json`, and
`estimand_distribution_confirmatory100_20260811/confirmatory_analysis/`.

## Independent H6 observation-design replication, 100 fresh replicates (2026-08-11)

The original 99-pair H6 result remains confirmatorily inconclusive.  A separate
preregistered fixed-N replication used only replicates 121--220 and did not
pool, replace, or impute any earlier row.  All 100 native-precision Python JMLE
pairs were finite.  Planned-connected minus complete threshold RMSE was
`+0.075464` logits (95% CI `[+0.065075, +0.085853]`, one-sided
`p=2.36e-26`), and its interval half-width `0.010389` passed the registered
`0.015` precision target.  The effect therefore replicated for this normal-
Person, heterogeneous-PCM DGM; it is not a universal missing-data-bias claim.

The new resilient execution layer separately passed 200/200 Python evidence
and 200/200 FACETS 4.5.0 calibration gates with zero retries.  Worst-run direct
differences were `0.012475` logits for main effects and `0.010601` for PCM
thresholds; minimum within-facet Spearman was 1.0.  A 200/200 resume replay
performed zero refits.  Base-R 4.5.1 rebuilt the endpoint from unrounded Python
threshold rows and agreed within `1.78e-15`.  The retained v1--v3 preparations
document two result-blind preexecution-audit corrections and the legacy FACETS
path-length failure that motivated the 220-character worst-case path gate.

See `H6_DESIGN_REPLICATION100_20260811.md`,
`h6_design_replication_plan_20260811.json`,
`h6_design_replication100_assessment_20260811.json`, and
`h6_design_replication100_v4_20260811/aggregate/`.

## Fixed-density informative-assignment confirmation (2026-08-11)

The next registered layer separates response density from assignment
dependence.  Both sparse designs retain 960 rows, two Raters per Person, 40
Persons and 240 rows per Rater, one Person-Rater component, constrained PCM
nullity zero, and the same generated complete responses.  The stress design
alone aligns ascending latent-ability quartiles with ascending mean Rater
severity.  This is a deliberately constructed mechanism, not an empirical
MAR/MNAR diagnosis.

After a fixed 20-replicate screen, the fresh confirmation used only replicates
241--340.  All 800/800 attempt units were evidence-ready: 200/200 resilient
FACETS/Python JMLE pairs, 400/400 fixed/free-SD Q31 MML fits, and 200/200 exact-
CMLE fits.  The primary free-SD MML Rater-RMSE contrast (aligned minus planned)
was `+0.098634` logits (95% CI `[+0.084322, +0.112946]`, one-sided
`p=7.70e-25`); its half-width `0.014312` passed the registered `0.015` target.
All four secondary directions passed the primary gate and Holm family-wise
control: fixed-SD Rater RMSE `+0.068506`, free-MML population SD `-0.118394`,
and free/fixed Rater-compression slopes `-0.489320`/`-0.383835`.

Fixing the correct Person SD reduced but did not remove the Rater-recovery
effect.  JMLE and exact-CMLE Rater-RMSE contrasts were small descriptive
context only; they do not establish estimator superiority.  FACETS remained a
same-estimand JMLE calibrator (200/200 ready, zero retries, minimum within-facet
Spearman 1.0), not a gold standard for marginal or conditional estimands.  R
4.5.1 independently rebuilt all 500 endpoint contrasts, t statistics,
intervals, Holm adjustments, and logical gates within `3.55e-15` of Python.

See `INFORMATIVE_ASSIGNMENT_SCREENING20_20260811.md`,
`INFORMATIVE_ASSIGNMENT_CONFIRMATORY100_20260811.md`,
`informative_assignment_confirmatory_plan_20260811.json`, and
`informative_assignment_confirmatory100_20260811/aggregate/`.

### Product integration boundary

The Streamlit result path now exposes this evidence without converting it into
an empirical MNAR classifier.  `mfrm_app/design_assignment.py` builds an
outcome-blind audit from unique Person-Rater assignments, reports exposure,
density, co-rating, and the direct Rater-overlap graph, and always labels the
assignment mechanism as unidentified from observed assignments.  Score,
estimated ability, estimated severity, residuals, and fit statistics are not
inputs to that audit.  A separate assignment sensitivity contract is shown
as eligible/not run.  For strict eligible JMLE/MML RSM/PCM fits, the
Streamlit Prediction/Simulation path first uses a paired parametric V1: it
constructs a score-free connected degree-preserving 2-switch counterfactual,
generates responses with common row-level uniforms under each assignment, and
refits only the current method.  Contrasts require both refits to be
inference-ready and retain explicit fit/failure denominators.  It also maps
requested alignment doses to achievable points on the discrete greedy switch
path and reports each dose minus the observed assignment using the same random
draw within replicate.  Dose is normalized path-objective progress, not an
estimated propensity or a global worst-case guarantee.  The external
fresh-100 effects remain rationale for the workflow, not correction factors
for a user's data.

When Person-Rater blocks have unequal Task/Criterion signatures, the strict
path remains blocked and the product audits a separate binary MILP fallback in
`mfrm_app/assignment_context_milp.py`.  Each source block remains indivisible;
Person degree, unique Person-target-Rater assignment, per-Rater Person and row
exposure, and every Rater-by-context-cell row count are constraints.  A
deterministic spanning tree of the observed direct-overlap graph supplies
locked witness-Person constraints, so connectivity is guaranteed inside the
optimization.  The observed assignment must reproduce the complete constraint
vector before HiGHS is called.  An endpoint is exposed only for a certified
global optimum with zero reported MIP gap, negligible constraint/integrality
residuals, a changed edge set, and positive direction-adjusted gain.  This is a
constrained endpoint—not a min-cost-flow dose curve, unconstrained worst case,
or assignment-mechanism identification.

MML now has two explicitly separate Person-response generators.  The default
uses the source posterior/EAP measures and labels them as reference coordinates,
not truth.  The optional `mml_population_rank_preserving` mode draws iid normal
values at the source fit's fixed or estimated population SD, sorts them, and
maps the order statistics to the strictly ordered source EAP ranks.  One draw is
shared by every assignment scenario within a replicate, independently of the
common response-uniform stream.  This preserves the fitted rank-to-assignment
relation needed by the retained assignment contrast while varying population
magnitudes.  It is deliberately labelled conditional on that ordering and not
an unconditional new-person sample, assignment-mechanism estimate, or bias
correction.  Tied EAP ranks, JMLE, invalid SD, and latent regression fail closed.
Individual generated coordinates remain private; aggregate generation summaries
and the generator contract can enter a public export.

The same result path now places an estimand contract before technical result
tables.  It distinguishes fixed-Person joint JMLE, Gaussian-population marginal
MML, and Person-total conditional exact CMLE, prohibits cross-basis
likelihood/AIC/BIC ranking, and states that FACETS 4.5.0 is an external
same-estimand comparator only for qualified JMLE scopes.  Its FACETS precision
boundary also prohibits recomputing or adjudicating fit from two-decimal
display values.  The estimand contract, assignment audit, sensitivity plan, and
validation evidence register, selected-runner gates, block profiles, and the
fallback solver/context-margin/witness audits are included in the one-click
CSV/Excel bundles.  Completed session evidence is included in the matching
complete/private Downloads bundle; public filtering removes Person-level
assignment maps, block profiles, switch ledgers containing Person IDs, and
connectivity witnesses.  See
`tests/test_design_assignment.py`, `tests/test_design_assignment_integration.py`,
`tests/test_assignment_sensitivity.py`,
`tests/test_assignment_context_milp.py`,
`tests/test_assignment_generator.py`,
`tests/test_assignment_sensitivity_integration.py`, and
`tests/test_assignment_sensitivity_ui.py`.  The product decision and current
V1 boundary are summarized in
`DESIGN_ASSIGNMENT_PRODUCT_INTEGRATION_20260811.md`.

## Current repository-only stress pilot

The 2026-08-09 Python–R stress pilot is a new, reproducible repository-only
exception to the otherwise archived material in this directory:

- protocol: `CROSS_ENGINE_STRESS_PROTOCOL_20260809.md`;
- results and retained CSV evidence: `cross_engine_stress_20260809/`;
- standalone Python generator/fitter/summarizer: `cross_engine_stress.py`;
- matched mfrmr/TAM/immer/sirt runner: `cross_engine_stress.R`.

It compares RSM/PCM cumulative-difficulty surfaces under explicit estimator,
quadrature, category-map, and identification labels. It is a deterministic
feasibility and failure-mode pilot, not a CI release gate, coverage study,
package ranking, GPCM-equivalence claim, or sample-size rule.

## Decision stability and sparse-design stress

The same date also has a Python-only numerical/presentation stress matrix:

- runner: `decision_stability_stress.py`;
- retained evidence and figures: `decision_stability_20260809/`; and
- core threshold contract: `../mfrm_app/decision_stability.py`.

It exercises floating-point neighbours and display-rounding bands around fit
MNSQ 0.50/1.50/2.00, bias alpha=.05 and |bias|=.50, nine balanced/planned-
missing/sparse/zero-category designs, and anchor coverage from 0% to 100%.
It is a deterministic decision-sensitivity and failure-mode check. It is not
a power, coverage, sample-size, universal anchor-share, or model-validity
claim. Reproduce it with:

```bash
MPLCONFIGDIR=/tmp/mfrm_matplotlib_cache \
  python3 validation/decision_stability_stress.py
```

## Repeated-simulation operating-characteristics lane

The first Python JMLE repeated-simulation scaffold is governed by
`OPERATING_CHARACTERISTICS_PROTOCOL_20260809.md` and reproduced with
`operating_characteristics_pilot.py`. Its retained smoke artifacts are under
`operating_characteristics_20260809/`, including attempted-run/failure
accounting, decision denominators, Monte Carlo uncertainty, recovery and
coverage by identification scale, raw/display conclusion sensitivity, a
UI-ready first-read table, and three figures.

The null/alternative pairs use recorded common random numbers. Correct and
+0.25-logit contaminated anchor conditions also share generated data so the
anchor sensitivity is paired. The sparse design is an intentional ineligible-
focal-cell control. The two-replicate smoke profile is orchestration evidence
only; it is not a false-positive, power, coverage, sample-size, package-parity,
or anchor-share claim. Public UI exposure remains withheld pending the frozen
study-depth and remaining matched-estimand gates in the protocol.

The prospective precision/stopping contract is now frozen in
`operating_characteristics_precision_plan_20260809.json`. The runner validates
and copies that exact file into each output bundle, records code and plan
hashes, and can verify that a larger manifest preserves every prior RunId,
condition definition, and seed. The separate 20-replicate output is retained
under `operating_characteristics_pilot20_20260809/`; its parent-manifest audit
confirms that all 16 smoke rows are unchanged and only 144 consecutive rows
were added.

The pilot is a stop signal, not a study-depth result. All 160 fits returned in
82.1 recorded fit-seconds, 156 carried the optimizer convergence flag, and 124
met the registered app bias-decision eligibility rule. However, reconstructed
terminal-gradient sup norms ranged from `0.00013306` to `0.102327`, so 0/160
passed the separately labelled post-pilot `1e-4` numerical-readiness
sensitivity. The two sparse conditions each yielded only 4/20 eligible focal
decisions; Wilson-lower-bound inflation would require 6,200 attempts per
condition, beyond the registered 2,500 cap. Fourteen fit statistics fell in a
3-decimal display boundary and five would be classified differently if the
rounded display, rather than the raw value, drove the label. Paired hard-anchor
contamination moved fixed Raters by exactly +0.25 logits on average while
unanchored Raters averaged zero shift. See
`operating_characteristics_pilot20_20260809/PILOT20_ASSESSMENT.md` and its four
review figures. Every pilot first-read public-surface flag remains false.

### Strict JMLE numerical qualification and structural gate

The prospectively frozen follow-up did not retroactively alter the pilot. Stage
A reproduced 16 pilot fits and confirmed the analytical gradient against both
the optimizer Jacobian and finite differences, while all 16 original terminal
vectors still failed the raw `1e-4` gradient gate. Stage A2 compared two
predeclared continuation candidates. Its priority L-BFGS-B precision candidate
passed all 16 numerical contracts; BFGS reached acceptable gradients but
reported precision loss in 14/16 fits and made a much larger flat-direction
move, so it was not selected. Stage B2 applied only the selected candidate to
the frozen 160-run manifest. All 160 passed the raw `1e-4` gate and all other
frozen numerical checks; 126 passed at `1e-5` and 43 at `1e-6`. No run changed
its gate decision when the raw gradient was displayed to three significant
digits.

Numerical success was not treated as inferential success. A separate movement
audit found a maximum free-coordinate change of 57.02 logits despite negligible
likelihood gains. The exact free-coordinate eta-design audit localized the
problem: every one of the 40 sparse runs had rank deficiency/nullity 7 and
eight disconnected Person-Rater components. Null-space energy was concentrated
in Person (72.69%) and Rater (27.31%) coordinates. Balanced and anchor
conditions had nullity zero and a connected Person-Rater graph. These results
are retained respectively under:

- `operating_characteristics_strict_jmle_smoke_20260809/`;
- `operating_characteristics_strict_jmle_a2_smoke_20260809/`;
- `operating_characteristics_strict_jmle_b2_20260809/`;
- `operating_characteristics_jmle_movement_20260809/`; and
- `operating_characteristics_identifiability_20260809/`.

The application now computes the same eta-rank audit before JMLE optimization.
It exposes `Converged` and `InferenceReady` separately, withholds conditional
bias for a rank-deficient eta design, and exports summary, connectivity,
null-space-energy, and coordinate-weight tables. The post-change integration
bridge in `jmle_identifiability_integration_20260809/` matched the frozen rank,
nullity, connectivity, scope, and readiness contract for 160/160 RunIds. It
changed neither optimizer controls nor estimates and does not automatically
substitute MML. The guard is deliberately scoped to the Person/facet eta
blocks; it does not by itself qualify PCM/GPCM steps or GPCM slopes.

The old Stage A/A2/B2 adapters retain pre-change source hashes and intentionally
refuse a post-change rerun. The supported bridge from frozen evidence to the
current core is reproduced with:

```bash
MPLCONFIGDIR=/tmp/mfrm_matplotlib_cache \
  python3 validation/jmle_identifiability_integration_check.py
```

Strict precision-result reconstruction and the choice among finite
extreme-score corrections remain open for rank-full designs. The first
output-only boundary policy is now
implemented and retained under `jmle_extreme_score_integration_20260809/`.
It matched all 8,320 Person/run rows in the frozen evidence. The 120
structurally identified runs contained four all-minimum Persons; all four
rank-full movements of at least 1 logit were confined to those theta
coordinates, while maximum non-theta movement was 0.000293 logits. The app now
retains their existing `Estimate` only as a technical optimizer/constraint
value and withholds `ReportableEstimate`. Exact integer score patterns—not an
estimate magnitude or rounded display—drive the boundary decision. The same
guard identified 53 extreme Person rows in sparse runs, whose entire Person
output remains additionally unready because the eta design is rank deficient.

Reproduce the current-source bridge with:

```bash
MPLCONFIGDIR=/tmp/mfrm_matplotlib_cache \
  python3 validation/jmle_extreme_score_integration_check.py
```

No finite extreme-score correction has entered the application, and precision
polish has not entered the application optimizer. The first repository-only
comparison is now retained under `fixed_calibration_wle_20260809/` and
`fixed_calibration_wle_downstream_20260809/`. The generic fixed-calibration
Warm WLE scorer matched TAM theta and conditional-information SEs in all 16
frozen RSM/PCM/GPCM fixtures, including missing and exact-extreme response
patterns. The original TAM function hash recorded in the plan could not be
reproduced; it was not overwritten. A pre-result amendment froze an explicit
deparse hash and a secondary formals/body digest, both of which matched the
loaded TAM 4.3-25 code.

The downstream replay held polished facets and steps fixed across the 120
rank-full runs. It compared 7,600 Persons, found four exact-extreme shifts of
23.72--25.19 logits and a median absolute interior shift of about 0.009 logits,
and retained 51 Person-fit-zone changes. Displaying fit at `.3g` would change
31 JMLE and 23 WLE raw-zone decisions. Focal Holm and combined strong-bias
decisions were unchanged, but one practical `|bias| >= 0.50` decision moved
from 0.500839 under JMLE to 0.498068 under WLE. All 40 sparse runs remain
withheld because WLE Person scoring cannot identify fixed facets when the eta
design is rank deficient.

Reproduce both repository-only stages with:

```bash
python3 validation/fixed_calibration_wle_fixture.py \
  --output validation/fixed_calibration_wle_20260809
Rscript validation/fixed_calibration_wle_tam.R \
  --input validation/fixed_calibration_wle_20260809 \
  --output validation/fixed_calibration_wle_20260809 --repo .
python3 validation/fixed_calibration_wle_compare.py \
  --input validation/fixed_calibration_wle_20260809 \
  --output validation/fixed_calibration_wle_20260809
MPLCONFIGDIR=/tmp/mfrm_matplotlib_cache \
  python3 validation/fixed_calibration_wle_downstream.py
```

Neither result authorizes a Streamlit estimator option or automatic fallback.
The follow-on `cmle_wle_bridge_20260810/` evidence connects an inference-ready
exact CMLE RSM/PCM result to the same generic WLE core. Across 14 deterministic
Person cases, including two exact extremes, bridge/direct theta and SE values
and CMLE category-surface intercepts agreed exactly; coefficient order was
irrelevant and a missing coefficient failed closed. Reproduce it with
`python3 validation/cmle_wle_bridge_check.py`. This is fit-sample scoring with
calibration treated as fixed, not new-Person prediction or calibration-
uncertainty propagation.

The prospectively frozen follow-on under
`cmle_wle_calibration_sensitivity_20260810/` samples the exact-CMLE asymptotic
free-coordinate covariance and re-scores the same fit-sample Persons. With
2,000 draws for each of RSM and PCM, all 32,000 draw-Person scores returned.
Model-median calibration-draw SD was `0.0616`--`0.0801` logits; the maximum was
`0.518` logits, or `0.427` of the corresponding conditional WLE SE. Exact
extremes were the most sensitive cases. Zero covariance scale returned
numerical zero and identical seeds reproduced exactly; increasing covariance
scale increased the median sensitivity for both models. The draw quantiles are
not confidence intervals, and `QuadratureSensitivitySE` is not inference-
qualified because the CMLE and WLE stages reuse responses and their dependence
is absent from the independent normal coefficient draws. The retained
`cmle_wle_first_read_projection.csv` separates research-ready computation,
cautions, and withheld inference/UI claims. Reproduce it with:

```bash
MPLCONFIGDIR=/tmp/mfrm_cmle_wle_matplotlib \
  python3 validation/cmle_wle_calibration_sensitivity.py
```

The original prospective plan is unchanged. Because it hashed the design
record itself, documenting the observed result changed that non-executable
identity; the explicitly post-result, documentation-only
`cmle_wle_calibration_sensitivity_documentation_amendment_20260810.json`
records both hashes and states that no numerical contract or gate changed.

The next stage was prospectively defined before implementation in
`cmle_wle_bootstrap_plan_20260810.json` and
`../docs/cmle_wle_bootstrap_design.md`. It keeps two estimands separate:
fixed-score conditional-pattern resampling diagnoses calibration/WLE coupling
at the observed Person total, while joint plug-in response resampling allows
Person totals to vary under fitted CMLE + WLE parameters. The planned 200
replicates per RSM/PCM lane are an implementation pilot only; no success-rate
promotion threshold or interval claim has been registered. The initial core is
now in `mfrm_app/cmle_wle_bootstrap.py`; its six sampler/refit/resource tests
pass. The frozen pilot under `cmle_wle_bootstrap_pilot_20260810/` completed all
800 attempted refits (200 per RSM/PCM lane) with full conditional rank, CMLE
inference readiness, and WLE availability. Fixed-score totals had zero
mismatches; same-seed replay was exact; and joint category frequencies were at
most `1.57` Monte Carlo SD from expectation. Median Person bootstrap SD was
`0.0602` logits for the fixed-score lane and `0.615` for the joint lane;
maximum SD was `0.744`. The joint lane produced 321 extreme-to-interior and 120
interior-to-extreme transitions. The result is
`pilot_complete_promotion_withheld`: repeated-truth coverage remains
unavailable, percentile output is not a confidence interval, and lane variance
differences are not a decomposition. Reproduce the immutable original pilot
with:

```bash
MPLCONFIGDIR=/tmp/mfrm_cmle_bootstrap_matplotlib \
  python3 validation/cmle_wle_bootstrap_pilot.py
```

The post-pilot Person-fit gate is registered in
`cmle_wle_person_fit_plan_20260810.json`. The independent Python kernel in
`../mfrm_app/cmle_wle_fit.py` evaluates untrimmed fixed-theta Person Infit and
Outfit, fails closed on non-positive/non-finite fitted variance, classifies from
raw values, and leaves ZSTD/p-values explicitly unavailable. The frozen
`cmle_wle_person_fit_20260810/` comparison passed for 32 RSM/PCM Person rows and
162 administered observations, including four exact-extreme Persons and 15
missing rectangular cells per model. Python and `sirt::pcm.fit` differed by at
most `4.89e-15` for Infit and `4.00e-15` for Outfit; row-level expectations,
variances, fourth moments, and squared-standardized residuals also matched
within `4.89e-15`. `TAM::tam.jml.fit` shares the Infit formula and pre-trim
Outfit contributions but trims unusually large Outfit contributions by
default, so it is retained as a formula audit rather than an unconditional
Outfit parity target. Reproduce the cross-engine gate with:

```bash
python3 validation/cmle_wle_person_fit_validate.py
```

The prospectively registered bootstrap follow-on is in
`cmle_wle_bootstrap_fit_extension_plan_20260810.json`. The first run under
`cmle_wle_bootstrap_fit_extension_20260810/` remains a failed result: all 800
Person-fit refits were available, but an asymmetric comparison of the parsed
original CSV with unpersisted extension floats produced a one-ULP
`4.44e-16` WLE identity difference against an exact-zero gate. The gate and
tolerance were not relaxed. Before correction,
`cmle_wle_bootstrap_fit_extension_amendment_20260810.json` froze the failure
identity and required a complete rerun with symmetric `%.17g` persisted
readback. The corrected evidence under
`cmle_wle_bootstrap_fit_extension_corrected_20260810/` passed: 800 attempts and
12,800 Person-replicates were fit-ready, the minimum per-cell Wilson 95% lower
bound was `0.981`, fixed-score total mismatches and WLE identity differences
were both zero, and neither the original nor failed inventory changed. Raw
classes changed 2,372 times for Infit and 2,454 for Outfit; 2,554
Person-replicates changed at least one. Twenty-four values were within the
three-decimal display boundary and nine displayed values implied a different
class, while numerical-ULP boundary counts were zero. The boundary-focused
figure enlarges those rare but decision-relevant rows rather than hiding them
under 12,000+ stable values. Reproduce the corrected extension and descriptive
visual with:

```bash
MPLCONFIGDIR=/tmp/mfrm_matplotlib_cache \
  python3 validation/cmle_wle_bootstrap_fit_extension.py
MPLCONFIGDIR=/tmp/mfrm_matplotlib_cache \
  python3 validation/cmle_wle_bootstrap_fit_visualize.py
```

`CMLE_WLE_BOOTSTRAP_FIT_CRITICAL_REVIEW.md` further separates transition
direction, baseline classes, Person concentration, and baseline-extreme
behavior. Per-statistic transition shares were 14.09%--14.66% in the
fixed-score lane and 22.47%--24.47% in the joint lane; baseline exact extremes
had 0/800 transitions when their totals were fixed and 206/800 when joint
resampling could change their totals/status. The largest Person-specific share
was 40.5%, reinforcing that the aggregate is not a misfit-prevalence estimate.

A separate post-result audit freezes a 125-triplet grid around
0.50/1.50/2.00 and classifies only retained raw MnSq. It reproduces all
canonical decisions, while counterfactual rounding produces nine class
mismatches and five either-statistic transition disagreements at three
decimals and none at six decimals. Full-grid aggregate transition ranges are
1,789--3,461 (Infit), 1,877--3,550 (Outfit), and 1,975--3,675 (either). Run it
with:

```bash
MPLCONFIGDIR=/tmp/mfrm_matplotlib_cache \
  python3 validation/cmle_wle_fit_threshold_surface.py
```

After the noisy/distorting one-at-a-time binary profile was observed to be
flat, a timestamped diagnostic amendment froze a complete nominal-class
decomposition before those cell counts were computed. Its 16-cell matrices
reproduce all parent off-diagonal counts and the stored canonical matrix.
Moving only the noisy/distorting boundary from 1.90 to 2.10 leaves aggregate
Infit/Outfit/either transitions exactly at 2,372/2,454/2,554, but reallocates
37 Infit and 51 Outfit replicate labels from `distorting` to `noisy`. Reproduce
the addendum with:

```bash
MPLCONFIGDIR=/tmp/mfrm_matplotlib_cache \
  python3 validation/cmle_wle_fit_threshold_severity.py
```

These are descriptive lane-specific transitions and post-result sensitivity
surfaces, not threshold validation.

The next Phase-A gate was registered before its known-truth outcomes. Its eight
paired fixed-calibration conditions vary observations per Person, categories,
uniform-random contamination, local response copying, and Person-specific
threshold heterogeneity. All 160 attempts returned, retaining 624,000 response
rows and 32,000 Person rows for each Person-measure source. The raw canonical
either-upper clean rate was 3.4%--3.7% with 24 observations and 10.1%--10.8%
with six; mechanism-specific affected detection ranged from 3.8% to 46.6%.
Run it with:

```bash
MPLCONFIGDIR=/tmp/mfrm_matplotlib_cache \
  python3 validation/cmle_wle_fit_known_truth_pilot.py
```

A transparently registered post-result addendum reports WLE recovery and
effective threshold dimensions without changing the parent evidence. It
retains 163 exact-extreme rows, reports non-extreme RMSE of 0.263--0.720
logits, and finds 2,333/32,000 canonical WLE versus generating-theta flag
disagreements. The nominal 125 triplets reduce to 5 (`either_upper`), 5
(`either_distorting`), 25 (`either_nonacceptable`), and 5 (`either_overfit`)
effective inputs. Reproduce it with:

```bash
MPLCONFIGDIR=/tmp/mfrm_matplotlib_cache \
  python3 validation/cmle_wle_fit_known_truth_pilot_addendum.py
```

The critical review documents why sparse clean rates, mechanism-specific
detection, and fit/error reversals block a universal traffic light. These are
20-replicate fixed-calibration design checks, not confirmatory operating
characteristics. ZSTD, p-values, coverage-qualified intervals, total SE,
estimated-CMLE connectivity, native anchor contamination, matched repeated
cross-engine fits, and Streamlit integration remain withheld.

The native exact-CMLE Phase-B increment is registered in
`cmle_native_hard_anchor_plan_20260810.json`; the smoke-specific retained-run
and anchor-level choices are frozen in
`cmle_native_hard_anchor_smoke_amendment_20260810.json`. It implements affine
hard facet anchors across the conditional likelihood, category surfaces,
fit-sample WLE, Person fit, bootstrap kernels/refits, and calibration draws.
The 24,000-row, 8-RunId, 32-fit smoke passed every engineering contract and is
reproduced with:

```bash
MPLCONFIGDIR=/tmp/mfrm_matplotlib_cache \
  python3 validation/cmle_native_hard_anchor_smoke.py
```

After inspecting that smoke, a separate post-result diagnostic was registered
before computation. It confirms the algebraic common-origin invariance: a
uniform +0.25 shift in all selected Rater anchors changes the Rater/WLE origin
but cancels from the fixed-total conditional likelihood and Person fit. The
maximum conditional-log-likelihood difference was `1.82e-12`, maximum
Infit/Outfit difference was `7.37e-08`, and raw either-upper mismatches were
zero. Reproduce it with:

```bash
python3 validation/cmle_native_hard_anchor_common_shift_diagnostic.py
```

The next differential-anchor stress was separately frozen before computation
in `cmle_native_anchor_differential_stress_plan_20260810.json`. It reused
120,000 retained response rows over 40 RunIds and crossed four dense/sparse,
clean/random conditions with eight anchor scenarios. Reproduce all 320 fits
with:

```bash
MPLCONFIGDIR=/tmp/mfrm_matplotlib_cache \
  python3 validation/cmle_native_anchor_differential_stress.py
```

All 320 fits returned inference-ready and passed input, anchor, dimension,
rank, and raw-fit recomputation gates. That numerical success did not validate
the anchors. Correct nested anchor sets did not improve WLE RMSE monotonically;
share was confounded with which deterministic Raters were fixed and with which
free Raters remained in the recovery denominator. A common +0.25 error changed
no raw fit flags, whereas a zero-mean differential error changed 2.4%--15.8%
of paired raw either-upper decisions by condition/group and moved some Person
WLEs by as much as `0.287525` logits. Wrong fixed values occasionally achieved
higher in-sample conditional likelihood in sparse samples. Across 64,000
Person-scenario rows, 3-decimal display-driven decisions disagreed with raw
decisions 12 times; 6-decimal values disagreed zero times in this dataset.
The report and post-result critical review are under
`cmle_native_anchor_differential_stress_20260810/`.

The next prospective plan, `cmle_native_anchor_connectivity_plan_20260810.json`,
randomizes correct anchor content by replicate and directly manipulates
within-Person Rater bridge topology. Reproduce 70 unique datasets and 350
topology-anchor audits with:

```bash
MPLCONFIGDIR=/tmp/mfrm_matplotlib_cache \
  python3 validation/cmle_native_anchor_connectivity_stress.py
```

The realized non-extreme-Person graph prediction agreed with exact prefit
eligibility in 350/350 cases. Partial anchors did not identify a fully
disconnected design. In the strong two-component design, two randomized
correct anchors were eligible in 7/10 replicates—exactly those in which both
components were covered. All 250 connected cases passed prefit, but the
one-Person-per-edge minimal chain yielded one non-finite trial-point failure
and 12 returned fits with deficient final information rank. The other 200
weak/strong chain/hub fits were inference-ready. Across 82,200 inference-ready
Person rows, 3-decimal MnSq changed 13 raw either-upper decisions and 6-decimal
MnSq changed none in this dataset. The generated report, figures, ledgers, and
post-result critical review are under
`cmle_native_anchor_connectivity_20260810/`.

The finite-domain optimizer remediation was separately frozen before changing
the core in `cmle_finite_domain_optimizer_remediation_plan_20260810.json`.
Reproduce its same-byte 287-fit replay and selected tests with:

```bash
MPLCONFIGDIR=/tmp/mfrm_matplotlib_cache \
  python3 validation/cmle_finite_domain_optimizer_remediation.py
```

All 287 prefit-eligible fits now returned and exactly the same 274 were
inference-ready. The former exception recorded 12 rejected non-finite BFGS
trials and returned with rank 8/10, nullity 2, gradient sup norm `0.00296`, and
inference withheld. The 286 previously returned fits had zero readiness
mismatches; all 82,200 shared Person raw/rounded flags agreed. Maximum absolute
differences were `3.55e-15` for structural estimates and `2.27e-13` for
conditional log likelihood. The public exact objective still raises on a
direct non-finite evaluation, and the deliberately tested best-finite fallback
forces optimizer success false. See the report and critical review under
`cmle_finite_domain_optimizer_remediation_20260810/`.

Finite-CMLE existence was then frozen independently in
`cmle_finite_mle_existence_plan_20260810.json`. The first research audit
enumerates every attainable fixed-score sufficient statistic, constructs the
conditional convex-support cone, and searches for a nonzero supporting
direction across a three-value tolerance grid. Reproduce it with:

```bash
MPLCONFIGDIR=/tmp/mfrm_matplotlib_cache \
  python3 validation/cmle_finite_mle_existence_stress.py --overwrite
```

All 93 prospective structural/boundary/interior controls matched, including
binary complete/reverse separation, near-boundary minority responses,
multi-Rater graph cases, RSM/PCM, unused declared categories, and correct or
contaminated anchors. The deliberately low 50,000-configuration retained cap
classified 63 structural, 13 boundary, and 47 interior cases and retained 227
as typed unavailable; six 192,745--217,093-configuration high-cap sentinels
were all interior. That result validates the criterion but not scalability.

The scalable replacement was separately frozen in
`cmle_finite_mle_oracle_plan_20260810.json`. It uses exact fixed-score dynamic
programming as a support-function oracle and adds only violated support rows to
each bounded coordinate LP. Reproduce it with:

```bash
MPLCONFIGDIR=/tmp/mfrm_matplotlib_cache \
  python3 validation/cmle_finite_mle_oracle_stress.py --overwrite
```

The final evidence matched 1,184/1,184 exhaustive support maxima, 93/93 fixture
statuses, 123/123 completed retained exhaustive statuses, and 6/6 high-cap
sentinels. The same-byte 350-case replay resolved 63 structural, 13 boundary,
and 274 interior cases with no unavailable or tolerance-unstable result. All
274 current inference-ready fits were interior; the 13 boundary cases were
already non-ready. Median theoretical configurations were 97,040 versus 163
generated constraints. Eligible-case P95/max time was 0.521/0.568 seconds and
the 350-case total was 109.216 seconds, passing the frozen performance gate.
The initial 156.290-second implementation failed that gate before a value-only
DP avoided unnecessary maximizing-response reconstruction; no cases or
tolerances were removed. Reports, ledgers, cut histories, and the critical
review are under `cmle_finite_mle_oracle_20260810/`.

The readiness connection was frozen separately in
`cmle_finite_mle_readiness_integration_plan_20260810.json`. Reproduce the full
287-fit plus 82,200-Person same-byte replay with:

```bash
MPLCONFIGDIR=/tmp/mfrm_matplotlib_cache \
  python3 validation/cmle_finite_mle_readiness_integration.py --overwrite
```

All 287 fits returned. The 274 oracle-interior fits were inference-ready; all
13 boundary fits were non-ready and carried
`finite_mle_boundary_no_finite_cmle`. Independent oracle status and the frozen
integrated-readiness rule matched 287/287. Structural estimates, SE,
conditional likelihood, WLE, Infit, and Outfit differed from the prior
finite-domain evidence by no more than `2.27e-13`; all raw and rounded-
counterfactual flags matched over 82,200 Person rows. The selected suite passed
102/102 tests. See the report and critical review under
`cmle_finite_mle_readiness_integration_20260810/`.

The user-facing early-stop contract was frozen next in
`cmle_structured_workflow_plan_20260810.json`. The Streamlit-free v1 workflow
always returns input/design, finite-existence, structural-optimization, and
Person-scoring-availability stages with Japanese and English text. Reproduce
the fixture and same-byte retained replay with:

```bash
MPLCONFIGDIR=/tmp/mfrm_matplotlib_cache \
  python3 validation/cmle_structured_workflow_stress.py --overwrite
```

All 93 fixture and 350 retained terminal states matched. The workflow stopped
63 structural and 13 boundary cases before optimization and optimized only 274
interior cases; forbidden optimizer attempts were zero. Structural estimates
and SE matched the direct integrated path within `4.44e-16`, likelihood within
`2.27e-13`, and anchor-exact mismatches were zero. P95 times were
0.068/0.306/1.084 seconds for structural block, boundary stop, and ready
interior workflow. Person scoring is available/not-run in v1, and nonempty
bilingual messages are not comprehension evidence. See the report and review
under `cmle_structured_workflow_20260810/`.

Cross-engine boundary semantics were then frozen in
`cmle_cross_engine_boundary_plan_20260810.json` with two result-blind source-
identity amendments necessitated by concurrent changes in the external mfrmr
development tree. The second amendment pins a build-minimal mfrmr 0.2.3
snapshot whose pre-copy, copied, and post-copy relative content hashes all
equal `06870f2f...21d3`. Reproduce the registered 13-case run while that
temporary snapshot remains available with:

```bash
MPLCONFIGDIR=/tmp/mfrm_matplotlib_cache \
  python3 validation/cmle_cross_engine_boundary.py
```

The contract passed. For the five oracle-interior additive-RSM controls,
native Python and matched `immer::immer_cml()` free-coordinate estimates
differed by at most `1.78e-9` and conditional log likelihood by at most
`5.33e-15`. For seven oracle-boundary controls, `immer` at `maxit=2000`
returned code zero in 3/7 and finite coefficient/SE vectors in 4/7; increasing
the cap to 10,000 raised code-zero returns to 4/7, but the oracle-qualified
final readiness remained 0/7. The boundary ledger retains raw gradients and
classifications at `1e-4`, `1e-5`, `1e-6`, and `1e-8`; no printed rounding or
MnSq participates in finite-existence classification.

The same rows were also passed to deliberately non-parity sensitivity lanes.
The frozen mfrmr JML snapshot returned and reported convergence for all 7/7
Python-CMLE boundary cases, with category-support warnings and maximum
absolute estimates up to `103.75`. TAM and sirt MML each returned one of seven
boundary fits; unsupported declared category support and other failures remain
typed. sirt is isolated in one child R process per case after an adapter-
development input-shape error demonstrated that a native crash could otherwise
erase all engine evidence. These results do not compare JML/MML numbers with
CMLE. The decision, 61 recursively hashed outputs, detailed maxit table,
engine-by-status table, and critical review are under
`cmle_cross_engine_boundary_20260810/`.

These are engineering and ten-replicate debugging diagnostics, not performance
results. Universal anchor share, drift distributions, group anchors, matched
PCM/hard-anchor cross-engine behavior, downstream workflow Person scoring, and
UI comprehension remain withheld.

Sparse designs require prospective redesign, anchoring that actually connects
the estimand, or a separately justified and qualified estimator; brute-force
replication cannot resolve structural nonidentification.

For cross-engine handoff, the Python runner retains the exact ratings, facet
truth, and anchor rows with per-run identities and SHA-256 file hashes. The
companion `operating_characteristics_bridge.R` validates those same bytes and
imports them without recreating R-side random draws. The retained bridge smoke
passed all 16 R checks for 16 RunIds and 18,526 rating rows using an isolated
dirty-development snapshot of mfrmr 0.2.3 at source commit `7ee1fd5`, TAM
4.3-25, immer 1.5-13, and sirt 4.2-133. This is input readiness; the four R fit
adapters do not become equivalent merely by passing it.

The first actual adapter now fits mfrmr JMLE in separate Python-control and
strict-gradient modes (`operating_characteristics_mfrmr.R`), then normalizes
the evidence with `operating_characteristics_compare.py`. Strict mfrmr results
were inference-ready for all 12 identified non-sparse RunIds and rejected all
four rank-deficient sparse RunIds before optimization. Across the 108 facet
estimates jointly included by Python and strict mfrmr, the mean absolute
difference was 0.0000124 logits and the maximum grouped absolute difference was
0.000146 logits. All 12 jointly eligible focal strong-bias decisions agreed;
focal bias estimates differed by 0.0000287 logits on average.

In contrast, the looser Python-control mfrmr mode returned the 12 identified
fits with optimizer code 0 but withheld inference readiness for all of them
because terminal gradients ranged from 0.0283 to 0.1784. This is retained as a
stopping-rule sensitivity finding, not pooled with the strict mode. The
two-replicate result remains a numerical/readiness smoke comparison, not an
operating-characteristic or package-superiority claim. TAM is retained in the
separate MML lane described below.

The second actual adapter lane compares the native exact Python CMLE core with
`immer::immer_cml()` on the same byte-validated ratings. It retains 12
structurally eligible RunIds (8 unique eligible rating datasets). Independent
Python exact-moment and R/immer conditional-information preflights agree on
rank/nullity for 16/16 RunIds, including rejection of four sparse RunIds at
rank 5/12. Across all 96 matched returned free coordinates, the maximum estimate
difference is 0.0000000563 logits and the maximum conditional-loglikelihood
difference is 0.00000000000477. Exact-information and `immer` Hessian SEs also
agree to a maximum absolute difference of 0.00000000641.

This close likelihood/coordinate agreement does not erase the numerical
readiness distinction. All 12 `immer` fits returned optimizer code 0, but its
terminal gradients pass 12/12 RunIds at `1e-4`, 7/12 at the primary `1e-5`,
and 0/12 at `1e-6`; Python's exact-Newton polish passes its declared threshold
for all 12. The comparison therefore preserves all returned numerical pairs
while limiting primary inference-ready agreement to 7 RunIds. See
`operating_characteristics_20260809/CMLE_IMMER_RESULTS.md` and the paired
coefficient/readiness figures.

Neither CMLE adapter consumes anchors or estimates the focal Rater x Task
interaction. The eight anchor-requesting RunIds are explicitly unanchored
parity references and ineligible for an anchored-condition claim. Clean and
contaminated anchor conditions share rating bytes, so their identical CMLE
results demonstrate ignored anchor inputs—not resistance to anchor bias.
The paired null/+0.60 generator comparison is retained separately as omitted-
interaction leakage: in `balanced_small`, additive CMLE free coordinates moved
by 0.0612 logits on average and at most 0.172 logits. This is not a CMLE local-
bias estimate or a detection-rate result.

The third adapter lane fits `sirt::rm.facets()` with 30- and 61-point fixed
quadrature grids. Task x Criterion combinations become six virtual PCM items;
Rater severity is retained, and Person ability follows an estimated common
normal population distribution. That is MML, not the Python additive-RSM JMLE
or exact-CMLE estimand, so `operating_characteristics_sirt_compare.py` labels
the Python comparison as cross-estimator sensitivity rather than parity.

All 32 q30/q61 fits returned and 30 were analysis eligible. The same sparse
`+0.60` RunId reached the 1,200-iteration cap under both grids and is retained
but excluded from eligible summaries. In balanced-small data, jointly included
Python/sirt Rater estimates differed by 0.01232 logits on average and 0.03481
at most, preserving all four within-run Rater orderings. In sparse data, sirt
was eligible for 3/4 RunIds while exact CMLE and strict mfrmr rejected all four
structurally; this is an assumption-based bridge through the common Person
distribution, not observed-design connectivity.

Quadrature results distinguish raw from identified comparisons: among eligible
runs the maximum q61-q30 difference was 0.001662 logits for Rater severity,
0.001387 for mean-centered virtual-item location, 0.05011 for raw virtual-item
location, and 0.07862 for Person EAP. A +0.25 contaminated-anchor input moved
sirt Rater estimates by essentially +0.25, demonstrating transmission rather
than robustness. See `operating_characteristics_20260809/SIRT_RESULTS.md`; the
paired figures mark returned-but-excluded values separately. Information
criteria are withheld because the inspected sirt 4.2-133 parameter-count path
uses the numeric Rater-centering mode in the count. This remains a smoke
sensitivity result, not an operating-characteristic or superiority claim. The
ordered `sirt_first_read_summary.csv` is the future one-click UI projection,
but every row retains `PublicSurfaceEnabled=False` and its final status is
`Withheld`.

The fourth adapter lane fits `TAM::tam.mml.mfr()` under the exact additive RSM
surface: Criterion item + Rater + Task + common step. This is closer to the
Python response-surface specification than the sirt virtual-item mapping, but
TAM still integrates over an estimated normal Person population and is not
JMLE/CMLE parity. The primary grid uses 61 fixed nodes on -8--8; a 21-node mode
is retained separately as quadrature sensitivity.

All 32 TAM fits returned and met the retained loop/progress convergence audit;
30 were analysis eligible. Two q21 sparse fits stored Person variance as
`0.0010000001`, effectively TAM's configured `0.001` lower boundary. A naive
exact `<=0.001` comparison missed both, while the registered `1e-8` numerical
band caught them. All 16 primary q61 fits were eligible. This is a concrete
example of floating-point policy changing readiness even when displayed values
look identical.

In balanced-small data, Python-JMLE/TAM-MML Rater differences averaged 0.01383
logits and reached 0.03989. TAM and sirt Rater estimates were closer there
(MAE 0.001513), but that agreement did not survive contaminated anchors. TAM's
Rater sum-zero constraint gave the two fixed anchors their injected +0.25 shift
and forced the two unanchored Raters to compensate by about -0.25; sirt instead
transmitted the location shift differently. Exact anchor reproduction is not
anchor robustness.

Quadrature sensitivity was output-dependent. Among jointly eligible runs, the
maximum q61-q21 difference was 0.02224 logits for Rater, 0.3814 for Person EAP,
and 0.6990 for cumulative response-surface difficulty. Primary TAM was eligible
for all four sparse RunIds while exact CMLE and strict mfrmr rejected all four
by observed-design rank. This is assumption-based MML identification, not new
observed connectivity. See `operating_characteristics_20260809/TAM_RESULTS.md`,
the variance-boundary audit, and the three paired figures.

`tam_first_read_summary.csv` is the future one-click projection and remains
fully disabled for public use. Derived last-level TAM coverage is withheld
because the inspected 4.3-25 expansion uses a diagonal matrix of free-xi SEs
without a retained full covariance; information criteria are retained only for
within-TAM audit.

## Native exact CMLE structural calibration

The repository also contains a separate 2026-08-09 deterministic CMLE
implementation check:

- numerical/design contract: `../docs/cmle_phase0_design.md`;
- pure Python core: `../mfrm_app/cmle.py`;
- reused 20-dataset stress runner: `cmle_stress.py`;
- exact-enumeration, negative-control, and immer parity tests:
  `../tests/test_cmle.py`;
- R-side matched-design runner: `../tests/data/run_immer_cmle_parity.R`; and
- retained result summary: `cmle_phase0_20260809/RESULTS.md`.

The separate scaling envelope is reproduced by `cmle_scaling.py` and retained
under `cmle_scaling_20260809/`. It varies Person count, complete response-unit
count, PCM parameter count, and unique missingness patterns in fresh processes,
and includes negative controls for the default work and 512 MiB memory guards.

This research core is not called by Streamlit and does not add a public
estimator option. It supports RSM/PCM structural calibration and native hard
facet anchors only; Step/group/soft anchors remain unsupported. It blocks
GPCM, non-unit row weights, duplicate response units, malformed anchors, and
deficient conditional-information rank. Its rank and covariance evidence use exact
conditional sufficient-statistic moments, not a finite-difference Hessian.
Raw BFGS status and any safeguarded exact-Newton polishing remain separately
auditable.

The 2026-08-10 PCM hard-anchor cross-engine increment is governed by
`cmle_pcm_anchor_cross_engine_plan_20260810.json`, its two transparent identity
remediation records, and retained evidence in
`cmle_pcm_anchor_cross_engine_20260810/`. The final contract passed. Native
Python exact PCM CMLE and a matched `immer::immer_cml()` `W`/`b_const`
parameterization agreed across six oracle-interior cases within `2.45e-9` for
free coordinates and `3.55e-14` for conditional log likelihood; zero-point
objective and gradient differences were at most `4.27e-14`. All six
separation/unused-declared-support conditions remained blocked by the exact
support oracle at each retained optimizer cap. The mfrmr JML and TAM/sirt MML
tables are different-estimand sensitivity ledgers, not parity evidence, and
TAM/sirt hard-anchor cases are explicitly unsupported. The frozen mfrmr 0.2.3
source is dirty and temporary, so a clean-release replay remains mandatory.
The critical review and raw full-precision CSV files are authoritative; public
UI remains withheld.

The subsequent one-click view-model contract is frozen in
`cmle_one_click_result_contract_plan_20260810.json` and retained under
`cmle_one_click_result_contract_20260810/`. Its five analytical states all
returned six ordered bilingual cards. Exactly two ready calibrations ran
fixed-calibration WLE/MnSq; input-invalid, structurally unidentified, and
finite-support-boundary cases did not call Person scoring. Across 14 ready
Persons, WLE estimate, conditional SE, Infit, and Outfit exactly reproduced a
direct call on the same fitted calibration. A constructed `1.5004` Infit
remained raw-classified as `noisy` while its displayed `1.500` counterfactual
was `acceptable`. This is deterministic orchestration and presentation-data
evidence, not a task-based usability/comprehension study. All public-surface
flags remain false.

Current-schema private save/reload/replay identity is frozen next in
`cmle_one_click_archive_identity_plan_20260810.json` and retained under
`cmle_one_click_archive_identity_20260810/`. Three deterministic ZIPs passed
21/21 replay checks. Repeated builds were byte-identical, all 56 archived
WLE/conditional-SE/Infit/Outfit values reloaded with the same binary64 hex
representation, and all six bilingual cards survived round-trip. Ordered and
semantic input hashes are separate: a row permutation retains the latter but
changes AnalysisID. Reversing hard-anchor rows is normalized and leaves the
configuration identity unchanged. Changed, missing, undeclared, and duplicate
entries are rejected before loading. The ZIPs contain raw rows and identifiers;
public-mode construction is deliberately rejected. They are not encrypted,
de-identified, anonymized, or approved for sharing, and cross-version/platform
replay remains unqualified.

Non-public bilingual instrument readiness is frozen in
`cmle_one_click_comprehension_readiness_plan_20260810.json` and retained under
`cmle_one_click_comprehension_readiness_20260810/`. Five machine states produce
10 English/Japanese static previews and 100 participant tasks with a separate
100-row scoring key. Six critical misconception domains cannot be offset by a
higher total score. Thirty synthetic packets exercised correct, all-wrong,
false-green, duplicate, missing, unknown-option, and invalid-time behavior;
the automated contract and 40 selected tests passed. These are scoring and
static-structure checks only. `human_gate_status.csv` records zero participants,
no comprehension rate, no language-equivalence result, and a disabled public
surface. Real-browser rendering and assistive-technology testing were not
available in this run and are not claimed.

The next private operations-only increment is registered in
`cmle_one_click_cognitive_interview_operations_plan_20260810.json` and retained
under `cmle_one_click_cognitive_interview_operations_20260810/`. Within each
language, 10 planned session slots exhaust all ten unordered pairs of the five
analytical states while balancing case order and novice/experienced target
strata. The kit contains 20 answer-free participant packet files (400 tasks),
a separate 100-row private moderator guide, a 400-row blank record template,
English/Japanese scripts, and a data dictionary. Fourteen synthetic validator
cases and 47 selected tests passed. Exact-column validation rejects direct or
contact identifier fields and typed record failures, but manual `PIIReviewed`
attestation is not automated de-identification. The evidence grants no ethics
approval or recruitment authority, contains no human data, and leaves
comprehension, language equivalence, accessibility, and public UI unassessed.

Frozen version identity and directional confirmatory mechanics are registered
in `cmle_one_click_versioned_confirmatory_gate_plan_20260810.json` and retained
under `cmle_one_click_versioned_confirmatory_gate_20260810/`. Instrument
identity covers exact task/key bytes, all 10 preview identities, and the six
directional rules. Any task, key, preview, or rule mutation changes the
composite identity; content aliases and cross-content pooling fail closed.
Dangerous-error numerators use one eligible primary case role per slot and are
not synonymous with every wrong answer. English and Japanese remain 12
separate cells. The one-sided 95% Wilson zero-error upper bound is `0.101310`
at `n=24` and `0.097654` at `n=25`; the 12-cell Bonferroni sensitivity remains
`0.217782` at `n=25`. Ten synthetic gate scenarios and 55 selected tests
passed. These are arithmetic and validator results, not a sample-size
recommendation, power calculation, pilot result, confirmatory result, or
language-equivalence claim. The minimum confirmatory n and public UI remain
unregistered/withheld.

The prospective no-human-data sensitivity plan is registered in
`cmle_one_click_confirmatory_sample_size_sensitivity_plan_20260810.json`, with
a timestamped amendment documenting that preliminary arithmetic revealed a
nonmonotone sawtooth. Evidence is retained under
`cmle_one_click_confirmatory_sample_size_sensitivity_20260810/`. At assumed
true dangerous-error probability `0.05`, the exact single-cell 90% pass-
probability first crossing is `n=224`, but the first n sustained through the
registered `n=5000` search is `n=260`. The conservative 12-cell union-bound
counterparts are `n=422` and `n=456`. Attrition assurance, cluster design-
effect heuristics, and blocked-mechanism exposure are retained as separate
sensitivity tables rather than combined into a false precision claim. Two
figures were visually inspected, and 63 selected tests passed. The search
maximum is finite, retention assumes independent homogeneous Bernoulli loss,
and the cluster calculation is not an exact Wilson correction. No minimum
valid n or recruitment total is registered; human participants remain zero
and the public surface stays withheld.

The subsequent registered synthetic dependence/MNAR stress is retained under
`cmle_one_click_confirmatory_dependence_stress_20260810/`. A Gaussian-threshold
Bernoulli generator preserves registered marginal dangerous-error
probabilities while varying participant latent dependence, consecutive fixed-
size moderator/site-like clusters, domain or blocked-mechanism heterogeneity,
and outcome-dependent retention. Eight scenarios by four planned-slot values
and 2,000 replicates yielded 64,000 joint-gate replicates, 384 language-domain
rows, and 384 mechanism rows. At n=500, the danger-under-retained MNAR scenario
had complete-data truth `0.099894`, observed valid-record risk `0.060315`, and
false reassurance `0.4035` (Monte Carlo Wilson 95% `0.3822`–`0.4252`). The
blocked-mechanism hotspot and combined adversarial false-reassuring rates were
`0.4955` and `0.9675`. All registration, deterministic replay, marginal,
exact-binomial, retention, logic, mechanism-allocation, and Monte Carlo
precision audits plus 75 selected tests passed. Latent rho is not observed
binary ICC, the MNAR models are not identified from data, and simulated equal
language parameters are not equivalence evidence. No n or recruitment total
is selected; human participants remain zero and public UI remains withheld.

The fail-closed private pre-recruitment workflow is registered in
`cmle_one_click_confirmatory_protocol_preflight_plan_20260810.json` and retained
under `cmle_one_click_confirmatory_protocol_preflight_20260810/`. Its decision
register has five fixed rows and 13 unresolved blocking rows. The exact-schema
zero-row attempt ledger, 22-field dictionary, and 15-code categorical
invalidity vocabulary reject direct identifier/free-text columns, slot
coverage drift, unregistered reasons, outcome-aware exclusion, and incomplete
review attestations. The 72-row partial-identification surface reports both
all-invalid-safe and all-invalid-dangerous bounds. With valid n=500 and
observed 5% error, the registered invalid/valid ratio 0.025 remains worst-case
robust, while 0.05 does not (`WilsonUpper95=0.118434`). The deterministic
nine-member private ZIP reproduces byte-for-byte and fails hash plus first-read
checks after a synthetic `BLOCKED` to `READY` mutation. Eight ledger attack
fixtures and 91 selected tests passed. The human-study state nevertheless
remains not started: ethics approval, recruitment authority, n selection,
human data, confirmatory results, language equivalence, and public UI are all
false/absent. Schema and PII attestations do not prove de-identification or
regulatory compliance.

The non-recommending decision workbench is registered in
`cmle_one_click_confirmatory_decision_workbench_plan_20260810.json`, with an
implementation-preceding broad-stage/topological-order clarification, and is
retained under `cmle_one_click_confirmatory_decision_workbench_20260810/`.
Thirteen blockers have 39 exact options and 32 acyclic prerequisite edges.
Every option retains strength, primary risk, evidence, and qualitative impact
fields; recommendation, rank, default, and auto-selection flags remain false.
SCI-02 numeric sample size has seven prerequisites. Seven synthetic selection
attacks and 105 selected tests passed. The self-contained bilingual HTML has no
form, script, or external resource, and the deterministic seven-member ZIP
rejects a `BLOCKED` to `READY` mutation. OptionsSelected remains 0,
SelectionsComplete and SubstantiveEvidenceVerified remain false, and the
workbench cannot set RecruitmentReady or PublicSurfaceEnabled true.

SCI-01 estimand decision support is separately registered in
`cmle_one_click_confirmatory_sci01_estimand_memo_plan_20260810.json` and
retained under `cmle_one_click_confirmatory_sci01_estimand_memo_20260810/`.
The three unranked alternatives are defined by target population, formula,
invalid-outcome handling, identification requirements, allowed/forbidden
claims, and required prospective material. Applying all three views to the 72
frozen integer-count scenarios yields 216 projection rows. Thirty-six source
scenarios have the joint pattern `observed-valid pass`, `all-invalid-dangerous
proxy fail`, and `dual observed-pass-not-robust`; no label uses displayed
rounding. Because the source surface contains neither scheduled-slot totals nor
a prospective invalidity-code-to-composite map, the Option B column is
explicitly not the final scheduled-slot estimand. Six selection attacks and
118 selected regression tests passed. The HTML and seven-member private ZIP
are read-only, self-contained, deterministic, and blocked. This is static
offline decision support, not public browser/accessibility acceptance or a
user-facing analysis button. OptionsSelected, sample-size selection,
RecruitmentReady, human participants, confirmatory outcomes, and public UI all
remain zero/false.

The shared SCI-01/SCI-05 invalidity adjudication contract is registered in
`cmle_one_click_confirmatory_invalidity_adjudication_plan_20260810.json` and
retained under `cmle_one_click_confirmatory_invalidity_adjudication_20260810/`.
All 15 codes remain distinct; only `none` is fixed and the other 14 have blank
owner-controlled dispositions. Consent/withdrawal, accessibility, duplicate,
and amendment paths have dedicated fail-closed constraints. Five diagnostic
composite fractions crossed with 72 source scenarios produce 360 rows and 36
within-grid pass/fail transitions, always from integer allocations and raw
Wilson bounds. The fraction grid neither infers code frequencies nor calculates
a scheduled-slot estimand. Five attack cases, deterministic eight-member ZIP
validation, and 128 selected tests passed. Static offline HTML is not public
browser, assistive-technology, or comprehension acceptance. ResolvedCodes=0,
RecruitmentReady=false, and PublicSurfaceEnabled=false.

## Archived Validation Stance

The goal is implementation credibility, not unconditional numerical identity
with another package. FACETS, TAM, sirt, mirt, and this app can differ in
constraints, parameterization, quadrature, optimizer details, latent-variance
treatment, and omitted likelihood constants.

These records can explain historical implementation checks. They do not support
a current product claim of cross-package parity, and they are not required to
run or release the Python application.

## Frozen Fixture Inventory

The historical fixture exporter has been removed from the CLI, Makefile, and
GitHub workflow. Existing sanitized fixtures may be retained temporarily for
compatibility maintenance, but no public command regenerates them and they are
not release evidence. The former generated folder included:

- Python outputs for JMLE RSM, JMLE PCM, MML RSM, MML latent regression, and MML GPCM scenarios.
- sirt-specific person-rater response fixture files: `sirt_rater_facets_response.csv` and `sirt_rater_facets_items.csv`.
- `python_parity_manifest.csv`.
- `cross_package_validation_plan.csv`.
- `cross_package_parameterization_notes.csv`.
- `cross_package_tolerance_policy.csv`.
- `external_reference_documentation.csv`.
- `external_simulation_reference_inventory.csv`.
- `external_simulation_template_inventory.csv`.
- `external_validation_artifact_checklist.csv`.
- `external_validation_report_template.csv`.
- `mfrmr_015_migration_coverage.csv`.
- `r_crosscheck_scaffold.R`.
- `README_external_simulation_templates.md` plus sanitized Python/R/Julia
  template scripts for optional Simulation-style handoff checks.
- notes explaining why exact equality is not expected.

Generated files remain intentionally ignored by Git under
`validation/generated/`.

## Cross-Package Validation Matrix

| Python area | External reference | What to compare | What not to claim |
| --- | --- | --- | --- |
| Data preparation and category support | all packages | row counts, score map, category support, missingness pattern | parameter equality before score support matches |
| JMLE RSM/PCM | FACETS-like fixed-effect workflows and ConQuest-style calibration | centered estimates, ordering, step behavior, fit outliers | exact equality without matching constraints and extreme-score handling |
| MML RSM | TAM `tam.mml.mfr` style faceted MML | convergence, ordering, EAP direction, broad fit behavior | direct likelihood equality unless quadrature, priors, variance, constants, and constraints match |
| Latent regression | TAM latent regression; mirt `mixedmirt` latent-regression concepts | covariate coding, coefficient direction, group-level EAP shifts | TAM/mirt identity while this app uses fixed `population_prior_sd` behavior |
| GPCM | mirt item-level GPCM checks | positive slopes, slope ordering, EAP ordering, item-level response behavior | raw slope equality without a facet-to-item parameterization map |
| Rater facets | sirt `rm.facets` | rater severity ordering, convergence, fit-log availability | exact equality when slope options and constraints differ |
| Plausible values | TAM, mirt, and sirt plausible-value or factor-score routines | distribution means, variances, covariate trends | draw-by-draw plausible-value equality |
| Anchor/linking | FACETS/TAM/ConQuest-style anchor workflows | anchor count, connectedness, hard-constraint drift | linking claims without common-scale evidence |
| Strict marginal diagnostics | TAM/mirt residual and fit diagnostics where applicable | directional flags and sparse-cell warnings | proof of model truth |
| Archived simulation validation sweep | archived mfrmr/Python/Julia/FACETS artifacts | manifest status counts, full-reference summaries, runtime summaries, non-empty validation-input replicates | public runtime dependency on private validation artifacts or exact parity without parameterization notes |
| mfrmr 0.1.5 / 0.1.6 migration coverage | mfrmr package source and package documentation | feature-level support, boundaries, next validation action | one-to-one helper parity or runtime wrapping |

Historical fixtures may contain a machine-readable version of this plan in
`cross_package_validation_plan.csv`.

## Optional External R Handoff

The frozen `r_crosscheck_scaffold.R` was an optional validation helper, not an
app dependency. If a maintainer intentionally audits a retained historical
fixture, the old scaffold can be run from that fixture folder:

```bash
Rscript r_crosscheck_scaffold.R
```

It writes:

- `r_crosscheck_status.csv`: whether each optional package check ran, failed, or was skipped because the package was missing.
- `r_crosscheck_package_versions.csv`: R and package versions for TAM, sirt, mirt, tidyr, and dplyr.
- `r_crosscheck_report.md`: human-readable status plus the Python manifest.
- `r_crosscheck_file_manifest.csv`: files present after the R check.

Then fill `external_validation_report_template.csv` only for rows with actual
supporting evidence. Rows marked `missing`, `error`, or `Not run` should not be
used in public parity claims.

## Default Tolerance Policy

- Exact row counts, rating-category support, score recoding, and missingness flags should match.
- For comparable centered fixed-effect measures, start review at absolute differences greater than about 0.05 logits. Tighten or relax only after documenting package constraints and identification constants.
- For rank ordering of comparable Rasch-family effects, review rank correlations below 0.95.
- For item-level GPCM checks, review rank correlations below 0.90 because the comparison is not the same arbitrary-facet parameterization.
- For latent-regression coefficients, review sign reversals or standardized differences greater than about 0.10.
- For plausible-value distributions, review mean shifts greater than about 0.10 logits or variance ratios outside 0.80 to 1.25.
- Treat log-likelihood differences as interpretable only when category coding, quadrature, priors, constraints, and omitted constants are aligned.
- Do not compare slopes, latent variances, anchored constants, or package-specific nuisance terms as if they were identical parameters without a parameterization map.

## External Reference Roles

- TAM: faceted MML design, latent regression, EAP, and multifacet reference checks.
- mirt: GPCM, EAP/factor-score, plausible-value, and broader IRT diagnostic reference checks.
- sirt: rater-facet, hierarchical rater-model, and plausible-value reference checks.
- mfrmr: functional capability reference for migration completeness.
- Archived simulation artifacts: numerical validation evidence, including
  observed long-form data, repeated engine-refit manifests, runtime summaries,
  full-reference backup outputs, diagnostic spot checks, and validation-input
  replicates. Keep these artifacts out of the public repository unless they are
  sanitized, de-identified, and size-reviewed.

## Official Documentation Touchpoints

- TAM `tam.mml.mfr`: faceted MML design with `formulaA`, `formulaY`, facets, constraints, and variance controls. See the CRAN TAM manual: https://cran.r-project.org/web/packages/TAM/TAM.pdf
- mirt `mirt`: itemtype choices include Rasch/PCM-style and GPCM-style models; `quadpts`, `dentype`, and optimizer choices affect comparability. See the mirt reference: https://philchalmers.github.io/mirt/reference/mirt.html
- mirt `fscores`: EAP is the default factor-score method, and plausible-value support is exposed through `plausible.draws`. See: https://philchalmers.github.io/mirt/reference/fscores.html
- mirt `mixedmirt`: latent-regression inputs are modeled through `lr.fixed` / `lr.random` style arguments. See: https://philchalmers.github.io/mirt/docs/reference/mixedmirt.html
- sirt `rm.facets`: rater-facet models use person-rater rows, rater severity, optional item/rater slopes, EAP factor scores, and modelfit methods. See the sirt pkgdown reference: https://alexanderrobitzsch.github.io/sirt/reference/rm.facets.html
- sirt plausible-value tools: distributional checks should be preferred over draw-by-draw equality. See the CRAN sirt manual: https://cran.r-project.org/web/packages/sirt/sirt.pdf
- mfrmr 0.1.5 / 0.1.6: the package source is used as the migration coverage reference; the generated `mfrmr_015_migration_coverage.csv` and `mfrmr_016_migration_coverage.csv` separate Python support from one-to-one parity claims.

## Reporting Rule

When reporting a cross-package check, archive:

- the generated fixture folder;
- the external package versions;
- the exact R script used;
- the parameterization map;
- the tolerance table;
- the final comparison table; and
- a short note explaining which differences are expected.

Do not report exact cross-package parity unless those artifacts are included.

## Latest FACETS Pilot Status

The current FACETS 4.5.0 qualification bundle is
`facets_450_pilot20_table8_table13_v2_20260811/`. It replays the immutable,
independently generated 20-replicate pilot bytes through FACETS and records the
current Python comparison results under separate hashes. It archives Table 7
parameter agreement, Table 8 category/threshold mappings, and
the canonical arranged-by-N Table 13 Rater-by-Task bias mapping. Read
`FACETS_COMPLEMENTARY_WORKBENCH_AUDIT_20260811.md` and
`facets_compatibility_matrix.json` before making a compatibility claim. The
machine-readable disclosure is `facets_pilot20_amendment_20260811.json`.

The Table 8 adapter deliberately uses a second `Umean=0,1,2` report for the
one-decimal category Outfit field because the six-decimal fixed-width rendering
truncates that token. The six-decimal primary report remains authoritative for
measures and thresholds. Missing Table 13 cells are retained in the expected-cell
audit as unavailable, never converted to negative findings. The pilot includes
nondegenerate Table 13 flags, but remains parser and pipeline qualification; it
does not provide confirmatory false-positive rate, power, coverage, or estimator
superiority evidence.

See `R_CROSSCHECK_STATUS.md` for the latest archived smoke run. That file
records whether the generated fixture could be fitted by the installed TAM,
mirt, and sirt packages, while keeping the same non-parity interpretation rule
used throughout this validation directory.

See `SIMULATION_REFERENCE_STATUS.md` for the archived validation-artifact
inventory that should guide external-data numerical validation without bundling
private or large datasets into this public app repository.

Historical fixtures may still contain sanitized Python/R/Julia templates.
Current Downloads and demo archives do not expose them; Python-native
reproduction assets are the supported route.

## Known assignment-mechanism oracle

`KNOWN_ASSIGNMENT_MECHANISM_PILOT_20260811.md` documents the frozen
finite-state qualification of the score-blind, degree-conditioned exponential
assignment generator. Its retained bundle contains all 90 exact states at
three `gamma` values, independent-chain diagnostics, state frequencies,
oracle comparisons, materialization audits, hashes, and a machine-readable
PASS assessment.

This evidence applies only to the assignment generator and equal-context row
materializer. It contains no FACETS run and no estimator-bias result because
responses are intentionally downstream of assignment. A later response study
must preserve this separation, use FACETS 4.5 only for same-estimand JMLE
calibration, and keep the ordinary two-decimal FACETS fit display out of raw
fit calculations.

## Known-assignment large-design and response screening

`KNOWN_ASSIGNMENT_LARGE_DESIGN_CALIBRATION_20260811.md` retains the failed
edge-pair Metropolis and Person-pair heat-bath finite-chain settings, then
documents the qualified four-Rater exact-DP route and outcome-blind selection
of `|gamma|=0.8`.

`KNOWN_ASSIGNMENT_RESPONSE_SCREENING10_20260811.md` connects that generator to
one fixed 80-Person heterogeneous-threshold PCM. The final v2 aggregate passed
all registered denominator, constraint, FACETS/Python JMLE, and endpoint-shape
gates. Base R 4.5.1 independently reproduced all ten estimator-by-dose Rater
RMSE contrasts to `1.39e-17` maximum difference. The only screening intervals
excluding zero were normal-person MML Rater RMSE/MAE contrasts. This evidence
is conditional on one fixed Person vector and cannot support estimator ranking,
real missingness-mechanism inference, or confirmation.

FACETS measure/threshold recovery uses the high-precision primary report;
ordinary rounded fit fields remain a separate evidence layer and are not raw
inputs to these recovery statistics.

`KNOWN_ASSIGNMENT_DENSE_DP_QUALIFICATION_20260811.md` records the next
implementation-only gate. The vectorized dense-margin partition function
matched both the 90-state oracle and the retained recursive 80-by-4 exact-DP
results, with maximum log-normalizer difference `2.84e-14`. Its maximum
partition construction time was 0.202 seconds and its storage boundary is an
explicit fail-closed cell cap. No responses or FACETS outputs entered this
qualification; it enables, but does not replace, a fresh multi-Person-vector
preflight before confirmatory registration.

`KNOWN_ASSIGNMENT_MULTIVECTOR_PREFLIGHT_STATUS_20260811.md` records the frozen
four-vector preflight and its non-destructive recovery. All 48 original attempt
units were preserved. A clean short-path, visible-mode FACETS supplement passed
12/12 same-input calibration gates with zero retries, and original versus
supplemental Python JMLE replay differed by at most `4.44e-16`. The derivative
aggregate passed every operational gate; PF1--PF4 all met their registered n=4
direction rules, and R 4.5.1 reproduced all values within `1.67e-16`. This is a
nonconfirmatory advancement signal only. The four vectors cannot enter the
separate confirmation, estimator ranking remains prohibited, and ordinary
two-decimal FACETS fit displays were not used as raw inputs.
