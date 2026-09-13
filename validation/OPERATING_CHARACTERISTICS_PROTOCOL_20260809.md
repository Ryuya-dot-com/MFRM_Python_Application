# MFRM operating-characteristics protocol (2026-08-09)

## Status and decision purpose

This protocol governs the repeated-simulation lane for the standalone Python
application and future matched adapters for the in-development `mfrmr` 0.2.3,
TAM, immer, and sirt implementations. The retained `smoke` profile verifies
orchestration, denominators, negative controls, and exports. It is not evidence
of stable false-positive rate, power, coverage, sample-size sufficiency, or
cross-package equivalence.

The study is deliberately repository-only. No simulation rate or estimator is
shown in the public Streamlit application until the promotion gates below are
met. In particular, optimizer completion is not interchangeable with
convergence, inference readiness, or an eligible decision.

## ADEMP specification

### Aims

1. Estimate false-positive rate and power of prespecified conditional-bias
   decisions without converting failed, sparse, or unavailable analyses into
   negative decisions.
2. Estimate parameter bias, RMSE, MAE, SE availability, and conditional-Wald
   coverage on the correct identification scale.
3. quantify sensitivity to planned sparsity, missing observations, anchor
   share, contaminated anchors, and floating-point/display boundaries.
4. Compare matched estimands across Python JMLE and, where structurally
   eligible, Python CMLE, mfrmr, TAM, immer, and sirt. Package rankings and
   equality claims are outside scope.

### Data-generating mechanisms

The first matrix uses a polytomous RSM with Person, Rater, Task, and Criterion
effects, four score categories, and a focal `R01 x T01` local interaction.
Null and +0.60-logit alternatives share latent draws and uniforms through a
recorded common-random-number seed group. The following design stresses are
paired by truth state:

| Design | Persons | Raters | Raters/person | Missingness | Anchors | Purpose |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| `balanced_small` | 30 | 4 | 4 | 0% | 0% | small balanced baseline |
| `balanced_large_anchors` | 80 | 4 | 4 | 0% | 50% of Raters | correct-anchor control |
| `sparse_missing` | 18 | 8 | 1 | 35% random | 0% | low-count and ineligible-cell control |
| `anchor_drift` | 80 | 4 | 4 | 0% | 50%, shifted +0.25 | paired contaminated-anchor control |

The two large anchor designs share generated ratings and true parameters; only
the supplied anchor targets differ. Anchor percentage is descriptive, not a
universal adequacy threshold. Later registered extensions should vary graph
connectivity, score-dependent missingness, category rarity, extreme Persons,
facet variance, PCM steps, and effect magnitudes without overwriting this
initial matrix.

### Estimands

- Decision operating characteristics: false-positive rate under `TruthBias=0`
  and power under `TruthBias!=0`, separately for Holm, practical-magnitude,
  conjunctive strong, and any non-sparse flag rules.
- Failure characteristics: fit-return, convergence, inference-readiness, and
  analysis-eligibility rates plus failure stage and normalized reason.
- Recovery: bias, RMSE, MAE, and Monte Carlo SE for each facet.
- Uncertainty: SE availability and nominal 95% conditional-Wald coverage among
  included returned estimates with finite positive SE.
- Presentation sensitivity: conclusion changes between decisions based on raw
  floating-point values and displayed rounded values. Raw decisions are
  authoritative.

For unidentified facet locations, recovery errors are mean-aligned within
facet and replicate. When fixed anchors identify a facet's absolute location,
absolute errors are used; re-centering would hide anchor contamination. These
comparison scales must remain separate in every table and figure.

### Methods

The baseline adapter fits the Python application's RSM JMLE and routes its
existing fit, DFF conditional-bias, sparse-cell, and decision-stability
contracts into an engine-neutral schema. The implemented CMLE increment now
applies the same frozen generated datasets and eligibility fields to:

- Python exact CMLE and `immer::immer_cml()` only for their common eligible
  unanchored RSM structural-calibration estimand;
- mfrmr 0.2.3 under explicit JML mode and identification labels;
- TAM MML under the additive Criterion-item + Rater + Task + common-step RSM
  surface, with an estimated normal Person population and separate quadrature
  modes;
- sirt MML under a separately labelled virtual-item PCM mapping and common
  estimated Person population; comparisons with Python JMLE are estimator
  sensitivity, not parity.

Unsupported specifications remain typed `unsupported`; they are not errors and
are not imputed. Likelihoods, raw coefficients, Person estimates, or slope
parameters are not compared until their parameterization and conditioning
contracts match.

### Performance measures

Every binary rate includes its eligible denominator, unavailable/ineligible
count, Monte Carlo SE, and Wilson 95% interval. Recovery summaries include the
number attempted, included, unavailable, and coverage-eligible. Results are
stratified by condition, truth, engine, estimator, parameter type, and
comparison scale; null and alternative conditions are never pooled.

## Replication profiles and precision

- `smoke`: 2 replicates/condition. Orchestration only.
- `pilot`: 20 replicates/condition. Runtime and design debugging only.
- `study`: 100 replicates/condition. Initial screening depth, not automatic
  validation.

The final budget should be chosen from the registered precision target and
observed failure rate, not from these labels alone. As a reference, about 500
eligible null replicates give a binomial Monte Carlo SE near .01 when the true
false-positive rate is .05; more are needed when subgrouping, failures, or
tail precision matter. Budget expansions must preserve the existing manifest
identity and append replicates rather than silently replace prior seeds.

The prospective details are machine-readable in
`operating_characteristics_precision_plan_20260809.json`. It fixes the raw
decision contract, 500 eligible-result confirmatory reference target,
condition-specific Wilson-lower-bound attempt inflation, 2,500-attempt cap,
no outcome peeking, pilot/confirmatory non-pooling, and a public status of
`Withheld`. `audit_manifest_extension()` operationalizes extension as exact
preservation of every prior RunId/row/seed plus consecutive higher replicate
indices; a condition-major expanded table is not incorrectly treated as a
literal row prefix.

### Post-registration 20-replicate pilot audit

The isolated Python pilot under
`operating_characteristics_pilot20_20260809/` retained 160 attempted fits and
passed all 10 artifact-integrity checks. It is not promoted to performance
evidence. All fits returned, 156 optimizer flags reported convergence, and 124
focal decisions met the registered app eligibility rule. The reconstructed
terminal-gradient sup norm nevertheless exceeded `1e-4` in all 160 runs
(range `0.00013306`--`0.102327`). This post-pilot threshold is retained as a
blocker sensitivity, not retroactively substituted as the primary rule. A
separately labelled, prospectively specified strict Python JMLE mode must be
qualified before a confirmatory manifest is fixed.

The sparse null and alternative conditions each yielded 4/20 eligible focal
decisions. Their registered Wilson-lower-bound calculation requires 6,200
attempts to target 500 eligible decisions, so both fail the 2,500-attempt cap
and require design/estimand revision rather than brute-force replication.
Among 39,680 row-level fit-statistic audits, 14 lay inside a 3-decimal display
boundary and five raw/display labels disagreed. Raw values remain
authoritative. The frozen v1 plan also contains an unresolved planning
conflict: its continuous-bias section invokes a pilot sample variance while
its budget rule forbids pilot bias outcomes. The v1 file remains unchanged;
continuous-outcome precision requires a prospective amendment or independent
variance-planning source.

### Strict JMLE follow-up and structural identifiability gate

The pilot's optimizer flag and reconstructed-gradient conflict was resolved
prospectively, not by changing the registered pilot rule. Stage A retained the
original 16 terminal vectors, verified the analytical gradient, and confirmed
that none met the raw `1e-4` sensitivity. Stage A2 froze two continuation
candidates before comparison. The priority L-BFGS-B precision candidate passed
all 16 contracts; BFGS was retained as sensitivity evidence but was not selected
because 14/16 fits ended in precision loss and its coordinate movement was much
larger. Revised Stage B2 then applied only the selected candidate to the frozen
160-run manifest. It passed 160/160 numerical contracts at `1e-4`, 126/160 at
`1e-5`, and 43/160 at `1e-6`, with no raw-versus-three-significant-digit gate
disagreement.

That numerical result is necessary but not sufficient. The separately frozen
movement audit found up to 57.02 logits of free-coordinate movement for tiny
likelihood gains. The exact eta-design-rank audit found nullity 7 and eight
Person-Rater components in all 40 sparse runs, while all 120 balanced/anchor
runs had nullity zero and one component. The sparse null space was concentrated
in Person and Rater coordinates. Accordingly, sparse JMLE remains structurally
nonidentified regardless of optimizer success or replication count.

The application-level integration contract authorizes only a pre-optimizer
audit of the free Person/facet eta design, machine-readable evidence, separate
`Converged` and `InferenceReady` status, and withholding of conditional-bias
output when rank deficient. It does not authorize optimizer changes, estimate
replacement, automatic JMLE-to-MML switching, or a public performance claim.
The post-change bridge matched rank, nullity, Person-Rater connectivity, audit
scope, and readiness on all 160 retained RunIds. This audit does not qualify
PCM/GPCM step or GPCM slope identification, which remain governed separately.

The retained decision trail is:

- `operating_characteristics_strict_jmle_plan_20260809.json` and
  `operating_characteristics_strict_jmle_smoke_20260809/`;
- `operating_characteristics_strict_jmle_a2_plan_20260809.json` and
  `operating_characteristics_strict_jmle_a2_smoke_20260809/`;
- `operating_characteristics_strict_jmle_b2_plan_20260809.json` and
  `operating_characteristics_strict_jmle_b2_20260809/`;
- `operating_characteristics_jmle_movement_plan_20260809.json` and
  `operating_characteristics_jmle_movement_20260809/`;
- `operating_characteristics_identifiability_plan_20260809.json` and
  `operating_characteristics_identifiability_20260809/`; and
- `jmle_identifiability_integration_plan_20260809.json` and
  `jmle_identifiability_integration_20260809/`.

### JMLE Person score-boundary output gate

The movement outcome motivated a separate post-hoc evidence audit and
prospective application-change contract in
`jmle_extreme_score_integration_plan_20260809.json`. It defines a boundary
Person solely by exact retained all-minimum or all-maximum integer scores. It
does not use the magnitude or displayed rounding of a terminal theta value.
The current-source bridge matched 8,320/8,320 retained Person rows. Four
boundary Persons occurred among the 120 structurally identified runs and
accounted for every >=1-logit movement in that subset; its maximum non-theta
movement was 0.000293 logits. The 40 structurally deficient sparse runs had 53
additional boundary rows, but their interior Person rows also remain unready
under the preceding rank gate.

The authorized output change retains the existing optimizer/constraint
`Estimate` for reproduction, sets `FiniteJMLEEstimate=False` and
`PersonInferenceReady=False`, and withholds `ReportableEstimate` for a boundary
Person. Fit-level structural readiness and Person-measure readiness remain
separate. Bilingual UI warnings and quick/full/publication evidence tables
carry the same audit. No optimizer controls or estimates changed, and no finite
extreme-score correction, estimator switch, precision-polish integration, or
public performance claim is authorized. Retained evidence is under
`jmle_extreme_score_integration_20260809/`.

### Fixed-calibration Warm WLE comparison

The first finite-score comparison is prospectively governed by
`fixed_calibration_wle_plan_20260809.json` and the transparent pre-result hash
amendment `fixed_calibration_wle_plan_amendment_20260809.json`. It implements
the one-dimensional Warm adjusted score
`U(theta) + 0.5 I'(theta) / I(theta)` for generic category intercepts and theta
coefficients. This is a fixed-calibration WLE Person estimand, not a finite
JMLE MLE and not MML/EAP. The Python core matched `TAM::tam.mml.wle2` for all
16 frozen RSM/PCM/GPCM Person fixtures; maximum absolute theta and SE
differences were `1.68e-11` and `1.18e-11`.

The separately registered downstream plan replayed only the 120 structurally
identified Stage-B2 runs and held polished facets and steps fixed. It compared
7,600 Person scores, including four exact-extreme patterns. Median absolute
interior movement was about 0.009 logits, whereas exact-extreme movement was
23.72--25.19 logits. Raw fit-zone decisions changed in 51/17,360 element-metric
comparisons, all on the Person facet. At `.3g`, 31 JMLE and 23 WLE fit decisions
disagreed with the unrounded decision. The focal Holm and combined strong-bias
rules had zero switches; the practical `|bias| >= 0.50` rule switched once.
The 40 rank-deficient sparse runs were explicitly withheld. These results
qualify a numerical research implementation and sensitivity evidence only;
they do not authorize automatic correction, estimator switching, facet-
uncertainty claims, or a public WLE control.

The follow-on exact-CMLE-to-WLE research bridge is further governed by
`cmle_wle_calibration_sensitivity_plan_20260810.json`. At primary covariance
scale, 2,000 seeded free-coordinate draws per RSM/PCM model returned all 32,000
draw-Person scores. Median calibration-draw SD was `0.0616`--`0.0801` logits by
model, with a maximum `0.518` logits and a maximum ratio `0.427` to the
conditional WLE SE. Zero-scale, same-seed, covariance-scale monotonicity,
exact-extreme, PSD, symmetry, and work-limit checks are retained. These are
sensitivity outputs only: same-sample calibration/Person dependence is not
represented, draw quantiles are not confidence intervals, and quadrature is
not an inference-qualified total SE. The first-read projection therefore
keeps numerical research readiness separate from caution and withheld states.

## Failure, multiplicity, and floating-point rules

1. Every attempted condition-replicate-engine fit receives a row, including
   exceptions and unsupported estimands.
2. Failed, non-converged, inference-unready, missing-SE, and sparse focal-cell
   results are visible in accounting and never counted as negative decisions.
3. The Holm-adjusted probability is computed over the full retained bias-cell
   family before selecting the focal cell.
4. Comparisons use finite unrounded numbers. Display rounding can annotate but
   cannot determine the primary conclusion.
5. Exact endpoint semantics and numerical/display boundary bands come from
   `mfrm_app.decision_stability`; conclusion mismatches are retained as an
   operating characteristic.

## Reproducibility and retained evidence

Each run retains the condition table, condition-replicate manifest, SHA-256
manifest fingerprint, base seed, seed-coupling label, engine and estimator,
raw run accounting, failure reasons, parameter rows, decision rows, aggregate
tables, first-read summary, and figures. The Python smoke profile is reproduced
from the repository root with:

```bash
MPLCONFIGDIR=/tmp/mfrm_matplotlib_cache \
  python3 validation/operating_characteristics_pilot.py --profile smoke
```

The runner has an explicit maximum-run guard. Generated datasets and adapter
outputs should be content-fingerprinted before cross-engine results are joined.

The Python runner now publishes `generated_ratings.csv`,
`generated_facet_truth.csv`, and `generated_anchors.csv` through atomic staging,
plus per-run `DataId`/`FitInputId` identities and byte-level SHA-256 file hashes.
The R bridge must validate those bytes rather than recreate the observations
from the recorded seed, because NumPy and R random-number implementations are
not assumed to produce identical draws. A regenerated Python bundle removes
derived R bridge outputs so stale validation cannot survive a changed input.

After installing the intended mfrmr development snapshot into an isolated
temporary library, validate and import the exact bundle with:

```bash
Rscript validation/operating_characteristics_bridge.R \
  --input validation/operating_characteristics_20260809 \
  --output validation/operating_characteristics_20260809 \
  --mfrmr-lib <isolated-mfrmr-0.2.3-library> \
  --mfrmr-git-head <source-commit> \
  --mfrmr-source-state <snapshot-state>
```

`BRIDGE_RESULTS.md`, `bridge_validation_r.csv`, and
`bridge_engine_availability_r.csv` establish input readiness and package/code
identity only. They do not contain R estimates and cannot support an agreement
claim.

The first fit adapter uses mfrmr JMLE in two non-pooled modes:

- `MFRMR_JML_MATCHED_CONTROL`: `maxit=160`, `reltol=1e-6`, matching the
  Python smoke runner's requested numerical controls. This is a stopping-rule
  sensitivity condition, not the primary mfrmr readiness mode.
- `MFRMR_JML_STRICT`: `maxit=500`, `reltol=1e-9`, requiring mfrmr's terminal
  gradient and inference-readiness contract to pass.

Both modes use the exact retained long data, `Person` plus the additive
`Rater`/`Task`/`Criterion` RSM structure, declared 0-3 category support, and
the exact retained hard anchors. The adapter checks anchored estimates against
their supplied targets to `1e-10`, retains all pre-fit rank rejections, and
computes the same focal `R01 x T01` conditional plug-in bias family with Holm,
BH-plus-|t|, practical 0.50-logit, and sparse-cell rules. Only inference-ready,
anchor-valid, bias-ready, non-sparse rows are analysis eligible.

Run it after the bridge with:

```bash
Rscript validation/operating_characteristics_mfrmr.R \
  --input validation/operating_characteristics_20260809 \
  --output validation/operating_characteristics_20260809 \
  --mfrmr-lib <isolated-mfrmr-0.2.3-library> \
  --mfrmr-git-head <source-commit> \
  --mfrmr-source-state <snapshot-state>

MPLCONFIGDIR=/tmp/mfrm_matplotlib_cache \
  python3 validation/operating_characteristics_compare.py
```

Python bundle regeneration invalidates both bridge-only and derived mfrmr
comparison artifacts. Numerical agreement is summarized only after readiness;
finite but withheld matched-control fits remain visible and are never pooled
with strict results.

The exact-CMLE comparison is reproduced after the R bridge with:

```bash
python3 validation/operating_characteristics_cmle.py \
  --input validation/operating_characteristics_20260809 \
  --output validation/operating_characteristics_20260809

Rscript validation/operating_characteristics_immer_cmle.R \
  --input validation/operating_characteristics_20260809 \
  --output validation/operating_characteristics_20260809

MPLCONFIGDIR=/tmp/mfrm_matplotlib_cache \
  python3 validation/operating_characteristics_cmle_compare.py \
  --input validation/operating_characteristics_20260809 \
  --output validation/operating_characteristics_20260809
```

Both adapters estimate the same unanchored additive RSM conditional
likelihood, with Person conditioned out and identical sum-to-zero free
coordinates. A zero-weight `immer` sentinel declares every virtual unit's
0--3 support without contributing a sufficient statistic or score frequency;
`nullcats="keep"` prevents observed zero-count categories from being assigned
zero probability. Python exact moments and an independently reconstructed
R/immer conditional-information matrix must agree on the pre-fit rank gate.
Anchor-requesting runs remain useful unanchored
numerical-parity references but are ineligible for the requested anchored
condition; local Rater x Task bias is outside this additive CMLE estimand.

The smoke result retained in `CMLE_IMMER_RESULTS.md` found 12 jointly returned
RunIds (8 unique eligible rating datasets), 96 matched free coordinates, a
maximum estimate difference of `5.63e-8` logits, and a maximum conditional-
loglikelihood difference of `4.77e-12`. However, `immer` optimizer code 0 did
not imply the prespecified inference readiness: terminal-gradient readiness
was 12/12 at `1e-4`, 7/12 at the primary `1e-5`, and 0/12 at `1e-6`.
These sensitivity counts are retained rather than selecting a favorable
threshold after seeing the fit. The reconstructed analytical `immer` gradient
was also checked against centered finite differences: the largest component
discrepancy was `3.88e-7`, smaller than every returned RunId's distance to the
primary boundary. The two pre-fit rank audits agreed for 16/16 RunIds; four
sparse RunIds were rejected before optimization at conditional-information
rank 5/12 (nullity 7).
Because additive CMLE has no local Rater x Task interaction parameter, the
paired null/+0.60 data are also summarized as coordinate leakage rather than a
bias decision. The balanced-small smoke change was 0.0612 logits on average
and at most 0.172 logits across free coordinates; these descriptive values are
not false-positive, power, or interaction-recovery estimates.

The TAM MML sensitivity adapter is reproduced after the R bridge with:

```bash
Rscript validation/operating_characteristics_tam.R \
  --input validation/operating_characteristics_20260809 \
  --output validation/operating_characteristics_20260809

MPLCONFIGDIR=/tmp/mfrm_matplotlib_cache \
  python3 validation/operating_characteristics_tam_compare.py \
  --input validation/operating_characteristics_20260809 \
  --output validation/operating_characteristics_20260809
```

`TAM::tam.mml.mfr()` uses Criterion as the item dimension and Rater and Task as
facets under `~ item + rater + task + step`. This reproduces the registered
additive RSM surface directly, but Person parameters are integrated under an
estimated normal population. Python-JMLE/TAM comparisons are therefore
cross-estimator sensitivity, not likelihood or coefficient parity. Requested
Rater anchors are fixed through free `xsi` indices discovered by a constrained
design preflight; a requested derived last-level anchor fails closed.

The primary mode uses 61 equally spaced nodes on -8--8 and the sensitivity mode
uses 21. Both use `convD=1e-6`, `conv=1e-5`, `convM=1e-6`, 20 M steps, and a
1,200-iteration cap. TAM's source loop keeps variance/regression convergence as
persistent latches, so exit before the cap alone is insufficient. The adapter
also retains the progress-formatter function identity, parses terminal item,
regression, and variance changes, and combines these with the exact deviance
history. Parsed progress values have text-format precision; the primary source
loop remains authoritative for its item/deviance gate, and alternative
threshold counts are labelled sensitivity evidence rather than gradients.

All 32 q21/q61 fits returned and passed that convergence audit; 30 were analysis
eligible. Two q21 sparse fits landed at stored variance `0.0010000001` against
the configured `0.001` lower bound. The variance gate uses a frozen `1e-8`
numerical band, so a naive exact comparison cannot convert the boundary fit to
ready evidence. All 16 primary q61 fits were eligible. At parsed terminal-
progress thresholds, 32/32 pass `1e-4`, 32/32 pass the primary `1e-5`, and 0/32
pass `1e-6`; the latter counts do not claim more precision than the formatter.

Among jointly eligible q21/q61 rows, maximum differences were 0.02224 logits
for Rater, 0.3814 for Person EAP, and 0.6990 for cumulative response-surface
difficulty. Returned variance-boundary rows remain under a separate all-
returned scope. Balanced-small Python-JMLE/TAM-MML Rater MAE was 0.01383
logits; TAM/sirt MML Rater MAE was 0.001513. These small descriptive differences
do not erase Person-treatment or item-threshold differences.

The paired contaminated anchors expose constraint mechanics. The two fixed TAM
Raters receive +0.25 exactly, while the two unanchored Raters compensate by
about -0.25 under the Rater sum-zero constraint. The local Rater x Task
interaction remains omitted and is summarized only as facet leakage. Derived
last-level coverage is withheld because TAM 4.3-25 constructs the expanded SE
from a diagonal matrix of free-xi SEs without retaining full covariance.
Information criteria are retained for within-TAM audit only.

`tam_first_read_summary.csv` orders scope, run accounting, estimator boundary,
MML structure, quadrature, floating variance boundary, anchor contamination,
sparse identification, and public status. It is a future UI projection;
`PublicSurfaceEnabled` remains false.

The sirt MML sensitivity adapter is reproduced after the R bridge with:

```bash
Rscript validation/operating_characteristics_sirt.R \
  --input validation/operating_characteristics_20260809 \
  --output validation/operating_characteristics_20260809

MPLCONFIGDIR=/tmp/mfrm_matplotlib_cache \
  python3 validation/operating_characteristics_sirt_compare.py \
  --input validation/operating_characteristics_20260809 \
  --output validation/operating_characteristics_20260809
```

`sirt::rm.facets()` is not treated as a JMLE/CMLE parity target. The adapter
maps each Task x Criterion combination to a virtual PCM item, retains Rater
severity, fixes requested Rater anchors through `b.rater.fixed`, and estimates
a normal Person population on a fixed -8--8 grid. The primary 61-point grid
and a 30-point sensitivity grid are never pooled. A fit is analysis eligible
only when it returns, exits before the 1,200-iteration cap, retains finite
Rater/item-threshold uncertainty and Person EAPs, passes the quadrature-tail
mass check, and reproduces hard anchors exactly. `rm.facets()` does not return
its final global-parameter change, so reaching the cap remains unresolved even
when the returned deviance change is small.

The retained smoke evidence has 32/32 returned fits and 30/32 eligible fits;
one sparse +0.60 RunId reached the cap in both grids. Quadrature sensitivity is
reported separately for all jointly returned and jointly eligible rows. Among
eligible runs, maximum q61-q30 differences were 0.001662 logits for Rater
severity, 0.001387 for mean-centered virtual-item location, 0.05011 for raw
virtual-item location, 0.07862 for Person EAP, and 0.003539 for population SD.
Threshold counts always use unrounded differences and publish their exact
denominator. `sirt_first_read_summary.csv` orders the same evidence as scope,
run accounting, estimator boundary, quadrature, sparse identification, anchor
contamination, and public-surface status. It is a future UI projection, not an
enabled UI route; `PublicSurfaceEnabled` remains false.

The Python comparison is intentionally asymmetric. Balanced-small Rater
estimates were close descriptively, but sparse sirt fits can be identified by
the estimated common Person distribution when the observed Rater graph fails
the strict mfrmr and exact-CMLE rank gates. This is assumption-based MML
identification, not evidence that the observed design became connected. The
paired anchor-drift input also shifts sirt Rater estimates by essentially the
injected +0.25 logits, which is contamination transmission rather than anchor
robustness. Local Rater x Task bias is omitted by this specification and is
reported only as leakage into Rater main effects. Fixed-anchor SEs and sirt
4.2-133 information criteria are not interpreted; the latter are withheld
because the inspected internal parameter-count path subtracts the numeric
Rater-centering mode from the count.

## Promotion gates

Public UI exposure requires all of the following, reviewed together:

1. deterministic smoke reproduction and passing contract tests on supported
   Python platforms;
2. negative controls that visibly reduce eligibility or worsen recovery in the
   intended direction without manufacturing false precision;
3. a frozen study manifest with condition-specific Monte Carlo precision and
   complete failure accounting;
4. matched-estimand or matched-response-surface cross-engine results for
   mfrmr/TAM/immer plus separately registered sirt MML sensitivity, with
   unsupported or non-comparable cells explicitly retained;
5. sensitivity analyses for floating-point thresholds, missingness, sparsity,
   anchor share and contamination, and multiplicity;
6. an independent review of labels, downloads, accessibility, and the rule
   that no automatic sample-size or anchor-share recommendation is inferred;
7. a one-click result surface that leads with readiness, denominators, and
   limitations, while keeping full tables and provenance downloadable.

Until then, `first_read_summary.csv` must report the public application surface
as `Withheld`. Passing a replicate-count threshold alone never changes that
status.
