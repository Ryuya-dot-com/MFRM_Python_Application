# Anchor User Guidelines (MFRM Streamlit App)

## 1) Purpose
Use anchors when you need stable scale linkage across runs, cohorts, prompts, or administrations.

Fixed anchors are one route to a shared scale. Joint calibration can instead
combine datasets through shared observations, without supplying fixed measures.
See the [two-school lesson](#two-school-recordings) below before choosing a route.

Typical use cases:
- Reusing previously calibrated rater/task measures.
- Linking two administrations with only partial overlap.
- Stabilizing known reference elements.

## 2) Which option to use
- Anchor table: fix specific element measures to exact values.
- Group anchor table: fix a group's mean measure; its individual levels may differ. A stable group mean must be justified by the linking design.
- Dummy facet: fix all levels of a facet to 0. Differences between these levels are not estimated or adjusted for.
- Non-centered facet: omit sum-to-zero centering for the selected facet. The full model still needs identifying constraints; inspect the identification audit.

## 3) Input formats
### 3.1 Anchor table
Required columns:
- Facet
- Level
- Anchor

Example:
```csv
Facet,Level,Anchor
Rater,R01,0.50
Rater,R02,-0.20
Task,T01,0.30
```

### 3.2 Group anchor table
Required columns:
- Facet
- Level
- Group
- GroupValue

Interpretation:
- For each group, the app enforces:
  - sum(level measures in group) = GroupValue * group size

Example:
```csv
Facet,Level,Group,GroupValue
Rater,R01,Core,0.00
Rater,R02,Core,0.00
Rater,R03,Core,0.00
Rater,R04,New,0.20
Rater,R05,New,0.20
```

## 4) Matching rules (critical)
- Facet and Level must exactly match the values in your analysis data.
- Matching is string-based.
- Leading/trailing spaces can break matching.
- Case mismatches can break matching.

Rows with an unknown facet or level are excluded from the constraints and
recorded in the anchor audit. Missing or non-finite anchor values are also
excluded. The default `warn` policy retains these issues for review; `error`
stops estimation when the audit contains warnings or errors. Always inspect
the audit and applied-anchor counts before interpreting a linked scale.

## 5) Practical workflow
1. Run once without anchors and inspect fit and counts.
2. Pick stable reference elements (reasonable fit, enough observations).
3. Build anchor table or group anchor table.
4. Rerun with anchors.
5. Confirm results:
- The anchor/linking audit: issues and applied-anchor counts by facet.
- Exported `anchor_audit_summary.csv` and the report's Anchor status column.

## 6) Anchor status in report
- A: directly anchored.
- G: group-anchor constrained.
- X: group-anchor constrained and extreme-only element.

## 7) Recommended guardrails
- Start with a small number of high-quality anchors.
- Avoid anchoring many weak or misfitting elements.
- Use a group anchor only when the shared scale and group mean are defensible; grouping does not establish stability or equivalence by itself.
- Recheck fit after applying anchors.
- Keep a versioned record of anchor files used in each run.

## 8) Common mistakes
- Wrong Facet name (for example, `Raters` vs `Rater`).
- Level label mismatch (`R1` vs `R01`).
- Hidden spaces in pasted text.
- Mixing scales from incompatible runs without common elements.

## 9) Notes for estimation mode
- In JMLE, person parameters are explicitly estimated with constraints.
- In ordinary MML engines (EM, direct, hybrid, auto), fixed and group-mean
  constraints are implemented for non-Person facet locations, such as Rater,
  Task and Criterion. This does not provide a general interface for fixing all
  thresholds or GPCM slopes. Match the population model and scale origin as
  well as the anchor values when comparing calibrations.
- In MML, person estimates are EAP scores under a population model. Fixed
  Person anchors, Person group anchors and a dummy Person facet are unsupported
  and stop estimation. Earlier versions could record these inputs without
  enforcing them in the marginal likelihood; rerun an affected analysis with
  supported constraints before interpreting its linking claims.
- The explicit joint free-population-SD engine excludes all fixed/group anchors
  and dummy facets; free-SD inference remains withheld.
- Combining ratings under shared Person IDs is joint calibration and does not
  require fixed Person anchors. MML can integrate each shared person's ratings
  together. Fixing a common recording's Person value is a different operation,
  supported by this app's JMLE route, not its MML route.

Anchoring is not inherently restricted to JML. For example,
[TAM's MML documentation](https://alexanderrobitzsch.github.io/TAM/reference/tam.mml.html)
provides `xsi.fixed` in `tam.mml` and `tam.mml.mfr` to fix structural parameters.
This establishes method support, not numerical equivalence with this app.

### DFF when reviewing an anchor design

The app's **Bias / Interaction** view provides conditional facet-pair estimates,
an interactive heatmap, descriptive distributions and screening CSVs. DFF means
that a particular combination departs from the additive model; overall rater
severity alone is not DFF. The plotted bias sign follows the score direction:
positive means higher scores than predicted, negative means lower scores.

Start with the heatmap and its four-column table (two facet levels, rating count
and bias). Open the reading guide for five questions on overall severity, signs,
missing evidence, inference and school assignments. A blank is not zero bias;
counts refer to ratings, not necessarily distinct persons. Color scales can
change between plots, so compare numeric values rather than color intensity.
Detailed calculations and CSV downloads are available in a collapsed section.
After fitting the two-school example, the lesson's DFF button opens this view;
Rater × Criterion is the initial pair. Changing the display does not refit the
model. If a different pair was computed, select the required pair in the analysis
settings and explicitly rerun.

The current calculation fixes estimated calibration values and, for MML, also
fixes person EAP scores. It is not a joint marginal interaction fit. Error rates,
coverage and calibration uncertainty remain unqualified; adjusted probabilities
and flags do not establish DFF or absence of bias. The separate classical DIF
screen uses an external Person-group file and observed-score matching. It is
not interchangeable with this MFRM interaction analysis and currently provides
tables rather than group-specific response-curve plots.

In the two-school example, Rater × Criterion cells are observed, but each local
student group is scored only by its own school's raters. Shared reference
recordings link the panels; they do not create the missing cross-school student
ratings needed to assess differential treatment of A versus B students by every
rater. Plan that overlap before making a Rater × School DFF claim.
See the FACETS 64-bit manual's [DFF overview](https://www.winsteps.com/facetman64/introduction.htm)
and [bias contrasts](https://www.winsteps.com/facetman64/table14.htm) for the
concept; this app's inferential qualification is assessed separately.

### Group information and observed assignments

The two-school lesson includes a Person/School CSV. A01–A20 belong to A and
B01–B20 to B. C01–C10 have a blank School: they link the calibration but are not
members of either student group. Do not assign references to an arbitrary
school or duplicate them in both groups. The lesson's count table is therefore:

| Rater | School A students | School B students |
|---|---:|---:|
| A_R1 | 20 | 0 |
| A_R2 | 20 | 0 |
| B_R1 | 0 | 20 |
| B_R2 | 0 | 20 |

Each cell counts distinct persons, not the three criterion ratings per person.
Shared references connect the panels but do not supply those missing student
ratings. Depending on the intended contrast, plan cross-school ratings of
performances spanning the relevant ability range and criteria. Complete cells
alone do not establish identification, adequate precision or fairness.

The classical DIF screen accepts an explicit group file and shows the same
kind of assignment table before its Run button. It preserves text IDs including
leading zeros, rejects duplicate/blank Person IDs, and excludes blank groups
without inventing an additional group. It reports unmatched persons and extra
IDs; spelling and whitespace in IDs must match the response data exactly.
The coverage CSV lists observed combinations as FacetLevel, Group and Persons,
with Facet and GroupVariable identifying the selected columns.
Zero-count combinations are omitted from that CSV. Dense matrices above 50,000
cells are not allocated; the observed-count CSV remains available.

These checks do not run MFRM or DIF estimation. The classical screen uses
person-by-item mean scores and observed-score matching. Averaging ratings does
**not** adjust rater severity. Different rater assignments between groups can
confound comparisons, even if common references link an earlier calibration.

<a id="external-group-dff-plan"></a>
### External-group MFRM DFF: specification draft, not implemented

Keep the following decisions separate from the existing conditional heatmap
and classical DIF screen. No new statistical method or group response curve is
released by this draft.

1. **Define the group and the comparison.** School of the student, school of the
   rater and a label shown to a rater are different variables. The current group
   file represents a fixed person attribute, not a label that changes across
   ratings of the same performance. Specify the rater/group contrasts before
   examining estimates. Common reference recordings can remain calibration
   observations while being outside the student-group contrast; their modeling
   assumptions must be stated explicitly.
2. **Start with a uniform interaction.** A candidate adjacent-category model is
   `log[P(Y_pri=k)/P(Y_pri=k-1)] = theta_p - rho_r - delta_i - tau_ik + gamma_r,g(p)`.
   `gamma` follows the score direction: positive means higher scores than the
   additive prediction. RSM shares `tau_k`; a PCM decomposition needs explicit
   item-location/step constraints. This is a proposed model, not the current
   fixed-calibration cell calculation. Criterion-by-group DIF would use a
   different interaction term and a separately defined estimand.
3. **Identify the scale and interaction together.** For the reference-relative
   model checked below, reference recordings define the rater baseline and have
   no interaction term. Set the sum of `gamma_r,g` over raters to zero separately
   for each student group; group-average performance remains in the free person
   locations. Do not also impose `gamma_r,A + gamma_r,B = 0` without substantive
   justification. With the reference baseline, this extra symmetry assumption
   restricts the model; it is not simply a choice of origin. Check the rank and
   intended contrasts of the actual constrained design. Do not add School as an
   unrestricted main effect nested in Person, or force school ability means to
   be equal merely to identify a desired contrast.
4. **Specify JML and MML targets separately.** JML estimates person effects;
   a joint MML implementation must integrate the interaction model over the
   person distribution and state how group means, variances and references are
   modeled. Fixing MML EAPs is a conditional diagnostic, not that joint marginal
   likelihood. Fixed anchors require common scale constraints; their external
   uncertainty is a separate issue. TAM's group and design-matrix interfaces
   are reference capabilities, not proof of numerical parity with this app.
5. **Specify contrast uncertainty before testing.** For a planned contrast
   `L beta`, the covariance target is `L Var(beta) L'`, including shared fitted
   parameters. Combining individual SEs without covariance is insufficient.
   If a bootstrap is chosen, specify the generating population, rating design,
   repeated-person structure, external-anchor treatment and full-refit rule
   first. Bootstrap is not selected automatically by this draft. Predefine the
   multiplicity family; corrections do not repair an unidentified comparison.
6. **Separate curves from qualification.** Future group curves must use the
   same ability scale and stated rater/criterion conditions, show observed
   support and distinguish estimates from extrapolation. A nonuniform effect
   adds ability-dependent differences and needs a separately identified model;
   it is not merely another color in the present heatmap.

<a id="external-group-dff-identification"></a>
#### Design checks: three cases completed, no estimator released

The 2026-10-16 check covers **3/3 planned designs**, each under two explicitly
different model specifications. It uses two students per school, two common
reference recordings, four raters, three criteria and a five-category RSM.
Person locations are free; rater, criterion and four step locations each sum
to zero. No equality of school ability means is imposed. References are outside
the A/B groups, with no group-interaction contribution.

The reference-relative model has 20 free parameters: six person locations,
three rater coordinates, two criterion coordinates, three step coordinates
and six group-interaction coordinates. The additional A/B symmetry restriction
reduces the interaction coordinates from six to three, giving 17 parameters.
The table gives **rank / number of free parameters**, including RSM steps.

| Student assignment, with shared references in every case | Reference-relative effects | Additional A/B symmetry |
|---|---:|---:|
| Every student rated by all four raters | 20 / 20 | 17 / 17 |
| Students rated only by their two local raters | 16 / 20 | 16 / 17 |
| Crossed design except no B student rated by A_R1 | 19 / 20 | 17 / 17 |

*Note.* Full rank concerns structural identification under the stated model.
It does not establish finite estimates, adequate precision, model fit or
fairness. The symmetry column is an assumption-sensitivity comparison, not a
recommended way to repair missing observations. These are this project's
calculations, not results supplied by FACETS or TAM.

For each design, floating point SVD of the app's free-coordinate expansion was
checked against exact rational elimination of the full design plus explicit
constraint equations. Each contrast was checked both against the numerical
null space and for exact membership in the constrained design's row space.
The [reproducible check](../validation/external_group_dff_identification_check.py)
and [retained results](../validation/generated/external_group_dff_identification_20261016/)
record ten predeclared contrasts per specification:

- Four individual rater differences `gamma_r,A - gamma_r,B`.
- Four differences of those differences, comparing A_R1/A_R2, B_R1/B_R2,
  A_R1/B_R1 and A_R2/B_R1.
- Two within-school departures, A_R1 versus A_R2 in A and B_R1 versus B_R2 in B.
  These last two are not between-school DFF contrasts.

**What the results allow.** All ten contrasts are structurally estimable in
the crossed design. Under the reference-relative model, the nested design
identifies only the two within-school departures. In the missing-cell design,
the two differences of differences involving only fully observed raters
(B_R1/B_R2 and A_R2/B_R1) and the two within-school departures remain estimable.
The four individual rater A-minus-B effects do not: even the observed rater
effects share an unidentified group location when one cell is absent. Thus a
missing cell neither invalidates every contrast nor licenses filling that cell.

**Why shared references are insufficient for DFF.** In the nested design,
decrease every A student's location by `t` and add `(t, t, -t, -t)` to the four
A-group interaction effects. Leave all B effects, reference locations, rater
severities, criteria and steps unchanged. Every observed adjacent logit and
category probability stays the same, while the individual A-minus-B effects
change. Four independent such directions are verified in the tests. Even fixing
the reference person values or rater severities would leave these directions.
The actual 50-person lesson has the same four-dimensional deficiency when this
interaction is added; its additive model remains structurally identified.

The initial double-centering proposal hides some of this ambiguity by requiring
opposite A/B departures. It still leaves one unidentified direction in the
nested design, but makes all ten selected contrasts estimable in the missing-cell
case. The data did not recover the missing cell: the added assumption did work.
Adding one cross-rated student removes that last direction only in the symmetric
model; the reference-relative model still has three. This is not a sample-size
or precision recommendation.

**Next implementation gate.** State the intended contrasts and reference
assumptions first. Plan overlapping ratings of students from both schools across
the relevant performance range, then check identification and precision for
that actual design. Keep partial identification separate from a fully identified
fit. Only then specify joint covariance, uncertainty and supported group curves.
The checks above concern JML/RSM structure, not PCM/GPCM or a marginal likelihood.
MML needs an explicit group-distribution model and its own identification check;
fixing equal group distributions would introduce assumptions, not new overlap.
Error rates, power, coverage and cross-engine acceptance still require separate
protocols; formal G3 and subsequent MML studies remain paused.

Primary references: [FACETS nested-design guidance](https://www.winsteps.com/facetman64/nested.htm)
illustrates how imposed group constraints change the interpretation of effects;
[FACETS bias contrasts](https://www.winsteps.com/facetman64/table14.htm) distinguishes
overall measures and local departures; [TAM MML documentation](https://alexanderrobitzsch.github.io/TAM/reference/tam.mml.html)
describes structural design, groups and fixed parameters. The proposed acceptance
criteria above are this project's decisions, not endorsements by those sources.

<a id="external-group-dff-output-contract"></a>
#### Output specification: name the comparison before assigning a status

This specification has **3/3 deterministic acceptance cases**. It defines the
future external-group output; it does not release a new estimator, confidence
interval, test or response curve. Current application changes are explanatory
Help only. The validation script emits design-review JSON/CSV fixtures, not
analysis results for uploaded data.

Use two separate facts: whether the model is structurally identified, and
whether **each requested comparison** is identified. A partially identified
model does not become fully identified when the user selects only estimable
contrasts. A count table, successful optimization, or an arbitrary generalized
inverse cannot substitute for checking the actual constrained design and target.
The row-space/null-space criterion is also described by the
[author's documentation for R's estimability package](https://search.r-project.org/CRAN/refmans/estimability/html/estimability-package.html).
That reference supports the linear-algebra definition, not this app's nonlinear
estimator, its intervals or numerical parity with R/TAM.

For these fixtures, the target family is the **eight between-school contrasts**
listed above. The two within-school comparisons remain outside that denominator.
Define comparisons before inspecting outcomes; retain every requested target,
including those that cannot be estimated.

| State | First user-facing message | Design-review content | Next action |
|---|---|---|---|
| Model identified | The model effects can be separated under the stated assumptions. Estimation and precision still need checking. | Named targets, observed support and assumptions. Full rank is not a reporting clearance. | Use a validated joint estimator; check finite estimates, convergence, fit and precision. |
| Partially identified, with an estimable requested contrast | Only the listed comparisons can be separated. | A status for every requested contrast; show which ones remain unidentified. Keep the model limitation even if all selected contrasts are estimable. | Narrow the scientific question transparently, or improve the assignment. Validate an identified parameterization before estimating supported contrasts. |
| Requested targets unidentified | This assignment cannot answer the requested comparison. | Requested targets and reasons; descriptive counts remain available. | Review cross-school student ratings and reference assumptions. Do not translate this state into “no DFF.” |

Missing, stale, inconsistent or out-of-scope evidence is **not assessed**, not
one of these three conclusions. Before any future automatic classification,
bind the audit to the response design, group mapping, weights/row exclusions,
model, estimator, constraints and anchors. Revalidate after any change. The
bounded fixture helper rejects incomplete evidence or an empty target list;
it is not a production auditor. The current count-table interface continues
to describe observed assignments without assigning these model states.

**Concrete acceptance cases.** All cases use the reference-relative JML/RSM
model, with independent rank and contrast checks from the preceding section.

| Assignment | Model state | Estimable requested school contrasts |
|---|---|---:|
| Crossed students, with common references | Model identified | 8 / 8 |
| Missing A_R1 ratings of school B students | Partially identified | 2 / 8 |
| Local raters only, with common references | Requested targets unidentified | 0 / 8 |

The two estimable between-school contrasts in the second row compare school
differences for B_R1/B_R2 and A_R2/B_R1. If those alone are selected, the count is
2/2 but the model is still partially identified. If only the four individual
rater A/B effects are requested, the same design cannot answer any of them
(0/4). These counts must not be increased by adding estimable within-school
comparisons. Changing to the symmetric interaction model changes the assumptions
and the audit; it must not silently replace the intended model.

**Tables and downloads.** A future public comparison table should use a short
title, such as *Planned School Comparisons*, and a separate Note identifying
the target, model, reference baseline, sign, conditions and uncertainty method.
The initial view needs only Comparison, Design status and Observed support.
Keep rank, null-space coordinates and development identifiers in technical
details. Keep unsupported requested comparisons visible with a plain-language
reason, not a zero estimate or a reassuring green badge. In a future numerical
CSV, their estimate/SE/interval/p-value fields must be missing with an explicit
reason, never zero. The current fixture CSV contains no numerical effect fields.

**Intervals and plots have separate requirements.**

- A structurally estimable contrast does not establish finite estimates or a
  trustworthy SE. For a qualified joint fit, target `L Var(beta) L'`, with the
  covariance of jointly fitted parameters. Do not combine individual SEs as
  though the estimates were independent. In a deficient model, validate the
  reduced parameterization and invariance of supported contrasts before using
  any generalized-inverse covariance; do not give unidentified effects intervals.
- Before a future group response curve is drawn, check that its predictions
  are identified under all remaining parameter ambiguities. An estimable
  difference between raters does not by itself identify each rater's curve.
  State the common ability scale, criterion/rater conditions and reference
  assumptions; distinguish observed support from assumption-based extrapolation.
- No inference is released by these fixtures: point estimates, intervals,
  p-values and group curves remain withheld in **all three states**. Existing
  conditional bias plots and classical DIF outputs retain their own scope and
  cannot be substituted for external-group MFRM results. An unidentified target
  also remains in the recorded planned comparison family; omitting it must not
  silently redefine multiplicity or the scientific question.

The bilingual Help section **Planning school-group comparisons** is reachable
through the global Help button: select the differential-interaction topic and
expand further details. The heatmap popover describes that route. It explains
these distinctions without internal status codes. Human comprehension and
browser acceptance remain separate from automated navigation checks.
Formal G3 and later MML studies remain paused. The next implementation decision
is the joint-estimation and uncertainty specification; no bootstrap method or
default interval procedure is selected by the design checks alone.

## 10) Minimal checklist before publishing results
- Anchor files archived with the analysis run.
- Anchor summary non-zero where expected.
- Report shows intended A/G/X flags.
- Key facet ranges and ordering are stable across reruns.
- Shared content, connectedness, drift, and uncertainty support the intended comparison; counts and A/G/X flags alone do not establish valid linking.

<a id="two-school-recordings"></a>

## 11) Worked lesson: two schools and common recordings

Select **Built-in example → Two schools, common recordings (360 ratings)**.
Open **Lesson: compare two schools using common recordings** above the input
preview. The app offers the lesson in English and Japanese. Download its CSV
from **Example details, references & CSV** in the sidebar.

| Performances | Who rates them? | Ratings on three criteria |
|---|---|---:|
| A01–A20: school A students | A_R1 and A_R2 | 120 |
| B01–B20: school B students | B_R1 and B_R2 | 120 |
| C01–C10: common references | All four raters | 120 |

Each person supplies one recording on the same speaking task. One row is one
rater's score on one criterion, from 0 to 4. The generator uses RSM, person
locations drawn from N(0, 1) for both schools and references, rater severities
−0.6, −0.2, 0.2, 0.6, criterion difficulties −0.4, 0, 0.4 and thresholds
−1.5, −0.5, 0.5, 1.5. These values are illustrative, not published estimates.
Equal generating distributions do not force equal realized school means.

1. **Prepare the overlap.** Use the same C ID for the identical recording at
   both schools. Preserve individual ratings and use unique rater IDs. Different
   performances must not share an ID. Reference performances are excluded from
   school student summaries; in real data, distinguish reference membership from
   school membership in a separate roster.
2. **Fit jointly.** Guided defaults use RSM/JMLE with Person, Rater and Criterion.
   The singleton Task column is excluded automatically. Leave fixed/group anchors
   empty. The model constrains facet origins and estimates persons on that common
   scale; it does not force the two school means to be equal. School is grouping
   metadata, not an additional freely estimated facet nested within Person/Rater.
3. **Inspect the comparison.** Check interpretation holds, Person–Rater
   connectivity, fit, category use, rubric compatibility and reference coverage.
   Removing or giving school-specific IDs to all common recordings disconnects
   this design despite identical criterion labels. Connectivity is necessary here
   but does not establish precision or invariance. Ten recordings are an example,
   not a validated minimum. No school-specific reference drift is planted or
   empirically assessed by this teaching example.
4. **Read and report.** Higher Person logits indicate higher modeled performance;
   higher Rater logits indicate greater severity. Present qualified estimates and
   fit, a Wright map when useful, the rating assignment and common-reference
   selection. A/B student means can describe the samples after excluding C.
   The app does not supply a validated school-mean contrast interval or a causal
   school-effect estimate. Individual SEs alone cannot supply the contrast SE.

### Later administrations: fixed anchors require another review

Reuse vetted reference measures only after checking labels, common content,
scale origin, task/criterion compatibility and reference stability. In JMLE,
Person anchors refer to the same modeled performances; an MML EAP score is not
an independently fixed person parameter. Do not reuse all estimated anchors
indiscriminately or fix A and B means to the same value to create comparability.

The app's **Hard-anchor consistency** check compares constrained estimates with
their imposed values. Zero deviation is expected under a hard constraint and is
not an empirical drift test. Independent calibration comparisons, with scale
alignment, uncertainty and content review, remain separate work. Current SEs
do not propagate uncertainty in externally estimated anchor values.

### Sources and scope

Reviewed 2026-10-16: FACETS **64-bit** official manual,
[Equating tests](https://www.winsteps.com/facetman64/equating.htm),
[Labels and anchors](https://www.winsteps.com/facetman64/labels.htm), and
[Nested facet designs](https://www.winsteps.com/facetman64/nested.htm).
The equating page recommends examining common elements across independent
calibrations before anchoring or combining datasets. Its common-element count
is explicitly a rule of thumb, not an acceptance test. This synthetic lesson
demonstrates a connected joint-calibration workflow; it does not replace that
empirical review or establish numerical agreement with FACETS/TAM.
