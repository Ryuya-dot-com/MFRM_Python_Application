# Changelog

All notable changes to this standalone Streamlit distribution should be recorded here.

## Unreleased — Windows desktop correction, 2026-10-29

- Preserve observation indices when assembling score probabilities. This fixes
  a crash after missing-score rows are excluded and prevents probabilities
  from being attached to different observations after row relabeling.
  Add regressions for permuted, gapped, duplicate and text indices, plus
  desktop save/restore with an excluded score row.
- Read FACETS `PtMea` correctly and keep `PtBis`/`PtBisI` in separate fields.
  Document the application's person-measure correlation definition and the
  limited FACETS comparison in `validation/FACETS_WINDOWS_DESKTOP_REVIEW_20261029.md`.
- Handle Windows pipe EOF during worker polling as well as receiving. Retain
  the completion/failure message until the child exits, then check the exit
  code and saved result. This fixes `WinError 109` after successful analysis
  and restores cancellation, failure reporting and restart.
- Add real spawned-process regressions for cancellation, invalid input,
  successful retry, project restoration and a child that exits without a
  result. Add a Windows/Python 3.12 desktop CI job.
- Replace deprecated pandas in-place renaming and implicit numeric-to-string
  fills in publication tables. Preserve blank held values, original numeric
  calculations and interpretation warnings; do not suppress warnings globally.
- Verify 125 related tests with warnings treated as errors and all 56 engine
  self-tests on Windows 11/Python 3.12. Build and run a local Windows executable,
  including worker/save/restore checks. See the bounded GUI and packaging
  acceptance record in `docs/desktop_prototype.md`; this is not a public
  desktop binary release or clean-machine certification.

## 0.2.19-beta — 2026-10-28

This version consolidates the changes recorded below since 0.2.18-beta.
Required Linux/Python 3.11 and 3.12 CI passed for source candidate
`f5139bb7ec02fbfa2f1079ef6e381acff468d683`
([run](https://github.com/Ryuya-dot-com/MFRM_Python_Application/actions/runs/40687727467)).
The final main CI and hosted acceptance are recorded with the GitHub prerelease;
source publication alone does not establish hosted acceptance. Later documentation-only
preparation does not constitute a new CI run. The Tkinter sources remain a
prototype; this Web update does not release desktop binaries.

- Provide selected fixed-SD RSM/PCM facet levels and prespecified differences
  with pointwise asymptotic normal intervals when termination, identification,
  boundary, full observed-information and quadrature-sensitivity checks pass.
  Support the default centering and grouped 1,200/3,000-rating inputs. Preserve
  full covariance for differences and withhold failed intervals without refitting.
  Share the interval record across bilingual tables, figures, CSV/JSON/HTML and
  publication tables. Reject stale or inconsistent saved interval records.
  Finite-sample coverage, Person-score inference, free-SD inference and automatic
  equivalence conclusions remain unqualified; older study decisions are unchanged.
- Keep subsequently computed or loaded intervals in the opted-in history snapshot
  for that exact run, including a restored run. Preserve metadata-only histories.
- Make the selected-interval report table identify its facet and logit units,
  refer to its actual qualification label, and keep that label legible in Word/PDF.
- Align the candidate version and software citation. Leave its release date
  unset while unpublished; show an undated citation instead of inventing a date.

### Directory review corrections and earlier development changes

- Carry person-level JMLE holds into shared Web/desktop reporting and reopened
  projects. Blank endpoint-only person estimates, fit and uncertainty in report
  tables; withhold person reliability without silently dropping people or
  suppressing valid rater/task results. Keep raw calculations for audit and
  label retained extreme-person figure locations as terminal calculation values.
- Preserve the opt-in history snapshot preference across navigation and input
  changes. Keep history reachable when current inputs make the last result stale,
  so users can restore its saved result without another fit.
- Use `MFRM` as the browser title and retain the current interpretation-hold
  message in saved basic HTML reports. Verify the source candidate in a separate
  Python 3.12 environment, keeping local research artifacts outside the bundle.
- Recheck retained R comparisons with a unit-scale ULP allowance only for
  recomputed metrics, propagating operand-scale rounding for truth-point finite
  values reevaluated across numerical environments. Keep source/returned values and decisions exact; reject
  changed raw R values even when their difference is one ULP. Statistical
  thresholds and inference qualification remain unchanged.

- Repair release self-tests after the shared-engine extraction and UI
  simplification: inject the forced MML fallback into the actual engine module,
  accept concise nonempty chart headings, and retain accessible labels for the
  Basic/Precision result navigation. Estimation and inference holds are unchanged.

- Keep interval levels in their original fractional units in table CSVs while
  displaying percentages on screen. Verify native table saving by reading back
  saved files, distinguishing cancelled headless file-picker requests from an
  application export failure.
- Lead guided calculation details with stopping/hold status and a return to
  existing data/settings. Fold technical records and remove only duplicate
  alerts covered by the open panel; preserve other warnings and saved results.

- Correct the beginner result view after its disclosure review: distinguish
  withheld and all-minimum/all-maximum person rows, preserve readable facet
  estimates, and offer Basic/Precision tables with only reader-facing columns.
  Keep technical source values in calculation records; do not fabricate missing
  intervals. Show fit rating counts with their weight-sum meaning and precision.
  Route unfinished runs directly to opened calculation details instead of document
  generation. Preserve existing inference holds and saved numerical results.

- Simplify the beginner Streamlit workflow: fold analysis settings, show one
  first-read summary and evidence disclosure, and keep result pages focused on
  their selected content. Show measures and fit before technical columns;
  retain the complete tables, qualifications and exports. Move convergence
  details and explanatory notes behind disclosures while keeping inference
  and extreme-person holds visible. Repair sticky navigation after scrolling.

- Prevent history and saved-result comparisons from taking IC differences
  across unknown/different conditions or unfinished/nonfinite results. Match
  normalized inputs, constraints, population settings and quadrature; keep
  fixed versus estimated SD distinct. Omit pooled criterion summaries when
  records do not qualify for a same-condition numerical difference.
  Explain omissions in Japanese/English and CSV while retaining inference holds.
  Use the same qualified BIC in compact and detailed history, preserve legacy
  saved values, and retain person-based BIC metadata in new compact entries.
  Download history records without rerunning the page.

- Offer an explicit finer-quadrature refit after the optional fixed-SD MML
  sensitivity check. Reuse the constrained direct optimizer, reevaluate its
  coordinates at a still finer rule, and retain every outcome in a separate
  CSV. Exclude unfinished or invalid candidates from pair differences; preserve
  the original analysis, normal comparison export, defaults and inference holds.
  Verify three retained synthetic candidates against independent references and
  saved ConQuest RSM/PCM results; check bilingual downloads and navigation.
  Identify the retained refit table after navigation and explain incomplete
  outcomes only for affected candidates. Keep translated statuses/reasons in
  CSV and name candidates near bounds; do not invent missing reasons in older
  records or trigger refits while displaying them. Numerical behavior is unchanged.

- Preserve EM's actual stopping tolerance and tested relative log-likelihood
  change, together with effective inner optimizer controls, in convergence
  records and exports. Do not interpret the duplicated final trace value as a
  zero stopping change or reconstruct missing controls for old results.
  Clarify that a tighter stopping tolerance does not change quadrature accuracy;
  estimation updates and defaults are unchanged.

- Add an on-demand numerical check under model-comparison calculation details.
  For small fixed-SD MML comparisons, retain each candidate's fitted coordinates
  and evaluate finer integration rules without refitting. Export likelihood/IC
  sensitivity, gradients and provenance; preserve the original comparison and
  inference holds. Cache the check across navigation and language changes, and
  download its bilingual record without rerunning the page.

- Validate AIC/BIC on a retained 32-person input using independent likelihood
  integration, R's information criteria and matched ConQuest RSM/PCM results.
  Reject model comparisons with known person/facet confounding, including old
  cached comparisons. Keep numeric-integration and model-selection qualifications
  separate; estimation defaults and inference holds are unchanged.

- Save model-comparison tables with Japanese or English headings, calculation
  status, selected facets, source analysis, estimation method and interpretation
  limits. Preserve full numeric precision in CSV; show compact numbers onscreen.
  Keep the technical record under calculation details, and download either file
  without rerunning the page. Verify native downloads and retained selections in
  an isolated local test browser; public deployment remains pending.

- Select step and slope facets explicitly for RSM/PCM/GPCM comparisons.
  Reuse a current fit only when its structural specification matches; keep
  separate cached results for each choice without fitting on selection changes.
  Record the source analysis and chosen facets in CSV, exclude unfinished fits
  from differences/rankings, and keep model-selection inference withheld.
  Use a compact bilingual table with qualified exports instead of a grid's
  context-free CSV; place calculation details and interpretation limits in
  collapsed sections. Estimation algorithms and defaults are unchanged.

- Preserve non-default prediction/simulation controls and uploaded inputs across
  panel navigation within the same analysis. Show retained filenames and an
  explicit removal control; distinguish widget cleanup from deleting a file.
  Clearing a DIF group file also clears its old screening result.
- Connect the app-engine runner's explicit bias-pair argument to execution,
  including an empty list meaning no pairs. Bind original input mappings and
  iteration settings in script generators to the requested analysis; foreign
  or unowned session snapshots no longer override it. Numerical algorithms and
  public/legacy route boundaries are unchanged.

- Bind optional-analysis caches and prediction/refit exports to the complete
  existing analysis identity. Exclude foreign or unowned cached results after
  history restore and export the settings attached to the actual computation.
- Retain classical DIF inputs, output, audit and caveats across language changes
  and panel navigation; invalidate screens when response/group/settings inputs
  change. CSV downloads do not trigger recalculation.
- Retain bias screening thresholds per facet pair and use them consistently in
  local screens, first-read summaries, CSV, PDF/Word/HTML and the app-engine
  reproduction runner. Record the selected criteria without changing estimates
  or releasing inferential holds; invalidate earlier document caches.

- Make first-read summaries follow each facet's reliability qualification, retain
  unavailable PCA and missing stability evidence, and preserve bias inference
  holds and unassessed pairs even when screening flags are zero. Localize the
  explanations in English/Japanese and carry them into basic/detailed PDF,
  Word and HTML reports. Preserve valid neighboring facets and saved numbers;
  invalidate cached reports with the earlier summaries. Numerical engines and
  inference qualifications are unchanged.

- Carry recorded optimizer status and first-read cautions into basic and detailed
  PDF, Word and HTML reports, including residual-PCA stability, interaction review
  flags and skipped checks. Reuse the screen's diagnostic rules; preserve inference
  holds and existing table values. Invalidate older document caches and keep PDF
  section headings with their following content.
- Qualify newly fitted summary tables only after attaching the optimizer record
  and identification audit. Successful termination no longer becomes an artificial
  "termination unverified" hold that blocks the sample guide. Guide progress uses
  explicit Boolean optimizer evidence. Failed attempts retain their holds and can
  open existing diagnostics/settings for an explicit retry; bounded exception
  descriptions survive language changes. Numerical routines, stopping criteria
  and method-specific inference restrictions are unchanged.
- Clarify optional-tool setup in Japanese: translate recovery inputs and
  network entry controls, explain replicate counts/seed/model/estimator, and
  preview the planned synthetic rating count and score range. Align setup
  controls when labels wrap. Keep numerical definitions, cached settings and
  inference holds unchanged; the default ten replicates is not a precision
  guarantee.
- Separate optional model comparison, G/D planning, network tools, recovery
  simulation and fixed-SD numerical-check execution from current-result Tables.
  Open one purpose-led tool at a time from Guided supporting material or the
  advanced content selector. Keep existing check results and holds in Tables;
  retain simulation settings and cached output when leaving the panel. No
  estimator, qualification rule or technical export schema changes.
- Group Guided report checks by the question they address and localize the
  questions, current statuses and review guidance in English/Japanese. Keep
  every distinct status, original finding and action; source details remain
  available separately, including English technical text. Keep calibration
  holds visible and optional analyses optional. Estimation, qualification
  rules and exported evidence schemas are unchanged.
- Align input-readiness messages with Guided and Advanced run controls. Keep
  sample counts, missingness and score-pattern screens advisory; block known
  input errors, including unusable retained rows and constant retained scores.
  Preserve estimator validation and post-fit holds. Remove count-based precision
  promises, clarify bilingual actions and display the next action as wrapping text.
- Replace nested advanced-report tabs with one content selector. Generate only
  the selected report and preserve its selection across language changes and
  leaving/returning to the view. Separate the first reporting journey from
  optional technical checks in the README.
- Restrict G/D ANOVA to complete balanced unit-weight designs with one or
  two non-object random facets, at least two levels per factor and one rating
  per cell. Reject missing/repeated cells and invalid scores explicitly;
  retain signed variance estimates beside zero-clipped calculation values and
  always include the observed design in D-study output. Withhold unqualified
  G/Phi intervals in the bilingual panel and CSV; record technical bootstrap
  attempts, per-coefficient failures and current inference status. Export
  calculation conditions, source identity and reasons with the results.
  Save the dedicated G/D CSV without rerunning the report page.
- Apply one current optimizer-termination rule to live interval plots, technical
  CSVs, reliability, scoring prose, manuscript guidance, and saved-history/comparison
  displays. Failed, absent or non-Boolean success cannot inherit stale Ready labels.
  Retain technical values and original snapshots; replace affirmative stale
  cautions with a stable hold reason and invalidate document caches.
- Complete a separate bounded writing-robustness study: 96 fixed-SD RSM fits and
  eight independent integral checks passed their numerical criteria. Ability-related
  assignment under the constructed weak link produced a mean rater-block contrast
  of 0.843 despite a generating contrast of zero. Report paired Monte Carlo
  uncertainty, full covariance and generating-contrast containment; no inference
  hold is released. Add boundary/rounding analysis of the prior retained cohort,
  reproducible tables and figures, and practical reporting guidance.
- Add result-aware applied-report tables and manuscript drafts: observed design,
  workload and shared-person rater overlap, category counts including unused
  categories, separate likelihood weights, and recorded estimation settings.
  Recompute retained counts instead of copying stale display summaries; distinguish
  optimizer termination from numerical qualification. Apply current holds to
  direct/reopened manuscript exports and publication SE/CI, reliability and figure
  notes when optimizer termination is failed or unrecorded. Preserve the concise UI.
- Add the separately authorized medium writing-assessment study and its retained
  protocol: 200/500 persons, matched rating workloads, full-covariance contrasts,
  finite-GH stationarity, matched-accuracy performance, and independent references.
  Complete 192 datasets / 384 finite-GH fits, eight independent integral checks,
  and four TAM-point comparisons. Retain timeout censoring in the matched benchmark
  and Monte Carlo uncertainty in exploratory coverage tables and figures.
  Preserve the old frozen cohorts and all inference holds; results concern the
  prespecified fixed-SD RSM conditions, not all application estimators.
- Rewrite the main UI around grouped controls and progressive disclosure.
  Remove repeated headings, secondary labels and instructions; place start
  choices side by side, combine input previews with their counts, lead results
  with the measures table, and group document options with the download action.
  Consolidate version, scope, license and citation in About MFRM. Preserve
  visible input labels, accessible navigation names, privacy notices and
  inference holds; keep bilingual selection and session restoration intact.
- Complete the user-authorized bounded pilot: 48 synthetic MFRM datasets,
  96 frozen-source Python fits, four predetermined TAM comparisons and separate
  normal-model shrinkage controls. Retain all attempts and prespecified criteria;
  distinguish native stopping success from stationarity eligibility, covariance
  effects and point-estimation error from interval coverage. No MFRM coverage
  claim or resumption of paused formal G3/MML cohorts.
  Fix live forest plots bypassing inference holds, reliability calculations
  mixing unmatched estimate/SE inputs, and saved reliability qualification in
  publications and exports. Clarify discrete residual Q–Q interpretation and
  the fixed-zero shrinkage prior; refresh document caches and bilingual text.
  Related tests: 265 passed; built-in self-tests: 56 passed.
- Carry current uncertainty holds through fresh and restored step/measure tables,
  FACETS-style tables, assumption audits, scoring-quality decisions and checklists.
  Remove duplicate unqualified APA conclusions and EB confidence bars; retain
  technical point values with explicit plug-in uncertainty metadata, including
  the zero-variance boundary. Preserve the EAP reliability formula when its
  inference status is withheld. Refresh cached documents and bilingual guidance.
  No new estimator, interval coverage qualification or G3 study is introduced.
- Close the advanced report's facet-equivalence bypass: remove unqualified
  TOST/BIC/ROPE calculations, interval plots and combined conclusions that used
  SEs without joint covariance or current inference holds. Show bilingual hold
  guidance for every estimator and reopened results, and export a stable
  qualification reason; invalidate document caches from before the hold.
  Preserve the old renderer as local audit evidence.
  Document normal-reference controls for equivalence and multiplicity; no new
  inferential method or estimator qualification is released.
- Check equivalent PCM parameterizations at three retained points, including
  reference constraints, sign reversal and latent-origin translation. Verify
  native likelihood/derivative agreement and full-covariance contrast
  invariance; distinguish covariance omission, changed estimands and expected
  singularity after constrained expansion. No production or inference changes.
- Audit fixed-SD PCM information and covariance at two retained original fits
  and six saved refit endpoints. Separate derivative accuracy, quadrature
  error and stationarity; show differing effects on level and contrast SEs.
  Retain positive-definiteness, reference refinements and all starting points.
  No production SE, default quadrature or inferential qualification changes.
- Extend the uncertainty review with analytic controls distinguishing EAP
  sampling variability, posterior SD and interval coverage; shared calibration
  uncertainty and contrasts; retained TAM SE definitions; and eigenvalue-floor
  artifacts. Clarify 95% normal-approximation versus 90% grid-quantile intervals
  and the fixed/free-SD inference hold in reporting guidance. No production
  uncertainty algorithm or inferential qualification changes.
- Assess a development-only continuous posterior CDF/quantile prototype at
  fixed calibration. Retain 44 independent integration comparisons, native
  app replays, known-normal and affine-transform controls, and explicit input
  and integration-budget failures. Document interpolation artifacts, posterior
  collapse on coarse GH rules, decision losses and requirements for consistent
  downstream scoring. Production algorithms and inference holds are unchanged.
- Extend the numerical review to decision boundaries. Retain fixed-calibration
  controls where accurate EAP/SD coexist with threshold reversals in grid
  quantiles, native TAM posterior replays, and PCM integration/recalibration
  decompositions. Distinguish interval definitions, display proximity and
  parameter boundaries. No estimator, scoring default or inference release changes.
- Extend the TAM/Python review using installed TAM 4.3.25 function bodies and
  bounded numerical controls. Separate EM updates, integration rules, free-SD
  objectives, score/fit definitions and decision uncertainty. Retain a stopped
  TAM fit's returned-point discrepancy and a long-response posterior underflow
  control alongside their successful controls. Add interpretation guidance;
  no estimator, numerical default or inference qualification changes.
- Correct two findings from the October 23 TAM review: reject unsuccessful,
  nonfinite, infeasible or worsening EM M-steps and retain their failure reason;
  withhold fixed-SD MML inference as well as free-SD inference in current and
  reopened results, reports and assignment-sensitivity preflight. Retain valid
  estimates and optimizer termination separately from inferential qualification.
  Record the independent continuous-integral check of 24 fixed-slope/fixed-SD
  GPCM endpoints and distinguish TAM coordinate-curvature SE from joint SE.
  Document backend tradeoffs; no R runtime, new estimator or inference release.
- Define three external-group DFF output states using the requested comparisons,
  with executable design-only fixtures: 8/8, 2/8 and 0/8 between-school contrasts
  for the crossed, missing-cell and nested cases. Exclude within-school effects
  from that denominator, retain partial identification after target selection,
  and reject absent or inconsistent audit evidence. All numerical releases stay
  withheld. Add bilingual school-comparison guidance to the existing Help topic;
  no uploaded-data auditor, estimator or new contextual navigation edge is added.
  Preserve the heatmap metric through Help and resolve diagnostic menu labels
  at render time so widget value serialization uses the displayed language.
- Check three proposed external-group DFF designs with both SVD and independent
  rational constraint calculations. Distinguish reference-relative group effects
  from an additional A/B symmetry assumption: the latter can make missing-cell
  contrasts identifiable by assumption. Shared references link the additive
  school lesson but leave DFF ambiguities. Add bilingual next-step guidance;
  no external-group estimator or inferential qualification is released.
- Add distinct-person group assignment counts to the two-school lesson and
  classical DIF preflight, with ungrouped reference recordings and a lesson CSV.
  Reject duplicate/blank Person IDs, preserve text IDs on CSV/Excel upload and
  exclude blank groups. Avoid allocating dense coverage matrices over 50,000
  cells while retaining observed-count downloads. Correct the claim that mean
  scores remove rater severity. Document a separate external-group MFRM-DFF
  specification draft; no new interaction estimator or inference is released.
- Connect the two-school lesson to DFF interpretation: default its selected pair
  to Rater × Criterion without overriding user choices, and add a navigation
  button that does not refit. Put the heatmap and a four-column reading table
  before collapsed technical details. Add bilingual guidance and five questions
  on severity, signs, missing cells, inference and cross-school assignment.
  Preserve missing-only heatmap rows, distinguish genuine zero values, and
  localize chart labels. Statistical qualification remains unchanged.
- Reject Person fixed/group/dummy constraints in ordinary MML before optimization:
  they were audited but did not constrain the marginal likelihood or EAP scores.
  Preserve non-Person facet anchors, explain the engine-specific scope in both
  languages, and distinguish shared-person joint calibration from Person anchoring.
  Clarify DFF definitions, score-direction signs and the fixed-calibration/EAP
  boundary beside the existing heatmap; no new inferential claim is qualified.
- Add a two-school synthetic example: 40 students, ten shared reference
  recordings and 360 ratings, with a bilingual joint-calibration lesson and CSV.
  Explain common IDs, reference exclusions, model directions, reporting limits
  and later fixed anchoring. Separate anchor-count thresholds from validated
  links, and hard-constraint agreement from empirical reference stability.
  Unanchored joint calibration is no longer labeled "Not linked" by the anchor
  checklist; dummy constraints do not count as evidence of a reference link.
- Append an English Q&A before References in basic and detailed Word, PDF and
  HTML reports, and offer the same guidance in Japanese or English beside the
  download. Explain overfit, PCM versus RSM, facet directions and
  bias interactions with definitions, examples and next checks. Include the
  cited sources once, distinguish undated manuals from retrieval dates, retain
  interpretation holds and invalidate cached reports without the Q&A.
- Define Strata, underfit and overfit next to results and in report Notes, with
  bilingual formulas, examples and decision limits. Distinguish model expectation
  from review thresholds. Reuse the qualified reliability table in the reporting
  view, removing a duplicate SE-based calculation that contradicted EAP summaries.
  Expand the planned comprehension review to six tasks; human acceptance is pending.

- Preserve emphasis, code text, nested bullets and tables in methods/writing
  appendices across Word, PDF and HTML; escape literal markup and keep headings
  with their explanations. Replace contradictory population-scale wording and
  internal state labels in the writing guide. Propagate calibration-wide holds
  to dependent claim-guide sections so reliability or targeting is not marked
  ready while final interpretation is withheld.
- Explain what to add to the basic report and link directly to existing table
  downloads, interpretation checks and manuscript templates. Remove a shared
  instruction that pointed to a control absent from the detailed view, and
  distinguish templates from completed drafts. Specify formative review tasks
  and scoring under the current five-novice/two-experienced target; human
  acceptance remains untested.
- Start report export with a basic set of two tables and two figures where
  available; make person rows, additional diagnostics and methods/writing
  appendices optional. Add bilingual reading steps and examples for language,
  piano and skating assessments; preserve all interpretation holds.
- Use separate table numbers, brief italic titles and explanatory Notes in
  Word, PDF and HTML. Remove colored grids, retain horizontal rules, combine
  interval bounds, and move constant qualifications into Notes while keeping
  mixed uncertainty bases by row. Keep results before detailed appendices and
  include report scope in the download cache key.
- Guide CSV selection through 14 bilingual table names and reading notes;
  retain access to every technical table and the existing privacy controls.
  Include category diagnostics and step ordering in normal download bundles.
- Correct the Matplotlib category-curve fallback, which treated long-form
  column names as separate curves. Check plotted probabilities for RSM, PCM
  and GPCM against the shared curve data, including category labels and sums.
- Include Infit/Outfit alongside estimates and uncertainty in publication
  tables. Keep interpretation holds in figure captions, replace internal bias
  status codes with readable labels, and improve Word/PDF column widths and
  caption pagination. Add a Japanese writing-assessment guide for three tables
  and two figures; broader reporting and novice-user acceptance remain open.
- Refine the active roadmap around graduate students learning MFRM: audit eight
  reporting areas, complete a language-assessment journey, then develop piano
  and skating examples. Add checkbox milestones with separate denominators for
  output acceptance, learning scenarios, human understanding and performance.
  Review primary FACETS, research and competition sources; distinguish existing
  SE/posterior-SD methods from proposed bootstrap work and retain the pause on
  formal uncertainty studies. The roadmap changes do not qualify statistical
  inference.
- Consolidate the bilingual estimation tutorial into one current explanation.
  Rewrite the full Help guides around user actions, remove development-only
  state names from guidance, and keep detailed workflow checks in collapsed
  disclosures. Distinguish convergence, uncertainty methods and interpretation
  holds; remove automatic fit/category/reliability conclusions from reviewed
  guidance and chart explanations.
- Explain model-comparison and bias holds in ordinary language, while keeping
  technical reason codes in exports. Exercise all eight Help panels in both
  languages without estimation and verify that a cached model recommendation
  remains withheld after the display changes.
- Hold JMLE interpretation for extreme-score free facet levels and preserve the
  screening reason in summaries, guidance and exports.
- Carry required calibration holds through publication decisions, APA sentences,
  report summaries and measure uncertainty labels. Preserve rerun/hold states
  when summarizing them again; unknown statuses do not become ready.
- Treat zero as a valid reportable Person coordinate in assignment-sensitivity
  checks, and reject non-finite reportable coordinates.
- Guard dense rank-audit allocations before constructing matrices; reuse MML
  posterior summaries and omit exactly-zero quadrature weights in the joint kernel.
- Reuse readiness and report-summary tables within a single export call, with
  complete input fingerprints, copy isolation and scope expiry.
- Localize first-read guidance, progress and sample explanations; correct stale
  testlet navigation and unsupported claims in synthetic example descriptions.
- Put installation and first-use instructions first, and distinguish local-only
  retained evidence from repository documentation.
- Distinguish corrected contrast calculations from the still-pending statistical
  validation in the README and replace internal validation-stage labels in
  manuscript guidance with the checks readers need to understand.
- Restrict copy-hygiene checks to maintained surfaces, preserving historical
  research records; update older cache tests to supply the required report tables.
- Revise 382 Japanese interface strings across setup, data import, first-read
  guidance and reporting, replacing literal translations and implementation
  jargon with concrete instructions. Rewrite the Japanese estimation tutorial.
- Correct both languages' setting descriptions: group anchors fix group means,
  Dummy facets fix all levels at zero, Xtreme affects report measures, and
  Umean/Uscale apply a linear display transformation. Preserve interpretation
  holds and distinguish uncertainty intervals from universally valid confidence
  intervals. No numerical calculation changes are made in this copy revision.
- Translate the wide-table detection explanation and describe it as a heuristic;
  ask users to check the proposed score and ID columns.
- Separate current user/developer instructions from historical v0.2.0 model
  guides and v0.1.x publication actions. Add a documentation map and require
  workflow changes to update the instructions users actually encounter.
- Correct the README's fixed-only population-SD description and the bundled
  anchor guide's silent-ignore, dummy-facet and group-stability claims against
  the current implementation. Identify local-only evidence links and preserve
  dated research findings without treating them as current validation.
- Carry current JMLE and free-population-SD MML interpretation holds into each
  measure row in the UI, CSV bundles, APA tables and scaled FACETS-style tables,
  including saved results with older positive uncertainty flags. Preserve
  original numerical values and uncertainty methods in technical exports.
- Give Word, PDF and HTML measure tables an explicit uncertainty basis and
  nearby explanation distinguishing conditional SEs from EAP posterior SDs.
  Leave held SE/interval and reliability/separation precision cells blank in
  publication tables, with their scope and the technical CSV source stated.
- Wrap PDF table cells within the page width and make Word tables more compact,
  retaining captions with tables and repeating headers across pages.
- Refresh the active roadmap around existing end-to-end journeys, user
  understanding, target-specific mathematical evidence and measured wait time.
  Separate the current queue from retained research chronology, add a bounded
  UX1 acceptance protocol and link it from the release checklist. Keep formal
  G3/MML studies paused and distinguish local fixes from public acceptance.
- Replace EM's function-local result class with SciPy's serializable
  `OptimizeResult`, enabling the existing report and document caches. Preserve
  the optimization, numerical outputs and inference holds; strengthen the
  fixed-SD adapter's unchanged-input check to include the complete result.
- Translate stale/restored-result notices and rerun controls, explaining which
  analysis owns the displayed output without raw internal setting identifiers.
- Exercise four real input-to-export AppTest journeys, bilingual score-mapping
  recovery and a restored free-SD snapshot with older positive uncertainty
  flags. Record limited timing evidence and remaining UX/J0 tasks in the audit
  and active roadmap.
- Clarify the TAM reference role against its versioned manual and retained
  numerical evidence. Distinguish fixed-point PCM agreement from refitting,
  quadrature accuracy and joint uncertainty; retain the historical EAP failure
  and prioritize the remaining N1/N2 questions in the active roadmap.
- Separate stopping error from integration differences in the two retained
  TAM comparison inputs using the existing fixed-SD checker. Clarify in both
  bundled guides that small finite-rule gradients and stable EAPs do not
  establish posterior-SD accuracy or stability after recalibration.
- Distinguish cross-language agreement from quadrature accuracy in the current
  comparison record and bilingual guides. A fixed-point R control using the
  same nodes and masses reproduced both Python's finite-rule results and the
  remaining GH31 differences from the continuous reference.
- Specify the next TAM comparisons by model: current-code RSM checks first,
  then GPCM parameter/constraint mapping before any refit-equivalence claim.
  Separate historical native TAM evidence from independent R controls, and
  preserve the pause on formal uncertainty and coverage studies.
- Recheck the current RSM mapping at all six retained TAM endpoints across
  regular, sparse and unbalanced inputs. The affine design matches exactly;
  fresh probabilities, likelihoods, scores and structural derivatives agree
  within the declared tolerances. Preserve the historical cached-EAP failure
  and distinguish this check from refitting, integration accuracy or joint
  uncertainty. Move the next model comparison to GPCM constraint mapping.
- Map current GPCM probabilities to TAM at six retained points and three
  unit-slope PCM controls. Check likelihoods, posterior scores and derivatives
  including log-slopes; detect omitted intercept scaling and frozen-design
  errors. Document the remaining nonlinear shared-effect and scale constraints
  before claiming equivalent refits. Preserve the coarse-quadrature failures.
- Compare twelve fixed-slope/fixed-SD GPCM refits per engine across retained
  inputs, two starts and two matched grids. Verify TAM stopping conditions and
  frozen parameters throughout its iterations; cross-evaluate all returned
  points and retain cached scores separately. Correct a grid-width omission
  in an auxiliary likelihood diagnostic without changing the refits or their
  primary checks. Keep continuous accuracy and full-model inference separate.
- Reprioritize the roadmap around UX1 and J0 on the same analysis-to-export
  journey. Define bounded performance comparisons, distinguish existing reuse
  and grouped calculations from proposed changes, and prioritize history copy,
  browser operation and user understanding. Limit the next numerical check to
  retained GPCM endpoints; prioritize defects affecting current outputs and
  keep formal studies paused. Align active planning and release references.
- Reuse manuscript claim/evidence frames and shared-context agreement summaries
  within the existing bounded export scope. Preserve the first candidate's
  missed performance target; the revised candidate passes five-pair comparisons
  on 12 and 960 ratings with unchanged public tables, ZIP/HTML content, workbook
  cells, analysis identity and holds. Record the remaining performance scope.
- Translate history controls and table labels in both languages, explain the
  model-selection hold without internal gate labels, preserve technical CSV
  reason codes, and check restore/delete/cancel behavior without refitting.
- Remove correlation-based claims of method validity from saved-run comparisons.
  Summarize finite, uniquely matched estimates by facet; retain element counts,
  reject ambiguous identifiers and flag missing or changed input records.
  Reuse publication holds for reliability cells without modifying saved values.
- Translate comparison controls, tables, plots and explanatory copy; distinguish
  identical history labels with entry numbers. Verify four fitted cases and
  bilingual selection without presentation-triggered estimation.
- Extend report-reuse timing evidence to five fresh-process pairs with a real
  12-rating fit and a changed-input refit. Retain the predeclared targets,
  complete outputs and limits; no new performance backend is introduced.
- Scientific operating-characteristic qualification remains incomplete.

Initial validation on Python 3.12: the native test selection was checked in split runs
and targeted reruns after corrections (3,023 passed, two optional `immer`
comparisons skipped). The environment doctor, release check, 56 self-tests,
quick benchmark and demo-report export also completed. Browser visual acceptance
and formal operating-characteristic studies were not run. The facet boundary
screen is conservative and does not certify general finite-MLE existence.

The follow-up reporting, status-normalization, and coordinate checks passed
595 distinct related tests across focused runs and corrected reruns, including
English/Japanese AppTest checks with an actual held JMLE fit. The 56 self-tests,
environment doctor, compile and release checks also passed. This follow-up did
not repeat the entire native suite or perform browser visual acceptance.

The Japanese-language follow-up passed 464 distinct related tests, including
English/Japanese AppTest checks. After the final data-layout wording change,
30 affected localization, table-conversion and help-integration tests passed
again. Additional AppTest smoke checks exercised Japanese advanced controls
and wide/long table switching with score preservation. The environment doctor,
syntax compilation and whitespace checks passed. The full native suite and
browser visual acceptance were not repeated for this copy revision.

The documentation-reference review checked 30 tracked Markdown documents
outside the validation/archive trees: 579 local link occurrences resolved in
this workspace, while 177 distinct targets were local-only rather than tracked
distribution files. Their availability limits are now stated at the documentation
entry points and research records. Relevant scope, anchor and population-SD
tests passed (48 tests); a direct check confirmed invalid anchor rows are
excluded and audited and the app loads the corrected bundled guide. The doctor
and whitespace checks passed. The two deployment references were checked
against Streamlit's official documentation; other external references and
historical study results were not comprehensively revalidated.

The result/export consistency review passed 155 related tests, including real
ordinary JMLE, facet-boundary JMLE, fixed-SD MML and free-SD MML fits, CSV numeric
round trips, actual Word/PDF/HTML construction and English/Japanese AppTest
report-table checks. Changed Word/PDF table pages were rendered and visually
inspected; 43 affected tests passed again after the Word formatting adjustment.
The environment doctor, release check, syntax compilation and whitespace checks
also passed. These checks do not qualify the estimators' sampling uncertainty or
replace browser/user acceptance; formal validation studies remain paused.

The roadmap-only follow-up checked 30 local link occurrences in the active
plan and updated entry points, preserved all 74 existing roadmap anchors,
and verified that the retained research text and detailed acceptance conditions
were unchanged. Whitespace checks passed. It did not rerun the application,
conduct user reviews, resume formal studies or change deployment status.

The journey follow-up passed 127 distinct related tests: 49 reporting/reuse/
locale checks, 71 engine/quadrature/adapter checks and seven AppTest journeys.
Six pre/post-change EM fits (RSM/PCM/GPCM × fixed/free SD) retained exactly equal
numeric results, traces, termination states and diagnostic tables, excluding
elapsed time. The doctor, release check, compile and whitespace checks passed.
One unmocked 12-row free-SD journey measured repeated PDF display at 2.596 →
0.712 seconds and measures-CSV selection at 23.247 → 6.859 seconds. Word's first
build did not improve; first-render and total times include initialization/cache
effects. See the [scope and timing record](docs/ux_adversarial_audit.md#local-execution-record--2026-10-14).
The full native suite, browser accessibility/visual acceptance, actual-user
reviews and formal statistical studies were not run; no deployment was made.

The TAM follow-up ran two fixed-calibration PCM comparisons and the existing
two-item SE-definition control with TAM 4.3.25 / R 4.6.1. Same-grid PCM NLL,
EAP and posterior SD agreed within the predeclared 1e-8 tolerances; the control
reproduced the distinction between diagonal-curvature and joint-information
SEs. Separately, 51 related Python tests passed and 81 historical evidence-file
hashes matched their manifests. No app implementation or estimator settings
were changed. GH31/grid differences and nonzero returned gradients remain
unresolved; these development checks do not qualify inference or resume the
paused studies. See the [conditions and results](docs/mml_quadrature_cross_engine_review.md#tam-review-20261014).

The fixed-SD follow-up evaluated all 12 planned trajectories (two inputs,
three GH orders, two starts) and 14 points under two continuous-integration
settings. All finite-rule final gradients were below 9.59e-10, while GH31
posterior-SD differences remained about 0.02 logit. GH121 reduced these
differences to about 1.20e-6 without establishing general accuracy. The existing
reference checks passed at all 14 points; 24 related tests and two bilingual
UI/bundled-guide checks passed. Inputs,
sources, every trajectory and hashes were retained. App numerical settings,
original results and inference holds were unchanged, and no new TAM run or
formal study was started. See the [decomposition and limits](docs/mml_quadrature_cross_engine_review.md#tam-quadrature-followup-20261014).

The language-control follow-up evaluated two retained PCM points under GH31,
GH61, GH121 and a 401-point grid. All eight R/Python comparisons met the
predeclared 1e-10 tolerance: maximum absolute differences were 4.27e-13 for
summed NLL, 1.03e-15 for EAP and 8.89e-16 for posterior SD. Prior-only moment
controls also passed. These were direct R evaluations with shared quadrature
arrays, not native TAM refits or independent node-generation checks. The
approximately 0.02-logit GH31 posterior-SD difference persisted in R. See the
[comparison conditions](docs/mml_quadrature_cross_engine_review.md#language-control-20261014).
Two bilingual UI and bundled-guide tests passed after the explanation update.

## 0.2.18-beta — 2026-09-15

- Rewrite Japanese guidance throughout setup, results, diagnostics, and export;
  keep detailed numerical explanations behind the relevant controls.
- Add a Home action that preserves the current data, settings, and results.
  Display the MIT license and provide APA, BibTeX, and CFF software citations
  for Ryuya Komuro. Report citations follow the version recorded in the fit.
- Add current-fit comparison exports and optional offline PCM scripts for TAM
  and ConQuest. The Streamlit application does not launch external engines.
- Offer an on-demand, read-only fixed-SD PCM returned-point check. Optimizer
  termination, gradient stationarity, and quadrature sensitivity remain
  separate judgments; checking a result does not replace its estimates.
- Retain development numerical audits and stable likelihood-difference
  optimization helpers. These do not qualify free-SD uncertainty, establish
  universal quadrature-order recommendations, or change the default estimator.
  Hash-bound study outputs and source archives remain local.

## 0.2.17-beta — 2026-09-13

- Publish the reviewed unified development application: simpler onboarding and
  sidebar, focused result navigation, direct PCA/bias checks, and task-based
  Report & Export with an English APA manuscript template.
- Keep the display-density control label accessible inside its expander and
  refresh the README preview with the built-in synthetic sample.
- Preserve statistical output guards, identified optimizer coordinates, and
  the distinction between optimizer termination and inference readiness.
  Free-population-SD MML uncertainty and inference remain withheld pending
  qualification. Development JAX/NumPyro plans and exact CMLE research do not
  become public estimator capabilities.
- Make release checks reproducible from a clean checkout: generate synthetic
  instrument inputs during tests, use absolute AppTest paths, and reserve
  local study replay for the explicit `retained_evidence` test selection.
  Update stale UI assertions and distinguish methodological references from
  executable external-engine calls in the static boundary check.
- Retain both development and prior public histories without a force push.
  The reviewed unified tree is the release source; remaining public-branch
  performance work is not represented as integrated or validated here.

## Unreleased

- **Known probabilistic assignment validation (repository-only).** Added a
  fixed-margin exponential Person-Rater assignment family, a 90-state exact
  oracle, score-free equal-context row materialization, finite-chain
  diagnostics, and a four-Rater exact dynamic-programming sampler. Two
  prospectively gated 80×4 MCMC settings were retained as failures because
  Statistic ESS was below 300; exact-DP qualification passed and selected
  `|gamma|=0.8` without responses. A separately frozen 10-replicate PCM
  screening completed 120/120 attempts, including 30/30 FACETS 4.5/Python JMLE
  calibration pairs. Only fixed/free normal-person MML Rater RMSE/MAE showed
  non-zero screening intervals under both assignment directions; this remains
  fixed-Person exploratory evidence, not a robustness or estimator-ranking
  claim. Base R 4.5.1 independently reconstructed the Rater RMSE contrasts.
  FACETS rounded fit display values were not used as raw recovery inputs.
  A separately frozen dense-margin exact-DP implementation then matched the
  90-state oracle and the retained recursive 80-by-4 partition functions
  (`2.84e-14` maximum difference).  Its largest partition build was 0.202 s,
  about 466 times faster than the retained same-machine +0.8 recursive
  diagnostic, while preserving an explicit 6,000,000-cell fail-closed cap.
  This qualifies only an execution route; a fresh multi-Person-vector
  preflight remains required before confirmatory endpoint registration.
  That four-vector preflight subsequently prepared 12 exact-assignment
  datasets and completed 48/48 attempt units with all Python JMLE, fixed/free
  MML, and exact-CMLE evidence ready.  The initial FACETS lane failed because
  it used a hidden Windows launch state: direct CMD `BATCH=YES`, Python
  `BATCH=YES`, and even `BATCH=NO` with a hidden `STARTUPINFO` window all
  reproduced `0xC0000005`, while visible `BATCH=NO` completed normally.  A new
  output-gated visible launcher waits for a nonempty stable report and all
  expected Scorefiles, then posts `WM_CLOSE`; its official and historical
  controls passed without forced termination.  A second operational defect
  was then isolated: the supplemental auxiliary Scorefile reached the legacy
  260-character Windows boundary and stalled after Table 5.  The shared helper
  now fails before launch above a qualified 220-character budget, and the clean
  short-path supplement passed 12/12 calibration pairs with zero retries.
  Worst direct differences were 0.012604 logits for main effects and 0.009033
  for PCM thresholds; minimum within-facet Spearman was 1.0.  A non-destructive
  evidence join preserved every original failure and marker, replayed original
  versus supplemental Python JMLE within `4.44e-16`, and created a separate
  preflight aggregate.  All four registered n=4 sign-stability diagnostics met
  their advancement directions; R 4.5.1 rebuilt all 16 vector values and four
  summaries within `1.67e-16`.  The result remains nonconfirmatory and does not
  rank estimators.  FACETS' ordinary two-decimal fit display was never treated
  as unrounded input.  Earlier dump evidence is retained as diagnosis of the
  hidden-window crash path, not evidence that interactive FACETS 4.5.0 is
  generally unusable.

- **Native exact-CMLE hard facet anchors.** Added prospectively specified
  affine hard-anchor coordinates to the repository-only RSM/PCM conditional
  likelihood and propagated fixed row/category offsets through structural
  surfaces, fit-sample WLE, Person MnSq, bootstrap generation/refits, and
  calibration-draw sensitivity. Invalid/duplicate/unknown/non-finite, Person,
  and Step anchors fail closed; anchored estimates are exact with SE zero and
  remaining free-coordinate rank is re-audited. Forty-six selected CMLE/WLE
  tests pass. The retained 24,000-response hard-anchor smoke completed 32/32
  inference-ready same-byte fits and passed anchor, dimension, persistence,
  and raw-fit gates. A transparently post-result diagnostic confirmed that a
  common +0.25 Rater-anchor shift moves WLE by +0.25 but cancels from the
  conditional likelihood and Person fit (maximum log-likelihood difference
  `1.82e-12`; zero raw flag mismatches). This is an explicit warning, not an
  anchor-robustness claim. Added a separately registered 320-fit differential-
  anchor stress over 120,000 retained response rows. All fits were inference-
  ready, but correct nested anchor sets did not improve WLE recovery
  monotonically; anchor share remained confounded with anchor content. A zero-
  mean differential error changed 2.4%--15.8% of paired raw either-upper fit
  decisions by condition/group and moved individual WLEs by up to `0.287525`
  logits. Three-decimal display-driven decisions disagreed with raw MnSq
  decisions 12 times across 64,000 Person-scenario rows; six decimals had zero
  disagreements in this dataset, without becoming a universal precision
  guarantee. Added a prospective randomized-anchor connectivity stress with 70
  unique datasets and 350 exact prefit audits. The realized informative-Person
  graph rule matched exact eligibility in all 350 cases; two random anchors
  covered a two-component design in 7/10 assignments, while partial anchors
  never identified the fully disconnected design. Of 287 prefit-eligible
  cases, 286 returned and 274 were inference-ready. Every shortfall occurred
  in the one-Person-per-edge minimal chain, including one non-finite
  conditional-normalizer trial and 12 final-information rank failures. The
  remaining 200 weak/strong chain/hub fits were inference-ready. Thirteen of
  82,200 Person decisions disagreed when 3-decimal display MnSq replaced raw
  values; six decimals disagreed zero times in this dataset. Added a
  prospectively registered finite-domain optimizer adapter that records and
  rejects non-finite BFGS trials while retaining the exact public objective and
  all readiness gates. Same-byte returnability rose from 286/287 to 287/287;
  readiness remained 274/287. The target recorded 12 invalid trials and
  returned at rank 8/10, nullity 2, with inference withheld. All baseline
  readiness states and 82,200 shared Person flags were unchanged; maximum
  structural and loglikelihood differences were `3.55e-15` and `2.27e-13`.
  Added focused invalid-trial and best-finite-fallback tests; 65 selected tests
  pass. Conditional separation/finite-MLE diagnosis, drift distributions,
  cross-engine replay, group anchors, and public UI remain withheld.

- **Finite-CMLE existence and scalable support oracle.** Added a separately
  prospective convex-support boundary audit that keeps structural rank,
  finite-maximizer existence, optimization, fit, and anchor validity as
  distinct gates. The initial exhaustive implementation matched all 93
  registered binary, multi-Rater, RSM, PCM, unused-category, and anchor
  controls; its deliberately low retained-replay cap classified 63 structural,
  13 boundary, and 47 interior cases while retaining 227 as unavailable. A
  second prospective implementation now uses an exact fixed-score dynamic-
  programming support oracle with fail-closed LP constraint generation. It
  matched 1,184/1,184 exhaustive support maxima, 93/93 fixture statuses,
  123/123 completed retained exhaustive statuses, and 6/6 high-cap sentinels.
  The same-byte 350-case replay resolved 63 structural, 13 boundary, and 274
  interior cases with zero unavailable/tolerance-unstable results and zero
  boundary-direction trace failures. All 274 current inference-ready cases
  were interior. Median theoretical configurations were 97,040 versus 163
  generated cuts; the eligible-case P95/max were 0.521/0.568 seconds and the
  350-case total was 109.216 seconds under the registered 120-second gate.
  A third prospective integration now enables the oracle by default in the
  repository-only `fit_cmle` readiness path. All 287 eligible same-byte fits
  returned; 274 interior cases remained ready and 13 boundary cases were
  non-ready with a typed reason. Oracle status/readiness matched 287/287,
  82,200 Person rows retained every raw/rounded flag, and numerical differences
  were at most `2.27e-13`. The selected integration suite passes 102/102 tests.
  Structured pre-optimizer boundary stopping, matched R boundary behavior, and
  public one-click UX remain withheld.

- **Structured pre-optimization CMLE workflow.** Added a prospectively frozen,
  Streamlit-free bilingual orchestration layer with stable terminal states and
  four ordered stages for input/design, finite-MLE existence, structural
  optimization, and downstream Person-scoring availability. All 93 registered
  fixture states and all 350 retained states matched. The workflow stopped 63
  structurally unidentified and 13 boundary cases before optimization and
  optimized only 274 interior cases; forbidden early-stop optimizer calls were
  zero. Ready structural estimates/SE differed from the direct integrated path
  by at most `4.44e-16`, conditional log likelihood by `2.27e-13`, and anchor-
  exact mismatches were zero. P95 times were 0.068/0.306/1.084 seconds for
  design-block/boundary/ready states. The v1 workflow does not perform Person
  scoring, and machine-checked bilingual text is not a comprehension study.
  Cross-engine boundary replay and public UI remain withheld.

- **Cross-engine finite-CMLE boundary semantics.** Added a prospectively
  registered 13-case additive-RSM replay with five oracle-interior, seven
  oracle-boundary, and one structurally unidentified control. Matched native
  Python/`immer::immer_cml()` interior free coordinates agreed within
  `1.78e-9`, and conditional log likelihoods within `5.33e-15`. On boundary
  controls at `maxit=2000`, `immer` nevertheless returned optimizer code zero
  in 3/7 and finite coefficient/SE vectors in 4/7; the exact support-oracle
  gate kept readiness at 0/7. A content-frozen dirty mfrmr 0.2.3 JML snapshot
  returned and labelled all 7/7 boundary fits converged, with category-support
  warnings and maximum absolute estimates up to `103.75`. TAM and sirt MML
  outcomes are retained as estimator-mismatched sensitivity ledgers only;
  sirt runs are isolated by case so a native-code failure cannot erase other
  engine evidence. Full-precision gradient classifications at `1e-4` through
  `1e-8`, maxit sensitivity, typed support failures, package/function hashes,
  and child-process logs are retained. PCM and hard-anchor mapping is handled
  in the separately registered study below; comprehension testing and public
  UI remain withheld.

- **PCM hard-anchor cross-engine mapping.** Added a prospectively registered
  12-case PCM study covering unanchored, Rater-anchor, Criterion-anchor,
  contaminated-anchor, differential-anchor, separation, and unused-declared-
  category conditions. For six exact-support-oracle interior cases, native
  Python CMLE and matched `immer::immer_cml()` `W`/`b_const` coordinates agreed
  within `2.45e-9`; conditional log likelihood agreed within `3.55e-14`, and
  zero-point objective/gradient agreement was within `4.27e-14`. All six
  boundary cases remained non-ready at three retained iteration caps. Exact
  Python hard anchors returned their supplied values with SE zero, but this is
  an implementation check rather than evidence that the anchors are unbiased.
  mfrmr 0.2.3 JML and TAM/sirt MML remain explicitly different-estimand
  sensitivity lanes; TAM/sirt anchor cases are typed unsupported. The first
  replay exposed path-dependent R source attributes in the function identity;
  a result-transparent remediation removed source references. A second failed
  replay then exposed a registered per-function versus final-vector hash
  mismatch; a frozen amendment corrected only that aggregation bookkeeping,
  and the exact replay passed all gates plus 91 selected tests. High-dimensional
  sparse PCM, clean released mfrmr 0.2.3 reproduction, step/group anchors,
  and public UI remain withheld. Downstream scoring/rounding cards are handled
  by the separate repository-only contract below; comprehension testing still
  remains open.

- **Repository-only one-click CMLE result contract.** Added a guarded single-
  call orchestration layer that returns six stable bilingual cards for input
  and design, finite-MLE existence, structural calibration, fixed-calibration
  Person scoring, Person fit/rounding, and uncertainty/scope. Five registered
  analytical states returned the same 30-card schema: two ready calibrations
  alone reached WLE/MnSq, while input-invalid, structurally unidentified, and
  finite-support-boundary cases made zero downstream attempts. Fourteen ready
  Person handoffs reproduced direct WLE estimates, conditional SEs, Infit, and
  Outfit with maximum absolute difference zero. The anchored PCM case retained
  two supplied main-effect anchors exactly with SE zero and a persistent
  validity warning. A registered `1.5004` raw Infit stayed `noisy` although
  three-decimal display yielded `1.500` and the counterfactual `acceptable`
  label. All card decisions use unrounded values. The final uncertainty card,
  ZSTD/p-values/intervals/total SE, new-Person scoring, automatic fallback,
  task-comprehension claims, and Streamlit exposure remain withheld; 35
  selected tests pass.

- **Deterministic private CMLE result archives.** Added a current-schema,
  Streamlit-free private reproduction ZIP for one-click CMLE results. Ordered
  input, semantic row-order input, canonical settings, result tables, archive
  content, and deterministic ZIP transport receive distinct SHA-256
  identities. Sorted entries use fixed timestamps, permissions, and compression;
  elapsed timing is excluded from result identity, while numeric CSV uses 17
  significant digits. Three registered ready/anchored/boundary archives were
  byte-identical across repeated construction and passed 21/21 exact replay
  checks. All 56 retained WLE/SE/Infit/Outfit values reloaded with identical
  binary64 hex values. A row permutation preserved semantic input identity but
  changed ordered identity and AnalysisID; reversing anchor-row order preserved
  normalized anchor bytes and AnalysisID. Changed, deleted, undeclared, and
  duplicate ZIP entries were rejected before loading. Public construction
  fails closed because the archive contains raw rows and Person/facet IDs;
  encryption, de-identification, public disclosure, cross-version migration,
  cross-platform bitwise identity, and Streamlit upload/download remain
  unqualified. Ninety-five selected tests pass.

- **Non-public bilingual CMLE comprehension instrument readiness.** Added
  standalone escaped English/Japanese six-card HTML previews for five retained
  analytical states, a 100-row participant packet without answer fields, and
  a separately retained 100-row researcher key. Six critical misconception
  domains are non-compensatory, including raw `1.5004` versus displayed
  `1.500` MnSq classification, false-green calibration, conditional-SE scope,
  anchor validity, invented Person scores, and private/public sharing. Thirty
  synthetic packets confirmed that perfect responses pass, all-wrong and
  single false-green responses fail, and duplicate, missing, unknown-option,
  or invalid-time packets are retained as invalid. All automated instrument-
  readiness gates and 40 selected tests pass. No human responses were
  collected; comprehension, language equivalence, real-browser/assistive-
  technology accessibility, and public Streamlit exposure remain unqualified.

- **Private CMLE cognitive-interview operations kit.** Added a deterministic
  English/Japanese complete-pair schedule with 10 planned slots per language.
  Each of five analytical states appears four times per language, twice in
  each presentation position and twice in each prior-MFRM-experience target
  stratum. Generated 20 answer-free participant packet files (400 assigned
  tasks), a separate 100-row private moderator guide requiring response lock
  before probing, and a 400-row blank record template. Exact-column and value
  validation rejects direct/contact identifier fields, assignment drift,
  missing/duplicate items, unknown response/issue codes, inconsistent issue
  severity, unreviewed notes, premature probes, and invalid timing/confidence.
  Fourteen synthetic validator cases and 47 selected tests pass. This is not
  ethics approval or recruitment authorization; no human data were collected
  and the public surface remains disabled.

- **Versioned directional CMLE confirmatory-gate mechanics.** Added a composite
  frozen-instrument identity over participant task bytes, private key bytes,
  ten preview hashes, and directional dangerous-error rules; content aliases,
  post-response freezes, and permission to pool differing content fail closed.
  Distinguished any critical-item error from a registered unsafe response
  direction: `display_acceptable` is dangerous for raw `1.5004`, while
  `cannot_decide` remains incorrect but non-dangerous. Added a two-case-per-slot
  template contributing at most one primary observation per danger domain,
  per-language/item one-sided 95% Wilson gates, raw strict `< 0.10` decisions,
  and a 12-cell Bonferroni sensitivity. Zero-error bounds are `0.101310` at
  `n=24` and `0.097654` at `n=25`; the latter is not a sample-size
  recommendation and its familywise sensitivity remains `0.217782`. Ten
  synthetic gate scenarios and 55 selected tests pass. Human participants, a
  registered confirmatory minimum n, language equivalence, and public UI remain
  unavailable.

- **No-human-data confirmatory sample-size sensitivity.** Added exact binomial
  pass-probability surfaces over the frozen directional Wilson gate, with a
  prospective amendment requiring both first crossing and sustained-through-
  search reporting after preliminary arithmetic exposed discrete sawtooth
  behavior. At true dangerous-error probability `0.05`, a single cell first
  reaches 90% at `n=224` and remains there through `n=5000` from `n=260`; the
  conservative 12-cell union-bound counterparts are `n=422` and `n=456`.
  Added independent homogeneous-retention assurance, cluster design-effect
  heuristics, and blocked-mechanism exposure sensitivities plus two visually
  inspected figures. Sixty-three selected tests and all registered arithmetic
  audits pass. No confirmatory n was selected, no participants were recruited,
  and clustered inference, language equivalence, and public UI remain withheld.

- **Synthetic confirmatory dependence and MNAR stress.** Added a prospectively
  registered Gaussian-threshold generator that preserves cell marginal truth
  while varying participant latent dependence, moderator/site-like clusters,
  domain and blocked-mechanism heterogeneity, and outcome-dependent validity.
  Across eight scenarios, four planned-slot values, and 2,000 replicates per
  condition, 64,000 synthetic joint-gate replicates passed identity, marginal,
  exact-binomial, retention, logic, mechanism-allocation, precision, and
  determinism audits plus 75 selected tests. At n=500, dangerous-outcome
  under-retention changed complete-data truth near 0.10 to an observed rate of
  `0.060315` and a `0.4035` false-reassuring all-12 pass rate; a blocked-
  mechanism hotspot gave `0.4955`, and the combined adversarial condition
  `0.9675`. A floating-point endpoint guard now prevents a theoretical zero
  Monte Carlo Wilson lower bound from becoming a tiny positive value through
  cancellation. Results are constructed stress evidence, not human behavior,
  language equivalence, a selected sample size, or permission to expose UI.

- **Fail-closed private confirmatory protocol preflight.** Added an 18-row
  decision register with five frozen parent decisions and 13 unresolved
  scientific/external-authority blockers; a no-silent-drop, zero-row exact-
  schema denominator ledger; a 15-code invalidity vocabulary; an external
  ethics/governance evidence template bound to protocol content SHA-256; and a
  72-row all-invalid-safe/all-invalid-dangerous partial-identification surface.
  For valid n=500 and observed 5% danger errors, worst-case robustness remains
  at the registered invalid/valid ratio 0.025 but fails at 0.05, where the raw
  Wilson upper bound is `0.118434`. Added a deterministic nine-member private
  ZIP with fixed metadata, internal hashes, bilingual blocked first read, and
  tamper rejection. Eight synthetic ledger attacks and 91 selected tests pass.
  Contract passage means fail-closed software behavior only: ethics approval,
  recruitment authority, selected n, human data, confirmatory results, and
  public Streamlit exposure remain false/absent.

- **Private non-recommending protocol decision workbench.** Expanded the 13
  unresolved blockers into 39 explicitly non-ranked, non-default, non-
  recommended options with strengths, risks, required evidence, and qualitative
  impacts. Added a 32-edge acyclic dependency graph and an implementation-
  preceding amendment clarifying that broad display stages never permit
  prerequisite skipping. Numeric sample-size decision SCI-02 is guarded by
  seven scientific prerequisites. Added a blank 13-row selection template and
  validators for hidden/unknown/cross-decision options, premature child
  resolution, absent numeric assumptions, repository-generated external
  evidence, and outcome-informed decisions. A fully populated synthetic sheet
  may be selection-complete but never makes evidence verified or recruitment
  ready. Added self-contained bilingual read-only HTML and a deterministic
  seven-member private ZIP with tamper rejection. All gates and 105 selected
  tests pass; options selected, n selected, human data, recruitment readiness,
  and public UI remain zero/false.

- **Private SCI-01 estimand comparison memo.** Added a prospectively registered,
  non-recommending comparison of the observed-valid, scheduled-slot composite,
  and dual-with-bounds alternatives. The same 72 exact integer-count scenarios
  are projected into 216 rows using raw strict `< 0.10` decisions; 36 scenarios
  show an observed-valid pass together with an all-invalid-dangerous proxy fail
  and a non-robust dual result. The audit also identifies that the frozen
  surface has neither the complete scheduled-slot denominator nor the
  prospective invalidity-code map needed to calculate Option B. Its displayed
  value is therefore typed as a diagnostic proxy, never as the final
  scheduled-slot estimand. Added 36 unquantified downstream-impact rows, a
  fail-closed external-owner selection template, a bilingual read-only HTML,
  a deterministic seven-member private ZIP, a same-data figure, and tamper
  checks. All gates and 118 selected tests pass. SCI-01 remains unselected;
  evidence verification, sample size, recruitment, human data, and public UI
  remain false/absent.

- **Private SCI-01/SCI-05 invalidity adjudication.** Preserved the 15-code
  vocabulary while fixing only `none` and leaving 14 codes externally
  unresolved. Added mechanism families, owner-review responsibilities, five
  unranked disposition dimensions, and guards against automatic danger
  classification, withdrawal without ethics authority, dangerous duplicate
  records, silent accessibility pooling, unregistered-reason resolution, and
  post-outcome adjudication. A registered 72-by-5 diagnostic surface produces
  360 raw-decision rows; 36 scenarios cross the strict 0.10 gate as the
  fraction of invalid records treated as composite events changes. Added a
  read-only bilingual workbench, figure, deterministic eight-member private
  ZIP, and tamper validation. All gates and 128 selected tests pass; resolved
  codes remain zero and recruitment/public UI remain blocked.

- **Repeated-simulation operating-characteristics scaffold.** Added a
  Streamlit-free deterministic manifest and aggregation contract that preserves
  attempted, failed, unavailable, and eligible denominators; separates power
  from false-positive rate; and reports Wilson intervals, Monte Carlo SE,
  recovery, SE availability, conditional-Wald coverage, and raw/display
  conclusion sensitivity. Added paired null/alternative, sparse/missing, and
  correct/contaminated-anchor Python JMLE smoke conditions, retained tables and
  figures, a UI-ready first-read summary, and an ADEMP-style protocol for the
  subsequent Python CMLE and mfrmr/TAM/immer/sirt adapters. The smoke evidence and
  public application surface remain explicitly pilot-only/withheld.
  The runner also retains exact generated ratings, truth, and anchors through
  atomic staging with per-run identities and byte-level SHA-256 hashes. A new R
  bridge validates/imports the identical rows, invalidates stale bridge output
  after Python regeneration, and records package versions, source identity,
  and primary-function hashes before any R fit adapter is permitted to run.
  Added the first mfrmr 0.2.3 JMLE fit adapter with separately retained
  Python-control and strict-gradient modes, hard-anchor equality checks,
  pre-fit rank failures, focal bias decisions, parameter recovery, and
  code/bundle provenance. The paired smoke comparison exposes optimizer-code
  versus inference-readiness differences, near numerical agreement among 108
  jointly included facet estimates, complete focal strong-decision agreement,
  and mfrmr's pre-optimization rejection of four rank-deficient sparse designs.
  Added a separate matched unanchored native-Python-CMLE/`immer_cml` adapter.
  It uses the same validated rating bytes, an explicit zero-weight declared-
  support sentinel for immer, independently reproduced Python/R conditional-
  information preflights with 16/16 rank/nullity agreement,
  and source/function/script identities. All 12 structurally eligible RunIds
  returned; 96 free coordinates and their SEs agreed within 0.0000000563 and
  0.00000000641 logits respectively, and conditional log likelihoods agreed
  within 0.00000000000477. Readiness remains threshold-sensitive despite every
  immer optimizer code being zero: 12/12 pass at `1e-4`, 7/12 at the primary
  `1e-5`, and 0/12 at `1e-6`. Sparse rank failures and unsupported anchor/local-
  bias scopes remain explicit. A paired null/+0.60 leakage table shows how the
  omitted local interaction perturbs additive CMLE coordinates without
  relabelling that sensitivity as bias detection. No public CMLE option or
  performance claim was enabled.
  Added a separate `sirt::rm.facets()` MML sensitivity adapter with exact input
  and source identities, Task x Criterion virtual PCM items, hard-Rater-anchor
  checks, 30/61-point quadrature modes, convergence/tail-mass/uncertainty
  gates, and Python normalization. All 32 fits returned and 30 were eligible;
  the retained sparse cap failures are visibly excluded rather than converted
  to negative evidence. Raw and mean-centered item-location sensitivity are
  reported separately, as are Person EAP and population-distribution changes.
  The paired contaminated-anchor condition shows near-exact transmission of
  the injected +0.25 shift. Sparse sirt readiness is explicitly described as
  assumption-based MML identification through a common Person distribution,
  not observed connectivity or JMLE/CMLE parity. sirt 4.2-133 information
  criteria and fixed-anchor SE interpretations are withheld. No public sirt
  runtime or estimator-ranking claim was added. An ordered first-read CSV
  provides a future one-click UI projection while keeping every public-surface
  flag false and the final status `Withheld`.
  Added a separate `TAM::tam.mml.mfr()` MML sensitivity adapter using the
  registered additive Criterion-item + Rater + Task + common-step RSM surface.
  It retains 21/61-point grids, hard-anchor preflights, loop and parsed terminal-
  progress convergence evidence, posterior endpoint mass, free-xi uncertainty,
  generalized response surfaces, Person EAPs, and TAM/progress-function hashes.
  All 32 fits returned; 30 were eligible. Two q21 sparse fits stored variance as
  `0.0010000001` at the configured `0.001` lower bound and are excluded by a
  `1e-8` numerical band, exposing a decision that naive exact comparison would
  miss. The comparison separates Python-JMLE/TAM-MML Person treatment, TAM/sirt
  item-threshold and constraint mechanics, and all-returned versus eligible
  quadrature scopes. Contaminated fixed anchors shift +0.25 while unanchored
  Raters compensate by about -0.25 under TAM's sum-zero constraint. Derived
  last-level SE coverage and cross-estimator IC comparison are withheld. Three
  figures and an ordered disabled first-read projection were added; no public
  TAM runtime or cross-engine performance claim was enabled.
  Added a prospective machine-readable precision/stopping plan, exact
  seed-preserving manifest-extension audit, code/plan identities, and a
  separately retained 20-replicate Python pilot. The pilot completed 160/160
  returned fits in 82.1 recorded fit-seconds, with 156 optimizer convergence
  flags and 124 registered app-eligible focal decisions. It blocks further
  expansion: all 160 reconstructed terminal-gradient sup norms exceeded the
  post-pilot `1e-4` sensitivity, and both sparse conditions produced only 4/20
  eligible decisions, implying 6,200 attempts under the registered Wilson-
  lower-bound rule and therefore exceeding the 2,500 cap. Row-level auditing
  retained 39,680 fit-statistic evaluations, 14 display-boundary rows, and five
  raw/display label disagreements. Paired anchor contamination shifted fixed
  Raters exactly +0.25 on average. Four review figures and an ordered disabled
  first-read assessment were added. The frozen v1 continuous-outcome planning
  conflict is documented rather than silently repaired; strict Python JMLE
  qualification and sparse design/estimand revision remain required before a
  confirmatory manifest or public surface.
  Added prospective Stage A/A2/B2 strict-JMLE follow-ups without rewriting the
  frozen pilot: the original 16-run baseline confirmed correct analytical
  gradients but failed the raw terminal-gradient gate; the predeclared
  L-BFGS-B precision candidate passed all 16 smoke runs and then all 160 frozen
  numerical gates at `1e-4` (126 at `1e-5`, 43 at `1e-6`). A separate movement
  audit exposed flat-direction changes up to 57.02 logits. Exact free-coordinate
  eta-rank auditing then showed structural nullity 7 and eight Person-Rater
  components in all 40 sparse runs, with null-space energy concentrated in
  Person and Rater coordinates. Added a pre-optimizer JMLE guard that reports
  rank/nullity/connectivity, separates optimizer convergence from inference
  readiness, withholds conditional bias for rank-deficient fits, and exports
  the audit. Its bridge check matched all 160 retained RunIds. Optimizer
  controls and estimates remain unchanged, no automatic JMLE-to-MML switch was
  introduced, and the eta-only audit does not authorize a public performance
  claim or qualify threshold/slope identification.
  Added a separately registered JMLE Person score-boundary guard. It classifies
  all-minimum/all-maximum Persons from exact retained integer scores rather
  than terminal-estimate magnitude or display rounding; retains the existing
  optimizer/constraint value as technical `Estimate`; withholds
  `ReportableEstimate`; and adds finite-MLE, direction, Person-readiness, and
  value-role fields to results, diagnostics, bilingual UI warnings, quick/full
  bundles, demo exports, and APA tables. The frozen integration check matched
  8,320/8,320 Person rows. Among 120 structurally identified runs, four extreme
  Persons accounted for every >=1-logit movement, while maximum non-theta
  movement was 0.000293 logits. The 53 extreme rows in sparse runs remain
  additionally blocked by structural nonidentification. No estimates,
  optimizer controls, finite extreme-score correction, estimator switch, or
  precision-polish core path changed.
  Added a separately registered fixed-calibration Warm WLE research core with
  stable generic category logits, positive row weights, Person-specific
  missingness, analytical score/information/information-derivative moments,
  bracketed Brent roots, exact-extreme handling, and explicit no-root states.
  Python theta and conditional-information SEs matched `TAM::tam.mml.wle2` in
  16/16 frozen RSM/PCM/GPCM Person fixtures within `1.68e-11` and `1.18e-11`.
  The original unreproducible TAM function hash remains visible; a pre-result
  amendment records two reproducible source identities without changing any
  numerical gate. A downstream replay held polished facets and steps fixed for
  120 rank-full runs and 7,600 Persons while withholding all 40 structurally
  deficient sparse runs. Exact-extreme Person shifts reached 25.19 logits;
  51 Person-fit-zone decisions changed, and `.3g` display rounding disagreed
  with 54 raw fit decisions across the JMLE/WLE surfaces. Holm and combined
  strong focal-bias decisions were stable, but one practical `|bias| >= 0.50`
  decision moved from 0.500839 to 0.498068. Retained tables and two reviewed
  figures expose these differences. No Streamlit WLE option, automatic JMLE
  correction/fallback, or public parity claim was enabled.
  Added a repository-only exact-CMLE-to-WLE fit-sample bridge with an explicit
  combined-estimator label and fail-closed calibration/parameter-identity
  checks. Frozen RSM/PCM evidence covered 14 Persons including two exact
  extremes: bridge/direct WLE theta and SE differences and reconstructed
  CMLE category-surface differences were all zero, with maximum adjusted-score
  residual `6.18e-16`. Coefficient order was invariant and missing identity was
  rejected. WLE uncertainty remains conditional on fixed fitted CMLE
  calibration; new-Person scoring and UI integration remain withheld.
  Added a prospectively frozen asymptotic exact-CMLE calibration-sensitivity
  diagnostic for the WLE bridge. RSM and PCM each used 2,000 seeded
  free-coordinate covariance draws; all 32,000 draw-Person scores returned,
  exact extremes remained finite, and covariance scaling produced the required
  monotone response. Model-median calibration-draw SD was `0.0616`--`0.0801`
  logits, the maximum was `0.518`, and the largest ratio to conditional WLE SE
  was `0.427`. Material covariance non-symmetry and non-PSD inputs now fail
  closed. Draw quantiles remain explicitly non-confidence intervals, and a
  quadrature sensitivity value is not labelled a total inferential SE because
  same-sample CMLE/WLE dependence is omitted. A retained six-card first-read
  projection separates research-ready, caution, and withheld states; it does
  not enable the public UI.
  Added a prospective CMLE-WLE bootstrap design and immutable input-identity
  plan before sampler implementation. The plan separates fixed-score
  conditional-pattern refits from joint plug-in response refits, requires all
  failures and rank states to remain in the denominator, and prohibits calling
  either lane a coverage-qualified interval or total SE. Its 200-replicate
  RSM/PCM matrix is an implementation pilot with no success-rate promotion
  threshold; a later repeated-truth ADEMP protocol must be registered before
  interval claims.
  Implemented the first repository-only bootstrap core after that plan was
  frozen. It uses a scaled forward/backward coefficient recurrence for exact
  fixed-Person-total RSM/PCM response-pattern draws, a separate fitted CMLE +
  WLE joint category sampler, deterministic 64-bit child seeds without float
  coercion, pre-refit rank/nullity capture, and complete failure ledgers with
  work/output guards. Six tests pass, including brute-force conditional-
  normalizer agreement, 20,000-draw distribution agreement, exact RSM/PCM
  total preservation, joint total variation, seed reproducibility, and
  resource failure. The frozen 200-replicate RSM/PCM pilot then completed all
  800 attempted refits with full conditional rank, CMLE inference readiness,
  and WLE availability; the Wilson 95% lower bound is `0.981` per cell.
  Fixed-score totals had zero mismatches, same-seed replay was exact, and joint
  category frequencies were within `1.57` Monte Carlo SD of expectation.
  Median Person bootstrap SD was `0.0602` logits for fixed-score resampling and
  `0.615` for joint plug-in resampling; maximum SD was `0.744`. The joint lane
  produced 321 extreme-to-interior and 120 interior-to-extreme transitions;
  maximum bootstrap-mean displacement was `0.370` logits. Sampler status
  passed, but output-contract/public promotion remains withheld because
  repeated-truth coverage is unavailable.
  Added a separately registered fixed-theta CMLE-WLE Person-fit core and
  Python/`sirt::pcm.fit` adapter. RSM/PCM evidence compared 32 Persons and 162
  administered observations, including exact extremes and non-rectangular
  missingness; maximum Infit and Outfit differences were `4.89e-15` and
  `4.00e-15`, with row-level probability moments agreeing at binary64 noise.
  Non-positive/non-finite fitted variance fails closed, unrounded MnSq drives
  the existing 0.50/1.50/2.00 decision contract, and ZSTD/p-values remain
  explicitly unavailable. TAM shares the Infit and pre-trim Outfit formulas
  but its default JML path trims large Outfit contributions, so default Outfit
  equivalence is not claimed.
  Added a prospectively frozen 800-replicate bootstrap fit extension without
  modifying the original pilot. Every replicate was Person-fit ready; 12,800
  Person-replicates produced 2,372 Infit transitions, 2,454 Outfit transitions,
  and 2,554 transitions in either statistic. The raw/display audit retained 24
  display boundaries and nine classification mismatches; all transition
  decisions used unrounded values. The first run remains visibly failed after
  an asymmetric CSV-versus-memory identity comparison introduced a one-ULP
  difference against the frozen exact-zero gate. A pre-correction amendment
  kept the gate/tolerance unchanged and required a full rerun with symmetric
  `%.17g` persisted readback. The corrected run passed with zero WLE difference,
  zero fixed-score total mismatches, and unchanged original/failed inventories.
  Added a boundary-focused reviewed figure so the 24 rare rows are not hidden
  beneath 12,000+ stable rows. Threshold optimality, ZSTD, p-values, coverage,
  total SE, sparse/anchor/misspecification stress, and Streamlit integration
  remain withheld.

- **CMLE-WLE MnSq rounding and threshold sensitivity.** Added a research-only
  raw-value threshold classifier, complete 125-triplet post-result sensitivity
  surface, one-through-six-decimal counterfactual display audit, threshold
  proximity/concentration evidence, and reviewed figures. Canonical raw
  classifications and transition counts reproduced exactly. At three display
  decimals, nine replicate-statistic classes and five either-statistic
  transition indicators disagreed with raw decisions; six decimals had no
  disagreements. Aggregate transitions spanned 1,789--3,461 for Infit,
  1,877--3,550 for Outfit, and 1,975--3,675 for either statistic across the
  frozen grid. A pre-computation diagnostic amendment added complete 16-cell
  baseline-by-replicate class matrices after a flat noisy-boundary profile was
  observed. Moving that boundary from 1.90 to 2.10 preserved the binary totals
  2,372/2,454/2,554 but relabelled 37 Infit and 51 Outfit replicate statistics
  from `distorting` to `noisy`. Parent evidence remains immutable; threshold
  validity, automatic selection, known-truth operating characteristics, and UI
  integration remain withheld.

- **Known-truth CMLE-WLE Person-fit design pilot.** Prospectively registered
  and ran eight paired fixed-calibration conditions over 20 replicates each,
  retaining all 160 attempts, 624,000 generated responses, and 32,000 Person
  rows per WLE/generating-theta source. The canonical raw either-upper clean
  rate was 3.4%--3.7% with 24 observations and 10.1%--10.8% with six;
  affected detection ranged from 3.8% to 46.6% across threshold-heterogeneity,
  local-copy, and random-response mechanisms. A post-result registered
  addendum reports unrecentered WLE bias/MAE/RMSE, retains 163 exact-extreme
  rows separately, and records 2,333/32,000 WLE-versus-generating-theta flag
  disagreements. It also verifies that the nominal 125-grid has only 5/5/25/5
  effective input configurations for the four rules. Same-seed replay, truth
  identity, raw-class reproduction, support, and irrelevant-threshold
  invariance passed with zero mismatches. This remains a 20-replicate
  fixed-calibration design pilot; anchor/connectivity, estimated calibration,
  confirmatory precision, repeated cross-engine fits, and public UI remain
  withheld.

- **Threshold-stable fit/bias decisions and design stress.** Added a pure
  floating-point decision contract that classifies from unrounded values,
  fixes exact MNSQ endpoint semantics, and exports numerical-boundary,
  display-rounding-boundary, and raw/display mismatch evidence. Unified fit
  table colours, summaries, casebooks, rater dashboards, report prose, and fit
  scatter annotations on the same contract; conditional bias audits now
  expose threshold-sensitive cells, while missing fit statistics can no longer
  produce a `Ready` stability summary. Quick results and full table downloads
  include detailed fit/bias decision-stability tables. Anchor audits report unique coverage counts,
  unanchored counts, descriptive share, and a stacked chart without imposing a
  universal percentage cutoff. Added a deterministic nine-design sparse-data,
  anchor, fit-boundary, and bias-boundary stress matrix with retained CSV and
  PNG evidence. Non-finite group-anchor values are rejected with an explicit
  audit warning instead of being silently replaced by zero. Essential results now show a compact run-readiness snapshot
  before the one-click goal route.
- **Repository-only native exact CMLE structural calibration.** Added a Streamlit-free RSM/PCM
  conditional maximum-likelihood core with explicit rating support, scaled
  response-polynomial moment dynamic programming, an analytical gradient, exact
  conditional information, typed extreme-Person states, and exact sum-to-zero
  facet/threshold coordinates. The pre-fit audit blocks GPCM, non-unit row
  weights, unlabelled duplicate response units, no-information designs, and
  deficient conditional rank. Exact enumeration and matched
  `immer::immer_cml()` RSM/PCM fixtures are covered in `tests/test_cmle.py`.
  Rank and covariance now come directly from conditional sufficient-statistic
  covariance rather than finite differences; the reused 20-dataset pilot kept
  all readiness decisions while reducing recorded median fit time to 0.072 s.
  Added a conservative 512 MiB analytical-information memory guard alongside
  the work and parameter caps. A 16-case fresh-process scaling matrix covers
  up to 5,000 Persons, 64 response units, 30 fitted parameters, and 100
  missingness patterns, with explicit work/memory negative controls. Guarded
  exact-Newton polishing resolves above-threshold BFGS precision-loss endings
  without overwriting the raw optimizer status.
  This research core is not yet exposed in the Streamlit estimator selector;
  CMLE Person scoring, anchors, downstream diagnostics, scaling evidence, and
  public reporting remain gated.
- **S1 exact identified parameterization.** Replaced redundant centered step
  and GPCM log-slope optimizer blocks with `n-1` free coordinates whose final
  value is derived as the negative sum. RSM, PCM, and bounded GPCM now count
  identified optimizer dimensions in `KParams`; analytical gradients apply the
  exact expansion Jacobian. Bounded GPCM uses a linear SLSQP constraint so the
  derived final log-slope obeys the same configured bounds without narrowing
  the feasible sum-zero space. Expanded result tables remain compatible, and
  quick/full/APA/reproduction exports include `parameterization_audit.csv`.
  Parameter-recovery rows now also record `IncludedInSummary` and
  `SummaryBasis`, making convergence-based replicate inclusion auditable.
- **S0 statistical output protection.** Added a versioned, fail-closed output
  qualification contract with stable remediation reason codes. AIC/AICc/BIC
  and derived weights remain technical-only; automatic model recommendations,
  model-choice LR decisions, rank-deficient structural covariance claims, and
  pairwise bias local measures are withheld. Raw diagnostic values remain in
  explicitly qualified audit exports, while Report, run-history, comparison,
  quick-download, full-download, and APA export surfaces carry the same status.
- **Optional five-step sample guide.** Replaced the overlapping first-run
  banner, static three-step explanation, and always-rendered tutorial with one
  focused landing decision and a session-scoped five-node guide. The guide
  checks sample data roles before running the production estimator, requires a
  bounded interpretation answer, and ends with a non-scientific learning
  checkpoint. Skip, exit, resume, restart, locale, view, and Help events are
  governed by a pure `mfrm_app.guidance` reducer that cannot change scientific
  readiness. The ordinary workspace remains directly available, and starting
  with user data opens the privacy-aware paste preflight rather than forcing a
  sample tour.
- **Persistent result navigation without shortcut overhead.** Removed the
  sidebar keyboard-shortcut cheat sheet and its Streamlit maintenance commands.
  Essential and All-panels result selectors now live in a sticky navigation
  dock that remains available while long result sections scroll. On narrow
  screens, the Essential switcher stays on a compact touch-scrollable line
  instead of expanding into a tall fixed overlay.
- **Phase-aware UX hierarchy.** Added a pure presentation policy for data
  origin and workflow phase, contextual privacy severity, action-first
  onboarding, collapsed raw-row previews, and a unified pre-run workspace.
  Successful fits now clear repeated setup content, retain one response-row
  audit, and foreground one recommended result action while keeping alternate
  routes, progress, and claim boundaries in supporting detail. The adversarial
  baseline and long-term acceptance measures are recorded in
  `docs/ux_adversarial_audit.md`.
- **Workspace-first guided Run action.** Moved the Guided-defaults primary Run
  control from the technical sidebar into the main setup workspace, after
  column mapping and readiness evidence. Its bilingual user-facing label is
  `Run this analysis / この設定で分析する`, with the selected model, method,
  and analysis depth summarized immediately above it. Guided mode no longer
  renders a duplicate `FACETS-mode` sidebar action; Advanced controls retain
  the existing sidebar route. The estimator inputs, defaults, result identity,
  stale-result rerun behavior, and downloads are unchanged. English/Japanese
  placement checks, 27 focused localization/onboarding tests, and 17 full-fit,
  Help-return, and small-data AppTests passed.
- **Task-first two-level data source chooser.** Replaced the ten-choice flat
  source radio with four user-task classes (built-in example, generate, paste,
  upload) and an always-visible scenario selector when the example class is
  active. The analysis-facing `scenario:*`, `simulate`, `paste`, and `upload`
  IDs remain unchanged. Added migration from older `data_source_flat` session
  state, immediate privacy-severity projection from the new class, and stable
  restoration of the last selected sample after visiting another source.
  Thirty focused state/localization/UI tests and 24 full-fit, Help, simulation,
  and small-data AppTests pass.
- **First-run decisions separated from technical controls.** Guided defaults
  now shows only Person/Score/facet mapping, model, estimation method, and
  analysis coverage before the main Run action. Weighting, visualization CI
  and styling, facet regularization, latent regression, score-scale overrides,
  identification/optimizer settings, anchors, anchor audit policy, and report
  scaling are available in Advanced controls. Guided mode explicitly restores
  documented standard values when returning from an advanced run, so hidden
  custom settings cannot silently affect a new guided analysis. A detected
  weight column is excluded from facet suggestions and explained rather than
  being silently applied. Nine dedicated UI-contract AppTests and 47 broader
  scenario-fit, guide, small-data, locale, and contextual-Help tests pass; an
  additional end-to-end check confirms Guided and Advanced defaults reproduce
  the same parameter estimates to `1e-12` absolute tolerance.
- **Fail-closed browser accessibility acceptance contract.** Added a pure,
  versioned catalog that expands eight user tasks across English/Japanese and
  seven keyboard, narrow, touch, zoom, reduced-motion, and screen-reader
  profiles into 52 browser cases and 482 required evidence rows. Blank or
  incomplete evidence remains `NOT_ACCEPTED`; `PASS` requires every declared
  portable evidence type plus browser-build, UTC, observation, and reviewer
  metadata, and stale catalog fingerprints are rejected. Added export/evaluate
  CLI commands and a browser runbook. Static preflight also broadens visible
  focus styling to links, disclosures, and common ARIA roles, adds
  forced-colors/coarse-pointer handling, and gives compact HTML tables an
  accessible name with row/column header scopes. Eleven contract tests and 29
  combined accessibility/UI regression tests pass. No browser or screen-reader
  acceptance claim is made before the completed evidence gate passes.
- **Standalone Python evidence and release boundary.** Added deterministic
  AnalysisIdentity/EvidenceRecord contracts, canonical readiness states, and a
  fail-closed MML prior-SD sensitivity decision. Default downloads, demo
  archives, `--release-check`, `--self-test`, `make verify`, and GitHub CI now
  use only Python-native evidence and reproduction assets.
- **Legacy compatibility isolation.** Removed Posterior Viewer and Stan
  generation from public sidebar, Report, Help, and tutorial routes; removed
  the parity-fixture CLI/Make/CI route; and placed dormant compatibility tests
  behind the non-release `legacy_compat` marker. Public Yardstick and
  rating-scale recode bundles now take Python-only generator branches directly.
- **Demo handoff accuracy.** Demo checklists and manuscript handoff prose now
  name the actual `MFRM_Demo_*`, `method_appendix.md`, and
  `manuscript_template.md` artifacts produced by the deterministic exporter.
  The same demo-only artifact profile now covers generated tables, Markdown,
  HTML, visual evidence, and binder README files. The exporter also writes an
  `export_privacy_manifest.csv` that explicitly identifies the built-in data
  as synthetic and does not transfer that status to user-supplied data.
- **FACETS-primary fit standardisation.** RSM and PCM diagnostics now default
  to Wright-Masters fourth-moment d.f., Wilson-Hilferty ZSTD, and the FACETS
  absolute cap of 9. Primary `DF_Infit` / `DF_Outfit` and ZSTD columns follow
  FACETS while the previous engine convention is retained in `*_ENGINE`
  sidecars. Category diagnostics, FACETS-style report exports, and generated
  Python/R reproducibility scripts use the same profile. The Fit Details UI
  starts on FACETS-primary but retains explicit engine and engine-primary
  comparison choices. Bounded GPCM outputs are marked
  `facets_style_approximation_for_gpcm` and display a warning that exact
  FACETS equivalence is not claimed.

## 0.2.15-beta - 2026-06-17

### Added

- **README and version alignment for public distribution**. The app version is
  bumped to `0.2.15-beta`, and the README now documents the current beta label,
  release-readiness commands, standalone runtime boundary, demo/parity export
  smoke checks, and the new `mfrm_app/` helper package.
- **Guided UX and export-readiness documentation surface**. The README now
  reflects the Essential-view reading-order guidance, localized help popover
  overlays, modular helper boundaries, privacy/cache expectations, and
  reproducibility/export contracts added in the current public beta line.
- **Zotero-informed method reference audit**. A safe read-only audit of
  the local Zotero/BibTeX library was used to expand the in-app APA
  reference library with Rasch, fit, local-dependence, external-package,
  and recent MFRM extension anchors. New `build_method_reference_audit()`
  exports `method_reference_audit.csv`, mapping model, estimation,
  fit/person-fit, residual PCA, bias/local interaction, simulation,
  and external-validation surfaces to citation tokens, BibTeX/RIS export
  status, manuscript use, and claim boundaries. Current-run reference
  rows now feed the APA report, manuscript template, and
  `claim_to_evidence_matrix.csv`, so citations and statistical claim
  boundaries stay aligned with the analyses actually computed. Quick,
  download, demo bundles, and the manuscript binder include these tables.
- **APA sentence evidence audit**. Generated APA-style prose now has an
  exportable `apa_report_sentence_audit.csv` that maps each draft
  Method/Results sentence to the manuscript claim area, output evidence,
  evidence files, suggested citations, citation boundaries, and copy
  decision. The audit is visible in the APA Report and Claim Guide areas
  and is included in quick/download/demo/manuscript-binder exports.
- **Final-result manuscript handoff export**. Downloads now include
  `mfrm_manuscript_handoff.md` and
  `manuscript_handoff_checklist.csv`, a concise file-opening,
  download, submission-preflight, privacy, and archiving guide that
  links the Publication Document, table bundle, figure bundle, method
  appendix, manuscript template, config JSON, app-engine runner, and
  OSF package. The same assets are included in the OSF ZIP and
  synthetic demo report export.
- **Claim-to-evidence matrix and manuscript binder**. The Downloads
  tab and demo export now include `claim_to_evidence_matrix.csv`,
  which maps each manuscript area to claim status, gate status,
  exported tables, figures, readiness checks, reviewer questions,
  caveats, and archive files. A curated
  `MFRM_Manuscript_Binder.zip` / `MFRM_Demo_Manuscript_Binder.zip`
  collects the handoff, evidence matrix, action/gate tables, method
  appendix, manuscript template, config, and figure manifest for
  coauthor review and submission preparation.
- **Visual evidence binder**. Figure downloads now include
  `MFRM_Visual_Evidence_Binder.zip`, with the actual PNG/HTML figure
  files, `visual_evidence_map.csv`, `visual_caption_drafts.md`,
  `figure_manifest.csv`, visual interpretation and method-evidence
  tables, and the claim-to-evidence matrix. The demo report exports
  `MFRM_Demo_Visual_Evidence_Binder.zip` and
  `visual_evidence_map.csv` as well.
- **SE/CI coverage diagnostic for ADEMP recovery**.
  `evaluate_parameter_recovery(...)` now preserves `SE_Method`,
  `SE_Status`, `CI_Method`, `CI_Status`, and `UncertaintyCaution` in the
  recovery long table, and attaches a separate `coverage_report`.
  `build_se_ci_coverage_report(...)` and `evaluate_se_ci_coverage(...)`
  expose nominal coverage, empirical Coverage95, coverage error,
  SE/CI status summaries, interpretation, and recommended action for
  manuscript caveats. The Report tab displays this diagnostic under the
  parameter-recovery panel and offers `mfrm_se_ci_coverage_report.csv`.
- **Executable MML fixed-prior-SD sensitivity**.
  `evaluate_mml_prior_sd_sensitivity(...)` now refits the current MML
  result across selected `population_prior_sd` values, reusing the
  current result for the base SD and returning `summary`,
  `measure_deltas`, optional `population_deltas`, and an interpretation
  `report`. The Fit Details tab exposes this as an on-demand screen with
  a ZIP download so final reports can document whether person EAPs,
  non-person facet measures, ranks, and latent-regression coefficients
  are stable to the fixed population variance setting.
- **Conditional bias-inference audit**.
  `build_bias_inference_audit(...)` now summarizes every exported
  bias/local-interaction facet pair with cells screened, DFF flags,
  strong-review cells, sparse cells, Holm/BH counts, practical-logit
  counts, profile-CI status summaries, inference-tier summaries, and
  connected-scale status. DFF tables now preserve the conditional
  likelihood/SE/multiplicity basis from `estimate_bias_interaction(...)`
  and add cell-level claim status plus reporting cautions. The
  first-read guide, final-readiness gate, manuscript claim guide,
  APA-style narrative, quick result bundle, Downloads ZIP, and demo
  exports now use this audit instead of raw `|t| >= 2` counts.
- **Residual PCA stability sensitivity audit**.
  `compute_pca_bundle(...)` now augments the basic residual-matrix
  checks with leave-one-column-out EV1/share/loading stability and a
  deterministic row-bootstrap audit of EV1, EV1 CV, 5th/95th percentile
  EV1, loading correlations, and EV=2/3 threshold crossings.
  `collect_pca_stability_tables(...)` exports overall and per-facet
  stability rows together as `pca_stability_all_scopes.csv`; the
  Dimensionality panel now displays and downloads the stability audit
  alongside the scree plot and loadings.
- **MML covariance and SE/CI claim audit**.
  `compute_mml_parameter_covariance(...)` now attaches an `audit` table
  with Hessian rank, rank deficiency, eigenvalue floor tolerance,
  regularized eigenvalue count, condition number, covariance diagonal
  range, claim status, and recommended action. Normal diagnostics export
  the same table as `mml_covariance_audit.csv`, and measure rows now
  carry covariance status, claim status, condition number, rank, rank
  deficiency, and regularized-eigenvalue metadata. ADEMP SE/CI coverage
  reports now add `SEBasisRisk` and `CoverageClaimStatus` so coverage is
  not cited without the matching SE basis and covariance caveat.
- **Visualization preferences and visual QA preflight**. The sidebar now
  lets users choose figure theme, label-density policy, base font size,
  label truncation length, static export width, minimum figure height,
  and caption-detail level. These settings are saved in
  `visualization_settings.json` and `visualization_settings.csv`,
  included in config JSON inspection via `visualization_preferences`,
  propagated into `figure_manifest.csv`, `visual_evidence_map.csv`,
  and caption drafts, and checked in
  `visual_qa_preflight.csv` before manuscript use.
- **BibTeX and RIS export of the cited reference library**. The Report
  tab's APA-style report now includes a "Bibliography downloads"
  expander with two buttons that emit a `.bib` (LaTeX / Zotero /
  BibDesk) or `.ris` (EndNote / Mendeley / Zotero) file containing
  exactly the references cited in the generated narrative. The
  conversion is heuristic — `_parse_apa_entry` reads each APA prose
  string in `_APA_REFERENCE_LIBRARY` and classifies it as
  `@article` / `@book` / `@incollection`, extracting authors / year /
  title / journal-or-publisher / volume / number / pages / DOI from
  the canonical APA 7 layout. RIS uses the matching `TY  - JOUR / BOOK
  / CHAP` records with one `AU  -` line per author. New helpers:
  `apa_entry_to_bibtex(key, apa)`, `apa_entry_to_ris(key, apa)`,
  `build_bibtex_from_cited(text, always_include=...)`,
  `build_ris_from_cited(text, always_include=...)`,
  `_cited_reference_keys(text, always_include=...)`. The author
  parser handles single authors, 2–4 authors with `&` joiners,
  hyphenated initials (`K.-Y.`), two-initial given names (`A. A.`),
  and `Jr.` suffixes. Titles ending in `?` / `!` are preserved (the
  Linacre 2002 RMT short note is one such case). All 66 entries in
  the live library parse without falling back to `@misc` — a
  regression test pins this so a future APA addition whose format
  breaks the parser fails the suite immediately. 24 unit tests in
  `tests/test_bibtex_ris_export.py` cover author conversion, every
  citation type, the question-mark-title edge case, the
  Winsteps.com-publisher case, proceedings entries with no
  publisher, and parity between the BibTeX bundle and the APA
  reference list. New i18n namespace `apa_report` in both
  `locales/{en,ja}.json`.
- **Configurable visualization CI level**. The classical 95 % Wald CI
  band on forest plots and EB-shrinkage error bars is now a
  user-selectable level: 50 % / 66 % / 80 % / 89 % / 90 % / 95 % /
  99 %. Sidebar control labelled "Visualization CI level" lives in
  `run_facets_mode` next to Workflow mode; selection persists in
  `st.session_state["viz_ci_level"]` and propagates to
  `_draw_measures_forest_plotly` (outer CI band + caption + title)
  and `_make_eb_shrinkage_figure` (raw and shrunk SE error bars + legend +
  title). New helpers `get_viz_ci_level()`, `_ci_z(level)`, and
  `_ci_level_pct_label(level)` centralise the level → critical value
  conversion via `Phi^{-1}((1 + c) / 2)`. The numeric `CI_Lower` /
  `CI_Upper` columns in exported measure tables stay at 95 % for
  backwards compatibility with downstream parsers — only the rendered
  visualizations change. The three pre-existing per-feature CI
  selectboxes (DIMTEST p-value, G-study bootstrap, GPCM fair-average
  SE) also broaden from `[0.90, 0.95, 0.99]` to the full
  `VIZ_CI_LEVEL_OPTIONS` list. 22 unit tests in
  `tests/test_viz_ci_level.py` pin the `Phi^{-1}` values for all seven
  canonical levels, the soft-fallback to 95 % on NaN / out-of-range
  / non-numeric input, the percent-label trailing-zero formatting,
  and the sidebar selector defaults / clamping.
- **Smart column-role detection** for the upload UI. New
  `auto_detect_column_roles(df)` helper combines a multilingual keyword
  dictionary (English + Japanese: 学生 / 受験者 / 評価者 / 採点者 / 観点 /
  項目 / 得点 / 評定 / ...) with dtype and cardinality cues to pre-fill
  the Person / Score / Rater / Criterion / Weight pickers in
  `run_facets_mode`. Score detection prefers numeric small-cardinality
  columns (≤ 12 distinct values), Person prefers high-cardinality
  non-numeric columns, and a confidence caption surfaces under each
  selectbox when the auto-detector picked the active value. Falls back
  to the legacy keyword-only `guess_col` and `pick_default_facets`
  when no Japanese / English keywords match, so existing English-only
  CSV workflows are unchanged. 19 unit tests in
  `tests/test_column_role_detection.py` cover canonical English / JA
  headers, mixed-language headers, dtype/cardinality bias,
  Weight-only-when-named, no-double-assignment, and confidence-in-unit
  interval. New i18n keys in `locales/{en,ja}.json` keep both locales
  in parity for the auto-detection captions. Companion fix:
  `_column_is_numeric` now uses a direct ratio (`coerced / non_blank
  >= threshold`) instead of `int(threshold * len(...))`, which was off
  by one on small columns and let a 2-of-4-numeric column pass the
  70 % threshold.
- Bilingual i18n scaffold (English / Japanese) for the sidebar header, app
  title, language picker, and starter messages; translation dictionaries live
  in `locales/{en,ja}.json` and resolve via the new `t()` helper with English
  fallback. Persists the choice in `st.session_state["lang"]`.
- mfrmr 0.1.6 migration coverage table and exports, extending the prior
  0.1.5 map with EB shrinkage advisory, facet sample-size adequacy,
  nesting/crossing audit, intraclass-cluster-ICC design effect,
  information curves, and misfit/weighting audit support.
- Full publication mode now computes non-person facet EB shrinkage as a
  post-hoc reporting advisory; Standard mode leaves it optional/off by
  default.
- `LICENSE_NOTICE.md` clarifies MIT licensing, commercial-use permission,
  no-warranty/as-is status, and responsible-use boundaries.
- Teacher-facing paste-data guidance now appears in the sidebar, onboarding
  text, README, and Help quick start, with a downloadable classroom template.
- Upload parsing now accepts Excel `.xlsx/.xlsm`, Apache Parquet, and JSON /
  JSON-lines in addition to CSV/TSV/TXT for rating data, person data, anchors,
  and held-out prediction rows.
- Bias/interaction output now includes a DFF-style screening table with
  sparse-cell flags, Holm adjustment, BH/FDR adjustment, practical logit
  thresholds, interpretation text, and CSV/export-bundle coverage.
- Advanced Stan-model downloads now include model-specific data templates and
  clearer scope notes for local dependence and Mixture Rasch.
- Anchor/linking output now includes an anchor/equating workflow checklist for
  current-run linking evidence and export review.
- **GPCM fair-average kernel test suite** (`tests/test_gpcm_fair_average.py`,
  13 tests). Pins reduction to the PCM/RSM softmax reference at slope = 1,
  the worked example at slope = 1.2 (1.8381511), the analytic derivative
  identity `dE[X]/d eta = a * Var[X|eta]`, invariance under the GPCM
  identification rescaling, the NaN-on-degenerate-slope contract, binary
  and singleton boundary fixtures, and the slope-1 / per-level-slope wiring
  in `calc_facets_report_tbls()`. Includes a 27-case parity fixture
  generated from the mfrmr R kernel (mfrmr 0.2.0, R 4.5.2) at 1e-10
  agreement.
- **Naming-hygiene gate** (`tests/test_naming_hygiene.py`). Parameterised
  test that scans code, tests, locales, validation docs, project docs,
  CI workflows, and project metadata for internal milestone codenames
  (numeric, single-letter, and a finite vocabulary of multi-letter forms
  such as NATO phonetics and spelled-out numbers). Public commits must
  not ship internal sprint or phase labels.
- **MML observed-information covariance** (`compute_mml_parameter_covariance`).
  For an MML fit, the structural parameter covariance is computed as the
  inverse of the observed Fisher information, evaluated as the numerical
  Jacobian of the analytical gradient returned by
  `mfrm_loglik_mml_value_grad` (O(P) gradient calls, not O(P^2) function
  calls). Inversion uses an eigendecomposition with regularization for
  near-singular spectra; a `regularized` flag and a status field
  (`"ok"`, `"regularized"`, `"not_applicable"`, `"fallback"`) propagate
  to downstream consumers so they can label the resulting SE as
  regularization-aware.
- **Structural delta-method SE / CI for the GPCM Linacre fair-average**
  (`add_gpcm_fair_average_delta_se`, `fair_average_table(fair_se=True)`).
  For each non-Person row, the routine computes
  `SE = sqrt(grad^T Cov grad)` where `Cov` is the MML observed-information
  covariance and `grad` is the finite-difference gradient of the
  fair-average scalar function with respect to the optimisation parameter
  vector. Confidence intervals are reported at the requested `ci_level`
  (default 0.95) and clipped to the rating bounds. Person rows are marked
  `"not available"` because MML EAP person estimates are not part of the
  structural Hessian; non-GPCM fits are marked `"not_applicable"`.
  Reference: standard large-sample delta method (e.g. Cramer 1946);
  observed-information identity for MML estimators (Louis 1982).
- **FACETS-style tables tab now exposes the fair-average SE / CI inline.**
  A new toggle at the top of the tab lets the user request fair-average
  standard errors and confidence intervals; a confidence-level selector
  (90 % / 95 % / 99 %) sits next to it. The toggle is disabled for non
  GPCM-MML fits with a one-line explanation. When enabled, a spinner runs
  the Hessian + delta-method evaluation once per fit, the result is cached
  in `st.session_state` keyed on a content hash of the optimisation
  parameter vector so subsequent tab interactions are instantaneous, and
  the table columns are reordered so each `Fair(M)` / `Fair(Z)` value
  sits adjacent to its S.E., CI bounds, and status flag. Method and
  citation captions are rendered below the controls when the toggle is
  on. Strings are localised in both English and Japanese.

### Changed

- **Onboarding banner is now compressed and foldable.** The 3-step
  quickstart guide moved into a Streamlit expander (default expanded for
  first-time visitors, collapsible thereafter) so the Run / Dismiss
  buttons stay visible without paying a large vertical footprint for the
  explanatory text. The redundant tip-caption row (which advertised
  past v0.2.x features and pointed at the Tutorial expander that already
  sits directly below) was dropped; the Tutorial pointer is now folded
  into the bottom of the banner body.
- **Decorative emojis removed from sample-data scenario labels and from
  the onboarding banner.** Functional status indicators (red / amber /
  green readiness colours, Wright & Linacre fit-band swatches) are
  retained because they pair text labels with a redundant colour signal
  for accessibility (WCAG 1.4.1).
- **Default sample-data scenario (`writing_essay`) description expanded.**
  Added a "why it is the default" paragraph, a "what to look at first"
  walkthrough of the result tabs in their natural reading order, and a
  "learning points specific to this data set" section that notes the
  modest rater-severity spread and the PCA sample-size caveat at n = 30.
- **`mfrmr_020_migration_coverage_table()`** documents the migration
  surface for the mfrmr 0.2.0 R package release. It inherits every
  0.1.5 / 0.1.6 row (renaming the area column), overrides the bounded-
  GPCM row to reflect this release's slope-aware fair-average and
  structural delta-method SE delivery, and appends eleven
  0.2.0-specific rows: the slope-aware fair-average and SE pair
  (Ready), the MML observed-information covariance (Ready), and nine
  Planned rows tracking the slope-aware bias inference, the FACETS
  df / ZSTD reporting layer, the network analysis triple, the G/D-
  study planner, ADEMP parameter recovery, Snijders polytomous lz*
  correction, category-specific information curves, the model-choice
  review helper, and the kable / flextable / monochrome APA presets.
  The new table is surfaced in the Help tab's coverage panel, in the
  per-run download bundle, in the disk export, and in the
  external-validation report template; the public release-readiness
  check now reports against the 0.2.0 table.
- **Category-specific information curve decomposition.** A new helper
  `build_category_information_curve_data()` returns the per-category
  contribution to the total Fisher information at every theta on the
  curve grid: `I_k(theta) = a^2 * P_k(theta) * (k - E[X | theta])^2`
  (Muraki, 1993, Eq. 16). The contributions sum to the existing
  total information curve by construction; the identity is pinned at
  machine precision across a synthetic three-level GPCM fixture. A
  new optional expander section in **Visuals -> Information
  Curves** renders the per-category breakdown as one line per
  category using a Wong (2011) colour-blind-safe palette, with the
  total information curve overlaid as a dashed line so the
  decomposition-equals-total identity is visible at a glance. A
  per-category summary table reports the peak information, the theta
  at the peak, the integrated information under the curve, and the
  category's share of the integrated total. The data is also
  downloadable as CSV. Reference: Muraki (1993) Applied Psychological
  Measurement 17(4), Eq. 16.
- **Slope-aware GPCM bias inference: information identity, likelihood-
  ratio test, and profile-likelihood confidence interval.** Three
  related fixes / additions to `estimate_bias_interaction()`:
  * **Information identity.** Under GPCM the conditional Fisher
    information for the additive bias shift `b` is
    `I(b) = sum_i a_i^2 * Var[X_i | eta_i + b] * w_i`, where `a_i` is
    the per-observation discrimination. The previous implementation
    used `I(b) = sum_i Var[X_i | eta_i + b] * w_i`, which is the
    RSM / PCM (a = 1) form; the resulting GPCM bias `S.E.` column was
    systematically too large because the `a^2` factor was missing. The
    fix is gated on `model == "GPCM"`, so RSM and PCM fits are
    byte-identical to the previous release. Reference: Muraki (1993)
    Applied Psychological Measurement 17(4), Eqs. 7 and 16.
  * **Likelihood-ratio test.** Each GPCM bias cell now reports `LR
    ChiSq`, `LR d.f.` (always 1), and `LR Prob.` from the chi-square
    pivotal `Lambda = 2 * (loglik(bias_hat) - loglik(0)) ~ chi2_1`
    (Wilks, 1938). RSM and PCM cells leave these columns as NaN with
    an explanatory `Likelihood Basis` string instead of silently
    reusing the t-based screening tier.
  * **Profile-likelihood confidence interval.** The new
    `_profile_bias_ci` helper inverts the same chi-square pivotal to
    produce a `(1 - alpha)` confidence interval for the additive bias
    shift, by solving `2 * (NLL(b) - NLL(bias_hat)) = chi2_{1, 1-alpha}`
    on both sides of `bias_hat` via Brent root-finding inside
    `[-max_abs, max_abs]`. When the likelihood never falls far enough
    inside the bracket the endpoint is returned with `Profile CI
    Status = "limited by search range"` so the caller can render the
    truncation honestly. Status values are `"ok"`, `"limited by
    search range"`, or `"not available"`. Reference: Cox (1975)
    Biometrika 62(2).
  * **UI.** The bias / interaction tab inserts `LR ChiSq`, `LR Prob.`,
    `Profile CI Lower`, `Profile CI Upper`, and `Profile CI Status`
    columns immediately after the existing t-based block when at
    least one cell carries a finite `LR ChiSq` (i.e. for GPCM fits).
    A short caption below the table cites Wilks (1938) and Cox (1975)
    and reminds readers that theta, step thresholds, slopes, and
    other facet estimates are held fixed inside the conditional
    profile. RSM / PCM fits show the unavailable caption.
  * **mfrmr 0.2.0 coverage.** `mfrmr_020_migration_coverage_table()`
    has the `GPCM bias inference - slope-aware` row and the
    `Category-specific information curves` row both updated from
    `Planned` to `Ready`. The new `PythonEvidence` cites the
    slope-aware Fisher information identity (Muraki, 1993, Eqs. 7,
    16), the LR chi-square pivotal (Wilks, 1938), the
    profile-likelihood CI inversion (Cox, 1975), and the publication-
    document integration; the `Boundary` paragraph documents the
    still-uniformly-`screening` inference tier and the conditional
    profile semantics (theta, step thresholds, slopes, and other
    facet measures held fixed). The Bounded GPCM row's `Boundary` is
    also updated to remove the obsolete "bias inference still uses
    the non-slope-aware identity pending a follow-up" caveat; the
    paragraph now points directly at the dedicated 0.2.0 bias row
    for the slope-aware machinery. The
    `external_validation_report_template` claim row for
    `mfrmr 0.2.0 migration` is updated in lock-step so the
    public-wording paragraph names slope-aware GPCM bias inference
    among the shipped pieces rather than among the roadmap items.
  * **Publication-document integration.** The Word, PDF, and HTML
    manuscript exports now embed a bias / interaction table after the
    element-measures and reliability tables. For GPCM fits the
    table includes the new ``LR ChiSq``, ``LR Prob.``, ``Profile CI
    Lower``, ``Profile CI Upper``, and ``Profile CI Status`` columns;
    RSM and PCM fits show the t-based screening columns only. Two
    new helpers ``_bias_table_for_publication`` and
    ``_bias_table_caption`` are shared across the three builders so
    the column selection and the caption stay in lock-step.
    Captions cite Muraki (1993) for the slope-aware information
    identity, Wilks (1938) for the LR chi-square pivotal, and
    Cox (1975) for the profile-likelihood confidence interval; the
    APA reference library and the citation map gain entries for
    Wilks (1938), Cox (1975), Louis (1982), Cramer (1946), and
    Muraki (1993) so a Methods narrative that mentions these
    citations in ``(Author, Year)`` form auto-resolves to a
    properly alphabetised references list.
  * **R parity.** A new test
    (`test_bias_estimation_matches_r_reference_within_tolerance`)
    runs the same GPCM MML fit and bias estimation in Python and in
    the mfrmr 0.2.0 R reference (R 4.5.2) against a shared
    deterministic CSV (`tests/data/r_bias_parity_input.csv`,
    180 rows). The R-side output (`tests/data/r_bias_parity_output.json`,
    six per-cell reference values) is generated once via `Rscript` +
    `estimate_bias()` and committed to the repository. The parity test
    verifies that Python's per-cell `Bias Size`, `S.E.`, `LR ChiSq`,
    `LR Prob.`, and both profile-likelihood CI endpoints agree with
    the R reference within tolerances calibrated to the observed
    implementation-vs-implementation convergence variance (the two
    independent MML fits converge to log-likelihoods within `~0.02`
    of each other, and the downstream cells inherit that variance);
    all observed disagreements stay well below the `0.05`-logit
    "negligible" threshold cited in Linacre's FACETS Manual. The
    categorical fields (`Profile CI Status`, `InferenceTier`) must
    agree exactly.
- **Wide-format Excel / spreadsheet upload with automatic melt to
  long format.** Excel users typically have one row per
  ``(Person, Rater)`` and one column per item being scored (e.g. C1,
  C2, C3 for criteria); the MFRM likelihood pipeline expects long
  format (one row per observation), so the upload sidebar now offers
  an inline melt step. A new ``detect_wide_format_columns(df)``
  heuristic classifies the parsed table as long or wide by counting
  numeric vs non-numeric columns: 3+ numeric columns plus 1+ id
  column flips the auto-detect to "wide" and the layout expander
  defaults to open. A new ``apply_wide_to_long_pivot(df, id_cols,
  score_cols, new_facet_name, score_col_name)`` helper performs the
  melt, dropping blank / non-numeric score cells so the likelihood
  never sees zeros for missing observations. The sidebar UI surfaces
  ID-columns / score-columns multiselects, a new-facet name input
  (default ``Criterion``), and a long-format preview after pivot. Math
  contract pinned in ``tests/test_wide_to_long_pivot.py`` (11 tests:
  detection on canonical wide / long fixtures, pivot row-count
  identity, id-column preservation, blank cell handling, refusal on
  overlapping / missing / name-collision inputs, end-to-end fit
  through ``mfrm_estimate`` on pivoted data).
- **Self-tests for the seven 0.2.0 helpers** added in this session
  (FACETS d.f. / ZSTD alignment, model-choice AICc + Akaike weights,
  G/D-study Brennan 3-way identity, design network topology, rater
  severity higher-prop convention, APA preset Markdown / HTML
  round-trip, DIMTEST nonparametric smoke). The module-level
  ``_self_test_*`` functions run at import time via
  ``run_self_tests()`` (now 55 tests, was 48) and complement the
  pytest suite by exercising the helpers on a tiny shared fixture
  ``_make_self_test_3way_rating_data()``. Each test pins a small
  closed-form contract: the Welch-Satterthwaite denominator
  identity, the Hurvich-Tsai AICc closed form and weight-sum-to-one
  identity, Brennan Eq. 3.18 G computation, node-count parity in
  the design graph, the Rater1HigherProp + Rater2HigherProp +
  TieRate = 1 identity, HTML escaping in the APA preset, and the
  Z / p-value finiteness on a small DIMTEST run.
- **R parity for the Brennan 3-way G-study decomposition** against
  lme4 REML with explicit two-way interaction terms. The R fixture
  generator now fits

      Score ~ 1 + (1 | Person) + (1 | Rater) + (1 | Criterion) +
              (1 | Person:Rater) + (1 | Person:Criterion) +
              (1 | Rater:Criterion)

  via lme4 (the three-way interaction is intentionally omitted
  because it is confounded with residual under one observation per
  cell) and dumps the variance components plus the closed-form
  G / Phi to ``tests/data/r_helper_parity_output.json``. Four new
  parity tests verify the agreement: the variance-component source
  list matches across implementations; Person is a dominant
  variance source in both; the observed (n_p, n_r, n_c) sample
  sizes round-trip; and the headline G / Phi coefficients agree to
  ~7e-2 absolute. Individual variance components can disagree at the
  zero boundary because Henderson III MoM (Python) clamps a
  negative-MoM estimate at zero and leaves the residual unchanged,
  while lme4 REML re-solves the remaining variances under the
  non-negativity constraint (Brennan 2001 p. 81; Searle, Casella &
  McCulloch 1992 Ch. 4.6) — both are valid but allocate boundary
  variance differently, which is documented in the test docstrings.
- **D-study CI bands at every planned design point**. The cluster
  bootstrap for G / Phi now propagates through each
  ``(n_facet_a, n_facet_b)`` row of the D-study forecast grid by
  reusing the same B variance-component replicates: each replicate
  re-evaluates Brennan (2001) Eq. 3.18 at every grid point so the
  marginal cost per band is one closed-form arithmetic step rather
  than another bootstrap loop. With ``bootstrap_ci=True`` the
  ``d_study`` DataFrame gains ``G_lower``, ``G_upper``,
  ``Phi_lower``, ``Phi_upper`` columns; manuscripts can now quote
  "G = 0.85, 95% CI [0.79, 0.91] at n_rater = 4, n_criterion = 6"
  directly from a single table. The observed-design row of the
  grid matches the observed-design bootstrap CI to machine
  precision (same replicate set). The Report-tab UI shows the
  bands inline; a new caption explains the cost-sharing structure.
  Math contract pinned in ``tests/test_g_d_study.py`` (CI columns
  present only when bootstrap is enabled, ordered bounds inside
  [0, 1], observed-design row matches the observed-design CI, band
  width shrinks with design size on >= 70 % of comparable pairs,
  seed determinism).
- **Nonparametric DIMTEST for essential unidimensionality** (Stout,
  1987, Psychometrika 52; polytomous adaptation per Nandakumar & Yu,
  1996). A new ``compute_dimtest_nonparametric(res, diagnostics)``
  helper tests whether the criteria measure a single latent
  dimension after the fitted MFRM main effects are accounted for.
  The test averages the rater scores per (Person, Item), partitions
  the items into an Assessment subtest (AT) and a Partitioning
  subtest (PT), stratifies persons by their PT total score, and
  computes the pooled within-stratum covariance among AT items as
  the test statistic T_L. Under essential unidimensionality T_L is
  asymptotically zero; departures signal a residual secondary
  dimension. The standard error and p-value come from a cluster
  bootstrap on persons (Efron & Tibshirani, 1993, Ch. 19) rather
  than the asymptotic null distribution; this keeps the pipeline
  purely nonparametric and side-steps the parametric bias-
  correction term T_B from Stout (1987) that would require fitting
  an auxiliary unidimensional IRT model on the AT subtest.

  AT / PT subtest selection follows Roussos & Stout (1996): the
  default ``pc2_loading_sign`` auto-split uses the PC2 loadings on
  the item facet from the existing residual-PCA pipeline. The
  auto-split inherits a known selection-bias inflation of Type I
  error (PC2 already maximises an approximately-orthogonal
  covariance direction); for a rigorous confirmatory test, pass
  ``at_items=`` and ``pt_items=`` derived from substantive
  considerations or from a separate calibration sample. The
  ``split_method`` field on the returned bundle names the active
  path.

  The Dimensionality tab gains a DIMTEST panel below the residual-
  PCA tabs with controls for the item facet, the AT / PT split
  method (auto vs user-specified multi-select), bootstrap
  replicates, and confidence level. The result reports T_L, Z,
  p-value, person count, the bootstrap CI, the AT / PT split, and a
  Reject (p < 0.05) / Retain status banner. Math contract pinned in
  ``tests/test_dimtest_nonparametric.py`` (10 tests: refusal on
  invalid input, overlap detection, unidim-fails-to-reject on a-
  priori split, 2-dim-rejects on correct split, T_L magnitude
  ordering, bootstrap determinism, bundle field completeness, the
  selection-bias caveat appears on the auto-split path, and the
  closed-form two-sided p-value identity). Three new APA references
  (Stout, 1987; Nandakumar & Yu, 1996; Roussos & Stout, 1996) are
  added to the library and citation map.
- **Bootstrap CI for G / Phi coefficients** (Efron & Tibshirani, 1993,
  Ch. 19). ``compute_generalizability_study()`` now accepts
  ``bootstrap_ci=True`` and runs a cluster-bootstrap resampling on
  persons (the object of measurement): each replicate resamples
  persons with replacement, re-runs the variance-component
  decomposition on the bootstrapped sample, and recomputes G / Phi.
  The bundle's new ``bootstrap_ci`` field reports
  ``G_lower / G_upper``, ``Phi_lower / Phi_upper``,
  ``G_replicates / Phi_replicates`` (the full replicate arrays for
  downstream diagnostics), and the active settings
  (``confidence``, ``n_bootstrap``, ``n_success``, ``method =
  "cluster_bootstrap_on_persons"``). Person-level cluster
  resampling is the standard nonparametric bootstrap for G-theory
  because it respects the within-person correlation structure that
  drives both relative and absolute error. The Report tab gains a
  Bootstrap CI expander with a checkbox, replicate count, and
  confidence-level dropdown; when enabled the panel reports the
  CI in a single caption line below the observed G / Phi metrics.
  Math contract pinned in ``tests/test_g_d_study.py`` (default-off
  behaviour, finite-ordered bands when enabled, CI widening with
  confidence level, seed determinism, [0, 1] bounds on the
  replicate arrays). One new APA reference (Efron & Tibshirani,
  1993).
- **G/D-study: full Brennan (2001) 3-way decomposition with explicit
  two-way interaction terms** for the canonical balanced p x i x j
  design with one observation per cell.
  ``compute_generalizability_study()`` now detects whether the
  design admits the full random-effects ANOVA (object_facet plus
  exactly two random facets, balanced, one observation per cell)
  and dispatches to a new
  ``_crossed_anova_two_way_variance_components`` helper that returns
  seven variance components: the three main effects, the three
  two-way interactions (``object:facetA``, ``object:facetB``,
  ``facetA:facetB``), and a Residual term that is the three-way
  interaction confounded with error (the standard one-observation-
  per-cell confounding). The method-of-moments estimators follow
  Brennan (2001, Table 3.5):

      sigma2_p   = (MS_p   - MS_pi - MS_pj + MS_pij) / (n_i n_j)
      sigma2_i   = (MS_i   - MS_pi - MS_ij + MS_pij) / (n_p n_j)
      sigma2_j   = (MS_j   - MS_pj - MS_ij + MS_pij) / (n_p n_i)
      sigma2_pi  = (MS_pi  - MS_pij) / n_j
      sigma2_pj  = (MS_pj  - MS_pij) / n_i
      sigma2_ij  = (MS_ij  - MS_pij) / n_p
      sigma2_pij = MS_pij

  Negative MoM estimates are clamped at zero per the Henderson III
  convention (Brennan 2001 p. 81). The G / Phi formulas use Brennan
  (2001) Eq. 3.18 / 3.19 with explicit interaction variances:

      sigma2(delta) = sigma2_pi / n_i + sigma2_pj / n_j +
                      sigma2_pij / (n_i n_j)
      sigma2(Delta) = sigma2(delta) + sigma2_i / n_i +
                      sigma2_j / n_j + sigma2_ij / (n_i n_j)
      G   = sigma2_p / (sigma2_p + sigma2(delta))
      Phi = sigma2_p / (sigma2_p + sigma2(Delta))

  The bundle's ``observed_coefficients.details.decomposition``
  field names the active mode (``"full_3way_brennan_eq_3_18"``
  vs ``"main_effects_only_approximation"``) so downstream
  consumers and manuscripts can name the formula they used. For
  designs with one random facet, three or more random facets, or
  unbalanced cells, the helper falls back to the existing main-
  effects-only ANOVA approximation; the existing caveat is
  preserved for that path. Math contract pinned in
  ``tests/test_g_d_study.py``: the seven-source variance-component
  layout, Brennan Eq. 3.18 closed-form G / Phi identity, and the
  decomposition label. ``mfrmr_020_migration_coverage_table()``
  boundary text updated so manuscripts no longer see "main-effects-
  only" as the default characterisation; the new boundary names
  the dispatch logic explicitly.
- **Akaike / Schwarz / AICc weights in model-choice guidance** (Burnham
  & Anderson, 2002, Eq. 2.8). Each model-choice comparison row now
  also carries an ``AICWeight``, ``AICcWeight``, and ``BICWeight``
  column. The weights are computed as ``w_i = exp(-DeltaIC_i / 2) /
  sum_j exp(-DeltaIC_j / 2)`` and admit a probability interpretation:
  ``w_i`` is the relative likelihood that model ``i`` is the
  Kullback-Leibler best of the candidate set (Burnham & Anderson,
  2002, Section 2.9). The weights sum to 1 over the finite-DeltaIC
  rows by construction; manuscripts that quote a "best model
  probability" can read it directly off the table. Math contract
  pinned in ``tests/test_model_choice_guidance.py`` (sum-to-one,
  closed-form match, best-row carries the largest weight).
- **AICc finite-sample correction in model-choice guidance** (Hurvich &
  Tsai, 1989, Eq. 1). Each model-choice comparison row now carries
  ``AICc = AIC + 2 K (K + 1) / (N - K - 1)`` alongside the existing AIC
  and BIC, plus the corresponding ``DeltaAICc`` and
  ``AICcEvidenceRatio = exp(DeltaAICc / 2)`` columns. A new
  ``AICcRecommended`` flag fires when ``N / K < 40`` (Burnham &
  Anderson, 2002, p. 66), and the tiered recommendation logic switches
  the secondary criterion from AIC to AICc whenever any candidate
  satisfies the recommendation window — small-design users no longer
  get an over-confident AIC-based ranking. The N / K ratio is reported
  per row so the applicability is recoverable from the table itself.
  Math contract pinned in
  ``tests/test_model_choice_guidance.py`` (closed-form AICc identity,
  applicability flag, ``AICc >= AIC`` bound, DeltaAICc minimum-is-zero,
  caveat-mentions-AICc). One new APA reference (Hurvich & Tsai, 1989).
- **R parity fixtures for four 0.2.0 helpers** (FACETS d.f. / ZSTD
  alignment, rater severity directed network, rater halo
  correlation network, design network connectivity). The R-side
  reference is regenerated by ``tests/data/generate_r_helper_parity_fixture.R``
  against the same deterministic CSV the existing bias-inference
  parity test consumes; the JSON output
  (``tests/data/r_helper_parity_output.json``) carries the per-facet
  d.f. / ZSTD columns from ``dx$fit`` under ``fit_df_method = "both"``,
  the directed rater pair metrics from
  ``rater_network_analysis(mode = "severity_direction")``, the
  Spearman correlations from ``rater_halo_network_analysis(method =
  "spearman")``, and the graph-level summary + per-node metrics from
  ``mfrm_network_analysis()``. Both ecosystems run JML with
  ``maxit = 200`` and ``reltol = 1e-6`` so the
  implementation-vs-implementation drift stays inside the calibrated
  0.02-log-likelihood band the existing bias-parity test
  characterises. Five new parity tests pin the agreement:
  exact counts on N / higher-counts / articulation flags, ~5e-2
  absolute tolerance on d.f., ~1.5e-1 on ZSTD, machine precision on
  the rater severity proportions, ~5e-6 on Spearman rho, and
  byte-identical topology on the design graph.

  Two **convention alignments to mfrmr** went with the R parity
  fixtures so the numerical comparison is apples-to-apples:

  * ``zstd_from_mnsq`` now guards ``df >= 1`` (matching the mfrmr
    engine convention) instead of ``df > 0``. Cells with d.f. < 1
    now report NaN ZSTD rather than a sign-flipped value driven by
    the ``1 - 2/(9 df)`` centring term. The FACETS-convention
    helper ``zstd_from_mnsq_facets`` retains ``df > 0`` with the
    explicit cap (the Wright-Masters Welch-Satterthwaite d.f. is
    fractional by design).
  * ``_rater_severity_pair_metrics`` now reports
    ``Rater1HigherProp = Rater1HigherCount / N`` (the mfrmr /
    FACETS / Eckes (2011) "share of all shared contexts" convention)
    as the primary proportion. The previous "share of disagreements"
    conditional probability remains available as
    ``Rater1HigherCondProp`` for callers that want the
    sum-to-one normalisation on the disagreement subset, and the
    new ``TieCount`` column carries the agreement count so
    ``HigherCount1 + HigherCount2 + TieCount = N`` holds by
    construction. The R parity fixture's
    ``Rater1HigherProp`` line up with the Python column at
    machine precision after the alignment.

- **Rater severity / leniency directed network and rater halo
  cross-criterion screen.** Two further network screens
  complete the mfrmr 0.2.0 network-analysis triple alongside the
  earlier design network. ``compute_rater_severity_network(res)``
  builds the directed rater-leniency graph: for every pair of raters
  with shared scoring contexts, the edge ``A -> B`` carries weight
  ``P(A > B | A != B)`` from the directional higher-count diagnostic,
  and the per-rater ``LeniencyIndex = out_mass - in_mass`` plus
  ``SeverityRank`` give a Rasch-style severity ranking that recovers
  the data-generating order exactly on the synthetic fixture (R1
  generated with the smallest ``rater_eff`` ends up at rank 1).
  ``compute_rater_halo_network(res)`` builds the
  ``(Rater, Criterion)`` node graph with Spearman edge weights and
  Bonferroni-adjusted retention; per-rater ``ReviewStatus`` of
  ``warning`` / ``review`` / ``ok`` is keyed on the configurable
  ``halo_weight_review`` and ``halo_contrast_review`` thresholds. The
  warning tier upgrades when the halo signal dominates the retained
  edges (high same-rater cross-criterion mean correlation AND either
  a clear halo / non-halo contrast or no non-halo edges retained at
  all). The Report tab gains two new panels — Rater severity and
  Rater halo — each with summary, per-rater table, expander for the
  full pair / edge metrics, and per-table CSV download. Math contract
  pinned in ``tests/test_rater_networks.py`` (16 tests covering
  refusal on invalid input, the directional-count identity
  ``Rater1HigherProp + Rater2HigherProp = 1``, severity-rank recovery
  on synthetic data, ``sum(LeniencyIndex) = 0`` under a fully crossed
  design, halo refusal on coinciding facets, warning detection on
  halo-patterned data, ``ok`` retention on clean RSM data, the
  contrast identity ``HaloNonHaloContrast = HaloMeanWeight -
  NonHaloMeanWeight``, settings round-trip, and node-count parity).
  Two new APA references (Lai, Wolfe, & Vickers, 2015; Lamprianou,
  2025) are added. The
  ``mfrmr_020_migration_coverage_table()`` row for "Network analysis
  (mfrm / rater / rater-halo)" flips from "Ready (design network);
  rater / halo screens remain follow-up work" to plain "Ready". The
  Welch test contrasting halo vs non-halo edge weights is reported
  but described as descriptive only (edge weights are clustered by
  rater and node), so ``ReviewStatus = 'warning'`` should be read as
  a follow-up trigger rather than a causal halo diagnosis.
- **Design network analysis (connectivity, articulation points, bridges,
  betweenness centrality).** A new ``compute_design_network_analysis(res)``
  helper treats the rating design as an undirected weighted graph
  (nodes = ``(Facet, Level)`` pairs, edges = co-observed in at least
  one rating row, weight = co-observation count) and computes the
  canonical connectivity diagnostics from networkx (Hagberg, Schult,
  & Swart, 2008): graph-level summary (Nodes, Edges, Components,
  LargestComponentNodes, LargestComponentShare, Density, MeanDegree,
  MeanStrength, ArticulationPoints, Bridges, Connected, Diameter,
  MeanDistance), per-node metrics (Degree, Strength, Betweenness,
  Closeness, EigenvectorCentrality) with an ``IsArticulationPoint``
  flag, per-edge bridge flag, and per-facet aggregates. The Report
  tab gains a Design Network section with the summary table, the
  per-facet aggregate, articulation-point and bridge-edge callouts
  (green-success messages when none are found), a per-node metrics
  expander, and a CSV download for the combined summary + node +
  edge + facet tables. ``min_observations`` lets the user drop weak
  edges; ``include_graph=True`` round-trips the underlying
  ``networkx.Graph`` for downstream visualisation. Math contract
  pinned in ``tests/test_design_network.py`` (14 tests covering
  refusal on invalid input, node-count identity, fully-crossed
  designs having 1 component and no cut-points, disjoint designs
  reporting Components >= 2, facet-summary aggregation parity, cut-
  node / bridge-edge subset invariants, ``min_observations`` filter,
  ``include_graph`` round-trip, and centrality [0, 1] bounds). Two
  new APA references (Csardi & Nepusz, 2006; Hagberg, Schult, &
  Swart, 2008) are added. The
  ``mfrmr_020_migration_coverage_table()`` row for "Network analysis
  (mfrm / rater / rater-halo)" flips from Planned to
  "Ready (design network); rater / halo screens remain follow-up
  work" — the directed severity / halo screens require per-rater
  pairwise comparison machinery not yet in scope.
- **APA 7 table presets and monochrome figure palette.** Two new
  helpers, ``apa_table_to_markdown(df, ...)`` and
  ``apa_table_to_html(df, ...)``, re-emit any report DataFrame as a
  manuscript-ready APA 7 snippet (bold table number, italic caption,
  header / body grid with top + header + bottom rules, ``Note.`` block
  below). Float digits are user-configurable; boolean columns render
  as ``Yes`` / ``No``; missing values render as the empty string; HTML
  output is HTML-escaped end-to-end (cells, captions, notes) so
  ``<script>`` payloads cannot break out of the cell. A companion
  monochrome palette is exposed via ``get_monochrome_palette()`` (eight
  grayscale steps from near-black to very light grey) and
  ``get_monochrome_dash_patterns()`` (six Plotly dash patterns) so
  figure exports under a print-friendly preset stay distinguishable
  even at low DPI. The Report tab gains an APA-presets panel with a
  table picker, format radio (Markdown / HTML), digits / table-number
  / caption / note controls, an inline preview, and per-format
  download buttons; the monochrome palette is previewed with colour
  swatches and dash-pattern names. Math contract pinned in
  ``tests/test_apa_presets.py`` (18 tests). The
  ``mfrmr_020_migration_coverage_table()`` row for "APA table presets
  (kable / flextable / monochrome)" flips from Planned to Ready. The
  helpers complement the existing single-file Publication Document
  rather than replace it: use the APA preset for per-table snippets
  authors paste into a longer manuscript, and the Publication
  Document for the bundled Word / PDF / HTML deliverable.
- **Generalizability theory (G-study) and D-study forecast.** A new
  ``compute_generalizability_study(res)`` helper performs a method-of-
  moments crossed-ANOVA decomposition on the rating data and reports
  the canonical generalizability-theory outputs (Cronbach, Gleser,
  Nanda, & Rajaratnam, 1972; Brennan, 2001): per-source variance
  components with proportion-of-variance shares, G (relative decision)
  and Phi (absolute decision) coefficients at the observed design,
  and a D-study forecast grid that scales G / Phi over planned
  numbers of raters and criteria. The dispatch follows the standard
  closed form: ``G = sigma2_p / (sigma2_p + sigma2_e / n_total)`` and
  ``Phi = sigma2_p / (sigma2_p + sum(sigma2_facet / n_facet) +
  sigma2_e / n_total)``. The Report tab gains a G-study / D-study
  section with the variance-component table, an observed-design
  G / Phi metric strip, the D-study forecast table, and a CSV
  download; the panel runs every page load because the math is fast
  (a few ANOVA reductions on the rating data). Math contract pinned
  in ``tests/test_g_d_study.py`` (15 tests: refusal on invalid input,
  variance-component non-negativity, required source list, proportion-
  of-variance sum = 1, closed-form G / Phi identity, ``Phi <= G``,
  D-study monotonicity in each facet sample size, D-study grid
  contains the observed design, row count = product of facet grids,
  bounds in [0, 1], and Brennan / Cronbach citation hygiene). Two
  new APA references (Brennan, 2001; Cronbach, Gleser, Nanda, &
  Rajaratnam, 1972) are added. ``mfrmr_020_migration_coverage_table()``
  flips the "G/D-study planning (mfrm_d_study)" row from Planned to
  Ready. The implementation deliberately uses the one-observation-
  per-cell main-effects-only approximation (no explicit person-by-
  facet interaction terms), which is the standard simplification when
  ``lme4`` random-effects fitting is not available; the caveat
  surfaces in the bundle's ``caveat`` field for manuscript citation.
- **Parameter recovery (ADEMP) simulation.** A new
  ``evaluate_parameter_recovery(...)`` helper runs a Monte Carlo
  parameter-recovery study under one explicit data-generating setup
  and reports the performance measures of Morris, White, and Crowther
  (2019): bias, MCSE(bias), RMSE, MCSE(RMSE), MAE, raw / aligned
  errors, Pearson correlation, mean SE, SE-availability rate, and
  95 % coverage. The data-generating mechanism samples person
  measures, rater severities, and criterion difficulties from
  mean-zero independent normals; step thresholds are equally spaced
  over a chosen logit span (RSM uses shared thresholds; PCM and GPCM
  use per-criterion perturbed thresholds). Each replicate refits the
  requested model via ``mfrm_estimate``, runs a lightweight
  ``mfrm_diagnostics`` pass to populate SE on the measures table,
  mean-aligns the location-indeterminate blocks (Person, Rater,
  Criterion) against the truth, and computes per-row error +
  coverage95. The bundle carries six fields: ``recovery``
  (long table), ``recovery_summary`` (per-parameter-type aggregates),
  ``coverage_report`` (SE/CI coverage diagnostic),
  ``rep_overview`` (per-rep convergence and timing), ``ademp``
  (Aims / Data-generating mechanism / Estimands / Methods /
  Performance measures narrative), and ``settings`` (round-trip of
  the input arguments). The Report tab gains a Parameter-recovery
  (ADEMP) section with eight controls (model, fit_method,
  replicates, seed, n_person, n_rater, n_criterion, n_cat) and
  caches the bundle in session_state keyed on the input settings.
  Math contract pinned in ``tests/test_parameter_recovery.py``
  (18 tests covering refusal on bad inputs, three location-block
  ParameterType values, long-table row counts, the mean-alignment
  identity at machine precision, Bias = mean(ErrorAligned), RMSE /
  MAE closed forms, Coverage95 closed form, SE/CI status propagation,
  explicit coverage-report construction, fixed-seed determinism,
  ADEMP narrative completeness, and correlation positivity on a
  clean RSM fit). One new APA reference (Morris, White, & Crowther,
  2019) is added to the library and citation map. The
  ``mfrmr_020_migration_coverage_table()`` row for "Parameter
  recovery (ADEMP)" flips from Planned to Ready.
- **Model-choice guidance: RSM / PCM / GPCM side-by-side comparison.**
  A new ``compute_model_choice_comparison(res)`` helper refits the two
  non-current rating-scale models on the same data and returns a
  publication-style comparison bundle. Each row carries ``Model``,
  ``FitStatus``, ``N``, ``KParams``, ``LogLik``, ``AIC`` (Akaike,
  1974), ``DeltaAIC``, ``AICEvidenceRatio``, ``BIC`` (Schwarz, 1978),
  ``DeltaBIC``, and ``BICEvidenceRatio``; the evidence ratios use
  the Akaike / Schwarz closed form ``exp(DeltaIC / 2)`` (Burnham &
  Anderson, 2002, Eq. 2.10). A nested likelihood-ratio table reports
  ``Lambda = 2 * (LL_alt - LL_null) ~ chi2_{df_alt - df_null}`` for
  every pair where both fits succeeded (RSM in PCM, PCM in GPCM, RSM
  in GPCM; Wilks, 1938). A tiered recommendation (``strong`` /
  ``moderate`` / ``weak`` / ``tie``) is keyed on the BIC gap to the
  second-best candidate; the rationale string names the per-criterion
  margin so manuscripts can quote it directly. The Report tab gains
  a "Model-choice guidance" section that exposes the comparison behind
  a "Run model-choice comparison" button; the result is cached in
  session_state keyed on a content-addressed hash of the result so
  re-loading the same fit skips the refit. Runs with anchors,
  population / latent-regression terms, or facet regularization are
  refused with an explanatory reason because their likelihoods are
  not directly comparable on a common scale. Math contract pinned in
  ``tests/test_model_choice_guidance.py`` (13 tests: invalid-input
  refusal, regularization / population-formula refusal, three-row
  output, DeltaIC = 0 at the minimum, evidence-ratio closed form,
  LR chi-square = 2 (LL_alt - LL_null), scipy parity on p-values,
  preference recovery on RSM-generated and GPCM-heterogeneous-slope
  data, caveat citation, column-order stability, fit_times zero for
  current model). Three new APA references (Akaike, 1974; Schwarz,
  1978; Burnham & Anderson, 2002) are added to the library and
  citation map. The ``mfrmr_020_migration_coverage_table()`` row
  for "Model-choice guidance (RSM / PCM / GPCM)" flips from Planned
  to Ready.
- **FACETS d.f. / ZSTD reporting alignment.** Two d.f. conventions for
  the standardised Infit / Outfit fit statistic now ship side-by-side.
  The engine convention (the default when this alignment first shipped;
  superseded by the Unreleased FACETS-primary change above) uses
  ``DF_Infit = sum(Var * w)`` and ``DF_Outfit = sum(w)`` and reports
  the Wilson-Hilferty (1931) standardised value through
  ``zstd_from_mnsq``. The new FACETS / Wright-Masters convention
  (``Wright & Masters, 1982``, Eqs. 4.20 / 4.27) applies a Welch-
  Satterthwaite d.f. that captures variance heterogeneity through the
  per-observation fourth central moment
  ``M4_i = sum_k P_k * (k - E[X | theta])^4``:
  ``DF_Infit_FACETS = 2 * (sum Var * w)^2 / sum w * (M4 - Var^2)`` and
  ``DF_Outfit_FACETS = 2 * (sum w)^2 / sum w * (M4 / Var^2 - 1)``.
  ``compute_obs_table`` now adds a ``FourthCentralMoment`` column;
  ``calc_overall_fit`` and ``calc_facet_fit`` accept a ``fit_df_method``
  argument with values ``"engine"``, ``"facets"``, or ``"both"``; and
  ``zstd_from_mnsq_facets`` applies the FACETS ``+/- 9`` cap (with the
  cap exposed as a parameter so users can disable it). The Fit Details
  tab has a three-way radio that switches the on-screen fit table
  between conventions without re-fitting, plus a metadata caption that
  names the active method, ZSTD transform, and cap on every render.
  An Import-FACETS-fit-table expander accepts a CSV/TSV with ``Facet /
  Level / InfitZSTD / OutfitZSTD`` columns, joins on ``(Facet, Level)``,
  and renders a per-row delta table (Python minus imported) so users
  can verify the Python output against an external FACETS run. Math
  contract pinned in ``tests/test_facets_df_zstd_alignment.py`` (15
  tests covering the closed-form fourth moment, ``M4 >= Var^2`` bound,
  Wilson-Hilferty cap behaviour, Welch-Satterthwaite formula, and the
  three-way dispatch column shape including the engine-only backward-
  compatibility guarantee). The ``mfrmr_020_migration_coverage_table()``
  row for "FACETS df / ZSTD reporting alignment" flips from Planned to
  Ready. References: Wright & Masters (1982); Wilson-Hilferty cube-root
  transformation as cited in Linacre (2002).
- **Polytomous person-fit indices: lz (Drasgow et al., 1985) and lz*
  (Snijders, 2001).** A new ``compute_person_fit_indices(res)`` reports
  the standardised polytomous log-likelihood for every fitted person.
  Under JMLE the report carries the Snijders (2001, Eq. 16) correction:
  ``c_n = Cov[l, S] / I(theta)``, corrected variance
  ``Var[l] - Cov[l, S]^2 / I(theta)``, and ``lz* = (l - E[l] - c_n * S)
  / sqrt(corrected_var)``. Under MML / EAP the function reports the
  unadjusted ``lz`` with the status field
  ``"not_applicable_eap"`` so the unadjusted form is never silently
  treated as Snijders-corrected. The six per-observation diagnostic
  columns ``PrObserved``, ``ItemEntropy``, ``ItemVarLogP``,
  ``ItemLogPScoreCov``, ``ScoreInformation``, and
  ``ObservedScoreDerivative`` are added to ``compute_obs_table(res)``
  in a single pass over the per-observation probability matrix; the lz
  / lz* aggregation reads them by name. Reporting columns
  (``ReportIndex``, ``ReportValue``, ``ReportFlagLevel``,
  ``ReviewStatus``, ``ReviewReason``, ``ReportCaveat``) choose lz*
  when computed and fall back to lz otherwise, using the standard
  normal thresholds ``|z| > qnorm(0.975)`` (5 %) and
  ``|z| > qnorm(0.995)`` (1 %). The Measures tab surfaces a compact
  panel: a four-metric strip (persons measured, flagged at 5 %, flagged
  at 1 %, chosen report index), a flagged-rows table sorted by
  ``|ReportValue|`` when any person trips the practical threshold, the
  full per-person table behind an expander, a CSV download, and a
  method-aware caveat (``Snijders correction conditional on fitted
  non-person calibration`` vs ``uncorrected lz; interpret tail z-scores
  conservatively``). The lz / lz* indices also ship in the per-run
  downloads bundle as ``person_fit_indices``. Two new APA references
  (Drasgow, Levine & Williams, 1985; Snijders, 2001) are added to the
  library and citation map so a Methods narrative that mentions either
  in ``(Author, Year)`` form auto-resolves to a properly alphabetised
  references list. The math contract is pinned in
  ``tests/test_person_fit_indices.py`` (17 tests covering the closed-
  form lz identity, the closed-form lz* / c_n / corrected-variance
  identity, the threshold rule, status-field fallback for MML and for
  missing-column / degenerate inputs, the small JMLE end-to-end fit,
  and column-order stability). The ``mfrmr_020_migration_coverage_table()``
  row for "Polytomous person-fit correction (Snijders lz*)" flips from
  Planned to Ready with the new evidence narrative.
- **Remaining six sample-data scenarios** (large-scale writing, L2
  speaking, clinical OSCE, writing with missing, music peer-rating,
  reading testlet binary) now follow the same expanded template:
  a one-paragraph design summary, a "what this scenario teaches"
  paragraph, a tab-by-tab "what to look at first" walkthrough, and
  scenario-specific learning points (PCA at n = 120 for the severity-
  outlier scenario, PCM/GPCM step thresholds for the L2 speaking
  criterion contrast, low person separation as a feature for the
  clinical OSCE cohort, the missing-rate and connectivity readiness
  flags for the MAR scenario, the sparse Person × Rater graph for the
  round-robin peer rating, and the local-item-dependence / Stan
  TESTLET_RI hand-off for the binary reading testlet).
- **Decorative emojis removed across the app and i18n.** Help text,
  status messages, button labels, dict-keyed tab labels, and section
  comments no longer carry decorative emoji ornaments (document,
  chart, inbox, paste, microscope, abacus, test-tube, clock, shuffle,
  book, notepad, ruler, floppy, party-popper, alarm-bell, and other
  category badges). Status / severity markers are explicitly retained
  because they pair text labels with a redundant colour signal for
  WCAG 1.4.1 accessibility: red / amber / green readiness lights, the
  Wright & Linacre fit-band swatches, the success / warning / error
  severity icons, the no-severity default circle, and the structural
  yes / no markers used inside capability-matrix tables. Toast
  notifications now use consistent severity icons (`icon="✅"` for
  success, `icon="❌"` for failure) instead of celebratory or alarm-
  bell glyphs.
- Design evaluation and report/download bundles now include the new
  0.1.6-oriented audit tables so Essential-view users see warnings in the UI and
  expert users can export full CSV/Excel evidence.
- Config JSON import contract is pinned at 25 replayable settings so restored
  runs preserve the current analysis-depth and diagnostic controls.
- Public-facing README, Help, release-check, and validation wording now avoid
  machine-specific paths and internal workspace terminology.
- Reproducible script downloads now include a local batch workflow note so
  repeated refits stay outside hosted Streamlit.

### Fixed

- **GPCM Linacre fair-average on non slope-facet rows** (Person, and any
  non-step facet such as Rater or Task under a `step_facet = Criterion`
  fit) is now evaluated at slope = 1, the geometric-mean discrimination
  under the sum-to-zero log-slope identification. The previous code
  substituted the arithmetic mean of the estimated slopes, which by the
  arithmetic-mean / geometric-mean inequality is always >= 1 and so
  reported a systematically over-discriminated fair average there. The
  step-facet's own rows continue to use that level's own discrimination
  and that level's own step thresholds. Reported FairM / FairZ values
  in the FACETS-style report table will change for GPCM fits; RSM and
  PCM fits, and non slope-facet rows under PCM, are unaffected.
  Reference: Muraki (1992), Eqs. 2-3 and 10; Muraki (1993), Eqs. 7
  and 16; Linacre, FACETS Manual section "Fair Average". Cross-checked
  against the mfrmr R reference (0.2.0) at 1e-10 across 27 fixture
  cases.
- `expected_score_from_eta()` now returns `np.nan` for a non-finite or
  non-positive slope instead of silently falling back to slope = 1, so
  a degenerate GPCM fit can no longer masquerade as a healthy one
  through the fair-average path.
- Data source radio now shows the 7 built-in samples plus Paste and Upload
  exactly once; Paste/Upload are no longer repeated after every sample option.

## 0.2.14-beta - 2026-04-18

User feedback after v0.2.13: documentation said "four scenarios" in
several places after a fifth was added; the readiness panel blocked
the binary testlet scenario with a hard [ISSUE] for its single-scorer
facet (working as designed but UX-hostile); the Run history area
needed a FACETS-style one-click bundle download; and two more
scenarios (missing data + peer-rating) were requested. This release
addresses all five points and — once again — the v0.2.12 end-to-end
gate caught two real bugs in CI before users would have hit them.

### Added

- **Two new sample scenarios** (registry now has 7 total):
  - **Writing with missing** (80×4×2×3, ~1,632 obs after 15 %
    MAR drop) — exercises MFRM's tolerance of missing-at-random
    observations, the readiness-panel missing-rate flag, and the
    connectivity check. (Little & Rubin, 2002)
  - **Music peer-rating** (120 musicians × 2 cyclic raters × 2
    pieces × 4 criteria, 1,920 obs) — sparse round-robin design
    where every musician is both an examinee and a rater. Person_i
    rated by Person_(i+1) and Person_(i+2). 240 unique pairs out
    of 14,400 possible. (Myford & Wolfe, 2004; Linacre, 2007)
- **`render_quick_results_download()`** — FACETS-style one-click
  results bundle surfaced right below the Run history panel.
  Offers ZIP (CSVs) and Excel (multi-sheet) variants. Users
  can download the full results bundle without discovering the
  Downloads sub-tab. Built via new helper
  `build_result_bundle_frames()`.
- **`_generate_mfrm_peer_rating_data()`** — sparse cyclic
  generator for round-robin peer-assessment data. The first
  generator in the registry that does NOT assume a fully-crossed
  design.
- **`missing_rate` parameter support** in the scenario registry —
  any scenario can declare a stochastic MAR drop applied after
  generation (separate RNG stream so the response seed is
  unaffected).
- **`tests/test_v0214_scenarios.py`** (12 tests) covering the
  missing-rate drop, the peer-rating generator's structural
  invariants (no self-rating, exactly 2 raters per examinee,
  shared ID pool), the result-bundle helper, and the singleton-
  facet readiness downgrade.

### Changed

- **Readiness check for singleton facets downgraded** from
  [ISSUE] to [CAUTION]. MFRM centers a singleton facet to 0
  without error; a hard block was overkill for legitimate single-
  scorer / single-form designs (e.g., binary reading testlets).
  The detail message now explains the model behaviour.
- **`pick_default_facets` post-filter**: columns with < 2 unique
  levels are excluded from the auto-selected facet set so users
  do not have to deselect them manually for every scenario.
- **Onboarding banner is now dynamic** — scenario list builds from
  `SAMPLE_DATA_SCENARIOS`, no more hardcoded "four built-in
  scenarios" text. Adding a scenario will not desynchronise docs.
- **Sidebar Data source caption** uses the same dynamic count.
- **Inline `_self_test_sample_data_scenarios`** allows ±5 %
  tolerance on row count for scenarios with `missing_rate > 0`.
- **README.md** scenario list expanded to mention all five built-
  ins (was four).

### Fixed (caught by the v0.2.12 end-to-end gate before release)

- **`calc_interrater_agreement` `observed=False` groupby
  expansion** crashed `pd.insert` with
  `ValueError: Length of values (960) does not match length of
  index (57600)` on the music peer-rating scenario. The
  Categorical-grouped agg was materialising the full 120 ×
  N_context Cartesian even though only 240 cells had data. Fix:
  `observed=True` on both group-by sites so empty cells are
  skipped — correct behaviour for sparse designs and identical for
  fully-crossed ones.
- **`guess_col` greedy substring match** (carryover from v0.2.13):
  reading-testlet column auto-detect picked `Scorer` for the
  `score` pattern, which collided after rename. Fix landed in
  v0.2.13.

### Verification

- 41 inline self-tests pass.
- **226 pytest tests pass** (was 200) — +26: 12 v0.2.14 tests, 7
  binary/testlet (carried), and 7 across the existing files for
  the new scenarios.
- End-to-end AppTest covers all 7 scenarios × Run without any
  exception (≈ 3 min CI runtime).

## 0.2.13-beta - 2026-04-18

Binary / testlet data support, plus two bugs the new end-to-end
AppTest gate (introduced v0.2.12) caught before release. Demonstrates
the gate working as intended: issues surface in CI, not production.

### Added

- **`reading_testlet_binary` sample scenario** — 100 examinees × 1
  scorer × 6 passages × 4 items = 2,400 binary (0/1) observations.
  Grounded in the testlet literature (Wainer & Kiely, 1987;
  Bradlow, Wainer & Wang, 1999). Columns are
  `Person / Scorer / Text / Item / Score`, the tidy-long shape that
  reading / listening / SJT tests produce.
- **Help → Model Capability → "Testlet-format Data"** section
  (~120 lines). Compares MFRM to classical testlet models
  (Bradlow+ 1999) and the bi-factor equivalent (Rijmen, 2010;
  DeMars, 2006). Includes a recommendation flow for when to use
  MFRM's fixed-effect Text facet vs. the **TESTLET_RI** /
  **TESTLET_BIFACTOR** Stan downloads for random-effects
  modelling of local item dependence.
- **Five new APA 7 references**:
  Bradlow, Wainer & Wang (1999), DeMars (2006), Rijmen (2010),
  Wainer & Kiely (1987), Wang, Bradlow & Wainer (2002). All
  wired into `_APA_REFERENCE_LIBRARY` + `_CITATION_TO_KEY`.
- **`tests/test_binary_testlet.py`** (7 tests) — schema, JMLE
  convergence on 2,400-row binary data, citation resolution,
  Help-doc content pin, generator determinism.

### Fixed (bugs the new end-to-end gate caught before release)

- **`guess_col` greedy substring match.** Column auto-detection on
  a testlet dataset (`Person, Scorer, Text, Item, Score`) picked
  `Scorer` as the Score column because `"score" in "scorer"` is
  True. After rename-collision, `pd.to_numeric(df["Score"])` hit
  a DataFrame with two "Score" columns and raised `arg must be a
  list, tuple, 1-d array, or Series`. Fix: `guess_col` now makes
  two passes — exact case-insensitive match first, substring match
  only if nothing exact was found.
- **`prepare_mfrm_data` silent rename-collision** — the same class
  of bug could still surface if a user-uploaded dataset had a
  facet literally named "Score" or "Person". Added a
  `df.columns.duplicated()` check with a clear error message.
- **`np.nanvar(..., ddof=1)` runtime warning** on singleton-facet
  groups (the new binary testlet has a single scorer). Pre-check
  `len(finite_est) >= 2` to avoid the numpy warning on
  zero-variance slices.
- **APA citation regex** rejected 3-author forms like
  `(Bradlow, Wainer & Wang, 1999)`. Loosened to accept any
  `(Surname…, YYYY)` shape.

### Unchanged on purpose

- Binary support was *already* present in the core
  RSM/PCM/GPCM estimator — a Score column with 0/1 values and
  `n_cat = 2` produces a single Rasch–Andrich threshold, which
  is mathematically equivalent to the 1PL Rasch model. This
  release just documents and tests that pathway and adds the
  testlet-shaped scenario that exercises it.
- For discrimination-aware binary analysis (2PL), the existing
  sidebar **Advanced models → IRT_2PL_BINARY** Stan generator
  remains the recommended download-only path.

### Verification

- 41 inline self-tests pass.
- **200 pytest tests** pass (was 186) — +14 including 7 new
  binary/testlet tests and the updated regression guards.
- End-to-end AppTest covers all 5 scenarios × Run without
  exception (≈ 2 min CI runtime).

## 0.2.12-beta - 2026-04-18

Process hotfix responding to a user question after the v0.2.11
widget-crash ship-and-patch cycle: *"Why did this happen and aren't
there other bugs like it?"*

Honest root-cause: v0.2.7 added four sample scenarios including
Clinical OSCE (3 competencies) and the Writing Essay default (4
raters). Tests verified the generator functions but never ran the
generated data through the full UI pipeline. This release installs
the CI gate that would have caught the v0.2.10 bug and completes the
systematic audit of related hazards.

### Added

- **`tests/test_end_to_end_scenarios.py`** (5 tests). For every key
  in `SAMPLE_DATA_SCENARIOS`, boot the app via
  `streamlit.testing.v1.AppTest`, switch the data-source radio to
  that scenario, click **Run FACETS-mode estimation**, and assert
  `at.exception` stays empty across the full post-run tab render.
  This would have caught the v0.2.10 `StreamlitValueBelowMinError`
  before release. Total end-to-end runtime ≈ 2 min on a CI runner.
- **Three additional boundary-case tests** in
  `tests/test_small_data_widgets.py`: 1-row forest plot, 1-row
  misfit ranking, 2×2 bias heatmap. These cover user-provided data
  shapes that do not match any built-in scenario but are plausible
  in real research (e.g., a pilot design with 2 raters).

### Audited (no additional fixes required)

A 5-pattern systematic sweep of the whole monolith found no further
reachable crashes:

1. **Data-dependent widget bounds** — 18 call sites scanned. The 3
   unsafe ones were already fixed in v0.2.11; the remaining 15 use
   literal bounds only.
2. **Statistical calls on small samples** (`std(ddof=1)`,
   `quantile`, `percentile`) — 7 call sites. All either guarded by
   `if n > 1` / `if levels > 1` or downstream of a `len(elements) <
   2` return (Facet Equivalence, Rater Dashboard, Q-Q plot).
3. **Division by data length** — 8 call sites. All guarded by
   `if n_ok == 0 / len(data) > 0 / max(total, 1)` or similar.
4. **Facet Equivalence chi-square** — guarded by
   `len(elements) < 2: return` at L16948 so `df_chi = n_elem - 1 ≥
   1`.
5. **Q-Q residual plot `.std(ddof=1)`** — guarded by `residuals.size
   < 20: return` at L18763.

### Process improvement

- Scenario additions now require an end-to-end AppTest gate before
  release. Unit tests for the generator alone are not sufficient:
  the bug class that shipped was "valid DataFrame → widget
  constraint violation", which only the full pipeline exercises.
- Total pytest count: 178 → **186** (+5 end-to-end, +3 boundary).
- 41 inline self-tests continue to pass.

## 0.2.11-beta - 2026-04-18

Widget small-data hotfix. A user reported a hard crash with
`StreamlitValueBelowMinError` on the Visuals tab when running against
a facet with fewer than 5 elements (e.g. the Clinical OSCE scenario
has 3 competencies). The forest-plot widget's `min_value=5` was
incompatible with `value=min(40, len(sub))` when `len(sub) < 5`. An
audit found two sibling widgets with the same data-dependent
min/value hazard.

### Fixed

- **`_draw_measures_forest_plotly`** — when the selected facet has
  fewer than 6 rows, skip the "Show top N" number input entirely
  and plot every row. When it has ≥ 6, clamp the default to the
  widget's own range so `value >= min_value` always holds.
- **`_draw_misfit_ranking`** — the `st.slider("Number of elements
  to show")` had `min_value=5, max_value=len(df)`; on fit tables
  with `len(df) < 5` this flipped to max<min and triggered
  `StreamlitAPIException`. Now auto-skips the widget for tiny
  tables (< 6 rows).
- **Bias heatmap axis-label slider** — `max_axis_default` could
  exceed `max_value` when the pivot was narrow. Now all three of
  (min, max, default) are clamped to the pivot's actual extent,
  and the slider is skipped entirely for ≤ 8-label pivots.

### Added

- **`tests/test_small_data_widgets.py`** (3 tests): directly invokes
  `_draw_measures_forest_plotly` and `_draw_misfit_ranking` on
  3-row fixtures, and runs the full AppTest with the Clinical OSCE
  scenario. Pass condition = no `StreamlitValueBelowMinError` /
  `ValueAboveMaxError` / `APIException` raised. These are the
  cheapest guards against the class of bug that just shipped.

### Unchanged on purpose

- All v0.2.10 one-click PDF / Word / HTML downloads and matplotlib
  figure fallbacks continue to work unchanged.
- Non-affected widgets (maxit, reltol, rating_min, rating_max,
  number-of-latent-classes, minimum-common-anchors, etc.) all have
  literal min/max/value bounds and are not at risk.

## 0.2.10-beta - 2026-04-18

Publication-Document hotfix. On Streamlit Community Cloud, the
Publication PDF shipped without any plots because kaleido 1.x
requires a Chrome binary that Cloud does not install — every figure
call silently returned `None` and the PDF ended up as a
tables-and-text document. This release adds a matplotlib fallback
that never needs Chrome, and streamlines the download UI so users
can grab the PDF in one click.

### Added

- **Matplotlib fallback for publication figures.** Four new helpers
  render the core figures directly from the result / diagnostics
  dicts, bypassing the Plotly → kaleido → Chrome chain entirely:
  - `_mpl_wright_map_png()`
  - `_mpl_fit_scatter_png()`
  - `_mpl_category_probability_png()`
  - `_mpl_facet_distribution_png()`
- **Hybrid selector** `_plotly_or_matplotlib_png()` tries the Plotly
  export first (for parity with the on-screen look when Chrome is
  available) and falls back to matplotlib when it returns `None`. So
  the Publication Document always ships with its core figures,
  regardless of host environment.
- **`matplotlib>=3.8`** added to `requirements.txt`.
- **New pytest file** `tests/test_publication_figures.py` (10 tests)
  pins: every matplotlib helper produces valid PNG bytes, the
  payload integration function emits all 4 figures, the full PDF
  builder's output contains compressed image streams, and empty
  inputs never raise.
- **Strengthened `_self_test_publication_document_pdf`** — the test
  now enforces `len(pdf) > 100 KB` and `b"FlateDecode" in pdf` so
  the regression where figures silently vanished will fire the next
  time it happens.

### Changed

- **One-click Publication Document downloads.** The previous
  two-step "Generate → wait → Download" flow is replaced by a
  single prominent download button per format. The bytes are built
  eagerly on tab entry with a visible spinner, and cached in
  `session_state` keyed by a run-level fingerprint so subsequent
  sidebar reruns reuse the cache rather than rebuilding.
- **PDF button promoted to primary action** (`type="primary"`) and
  relabelled "Download PDF (with plots)" to make its content
  explicit after the figure-regression fix.
- Consistent button labels across formats: "Download Word (.docx)",
  "Download PDF (with plots)", "Download HTML".

### Unchanged on purpose

- The Plotly export path is still tried first when kaleido + Chrome
  are available; it produces PNGs that match the on-screen look
  exactly. The matplotlib fallback only fires when that path
  returns `None`.
- Figure content, captions, and ordering remain identical to v0.2.9.
- Word / HTML builders unchanged apart from sharing the same
  one-click UI pattern; their existing embedding paths already
  benefited from the shared `_publication_figure_payloads()` fix.

## 0.2.9-beta - 2026-04-17

Quality hotfix responding to a UX re-audit. Three gaps identified in
the v0.2.8 scoring review are addressed in this release:
verbose captions that over-explained features, missing ingestion-time
outlier detection, and thin pytest coverage.

### Added

- **Ingestion-time outlier detector** (`detect_data_outliers`).
  Screens uploaded / pasted / sample data for six anomaly types
  before the user clicks Run, and appends findings to the
  readiness panel with appropriate severity:
  1. Out-of-range scores (when rating_min / rating_max are set)
  2. Negative score values (rating scales are conventionally ≥ 0)
  3. Zero-variance persons (all-identical responses)
  4. Zero-variance facet elements (rater / task / item with
     constant scores — catches classic central-tendency raters)
  5. Extreme-frequency persons via Tukey fence (k=3)
  6. Ceiling / floor saturation (> 90 % of scores at scale extremes)
- **Six new pytest files** covering previously untested surfaces.
  Total pytest count: 10 → 165.
  - `tests/test_sample_scenarios.py` (29 tests) — every scenario's
    shape, determinism, citation resolution, PCA-readiness
  - `tests/test_data_outliers.py` (12 tests) — every anomaly type,
    clean-data-zero-findings invariant, defensive empty-data handling
  - `tests/test_readiness_report.py` (9 tests) — severity logic,
    outlier integration, overall-max invariant
  - `tests/test_apa_references.py` (12 tests) — library integrity,
    citation-map round-trip, key-naming convention, sorted output
  - `tests/test_help_popover_library.py` (58 tests) — 18 topics ×
    required-field check, v0.2.6 additions, unknown-key no-op
  - `tests/test_response_model_generators.py` (29 tests) — seven
    Stan generators + dispatcher, brace balance, required blocks
  - `tests/test_config_whitelist.py` (6 tests) — 23-key contract,
    filter drop-through, critical-keys always present
- **`_self_test_data_outlier_detection`** inline self-test (#41)
  pins the six detection types against a synthetic dirty fixture.

### Changed

- **Verbose sidebar captions trimmed.** Twelve captions and `help=`
  strings shortened to 1-sentence form. Full explanations moved to
  Help → Interpretation Guide rather than repeated at every
  widget. Affected widgets: Column mapping caption, Facet selector
  help, Workflow mode help, Model radio help, Advanced-models
  caption, population_formula help, standardize-numeric help,
  categorical-terms help, noncenter_facet help, dummy_facets help,
  positive_facets help, Analysis depth help, residual-PCA help,
  strict-marginal help, anchor-templates caption.

### Unchanged on purpose

- The readiness panel's existing status UI absorbs outlier
  findings without any rendering changes — it was already designed
  for extensible check lists.
- Existing self-tests, sample scenarios, popover library, and
  Essential-mode tab filters all unchanged.

## 0.2.8-beta - 2026-04-17

Sample-scenario UX hotfix. The v0.2.7 scenario registry hid the
switcher behind a two-step flow: users had to first select
"Sample data (built-in)" on a radio, then spot a selectbox below
labelled "Sample scenario". Early testers missed that the selectbox
existed and assumed the app shipped a single sample. This release
pulls the scenario choice out of hiding.

### Changed

- **Flat single-radio data-source picker.** The sidebar now shows
  one radio with six options: the four built-in scenarios (with
  label + observation count) followed by Paste CSV/TSV and Upload
  file. No more two-step "Sample data → Sample scenario" flow.
- **Inline info card** replaces the collapsed expander as the
  primary surface for scenario metadata. The selected scenario's
  name, design dimensions, observation count, and one-line
  diagnostic hint are always visible in the sidebar directly under
  the radio.
- **"Try another scenario" quick-switch buttons.** When any scenario
  is loaded, the sidebar now shows one-click buttons for each of the
  other three scenarios so users can hop between designs without
  scrolling back to the radio. Each button is tooltipped with the
  scenario's short description.

### Added

- **Main-area data banner** (`render_loaded_data_banner`). When a
  built-in scenario is loaded, a one-line callout appears above the
  readiness panel naming the scenario, its dimensions, and pointing
  at the sidebar switcher. Clears automatically when the user
  switches to Paste or Upload so the banner never implies a
  user-provided file came from the literature.
- **Onboarding banner updated** to list all four built-in
  scenarios (Writing essay / Large-scale writing / L2
  speaking / Clinical OSCE) with observation counts, so the
  quickstart flow advertises the scenarios before the user even
  reaches the sidebar.

### Unchanged on purpose

- The "Run with sample data" onboarding quickstart button still
  fires the default writing-essay scenario — first-timers still get
  a small dataset before they meet the switcher.
- All four scenarios, their parameter values, and their APA 7
  references are identical to v0.2.7. Only the UI wrapping changed.

## 0.2.7-beta - 2026-04-17

Sample-data enrichment. Replaces the single hardcoded 960-observation
demo with a registry of four literature-grounded synthetic scenarios,
each sized and parameterised for a different diagnostic emphasis.
Users who want to see residual PCA, bias heatmaps, or the rater-
severity Wright-map callout at publication-realistic scale can now
pick a scenario that actually exercises those paths.

### Added

- **Four built-in sample-data scenarios** selectable from a new
  `Sample scenario` sidebar selectbox:
  - **Writing essay** (30×4×2×4, 960 obs) — the v0.1+ default,
    retained unchanged. (Eckes, 2011; McNamara, 1996; Linacre, 1989)
  - **Large-scale writing** (120×4×2×3, 2,880 obs) — sized for
    stable residual-PCA; one rater deliberately injected at +1.6
    logits of severity so the bias heatmap and misfit ranking have a
    legible signal. (Myford & Wolfe, 2003, 2004; Engelhard, 1994,
    2013; Smith, 2002)
  - **L2 speaking** (80×3×3×5, 3,600 obs) — analytic rubric with
    five criteria; Pronunciation is hardest by design. Tighter rater
    severity spread typical of trained panels.
    (McNamara, 1996; Bachman & Palmer, 1996; Luoma, 2004)
  - **Clinical OSCE** (60×4×5×3, 3,600 obs) — five stations ×
    three competencies on a compact 4-point scale; station difficulty
    is the dominant source of variation. (Downing & Yudkowsky, 2009;
    Tavakol & Dennick, 2011; Wolfe & Song, 2015)
- **"About this dataset" sidebar expander** renders the chosen
  scenario's description, design dimensions, observation count, and
  APA 7 reference list built from the same `_APA_REFERENCE_LIBRARY`
  the Publication Document uses.
- **Five new APA 7 references** added to support the new scenarios:
  Bachman & Palmer (1996), Downing & Yudkowsky (2009), Engelhard
  (2013), Luoma (2004), Tavakol & Dennick (2011). All integrated
  into `_CITATION_TO_KEY` so they are pickupable by the narrative
  citation scanner.
- **`_self_test_sample_data_scenarios`** (test #40) pins every
  scenario's row-count / score-bounds / citation resolution /
  determinism contract and enforces backward compatibility of the
  legacy `sample_mfrm_data()` signature.

### Changed

- `sample_mfrm_data(seed)` is now a thin wrapper over the new
  `sample_mfrm_data_by_key(DEFAULT_SAMPLE_SCENARIO_KEY, seed)`. The
  default scenario still produces the same 960-row writing-essay
  demo as v0.1–v0.2.6, so existing screenshots, saved configs, and
  onboarding tours continue to match.
- New module-level helper `_generate_mfrm_rsm_from_params(params, seed)`
  is the single RSM-family generator every scenario shares; each
  scenario only owns its parameter dict, not the sampling loop.
- Sample-data CSV download filename now encodes the chosen scenario
  (`mfrm_sample_<key>.csv`) so users can re-download each scenario
  without overwriting the previous one.

### Unchanged on purpose

- `cached_sample_mfrm_data(seed)` is retained with the same
  signature for existing call sites; `cached_sample_mfrm_data_by_key`
  is the new per-scenario cache.
- The "Run with sample data" onboarding banner still fires the
  default writing-essay scenario — first-time users get the familiar
  small dataset before they discover the selectbox.

## 0.2.6-beta - 2026-04-17

Progressive-disclosure polish + self-test expansion. Continues the
cognitive-load reduction started in v0.2.5: every major diagnostic
chart now has a just-in-time help popover, the Essential view-density
mode is extended to Help and Report (not just Visuals), and four new
contract-level self-tests pin previously manually-tested code paths.

### Added

- **Help popover library expanded from 13 → 18 topics.** New entries:
  `misfit_ranking`, `facet_distribution`, `rater_agreement`,
  `threshold_map`, `zstd_distribution`. Each follows the existing
  What it shows / How to read / Watch for schema.
- **Five more charts wired to popovers:** ZSTD distribution histogram,
  Top-misfit ranking, Facet element distribution, Step / Threshold
  Ordering, and Inter-Element Agreement. Combined with v0.2.5's
  Wright Map / Scree / Fit-scatter / Category-probability /
  Bias-heatmap wiring, every chart in the core results tabs now
  has a one-click "How to read" affordance next to its subheader.
- **Essential mode extended to Help and Report sub-tabs.** Help drops
  from 10 → 7 tabs in Essential (hides Rating Scale Guide, Model
  Capability, Public Beta). Report drops from 10 → 6 sub-tabs in
  Essential (hides Manuscript Template, Method Appendix, Facet
  Equivalence, Stan Code). Full mode (sidebar toggle) restores all
  sub-tabs for publication-depth work.
- **Four new inline self-tests** raise the total from 35 → 39:
  `_self_test_posterior_load_netcdf` (synthetic 2-chain × 100-draw
  round-trip via `az.from_dict` → `.to_netcdf()` → uploader wrapper,
  portable across arviz 0.x / 1.x),
  `_self_test_posterior_load_cmdstan_csvs` (two minimal CmdStan CSVs
  → `_posterior_load_cmdstan_csvs` → shape assertions),
  `_self_test_config_json_import_whitelist` (23-key frozenset +
  critical-key check + filter drop-through),
  `_self_test_run_history_clear_confirmation` (push N+2 runs, verify
  cap, verify clear + confirm flag independence).
- **`_self_test_essential_mode_tab_filters`** pins the exact
  Essential-visible label sets for the Help and Report tabs so a
  future tab rename or reorder surfaces as a failing test rather
  than a silently-hidden sub-tab.

### Changed

- `_CONFIG_JSON_IMPORT_WHITELIST` extracted from an inline `set` inside
  the Report tab's import handler to a module-level `frozenset` near
  `_RUN_HISTORY_KEY`, so the new whitelist test pins the same object
  the handler uses. No behavioural change — same 23 keys, same filter.
- `show_help_section` and `show_report_section` now build their tab
  lists dynamically based on the View density setting, with a
  `{label: tab}` dict replacing positional `help_tabs[4]` indexing so
  hiding a tab in the middle no longer shifts downstream indices.

### Unchanged on purpose

- All 10 Help tabs and 10 Report sub-tabs remain reachable in Full
  mode — Essential only hides the *default* view.
- Popover library topics from v0.2.5 are untouched; v0.2.6 only
  appends new entries.

## 0.2.5-beta - 2026-04-17

Cognitive-load hotfix. The v0.2.0 → v0.2.4 feature stack added 23
diagnostic visualisations and roughly 40 new widgets, tabs, and
expanders, crossing the threshold where the default view started to
overwhelm first-time users. This release layers three fixes on top
without removing any existing feature.

### Added

- **Just-in-time help popovers** backed by `_HELP_POPOVER_LIBRARY`
  with 13 three-section entries (What it shows / How to read /
  Watch for). The new helper `render_help_popover(topic_key)`
  mounts an `st.popover` button next to chart titles so users get
  focused guidance without leaving the current tab. Wired into the
  Pathway Map, Forest plot, Q-Q plot, ECDF, Observation-coverage
  heatmap, Category-usage bar, Posterior trace, and Posterior
  Rhat/ESS bars. On Streamlit < 1.32 the helper falls back to an
  expander so help is still reachable.
- **View density toggle** in the sidebar (radio: Essential / Full,
  default Essential). Essential hides the three advanced Visuals
  sub-tabs (Forest / Q-Q / ECDF) and shows a compact caption
  pointing users at the Full mode when they need publication-depth
  diagnostics. State survives reruns via
  `st.session_state["app_view_density"]`.
- Inline self-test `_self_test_help_popover_library` (#34) pins the
  13 required topic keys and enforces minimum prose length so the
  library cannot silently empty out.

### Changed

- **First-read guide is collapsed by default** (was unfolded on
  every rerun). Returning users see the result tabs immediately;
  first-time users still see the affordance at the top.
- **Visual interpretation roadmap expander** in the Visuals tab
  defaults to collapsed (was `expanded=True`), so the four core
  diagnostic sub-tabs are the first thing users see.

### Unchanged on purpose

- Readiness panel in Warning / Issue state still renders expanded
  (alerts must stay visible).
- Anchor-issues, submission action plan "first actions", and
  copy-edit wording repair expanders stay expanded (all three are
  actionable blockers before final reporting).

## 0.2.4-beta - 2026-04-17

Visualization-coverage hotfix. Six new diagnostic plots land across
the Data tab, Visuals tab, and Posterior Viewer, closing the gaps
identified by the post-v0.2.3 chart-coverage audit. All new plots
are Plotly-based, respect the v0.2.2 WCAG palette (primary teal
`#0d7a5a`, warn `#c24e00`, neutral `#666666`), and degrade to an
st.info notice when their required data is missing rather than
rendering broken canvases.

### Added — Data tab

- **Observation coverage heatmap** (`_draw_data_coverage_heatmap`).
  Person × first-facet grid, blue cells = ≥ 1 observation, white =
  none. Capped at 80 persons for readability; fires a warning when
  > 40 % of displayed cells are empty (connectivity risk).
- **Category usage bar chart** (`_draw_category_usage_bar`). Raw
  observed frequency per score category, with an optional "split
  by facet" stack. Categories below `max(5, 1%)` of the total are
  called out in a warning so under-used rating categories are
  caught before they cascade into disordered thresholds downstream
  (Linacre 2004).

### Added — Visuals tab (three new sub-tabs)

- **Forest plot** (`_draw_measures_forest_plotly`). Per-facet
  forest of element measures with dot = Estimate, thick bar =
  50 % CI (± 0.67 · SE), thin bar = 95 % CI (± 1.96 · SE). The
  frequentist analogue of the Bayesian forest already in
  Posterior Viewer; covers JMLE / MML users who do not cross into
  Stan.
- **Q-Q plot** (`_draw_residual_qq_plotly`). Standardized residuals
  vs N(0, 1) theoretical quantiles with 45° reference. Requires
  ≥ 20 residuals. Heavy tails flag extreme-score misfit; S-shape
  flags possible residual structure.
- **ECDF of measures** (`_draw_measure_ecdf_plotly`). Empirical
  cumulative distribution of person and facet-element measures on
  a shared logit axis. Flat stretches = measurement gaps; steep
  jumps = clusters. Complements the Wright Map binned histogram.

### Added — Posterior Viewer

- **Rhat / ESS bar chart tab** (`_posterior_rhat_ess_figure`).
  Two-panel horizontal bar chart with dashed reference at
  Rhat = 1.01 and ESS = 400. Bars turn red (#c24e00) when either
  threshold is crossed so convergence issues jump out without
  reading the full summary table row by row. Reuses
  `_posterior_compute_summary` so numbers stay in sync.

## 0.2.3-beta - 2026-04-17

Third-pass hotfix landing the remaining findings from a second,
post-v0.2.2 comprehensive gap audit. Focus: reconciling in-app text
and metadata tables with the feature surface that actually ships.

### Fixed

- `public_release_readiness_table` row "Simulation validation
  templates" had Status predicate reading from a 4-key scripts
  dict while Evidence sentence quoted a 3-row inventory. Aligned
  both on the inventory so "X rows are documented" and the
  Ready / Review flag cannot drift apart.
- `_render_reporting_checklist` two cells were hardcoded to
  unchecked (PCA of residuals, Average measures per category)
  regardless of whether the diagnostic was actually available.
  Wire them to the real data.
- 11 stale tab-reference strings updated to match the v0.2.0 Report
  meta-grouping: `Report → Stan Code` → `Report → Exports → Stan
  Code`, `Report → Tables` → `Report → Tables & checks → Tables`,
  `Report → APA Report` → `Report → Reports → APA Report`, etc.

### Added

- Help tab Quick Start + Troubleshooting documenting every v0.2.x
  feature and diagnostic (shipped in a separate commit on the same
  release branch).
- `README.md` "What's new in v0.2.0 → v0.2.2" covers the four
  feature tracks plus both post-release hotfix passes.
- `SECURITY.md` Posterior Viewer Uploads and Advanced-Model Stan
  Code Downloads sections document the two v0.2.x file-upload /
  file-download surfaces distinct from rating-data handling.
- `public_beta_limitations_table` rows for Advanced models (Stan,
  download-only), Posterior Viewer (upload), and Publication
  Document (Word / PDF / HTML).
- `mfrmr_015_migration_coverage_table` rows for Advanced response
  models (Stan generators only) and Posterior Viewer (diagnostic).
- `visual_method_evidence_table` row for Posterior trace / ridge /
  pair / forest diagnostics with primary references.
- `external_reference_documentation_table` rows for arviz and
  cmdstan / cmdstanpy / rstan so reviewers have explicit hooks.
- Dimensionality readiness row now concatenates
  `diagnostics['pca_reason']` into its Evidence sentence, so the
  readiness table, manuscript claim guide, publication gate, and
  submission action plan all explain WHY PCA was skipped (small
  data / single item-combination / disabled by Analysis depth /
  >95% NaN).
- Onboarding-banner tip line carries a one-line v0.2.x feature
  catalogue (Publication Document / Posterior Viewer / Advanced
  models / Run history) so new users discover them before clicking
  "Got it".
- Self-tests `_self_test_publication_document_html` (#32) and
  `_self_test_diagnose_estimation_error_patterns` (#33).

## 0.2.2-beta - 2026-04-17

Second-pass hotfix completing the remaining MAJOR and MINOR items
from the pre-redeploy 8-agent UX audit.

### Added

- Run breadcrumb banner (`_render_run_breadcrumb`) above the critical-
  warnings row: one compact line showing model / method / convergence
  status / iterations / obs × person × facet counts / 8-char run
  fingerprint. Users flipping between history snapshots finally have
  a persistent "where am I?" signal.
- MFRM / Bayesian quick-reference glossary (`_MFRM_GLOSSARY` +
  `render_glossary_expander`) with 20 curated entries (logit, Infit,
  Outfit, MnSq, ZSTD, step facet, anchor, separation, strata, RSM,
  PCM, GPCM, SE, CI, Rhat, ESS, divergent, E-BFMI, …). Wired into
  the existing Help → Glossary tab so both the quick reference and
  the deep 50+ term glossary are one click away.
- `reorder_measure_columns(df)` + module-level priority list so every
  measure table surfaces Facet / Level / Estimate / SE / CI /
  Infit / Outfit first; Anchor, N, ReliabilityNote, etc. move to the
  tail. Applied to the Combined measures + fit table in the Measures
  tab; callers can opt in one table at a time.
- `style_fit_columns(df)` returns a pandas Styler that colour-codes
  Infit / Outfit mean-square cells by the Wright & Linacre (1994)
  bands (0.5–1.5 inclusive acceptable · >1.5–2.0 noisy · >2.0
  distorting · < 0.5 over-fit). ZSTD columns are excluded (they
  use a different interpretation scale). Applied to the Combined
  measures + fit table.
- Scree plot now draws three reference lines — EV = 1 (expected
  null), EV = 2 (caution), EV = 3 (strong secondary residual signal) —
  matching the thresholds in the chart guide and narrative.
  Previously the EV = 3 threshold was documented but invisible on
  the plot.
- Inline self-tests: `_self_test_style_fit_columns` (#29),
  `_self_test_mfrm_glossary` (#30), `_self_test_reorder_measure_columns`
  (#31).

### Changed

- Stan downloads now advertise `application/x-stan` instead of
  `text/plain`; browsers trigger a save dialog instead of inline
  preview, and Stan-aware editors pick up syntax highlighting.
- Scree plot neutral-threshold line shifted from `#999999` to the
  stronger `#666666` grey to satisfy WCAG AA contrast on white
  (was 2.85 : 1, now 5.7 : 1).
- Posterior Viewer sidebar carries a reassurance caption when a
  FACETS-mode run is in session_state, so users know switching to
  the viewer does not discard their live run.

### Fixed

- `reorder_measure_columns` short-circuited on `df.empty` (True for
  any rowless DataFrame), leaving placeholder / pre-estimation
  tables with the wrong column order. Now guarded on
  `len(df.columns) == 0` so only genuinely schema-less frames
  pass through untouched.

## 0.2.1-beta - 2026-04-17

Pre-redeploy hotfix landing the 8 highest-impact findings from a
parallel UI / UX audit (downloads, tables, charts, information
architecture, microcopy, onboarding, accessibility, performance).

### Added

- Downloadable built-in sample dataset (Download sample CSV) in the
  sidebar when "Sample data (built-in)" is selected, so users can
  inspect the exact column structure the estimator expects before
  preparing their own upload.
- File-size preflight on CSV uploads: ≥50 MB shows a st.warning,
  ≥100 MB shows a st.error, before pandas parsing consumes memory.
- Clear-history confirmation: the Clear history button now uses a
  two-step flow (first click surfaces a warning with explicit
  "Yes, delete all N snapshots" / "Cancel" buttons), preventing
  accidental one-click loss of the run stack.
- Publication Document signpost at the top of the Downloads tab, so
  users expecting a Word/PDF manuscript file are pointed at
  Report → Exports → Publication Document instead of bouncing.
- Read-only companion import for `mfrm_config.json`: upload a
  previously downloaded config file to inspect its settings as a
  table. Does not overwrite sidebar state (most fields are
  run-specific metadata rather than replayable inputs).
- Canonical rounding helper `format_measure_table()` with `3` decimals
  for measure-like columns (Estimate / SE / CI / Infit / Outfit /
  RMSE / Separation / Eigenvalue …) and `4` decimals for probability-
  like columns (Probability / Proportion / Percent / Share / Rate).
- Inline self-tests `_self_test_format_measure_table` (case #28).

### Changed

- `render_chart_guide()` (v0.2.0-beta) was dead code — the
  11-entry library was registered but never rendered. Now called
  after the scree plot, yardstick Wright map, category probability
  curves, pathway map, facet distribution, posterior trace, and
  posterior ridge, so every diagnostic plot carries a consistent
  "How to read this" expander backed by the library.
- Readiness-panel states render as textual [OK] / [CAUTION] / [ISSUE]
  labels at every render site, so users can distinguish severity
  without hue alone.
- Run-history snapshot cap lowered from 10 → 5 to cut
  session-state memory pressure on Streamlit Community Cloud
  (each deep-copy of facets_mode_output can reach 50+ MB on
  realistic datasets).
- Four most prominent empty-state messages in the Measures tab
  ("No person measures available" → "No person measures yet. Click
  Run FACETS-mode estimation …") rewritten to be action-oriented
  and explain when the table becomes available.

### Accessibility

- Global CSS additions via `_inject_desktop_readability_css`:
  - `:focus-visible` outline 3 px #0066cc at 2 px offset for
    keyboard users (WCAG 2.4.7 Focus Visible).
  - `prefers-reduced-motion` suppresses all animations /
    transitions / smooth scrolling for vestibular-sensitive users.
  - `@media (max-width: 899px)` wraps tabs and tightens the
    container padding so the wide layout does not crush tablet /
    iPad content.

## 0.2.0-beta - 2026-04-17

### Added — advanced model Stan code generators (download-only)

- Seven advanced response models registered in `_ADVANCED_RESPONSE_MODELS`:
  `DINA` (CDM; de la Torre 2009), `HRM` (Patz et al. 2002), `TESTLET_RI`
  and `TESTLET_BIFACTOR` (Bradlow 1999 / DeMars 2006), `MIXTURE_RASCH`
  (Rost 1990), `IRT_2PL_BINARY` (Birnbaum 1968), `PAIRWISE_BTL`
  (Bradley & Terry 1952). Each entry carries `family`, `binary`,
  `needs_q_matrix`, `needs_testlet_column`, and `needs_class_count`
  metadata so the sidebar can surface the right extra widgets.
- Per-model Stan code generators:
  - `generate_dina_stan_code(n_items, n_attributes)` — Bernoulli
    likelihood with slip / guess Beta(2, 10) priors, enumerated
    profile class, `log_lik`, and `profile_class` in
    `generated quantities`.
  - `generate_hrm_stan_code(n_categories)` — signal-detection
    parameters phi (accuracy, lognormal) and eta (bias), ordered
    kappa thresholds, `y_rep` posterior predictive.
  - `generate_testlet_stan_code(n_categories, bifactor)` — random
    intercept by default; bifactor adds a `theta_testlet_general`
    latent factor per testlet.
  - `generate_mixture_rasch_stan_code(n_classes)` — latent classes
    with class-specific item difficulties, Dirichlet prior on class
    probabilities.
  - `generate_2pl_binary_stan_code()` — free discrimination
    (lognormal) + difficulty, `y_rep` PPC.
  - `generate_pairwise_btl_stan_code()` — ability vector + pairwise
    Bernoulli logit.
  - `generate_advanced_model_stan_code(name, ...)` — dispatcher.
- `validate_q_matrix(df)` validates DINA Q-matrices (0/1 values,
  attribute coverage, item coverage, per-item / per-attribute row
  and column sums) and returns a structured verdict with messages.
- New sidebar section "Advanced models (Stan, download only)"
  (collapsed by default) with:
  - enable checkbox
  - model family picker (labels from the registry)
  - Q-matrix uploader (DINA), class-count input (Mixture Rasch)
  - "Generate Stan code" button that serialises the Stan program
    into session_state and exposes a Download .stan button.
- Inline self-test `_self_test_advanced_model_generators` walks every
  registered advanced model, verifies that the generated Stan code
  contains each of the four canonical blocks, has balanced braces,
  carries model-specific keyword markers (bernoulli_lpmf for DINA,
  phi for HRM, gamma_testlet for testlets, simplex + dirichlet for
  mixture, alpha_item for 2PL, ability for BTL), and exercises
  `validate_q_matrix` on well-formed and broken inputs.

### Added — UX tweaks

- Toast notification (`st.toast`) fires in addition to the persistent
  `st.status` accordion whenever estimation completes; message variants
  communicate convergence, non-convergence, and failure.
  Users who scroll away from the sidebar still get a completion signal.
- Report tab regrouped into three meta-categories:
  Reports (APA Report, Manuscript Template, Method Appendix, Claim Guide),
  Tables & checks (Tables, Reporting Checklist, Facet Equivalence, Readiness),
  Exports (Stan Code, Publication Document).
  Individual sub-section renderers are unchanged; only the surface shape is.
- Unified chart interpretation helper `render_chart_guide(chart_name)`
  backed by `_CHART_GUIDE_LIBRARY`. Every diagnostic plot can now drop
  a consistent "How to read this" expander with headline + body text
  sourced from a single library (Wright map, pathway map, category
  probability curves, threshold map, ICC, scree, facet distribution,
  reliability, rater agreement, posterior trace, posterior ridge).
- Keyboard shortcut cheat sheet (`render_keyboard_shortcuts_help`) lives
  in a collapsed sidebar expander so the shortcut surface is documented
  once and kept up to date (R rerun, C clear cache, Esc close, ? cheat
  sheet, Tab focus, Enter activate, Ctrl/Cmd+F search). This historical
  surface was removed in the current Unreleased work in favor of persistent
  result navigation and native control behavior.
- Cache-stale detection banner now exposes a one-click **Rerun now**
  button that sets `_facets_mode_force_rerun` and triggers
  `st.rerun()`, so users don't have to scroll back to the sidebar.
- Inline self-test `_self_test_chart_guide_library` pins the chart-guide
  library key set and enforces a minimum text budget per entry so the
  explanatory text doesn't silently empty out over time.

### Added — Posterior Viewer

- New top-level app mode "Posterior Viewer (upload)" exposed via the
  sidebar radio, for inspecting externally-produced posterior draws
  without leaving the browser. Estimation itself still happens locally
  (Streamlit Cloud is used only as a visualisation frontend).
- Upload pipeline supports three formats with auto-detection by file
  extension:
  - CmdStan CSV (one or more per-chain outputs, loaded via
    `arviz.from_cmdstan`).
  - Apache Parquet (any DataFrame with `chain`, `draw`, and one
    column per parameter).
  - ArviZ NetCDF `.nc` (round-tripped via `arviz.from_netcdf`).
- Summary table with mean / sd / median / 5% / 95% / n_eff / Rhat
  (the ESS and Rhat columns come from `arviz.ess` / `arviz.rhat`) plus
  a one-click CSV download.
- HMC transition diagnostics banner with coloured metrics: divergent
  transitions count + %, max-treedepth hits, acceptance rate mean,
  step-size mean, and E-BFMI min across chains. Inline guidance triggers
  when thresholds are crossed (e.g. `adapt_delta 0.95 → 0.99` for
  divergences, `max_treedepth 10 → 12` for treedepth hits, or
  re-parameterisation for E-BFMI < 0.3).
- Plot suite (Plotly, matches the rest of the app):
  - Trace — chain-coloured trace per parameter
  - Ridge — offset KDE ridge across selected parameters
  - Pair  — scatter-matrix with divergence highlights
  - Forest — posterior mean + 50% and 95% CIs
- Parameter multi-select (defaults to the first 6 parameters) drives
  every plot and the summary table.
- Inline self-test `_self_test_posterior_viewer_loaders` round-trips
  synthetic parquet draws through the loader, summary, HMC diagnostics,
  and plot builders. Gracefully skipped when `arviz` / `pyarrow` are
  not yet installed.

### Added — publication export

- Word (.docx) publication document builder with auto-generated abstract,
  exhaustive Methods (estimator / priors / anchor policy / convergence
  criteria / unused-category handling), manuscript-style Results,
  embedded figures (Wright map, fit scatter, category probability curves,
  facet distribution), core tables (element measures, reliability), and
  an APA 7 reference list. Runs lazily via `python-docx` (new dependency).
- PDF publication document builder with the same structure rendered as a
  letter-size reportlab platypus flowable stack (hanging-indent
  references, 6-inch figure embeds, 8 pt table styling). The PDF
  excludes the reproducibility code appendix per the monolith's public
  beta convention; use the Word export for code-bearing variants.
- HTML publication document builder (`build_publication_html_bytes`):
  self-contained single-file HTML with inline CSS, base64 figure embeds,
  and a minimal markdown → HTML converter. Works offline with no external
  assets.
- APA 7 reference-list generator (`collect_cited_references`,
  `build_apa_reference_list`) keyed by the narrative's `(Author, Year)`
  citations, with an always-include anchor set so every publication
  document carries the foundational Rasch references (Andrich 1978,
  Masters 1982, Linacre 1989 / 2024, Myford & Wolfe 2003 / 2004,
  Smith 2002, Wright & Masters 1982).
- New Report tab sub-tab "Publication Document" exposing three
  download paths (Word / PDF / HTML) from identical narrative sources.
- Inline self-tests: `_self_test_apa_reference_list`,
  `_self_test_publication_document_word`,
  `_self_test_publication_document_pdf`.

### Added — PCA + bias/reliability silent-fail diagnostics

- `diagnose_pca_skip_reason` surfaces exactly why residual PCA was
  skipped (small person count, too few item-combinations, all-NaN
  columns, NaN share above 95%, exception, or an Analysis depth that
  disabled PCA) instead of the previous silent "PCA is not available"
  fallback. The Dimensionality tab reads the stored reason and shows
  it with concrete fix suggestions.
- `estimate_bias_interaction` now returns a `_skip_reason`-carrying
  dict on every early-return path (missing result / diagnostics, empty
  observation table, same-facet pair, unknown facet names, no valid
  cells after filtering). `show_bias_section` surfaces the reason
  instead of crashing on the missing `table` key.
- `calc_reliability` records a `ReliabilityNote` column explaining
  why Separation / Reliability is NaN for a given facet ("only N
  level", "all levels have identical measures", "SE zero or NaN",
  "measurement variance ≤ error variance").

### Added — UX layer (st.status, error remedy, readiness, history, comparison)

- Replace the plain `st.spinner` with an expandable `st.status`
  accordion that writes step-by-step progress (Estimate → Diagnostics
  → Report tables → Bias interactions) and the final elapsed seconds.
- Replace the generic estimation "common causes" checklist with a
  pattern-matched diagnosis → action guide
  (`diagnose_estimation_error`): singular matrix → non-centered facet,
  maxit reached → `maxit=1000`, rating scale error → Score column
  check, memory error → reduce facet cardinality, anchor error →
  fix anchor tables, etc. Falls back to the generic list for
  unrecognised failures.
- Pre-estimation data-quality readiness panel above the Run button
  (`build_readiness_report`, `render_readiness_panel`) with a
  [OK] / [CAUTION] / [ISSUE] readiness summary across eight checks (observation
  count, person count, score column dtype, facet count, per-facet
  level count, per-person coverage, column-role overlap).
- Dismissible onboarding banner with a "Run with sample data"
  one-click quickstart that fires estimation with the built-in
  sample dataset and default column mapping.
- Session-scoped run history (`record_run_in_history`,
  `render_run_history_panel`) caps the last 10 runs and offers a
  Restore button per entry plus a "Clear history" control.
- Two-run comparison panel (`render_comparison_panel`,
  `render_comparison_selector`) on top of the history stack:
  convergence / iterations / elapsed delta, element-level Pearson r
  + RMSE scatter with a 45° reference line, facet-level reliability
  diff table, and an automatic interpretation ribbon
  (near-perfect / strong / moderate / weak agreement).

### Dependencies

- Added `python-docx >= 1.0` (Word publication export).
- Added `reportlab >= 4.0` (PDF publication export).
- Added `arviz >= 0.17` (Posterior Viewer — Rhat / ESS / InferenceData I/O).
- Added `netcdf4 >= 1.6` (Posterior Viewer — NetCDF file round-trip).
- Added `pyarrow >= 15.0` (Posterior Viewer — Apache Parquet draws).

### Changed

- App version label: `0.1.2-beta` → `0.2.0-beta`.
- Report tab restructured into three meta-categories (Reports /
  Tables & checks / Exports); individual sub-section renderers
  unchanged.
- Residual-PCA skip-path surfaces a structured reason instead of the
  generic "not available" fallback.
- Estimation failure UI uses pattern-matched specific remedies instead
  of a generic four-bullet checklist.

## Unreleased

No unreleased changes yet.

## 0.1.2-beta - 2026-04-12

### Added

- Guided visual interpretation checklist in the app and downloadable report bundle.
- PCM/GPCM category probability curves can now be inspected by selected step-facet level, not only as averaged curves.
- Category probability curve exports now use the same curve builder as the Visuals tab and include long-form data for all curve scopes.
- Synthetic guided demo report export through `--export-demo-report`, including report tables, method appendix, and interactive category curves.
- Visual method evidence table plus readability safeguards for dense Wright maps, yardsticks, marginal heatmaps, and bias heatmaps.
- Latent-regression covariate type preview and export table for numeric, categorical, and integer-code review decisions.
- Public-beta limitation and release-readiness tables in the app, demo report export, and CLI release check.
- mfrmr 0.1.5 migration coverage table plus external validation artifact checklist, official-reference table, and reviewer template in the parity fixture.
- README preview screenshot generated from the real Streamlit app using built-in synthetic sample data.
- Archived local R cross-check smoke status for generated TAM, mirt, and sirt handoff fixtures.
- Result-aware manuscript claim guide in the Report tab, table downloads, and demo report export.
- Result-aware Markdown manuscript template for Methods, Results, limitations, reviewer preflight checks, and OSF/demo report exports.
- Publication gate summary that aligns APA Report conclusions with readiness checks and manuscript claim guardrails.
- Case-specific interpretation guidance for sparse categories, dimensionality, bias screens, MML marginal checks, rater reliability, and linking claims.
- Result-specific avoid/safer wording repairs in APA Report guidance, manuscript templates, and case-interpretation guidance exports.
- Submission action plan that combines publication gates, readiness checks, claim guardrails, and wording repairs into a prioritized first-read table.
- Desktop readability refinements for result-tab wrapping, wrapped guide tables, and dense or tightly spaced Wright map / yardstick labels.
- Publication-styled figure export bundle with 300 DPI PNG target, matching HTML figures, and a figure-use manifest.
- Archived simulation-validation inventory for external numerical validation planning without bundling private data or machine-specific absolute paths.
- Sanitized Python/R/Julia external Simulation validation templates, downloadable from the app, OSF/demo report exports, and parity fixture exports.

### Changed

- App version label now marks this public beta release as `0.1.2-beta`.
- Plotly/Kaleido floors now target Kaleido 1.x PNG export behavior, with HTML figure fallback when Chrome/Chromium is unavailable.
- Final-report readiness, first-read guidance, and report text now use the same residual, reliability, and PCA thresholds.
- APA Report summary wording now respects claim-guide caveats and avoids broad no-bias claims outside computed screens.
- New/held-out prediction session-state invalidation now fingerprints uploaded file content, not only file name and size.
- Uploaded CSV/TSV reads are now non-destructive, so sidebar previews do not consume the file stream before estimation.
- `make verify` now includes `python streamlit_app.py --release-check`.
- The optional R cross-check scaffold now records package versions and output file manifests for archived validation evidence.

### Validation

- `make verify` passed with 21 built-in self-tests, 10 pytest checks, release check, benchmark smoke, demo report export, and parity fixture export.
- GitHub Actions passed on Python 3.11 and 3.12 for the final pre-release template commit before version finalization.

## 0.1.1-beta - 2026-04-12

### Added

- Deployment guide for local and Streamlit Community Cloud use.
- Cross-package validation contract for TAM, sirt, mirt, FACETS-like, and ConQuest-style comparisons.

### Validation

- `make verify` passed with 18 built-in self-tests, 2 pytest checks, benchmark smoke, and parity fixture export.
- GitHub Actions passed on Python 3.11 and 3.12 for the validation-protocol commit.

## 0.1.0-beta - 2026-04-12

Initial public-beta repository setup.

### Added

- Standalone single-file Streamlit app for Many-Facet Rasch Model workflows.
- RSM, PCM, bounded GPCM, JMLE, MML EM/Direct/Hybrid/Auto, latent regression, plausible values, prediction, simulation, strict marginal diagnostics, residual PCA, anchor audit, and linking review workflows.
- Data privacy warning in the app and README.
- `st.cache_data` / `st.cache_resource` usage for built-in sample data, static bundled assets, read-only core namespace, and bounded export-byte caches.
- Reproducibility/config fingerprints in exported analysis metadata.
- `--doctor`, `--self-test`, `--benchmark-quick`, and `--export-parity-fixture` CLI paths.
- CI smoke workflow for Python 3.11 and 3.12.
- Validation scaffold for optional TAM, mirt, and sirt external checks.
- Anchor templates and a user guideline.
- MIT license file.

### Validation

- Local compile, doctor, self-test, AppTest, benchmark smoke, parity fixture export, optional R scaffold, fresh virtual environment install/self-test, and Streamlit health checks were run during repository preparation.

### Caveats

- Release status is public beta / research preview.
- Exact numerical parity with FACETS, TAM, sirt, mirt, or `mfrmr` is not claimed.
- Licensed under MIT.

## Research status record moved from README (2026-10-14)

The following text preserves the earlier README's research and verification
record. It is historical context, not an installation guide or a statement of
current qualification. Paths marked **local retained evidence, not distributed**
refer to artifacts absent from a fresh clone; their contents and recorded test
counts have not been independently reverified by moving this text. See
[README](README.md) for current usage and [ROADMAP](ROADMAP.md) for the active plan.

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
`critical review` (local retained evidence, not distributed; `validation/cmle_pcm_anchor_cross_engine_20260810/CMLE_PCM_ANCHOR_CROSS_ENGINE_CRITICAL_REVIEW.md`).

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
`one-click critical review` (local retained evidence, not distributed; `validation/cmle_one_click_result_contract_20260810/CMLE_ONE_CLICK_RESULT_CONTRACT_CRITICAL_REVIEW.md`).

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
`private-archive critical review` (local retained evidence, not distributed; `validation/cmle_one_click_archive_identity_20260810/CMLE_ONE_CLICK_PRIVATE_ARCHIVE_CRITICAL_REVIEW.md`).

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
`comprehension-readiness critical review` (local retained evidence, not distributed; `validation/cmle_one_click_comprehension_readiness_20260810/HTML_PREVIEW_INSTRUMENT_CRITICAL_REVIEW.md`).

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
`operations critical review` (local retained evidence, not distributed; `validation/cmle_one_click_cognitive_interview_operations_20260810/COGNITIVE_INTERVIEW_OPERATIONS_CRITICAL_REVIEW.md`).

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
`versioned-gate critical review` (local retained evidence, not distributed; `validation/cmle_one_click_versioned_confirmatory_gate_20260810/VERSIONED_CONFIRMATORY_GATE_CRITICAL_REVIEW.md`).

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
`sample-size sensitivity critical review` (local retained evidence, not distributed; `validation/cmle_one_click_confirmatory_sample_size_sensitivity_20260810/CONFIRMATORY_SAMPLE_SIZE_SENSITIVITY_CRITICAL_REVIEW.md`).

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
`dependence/MNAR critical review` (local retained evidence, not distributed; `validation/cmle_one_click_confirmatory_dependence_stress_20260810/CONFIRMATORY_DEPENDENCE_STRESS_CRITICAL_REVIEW.md`).

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
`protocol-preflight critical review` (local retained evidence, not distributed; `validation/cmle_one_click_confirmatory_protocol_preflight_20260810/CONFIRMATORY_PROTOCOL_PREFLIGHT_CRITICAL_REVIEW.md`)
and the
`private preflight ZIP` (local retained evidence, not distributed; `validation/cmle_one_click_confirmatory_protocol_preflight_20260810/private_confirmatory_protocol_preflight.zip`).

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
`offline workbench` (local retained evidence, not distributed; `validation/cmle_one_click_confirmatory_decision_workbench_20260810/protocol_decision_workbench.html`),
`critical review` (local retained evidence, not distributed; `validation/cmle_one_click_confirmatory_decision_workbench_20260810/CONFIRMATORY_DECISION_WORKBENCH_CRITICAL_REVIEW.md`),
and
`private workbench ZIP` (local retained evidence, not distributed; `validation/cmle_one_click_confirmatory_decision_workbench_20260810/private_confirmatory_decision_workbench.zip`).

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
`SCI-01 offline memo` (local retained evidence, not distributed; `validation/cmle_one_click_confirmatory_sci01_estimand_memo_20260810/sci01_estimand_comparison.html`),
`critical review` (local retained evidence, not distributed; `validation/cmle_one_click_confirmatory_sci01_estimand_memo_20260810/SCI01_ESTIMAND_MEMO_CRITICAL_REVIEW.md`),
and
`same-data figure` (local retained evidence, not distributed; `validation/cmle_one_click_confirmatory_sci01_estimand_memo_20260810/sci01_estimand_projection.png`).

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
`private invalidity workbench` (local retained evidence, not distributed; `validation/cmle_one_click_confirmatory_invalidity_adjudication_20260810/invalidity_adjudication_workbench.html`)
remains unselected and recruitment-blocking.

A prospective machine-readable precision plan and a separate 20-replicate
Python pilot are now retained under `validation/`. The pilot completed 160/160
fits in about 82 recorded fit-seconds, but it blocks rather than promotes the
study: only 124/160 focal decisions met the registered app eligibility rule,
the two sparse conditions each produced 4/20 eligible decisions, and 0/160
fits passed the post-pilot `1e-4` terminal-gradient sup-norm sensitivity.
Row-level floating-point auditing also found five fit classifications that
would change if 3-decimal display values drove the decision. See
`PILOT20_ASSESSMENT.md` (local retained evidence, not distributed; `validation/operating_characteristics_pilot20_20260809/PILOT20_ASSESSMENT.md`).
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
`FIXED_CALIBRATION_WLE_DOWNSTREAM_RESULTS.md` (local retained evidence, not distributed; `validation/fixed_calibration_wle_downstream_20260809/FIXED_CALIBRATION_WLE_DOWNSTREAM_RESULTS.md`).
The first explicit exact-CMLE-to-WLE fit-sample bridge also passes its frozen
RSM/PCM handoff checks with zero bridge/direct-score and category-surface
differences; it remains repository-only and its WLE SE treats CMLE calibration
as fixed. See
`CMLE_WLE_BRIDGE_RESULTS.md` (local retained evidence, not distributed; `validation/cmle_wle_bridge_20260810/CMLE_WLE_BRIDGE_RESULTS.md`).
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
`CMLE_WLE_CALIBRATION_SENSITIVITY_RESULTS.md` (local retained evidence, not distributed; `validation/cmle_wle_calibration_sensitivity_20260810/CMLE_WLE_CALIBRATION_SENSITIVITY_RESULTS.md`).
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
`CMLE_WLE_BOOTSTRAP_PILOT_RESULTS.md` (local retained evidence, not distributed; `validation/cmle_wle_bootstrap_pilot_20260810/CMLE_WLE_BOOTSTRAP_PILOT_RESULTS.md`).

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
`CMLE_WLE_PERSON_FIT_RESULTS.md` (local retained evidence, not distributed; `validation/cmle_wle_person_fit_20260810/CMLE_WLE_PERSON_FIT_RESULTS.md`)
and
`CMLE_WLE_BOOTSTRAP_FIT_EXTENSION_RESULTS.md` (local retained evidence, not distributed; `validation/cmle_wle_bootstrap_fit_extension_corrected_20260810/CMLE_WLE_BOOTSTRAP_FIT_EXTENSION_RESULTS.md`).
The accompanying
`critical decomposition` (local retained evidence, not distributed; `validation/cmle_wle_bootstrap_fit_extension_corrected_20260810/CMLE_WLE_BOOTSTRAP_FIT_CRITICAL_REVIEW.md`)
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
`threshold surface` (local retained evidence, not distributed; `validation/cmle_wle_fit_threshold_surface_20260810/CMLE_WLE_FIT_THRESHOLD_SURFACE_RESULTS.md`)
and
`class-decomposition addendum` (local retained evidence, not distributed; `validation/cmle_wle_fit_threshold_surface_severity_20260810/CMLE_WLE_FIT_THRESHOLD_SEVERITY_ADDENDUM.md`).

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
`known-truth pilot` (local retained evidence, not distributed; `validation/cmle_wle_fit_known_truth_pilot_20260810/CMLE_WLE_FIT_KNOWN_TRUTH_PILOT_RESULTS.md`),
`bias/dimension addendum` (local retained evidence, not distributed; `validation/cmle_wle_fit_known_truth_pilot_addendum_20260810/CMLE_WLE_FIT_KNOWN_TRUTH_PILOT_ADDENDUM.md`),
and
`critical review` (local retained evidence, not distributed; `validation/cmle_wle_fit_known_truth_pilot_addendum_20260810/CMLE_WLE_FIT_KNOWN_TRUTH_CRITICAL_REVIEW.md`).

The next repository-only hard-anchor increment now implements native hard facet
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
See the `hard-anchor smoke` (local retained evidence, not distributed; `validation/cmle_native_hard_anchor_smoke_20260810/CMLE_NATIVE_HARD_ANCHOR_SMOKE.md`)
and `common-shift diagnostic` (local retained evidence, not distributed; `validation/cmle_native_hard_anchor_common_shift_20260810/CMLE_HARD_ANCHOR_COMMON_SHIFT_DIAGNOSTIC.md`).

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
`registered stress report` (local retained evidence, not distributed; `validation/cmle_native_anchor_differential_stress_20260810/CMLE_NATIVE_ANCHOR_DIFFERENTIAL_STRESS.md`)
and its
`critical review` (local retained evidence, not distributed; `validation/cmle_native_anchor_differential_stress_20260810/CMLE_NATIVE_ANCHOR_DIFFERENTIAL_CRITICAL_REVIEW.md`).

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
`connectivity report` (local retained evidence, not distributed; `validation/cmle_native_anchor_connectivity_20260810/CMLE_NATIVE_ANCHOR_CONNECTIVITY_STRESS.md`)
and
`critical review` (local retained evidence, not distributed; `validation/cmle_native_anchor_connectivity_20260810/CMLE_NATIVE_ANCHOR_CONNECTIVITY_CRITICAL_REVIEW.md`).

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
`remediation report` (local retained evidence, not distributed; `validation/cmle_finite_domain_optimizer_remediation_20260810/CMLE_FINITE_DOMAIN_OPTIMIZER_REMEDIATION.md`)
and
`critical review` (local retained evidence, not distributed; `validation/cmle_finite_domain_optimizer_remediation_20260810/CMLE_FINITE_DOMAIN_OPTIMIZER_REMEDIATION_CRITICAL_REVIEW.md`).

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
`oracle stress report` (local retained evidence, not distributed; `validation/cmle_finite_mle_oracle_20260810/CMLE_FINITE_MLE_SUPPORT_ORACLE_STRESS.md`)
and
`critical review` (local retained evidence, not distributed; `validation/cmle_finite_mle_oracle_20260810/CMLE_FINITE_MLE_SUPPORT_ORACLE_CRITICAL_REVIEW.md`).
The oracle is now connected to the repository-only `fit_cmle` readiness path
under a third prospective contract. All 287 same-byte eligible fits returned;
274 oracle-interior fits remained ready and all 13 boundary fits remained
non-ready with a dedicated reason. Independent oracle status matched 287/287,
102 selected tests passed, and 82,200 Person rows retained every raw and
rounded-counterfactual fit decision. Estimate, SE, likelihood, WLE, Infit, and
Outfit changes were at most `2.27e-13`. The core still retains technical
optimizer outputs on a boundary for research reproducibility; a future UI
orchestrator should stop before optimization and explain the boundary. See the
`integration report` (local retained evidence, not distributed; `validation/cmle_finite_mle_readiness_integration_20260810/CMLE_FINITE_MLE_READINESS_INTEGRATION.md`)
and
`integration critical review` (local retained evidence, not distributed; `validation/cmle_finite_mle_readiness_integration_20260810/CMLE_FINITE_MLE_READINESS_INTEGRATION_CRITICAL_REVIEW.md`).
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
`workflow report` (local retained evidence, not distributed; `validation/cmle_structured_workflow_20260810/CMLE_STRUCTURED_WORKFLOW_STRESS.md`)
and
`workflow critical review` (local retained evidence, not distributed; `validation/cmle_structured_workflow_20260810/CMLE_STRUCTURED_WORKFLOW_CRITICAL_REVIEW.md`).

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
`cross-engine decision` (local retained evidence, not distributed; `validation/cmle_cross_engine_boundary_20260810/decision.json`)
and
`critical review` (local retained evidence, not distributed; `validation/cmle_cross_engine_boundary_20260810/CMLE_CROSS_ENGINE_BOUNDARY_CRITICAL_REVIEW.md`).

Before using results for high-stakes scoring, placement, certification, employment, or institutional decisions, cross-check the analysis with an established workflow and document the model assumptions.


## What's new in the current beta line

The current app label is v0.2.18-beta. This beta line adds the v0.2.14
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
