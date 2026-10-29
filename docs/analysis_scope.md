# Analyses, interpretation and reporting limits

[Start here](../README.md) · [Windows guide](desktop_user_guide.en.md) · [日本語の操作・結果の読み方](desktop_user_guide.ja.md#分析を選ぶ)

Use this page to choose an analysis and understand what its results can support.
It describes the Web application and the desktop preview documented in this
source branch. The desktop preview includes RSM/PCM, bias diagnostics, selected
intervals, joint interactions and G/D studies; it is not yet a general desktop
release. The published Web release and this desktop preview do not have
identical controls or export formats.

## Choose the analysis for your question

| Question | Analysis | What to distinguish |
|---|---|---|
| How do person measures, rater severity and task difficulty compare? | Many-facet RSM or PCM | RSM shares category steps. PCM allows separate steps for levels of the chosen step facet. A main rater effect is not a rater-by-task interaction. |
| Does a particular combination of facet levels differ from the main-effects pattern? | Conditional bias diagnostics | These hold the original calibration fixed. They identify patterns to investigate, without confirming the presence or absence of bias. |
| What is the uncertainty in a facet level or a planned difference? | Specified structural intervals for fixed-SD MML | Only targets passing the checks receive intervals. A difference uses the full covariance, not the sum of marginal SEs. |
| What happens when an interaction is fitted together with main effects and steps? | Desktop joint interaction analysis | This refits a separate model. It does not turn the fixed-calibration diagnostic into a joint test. |
| How would numbers of raters or tasks affect score reliability? | G study and D study | These use observed scores and a declared random-facet design, not the MFRM logit scale. |

The Web interface also provides advanced workflows such as bounded GPCM,
latent regression, anchors and linking review, residual PCA, prediction and
simulation. Their availability does not establish that every displayed
statistic can support inference. Read the explanation for the selected output.

## Estimation and standard errors

JMLE jointly estimates person and facet parameters. Desktop MML uses EM with
a fixed normal population SD; the SD is a setting, not an estimated result.
The Web interface offers additional MML controls, including optional population
SD estimation, whose uncertainty and model-ranking uses remain unavailable.

A standard error must be interpreted for its stated target and method.
Conditional-information errors hold other parameters fixed. Structural
observed-information errors use a covariance calculation for the model
parameters. An MML person's EAP posterior SD describes a posterior distribution
with calibration fixed; it is not the ML standard error of a structural effect.
Do not treat these quantities as interchangeable because each is labelled SE.

The application does not support inferential use of general MML tables,
Person-score uncertainty, estimated population-SD uncertainty or automatic
equivalence conclusions. Selected structural intervals have the separate,
limited scope below. An available selected interval does not remove the
restrictions on other outputs.

## Specified structural intervals

For supported unit-weight, fixed-SD RSM/PCM analyses, this output provides
pointwise asymptotic normal intervals for selected levels of a non-person
facet and a prespecified difference. It uses the full structural covariance.

The calculation checks termination, identification, parameter boundaries,
observed information and sensitivity to quadrature settings. These checks
assess numerical suitability. They do not prove that the response model fits
the population, give a strict quadrature error bound, or guarantee
finite-sample coverage.

The approximation assumes the specified response model and fixed normal
population, independent persons, ignorable assignment/missingness, an
identified interior parameter and growing information. It does not provide
simultaneous intervals, Person-score intervals or a test of equivalence.

The [English](fixed_sd_quadrature_check.en.md) and
[Japanese](fixed_sd_quadrature_check.ja.md) Web guides explain selecting,
saving and rechecking these intervals. In Desktop, register the requests in
**Settings → Bias / intervals** before running the analysis.

## Bias diagnostics and joint interactions

Conditional bias diagnostics describe cells with the original calibration
held fixed. Sparse cells, extreme scores and disconnected rating designs can
limit what can be estimated or compared. Displayed diagnostic values do not
establish statistical significance, absence of bias or assessment fairness.

The desktop joint-interaction option fits main effects, steps and one selected
structural-facet interaction together in a separate RSM/PCM model. It uses
unit weights and a normal person population fixed at mean 0 and SD 1.
Eligible requested cells and prespecified contrasts receive pointwise
approximate Wald intervals and p values. Holm adjustment applies to the
recorded family of requested tests; it does not make the intervals simultaneous.

This model assumes independent persons, the stated response/population model,
ignorable assignment/missingness and an identified interior solution with
growing information. Numerical checks are not a finite-sample error-rate
guarantee. An interaction is evidence about the fitted model, not a verdict
about fairness or a person. The original main-effects results remain separate.

## G studies and D studies

The supported design is completely crossed and balanced, with one or two
random facets besides the object of measurement, at least two levels of each
factor, one observation per cell and unit observation weights. Missing or
repeated cells and unsupported designs produce a reason rather than an
unsupported coefficient.

G describes relative decisions and Phi describes absolute decisions. The D
study changes planned numbers of observations under the stated design and
population assumptions. Negative estimated variance components are shown
alongside the zero-clipped components used for the coefficients. A coefficient
alone does not establish an acceptable assessment or validate a future design.

The desktop G/D entry uses the original observed scores without requiring an
MFRM fit. Multivariate G-theory and uncertainty intervals for G/Phi are not
implemented. These restrictions also apply when G/D is requested alongside MFRM.

## When a result cannot support inference

Read the message for the particular table or figure.

- **Input or design problem:** correct the indicated column mapping, score
  range, missingness or design problem. For example, unrelated rating groups
  may lack a common scale. Keep the original data and document any change.
- **Numerical or boundary problem:** review the stopping reason, identification
  and extreme scores. A finite displayed estimate or an optimizer stopping
  flag does not establish valid inference.
- **Unsupported uncertainty method or use:** the application has not
  established that use for this output. Repeatedly changing the data or
  rerunning the same settings does not resolve that limitation.
- **Changed settings or a saved result from another implementation:** check
  which analysis the displayed output belongs to. Reopened results retain the
  saved values but are checked against current output-specific rules.

A blank SE or interval means it is unavailable, not zero. Do not infer equality,
absence of bias or a reliable ranking from missing uncertainty or a
non-significant result. Use an eligible dedicated output only within its
stated scope.

## Comparing programs and reporting results

Match the model, data, weights, constraints, population distribution, scoring
method and numerical settings before comparing software.

In particular, this application's non-person PTMEA/PtMeaCorr correlates
observed scores with person measures. FACETS PtMea uses the signed sum of
measures for each observation. They are different diagnostics. The application's
RSM/PCM fit ZSTD follows a FACETS-style fourth-moment/Wilson–Hilferty convention
with an absolute cap of 9; the GPCM adaptation is approximate. Agreement in one
diagnostic does not establish equivalence of the programs.

AIC/AICc/BIC, model-choice LR values and automatic model recommendations remain
restricted to technical review. Residual PCA is a diagnostic requiring attention
to sparse observations, overlap and stability, not an automatic dimensionality
decision. Simulations describe their specified data-generating conditions;
they do not validate a different empirical design.

When reporting, retain the model, estimator, constraints, scoring rubric,
population-SD treatment, analysis version, uncertainty method and relevant
limitations. Preserve captions and notes with exported figures and tables.
Software citation does not replace citations for the statistical methods.

Implementation references for technical readers:
[output rules](../mfrm_app/output_qualification.py),
[structural intervals](../mfrm_app/structural_intervals.py),
[joint interactions](../mfrm_app/joint_interaction.py),
[G/D calculations](../mfrm_app/generalizability.py).
The [FACETS comparison](../validation/FACETS_WINDOWS_DESKTOP_REVIEW_20261029.md)
states the conditions and limits of its numerical checks.
