#' mfrmr: Flexible Many-Facet Rasch Modeling in R
#'
#' @description
#' `mfrmr` provides estimation, diagnostics, and reporting utilities for
#' many-facet Rasch models (MFRM) without relying on FACETS/TAM/sirt backends.
#'
#' @details
#' Recommended workflow:
#'
#' 1. Fit model with [fit_mfrm()]
#' 2. Compute diagnostics with [diagnose_mfrm()]
#' 3. Run residual PCA with [analyze_residual_pca()] if needed
#' 4. Estimate interactions with [estimate_bias()]
#' 5. Build narrative/report outputs with [build_apa_outputs()] and [build_visual_summaries()]
#'
#' Function families:
#' - Model fitting: [fit_mfrm()], [summary.mfrm_fit()], [plot.mfrm_fit()]
#' - FACETS-style workflow wrapper: [run_mfrm_facets()], [mfrmRFacets()]
#' - Diagnostics: [diagnose_mfrm()], `summary(diag)`,
#'   [analyze_residual_pca()], [plot_residual_pca()]
#' - Bias and interaction: [estimate_bias()], `summary(bias)`,
#'   [bias_interaction_report()], [plot_bias_interaction()]
#' - Reporting: [build_apa_outputs()], [build_visual_summaries()], [apa_table()]
#' - Data and anchors: [describe_mfrm_data()], [audit_mfrm_anchors()],
#'   [make_anchor_table()], [load_mfrmr_data()]
#'
#' Data interface:
#' - Input analysis data is long format (one row per observed rating).
#' - Packaged simulation data is available via [load_mfrmr_data()] or `data()`.
#'
#' @section Interpreting output:
#' Core object classes are:
#' - `mfrm_fit`: fitted model parameters and metadata.
#' - `mfrm_diagnostics`: fit/reliability/flag diagnostics.
#' - `mfrm_bias`: interaction bias estimates.
#' - `mfrm_bundle` families: FACETS-style summary/report tables.
#'
#' @section Typical workflow:
#' 1. Prepare long-format data.
#' 2. Fit with [fit_mfrm()].
#' 3. Diagnose with [diagnose_mfrm()].
#' 4. Report with [build_apa_outputs()] and [build_visual_summaries()].
#'
#' @section Statistical background:
#' Key statistics reported throughout the package:
#'
#' **Infit (Information-Weighted Mean Square)**
#'
#' Weighted average of squared standardized residuals, where weights are the
#' model-based variance of each observation.
#' `Infit = sum(Z^2 * Var_i * w_i) / sum(Var_i * w_i)`.
#' Expected value is 1.0 under model fit.  Values below 0.5 suggest overfit
#' (Mead-style responses); values above 1.5 suggest underfit (noise or
#' misfit).  Infit is most sensitive to unexpected patterns among on-target
#' observations (Wright & Masters, 1982).
#'
#' **Outfit (Unweighted Mean Square)**
#'
#' Simple average of squared standardized residuals.
#' `Outfit = sum(Z^2 * w_i) / sum(w_i)`.
#' Same expected value and flagging thresholds as Infit, but more sensitive
#' to extreme off-target outliers (e.g., a high-ability person scoring the
#' lowest category).
#'
#' **ZSTD (Standardized Fit Statistic)**
#'
#' Wilson-Hilferty cube-root transformation that converts the mean-square
#' chi-square ratio to an approximate standard normal deviate:
#' `ZSTD = (MnSq^(1/3) - (1 - 2/(9*df))) / sqrt(2/(9*df))`.
#' Values near 0 indicate expected fit; `|ZSTD| > 2` flags potential misfit.
#' ZSTD is reported alongside every Infit and Outfit value.
#'
#' **PTMEA (Point-Measure Correlation)**
#'
#' Pearson correlation between observed scores and estimated person measures
#' within each facet level.  Positive values indicate that scoring aligns
#' with the latent trait dimension; negative values suggest reversed
#' orientation or scoring errors.
#'
#' **Separation**
#'
#' Ratio of the adjusted standard deviation of element measures to their
#' root-mean-square standard error:
#' `Separation = SD(measures) / RMSE(SE)`.
#' Higher values indicate the facet discriminates more statistically distinct
#' levels along the measured variable.
#'
#' **Reliability**
#'
#' `Reliability = Separation^2 / (1 + Separation^2)`.
#' Analogous to Cronbach's alpha or KR-20 for the reproducibility of element
#' ordering.  Values above 0.8 are generally considered adequate for
#' distinguishing among elements (Linacre, 2025).
#'
#' **Strata**
#'
#' Number of statistically distinguishable groups of elements:
#' `Strata = (4 * Separation + 1) / 3`.
#' Three or more strata are recommended (Wright & Masters, 1982).
#'
#' @section Model selection:
#' **RSM vs PCM**
#'
#' The Rating Scale Model (RSM; Andrich, 1978) assumes all levels of the
#' step facet share identical threshold parameters.  The Partial Credit
#' Model (PCM; Masters, 1982) allows each level of the `step_facet` to have
#' its own set of thresholds.  Use RSM when the rating rubric is identical
#' across all items/criteria; use PCM when different items have different
#' numbers or types of response categories, or when there is reason to
#' believe category boundaries vary by item.
#'
#' **MML vs JML**
#'
#' Marginal Maximum Likelihood (MML) integrates over the person ability
#' distribution using Gauss-Hermite quadrature and does not directly estimate
#' person parameters; person estimates are computed post-hoc via Expected A
#' Posteriori (EAP).  Joint Maximum Likelihood (JML/JMLE) estimates all
#' person and facet parameters simultaneously as fixed effects.
#'
#' MML is generally preferred for smaller samples because it avoids the
#' incidental-parameter problem of JML.  JML is faster and does not assume
#' a normal person distribution, which may be an advantage when the
#' population shape is strongly non-normal.
#'
#' See [fit_mfrm()] for usage.
#'
#' @examples
#' mfrm_threshold_profiles()
#' list_mfrmr_data()
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import purrr
#' @import stringr
#' @import psych
#' @import utils
#' @importFrom rlang .data
#' @importFrom stats na.omit optim rnorm sd setNames uniroot
#'
#' @name mfrmr-package
"_PACKAGE"
