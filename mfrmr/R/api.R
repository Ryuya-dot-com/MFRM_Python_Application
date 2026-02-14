#' Fit a many-facet Rasch model with a flexible number of facets
#'
#' This is the package entry point. It wraps `mfrm_estimate()` and defaults to
#' `method = "MML"`. Any number of facet columns can be supplied via `facets`.
#'
#' @param data A data.frame with one row per observation.
#' @param person Column name for the person (character scalar).
#' @param facets Character vector of facet column names.
#' @param score Column name for observed category score.
#' @param rating_min Optional minimum category value.
#' @param rating_max Optional maximum category value.
#' @param weight Optional weight column name.
#' @param keep_original Keep original category values.
#' @param model `"RSM"` or `"PCM"`.
#' @param method `"MML"` (default) or `"JML"` / `"JMLE"`.
#' @param step_facet Step facet for PCM.
#' @param anchors Optional anchor table.
#' @param group_anchors Optional group-anchor table.
#' @param noncenter_facet One facet to leave non-centered.
#' @param dummy_facets Facets to fix at zero.
#' @param positive_facets Facets with positive orientation.
#' @param quad_points Quadrature points for MML.
#' @param maxit Maximum optimizer iterations.
#' @param reltol Optimization tolerance.
#'
#' @details
#' Data must be in **long format** (one row per observed rating event).
#'
#' Supported model/estimation combinations:
#' - `model = "RSM"` with `method = "MML"` or `"JML"/"JMLE"`
#' - `model = "PCM"` with a designated `step_facet` (defaults to first facet)
#'
#' Anchor inputs are optional:
#' - `anchors` should contain facet/level/fixed-value information.
#' - `group_anchors` should contain facet/level/group/group-value information.
#' Both are normalized internally, so column names can be flexible
#' (`facet`, `level`, `anchor`, `group`, `groupvalue`, etc.).
#'
#' Facet sign orientation:
#' - facets listed in `positive_facets` are treated as `+1`
#' - all other facets are treated as `-1`
#' This affects interpretation of reported facet measures.
#'
#' @return
#' An object of class `mfrm_fit` (named list) with:
#' - `summary`: one-row model summary (`LogLik`, `AIC`, `BIC`, convergence)
#' - `facets$person`: person estimates (`Estimate`; plus `SD` for MML)
#' - `facets$others`: facet-level estimates for each facet
#' - `steps`: estimated threshold/step parameters
#' - `config`: resolved model configuration used for estimation
#' - `prep`: preprocessed data/level metadata
#' - `opt`: raw optimizer result from [stats::optim()]
#'
#' @seealso [diagnose_mfrm()], [estimate_bias()], [build_apa_outputs()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#'
#' fit <- fit_mfrm(
#'   data = toy,
#'   person = "Person",
#'   facets = c("Rater", "Criterion"),
#'   score = "Score",
#'   method = "JML",
#'   model = "RSM",
#'   maxit = 25
#' )
#' fit$summary
fit_mfrm <- function(data,
                     person,
                     facets,
                     score,
                     rating_min = NULL,
                     rating_max = NULL,
                     weight = NULL,
                     keep_original = FALSE,
                     model = c("RSM", "PCM"),
                     method = c("MML", "JML", "JMLE"),
                     step_facet = NULL,
                     anchors = NULL,
                     group_anchors = NULL,
                     noncenter_facet = "Person",
                     dummy_facets = NULL,
                     positive_facets = NULL,
                     quad_points = 15,
                     maxit = 400,
                     reltol = 1e-6) {
  model <- toupper(match.arg(model))
  method <- toupper(match.arg(method))
  method <- ifelse(method == "JML", "JMLE", method)

  fit <- mfrm_estimate(
    data = data,
    person_col = person,
    facet_cols = facets,
    score_col = score,
    rating_min = rating_min,
    rating_max = rating_max,
    weight_col = weight,
    keep_original = keep_original,
    model = model,
    method = method,
    step_facet = step_facet,
    anchor_df = anchors,
    group_anchor_df = group_anchors,
    noncenter_facet = noncenter_facet,
    dummy_facets = dummy_facets,
    positive_facets = positive_facets,
    quad_points = quad_points,
    maxit = maxit,
    reltol = reltol
  )

  class(fit) <- c("mfrm_fit", class(fit))
  fit
}

#' Compute diagnostics for an `mfrm_fit` object
#'
#' @param fit Output from [fit_mfrm()].
#' @param interaction_pairs Optional list of facet pairs.
#' @param top_n_interactions Number of top interactions.
#' @param whexact Use exact ZSTD transformation.
#' @param residual_pca Residual PCA mode: `"none"`, `"overall"`, `"facet"`, or `"both"`.
#' @param pca_max_factors Maximum number of PCA factors to retain per matrix.
#'
#' @details
#' This function computes a diagnostic bundle used by downstream reporting.
#'
#' `interaction_pairs` controls which facet interactions are summarized.
#' Each element can be:
#' - a length-2 character vector such as `c("Rater", "Criterion")`, or
#' - omitted (`NULL`) to let the function select top interactions automatically.
#'
#' Residual PCA behavior:
#' - `"none"`: skip PCA
#' - `"overall"`: compute only overall residual PCA
#' - `"facet"`: compute only facet-specific residual PCA
#' - `"both"`: compute both sets
#'
#' @return
#' A named list (diagnostics object) including:
#' - `obs`: observed/expected/residual-level table
#' - `measures`: facet/person fit table (`Infit`, `Outfit`, `ZSTD`, `PTMEA`)
#' - `overall_fit`: overall fit summary
#' - `fit`: element-level fit diagnostics
#' - `reliability`: separation/reliability by facet
#' - `interactions`: top interaction diagnostics
#' - `residual_pca_overall`: optional overall PCA object
#' - `residual_pca_by_facet`: optional facet PCA objects
#'
#' @seealso [fit_mfrm()], [analyze_residual_pca()], [build_visual_summaries()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#' names(diag)
diagnose_mfrm <- function(fit,
                          interaction_pairs = NULL,
                          top_n_interactions = 20,
                          whexact = FALSE,
                          residual_pca = c("none", "overall", "facet", "both"),
                          pca_max_factors = 10L) {
  residual_pca <- match.arg(tolower(residual_pca), c("none", "overall", "facet", "both"))

  mfrm_diagnostics(
    fit,
    interaction_pairs = interaction_pairs,
    top_n_interactions = top_n_interactions,
    whexact = whexact,
    residual_pca = residual_pca,
    pca_max_factors = pca_max_factors
  )
}

extract_pca_eigenvalues <- function(pca_bundle) {
  if (is.null(pca_bundle)) return(numeric(0))

  eig <- numeric(0)
  if (!is.null(pca_bundle$pca) && "values" %in% names(pca_bundle$pca)) {
    eig <- suppressWarnings(as.numeric(pca_bundle$pca$values))
  }
  if (length(eig) == 0 && !is.null(pca_bundle$cor_matrix)) {
    eig <- tryCatch(
      suppressWarnings(as.numeric(eigen(pca_bundle$cor_matrix, symmetric = TRUE, only.values = TRUE)$values)),
      error = function(e) numeric(0)
    )
  }

  eig[is.finite(eig)]
}

build_pca_variance_table <- function(pca_bundle, facet = NULL) {
  eig <- extract_pca_eigenvalues(pca_bundle)
  if (length(eig) == 0) return(data.frame())

  total <- sum(eig, na.rm = TRUE)
  prop <- if (is.finite(total) && total > 0) eig / total else rep(NA_real_, length(eig))
  out <- data.frame(
    Component = seq_along(eig),
    Eigenvalue = eig,
    Proportion = prop,
    Cumulative = cumsum(prop),
    stringsAsFactors = FALSE
  )
  if (!is.null(facet)) out$Facet <- facet
  out
}

infer_facet_names <- function(diagnostics) {
  if (!is.null(diagnostics$facet_names) && length(diagnostics$facet_names) > 0) {
    return(unique(as.character(diagnostics$facet_names)))
  }

  if (!is.null(diagnostics$measures) && "Facet" %in% names(diagnostics$measures)) {
    f <- unique(as.character(diagnostics$measures$Facet))
    f <- setdiff(f, "Person")
    if (length(f) > 0) return(f)
  }

  if (!is.null(diagnostics$obs)) {
    nm <- names(diagnostics$obs)
    skip <- c(
      "Person", "Score", "Weight", "score_k", "PersonMeasure",
      "Observed", "Expected", "Var", "Residual", "StdResidual", "StdSq"
    )
    f <- setdiff(nm, skip)
    if (length(f) > 0) return(f)
  }

  character(0)
}

#' Run residual PCA for unidimensionality checks
#'
#' FACETS-style residual diagnostics can be inspected in two ways:
#' 1) overall residual PCA on the person x combined-facet matrix
#' 2) facet-specific residual PCA on person x facet-level matrices
#'
#' @param diagnostics Output from [diagnose_mfrm()].
#' @param mode `"overall"`, `"facet"`, or `"both"`.
#' @param facets Optional subset of facets for facet-specific PCA.
#' @param pca_max_factors Maximum number of retained components.
#'
#' @details
#' The function works on standardized residual structures derived from
#' [diagnose_mfrm()].
#'
#' Output tables use:
#' - `Component`: principal-component index (1, 2, ...)
#' - `Eigenvalue`: eigenvalue for each component
#' - `Proportion`: component variance proportion
#' - `Cumulative`: cumulative variance proportion
#'
#' For `mode = "facet"` or `"both"`, `by_facet_table` additionally includes
#' a `Facet` column.
#'
#' @return
#' A named list with:
#' - `mode`: resolved mode used for computation
#' - `facet_names`: facets analyzed
#' - `overall`: overall PCA bundle (or `NULL`)
#' - `by_facet`: named list of facet PCA bundles
#' - `overall_table`: variance table for overall PCA
#' - `by_facet_table`: stacked variance table across facets
#'
#' @seealso [diagnose_mfrm()], [plot_residual_pca()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' diag <- diagnose_mfrm(fit, residual_pca = "both")
#' pca <- analyze_residual_pca(diag, mode = "both")
#' head(pca$overall_table)
analyze_residual_pca <- function(diagnostics,
                                 mode = c("overall", "facet", "both"),
                                 facets = NULL,
                                 pca_max_factors = 10L) {
  mode <- match.arg(tolower(mode), c("overall", "facet", "both"))

  if (is.null(diagnostics$obs) || nrow(diagnostics$obs) == 0) {
    stop("diagnostics$obs is empty. Run diagnose_mfrm() first.")
  }

  facet_names <- infer_facet_names(diagnostics)
  if (!is.null(facets)) {
    facet_names <- intersect(facet_names, as.character(facets))
    if (length(facet_names) == 0) {
      stop("No matching facets found in diagnostics for `facets`.")
    }
  }

  out_overall <- NULL
  out_by_facet <- list()

  if (mode %in% c("overall", "both")) {
    can_reuse <- is.null(facets) && !is.null(diagnostics$residual_pca_overall)
    out_overall <- if (can_reuse) {
      diagnostics$residual_pca_overall
    } else {
      compute_pca_overall(
        obs_df = diagnostics$obs,
        facet_names = facet_names,
        max_factors = pca_max_factors
      )
    }
  }

  if (mode %in% c("facet", "both")) {
    can_reuse <- is.null(facets) && !is.null(diagnostics$residual_pca_by_facet)
    out_by_facet <- if (can_reuse) {
      diagnostics$residual_pca_by_facet
    } else {
      compute_pca_by_facet(
        obs_df = diagnostics$obs,
        facet_names = facet_names,
        max_factors = pca_max_factors
      )
    }

    if (!is.null(facets) && length(out_by_facet) > 0) {
      out_by_facet <- out_by_facet[intersect(names(out_by_facet), facet_names)]
    }
  }

  overall_table <- build_pca_variance_table(out_overall)
  by_facet_table <- if (length(out_by_facet) == 0) {
    data.frame()
  } else {
    tbls <- lapply(names(out_by_facet), function(f) {
      build_pca_variance_table(out_by_facet[[f]], facet = f)
    })
    tbls <- tbls[vapply(tbls, nrow, integer(1)) > 0]
    if (length(tbls) == 0) data.frame() else dplyr::bind_rows(tbls)
  }

  list(
    mode = mode,
    facet_names = facet_names,
    overall = out_overall,
    by_facet = out_by_facet,
    overall_table = overall_table,
    by_facet_table = by_facet_table
  )
}

resolve_pca_input <- function(x) {
  if (is.null(x)) stop("Input cannot be NULL.")
  if (!is.null(x$overall_table) || !is.null(x$by_facet_table)) return(x)
  if (!is.null(x$obs)) return(analyze_residual_pca(x, mode = "both"))
  stop("Input must be diagnostics from diagnose_mfrm() or result from analyze_residual_pca().")
}

extract_loading_table <- function(pca_bundle, component = 1L, top_n = 20L) {
  if (is.null(pca_bundle) || is.null(pca_bundle$pca) || is.null(pca_bundle$pca$loadings)) {
    return(data.frame())
  }

  loads <- tryCatch(as.matrix(unclass(pca_bundle$pca$loadings)), error = function(e) NULL)
  if (is.null(loads) || nrow(loads) == 0) return(data.frame())
  if (component > ncol(loads) || component < 1) return(data.frame())

  vals <- suppressWarnings(as.numeric(loads[, component]))
  vars <- rownames(loads)
  if (is.null(vars)) vars <- paste0("V", seq_along(vals))

  ok <- is.finite(vals)
  if (!any(ok)) return(data.frame())

  tbl <- data.frame(
    Variable = vars[ok],
    Loading = vals[ok],
    stringsAsFactors = FALSE
  )
  tbl <- tbl[order(abs(tbl$Loading), decreasing = TRUE), , drop = FALSE]
  top_n <- max(1L, as.integer(top_n))
  head(tbl, n = min(nrow(tbl), top_n))
}

#' Visualize residual PCA results
#'
#' @param x Output from [analyze_residual_pca()] or [diagnose_mfrm()].
#' @param mode `"overall"` or `"facet"`.
#' @param facet Facet name for `mode = "facet"`.
#' @param plot_type `"scree"` or `"loadings"`.
#' @param component Component index for loadings plot.
#' @param top_n Maximum number of variables shown in loadings plot.
#'
#' @details
#' `x` can be either:
#' - output of [analyze_residual_pca()], or
#' - a diagnostics object from [diagnose_mfrm()] (PCA is computed internally).
#'
#' Plot types:
#' - `"scree"`: component vs eigenvalue line plot
#' - `"loadings"`: horizontal bar chart of top absolute loadings
#'
#' For `mode = "facet"` and `facet = NULL`, the first available facet is used.
#'
#' @return A `plotly` htmlwidget object.
#'
#' @seealso [analyze_residual_pca()], [diagnose_mfrm()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' diag <- diagnose_mfrm(fit, residual_pca = "both")
#' pca <- analyze_residual_pca(diag, mode = "both")
#' plt <- plot_residual_pca(pca, mode = "overall", plot_type = "scree")
#' class(plt)
plot_residual_pca <- function(x,
                              mode = c("overall", "facet"),
                              facet = NULL,
                              plot_type = c("scree", "loadings"),
                              component = 1L,
                              top_n = 20L) {
  mode <- match.arg(tolower(mode), c("overall", "facet"))
  plot_type <- match.arg(tolower(plot_type), c("scree", "loadings"))
  pca_obj <- resolve_pca_input(x)

  pca_bundle <- NULL
  title_suffix <- ""

  if (mode == "overall") {
    pca_bundle <- pca_obj$overall
    title_suffix <- "Overall Residual PCA"
  } else {
    if (is.null(pca_obj$by_facet) || length(pca_obj$by_facet) == 0) {
      stop("No facet-level PCA results available.")
    }
    if (is.null(facet)) facet <- names(pca_obj$by_facet)[1]
    if (!facet %in% names(pca_obj$by_facet)) {
      stop("Requested facet not found in PCA results.")
    }
    pca_bundle <- pca_obj$by_facet[[facet]]
    title_suffix <- paste0("Residual PCA - ", facet)
  }

  if (plot_type == "scree") {
    tbl <- build_pca_variance_table(pca_bundle)
    if (nrow(tbl) == 0) stop("No eigenvalues available for scree plot.")
    p <- plotly::plot_ly(
      data = tbl,
      x = ~Component,
      y = ~Eigenvalue,
      type = "scatter",
      mode = "lines+markers",
      hovertemplate = "Component %{x}<br>Eigenvalue %{y:.3f}<extra></extra>"
    )
    return(plotly::layout(
      p,
      title = paste0(title_suffix, " (Scree)"),
      xaxis = list(title = "Component"),
      yaxis = list(title = "Eigenvalue")
    ))
  }

  load_tbl <- extract_loading_table(
    pca_bundle = pca_bundle,
    component = as.integer(component),
    top_n = as.integer(top_n)
  )
  if (nrow(load_tbl) == 0) stop("No loadings available for the requested component.")

  load_tbl$Direction <- ifelse(load_tbl$Loading >= 0, "Positive", "Negative")
  load_tbl$Variable <- factor(load_tbl$Variable, levels = rev(load_tbl$Variable))

  p <- plotly::plot_ly(
    data = load_tbl,
    x = ~Loading,
    y = ~Variable,
    type = "bar",
    orientation = "h",
    color = ~Direction,
    colors = c("Negative" = "#d95f02", "Positive" = "#1b9e77"),
    hovertemplate = "%{y}<br>Loading %{x:.3f}<extra></extra>"
  )
  plotly::layout(
    p,
    title = paste0(title_suffix, " (Loadings: PC", as.integer(component), ")"),
    xaxis = list(title = "Loading"),
    yaxis = list(title = "")
  )
}

#' Estimate FACETS-style bias/interaction terms iteratively
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Output from [diagnose_mfrm()].
#' @param facet_a First facet name.
#' @param facet_b Second facet name.
#' @param max_abs Bound for absolute bias size.
#' @param omit_extreme Omit extreme-only elements.
#' @param max_iter Iteration cap.
#' @param tol Convergence tolerance.
#'
#' @details
#' The function estimates interaction contrasts for `facet_a x facet_b`
#' with iterative recalibration in a FACETS-like style.
#'
#' Typical usage:
#' 1. fit model via [fit_mfrm()]
#' 2. compute diagnostics via [diagnose_mfrm()]
#' 3. call `estimate_bias()` for one facet pair
#' 4. format output with [build_fixed_reports()]
#'
#' @return
#' A named list with:
#' - `table`: interaction rows with effect size, SE, t, p, fit columns
#' - `summary`: compact summary statistics
#' - `chi_sq`: fixed-effect chi-square style summary
#' - `facet_a`, `facet_b`: analyzed facet names
#' - `iterations`: iteration history/metadata
#'
#' @seealso [build_fixed_reports()], [build_apa_outputs()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#' bias <- estimate_bias(fit, diag, facet_a = "Rater", facet_b = "Criterion", max_iter = 2)
#' names(bias)
estimate_bias <- function(fit,
                          diagnostics,
                          facet_a,
                          facet_b,
                          max_abs = 10,
                          omit_extreme = TRUE,
                          max_iter = 4,
                          tol = 1e-3) {
  estimate_bias_interaction(
    res = fit,
    diagnostics = diagnostics,
    facet_a = facet_a,
    facet_b = facet_b,
    max_abs = max_abs,
    omit_extreme = omit_extreme,
    max_iter = max_iter,
    tol = tol
  )
}

#' Build FACETS-style fixed-width text reports
#'
#' @param bias_results Output from [estimate_bias()].
#' @param target_facet Optional target facet for pairwise contrast table.
#'
#' @details
#' This function generates plain-text, fixed-width output intended to be read in
#' console/log environments or exported into text reports.
#'
#' The pairwise section compares levels inside `target_facet` across contexts
#' defined by the opposite facet.
#'
#' @return
#' A named list with:
#' - `bias_fixed`: fixed-width interaction table text
#' - `pairwise_fixed`: fixed-width pairwise contrast text
#' - `pairwise_table`: underlying pairwise data.frame
#'
#' @seealso [estimate_bias()], [build_apa_outputs()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#' bias <- estimate_bias(fit, diag, facet_a = "Rater", facet_b = "Criterion", max_iter = 2)
#' fixed <- build_fixed_reports(bias)
#' names(fixed)
build_fixed_reports <- function(bias_results, target_facet = NULL) {
  if (is.null(bias_results) || is.null(bias_results$table) || nrow(bias_results$table) == 0) {
    return(list(bias_fixed = "No bias data", pairwise_fixed = "No pairwise data"))
  }

  facet_a <- bias_results$facet_a
  facet_b <- bias_results$facet_b

  tbl_display <- bias_results$table |>
    dplyr::select(
      Sq,
      `Observd Score`,
      `Expctd Score`,
      `Observd Count`,
      `Obs-Exp Average`,
      `Bias Size`,
      `S.E.`,
      t,
      `d.f.`,
      `Prob.`,
      Infit,
      Outfit,
      FacetA_Index,
      FacetA_Level,
      FacetA_Measure,
      FacetB_Index,
      FacetB_Level,
      FacetB_Measure
    ) |>
    dplyr::rename(
      `Model S.E.` = `S.E.`,
      `Infit MnSq` = Infit,
      `Outfit MnSq` = Outfit
    )

  names(tbl_display)[names(tbl_display) == "FacetA_Index"] <- paste0(facet_a, " N")
  names(tbl_display)[names(tbl_display) == "FacetA_Level"] <- facet_a
  names(tbl_display)[names(tbl_display) == "FacetA_Measure"] <- paste0(facet_a, " measr")
  names(tbl_display)[names(tbl_display) == "FacetB_Index"] <- paste0(facet_b, " N")
  names(tbl_display)[names(tbl_display) == "FacetB_Level"] <- facet_b
  names(tbl_display)[names(tbl_display) == "FacetB_Measure"] <- paste0(facet_b, " measr")

  bias_cols <- c(
    "Sq",
    "Observd Score",
    "Expctd Score",
    "Observd Count",
    "Obs-Exp Average",
    "Bias Size",
    "Model S.E.",
    "t",
    "d.f.",
    "Prob.",
    "Infit MnSq",
    "Outfit MnSq",
    paste0(facet_a, " N"),
    facet_a,
    paste0(facet_a, " measr"),
    paste0(facet_b, " N"),
    facet_b,
    paste0(facet_b, " measr")
  )

  bias_formats <- list(
    Sq = "{}",
    `Observd Score` = "{:.2f}",
    `Expctd Score` = "{:.2f}",
    `Observd Count` = "{:.0f}",
    `Obs-Exp Average` = "{:.2f}",
    `Bias Size` = "{:.2f}",
    `Model S.E.` = "{:.2f}",
    t = "{:.2f}",
    `d.f.` = "{:.0f}",
    `Prob.` = "{:.4f}",
    `Infit MnSq` = "{:.2f}",
    `Outfit MnSq` = "{:.2f}"
  )
  bias_formats[[paste0(facet_a, " N")]] <- "{:.0f}"
  bias_formats[[paste0(facet_a, " measr")]] <- "{:.2f}"
  bias_formats[[paste0(facet_b, " N")]] <- "{:.0f}"
  bias_formats[[paste0(facet_b, " measr")]] <- "{:.2f}"

  bias_fixed <- build_bias_fixed_text(
    table_df = tbl_display,
    summary_df = bias_results$summary,
    chi_df = bias_results$chi_sq,
    facet_a = facet_a,
    facet_b = facet_b,
    columns = bias_cols,
    formats = bias_formats
  )

  if (is.null(target_facet)) target_facet <- facet_a
  context_facet <- ifelse(target_facet == facet_a, facet_b, facet_a)
  pairwise_tbl <- calc_bias_pairwise(bias_results$table, target_facet, context_facet)

  pairwise_fixed <- if (nrow(pairwise_tbl) == 0) {
    "No pairwise data"
  } else {
    pair_display <- pairwise_tbl |>
      dplyr::select(
        `Target N`,
        Target,
        `Target Measure`,
        `Target S.E.`,
        `Context1 N`,
        Context1,
        `Local Measure1`,
        SE1,
        `Obs-Exp Avg1`,
        Count1,
        `Context2 N`,
        Context2,
        `Local Measure2`,
        SE2,
        `Obs-Exp Avg2`,
        Count2,
        Contrast,
        SE,
        t,
        `d.f.`,
        `Prob.`
      ) |>
      dplyr::rename(
        `Target Measr` = `Target Measure`,
        `Context1 Measr` = `Local Measure1`,
        `Context1 S.E.` = SE1,
        `Context2 Measr` = `Local Measure2`,
        `Context2 S.E.` = SE2
      )

    pair_cols <- c(
      "Target N", "Target", "Target Measr", "Target S.E.",
      "Context1 N", "Context1", "Context1 Measr", "Context1 S.E.",
      "Obs-Exp Avg1", "Count1", "Context2 N", "Context2", "Context2 Measr", "Context2 S.E.",
      "Obs-Exp Avg2", "Count2", "Contrast", "SE", "t", "d.f.", "Prob."
    )

    pair_formats <- list(
      `Target N` = "{:.0f}",
      `Target Measr` = "{:.2f}",
      `Target S.E.` = "{:.2f}",
      `Context1 N` = "{:.0f}",
      `Context1 Measr` = "{:.2f}",
      `Context1 S.E.` = "{:.2f}",
      `Obs-Exp Avg1` = "{:.2f}",
      Count1 = "{:.0f}",
      `Context2 N` = "{:.0f}",
      `Context2 Measr` = "{:.2f}",
      `Context2 S.E.` = "{:.2f}",
      `Obs-Exp Avg2` = "{:.2f}",
      Count2 = "{:.0f}",
      Contrast = "{:.2f}",
      SE = "{:.2f}",
      t = "{:.2f}",
      `d.f.` = "{:.0f}",
      `Prob.` = "{:.4f}"
    )

    build_pairwise_fixed_text(
      pair_df = pair_display,
      target_facet = target_facet,
      context_facet = context_facet,
      columns = pair_cols,
      formats = pair_formats
    )
  }

  list(
    bias_fixed = bias_fixed,
    pairwise_fixed = pairwise_fixed,
    pairwise_table = pairwise_tbl
  )
}

#' Build APA text outputs from model results
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Output from [diagnose_mfrm()].
#' @param bias_results Optional output from [estimate_bias()].
#' @param context Optional named list for report context.
#' @param whexact Use exact ZSTD transformation.
#'
#' @details
#' `context` is an optional named list for narrative customization.
#' Frequently used fields include:
#' - `assessment`, `setting`, `scale_desc`
#' - `rater_training`, `raters_per_response`
#' - `rater_facet` (used for targeted reliability note text)
#'
#' Output text includes residual PCA interpretation if PCA diagnostics are
#' available in `diagnostics`.
#'
#' @return
#' A named list with:
#' - `report_text`: APA-style Method/Results prose
#' - `table_figure_notes`: consolidated notes for tables/figures
#' - `table_figure_captions`: caption-ready labels
#'
#' @seealso [build_visual_summaries()], [estimate_bias()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' diag <- diagnose_mfrm(fit, residual_pca = "both")
#' apa <- build_apa_outputs(fit, diag)
#' names(apa)
build_apa_outputs <- function(fit,
                              diagnostics,
                              bias_results = NULL,
                              context = list(),
                              whexact = FALSE) {
  list(
    report_text = build_apa_report_text(
      res = fit,
      diagnostics = diagnostics,
      bias_results = bias_results,
      context = context,
      whexact = whexact
    ),
    table_figure_notes = build_apa_table_figure_notes(
      res = fit,
      diagnostics = diagnostics,
      bias_results = bias_results,
      context = context,
      whexact = whexact
    ),
    table_figure_captions = build_apa_table_figure_captions(
      res = fit,
      diagnostics = diagnostics,
      bias_results = bias_results,
      context = context
    )
  )
}

#' List literature-based warning threshold profiles
#'
#' @return A named list with `profiles` (`strict`, `standard`, `lenient`)
#'   and `pca_reference_bands`.
#' @details
#' Use this function to inspect available profile presets before calling
#' [build_visual_summaries()].
#'
#' `profiles` contains thresholds used by warning logic
#' (sample size, fit ratios, PCA cutoffs, etc.).
#' `pca_reference_bands` contains literature-oriented descriptive bands used in
#' summary text.
#'
#' @seealso [build_visual_summaries()]
#' @examples
#' mfrm_threshold_profiles()
mfrm_threshold_profiles <- function() {
  warning_threshold_profiles()
}

#' Build warning and narrative summaries for visual outputs
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Output from [diagnose_mfrm()].
#' @param threshold_profile Threshold profile name (`strict`, `standard`, `lenient`).
#' @param thresholds Optional named overrides for profile thresholds.
#' @param summary_options Summary options for `build_visual_summary_map()`.
#' @param whexact Use exact ZSTD transformation.
#'
#' @details
#' This function returns figure-indexed text maps (`figure1` ... `figure11`)
#' to support dashboard/report rendering without hard-coding narrative strings
#' in UI code.
#'
#' `thresholds` can override any profile field by name. Common overrides:
#' - `n_obs_min`, `n_person_min`
#' - `misfit_ratio_warn`, `zstd2_ratio_warn`, `zstd3_ratio_warn`
#' - `pca_first_eigen_warn`, `pca_first_prop_warn`
#'
#' `summary_options` supports:
#' - `detail`: `"standard"` or `"detailed"`
#' - `max_facet_ranges`: max facet-range snippets shown in figure summaries
#' - `top_misfit_n`: number of top misfit entries included
#'
#' @return
#' A named list with:
#' - `warning_map`: figure-level warning text vectors
#' - `summary_map`: figure-level descriptive text vectors
#'
#' @seealso [mfrm_threshold_profiles()], [build_apa_outputs()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' diag <- diagnose_mfrm(fit, residual_pca = "both")
#' vis <- build_visual_summaries(fit, diag, threshold_profile = "strict")
#' names(vis$warning_map)
build_visual_summaries <- function(fit,
                                   diagnostics,
                                   threshold_profile = "standard",
                                   thresholds = NULL,
                                   summary_options = NULL,
                                   whexact = FALSE) {
  list(
    warning_map = build_visual_warning_map(
      res = fit,
      diagnostics = diagnostics,
      whexact = whexact,
      thresholds = thresholds,
      threshold_profile = threshold_profile
    ),
    summary_map = build_visual_summary_map(
      res = fit,
      diagnostics = diagnostics,
      whexact = whexact,
      options = summary_options,
      thresholds = thresholds,
      threshold_profile = threshold_profile
    )
  )
}

#' @export
print.mfrm_fit <- function(x, ...) {
  if (is.list(x) && !is.null(x$summary) && nrow(x$summary) > 0) {
    cat("mfrm_fit object\n")
    print(x$summary)
  } else {
    cat("mfrm_fit object (empty summary)\n")
  }
  invisible(x)
}
