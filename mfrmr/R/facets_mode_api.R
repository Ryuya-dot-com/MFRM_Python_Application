# FACETS-mode API wrappers (defaults to public-spec: RSM + JML/JMLE path)

normalize_facets_mode_data <- function(data) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame. Got: ", class(data)[1], ". ",
         "Convert with as.data.frame() if needed.", call. = FALSE)
  }
  out <- as.data.frame(data, stringsAsFactors = FALSE)
  if (nrow(out) == 0) {
    stop("Input data is empty. ",
         "Supply a data.frame with at least one row of observations.",
         call. = FALSE)
  }
  out
}

normalize_spec_input <- function(x, arg_name) {
  if (is.null(x)) return(NULL)
  if (!is.data.frame(x)) stop("`", arg_name, "` must be NULL or a data.frame.", call. = FALSE)
  as.data.frame(x, stringsAsFactors = FALSE)
}

infer_facets_mode_mapping <- function(dat, person = NULL, facets = NULL, score = NULL, weight = NULL) {
  cols <- names(dat)
  if (length(cols) < 3) {
    stop("Input data must have at least 3 columns ",
         "(person, score, and at least one facet). Got: ", length(cols), ".",
         call. = FALSE)
  }

  person_col <- if (!is.null(person) && nzchar(person)) {
    person
  } else {
    guess_col(cols, c("person", "participant", "student"), fallback = 1)
  }
  score_col <- if (!is.null(score) && nzchar(score)) {
    score
  } else {
    guess_col(cols, c("score", "rating", "mark"), fallback = min(2L, length(cols)))
  }

  if (!person_col %in% cols) stop("Person column not found: ", person_col, call. = FALSE)
  if (!score_col %in% cols) stop("Score column not found: ", score_col, call. = FALSE)

  weight_col <- NULL
  if (!is.null(weight) && nzchar(weight)) {
    if (!weight %in% cols) stop("Weight column not found: ", weight, call. = FALSE)
    weight_col <- weight
  }

  facet_cols <- facets
  if (is.null(facet_cols) || length(facet_cols) == 0) {
    blocked <- unique(c(person_col, score_col, weight_col))
    candidates <- setdiff(cols, blocked)
    preferred <- candidates[stringr::str_detect(tolower(candidates), "rater|task|criterion|criteria|facet")]
    facet_cols <- if (length(preferred) > 0) preferred else candidates
  }
  facet_cols <- as.character(facet_cols)
  if (length(facet_cols) == 0) {
    stop("No facet columns detected. Supply `facets` explicitly.", call. = FALSE)
  }
  missing_facets <- setdiff(facet_cols, cols)
  if (length(missing_facets) > 0) {
    stop("Facet column(s) not found: ", paste(missing_facets, collapse = ", "), call. = FALSE)
  }

  list(person = person_col, facets = facet_cols, score = score_col, weight = weight_col)
}

#' Run FACETS-mode estimation workflow (public-spec mode)
#'
#' This helper mirrors `mfrmRFacets.R` behavior as a package API and keeps
#' FACETS-like defaults (`model = "RSM"`, `method = "JML"`), while allowing
#' users to choose compatible estimation options.
#'
#' `run_mfrm_facets()` is intended as a one-shot workflow helper:
#' fit -> diagnostics -> key report tables.
#' Returned objects can be inspected with `summary()` and `plot()`.
#'
#' @param data A data.frame in long format.
#' @param person Optional person column name. If `NULL`, guessed from names.
#' @param facets Optional facet column names. If `NULL`, inferred from remaining
#'   columns after person/score/weight mapping.
#' @param score Optional score column name. If `NULL`, guessed from names.
#' @param weight Optional weight column name.
#' @param keep_original Passed to [fit_mfrm()].
#' @param model MFRM model (`"RSM"` default, or `"PCM"`).
#' @param method Estimation method (`"JML"` default; `"JMLE"` and `"MML"` also supported).
#' @param step_facet Step facet for PCM mode; passed to [fit_mfrm()].
#' @param anchors Optional anchor table (data.frame).
#' @param group_anchors Optional group-anchor table (data.frame).
#' @param noncenter_facet Non-centered facet passed to [fit_mfrm()].
#' @param dummy_facets Optional dummy facets fixed at zero.
#' @param positive_facets Optional facets with positive orientation.
#' @param quad_points Quadrature points for MML; passed to [fit_mfrm()].
#' @param maxit Maximum optimizer iterations.
#' @param reltol Optimization tolerance.
#' @param top_n_interactions Number of rows for interaction diagnostics.
#'
#' @return A list with components:
#' - `fit`: [fit_mfrm()] result
#' - `diagnostics`: [diagnose_mfrm()] result
#' - `iteration`: [estimation_iteration_report()] result
#' - `fair_average`: [fair_average_table()] result
#' - `rating_scale`: [rating_scale_table()] result
#' - `run_info`: run metadata table
#' - `mapping`: resolved column mapping
#'
#' @section Estimation-method notes:
#' - `method = "JML"` (default): FACETS-style joint estimation route.
#' - `method = "JMLE"`: explicit JMLE label; internally equivalent to JML route.
#' - `method = "MML"`: marginal maximum likelihood route using `quad_points`.
#'
#' `model = "PCM"` is supported; set `step_facet` when facet-specific step
#' structure is needed.
#'
#' @section Visualization:
#' - `plot(out, type = "fit")` delegates to [plot.mfrm_fit()] and returns
#'   fit-level visual bundles (e.g., Wright/pathway/CCC).
#' - `plot(out, type = "qc")` delegates to [plot_qc_dashboard()] and returns
#'   a QC dashboard plot object.
#'
#' @section Interpreting output:
#' Start with `summary(out)`:
#' - check convergence and iteration count in `overview`.
#' - confirm resolved columns in `mapping`.
#'
#' Then inspect:
#' - `out$rating_scale` for category/threshold behavior.
#' - `out$fair_average` for observed-vs-model scoring tendencies.
#' - `out$diagnostics` for misfit/reliability/interactions.
#'
#' @section Typical workflow:
#' 1. Run `run_mfrm_facets()` with explicit column mapping.
#' 2. Check `summary(out)` and `summary(out$diagnostics)`.
#' 3. Visualize with `plot(out, type = "fit")` and `plot(out, type = "qc")`.
#' 4. Export selected tables for reporting (`out$rating_scale`, `out$fair_average`).
#'
#' @seealso [fit_mfrm()], [diagnose_mfrm()], [estimation_iteration_report()],
#'   [fair_average_table()], [rating_scale_table()]
#'
#' @examples
#' toy <- mfrmr:::sample_mfrm_data(seed = 123)
#'
#' # FACETS-like default: RSM + JML
#' out <- run_mfrm_facets(
#'   data = toy,
#'   person = "Person",
#'   facets = c("Rater", "Task", "Criterion"),
#'   score = "Score",
#'   maxit = 15
#' )
#' names(out)
#' out$fit$summary[, c("Model", "Method")]
#' s <- summary(out)
#' s$overview[, c("Model", "Method", "Converged")]
#' p_fit <- plot(out, type = "fit", draw = FALSE)
#' class(p_fit)
#' p_qc <- plot(out, type = "qc", draw = FALSE)
#' class(p_qc)
#'
#' # Optional: MML route
#' out_mml <- run_mfrm_facets(
#'   data = toy,
#'   person = "Person",
#'   facets = c("Rater", "Task", "Criterion"),
#'   score = "Score",
#'   method = "MML",
#'   quad_points = 7,
#'   maxit = 15
#' )
#' out_mml$fit$summary[, c("Model", "Method")]
#' @export
run_mfrm_facets <- function(data,
                            person = NULL,
                            facets = NULL,
                            score = NULL,
                            weight = NULL,
                            keep_original = FALSE,
                            model = c("RSM", "PCM"),
                            method = c("JML", "JMLE", "MML"),
                            step_facet = NULL,
                            anchors = NULL,
                            group_anchors = NULL,
                            noncenter_facet = "Person",
                            dummy_facets = NULL,
                            positive_facets = NULL,
                            quad_points = 15,
                            maxit = 400,
                            reltol = 1e-6,
                            top_n_interactions = 20L) {
  model <- toupper(match.arg(model))
  method <- toupper(match.arg(method))

  dat <- normalize_facets_mode_data(data)
  mapping <- infer_facets_mode_mapping(
    dat = dat,
    person = person,
    facets = facets,
    score = score,
    weight = weight
  )

  anchor_df <- normalize_spec_input(anchors, "anchors")
  group_anchor_df <- normalize_spec_input(group_anchors, "group_anchors")

  fit <- fit_mfrm(
    data = dat,
    person = mapping$person,
    facets = mapping$facets,
    score = mapping$score,
    weight = mapping$weight,
    keep_original = keep_original,
    model = model,
    method = method,
    step_facet = step_facet,
    anchors = anchor_df,
    group_anchors = group_anchor_df,
    noncenter_facet = noncenter_facet,
    dummy_facets = dummy_facets,
    positive_facets = positive_facets,
    quad_points = as.integer(quad_points),
    maxit = as.integer(maxit),
    reltol = as.numeric(reltol)
  )

  diagnostics <- diagnose_mfrm(fit, top_n_interactions = as.integer(top_n_interactions))
  iter <- estimation_iteration_report(fit)
  fair <- fair_average_table(fit, diagnostics = diagnostics)
  rating <- rating_scale_table(fit, diagnostics = diagnostics)

  run_info <- data.frame(
    key = c(
      "person", "score", "facets", "weight",
      "model", "method_input", "method_used", "step_facet", "noncenter_facet",
      "dummy_facets", "positive_facets", "keep_original", "quad_points",
      "maxit", "reltol", "top_n_interactions"
    ),
    value = c(
      mapping$person,
      mapping$score,
      paste(mapping$facets, collapse = ","),
      if (is.null(mapping$weight)) "" else mapping$weight,
      model,
      method,
      as.character(fit$summary$Method[[1]]),
      if (is.null(step_facet)) "" else as.character(step_facet),
      noncenter_facet,
      paste(if (is.null(dummy_facets)) character(0) else dummy_facets, collapse = ","),
      paste(if (is.null(positive_facets)) character(0) else positive_facets, collapse = ","),
      as.character(isTRUE(keep_original)),
      as.character(as.integer(quad_points)),
      as.character(as.integer(maxit)),
      as.character(as.numeric(reltol)),
      as.character(as.integer(top_n_interactions))
    ),
    stringsAsFactors = FALSE
  )

  out <- list(
    fit = fit,
    diagnostics = diagnostics,
    iteration = iter,
    fair_average = fair,
    rating_scale = rating,
    run_info = run_info,
    mapping = mapping
  )
  class(out) <- c("mfrm_facets_run", class(out))

  out
}

#' Compatibility alias for FACETS-mode workflow
#'
#' @rdname run_mfrm_facets
#' @export
mfrmRFacets <- function(data,
                        person = NULL,
                        facets = NULL,
                        score = NULL,
                        weight = NULL,
                        keep_original = FALSE,
                        model = c("RSM", "PCM"),
                        method = c("JML", "JMLE", "MML"),
                        step_facet = NULL,
                        anchors = NULL,
                        group_anchors = NULL,
                        noncenter_facet = "Person",
                        dummy_facets = NULL,
                        positive_facets = NULL,
                        quad_points = 15,
                        maxit = 400,
                        reltol = 1e-6,
                        top_n_interactions = 20L) {
  run_mfrm_facets(
    data = data,
    person = person,
    facets = facets,
    score = score,
    weight = weight,
    keep_original = keep_original,
    model = model,
    method = method,
    step_facet = step_facet,
    anchors = anchors,
    group_anchors = group_anchors,
    noncenter_facet = noncenter_facet,
    dummy_facets = dummy_facets,
    positive_facets = positive_facets,
    quad_points = quad_points,
    maxit = maxit,
    reltol = reltol,
    top_n_interactions = top_n_interactions
  )
}
