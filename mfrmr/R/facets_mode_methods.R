# S3 methods for run_mfrm_facets() workflow objects

round_numeric_frame <- function(df, digits = 3L) {
  if (!is.data.frame(df) || nrow(df) == 0) return(df)
  out <- df
  num_cols <- vapply(out, is.numeric, logical(1))
  out[num_cols] <- lapply(out[num_cols], round, digits = digits)
  out
}

#' Summarize a FACETS-mode workflow run
#'
#' @param object Output from [run_mfrm_facets()].
#' @param digits Number of digits for numeric rounding in summaries.
#' @param top_n Maximum rows shown in nested preview tables.
#' @param ... Passed through to nested summary methods.
#'
#' @details
#' This method returns a compact cross-object summary that combines:
#' - model overview (`object$fit$summary`)
#' - resolved column mapping
#' - run settings (`run_info`)
#' - nested summaries of `fit` and `diagnostics`
#'
#' @section Interpreting output:
#' - `overview`: convergence, information criteria, and scale size.
#' - `mapping`: sanity check for auto/explicit column mapping.
#' - `fit` / `diagnostics`: drill-down summaries for reporting decisions.
#'
#' @section Typical workflow:
#' 1. Run [run_mfrm_facets()] to execute a one-shot pipeline.
#' 2. Inspect with `summary(out)` for mapping and convergence checks.
#' 3. Review nested objects (`out$fit`, `out$diagnostics`) as needed.
#'
#' @return An object of class `summary.mfrm_facets_run`.
#'
#' @seealso [run_mfrm_facets()], [summary.mfrm_fit()], `summary()`
#'
#' @examples
#' toy <- mfrmr:::sample_mfrm_data(seed = 123)
#' out <- run_mfrm_facets(
#'   data = toy,
#'   person = "Person",
#'   facets = c("Rater", "Task", "Criterion"),
#'   score = "Score",
#'   maxit = 15
#' )
#' s <- summary(out)
#' s$overview[, c("Model", "Method", "Converged")]
#' s$mapping
#' @export
summary.mfrm_facets_run <- function(object, digits = 3, top_n = 10, ...) {
  if (!inherits(object, "mfrm_facets_run")) {
    stop("`object` must be an mfrm_facets_run object from run_mfrm_facets().", call. = FALSE)
  }
  digits <- max(0L, as.integer(digits))
  top_n <- max(1L, as.integer(top_n))

  mapping_tbl <- tibble::tibble(
    Key = c("Person", "Score", "Facets", "Weight"),
    Value = c(
      object$mapping$person,
      object$mapping$score,
      paste(object$mapping$facets, collapse = ", "),
      if (is.null(object$mapping$weight)) "" else object$mapping$weight
    )
  )

  out <- list(
    overview = tibble::as_tibble(object$fit$summary),
    mapping = mapping_tbl,
    run_info = tibble::as_tibble(object$run_info),
    fit = summary(object$fit, digits = digits, top_n = top_n, ...),
    diagnostics = summary(object$diagnostics, digits = digits, top_n = top_n, ...),
    digits = digits
  )
  class(out) <- "summary.mfrm_facets_run"
  out
}

#' @export
print.summary.mfrm_facets_run <- function(x, ...) {
  digits <- x$digits
  if (is.null(digits) || !is.finite(digits)) digits <- 3L

  cat("FACETS-mode Workflow Summary\n")
  if (!is.null(x$overview) && nrow(x$overview) > 0) {
    ov <- round_numeric_frame(as.data.frame(x$overview), digits = digits)[1, , drop = FALSE]
    cat(sprintf("  Model: %s | Method: %s\n", ov$Model, ov$Method))
    cat(sprintf("  N: %s | Persons: %s | Facets: %s | Categories: %s\n", ov$N, ov$Persons, ov$Facets, ov$Categories))
    cat(sprintf("  LogLik: %s | AIC: %s | BIC: %s\n", ov$LogLik, ov$AIC, ov$BIC))
    cat(sprintf("  Converged: %s | Iterations: %s\n", ifelse(isTRUE(ov$Converged), "Yes", "No"), ov$Iterations))
  }

  if (!is.null(x$mapping) && nrow(x$mapping) > 0) {
    cat("\nColumn mapping\n")
    print(as.data.frame(x$mapping), row.names = FALSE)
  }

  cat("\nDetailed objects:\n")
  cat(" - summary(out$fit)\n")
  cat(" - summary(out$diagnostics)\n")

  invisible(x)
}

#' @export
print.mfrm_facets_run <- function(x, ...) {
  print(summary(x, ...), ...)
  invisible(x)
}

#' Plot outputs from a FACETS-mode workflow run
#'
#' @param x A `mfrm_facets_run` object from [run_mfrm_facets()].
#' @param y Unused.
#' @param type Plot route: `"fit"` delegates to [plot.mfrm_fit()] and `"qc"`
#'   delegates to [plot_qc_dashboard()].
#' @param ... Additional arguments passed to the selected plot function.
#'
#' @details
#' This method is a router for fast visualization from a one-shot workflow
#' result:
#' - `type = "fit"` for model-level displays.
#' - `type = "qc"` for multi-panel quality-control diagnostics.
#'
#' @section Interpreting output:
#' Returns the plotting object produced by the delegated route:
#' [plot.mfrm_fit()] for `"fit"` and [plot_qc_dashboard()] for `"qc"`.
#'
#' @section Typical workflow:
#' 1. Run [run_mfrm_facets()].
#' 2. Start with `plot(out, type = "fit", draw = FALSE)`.
#' 3. Continue with `plot(out, type = "qc", draw = FALSE)` for diagnostics.
#'
#' @return A plotting object from the delegated plot route.
#'
#' @seealso [run_mfrm_facets()], [plot.mfrm_fit()], [plot_qc_dashboard()]
#'
#' @examples
#' toy <- mfrmr:::sample_mfrm_data(seed = 123)
#' out <- run_mfrm_facets(
#'   data = toy,
#'   person = "Person",
#'   facets = c("Rater", "Task", "Criterion"),
#'   score = "Score",
#'   maxit = 15
#' )
#' p_fit <- plot(out, type = "fit", draw = FALSE)
#' class(p_fit)
#' p_qc <- plot(out, type = "qc", draw = FALSE)
#' class(p_qc)
#'
#' @export
plot.mfrm_facets_run <- function(x, y = NULL, type = c("fit", "qc"), ...) {
  if (!inherits(x, "mfrm_facets_run")) {
    stop("`x` must be an mfrm_facets_run object from run_mfrm_facets().", call. = FALSE)
  }
  type <- match.arg(type)
  if (identical(type, "fit")) {
    return(plot(x$fit, ...))
  }
  plot_qc_dashboard(fit = x$fit, diagnostics = x$diagnostics, ...)
}
