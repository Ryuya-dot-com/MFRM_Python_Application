#!/usr/bin/env Rscript

# One-case sirt worker. Native-code failures are isolated by the Python parent.

parse_args <- function(values) {
  output <- list()
  index <- 1L
  while (index <= length(values)) {
    key <- sub("^--", "", values[[index]])
    if (index == length(values)) stop("Missing value for --", key, call. = FALSE)
    output[[key]] <- values[[index + 1L]]
    index <- index + 2L
  }
  output
}

`%||%` <- function(value, replacement) {
  if (is.null(value) || length(value) == 0L) replacement else value
}

safe_max_abs <- function(values) {
  values <- as.numeric(values)
  if (!length(values) || !any(is.finite(values))) return(NA_real_)
  max(abs(values[is.finite(values)]))
}

prepare_sirt_wide <- function(data) {
  data$Item <- as.character(data$Unit)
  key <- paste(data$Person, data$Rater, data$Item, sep = "\r")
  if (anyDuplicated(key)) stop("Duplicate Person x Rater x response-unit cells.")
  grid <- unique(data[c("Person", "Rater")])
  grid <- grid[order(grid$Person, grid$Rater), , drop = FALSE]
  items <- sort(unique(as.character(data$Item)))
  response <- matrix(
    NA_real_, nrow = nrow(grid), ncol = length(items),
    dimnames = list(NULL, items)
  )
  row_index <- match(
    paste(data$Person, data$Rater), paste(grid$Person, grid$Rater)
  )
  column_index <- match(data$Item, items)
  response[cbind(row_index, column_index)] <- as.numeric(data$Score)
  list(response = response, grid = grid)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
for (name in c("input", "output", "case-id")) {
  if (!nzchar(args[[name]] %||% "")) stop("Missing --", name, call. = FALSE)
}
if (!requireNamespace("sirt", quietly = TRUE)) stop("sirt unavailable.")

manifest <- utils::read.csv(
  file.path(args$input, "manifest.csv"), stringsAsFactors = FALSE,
  check.names = FALSE
)
ratings <- utils::read.csv(
  file.path(args$input, "ratings.csv"), stringsAsFactors = FALSE,
  check.names = FALSE
)
python_runs <- utils::read.csv(
  file.path(args$input, "python_runs.csv"), stringsAsFactors = FALSE,
  check.names = FALSE
)
case_id <- as.character(args[["case-id"]])
manifest_row <- manifest[manifest$CaseId == case_id, , drop = FALSE]
python_row <- python_runs[python_runs$CaseId == case_id, , drop = FALSE]
data <- ratings[ratings$CaseId == case_id, , drop = FALSE]
if (nrow(manifest_row) != 1L || nrow(python_row) != 1L || !nrow(data)) {
  stop("Worker case identity missing or duplicated.")
}

row <- data.frame(
  CaseId = case_id,
  Family = as.character(manifest_row$Family),
  Mechanism = as.character(manifest_row$Mechanism),
  ExpectedStatus = as.character(manifest_row$ExpectedStatus),
  PythonWorkflowStatus = as.character(python_row$WorkflowStatus),
  Engine = "sirt",
  Estimator = "MML",
  Mode = "rm.facets_isolated_worker",
  FitAttempted = FALSE,
  FitReturned = FALSE,
  Converged = FALSE,
  Iterations = NA_real_,
  LogLik = NA_real_,
  MaxAbsEstimate = NA_real_,
  MaxAbsSE = NA_real_,
  NonfiniteEstimateCount = NA_real_,
  NonfiniteSECount = NA_real_,
  ObservedSupportMax = max(as.numeric(data$Score)),
  DeclaredRatingMax = as.integer(manifest_row$RatingMax),
  SupportMismatch = FALSE,
  Warnings = "",
  Error = "",
  ElapsedSeconds = 0,
  stringsAsFactors = FALSE
)
row$SupportMismatch <- row$ObservedSupportMax < row$DeclaredRatingMax
if (row$SupportMismatch) {
  row$Error <- paste(
    "category_support: sirt infers response-unit support from observed maxima;",
    "declared but unused top category not injected"
  )
  utils::write.csv(row, args$output, row.names = FALSE)
  quit(save = "no", status = 0L)
}

wide <- tryCatch(prepare_sirt_wide(data), error = function(error) error)
if (inherits(wide, "error")) {
  row$Error <- paste0("input: ", conditionMessage(wide))
  utils::write.csv(row, args$output, row.names = FALSE)
  quit(save = "no", status = 0L)
}

warnings <- character()
started <- proc.time()[["elapsed"]]
row$FitAttempted <- TRUE
fit <- tryCatch({
  invisible(utils::capture.output(
    value <- withCallingHandlers(
      sirt::rm.facets(
        dat = wide$response,
        pid = wide$grid$Person,
        rater = wide$grid$Rater,
        theta.k = seq(-6, 6, length.out = 41L),
        est.b.rater = TRUE,
        est.a.item = FALSE,
        est.a.rater = FALSE,
        rater_item_int = FALSE,
        est.mean = FALSE,
        b.rater.center = 2,
        maxdevchange = 0.01,
        globconv = 0.001,
        maxiter = 400L
      ),
      warning = function(warning) {
        warnings <<- c(warnings, conditionMessage(warning))
        invokeRestart("muffleWarning")
      }
    )
  ))
  value
}, error = function(error) error)
row$ElapsedSeconds <- proc.time()[["elapsed"]] - started
row$Warnings <- paste(unique(warnings), collapse = " | ")
if (inherits(fit, "error")) {
  row$Error <- paste0(class(fit)[[1L]], ": ", conditionMessage(fit))
  utils::write.csv(row, args$output, row.names = FALSE)
  quit(save = "no", status = 0L)
}

row$FitReturned <- TRUE
row$Iterations <- as.numeric(fit$iter %||% NA_real_)
row$Converged <- is.finite(row$Iterations) && row$Iterations < 400L
row$LogLik <- as.numeric(
  tryCatch(stats::logLik(fit), error = function(error) NA_real_)
)
estimates <- c(
  as.numeric(fit$b.rater %||% numeric()),
  as.numeric(fit$tau.item %||% numeric())
)
standard_errors <- c(
  as.numeric(fit$se.b.rater %||% numeric()),
  as.numeric(fit$se.tau.item %||% numeric())
)
row$MaxAbsEstimate <- safe_max_abs(estimates)
row$MaxAbsSE <- safe_max_abs(standard_errors)
row$NonfiniteEstimateCount <- sum(!is.finite(estimates))
row$NonfiniteSECount <- sum(!is.finite(standard_errors))
utils::write.csv(row, args$output, row.names = FALSE)
