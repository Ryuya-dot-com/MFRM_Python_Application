#!/usr/bin/env Rscript

# Independent base-R evaluator/refitter for the registered free-SD PCM MML lane.
# It intentionally does not source Python or use a Python-derived R package.

parse_args <- function(values) {
  result <- list()
  index <- 1L
  while (index <= length(values)) {
    key <- sub("^--", "", values[[index]])
    if (index == length(values)) stop("Missing value for --", key, call. = FALSE)
    result[[key]] <- values[[index + 1L]]
    index <- index + 2L
  }
  result
}

`%||%` <- function(value, replacement) {
  if (is.null(value) || !length(value)) replacement else value
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("study", "output", "vectors")
missing <- required[!vapply(required, function(key) nzchar(args[[key]] %||% ""), logical(1))]
if (length(missing)) stop("Missing arguments: ", paste(missing, collapse = ", "), call. = FALSE)

study_dir <- normalizePath(args$study, mustWork = TRUE)
output_dir <- normalizePath(args$output, mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
vector_ids <- as.integer(strsplit(args$vectors, ",", fixed = TRUE)[[1L]])
if (anyNA(vector_ids) || !length(vector_ids)) stop("--vectors must contain integers", call. = FALSE)
maxit <- as.integer(args$maxit %||% "400")
if (!is.finite(maxit) || maxit < 1L) stop("--maxit must be positive", call. = FALSE)

input_dir <- file.path(study_dir, "retained_input")
manifest <- utils::read.csv(file.path(input_dir, "manifest.csv"), stringsAsFactors = FALSE)
ratings_all <- utils::read.csv(file.path(input_dir, "generated_ratings.csv"), stringsAsFactors = FALSE)
attempts <- utils::read.csv(file.path(input_dir, "attempt_manifest.csv"), stringsAsFactors = FALSE)
attempts <- attempts[
  attempts$AttemptType == "PYTHON_MML_FREE_SD_Q31_PCM" &
    attempts$PersonVector %in% vector_ids,
  , drop = FALSE
]
if (nrow(attempts) != 3L * length(unique(vector_ids))) {
  stop("Expected three free-MML attempts per selected Person vector", call. = FALSE)
}

logsumexp <- function(values) {
  maximum <- max(values)
  maximum + log(sum(exp(values - maximum)))
}

normal_quadrature <- function(n, sigma) {
  # Golub-Welsch construction for physicists' Hermite quadrature.  The first
  # eigenvector row squared is already normalized to N(0,1) weights after the
  # sqrt(2) node transformation.
  off <- sqrt((1L:(n - 1L)) / 2)
  jacobi <- matrix(0, nrow = n, ncol = n)
  jacobi[cbind(1L:(n - 1L), 2L:n)] <- off
  jacobi[cbind(2L:n, 1L:(n - 1L))] <- off
  eig <- eigen(jacobi, symmetric = TRUE)
  order_index <- order(eig$values)
  list(
    nodes = sqrt(2) * eig$values[order_index] * sigma,
    weights = eig$vectors[1L, order_index]^2
  )
}

expand_parameters <- function(par) {
  if (length(par) != 12L) stop("Expected 12 coordinates", call. = FALSE)
  list(
    Rater = c(par[1:3], -sum(par[1:3])),
    Task = c(par[4:5], -sum(par[4:5])),
    Criterion = par[6:7],
    steps = rbind(
      C01 = c(par[8:9], -sum(par[8:9])),
      C02 = c(par[10:11], -sum(par[10:11]))
    ),
    sigma = exp(par[12])
  )
}

prepare_data <- function(data) {
  data$Person <- factor(data$Person, levels = sort(unique(data$Person)))
  data$Rater <- factor(data$Rater, levels = c("R01", "R02", "R03", "R04"))
  data$Task <- factor(data$Task, levels = c("T01", "T02", "T03"))
  data$Criterion <- factor(data$Criterion, levels = c("C01", "C02"))
  data$Score <- as.integer(data$Score)
  if (anyNA(data[, c("Person", "Rater", "Task", "Criterion", "Score")])) {
    stop("Cross-fit data contains unmapped values", call. = FALSE)
  }
  list(
    data = data,
    person = as.integer(data$Person),
    rater = as.integer(data$Rater),
    task = as.integer(data$Task),
    criterion = as.integer(data$Criterion),
    score = data$Score,
    n_person = nlevels(data$Person)
  )
}

marginal_components <- function(par, prepared, q = 31L, retain_posterior = FALSE) {
  expanded <- expand_parameters(par)
  quadrature <- normal_quadrature(q, expanded$sigma)
  log_weights <- log(quadrature$weights)
  base_eta <- -expanded$Rater[prepared$rater] -
    expanded$Task[prepared$task] -
    expanded$Criterion[prepared$criterion]
  cumulative <- t(apply(expanded$steps, 1L, function(value) c(0, cumsum(value))))
  threshold_by_row <- cumulative[prepared$criterion, , drop = FALSE]
  categories <- 0:3
  person_node_ll <- matrix(0, nrow = prepared$n_person, ncol = q)
  row_index <- seq_along(prepared$score)
  for (node_index in seq_len(q)) {
    eta <- quadrature$nodes[[node_index]] + base_eta
    logits <- outer(eta, categories) - threshold_by_row
    maximum <- pmax(logits[, 1L], logits[, 2L], logits[, 3L], logits[, 4L])
    denominator <- maximum + log(rowSums(exp(logits - maximum)))
    observed <- logits[cbind(row_index, prepared$score + 1L)]
    by_person <- rowsum(observed - denominator, prepared$person, reorder = FALSE)
    person_node_ll[, node_index] <- by_person[, 1L]
  }
  joint <- sweep(person_node_ll, 2L, log_weights, "+")
  maximum <- apply(joint, 1L, max)
  person_values <- maximum + log(rowSums(exp(joint - maximum)))
  posterior <- if (retain_posterior) exp(joint - person_values) else NULL
  list(
    loglik = sum(person_values),
    posterior = posterior,
    nodes = quadrature$nodes,
    weights = quadrature$weights,
    expanded = expanded
  )
}

objective <- function(par, prepared, q = 31L) {
  value <- tryCatch(-marginal_components(par, prepared, q = q)$loglik, error = function(error) Inf)
  if (is.finite(value)) value else .Machine$double.xmax / 100
}

finite_difference_gradient <- function(par, prepared, q = 31L, step = 1e-5) {
  gradient <- numeric(length(par))
  for (index in seq_along(par)) {
    local_step <- step * max(1, abs(par[[index]]))
    plus <- minus <- par
    plus[[index]] <- plus[[index]] + local_step
    minus[[index]] <- minus[[index]] - local_step
    gradient[[index]] <- (objective(plus, prepared, q) - objective(minus, prepared, q)) /
      (2 * local_step)
  }
  gradient
}

python_start <- function(recovery, thresholds, run) {
  facet <- function(name, levels) {
    rows <- recovery[recovery$Facet == name, , drop = FALSE]
    values <- rows$Estimate[match(levels, rows$Level)]
    if (anyNA(values)) stop("Missing Python facet coordinate: ", name, call. = FALSE)
    as.numeric(values)
  }
  rater <- facet("Rater", c("R01", "R02", "R03", "R04"))
  task <- facet("Task", c("T01", "T02", "T03"))
  criterion <- facet("Criterion", c("C01", "C02"))
  steps <- lapply(c("C01", "C02"), function(level) {
    rows <- thresholds[thresholds$StepFacetLevel == level, , drop = FALSE]
    rows <- rows[order(rows$Category), , drop = FALSE]
    as.numeric(rows$Estimate)
  })
  if (any(vapply(steps, length, integer(1)) != 3L)) stop("Missing Python PCM steps", call. = FALSE)
  sigma <- as.numeric(run$EstimatedPopulationSD[[1L]])
  if (!is.finite(sigma) || sigma <= 0) stop("Invalid Python population SD", call. = FALSE)
  c(rater[1:3], task[1:2], criterion, steps[[1L]][1:2], steps[[2L]][1:2], log(sigma))
}

run_rows <- list()
parameter_rows <- list()
gradient_rows <- list()
run_index <- 0L
parameter_index <- 0L
gradient_index <- 0L
coordinate_names <- c(
  "Rater::R01", "Rater::R02", "Rater::R03",
  "Task::T01", "Task::T02", "Criterion::C01", "Criterion::C02",
  "Step::C01::1", "Step::C01::2", "Step::C02::1", "Step::C02::2",
  "LogPopulationSD"
)

for (attempt_index in seq_len(nrow(attempts))) {
  attempt <- attempts[attempt_index, , drop = FALSE]
  run_id <- as.character(attempt$RunId[[1L]])
  ordinal <- as.integer(attempt$AttemptOrdinal[[1L]])
  artifact <- file.path(study_dir, "work", sprintf("%05d", ordinal))
  run <- utils::read.csv(file.path(artifact, "run_ledger.csv"), stringsAsFactors = FALSE)
  recovery <- utils::read.csv(file.path(artifact, "recovery.csv"), stringsAsFactors = FALSE)
  thresholds <- utils::read.csv(file.path(artifact, "thresholds.csv"), stringsAsFactors = FALSE)
  data <- ratings_all[ratings_all$RunId == run_id, c("Person", "Rater", "Task", "Criterion", "Score")]
  prepared <- prepare_data(data)
  start <- python_start(recovery, thresholds, run)
  python_q31 <- marginal_components(start, prepared, q = 31L, retain_posterior = TRUE)
  python_q61 <- marginal_components(start, prepared, q = 61L)
  python_gradient <- finite_difference_gradient(start, prepared, q = 31L)
  sigma_fixed_point <- sqrt(sum(python_q31$posterior * rep(python_q31$nodes^2, each = nrow(python_q31$posterior))) /
    nrow(python_q31$posterior))

  fit <- optim(
    par = start,
    fn = objective,
    gr = finite_difference_gradient,
    prepared = prepared,
    q = 31L,
    method = "BFGS",
    control = list(maxit = maxit, reltol = 1e-12)
  )
  fit_gradient <- finite_difference_gradient(fit$par, prepared, q = 31L)
  fit_q61 <- marginal_components(fit$par, prepared, q = 61L)
  fit61 <- optim(
    par = fit$par,
    fn = objective,
    gr = finite_difference_gradient,
    prepared = prepared,
    q = 61L,
    method = "BFGS",
    control = list(maxit = maxit, reltol = 1e-12)
  )
  fit61_gradient <- finite_difference_gradient(fit61$par, prepared, q = 61L)
  expanded_python <- expand_parameters(start)
  expanded_r <- expand_parameters(fit$par)
  expanded_r61 <- expand_parameters(fit61$par)

  run_index <- run_index + 1L
  run_rows[[run_index]] <- data.frame(
    RunId = run_id,
    PersonVector = as.integer(attempt$PersonVector[[1L]]),
    Gamma = as.numeric(attempt$Gamma[[1L]]),
    PythonRecordedLogLik = as.numeric(run$LogLik[[1L]]),
    PythonSolutionRLogLikQ31 = python_q31$loglik,
    PythonSolutionRLogLikQ61 = python_q61$loglik,
    PythonSolutionRGradientSupNormQ31 = max(abs(python_gradient)),
    PythonSigma = expanded_python$sigma,
    PythonSigmaFixedPoint = sigma_fixed_point,
    PythonSigmaFixedPointResidual = sigma_fixed_point - expanded_python$sigma,
    RConvergenceCode = as.integer(fit$convergence),
    RMessage = as.character(fit$message %||% ""),
    RLogLikQ31 = -as.numeric(fit$value),
    RQ31SolutionRLogLikQ61 = fit_q61$loglik,
    RGradientSupNormQ31 = max(abs(fit_gradient)),
    RSigma = expanded_r$sigma,
    RFunctionEvaluations = as.integer(fit$counts[["function"]]),
    RGradientEvaluations = as.integer(fit$counts[["gradient"]]),
    RQ61ConvergenceCode = as.integer(fit61$convergence),
    RQ61Message = as.character(fit61$message %||% ""),
    ROptimizedLogLikQ61 = -as.numeric(fit61$value),
    RGradientSupNormQ61 = max(abs(fit61_gradient)),
    RQ61Sigma = expanded_r61$sigma,
    RQ61FunctionEvaluations = as.integer(fit61$counts[["function"]]),
    RQ61GradientEvaluations = as.integer(fit61$counts[["gradient"]]),
    stringsAsFactors = FALSE
  )

  for (coordinate_index in seq_along(coordinate_names)) {
    gradient_index <- gradient_index + 1L
    gradient_rows[[gradient_index]] <- data.frame(
      RunId = run_id,
      Coordinate = coordinate_names[[coordinate_index]],
      PythonSolutionRGradientQ31 = python_gradient[[coordinate_index]],
      RQ31SolutionRGradientQ31 = fit_gradient[[coordinate_index]],
      RQ61SolutionRGradientQ61 = fit61_gradient[[coordinate_index]],
      stringsAsFactors = FALSE
    )
  }

  blocks <- list(
    Rater = c("R01", "R02", "R03", "R04"),
    Task = c("T01", "T02", "T03"),
    Criterion = c("C01", "C02")
  )
  for (block in names(blocks)) {
    levels <- blocks[[block]]
    for (level_index in seq_along(levels)) {
      parameter_index <- parameter_index + 1L
      parameter_rows[[parameter_index]] <- data.frame(
        RunId = run_id, Block = block, Level = levels[[level_index]],
        PythonEstimate = expanded_python[[block]][[level_index]],
        REstimateQ31 = expanded_r[[block]][[level_index]],
        REstimateQ61 = expanded_r61[[block]][[level_index]], stringsAsFactors = FALSE
      )
    }
  }
  for (criterion_index in seq_len(2L)) {
    for (step_index in seq_len(3L)) {
      parameter_index <- parameter_index + 1L
      parameter_rows[[parameter_index]] <- data.frame(
        RunId = run_id, Block = "Step",
        Level = paste0(c("C01", "C02")[[criterion_index]], "::", step_index),
        PythonEstimate = expanded_python$steps[criterion_index, step_index],
        REstimateQ31 = expanded_r$steps[criterion_index, step_index],
        REstimateQ61 = expanded_r61$steps[criterion_index, step_index], stringsAsFactors = FALSE
      )
    }
  }
}

runs <- do.call(rbind, run_rows)
parameters <- do.call(rbind, parameter_rows)
gradients <- do.call(rbind, gradient_rows)
parameters$DifferenceRQ31MinusPython <- parameters$REstimateQ31 - parameters$PythonEstimate
parameters$DifferenceRQ61MinusRQ31 <- parameters$REstimateQ61 - parameters$REstimateQ31
utils::write.csv(runs, file.path(output_dir, "crossfit_runs.csv"), row.names = FALSE)
utils::write.csv(parameters, file.path(output_dir, "crossfit_parameters.csv"), row.names = FALSE)
utils::write.csv(gradients, file.path(output_dir, "crossfit_gradients.csv"), row.names = FALSE)
software_versions <- extSoftVersion()
software_value <- function(name) {
  if (name %in% names(software_versions)) software_versions[[name]] else NA_character_
}
utils::write.csv(
  data.frame(
    RVersion = as.character(getRversion()),
    Platform = R.version$platform,
    BLAS = software_value("BLAS"),
    LAPACK = software_value("LAPACK"),
    SelectedVectors = paste(vector_ids, collapse = ","),
    Datasets = nrow(runs),
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, "runtime_identity.csv"),
  row.names = FALSE
)
