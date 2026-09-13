#!/usr/bin/env Rscript

# Cross-engine companion runner for validation/cross_engine_stress.py.
# Fits the exact generated rows with mfrmr 0.2.3, TAM, immer, and sirt.

parse_args <- function(x) {
  out <- list()
  i <- 1L
  while (i <= length(x)) {
    key <- sub("^--", "", x[[i]])
    if (i == length(x)) stop("Missing value for --", key, call. = FALSE)
    out[[key]] <- x[[i + 1L]]
    i <- i + 2L
  }
  out
}

`%||%` <- function(value, replacement) {
  if (is.null(value) || length(value) == 0L) replacement else value
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("input", "output", "mfrmr-lib")
missing_args <- required[!vapply(required, function(x) nzchar(args[[x]] %||% ""), logical(1))]
if (length(missing_args) > 0L) {
  stop("Missing required arguments: ", paste(missing_args, collapse = ", "), call. = FALSE)
}

.libPaths(c(args[["mfrmr-lib"]], .libPaths()))
for (pkg in c("mfrmr", "TAM", "immer", "sirt", "digest")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required package unavailable: ", pkg, call. = FALSE)
  }
}
suppressPackageStartupMessages(library(mfrmr))

input_dir <- normalizePath(args$input, mustWork = TRUE)
output_dir <- normalizePath(args$output, mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

manifest <- utils::read.csv(file.path(input_dir, "manifest.csv"), stringsAsFactors = FALSE)
ratings <- utils::read.csv(file.path(input_dir, "simulated_ratings.csv"), stringsAsFactors = FALSE)

as_flag <- function(x) isTRUE(as.logical(x))
as_scalar <- function(x, default = NA_real_) {
  if (is.null(x) || length(x) == 0L) return(default)
  value <- suppressWarnings(as.numeric(x[[1L]]))
  if (is.finite(value)) value else default
}

run_rows <- list()
surface_rows <- list()
person_rows <- list()
run_index <- 0L
surface_index <- 0L
person_index <- 0L

append_run <- function(row) {
  run_index <<- run_index + 1L
  run_rows[[run_index]] <<- row
  invisible(NULL)
}

append_surface <- function(surface) {
  if (is.null(surface) || nrow(surface) == 0L) return(invisible(NULL))
  surface_index <<- surface_index + 1L
  surface_rows[[surface_index]] <<- surface
  invisible(NULL)
}

append_person <- function(person) {
  if (is.null(person) || nrow(person) == 0L) return(invisible(NULL))
  person_index <<- person_index + 1L
  person_rows[[person_index]] <<- person
  invisible(NULL)
}

formula_for <- function(model) {
  if (identical(model, "RSM")) {
    stats::as.formula("~ item + rater + step")
  } else {
    stats::as.formula("~ item + rater + item:step")
  }
}

prepare_wide <- function(data) {
  criteria <- sort(unique(as.character(data$Criterion)))
  grid <- unique(data[c("Person", "Rater")])
  grid <- grid[order(grid$Person, grid$Rater), , drop = FALSE]
  key <- paste(data$Person, data$Rater, data$Criterion, sep = "\r")
  response <- vapply(criteria, function(criterion) {
    target <- paste(grid$Person, grid$Rater, criterion, sep = "\r")
    index <- match(target, key)
    value <- rep(NA_integer_, nrow(grid))
    retained <- !is.na(index)
    value[retained] <- as.integer(data$Score[index[retained]])
    value
  }, integer(nrow(grid)))
  response <- as.data.frame(response, stringsAsFactors = FALSE)
  names(response) <- criteria
  retained <- rowSums(!is.na(response)) > 0L
  list(grid = grid[retained, , drop = FALSE], resp = response[retained, , drop = FALSE])
}

prepare_jml_design <- function(data, model) {
  wide <- prepare_wide(data)
  output <- utils::capture.output(
    design <- suppressWarnings(TAM::tam.mml.mfr(
      resp = wide$resp,
      facets = data.frame(rater = wide$grid$Rater, stringsAsFactors = FALSE),
      pid = wide$grid$Person,
      formulaA = formula_for(model),
      constraint = "cases",
      control = list(maxiter = 2L, progress = FALSE),
      verbose = FALSE
    ))
  )
  item <- rownames(design$A)
  item_map <- data.frame(
    Item = item,
    Criterion = sub("-rater.*$", "", item),
    Rater = sub("^.*-rater", "", item),
    stringsAsFactors = FALSE
  )
  list(
    resp = design$resp,
    pid = as.character(design$pid),
    A_tam = design$A,
    A_immer = design$A[, -1L, , drop = FALSE],
    item_map = item_map,
    output = output
  )
}

surface_from_matrix <- function(matrix, item_map, dataset_id, engine, mode) {
  matrix <- as.matrix(matrix)
  if (is.null(rownames(matrix))) rownames(matrix) <- item_map$Item
  index <- match(item_map$Item, rownames(matrix))
  if (anyNA(index)) stop("Cumulative-difficulty rows do not align.", call. = FALSE)
  matrix <- matrix[index, , drop = FALSE]
  rows <- lapply(seq_len(ncol(matrix)), function(k) {
    data.frame(
      DatasetId = dataset_id,
      Engine = engine,
      Mode = mode,
      Rater = item_map$Rater,
      Criterion = item_map$Criterion,
      Category = k,
      Estimate = as.numeric(matrix[, k]),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

step_values <- function(table, criterion, model) {
  table <- as.data.frame(table, stringsAsFactors = FALSE)
  if (identical(model, "PCM")) {
    table <- table[as.character(table$StepFacet) == criterion, , drop = FALSE]
  }
  index <- suppressWarnings(as.integer(sub(".*_", "", table$Step)))
  as.numeric(table$Estimate[order(index)])
}

surface_from_mfrmr <- function(fit, dataset_id, mode, model) {
  facets <- as.data.frame(fit$facets$others, stringsAsFactors = FALSE)
  rater <- facets[facets$Facet == "Rater", , drop = FALSE]
  criterion <- facets[facets$Facet == "Criterion", , drop = FALSE]
  rater_values <- stats::setNames(as.numeric(rater$Estimate), as.character(rater$Level))
  criterion_values <- stats::setNames(as.numeric(criterion$Estimate), as.character(criterion$Level))
  rows <- list()
  index <- 0L
  for (rr in names(rater_values)) {
    for (cc in names(criterion_values)) {
      steps <- step_values(fit$steps, cc, model)
      for (k in seq_along(steps)) {
        index <- index + 1L
        rows[[index]] <- data.frame(
          DatasetId = dataset_id,
          Engine = "mfrmr",
          Mode = mode,
          Rater = rr,
          Criterion = cc,
          Category = k,
          Estimate = k * (rater_values[[rr]] + criterion_values[[cc]]) + sum(steps[seq_len(k)]),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, rows)
}

persons_from_mfrmr <- function(fit, dataset_id, mode) {
  table <- as.data.frame(fit$facets$person, stringsAsFactors = FALSE)
  estimate <- if ("OptimizerEstimate" %in% names(table)) table$OptimizerEstimate else table$Estimate
  if (!any(is.finite(suppressWarnings(as.numeric(estimate))))) estimate <- table$Estimate
  status <- if ("ParameterStatus" %in% names(table)) table$ParameterStatus else "not_reported"
  data.frame(
    DatasetId = dataset_id,
    Engine = "mfrmr",
    Mode = mode,
    Person = as.character(table$Person),
    Estimate = as.numeric(estimate),
    ParameterStatus = as.character(status),
    stringsAsFactors = FALSE
  )
}

base_run <- function(manifest_row, engine, mode, estimator, quadrature = NA_real_,
                     population_sd_mode = "not_applicable") {
  data.frame(
    DatasetId = as.character(manifest_row$DatasetId),
    Engine = engine,
    Mode = mode,
    Estimator = estimator,
    Model = as.character(manifest_row$Model),
    Scenario = as.character(manifest_row$Scenario),
    Replicate = as.integer(manifest_row$Replicate),
    Quadrature = as.numeric(quadrature),
    PopulationSDMode = population_sd_mode,
    PopulationSD = NA_real_,
    FitReturned = FALSE,
    Converged = FALSE,
    ConvergenceBasis = "exception",
    InferenceReady = FALSE,
    LogLik = NA_real_,
    Npar = NA_real_,
    Iterations = NA_real_,
    GradientNorm = NA_real_,
    ElapsedSeconds = NA_real_,
    Error = "",
    Warnings = "",
    ExpectedConnected = as_flag(manifest_row$ExpectedConnected),
    stringsAsFactors = FALSE
  )
}

record_failure <- function(manifest_row, engine, mode, estimator, elapsed, error,
                           quadrature = NA_real_, population_sd_mode = "not_applicable") {
  row <- base_run(manifest_row, engine, mode, estimator, quadrature, population_sd_mode)
  row$ElapsedSeconds <- elapsed
  row$Error <- conditionMessage(error)
  append_run(row)
}

fit_mfrmr_mode <- function(data, manifest_row, mode, estimator, q = NA_integer_, free_sd = FALSE) {
  started <- proc.time()[["elapsed"]]
  value <- tryCatch({
    call <- list(
      data = data,
      person = "Person",
      facets = c("Rater", "Criterion"),
      score = "Score",
      rating_min = 0,
      rating_max = 3,
      model = as.character(manifest_row$Model),
      method = estimator,
      maxit = if (identical(estimator, "JML")) 500L else 600L,
      reltol = if (identical(estimator, "JML")) 1e-8 else 1e-7,
      optimizer = "BFGS",
      min_obs_per_element = 1L,
      min_obs_per_category = 1L
    )
    if (identical(as.character(manifest_row$Model), "PCM")) call$step_facet <- "Criterion"
    if (identical(estimator, "MML")) {
      call$quad_points <- as.integer(q)
      call$mml_engine <- "direct"
      if (isTRUE(free_sd)) {
        call$population_formula <- stats::as.formula("~ 1")
        call$person_data <- data.frame(Person = sort(unique(as.character(data$Person))))
        call$person_id <- "Person"
      }
    }
    suppressMessages(suppressWarnings(do.call(mfrmr::fit_mfrm, call)))
  }, error = function(e) e)
  elapsed <- proc.time()[["elapsed"]] - started
  if (inherits(value, "error")) {
    record_failure(
      manifest_row, "mfrmr", mode, estimator, elapsed, value,
      if (identical(estimator, "MML")) q else NA_real_,
      if (identical(estimator, "MML")) if (free_sd) "estimated" else "fixed" else "not_applicable"
    )
    return(invisible(NULL))
  }
  fit <- value
  summary <- as.data.frame(fit$summary, stringsAsFactors = FALSE)[1L, , drop = FALSE]
  run <- base_run(
    manifest_row, "mfrmr", mode, estimator,
    if (identical(estimator, "MML")) q else NA_real_,
    if (identical(estimator, "MML")) if (free_sd) "estimated" else "fixed" else "not_applicable"
  )
  run$PopulationSD <- if (identical(estimator, "MML")) {
    if (free_sd) sqrt(as_scalar(fit$population$sigma2)) else 1
  } else NA_real_
  run$FitReturned <- TRUE
  run$Converged <- isTRUE(summary$Converged)
  run$ConvergenceBasis <- as.character(summary$ConvergenceBasis %||% "mfrmr_summary")
  run$InferenceReady <- isTRUE(summary$InferenceReady)
  run$LogLik <- as_scalar(summary$LogLik)
  run$Npar <- as_scalar(summary$Npar)
  run$Iterations <- as_scalar(summary$Iterations)
  run$GradientNorm <- as_scalar(summary$TerminalGradientSupNorm)
  run$ElapsedSeconds <- elapsed
  append_run(run)
  append_surface(surface_from_mfrmr(fit, as.character(manifest_row$DatasetId), mode, as.character(manifest_row$Model)))
  append_person(persons_from_mfrmr(fit, as.character(manifest_row$DatasetId), mode))
  invisible(NULL)
}

fit_tam_jml <- function(prepared, manifest_row, mode, adj) {
  started <- proc.time()[["elapsed"]]
  value <- tryCatch({
    output <- utils::capture.output(
      fit <- suppressWarnings(TAM::tam.jml(
        resp = prepared$resp,
        A = prepared$A_tam,
        adj = adj,
        bias = FALSE,
        constraint = "cases",
        verbose = FALSE,
        control = list(maxiter = 600L, Msteps = 10L, conv = 1e-8, progress = FALSE)
      ))
    )
    fit
  }, error = function(e) e)
  elapsed <- proc.time()[["elapsed"]] - started
  if (inherits(value, "error")) {
    record_failure(manifest_row, "TAM", mode, "JML", elapsed, value)
    return(invisible(NULL))
  }
  fit <- value
  run <- base_run(manifest_row, "TAM", mode, "JML")
  run$FitReturned <- TRUE
  run$Converged <- as_scalar(fit$iter) < 600
  run$ConvergenceBasis <- "iteration_cap_proxy"
  run$InferenceReady <- NA
  run$LogLik <- -0.5 * as_scalar(fit$deviance)
  run$Npar <- length(fit$xsi) + length(fit$theta)
  run$Iterations <- as_scalar(fit$iter)
  run$ElapsedSeconds <- elapsed
  append_run(run)
  append_surface(surface_from_matrix(
    -as.matrix(fit$AXsi[, -1L, drop = FALSE]), prepared$item_map,
    as.character(manifest_row$DatasetId), "TAM", mode
  ))
  append_person(data.frame(
    DatasetId = as.character(manifest_row$DatasetId), Engine = "TAM", Mode = mode,
    Person = prepared$pid, Estimate = as.numeric(fit$theta),
    ParameterStatus = ifelse(is.finite(fit$theta), "finite_trace", "nonfinite"),
    stringsAsFactors = FALSE
  ))
  invisible(NULL)
}

fit_immer_jml <- function(prepared, manifest_row, mode, method) {
  started <- proc.time()[["elapsed"]]
  value <- tryCatch(suppressWarnings(immer::immer_jml(
    dat = prepared$resp,
    A = prepared$A_immer,
    est_method = method,
    eps = 0.3,
    center_theta = TRUE,
    maxiter = 1000L,
    conv = 1e-8,
    verbose = FALSE,
    use_Rcpp = TRUE,
    shortcut = TRUE
  )), error = function(e) e)
  elapsed <- proc.time()[["elapsed"]] - started
  if (inherits(value, "error")) {
    record_failure(manifest_row, "immer", mode, "JML", elapsed, value)
    return(invisible(NULL))
  }
  fit <- value
  run <- base_run(manifest_row, "immer", mode, "JML")
  run$FitReturned <- TRUE
  run$Converged <- as_scalar(fit$iter) < 1000
  run$ConvergenceBasis <- "iteration_cap_proxy"
  run$InferenceReady <- NA
  run$LogLik <- as_scalar(fit$loglike)
  run$Npar <- length(fit$xsi) + length(fit$theta)
  run$Iterations <- as_scalar(fit$iter)
  run$ElapsedSeconds <- elapsed
  append_run(run)
  append_surface(surface_from_matrix(
    as.matrix(fit$b), prepared$item_map,
    as.character(manifest_row$DatasetId), "immer", mode
  ))
  append_person(data.frame(
    DatasetId = as.character(manifest_row$DatasetId), Engine = "immer", Mode = mode,
    Person = prepared$pid, Estimate = as.numeric(fit$theta),
    ParameterStatus = ifelse(is.finite(fit$theta), "finite_trace", "nonfinite"),
    stringsAsFactors = FALSE
  ))
  invisible(NULL)
}

fit_tam_mml <- function(data, manifest_row, q) {
  mode <- paste0("TAM_MML_Q", q)
  wide <- prepare_wide(data)
  started <- proc.time()[["elapsed"]]
  value <- tryCatch({
    output <- utils::capture.output(
      fit <- suppressWarnings(TAM::tam.mml.mfr(
        resp = wide$resp,
        facets = data.frame(rater = wide$grid$Rater, stringsAsFactors = FALSE),
        pid = wide$grid$Person,
        formulaA = formula_for(as.character(manifest_row$Model)),
        constraint = "cases",
        control = list(
          nodes = seq(-8, 8, length.out = as.integer(q)),
          maxiter = 1000L, progress = FALSE
        ),
        verbose = FALSE
      ))
    )
    fit
  }, error = function(e) e)
  elapsed <- proc.time()[["elapsed"]] - started
  if (inherits(value, "error")) {
    record_failure(manifest_row, "TAM", mode, "MML", elapsed, value, q, "estimated")
    return(invisible(NULL))
  }
  fit <- value
  item <- rownames(fit$A)
  item_map <- data.frame(
    Item = item,
    Criterion = sub("-rater.*$", "", item),
    Rater = sub("^.*-rater", "", item),
    stringsAsFactors = FALSE
  )
  run <- base_run(manifest_row, "TAM", mode, "MML", q, "estimated")
  run$PopulationSD <- sqrt(as_scalar(fit$variance[1L, 1L]))
  run$FitReturned <- TRUE
  run$Converged <- as_scalar(fit$iter) < 1000
  run$ConvergenceBasis <- "iteration_cap_proxy"
  run$InferenceReady <- NA
  run$LogLik <- -0.5 * as_scalar(fit$deviance)
  run$Npar <- as_scalar(fit$ic$np)
  run$Iterations <- as_scalar(fit$iter)
  run$ElapsedSeconds <- elapsed
  append_run(run)
  append_surface(surface_from_matrix(
    -as.matrix(fit$AXsi[, -1L, drop = FALSE]), item_map,
    as.character(manifest_row$DatasetId), "TAM", mode
  ))
  append_person(data.frame(
    DatasetId = as.character(manifest_row$DatasetId), Engine = "TAM", Mode = mode,
    Person = as.character(fit$person$pid), Estimate = as.numeric(fit$person$EAP),
    ParameterStatus = "eap", stringsAsFactors = FALSE
  ))
  invisible(NULL)
}

fit_sirt_mml <- function(data, manifest_row, q) {
  mode <- paste0("SIRT_MML_Q", q)
  wide <- prepare_wide(data)
  maxiter <- 1200L
  started <- proc.time()[["elapsed"]]
  value <- tryCatch({
    output <- utils::capture.output(
      fit <- suppressWarnings(sirt::rm.facets(
        dat = wide$resp,
        pid = wide$grid$Person,
        rater = wide$grid$Rater,
        theta.k = seq(-8, 8, length.out = as.integer(q)),
        est.b.rater = TRUE,
        est.a.item = FALSE,
        est.a.rater = FALSE,
        est.mean = FALSE,
        globconv = 1e-4,
        maxiter = maxiter
      ))
    )
    fit
  }, error = function(e) e)
  elapsed <- proc.time()[["elapsed"]] - started
  if (inherits(value, "error")) {
    record_failure(manifest_row, "sirt", mode, "MML", elapsed, value, q, "estimated")
    return(invisible(NULL))
  }
  fit <- value
  raters <- sort(unique(as.character(wide$grid$Rater)))
  criteria <- colnames(wide$resp)
  b <- stats::setNames(as.numeric(fit$b.rater), raters)
  tau <- as.matrix(fit$tau.item)
  rownames(tau) <- criteria
  rows <- list()
  index <- 0L
  for (rr in raters) {
    for (cc in criteria) {
      for (k in seq_len(ncol(tau))) {
        index <- index + 1L
        rows[[index]] <- data.frame(
          DatasetId = as.character(manifest_row$DatasetId),
          Engine = "sirt", Mode = mode, Rater = rr, Criterion = cc,
          Category = k, Estimate = as.numeric(tau[cc, k]) + k * b[[rr]],
          stringsAsFactors = FALSE
        )
      }
    }
  }
  run <- base_run(manifest_row, "sirt", mode, "MML", q, "estimated")
  run$PopulationSD <- as_scalar(fit$sigma)
  run$FitReturned <- TRUE
  run$Converged <- as_scalar(fit$iter) < maxiter
  run$ConvergenceBasis <- "iteration_cap_proxy"
  run$InferenceReady <- NA
  run$LogLik <- as_scalar(stats::logLik(fit))
  run$Npar <- as_scalar(fit$ic$np)
  run$Iterations <- as_scalar(fit$iter)
  run$ElapsedSeconds <- elapsed
  append_run(run)
  append_surface(do.call(rbind, rows))
  append_person(data.frame(
    DatasetId = as.character(manifest_row$DatasetId), Engine = "sirt", Mode = mode,
    Person = as.character(fit$person$pid), Estimate = as.numeric(fit$person$EAP),
    ParameterStatus = "eap", stringsAsFactors = FALSE
  ))
  invisible(NULL)
}

for (i in seq_len(nrow(manifest))) {
  row <- manifest[i, , drop = FALSE]
  data <- ratings[ratings$DatasetId == row$DatasetId, c("Person", "Rater", "Criterion", "Score")]
  if (as_flag(row$RunJML)) {
    fit_mfrmr_mode(data, row, "MFRMR_JML_RAW", "JML")
    prepared <- tryCatch(prepare_jml_design(data, as.character(row$Model)), error = function(e) e)
    if (inherits(prepared, "error")) {
      for (spec in list(
        c("TAM", "TAM_JML_RAW"), c("TAM", "TAM_JML_ADJ"),
        c("immer", "IMMER_JML_PERSON_EPS"), c("immer", "IMMER_EPS_ADJ")
      )) {
        record_failure(row, spec[[1L]], spec[[2L]], "JML", 0, prepared)
      }
    } else {
      fit_tam_jml(prepared, row, "TAM_JML_RAW", 0)
      fit_tam_jml(prepared, row, "TAM_JML_ADJ", 0.3)
      fit_immer_jml(prepared, row, "IMMER_JML_PERSON_EPS", "jml")
      fit_immer_jml(prepared, row, "IMMER_EPS_ADJ", "eps_adj")
    }
  }
  if (as_flag(row$RunMML)) {
    fit_mfrmr_mode(data, row, "MFRMR_MML_FREE_Q31", "MML", 31L, TRUE)
    fit_tam_mml(data, row, 21L)
    if (as_flag(row$RunSirt) && identical(as.character(row$Model), "PCM")) {
      fit_sirt_mml(data, row, 30L)
    }
    if (as_flag(row$RunQ61)) {
      fit_mfrmr_mode(data, row, "MFRMR_MML_FREE_Q61", "MML", 61L, TRUE)
      fit_tam_mml(data, row, 61L)
      if (as_flag(row$RunSirt) && identical(as.character(row$Model), "PCM")) {
        fit_sirt_mml(data, row, 61L)
      }
    }
    if (identical(as.character(row$Scenario), "balanced")) {
      fit_mfrmr_mode(data, row, "MFRMR_MML_FIXED_Q31_PERSON", "MML", 31L, FALSE)
      fit_mfrmr_mode(data, row, "MFRMR_MML_FIXED_Q15_PERSON", "MML", 15L, FALSE)
    }
  }
}

runs <- if (length(run_rows)) do.call(rbind, run_rows) else data.frame()
surfaces <- if (length(surface_rows)) do.call(rbind, surface_rows) else data.frame()
persons <- if (length(person_rows)) do.call(rbind, person_rows) else data.frame()
utils::write.csv(runs, file.path(output_dir, "fit_runs_r.csv"), row.names = FALSE)
utils::write.csv(surfaces, file.path(output_dir, "surfaces_r.csv"), row.names = FALSE)
utils::write.csv(persons, file.path(output_dir, "persons_r.csv"), row.names = FALSE)

function_hash <- function(package, name) {
  fn <- get(name, envir = asNamespace(package), inherits = FALSE)
  digest::digest(list(formals = formals(fn), body = body(fn)), algo = "sha256", serialize = TRUE)
}

identity <- data.frame(
  Engine = c("mfrmr", "TAM", "immer", "sirt"),
  Version = vapply(c("mfrmr", "TAM", "immer", "sirt"), function(pkg) as.character(utils::packageVersion(pkg)), character(1)),
  PrimaryFunction = c("fit_mfrm", "tam.mml.mfr/tam.jml", "immer_jml", "rm.facets"),
  FunctionSHA256 = c(
    function_hash("mfrmr", "fit_mfrm"),
    digest::digest(c(function_hash("TAM", "tam.mml.mfr"), function_hash("TAM", "tam.jml")), algo = "sha256"),
    function_hash("immer", "immer_jml"),
    function_hash("sirt", "rm.facets")
  ),
  SourceGitHead = c(args[["mfrmr-git-head"]] %||% NA_character_, NA, NA, NA),
  SourceState = c(args[["mfrmr-source-state"]] %||% "development snapshot", "CRAN installed", "CRAN installed", "CRAN installed"),
  stringsAsFactors = FALSE
)
utils::write.csv(identity, file.path(output_dir, "r_runtime_identity.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "r_session_info.txt"))
