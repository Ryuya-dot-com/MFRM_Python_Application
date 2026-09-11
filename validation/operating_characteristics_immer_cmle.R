#!/usr/bin/env Rscript

# Fit immer::immer_cml() to the same byte-validated OC rating rows used by the
# native Python CMLE adapter.  Python exact moments and an independently
# reconstructed R/immer conditional-information matrix must agree at pre-fit;
# structurally deficient runs remain in the ledger but reach neither optimizer.

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
required <- c("input", "output")
missing <- required[!vapply(required, function(name) nzchar(args[[name]] %||% ""), logical(1))]
if (length(missing)) stop("Missing required arguments: ", paste(missing, collapse = ", "), call. = FALSE)
for (package in c("immer", "digest", "psychotools")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Required adapter dependency unavailable: ", package, call. = FALSE)
  }
}

input_dir <- normalizePath(args$input, mustWork = TRUE)
output_dir <- normalizePath(args$output, mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
read_input <- function(filename) {
  path <- file.path(input_dir, filename)
  if (!file.exists(path)) stop("Missing adapter input: ", filename, call. = FALSE)
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

manifest <- read_input("manifest.csv")
ratings <- read_input("generated_ratings.csv")
anchors <- read_input("generated_anchors.csv")
inventory <- read_input("generated_bundle_files.csv")
bridge_validation <- read_input("bridge_validation_r.csv")
bridge_engines <- read_input("bridge_engine_availability_r.csv")
python_runs <- read_input("python_cmle_runs.csv")

bundle_inventory_sha256 <- digest::digest(
  file = file.path(input_dir, "generated_bundle_files.csv"),
  algo = "sha256",
  serialize = FALSE
)
if (!all(as.logical(bridge_validation$Passed)) ||
    !all(as.character(bridge_validation$BundleInventorySHA256) == bundle_inventory_sha256)) {
  stop("R bridge validation failed or is stale for the current bundle.", call. = FALSE)
}
hash_ok <- vapply(seq_len(nrow(inventory)), function(i) {
  path <- file.path(input_dir, as.character(inventory$File[[i]]))
  file.exists(path) && identical(
    digest::digest(file = path, algo = "sha256", serialize = FALSE),
    as.character(inventory$SHA256[[i]])
  )
}, logical(1))
if (!all(hash_ok)) stop("Generated bundle byte hash changed after bridge validation.", call. = FALSE)
if (anyDuplicated(manifest$RunId) || anyDuplicated(python_runs$RunId) ||
    !setequal(as.character(manifest$RunId), as.character(python_runs$RunId))) {
  stop("Python CMLE run ledger and manifest identities do not match.", call. = FALSE)
}
if (!all(as.character(python_runs$BundleInventorySHA256) == bundle_inventory_sha256)) {
  stop("Python CMLE adapter is stale for the current generated bundle.", call. = FALSE)
}

function_hash <- function(package, functions) {
  namespace <- asNamespace(package)
  hashes <- vapply(functions, function(name) {
    fn <- get(name, envir = namespace, inherits = FALSE)
    digest::digest(list(formals = formals(fn), body = body(fn)), algo = "sha256", serialize = TRUE)
  }, character(1))
  digest::digest(hashes, algo = "sha256", serialize = TRUE)
}
immer_function_hash <- function_hash("immer", c("immer_jml", "immer_cml"))
registered <- bridge_engines[bridge_engines$Engine == "immer", , drop = FALSE]
if (nrow(registered) != 1L || !isTRUE(as.logical(registered$Available[[1L]])) ||
    !identical(as.character(registered$Version[[1L]]), as.character(utils::packageVersion("immer"))) ||
    !identical(as.character(registered$FunctionSHA256[[1L]]), immer_function_hash)) {
  stop("Loaded immer code identity differs from the validated R bridge.", call. = FALSE)
}

sum_zero <- function(n) {
  if (n <= 1L) return(matrix(0, nrow = n, ncol = 0L))
  out <- matrix(0, nrow = n, ncol = n - 1L)
  out[seq_len(n - 1L), ] <- diag(n - 1L)
  out[n, ] <- -1
  out
}

make_wide <- function(data) {
  data$VirtualUnit <- paste(data$Rater, data$Task, data$Criterion, sep = "__")
  persons <- sort(unique(as.character(data$Person)))
  units <- sort(unique(as.character(data$VirtualUnit)))
  response <- matrix(NA_real_, nrow = length(persons), ncol = length(units),
                     dimnames = list(persons, units))
  response[cbind(match(data$Person, persons), match(data$VirtualUnit, units))] <- as.numeric(data$Score)
  mapping <- unique(data[c("VirtualUnit", "Rater", "Task", "Criterion")])
  mapping <- mapping[match(units, mapping$VirtualUnit), , drop = FALSE]
  list(response = response, mapping = mapping)
}

build_design <- function(response, mapping, rating_max) {
  # immer infers each item's support from its observed maximum.  A zero-weight
  # complete sentinel declares the Python 0..rating_max support without adding
  # any observed sufficient statistic or score frequency.
  augmented <- rbind(response, rep(as.numeric(rating_max), ncol(response)))
  rownames(augmented)[nrow(augmented)] <- "__declared_support_zero_weight__"
  weights <- c(rep(1, nrow(response)), 0)
  prep <- immer:::lpcm_data_prep(augmented, weights = weights, a = NULL)
  if (!all(prep$maxK == rating_max)) stop("Declared support sentinel did not reach every virtual unit.")
  pars_info <- prep$pars_info
  facet_names <- c("Rater", "Task", "Criterion")
  facet_levels <- lapply(facet_names, function(name) sort(unique(as.character(mapping[[name]]))))
  names(facet_levels) <- facet_names
  contrasts <- lapply(facet_levels, function(levels) sum_zero(length(levels)))
  step_contrast <- sum_zero(as.integer(rating_max))
  parameter_names <- unlist(c(
    lapply(facet_names, function(name) {
      contrast <- contrasts[[name]]
      if (!ncol(contrast)) return(character())
      paste0("facet:", name, ":free:", facet_levels[[name]][seq_len(ncol(contrast))])
    }),
    list(if (ncol(step_contrast)) paste0("step:__shared__:free:", seq_len(ncol(step_contrast))) else character())
  ), use.names = FALSE)
  W <- matrix(0, nrow = nrow(pars_info), ncol = length(parameter_names),
              dimnames = list(rownames(pars_info), parameter_names))
  offsets <- integer(length(facet_names))
  cursor <- 0L
  for (j in seq_along(facet_names)) {
    offsets[[j]] <- cursor
    cursor <- cursor + ncol(contrasts[[facet_names[[j]]]])
  }
  step_offset <- cursor
  for (row_index in seq_len(nrow(pars_info))) {
    unit <- as.character(pars_info$item[[row_index]])
    category <- as.integer(pars_info$cat[[row_index]])
    meta <- mapping[mapping$VirtualUnit == unit, , drop = FALSE]
    if (nrow(meta) != 1L) stop("Virtual-unit mapping is missing or duplicated: ", unit)
    for (j in seq_along(facet_names)) {
      name <- facet_names[[j]]
      contrast <- contrasts[[name]]
      if (ncol(contrast)) {
        level_index <- match(as.character(meta[[name]][[1L]]), facet_levels[[name]])
        cols <- offsets[[j]] + seq_len(ncol(contrast))
        W[row_index, cols] <- category * contrast[level_index, ]
      }
    }
    if (ncol(step_contrast)) {
      cols <- step_offset + seq_len(ncol(step_contrast))
      W[row_index, cols] <- colSums(step_contrast[seq_len(category), , drop = FALSE])
    }
  }
  list(
    response = augmented,
    weights = weights,
    W = W,
    parameters = parameter_names,
    prepared = prep,
    b_const = rep(0, nrow(W))
  )
}

conditional_gradient <- function(object, par) {
  esf_par0 <- as.numeric(object$W %*% as.numeric(par) + object$b_const)
  pieces <- lapply(seq_len(object$NP), function(pp) {
    # Recreate immer's grouping exactly from its prepared item ids.
    values <- split(
      esf_par0[object$parm_index[[pp]]],
      object$pars_info$itemid[object$parm_index[[pp]]]
    )
    esf <- psychotools::elementary_symmetric_functions(par = values, order = 1, diff = FALSE)
    gamma0 <- esf[[1L]]
    gamma1 <- esf[[2L]]
    W1 <- object$W[object$parm_index[[pp]], , drop = FALSE]
    observed <- object$suffstat[[pp]] %*% W1
    expected <- -colSums((object$score_freq[[pp]] * (gamma1 / gamma0)) %*% W1)
    as.numeric(observed) + expected
  })
  Reduce(`+`, pieces)
}

conditional_objective <- function(object, par) {
  esf_par0 <- as.numeric(object$W %*% as.numeric(par) + object$b_const)
  conditional_parts <- vapply(seq_len(object$NP), function(pp) {
    indices <- object$parm_index[[pp]]
    b <- esf_par0[indices]
    values <- split(b, object$pars_info$itemid[indices])
    esf0 <- psychotools::elementary_symmetric_functions(
      par = values, order = 0, diff = FALSE
    )[[1L]]
    -sum(object$suffstat[[pp]] * b) -
      sum(object$score_freq[[pp]] * log(esf0))
  }, numeric(1))
  -sum(conditional_parts)
}

finite_difference_gradient <- function(object, par) {
  vapply(seq_along(par), function(index) {
    # The objective is O(10^3), so 1e-6 produces subtraction near floating-
    # point resolution.  A centered 1e-4 step gives a more stable audit while
    # retaining O(h^2) truncation error for this smooth conditional likelihood.
    step <- 1e-4 * max(1, abs(par[[index]]))
    plus <- par
    minus <- par
    plus[[index]] <- plus[[index]] + step
    minus[[index]] <- minus[[index]] - step
    (conditional_objective(object, plus) - conditional_objective(object, minus)) / (2 * step)
  }, numeric(1))
}

finite_difference_information <- function(object, par) {
  hessian <- vapply(seq_along(par), function(index) {
    step <- 1e-4 * max(1, abs(par[[index]]))
    plus <- par
    minus <- par
    plus[[index]] <- plus[[index]] + step
    minus[[index]] <- minus[[index]] - step
    (conditional_gradient(object, plus) - conditional_gradient(object, minus)) / (2 * step)
  }, numeric(length(par)))
  0.5 * (hessian + t(hessian))
}

rank_evidence <- function(information) {
  singular <- svd(information, nu = 0, nv = 0)$d
  tolerance <- max(1e-9, max(singular) * ncol(information) * 1e-8)
  list(
    rank = sum(singular > tolerance),
    nullity = ncol(information) - sum(singular > tolerance),
    condition = if (min(singular) > tolerance) max(singular) / min(singular) else Inf,
    minimum_eigenvalue = min(eigen(information, symmetric = TRUE, only.values = TRUE)$values),
    tolerance = tolerance
  )
}

base_run <- function(manifest_row, python_row, rows, anchor_rows) {
  data.frame(
    SchemaVersion = as.character(manifest_row$SchemaVersion),
    RunId = as.character(manifest_row$RunId),
    ConditionId = as.character(manifest_row$ConditionId),
    Design = as.character(manifest_row$Design),
    TruthBias = as.numeric(manifest_row$TruthBias),
    TruthPositive = as.logical(manifest_row$TruthPositive),
    Replicate = as.integer(manifest_row$Replicate),
    Seed = as.numeric(manifest_row$Seed),
    Engine = "immer",
    Estimator = "CMLE",
    Mode = "IMMER_CML_MATCHED_UNANCHORED",
    Model = "RSM",
    Estimand = "unanchored additive structural parameters after conditioning out Person",
    Rows = as.integer(rows),
    RequestedAnchorRows = as.integer(anchor_rows),
    RequestedAnchorSupported = anchor_rows == 0L,
    BiasEstimandSupported = FALSE,
    SharedConditionalDesignEligible = as.logical(python_row$ConditionalDesignEligible),
    FitAttempted = FALSE,
    FitReturned = FALSE,
    Converged = FALSE,
    InferenceReady = FALSE,
    ParityEligible = FALSE,
    RequestedConditionEligible = FALSE,
    FailureStage = "prefit",
    FailureReason = "",
    ConditionalRankFromPythonAudit = as.numeric(python_row$ConditionalRank),
    ConditionalNullityFromPythonAudit = as.numeric(python_row$ConditionalNullity),
    RConditionalRankAtZero = NA_real_,
    RConditionalNullityAtZero = NA_real_,
    RConditionalDesignEligible = FALSE,
    PrefitRankAgreement = FALSE,
    KParams = as.numeric(python_row$KParams),
    ConditionalLogLik = NA_real_,
    ConditionalDeviance = NA_real_,
    ConvergenceCode = NA_integer_,
    FunctionEvaluations = NA_real_,
    GradientEvaluations = NA_real_,
    GradientSupNorm = NA_real_,
    GradientFiniteDifferenceMaxAbsError = NA_real_,
    GradientPrimaryThresholdMargin = NA_real_,
    GradientCheckSupportsPrimaryClassification = FALSE,
    StationarityTolerance = 1e-5,
    GradientToleranceRatio = NA_real_,
    GradientReadyAt1e4 = FALSE,
    GradientReadyAt1e5 = FALSE,
    GradientReadyAt1e6 = FALSE,
    InformationRank = NA_real_,
    InformationNullity = NA_real_,
    InformationConditionNumber = NA_real_,
    InformationMinEigenvalue = NA_real_,
    DeclaredSupportSentinelRows = 1L,
    DeclaredSupportSentinelWeight = 0,
    Warnings = "",
    ElapsedSeconds = NA_real_,
    BundleInventorySHA256 = bundle_inventory_sha256,
    ImmerVersion = as.character(utils::packageVersion("immer")),
    ImmerFunctionSHA256 = immer_function_hash,
    stringsAsFactors = FALSE
  )
}

run_rows <- list()
coefficient_rows <- list()
for (manifest_index in seq_len(nrow(manifest))) {
  manifest_row <- manifest[manifest_index, , drop = FALSE]
  run_id <- as.character(manifest_row$RunId)
  data <- ratings[ratings$RunId == run_id, c("Person", "Rater", "Task", "Criterion", "Score"), drop = FALSE]
  anchor_rows <- sum(anchors$RunId == run_id)
  python_row <- python_runs[python_runs$RunId == run_id, , drop = FALSE]
  run <- base_run(manifest_row, python_row, nrow(data), anchor_rows)
  prefit_started <- proc.time()[["elapsed"]]
  wide <- make_wide(data)
  design <- tryCatch(
    build_design(wide$response, wide$mapping, as.integer(manifest_row$Categories) - 1L),
    error = function(error) error
  )
  if (inherits(design, "error")) {
    run$FailureReason <- paste0("R matched-design construction failed: ", conditionMessage(design))
    run$ElapsedSeconds <- proc.time()[["elapsed"]] - prefit_started
    run_rows[[length(run_rows) + 1L]] <- run
    next
  }
  prefit_object <- design$prepared
  prefit_object$W <- design$W
  prefit_object$b_const <- design$b_const
  information_zero <- finite_difference_information(
    prefit_object, rep(0, ncol(design$W))
  )
  r_prefit <- rank_evidence(information_zero)
  run$RConditionalRankAtZero <- r_prefit$rank
  run$RConditionalNullityAtZero <- r_prefit$nullity
  run$RConditionalDesignEligible <- r_prefit$nullity == 0L &&
    r_prefit$minimum_eigenvalue > r_prefit$tolerance
  run$PrefitRankAgreement <-
    run$ConditionalRankFromPythonAudit == run$RConditionalRankAtZero &&
    run$ConditionalNullityFromPythonAudit == run$RConditionalNullityAtZero
  if (!isTRUE(run$PrefitRankAgreement)) {
    run$FailureReason <- paste0(
      "Python/R conditional prefit rank disagreement: Python ",
      run$ConditionalRankFromPythonAudit, "/", run$KParams,
      "; R ", run$RConditionalRankAtZero, "/", run$KParams
    )
    run$ElapsedSeconds <- proc.time()[["elapsed"]] - prefit_started
    run_rows[[length(run_rows) + 1L]] <- run
    next
  }
  if (!isTRUE(as.logical(python_row$ConditionalDesignEligible)) ||
      !isTRUE(run$RConditionalDesignEligible)) {
    run$FailureReason <- paste0(
      "independent Python/R conditional-information prefit rejection: rank ",
      run$RConditionalRankAtZero, "/", run$KParams,
      " (nullity ", run$RConditionalNullityAtZero, "); ",
      as.character(python_row$FailureReason)
    )
    run$ElapsedSeconds <- proc.time()[["elapsed"]] - prefit_started
    run_rows[[length(run_rows) + 1L]] <- run
    next
  }

  warning_messages <- character()
  started <- proc.time()[["elapsed"]]
  run$FitAttempted <- TRUE
  fit <- tryCatch(
    withCallingHandlers(
      immer::immer_cml(
        design$response,
        weights = design$weights,
        W = design$W,
        par_init = rep(0, ncol(design$W)),
        nullcats = "keep",
        use_rcpp = FALSE,
        control = list(maxit = 2000L, reltol = 1e-12)
      ),
      warning = function(warning) {
        warning_messages <<- c(warning_messages, conditionMessage(warning))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(error) error
  )
  run$ElapsedSeconds <- proc.time()[["elapsed"]] - started
  run$Warnings <- paste(unique(warning_messages), collapse = " | ")
  if (inherits(fit, "error")) {
    run$FailureStage <- "fit"
    run$FailureReason <- paste0(class(fit)[[1L]], ": ", conditionMessage(fit))
    run_rows[[length(run_rows) + 1L]] <- run
    next
  }

  run$FitReturned <- TRUE
  run$ConditionalLogLik <- as.numeric(fit$loglike)
  run$ConditionalDeviance <- as.numeric(fit$deviance)
  run$ConvergenceCode <- as.integer(fit$result_optim$convergence)
  run$FunctionEvaluations <- as.numeric(fit$result_optim$counts[["function"]])
  run$GradientEvaluations <- as.numeric(fit$result_optim$counts[["gradient"]])
  gradient <- conditional_gradient(fit, as.numeric(fit$coefficients))
  numeric_gradient <- finite_difference_gradient(fit, as.numeric(fit$coefficients))
  run$GradientSupNorm <- max(abs(gradient))
  run$GradientFiniteDifferenceMaxAbsError <- max(abs(gradient - numeric_gradient))
  run$GradientPrimaryThresholdMargin <- abs(run$GradientSupNorm - run$StationarityTolerance)
  run$GradientCheckSupportsPrimaryClassification <-
    run$GradientFiniteDifferenceMaxAbsError < run$GradientPrimaryThresholdMargin
  run$GradientToleranceRatio <- run$GradientSupNorm / run$StationarityTolerance
  run$GradientReadyAt1e4 <- is.finite(run$GradientSupNorm) && run$GradientSupNorm <= 1e-4
  run$GradientReadyAt1e5 <- is.finite(run$GradientSupNorm) && run$GradientSupNorm <= 1e-5
  run$GradientReadyAt1e6 <- is.finite(run$GradientSupNorm) && run$GradientSupNorm <= 1e-6
  information <- 0.5 * (fit$result_optim$hessian + t(fit$result_optim$hessian))
  fitted_rank <- rank_evidence(information)
  information_rank <- fitted_rank$rank
  run$InformationRank <- fitted_rank$rank
  run$InformationNullity <- fitted_rank$nullity
  run$InformationConditionNumber <- fitted_rank$condition
  run$InformationMinEigenvalue <- fitted_rank$minimum_eigenvalue
  run$Converged <- run$ConvergenceCode == 0L && is.finite(run$ConditionalLogLik)
  run$InferenceReady <- run$Converged && is.finite(run$GradientSupNorm) &&
    run$GradientSupNorm <= run$StationarityTolerance &&
    information_rank == ncol(information) &&
    fitted_rank$minimum_eigenvalue > fitted_rank$tolerance &&
    all(is.finite(fit$vcov))
  run$ParityEligible <- run$InferenceReady
  run$RequestedConditionEligible <- run$ParityEligible && anchor_rows == 0L
  if (!run$InferenceReady) {
    run$FailureStage <- "readiness"
    run$FailureReason <- paste0(
      "immer CMLE readiness withheld: code=", run$ConvergenceCode,
      "; gradient_sup=", format(run$GradientSupNorm, digits = 8),
      "; information_rank=", information_rank, "/", ncol(information)
    )
  } else if (anchor_rows > 0L) {
    run$FailureStage <- "requested_scope"
    run$FailureReason <- "anchors requested by condition but unsupported in matched CMLE estimand; unanchored parity fit retained"
  } else {
    run$FailureStage <- ""
    run$FailureReason <- ""
  }

  coefficients <- data.frame(
    RunId = run_id,
    ConditionId = as.character(manifest_row$ConditionId),
    Engine = "immer",
    Estimator = "CMLE",
    Mode = "IMMER_CML_MATCHED_UNANCHORED",
    Parameter = colnames(fit$W),
    Estimate = as.numeric(fit$coefficients),
    Gradient = as.numeric(gradient),
    SE = sqrt(pmax(diag(fit$vcov), 0)),
    IncludedInParity = isTRUE(run$ParityEligible),
    stringsAsFactors = FALSE
  )
  coefficient_rows[[length(coefficient_rows) + 1L]] <- coefficients
  run_rows[[length(run_rows) + 1L]] <- run
}

runs <- do.call(rbind, run_rows)
coefficients <- if (length(coefficient_rows)) do.call(rbind, coefficient_rows) else data.frame()
utils::write.csv(runs, file.path(output_dir, "immer_cmle_runs.csv"), row.names = FALSE)
utils::write.csv(coefficients, file.path(output_dir, "immer_cmle_coefficients.csv"), row.names = FALSE)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1L]]) else NA_character_
identity <- data.frame(
  BundleInventorySHA256 = bundle_inventory_sha256,
  ImmerVersion = as.character(utils::packageVersion("immer")),
  ImmerFunctionSHA256 = immer_function_hash,
  AdapterScriptSHA256 = if (!is.na(script_path)) {
    digest::digest(file = normalizePath(script_path), algo = "sha256", serialize = FALSE)
  } else NA_character_,
  Mode = "IMMER_CML_MATCHED_UNANCHORED",
  FacetOrder = "Rater|Task|Criterion",
  RatingSupport = "zero-weight sentinel declares 0..Categories-1",
  NullCategoryPolicy = "keep",
  AnchorScope = "unsupported; unanchored parity fit only",
  BiasScope = "unsupported by additive CMLE estimand",
  stringsAsFactors = FALSE
)
utils::write.csv(identity, file.path(output_dir, "immer_cmle_adapter_identity.csv"), row.names = FALSE)
