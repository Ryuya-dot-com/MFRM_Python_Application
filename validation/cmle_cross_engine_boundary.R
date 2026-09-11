#!/usr/bin/env Rscript

# Registered boundary-semantics comparison for native Python CMLE, matched
# immer CMLE, and estimator-mismatched mfrmr/TAM/sirt sensitivity fits.

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

args <- parse_args(commandArgs(trailingOnly = TRUE))
for (name in c("input", "output", "plan-sha256", "mfrmr-git-head",
               "mfrmr-status-sha256", "mfrmr-content-sha256")) {
  if (!nzchar(args[[name]] %||% "")) stop("Missing required argument --", name, call. = FALSE)
}
for (package in c("digest", "immer", "psychotools", "mfrmr", "TAM", "sirt")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Required package unavailable: ", package, call. = FALSE)
  }
}

input_dir <- normalizePath(args$input, mustWork = TRUE)
output_dir <- normalizePath(args$output, mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
manifest <- utils::read.csv(
  file.path(input_dir, "manifest.csv"), stringsAsFactors = FALSE, check.names = FALSE
)
ratings <- utils::read.csv(
  file.path(input_dir, "ratings.csv"), stringsAsFactors = FALSE, check.names = FALSE
)
python_runs <- utils::read.csv(
  file.path(input_dir, "python_runs.csv"), stringsAsFactors = FALSE, check.names = FALSE
)
inventory <- utils::read.csv(
  file.path(input_dir, "bundle_inventory.csv"), stringsAsFactors = FALSE, check.names = FALSE
)

hash_file <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

hash_ok <- vapply(seq_len(nrow(inventory)), function(index) {
  path <- file.path(input_dir, as.character(inventory$File[[index]]))
  file.exists(path) && identical(hash_file(path), as.character(inventory$SHA256[[index]]))
}, logical(1))
if (!all(hash_ok)) stop("Input bundle byte identity failed.", call. = FALSE)
if (!identical(sort(as.character(manifest$CaseId)), sort(as.character(python_runs$CaseId)))) {
  stop("Manifest/Python case identities differ.", call. = FALSE)
}

function_hash <- function(package, functions) {
  namespace <- asNamespace(package)
  hashes <- vapply(functions, function(name) {
    fn <- get(name, envir = namespace, inherits = FALSE)
    digest::digest(
      list(formals = formals(fn), body = body(fn)),
      algo = "sha256", serialize = TRUE
    )
  }, character(1))
  digest::digest(hashes, algo = "sha256", serialize = TRUE)
}

sum_zero <- function(n) {
  if (n <= 1L) return(matrix(0, nrow = n, ncol = 0L))
  output <- matrix(0, nrow = n, ncol = n - 1L)
  output[seq_len(n - 1L), ] <- diag(n - 1L)
  output[n, ] <- -1
  output
}

make_wide <- function(data) {
  data$VirtualUnit <- paste(data$Rater, data$Unit, sep = "__")
  persons <- sort(unique(as.character(data$Person)))
  units <- sort(unique(as.character(data$VirtualUnit)))
  response <- matrix(
    NA_real_, nrow = length(persons), ncol = length(units),
    dimnames = list(persons, units)
  )
  indices <- cbind(match(data$Person, persons), match(data$VirtualUnit, units))
  if (anyDuplicated(data.frame(indices))) {
    stop("Person/virtual-unit response is duplicated.", call. = FALSE)
  }
  response[indices] <- as.numeric(data$Score)
  mapping <- unique(data[c("VirtualUnit", "Rater", "Unit")])
  mapping <- mapping[match(units, mapping$VirtualUnit), , drop = FALSE]
  list(response = response, mapping = mapping, persons = persons)
}

make_person_rater_wide <- function(data) {
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
  list(response = as.data.frame(response), grid = grid, items = items)
}

build_immer_design <- function(response, mapping, rating_max) {
  augmented <- rbind(response, rep(as.numeric(rating_max), ncol(response)))
  rownames(augmented)[nrow(augmented)] <- "__declared_support_zero_weight__"
  weights <- c(rep(1, nrow(response)), 0)
  prepared <- immer:::lpcm_data_prep(augmented, weights = weights, a = NULL)
  if (!all(prepared$maxK == rating_max)) {
    stop("immer support sentinel failed to declare 0..rating_max.")
  }
  rater_levels <- sort(unique(as.character(mapping$Rater)))
  rater_contrast <- sum_zero(length(rater_levels))
  step_contrast <- sum_zero(as.integer(rating_max))
  parameter_names <- c(
    if (ncol(rater_contrast)) {
      paste0("facet:Rater:free:", rater_levels[seq_len(ncol(rater_contrast))])
    } else character(),
    if (ncol(step_contrast)) {
      paste0("step:__shared__:free:", seq_len(ncol(step_contrast)))
    } else character()
  )
  W <- matrix(
    0, nrow = nrow(prepared$pars_info), ncol = length(parameter_names),
    dimnames = list(rownames(prepared$pars_info), parameter_names)
  )
  for (row_index in seq_len(nrow(prepared$pars_info))) {
    unit <- as.character(prepared$pars_info$item[[row_index]])
    category <- as.integer(prepared$pars_info$cat[[row_index]])
    meta <- mapping[mapping$VirtualUnit == unit, , drop = FALSE]
    if (nrow(meta) != 1L) stop("Missing/duplicate virtual-unit metadata: ", unit)
    if (ncol(rater_contrast)) {
      level_index <- match(as.character(meta$Rater[[1L]]), rater_levels)
      W[row_index, seq_len(ncol(rater_contrast))] <-
        category * rater_contrast[level_index, ]
    }
    if (ncol(step_contrast) && category > 0L) {
      columns <- ncol(rater_contrast) + seq_len(ncol(step_contrast))
      W[row_index, columns] <- colSums(
        step_contrast[seq_len(category), , drop = FALSE]
      )
    }
  }
  list(
    response = augmented, weights = weights, W = W,
    prepared = prepared, b_const = rep(0, nrow(W))
  )
}

conditional_gradient <- function(object, par) {
  linear <- as.numeric(object$W %*% as.numeric(par) + object$b_const)
  pieces <- lapply(seq_len(object$NP), function(person_pattern) {
    indices <- object$parm_index[[person_pattern]]
    values <- split(linear[indices], object$pars_info$itemid[indices])
    esf <- psychotools::elementary_symmetric_functions(
      par = values, order = 1, diff = FALSE
    )
    W1 <- object$W[indices, , drop = FALSE]
    observed <- object$suffstat[[person_pattern]] %*% W1
    expected <- -colSums(
      (object$score_freq[[person_pattern]] * (esf[[2L]] / esf[[1L]])) %*% W1
    )
    as.numeric(observed) + expected
  })
  Reduce(`+`, pieces)
}

conditional_objective <- function(object, par) {
  linear <- as.numeric(object$W %*% as.numeric(par) + object$b_const)
  parts <- vapply(seq_len(object$NP), function(person_pattern) {
    indices <- object$parm_index[[person_pattern]]
    values <- split(linear[indices], object$pars_info$itemid[indices])
    esf0 <- psychotools::elementary_symmetric_functions(
      par = values, order = 0, diff = FALSE
    )[[1L]]
    -sum(object$suffstat[[person_pattern]] * linear[indices]) -
      sum(object$score_freq[[person_pattern]] * log(esf0))
  }, numeric(1))
  -sum(parts)
}

finite_difference_gradient <- function(object, par) {
  vapply(seq_along(par), function(index) {
    step <- 1e-4 * max(1, abs(par[[index]]))
    plus <- minus <- par
    plus[[index]] <- plus[[index]] + step
    minus[[index]] <- minus[[index]] - step
    (conditional_objective(object, plus) - conditional_objective(object, minus)) /
      (2 * step)
  }, numeric(1))
}

finite_difference_information <- function(object, par) {
  hessian <- vapply(seq_along(par), function(index) {
    step <- 1e-4 * max(1, abs(par[[index]]))
    plus <- minus <- par
    plus[[index]] <- plus[[index]] + step
    minus[[index]] <- minus[[index]] - step
    (conditional_gradient(object, plus) - conditional_gradient(object, minus)) /
      (2 * step)
  }, numeric(length(par)))
  0.5 * (hessian + t(hessian))
}

rank_evidence <- function(information) {
  singular <- svd(information, nu = 0, nv = 0)$d
  tolerance <- max(1e-9, max(singular) * ncol(information) * 1e-8)
  retained_rank <- sum(singular > tolerance)
  eigenvalues <- eigen(information, symmetric = TRUE, only.values = TRUE)$values
  list(
    rank = retained_rank,
    nullity = ncol(information) - retained_rank,
    condition = if (min(singular) > tolerance) max(singular) / min(singular) else Inf,
    minimum_eigenvalue = min(eigenvalues),
    tolerance = tolerance
  )
}

safe_max_abs <- function(values) {
  values <- as.numeric(values)
  if (!length(values) || !any(is.finite(values))) return(NA_real_)
  max(abs(values[is.finite(values)]))
}

base_identity <- function(case_id, engine, estimator, mode) {
  manifest_row <- manifest[manifest$CaseId == case_id, , drop = FALSE]
  python_row <- python_runs[python_runs$CaseId == case_id, , drop = FALSE]
  data.frame(
    CaseId = case_id,
    Family = as.character(manifest_row$Family),
    Mechanism = as.character(manifest_row$Mechanism),
    ExpectedStatus = as.character(manifest_row$ExpectedStatus),
    PythonWorkflowStatus = as.character(python_row$WorkflowStatus),
    Engine = engine,
    Estimator = estimator,
    Mode = mode,
    stringsAsFactors = FALSE
  )
}

immer_rows <- list()
immer_coefficients <- list()
for (case_id in as.character(manifest$CaseId)) {
  data <- ratings[ratings$CaseId == case_id, , drop = FALSE]
  manifest_row <- manifest[manifest$CaseId == case_id, , drop = FALSE]
  python_row <- python_runs[python_runs$CaseId == case_id, , drop = FALSE]
  wide <- make_wide(data)
  design <- build_immer_design(
    wide$response, wide$mapping, as.integer(manifest_row$RatingMax)
  )
  prefit <- design$prepared
  prefit$W <- design$W
  prefit$b_const <- design$b_const
  information_zero <- finite_difference_information(
    prefit, rep(0, ncol(design$W))
  )
  prefit_rank <- rank_evidence(information_zero)
  rank_agreement <-
    as.integer(python_row$PrefitRank) == prefit_rank$rank &&
    as.integer(python_row$PrefitNullity) == prefit_rank$nullity

  for (maxit in c(200L, 2000L, 10000L)) {
    row <- base_identity(case_id, "immer", "CMLE", paste0("maxit_", maxit))
    row$MaxIt <- maxit
    row$KParams <- ncol(design$W)
    row$PythonPrefitRank <- as.integer(python_row$PrefitRank)
    row$PythonPrefitNullity <- as.integer(python_row$PrefitNullity)
    row$RPrefitRank <- prefit_rank$rank
    row$RPrefitNullity <- prefit_rank$nullity
    row$PrefitRankAgreement <- rank_agreement
    row$FitAttempted <- FALSE
    row$FitReturned <- FALSE
    row$ConvergenceCode <- NA_integer_
    row$ConditionalLogLik <- NA_real_
    row$FunctionEvaluations <- NA_real_
    row$GradientEvaluations <- NA_real_
    row$GradientSupNorm <- NA_real_
    row$GradientFiniteDifferenceMaxAbsError <- NA_real_
    row$GradientReadyAt1e4 <- FALSE
    row$GradientReadyAt1e5 <- FALSE
    row$GradientReadyAt1e6 <- FALSE
    row$GradientReadyAt1e8 <- FALSE
    row$InformationRank <- NA_real_
    row$InformationNullity <- NA_real_
    row$InformationConditionNumber <- NA_real_
    row$InformationMinEigenvalue <- NA_real_
    row$AllCoefficientsFinite <- FALSE
    row$AllSEFinite <- FALSE
    row$MaxAbsEstimate <- NA_real_
    row$MaxAbsSE <- NA_real_
    row$Warnings <- ""
    row$Error <- ""
    row$ElapsedSeconds <- 0
    if (!rank_agreement || prefit_rank$nullity != 0L) {
      row$Error <- if (!rank_agreement) {
        "Python/R prefit-rank disagreement"
      } else "structural nonidentification; optimizer skipped"
      immer_rows[[length(immer_rows) + 1L]] <- row
      next
    }

    warnings <- character()
    started <- proc.time()[["elapsed"]]
    row$FitAttempted <- TRUE
    fit <- tryCatch(
      withCallingHandlers(
        immer::immer_cml(
          design$response,
          weights = design$weights,
          W = design$W,
          par_init = rep(0, ncol(design$W)),
          nullcats = "keep",
          use_rcpp = FALSE,
          control = list(maxit = maxit, reltol = 1e-12)
        ),
        warning = function(warning) {
          warnings <<- c(warnings, conditionMessage(warning))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(error) error
    )
    row$ElapsedSeconds <- proc.time()[["elapsed"]] - started
    row$Warnings <- paste(unique(warnings), collapse = " | ")
    if (inherits(fit, "error")) {
      row$Error <- paste0(class(fit)[[1L]], ": ", conditionMessage(fit))
      immer_rows[[length(immer_rows) + 1L]] <- row
      next
    }

    row$FitReturned <- TRUE
    coefficients <- as.numeric(fit$coefficients)
    names(coefficients) <- colnames(fit$W)
    standard_errors <- if (!is.null(fit$vcov) &&
      all(dim(fit$vcov) == c(length(coefficients), length(coefficients)))) {
      sqrt(pmax(diag(fit$vcov), 0))
    } else rep(NA_real_, length(coefficients))
    gradient <- conditional_gradient(fit, coefficients)
    numeric_gradient <- finite_difference_gradient(fit, coefficients)
    information <- 0.5 * (fit$result_optim$hessian + t(fit$result_optim$hessian))
    fitted_rank <- rank_evidence(information)
    row$ConvergenceCode <- as.integer(fit$result_optim$convergence)
    row$ConditionalLogLik <- as.numeric(fit$loglike)
    row$FunctionEvaluations <- as.numeric(fit$result_optim$counts[["function"]])
    row$GradientEvaluations <- as.numeric(fit$result_optim$counts[["gradient"]])
    row$GradientSupNorm <- max(abs(gradient))
    row$GradientFiniteDifferenceMaxAbsError <- max(abs(gradient - numeric_gradient))
    row$GradientReadyAt1e4 <- row$GradientSupNorm <= 1e-4
    row$GradientReadyAt1e5 <- row$GradientSupNorm <= 1e-5
    row$GradientReadyAt1e6 <- row$GradientSupNorm <= 1e-6
    row$GradientReadyAt1e8 <- row$GradientSupNorm <= 1e-8
    row$InformationRank <- fitted_rank$rank
    row$InformationNullity <- fitted_rank$nullity
    row$InformationConditionNumber <- fitted_rank$condition
    row$InformationMinEigenvalue <- fitted_rank$minimum_eigenvalue
    row$AllCoefficientsFinite <- all(is.finite(coefficients))
    row$AllSEFinite <- all(is.finite(standard_errors))
    row$MaxAbsEstimate <- safe_max_abs(coefficients)
    row$MaxAbsSE <- safe_max_abs(standard_errors)
    coefficient_table <- data.frame(
      CaseId = case_id,
      MaxIt = maxit,
      Parameter = names(coefficients),
      Estimate = coefficients,
      SE = standard_errors,
      Gradient = as.numeric(gradient),
      stringsAsFactors = FALSE
    )
    immer_coefficients[[length(immer_coefficients) + 1L]] <- coefficient_table
    immer_rows[[length(immer_rows) + 1L]] <- row
  }
}

run_mfrmr <- function(case_id, data, manifest_row) {
  row <- base_identity(case_id, "mfrmr", "JML", "fit_mfrm_JML")
  row$FitAttempted <- TRUE
  row$FitReturned <- FALSE
  row$Converged <- FALSE
  row$InferenceReady <- FALSE
  row$Iterations <- NA_real_
  row$ConvergenceCode <- NA_real_
  row$LogLik <- NA_real_
  row$MaxAbsEstimate <- NA_real_
  row$MaxAbsSE <- NA_real_
  row$NonfiniteEstimateCount <- NA_real_
  row$NonfiniteSECount <- NA_real_
  row$Warnings <- ""
  row$Error <- ""
  warnings <- character()
  started <- proc.time()[["elapsed"]]
  fit <- tryCatch(
    withCallingHandlers(
      mfrmr::fit_mfrm(
        data = data[c("Person", "Rater", "Score")],
        person = "Person", facets = "Rater", score = "Score",
        rating_min = 0, rating_max = as.integer(manifest_row$RatingMax),
        keep_original = TRUE, model = "RSM", method = "JML",
        noncenter_facet = "Person", min_obs_per_element = 1L,
        min_obs_per_category = 1L, maxit = 800L, reltol = 1e-10,
        optimizer = "BFGS"
      ),
      warning = function(warning) {
        warnings <<- c(warnings, conditionMessage(warning))
        invokeRestart("muffleWarning")
      }
    ), error = function(error) error
  )
  row$ElapsedSeconds <- proc.time()[["elapsed"]] - started
  row$Warnings <- paste(unique(warnings), collapse = " | ")
  if (inherits(fit, "error")) {
    row$Error <- paste0(class(fit)[[1L]], ": ", conditionMessage(fit))
    return(row)
  }
  row$FitReturned <- TRUE
  summary <- as.data.frame(fit$summary, stringsAsFactors = FALSE)[1L, , drop = FALSE]
  value <- function(name, default = NA) if (name %in% names(summary)) summary[[name]][[1L]] else default
  row$Converged <- isTRUE(as.logical(value("Converged", FALSE)))
  row$InferenceReady <- isTRUE(as.logical(value("InferenceReady", FALSE)))
  row$Iterations <- as.numeric(value("Iterations", NA_real_))
  row$ConvergenceCode <- as.numeric(value("ConvergenceCode", NA_real_))
  row$LogLik <- as.numeric(value("LogLik", NA_real_))
  tables <- Filter(Negate(is.null), list(fit$facets$others, fit$steps))
  estimates <- unlist(lapply(tables, function(table) {
    table <- as.data.frame(table)
    if ("Estimate" %in% names(table)) table$Estimate else numeric()
  }), use.names = FALSE)
  standard_errors <- unlist(lapply(tables, function(table) {
    table <- as.data.frame(table)
    if ("SE" %in% names(table)) table$SE else numeric()
  }), use.names = FALSE)
  row$MaxAbsEstimate <- safe_max_abs(estimates)
  row$MaxAbsSE <- safe_max_abs(standard_errors)
  row$NonfiniteEstimateCount <- sum(!is.finite(as.numeric(estimates)))
  row$NonfiniteSECount <- sum(!is.finite(as.numeric(standard_errors)))
  row
}

run_tam <- function(case_id, data, manifest_row) {
  row <- base_identity(case_id, "TAM", "MML", "tam.mml.mfr")
  row$FitAttempted <- TRUE
  row$FitReturned <- FALSE
  row$Converged <- FALSE
  row$Iterations <- NA_real_
  row$LogLik <- NA_real_
  row$MaxAbsEstimate <- NA_real_
  row$MaxAbsSE <- NA_real_
  row$NonfiniteEstimateCount <- NA_real_
  row$NonfiniteSECount <- NA_real_
  row$ObservedSupportMax <- max(as.numeric(data$Score))
  row$DeclaredRatingMax <- as.integer(manifest_row$RatingMax)
  row$SupportMismatch <- row$ObservedSupportMax < row$DeclaredRatingMax
  row$Warnings <- ""
  row$Error <- ""
  if (row$SupportMismatch) {
    row$FitAttempted <- FALSE
    row$Error <- paste(
      "category_support: TAM infers response-unit support from observed maxima;",
      "declared but unused top category not injected"
    )
    return(row)
  }
  wide <- make_person_rater_wide(data)
  facet_formula <- if (as.integer(manifest_row$RatingMax) > 1L) {
    ~ rater + step
  } else {
    ~ rater
  }
  warnings <- character()
  started <- proc.time()[["elapsed"]]
  fit <- tryCatch({
    invisible(utils::capture.output(
      value <- withCallingHandlers(
        TAM::tam.mml.mfr(
          resp = wide$response,
          facets = data.frame(rater = wide$grid$Rater, stringsAsFactors = FALSE),
          pid = wide$grid$Person,
          formulaA = facet_formula,
          constraint = "cases",
          control = list(
            nodes = seq(-6, 6, length.out = 41L), maxiter = 400L,
            convD = 1e-4, conv = 1e-4, convM = 1e-4,
            Msteps = 4L, dev_crit = "absolute", progress = FALSE
          ), verbose = FALSE
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
    return(row)
  }
  row$FitReturned <- TRUE
  row$Iterations <- as.numeric(fit$iter %||% NA_real_)
  row$Converged <- is.finite(row$Iterations) && row$Iterations < 400L
  row$LogLik <- -0.5 * as.numeric(fit$deviance %||% NA_real_)
  parameters <- as.data.frame(fit$xsi.facets %||% data.frame(), stringsAsFactors = FALSE)
  estimates <- if ("xsi" %in% names(parameters)) parameters$xsi else numeric()
  standard_errors <- if ("se.xsi" %in% names(parameters)) parameters$se.xsi else numeric()
  row$MaxAbsEstimate <- safe_max_abs(estimates)
  row$MaxAbsSE <- safe_max_abs(standard_errors)
  row$NonfiniteEstimateCount <- sum(!is.finite(as.numeric(estimates)))
  row$NonfiniteSECount <- sum(!is.finite(as.numeric(standard_errors)))
  row
}

run_sirt <- function(case_id, data, manifest_row) {
  row <- base_identity(case_id, "sirt", "MML", "rm.facets")
  row$FitAttempted <- TRUE
  row$FitReturned <- FALSE
  row$Converged <- FALSE
  row$Iterations <- NA_real_
  row$LogLik <- NA_real_
  row$MaxAbsEstimate <- NA_real_
  row$MaxAbsSE <- NA_real_
  row$NonfiniteEstimateCount <- NA_real_
  row$NonfiniteSECount <- NA_real_
  row$ObservedSupportMax <- max(as.numeric(data$Score))
  row$DeclaredRatingMax <- as.integer(manifest_row$RatingMax)
  row$SupportMismatch <- row$ObservedSupportMax < row$DeclaredRatingMax
  row$Warnings <- ""
  row$Error <- ""
  wide <- make_wide(data)
  warnings <- character()
  started <- proc.time()[["elapsed"]]
  fit <- tryCatch({
    invisible(utils::capture.output(
      value <- withCallingHandlers(
        sirt::rm.facets(
          dat = wide$response, pid = wide$persons,
          rater = wide$mapping$Rater,
          theta.k = seq(-6, 6, length.out = 41L),
          est.b.rater = TRUE, est.a.item = FALSE, est.a.rater = FALSE,
          rater_item_int = FALSE, est.mean = FALSE,
          b.rater.center = 2, maxdevchange = 0.01,
          globconv = 0.001, maxiter = 400L
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
    return(row)
  }
  row$FitReturned <- TRUE
  row$Iterations <- as.numeric(fit$iter %||% NA_real_)
  row$Converged <- is.finite(row$Iterations) && row$Iterations < 400L
  row$LogLik <- as.numeric(tryCatch(stats::logLik(fit), error = function(error) NA_real_))
  estimates <- c(as.numeric(fit$b.rater %||% numeric()), as.numeric(fit$tau.item %||% numeric()))
  standard_errors <- c(
    as.numeric(fit$se.b.rater %||% numeric()),
    as.numeric(fit$se.tau.item %||% numeric())
  )
  row$MaxAbsEstimate <- safe_max_abs(estimates)
  row$MaxAbsSE <- safe_max_abs(standard_errors)
  row$NonfiniteEstimateCount <- sum(!is.finite(estimates))
  row$NonfiniteSECount <- sum(!is.finite(standard_errors))
  row
}

secondary_rows <- list()
for (case_id in as.character(manifest$CaseId)) {
  data <- ratings[ratings$CaseId == case_id, , drop = FALSE]
  manifest_row <- manifest[manifest$CaseId == case_id, , drop = FALSE]
  secondary_rows[[length(secondary_rows) + 1L]] <- run_mfrmr(case_id, data, manifest_row)
  secondary_rows[[length(secondary_rows) + 1L]] <- run_tam(case_id, data, manifest_row)
}

rbind_fill <- function(frames) {
  columns <- unique(unlist(lapply(frames, names), use.names = FALSE))
  normalized <- lapply(frames, function(frame) {
    missing <- setdiff(columns, names(frame))
    for (name in missing) frame[[name]] <- NA
    frame[columns]
  })
  do.call(rbind, normalized)
}

immer_runs <- do.call(rbind, immer_rows)
immer_coefficient_table <- if (length(immer_coefficients)) {
  do.call(rbind, immer_coefficients)
} else data.frame()
secondary_runs <- rbind_fill(secondary_rows)
utils::write.csv(immer_runs, file.path(output_dir, "immer_runs.csv"), row.names = FALSE)
utils::write.csv(
  immer_coefficient_table,
  file.path(output_dir, "immer_coefficients.csv"), row.names = FALSE
)
utils::write.csv(
  secondary_runs, file.path(output_dir, "secondary_engine_runs.csv"), row.names = FALSE
)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1L]]) else NA_character_
identity <- data.frame(
  PlanSHA256 = as.character(args[["plan-sha256"]]),
  AdapterSHA256 = if (!is.na(script_path)) hash_file(normalizePath(script_path)) else NA_character_,
  RVersion = as.character(getRversion()),
  ImmerVersion = as.character(utils::packageVersion("immer")),
  ImmerFunctionSHA256 = function_hash("immer", c("immer_jml", "immer_cml")),
  MfrmrVersion = as.character(utils::packageVersion("mfrmr")),
  MfrmrFunctionSHA256 = function_hash("mfrmr", "fit_mfrm"),
  TAMVersion = as.character(utils::packageVersion("TAM")),
  TAMFunctionSHA256 = function_hash("TAM", "tam.mml.mfr"),
  SirtVersion = as.character(utils::packageVersion("sirt")),
  SirtFunctionSHA256 = function_hash("sirt", "rm.facets"),
  MfrmrSourceGitHead = as.character(args[["mfrmr-git-head"]]),
  MfrmrSourceStatusSHA256 = as.character(args[["mfrmr-status-sha256"]]),
  MfrmrSourceContentSHA256 = as.character(args[["mfrmr-content-sha256"]]),
  stringsAsFactors = FALSE
)
utils::write.csv(identity, file.path(output_dir, "r_engine_identity.csv"), row.names = FALSE)
