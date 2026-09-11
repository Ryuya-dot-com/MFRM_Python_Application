#!/usr/bin/env Rscript

# Matched PCM hard-anchor CMLE adapter plus different-estimator sensitivities.

parse_args <- function(values) {
  output <- list(); index <- 1L
  while (index <= length(values)) {
    key <- sub("^--", "", values[[index]])
    if (index == length(values)) stop("Missing --", key, call. = FALSE)
    output[[key]] <- values[[index + 1L]]; index <- index + 2L
  }
  output
}
`%||%` <- function(value, replacement) {
  if (is.null(value) || length(value) == 0L) replacement else value
}
args <- parse_args(commandArgs(trailingOnly = TRUE))
for (name in c("input", "output", "plan-sha256", "mfrmr-source-sha256")) {
  if (!nzchar(args[[name]] %||% "")) stop("Missing --", name, call. = FALSE)
}
for (package in c("digest", "immer", "psychotools", "mfrmr", "TAM", "sirt")) {
  if (!requireNamespace(package, quietly = TRUE)) stop("Missing package: ", package)
}

input_dir <- normalizePath(args$input, mustWork = TRUE)
output_dir <- normalizePath(args$output, mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
read_input <- function(name) utils::read.csv(
  file.path(input_dir, name), stringsAsFactors = FALSE, check.names = FALSE
)
manifest <- read_input("manifest.csv")
ratings <- read_input("ratings.csv")
anchors <- read_input("anchors.csv")
python_runs <- read_input("python_runs.csv")
inventory <- read_input("bundle_inventory.csv")
hash_file <- function(path) digest::digest(
  file = path, algo = "sha256", serialize = FALSE
)
if (!all(vapply(seq_len(nrow(inventory)), function(index) {
  path <- file.path(input_dir, inventory$File[[index]])
  file.exists(path) && identical(hash_file(path), inventory$SHA256[[index]])
}, logical(1)))) stop("Bundle identity failed.")

function_hash <- function(package, functions, remove_source = TRUE) {
  namespace <- asNamespace(package)
  hashes <- vapply(functions, function(name) {
    fn <- get(name, envir = namespace, inherits = FALSE)
    if (isTRUE(remove_source)) fn <- utils::removeSource(fn)
    digest::digest(list(formals = formals(fn), body = body(fn)),
                   algo = "sha256", serialize = TRUE)
  }, character(1))
  digest::digest(hashes, algo = "sha256", serialize = TRUE)
}
sum_zero <- function(n) {
  if (n <= 1L) return(matrix(0, nrow = n, ncol = 0L))
  output <- matrix(0, nrow = n, ncol = n - 1L)
  output[seq_len(n - 1L), ] <- diag(n - 1L); output[n, ] <- -1
  output
}
make_wide <- function(data) {
  data$VirtualUnit <- paste(data$Rater, data$Criterion, sep = "__")
  persons <- sort(unique(data$Person)); units <- sort(unique(data$VirtualUnit))
  response <- matrix(NA_real_, length(persons), length(units),
                     dimnames = list(persons, units))
  indices <- cbind(match(data$Person, persons), match(data$VirtualUnit, units))
  if (anyDuplicated(data.frame(indices))) stop("Duplicate response unit.")
  response[indices] <- data$Score
  mapping <- unique(data[c("VirtualUnit", "Rater", "Criterion")])
  mapping <- mapping[match(units, mapping$VirtualUnit), , drop = FALSE]
  list(response = response, mapping = mapping, persons = persons)
}
make_person_rater_wide <- function(data) {
  key <- paste(data$Person, data$Rater, data$Criterion, sep = "\r")
  if (anyDuplicated(key)) stop("Duplicate Person x Rater x Criterion cell.")
  grid <- unique(data[c("Person", "Rater")])
  grid <- grid[order(grid$Person, grid$Rater), , drop = FALSE]
  items <- sort(unique(data$Criterion))
  response <- matrix(NA_real_, nrow(grid), length(items), dimnames = list(NULL, items))
  response[cbind(
    match(paste(data$Person, data$Rater), paste(grid$Person, grid$Rater)),
    match(data$Criterion, items)
  )] <- data$Score
  list(response = as.data.frame(response), grid = grid, items = items)
}
case_anchors <- function(case_id) {
  anchors[anchors$CaseId == case_id, c("Facet", "Level", "Value"), drop = FALSE]
}

build_immer_design <- function(response, mapping, rating_max, anchor_table) {
  augmented <- rbind(response, rep(rating_max, ncol(response)))
  rownames(augmented)[nrow(augmented)] <- "__support_zero_weight__"
  weights <- c(rep(1, nrow(response)), 0)
  prepared <- immer:::lpcm_data_prep(augmented, weights = weights, a = NULL)
  if (!all(prepared$maxK == rating_max)) stop("Support sentinel failed.")
  facet_names <- c("Rater", "Criterion")
  levels <- lapply(facet_names, function(name) sort(unique(mapping[[name]])))
  names(levels) <- facet_names
  contrasts <- list(); offsets <- list(); labels <- list()
  for (facet in facet_names) {
    facet_anchors <- anchor_table[anchor_table$Facet == facet, , drop = FALSE]
    if (nrow(facet_anchors)) {
      free <- setdiff(levels[[facet]], facet_anchors$Level)
      contrast <- matrix(0, length(levels[[facet]]), length(free))
      if (length(free)) {
        contrast[cbind(match(free, levels[[facet]]), seq_along(free))] <- 1
      }
      offset <- rep(0, length(levels[[facet]]))
      offset[match(facet_anchors$Level, levels[[facet]])] <- facet_anchors$Value
      contrasts[[facet]] <- contrast; offsets[[facet]] <- offset; labels[[facet]] <- free
    } else {
      contrasts[[facet]] <- sum_zero(length(levels[[facet]]))
      offsets[[facet]] <- rep(0, length(levels[[facet]]))
      labels[[facet]] <- levels[[facet]][seq_len(ncol(contrasts[[facet]]))]
    }
  }
  step_contrast <- sum_zero(rating_max)
  parameter_names <- c(
    unlist(lapply(facet_names, function(facet) {
      paste0("facet:", facet, ":free:", labels[[facet]])
    }), use.names = FALSE),
    unlist(lapply(levels$Criterion, function(level) {
      paste0("step:", level, ":free:", seq_len(ncol(step_contrast)))
    }), use.names = FALSE)
  )
  W <- matrix(0, nrow(prepared$pars_info), length(parameter_names),
              dimnames = list(rownames(prepared$pars_info), parameter_names))
  b_const <- rep(0, nrow(W))
  rater_offset <- 0L
  criterion_offset <- ncol(contrasts$Rater)
  step_offset <- criterion_offset + ncol(contrasts$Criterion)
  for (row_index in seq_len(nrow(prepared$pars_info))) {
    unit <- as.character(prepared$pars_info$item[[row_index]])
    category <- as.integer(prepared$pars_info$cat[[row_index]])
    meta <- mapping[mapping$VirtualUnit == unit, , drop = FALSE]
    rater_index <- match(meta$Rater[[1]], levels$Rater)
    criterion_index <- match(meta$Criterion[[1]], levels$Criterion)
    if (ncol(contrasts$Rater)) {
      columns <- rater_offset + seq_len(ncol(contrasts$Rater))
      W[row_index, columns] <- category * contrasts$Rater[rater_index, ]
    }
    if (ncol(contrasts$Criterion)) {
      columns <- criterion_offset + seq_len(ncol(contrasts$Criterion))
      W[row_index, columns] <- category * contrasts$Criterion[criterion_index, ]
    }
    b_const[[row_index]] <- category * (
      offsets$Rater[[rater_index]] + offsets$Criterion[[criterion_index]]
    )
    if (ncol(step_contrast) && category > 0L) {
      columns <- step_offset + (criterion_index - 1L) * ncol(step_contrast) +
        seq_len(ncol(step_contrast))
      W[row_index, columns] <- colSums(step_contrast[seq_len(category), , drop = FALSE])
    }
  }
  list(response = augmented, weights = weights, prepared = prepared,
       W = W, b_const = b_const)
}

conditional_gradient <- function(object, par) {
  linear <- as.numeric(object$W %*% par + object$b_const)
  pieces <- lapply(seq_len(object$NP), function(pattern) {
    indices <- object$parm_index[[pattern]]
    values <- split(linear[indices], object$pars_info$itemid[indices])
    esf <- psychotools::elementary_symmetric_functions(par = values, order = 1, diff = FALSE)
    W1 <- object$W[indices, , drop = FALSE]
    observed <- object$suffstat[[pattern]] %*% W1
    expected <- -colSums((object$score_freq[[pattern]] * (esf[[2L]] / esf[[1L]])) %*% W1)
    as.numeric(observed) + expected
  })
  Reduce(`+`, pieces)
}
conditional_objective <- function(object, par) {
  linear <- as.numeric(object$W %*% par + object$b_const)
  parts <- vapply(seq_len(object$NP), function(pattern) {
    indices <- object$parm_index[[pattern]]
    values <- split(linear[indices], object$pars_info$itemid[indices])
    esf0 <- psychotools::elementary_symmetric_functions(par = values, order = 0, diff = FALSE)[[1L]]
    -sum(object$suffstat[[pattern]] * linear[indices]) -
      sum(object$score_freq[[pattern]] * log(esf0))
  }, numeric(1))
  -sum(parts)
}
finite_difference_information <- function(object, par) {
  hessian <- vapply(seq_along(par), function(index) {
    step <- 1e-4 * max(1, abs(par[[index]])); plus <- minus <- par
    plus[[index]] <- plus[[index]] + step; minus[[index]] <- minus[[index]] - step
    (conditional_gradient(object, plus) - conditional_gradient(object, minus)) / (2 * step)
  }, numeric(length(par)))
  0.5 * (hessian + t(hessian))
}
finite_difference_gradient <- function(object, par) {
  vapply(seq_along(par), function(index) {
    step <- 1e-4 * max(1, abs(par[[index]])); plus <- minus <- par
    plus[[index]] <- plus[[index]] + step; minus[[index]] <- minus[[index]] - step
    (conditional_objective(object, plus) - conditional_objective(object, minus)) / (2 * step)
  }, numeric(1))
}
rank_evidence <- function(information) {
  singular <- svd(information, nu = 0, nv = 0)$d
  tolerance <- max(1e-9, max(singular) * ncol(information) * 1e-8)
  rank <- sum(singular > tolerance)
  list(rank = rank, nullity = ncol(information) - rank,
       min_eigen = min(eigen(information, symmetric = TRUE, only.values = TRUE)$values),
       tolerance = tolerance)
}
safe_max_abs <- function(values) {
  values <- as.numeric(values); values <- values[is.finite(values)]
  if (!length(values)) NA_real_ else max(abs(values))
}
base_row <- function(case_id, engine, estimator, mode) {
  m <- manifest[manifest$CaseId == case_id, , drop = FALSE]
  p <- python_runs[python_runs$CaseId == case_id, , drop = FALSE]
  data.frame(CaseId = case_id, Mechanism = m$Mechanism,
             ExpectedStatus = m$ExpectedStatus, AnchorLabel = m$AnchorLabel,
             AnchorCount = m$AnchorCount, PythonWorkflowStatus = p$WorkflowStatus,
             Engine = engine, Estimator = estimator, Mode = mode,
             stringsAsFactors = FALSE)
}
rbind_fill <- function(frames) {
  columns <- unique(unlist(lapply(frames, names), use.names = FALSE))
  do.call(rbind, lapply(frames, function(frame) {
    for (name in setdiff(columns, names(frame))) frame[[name]] <- NA
    frame[columns]
  }))
}

immer_rows <- list(); coefficient_rows <- list(); gradient_rows <- list()
for (case_id in manifest$CaseId) {
  data <- ratings[ratings$CaseId == case_id, , drop = FALSE]
  p <- python_runs[python_runs$CaseId == case_id, , drop = FALSE]
  wide <- make_wide(data)
  design <- build_immer_design(wide$response, wide$mapping, 2L, case_anchors(case_id))
  prefit <- design$prepared; prefit$W <- design$W; prefit$b_const <- design$b_const
  zero <- rep(0, ncol(design$W))
  gradient_zero <- conditional_gradient(prefit, zero)
  names(gradient_zero) <- colnames(design$W)
  information_zero <- finite_difference_information(prefit, zero)
  rank_zero <- rank_evidence(information_zero)
  gradient_rows[[length(gradient_rows) + 1L]] <- data.frame(
    CaseId = case_id, Parameter = names(gradient_zero),
    RGradientZero = gradient_zero, RNLLZero = conditional_objective(prefit, zero),
    stringsAsFactors = FALSE
  )
  for (maxit in c(200L, 2000L, 10000L)) {
    row <- base_row(case_id, "immer", "CMLE", paste0("maxit_", maxit))
    row$MaxIt <- maxit; row$KParams <- ncol(design$W)
    row$ParameterNames <- paste(colnames(design$W), collapse = "|")
    row$RPrefitRank <- rank_zero$rank; row$RPrefitNullity <- rank_zero$nullity
    row$RNLLZero <- conditional_objective(prefit, zero)
    row$BConstMaxAbs <- safe_max_abs(design$b_const)
    row$FitAttempted <- TRUE; row$FitReturned <- FALSE
    row$ConvergenceCode <- NA_integer_; row$ConditionalLogLik <- NA_real_
    row$GradientSupNorm <- NA_real_; row$InformationRank <- NA_real_
    row$InformationNullity <- NA_real_; row$AllCoefficientsFinite <- FALSE
    row$AllSEFinite <- FALSE; row$MaxAbsEstimate <- NA_real_; row$MaxAbsSE <- NA_real_
    row$Warnings <- ""; row$Error <- ""
    warnings <- character(); started <- proc.time()[["elapsed"]]
    fit <- tryCatch(withCallingHandlers(
      immer::immer_cml(design$response, weights = design$weights, W = design$W,
                      b_const = design$b_const, par_init = zero, nullcats = "keep",
                      use_rcpp = FALSE,
                      control = list(maxit = maxit, reltol = 1e-12)),
      warning = function(warning) {
        warnings <<- c(warnings, conditionMessage(warning)); invokeRestart("muffleWarning")
      }), error = function(error) error)
    row$ElapsedSeconds <- proc.time()[["elapsed"]] - started
    row$Warnings <- paste(unique(warnings), collapse = " | ")
    if (inherits(fit, "error")) {
      row$Error <- paste0(class(fit)[[1]], ": ", conditionMessage(fit))
      immer_rows[[length(immer_rows) + 1L]] <- row; next
    }
    row$FitReturned <- TRUE
    estimate <- as.numeric(fit$coefficients); names(estimate) <- colnames(fit$W)
    se <- if (!is.null(fit$vcov)) sqrt(pmax(diag(fit$vcov), 0)) else rep(NA_real_, length(estimate))
    gradient <- conditional_gradient(fit, estimate)
    numeric_gradient <- finite_difference_gradient(fit, estimate)
    information <- 0.5 * (fit$result_optim$hessian + t(fit$result_optim$hessian))
    rank <- rank_evidence(information)
    row$ConvergenceCode <- fit$result_optim$convergence
    row$ConditionalLogLik <- fit$loglike; row$GradientSupNorm <- max(abs(gradient))
    row$GradientFiniteDifferenceMaxAbsError <- max(abs(gradient - numeric_gradient))
    row$InformationRank <- rank$rank; row$InformationNullity <- rank$nullity
    row$AllCoefficientsFinite <- all(is.finite(estimate)); row$AllSEFinite <- all(is.finite(se))
    row$MaxAbsEstimate <- safe_max_abs(estimate); row$MaxAbsSE <- safe_max_abs(se)
    coefficient_rows[[length(coefficient_rows) + 1L]] <- data.frame(
      CaseId = case_id, MaxIt = maxit, Parameter = names(estimate),
      Estimate = estimate, SE = se, Gradient = gradient, stringsAsFactors = FALSE
    )
    immer_rows[[length(immer_rows) + 1L]] <- row
  }
}

run_mfrmr <- function(case_id) {
  data <- ratings[ratings$CaseId == case_id, c("Person", "Rater", "Criterion", "Score")]
  anchor_table <- case_anchors(case_id)
  anchor_arg <- if (nrow(anchor_table)) {
    transform(anchor_table, Anchor = Value)[c("Facet", "Level", "Anchor")]
  } else NULL
  row <- base_row(case_id, "mfrmr", "JML", "fit_mfrm_PCM_JML")
  row$FitAttempted <- TRUE; row$FitReturned <- FALSE; row$Converged <- FALSE
  row$InferenceReady <- FALSE; row$AnchorContractPassed <- nrow(anchor_table) == 0L
  row$MaxAbsEstimate <- NA_real_; row$MaxAbsSE <- NA_real_; row$Warnings <- ""; row$Error <- ""
  warnings <- character(); started <- proc.time()[["elapsed"]]
  fit <- tryCatch(withCallingHandlers(
    mfrmr::fit_mfrm(data = data, person = "Person", facets = c("Rater", "Criterion"),
      score = "Score", rating_min = 0, rating_max = 2, keep_original = TRUE,
      model = "PCM", method = "JML", step_facet = "Criterion", anchors = anchor_arg,
      noncenter_facet = "Person", anchor_policy = "warn", min_common_anchors = 1L,
      min_obs_per_element = 1L, min_obs_per_category = 1L,
      maxit = 800L, reltol = 1e-10, optimizer = "BFGS"),
    warning = function(warning) {
      warnings <<- c(warnings, conditionMessage(warning)); invokeRestart("muffleWarning")
    }), error = function(error) error)
  row$ElapsedSeconds <- proc.time()[["elapsed"]] - started
  row$Warnings <- paste(unique(warnings), collapse = " | ")
  if (inherits(fit, "error")) { row$Error <- paste0(class(fit)[[1]], ": ", conditionMessage(fit)); return(row) }
  row$FitReturned <- TRUE
  summary <- as.data.frame(fit$summary)[1, , drop = FALSE]
  get_value <- function(name, default = NA) if (name %in% names(summary)) summary[[name]][[1]] else default
  row$Converged <- isTRUE(as.logical(get_value("Converged", FALSE)))
  row$InferenceReady <- isTRUE(as.logical(get_value("InferenceReady", FALSE)))
  tables <- Filter(Negate(is.null), list(fit$facets$others, fit$steps))
  estimates <- unlist(lapply(tables, function(x) if ("Estimate" %in% names(x)) x$Estimate else numeric()))
  ses <- unlist(lapply(tables, function(x) if ("SE" %in% names(x)) x$SE else numeric()))
  row$MaxAbsEstimate <- safe_max_abs(estimates); row$MaxAbsSE <- safe_max_abs(ses)
  if (nrow(anchor_table)) {
    facets <- as.data.frame(fit$facets$others)
    check <- merge(anchor_table, facets[c("Facet", "Level", "Estimate")],
                   by = c("Facet", "Level"), all.x = TRUE, sort = FALSE)
    row$AnchorMaxAbsDeviation <- max(abs(check$Value - check$Estimate))
    row$AnchorContractPassed <- all(is.finite(check$Estimate)) &&
      all(check$Estimate == check$Value)
  } else row$AnchorMaxAbsDeviation <- 0
  row
}

run_tam <- function(case_id) {
  row <- base_row(case_id, "TAM", "MML", "tam.mml.mfr_PCM")
  anchor_table <- case_anchors(case_id)
  row$FitAttempted <- FALSE; row$FitReturned <- FALSE; row$Converged <- FALSE
  row$MaxAbsEstimate <- NA_real_; row$MaxAbsSE <- NA_real_; row$Warnings <- ""; row$Error <- ""
  if (nrow(anchor_table)) { row$Error <- "anchor_scope_unsupported: equivalent PCM fixed-offset mapping not registered"; return(row) }
  data <- ratings[ratings$CaseId == case_id, , drop = FALSE]
  if (max(data$Score) < 2) { row$Error <- "category_support: declared top category unavailable to TAM"; return(row) }
  wide <- make_person_rater_wide(data); row$FitAttempted <- TRUE
  warnings <- character(); started <- proc.time()[["elapsed"]]
  fit <- tryCatch({
    invisible(utils::capture.output(value <- withCallingHandlers(
      TAM::tam.mml.mfr(resp = wide$response,
        facets = data.frame(rater = wide$grid$Rater), pid = wide$grid$Person,
        formulaA = ~ item + rater + item:step, constraint = "cases",
        control = list(nodes = seq(-6, 6, length.out = 41), maxiter = 400L,
                       convD = 1e-4, conv = 1e-4, convM = 1e-4,
                       Msteps = 4L, progress = FALSE), verbose = FALSE),
      warning = function(warning) {
        warnings <<- c(warnings, conditionMessage(warning)); invokeRestart("muffleWarning")
      }))); value
  }, error = function(error) error)
  row$ElapsedSeconds <- proc.time()[["elapsed"]] - started
  row$Warnings <- paste(unique(warnings), collapse = " | ")
  if (inherits(fit, "error")) { row$Error <- paste0(class(fit)[[1]], ": ", conditionMessage(fit)); return(row) }
  row$FitReturned <- TRUE; row$Converged <- is.finite(fit$iter) && fit$iter < 400L
  parameters <- as.data.frame(fit$xsi.facets %||% data.frame())
  row$MaxAbsEstimate <- safe_max_abs(parameters$xsi %||% numeric())
  row$MaxAbsSE <- safe_max_abs(parameters$se.xsi %||% numeric()); row
}

secondary <- list()
for (case_id in manifest$CaseId) {
  secondary[[length(secondary) + 1L]] <- run_mfrmr(case_id)
  secondary[[length(secondary) + 1L]] <- run_tam(case_id)
}
utils::write.csv(do.call(rbind, immer_rows), file.path(output_dir, "immer_runs.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, coefficient_rows), file.path(output_dir, "immer_coefficients.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, gradient_rows), file.path(output_dir, "immer_prefit_gradient.csv"), row.names = FALSE)
utils::write.csv(rbind_fill(secondary), file.path(output_dir, "secondary_engine_runs.csv"), row.names = FALSE)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else NA_character_
identity <- data.frame(
  PlanSHA256 = args[["plan-sha256"]], AdapterSHA256 = hash_file(normalizePath(script_path)),
  RVersion = as.character(getRversion()), ImmerVersion = as.character(packageVersion("immer")),
  ImmerFunctionSHA256 = function_hash("immer", c("immer_jml", "immer_cml")),
  ImmerLegacyFunctionSHA256 = function_hash("immer", c("immer_jml", "immer_cml"), FALSE),
  MfrmrVersion = as.character(packageVersion("mfrmr")),
  MfrmrFunctionSHA256 = function_hash("mfrmr", "fit_mfrm"),
  MfrmrLegacyFunctionSHA256 = function_hash("mfrmr", "fit_mfrm", FALSE),
  TAMVersion = as.character(packageVersion("TAM")),
  TAMFunctionSHA256 = function_hash("TAM", "tam.mml.mfr"),
  TAMLegacyFunctionSHA256 = function_hash("TAM", "tam.mml.mfr", FALSE),
  SirtVersion = as.character(packageVersion("sirt")),
  SirtFunctionSHA256 = function_hash("sirt", "rm.facets"),
  SirtLegacyFunctionSHA256 = function_hash("sirt", "rm.facets", FALSE),
  MfrmrSourceSHA256 = args[["mfrmr-source-sha256"]], stringsAsFactors = FALSE
)
utils::write.csv(identity, file.path(output_dir, "r_engine_identity.csv"), row.names = FALSE)
