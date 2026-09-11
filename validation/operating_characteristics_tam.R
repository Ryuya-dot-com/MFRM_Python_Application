#!/usr/bin/env Rscript

# Fail-closed TAM MML adapter for the frozen operating-characteristics bundle.

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
missing <- required[!vapply(required, function(name) {
  nzchar(args[[name]] %||% "")
}, logical(1))]
if (length(missing)) {
  stop("Missing required arguments: ", paste(missing, collapse = ", "), call. = FALSE)
}

for (package in c("TAM", "digest")) {
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
truth <- read_input("generated_facet_truth.csv")
anchors <- read_input("generated_anchors.csv")
inventory <- read_input("generated_bundle_files.csv")
bridge_validation <- read_input("bridge_validation_r.csv")
bridge_engines <- read_input("bridge_engine_availability_r.csv")

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
if (!all(hash_ok)) {
  stop("Generated bundle byte hash changed after bridge validation.", call. = FALSE)
}

function_hash <- function(package, functions) {
  namespace <- asNamespace(package)
  hashes <- vapply(functions, function(name) {
    fn <- get(name, envir = namespace, inherits = FALSE)
    digest::digest(
      list(formals = formals(fn), body = body(fn)),
      algo = "sha256",
      serialize = TRUE
    )
  }, character(1))
  digest::digest(hashes, algo = "sha256", serialize = TRUE)
}

tam_namespace <- asNamespace("TAM")
tam_function_hash <- function_hash("TAM", c("tam.mml.mfr", "tam.jml"))
tam_mfr_body_hash <- digest::digest(
  list(
    formals = formals(TAM::tam.mml.mfr),
    body = body(TAM::tam.mml.mfr)
  ),
  algo = "sha256",
  serialize = TRUE
)
progress_function <- get("tam_mml_progress_em", envir = tam_namespace, inherits = FALSE)
progress_body_hash <- digest::digest(
  list(formals = formals(progress_function), body = body(progress_function)),
  algo = "sha256",
  serialize = TRUE
)
registered <- bridge_engines[bridge_engines$Engine == "TAM", , drop = FALSE]
if (nrow(registered) != 1L ||
    !isTRUE(as.logical(registered$Available[[1L]])) ||
    !identical(
      as.character(registered$Version[[1L]]),
      as.character(utils::packageVersion("TAM"))
    ) ||
    !identical(as.character(registered$FunctionSHA256[[1L]]), tam_function_hash)) {
  stop("Loaded TAM code identity differs from the validated R bridge.", call. = FALSE)
}

modes <- data.frame(
  Mode = c("TAM_MML_Q21_SENSITIVITY", "TAM_MML_Q61_PRIMARY"),
  QuadraturePoints = c(21L, 61L),
  PrimaryMode = c(FALSE, TRUE),
  stringsAsFactors = FALSE
)
maxiter <- 1200L
conv_deviance <- 1e-6
conv_parameter <- 1e-5
conv_mstep <- 1e-6
msteps <- 20L
tail_mass_tolerance <- 1e-6
minimum_variance <- 0.001
variance_boundary_tolerance <- 1e-8

prepare_wide <- function(data, rating_max) {
  key <- paste(data$Person, data$Rater, data$Task, data$Criterion, sep = "\r")
  if (anyDuplicated(key)) stop("Duplicate Person x Rater x Task x Criterion cells.")
  grid <- unique(data[c("Person", "Rater", "Task")])
  grid <- grid[order(grid$Person, grid$Rater, grid$Task), , drop = FALSE]
  items <- sort(unique(as.character(data$Criterion)))
  response <- matrix(
    NA_integer_,
    nrow = nrow(grid),
    ncol = length(items),
    dimnames = list(NULL, items)
  )
  row_index <- match(
    paste(data$Person, data$Rater, data$Task),
    paste(grid$Person, grid$Rater, grid$Task)
  )
  column_index <- match(data$Criterion, items)
  response[cbind(row_index, column_index)] <- as.integer(data$Score)
  retained <- rowSums(!is.na(response)) > 0L
  response <- response[retained, , drop = FALSE]
  grid <- grid[retained, , drop = FALSE]
  item_counts <- table(
    factor(data$Criterion, levels = items),
    factor(data$Score, levels = 0:rating_max)
  )
  generalized_counts <- table(
    factor(data$Rater, levels = sort(unique(data$Rater))),
    factor(data$Task, levels = sort(unique(data$Task))),
    factor(data$Criterion, levels = items)
  )
  list(
    response = as.data.frame(response, stringsAsFactors = FALSE),
    grid = grid,
    items = items,
    declared_support_complete = all(apply(response, 2L, max, na.rm = TRUE) == rating_max),
    zero_item_categories = sum(item_counts == 0),
    minimum_item_category_count = min(item_counts),
    zero_generalized_cells = sum(generalized_counts == 0),
    minimum_generalized_cell_count = min(generalized_counts),
    minimum_person_rows = min(table(grid$Person)),
    minimum_rater_rows = min(table(grid$Rater))
  )
}

last_deviance_change <- function(history) {
  history <- as.matrix(history)
  if (nrow(history) < 2L) return(NA_real_)
  values <- as.numeric(history[, 2L])
  values <- values[is.finite(values)]
  if (length(values) < 2L) return(NA_real_)
  abs(values[[length(values)]] - values[[length(values) - 1L]])
}

last_progress_value <- function(output, label) {
  matches <- grep(paste0("^  ", label, ":"), output, value = TRUE)
  if (!length(matches)) return(NA_real_)
  value <- sub(paste0("^  ", label, ":\\s*"), "", matches[[length(matches)]])
  suppressWarnings(as.numeric(value))
}

build_anchor_fixed <- function(wide, run_anchors) {
  if (!nrow(run_anchors)) return(NULL)
  if (!all(run_anchors$Facet == "Rater")) {
    stop("TAM adapter supports Rater anchors only.")
  }
  output <- utils::capture.output(
    design <- suppressWarnings(TAM::tam.mml.mfr(
      resp = wide$response,
      facets = data.frame(
        rater = wide$grid$Rater,
        task = wide$grid$Task,
        stringsAsFactors = FALSE
      ),
      pid = wide$grid$Person,
      formulaA = ~ item + rater + task + step,
      constraint = "cases",
      control = list(maxiter = 2L, progress = FALSE),
      verbose = FALSE
    ))
  )
  labels <- rownames(design$xsi)
  target <- paste0("rater", as.character(run_anchors$Level))
  index <- match(target, labels)
  if (anyNA(index)) {
    stop("One or more requested TAM anchors are derived or absent from the free xi vector.")
  }
  fixed <- cbind(index, as.numeric(run_anchors$Anchor))
  storage.mode(fixed) <- "double"
  fixed
}

base_run <- function(manifest_row, mode_row, data, wide, run_anchors) {
  data.frame(
    SchemaVersion = as.character(manifest_row$SchemaVersion),
    RunId = as.character(manifest_row$RunId),
    ConditionId = as.character(manifest_row$ConditionId),
    Design = as.character(manifest_row$Design),
    TruthBias = as.numeric(manifest_row$TruthBias),
    TruthPositive = as.logical(manifest_row$TruthPositive),
    Replicate = as.integer(manifest_row$Replicate),
    Seed = as.numeric(manifest_row$Seed),
    Engine = "TAM",
    Estimator = "MML",
    Mode = as.character(mode_row$Mode),
    PrimaryMode = as.logical(mode_row$PrimaryMode),
    QuadraturePoints = as.integer(mode_row$QuadraturePoints),
    ThetaMinimum = -8,
    ThetaMaximum = 8,
    Model = "additive RSM: Criterion item + Rater + Task + common step",
    Estimand = paste(
      "MML additive RSM facet severities and common steps;",
      "estimated normal Person distribution"
    ),
    Rows = nrow(data),
    PersonTaskRaterRows = nrow(wide$grid),
    Persons = length(unique(data$Person)),
    Raters = length(unique(data$Rater)),
    Tasks = length(unique(data$Task)),
    Criteria = length(unique(data$Criterion)),
    RequestedAnchorRows = nrow(run_anchors),
    AnchorPreflightReturned = !nrow(run_anchors),
    AnchorContractPassed = NA,
    AnchorMaxAbsDeviation = NA_real_,
    BiasEstimandSupported = FALSE,
    DeclaredSupportComplete = wide$declared_support_complete,
    ZeroItemCategories = wide$zero_item_categories,
    MinimumItemCategoryCount = wide$minimum_item_category_count,
    ZeroGeneralizedCells = wide$zero_generalized_cells,
    MinimumGeneralizedCellCount = wide$minimum_generalized_cell_count,
    MinimumPersonRows = wide$minimum_person_rows,
    MinimumRaterRows = wide$minimum_rater_rows,
    FitReturned = FALSE,
    LoopExitedBeforeCap = FALSE,
    TerminalProgressParsed = FALSE,
    Converged = FALSE,
    InferenceReady = FALSE,
    AnalysisEligible = FALSE,
    FailureStage = "fit",
    FailureReason = "",
    Warnings = "",
    Iterations = NA_real_,
    FinalAbsoluteDevianceChange = NA_real_,
    TerminalItemParameterChange = NA_real_,
    TerminalRegressionParameterChange = NA_real_,
    TerminalVarianceParameterChange = NA_real_,
    TerminalMaxParameterChange = NA_real_,
    DevianceTolerance = conv_deviance,
    ParameterTolerance = conv_parameter,
    MStepTolerance = conv_mstep,
    MSteps = msteps,
    ConvergenceBasis = paste(
      "TAM loop exit plus parsed terminal variance/regression progress;",
      "item/deviance loop criteria are also checked from retained output/history"
    ),
    ProgressPrecisionBoundary = paste(
      "terminal parameter changes are parsed from TAM progress text;",
      "the source loop supplies the authoritative primary item/deviance decision"
    ),
    Deviance = NA_real_,
    LogLik = NA_real_,
    PopulationMean = 0,
    PopulationVariance = NA_real_,
    PopulationSD = NA_real_,
    MinimumPopulationVariance = minimum_variance,
    PopulationVarianceBoundaryTolerance = variance_boundary_tolerance,
    PopulationVarianceAtLowerBound = FALSE,
    QuadratureTailMass = NA_real_,
    QuadratureTailMassTolerance = tail_mass_tolerance,
    FacetSEReady = FALSE,
    StepSEReady = FALSE,
    EAPAvailable = FALSE,
    EAPReliability = NA_real_,
    Unestimable99Detected = FALSE,
    ReportedNParameters = NA_real_,
    AIC = NA_real_,
    BIC = NA_real_,
    InformationCriteriaComparableAcrossEstimators = FALSE,
    InformationCriteriaBoundary = paste(
      "retained for within-TAM audit only; MML likelihood and parameter count",
      "are not compared with Python JMLE, CMLE, or sirt"
    ),
    DerivedFacetSEBoundary = paste(
      "TAM 4.3-25 expands case-constrained SEs from a diagonal matrix of free-xi SEs;",
      "coverage is withheld for derived last levels"
    ),
    ElapsedSeconds = NA_real_,
    BundleInventorySHA256 = bundle_inventory_sha256,
    TAMVersion = as.character(utils::packageVersion("TAM")),
    TAMFunctionSHA256 = tam_function_hash,
    TAMMFRFunctionBodySHA256 = tam_mfr_body_hash,
    TAMProgressFunctionBodySHA256 = progress_body_hash,
    stringsAsFactors = FALSE
  )
}

run_rows <- list()
facet_rows <- list()
step_rows <- list()
person_rows <- list()
surface_rows <- list()

for (manifest_index in seq_len(nrow(manifest))) {
  manifest_row <- manifest[manifest_index, , drop = FALSE]
  run_id <- as.character(manifest_row$RunId)
  data <- ratings[
    ratings$RunId == run_id,
    c("Person", "Rater", "Task", "Criterion", "Score"),
    drop = FALSE
  ]
  run_truth <- truth[
    truth$RunId == run_id,
    c("Facet", "Level", "Truth"),
    drop = FALSE
  ]
  run_anchors <- anchors[
    anchors$RunId == run_id,
    c("Facet", "Level", "Anchor"),
    drop = FALSE
  ]
  wide <- tryCatch(
    prepare_wide(data, as.integer(manifest_row$Categories) - 1L),
    error = function(error) error
  )

  for (mode_index in seq_len(nrow(modes))) {
    mode_row <- modes[mode_index, , drop = FALSE]
    if (inherits(wide, "error")) {
      placeholder <- list(
        grid = data.frame(), declared_support_complete = FALSE,
        zero_item_categories = NA, minimum_item_category_count = NA,
        zero_generalized_cells = NA, minimum_generalized_cell_count = NA,
        minimum_person_rows = NA, minimum_rater_rows = NA
      )
      run <- base_run(manifest_row, mode_row, data, placeholder, run_anchors)
      run$FailureStage <- "input"
      run$FailureReason <- conditionMessage(wide)
      run_rows[[length(run_rows) + 1L]] <- run
      next
    }
    run <- base_run(manifest_row, mode_row, data, wide, run_anchors)
    if (!isTRUE(wide$declared_support_complete)) {
      run$FailureStage <- "category_support"
      run$FailureReason <- paste(
        "TAM infers Criterion support from observed maxima and this adapter",
        "does not alter the registered category map"
      )
      run_rows[[length(run_rows) + 1L]] <- run
      next
    }

    fixed <- tryCatch(
      build_anchor_fixed(wide, run_anchors),
      error = function(error) error
    )
    if (inherits(fixed, "error")) {
      run$FailureStage <- "anchors"
      run$FailureReason <- conditionMessage(fixed)
      run_rows[[length(run_rows) + 1L]] <- run
      next
    }
    run$AnchorPreflightReturned <- TRUE

    warning_messages <- character()
    console_output <- character()
    started <- proc.time()[["elapsed"]]
    fit <- tryCatch({
      console_output <- utils::capture.output(
        value <- withCallingHandlers(
          TAM::tam.mml.mfr(
            resp = wide$response,
            facets = data.frame(
              rater = wide$grid$Rater,
              task = wide$grid$Task,
              stringsAsFactors = FALSE
            ),
            pid = wide$grid$Person,
            formulaA = ~ item + rater + task + step,
            constraint = "cases",
            xsi.fixed = fixed,
            control = list(
              nodes = seq(-8, 8, length.out = as.integer(mode_row$QuadraturePoints)),
              maxiter = maxiter,
              convD = conv_deviance,
              conv = conv_parameter,
              convM = conv_mstep,
              Msteps = msteps,
              min.variance = minimum_variance,
              dev_crit = "absolute",
              progress = TRUE
            ),
            verbose = FALSE
          ),
          warning = function(warning) {
            warning_messages <<- c(warning_messages, conditionMessage(warning))
            invokeRestart("muffleWarning")
          }
        )
      )
      value
    }, error = function(error) error)
    run$ElapsedSeconds <- proc.time()[["elapsed"]] - started
    run$Warnings <- paste(unique(warning_messages), collapse = " | ")
    if (inherits(fit, "error")) {
      run$FailureReason <- paste0(class(fit)[[1L]], ": ", conditionMessage(fit))
      run_rows[[length(run_rows) + 1L]] <- run
      next
    }

    run$FitReturned <- TRUE
    run$Iterations <- as.numeric(fit$iter)
    run$LoopExitedBeforeCap <- is.finite(run$Iterations) && run$Iterations < maxiter
    run$FinalAbsoluteDevianceChange <- last_deviance_change(fit$deviance.history)
    run$TerminalItemParameterChange <- last_progress_value(
      console_output,
      "Maximum item intercept parameter change"
    )
    run$TerminalRegressionParameterChange <- last_progress_value(
      console_output,
      "Maximum regression parameter change"
    )
    run$TerminalVarianceParameterChange <- last_progress_value(
      console_output,
      "Maximum variance parameter change"
    )
    terminal_values <- c(
      run$TerminalItemParameterChange,
      run$TerminalRegressionParameterChange,
      run$TerminalVarianceParameterChange
    )
    run$TerminalProgressParsed <- all(is.finite(terminal_values))
    run$TerminalMaxParameterChange <- if (run$TerminalProgressParsed) {
      max(terminal_values)
    } else NA_real_
    run$Deviance <- as.numeric(fit$deviance)
    run$LogLik <- -0.5 * run$Deviance
    run$PopulationVariance <- as.numeric(fit$variance[1L, 1L])
    run$PopulationSD <- sqrt(run$PopulationVariance)
    run$PopulationVarianceAtLowerBound <- is.finite(run$PopulationVariance) &&
      run$PopulationVariance <= minimum_variance + variance_boundary_tolerance
    post <- as.matrix(fit$post)
    run$QuadratureTailMass <- mean(rowSums(post[, c(1L, ncol(post)), drop = FALSE]))
    run$EAPAvailable <- nrow(as.data.frame(fit$person)) > 0L &&
      all(is.finite(as.numeric(fit$person$EAP)))
    run$EAPReliability <- as.numeric(fit$EAP.rel)
    run$ReportedNParameters <- as.numeric(fit$ic$np)
    run$AIC <- as.numeric(fit$ic$AIC)
    run$BIC <- as.numeric(fit$ic$BIC)

    facet_table <- as.data.frame(fit$xsi.facets, stringsAsFactors = FALSE)
    constraint_table <- as.data.frame(fit$xsi.constr$xsi.table, stringsAsFactors = FALSE)
    facet_table$DerivedConstraint <- as.logical(
      constraint_table$constraint[match(facet_table$parameter, constraint_table$parameter)]
    )
    fixed_labels <- if (is.null(fixed)) character() else rownames(fit$xsi)[fixed[, 1L]]
    facet_table$Fixed <- facet_table$parameter %in% fixed_labels
    run$Unestimable99Detected <- any(
      is.finite(facet_table$xsi) & abs(facet_table$xsi - 99) <= 1e-8
    )
    free_facet <- facet_table$facet %in% c("item", "rater", "task") &
      !facet_table$Fixed & !facet_table$DerivedConstraint
    free_step <- facet_table$facet == "step" &
      !facet_table$Fixed & !facet_table$DerivedConstraint
    run$FacetSEReady <- any(free_facet) && all(
      is.finite(facet_table$se.xsi[free_facet]) & facet_table$se.xsi[free_facet] > 0
    )
    run$StepSEReady <- any(free_step) && all(
      is.finite(facet_table$se.xsi[free_step]) & facet_table$se.xsi[free_step] > 0
    )

    if (nrow(run_anchors)) {
      estimated_anchor <- facet_table$xsi[
        match(paste0("rater", run_anchors$Level), facet_table$parameter)
      ]
      deviations <- abs(estimated_anchor - as.numeric(run_anchors$Anchor))
      run$AnchorMaxAbsDeviation <- max(deviations)
      run$AnchorContractPassed <- all(is.finite(deviations)) && all(deviations <= 1e-10)
    } else {
      run$AnchorMaxAbsDeviation <- 0
      run$AnchorContractPassed <- TRUE
    }

    run$Converged <- run$LoopExitedBeforeCap &&
      run$TerminalProgressParsed &&
      is.finite(run$FinalAbsoluteDevianceChange) &&
      run$FinalAbsoluteDevianceChange <= conv_deviance &&
      run$TerminalItemParameterChange <= conv_parameter &&
      run$TerminalRegressionParameterChange <= conv_parameter &&
      run$TerminalVarianceParameterChange <= conv_parameter
    run$InferenceReady <- run$Converged &&
      run$FacetSEReady && run$StepSEReady && run$EAPAvailable &&
      !run$Unestimable99Detected &&
      !run$PopulationVarianceAtLowerBound &&
      is.finite(run$QuadratureTailMass) &&
      run$QuadratureTailMass <= tail_mass_tolerance
    run$AnalysisEligible <- run$InferenceReady && run$AnchorContractPassed
    if (run$AnalysisEligible) {
      run$FailureStage <- ""
      run$FailureReason <- ""
    } else if (!run$Converged) {
      run$FailureStage <- "convergence"
      run$FailureReason <- paste(
        "TAM loop exit and independently retained terminal progress",
        "did not satisfy the adapter convergence contract"
      )
    } else if (run$Unestimable99Detected) {
      run$FailureStage <- "estimability"
      run$FailureReason <- "TAM returned one or more documented unestimable facet values at 99"
    } else if (run$PopulationVarianceAtLowerBound) {
      run$FailureStage <- "population_variance"
      run$FailureReason <- "estimated Person variance is on the configured TAM lower bound"
    } else if (!run$FacetSEReady || !run$StepSEReady) {
      run$FailureStage <- "uncertainty"
      run$FailureReason <- "one or more free facet/common-step SEs are unavailable"
    } else if (!run$EAPAvailable) {
      run$FailureStage <- "person_scores"
      run$FailureReason <- "finite EAP scores unavailable"
    } else if (run$QuadratureTailMass > tail_mass_tolerance) {
      run$FailureStage <- "quadrature"
      run$FailureReason <- "mean posterior endpoint mass exceeds review tolerance"
    } else if (!run$AnchorContractPassed) {
      run$FailureStage <- "anchors"
      run$FailureReason <- "fixed Rater estimates differ from supplied anchors"
    }

    facet_map <- list(item = "Criterion", rater = "Rater", task = "Task")
    for (tam_facet in names(facet_map)) {
      output <- facet_table[facet_table$facet == tam_facet, , drop = FALSE]
      output$Facet <- facet_map[[tam_facet]]
      output$Level <- if (tam_facet == "item") {
        as.character(output$parameter)
      } else {
        sub(paste0("^", tam_facet), "", as.character(output$parameter))
      }
      facet_truth <- run_truth[
        run_truth$Facet == facet_map[[tam_facet]],
        c("Level", "Truth"),
        drop = FALSE
      ]
      output <- merge(facet_truth, output, by = "Level", all.x = TRUE, sort = FALSE)
      output$Estimate <- as.numeric(output$xsi)
      output$SE <- as.numeric(output$se.xsi)
      output$Anchored <- as.logical(output$Fixed)
      output$EstimateAligned <- output$Estimate
      output$TruthAligned <- as.numeric(output$Truth)
      output$ComparisonScale <- "anchor_identified_absolute"
      if (!(tam_facet == "rater" && nrow(run_anchors))) {
        output$EstimateAligned <- output$Estimate - mean(output$Estimate)
        output$TruthAligned <- output$Truth - mean(output$Truth)
        output$ComparisonScale <- "mean_aligned_location"
      }
      output$ErrorAligned <- output$EstimateAligned - output$TruthAligned
      output$IncludedInSummary <- run$AnalysisEligible
      output$CoverageEligible <- run$AnalysisEligible &
        !output$Anchored & !output$DerivedConstraint &
        is.finite(output$SE) & output$SE > 0
      output$RunId <- run_id
      output$ConditionId <- as.character(manifest_row$ConditionId)
      output$Design <- as.character(manifest_row$Design)
      output$TruthBias <- as.numeric(manifest_row$TruthBias)
      output$Replicate <- as.integer(manifest_row$Replicate)
      output$Seed <- as.numeric(manifest_row$Seed)
      output$Engine <- "TAM"
      output$Estimator <- "MML"
      output$Mode <- as.character(mode_row$Mode)
      facet_rows[[length(facet_rows) + 1L]] <- output[c(
        "RunId", "ConditionId", "Engine", "Estimator", "Mode", "Design",
        "TruthBias", "Replicate", "Seed", "Facet", "Level", "Truth",
        "Estimate", "SE", "EstimateAligned", "TruthAligned", "ErrorAligned",
        "ComparisonScale", "Anchored", "DerivedConstraint",
        "IncludedInSummary", "CoverageEligible"
      )]
    }

    step_output <- facet_table[facet_table$facet == "step", , drop = FALSE]
    step_output$Step <- as.integer(sub("^step", "", step_output$parameter))
    step_output$Estimate <- as.numeric(step_output$xsi)
    step_output$SE <- as.numeric(step_output$se.xsi)
    step_output$RunId <- run_id
    step_output$ConditionId <- as.character(manifest_row$ConditionId)
    step_output$Design <- as.character(manifest_row$Design)
    step_output$TruthBias <- as.numeric(manifest_row$TruthBias)
    step_output$Replicate <- as.integer(manifest_row$Replicate)
    step_output$Seed <- as.numeric(manifest_row$Seed)
    step_output$Engine <- "TAM"
    step_output$Estimator <- "MML"
    step_output$Mode <- as.character(mode_row$Mode)
    step_output$IncludedInSummary <- run$AnalysisEligible
    step_output$CoverageEligible <- FALSE
    step_output$UncertaintyBoundary <- paste(
      "common-step truth is not registered in the current recovery bundle;",
      "derived last-step covariance is unavailable"
    )
    step_rows[[length(step_rows) + 1L]] <- step_output[c(
      "RunId", "ConditionId", "Engine", "Estimator", "Mode", "Design",
      "TruthBias", "Replicate", "Seed", "Step", "Estimate", "SE",
      "DerivedConstraint", "IncludedInSummary", "CoverageEligible",
      "UncertaintyBoundary"
    )]

    persons <- as.data.frame(fit$person, stringsAsFactors = FALSE)
    person_rows[[length(person_rows) + 1L]] <- data.frame(
      RunId = run_id,
      ConditionId = as.character(manifest_row$ConditionId),
      Engine = "TAM",
      Estimator = "MML",
      Mode = as.character(mode_row$Mode),
      Design = as.character(manifest_row$Design),
      TruthBias = as.numeric(manifest_row$TruthBias),
      Replicate = as.integer(manifest_row$Replicate),
      Seed = as.numeric(manifest_row$Seed),
      Person = as.character(persons$pid),
      Estimate = as.numeric(persons$EAP),
      SE = as.numeric(persons$SD.EAP),
      ParameterStatus = "EAP under estimated TAM normal Person distribution",
      IncludedInSummary = run$AnalysisEligible,
      stringsAsFactors = FALSE
    )

    item_names <- rownames(fit$A)
    if (is.null(item_names) || length(item_names) != nrow(fit$AXsi_)) {
      stop("TAM generalized-item labels do not align with AXsi_ rows.", call. = FALSE)
    }
    item_map <- data.frame(
      GeneralizedItem = item_names,
      Criterion = sub("-rater.*$", "", item_names),
      Rater = sub("^.*-rater([^ -]+)-task.*$", "\\1", item_names),
      Task = sub("^.*-task", "", item_names),
      stringsAsFactors = FALSE
    )
    for (category in seq_len(ncol(fit$AXsi_) - 1L)) {
      surface_rows[[length(surface_rows) + 1L]] <- data.frame(
        RunId = run_id,
        ConditionId = as.character(manifest_row$ConditionId),
        Engine = "TAM",
        Estimator = "MML",
        Mode = as.character(mode_row$Mode),
        Design = as.character(manifest_row$Design),
        TruthBias = as.numeric(manifest_row$TruthBias),
        Replicate = as.integer(manifest_row$Replicate),
        Seed = as.numeric(manifest_row$Seed),
        GeneralizedItem = item_map$GeneralizedItem,
        Rater = item_map$Rater,
        Task = item_map$Task,
        Criterion = item_map$Criterion,
        Category = category,
        CumulativeDifficulty = as.numeric(fit$AXsi_[, category + 1L]),
        IncludedInSummary = run$AnalysisEligible,
        stringsAsFactors = FALSE
      )
    }
    run_rows[[length(run_rows) + 1L]] <- run
  }
}

runs <- do.call(rbind, run_rows)
facet_recovery <- do.call(rbind, facet_rows)
step_estimates <- do.call(rbind, step_rows)
person_estimates <- do.call(rbind, person_rows)
surface_estimates <- do.call(rbind, surface_rows)
utils::write.csv(runs, file.path(output_dir, "tam_runs.csv"), row.names = FALSE)
utils::write.csv(
  facet_recovery,
  file.path(output_dir, "tam_facet_recovery.csv"),
  row.names = FALSE
)
utils::write.csv(
  step_estimates,
  file.path(output_dir, "tam_step_estimates.csv"),
  row.names = FALSE
)
utils::write.csv(
  person_estimates,
  file.path(output_dir, "tam_person_estimates.csv"),
  row.names = FALSE
)
utils::write.csv(
  surface_estimates,
  file.path(output_dir, "tam_surface_estimates.csv"),
  row.names = FALSE
)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) {
  sub("^--file=", "", script_arg[[1L]])
} else NA_character_
identity <- data.frame(
  BundleInventorySHA256 = bundle_inventory_sha256,
  TAMVersion = as.character(utils::packageVersion("TAM")),
  TAMFunctionSHA256 = tam_function_hash,
  TAMMFRFunctionBodySHA256 = tam_mfr_body_hash,
  TAMProgressFunctionBodySHA256 = progress_body_hash,
  AdapterScriptSHA256 = if (!is.na(script_path)) {
    digest::digest(
      file = normalizePath(script_path),
      algo = "sha256",
      serialize = FALSE
    )
  } else NA_character_,
  Modes = paste(modes$Mode, collapse = "|"),
  Model = "additive RSM: Criterion item + Rater + Task + common step",
  PersonDistribution = "estimated normal distribution on fixed -8..8 grid",
  AnchorScope = paste(
    "Rater free-xi values fixed; requested derived anchor levels fail closed;",
    "derived constraint SE not used for coverage"
  ),
  BiasScope = "local Rater x Task interaction unsupported",
  ConvergenceScope = paste(
    "TAM loop exit supplemented by parsed terminal variance/regression progress;",
    "progress formatter identity retained"
  ),
  InformationCriteriaScope = "within-TAM audit only; no cross-estimator comparison",
  stringsAsFactors = FALSE
)
utils::write.csv(
  identity,
  file.path(output_dir, "tam_adapter_identity.csv"),
  row.names = FALSE
)
