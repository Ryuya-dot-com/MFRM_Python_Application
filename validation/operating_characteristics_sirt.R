#!/usr/bin/env Rscript

# Fit sirt::rm.facets() to the byte-validated OC bundle.  This is an MML
# sensitivity lane, not JMLE parity: Task x Criterion combinations become six
# virtual items with item-specific PCM thresholds and the Person distribution
# identifies otherwise disconnected rater groups.

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
for (package in c("sirt", "digest")) {
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
if (!all(hash_ok)) stop("Generated bundle byte hash changed after bridge validation.", call. = FALSE)

function_body_hash <- digest::digest(
  list(
    formals = formals(sirt::rm.facets),
    body = body(sirt::rm.facets)
  ),
  algo = "sha256",
  serialize = TRUE
)
function_hash <- digest::digest(
  c(rm.facets = function_body_hash),
  algo = "sha256",
  serialize = TRUE
)
registered <- bridge_engines[bridge_engines$Engine == "sirt", , drop = FALSE]
if (nrow(registered) != 1L || !isTRUE(as.logical(registered$Available[[1L]])) ||
    !identical(as.character(registered$Version[[1L]]), as.character(utils::packageVersion("sirt"))) ||
    !identical(as.character(registered$FunctionSHA256[[1L]]), function_hash)) {
  stop("Loaded sirt code identity differs from the validated R bridge.", call. = FALSE)
}

modes <- data.frame(
  Mode = c("SIRT_MML_Q30_SENSITIVITY", "SIRT_MML_Q61_PRIMARY"),
  QuadraturePoints = c(30L, 61L),
  PrimaryMode = c(FALSE, TRUE),
  stringsAsFactors = FALSE
)
maxiter <- 1200L
globconv <- 1e-4
maxdevchange <- 1e-6
tail_mass_tolerance <- 1e-6

prepare_wide <- function(data, rating_max) {
  data$Item <- paste(data$Task, data$Criterion, sep = "__")
  response_key <- paste(data$Person, data$Rater, data$Item, sep = "\r")
  if (anyDuplicated(response_key)) stop("Duplicate Person x Rater x virtual-item cells.")
  grid <- unique(data[c("Person", "Rater")])
  grid <- grid[order(grid$Person, grid$Rater), , drop = FALSE]
  items <- sort(unique(as.character(data$Item)))
  response <- matrix(
    NA_real_, nrow = nrow(grid), ncol = length(items),
    dimnames = list(NULL, items)
  )
  row_index <- match(paste(data$Person, data$Rater), paste(grid$Person, grid$Rater))
  column_index <- match(data$Item, items)
  response[cbind(row_index, column_index)] <- as.numeric(data$Score)

  counts <- table(
    factor(data$Item, levels = items),
    factor(data$Score, levels = 0:rating_max)
  )
  item_maxima <- apply(response, 2L, max, na.rm = TRUE)
  person_rater_counts <- table(factor(grid$Rater, levels = sort(unique(grid$Rater))))
  item_map <- unique(data[c("Item", "Task", "Criterion")])
  item_map <- item_map[match(items, item_map$Item), , drop = FALSE]
  list(
    response = response,
    grid = grid,
    items = items,
    item_map = item_map,
    item_category_counts = counts,
    declared_support_complete = all(item_maxima == rating_max),
    zero_item_categories = sum(counts == 0),
    minimum_item_category_count = min(counts),
    minimum_person_rater_rows = min(person_rater_counts)
  )
}

last_relative_deviance_change <- function(history) {
  values <- as.numeric(history[is.finite(history)])
  if (length(values) < 2L) return(NA_real_)
  abs((values[[length(values)]] - values[[length(values) - 1L]]) /
        values[[length(values) - 1L]])
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
    Engine = "sirt",
    Estimator = "MML",
    Mode = as.character(mode_row$Mode),
    PrimaryMode = as.logical(mode_row$PrimaryMode),
    QuadraturePoints = as.integer(mode_row$QuadraturePoints),
    ThetaMinimum = -8,
    ThetaMaximum = 8,
    Model = "virtual-item PCM with additive rater severity",
    Estimand = paste(
      "MML rater severity and Task x Criterion item-specific thresholds;",
      "estimated normal Person distribution"
    ),
    Rows = nrow(data),
    PersonRaterRows = nrow(wide$grid),
    Persons = length(unique(data$Person)),
    Raters = length(unique(data$Rater)),
    VirtualItems = length(wide$items),
    RequestedAnchorRows = nrow(run_anchors),
    AnchorContractPassed = NA,
    AnchorMaxAbsDeviation = NA_real_,
    BiasEstimandSupported = FALSE,
    DeclaredSupportComplete = wide$declared_support_complete,
    ZeroItemCategories = wide$zero_item_categories,
    MinimumItemCategoryCount = wide$minimum_item_category_count,
    MinimumPersonRaterRows = wide$minimum_person_rater_rows,
    FitReturned = FALSE,
    Converged = FALSE,
    InferenceReady = FALSE,
    AnalysisEligible = FALSE,
    FailureStage = "fit",
    FailureReason = "",
    Warnings = "",
    Iterations = NA_real_,
    FinalRelativeDevianceChange = NA_real_,
    GlobalParameterTolerance = globconv,
    RelativeDevianceTolerance = maxdevchange,
    ConvergenceBasis = paste(
      "rm.facets source-loop exit before maxiter implies both internal",
      "global parameter and relative-deviance criteria; only deviance history is returned"
    ),
    Deviance = NA_real_,
    LogLik = NA_real_,
    PopulationMean = NA_real_,
    PopulationSD = NA_real_,
    QuadratureTailMass = NA_real_,
    QuadratureTailMassTolerance = tail_mass_tolerance,
    RaterSEReady = FALSE,
    ItemThresholdSEReady = FALSE,
    EAPAvailable = FALSE,
    EAPReliability = NA_real_,
    ReportedNParameters = NA_real_,
    InformationCriteriaComparable = FALSE,
    InformationCriteriaBoundary = paste(
      "not compared: sirt 4.2-133 rm_facets_ic subtracts the numeric",
      "b.rater.center code from RR; unanchored center mode 2 reports RR-2"
    ),
    ElapsedSeconds = NA_real_,
    BundleInventorySHA256 = bundle_inventory_sha256,
    SirtVersion = as.character(utils::packageVersion("sirt")),
    SirtFunctionSHA256 = function_hash,
    SirtFunctionBodySHA256 = function_body_hash,
    stringsAsFactors = FALSE
  )
}

run_rows <- list()
rater_rows <- list()
item_rows <- list()
person_rows <- list()
for (manifest_index in seq_len(nrow(manifest))) {
  manifest_row <- manifest[manifest_index, , drop = FALSE]
  run_id <- as.character(manifest_row$RunId)
  data <- ratings[
    ratings$RunId == run_id,
    c("Person", "Rater", "Task", "Criterion", "Score"),
    drop = FALSE
  ]
  run_truth <- truth[truth$RunId == run_id, c("Facet", "Level", "Truth"), drop = FALSE]
  run_anchors <- anchors[anchors$RunId == run_id, c("Facet", "Level", "Anchor"), drop = FALSE]
  wide <- tryCatch(
    prepare_wide(data, as.integer(manifest_row$Categories) - 1L),
    error = function(error) error
  )

  for (mode_index in seq_len(nrow(modes))) {
    mode_row <- modes[mode_index, , drop = FALSE]
    if (inherits(wide, "error")) {
      placeholder <- list(
        grid = data.frame(), items = character(), declared_support_complete = FALSE,
        zero_item_categories = NA, minimum_item_category_count = NA,
        minimum_person_rater_rows = NA
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
        "sirt infers virtual-item support from observed maxima and cannot",
        "receive a zero-weight declared-support sentinel"
      )
      run_rows[[length(run_rows) + 1L]] <- run
      next
    }

    rater_levels <- sort(unique(as.character(wide$grid$Rater)))
    fixed <- NULL
    if (nrow(run_anchors)) {
      if (!all(run_anchors$Facet == "Rater")) {
        run$FailureStage <- "anchors"
        run$FailureReason <- "sirt adapter supports Rater anchors only"
        run_rows[[length(run_rows) + 1L]] <- run
        next
      }
      fixed <- rep(NA_real_, length(rater_levels))
      anchor_index <- match(as.character(run_anchors$Level), rater_levels)
      if (anyNA(anchor_index)) {
        run$FailureStage <- "anchors"
        run$FailureReason <- "anchor level not present in sirt rater levels"
        run_rows[[length(run_rows) + 1L]] <- run
        next
      }
      fixed[anchor_index] <- as.numeric(run_anchors$Anchor)
    }

    warning_messages <- character()
    console_output <- character()
    started <- proc.time()[["elapsed"]]
    fit <- tryCatch(
      {
        console_output <- utils::capture.output(
          value <- withCallingHandlers(
            sirt::rm.facets(
              dat = wide$response,
              pid = wide$grid$Person,
              rater = wide$grid$Rater,
              theta.k = seq(-8, 8, length.out = as.integer(mode_row$QuadraturePoints)),
              est.b.rater = TRUE,
              est.a.item = FALSE,
              est.a.rater = FALSE,
              rater_item_int = FALSE,
              est.mean = FALSE,
              b.rater.fixed = fixed,
              b.rater.center = 2,
              globconv = globconv,
              maxdevchange = maxdevchange,
              maxiter = maxiter
            ),
            warning = function(warning) {
              warning_messages <<- c(warning_messages, conditionMessage(warning))
              invokeRestart("muffleWarning")
            }
          )
        )
        value
      },
      error = function(error) error
    )
    run$ElapsedSeconds <- proc.time()[["elapsed"]] - started
    run$Warnings <- paste(unique(warning_messages), collapse = " | ")
    if (inherits(fit, "error")) {
      run$FailureReason <- paste0(class(fit)[[1L]], ": ", conditionMessage(fit))
      run_rows[[length(run_rows) + 1L]] <- run
      next
    }

    run$FitReturned <- TRUE
    run$Iterations <- as.numeric(fit$iter)
    run$FinalRelativeDevianceChange <- last_relative_deviance_change(fit$deviance.history)
    run$Deviance <- as.numeric(fit$deviance)
    run$LogLik <- as.numeric(stats::logLik(fit))
    run$PopulationMean <- as.numeric(fit$mu)
    run$PopulationSD <- as.numeric(fit$sigma)
    run$QuadratureTailMass <- sum(as.numeric(fit$pi.k[c(1L, length(fit$pi.k))]))
    run$EAPAvailable <- nrow(as.data.frame(fit$person)) > 0L &&
      all(is.finite(as.numeric(fit$person$EAP)))
    run$EAPReliability <- as.numeric(fit$EAP.rel)
    run$ReportedNParameters <- as.numeric(fit$ic$np)
    run$Converged <- is.finite(run$Deviance) && run$Iterations < maxiter &&
      is.finite(run$FinalRelativeDevianceChange) &&
      run$FinalRelativeDevianceChange <= maxdevchange

    rater_table <- as.data.frame(fit$rater, stringsAsFactors = FALSE)
    rater_table <- rater_table[match(rater_levels, as.character(rater_table$rater)), , drop = FALSE]
    fixed_flag <- !is.na(fixed %||% rep(NA_real_, length(rater_levels)))
    nonfixed <- which(!fixed_flag)
    run$RaterSEReady <- length(nonfixed) == 0L || all(
      is.finite(as.numeric(fit$se.b.rater[nonfixed])) &
        as.numeric(fit$se.b.rater[nonfixed]) > 0
    )
    threshold_se <- as.numeric(fit$se.tau.item)
    run$ItemThresholdSEReady <- length(threshold_se) > 0L &&
      all(is.finite(threshold_se) & threshold_se > 0)

    if (nrow(run_anchors)) {
      estimated_anchor <- as.numeric(fit$b.rater[match(run_anchors$Level, rater_levels)])
      deviations <- abs(estimated_anchor - as.numeric(run_anchors$Anchor))
      run$AnchorMaxAbsDeviation <- max(deviations)
      run$AnchorContractPassed <- all(is.finite(deviations)) && all(deviations <= 1e-10)
    } else {
      run$AnchorMaxAbsDeviation <- 0
      run$AnchorContractPassed <- TRUE
    }
    run$InferenceReady <- run$Converged && run$RaterSEReady &&
      run$ItemThresholdSEReady && run$EAPAvailable &&
      is.finite(run$QuadratureTailMass) &&
      run$QuadratureTailMass <= tail_mass_tolerance
    run$AnalysisEligible <- run$InferenceReady && run$AnchorContractPassed
    if (run$AnalysisEligible) {
      run$FailureStage <- ""
      run$FailureReason <- ""
    } else if (!run$Converged) {
      run$FailureStage <- "convergence"
      run$FailureReason <- "rm.facets did not satisfy the retained loop-exit evidence"
    } else if (!run$RaterSEReady || !run$ItemThresholdSEReady) {
      run$FailureStage <- "uncertainty"
      run$FailureReason <- "one or more nonfixed rater/item-threshold SEs are unavailable"
    } else if (!run$EAPAvailable) {
      run$FailureStage <- "person_scores"
      run$FailureReason <- "finite EAP scores unavailable"
    } else if (run$QuadratureTailMass > tail_mass_tolerance) {
      run$FailureStage <- "quadrature"
      run$FailureReason <- "quadrature endpoint mass exceeds review tolerance"
    } else if (!run$AnchorContractPassed) {
      run$FailureStage <- "anchors"
      run$FailureReason <- "fixed rater estimates differ from supplied anchors"
    }

    rater_truth <- run_truth[run_truth$Facet == "Rater", c("Level", "Truth")]
    rater_output <- merge(
      rater_truth,
      data.frame(
        Level = rater_levels,
        Estimate = as.numeric(fit$b.rater),
        SE = as.numeric(fit$se.b.rater),
        Anchored = fixed_flag,
        stringsAsFactors = FALSE
      ),
      by = "Level",
      all.x = TRUE,
      sort = FALSE
    )
    rater_output$EstimateAligned <- as.numeric(rater_output$Estimate)
    rater_output$TruthAligned <- as.numeric(rater_output$Truth)
    rater_output$ComparisonScale <- "anchor_identified_absolute"
    if (!nrow(run_anchors)) {
      rater_output$EstimateAligned <- rater_output$Estimate - mean(rater_output$Estimate)
      rater_output$TruthAligned <- rater_output$Truth - mean(rater_output$Truth)
      rater_output$ComparisonScale <- "mean_aligned_location"
    }
    rater_output$ErrorAligned <- rater_output$EstimateAligned - rater_output$TruthAligned
    rater_output$IncludedInSummary <- run$AnalysisEligible
    rater_output$CoverageEligible <- run$AnalysisEligible & !rater_output$Anchored &
      is.finite(rater_output$SE) & rater_output$SE > 0
    rater_output$RunId <- run_id
    rater_output$ConditionId <- as.character(manifest_row$ConditionId)
    rater_output$Design <- as.character(manifest_row$Design)
    rater_output$TruthBias <- as.numeric(manifest_row$TruthBias)
    rater_output$Replicate <- as.integer(manifest_row$Replicate)
    rater_output$Seed <- as.numeric(manifest_row$Seed)
    rater_output$Engine <- "sirt"
    rater_output$Estimator <- "MML"
    rater_output$Mode <- as.character(mode_row$Mode)
    rater_rows[[length(rater_rows) + 1L]] <- rater_output[c(
      "RunId", "ConditionId", "Engine", "Estimator", "Mode", "Design",
      "TruthBias", "Replicate", "Seed", "Level", "Truth", "Estimate", "SE",
      "EstimateAligned", "TruthAligned", "ErrorAligned", "ComparisonScale",
      "Anchored", "IncludedInSummary", "CoverageEligible"
    )]

    task_truth <- stats::setNames(
      run_truth$Truth[run_truth$Facet == "Task"],
      run_truth$Level[run_truth$Facet == "Task"]
    )
    criterion_truth <- stats::setNames(
      run_truth$Truth[run_truth$Facet == "Criterion"],
      run_truth$Level[run_truth$Facet == "Criterion"]
    )
    item_output <- wide$item_map
    item_output$Truth <- as.numeric(task_truth[item_output$Task]) +
      as.numeric(criterion_truth[item_output$Criterion])
    item_output$Estimate <- as.numeric(fit$delta.item)
    item_output$EstimateAligned <- item_output$Estimate - mean(item_output$Estimate)
    item_output$TruthAligned <- item_output$Truth - mean(item_output$Truth)
    item_output$ErrorAligned <- item_output$EstimateAligned - item_output$TruthAligned
    item_output$IncludedInSummary <- run$AnalysisEligible
    item_output$RunId <- run_id
    item_output$ConditionId <- as.character(manifest_row$ConditionId)
    item_output$Design <- as.character(manifest_row$Design)
    item_output$TruthBias <- as.numeric(manifest_row$TruthBias)
    item_output$Replicate <- as.integer(manifest_row$Replicate)
    item_output$Seed <- as.numeric(manifest_row$Seed)
    item_output$Engine <- "sirt"
    item_output$Estimator <- "MML"
    item_output$Mode <- as.character(mode_row$Mode)
    item_output$ComparisonScale <- "mean_aligned_virtual_item_location"
    item_output$SE <- NA_real_
    item_output$UncertaintyBoundary <- "delta.item covariance not returned by rm.facets"
    item_rows[[length(item_rows) + 1L]] <- item_output[c(
      "RunId", "ConditionId", "Engine", "Estimator", "Mode", "Design",
      "TruthBias", "Replicate", "Seed", "Item", "Task", "Criterion",
      "Truth", "Estimate", "SE", "EstimateAligned", "TruthAligned",
      "ErrorAligned", "ComparisonScale", "IncludedInSummary",
      "UncertaintyBoundary"
    )]

    persons <- as.data.frame(fit$person, stringsAsFactors = FALSE)
    person_output <- data.frame(
      RunId = run_id,
      ConditionId = as.character(manifest_row$ConditionId),
      Engine = "sirt",
      Estimator = "MML",
      Mode = as.character(mode_row$Mode),
      Design = as.character(manifest_row$Design),
      TruthBias = as.numeric(manifest_row$TruthBias),
      Replicate = as.integer(manifest_row$Replicate),
      Seed = as.numeric(manifest_row$Seed),
      Person = as.character(persons$pid),
      Estimate = as.numeric(persons$EAP),
      SE = as.numeric(persons$SE),
      ParameterStatus = "EAP under estimated sirt population distribution",
      IncludedInSummary = run$AnalysisEligible,
      stringsAsFactors = FALSE
    )
    person_rows[[length(person_rows) + 1L]] <- person_output
    run_rows[[length(run_rows) + 1L]] <- run
  }
}

runs <- do.call(rbind, run_rows)
rater_recovery <- do.call(rbind, rater_rows)
item_recovery <- do.call(rbind, item_rows)
person_estimates <- do.call(rbind, person_rows)
utils::write.csv(runs, file.path(output_dir, "sirt_runs.csv"), row.names = FALSE)
utils::write.csv(rater_recovery, file.path(output_dir, "sirt_rater_recovery.csv"), row.names = FALSE)
utils::write.csv(item_recovery, file.path(output_dir, "sirt_item_recovery.csv"), row.names = FALSE)
utils::write.csv(person_estimates, file.path(output_dir, "sirt_person_estimates.csv"), row.names = FALSE)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1L]]) else NA_character_
identity <- data.frame(
  BundleInventorySHA256 = bundle_inventory_sha256,
  SirtVersion = as.character(utils::packageVersion("sirt")),
  SirtFunctionSHA256 = function_hash,
  SirtFunctionBodySHA256 = function_body_hash,
  AdapterScriptSHA256 = if (!is.na(script_path)) {
    digest::digest(file = normalizePath(script_path), algo = "sha256", serialize = FALSE)
  } else NA_character_,
  Modes = paste(modes$Mode, collapse = "|"),
  Model = "Task x Criterion virtual-item PCM plus rater severity",
  PersonDistribution = "estimated normal distribution on fixed -8..8 grid",
  AnchorScope = "Rater fixed values supported; fixed-parameter SE not interpreted",
  BiasScope = "local Rater x Task interaction unsupported",
  InformationCriteriaScope = "withheld due unanchored center-mode parameter-count convention",
  stringsAsFactors = FALSE
)
utils::write.csv(identity, file.path(output_dir, "sirt_adapter_identity.csv"), row.names = FALSE)
