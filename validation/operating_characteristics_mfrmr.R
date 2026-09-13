#!/usr/bin/env Rscript

# Fit mfrmr JMLE to the byte-validated operating-characteristics bundle.
# The matched-control and strict numerical modes remain separate estimands.

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

as_flag <- function(value) isTRUE(as.logical(value[[1L]]))
as_number <- function(value, default = NA_real_) {
  if (is.null(value) || length(value) == 0L) return(default)
  out <- suppressWarnings(as.numeric(value[[1L]]))
  if (is.finite(out)) out else default
}
as_text <- function(value, default = "") {
  if (is.null(value) || length(value) == 0L || is.na(value[[1L]])) return(default)
  as.character(value[[1L]])
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("input", "output", "mfrmr-lib", "mfrmr-git-head", "mfrmr-source-state")
missing <- required[!vapply(required, function(name) nzchar(args[[name]] %||% ""), logical(1))]
if (length(missing)) {
  stop("Missing required arguments: ", paste(missing, collapse = ", "), call. = FALSE)
}

.libPaths(c(normalizePath(args[["mfrmr-lib"]], mustWork = TRUE), .libPaths()))
for (package in c("mfrmr", "digest")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Required adapter dependency unavailable: ", package, call. = FALSE)
  }
}
if (!identical(as.character(utils::packageVersion("mfrmr")), "0.2.3")) {
  stop("mfrmr adapter requires the isolated 0.2.3 development snapshot.", call. = FALSE)
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
if (!all(as.logical(bridge_validation$Passed))) {
  stop("R bridge validation contains a failed check.", call. = FALSE)
}
if (!all(as.character(bridge_validation$BundleInventorySHA256) == bundle_inventory_sha256)) {
  stop("R bridge validation is stale for the current generated bundle.", call. = FALSE)
}

hash_ok <- vapply(seq_len(nrow(inventory)), function(i) {
  path <- file.path(input_dir, as.character(inventory$File[[i]]))
  file.exists(path) && identical(
    digest::digest(file = path, algo = "sha256", serialize = FALSE),
    as.character(inventory$SHA256[[i]])
  )
}, logical(1))
if (!all(hash_ok)) stop("Generated bundle byte hash changed after bridge validation.", call. = FALSE)

fit_body_hash <- digest::digest(
  list(
    formals = formals(mfrmr::fit_mfrm),
    body = body(mfrmr::fit_mfrm)
  ),
  algo = "sha256",
  serialize = TRUE
)
fit_function_hash <- digest::digest(
  c(fit_mfrm = fit_body_hash),
  algo = "sha256",
  serialize = TRUE
)
registered_mfrmr <- bridge_engines[bridge_engines$Engine == "mfrmr", , drop = FALSE]
if (nrow(registered_mfrmr) != 1L ||
    !identical(as.character(registered_mfrmr$Version[[1L]]), "0.2.3") ||
    !identical(as.character(registered_mfrmr$FunctionSHA256[[1L]]), fit_function_hash)) {
  stop("Loaded mfrmr code identity differs from the validated R bridge.", call. = FALSE)
}

modes <- data.frame(
  Mode = c("MFRMR_JML_MATCHED_CONTROL", "MFRMR_JML_STRICT"),
  MaxIt = c(160L, 500L),
  RelTol = c(1e-6, 1e-9),
  NumericalRole = c(
    "Python-control sensitivity; not pooled with strict mode",
    "mfrmr inference-readiness mode"
  ),
  stringsAsFactors = FALSE
)

run_rows <- list()
parameter_rows <- list()
bias_rows <- list()
run_index <- 0L
parameter_index <- 0L
bias_index <- 0L

append_run <- function(row) {
  run_index <<- run_index + 1L
  run_rows[[run_index]] <<- row
}
append_parameters <- function(row) {
  if (is.null(row) || !nrow(row)) return(invisible(NULL))
  parameter_index <<- parameter_index + 1L
  parameter_rows[[parameter_index]] <<- row
}
append_bias <- function(row) {
  bias_index <<- bias_index + 1L
  bias_rows[[bias_index]] <<- row
}

base_run <- function(manifest_row, mode_row, rows) {
  data.frame(
    SchemaVersion = as.character(manifest_row$SchemaVersion),
    RunId = as.character(manifest_row$RunId),
    ConditionId = as.character(manifest_row$ConditionId),
    Design = as.character(manifest_row$Design),
    TruthBias = as.numeric(manifest_row$TruthBias),
    TruthPositive = as.logical(manifest_row$TruthPositive),
    Replicate = as.integer(manifest_row$Replicate),
    Seed = as.numeric(manifest_row$Seed),
    Engine = "mfrmr",
    Estimator = "JMLE",
    Mode = as.character(mode_row$Mode),
    NumericalRole = as.character(mode_row$NumericalRole),
    MaxIt = as.integer(mode_row$MaxIt),
    RelTol = as.numeric(mode_row$RelTol),
    FitReturned = FALSE,
    Converged = FALSE,
    InferenceReady = FALSE,
    BiasAvailable = FALSE,
    BiasOptimizationReady = FALSE,
    FocalCellSparse = TRUE,
    AnchorContractPassed = NA,
    AnchorMaxAbsDeviation = NA_real_,
    AnalysisEligible = FALSE,
    FailureStage = "fit",
    FailureReason = "",
    Warnings = "",
    Rows = as.integer(rows),
    Iterations = NA_real_,
    GradientNorm = NA_real_,
    GradientReviewTolerance = NA_real_,
    ConvergenceCode = NA_real_,
    ConvergenceStatus = "",
    ReadinessReasonCodes = "",
    LogLik = NA_real_,
    ElapsedSeconds = NA_real_,
    BundleInventorySHA256 = bundle_inventory_sha256,
    MfrmrVersion = as.character(utils::packageVersion("mfrmr")),
    MfrmrFunctionSHA256 = fit_function_hash,
    MfrmrFitBodySHA256 = fit_body_hash,
    SourceGitHead = as.character(args[["mfrmr-git-head"]]),
    SourceState = as.character(args[["mfrmr-source-state"]]),
    stringsAsFactors = FALSE
  )
}

base_bias <- function(manifest_row, mode_row) {
  data.frame(
    RunId = as.character(manifest_row$RunId),
    ConditionId = as.character(manifest_row$ConditionId),
    Design = as.character(manifest_row$Design),
    TruthBias = as.numeric(manifest_row$TruthBias),
    TruthPositive = as.logical(manifest_row$TruthPositive),
    Replicate = as.integer(manifest_row$Replicate),
    Seed = as.numeric(manifest_row$Seed),
    Engine = "mfrmr",
    Estimator = "JMLE",
    Mode = as.character(mode_row$Mode),
    AnalysisEligible = FALSE,
    BiasEstimate = NA_real_,
    BiasSE = NA_real_,
    p = NA_real_,
    p_holm = NA_real_,
    p_bh = NA_real_,
    t = NA_real_,
    ObsN = NA_real_,
    AbsBias = NA_real_,
    SparseCell = TRUE,
    BiasOptimizationReady = FALSE,
    DecisionHolmRaw = NA,
    DecisionHolmDisplayed = NA,
    DecisionPracticalRaw = NA,
    DecisionPracticalDisplayed = NA,
    DecisionStrongRaw = NA,
    DecisionStrongDisplayed = NA,
    DecisionAnyNonSparseFlag = NA,
    InferenceTier = "screening",
    InterpretationBoundary = paste(
      "Conditional plug-in bias screen; raw probabilities and magnitudes govern decisions;",
      "failed, inference-unready, or sparse fits remain ineligible."
    ),
    stringsAsFactors = FALSE
  )
}

build_parameter_rows <- function(manifest_row, mode_row, diagnostics, ready, run_anchors) {
  measures <- as.data.frame(diagnostics$measures, stringsAsFactors = FALSE)
  measures <- measures[measures$Facet %in% c("Rater", "Task", "Criterion"), , drop = FALSE]
  if (!nrow(measures)) return(data.frame())
  measure_columns <- c("Facet", "Level", "Estimate", "SE")
  if (!all(measure_columns %in% names(measures))) return(data.frame())
  run_truth <- truth[truth$RunId == manifest_row$RunId, c("Facet", "Level", "Truth")]
  merged <- merge(
    run_truth,
    measures[measure_columns],
    by = c("Facet", "Level"),
    all.x = TRUE,
    sort = FALSE
  )
  merged$RawError <- as.numeric(merged$Estimate) - as.numeric(merged$Truth)
  merged$EstimateAligned <- NA_real_
  merged$TruthAligned <- as.numeric(merged$Truth)
  merged$ErrorAligned <- NA_real_
  merged$ComparisonScale <- "mean_aligned_location"
  anchored_facets <- unique(as.character(run_anchors$Facet))
  for (facet in unique(as.character(merged$Facet))) {
    index <- which(as.character(merged$Facet) == facet)
    if (facet %in% anchored_facets) {
      merged$EstimateAligned[index] <- as.numeric(merged$Estimate[index])
      merged$ErrorAligned[index] <- merged$RawError[index]
      merged$ComparisonScale[index] <- "anchor_identified_absolute"
    } else {
      shift <- mean(merged$RawError[index], na.rm = TRUE)
      merged$EstimateAligned[index] <- as.numeric(merged$Estimate[index]) - shift
      merged$ErrorAligned[index] <- merged$EstimateAligned[index] - merged$TruthAligned[index]
    }
  }
  anchor_key <- paste(run_anchors$Facet, run_anchors$Level, sep = "\r")
  merged_key <- paste(merged$Facet, merged$Level, sep = "\r")
  merged$Anchored <- merged_key %in% anchor_key
  merged$IncludedInSummary <- isTRUE(ready) & is.finite(merged$ErrorAligned)
  merged$RunId <- as.character(manifest_row$RunId)
  merged$ConditionId <- as.character(manifest_row$ConditionId)
  merged$Engine <- "mfrmr"
  merged$Estimator <- "JMLE"
  merged$Mode <- as.character(mode_row$Mode)
  merged$Design <- as.character(manifest_row$Design)
  merged$TruthBias <- as.numeric(manifest_row$TruthBias)
  merged$Replicate <- as.integer(manifest_row$Replicate)
  merged$Seed <- as.numeric(manifest_row$Seed)
  merged$ParameterType <- merged$Facet
  merged$CoverageMethod <- ifelse(
    merged$ComparisonScale == "anchor_identified_absolute",
    "anchor-identified absolute conditional-Wald diagnostic",
    "mean-aligned conditional-Wald diagnostic"
  )
  merged[c(
    "RunId", "ConditionId", "Engine", "Estimator", "Mode", "Design",
    "TruthBias", "Replicate", "Seed", "Facet", "Level", "Truth", "Estimate",
    "SE", "RawError", "EstimateAligned", "TruthAligned", "ErrorAligned",
    "ComparisonScale", "Anchored", "IncludedInSummary", "ParameterType",
    "CoverageMethod"
  )]
}

for (manifest_index in seq_len(nrow(manifest))) {
  manifest_row <- manifest[manifest_index, , drop = FALSE]
  data <- ratings[
    ratings$RunId == manifest_row$RunId,
    c("Person", "Rater", "Task", "Criterion", "Score"),
    drop = FALSE
  ]
  run_anchors <- anchors[
    anchors$RunId == manifest_row$RunId,
    c("Facet", "Level", "Anchor"),
    drop = FALSE
  ]
  anchor_arg <- if (nrow(run_anchors)) run_anchors else NULL

  for (mode_index in seq_len(nrow(modes))) {
    mode_row <- modes[mode_index, , drop = FALSE]
    run <- base_run(manifest_row, mode_row, nrow(data))
    bias_row <- base_bias(manifest_row, mode_row)
    warning_messages <- character()
    started <- proc.time()[["elapsed"]]
    fit <- tryCatch(
      withCallingHandlers(
        mfrmr::fit_mfrm(
          data = data,
          person = "Person",
          facets = c("Rater", "Task", "Criterion"),
          score = "Score",
          rating_min = 0,
          rating_max = as.integer(manifest_row$Categories) - 1L,
          keep_original = TRUE,
          model = "RSM",
          method = "JML",
          anchors = anchor_arg,
          noncenter_facet = "Person",
          anchor_policy = "warn",
          min_common_anchors = 2L,
          min_obs_per_element = 1L,
          min_obs_per_category = 1L,
          maxit = as.integer(mode_row$MaxIt),
          reltol = as.numeric(mode_row$RelTol),
          optimizer = "BFGS"
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
      run$FailureReason <- paste0(class(fit)[[1L]], ": ", conditionMessage(fit))
      append_run(run)
      append_bias(bias_row)
      next
    }

    run$FitReturned <- TRUE
    summary <- as.data.frame(fit$summary, stringsAsFactors = FALSE)[1L, , drop = FALSE]
    run$Converged <- as_flag(summary$Converged)
    run$InferenceReady <- as_flag(summary$InferenceReady)
    run$Iterations <- as_number(summary$Iterations)
    run$GradientNorm <- as_number(summary$TerminalGradientSupNorm)
    run$GradientReviewTolerance <- as_number(summary$GradientReviewTolerance)
    run$ConvergenceCode <- as_number(summary$ConvergenceCode)
    run$ConvergenceStatus <- as_text(summary$ConvergenceStatus)
    run$ReadinessReasonCodes <- as_text(summary$ReadinessReasonCodes)
    run$LogLik <- as_number(summary$LogLik)

    facet_estimates <- as.data.frame(fit$facets$others, stringsAsFactors = FALSE)
    if (nrow(run_anchors)) {
      anchor_check <- merge(
        run_anchors,
        facet_estimates[c("Facet", "Level", "Estimate")],
        by = c("Facet", "Level"),
        all.x = TRUE,
        sort = FALSE
      )
      deviation <- abs(as.numeric(anchor_check$Estimate) - as.numeric(anchor_check$Anchor))
      run$AnchorMaxAbsDeviation <- if (length(deviation)) max(deviation, na.rm = TRUE) else NA_real_
      run$AnchorContractPassed <- length(deviation) == nrow(run_anchors) &&
        all(is.finite(deviation)) && all(deviation <= 1e-10)
    } else {
      run$AnchorContractPassed <- TRUE
      run$AnchorMaxAbsDeviation <- 0
    }

    diagnostics <- tryCatch(
      mfrmr::diagnose_mfrm(
        fit,
        diagnostic_mode = "legacy",
        residual_pca = "none"
      ),
      error = function(error) error
    )
    if (inherits(diagnostics, "error")) {
      run$FailureStage <- "diagnostics"
      run$FailureReason <- conditionMessage(diagnostics)
      append_run(run)
      append_bias(bias_row)
      next
    }

    bias <- tryCatch(
      mfrmr::estimate_bias(
        fit,
        diagnostics,
        facet_a = "Rater",
        facet_b = "Task",
        omit_extreme = FALSE
      ),
      error = function(error) error
    )
    if (!inherits(bias, "error")) {
      table <- as.data.frame(bias$table, stringsAsFactors = FALSE)
      if (nrow(table)) {
        table$p_holm <- stats::p.adjust(as.numeric(table[["Prob."]]), method = "holm")
        table$p_bh <- stats::p.adjust(as.numeric(table[["Prob."]]), method = "BH")
        focal <- table[
          as.character(table$FacetA_Level) == "R01" &
            as.character(table$FacetB_Level) == "T01",
          ,
          drop = FALSE
        ]
        if (nrow(focal) == 1L) {
          value <- as.numeric(focal[["Bias Size"]])
          se <- as.numeric(focal[["S.E."]])
          probability <- as.numeric(focal[["Prob."]])
          p_holm <- as.numeric(focal$p_holm)
          p_bh <- as.numeric(focal$p_bh)
          statistic <- as.numeric(focal$t)
          obs_n <- as.numeric(focal$ObsN)
          sparse <- !is.finite(obs_n) || obs_n < 5
          optimization_ready <- identical(as.character(focal$OptimizationStatus), "ok")
          available <- all(is.finite(c(value, se, probability, p_holm, p_bh, statistic)))
          holm_raw <- is.finite(p_holm) && p_holm < 0.05
          practical_raw <- is.finite(value) && abs(value) >= 0.50
          bh_t_raw <- is.finite(p_bh) && p_bh < 0.05 && is.finite(statistic) && abs(statistic) >= 2
          bias_row$BiasEstimate <- value
          bias_row$BiasSE <- se
          bias_row$p <- probability
          bias_row$p_holm <- p_holm
          bias_row$p_bh <- p_bh
          bias_row$t <- statistic
          bias_row$ObsN <- obs_n
          bias_row$AbsBias <- abs(value)
          bias_row$SparseCell <- sparse
          bias_row$BiasOptimizationReady <- optimization_ready
          bias_row$DecisionHolmRaw <- holm_raw
          bias_row$DecisionHolmDisplayed <- round(p_holm, 4) < 0.05
          bias_row$DecisionPracticalRaw <- practical_raw
          bias_row$DecisionPracticalDisplayed <- round(abs(value), 4) >= 0.50
          bias_row$DecisionStrongRaw <- holm_raw && practical_raw && !sparse
          bias_row$DecisionStrongDisplayed <-
            bias_row$DecisionHolmDisplayed && bias_row$DecisionPracticalDisplayed && !sparse
          bias_row$DecisionAnyNonSparseFlag <-
            !sparse && (holm_raw || bh_t_raw || practical_raw)
          run$BiasAvailable <- available
          run$BiasOptimizationReady <- optimization_ready
          run$FocalCellSparse <- sparse
        } else {
          run$FailureStage <- "bias"
          run$FailureReason <- "focal R01 x T01 cell missing or duplicated"
        }
      }
    } else {
      run$FailureStage <- "bias"
      run$FailureReason <- conditionMessage(bias)
    }

    eligible <- isTRUE(run$FitReturned) && isTRUE(run$Converged) &&
      isTRUE(run$InferenceReady) && isTRUE(run$AnchorContractPassed) &&
      isTRUE(run$BiasAvailable) && isTRUE(run$BiasOptimizationReady) &&
      !isTRUE(run$FocalCellSparse)
    run$AnalysisEligible <- eligible
    bias_row$AnalysisEligible <- eligible
    if (eligible) {
      run$FailureStage <- ""
      run$FailureReason <- ""
    } else if (!nzchar(run$FailureReason)) {
      if (!run$Converged) {
        run$FailureStage <- "convergence"
        run$FailureReason <- "fit returned without convergence"
      } else if (!run$InferenceReady) {
        run$FailureStage <- "readiness"
        run$FailureReason <- paste0(
          "inference readiness withheld: ",
          run$ReadinessReasonCodes %||% run$ConvergenceStatus
        )
      } else if (!run$AnchorContractPassed) {
        run$FailureStage <- "anchors"
        run$FailureReason <- "fixed-anchor estimates do not equal supplied targets"
      } else if (!run$BiasOptimizationReady) {
        run$FailureStage <- "bias"
        run$FailureReason <- "focal bias optimizer not ready"
      } else if (run$FocalCellSparse) {
        run$FailureStage <- "bias"
        run$FailureReason <- "focal bias cell below min_n"
      } else {
        run$FailureStage <- "bias"
        run$FailureReason <- "analysis eligibility not satisfied"
      }
    }

    append_parameters(build_parameter_rows(
      manifest_row,
      mode_row,
      diagnostics,
      ready = run$InferenceReady && run$AnchorContractPassed,
      run_anchors = run_anchors
    ))
    append_run(run)
    append_bias(bias_row)
  }
}

runs <- do.call(rbind, run_rows)
parameters <- if (length(parameter_rows)) do.call(rbind, parameter_rows) else data.frame()
bias_decisions <- do.call(rbind, bias_rows)
utils::write.csv(runs, file.path(output_dir, "mfrmr_runs.csv"), row.names = FALSE)
utils::write.csv(parameters, file.path(output_dir, "mfrmr_parameter_recovery.csv"), row.names = FALSE)
utils::write.csv(bias_decisions, file.path(output_dir, "mfrmr_bias_decisions.csv"), row.names = FALSE)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1L]]) else NA_character_
adapter_identity <- data.frame(
  BundleInventorySHA256 = bundle_inventory_sha256,
  MfrmrVersion = as.character(utils::packageVersion("mfrmr")),
  MfrmrFunctionSHA256 = fit_function_hash,
  MfrmrFitBodySHA256 = fit_body_hash,
  SourceGitHead = as.character(args[["mfrmr-git-head"]]),
  SourceState = as.character(args[["mfrmr-source-state"]]),
  AdapterScriptSHA256 = if (!is.na(script_path)) {
    digest::digest(file = normalizePath(script_path), algo = "sha256", serialize = FALSE)
  } else NA_character_,
  Runs = nrow(runs),
  stringsAsFactors = FALSE
)
utils::write.csv(
  adapter_identity,
  file.path(output_dir, "mfrmr_adapter_identity.csv"),
  row.names = FALSE
)

message(
  "Completed ", nrow(runs), " mfrmr JMLE fit attempts across ",
  nrow(manifest), " bundle RunIds and ", nrow(modes), " numerical modes."
)
