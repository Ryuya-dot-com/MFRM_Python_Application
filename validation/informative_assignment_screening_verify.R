#!/usr/bin/env Rscript

# Independent base-R verification of the corrected informative-assignment
# screening aggregate.  Scientific contrasts are rebuilt from native-precision
# Python truth-error rows.  FACETS display-rounded values are not used to
# reconstruct recovery statistics.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("usage: Rscript informative_assignment_screening_verify.R AGGREGATE_DIR OUTPUT")
}

aggregate_dir <- normalizePath(args[[1L]], mustWork = TRUE)
output_path <- args[[2L]]

read_artifact <- function(filename) {
  read.csv(
    file.path(aggregate_dir, filename),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

as_flag <- function(values) {
  tolower(as.character(values)) == "true"
}

modes <- c(
  "PYTHON_JMLE",
  "PYTHON_MML_FIXED_SD08_Q31",
  "PYTHON_MML_FREE_SD_Q31",
  "PYTHON_EXACT_CMLE"
)

recovery <- read_artifact("recovery.csv")
thresholds <- read_artifact("thresholds.csv")
runs <- read_artifact("run_ledger.csv")
python_summary <- read_artifact("screening_contrast_summary.csv")
python_rater <- read_artifact("rater_level_summary.csv")
python_free_sd <- read_artifact("free_sd_paired.csv")

eligible_recovery <- recovery[
  recovery$EstimatorMode %in% modes &
    as_flag(recovery$IncludedInStudy) &
    recovery$Facet %in% c("Rater", "Task", "Criterion"),
]
eligible_recovery$ErrorAligned <- as.numeric(eligible_recovery$ErrorAligned)

facet_keys <- unique(eligible_recovery[
  , c("RunId", "Replicate", "Design", "EstimatorMode", "Facet")
])
facet_rows <- lapply(seq_len(nrow(facet_keys)), function(index) {
  key <- facet_keys[index,]
  selected <- eligible_recovery[
    eligible_recovery$RunId == key$RunId &
      eligible_recovery$EstimatorMode == key$EstimatorMode &
      eligible_recovery$Facet == key$Facet,
  ]
  errors <- selected$ErrorAligned
  data.frame(
    Replicate = as.integer(key$Replicate),
    Design = key$Design,
    EstimatorMode = key$EstimatorMode,
    RecoveryDomain = paste0("Facet:", key$Facet),
    RMSE = sqrt(mean(errors^2)),
    MAE = mean(abs(errors)),
    stringsAsFactors = FALSE
  )
})
facet_loss <- do.call(rbind, facet_rows)

eligible_thresholds <- thresholds[
  thresholds$EstimatorMode %in% modes & as_flag(thresholds$IncludedInStudy),
]
eligible_thresholds$TruthError <- as.numeric(eligible_thresholds$TruthError)
threshold_keys <- unique(eligible_thresholds[
  , c("RunId", "Replicate", "Design", "EstimatorMode")
])
threshold_rows <- lapply(seq_len(nrow(threshold_keys)), function(index) {
  key <- threshold_keys[index,]
  selected <- eligible_thresholds[
    eligible_thresholds$RunId == key$RunId &
      eligible_thresholds$EstimatorMode == key$EstimatorMode,
  ]
  errors <- selected$TruthError
  if (length(errors) != 6L || any(!is.finite(errors))) {
    stop(sprintf("invalid threshold group: %s / %s", key$RunId, key$EstimatorMode))
  }
  data.frame(
    Replicate = as.integer(key$Replicate),
    Design = key$Design,
    EstimatorMode = key$EstimatorMode,
    RecoveryDomain = "Threshold",
    RMSE = sqrt(mean(errors^2)),
    MAE = mean(abs(errors)),
    stringsAsFactors = FALSE
  )
})
threshold_loss <- do.call(rbind, threshold_rows)

wide_loss <- rbind(facet_loss, threshold_loss)
long_loss <- rbind(
  data.frame(
    wide_loss[, c("Replicate", "Design", "EstimatorMode", "RecoveryDomain")],
    Metric = "RMSE",
    Loss = wide_loss$RMSE,
    stringsAsFactors = FALSE
  ),
  data.frame(
    wide_loss[, c("Replicate", "Design", "EstimatorMode", "RecoveryDomain")],
    Metric = "MAE",
    Loss = wide_loss$MAE,
    stringsAsFactors = FALSE
  )
)

aligned <- long_loss[
  long_loss$Design == "ability_severity_aligned_connected",
  c("Replicate", "EstimatorMode", "RecoveryDomain", "Metric", "Loss")
]
planned <- long_loss[
  long_loss$Design == "planned_connected",
  c("Replicate", "EstimatorMode", "RecoveryDomain", "Metric", "Loss")
]
names(aligned)[[5L]] <- "LossAligned"
names(planned)[[5L]] <- "LossPlanned"
contrasts <- merge(
  aligned,
  planned,
  by = c("Replicate", "EstimatorMode", "RecoveryDomain", "Metric"),
  all = FALSE,
  sort = TRUE
)
contrasts$Contrast <- contrasts$LossAligned - contrasts$LossPlanned

summary_keys <- unique(contrasts[
  , c("EstimatorMode", "RecoveryDomain", "Metric")
])
r_summary_rows <- lapply(seq_len(nrow(summary_keys)), function(index) {
  key <- summary_keys[index,]
  values <- contrasts[
    contrasts$EstimatorMode == key$EstimatorMode &
      contrasts$RecoveryDomain == key$RecoveryDomain &
      contrasts$Metric == key$Metric,
    "Contrast"
  ]
  values <- values[is.finite(values)]
  n <- length(values)
  mean_value <- mean(values)
  sd_value <- stats::sd(values)
  se_value <- sd_value / sqrt(n)
  half_width <- stats::qt(0.975, df = n - 1L) * se_value
  data.frame(
    EstimatorMode = key$EstimatorMode,
    RecoveryDomain = key$RecoveryDomain,
    Metric = key$Metric,
    FinitePairs = n,
    MeanContrastAlignedMinusPlanned = mean_value,
    SDContrast = sd_value,
    SEContrast = se_value,
    ScreeningCI95Low = mean_value - half_width,
    ScreeningCI95High = mean_value + half_width,
    ScreeningCI95HalfWidth = half_width,
    IntervalExcludesZero = mean_value - half_width > 0 || mean_value + half_width < 0,
    stringsAsFactors = FALSE
  )
})
r_summary <- do.call(rbind, r_summary_rows)

summary_join <- merge(
  r_summary,
  python_summary,
  by = c("EstimatorMode", "RecoveryDomain", "Metric"),
  suffixes = c("R", "Python"),
  all = TRUE,
  sort = TRUE
)
summary_fields <- c(
  "FinitePairs",
  "MeanContrastAlignedMinusPlanned",
  "SDContrast",
  "SEContrast",
  "ScreeningCI95Low",
  "ScreeningCI95High",
  "ScreeningCI95HalfWidth"
)
summary_differences <- unlist(lapply(summary_fields, function(field) {
  abs(
    as.numeric(summary_join[[paste0(field, "R")]]) -
      as.numeric(summary_join[[paste0(field, "Python")]])
  )
}))
summary_logic_pass <- all(
  summary_join$IntervalExcludesZeroR == as_flag(summary_join$IntervalExcludesZeroPython)
)

rater <- eligible_recovery[eligible_recovery$Facet == "Rater",]
rater_aligned <- rater[
  rater$Design == "ability_severity_aligned_connected",
  c("Replicate", "EstimatorMode", "Level", "ErrorAligned")
]
rater_planned <- rater[
  rater$Design == "planned_connected",
  c("Replicate", "EstimatorMode", "Level", "ErrorAligned")
]
names(rater_aligned)[[4L]] <- "ErrorAlignedDesign"
names(rater_planned)[[4L]] <- "ErrorPlannedDesign"
rater_contrasts <- merge(
  rater_aligned,
  rater_planned,
  by = c("Replicate", "EstimatorMode", "Level"),
  all = FALSE,
  sort = TRUE
)
rater_contrasts$Contrast <- (
  rater_contrasts$ErrorAlignedDesign - rater_contrasts$ErrorPlannedDesign
)
rater_keys <- unique(rater_contrasts[, c("EstimatorMode", "Level")])
r_rater_rows <- lapply(seq_len(nrow(rater_keys)), function(index) {
  key <- rater_keys[index,]
  values <- rater_contrasts[
    rater_contrasts$EstimatorMode == key$EstimatorMode &
      rater_contrasts$Level == key$Level,
    "Contrast"
  ]
  data.frame(
    EstimatorMode = key$EstimatorMode,
    Level = key$Level,
    N = length(values),
    Mean = mean(values),
    SD = stats::sd(values),
    stringsAsFactors = FALSE
  )
})
r_rater <- do.call(rbind, r_rater_rows)
rater_join <- merge(
  r_rater,
  python_rater,
  by = c("EstimatorMode", "Level"),
  suffixes = c("R", "Python"),
  all = TRUE,
  sort = TRUE
)
rater_differences <- c(
  abs(as.numeric(rater_join$NR) - as.numeric(rater_join$NPython)),
  abs(as.numeric(rater_join$MeanR) - as.numeric(rater_join$MeanPython)),
  abs(as.numeric(rater_join$SDR) - as.numeric(rater_join$SDPython))
)

free_runs <- runs[
  runs$EstimatorMode == "PYTHON_MML_FREE_SD_Q31",
  c("Replicate", "Design", "EstimatedPopulationSD")
]
free_aligned <- free_runs[
  free_runs$Design == "ability_severity_aligned_connected",
  c("Replicate", "EstimatedPopulationSD")
]
free_complete <- free_runs[
  free_runs$Design == "complete",
  c("Replicate", "EstimatedPopulationSD")
]
free_planned <- free_runs[
  free_runs$Design == "planned_connected",
  c("Replicate", "EstimatedPopulationSD")
]
names(free_aligned)[[2L]] <- "ability_severity_aligned_connected"
names(free_complete)[[2L]] <- "complete"
names(free_planned)[[2L]] <- "planned_connected"
r_free_sd <- Reduce(
  function(left, right) merge(left, right, by = "Replicate", all = FALSE, sort = TRUE),
  list(free_aligned, free_complete, free_planned)
)
r_free_sd$ContrastAlignedMinusPlanned <- (
  r_free_sd$ability_severity_aligned_connected - r_free_sd$planned_connected
)
free_join <- merge(
  r_free_sd,
  python_free_sd,
  by = "Replicate",
  suffixes = c("R", "Python"),
  all = TRUE,
  sort = TRUE
)
free_fields <- c(
  "ability_severity_aligned_connected",
  "complete",
  "planned_connected",
  "ContrastAlignedMinusPlanned"
)
free_differences <- unlist(lapply(free_fields, function(field) {
  abs(
    as.numeric(free_join[[paste0(field, "R")]]) -
      as.numeric(free_join[[paste0(field, "Python")]])
  )
}))

max_or_inf <- function(values) {
  if (length(values) == 0L || any(!is.finite(values))) Inf else max(values)
}
summary_max <- max_or_inf(summary_differences)
rater_max <- max_or_inf(rater_differences)
free_max <- max_or_inf(free_differences)
overall_max <- max(summary_max, rater_max, free_max)
structural_pass <- (
  nrow(r_summary) == 32L &&
    nrow(contrasts) == 640L &&
    all(r_summary$FinitePairs == 20L) &&
    nrow(r_rater) == 16L &&
    nrow(rater_contrasts) == 320L &&
    nrow(r_free_sd) == 20L &&
    nrow(summary_join) == 32L &&
    nrow(rater_join) == 16L &&
    nrow(free_join) == 20L
)
numeric_pass <- is.finite(overall_max) && overall_max <= 1e-12
verification_pass <- structural_pass && numeric_pass && summary_logic_pass

free_contrast <- r_free_sd$ContrastAlignedMinusPlanned
verified <- data.frame(
  ScreeningCellsR = nrow(r_summary),
  PairedMetricRowsR = nrow(contrasts),
  RaterLevelCellsR = nrow(r_rater),
  RaterLevelPairsR = nrow(rater_contrasts),
  FreeSDPairsR = nrow(r_free_sd),
  MaximumScreeningSummaryDifference = summary_max,
  MaximumRaterLevelDifference = rater_max,
  MaximumFreeSDPairDifference = free_max,
  MaximumAbsolutePythonDifference = overall_max,
  ScreeningLogicalAgreementPass = summary_logic_pass,
  StructuralAgreementPass = structural_pass,
  PythonNumericAgreementPass = numeric_pass,
  FreeSDAlignedMeanR = mean(r_free_sd$ability_severity_aligned_connected),
  FreeSDPlannedMeanR = mean(r_free_sd$planned_connected),
  FreeSDMeanContrastR = mean(free_contrast),
  FreeSDContrastCI95LowR = mean(free_contrast) -
    stats::qt(0.975, df = length(free_contrast) - 1L) *
      stats::sd(free_contrast) / sqrt(length(free_contrast)),
  FreeSDContrastCI95HighR = mean(free_contrast) +
    stats::qt(0.975, df = length(free_contrast) - 1L) *
      stats::sd(free_contrast) / sqrt(length(free_contrast)),
  VerificationPass = verification_pass,
  stringsAsFactors = FALSE
)
write.csv(verified, output_path, row.names = FALSE, na = "")

cat(sprintf("screening_cells=%d\n", nrow(r_summary)))
cat(sprintf("paired_metric_rows=%d\n", nrow(contrasts)))
cat(sprintf("rater_level_cells=%d\n", nrow(r_rater)))
cat(sprintf("free_sd_pairs=%d\n", nrow(r_free_sd)))
cat(sprintf("max_abs_python_difference=%.17g\n", overall_max))
cat(sprintf("free_sd_mean_contrast=%.17g\n", mean(free_contrast)))
cat(sprintf("verification_pass=%s\n", tolower(as.character(verification_pass))))
if (!verification_pass) {
  stop("R/Python informative-assignment screening verification failed")
}
