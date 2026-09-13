#!/usr/bin/env Rscript

# Independent base-R reconstruction of the five registered informative-
# assignment confirmation endpoints.  Only native-precision Python recovery
# rows are used; FACETS display-rounded values never enter these statistics.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("usage: Rscript informative_assignment_confirmatory_verify.R AGGREGATE_DIR OUTPUT")
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
as_flag <- function(values) tolower(as.character(values)) == "true"

free_mode <- "PYTHON_MML_FREE_SD_Q31"
fixed_mode <- "PYTHON_MML_FIXED_SD08_Q31"
primary_id <- "IA1_FREE_MML_RATER_RMSE_ALIGNED_GT_PLANNED"
secondary_ids <- c(
  "IA2_FIXED_MML_RATER_RMSE_ALIGNED_GT_PLANNED",
  "IA3_FREE_MML_POPULATION_SD_ALIGNED_LT_PLANNED",
  "IA4_FREE_MML_RATER_COMPRESSION_SLOPE_LT_ZERO",
  "IA5_FIXED_MML_RATER_COMPRESSION_SLOPE_LT_ZERO"
)
severity <- c(R01 = -0.45, R02 = -0.15, R03 = 0.15, R04 = 0.45)

recovery <- read_artifact("recovery.csv")
runs <- read_artifact("run_ledger.csv")
python_contrasts <- read_artifact("endpoint_contrasts.csv")
python_primary <- read_artifact("primary_result.csv")
python_secondary <- read_artifact("secondary_results.csv")

rater <- recovery[
  recovery$EstimatorMode %in% c(free_mode, fixed_mode) &
    recovery$Facet == "Rater" &
    as_flag(recovery$IncludedInStudy),
]
rater$ErrorAligned <- as.numeric(rater$ErrorAligned)

run_keys <- unique(rater[, c("Replicate", "Design", "EstimatorMode")])
rmse_rows <- lapply(seq_len(nrow(run_keys)), function(index) {
  key <- run_keys[index,]
  selected <- rater[
    rater$Replicate == key$Replicate &
      rater$Design == key$Design &
      rater$EstimatorMode == key$EstimatorMode,
  ]
  if (nrow(selected) != 4L || any(!is.finite(selected$ErrorAligned))) {
    stop("invalid Rater recovery group")
  }
  data.frame(
    Replicate = as.integer(key$Replicate),
    Design = key$Design,
    EstimatorMode = key$EstimatorMode,
    RaterRMSE = sqrt(mean(selected$ErrorAligned^2)),
    stringsAsFactors = FALSE
  )
})
rmse <- do.call(rbind, rmse_rows)
aligned_rmse <- rmse[
  rmse$Design == "ability_severity_aligned_connected",
  c("Replicate", "EstimatorMode", "RaterRMSE")
]
planned_rmse <- rmse[
  rmse$Design == "planned_connected",
  c("Replicate", "EstimatorMode", "RaterRMSE")
]
names(aligned_rmse)[[3L]] <- "AlignedRMSE"
names(planned_rmse)[[3L]] <- "PlannedRMSE"
paired_rmse <- merge(
  aligned_rmse,
  planned_rmse,
  by = c("Replicate", "EstimatorMode"),
  all = FALSE,
  sort = TRUE
)
paired_rmse$Contrast <- paired_rmse$AlignedRMSE - paired_rmse$PlannedRMSE

make_endpoint <- function(frame, endpoint_id) {
  data.frame(
    EndpointId = endpoint_id,
    Replicate = as.integer(frame$Replicate),
    Contrast = as.numeric(frame$Contrast),
    stringsAsFactors = FALSE
  )
}
r_contrasts <- rbind(
  make_endpoint(paired_rmse[paired_rmse$EstimatorMode == free_mode,], primary_id),
  make_endpoint(
    paired_rmse[paired_rmse$EstimatorMode == fixed_mode,], secondary_ids[[1L]]
  )
)

free_runs <- runs[
  runs$EstimatorMode == free_mode,
  c("Replicate", "Design", "EstimatedPopulationSD")
]
free_aligned <- free_runs[
  free_runs$Design == "ability_severity_aligned_connected",
  c("Replicate", "EstimatedPopulationSD")
]
free_planned <- free_runs[
  free_runs$Design == "planned_connected",
  c("Replicate", "EstimatedPopulationSD")
]
names(free_aligned)[[2L]] <- "AlignedSD"
names(free_planned)[[2L]] <- "PlannedSD"
paired_sd <- merge(free_aligned, free_planned, by = "Replicate", all = FALSE, sort = TRUE)
paired_sd$Contrast <- as.numeric(paired_sd$AlignedSD) - as.numeric(paired_sd$PlannedSD)
r_contrasts <- rbind(r_contrasts, make_endpoint(paired_sd, secondary_ids[[2L]]))

aligned_level <- rater[
  rater$Design == "ability_severity_aligned_connected",
  c("Replicate", "EstimatorMode", "Level", "ErrorAligned")
]
planned_level <- rater[
  rater$Design == "planned_connected",
  c("Replicate", "EstimatorMode", "Level", "ErrorAligned")
]
names(aligned_level)[[4L]] <- "AlignedError"
names(planned_level)[[4L]] <- "PlannedError"
level_pairs <- merge(
  aligned_level,
  planned_level,
  by = c("Replicate", "EstimatorMode", "Level"),
  all = FALSE,
  sort = TRUE
)
level_pairs$ErrorContrast <- level_pairs$AlignedError - level_pairs$PlannedError
level_pairs$TrueSeverity <- severity[level_pairs$Level]
slope_keys <- unique(level_pairs[, c("Replicate", "EstimatorMode")])
slope_rows <- lapply(seq_len(nrow(slope_keys)), function(index) {
  key <- slope_keys[index,]
  selected <- level_pairs[
    level_pairs$Replicate == key$Replicate &
      level_pairs$EstimatorMode == key$EstimatorMode,
  ]
  if (nrow(selected) != 4L) stop("invalid Rater compression group")
  data.frame(
    Replicate = as.integer(key$Replicate),
    EstimatorMode = key$EstimatorMode,
    Contrast = sum(selected$TrueSeverity * selected$ErrorContrast) /
      sum(severity^2),
    stringsAsFactors = FALSE
  )
})
slopes <- do.call(rbind, slope_rows)
r_contrasts <- rbind(
  r_contrasts,
  make_endpoint(slopes[slopes$EstimatorMode == free_mode,], secondary_ids[[3L]]),
  make_endpoint(slopes[slopes$EstimatorMode == fixed_mode,], secondary_ids[[4L]])
)

r_contrasts <- r_contrasts[order(r_contrasts$EndpointId, r_contrasts$Replicate),]
python_contrasts <- python_contrasts[
  order(python_contrasts$EndpointId, python_contrasts$Replicate),
]
contrast_keys_pass <- identical(
  paste(r_contrasts$EndpointId, r_contrasts$Replicate),
  paste(python_contrasts$EndpointId, python_contrasts$Replicate)
)
contrast_max <- if (contrast_keys_pass) {
  max(abs(r_contrasts$Contrast - as.numeric(python_contrasts$Contrast)))
} else {
  Inf
}

summarize_endpoint <- function(endpoint_id, alternative) {
  values <- r_contrasts$Contrast[r_contrasts$EndpointId == endpoint_id]
  values <- values[is.finite(values)]
  n <- length(values)
  mean_value <- mean(values)
  sd_value <- stats::sd(values)
  se_value <- sd_value / sqrt(n)
  statistic <- mean_value / se_value
  raw_p <- if (alternative == "less") {
    stats::pt(statistic, df = n - 1L)
  } else {
    stats::pt(statistic, df = n - 1L, lower.tail = FALSE)
  }
  half <- stats::qt(0.975, df = n - 1L) * se_value
  data.frame(
    EndpointId = endpoint_id,
    Alternative = alternative,
    FinitePairedReplicates = n,
    RequiredPairedReplicates = 100L,
    FullPairGate = n == 100L,
    MeanContrast = mean_value,
    MonteCarloSD = sd_value,
    MonteCarloSE = se_value,
    TStatistic = statistic,
    DegreesOfFreedom = n - 1L,
    RawOneSidedP = raw_p,
    Lower95 = mean_value - half,
    Upper95 = mean_value + half,
    TwoSided95HalfWidth = half,
    DirectionPass = if (alternative == "less") mean_value < 0 else mean_value > 0,
    stringsAsFactors = FALSE
  )
}

r_primary <- summarize_endpoint(primary_id, "greater")
r_primary$AdjustedP <- r_primary$RawOneSidedP
r_primary$PrecisionPass <- r_primary$TwoSided95HalfWidth <= 0.015
r_primary$DirectionConfirmed <- (
  r_primary$FullPairGate &
    r_primary$DirectionPass &
    r_primary$RawOneSidedP <= 0.05
)
r_secondary <- rbind(
  summarize_endpoint(secondary_ids[[1L]], "greater"),
  summarize_endpoint(secondary_ids[[2L]], "less"),
  summarize_endpoint(secondary_ids[[3L]], "less"),
  summarize_endpoint(secondary_ids[[4L]], "less")
)
r_secondary$HolmAdjustedP <- stats::p.adjust(r_secondary$RawOneSidedP, method = "holm")
r_secondary$PrimaryGatePass <- r_primary$DirectionConfirmed[[1L]]
r_secondary$DirectionConfirmed <- (
  r_secondary$PrimaryGatePass &
    r_secondary$FullPairGate &
    r_secondary$DirectionPass &
    r_secondary$HolmAdjustedP <= 0.05
)

numeric_fields <- c(
  "FinitePairedReplicates",
  "RequiredPairedReplicates",
  "MeanContrast",
  "MonteCarloSD",
  "MonteCarloSE",
  "TStatistic",
  "DegreesOfFreedom",
  "RawOneSidedP",
  "Lower95",
  "Upper95",
  "TwoSided95HalfWidth"
)
compare_numeric <- function(r_frame, python_frame, extra_fields = character()) {
  fields <- c(numeric_fields, extra_fields)
  max(unlist(lapply(fields, function(field) {
    abs(as.numeric(r_frame[[field]]) - as.numeric(python_frame[[field]]))
  })))
}
primary_max <- compare_numeric(r_primary, python_primary, "AdjustedP")
secondary_max <- compare_numeric(r_secondary, python_secondary, "HolmAdjustedP")
primary_logic_pass <- (
  r_primary$FullPairGate == as_flag(python_primary$FullPairGate) &
    r_primary$DirectionPass == as_flag(python_primary$DirectionPass) &
    r_primary$PrecisionPass == as_flag(python_primary$PrecisionPass) &
    r_primary$DirectionConfirmed == as_flag(python_primary$DirectionConfirmed)
)
secondary_logic_pass <- all(
  r_secondary$FullPairGate == as_flag(python_secondary$FullPairGate) &
    r_secondary$DirectionPass == as_flag(python_secondary$DirectionPass) &
    r_secondary$PrimaryGatePass == as_flag(python_secondary$PrimaryGatePass) &
    r_secondary$DirectionConfirmed == as_flag(python_secondary$DirectionConfirmed)
)

overall_max <- max(contrast_max, primary_max, secondary_max)
structural_pass <- (
  nrow(r_contrasts) == 500L &
    all(table(r_contrasts$EndpointId) == 100L) &
    contrast_keys_pass &
    nrow(r_primary) == 1L &
    nrow(r_secondary) == 4L
)
numeric_pass <- is.finite(overall_max) && overall_max <= 1e-12
verification_pass <- (
  structural_pass &
    numeric_pass &
    primary_logic_pass &
    secondary_logic_pass
)

verified <- data.frame(
  EndpointContrastRowsR = nrow(r_contrasts),
  EndpointsR = length(unique(r_contrasts$EndpointId)),
  PairsPerEndpointR = min(table(r_contrasts$EndpointId)),
  MaximumEndpointContrastDifference = contrast_max,
  MaximumPrimarySummaryDifference = primary_max,
  MaximumSecondarySummaryDifference = secondary_max,
  MaximumAbsolutePythonDifference = overall_max,
  StructuralAgreementPass = structural_pass,
  PythonNumericAgreementPass = numeric_pass,
  PrimaryLogicalAgreementPass = primary_logic_pass,
  SecondaryLogicalAgreementPass = secondary_logic_pass,
  PrimaryMeanContrastR = r_primary$MeanContrast,
  PrimaryRawOneSidedPR = r_primary$RawOneSidedP,
  PrimaryLower95R = r_primary$Lower95,
  PrimaryUpper95R = r_primary$Upper95,
  PrimaryHalfWidthR = r_primary$TwoSided95HalfWidth,
  PrimaryDirectionConfirmedR = r_primary$DirectionConfirmed,
  PrimaryPrecisionPassR = r_primary$PrecisionPass,
  SecondaryConfirmedR = sum(r_secondary$DirectionConfirmed),
  VerificationPass = verification_pass,
  stringsAsFactors = FALSE
)
write.csv(verified, output_path, row.names = FALSE, na = "")

cat(sprintf("endpoint_rows=%d\n", nrow(r_contrasts)))
cat(sprintf("primary_mean=%.17g\n", r_primary$MeanContrast))
cat(sprintf("primary_one_sided_p=%.17g\n", r_primary$RawOneSidedP))
cat(sprintf("primary_half_width=%.17g\n", r_primary$TwoSided95HalfWidth))
cat(sprintf("secondary_confirmed=%d\n", sum(r_secondary$DirectionConfirmed)))
cat(sprintf("max_abs_python_difference=%.17g\n", overall_max))
cat(sprintf("verification_pass=%s\n", tolower(as.character(verification_pass))))
if (!verification_pass) stop("R/Python informative-assignment confirmation failed")
