#!/usr/bin/env Rscript

# Independent base-R verification of the frozen H6 replication aggregate.
# The registered endpoint is rebuilt from native-precision Python threshold
# rows; FACETS display values are never used to reconstruct the statistic.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: Rscript h6_design_replication_verify.R THRESHOLDS PYTHON_PRIMARY OUTPUT")
}

thresholds <- read.csv(args[[1L]], stringsAsFactors = FALSE, check.names = FALSE)
python <- read.csv(args[[2L]], stringsAsFactors = FALSE, check.names = FALSE)

included <- tolower(as.character(thresholds$IncludedInStudy)) == "true"
eligible <- thresholds[
  thresholds$EstimatorMode == "PYTHON_JMLE" & included,
  c("RunId", "Replicate", "Design", "TruthError")
]
eligible$TruthError <- as.numeric(eligible$TruthError)

groups <- split(eligible, list(eligible$RunId, drop = TRUE))
run_rows <- lapply(groups, function(frame) {
  if (nrow(frame) != 6L || any(!is.finite(frame$TruthError))) {
    return(NULL)
  }
  data.frame(
    Replicate = as.integer(frame$Replicate[[1L]]),
    Design = frame$Design[[1L]],
    ThresholdRMSE = sqrt(mean(frame$TruthError^2)),
    stringsAsFactors = FALSE
  )
})
run_rmse <- do.call(rbind, run_rows)
complete <- run_rmse[run_rmse$Design == "complete", c("Replicate", "ThresholdRMSE")]
planned <- run_rmse[
  run_rmse$Design == "planned_connected", c("Replicate", "ThresholdRMSE")
]
names(complete)[[2L]] <- "CompleteRMSE"
names(planned)[[2L]] <- "PlannedRMSE"
paired <- merge(complete, planned, by = "Replicate", all = FALSE, sort = TRUE)
contrast <- paired$PlannedRMSE - paired$CompleteRMSE

n <- length(contrast)
estimate <- mean(contrast)
sd_value <- stats::sd(contrast)
se_value <- sd_value / sqrt(n)
statistic <- estimate / se_value
p_one_sided <- stats::pt(statistic, df = n - 1L, lower.tail = FALSE)
half_width <- stats::qt(0.975, df = n - 1L) * se_value
lower <- estimate - half_width
upper <- estimate + half_width
replicated <- n == 100L && estimate > 0 && p_one_sided <= 0.05
precision <- n == 100L && half_width <= 0.015

comparisons <- c(
  FinitePairs = abs(n - python$FinitePairs[[1L]]),
  MeanContrast = abs(estimate - python$MeanContrast[[1L]]),
  SDContrast = abs(sd_value - python$SDContrast[[1L]]),
  SEContrast = abs(se_value - python$SEContrast[[1L]]),
  TStatistic = abs(statistic - python$TStatistic[[1L]]),
  OneSidedP = abs(p_one_sided - python$OneSidedP[[1L]]),
  CI95Low = abs(lower - python$CI95Low[[1L]]),
  CI95High = abs(upper - python$CI95High[[1L]]),
  CI95HalfWidth = abs(half_width - python$CI95HalfWidth[[1L]])
)
numeric_pass <- all(is.finite(comparisons)) && max(comparisons) <= 1e-12
logical_pass <- (
  identical(replicated, tolower(as.character(python$ReplicationDecision[[1L]])) == "true") &&
  identical(precision, tolower(as.character(python$PrecisionQualification[[1L]])) == "true") &&
  n == 100L
)

verified <- data.frame(
  EndpointId = "H6R_JMLE_NORMAL_THRESHOLD_RMSE_PLANNED_GT_COMPLETE",
  FinitePairsR = n,
  MeanContrastR = estimate,
  SDContrastR = sd_value,
  SEContrastR = se_value,
  TStatisticR = statistic,
  OneSidedPR = p_one_sided,
  CI95LowR = lower,
  CI95HighR = upper,
  CI95HalfWidthR = half_width,
  ReplicationDecisionR = replicated,
  PrecisionQualificationR = precision,
  MaximumAbsolutePythonDifference = max(comparisons),
  PythonNumericAgreementPass = numeric_pass,
  PythonLogicalAgreementPass = logical_pass,
  stringsAsFactors = FALSE
)
write.csv(verified, args[[3L]], row.names = FALSE, na = "")

cat(sprintf("finite_pairs=%d\n", n))
cat(sprintf("mean_contrast=%.17g\n", estimate))
cat(sprintf("one_sided_p=%.17g\n", p_one_sided))
cat(sprintf("ci95_half_width=%.17g\n", half_width))
cat(sprintf("max_abs_python_difference=%.17g\n", max(comparisons)))
cat(sprintf("numeric_agreement_pass=%s\n", tolower(as.character(numeric_pass))))
cat(sprintf("logical_agreement_pass=%s\n", tolower(as.character(logical_pass))))
if (!numeric_pass || !logical_pass) {
  stop("R/Python H6 replication verification failed")
}

