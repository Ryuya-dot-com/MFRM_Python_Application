#!/usr/bin/env Rscript

# Independent base-R verification of the frozen Python confirmatory statistics.
# This script verifies paired means, t intervals/tests, Holm adjustment, and the
# registered 100-pair gate. It does not create any new endpoint or inference.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: Rscript estimand_distribution_confirmatory_verify.R CONTRASTS PYTHON_RESULTS OUTPUT")
}

contrasts <- read.csv(args[[1L]], stringsAsFactors = FALSE, check.names = FALSE)
python <- read.csv(args[[2L]], stringsAsFactors = FALSE, check.names = FALSE)

registered <- data.frame(
  EndpointId = c(
    "H1_FREE_SD_RIGHT_SKEW_PLANNED_LT_NORMAL",
    "H2_FREE_SD_HEAVY_TAIL_PLANNED_LT_NORMAL",
    "H3_JMLE_MIXTURE_PLANNED_C02_STEP2_LT_NORMAL",
    "H4_FREE_MML_MIXTURE_PLANNED_C02_STEP2_LT_NORMAL",
    "H5_EXACT_CMLE_MIXTURE_PLANNED_C02_STEP2_LT_NORMAL",
    "H6_JMLE_NORMAL_THRESHOLD_RMSE_PLANNED_GT_COMPLETE"
  ),
  Alternative = c("less", "less", "less", "less", "less", "greater"),
  stringsAsFactors = FALSE
)

rows <- vector("list", nrow(registered))
for (index in seq_len(nrow(registered))) {
  endpoint <- registered$EndpointId[[index]]
  alternative <- registered$Alternative[[index]]
  values <- contrasts$Contrast[contrasts$EndpointId == endpoint]
  values <- values[is.finite(values)]
  n <- length(values)
  estimate <- mean(values)
  monte_carlo_sd <- stats::sd(values)
  monte_carlo_se <- monte_carlo_sd / sqrt(n)
  statistic <- estimate / monte_carlo_se
  raw_p <- if (alternative == "less") {
    stats::pt(statistic, df = n - 1L)
  } else {
    stats::pt(statistic, df = n - 1L, lower.tail = FALSE)
  }
  half_width <- stats::qt(0.975, df = n - 1L) * monte_carlo_se
  rows[[index]] <- data.frame(
    EndpointId = endpoint,
    Alternative = alternative,
    FinitePairedReplicatesR = n,
    MeanContrastR = estimate,
    MonteCarloSDR = monte_carlo_sd,
    MonteCarloSER = monte_carlo_se,
    TStatisticR = statistic,
    Lower95R = estimate - half_width,
    Upper95R = estimate + half_width,
    TwoSided95HalfWidthR = half_width,
    RawOneSidedPR = raw_p,
    stringsAsFactors = FALSE
  )
}

verified <- do.call(rbind, rows)
verified$HolmAdjustedPR <- stats::p.adjust(verified$RawOneSidedPR, method = "holm")
verified$FullPairGateR <- verified$FinitePairedReplicatesR == 100L
verified$DirectionPassR <- ifelse(
  verified$Alternative == "less",
  verified$MeanContrastR < 0,
  verified$MeanContrastR > 0
)
verified$DirectionConfirmedR <- (
  verified$FullPairGateR & verified$DirectionPassR & verified$HolmAdjustedPR <= 0.05
)

python <- python[match(verified$EndpointId, python$EndpointId), ]
if (any(is.na(python$EndpointId)) || !identical(verified$EndpointId, python$EndpointId)) {
  stop("Python result endpoint identities do not match the registered order")
}

comparisons <- c(
  MeanContrast = max(abs(verified$MeanContrastR - python$MeanContrast)),
  MonteCarloSD = max(abs(verified$MonteCarloSDR - python$MonteCarloSD)),
  MonteCarloSE = max(abs(verified$MonteCarloSER - python$MonteCarloSE)),
  TStatistic = max(abs(verified$TStatisticR - python$TStatistic)),
  Lower95 = max(abs(verified$Lower95R - python$Lower95)),
  Upper95 = max(abs(verified$Upper95R - python$Upper95)),
  TwoSided95HalfWidth = max(abs(verified$TwoSided95HalfWidthR - python$TwoSided95HalfWidth)),
  RawOneSidedP = max(abs(verified$RawOneSidedPR - python$RawOneSidedP)),
  HolmAdjustedP = max(abs(verified$HolmAdjustedPR - python$HolmAdjustedP))
)
logical_pass <- (
  identical(as.integer(verified$FinitePairedReplicatesR), as.integer(python$FinitePairedReplicates)) &&
  identical(as.logical(verified$FullPairGateR), as.logical(python$FullPairGate)) &&
  identical(as.logical(verified$DirectionPassR), as.logical(python$DirectionPass)) &&
  identical(as.logical(verified$DirectionConfirmedR), as.logical(python$DirectionConfirmed))
)
numeric_pass <- all(is.finite(comparisons)) && max(comparisons) <= 1e-12
verified$PythonNumericAgreementPass <- numeric_pass
verified$PythonLogicalAgreementPass <- logical_pass
write.csv(verified, args[[3L]], row.names = FALSE, na = "")

cat(sprintf("endpoints=%d\n", nrow(verified)))
cat(sprintf("max_abs_numeric_difference=%.17g\n", max(comparisons)))
cat(sprintf("numeric_agreement_pass=%s\n", tolower(as.character(numeric_pass))))
cat(sprintf("logical_agreement_pass=%s\n", tolower(as.character(logical_pass))))
if (!numeric_pass || !logical_pass) {
  stop("R/Python confirmatory statistic verification failed")
}
