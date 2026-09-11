#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript known_assignment_response_verify.R <aggregate_dir> <output_dir>")
}

aggregate_dir <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
output_dir <- args[[2]]
if (dir.exists(output_dir)) {
  stop(paste("Refusing to overwrite output directory:", output_dir))
}
dir.create(output_dir, recursive = FALSE)

recovery <- read.csv(
  file.path(aggregate_dir, "recovery.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
python_summary <- read.csv(
  file.path(aggregate_dir, "gamma_contrast_summary.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

included <- recovery$IncludedInStudy %in% c(TRUE, "True", "TRUE", 1, "1")
rater <- recovery[
  included & recovery$Facet == "Rater" & is.finite(recovery$ErrorAligned),
  c("RunId", "Gamma", "Replicate", "EstimatorMode", "ErrorAligned")
]

keys <- unique(rater[, c("RunId", "Gamma", "Replicate", "EstimatorMode")])
loss_rows <- vector("list", nrow(keys))
for (index in seq_len(nrow(keys))) {
  key <- keys[index, ]
  selected <-
    rater$RunId == key$RunId &
    rater$EstimatorMode == key$EstimatorMode
  values <- rater$ErrorAligned[selected]
  loss_rows[[index]] <- data.frame(
    RunId = key$RunId,
    Gamma = key$Gamma,
    Replicate = key$Replicate,
    EstimatorMode = key$EstimatorMode,
    RMSE = sqrt(mean(values ^ 2)),
    stringsAsFactors = FALSE
  )
}
loss <- do.call(rbind, loss_rows)
neutral <- loss[loss$Gamma == 0, c("Replicate", "EstimatorMode", "RMSE")]
names(neutral)[names(neutral) == "RMSE"] <- "NeutralRMSE"

contrast_parts <- list()
part_index <- 1L
for (stress_gamma in c(-0.8, 0.8)) {
  stress <- loss[loss$Gamma == stress_gamma, c("Replicate", "EstimatorMode", "RMSE")]
  names(stress)[names(stress) == "RMSE"] <- "StressRMSE"
  paired <- merge(
    stress,
    neutral,
    by = c("Replicate", "EstimatorMode"),
    all = FALSE,
    sort = TRUE
  )
  paired$StressGamma <- stress_gamma
  paired$ContrastStressMinusNeutral <- paired$StressRMSE - paired$NeutralRMSE
  contrast_parts[[part_index]] <- paired
  part_index <- part_index + 1L
}
contrasts <- do.call(rbind, contrast_parts)

summary_rows <- list()
row_index <- 1L
for (mode in sort(unique(contrasts$EstimatorMode))) {
  for (stress_gamma in c(-0.8, 0.8)) {
    values <- contrasts$ContrastStressMinusNeutral[
      contrasts$EstimatorMode == mode & contrasts$StressGamma == stress_gamma
    ]
    summary_rows[[row_index]] <- data.frame(
      EstimatorMode = mode,
      StressGamma = stress_gamma,
      N = length(values),
      RMeanContrast = mean(values),
      RSDContrast = sd(values),
      stringsAsFactors = FALSE
    )
    row_index <- row_index + 1L
  }
}
r_summary <- do.call(rbind, summary_rows)
python_rater <- python_summary[
  python_summary$RecoveryDomain == "Facet:Rater" & python_summary$Metric == "RMSE",
  c("EstimatorMode", "StressGamma", "N", "MeanContrastStressMinusNeutral", "SDContrastStressMinusNeutral")
]
comparison <- merge(
  r_summary,
  python_rater,
  by = c("EstimatorMode", "StressGamma"),
  all = TRUE,
  sort = TRUE
)
comparison$MeanDifference <- comparison$RMeanContrast - comparison$MeanContrastStressMinusNeutral
comparison$SDDifference <- comparison$RSDContrast - comparison$SDContrastStressMinusNeutral
comparison$NMatch <- comparison$N.x == comparison$N.y

maximum_mean_difference <- max(abs(comparison$MeanDifference))
maximum_sd_difference <- max(abs(comparison$SDDifference))
passed <-
  nrow(comparison) == 10L &&
  all(comparison$NMatch) &&
  maximum_mean_difference <= 1e-12 &&
  maximum_sd_difference <= 1e-12

write.csv(loss, file.path(output_dir, "r_rater_rmse_by_run.csv"), row.names = FALSE)
write.csv(contrasts, file.path(output_dir, "r_paired_rater_rmse_contrasts.csv"), row.names = FALSE)
write.csv(comparison, file.path(output_dir, "r_python_comparison.csv"), row.names = FALSE)
writeLines(
  c(
    paste0("decision=", ifelse(passed, "pass", "fail")),
    paste0("comparison_rows=", nrow(comparison)),
    paste0("maximum_mean_difference=", format(maximum_mean_difference, digits = 17)),
    paste0("maximum_sd_difference=", format(maximum_sd_difference, digits = 17)),
    "scope=Independent base-R reconstruction of within-estimator Rater RMSE gamma contrasts only.",
    "confirmatory_claim=false",
    paste0("R.version=", R.version.string)
  ),
  file.path(output_dir, "assessment.txt")
)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))

if (!passed) {
  quit(status = 1L)
}
