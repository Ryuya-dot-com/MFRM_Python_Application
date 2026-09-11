args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript known_assignment_multivector_verify.R <aggregate_dir> <output_dir>")
}

aggregate_dir <- normalizePath(args[[1L]], mustWork = TRUE)
output_dir <- normalizePath(args[[2L]], mustWork = FALSE)
if (dir.exists(output_dir)) {
  stop(paste("Refusing to overwrite R verification output:", output_dir))
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_input <- function(filename) {
  read.csv(
    file.path(aggregate_dir, filename),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

as_flag <- function(value) {
  tolower(as.character(value)) %in% c("true", "1")
}

at_gamma <- function(frame, gamma) {
  frame[abs(as.numeric(frame$Gamma) - gamma) < 1e-12, , drop = FALSE]
}

one_value <- function(frame, column, context) {
  if (nrow(frame) != 1L) {
    stop(paste("Expected one row for", context, "but found", nrow(frame)))
  }
  value <- as.numeric(frame[[column]][[1L]])
  if (!is.finite(value)) {
    stop(paste("Non-finite value for", context))
  }
  value
}

recovery <- read_input("recovery.csv")
runs <- read_input("run_ledger.csv")
python_diagnostics <- read_input("registered_preflight_diagnostics.csv")
python_summary <- read_input("registered_preflight_summary.csv")

free_mode <- "PYTHON_MML_FREE_SD_Q31"
fixed_mode <- "PYTHON_MML_FIXED_SD08_Q31"
severity <- c(R01 = -0.45, R02 = -0.15, R03 = 0.15, R04 = 0.45)

rater <- recovery[
  as_flag(recovery$IncludedInStudy) &
    recovery$Facet == "Rater" &
    recovery$EstimatorMode %in% c(free_mode, fixed_mode),
  ,
  drop = FALSE
]
rater$ErrorAligned <- as.numeric(rater$ErrorAligned)
if (nrow(rater) != 96L || any(!is.finite(rater$ErrorAligned))) {
  stop("Expected 96 finite MML Rater error rows")
}

rmse_rows <- list()
slopes <- list()
rmse_index <- 0L
slope_index <- 0L
for (vector in sort(unique(as.integer(rater$PersonVector)))) {
  for (mode in c(free_mode, fixed_mode)) {
    for (gamma in c(-0.8, 0.0, 0.8)) {
      group <- rater[
        as.integer(rater$PersonVector) == vector &
          rater$EstimatorMode == mode &
          abs(as.numeric(rater$Gamma) - gamma) < 1e-12,
        ,
        drop = FALSE
      ]
      if (nrow(group) != 4L) {
        stop(paste("Expected four Rater rows:", vector, mode, gamma))
      }
      rmse_index <- rmse_index + 1L
      rmse_rows[[rmse_index]] <- data.frame(
        PersonVector = vector,
        EstimatorMode = mode,
        Gamma = gamma,
        RaterRMSE = sqrt(mean(group$ErrorAligned ^ 2)),
        stringsAsFactors = FALSE
      )
      group$Severity <- unname(severity[group$Level])
      if (any(!is.finite(group$Severity))) {
        stop("Unknown Rater severity level")
      }
      slope_index <- slope_index + 1L
      slopes[[slope_index]] <- data.frame(
        PersonVector = vector,
        EstimatorMode = mode,
        Gamma = gamma,
        Slope = unname(coef(lm(ErrorAligned ~ Severity, data = group))[[2L]]),
        stringsAsFactors = FALSE
      )
    }
  }
}
rmse <- do.call(rbind, rmse_rows)
slope_table <- do.call(rbind, slopes)

diagnostics <- list()
diagnostic_index <- 0L
add_diagnostic <- function(id, vector, value, direction, required) {
  diagnostic_index <<- diagnostic_index + 1L
  diagnostics[[diagnostic_index]] <<- data.frame(
    DiagnosticId = id,
    PersonVector = vector,
    Value = value,
    RegisteredDirection = direction,
    RequiredDirectionalCount = required,
    ConfirmatoryClaimAllowed = FALSE,
    stringsAsFactors = FALSE
  )
}

vectors <- sort(unique(as.integer(rater$PersonVector)))
for (vector in vectors) {
  for (specification in list(
    list(
      id = "PF1_FREE_MML_SYMMETRIC_STRESS_RATER_RMSE",
      mode = free_mode,
      direction = "positive",
      required = 3L
    ),
    list(
      id = "PF2_FIXED_MML_SYMMETRIC_STRESS_RATER_RMSE",
      mode = fixed_mode,
      direction = "positive",
      required = 3L
    )
  )) {
    subset <- rmse[
      rmse$PersonVector == vector & rmse$EstimatorMode == specification$mode,
      ,
      drop = FALSE
    ]
    value <- 0.5 * (
      one_value(at_gamma(subset, -0.8), "RaterRMSE", "negative-gamma RMSE") +
        one_value(at_gamma(subset, 0.8), "RaterRMSE", "positive-gamma RMSE")
    ) - one_value(at_gamma(subset, 0.0), "RaterRMSE", "zero-gamma RMSE")
    add_diagnostic(
      specification$id, vector, value,
      specification$direction, specification$required
    )
  }

  free_runs <- runs[
    as_flag(runs$IncludedInStudy) &
      runs$EstimatorMode == free_mode &
      as.integer(runs$PersonVector) == vector,
    ,
    drop = FALSE
  ]
  sd_shift <- 0.5 * (
    one_value(at_gamma(free_runs, -0.8), "EstimatedPopulationSD", "negative-gamma SD") +
      one_value(at_gamma(free_runs, 0.8), "EstimatedPopulationSD", "positive-gamma SD")
  ) - one_value(at_gamma(free_runs, 0.0), "EstimatedPopulationSD", "zero-gamma SD")
  add_diagnostic(
    "PF3_FREE_MML_SYMMETRIC_STRESS_SD_SHIFT",
    vector, sd_shift, "negative", 3L
  )

  free_slopes <- slope_table[
    slope_table$PersonVector == vector & slope_table$EstimatorMode == free_mode,
    ,
    drop = FALSE
  ]
  slope_half_difference <- 0.5 * (
    one_value(at_gamma(free_slopes, -0.8), "Slope", "negative-gamma slope") -
      one_value(at_gamma(free_slopes, 0.8), "Slope", "positive-gamma slope")
  )
  add_diagnostic(
    "PF4_FREE_MML_DIRECTION_ALIGNED_RATER_SLOPE",
    vector, slope_half_difference, "positive", 4L
  )
}

r_diagnostics <- do.call(rbind, diagnostics)
r_diagnostics <- r_diagnostics[
  order(r_diagnostics$DiagnosticId, r_diagnostics$PersonVector),
  ,
  drop = FALSE
]
rownames(r_diagnostics) <- NULL

summary_rows <- list()
summary_index <- 0L
for (id in sort(unique(r_diagnostics$DiagnosticId))) {
  group <- r_diagnostics[r_diagnostics$DiagnosticId == id, , drop = FALSE]
  direction <- group$RegisteredDirection[[1L]]
  required <- as.integer(group$RequiredDirectionalCount[[1L]])
  directional_count <- if (direction == "positive") {
    sum(group$Value > 0)
  } else {
    sum(group$Value < 0)
  }
  mean_value <- mean(group$Value)
  mean_direction <- if (direction == "positive") mean_value > 0 else mean_value < 0
  summary_index <- summary_index + 1L
  summary_rows[[summary_index]] <- data.frame(
    DiagnosticId = id,
    N = nrow(group),
    Mean = mean_value,
    Minimum = min(group$Value),
    Maximum = max(group$Value),
    RegisteredDirection = direction,
    DirectionalCount = directional_count,
    RequiredDirectionalCount = required,
    AdvancementSignal = nrow(group) == 4L && mean_direction && directional_count >= required,
    PValueComputed = FALSE,
    ConfirmatoryClaimAllowed = FALSE,
    stringsAsFactors = FALSE
  )
}
r_summary <- do.call(rbind, summary_rows)

diagnostic_comparison <- merge(
  r_diagnostics,
  python_diagnostics,
  by = c("DiagnosticId", "PersonVector"),
  suffixes = c("R", "Python"),
  all = TRUE,
  sort = TRUE
)
if (nrow(diagnostic_comparison) != 16L) {
  stop("Expected 16 joined diagnostic rows")
}
diagnostic_difference <- max(abs(
  as.numeric(diagnostic_comparison$ValueR) -
    as.numeric(diagnostic_comparison$ValuePython)
))
diagnostic_metadata_match <- all(
  diagnostic_comparison$RegisteredDirectionR ==
    diagnostic_comparison$RegisteredDirectionPython
) && all(
  as.integer(diagnostic_comparison$RequiredDirectionalCountR) ==
    as.integer(diagnostic_comparison$RequiredDirectionalCountPython)
)

summary_comparison <- merge(
  r_summary,
  python_summary,
  by = "DiagnosticId",
  suffixes = c("R", "Python"),
  all = TRUE,
  sort = TRUE
)
if (nrow(summary_comparison) != 4L) {
  stop("Expected four joined summary rows")
}
summary_difference <- max(abs(c(
  as.numeric(summary_comparison$MeanR) - as.numeric(summary_comparison$MeanPython),
  as.numeric(summary_comparison$MinimumR) - as.numeric(summary_comparison$MinimumPython),
  as.numeric(summary_comparison$MaximumR) - as.numeric(summary_comparison$MaximumPython)
)))
summary_metadata_match <- all(
  as.integer(summary_comparison$DirectionalCountR) ==
    as.integer(summary_comparison$DirectionalCountPython)
) && all(
  as_flag(summary_comparison$AdvancementSignalR) ==
    as_flag(summary_comparison$AdvancementSignalPython)
)

tolerance <- 1e-12
passed <- is.finite(diagnostic_difference) &&
  is.finite(summary_difference) &&
  diagnostic_difference <= tolerance &&
  summary_difference <= tolerance &&
  diagnostic_metadata_match &&
  summary_metadata_match

write.csv(
  r_diagnostics,
  file.path(output_dir, "r_recomputed_diagnostics.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  r_summary,
  file.path(output_dir, "r_recomputed_summary.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  data.frame(
    Pass = passed,
    DiagnosticRows = nrow(diagnostic_comparison),
    SummaryRows = nrow(summary_comparison),
    MaximumDiagnosticAbsoluteDifference = diagnostic_difference,
    MaximumSummaryAbsoluteDifference = summary_difference,
    DiagnosticMetadataMatch = diagnostic_metadata_match,
    SummaryMetadataMatch = summary_metadata_match,
    Tolerance = tolerance,
    PValuesComputed = FALSE,
    ConfirmatoryClaimAllowed = FALSE
  ),
  file.path(output_dir, "verification.csv"),
  row.names = FALSE,
  na = ""
)

json <- sprintf(
  paste0(
    "{\n",
    "  \"schema_version\": \"known_assignment_multivector_r_verification_v1\",\n",
    "  \"pass\": %s,\n",
    "  \"diagnostic_rows\": %d,\n",
    "  \"summary_rows\": %d,\n",
    "  \"maximum_diagnostic_absolute_difference\": %.17g,\n",
    "  \"maximum_summary_absolute_difference\": %.17g,\n",
    "  \"diagnostic_metadata_match\": %s,\n",
    "  \"summary_metadata_match\": %s,\n",
    "  \"tolerance\": %.17g,\n",
    "  \"p_values_computed\": false,\n",
    "  \"confirmatory_claim_allowed\": false,\n",
    "  \"timing\": \"post-Python-endpoint implementation verification\"\n",
    "}\n"
  ),
  tolower(as.character(passed)),
  nrow(diagnostic_comparison),
  nrow(summary_comparison),
  diagnostic_difference,
  summary_difference,
  tolower(as.character(diagnostic_metadata_match)),
  tolower(as.character(summary_metadata_match)),
  tolerance
)
writeLines(json, file.path(output_dir, "assessment.json"), useBytes = TRUE)
cat(json)
if (!passed) {
  quit(save = "no", status = 2L)
}
