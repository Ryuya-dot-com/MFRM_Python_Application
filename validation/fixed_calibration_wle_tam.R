#!/usr/bin/env Rscript

# TAM fixed-calibration Warm-WLE adapter for the frozen parity fixtures.

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

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (is.null(args$input) || is.null(args$output) || is.null(args$repo)) {
  stop("Required arguments: --input, --output, --repo", call. = FALSE)
}
for (package in c("TAM", "digest", "jsonlite")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Required adapter dependency unavailable: ", package, call. = FALSE)
  }
}

input_dir <- normalizePath(args$input, mustWork = TRUE)
repo_dir <- normalizePath(args$repo, mustWork = TRUE)
output_dir <- normalizePath(args$output, mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

manifest_path <- file.path(input_dir, "fixture_manifest.json")
manifest <- jsonlite::fromJSON(manifest_path, simplifyVector = FALSE)
for (entry in manifest$files) {
  path <- file.path(repo_dir, entry$path)
  if (!file.exists(path)) stop("Manifest file is absent: ", entry$path, call. = FALSE)
  actual <- digest::digest(file = path, algo = "sha256", serialize = FALSE)
  if (!identical(actual, entry$sha256)) {
    stop("Manifest hash mismatch: ", entry$path, call. = FALSE)
  }
}

amendment <- jsonlite::fromJSON(
  file.path(repo_dir, "validation/fixed_calibration_wle_plan_amendment_20260809.json"),
  simplifyVector = FALSE
)
tam_version <- as.character(utils::packageVersion("TAM"))
tam_source <- paste(deparse(TAM::tam.mml.wle2), collapse = "\n")
tam_source_hash <- digest::digest(tam_source, algo = "sha256", serialize = FALSE)
tam_body_hash <- digest::digest(
  list(formals = formals(TAM::tam.mml.wle2), body = body(TAM::tam.mml.wle2)),
  algo = "sha256",
  serialize = TRUE
)
identity_passed <- identical(tam_version, amendment$reproducible_source_identity$version) &&
  identical(tam_source_hash, amendment$reproducible_source_identity$sha256) &&
  identical(tam_body_hash, amendment$reproducible_source_identity$secondary_sha256)
if (!identity_passed) {
  stop("Loaded TAM WLE code identity differs from the frozen amendment.", call. = FALSE)
}

calibration <- utils::read.csv(
  file.path(input_dir, "fixture_calibration.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
responses <- utils::read.csv(
  file.path(input_dir, "fixture_responses.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA")
)

expected_cases <- unlist(manifest$case_names, use.names = FALSE)
if (!identical(unique(calibration$Case), expected_cases) ||
    !identical(unique(responses$Case), expected_cases)) {
  stop("Fixture cases or case ordering differ from the frozen manifest.", call. = FALSE)
}

result_rows <- list()
for (case_name in expected_cases) {
  case_calibration <- calibration[calibration$Case == case_name, , drop = FALSE]
  case_responses <- responses[responses$Case == case_name, , drop = FALSE]
  items <- unique(case_calibration$Item)
  persons <- unique(case_responses$Person)
  categories <- sort(unique(as.integer(case_calibration$Category)))
  if (!identical(categories, 0:max(categories))) {
    stop("Categories must be consecutive from zero: ", case_name, call. = FALSE)
  }
  nitems <- length(items)
  maxK <- length(categories)
  AXsi <- matrix(NA_real_, nrow = nitems, ncol = maxK)
  B <- array(NA_real_, dim = c(nitems, maxK, 1L))
  for (row in seq_len(nrow(case_calibration))) {
    i <- match(case_calibration$Item[[row]], items)
    k <- as.integer(case_calibration$Category[[row]]) + 1L
    AXsi[i, k] <- as.numeric(case_calibration$Intercept[[row]])
    B[i, k, 1L] <- as.numeric(case_calibration$Slope[[row]])
  }
  if (any(!is.finite(AXsi)) || any(!is.finite(B))) {
    stop("Incomplete or non-finite calibration: ", case_name, call. = FALSE)
  }

  resp <- matrix(
    NA_integer_,
    nrow = length(persons),
    ncol = nitems,
    dimnames = list(persons, items)
  )
  keys <- paste(case_responses$Person, case_responses$Item, sep = "\r")
  if (anyDuplicated(keys)) stop("Duplicate Person-item response: ", case_name, call. = FALSE)
  for (row in seq_len(nrow(case_responses))) {
    p <- match(case_responses$Person[[row]], persons)
    i <- match(case_responses$Item[[row]], items)
    value <- case_responses$Observed[[row]]
    if (!is.na(value)) resp[p, i] <- as.integer(value)
  }
  if (any(rowSums(!is.na(resp)) == 0L)) {
    stop("A fixture Person has no observed response: ", case_name, call. = FALSE)
  }

  tam_object <- list(
    B = B,
    A = NULL,
    nitems = nitems,
    xsi = NULL,
    AXsi = AXsi,
    resp = resp,
    resp.ind = 1L - is.na(resp),
    pweights = rep(1, nrow(resp)),
    pid = persons
  )
  fit <- suppressWarnings(TAM::tam.mml.wle2(
    tam_object,
    WLE = TRUE,
    Msteps = 200L,
    convM = 1e-12,
    progress = FALSE,
    output.prob = FALSE
  ))
  result_rows[[case_name]] <- data.frame(
    Case = case_name,
    Person = as.character(fit$pid),
    EstimateTAM = as.numeric(fit$theta),
    StandardErrorTAM = as.numeric(fit$error),
    ObservedRowsTAM = as.integer(fit$N.items),
    PersonScoreTAM = as.numeric(fit$PersonScores),
    PersonMaximumTAM = as.numeric(fit$PersonMax),
    stringsAsFactors = FALSE
  )
}

tam_output <- do.call(rbind, result_rows)
rownames(tam_output) <- NULL
utils::write.csv(
  tam_output,
  file.path(output_dir, "tam_wle.csv"),
  row.names = FALSE,
  na = "",
  quote = TRUE
)
identity <- data.frame(
  Package = "TAM",
  Version = tam_version,
  Function = "TAM::tam.mml.wle2",
  DeparseSHA256 = tam_source_hash,
  FormalsBodySHA256 = tam_body_hash,
  FrozenAmendmentIdentityPassed = identity_passed,
  OriginalPlanHashReproducible = FALSE,
  stringsAsFactors = FALSE
)
utils::write.csv(
  identity,
  file.path(output_dir, "tam_identity.csv"),
  row.names = FALSE,
  quote = TRUE
)
