#!/usr/bin/env Rscript

# Frozen fixed-theta Person MnSq adapter for sirt::pcm.fit.

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
for (package in c("TAM", "sirt", "digest", "jsonlite")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Required adapter dependency unavailable: ", package, call. = FALSE)
  }
}

input_dir <- normalizePath(args$input, mustWork = TRUE)
repo_dir <- normalizePath(args$repo, mustWork = TRUE)
output_dir <- normalizePath(args$output, mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

plan <- jsonlite::fromJSON(
  file.path(repo_dir, "validation/cmle_wle_person_fit_plan_20260810.json"),
  simplifyVector = FALSE
)
hash_function <- function(fun) {
  source <- paste(deparse(fun), collapse = "\n")
  digest::digest(source, algo = "sha256", serialize = FALSE)
}
identity <- data.frame(
  Package = c("TAM", "TAM", "sirt"),
  Version = c(
    as.character(utils::packageVersion("TAM")),
    as.character(utils::packageVersion("TAM")),
    as.character(utils::packageVersion("sirt"))
  ),
  Function = c("TAM::tam.personfit", "TAM::tam.jml.fit", "sirt::pcm.fit"),
  DeparseSHA256 = c(
    hash_function(TAM::tam.personfit),
    hash_function(TAM::tam.jml.fit),
    hash_function(sirt::pcm.fit)
  ),
  ExpectedSHA256 = c(
    plan$reference_identity$TAM_tam_personfit_deparsed_sha256,
    plan$reference_identity$TAM_tam_jml_fit_deparsed_sha256,
    plan$reference_identity$sirt_pcm_fit_deparsed_sha256
  ),
  stringsAsFactors = FALSE
)
identity$VersionPassed <- c(
  identity$Version[[1L]] == plan$reference_identity$TAM_version,
  identity$Version[[2L]] == plan$reference_identity$TAM_version,
  identity$Version[[3L]] == plan$reference_identity$sirt_version
)
identity$SourcePassed <- identity$DeparseSHA256 == identity$ExpectedSHA256
identity$IdentityPassed <- identity$VersionPassed & identity$SourcePassed
if (!all(identity$IdentityPassed)) {
  stop("Loaded TAM/sirt function identity differs from the frozen plan.", call. = FALSE)
}

manifest <- jsonlite::fromJSON(
  file.path(input_dir, "fixture_manifest.json"), simplifyVector = FALSE
)
for (entry in manifest$files) {
  path <- file.path(repo_dir, entry$path)
  if (!file.exists(path)) stop("Manifest file is absent: ", entry$path, call. = FALSE)
  actual <- digest::digest(file = path, algo = "sha256", serialize = FALSE)
  if (!identical(actual, entry$sha256)) {
    stop("Manifest hash mismatch: ", entry$path, call. = FALSE)
  }
}

calibration <- utils::read.csv(
  file.path(input_dir, "sirt_calibration.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
responses <- utils::read.csv(
  file.path(input_dir, "sirt_responses.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
theta_table <- utils::read.csv(
  file.path(input_dir, "sirt_theta.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

person_rows <- list()
observation_rows <- list()
for (model in unlist(manifest$models, use.names = FALSE)) {
  cal <- calibration[calibration$Model == model, , drop = FALSE]
  rsp <- responses[responses$Model == model, , drop = FALSE]
  th <- theta_table[theta_table$Model == model, , drop = FALSE]
  items <- unique(cal$VirtualUnit)
  persons <- th$Person
  categories <- sort(unique(as.integer(cal$InternalCategory)))
  if (!identical(categories, seq_len(max(categories)))) {
    stop("sirt b categories must be consecutive from one: ", model, call. = FALSE)
  }
  K <- max(categories)
  b <- matrix(
    NA_real_, nrow = length(items), ncol = K,
    dimnames = list(items, paste0("b", seq_len(K)))
  )
  for (row in seq_len(nrow(cal))) {
    b[match(cal$VirtualUnit[[row]], items), cal$InternalCategory[[row]]] <-
      as.numeric(cal$CumulativeDifficulty[[row]])
  }
  if (any(!is.finite(b))) stop("Incomplete sirt calibration: ", model, call. = FALSE)

  dat <- matrix(
    NA_integer_, nrow = length(persons), ncol = length(items),
    dimnames = list(persons, items)
  )
  keys <- paste(rsp$Person, rsp$VirtualUnit, sep = "\r")
  if (anyDuplicated(keys)) stop("Duplicate Person-unit response: ", model, call. = FALSE)
  for (row in seq_len(nrow(rsp))) {
    dat[match(rsp$Person[[row]], persons), match(rsp$VirtualUnit[[row]], items)] <-
      as.integer(rsp$ObservedInternalCategory[[row]])
  }
  if (any(rowSums(!is.na(dat)) == 0L)) stop("Person with no response: ", model, call. = FALSE)

  fit <- sirt::pcm.fit(b = b, theta = as.numeric(th$WLEEstimate), dat = dat)
  person_rows[[model]] <- data.frame(
    Model = model,
    Person = persons,
    SirtOutfit = as.numeric(fit$personfit$outfit),
    SirtOutfitT = as.numeric(fit$personfit$outfit.t),
    SirtInfit = as.numeric(fit$personfit$infit),
    SirtInfitT = as.numeric(fit$personfit$infit.t),
    stringsAsFactors = FALSE
  )

  # Independent row-level reconstruction of the probability moments used by
  # pcm.fit.  This prevents aggregate MnSq agreement from hiding a category-
  # kernel mismatch.
  score_vec <- 0:K
  rows <- vector("list", nrow(rsp))
  for (row in seq_len(nrow(rsp))) {
    person_index <- match(rsp$Person[[row]], persons)
    item_index <- match(rsp$VirtualUnit[[row]], items)
    logits <- as.numeric(th$WLEEstimate[[person_index]]) * score_vec
    logits[-1L] <- logits[-1L] - b[item_index, ]
    probabilities <- exp(logits - max(logits))
    probabilities <- probabilities / sum(probabilities)
    expected <- sum(score_vec * probabilities)
    variance <- sum((score_vec - expected)^2 * probabilities)
    fourth <- sum((score_vec - expected)^4 * probabilities)
    observed <- as.integer(rsp$ObservedInternalCategory[[row]])
    residual <- observed - expected
    rows[[row]] <- data.frame(
      Model = model,
      Person = rsp$Person[[row]],
      VirtualUnit = rsp$VirtualUnit[[row]],
      ObservedInternalCategory = observed,
      SirtExpectedInternalCategory = expected,
      SirtVariance = variance,
      SirtFourthCentralMoment = fourth,
      SirtResidual = residual,
      SirtStandardizedSquaredResidual = residual^2 / variance,
      SirtObservedProbability = probabilities[[observed + 1L]],
      stringsAsFactors = FALSE
    )
  }
  observation_rows[[model]] <- do.call(rbind, rows)
}

person_output <- do.call(rbind, person_rows)
rownames(person_output) <- NULL
observation_output <- do.call(rbind, observation_rows)
rownames(observation_output) <- NULL
utils::write.csv(
  person_output,
  file.path(output_dir, "sirt_person_fit.csv"),
  row.names = FALSE,
  quote = TRUE,
  na = ""
)
utils::write.csv(
  observation_output,
  file.path(output_dir, "sirt_observation_moments.csv"),
  row.names = FALSE,
  quote = TRUE,
  na = ""
)
utils::write.csv(
  identity,
  file.path(output_dir, "reference_identity.csv"),
  row.names = FALSE,
  quote = TRUE
)
