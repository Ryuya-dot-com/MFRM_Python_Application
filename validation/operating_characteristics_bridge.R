#!/usr/bin/env Rscript

# Validate and import the exact Python-generated operating-characteristics
# bundle before any mfrmr/TAM/immer/sirt adapter is allowed to fit it.

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
if (!nzchar(args$input %||% "")) {
  stop("Missing required argument: --input", call. = FALSE)
}
if (nzchar(args[["mfrmr-lib"]] %||% "")) {
  .libPaths(c(normalizePath(args[["mfrmr-lib"]], mustWork = TRUE), .libPaths()))
}
if (!requireNamespace("digest", quietly = TRUE)) {
  stop("Required bridge dependency unavailable: digest", call. = FALSE)
}

input_dir <- normalizePath(args$input, mustWork = TRUE)
output_dir <- normalizePath(args$output %||% args$input, mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_bundle_csv <- function(filename) {
  path <- file.path(input_dir, filename)
  if (!file.exists(path)) stop("Missing bridge file: ", filename, call. = FALSE)
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

manifest <- read_bundle_csv("manifest.csv")
inventory <- read_bundle_csv("generated_bundle_files.csv")
identity <- read_bundle_csv("generated_data_identity.csv")
ratings <- read_bundle_csv("generated_ratings.csv")
facet_truth <- read_bundle_csv("generated_facet_truth.csv")
anchors <- read_bundle_csv("generated_anchors.csv")

checks <- list()
check_index <- 0L
record_check <- function(check, passed, evidence) {
  check_index <<- check_index + 1L
  checks[[check_index]] <<- data.frame(
    Check = check,
    Passed = isTRUE(passed),
    Evidence = as.character(evidence),
    stringsAsFactors = FALSE
  )
  invisible(NULL)
}

required_columns <- function(frame, columns, label) {
  missing <- setdiff(columns, names(frame))
  record_check(
    paste0(label, "_columns"),
    length(missing) == 0L,
    if (length(missing)) paste("missing", paste(missing, collapse = ";")) else "complete"
  )
}

required_columns(manifest, c("RunId", "Categories"), "manifest")
required_columns(
  identity,
  c(
    "RunId", "DataId", "FitInputId", "RatingsRows", "FacetTruthRows",
    "AnchorRows", "RatingsFingerprint", "FacetTruthFingerprint", "AnchorFingerprint"
  ),
  "identity"
)
required_columns(
  ratings,
  c("RunId", "Person", "Rater", "Task", "Criterion", "Score"),
  "ratings"
)
required_columns(facet_truth, c("RunId", "Facet", "Level", "Truth"), "facet_truth")
required_columns(anchors, c("RunId", "Facet", "Level", "Anchor"), "anchors")

hash_ok <- logical(nrow(inventory))
for (i in seq_len(nrow(inventory))) {
  path <- file.path(input_dir, as.character(inventory$File[[i]]))
  actual <- if (file.exists(path)) {
    digest::digest(file = path, algo = "sha256", serialize = FALSE)
  } else {
    "missing"
  }
  hash_ok[[i]] <- identical(tolower(actual), tolower(as.character(inventory$SHA256[[i]])))
}
record_check(
  "bundle_file_sha256",
  length(hash_ok) > 0L && all(hash_ok),
  paste(sum(hash_ok), "of", length(hash_ok), "byte-level hashes match")
)

manifest_ids <- sort(unique(as.character(manifest$RunId)))
identity_ids <- sort(unique(as.character(identity$RunId)))
record_check(
  "manifest_run_identity",
  identical(manifest_ids, identity_ids) && !anyDuplicated(identity$RunId),
  paste("manifest", length(manifest_ids), "identity", length(identity_ids))
)

check_counts <- function(frame, count_column, label) {
  tab <- table(as.character(frame$RunId))
  observed <- as.integer(tab[match(as.character(identity$RunId), names(tab))])
  observed[is.na(observed)] <- 0L
  expected <- as.integer(identity[[count_column]])
  record_check(
    paste0(label, "_run_counts"),
    identical(observed, expected),
    paste("rows", nrow(frame), "expected", sum(expected))
  )
}
check_counts(ratings, "RatingsRows", "ratings")
check_counts(facet_truth, "FacetTruthRows", "facet_truth")
check_counts(anchors, "AnchorRows", "anchors")

rating_key <- paste(
  ratings$RunId, ratings$Person, ratings$Rater, ratings$Task,
  ratings$Criterion, sep = "\r"
)
truth_key <- paste(facet_truth$RunId, facet_truth$Facet, facet_truth$Level, sep = "\r")
anchor_key <- paste(anchors$RunId, anchors$Facet, anchors$Level, sep = "\r")
record_check("ratings_response_unit_unique", !anyDuplicated(rating_key), "RunId x response unit")
record_check("facet_truth_unique", !anyDuplicated(truth_key), "RunId x Facet x Level")
record_check("anchors_unique", !anyDuplicated(anchor_key), "RunId x Facet x Level")

score_contract <- merge(
  ratings[c("RunId", "Score")],
  manifest[c("RunId", "Categories")],
  by = "RunId",
  all.x = TRUE,
  sort = FALSE
)
score_numeric <- suppressWarnings(as.numeric(score_contract$Score))
category_numeric <- suppressWarnings(as.numeric(score_contract$Categories))
score_ok <- is.finite(score_numeric) & score_numeric == floor(score_numeric) &
  score_numeric >= 0 & score_numeric < category_numeric
record_check(
  "rating_category_contract",
  nrow(score_contract) == nrow(ratings) && all(score_ok),
  paste(sum(score_ok), "of", length(score_ok), "scores within registered support")
)

clean <- identity[identity$ConditionId %in% grep(
  "^balanced_large_anchors__", identity$ConditionId, value = TRUE
), , drop = FALSE]
drift <- identity[identity$ConditionId %in% grep(
  "^anchor_drift__", identity$ConditionId, value = TRUE
), , drop = FALSE]
paired <- merge(clean, drift, by = c("TruthBias", "Replicate"), suffixes = c("_clean", "_drift"))
paired_data_ok <- nrow(paired) > 0L &&
  all(paired$DataId_clean == paired$DataId_drift) &&
  all(paired$RatingsFingerprint_clean == paired$RatingsFingerprint_drift) &&
  all(paired$FacetTruthFingerprint_clean == paired$FacetTruthFingerprint_drift)
paired_anchor_ok <- nrow(paired) > 0L &&
  all(paired$AnchorFingerprint_clean != paired$AnchorFingerprint_drift) &&
  all(paired$FitInputId_clean != paired$FitInputId_drift)
record_check(
  "paired_anchor_generated_data",
  paired_data_ok,
  paste(nrow(paired), "clean/drift pairs share ratings and truth identity")
)
record_check(
  "paired_anchor_fit_input_difference",
  paired_anchor_ok,
  paste(nrow(paired), "clean/drift pairs retain different anchor input identity")
)

validation <- do.call(rbind, checks)
bundle_inventory_sha256 <- digest::digest(
  file = file.path(input_dir, "generated_bundle_files.csv"),
  algo = "sha256",
  serialize = FALSE
)
validation$BundleInventorySHA256 <- bundle_inventory_sha256
utils::write.csv(
  validation,
  file.path(output_dir, "bridge_validation_r.csv"),
  row.names = FALSE
)

engines <- c("mfrmr", "TAM", "immer", "sirt")
primary_functions <- list(
  mfrmr = c("fit_mfrm"),
  TAM = c("tam.mml.mfr", "tam.jml"),
  immer = c("immer_jml", "immer_cml"),
  sirt = c("rm.facets")
)
function_hash <- function(package, functions) {
  if (!requireNamespace(package, quietly = TRUE)) return(NA_character_)
  namespace <- asNamespace(package)
  hashes <- vapply(functions, function(name) {
    if (!exists(name, envir = namespace, inherits = FALSE)) return(NA_character_)
    fn <- get(name, envir = namespace, inherits = FALSE)
    digest::digest(list(formals = formals(fn), body = body(fn)), algo = "sha256", serialize = TRUE)
  }, character(1))
  digest::digest(hashes, algo = "sha256", serialize = TRUE)
}
namespace_hash <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) return(NA_character_)
  namespace <- asNamespace(package)
  names <- sort(ls(namespace, all.names = TRUE))
  hashes <- vapply(names, function(name) {
    value <- get(name, envir = namespace, inherits = FALSE)
    if (!is.function(value)) return(NA_character_)
    digest::digest(list(formals = formals(value), body = body(value)), algo = "sha256", serialize = TRUE)
  }, character(1))
  digest::digest(hashes[!is.na(hashes)], algo = "sha256", serialize = TRUE)
}
native_library_hash <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) return(NA_character_)
  library_dir <- system.file("libs", package = package)
  if (!nzchar(library_dir) || !dir.exists(library_dir)) return(NA_character_)
  files <- sort(list.files(library_dir, recursive = TRUE, full.names = TRUE))
  files <- files[file.info(files)$isdir %in% FALSE]
  if (!length(files)) return(NA_character_)
  hashes <- vapply(files, function(path) {
    digest::digest(file = path, algo = "sha256", serialize = FALSE)
  }, character(1))
  digest::digest(hashes, algo = "sha256", serialize = TRUE)
}
availability <- do.call(rbind, lapply(engines, function(engine) {
  available <- requireNamespace(engine, quietly = TRUE)
  data.frame(
    Engine = engine,
    Available = available,
    Version = if (available) as.character(utils::packageVersion(engine)) else NA_character_,
    PrimaryFunctions = paste(primary_functions[[engine]], collapse = "/"),
    FunctionSHA256 = function_hash(engine, primary_functions[[engine]]),
    NamespaceFunctionSHA256 = namespace_hash(engine),
    NativeLibrarySHA256 = native_library_hash(engine),
    SourceGitHead = if (engine == "mfrmr") args[["mfrmr-git-head"]] %||% NA_character_ else NA_character_,
    SourceState = if (engine == "mfrmr") args[["mfrmr-source-state"]] %||% "installed-library" else "installed-library",
    BundleInventorySHA256 = bundle_inventory_sha256,
    BundleValidated = all(validation$Passed),
    AdapterStatus = if (all(validation$Passed)) "input_ready_fit_adapter_pending" else "blocked_invalid_input",
    stringsAsFactors = FALSE
  )
}))
utils::write.csv(
  availability,
  file.path(output_dir, "bridge_engine_availability_r.csv"),
  row.names = FALSE
)

engine_lines <- vapply(seq_len(nrow(availability)), function(i) {
  paste0(
    "- ", availability$Engine[[i]], " ", availability$Version[[i]],
    ": ", availability$AdapterStatus[[i]], "."
  )
}, character(1))
bridge_results <- c(
  "# Cross-engine operating-characteristics input bridge",
  "",
  paste0("- Validation checks: ", sum(validation$Passed), "/", nrow(validation), " passed."),
  paste0("- Imported RunIds: ", nrow(identity), "."),
  paste0("- Imported rating rows: ", nrow(ratings), "."),
  paste0("- Paired clean/contaminated-anchor inputs: ", nrow(paired), "."),
  paste0("- Bundle inventory SHA-256: `", bundle_inventory_sha256, "`."),
  "",
  "## Detected adapter targets",
  "",
  engine_lines,
  "",
  "## Boundary",
  "",
  "The exact generated inputs are byte-validated and ready for adapters, but no R fit",
  "is executed by this bridge. Package availability is not cross-engine agreement,",
  "and this smoke bundle is not operating-characteristic performance evidence."
)
writeLines(bridge_results, file.path(output_dir, "BRIDGE_RESULTS.md"), useBytes = TRUE)

run_counts <- merge(
  identity[c("RunId", "ConditionId", "DataId", "FitInputId", "RatingsRows", "FacetTruthRows", "AnchorRows")],
  manifest[c("RunId", "Design", "TruthBias", "Replicate", "Seed")],
  by = "RunId",
  all.x = TRUE,
  sort = FALSE
)
utils::write.csv(
  run_counts,
  file.path(output_dir, "bridge_import_summary_r.csv"),
  row.names = FALSE
)

if (!all(validation$Passed)) {
  failed <- validation$Check[!validation$Passed]
  stop("Bridge validation failed: ", paste(failed, collapse = ", "), call. = FALSE)
}

message(
  "Validated ", nrow(identity), " RunIds and imported ", nrow(ratings),
  " rating rows for matched mfrmr/TAM/immer/sirt adapters."
)
