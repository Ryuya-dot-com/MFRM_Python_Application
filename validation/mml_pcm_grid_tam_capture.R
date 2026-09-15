#!/usr/bin/env Rscript
# Matched-start native TAM and independent fixed-grid/continuous PCM arithmetic.
args <- commandArgs(TRUE)
stopifnot(length(args) == 3L)
spec <- jsonlite::read_json(args[1], simplifyVector = TRUE)
out <- args[2]
mode <- args[3]
stopifnot(!spec$scientific_inference_ready, !spec$qualification_eligible,
          mode %in% c("native", "evaluate"))
write_new <- function(value, path) {
  stopifnot(!file.exists(path))
  jsonlite::write_json(value, path, digits = NA, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null")
}
run_native <- function() {
  stopifnot(as.character(packageVersion("TAM")) == spec$tam$version)
  namespace <- asNamespace("TAM")
  names <- c("tam_mml_progress_em0", "tam_mml_progress_em")
  originals <- lapply(names, get, envir = namespace)
  source_path <- file.path(out, "tam_capture_installed_source.txt")
  stopifnot(!file.exists(source_path))
  writeLines(c(R.version.string, as.character(packageVersion("TAM")),
    deparse(TAM::tam.mml.mfr), deparse(get("tam_mml_mstep_regression", namespace))), source_path)
  on.exit({for (name in names) untrace(name, where = namespace)}, add = TRUE)
  trace(names[1], where = namespace, print = FALSE, tracer = quote({
    frames <- Filter(function(f) exists("tam_fct", f, inherits = FALSE) &&
      identical(f$tam_fct, "tam.mml.mfr"), sys.frames())
    stopifnot(length(frames) == 1L)
    e <- frames[[1]]
    .GlobalEnv$.tam_before[[iter]] <- list(iter = iter, xsi = e$xsi, variance = e$variance, beta = e$beta)
    if (iter == 1L) .GlobalEnv$.tam_A <- e$A
  }))
  trace(names[2], where = namespace, print = FALSE, tracer = quote({
    frames <- Filter(function(f) exists("tam_fct", f, inherits = FALSE) &&
      identical(f$tam_fct, "tam.mml.mfr"), sys.frames())
    stopifnot(length(frames) == 1L)
    e <- frames[[1]]
    .GlobalEnv$.tam_after[[iter]] <- list(iter = iter, xsi = e$xsi, variance = e$variance, beta = e$beta,
      deviance = deviance, xsi_change = xsi_change, beta_change = beta_change,
      variance_change = variance_change, deviance_change = deviance_change)
  }))
  for (dataset in spec$datasets) {
    folder <- file.path(out, dataset)
    original <- jsonlite::read_json(file.path(folder, "input.json"), simplifyVector = TRUE)
    d <- original$data
    grid <- unique(d[c("Person", "Rater")]); grid <- grid[order(grid$Person, grid$Rater), ]
    key <- paste(d$Person, d$Rater, d$Criterion)
    stopifnot(nrow(d) == 1600L, !anyDuplicated(key))
    response <- sapply(paste0("C", 1:5), function(c)
      d$Score[match(paste(grid$Person, grid$Rater, c), key)])
    stopifnot(identical(dim(response), c(320L,5L)), !anyNA(response))
    for (bound in c(12,20)) {
      .GlobalEnv$.tam_before <- list(); .GlobalEnv$.tam_after <- list()
      warnings <- character()
      fit <- withCallingHandlers(TAM::tam.mml.mfr(resp = response,
        facets = data.frame(rater = grid$Rater), pid = grid$Person,
        formulaA = ~ item + rater + item:step, constraint = "cases", est.variance = TRUE,
        xsi.inits = cbind(1:23, rep(0,23)), variance.inits = matrix(1,1,1),
        control = list(nodes = seq(-bound,bound,length.out = bound*20+1), snodes = 0L,
          maxiter = spec$tam$maxiter, conv = spec$tam$conv, convD = spec$tam$convD,
          convM = spec$tam$convM, Msteps = spec$tam$Msteps, progress = FALSE), verbose = FALSE),
        warning = function(w) {warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")})
      stopifnot(identical(rownames(fit$xsi),dimnames(fit$A)[[3]]),
        length(.tam_before) == fit$iter, length(.tam_after) == fit$iter,
        length(.tam_before[[1]]$xsi) == 23L, all(.tam_before[[1]]$xsi == 0),
        .tam_before[[1]]$variance == 1, .tam_before[[1]]$beta == 0, identical(.tam_A, fit$A))
      path <- file.path(folder, paste0("tam_b",bound))
      stopifnot(!file.exists(paste0(path,".rds")))
      saveRDS(fit, paste0(path,".rds"))
      write_new(list(bound = bound, spacing = .1, warnings = warnings,
        version = as.character(packageVersion("TAM")), initial = .tam_before[[1]],
        before = .tam_before, after = .tam_after, iter = fit$iter, A = fit$A, B = fit$B,
        AXsi = fit$AXsi, xsi = fit$xsi, returned_xsi = fit$xsi$xsi,
        variance = fit$variance, beta = fit$beta, person = fit$person,
        item = fit$item, deviance = fit$deviance, control = fit$control,
        reached_iteration_cap = fit$iter >= spec$tam$maxiter), paste0(path,".json"))
      cat(dataset, "TAM", bound, "iterations",fit$iter,"SD",sqrt(fit$variance[1]),"\n")
    }
  }
  for (name in names) untrace(name, where = namespace)
  stopifnot(all(vapply(seq_along(names), function(i) identical(get(names[i],namespace), originals[[i]]), logical(1))))
  on.exit(NULL)
}

stopifnot(mode == "native")
run_native()
