#!/usr/bin/env Rscript
# Trace the installed TAM updates without changing any estimation function.
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
spec <- jsonlite::read_json(args[1], simplifyVector = TRUE)
stopifnot(!spec$scientific_inference_ready, !spec$qualification_eligible)
data <- spec$data
grid <- unique(data[c("Person", "Rater")])
grid <- grid[order(grid$Person, grid$Rater), ]
key <- paste(data$Person, data$Rater, data$Criterion)
response <- sapply(paste0("C", 1:5), function(c)
  data$Score[match(paste(grid$Person, grid$Rater, c), key)])
stopifnot(nrow(data) == 480L, !anyNA(response), identical(dim(response), c(96L, 5L)))
namespace <- asNamespace("TAM")
original_progress <- get("tam_mml_progress_em", namespace)
trace("tam_mml_progress_em", where = namespace, print = FALSE, tracer = quote({
  .GlobalEnv$.tam_updates[[length(.GlobalEnv$.tam_updates) + 1L]] <- c(
    iter = iter, deviance = deviance, deviance_change = deviance_change,
    xsi_change = xsi_change, beta_change = beta_change, variance_change = variance_change)
}))
for (i in seq_len(nrow(spec$tam))) {
  config <- spec$tam[i, ]
  path <- file.path(args[2], paste0(config$id, ".json"))
  stopifnot(!file.exists(path))
  .tam_updates <- list()
  warnings <- character()
  fit <- withCallingHandlers(TAM::tam.mml.mfr(resp = response,
    facets = data.frame(rater = grid$Rater), pid = grid$Person,
    formulaA = ~ item + rater + item:step, constraint = "cases", est.variance = TRUE,
    control = list(nodes = seq(-config$bound, config$bound, length.out = config$q),
      snodes = 0L, maxiter = 2000L, conv = config$conv, convD = config$convD,
      convM = config$convM, Msteps = 20L, progress = FALSE), verbose = FALSE),
    warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
  updates <- as.data.frame(do.call(rbind, .tam_updates))
  stopifnot(nrow(updates) == fit$iter, identical(as.integer(updates$iter), seq_len(fit$iter)))
  last <- tail(updates, 1)
  # These are the installed tam.mml.mfr loop's sticky population flags and
  # terminal xsi/deviance criteria; unit slopes leave a4 at its initial zero.
  guard <- c(beta_ever_small = any(updates$beta_change < config$conv),
    variance_ever_small = any(updates$variance_change < config$conv),
    terminal_xsi_small = last$xsi_change <= config$conv,
    terminal_deviance_small = last$deviance_change <= config$convD)
  saveRDS(fit, file.path(args[2], paste0(config$id, ".rds")))
  result <- list(config = as.list(config), warnings = warnings,
    runtime = list(R = R.version.string, TAM = as.character(packageVersion("TAM"))),
    xsi = fit$xsi, xsi_facets = fit$xsi.facets, item = fit$item,
    A = fit$A, B = fit$B, AXsi = fit$AXsi, beta = fit$beta, variance = fit$variance,
    deviance = fit$deviance, iter = fit$iter, deviance_history = fit$deviance.history,
    control = fit$control, person = fit$person, pid = fit$pid,
    updates = updates, terminal_loop_checks = as.list(guard),
    loop_stopping_checks_pass = all(guard), reached_iteration_cap = fit$iter >= 2000L)
  jsonlite::write_json(result, path, digits = NA, auto_unbox = TRUE, pretty = TRUE, null = "null")
  cat(config$id, "SD", sqrt(fit$variance[1]), "iterations", fit$iter,
      "guard", guard, "last changes", unlist(last[3:6]), "\n")
}
untrace("tam_mml_progress_em", where = namespace)
stopifnot(identical(get("tam_mml_progress_em", namespace), original_progress))
