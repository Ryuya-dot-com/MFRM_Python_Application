#!/usr/bin/env Rscript
# TAM fits of the same fully crossed observed PCM stress fixture.
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
spec <- jsonlite::read_json(args[1], simplifyVector = TRUE)
data <- spec$data
stopifnot(nrow(data) == 480L, !spec$scientific_inference_ready, !spec$qualification_eligible)
grid <- unique(data[c("Person", "Rater")])
grid <- grid[order(grid$Person, grid$Rater), ]
key <- paste(data$Person, data$Rater, data$Criterion)
response <- sapply(paste0("C", 1:5), function(c)
  data$Score[match(paste(grid$Person, grid$Rater, c), key)])
stopifnot(!anyNA(response), identical(dim(response), c(96L, 5L)))
for (i in seq_len(nrow(spec$tam))) {
  config <- spec$tam[i, ]
  path <- file.path(args[2], paste0(config$id, ".json"))
  stopifnot(!file.exists(path))
  warnings <- character()
  fit <- withCallingHandlers(TAM::tam.mml.mfr(resp = response,
    facets = data.frame(rater = grid$Rater), pid = grid$Person,
    formulaA = ~ item + rater + item:step, constraint = "cases", est.variance = TRUE,
    control = list(nodes = seq(-config$bound, config$bound, length.out = config$q),
      snodes = 0L, maxiter = 2000L, conv = 1e-12, convD = 1e-12, convM = 1e-10,
      Msteps = 20L, progress = FALSE), verbose = FALSE),
    warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
  saveRDS(fit, file.path(args[2], paste0(config$id, ".rds")))
  result <- list(config = as.list(config), warnings = warnings,
    runtime = list(R = R.version.string, TAM = as.character(packageVersion("TAM"))),
    xsi = fit$xsi, xsi_facets = fit$xsi.facets, item = fit$item,
    A = fit$A, B = fit$B, AXsi = fit$AXsi, beta = fit$beta, variance = fit$variance,
    deviance = fit$deviance, iter = fit$iter, deviance_history = fit$deviance.history,
    control = fit$control, person = fit$person, pid = fit$pid)
  jsonlite::write_json(result, path, digits = NA, auto_unbox = TRUE, pretty = TRUE, null = "null")
  cat(config$id, "variance", fit$variance[1], "deviance", fit$deviance, "iterations", fit$iter, "\n")
}
