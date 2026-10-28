#!/usr/bin/env Rscript
# Development control: distinguish TAM's coordinate curvature from joint covariance.
# This deterministic two-item example is excluded from qualification studies.
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L, !dir.exists(args[1]))
out <- args[1]
dir.create(out, recursive = TRUE)
stopifnot(dir.exists(out))
contract <- list(classification = "DEVELOPMENT_ONLY", scientific_inference_ready = FALSE,
  counts = c(20L, 10L, 20L, 50L), patterns = rbind(c(0,0), c(0,1), c(1,0), c(1,1)),
  nodes = seq(-8, 8, by = .1), fixed_population_mean = 0, fixed_population_variance = 1,
  est_variance = FALSE,
  maxiter = 1000L, conv = 1e-10, convD = 1e-10, convM = 1e-10,
  finite_difference_step = 1e-4, tam_se_step = 1e-3,
  tolerances = list(nll = 1e-8, hessian = 1e-5, diagonal_se = 1e-6,
                    stationary_score = 1e-5, minimum_joint_se_gap = 1e-3))
jsonlite::write_json(contract, file.path(out, "contract.json"), digits = NA,
  pretty = TRUE, auto_unbox = TRUE)
file.copy("validation/tam_se_definition_check.R", file.path(out, "source.R"))
jsonlite::write_json(list(R = R.version.string, TAM = as.character(packageVersion("TAM"))),
  file.path(out, "runtime.json"), pretty = TRUE, auto_unbox = TRUE)
source_names <- c("tam.se", "tam_mml_se_quick", "tam_mml_se_quick_difference_quotient")
sources <- lapply(source_names, function(n) deparse(get(n, asNamespace("TAM"))))
names(sources) <- source_names
jsonlite::write_json(sources, file.path(out, "installed_tam_sources.json"), pretty = TRUE)
resp <- contract$patterns[rep(seq_len(4), contract$counts), , drop = FALSE]
colnames(resp) <- c("Item1", "Item2")
write.csv(resp, file.path(out, "responses.csv"), row.names = FALSE)
fit <- TAM::tam.mml(resp = resp, constraint = "cases",
  est.variance = contract$est_variance, variance.inits = matrix(1, 1, 1),
  control = list(nodes = contract$nodes, snodes = 0L, maxiter = contract$maxiter,
    conv = contract$conv, convD = contract$convD, convM = contract$convM,
    progress = FALSE), verbose = FALSE)
saveRDS(fit, file.path(out, "fit.rds"))
stopifnot(fit$beta[1] == 0, fit$variance[1] == 1,
  identical(as.numeric(fit$theta), contract$nodes), all(fit$pweights == 1),
  length(fit$xsi$xsi) == 2L, all(fit$B[, 1, 1] == 0), all(fit$B[, 2, 1] == 1),
  max(abs(fit$AXsi[, 2] + fit$xsi$xsi)) < 1e-12,
  all(fit$AXsi[, 1] == 0))

# Explicit binary probabilities. No TAM probability, score or Hessian helper.
evaluate <- function(b) {
  logits <- outer(contract$nodes, b, "-")
  p <- plogis(logits)
  log_prior <- dnorm(contract$nodes, log = TRUE) + log(.1)
  value <- 0
  score <- numeric(2)
  information <- matrix(0, 2, 2)
  for (r in seq_len(nrow(contract$patterns))) {
    y <- contract$patterns[r, ]
    count <- contract$counts[r]
    lp <- log_prior + rowSums(sweep(log(p), 2, y, "*") +
      sweep(log1p(-p), 2, 1-y, "*"))
    centered <- exp(lp - max(lp))
    posterior <- centered / sum(centered)
    conditional_score <- sweep(-p, 2, y, "+")
    mean_score <- colSums(posterior * conditional_score)
    residual <- sweep(conditional_score, 2, mean_score, "-")
    complete <- diag(colSums(posterior * p * (1-p)))
    missing <- crossprod(residual, posterior * residual)
    value <- value - count * (max(lp) + log(sum(centered)))
    score <- score + count * mean_score
    information <- information + count * (complete - missing)
  }
  list(nll = value, score = score, information = information)
}
point <- fit$xsi$xsi
oracle <- evaluate(point)
h <- contract$finite_difference_step
fd_information <- sapply(seq_along(point), function(j) {
  d <- numeric(length(point)); d[j] <- h
  (evaluate(point + d)$score - evaluate(point - d)$score) / (2*h)
})
tam_se <- TAM::tam.se(fit, numdiff.parm = contract$tam_se_step)
info <- oracle$information
stopifnot(min(eigen(info, symmetric = TRUE)$values) > 0)
joint_covariance <- solve(info)
diagonal_se <- 1 / sqrt(diag(info))
joint_se <- sqrt(diag(joint_covariance))
differences <- list(nll = abs(oracle$nll - fit$deviance/2),
  hessian = max(abs(info - fd_information)),
  diagonal_se = max(abs(tam_se$xsi$se - diagonal_se)),
  stationary_score = max(abs(oracle$score)),
  joint_se_gap = max(abs(tam_se$xsi$se - joint_se)))
t <- contract$tolerances
checks <- list(nll = differences$nll < t$nll,
  hessian = differences$hessian < t$hessian,
  diagonal_se = differences$diagonal_se < t$diagonal_se,
  stationary_score = differences$stationary_score < t$stationary_score,
  joint_se_distinct = differences$joint_se_gap > t$minimum_joint_se_gap,
  terminated_below_cap = fit$iter < contract$maxiter)
summary <- list(classification = contract$classification, scientific_inference_ready = FALSE,
  checks = checks, all_checks_pass = all(unlist(checks)), differences = differences,
  point = point, information = info, joint_covariance = joint_covariance,
  diagonal_se = diagonal_se, joint_se = joint_se, tam_se = tam_se$xsi,
  native_item_se = fit$xsi, iterations = fit$iter,
  runtime = list(R = R.version.string, TAM = as.character(packageVersion("TAM"))),
  source_sha256 = digest::digest(file = "validation/tam_se_definition_check.R", algo = "sha256"),
  interpretation = paste("In this fixed-population control, tam.se equals inverse square-root",
    "coordinate-wise observed curvature. It is not sqrt(diag(inverse(full information))).",
    "No population-SD uncertainty or coverage qualification is asserted."))
jsonlite::write_json(summary, file.path(out, "summary.json"), digits = NA,
  pretty = TRUE, auto_unbox = TRUE)
files <- setdiff(list.files(out, full.names = TRUE), file.path(out, "manifest.json"))
hashes <- setNames(lapply(files, function(p) digest::digest(file = p, algo = "sha256")), basename(files))
jsonlite::write_json(hashes, file.path(out, "manifest.json"), pretty = TRUE, auto_unbox = TRUE)
print(data.frame(item = rownames(fit$xsi), tam_se = tam_se$xsi$se,
                 diagonal_se = diagonal_se, joint_se = joint_se))
print(checks)
if (!all(unlist(checks))) quit(status = 1L)
