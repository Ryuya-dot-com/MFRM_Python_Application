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
expand_parameters <- function(par) {
  stopifnot(length(par) == 24L, all(is.finite(par)))
  free_steps <- matrix(par[9:23], nrow = 5, byrow = TRUE)
  steps <- cbind(free_steps, -rowSums(free_steps))
  list(rater = c(par[1:3], -sum(par[1:3])), criterion = par[4:8],
       cumulative = t(apply(steps, 1, function(s) c(0, cumsum(s)))), sigma = exp(par[24]))
}
log_probabilities <- function(x, theta, rows = seq_len(nrow(data))) {
  base <- theta - x$rater[rater[rows]] - x$criterion[criterion[rows]]
  logits <- outer(base, 0:4) - x$cumulative[criterion[rows], , drop = FALSE]
  shifted <- logits - apply(logits, 1, max)
  shifted - log(rowSums(exp(shifted)))
}
continuous_summary <- function(x, settings) {
  bound <- settings$bound
  per_person <- lapply(rows_by_person, function(rows) {
    log_integrand <- function(z) vapply(z, function(v) {
      p <- log_probabilities(x, x$sigma * v, rows)
      sum(p[cbind(seq_along(rows), data$Score[rows] + 1L)]) + dnorm(v, log = TRUE)
    }, numeric(1))
    # Log PCM likelihood is concave in theta; adding log phi(z) makes it strictly concave.
    mode <- optimize(log_integrand, c(-bound, bound), maximum = TRUE, tol = 1e-10)
    stopifnot(abs(mode$maximum) < bound - 0.1)
    cuts <- sort(unique(c(-bound, 0, mode$maximum, bound)))
    integrals <- lapply(0:2, function(k) {
      f <- function(z) exp(log_integrand(z) - mode$objective) * z^k
      parts <- lapply(seq_len(length(cuts) - 1L), function(i) integrate(
        f, cuts[i], cuts[i + 1L], subdivisions = 500L,
        rel.tol = settings$rel_tol, abs.tol = settings$abs_tol))
      stopifnot(all(vapply(parts, function(p) identical(p$message, "OK"), logical(1))))
      c(value = sum(vapply(parts, `[[`, numeric(1), "value")),
        error = sum(vapply(parts, `[[`, numeric(1), "abs.error")))
    })
    values <- vapply(integrals, `[[`, numeric(1), "value")
    errors <- vapply(integrals, `[[`, numeric(1), "error")
    ez <- values[2] / values[1]
    ez2 <- values[3] / values[1]
    stopifnot(values[1] > 0, ez2 > ez^2)
    list(log_marginal = mode$objective + log(values[1]), eap = x$sigma * ez,
         sd = x$sigma * sqrt(ez2 - ez^2), ez2 = ez2, mode_z = mode$maximum,
         scaled_integrals = values, scaled_numeric_errors = errors, log_scale = mode$objective)
  })
  log_mass <- vapply(per_person, `[[`, numeric(1), "log_marginal")
  # Each unweighted response likelihood is <= 1, giving absolute normal-tail bounds.
  tails <- c(mass = 2 * pnorm(-bound), first_absolute_theta_moment = x$sigma * 2 * dnorm(bound),
             second_theta_moment = x$sigma^2 * 2 * (bound * dnorm(bound) + pnorm(-bound)))
  list(nll = -sum(log_mass), eap = vapply(per_person, `[[`, numeric(1), "eap"),
       sd = vapply(per_person, `[[`, numeric(1), "sd"),
       log_sigma_nll_score = -sum(vapply(per_person, `[[`, numeric(1), "ez2") - 1),
       settings = settings, per_person = per_person, absolute_tail_bounds = tails,
       numeric_relative_mass_error_sum = sum(vapply(per_person, function(p)
         p$scaled_numeric_errors[1] / p$scaled_integrals[1], numeric(1))),
       tail_relative_mass_bound_sum = sum(exp(log(tails[1]) - log_mass)))
}

fixed_summary <- function(x, config) {
  theta <- seq(-config$bound, config$bound, length.out = config$q)
  logw <- dnorm(theta, sd = x$sigma, log = TRUE) + log(config$spacing)
  ll <- vapply(theta, function(t) {
    prob <- log_probabilities(x, t)
    chosen <- prob[cbind(seq_len(nrow(data)), data$Score + 1L)]
    vapply(rows_by_person, function(rows) sum(chosen[rows]), numeric(1))
  }, numeric(n))
  joint <- sweep(ll, 2, logw, "+")
  center <- apply(joint, 1, max)
  mass <- rowSums(exp(joint - center))
  post <- exp(joint - center) / mass
  eap <- as.numeric(post %*% theta)
  sd <- sqrt(rowSums(post * (matrix(theta, n, length(theta), byrow = TRUE) - eap)^2))
  nll <- -sum(center + log(mass))
  prior_mass <- sum(exp(logw))
  list(nll = nll, eap = eap, sd = sd, raw_prior_mass = prior_mass,
       normalized_prior_nll = nll + n * log(prior_mass))
}

run_native <- function() {
  stopifnot(as.character(packageVersion("TAM")) == spec$tam$version)
  namespace <- asNamespace("TAM")
  names <- c("tam_mml_progress_em0", "tam_mml_progress_em")
  originals <- lapply(names, get, envir = namespace)
  source_path <- file.path(out, "tam_installed_source.txt")
  stopifnot(!file.exists(source_path))
  writeLines(c(R.version.string, as.character(packageVersion("TAM")),
    deparse(TAM::tam.mml.mfr), deparse(get("tam_mml_mstep_regression", namespace))), source_path)
  on.exit({for (name in names) untrace(name, where = namespace)}, add = TRUE)
  trace(names[1], where = namespace, print = FALSE, tracer = quote({
    e <- parent.frame()
    .GlobalEnv$.tam_before[[iter]] <- list(iter = iter, xsi = e$xsi, variance = e$variance, beta = e$beta)
    if (iter == 1L) .GlobalEnv$.tam_A <- e$A
  }))
  trace(names[2], where = namespace, print = FALSE, tracer = quote({
    e <- parent.frame()
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

if (mode == "native") {
  run_native()
} else {
  # One dataset per invocation, allowing independent R checks to run concurrently.
  data <- spec$data
  n <- length(unique(data$Person))
  stopifnot(n == 80L, nrow(data) == 20L*n, all(data$Score %in% 0:4))
  person <- match(data$Person, sprintf("P%02d",seq_len(n)))
  rater <- match(data$Rater,paste0("R",1:4)); criterion <- match(data$Criterion,paste0("C",1:5))
  rows_by_person <- lapply(seq_len(n),function(p) which(person == p))
  stopifnot(all(vapply(rows_by_person,function(r) length(r)==20L &&
    nrow(unique(data[r,c("Rater","Criterion")]))==20L,logical(1))))
  cases <- jsonlite::read_json(file.path(out,"r_cases.json"))
  for (name in names(cases)) {
    case <- cases[[name]]
    x <- expand_parameters(unlist(case$coordinates))
    finite <- lapply(case$grids, function(g) fixed_summary(x,g))
    names(finite) <- vapply(case$grids,function(g) g$id,character(1))
    continuous <- if (case$continuous) continuous_summary(x, as.list(spec$integration[2,])) else NULL
    write_new(list(finite = finite, continuous = continuous),file.path(out,paste0(name,"_r.json")))
    cat(basename(out),name,"R checked\n")
  }
}
