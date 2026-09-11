#!/usr/bin/env Rscript
# Bounded development check: 24 Persons, 2 Raters/Tasks/Criteria, scores 0:3.
# Reuse the independent R likelihood, leaving its historical runner unchanged.
args <- commandArgs(TRUE)
stopifnot(length(args) == 3L, !file.exists(args[3]))
input <- jsonlite::read_json(args[1], simplifyVector = TRUE)
stopifnot(identical(input$classification, "OBSERVED_DEVELOPMENT_ONLY"))
helpers <- c("marginal_components", "objective", "finite_difference_gradient")
for (e in parse(args[2])) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) &&
      as.character(e[[2]])[1] %in% helpers) eval(e)
}
stopifnot(all(vapply(helpers, exists, logical(1), envir = .GlobalEnv, inherits = FALSE)))
data <- input$data
stopifnot(nrow(data) == 192L, all(data$Score %in% 0:3),
          identical(sort(unique(data$Person)), sprintf("P%03d", 0:23)),
          all(data$Rater %in% c("R1", "R2")),
          all(data$Task %in% c("T1", "T2")),
          all(data$Criterion %in% c("C1", "C2")))
prepared <- list(person = match(data$Person, sort(unique(data$Person))),
                 rater = match(data$Rater, c("R1", "R2")),
                 task = match(data$Task, c("T1", "T2")),
                 criterion = match(data$Criterion, c("C1", "C2")),
                 score = data$Score, n_person = 24L)
stopifnot(!is.unsorted(prepared$person), all(tabulate(prepared$person) == 8L))

# Existing statmod uses a different eigensolver and retains positive tail weights.
rules <- lapply(c(31L, 61L, 121L, 181L), statmod::gauss.quad.prob, dist = "normal")
names(rules) <- c("31", "61", "121", "181")
stopifnot(all(vapply(rules, function(x) {
  all(is.finite(x$nodes)) && all(is.finite(x$weights)) &&
    all(x$weights > 0) && abs(sum(x$weights) - 1) < 1e-12
}, logical(1))))
normal_quadrature <- function(n, sigma) {
  stopifnot(as.character(n) %in% names(rules), is.finite(sigma), sigma > 0)
  x <- rules[[as.character(n)]]
  list(nodes = x$nodes * sigma, weights = x$weights)
}
expand_parameters <- function(par) {
  stopifnot(model %in% c("RSM", "PCM"), all(is.finite(par)),
            length(par) == if (model == "RSM") 7L else 9L)
  steps <- if (model == "RSM") {
    matrix(rep(c(par[5:6], -sum(par[5:6])), 2L), nrow = 2L, byrow = TRUE)
  } else rbind(c(par[5:6], -sum(par[5:6])), c(par[7:8], -sum(par[7:8])))
  list(Rater = c(par[1], -par[1]), Task = c(par[2], -par[2]),
       Criterion = par[3:4], steps = steps, sigma = exp(tail(par, 1L)))
}
probabilities <- function(par, theta) {
  x <- expand_parameters(par)
  base <- -x$Rater[prepared$rater] - x$Task[prepared$task] - x$Criterion[prepared$criterion]
  steps <- t(apply(x$steps, 1L, function(v) c(0, cumsum(v))))
  logits <- outer(theta + base, 0:3) - steps[prepared$criterion, , drop = FALSE]
  shifted <- logits - apply(logits, 1L, max)
  exp(shifted) / rowSums(exp(shifted))
}
score_summary <- function(par, q) {
  c <- marginal_components(par, prepared, q, retain_posterior = TRUE)
  eap <- as.numeric(c$posterior %*% c$nodes)
  centered <- matrix(c$nodes, 24L, q, byrow = TRUE) - eap
  list(nll = -c$loglik, eap = eap,
       sd = sqrt(rowSums(c$posterior * centered^2)))
}

# Independent adaptive integral in standard-normal z, split at the pattern mode.
# The finite [-12,12] interval has explicit absolute normal-tail moment bounds.
continuous_summary <- function(par) {
  sigma <- expand_parameters(par)$sigma
  per_person <- lapply(seq_len(24L), function(p) {
    rows <- which(prepared$person == p)
    log_integrand <- function(z) vapply(z, function(v) {
      probs <- probabilities(par, sigma * v)
      sum(log(probs[cbind(rows, prepared$score[rows] + 1L)])) + dnorm(v, log = TRUE)
    }, numeric(1))
    mode <- optimize(log_integrand, c(-12, 12), maximum = TRUE, tol = 1e-10)
    stopifnot(abs(mode$maximum) < 11.9)
    integral <- lapply(0:2, function(k) {
      f <- function(z) exp(log_integrand(z) - mode$objective) * z^k
      parts <- lapply(list(c(-12, mode$maximum), c(mode$maximum, 12)), function(b) {
        integrate(f, b[1], b[2], subdivisions = 300L, rel.tol = 1e-10, abs.tol = 1e-12)
      })
      c(value = sum(vapply(parts, `[[`, numeric(1), "value")),
        error = sum(vapply(parts, `[[`, numeric(1), "abs.error")))
    })
    mass <- integral[[1]][["value"]]
    mean <- sigma * integral[[2]][["value"]] / mass
    variance <- sigma^2 * integral[[3]][["value"]] / mass - mean^2
    stopifnot(mass > 0, variance > 0)
    list(log_marginal = mode$objective + log(mass), eap = mean, sd = sqrt(variance),
         mode_z = mode$maximum,
         scaled_integrals = vapply(integral, `[[`, numeric(1), "value"),
         scaled_numeric_errors = vapply(integral, `[[`, numeric(1), "error"),
         log_scale = mode$objective)
  })
  list(nll = -sum(vapply(per_person, `[[`, numeric(1), "log_marginal")),
       eap = vapply(per_person, `[[`, numeric(1), "eap"),
       sd = vapply(per_person, `[[`, numeric(1), "sd"), per_person = per_person,
       absolute_tail_bounds = c(mass = 2 * pnorm(-12),
         first_absolute_theta_moment = sigma * 2 * dnorm(12),
         second_theta_moment = sigma^2 * 2 * (12 * dnorm(12) + pnorm(-12))))
}
evaluate <- function(par, q) {
  c(score_summary(par, q), list(
    gradient_h = finite_difference_gradient(par, prepared, q, step = 1e-5),
    gradient_half = finite_difference_gradient(par, prepared, q, step = 5e-6)))
}
# Analytic limiting control: zero structural effects and nearly zero latent SD
# give probability 1/4 per response and NLL 192*log(4).
for (model in c("RSM", "PCM")) {
  point <- c(rep(0, if (model == "RSM") 6L else 8L), log(1e-10))
  stopifnot(abs(objective(point, prepared, 31L) - 192 * log(4)) < 1e-10,
            max(abs(probabilities(point, 0) - 0.25)) < 1e-14)
}
cases <- lapply(seq_len(nrow(input$cases)), function(i) {
  case <- input$cases[i, ]
  model <<- case$model
  par <- as.numeric(case$coordinates[[1]])
  c(list(id = case$id, model = model, q = case$q), evaluate(par, case$q),
    list(probabilities = if (case$q == 31L) lapply(c(-3, 0, 3),
      function(theta) probabilities(par, theta)[seq_len(8L), ]) else NULL))
})
refits <- lapply(seq_len(nrow(input$fits)), function(i) {
  case <- input$fits[i, ]
  model <<- case$model
  start <- as.numeric(case$start[[1]])
  fit <- optim(start, objective, gr = finite_difference_gradient,
               prepared = prepared, q = case$q, method = "L-BFGS-B",
               lower = c(rep(-Inf, length(start) - 1L), log(0.05)),
               upper = c(rep(Inf, length(start) - 1L), log(10)),
               control = list(maxit = 250L, factr = 1, pgtol = 1e-8))
  cat(model, "Q", case$q, "R refit:", fit$convergence, "\n")
  # Post-result diagnostic after the initial development run had two code-52
  # line-search terminations. Apply once to ALL eight cases; retain primary codes.
  restart <- optim(fit$par, objective, gr = finite_difference_gradient,
                   prepared = prepared, q = case$q, method = "BFGS",
                   control = list(maxit = 100L, reltol = 1e-13))
  # Continuous comparison holds the Python point fixed, separately from this refit.
  python_par <- as.numeric(case$python_coordinates[[1]])
  list(model = model, q = case$q, coordinates = fit$par, start = start,
       convergence = fit$convergence, message = fit$message, counts = fit$counts,
       at_r_fit = evaluate(fit$par, case$q),
       diagnostic_restart = list(convergence = restart$convergence,
         coordinates = restart$par, counts = restart$counts,
         displacement = max(abs(restart$par - fit$par)),
         nll_gain = fit$value - restart$value, at_point = evaluate(restart$par, case$q)),
       continuous_at_python_fit = continuous_summary(python_par))
})
jsonlite::write_json(list(classification = "OBSERVED_DEVELOPMENT_ONLY",
  scientific_inference_ready = FALSE, rules = rules, cases = cases, refits = refits,
  runtime = list(R = R.version.string, platform = R.version$platform,
    statmod = as.character(packageVersion("statmod")),
    jsonlite = as.character(packageVersion("jsonlite")), LAPACK = La_version(),
    BLAS = unname(extSoftVersion()["BLAS"]))), args[3], digits = NA,
  pretty = TRUE, auto_unbox = TRUE, na = "null", null = "null")
