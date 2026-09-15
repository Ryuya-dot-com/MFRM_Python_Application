#!/usr/bin/env Rscript
# Independent scalar PCM evaluation for the fixed 24 x 4 x 5, score 0:4 fixture.
# Adapt the earlier mode-split R integration to this design, preserving old files.
args <- commandArgs(TRUE)
stopifnot(length(args) == 3L)
input <- jsonlite::read_json(args[1], simplifyVector = TRUE)
cases <- jsonlite::read_json(args[2])
stopifnot(identical(input$classification, "OBSERVED_DEVELOPMENT_ONLY"),
          !input$scientific_inference_ready, !input$qualification_eligible)
data <- input$data
stopifnot(nrow(data) == 480L, all(data$Score %in% 0:4),
          identical(sort(unique(data$Person)), sprintf("P%02d", 1:24)),
          all(data$Rater %in% paste0("R", 1:4)), all(data$Criterion %in% paste0("C", 1:5)),
          all(data$Score[data$Person == "P01"] == 0),
          all(data$Score[data$Person == "P02"] == 4))
person <- match(data$Person, sprintf("P%02d", 1:24))
rater <- match(data$Rater, paste0("R", 1:4))
criterion <- match(data$Criterion, paste0("C", 1:5))
rows_by_person <- lapply(1:24, function(p) which(person == p))
stopifnot(all(vapply(rows_by_person, function(r) length(r) == 20L &&
  nrow(unique(data[r, c("Rater", "Criterion")])) == 20L, logical(1))))
rules <- lapply(input$orders, statmod::gauss.quad.prob, dist = "normal")
names(rules) <- as.character(input$orders)
stopifnot(all(vapply(rules, function(r) all(is.finite(r$nodes)) && all(r$weights > 0) &&
  abs(sum(r$weights) - 1) < 1e-12 && abs(sum(r$weights * r$nodes^2) - 1) < 1e-12,
  logical(1))))

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
finite_summary <- function(x, rule) {
  theta <- x$sigma * rule$nodes
  log_likelihood <- vapply(theta, function(t) {
    p <- log_probabilities(x, t)
    selected <- p[cbind(seq_len(nrow(data)), data$Score + 1L)]
    vapply(rows_by_person, function(rows) sum(selected[rows]), numeric(1))
  }, numeric(24))
  joint <- sweep(log_likelihood, 2, log(rule$weights), "+")
  center <- apply(joint, 1, max)
  mass <- rowSums(exp(joint - center))
  posterior <- exp(joint - center) / mass
  eap <- as.numeric(posterior %*% theta)
  deviations <- matrix(theta, 24, length(theta), byrow = TRUE) - eap
  list(nll = -sum(center + log(mass)), eap = eap,
       sd = sqrt(rowSums(posterior * deviations^2)))
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

# Runnable analytic limiting check: every category has probability 1/5 as SD -> 0.
zero <- expand_parameters(c(rep(0, 23), log(1e-10)))
stopifnot(max(abs(exp(log_probabilities(zero, 0)) - 0.2)) < 1e-14)
for (settings in input$integration |> split(seq_len(nrow(input$integration)))) {
  control <- continuous_summary(zero, as.list(settings))
  stopifnot(abs(control$nll - 480 * log(5)) < 1e-8,
            max(abs(control$sd - 1e-10)) < 1e-18,
            abs(control$log_sigma_nll_score) < 1e-8)
}
for (name in names(cases)) {
  output <- file.path(args[3], paste0(name, "_r.json"))
  stopifnot(!file.exists(output))
  x <- expand_parameters(unlist(cases[[name]]$coordinates))
  result <- list(id = name, finite = lapply(rules, function(r) finite_summary(x, r)),
    probabilities = lapply(c(-3, 0, 3), function(t) exp(log_probabilities(x, t))[1:20, ]),
    continuous = lapply(seq_len(nrow(input$integration)), function(i)
      continuous_summary(x, as.list(input$integration[i, ]))),
    runtime = list(R = R.version.string, statmod = as.character(packageVersion("statmod")),
                   jsonlite = as.character(packageVersion("jsonlite"))))
  jsonlite::write_json(result, output, digits = NA, pretty = TRUE, auto_unbox = TRUE,
                       na = "null", null = "null")
  cat(name, "continuous NLL", result$continuous[[2]]$nll, "\n")
}
