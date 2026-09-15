#!/usr/bin/env Rscript
# Independent scalar weighted-PCM checks on new paired designs.
# Adapted from the preserved one-case evaluator; its source remains unchanged.
args <- commandArgs(TRUE)
stopifnot(length(args) == 3L)
input <- jsonlite::read_json(args[1], simplifyVector = TRUE)
cases <- jsonlite::read_json(args[2])
stopifnot(identical(input$classification, "OBSERVED_DEVELOPMENT_ONLY"), !input$scientific_inference_ready)
idx <- input$indices
person <- idx$person + 1L
rater <- idx$facets$Rater + 1L
task <- idx$facets$Task + 1L
criterion <- idx$facets$Criterion + 1L
stopifnot(length(idx$score_k) > 0L, length(idx$score_k) == length(idx$weight), all(idx$score_k %in% 0:3),
          all(is.finite(idx$weight)), all(idx$weight >= 0), sum(idx$weight) > 0,
          all(rater %in% 1:2), all(task %in% 1:2), all(criterion %in% 1:2))
rows_by_person <- lapply(seq_len(input$n_person), function(p) which(person == p))

continuous <- function(par, settings) {
  stopifnot(length(par) == 10L, all(is.finite(par)))
  sigma <- exp(par[10])
  raters <- c(0.35, par[2])
  tasks <- c(par[3], 0.4-par[3])
  criteria <- par[4:5]
  cumulative <- rbind(c(0, par[6], par[6]+par[7], 0), c(0, par[8], par[8]+par[9], 0))
  bound <- settings$bound
  values <- lapply(seq_along(rows_by_person), function(p) {
    rows <- rows_by_person[[p]]
    mu <- input$x[p]*par[1]
    log_integrand <- function(z) vapply(z, function(v) {
      if (!length(rows)) return(dnorm(v, log = TRUE))
      eta <- mu+sigma*v+raters[rater[rows]]-tasks[task[rows]]-criteria[criterion[rows]]
      logits <- outer(eta, 0:3)-cumulative[criterion[rows], , drop = FALSE]
      shifted <- logits-apply(logits, 1, max)
      logp <- shifted-log(rowSums(exp(shifted)))
      sum(idx$weight[rows]*logp[cbind(seq_along(rows), idx$score_k[rows]+1L)])+dnorm(v, log = TRUE)
    }, numeric(1))
    mode <- optimize(log_integrand, c(-bound, bound), maximum = TRUE, tol = 1e-10)
    stopifnot(abs(mode$maximum) < bound-0.1)
    cuts <- sort(unique(c(-bound, 0, mode$maximum, bound)))
    integrals <- lapply(0:2, function(k) {
      f <- function(z) exp(log_integrand(z)-mode$objective)*z^k
      parts <- lapply(seq_len(length(cuts)-1L), function(i) integrate(f, cuts[i], cuts[i+1L],
        subdivisions = 500L, rel.tol = settings$rel_tol, abs.tol = settings$abs_tol))
      stopifnot(all(vapply(parts, function(x) identical(x$message, "OK"), logical(1))))
      c(value = sum(vapply(parts, `[[`, numeric(1), "value")),
        error = sum(vapply(parts, `[[`, numeric(1), "abs.error")))
    })
    a <- vapply(integrals, `[[`, numeric(1), "value")
    ez <- a[2]/a[1]
    ez2 <- a[3]/a[1]
    stopifnot(a[1] > 0, ez2 > ez^2)
    list(log_mass = mode$objective+log(a[1]), eap = mu+sigma*ez,
         sd = sigma*sqrt(ez2-ez^2), ez2 = ez2,
         numeric_relative_mass_error = integrals[[1]]["error"]/a[1])
  })
  mass <- vapply(values, `[[`, numeric(1), "log_mass")
  eap <- vapply(values, `[[`, numeric(1), "eap")
  sd <- vapply(values, `[[`, numeric(1), "sd")
  empty <- which(vapply(rows_by_person, function(rows) sum(idx$weight[rows]), numeric(1)) == 0)
  if (length(empty)) stopifnot(max(abs(eap[empty]-input$x[empty]*par[1])) < 1e-9, max(abs(sd[empty]-sigma)) < 1e-9)
  list(nll = -sum(mass), eap = eap, sd = sd,
       log_sigma_nll_score = -sum(vapply(values, `[[`, numeric(1), "ez2")-1),
       numeric_relative_mass_error_sum = sum(vapply(values, `[[`, numeric(1), "numeric_relative_mass_error")),
       tail_relative_mass_bound_sum = sum(exp(log(2*pnorm(-bound))-mass)))
}
for (name in names(cases)) {
  output <- file.path(args[3], paste0(name, "_r.json"))
  stopifnot(!file.exists(output))
  result <- list(continuous = lapply(seq_len(nrow(input$r_integration)), function(i)
      continuous(unlist(cases[[name]]$coordinates), as.list(input$r_integration[i, ]))),
    runtime = list(R = R.version.string, jsonlite = as.character(packageVersion("jsonlite"))))
  jsonlite::write_json(result, output, pretty = TRUE, auto_unbox = TRUE, digits = NA)
  cat(name, "NLL", result$continuous[[2]]$nll, "\n")
}
