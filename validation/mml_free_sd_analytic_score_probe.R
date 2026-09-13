#!/usr/bin/env Rscript
# Follow-up on the retained 24-Person fixture; not a qualification experiment.
# Run from the repository root with an UNUSED output JSON path.
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L, !file.exists(args[1]))
previous <- "validation/mml_free_sd_cross_language_20260912"
script <- "validation/mml_free_sd_analytic_score_probe.R"
sha256 <- function(path) digest::digest(file = path, algo = "sha256")
manifest <- jsonlite::read_json(file.path(previous, "summary.json"))
for (path in names(manifest$source_sha256)) {
  stopifnot(identical(sha256(path), manifest$source_sha256[[path]]))
}
for (path in names(manifest$artifact_sha256)) {
  stopifnot(identical(sha256(file.path(previous, path)), manifest$artifact_sha256[[path]]))
}
input <- jsonlite::read_json(file.path(previous, "input.json"), simplifyVector = TRUE)
python <- jsonlite::read_json(file.path(previous, "python.json"))
prior_r <- jsonlite::read_json(file.path(previous, "r.json"))
stopifnot(identical(input$classification, "OBSERVED_DEVELOPMENT_ONLY"),
          !manifest$scientific_inference_ready, !manifest$qualification_eligible)

# Reuse only named definitions; never run either historical experiment again.
imports <- list(
  "validation/known_assignment_mml_crossfit.R" =
    c("marginal_components", "objective", "finite_difference_gradient"),
  "validation/mml_free_sd_cross_language_probe.R" =
    c("data", "prepared", "rules", "normal_quadrature", "expand_parameters", "probabilities"))
for (path in names(imports)) {
  for (e in parse(path)) {
    if (is.call(e) && identical(e[[1]], as.name("<-")) &&
        as.character(e[[2]])[1] %in% imports[[path]]) eval(e)
  }
}
stopifnot(all(vapply(unlist(imports), exists, logical(1), envir = .GlobalEnv, inherits = FALSE)))
names(rules) <- c("31", "61", "121", "181")
stopifnot(nrow(data) == 192L, all(prepared$score %in% 0:3),
          !anyNA(unlist(prepared)), !is.unsorted(prepared$person),
          all(tabulate(prepared$person) == 8L), length(unique(data$Person)) == 24L,
          nrow(input$cases) == 32L, nrow(input$fits) == 8L,
          all(vapply(rules, function(x) all(x$weights > 0), logical(1))))

# d eta / d (r1,t1,c1,c2), with r2=-r1, t2=-t1; Criterion uncentered.
eta_design <- cbind(ifelse(prepared$rater == 1L, -1, 1),
                    ifelse(prepared$task == 1L, -1, 1),
                    -as.numeric(prepared$criterion == 1L),
                    -as.numeric(prepared$criterion == 2L))
# d logit_k / d (tau1,tau2), since tau3=-tau1-tau2.
step_design <- rbind(c(0, 0), c(-1, 0), c(-1, -1), c(0, 0))
analytic_gradient <- function(par, prepared, q) {
  bundle <- marginal_components(par, prepared, q, retain_posterior = TRUE)
  gradient <- numeric(length(par))
  for (j in seq_len(q)) {
    theta <- bundle$nodes[j]
    probs <- probabilities(par, theta)
    residual <- prepared$score - as.numeric(probs %*% (0:3))
    conditional_steps <- step_design[prepared$score + 1L, ] - probs %*% step_design
    if (model == "PCM") {
      conditional_steps <- cbind(conditional_steps * (prepared$criterion == 1L),
                                  conditional_steps * (prepared$criterion == 2L))
    }
    # GH weights are constant; theta=exp(log_sigma)*z moves with log_sigma.
    conditional <- cbind(eta_design * residual, conditional_steps, theta * residual)
    gradient <- gradient - colSums(conditional * bundle$posterior[prepared$person, j])
  }
  stopifnot(all(is.finite(gradient)))
  gradient
}
score_pair <- function(par, q) {
  bundle <- marginal_components(par, prepared, q, retain_posterior = TRUE)
  z <- rules[[as.character(q)]]$nodes
  list(moving_node_nll_gradient = tail(analytic_gradient(par, prepared, q), 1L),
       density_moment_nll_score = -sum(bundle$posterior %*% (z^2 - 1)))
}
local_audit <- function(par, q) {
  gradient <- analytic_gradient(par, prepared, q)
  hessian <- function(step) {
    jacobian <- vapply(seq_along(par), function(i) {
      h <- step * max(1, abs(par[i]))
      plus <- minus <- par
      plus[i] <- plus[i] + h
      minus[i] <- minus[i] - h
      (analytic_gradient(plus, prepared, q) - analytic_gradient(minus, prepared, q)) / (2 * h)
    }, numeric(length(par)))
    jacobian
  }
  h1 <- hessian(1e-4)
  h2 <- hessian(5e-5)
  symmetric <- (h2 + t(h2)) / 2
  eigenvalues <- eigen(symmetric, symmetric = TRUE, only.values = TRUE)$values
  correction <- if (min(eigenvalues) > 0) solve(symmetric, gradient) else NULL
  list(gradient = gradient, gradient_supnorm = max(abs(gradient)),
       hessian_step_difference = max(abs(h1 - h2)),
       hessian_asymmetry = max(abs(h2 - t(h2))), eigenvalues = eigenvalues,
       newton_correction = correction,
       predicted_nll_gain = if (is.null(correction)) NULL else sum(gradient * correction) / 2,
       nll_times_machine_epsilon = abs(objective(par, prepared, q)) * .Machine$double.eps)
}

# Thresholds inherited from the earlier implementation check; NOT inference gates.
limits <- list(gradient = 1e-6, coordinate = 1e-4, nll = 1e-8, fitted_score = 1e-4)
cases <- lapply(seq_len(nrow(input$cases)), function(i) {
  case <- input$cases[i, ]
  model <<- case$model
  par <- as.numeric(case$coordinates[[1]])
  gradient <- analytic_gradient(par, prepared, case$q)
  py <- unlist(python$cases[[case$id]]$gradient)
  fd <- finite_difference_gradient(par, prepared, case$q, step = 1e-5)
  half <- finite_difference_gradient(par, prepared, case$q, step = 5e-6)
  differences <- c(python = max(abs(gradient - py)),
                   fd_h = max(abs(gradient - fd)), fd_half = max(abs(gradient - half)))
  list(id = case$id, gradient = gradient, differences = as.list(differences),
       score_pair = score_pair(par, case$q), pass = all(differences < limits$gradient))
})
stopifnot(all(vapply(cases, `[[`, logical(1), "pass")))

# Fix this experiment before observing refit outcomes: same eight starts,
# objective, bounds and controls; only replace finite differences with gr.
control <- list(maxit = 250L, factr = 1, pgtol = 1e-8)
refits <- lapply(seq_len(nrow(input$fits)), function(i) {
  case <- input$fits[i, ]
  model <<- case$model
  old <- prior_r$refits[[i]]
  stopifnot(identical(old$model, model), old$q == case$q)
  old_par <- unlist(old$coordinates)
  old_audit <- local_audit(old_par, case$q)
  ladder <- lapply(c(1e-3, 1e-4, 1e-5, 5e-6, 1e-6, 1e-7), function(h) {
    fd <- finite_difference_gradient(old_par, prepared, case$q, step = h)
    list(step = h, gradient = fd, difference = max(abs(fd - old_audit$gradient)))
  })
  start <- as.numeric(case$start[[1]])
  fit <- optim(start, objective, gr = analytic_gradient, prepared = prepared, q = case$q,
               method = "L-BFGS-B", lower = c(rep(-Inf, length(start) - 1L), log(0.05)),
               upper = c(rep(Inf, length(start) - 1L), log(10)), control = control)
  audit <- local_audit(fit$par, case$q)
  py_par <- as.numeric(case$python_coordinates[[1]])
  coordinate_difference <- max(abs(fit$par - py_par))
  nll_difference <- objective(fit$par, prepared, case$q) - objective(py_par, prepared, case$q)
  cat(model, "Q", case$q, "analytic refit:", fit$convergence,
      "score:", format(audit$gradient_supnorm), "\n")
  list(model = model, q = case$q, start = start, coordinates = fit$par,
       convergence = fit$convergence, message = fit$message, counts = as.list(fit$counts),
       nll = fit$value, audit = audit,
       python_coordinate_difference = coordinate_difference, python_nll_difference = nll_difference,
       numeric_checks_pass = coordinate_difference < limits$coordinate &&
         abs(nll_difference) < limits$nll && audit$gradient_supnorm < limits$fitted_score,
       prior_r_convergence = old$convergence, prior_r_message = old$message,
       prior_r_audit = old_audit, prior_r_difference_ladder = ladder,
       score_pair_at_python_fit = score_pair(py_par, case$q))
})
numeric_pass <- all(vapply(refits, `[[`, logical(1), "numeric_checks_pass"))
status_pass <- all(vapply(refits, function(x) x$convergence == 0L, logical(1)))
jsonlite::write_json(list(classification = "OBSERVED_DEVELOPMENT_ONLY",
  qualification_eligible = FALSE, scientific_inference_ready = FALSE,
  source_sha256 = c(manifest$source_sha256, setNames(list(sha256(script)), script)),
  prior_artifact_sha256 = c(manifest$artifact_sha256,
    list("summary.json" = sha256(file.path(previous, "summary.json")))),
  limits = limits, control = control, cases = cases, refits = refits,
  gradient_checks_pass = TRUE, refit_numeric_checks_pass = numeric_pass,
  refit_status_checks_pass = status_pass, implementation_checks_pass = numeric_pass && status_pass,
  runtime = list(R = R.version.string, platform = R.version$platform,
    statmod = as.character(packageVersion("statmod")), digest = as.character(packageVersion("digest")),
    jsonlite = as.character(packageVersion("jsonlite")), LAPACK = La_version(),
    BLAS = unname(extSoftVersion()["BLAS"])),
  limits_of_evidence = c("One existing RSM-generated dataset, also fitted as PCM",
    "No SE, CI, coverage, general success-rate or global-optimum qualification",
    "Earlier finite-difference failures remain recorded; no optimizer restart in this follow-up")),
  args[1], digits = NA, pretty = TRUE, auto_unbox = TRUE, na = "null", null = "null")
quit(status = if (numeric_pass && status_pass) 0L else 1L)
