#!/usr/bin/env Rscript
# Independent probability formula and Richardson derivatives; shared finite GH rule.
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L, !file.exists(args[2]))
input <- jsonlite::read_json(args[1], simplifyVector = FALSE)
stopifnot(identical(input$classification, "DEVELOPMENT_ONLY"))
results <- lapply(input$cases, function(case) {
  cfg <- case$config
  facets <- unlist(cfg$facet_names)
  levels <- lapply(cfg$facet_levels, unlist)
  data <- do.call(rbind, lapply(case$data, as.data.frame, stringsAsFactors = FALSE))
  person <- match(data$Person, sort(unique(data$Person)))
  n <- nrow(data); q <- length(case$nodes); k <- cfg$n_cat
  nodes <- unlist(case$nodes); weights <- unlist(case$weights)
  stopifnot(cfg$model %in% c("RSM", "PCM", "GPCM"), cfg$noncenter_facet %in% facets,
    all(data$Score %in% seq(cfg$rating_min, cfg$rating_max)))
  if (!isTRUE(case$continuous_only)) stopifnot(
    all(is.finite(nodes)), all(weights > 0), abs(sum(weights) - 1) < 1e-12)
  indices <- lapply(facets, function(f) match(data[[f]], levels[[f]]))
  names(indices) <- facets
  if (isTRUE(case$group_rows)) {
    # Independent R reduction of repeated covariate cells. The original
    # rowwise objective remains available and is checked at the returned point.
    key <- do.call(paste, c(indices, sep=","))
    cell_index <- match(key,unique(key)); first <- match(seq_len(max(cell_index)),cell_index)
    nperson <- length(unique(person)); ncell <- length(first)
    observed_column <- (cell_index-1L)*k + data$Score-cfg$rating_min+1L
    counts <- matrix(tabulate(person+(observed_column-1L)*nperson, nbins=nperson*ncell*k),nperson,ncell*k)
  }
  point <- unlist(case$point)
  expand <- function(par) {
    cursor <- 0L; values <- list(); base <- numeric(n)
    for (f in facets) {
      count <- length(levels[[f]]) - as.integer(f != cfg$noncenter_facet)
      v <- par[cursor + seq_len(count)]; cursor <- cursor + count
      if (f != cfg$noncenter_facet) v <- c(v, -sum(v))
      values[[f]] <- v
      base <- base + cfg$facet_signs[[f]] * v[indices[[f]]]
    }
    nsteps <- if (cfg$model == "RSM") 1L else length(levels[[cfg$step_facet]])
    steps <- matrix(0, nsteps, k-1L)
    for (i in seq_len(nsteps)) {
      v <- par[cursor + seq_len(k-2L)]; cursor <- cursor + k-2L
      steps[i, ] <- c(v, -sum(v))
    }
    log_slopes <- numeric(0); slopes <- rep(1,n)
    if (cfg$model == "GPCM") {
      stopifnot(identical(cfg$slope_facet, cfg$step_facet))
      count <- length(levels[[cfg$slope_facet]])-1L
      v <- par[cursor + seq_len(count)]; cursor <- cursor+count
      log_slopes <- c(v,-sum(v))
      slopes <- exp(log_slopes)[indices[[cfg$slope_facet]]]
    }
    stopifnot(cursor + 1L == length(par))
    list(base = base, steps = steps, sigma = exp(par[cursor+1L]), facets = values,
         log_slopes=log_slopes, slopes=slopes)
  }
  if (isTRUE(case$continuous_only)) {
    x <- expand(point)
    cumulative <- t(apply(x$steps, 1L, function(s) c(0, cumsum(s))))
    row_steps <- if (cfg$model == "RSM") rep(1L, n) else indices[[cfg$step_facet]]
    continuous <- function(settings) {
      b <- settings$bound
      people <- lapply(sort(unique(person)), function(p) {
        rows <- which(person == p)
        log_joint <- function(z) vapply(z, function(v) {
          logits <- outer(x$base[rows] + x$sigma*v, 0:(k-1L)) - cumulative[row_steps[rows], , drop = FALSE]
          logits <- logits*x$slopes[rows]
          maximum <- apply(logits, 1L, max)
          lp <- logits - maximum - log(rowSums(exp(logits-maximum)))
          sum(lp[cbind(seq_along(rows), data$Score[rows]-cfg$rating_min+1L)]) + dnorm(v, log = TRUE)
        }, numeric(1))
        mode <- optimize(log_joint, c(-b, b), maximum = TRUE, tol = 1e-10)
        stopifnot(abs(mode$maximum) < b-.1)
        cuts <- sort(unique(c(-b, 0, mode$maximum, b)))
        integral <- lapply(0:2, function(j) {
          f <- function(z) exp(log_joint(z)-mode$objective)*z^j
          parts <- lapply(seq_len(length(cuts)-1L), function(i) integrate(f, cuts[i], cuts[i+1L],
              subdivisions = 500L, rel.tol = settings$rel_tol, abs.tol = settings$abs_tol))
          stopifnot(all(vapply(parts, function(v) identical(v$message, "OK"), logical(1))))
          c(value = sum(vapply(parts, `[[`, numeric(1), "value")),
            error = sum(vapply(parts, `[[`, numeric(1), "abs.error")))
        })
        mass <- integral[[1]][["value"]]
        ez <- integral[[2]][["value"]]/mass; ez2 <- integral[[3]][["value"]]/mass
        stopifnot(mass > 0, ez2 > ez^2)
        list(log_mass = mode$objective+log(mass), eap = x$sigma*ez,
             sd = x$sigma*sqrt(ez2-ez^2), ez2 = ez2, relative_error = integral[[1]][["error"]]/mass)
      })
      log_mass <- vapply(people, `[[`, numeric(1), "log_mass")
      list(nll = -sum(log_mass), eap = vapply(people, `[[`, numeric(1), "eap"),
           sd = vapply(people, `[[`, numeric(1), "sd"),
           density_log_sigma_score = -sum(vapply(people, `[[`, numeric(1), "ez2")-1),
           numeric_relative_mass_error_sum = sum(vapply(people, `[[`, numeric(1), "relative_error")),
           tail_relative_mass_bound_sum = sum(exp(log(2*pnorm(-b))-log_mass)))
    }
    values <- lapply(input$continuous_settings, continuous)
    cat(case$id, "continuous reference complete\n")
    return(list(id = case$id, continuous = values))
  }
  objective <- function(par, details=FALSE, grouped=isTRUE(case$group_rows)) {
    x <- expand(par)
    cumulative <- t(apply(x$steps, 1L, function(s) c(0, cumsum(s))))
    row_steps <- if (cfg$model == "RSM") rep(1L, n) else indices[[cfg$step_facet]]
    selected <- if (grouped) first else seq_len(n)
    nr <- length(selected)
    logits <- outer(rep(x$base[selected], q) + rep(x$sigma * nodes, each = nr), 0:(k-1L)) -
      cumulative[rep(row_steps[selected], q), , drop = FALSE]
    logits <- logits * rep(x$slopes[selected],q)
    maximum <- apply(logits, 1L, max)
    normalizer <- maximum + log(rowSums(exp(logits - maximum)))
    if (grouped) {
      category_log <- logits-normalizer
      cell_category_log <- matrix(aperm(array(category_log,c(nr,q,k)),c(3L,1L,2L)),nr*k,q)
      person_log <- counts %*% cell_category_log
    } else {
      observed <- logits[cbind(seq_len(n*q), rep(data$Score - cfg$rating_min + 1L, q))] - normalizer
      person_log <- rowsum(matrix(observed, n, q), person, reorder = TRUE)
    }
    lp <- sweep(person_log, 2L, log(weights), "+")
    maximum <- apply(lp, 1L, max)
    nll <- -sum(maximum + log(rowSums(exp(lp - maximum))))
    if (!details) return(nll)
    posterior <- exp(lp-maximum); posterior <- posterior/rowSums(posterior)
    ez <- drop(posterior %*% nodes)
    list(nll=nll, posterior=posterior, eap=x$sigma*ez,
         sd=x$sigma*sqrt(rowSums(posterior*(matrix(nodes,nrow(posterior),q,byrow=TRUE)-ez)^2)))
  }
  expanded <- function(par) {
    x <- expand(par)
    c(unlist(x$facets, use.names = FALSE), as.vector(t(x$steps)),
      x$log_slopes, if(length(x$log_slopes)) exp(x$log_slopes) else numeric(0), x$sigma)
  }
  # Differentiate at a zero-centered origin with an explicit absolute step.
  # Default relative steps collapse near zero-valued fitted coordinates and
  # amplify cancellation of O(N) NLL values in second differences.
  scale <- pmax(1, abs(point)); origin <- rep(0, length(point))
  shifted <- function(u) objective(point + scale * u)
  h <- numDeriv::hessian(shifted, origin, method.args = list(eps = .01, r = 4L)) / outer(scale, scale)
  h_half <- numDeriv::hessian(shifted, origin, method.args = list(eps = .005, r = 4L)) / outer(scale, scale)
  positive <- min(eigen(h_half, symmetric = TRUE, only.values = TRUE)$values) > 0
  covariance <- if (positive) solve(h_half) else NULL
  jacobian <- sweep(numDeriv::jacobian(function(u) expanded(point+scale*u), origin,
    method.args=list(eps=.001)),2L,scale,"/")
  cat(case$id, "reference complete\n")
  moments <- objective(point, details=TRUE)
  list(id = case$id, nll = moments$nll, eap=moments$eap, sd=moments$sd,
    posterior=moments$posterior,
    grouped_rows=isTRUE(case$group_rows), rowwise_nll_disagreement=abs(moments$nll-objective(point,grouped=FALSE)),
    gradient = numDeriv::grad(shifted, origin, method.args = list(eps = .001)) / scale,
    information = h_half,
    hessian_step_disagreement = max(abs(h - h_half)), positive_definite = positive,
    covariance = covariance, expanded = expanded(point), expanded_jacobian = jacobian,
    expanded_covariance = if (positive) jacobian %*% covariance %*% t(jacobian) else NULL)
})
jsonlite::write_json(list(classification = "DEVELOPMENT_ONLY", scientific_inference_ready = FALSE,
  runtime = list(R = R.version.string, numDeriv = as.character(packageVersion("numDeriv")),
    jsonlite = as.character(packageVersion("jsonlite"))), cases = results), args[2],
  digits = NA, pretty = TRUE, auto_unbox = TRUE, null = "null", na = "null")
