# MFRM estimation without FACETS or ordinal (RSM / PCM, JMLE / MML)

suppressPackageStartupMessages({
  library(shiny)
  library(tidyverse)
  library(gt)
  library(glue)
  library(plotly)
  library(psych)
})

# ---- math helpers ----
logsumexp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

weighted_mean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

get_weights <- function(df) {
  if (!is.null(df) && "Weight" %in% names(df)) {
    w <- suppressWarnings(as.numeric(df$Weight))
    w <- ifelse(is.finite(w) & w > 0, w, 0)
    return(w)
  }
  rep(1, nrow(df))
}

# Gauss-Hermite nodes/weights for standard normal integration
# Based on Golub-Welsch for Hermite polynomials
# Returns nodes and weights for phi(x) dx
#   integral f(x) phi(x) dx \approx sum w_i f(x_i)
gauss_hermite_normal <- function(n) {
  if (n < 1) stop("n must be >= 1")
  i <- seq_len(n - 1)
  a <- rep(0, n)
  b <- sqrt(i / 2)
  jmat <- diag(a) + diag(b, 1) + diag(b, -1)
  eig <- eigen(jmat, symmetric = TRUE)
  nodes <- eig$values
  weights <- sqrt(pi) * (eig$vectors[1, ]^2)
  # Convert from exp(-x^2) to standard normal
  nodes <- sqrt(2) * nodes
  weights <- weights / sqrt(pi)
  list(nodes = nodes, weights = weights)
}

center_sum_zero <- function(x) {
  if (length(x) == 0) return(x)
  x - mean(x)
}

expand_facet <- function(free, n_levels) {
  if (n_levels <= 1) return(rep(0, n_levels))
  c(free, -sum(free))
}

build_facet_constraint <- function(levels,
                                   anchors = NULL,
                                   groups = NULL,
                                   group_values = NULL,
                                   centered = TRUE) {
  lvl <- as.character(levels)
  anchors_vec <- rep(NA_real_, length(lvl))
  names(anchors_vec) <- lvl
  if (!is.null(anchors)) {
    anchors <- anchors[!is.na(names(anchors))]
    anchors_vec[names(anchors)] <- as.numeric(anchors)
  }

  groups_vec <- rep(NA_character_, length(lvl))
  names(groups_vec) <- lvl
  if (!is.null(groups)) {
    groups <- groups[!is.na(names(groups))]
    groups_vec[names(groups)] <- as.character(groups)
  }

  group_values_map <- numeric(0)
  if (!is.null(group_values)) {
    group_values <- group_values[!is.na(names(group_values))]
    group_values_map <- as.numeric(group_values)
    names(group_values_map) <- names(group_values)
  }

  spec <- list(
    levels = lvl,
    anchors = anchors_vec,
    groups = groups_vec,
    group_values = group_values_map,
    centered = centered
  )
  spec$n_params <- count_facet_params(spec)
  spec
}

count_facet_params <- function(spec) {
  anchors <- spec$anchors
  groups <- spec$groups
  free_idx <- which(is.na(anchors))
  if (length(free_idx) == 0) return(0)

  n_params <- 0
  group_ids <- unique(na.omit(groups[free_idx]))
  if (length(group_ids) > 0) {
    for (gid in group_ids) {
      group_levels <- which(groups == gid)
      free_in_group <- group_levels[is.na(anchors[group_levels])]
      k <- length(free_in_group)
      if (k > 1) n_params <- n_params + (k - 1)
    }
  }

  ungrouped_idx <- free_idx[is.na(groups[free_idx]) | groups[free_idx] == ""]
  m <- length(ungrouped_idx)
  if (m == 0) return(n_params)
  if (spec$centered) {
    n_params <- n_params + max(m - 1, 0)
  } else {
    n_params <- n_params + m
  }
  n_params
}

expand_facet_with_constraints <- function(free, spec) {
  out <- spec$anchors
  groups <- spec$groups
  group_values <- spec$group_values
  centered <- isTRUE(spec$centered)
  free_idx <- which(is.na(out))
  if (length(free_idx) == 0) return(out)

  used <- 0
  group_ids <- unique(na.omit(groups[free_idx]))
  if (length(group_ids) > 0) {
    for (gid in group_ids) {
      group_levels <- which(groups == gid)
      free_in_group <- group_levels[is.na(out[group_levels])]
      if (length(free_in_group) == 0) next
      group_value <- if (gid %in% names(group_values)) group_values[[gid]] else 0
      anchor_sum <- sum(out[group_levels], na.rm = TRUE)
      target_sum <- group_value * length(group_levels)
      k <- length(free_in_group)
      if (k == 1) {
        out[free_in_group] <- target_sum - anchor_sum
      } else {
        seg <- free[(used + 1):(used + k - 1)]
        used <- used + (k - 1)
        last_val <- target_sum - anchor_sum - sum(seg)
        out[free_in_group] <- c(seg, last_val)
      }
    }
  }

  ungrouped_idx <- free_idx[is.na(groups[free_idx]) | groups[free_idx] == ""]
  m <- length(ungrouped_idx)
  if (m == 0) return(out)
  if (centered) {
    if (m == 1) {
      out[ungrouped_idx] <- 0
    } else {
      seg <- free[(used + 1):(used + m - 1)]
      used <- used + (m - 1)
      out[ungrouped_idx] <- c(seg, -sum(seg))
    }
  } else {
    seg <- free[(used + 1):(used + m)]
    used <- used + m
    out[ungrouped_idx] <- seg
  }
  out
}

build_param_sizes <- function(config) {
  n_steps <- max(config$n_cat - 1, 0)
  sizes <- list(
    theta = if (config$method == "JMLE") config$theta_spec$n_params else 0
  )
  for (facet in config$facet_names) {
    sizes[[facet]] <- config$facet_specs[[facet]]$n_params
  }
  if (config$model == "RSM") {
    sizes$steps <- n_steps
  } else {
    if (is.null(config$step_facet) || !config$step_facet %in% config$facet_names) {
      stop("PCM requires a valid step facet.")
    }
    sizes$steps <- length(config$facet_levels[[config$step_facet]]) * n_steps
  }
  sizes
}

split_params <- function(par, sizes) {
  out <- list()
  idx <- 1
  for (nm in names(sizes)) {
    k <- sizes[[nm]]
    if (k == 0) {
      out[[nm]] <- numeric(0)
    } else {
      out[[nm]] <- par[idx:(idx + k - 1)]
      idx <- idx + k
    }
  }
  out
}

expand_params <- function(par, sizes, config) {
  parts <- split_params(par, sizes)
  theta <- if (config$method == "JMLE") {
    expand_facet_with_constraints(parts$theta, config$theta_spec)
  } else {
    numeric(0)
  }

  facets <- lapply(config$facet_names, function(facet) {
    expand_facet_with_constraints(parts[[facet]], config$facet_specs[[facet]])
  })
  names(facets) <- config$facet_names

  if (config$model == "RSM") {
    steps <- center_sum_zero(parts$steps)
    steps_mat <- NULL
  } else {
    n_levels <- length(config$facet_levels[[config$step_facet]])
    if (n_levels == 0 || config$n_cat <= 1) {
      steps_mat <- matrix(0, nrow = n_levels, ncol = max(config$n_cat - 1, 0))
    } else {
      steps_mat <- matrix(parts$steps, nrow = n_levels, byrow = TRUE)
      steps_mat <- t(apply(steps_mat, 1, center_sum_zero))
    }
    steps <- NULL
  }

  list(theta = theta, facets = facets, steps = steps, steps_mat = steps_mat)
}

# ---- data preparation ----
prepare_mfrm_data <- function(data, person_col, facet_cols, score_col,
                              rating_min = NULL, rating_max = NULL,
                              weight_col = NULL, keep_original = FALSE) {
  required <- c(person_col, facet_cols, score_col)
  if (!is.null(weight_col)) {
    required <- c(required, weight_col)
  }
  if (length(unique(required)) != length(required)) {
    stop("Person/score/facet columns must be distinct (no duplicates).")
  }
  if (!all(required %in% names(data))) {
    stop("Some selected columns are not in the data.")
  }
  if (any(duplicated(names(data)))) {
    dupes <- names(data)[duplicated(names(data))]
    if (any(required %in% dupes)) {
      stop("Selected columns include duplicate names in the data. Please rename columns to be unique.")
    }
  }
  if (length(facet_cols) == 0) {
    stop("Select at least one facet column.")
  }

  cols <- c(person_col, facet_cols, score_col)
  if (!is.null(weight_col)) {
    cols <- c(cols, weight_col)
  }
  df <- data |>
    select(all_of(cols)) |>
    rename(
      Person = all_of(person_col),
      Score = all_of(score_col)
    )
  if (!is.null(weight_col)) {
    df <- df |> rename(Weight = all_of(weight_col))
  }

  df <- df |>
    mutate(
      Person = as.character(Person),
      across(all_of(facet_cols), ~ as.character(.x)),
      Score = suppressWarnings(as.numeric(Score))
    )
  if (!"Weight" %in% names(df)) {
    df <- df |> mutate(Weight = 1)
  } else {
    df <- df |> mutate(Weight = suppressWarnings(as.numeric(Weight)))
  }

  df <- df |>
    tidyr::drop_na() |>
    filter(Weight > 0)

  df <- df |>
    mutate(Score = as.integer(Score))

  if (is.null(rating_min)) rating_min <- min(df$Score, na.rm = TRUE)
  if (is.null(rating_max)) rating_max <- max(df$Score, na.rm = TRUE)

  if (!isTRUE(keep_original)) {
    score_vals <- sort(unique(df$Score))
    expected_vals <- seq(rating_min, rating_max)
    if (!identical(score_vals, expected_vals)) {
      df <- df |>
        mutate(Score = match(Score, score_vals) + rating_min - 1)
      rating_max <- rating_min + length(score_vals) - 1
    }
  }

  df <- df |>
    mutate(score_k = Score - rating_min)

  df <- df |>
    mutate(
      Person = factor(Person),
      across(all_of(facet_cols), ~ factor(.x))
    )

  facet_names <- facet_cols
  facet_levels <- lapply(facet_names, function(f) levels(df[[f]]))
  names(facet_levels) <- facet_names

  list(
    data = df,
    rating_min = rating_min,
    rating_max = rating_max,
    facet_names = facet_names,
    levels = c(list(Person = levels(df$Person)), facet_levels),
    weight_col = if (!is.null(weight_col)) weight_col else NULL
  )
}

build_indices <- function(prep, step_facet = NULL) {
  df <- prep$data
  facets_idx <- lapply(prep$facet_names, function(f) as.integer(df[[f]]))
  names(facets_idx) <- prep$facet_names
  step_idx <- if (!is.null(step_facet)) {
    as.integer(df[[step_facet]])
  } else {
    NULL
  }
  list(
    person = as.integer(df$Person),
    facets = facets_idx,
    step_idx = step_idx,
    score_k = as.integer(df$score_k),
    weight = suppressWarnings(as.numeric(df$Weight))
  )
}

sample_mfrm_data <- function(seed = 20240131) {
  set.seed(seed)
  persons <- paste0("P", sprintf("%02d", 1:36))
  raters <- paste0("R", 1:3)
  tasks <- paste0("T", 1:4)
  criteria <- paste0("C", 1:3)
  df <- expand_grid(
    Person = persons,
    Rater = raters,
    Task = tasks,
    Criterion = criteria
  )
  ability <- rnorm(length(persons), 0, 1)
  rater_eff <- c(-0.4, 0, 0.4)
  task_eff <- seq(-0.5, 0.5, length.out = length(tasks))
  crit_eff <- c(-0.3, 0, 0.3)
  eta <- ability[match(df$Person, persons)] -
    rater_eff[match(df$Rater, raters)] -
    task_eff[match(df$Task, tasks)] -
    crit_eff[match(df$Criterion, criteria)]
  raw <- eta + rnorm(nrow(df), 0, 0.6)
  score <- as.integer(cut(
    raw,
    breaks = c(-Inf, -1.0, -0.3, 0.3, 1.0, Inf),
    labels = 1:5
  ))
  df$Score <- score
  df
}

format_tab_template <- function(df) {
  char_df <- df |> mutate(across(everything(), ~ replace_na(as.character(.x), "")))
  widths <- vapply(seq_along(char_df), function(i) {
    max(nchar(c(names(char_df)[i], char_df[[i]])), na.rm = TRUE)
  }, integer(1))
  format_row <- function(row_vec) {
    padded <- mapply(function(value, width) {
      value <- ifelse(is.na(value), "", value)
      stringr::str_pad(value, width = width, side = "right")
    }, row_vec, widths, SIMPLIFY = TRUE)
    paste(padded, collapse = "\t")
  }
  header <- format_row(names(char_df))
  rows <- apply(char_df, 1, format_row)
  paste(c(header, rows), collapse = "\n")
}

template_tab_source_demo <- sample_mfrm_data(seed = 20240131) |>
  slice_head(n = 24)
template_tab_source_toy <- sample_mfrm_data(seed = 20240131) |>
  slice_head(n = 8)
template_tab_text <- format_tab_template(template_tab_source_demo)
template_tab_text_toy <- format_tab_template(template_tab_source_toy)
template_header_text <- format_tab_template(template_tab_source_demo[0, ])
download_sample_data <- sample_mfrm_data(seed = 20240131)

guess_col <- function(cols, patterns, fallback = 1) {
  if (length(cols) == 0) return(character(0))
  hit <- which(stringr::str_detect(tolower(cols), paste(patterns, collapse = "|")))
  if (length(hit) > 0) return(cols[hit[1]])
  cols[min(fallback, length(cols))]
}

plotly_style_axes <- function(p,
                              margin = list(t = 60, r = 40, b = 60, l = 60),
                              base_size = 11) {
  p |>
    plotly::layout(
      font = list(size = base_size),
      margin = margin,
      xaxis = list(automargin = TRUE, title = list(standoff = 10)),
      yaxis = list(automargin = TRUE, title = list(standoff = 10))
    )
}

plotly_compact_legend <- function(p, legend_y = -0.25, bottom_margin = 90, base_size = 11) {
  p <- plotly_style_axes(
    p,
    margin = list(t = 50, r = 40, b = bottom_margin, l = 60),
    base_size = base_size
  )
  plotly::layout(
    p,
    legend = list(
      orientation = "h",
      x = 0.5,
      xanchor = "center",
      y = legend_y,
      yanchor = "top"
    )
  )
}

plotly_legend_top <- function(p, top_margin = 90, legend_y = 1.12, base_size = 11) {
  p <- plotly_style_axes(
    p,
    margin = list(t = top_margin, r = 40, b = 70, l = 60),
    base_size = base_size
  )
  plotly::layout(
    p,
    legend = list(
      orientation = "h",
      x = 0.5,
      xanchor = "center",
      y = legend_y,
      yanchor = "bottom"
    )
  )
}

plotly_legend_right <- function(p, right_margin = 140, base_size = 11) {
  p <- plotly_style_axes(
    p,
    margin = list(t = 60, r = right_margin, b = 70, l = 60),
    base_size = base_size
  )
  plotly::layout(
    p,
    legend = list(
      orientation = "v",
      x = 1.02,
      xanchor = "left",
      y = 1,
      yanchor = "top"
    )
  )
}

plotly_basic <- function(p, base_size = 11) {
  plotly_style_axes(p, base_size = base_size)
}

truncate_label <- function(x, width = 28) {
  stringr::str_trunc(as.character(x), width = width)
}

facet_report_id <- function(facet) {
  paste0("facet_report_", stringr::str_replace_all(as.character(facet), "[^A-Za-z0-9]", "_"))
}

# ---- likelihoods ----
loglik_rsm <- function(eta, score_k, step_cum, weight = NULL) {
  n <- length(eta)
  if (n == 0) return(0)
  k_cat <- length(step_cum)
  eta_mat <- outer(eta, 0:(k_cat - 1))
  log_num <- eta_mat - matrix(step_cum, nrow = n, ncol = k_cat, byrow = TRUE)
  row_max <- apply(log_num, 1, max)
  log_denom <- row_max + log(rowSums(exp(log_num - row_max)))
  log_num_obs <- log_num[cbind(seq_len(n), score_k + 1)]
  diff <- log_num_obs - log_denom
  if (is.null(weight)) {
    sum(diff)
  } else {
    sum(diff * weight)
  }
}

loglik_pcm <- function(eta, score_k, step_cum_mat, criterion_idx, weight = NULL) {
  n <- length(eta)
  if (n == 0) return(0)
  k_cat <- ncol(step_cum_mat)
  total <- 0
  for (c_idx in seq_len(nrow(step_cum_mat))) {
    rows <- which(criterion_idx == c_idx)
    if (length(rows) == 0) next
    eta_c <- eta[rows]
    step_cum <- step_cum_mat[c_idx, ]
    eta_mat <- outer(eta_c, 0:(k_cat - 1))
    log_num <- eta_mat - matrix(step_cum, nrow = length(rows), ncol = k_cat, byrow = TRUE)
    row_max <- apply(log_num, 1, max)
    log_denom <- row_max + log(rowSums(exp(log_num - row_max)))
    log_num_obs <- log_num[cbind(seq_len(length(rows)), score_k[rows] + 1)]
    diff <- log_num_obs - log_denom
    if (is.null(weight)) {
      total <- total + sum(diff)
    } else {
      total <- total + sum(diff * weight[rows])
    }
  }
  total
}

category_prob_rsm <- function(eta, step_cum) {
  n <- length(eta)
  if (n == 0) return(matrix(0, nrow = 0, ncol = length(step_cum)))
  k_cat <- length(step_cum)
  eta_mat <- outer(eta, 0:(k_cat - 1))
  log_num <- eta_mat - matrix(step_cum, nrow = n, ncol = k_cat, byrow = TRUE)
  row_max <- apply(log_num, 1, max)
  log_denom <- row_max + log(rowSums(exp(log_num - row_max)))
  exp(log_num - matrix(log_denom, nrow = n, ncol = k_cat))
}

category_prob_pcm <- function(eta, step_cum_mat, criterion_idx) {
  n <- length(eta)
  if (n == 0) return(matrix(0, nrow = 0, ncol = ncol(step_cum_mat)))
  k_cat <- ncol(step_cum_mat)
  probs <- matrix(0, nrow = n, ncol = k_cat)
  for (c_idx in seq_len(nrow(step_cum_mat))) {
    rows <- which(criterion_idx == c_idx)
    if (length(rows) == 0) next
    step_cum <- step_cum_mat[c_idx, ]
    eta_c <- eta[rows]
    eta_mat <- outer(eta_c, 0:(k_cat - 1))
    log_num <- eta_mat - matrix(step_cum, nrow = length(rows), ncol = k_cat, byrow = TRUE)
    row_max <- apply(log_num, 1, max)
    log_denom <- row_max + log(rowSums(exp(log_num - row_max)))
    probs[rows, ] <- exp(log_num - matrix(log_denom, nrow = length(rows), ncol = k_cat))
  }
  probs
}

zstd_from_mnsq <- function(mnsq, df, whexact = FALSE) {
  if (!is.finite(mnsq) || !is.finite(df) || df <= 0) return(NA_real_)
  if (isTRUE(whexact)) {
    return((mnsq - 1) * sqrt(df / 2))
  }
  (mnsq^(1/3) - (1 - 2 / (9 * df))) / sqrt(2 / (9 * df))
}

compute_base_eta <- function(idx, params, config) {
  eta <- rep(0, length(idx$score_k))
  facet_signs <- config$facet_signs
  for (facet in config$facet_names) {
    sign <- if (!is.null(facet_signs) && !is.null(facet_signs[[facet]])) {
      facet_signs[[facet]]
    } else {
      -1
    }
    eta <- eta + sign * params$facets[[facet]][idx$facets[[facet]]]
  }
  eta
}

compute_eta <- function(idx, params, config, theta_override = NULL) {
  theta <- if (!is.null(theta_override)) theta_override else params$theta
  eta <- if (length(theta) == 0) {
    rep(0, length(idx$score_k))
  } else {
    theta[idx$person]
  }
  eta + compute_base_eta(idx, params, config)
}

mfrm_loglik_jmle <- function(par, idx, config, sizes) {
  params <- expand_params(par, sizes, config)
  eta <- compute_eta(idx, params, config)

  if (config$model == "RSM") {
    step_cum <- c(0, cumsum(params$steps))
    ll <- loglik_rsm(eta, idx$score_k, step_cum, weight = idx$weight)
  } else {
    step_cum_mat <- t(apply(params$steps_mat, 1, function(x) c(0, cumsum(x))))
    ll <- loglik_pcm(eta, idx$score_k, step_cum_mat, idx$step_idx, weight = idx$weight)
  }
  -ll
}

mfrm_loglik_mml <- function(par, idx, config, sizes, quad) {
  params <- expand_params(par, sizes, config)
  n <- length(idx$score_k)
  if (n == 0) return(0)

  base_eta <- compute_base_eta(idx, params, config)
  rows_by_person <- split(seq_len(n), idx$person)
  log_w <- log(quad$weights)

  if (config$model == "RSM") {
    step_cum <- c(0, cumsum(params$steps))
    ll_person <- vapply(rows_by_person, function(rows) {
      base <- base_eta[rows]
      score_k <- idx$score_k[rows]
      w <- if (!is.null(idx$weight)) idx$weight[rows] else NULL
      ll_nodes <- vapply(quad$nodes, function(theta) {
        loglik_rsm(theta + base, score_k, step_cum, weight = w)
      }, numeric(1))
      logsumexp(log_w + ll_nodes)
    }, numeric(1))
  } else {
    step_cum_mat <- t(apply(params$steps_mat, 1, function(x) c(0, cumsum(x))))
    ll_person <- vapply(rows_by_person, function(rows) {
      base <- base_eta[rows]
      score_k <- idx$score_k[rows]
      crit <- idx$step_idx[rows]
      w <- if (!is.null(idx$weight)) idx$weight[rows] else NULL
      ll_nodes <- vapply(quad$nodes, function(theta) {
        loglik_pcm(theta + base, score_k, step_cum_mat, crit, weight = w)
      }, numeric(1))
      logsumexp(log_w + ll_nodes)
    }, numeric(1))
  }

  -sum(ll_person)
}

compute_person_eap <- function(idx, config, params, quad) {
  n <- length(idx$score_k)
  if (n == 0) {
    return(tibble(Person = character(0), Estimate = numeric(0), SD = numeric(0)))
  }
  base_eta <- compute_base_eta(idx, params, config)
  rows_by_person <- split(seq_len(n), idx$person)
  log_w <- log(quad$weights)

  if (config$model == "RSM") {
    step_cum <- c(0, cumsum(params$steps))
    estimates <- lapply(rows_by_person, function(rows) {
      base <- base_eta[rows]
      score_k <- idx$score_k[rows]
      w <- if (!is.null(idx$weight)) idx$weight[rows] else NULL
      ll_nodes <- vapply(quad$nodes, function(theta) {
        loglik_rsm(theta + base, score_k, step_cum, weight = w)
      }, numeric(1))
      log_post <- log_w + ll_nodes
      log_post <- log_post - logsumexp(log_post)
      post_w <- exp(log_post)
      eap <- sum(quad$nodes * post_w)
      sd <- sqrt(sum((quad$nodes - eap)^2 * post_w))
      c(eap = eap, sd = sd)
    })
  } else {
    step_cum_mat <- t(apply(params$steps_mat, 1, function(x) c(0, cumsum(x))))
    estimates <- lapply(rows_by_person, function(rows) {
      base <- base_eta[rows]
      score_k <- idx$score_k[rows]
      crit <- idx$step_idx[rows]
      w <- if (!is.null(idx$weight)) idx$weight[rows] else NULL
      ll_nodes <- vapply(quad$nodes, function(theta) {
        loglik_pcm(theta + base, score_k, step_cum_mat, crit, weight = w)
      }, numeric(1))
      log_post <- log_w + ll_nodes
      log_post <- log_post - logsumexp(log_post)
      post_w <- exp(log_post)
      eap <- sum(quad$nodes * post_w)
      sd <- sqrt(sum((quad$nodes - eap)^2 * post_w))
      c(eap = eap, sd = sd)
    })
  }

  est_mat <- do.call(rbind, estimates)
  tibble(Estimate = est_mat[, "eap"], SD = est_mat[, "sd"])
}

prepare_constraint_specs <- function(prep,
                                     anchor_df = NULL,
                                     group_anchor_df = NULL,
                                     noncenter_facet = "Person",
                                     dummy_facets = character(0)) {
  facet_names <- prep$facet_names
  all_facets <- c("Person", facet_names)

  anchor_df <- if (!is.null(anchor_df) && nrow(anchor_df) > 0) {
    anchor_df |>
      mutate(
        Facet = trimws(as.character(Facet)),
        Level = trimws(as.character(Level)),
        Anchor = as.numeric(Anchor)
      ) |>
      filter(Facet %in% all_facets, !is.na(Level))
  } else {
    tibble()
  }

  group_anchor_df <- if (!is.null(group_anchor_df) && nrow(group_anchor_df) > 0) {
    group_anchor_df |>
      mutate(
        Facet = trimws(as.character(Facet)),
        Level = trimws(as.character(Level)),
        Group = trimws(as.character(Group)),
        GroupValue = as.numeric(GroupValue)
      ) |>
      filter(Facet %in% all_facets, !is.na(Level), !is.na(Group))
  } else {
    tibble()
  }

  anchor_map <- setNames(vector("list", length(all_facets)), all_facets)
  group_map <- setNames(vector("list", length(all_facets)), all_facets)
  group_values <- setNames(vector("list", length(all_facets)), all_facets)

  for (facet in all_facets) {
    levels <- prep$levels[[facet]]
    if (!is.null(levels)) {
      if (facet %in% dummy_facets) {
        anchor_map[[facet]] <- setNames(rep(0, length(levels)), levels)
        group_map[[facet]] <- NULL
        group_values[[facet]] <- numeric(0)
        next
      }
      if (nrow(anchor_df) > 0) {
        df <- filter(anchor_df, Facet == facet)
        if (nrow(df) > 0) {
          anchors <- setNames(df$Anchor, df$Level)
          anchors <- anchors[names(anchors) %in% levels]
          if (length(anchors) > 0) anchor_map[[facet]] <- anchors
        }
      }

      if (nrow(group_anchor_df) > 0) {
        df <- filter(group_anchor_df, Facet == facet)
        if (nrow(df) > 0) {
          groups <- setNames(df$Group, df$Level)
          groups <- groups[names(groups) %in% levels]
          if (length(groups) > 0) group_map[[facet]] <- groups

          group_vals <- df |>
            select(Group, GroupValue) |>
            distinct() |>
            filter(!is.na(Group)) |>
            group_by(Group) |>
            summarize(
              GroupValue = {
                vals <- GroupValue[!is.na(GroupValue)]
                if (length(vals) == 0) NA_real_ else vals[1]
              },
              .groups = "drop"
            ) |>
            mutate(GroupValue = ifelse(is.na(GroupValue), 0, GroupValue))
          if (nrow(group_vals) > 0) {
            group_values[[facet]] <- setNames(group_vals$GroupValue, group_vals$Group)
          }
        }
      }
    }
  }

  theta_spec <- build_facet_constraint(
    levels = prep$levels$Person,
    anchors = anchor_map$Person,
    groups = group_map$Person,
    group_values = group_values$Person,
    centered = !(identical(noncenter_facet, "Person"))
  )

  facet_specs <- lapply(facet_names, function(facet) {
    build_facet_constraint(
      levels = prep$levels[[facet]],
      anchors = anchor_map[[facet]],
      groups = group_map[[facet]],
      group_values = group_values[[facet]],
      centered = !identical(noncenter_facet, facet)
    )
  })
  names(facet_specs) <- facet_names

  anchor_summary <- tibble(
    Facet = all_facets,
    AnchoredLevels = vapply(anchor_map, function(x) length(x), integer(1)),
    GroupAnchors = vapply(group_map, function(x) length(unique(x)), integer(1)),
    DummyFacet = all_facets %in% dummy_facets
  )

  list(
    theta_spec = theta_spec,
    facet_specs = facet_specs,
    anchor_summary = anchor_summary
  )
}

read_flexible_table <- function(text_value, file_input) {
  if (!is.null(file_input) && !is.null(file_input$datapath)) {
    sep <- ifelse(grepl("\\.tsv$|\\.txt$", file_input$name, ignore.case = TRUE), "\t", ",")
    return(read.csv(file_input$datapath,
                    sep = sep,
                    header = TRUE,
                    stringsAsFactors = FALSE,
                    check.names = FALSE))
  }
  if (is.null(text_value)) return(tibble())
  text_value <- trimws(text_value)
  if (!nzchar(text_value)) return(tibble())
  sep <- if (grepl("\t", text_value)) {
    "\t"
  } else if (grepl(";", text_value)) {
    ";"
  } else {
    ","
  }
  read.csv(text = text_value,
           sep = sep,
           header = TRUE,
           stringsAsFactors = FALSE,
           check.names = FALSE)
}

normalize_anchor_df <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(tibble())
  nm <- tolower(names(df))
  facet_col <- which(nm %in% c("facet", "facets"))
  level_col <- which(nm %in% c("level", "element", "label"))
  anchor_col <- which(nm %in% c("anchor", "value", "measure"))
  if (length(facet_col) == 0 || length(level_col) == 0 || length(anchor_col) == 0) {
    return(tibble())
  }
  tibble(
    Facet = as.character(df[[facet_col[1]]]),
    Level = as.character(df[[level_col[1]]]),
    Anchor = suppressWarnings(as.numeric(df[[anchor_col[1]]]))
  ) |>
    filter(!is.na(Facet), !is.na(Level))
}

normalize_group_anchor_df <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(tibble())
  nm <- tolower(names(df))
  facet_col <- which(nm %in% c("facet", "facets"))
  level_col <- which(nm %in% c("level", "element", "label"))
  group_col <- which(nm %in% c("group", "subset"))
  value_col <- which(nm %in% c("groupvalue", "value", "anchor"))
  if (length(facet_col) == 0 || length(level_col) == 0 || length(group_col) == 0 || length(value_col) == 0) {
    return(tibble())
  }
  tibble(
    Facet = as.character(df[[facet_col[1]]]),
    Level = as.character(df[[level_col[1]]]),
    Group = as.character(df[[group_col[1]]]),
    GroupValue = suppressWarnings(as.numeric(df[[value_col[1]]]))
  ) |>
    filter(!is.na(Facet), !is.na(Level), !is.na(Group))
}

# ---- estimation wrapper ----
mfrm_estimate <- function(data, person_col, facet_cols, score_col,
                          rating_min = NULL, rating_max = NULL,
                          weight_col = NULL, keep_original = FALSE,
                          model = c("RSM", "PCM"), method = c("JMLE", "MML"),
                          step_facet = NULL,
                          anchor_df = NULL,
                          group_anchor_df = NULL,
                          noncenter_facet = "Person",
                          dummy_facets = character(0),
                          positive_facets = character(0),
                          quad_points = 15, maxit = 400, reltol = 1e-6) {
  model <- match.arg(model)
  method <- match.arg(method)

  prep <- prepare_mfrm_data(
    data,
    person_col = person_col,
    facet_cols = facet_cols,
    score_col = score_col,
    rating_min = rating_min,
    rating_max = rating_max,
    weight_col = weight_col,
    keep_original = keep_original
  )
  if (model == "PCM") {
    if (is.null(step_facet)) {
      step_facet <- prep$facet_names[1]
    }
    if (!step_facet %in% prep$facet_names) {
      stop("Selected step facet is not in the facet list.")
    }
  } else {
    step_facet <- NULL
  }
  idx <- build_indices(prep, step_facet = step_facet)

  n_person <- length(prep$levels$Person)
  facet_levels <- prep$levels[prep$facet_names]
  n_cat <- prep$rating_max - prep$rating_min + 1

  if (is.null(noncenter_facet) || !noncenter_facet %in% c("Person", prep$facet_names)) {
    noncenter_facet <- "Person"
  }
  if (is.null(dummy_facets)) dummy_facets <- character(0)
  dummy_facets <- intersect(dummy_facets, c("Person", prep$facet_names))

  positive_facets <- intersect(positive_facets, prep$facet_names)
  facet_signs <- setNames(
    ifelse(prep$facet_names %in% positive_facets, 1, -1),
    prep$facet_names
  )

  config <- list(
    model = model,
    method = method,
    n_person = n_person,
    n_cat = n_cat,
    facet_names = prep$facet_names,
    facet_levels = facet_levels,
    step_facet = step_facet
  )
  config$weight_col <- if (!is.null(weight_col)) weight_col else NULL
  config$positive_facets <- positive_facets
  config$facet_signs <- facet_signs

  constraint_specs <- prepare_constraint_specs(
    prep = prep,
    anchor_df = anchor_df,
    group_anchor_df = group_anchor_df,
    noncenter_facet = noncenter_facet,
    dummy_facets = dummy_facets
  )
  config$theta_spec <- constraint_specs$theta_spec
  config$facet_specs <- constraint_specs$facet_specs
  config$noncenter_facet <- noncenter_facet
  config$dummy_facets <- dummy_facets
  config$anchor_summary <- constraint_specs$anchor_summary

  sizes <- build_param_sizes(config)

  step_init <- if (n_cat > 1) {
    seq(-1, 1, length.out = n_cat - 1)
  } else {
    numeric(0)
  }

  facet_starts <- unlist(lapply(config$facet_names, function(f) rep(0, sizes[[f]])))
  start <- c(
    rep(0, sizes$theta),
    facet_starts,
    if (model == "RSM") step_init else rep(step_init, length(facet_levels[[step_facet]]))
  )

  if (method == "JMLE") {
    opt <- optim(
      start,
      fn = mfrm_loglik_jmle,
      idx = idx,
      config = config,
      sizes = sizes,
      method = "BFGS",
      control = list(maxit = maxit, reltol = reltol)
    )
  } else {
    quad <- gauss_hermite_normal(quad_points)
    opt <- optim(
      start,
      fn = mfrm_loglik_mml,
      idx = idx,
      config = config,
      sizes = sizes,
      quad = quad,
      method = "BFGS",
      control = list(maxit = maxit, reltol = reltol)
    )
  }

  params <- expand_params(opt$par, sizes, config)

  if (method == "MML") {
    quad <- gauss_hermite_normal(quad_points)
    person_tbl <- compute_person_eap(idx, config, params, quad) |>
      mutate(Person = prep$levels$Person) |>
      select(Person, Estimate, SD)
  } else {
    person_tbl <- tibble(
      Person = prep$levels$Person,
      Estimate = params$theta
    )
  }

  facet_tbls <- lapply(config$facet_names, function(facet) {
    tibble(Level = prep$levels[[facet]], Estimate = params$facets[[facet]]) |>
      mutate(Facet = facet, .before = 1)
  })
  facet_tbl <- bind_rows(facet_tbls)

  step_tbl <- if (model == "RSM") {
    tibble(
      Step = paste0("Step_", seq_len(n_cat - 1)),
      Estimate = params$steps
    )
  } else {
    expand_grid(
      StepFacet = prep$levels[[step_facet]],
      Step = paste0("Step_", seq_len(n_cat - 1))
    ) |>
      mutate(Estimate = as.vector(t(params$steps_mat)))
  }

  k_params <- sum(unlist(sizes))
  loglik <- -opt$value
  n_obs <- if (!is.null(weight_col) && "Weight" %in% names(prep$data)) {
    sum(prep$data$Weight, na.rm = TRUE)
  } else {
    nrow(prep$data)
  }
  aic <- 2 * k_params - 2 * loglik
  bic <- log(n_obs) * k_params - 2 * loglik

  summary_tbl <- tibble(
    Model = model,
    Method = method,
    N = n_obs,
    Persons = n_person,
    Facets = length(config$facet_names),
    Categories = n_cat,
    LogLik = loglik,
    AIC = aic,
    BIC = bic,
    Converged = opt$convergence == 0,
    Iterations = opt$counts[["function"]]
  )

  list(
    summary = summary_tbl,
    facets = list(
      person = person_tbl,
      others = facet_tbl
    ),
    steps = step_tbl,
    config = config,
    prep = prep,
    opt = opt
  )
}

expected_score_table <- function(res) {
  prep <- res$prep
  idx <- build_indices(prep, step_facet = res$config$step_facet)
  config <- res$config
  sizes <- build_param_sizes(config)
  params <- expand_params(res$opt$par, sizes, config)
  theta_hat <- if (config$method == "JMLE") {
    params$theta
  } else {
    res$facets$person$Estimate
  }
  eta <- compute_eta(idx, params, config, theta_override = theta_hat)

  if (config$model == "RSM") {
    step_cum <- c(0, cumsum(params$steps))
    probs <- category_prob_rsm(eta, step_cum)
  } else {
    step_cum_mat <- t(apply(params$steps_mat, 1, function(x) c(0, cumsum(x))))
    probs <- category_prob_pcm(eta, step_cum_mat, idx$step_idx)
  }
  k_vals <- 0:(ncol(probs) - 1)
  expected_k <- as.vector(probs %*% k_vals)
  tibble(
    Observed = prep$data$Score,
    Expected = prep$rating_min + expected_k
  )
}

compute_obs_table <- function(res) {
  prep <- res$prep
  config <- res$config
  idx <- build_indices(prep, step_facet = config$step_facet)
  sizes <- build_param_sizes(config)
  params <- expand_params(res$opt$par, sizes, config)
  theta_hat <- if (config$method == "JMLE") {
    params$theta
  } else {
    res$facets$person$Estimate
  }
  person_levels <- prep$levels$Person
  person_measure <- if (config$method == "JMLE") {
    params$theta
  } else {
    res$facets$person$Estimate[match(person_levels, res$facets$person$Person)]
  }
  person_measure_by_row <- person_measure[idx$person]
  eta <- compute_eta(idx, params, config, theta_override = theta_hat)

  if (config$model == "RSM") {
    step_cum <- c(0, cumsum(params$steps))
    probs <- category_prob_rsm(eta, step_cum)
  } else {
    step_cum_mat <- t(apply(params$steps_mat, 1, function(x) c(0, cumsum(x))))
    probs <- category_prob_pcm(eta, step_cum_mat, idx$step_idx)
  }

  k_vals <- 0:(ncol(probs) - 1)
  expected_k <- as.vector(probs %*% k_vals)
  var_k <- as.vector(probs %*% (k_vals^2)) - expected_k^2
  var_k <- ifelse(var_k <= 1e-10, NA_real_, var_k)
  resid_k <- idx$score_k - expected_k
  std_sq <- resid_k^2 / var_k

  prep$data |>
    mutate(
      PersonMeasure = person_measure_by_row,
      Observed = prep$rating_min + idx$score_k,
      Expected = prep$rating_min + expected_k,
      Var = var_k,
      Residual = Observed - Expected,
      StdResidual = Residual / sqrt(Var),
      StdSq = std_sq
    )
}

compute_prob_matrix <- function(res) {
  prep <- res$prep
  config <- res$config
  idx <- build_indices(prep, step_facet = config$step_facet)
  sizes <- build_param_sizes(config)
  params <- expand_params(res$opt$par, sizes, config)
  theta_hat <- if (config$method == "JMLE") {
    params$theta
  } else {
    res$facets$person$Estimate
  }
  eta <- compute_eta(idx, params, config, theta_override = theta_hat)
  if (config$model == "RSM") {
    step_cum <- c(0, cumsum(params$steps))
    probs <- category_prob_rsm(eta, step_cum)
  } else {
    step_cum_mat <- t(apply(params$steps_mat, 1, function(x) c(0, cumsum(x))))
    probs <- category_prob_pcm(eta, step_cum_mat, idx$step_idx)
  }
  probs
}

compute_scorefile <- function(res) {
  obs <- compute_obs_table(res)
  if (nrow(obs) == 0) return(tibble())
  probs <- compute_prob_matrix(res)
  if (is.null(probs) || nrow(probs) != nrow(obs)) return(obs)
  cat_vals <- seq(res$prep$rating_min, res$prep$rating_max)
  prob_df <- as_tibble(probs)
  names(prob_df) <- paste0("P_", cat_vals)
  max_idx <- apply(probs, 1, which.max)
  max_prob <- probs[cbind(seq_len(nrow(probs)), max_idx)]
  most_likely <- cat_vals[max_idx]
  base_cols <- c("Person", res$config$facet_names)
  if ("Weight" %in% names(obs)) {
    base_cols <- c(base_cols, "Weight")
  }
  base_cols <- c(base_cols, "Score", "Observed", "Expected", "Residual", "StdResidual", "Var", "PersonMeasure")
  bind_cols(
    obs |>
      select(all_of(base_cols)),
    prob_df
  ) |>
    mutate(
      MostLikely = most_likely,
      MaxProb = max_prob
    )
}

compute_residual_file <- function(res) {
  obs <- compute_obs_table(res)
  if (nrow(obs) == 0) return(tibble())
  base_cols <- c("Person", res$config$facet_names)
  if ("Weight" %in% names(obs)) {
    base_cols <- c(base_cols, "Weight")
  }
  base_cols <- c(base_cols, "Score", "Observed", "Expected", "Residual", "StdResidual", "Var", "StdSq", "PersonMeasure")
  obs |>
    select(all_of(base_cols))
}

calc_overall_fit <- function(obs_df, whexact = FALSE) {
  w <- get_weights(obs_df)
  infit <- sum(obs_df$StdSq * obs_df$Var * w, na.rm = TRUE) / sum(obs_df$Var * w, na.rm = TRUE)
  outfit <- sum(obs_df$StdSq * w, na.rm = TRUE) / sum(w, na.rm = TRUE)
  df_infit <- sum(obs_df$Var * w, na.rm = TRUE)
  df_outfit <- sum(w, na.rm = TRUE)
  tibble(
    Infit = infit,
    Outfit = outfit,
    InfitZSTD = zstd_from_mnsq(infit, df_infit, whexact = whexact),
    OutfitZSTD = zstd_from_mnsq(outfit, df_outfit, whexact = whexact),
    DF_Infit = df_infit,
    DF_Outfit = df_outfit
  )
}

calc_facet_fit <- function(obs_df, facet_cols, whexact = FALSE) {
  obs_df <- obs_df |> mutate(.Weight = get_weights(obs_df))
  purrr::map_dfr(facet_cols, function(facet) {
    df <- obs_df |>
      group_by(.data[[facet]]) |>
      summarize(
        Infit = sum(StdSq * Var * .Weight, na.rm = TRUE) / sum(Var * .Weight, na.rm = TRUE),
        Outfit = sum(StdSq * .Weight, na.rm = TRUE) / sum(.Weight, na.rm = TRUE),
        DF_Infit = sum(Var * .Weight, na.rm = TRUE),
        DF_Outfit = sum(.Weight, na.rm = TRUE),
        N = sum(.Weight, na.rm = TRUE),
        .groups = "drop"
      )
    df |>
      mutate(
        InfitZSTD = zstd_from_mnsq(Infit, DF_Infit, whexact = whexact),
        OutfitZSTD = zstd_from_mnsq(Outfit, DF_Outfit, whexact = whexact)
      ) |>
      mutate(Facet = facet, Level = .data[[facet]]) |>
      select(Facet, Level, N, Infit, Outfit, InfitZSTD, OutfitZSTD, DF_Infit, DF_Outfit)
  })
}

calc_facet_se <- function(obs_df, facet_cols) {
  obs_df <- obs_df |> mutate(.Weight = get_weights(obs_df))
  purrr::map_dfr(facet_cols, function(facet) {
    obs_df |>
      group_by(.data[[facet]]) |>
      summarize(
        Info = sum(Var * .Weight, na.rm = TRUE),
        N = sum(.Weight, na.rm = TRUE),
        .groups = "drop"
      ) |>
      mutate(
        Facet = facet,
        Level = .data[[facet]],
        SE = ifelse(Info > 0, 1 / sqrt(Info), NA_real_)
      ) |>
      select(Facet, Level, N, SE)
  })
}

calc_bias_facet <- function(obs_df, facet_cols) {
  obs_df <- obs_df |> mutate(.Weight = get_weights(obs_df))
  purrr::map_dfr(facet_cols, function(facet) {
    obs_df |>
      group_by(.data[[facet]]) |>
      summarize(
        ObservedAverage = weighted_mean(Observed, .Weight),
        ExpectedAverage = weighted_mean(Expected, .Weight),
        MeanResidual = weighted_mean(Residual, .Weight),
        MeanStdResidual = weighted_mean(StdResidual, .Weight),
        MeanAbsStdResidual = weighted_mean(abs(StdResidual), .Weight),
        ChiSq = sum((StdResidual^2) * .Weight, na.rm = TRUE),
        SE_Residual = {
          n_w <- sum(.Weight, na.rm = TRUE)
          ifelse(n_w > 0, sqrt(sum(Var * .Weight, na.rm = TRUE)) / n_w, NA_real_)
        },
        SE_StdResidual = {
          n_w <- sum(.Weight, na.rm = TRUE)
          ifelse(n_w > 0, 1 / sqrt(n_w), NA_real_)
        },
        N = sum(.Weight, na.rm = TRUE),
        .groups = "drop"
      ) |>
      mutate(Facet = facet, Level = .data[[facet]]) |>
      mutate(
        Bias = MeanResidual,
        DF = ifelse(N > 1, N - 1, NA_real_),
        t_Residual = ifelse(is.finite(SE_Residual) & SE_Residual > 0, MeanResidual / SE_Residual, NA_real_),
        t_StdResidual = ifelse(is.finite(SE_StdResidual) & SE_StdResidual > 0, MeanStdResidual / SE_StdResidual, NA_real_),
        p_Residual = ifelse(is.finite(DF) & is.finite(t_Residual), 2 * stats::pt(-abs(t_Residual), df = DF), NA_real_),
        p_StdResidual = ifelse(is.finite(DF) & is.finite(t_StdResidual), 2 * stats::pt(-abs(t_StdResidual), df = DF), NA_real_),
        ChiDf = DF,
        ChiP = ifelse(is.finite(ChiSq) & is.finite(ChiDf) & ChiDf > 0, 1 - stats::pchisq(ChiSq, df = ChiDf), NA_real_)
      ) |>
      select(Facet, Level, N, ObservedAverage, ExpectedAverage, Bias, MeanResidual, MeanStdResidual,
             MeanAbsStdResidual, ChiSq, ChiDf, ChiP, SE_Residual, t_Residual, p_Residual,
             SE_StdResidual, t_StdResidual, p_StdResidual, DF)
  })
}

calc_bias_interactions <- function(obs_df, facet_cols, pairs = NULL, top_n = 20) {
  if (length(facet_cols) < 2) return(tibble())
  obs_df <- obs_df |> mutate(.Weight = get_weights(obs_df))
  if (is.null(pairs)) {
    combos <- combn(facet_cols, 2, simplify = FALSE)
  } else if (length(pairs) == 0) {
    return(tibble())
  } else {
    combos <- pairs
  }
  out <- purrr::map_dfr(combos, function(pair) {
    pair1 <- pair[1]
    pair2 <- pair[2]
    obs_df |>
      group_by(.data[[pair1]], .data[[pair2]]) |>
      summarize(
        ObservedAverage = weighted_mean(Observed, .Weight),
        ExpectedAverage = weighted_mean(Expected, .Weight),
        MeanResidual = weighted_mean(Residual, .Weight),
        MeanStdResidual = weighted_mean(StdResidual, .Weight),
        MeanAbsStdResidual = weighted_mean(abs(StdResidual), .Weight),
        ChiSq = sum((StdResidual^2) * .Weight, na.rm = TRUE),
        SE_Residual = {
          n_w <- sum(.Weight, na.rm = TRUE)
          ifelse(n_w > 0, sqrt(sum(Var * .Weight, na.rm = TRUE)) / n_w, NA_real_)
        },
        SE_StdResidual = {
          n_w <- sum(.Weight, na.rm = TRUE)
          ifelse(n_w > 0, 1 / sqrt(n_w), NA_real_)
        },
        N = sum(.Weight, na.rm = TRUE),
        .groups = "drop"
      ) |>
      mutate(
        Pair = paste(pair1, "x", pair2),
        Level = paste(.data[[pair1]], .data[[pair2]], sep = " | ")
      ) |>
      mutate(
        Bias = MeanResidual,
        DF = ifelse(N > 1, N - 1, NA_real_),
        t_Residual = ifelse(is.finite(SE_Residual) & SE_Residual > 0, MeanResidual / SE_Residual, NA_real_),
        t_StdResidual = ifelse(is.finite(SE_StdResidual) & SE_StdResidual > 0, MeanStdResidual / SE_StdResidual, NA_real_),
        p_Residual = ifelse(is.finite(DF) & is.finite(t_Residual), 2 * stats::pt(-abs(t_Residual), df = DF), NA_real_),
        p_StdResidual = ifelse(is.finite(DF) & is.finite(t_StdResidual), 2 * stats::pt(-abs(t_StdResidual), df = DF), NA_real_),
        ChiDf = DF,
        ChiP = ifelse(is.finite(ChiSq) & is.finite(ChiDf) & ChiDf > 0, 1 - stats::pchisq(ChiSq, df = ChiDf), NA_real_)
      ) |>
      select(Pair, Level, N, ObservedAverage, ExpectedAverage, Bias, MeanResidual, MeanStdResidual,
             MeanAbsStdResidual, ChiSq, ChiDf, ChiP, SE_Residual, t_Residual, p_Residual,
             SE_StdResidual, t_StdResidual, p_StdResidual, DF)
  })
  out |>
    mutate(AbsStd = abs(MeanStdResidual)) |>
    arrange(desc(AbsStd)) |>
    select(-AbsStd) |>
    slice_head(n = top_n)
}

safe_cor <- function(x, y, w = NULL) {
  ok <- is.finite(x) & is.finite(y)
  if (is.null(w)) {
    if (!any(ok)) return(NA_real_)
    x <- x[ok]
    y <- y[ok]
    if (length(unique(x)) < 2 || length(unique(y)) < 2) {
      return(NA_real_)
    }
    return(suppressWarnings(stats::cor(x, y, use = "complete.obs")))
  }
  ok <- ok & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  x <- x[ok]
  y <- y[ok]
  w <- w[ok]
  w_sum <- sum(w)
  if (w_sum <= 0) return(NA_real_)
  mx <- sum(w * x) / w_sum
  my <- sum(w * y) / w_sum
  vx <- sum(w * (x - mx)^2) / w_sum
  vy <- sum(w * (y - my)^2) / w_sum
  if (vx <= 0 || vy <= 0) return(NA_real_)
  cov <- sum(w * (x - mx) * (y - my)) / w_sum
  cov / sqrt(vx * vy)
}

weighted_mean_safe <- function(x, w) {
  weighted_mean(x, w)
}

calc_interrater_agreement <- function(obs_df, facet_cols, rater_facet, res = NULL) {
  if (is.null(obs_df) || nrow(obs_df) == 0) {
    return(list(summary = tibble(), pairs = tibble()))
  }
  if (is.null(rater_facet) || !rater_facet %in% facet_cols) {
    return(list(summary = tibble(), pairs = tibble()))
  }
  context_cols <- setdiff(facet_cols, rater_facet)
  if (length(context_cols) == 0) {
    return(list(summary = tibble(), pairs = tibble()))
  }

  df <- obs_df |>
    mutate(across(all_of(context_cols), as.character)) |>
    tidyr::unite(".context", all_of(context_cols), sep = "|", remove = FALSE) |>
    select(.context, !!rlang::sym(rater_facet), Observed, any_of("Weight")) |>
    mutate(.Weight = get_weights(.))

  df <- df |>
    group_by(.context, !!rlang::sym(rater_facet)) |>
    summarize(Score = weighted_mean(Observed, .Weight), .groups = "drop")

  if (nrow(df) == 0) {
    return(list(summary = tibble(), pairs = tibble()))
  }

  wide <- tryCatch(
    tidyr::pivot_wider(
      df,
      id_cols = .context,
      names_from = !!rlang::sym(rater_facet),
      values_from = Score
    ),
    error = function(e) NULL
  )
  if (is.null(wide)) {
    return(list(summary = tibble(), pairs = tibble()))
  }

  rater_cols <- setdiff(names(wide), ".context")
  if (length(rater_cols) < 2) {
    return(list(summary = tibble(), pairs = tibble()))
  }

  prob_map <- list()
  if (!is.null(res)) {
    probs <- compute_prob_matrix(res)
    if (!is.null(probs) && nrow(probs) == nrow(obs_df)) {
      prob_cols <- paste0(".p", seq_len(ncol(probs)))
      prob_df <- obs_df |>
        mutate(across(all_of(context_cols), as.character)) |>
        tidyr::unite(".context", all_of(context_cols), sep = "|", remove = FALSE) |>
        select(.context, !!rlang::sym(rater_facet), any_of("Weight"))
      prob_df[prob_cols] <- probs
      prob_df <- prob_df |>
        mutate(.Weight = get_weights(.))

      prob_avg <- prob_df |>
        group_by(.context, !!rlang::sym(rater_facet)) |>
        summarize(
          across(all_of(prob_cols), ~ weighted_mean(.x, .Weight)),
          .groups = "drop"
        )

      if (nrow(prob_avg) > 0) {
        for (i in seq_len(nrow(prob_avg))) {
          ctx <- prob_avg$.context[[i]]
          rater_val <- prob_avg[[rater_facet]][[i]]
          key <- paste(ctx, rater_val, sep = "||")
          prob_map[[key]] <- as.numeric(prob_avg[i, prob_cols, drop = TRUE])
        }
      }
    }
  }

  pairs <- combn(rater_cols, 2, simplify = FALSE)
  pair_tbl <- purrr::map_dfr(pairs, function(pair) {
    sub <- wide |>
      select(.context, !!rlang::sym(pair[1]), !!rlang::sym(pair[2])) |>
      filter(is.finite(.data[[pair[1]]]), is.finite(.data[[pair[2]]]))
    n_ok <- nrow(sub)
    if (n_ok == 0) {
      return(tibble(
        Rater1 = pair[1],
        Rater2 = pair[2],
        N = 0,
        Exact = NA_real_,
        ExpectedExact = NA_real_,
        Adjacent = NA_real_,
        MeanDiff = NA_real_,
        MAD = NA_real_,
        Corr = NA_real_
      ))
    }
    v1 <- sub[[pair[1]]]
    v2 <- sub[[pair[2]]]
    diff <- v1 - v2
    exact_count <- sum(diff == 0, na.rm = TRUE)

    exp_vals <- numeric(0)
    if (length(prob_map) > 0) {
      for (ctx in sub$.context) {
        key1 <- paste(ctx, pair[1], sep = "||")
        key2 <- paste(ctx, pair[2], sep = "||")
        p1 <- prob_map[[key1]]
        p2 <- prob_map[[key2]]
        if (is.null(p1) || is.null(p2)) next
        if (any(!is.finite(p1)) || any(!is.finite(p2))) next
        exp_vals <- c(exp_vals, sum(p1 * p2))
      }
    }
    exp_mean <- if (length(exp_vals) > 0) mean(exp_vals) else NA_real_

    tibble(
      Rater1 = pair[1],
      Rater2 = pair[2],
      N = n_ok,
      Exact = exact_count / n_ok,
      ExpectedExact = exp_mean,
      Adjacent = mean(abs(diff) <= 1, na.rm = TRUE),
      MeanDiff = mean(diff, na.rm = TRUE),
      MAD = mean(abs(diff), na.rm = TRUE),
      Corr = safe_cor(v1, v2)
    )
  })

  contexts_with_pairs <- sum(rowSums(!is.na(wide[rater_cols])) >= 2)
  total_pairs <- sum(pair_tbl$N, na.rm = TRUE)
  total_exact <- sum(pair_tbl$Exact * pair_tbl$N, na.rm = TRUE)
  expected_available <- any(is.finite(pair_tbl$ExpectedExact))
  total_expected <- if (expected_available) {
    sum(pair_tbl$ExpectedExact * pair_tbl$N, na.rm = TRUE)
  } else {
    NA_real_
  }
  summary_tbl <- tibble(
    RaterFacet = rater_facet,
    Raters = length(rater_cols),
    Pairs = nrow(pair_tbl),
    Contexts = contexts_with_pairs,
    TotalPairs = total_pairs,
    ExactAgreements = total_exact,
    ExpectedAgreements = ifelse(expected_available, total_expected, NA_real_),
    ExactAgreement = ifelse(total_pairs > 0, total_exact / total_pairs, NA_real_),
    ExpectedExactAgreement = ifelse(expected_available && total_pairs > 0, total_expected / total_pairs, NA_real_),
    AdjacentAgreement = weighted_mean_safe(pair_tbl$Adjacent, pair_tbl$N),
    MeanAbsDiff = weighted_mean_safe(pair_tbl$MAD, pair_tbl$N),
    MeanCorr = weighted_mean_safe(pair_tbl$Corr, pair_tbl$N)
  )

  list(summary = summary_tbl, pairs = pair_tbl)
}

calc_ptmea <- function(obs_df, facet_cols) {
  facet_cols <- setdiff(facet_cols, "Person")
  if (length(facet_cols) == 0) return(tibble())
  obs_df <- obs_df |> mutate(.Weight = get_weights(obs_df))
  purrr::map_dfr(facet_cols, function(facet) {
    obs_df |>
      group_by(.data[[facet]]) |>
      summarize(
        PTMEA = safe_cor(Observed, PersonMeasure, w = .Weight),
        N = sum(.Weight, na.rm = TRUE),
        .groups = "drop"
      ) |>
      mutate(Facet = facet, Level = .data[[facet]]) |>
      select(Facet, Level, PTMEA, N)
  })
}

expected_score_from_eta <- function(eta, step_cum, rating_min) {
  if (!is.finite(eta) || length(step_cum) == 0) return(NA_real_)
  probs <- category_prob_rsm(eta, step_cum)
  k_vals <- 0:(length(step_cum) - 1)
  rating_min + sum(probs * k_vals)
}

estimate_eta_from_target <- function(target, step_cum, rating_min, rating_max) {
  if (!is.finite(target) || length(step_cum) == 0) return(NA_real_)
  if (target <= rating_min) return(-Inf)
  if (target >= rating_max) return(Inf)
  f <- function(eta) expected_score_from_eta(eta, step_cum, rating_min) - target
  lower <- -10
  upper <- 10
  f_low <- f(lower)
  f_up <- f(upper)
  if (!is.finite(f_low) || !is.finite(f_up)) return(NA_real_)
  if (f_low * f_up > 0) {
    lower <- -20
    upper <- 20
    f_low <- f(lower)
    f_up <- f(upper)
    if (!is.finite(f_low) || !is.finite(f_up) || f_low * f_up > 0) return(NA_real_)
  }
  uniroot(f, lower = lower, upper = upper)$root
}

facet_anchor_status <- function(facet, levels, config) {
  spec <- if (facet == "Person") config$theta_spec else config$facet_specs[[facet]]
  if (is.null(spec)) return(rep("", length(levels)))
  anchors <- spec$anchors
  groups <- spec$groups
  status <- rep("", length(levels))
  if (!is.null(anchors)) {
    anchor_vals <- anchors[match(levels, names(anchors))]
    status[is.finite(anchor_vals)] <- "A"
  }
  if (!is.null(groups)) {
    group_vals <- groups[match(levels, names(groups))]
    status[status == "" & !is.na(group_vals) & group_vals != ""] <- "G"
  }
  status
}

calc_facets_report_tbls <- function(res,
                                    diagnostics,
                                    totalscore = TRUE,
                                    umean = 0,
                                    uscale = 1,
                                    udecimals = 2,
                                    omit_unobserved = FALSE,
                                    xtreme = 0) {
  if (is.null(res) || is.null(diagnostics)) return(list())
  obs_df <- diagnostics$obs
  measures <- diagnostics$measures
  if (nrow(obs_df) == 0 || nrow(measures) == 0) return(list())

  prep <- res$prep
  config <- res$config
  rating_min <- prep$rating_min
  rating_max <- prep$rating_max
  sizes <- build_param_sizes(config)
  params <- expand_params(res$opt$par, sizes, config)
  theta_hat <- if (config$method == "JMLE") {
    params$theta
  } else {
    res$facets$person$Estimate
  }
  theta_mean <- if (length(theta_hat) > 0) mean(theta_hat, na.rm = TRUE) else 0
  facet_means <- purrr::map_dbl(config$facet_names, function(f) {
    mean(params$facets[[f]], na.rm = TRUE)
  })
  names(facet_means) <- config$facet_names
  facet_signs <- config$facet_signs
  if (is.null(facet_signs) || length(facet_signs) == 0) {
    facet_signs <- setNames(rep(-1, length(config$facet_names)), config$facet_names)
  }

  if (config$model == "RSM") {
    step_cum_common <- c(0, cumsum(params$steps))
    step_cum_mean <- step_cum_common
  } else {
    step_mat <- params$steps_mat
    if (is.null(step_mat) || length(step_mat) == 0) {
      step_cum_common <- numeric(0)
      step_cum_mean <- numeric(0)
    } else {
      step_mean <- colMeans(step_mat, na.rm = TRUE)
      step_cum_common <- t(apply(step_mat, 1, function(x) c(0, cumsum(x))))
      step_cum_mean <- c(0, cumsum(step_mean))
    }
  }

  facet_names <- c("Person", config$facet_names)
  facet_levels_all <- lapply(facet_names, function(facet) {
    if (facet == "Person") {
      prep$levels$Person
    } else {
      prep$levels[[facet]]
    }
  })
  names(facet_levels_all) <- facet_names

  extreme_levels <- purrr::map(facet_names, function(facet) {
    if (!facet %in% names(obs_df)) return(character(0))
    obs_df |>
      group_by(.data[[facet]]) |>
      summarize(
        MinScore = min(Observed, na.rm = TRUE),
        MaxScore = max(Observed, na.rm = TRUE),
        .groups = "drop"
      ) |>
      filter(
        (MinScore == prep$rating_min & MaxScore == prep$rating_min) |
          (MinScore == prep$rating_max & MaxScore == prep$rating_max)
      ) |>
      mutate(Level = as.character(.data[[facet]])) |>
      pull(Level)
  })
  names(extreme_levels) <- facet_names

  extreme_flags <- list()
  for (facet in facet_names) {
    if (facet %in% names(obs_df)) {
      extreme_flags[[facet]] <- obs_df[[facet]] %in% extreme_levels[[facet]]
    }
  }
  if (length(extreme_flags) > 0) {
    extreme_flag_df <- as_tibble(extreme_flags)
    extreme_count <- rowSums(extreme_flag_df, na.rm = TRUE)
  } else {
    extreme_count <- rep(0, nrow(obs_df))
  }
  out <- purrr::map(facet_names, function(facet) {
    if (!facet %in% names(obs_df)) return(tibble())
    status_tbl <- obs_df |>
      group_by(.data[[facet]]) |>
      summarize(
        MinScore = min(Observed, na.rm = TRUE),
        MaxScore = max(Observed, na.rm = TRUE),
        TotalCountAll = n(),
        .groups = "drop"
      ) |>
      mutate(Level = as.character(.data[[facet]])) |>
      select(Level, MinScore, MaxScore, TotalCountAll)

    if (isTRUE(totalscore)) {
      score_source <- obs_df
    } else {
      flag <- extreme_flags[[facet]]
      if (is.null(flag)) {
        flag <- rep(FALSE, nrow(obs_df))
      }
      active_mask <- (extreme_count == 0) | ((extreme_count == 1) & flag)
      score_source <- obs_df[active_mask, , drop = FALSE]
    }
    score_tbl <- score_source |>
      group_by(.data[[facet]]) |>
      summarize(
        TotalScore = sum(Observed, na.rm = TRUE),
        TotalCount = n(),
        WeightdScore = sum(Observed * Weight, na.rm = TRUE),
        WeightdCount = sum(Weight, na.rm = TRUE),
        ObservedAverage = ifelse(sum(Weight, na.rm = TRUE) > 0,
                                 sum(Observed * Weight, na.rm = TRUE) / sum(Weight, na.rm = TRUE),
                                 NA_real_),
        .groups = "drop"
      ) |>
      mutate(Level = as.character(.data[[facet]])) |>
      select(Level, TotalScore, TotalCount, WeightdScore, WeightdCount, ObservedAverage)

    level_tbl <- tibble(Level = as.character(facet_levels_all[[facet]]))
    score_tbl <- level_tbl |>
      left_join(score_tbl, by = "Level") |>
      mutate(
        TotalScore = replace_na(TotalScore, 0),
        TotalCount = replace_na(TotalCount, 0),
        WeightdScore = replace_na(WeightdScore, 0),
        WeightdCount = replace_na(WeightdCount, 0),
        ObservedAverage = ifelse(WeightdCount > 0, ObservedAverage, NA_real_)
      )

    meas_tbl <- measures |>
      filter(Facet == facet) |>
      mutate(Level = as.character(Level)) |>
      select(Level, Estimate, SE, Infit, Outfit, InfitZSTD, OutfitZSTD, PTMEA)

    tbl <- level_tbl |>
      left_join(score_tbl, by = "Level") |>
      left_join(status_tbl, by = "Level") |>
      left_join(meas_tbl, by = "Level")

    anchor_status <- facet_anchor_status(facet, tbl$Level, config)
    status <- dplyr::case_when(
      tbl$TotalCountAll == 0 ~ "No data",
      tbl$MinScore == tbl$MaxScore & tbl$MinScore == rating_min ~ "Minimum",
      tbl$MinScore == tbl$MaxScore & tbl$MaxScore == rating_max ~ "Maximum",
      tbl$TotalCountAll == 1 ~ "One datum",
      TRUE ~ ""
    )

    sign_vec <- facet_signs[names(facet_means)]
    if (facet == "Person") {
      other_sum <- sum(sign_vec * facet_means, na.rm = TRUE)
      eta_m <- tbl$Estimate + other_sum
      eta_z <- tbl$Estimate
    } else {
      sign <- if (!is.null(facet_signs[[facet]])) facet_signs[[facet]] else -1
      other_sum <- sum(sign_vec[names(facet_means) != facet] *
                         facet_means[names(facet_means) != facet], na.rm = TRUE)
      eta_m <- theta_mean + other_sum + sign * tbl$Estimate
      eta_z <- sign * tbl$Estimate
    }

    if (config$model == "PCM" && !is.null(config$step_facet)) {
      step_levels <- prep$levels[[config$step_facet]]
      if (facet == config$step_facet && length(step_levels) > 0 && length(step_cum_common) > 0) {
        step_cum_list <- purrr::map(tbl$Level, function(lvl) {
          idx <- match(lvl, step_levels)
          if (is.na(idx) || idx < 1 || idx > nrow(step_cum_common)) {
            step_cum_mean
          } else {
            step_cum_common[idx, ]
          }
        })
      } else {
        step_cum_list <- rep(list(step_cum_mean), nrow(tbl))
      }
    } else {
      step_cum_list <- rep(list(step_cum_common), nrow(tbl))
    }

    fair_m <- purrr::map2_dbl(eta_m, step_cum_list, ~ expected_score_from_eta(.x, .y, rating_min))
    fair_z <- purrr::map2_dbl(eta_z, step_cum_list, ~ expected_score_from_eta(.x, .y, rating_min))

    xtreme_target <- ifelse(
      status == "Minimum", rating_min + xtreme,
      ifelse(status == "Maximum", rating_max - xtreme, NA_real_)
    )
    xtreme_eta <- purrr::map2_dbl(xtreme_target, step_cum_list, ~ {
      if (!is.finite(.x) || xtreme <= 0) return(NA_real_)
      estimate_eta_from_target(.x, .y, rating_min, rating_max)
    })
    measure_logit <- tbl$Estimate
    if (any(is.finite(xtreme_eta))) {
      if (facet == "Person") {
        measure_logit <- ifelse(
          is.finite(xtreme_eta),
          xtreme_eta - other_sum,
          measure_logit
        )
      } else {
        sign <- if (!is.null(facet_signs[[facet]])) facet_signs[[facet]] else -1
        measure_logit <- ifelse(
          is.finite(xtreme_eta),
          (xtreme_eta - theta_mean - other_sum) / sign,
          measure_logit
        )
      }
    }

    scale_factor <- ifelse(is.finite(uscale), uscale, 1)
    scale_origin <- ifelse(is.finite(umean), umean, 0)

    tbl <- tbl |>
      mutate(
        Anchor = anchor_status,
        Status = status,
        FairM = fair_m,
        FairZ = fair_z,
        Measure = ifelse(is.finite(measure_logit), measure_logit * scale_factor + scale_origin, NA_real_),
        ModelSE = ifelse(is.finite(SE), abs(scale_factor) * SE, NA_real_),
        RealSE = ifelse(is.finite(SE) & is.finite(Infit), abs(scale_factor) * SE * sqrt(pmax(Infit, 0)), NA_real_),
        Infit = ifelse(Status %in% c("Minimum", "Maximum"), NA_real_, Infit),
        Outfit = ifelse(Status %in% c("Minimum", "Maximum"), NA_real_, Outfit),
        InfitZSTD = ifelse(Status %in% c("Minimum", "Maximum"), NA_real_, InfitZSTD),
        OutfitZSTD = ifelse(Status %in% c("Minimum", "Maximum"), NA_real_, OutfitZSTD)
      ) |>
      transmute(
        TotalScore,
        TotalCount,
        WeightdScore,
        WeightdCount,
        ObservedAverage,
        FairM,
        FairZ,
        Measure,
        ModelSE,
        RealSE,
        InfitMnSq = Infit,
        InfitZStd = InfitZSTD,
        OutfitMnSq = Outfit,
        OutfitZStd = OutfitZSTD,
        PtMeaCorr = ifelse(facet == "Person", NA_real_, PTMEA),
        Anchor,
        Status,
        Level,
        TotalCountAll
      ) |>
      arrange(desc(Measure), desc(TotalCount))

    if (isTRUE(omit_unobserved)) {
      tbl <- tbl |> filter(TotalCountAll > 0)
    }

    tbl |>
      select(-TotalCountAll)
  })

  names(out) <- facet_names
  out
}

format_facets_report_gt <- function(tbl, facet, decimals = 2, totalscore = TRUE) {
  if (is.null(tbl) || nrow(tbl) == 0) {
    return(gt(tibble(Message = "No facet report available.")))
  }
  gt(tbl) |>
    cols_label(
      TotalScore = if (isTRUE(totalscore)) "Total Score" else "Obsvd Score",
      TotalCount = if (isTRUE(totalscore)) "Total Count" else "Obsvd Count",
      WeightdScore = "Weightd Score",
      WeightdCount = "Weightd Count",
      ObservedAverage = "Obsvd Average",
      FairM = "Fair(M) Average",
      FairZ = "Fair(Z) Average",
      Measure = "Measure",
      ModelSE = "Model S.E.",
      RealSE = "Real S.E.",
      InfitMnSq = "Infit MnSq",
      InfitZStd = "Infit ZStd",
      OutfitMnSq = "Outfit MnSq",
      OutfitZStd = "Outfit ZStd",
      PtMeaCorr = "PtMea Corr",
      Anchor = "Anch",
      Status = "Status",
      Level = "Element"
    ) |>
    cols_move_to_end(columns = Level) |>
    fmt_number(
      columns = c(TotalScore, TotalCount, WeightdScore, WeightdCount),
      decimals = 0
    ) |>
    fmt_number(
      columns = c(ObservedAverage, FairM, FairZ, Measure, ModelSE, RealSE,
                  InfitMnSq, InfitZStd, OutfitMnSq, OutfitZStd, PtMeaCorr),
      decimals = decimals
    ) |>
    tab_style(
      style = list(cell_fill(color = "#fff3cd")),
      locations = cells_body(rows = Status %in% c("Minimum", "Maximum", "One datum"))
    )
}

calc_expected_category_counts <- function(res) {
  if (is.null(res)) return(tibble())
  prep <- res$prep
  config <- res$config
  idx <- build_indices(prep, step_facet = config$step_facet)
  sizes <- build_param_sizes(config)
  params <- expand_params(res$opt$par, sizes, config)
  theta_hat <- if (config$method == "JMLE") {
    params$theta
  } else {
    res$facets$person$Estimate
  }
  eta <- compute_eta(idx, params, config, theta_override = theta_hat)
  if (config$model == "RSM") {
    step_cum <- c(0, cumsum(params$steps))
    probs <- category_prob_rsm(eta, step_cum)
  } else {
    step_cum_mat <- t(apply(params$steps_mat, 1, function(x) c(0, cumsum(x))))
    probs <- category_prob_pcm(eta, step_cum_mat, idx$step_idx)
  }
  if (length(probs) == 0) return(tibble())
  w <- idx$weight
  exp_counts <- if (is.null(w)) {
    colSums(probs, na.rm = TRUE)
  } else {
    colSums(probs * w, na.rm = TRUE)
  }
  total_exp <- sum(exp_counts)
  cat_vals <- seq(prep$rating_min, prep$rating_max)
  tibble(
    Category = cat_vals,
    ExpectedCount = exp_counts,
    ExpectedPercent = ifelse(total_exp > 0, 100 * exp_counts / total_exp, NA_real_)
  )
}

calc_category_stats <- function(obs_df, res = NULL, whexact = FALSE) {
  if (nrow(obs_df) == 0) return(tibble())
  obs_df <- obs_df |> mutate(.Weight = get_weights(obs_df))
  total_n <- sum(obs_df$.Weight, na.rm = TRUE)
  obs_summary <- obs_df |>
    group_by(Observed) |>
    summarize(
      Count = sum(.Weight, na.rm = TRUE),
      AvgPersonMeasure = weighted_mean(PersonMeasure, .Weight),
      ExpectedAverage = weighted_mean(Expected, .Weight),
      Infit = sum(StdSq * Var * .Weight, na.rm = TRUE) / sum(Var * .Weight, na.rm = TRUE),
      Outfit = sum(StdSq * .Weight, na.rm = TRUE) / sum(.Weight, na.rm = TRUE),
      MeanResidual = weighted_mean(Residual, .Weight),
      DF_Infit = sum(Var * .Weight, na.rm = TRUE),
      DF_Outfit = sum(.Weight, na.rm = TRUE),
      .groups = "drop"
    ) |>
    rename(Category = Observed)

  all_categories <- if (!is.null(res)) {
    seq(res$prep$rating_min, res$prep$rating_max)
  } else {
    sort(unique(obs_df$Observed))
  }

  cat_tbl <- tibble(Category = all_categories) |>
    left_join(obs_summary, by = "Category") |>
    mutate(
      Count = replace_na(Count, 0),
      Percent = ifelse(total_n > 0, 100 * Count / total_n, NA_real_),
      InfitZSTD = zstd_from_mnsq(Infit, DF_Infit, whexact = whexact),
      OutfitZSTD = zstd_from_mnsq(Outfit, DF_Outfit, whexact = whexact)
    )

  exp_tbl <- calc_expected_category_counts(res)
  if (nrow(exp_tbl) > 0) {
    cat_tbl <- cat_tbl |>
      left_join(exp_tbl, by = "Category") |>
      mutate(
        DiffCount = ifelse(is.finite(ExpectedCount), Count - ExpectedCount, NA_real_),
        DiffPercent = ifelse(is.finite(ExpectedPercent), Percent - ExpectedPercent, NA_real_)
      )
  }

  cat_tbl |>
    mutate(
      LowCount = Count < 10,
      InfitFlag = !is.na(Infit) & (Infit < 0.5 | Infit > 1.5),
      OutfitFlag = !is.na(Outfit) & (Outfit < 0.5 | Outfit > 1.5),
      ZSTDFlag = (is.finite(InfitZSTD) & abs(InfitZSTD) >= 2) |
        (is.finite(OutfitZSTD) & abs(OutfitZSTD) >= 2)
    ) |>
    arrange(Category)
}

make_union_find <- function(nodes) {
  parent <- setNames(nodes, nodes)
  find_root <- function(x) {
    px <- parent[[x]]
    if (is.null(px)) return(NA_character_)
    if (px != x) {
      parent[[x]] <<- find_root(px)
    }
    parent[[x]]
  }
  union_nodes <- function(a, b) {
    ra <- find_root(a)
    rb <- find_root(b)
    if (is.na(ra) || is.na(rb) || ra == rb) return(NULL)
    parent[[rb]] <<- ra
  }
  list(find = find_root, union = union_nodes)
}

calc_subsets <- function(obs_df, facet_cols) {
  if (is.null(obs_df) || nrow(obs_df) == 0 || length(facet_cols) == 0) {
    return(list(summary = tibble(), nodes = tibble()))
  }
  df <- obs_df |>
    select(all_of(facet_cols)) |>
    filter(if_all(everything(), ~ !is.na(.)))
  if (nrow(df) == 0) {
    return(list(summary = tibble(), nodes = tibble()))
  }

  nodes <- unique(unlist(lapply(facet_cols, function(facet) {
    paste0(facet, ":", as.character(df[[facet]]))
  })))
  if (length(nodes) == 0) {
    return(list(summary = tibble(), nodes = tibble()))
  }
  uf <- make_union_find(nodes)

  for (i in seq_len(nrow(df))) {
    row_vals <- unlist(df[i, facet_cols, drop = FALSE], use.names = FALSE)
    row_nodes <- paste0(facet_cols, ":", as.character(row_vals))
    row_nodes <- row_nodes[!is.na(row_nodes)]
    if (length(row_nodes) < 2) next
    base <- row_nodes[1]
    for (node in row_nodes[-1]) {
      uf$union(base, node)
    }
  }

  comp_ids <- vapply(nodes, uf$find, character(1))
  comp_levels <- unique(comp_ids)
  comp_index <- setNames(seq_along(comp_levels), comp_levels)

  node_tbl <- tibble(
    Node = nodes,
    Component = comp_ids,
    Subset = comp_index[comp_ids],
    Facet = stringr::str_extract(nodes, "^[^:]+"),
    Level = stringr::str_replace(nodes, "^[^:]+:", "")
  )

  facet_counts <- node_tbl |>
    group_by(Subset, Facet) |>
    summarize(Levels = n_distinct(Level), .groups = "drop") |>
    tidyr::pivot_wider(names_from = Facet, values_from = Levels, values_fill = 0)

  row_subset <- vapply(seq_len(nrow(df)), function(i) {
    node <- paste0(facet_cols[1], ":", as.character(df[[facet_cols[1]]][i]))
    comp_index[uf$find(node)]
  }, integer(1))
  obs_counts <- tibble(Subset = row_subset) |>
    count(Subset, name = "Observations")

  summary_tbl <- facet_counts |>
    left_join(obs_counts, by = "Subset") |>
    mutate(Observations = replace_na(Observations, 0)) |>
    arrange(desc(Observations))

  list(summary = summary_tbl, nodes = node_tbl)
}

calc_step_order <- function(step_tbl) {
  if (is.null(step_tbl) || nrow(step_tbl) == 0) return(tibble())
  step_tbl <- step_tbl |>
    mutate(StepIndex = suppressWarnings(as.integer(stringr::str_extract(Step, "\\d+"))))
  if (!"StepFacet" %in% names(step_tbl)) {
    step_tbl <- mutate(step_tbl, StepFacet = "Common")
  }
  step_tbl |>
    arrange(StepFacet, StepIndex) |>
    group_by(StepFacet) |>
    mutate(
      Spacing = Estimate - lag(Estimate),
      Ordered = ifelse(is.na(Spacing), NA, Spacing > 0)
    ) |>
    ungroup()
}

category_warnings_text <- function(cat_tbl, step_tbl = NULL) {
  if (is.null(cat_tbl) || nrow(cat_tbl) == 0) return("No category diagnostics available.")
  msgs <- character()
  unused <- cat_tbl |>
    filter(Count == 0)
  if (nrow(unused) > 0) {
    msgs <- c(msgs, paste0("Unused categories: ", paste(unused$Category, collapse = ", ")))
  }
  low_counts <- cat_tbl |>
    filter(Count < 10)
  if (nrow(low_counts) > 0) {
    msgs <- c(msgs, paste0("Low category counts (<10): ", paste(low_counts$Category, collapse = ", ")))
  }
  if ("DiffPercent" %in% names(cat_tbl)) {
    diff_bad <- cat_tbl |>
      filter(is.finite(DiffPercent), abs(DiffPercent) >= 5)
    if (nrow(diff_bad) > 0) {
      msgs <- c(msgs, paste0("Observed vs expected % differs by >= 5: ", paste(diff_bad$Category, collapse = ", ")))
    }
  }
  if ("InfitZSTD" %in% names(cat_tbl) && "OutfitZSTD" %in% names(cat_tbl)) {
    zstd_bad <- cat_tbl |>
      filter((is.finite(InfitZSTD) & abs(InfitZSTD) >= 2) | (is.finite(OutfitZSTD) & abs(OutfitZSTD) >= 2))
    if (nrow(zstd_bad) > 0) {
      msgs <- c(msgs, paste0("Category |ZSTD| >= 2: ", paste(zstd_bad$Category, collapse = ", ")))
    }
  }
  avg_tbl <- cat_tbl |>
    filter(is.finite(AvgPersonMeasure)) |>
    arrange(Category)
  if (nrow(avg_tbl) >= 3 && is.unsorted(avg_tbl$AvgPersonMeasure, strictly = FALSE)) {
    msgs <- c(msgs, "Category averages are not monotonic (Avg Measure by category).")
  }
  if (!is.null(step_tbl) && nrow(step_tbl) > 0) {
    disordered <- step_tbl |>
      filter(!is.na(Ordered), Ordered == FALSE)
    if (nrow(disordered) > 0) {
      bad <- disordered |>
        mutate(Label = paste0(StepFacet, ":", Step)) |>
        pull(Label)
      msgs <- c(msgs, paste0("Disordered thresholds detected: ", paste(bad, collapse = ", ")))
    }
  }
  if (length(msgs) == 0) {
    "No major category warnings detected."
  } else {
    paste(msgs, collapse = "\n")
  }
}

get_extreme_levels <- function(obs_df, facet_names, rating_min, rating_max) {
  out <- list()
  for (facet in facet_names) {
    if (!facet %in% names(obs_df)) {
      out[[facet]] <- character(0)
      next
    }
    stat <- obs_df |>
      group_by(.data[[facet]]) |>
      summarize(
        MinScore = min(Observed, na.rm = TRUE),
        MaxScore = max(Observed, na.rm = TRUE),
        .groups = "drop"
      )
    extreme <- stat |>
      filter(
        (MinScore == rating_min & MaxScore == rating_min) |
          (MinScore == rating_max & MaxScore == rating_max)
      )
    out[[facet]] <- as.character(extreme[[facet]])
  }
  out
}

estimate_bias_interaction <- function(res,
                                      diagnostics,
                                      facet_a,
                                      facet_b,
                                      max_abs = 10,
                                      omit_extreme = TRUE,
                                      max_iter = 4,
                                      tol = 1e-3) {
  if (is.null(res) || is.null(diagnostics)) return(list())
  obs_df <- diagnostics$obs
  if (is.null(obs_df) || nrow(obs_df) == 0) return(list())
  if (facet_a == facet_b) return(list())

  facet_names <- c("Person", res$config$facet_names)
  if (!facet_a %in% facet_names || !facet_b %in% facet_names) return(list())

  prep <- res$prep
  config <- res$config
  sizes <- build_param_sizes(config)
  params <- expand_params(res$opt$par, sizes, config)
  idx <- build_indices(prep, step_facet = config$step_facet)
  theta_hat <- if (config$method == "JMLE") params$theta else res$facets$person$Estimate
  eta_base <- compute_eta(idx, params, config, theta_override = theta_hat)
  score_k <- idx$score_k
  weight <- idx$weight
  step_idx <- idx$step_idx

  if (config$model == "RSM") {
    step_cum <- c(0, cumsum(params$steps))
    step_cum_mat <- NULL
  } else {
    step_cum <- NULL
    step_cum_mat <- t(apply(params$steps_mat, 1, function(x) c(0, cumsum(x))))
  }

  if (omit_extreme) {
    extreme_levels <- get_extreme_levels(
      obs_df,
      c(facet_a, facet_b),
      prep$rating_min,
      prep$rating_max
    )
  } else {
    extreme_levels <- list()
  }

  measures <- diagnostics$measures
  meas_map <- list()
  se_map <- list()
  if (!is.null(measures) && nrow(measures) > 0) {
    for (i in seq_len(nrow(measures))) {
      key <- paste(measures$Facet[i], as.character(measures$Level[i]), sep = "||")
      meas_map[[key]] <- measures$Estimate[i]
      se_map[[key]] <- measures$SE[i]
    }
  }

  level_map <- lapply(facet_names, function(facet) {
    if (facet == "Person") {
      as.character(prep$levels$Person)
    } else {
      as.character(prep$levels[[facet]])
    }
  })
  names(level_map) <- facet_names

  group_keys <- interaction(obs_df[[facet_a]], obs_df[[facet_b]], drop = TRUE, sep = "||")
  group_indices <- split(seq_len(nrow(obs_df)), group_keys)
  groups <- list()
  for (key in names(group_indices)) {
    lvl_vals <- strsplit(key, "\\|\\|")[[1]]
    lvl_a <- as.character(lvl_vals[1])
    lvl_b <- as.character(lvl_vals[2])
    if (isTRUE(omit_extreme)) {
      if (lvl_a %in% extreme_levels[[facet_a]] || lvl_b %in% extreme_levels[[facet_b]]) {
        next
      }
    }
    idx_rows <- group_indices[[key]]
    if (length(idx_rows) == 0) next
    groups[[key]] <- list(level_a = lvl_a, level_b = lvl_b, idx = idx_rows)
  }
  if (length(groups) == 0) return(list())

  estimate_bias_for_group <- function(idx_rows) {
    eta_sub <- eta_base[idx_rows]
    score_k_sub <- score_k[idx_rows]
    weight_sub <- if (!is.null(weight)) weight[idx_rows] else NULL
    step_idx_sub <- if (!is.null(step_idx)) step_idx[idx_rows] else NULL

    if (config$model == "RSM") {
      nll <- function(b) -loglik_rsm(eta_sub + b, score_k_sub, step_cum, weight = weight_sub)
    } else {
      nll <- function(b) -loglik_pcm(eta_sub + b, score_k_sub, step_cum_mat, step_idx_sub, weight = weight_sub)
    }
    opt <- tryCatch(stats::optimize(nll, interval = c(-max_abs, max_abs)), error = function(e) NULL)
    if (is.null(opt)) return(NA_real_)
    opt$minimum
  }

  iteration_metrics <- function(bias_map) {
    max_resid <- 0
    max_resid_pct <- NA_real_
    for (g in groups) {
      idx_rows <- g$idx
      key <- paste(g$level_a, g$level_b, sep = "||")
      bias <- bias_map[[key]]
      eta_sub <- eta_base[idx_rows] + ifelse(is.finite(bias), bias, 0)
      score_k_sub <- score_k[idx_rows]
      step_idx_sub <- if (!is.null(step_idx)) step_idx[idx_rows] else NULL
      probs <- if (config$model == "RSM") {
        category_prob_rsm(eta_sub, step_cum)
      } else {
        category_prob_pcm(eta_sub, step_cum_mat, step_idx_sub)
      }
      k_vals <- 0:(ncol(probs) - 1)
      expected_k <- as.vector(probs %*% k_vals)
      expected_score <- prep$rating_min + expected_k
      obs_score <- obs_df$Observed[idx_rows]
      if ("Weight" %in% names(obs_df)) {
        w <- obs_df$Weight[idx_rows]
        obs_score <- obs_score * w
        expected_score <- expected_score * w
      }
      obs_sum <- sum(obs_score, na.rm = TRUE)
      exp_sum <- sum(expected_score, na.rm = TRUE)
      resid_sum <- obs_sum - exp_sum
      if (abs(resid_sum) >= abs(max_resid)) {
        max_resid <- resid_sum
        max_resid_pct <- ifelse(exp_sum != 0, resid_sum / exp_sum * 100, NA_real_)
      }
    }
    list(
      max_resid = max_resid,
      max_resid_pct = max_resid_pct,
      max_resid_categories = NA_real_
    )
  }

  bias_map <- list()
  for (key in names(groups)) {
    bias_map[[key]] <- 0
  }

  iter_rows <- list()
  for (it in seq_len(max_iter)) {
    max_change_abs <- 0
    max_change_signed <- 0
    changes <- numeric(0)
    for (key in names(groups)) {
      g <- groups[[key]]
      bias_hat <- estimate_bias_for_group(g$idx)
      prev <- bias_map[[key]]
      if (is.finite(bias_hat) && is.finite(prev)) {
        delta <- bias_hat - prev
        if (abs(delta) >= max_change_abs) {
          max_change_abs <- abs(delta)
          max_change_signed <- delta
        }
        changes <- c(changes, abs(delta))
      } else {
        changes <- c(changes, NA_real_)
      }
      bias_map[[key]] <- bias_hat
    }
    resid_info <- iteration_metrics(bias_map)
    iter_rows[[it]] <- tibble(
      Iteration = it,
      MaxScoreResidual = resid_info$max_resid,
      MaxScoreResidualPct = resid_info$max_resid_pct,
      MaxScoreResidualCategories = resid_info$max_resid_categories,
      MaxLogitChange = max_change_signed,
      BiasCells = sum(changes > tol, na.rm = TRUE)
    )
    if (max_change_abs < tol) break
  }

  rows <- list()
  seq_id <- 1
  for (key in names(groups)) {
    g <- groups[[key]]
    lvl_a <- g$level_a
    lvl_b <- g$level_b
    idx_rows <- g$idx
    bias_hat <- bias_map[[key]]
    eta_sub <- eta_base[idx_rows]
    score_k_sub <- score_k[idx_rows]
    weight_sub <- if (!is.null(weight)) weight[idx_rows] else rep(1, length(idx_rows))
    step_idx_sub <- if (!is.null(step_idx)) step_idx[idx_rows] else NULL

    probs <- if (config$model == "RSM") {
      category_prob_rsm(eta_sub + ifelse(is.finite(bias_hat), bias_hat, 0), step_cum)
    } else {
      category_prob_pcm(eta_sub + ifelse(is.finite(bias_hat), bias_hat, 0), step_cum_mat, step_idx_sub)
    }
    k_vals <- 0:(ncol(probs) - 1)
    expected_k <- as.vector(probs %*% k_vals)
    var_k <- as.vector(probs %*% (k_vals^2)) - expected_k^2
    var_k <- ifelse(var_k <= 1e-10, NA_real_, var_k)
    resid_k <- score_k_sub - expected_k
    std_sq <- resid_k^2 / var_k

    info <- sum(var_k * weight_sub, na.rm = TRUE)
    se <- ifelse(is.finite(info) && info > 0, 1 / sqrt(info), NA_real_)
    infit <- ifelse(sum(var_k * weight_sub, na.rm = TRUE) > 0,
                    sum(std_sq * var_k * weight_sub, na.rm = TRUE) / sum(var_k * weight_sub, na.rm = TRUE),
                    NA_real_)
    outfit <- ifelse(sum(weight_sub, na.rm = TRUE) > 0,
                     sum(std_sq * weight_sub, na.rm = TRUE) / sum(weight_sub, na.rm = TRUE),
                     NA_real_)

    obs_slice <- obs_df[idx_rows, , drop = FALSE]
    w_obs <- if ("Weight" %in% names(obs_slice)) obs_slice$Weight else rep(1, nrow(obs_slice))
    obs_score <- sum(obs_slice$Observed * w_obs, na.rm = TRUE)
    exp_score <- sum((prep$rating_min + expected_k) * w_obs, na.rm = TRUE)
    obs_count <- sum(w_obs, na.rm = TRUE)
    obs_exp_avg <- ifelse(obs_count > 0, (obs_score - exp_score) / obs_count, NA_real_)

    n_obs <- nrow(obs_slice)
    df_t <- max(n_obs - 1, 0)
    t_val <- ifelse(is.finite(bias_hat) && is.finite(se) && se > 0, bias_hat / se, NA_real_)
    p_val <- ifelse(is.finite(t_val) && df_t > 0, 2 * stats::pt(-abs(t_val), df = df_t), NA_real_)

    idx_a <- match(lvl_a, level_map[[facet_a]])
    idx_b <- match(lvl_b, level_map[[facet_b]])
    meas_a <- meas_map[[paste(facet_a, lvl_a, sep = "||")]]
    meas_b <- meas_map[[paste(facet_b, lvl_b, sep = "||")]]
    se_a <- se_map[[paste(facet_a, lvl_a, sep = "||")]]
    se_b <- se_map[[paste(facet_b, lvl_b, sep = "||")]]

    rows[[seq_id]] <- tibble(
      Sq = seq_id,
      `Observd Score` = obs_score,
      `Expctd Score` = exp_score,
      `Observd Count` = obs_count,
      `Obs-Exp Average` = obs_exp_avg,
      `Bias Size` = bias_hat,
      `S.E.` = se,
      t = t_val,
      `d.f.` = df_t,
      `Prob.` = p_val,
      Infit = infit,
      Outfit = outfit,
      ObsN = n_obs,
      FacetA = facet_a,
      FacetA_Level = lvl_a,
      FacetA_Index = ifelse(is.na(idx_a), NA_real_, idx_a),
      FacetA_Measure = meas_a,
      FacetA_SE = se_a,
      FacetB = facet_b,
      FacetB_Level = lvl_b,
      FacetB_Index = ifelse(is.na(idx_b), NA_real_, idx_b),
      FacetB_Measure = meas_b,
      FacetB_SE = se_b
    )
    seq_id <- seq_id + 1
  }

  bias_tbl <- bind_rows(rows)
  if (nrow(bias_tbl) == 0) return(list())

  numeric_cols <- c("Observd Score", "Expctd Score", "Observd Count", "Obs-Exp Average", "Bias Size", "S.E.")
  pop_sd <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0) return(NA_real_)
    m <- mean(x)
    sqrt(mean((x - m)^2))
  }
  mean_row <- summarize(bias_tbl, across(all_of(numeric_cols), ~ mean(.x, na.rm = TRUE)))
  sd_pop_row <- summarize(bias_tbl, across(all_of(numeric_cols), ~ pop_sd(.x)))
  sd_sample_row <- summarize(bias_tbl, across(all_of(numeric_cols), ~ stats::sd(.x, na.rm = TRUE)))
  summary_tbl <- bind_rows(
    mean_row,
    sd_pop_row,
    sd_sample_row
  ) |>
    mutate(Statistic = c(paste0("Mean (Count: ", nrow(bias_tbl), ")"), "S.D. (Population)", "S.D. (Sample)"), .before = 1)

  se_bias <- bias_tbl$`S.E.`
  bias_vals <- bias_tbl$`Bias Size`
  w_chi <- ifelse(is.finite(se_bias) & se_bias > 0, 1 / (se_bias^2), NA_real_)
  ok <- is.finite(w_chi) & is.finite(bias_vals)
  fixed_chi <- if (sum(ok) >= 2) {
    sum(w_chi[ok] * bias_vals[ok]^2) - (sum(w_chi[ok] * bias_vals[ok])^2) / sum(w_chi[ok])
  } else {
    NA_real_
  }
  fixed_df <- max(nrow(bias_tbl) - 1, 0)
  fixed_prob <- ifelse(is.finite(fixed_chi) && fixed_df > 0, 1 - stats::pchisq(fixed_chi, df = fixed_df), NA_real_)
  chi_tbl <- tibble(
    FixedChiSq = fixed_chi,
    FixedDF = fixed_df,
    FixedProb = fixed_prob
  )

  list(
    facet_a = facet_a,
    facet_b = facet_b,
    table = bias_tbl,
    summary = summary_tbl,
    chi_sq = chi_tbl,
    iteration = bind_rows(iter_rows)
  )
}

calc_bias_pairwise <- function(bias_tbl, target_facet, context_facet) {
  if (is.null(bias_tbl) || nrow(bias_tbl) == 0) return(tibble())

  use_a <- bias_tbl$FacetA[1] == target_facet
  use_b <- bias_tbl$FacetB[1] == target_facet
  if (!(use_a || use_b)) return(tibble())

  target_prefix <- if (use_a) "FacetA" else "FacetB"
  context_prefix <- if (use_a) "FacetB" else "FacetA"

  sub <- bias_tbl |>
    filter(.data[[context_prefix]] %in% context_facet)
  if (nrow(sub) == 0) sub <- bias_tbl

  rows <- list()
  for (tgt_level in unique(sub[[paste0(target_prefix, "_Level")]])) {
    group_df <- sub |> filter(.data[[paste0(target_prefix, "_Level")]] == tgt_level)
    contexts <- unique(group_df[[paste0(context_prefix, "_Level")]])
    if (length(contexts) < 2) next
    ctx_pairs <- combn(contexts, 2, simplify = FALSE)
    for (pair in ctx_pairs) {
      c1 <- pair[1]
      c2 <- pair[2]
      r1 <- group_df |> filter(.data[[paste0(context_prefix, "_Level")]] == c1) |> slice(1)
      r2 <- group_df |> filter(.data[[paste0(context_prefix, "_Level")]] == c2) |> slice(1)

      tgt_measure <- r1[[paste0(target_prefix, "_Measure")]]
      tgt_se <- r1[[paste0(target_prefix, "_SE")]]
      tgt_index <- r1[[paste0(target_prefix, "_Index")]]

      bias1 <- r1$`Bias Size`
      bias2 <- r2$`Bias Size`
      bias_se1 <- r1$`S.E.`
      bias_se2 <- r2$`S.E.`

      local1 <- ifelse(is.finite(tgt_measure) & is.finite(bias1), tgt_measure + bias1, NA_real_)
      local2 <- ifelse(is.finite(tgt_measure) & is.finite(bias2), tgt_measure + bias2, NA_real_)
      se1 <- ifelse(is.finite(tgt_se) & is.finite(bias_se1), sqrt(tgt_se^2 + bias_se1^2), NA_real_)
      se2 <- ifelse(is.finite(tgt_se) & is.finite(bias_se2), sqrt(tgt_se^2 + bias_se2^2), NA_real_)

      contrast <- ifelse(is.finite(local1) & is.finite(local2), local1 - local2, NA_real_)
      se_contrast <- ifelse(is.finite(se1) & is.finite(se2), sqrt(se1^2 + se2^2), NA_real_)
      t_val <- ifelse(is.finite(contrast) & is.finite(se_contrast) & se_contrast > 0, contrast / se_contrast, NA_real_)

      n1 <- ifelse(is.finite(r1$ObsN), r1$ObsN, 0)
      n2 <- ifelse(is.finite(r2$ObsN), r2$ObsN, 0)
      df_num <- (se1^2 + se2^2)^2
      df_den <- 0
      if (is.finite(se1) && n1 > 1) df_den <- df_den + (se1^4) / (n1 - 1)
      if (is.finite(se2) && n2 > 1) df_den <- df_den + (se2^4) / (n2 - 1)
      df_t <- ifelse(df_den > 0, df_num / df_den, NA_real_)
      p_val <- ifelse(is.finite(t_val) & is.finite(df_t) & df_t > 0,
                      2 * stats::pt(-abs(t_val), df = df_t),
                      NA_real_)

      rows[[length(rows) + 1]] <- tibble(
        Target = tgt_level,
        `Target N` = tgt_index,
        `Target Measure` = tgt_measure,
        `Target S.E.` = tgt_se,
        Context1 = c1,
        `Context1 N` = r1[[paste0(context_prefix, "_Index")]],
        `Local Measure1` = local1,
        SE1 = se1,
        `Obs-Exp Avg1` = r1$`Obs-Exp Average`,
        Count1 = r1$`Observd Count`,
        Context2 = c2,
        `Context2 N` = r2[[paste0(context_prefix, "_Index")]],
        `Local Measure2` = local2,
        SE2 = se2,
        `Obs-Exp Avg2` = r2$`Obs-Exp Average`,
        Count2 = r2$`Observd Count`,
        Contrast = contrast,
        SE = se_contrast,
        t = t_val,
        `d.f.` = df_t,
        `Prob.` = p_val
      )
    }
  }
  bind_rows(rows)
}

extract_anchor_tables <- function(config) {
  facet_names <- c("Person", config$facet_names)
  anchor_tbl <- purrr::map_dfr(facet_names, function(facet) {
    spec <- if (facet == "Person") config$theta_spec else config$facet_specs[[facet]]
    anchors <- spec$anchors
    if (is.null(anchors)) return(tibble())
    df <- tibble(Facet = facet, Level = names(anchors), Anchor = as.numeric(anchors)) |>
      filter(is.finite(Anchor))
    if (nrow(df) == 0) return(tibble())
    df |>
      mutate(Source = ifelse(facet %in% config$dummy_facets, "Dummy facet", "Anchor"))
  })

  group_tbl <- purrr::map_dfr(facet_names, function(facet) {
    spec <- if (facet == "Person") config$theta_spec else config$facet_specs[[facet]]
    groups <- spec$groups
    if (is.null(groups)) return(tibble())
    df <- tibble(Facet = facet, Level = names(groups), Group = as.character(groups)) |>
      filter(!is.na(Group), Group != "")
    if (nrow(df) == 0) return(tibble())
    group_values <- spec$group_values
    df |>
      mutate(GroupValue = ifelse(Group %in% names(group_values), group_values[Group], NA_real_))
  })

  list(anchors = anchor_tbl, groups = group_tbl)
}

calc_reliability <- function(measure_df) {
  measure_df |>
    group_by(Facet) |>
    summarize(
      Levels = n(),
      SD = sd(Estimate, na.rm = TRUE),
      RMSE = sqrt(mean(SE^2, na.rm = TRUE)),
      Separation = ifelse(RMSE > 0, SD / RMSE, NA_real_),
      Strata = ifelse(RMSE > 0, (4 * (SD / RMSE) + 1) / 3, NA_real_),
      Reliability = ifelse(RMSE > 0, (Separation^2) / (1 + Separation^2), NA_real_),
      MeanInfit = mean(Infit, na.rm = TRUE),
      MeanOutfit = mean(Outfit, na.rm = TRUE),
      .groups = "drop"
    )
}

calc_facets_chisq <- function(measure_df) {
  if (is.null(measure_df) || nrow(measure_df) == 0) return(tibble())
  measure_df |>
    group_by(Facet) |>
    summarize(
      Levels = n(),
      MeanMeasure = mean(Estimate, na.rm = TRUE),
      SD = sd(Estimate, na.rm = TRUE),
      EstimateVec = list(Estimate),
      SEVec = list(SE),
      FixedChiSq = {
        w <- ifelse(is.finite(SE) & SE > 0, 1 / (SE^2), NA_real_)
        d <- Estimate
        ok <- is.finite(w) & is.finite(d)
        if (sum(ok) < 2) NA_real_ else sum(w[ok] * d[ok]^2) - (sum(w[ok] * d[ok])^2) / sum(w[ok])
      },
      FixedDF = pmax(Levels - 1, 0),
      RandomVar = {
        d <- Estimate
        se2 <- SE^2
        ok <- is.finite(d) & is.finite(se2)
        if (sum(ok) < 2) NA_real_ else (sum((d[ok] - mean(d[ok]))^2) / (sum(ok) - 1)) - (sum(se2[ok]) / sum(ok))
      },
      .groups = "drop"
    ) |>
    mutate(
      FixedProb = ifelse(is.finite(FixedChiSq) & FixedDF > 0,
                         1 - stats::pchisq(FixedChiSq, df = FixedDF),
                         NA_real_),
      RandomVar = ifelse(is.finite(RandomVar) & RandomVar > 0, RandomVar, NA_real_)
    ) |>
    rowwise() |>
    mutate(
      RandomChiSq = {
        if (!is.finite(RandomVar) || RandomVar <= 0) {
          NA_real_
        } else {
          d <- unlist(EstimateVec)
          se <- unlist(SEVec)
          w <- ifelse(is.finite(se) & se > 0, 1 / (RandomVar + se^2), NA_real_)
          ok <- is.finite(w) & is.finite(d)
          if (sum(ok) < 2) NA_real_ else sum(w[ok] * d[ok]^2) - (sum(w[ok] * d[ok])^2) / sum(w[ok])
        }
      },
      RandomDF = pmax(Levels - 2, 0),
      RandomProb = ifelse(is.finite(RandomChiSq) & RandomDF > 0,
                          1 - stats::pchisq(RandomChiSq, df = RandomDF),
                          NA_real_)
    ) |>
    ungroup() |>
    select(-EstimateVec, -SEVec)
}

ensure_positive_definite <- function(mat) {
  eig <- tryCatch(eigen(mat, symmetric = TRUE, only.values = TRUE)$values, error = function(e) NULL)
  if (is.null(eig)) return(mat)
  if (any(eig < .Machine$double.eps)) {
    smoothed <- tryCatch(suppressWarnings(psych::cor.smooth(mat)), error = function(e) NULL)
    if (!is.null(smoothed)) {
      if (is.list(smoothed) && !is.null(smoothed$R)) {
        mat <- smoothed$R
      } else {
        mat <- smoothed
      }
    }
  }
  mat
}

compute_pca_overall <- function(obs_df, facet_names) {
  if (length(facet_names) == 0) return(NULL)
  df_aug <- obs_df |>
    mutate(
      Person = as.character(Person),
      item_combination = paste(!!!rlang::syms(facet_names), sep = "_")
    ) |>
    select(Person, item_combination, StdResidual)

  residual_matrix_prep <- df_aug |>
    group_by(Person, item_combination) |>
    summarize(std_residual = mean(StdResidual, na.rm = TRUE), .groups = "drop")

  residual_matrix_wide <- tryCatch({
    residual_matrix_prep |>
      tidyr::pivot_wider(
        id_cols = Person,
        names_from = item_combination,
        values_from = std_residual,
        values_fill = list(std_residual = NA)
      ) |>
      tibble::column_to_rownames("Person")
  }, error = function(e) NULL)

  if (is.null(residual_matrix_wide) || nrow(residual_matrix_wide) < 2 || ncol(residual_matrix_wide) < 2) {
    return(NULL)
  }

  residual_matrix_clean <- residual_matrix_wide[, colSums(is.na(residual_matrix_wide)) < nrow(residual_matrix_wide), drop = FALSE]
  if (ncol(residual_matrix_clean) < 2) return(NULL)

  cor_matrix <- tryCatch(suppressWarnings(stats::cor(residual_matrix_clean, use = "pairwise.complete.obs")), error = function(e) NULL)
  if (is.null(cor_matrix)) return(NULL)
  cor_matrix[is.na(cor_matrix)] <- 0
  diag(cor_matrix) <- 1
  cor_matrix <- ensure_positive_definite(cor_matrix)

  n_factors <- max(1, min(10, ncol(cor_matrix) - 1, nrow(cor_matrix) - 1))
  pca_result <- tryCatch(psych::principal(cor_matrix, nfactors = n_factors, rotate = "none"), error = function(e) NULL)
  list(pca = pca_result, residual_matrix = residual_matrix_wide, cor_matrix = cor_matrix)
}

compute_pca_by_facet <- function(obs_df, facet_names) {
  out <- list()
  for (facet in facet_names) {
    facet_sym <- rlang::sym(facet)
    prep <- obs_df |>
      mutate(Person = as.character(Person), .Level = as.character(!!facet_sym)) |>
      select(Person, .Level, StdResidual) |>
      group_by(Person, .Level) |>
      summarize(std_residual = mean(StdResidual, na.rm = TRUE), .groups = "drop")

    wide <- tryCatch({
      tidyr::pivot_wider(
        prep,
        id_cols = Person,
        names_from = .Level,
        values_from = std_residual,
        values_fill = list(std_residual = NA)
      ) |>
        tibble::column_to_rownames("Person")
    }, error = function(e) NULL)

    if (is.null(wide) || nrow(wide) < 2 || ncol(wide) < 2) {
      out[[facet]] <- NULL
      next
    }

    keep <- colSums(is.na(wide)) < nrow(wide)
    wide <- wide[, keep, drop = FALSE]
    if (ncol(wide) < 2) {
      out[[facet]] <- NULL
      next
    }

    cor_mat <- tryCatch(suppressWarnings(stats::cor(wide, use = "pairwise.complete.obs")), error = function(e) NULL)
    if (is.null(cor_mat)) {
      out[[facet]] <- NULL
      next
    }
    cor_mat[is.na(cor_mat)] <- 0
    diag(cor_mat) <- 1
    cor_mat <- ensure_positive_definite(cor_mat)

    n_factors <- max(1, min(10, ncol(cor_mat) - 1, nrow(cor_mat) - 1))
    pca_obj <- tryCatch(psych::principal(cor_mat, nfactors = n_factors, rotate = "none"), error = function(e) NULL)
    out[[facet]] <- list(pca = pca_obj, cor_matrix = cor_mat, residual_matrix = wide)
  }
  out
}

mfrm_diagnostics <- function(res, interaction_pairs = NULL, top_n_interactions = 20, whexact = FALSE) {
  obs_df <- compute_obs_table(res)
  facet_cols <- c("Person", res$config$facet_names)
  overall_fit <- calc_overall_fit(obs_df, whexact = whexact)
  fit_tbl <- calc_facet_fit(obs_df, facet_cols, whexact = whexact)
  se_tbl <- calc_facet_se(obs_df, facet_cols)
  bias_tbl <- calc_bias_facet(obs_df, facet_cols)
  interaction_tbl <- calc_bias_interactions(obs_df, facet_cols, pairs = interaction_pairs,
                                             top_n = top_n_interactions)
  ptmea_tbl <- calc_ptmea(obs_df, facet_cols)
  subset_tbls <- calc_subsets(obs_df, facet_cols)

  person_tbl <- res$facets$person |>
    mutate(
      Facet = "Person",
      Level = Person,
      SE = if ("SD" %in% names(res$facets$person)) SD else NA_real_
    )
  facet_tbl <- res$facets$others |>
    mutate(Level = as.character(Level), SE = NA_real_)

  measures <- bind_rows(
    person_tbl |>
      select(Facet, Level, Estimate, SE),
    facet_tbl |>
      select(Facet, Level, Estimate, SE)
  ) |>
    left_join(se_tbl, by = c("Facet", "Level")) |>
    mutate(SE = dplyr::coalesce(SE.x, SE.y)) |>
    select(-SE.x, -SE.y) |>
    left_join(fit_tbl, by = c("Facet", "Level")) |>
    left_join(bias_tbl, by = c("Facet", "Level")) |>
    left_join(ptmea_tbl, by = c("Facet", "Level")) |>
    mutate(
      CI_Lower = ifelse(is.finite(SE), Estimate - 1.96 * SE, NA_real_),
      CI_Upper = ifelse(is.finite(SE), Estimate + 1.96 * SE, NA_real_)
    )

  if (!"N" %in% names(measures)) {
    if ("N.x" %in% names(measures) || "N.y" %in% names(measures)) {
      measures <- measures |>
        mutate(N = dplyr::coalesce(.data$N.x, .data$N.y))
    } else {
      measures <- measures |>
        mutate(N = NA_real_)
    }
  }

  reliability_tbl <- calc_reliability(measures)

  list(
    obs = obs_df,
    overall_fit = overall_fit,
    measures = measures,
    fit = fit_tbl,
    reliability = reliability_tbl,
    bias = bias_tbl,
    interactions = interaction_tbl,
    subsets = subset_tbls
  )
}

# ---- Shiny app ----
ui <- fluidPage(
  titlePanel("MFRM Estimation (No FACETS / No ordinal)"),
  sidebarLayout(
    sidebarPanel(
      radioButtons(
        "data_source",
        "Data source",
        choices = c("Sample data (built-in)" = "sample", "Paste table" = "paste", "Upload file" = "upload"),
        selected = "sample"
      ),
      conditionalPanel(
        "input.data_source == 'upload'",
        fileInput("file", "Upload CSV/TSV", accept = c(".csv", ".tsv", ".txt")),
        checkboxInput("header", "Header", TRUE),
        radioButtons("sep", "Separator",
                     choices = c(Comma = ",", Tab = "\t", Semicolon = ";"),
                     selected = ",")
      ),
      conditionalPanel(
        "input.data_source == 'paste'",
        helpText("Paste data into the template. Tab-delimited is recommended."),
        textAreaInput(
          "pasted_data",
          "Paste table data",
          value = template_tab_text,
          rows = 10
        ),
        checkboxInput("paste_header", "Header row", TRUE),
        radioButtons("paste_sep", "Separator",
                     choices = c(Tab = "\t", Comma = ",", Semicolon = ";"),
                     selected = "\t", inline = TRUE),
        tags$div(
          actionButton("paste_clear", "Clear table"),
          actionButton("paste_demo", "Load demo template"),
          actionButton("paste_toy", "Load toy template")
        )
      ),
      conditionalPanel(
        "input.data_source == 'sample'",
        helpText("Sample data is synthetic and included for quick demonstration.")
      ),
      conditionalPanel(
        "input.data_source != 'upload'",
        tags$div(
          downloadButton("download_sample_tsv", "Download sample TSV"),
          downloadButton("download_sample_csv", "Download sample CSV")
        )
      ),
      uiOutput("col_selectors"),
      uiOutput("facet_actions"),
      uiOutput("constraint_ui"),
      helpText("Select a Person column and 2+ facet columns (total facets >= 3)."),
      uiOutput("interaction_pair_selector"),
      numericInput("top_n_interactions", "Top bias interactions", value = 50, min = 5, step = 5),
      tags$h4("Bias/Interaction (Facets-style)"),
      checkboxInput("bias_run", "Estimate bias/interaction", value = FALSE),
      conditionalPanel(
        "input.bias_run == true",
        uiOutput("bias_pair_selector"),
        numericInput("bias_max_abs", "Max |bias| (logit)", value = 10, step = 1),
        checkboxInput("bias_omit_extreme", "Omit extreme elements", value = TRUE)
      ),
      tags$h4("Agreement"),
      uiOutput("agreement_rater_selector"),
      numericInput("top_n_misfit", "Top misfit levels (|ZSTD|)", value = 20, min = 5, step = 5),
      sliderInput("misfit_threshold", "Misfit ZSTD threshold", min = 0, max = 5, value = 2, step = 0.5),
      radioButtons(
        "misfit_compare",
        "Misfit threshold mode",
        choices = c(">=" = "gte", "<=" = "lte"),
        selected = "gte",
        inline = TRUE
      ),
      uiOutput("fit_facet_filter_ui"),
      uiOutput("fit_facet_quick_ui"),
      tags$hr(),
      tags$h4("Reporting (FACETS)"),
      radioButtons(
        "report_totalscore",
        "Total scores reported",
        choices = c("Yes (Total Score/Count)" = "yes", "No (Obsvd Score/Count)" = "no"),
        selected = "yes"
      ),
      radioButtons(
        "report_omit_unobserved",
        "Omit unobserved elements",
        choices = c("No" = "no", "Yes" = "yes"),
        selected = "no"
      ),
      numericInput("report_xtreme", "Xtreme correction (fraction)", value = 0, min = 0, step = 0.1),
      numericInput("report_umean", "Umean (origin)", value = 0, step = 0.1),
      numericInput("report_uscale", "Uscale (units per logit)", value = 1, step = 0.1),
      numericInput("report_udecimals", "Udecimals (report precision)", value = 2, min = 0, step = 1),
      numericInput("rating_min", "Rating min", value = 1, step = 1),
      numericInput("rating_max", "Rating max", value = 5, step = 1),
      radioButtons("model", "Model", choices = c("RSM", "PCM"), selected = "RSM"),
      radioButtons("method", "Estimation", choices = c("JMLE", "MML"), selected = "JMLE"),
      checkboxInput("keep_original", "Keep original category values (K)", value = FALSE),
      conditionalPanel(
        "input.model == 'PCM'",
        uiOutput("step_facet_selector"),
        uiOutput("pcm_curve_selector")
      ),
      conditionalPanel(
        "input.method == 'MML'",
        numericInput("quad_points", "Quadrature points (MML)", value = 15, min = 5, step = 2)
      ),
      numericInput("maxit", "Max iterations", value = 400, min = 50, step = 50),
      checkboxInput("whexact", "WHEXACT (Exact ZSTD)", value = FALSE),
      actionButton("run", "Run estimation", class = "btn-primary")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel(
          "Data",
          gt_output("data_gt"),
          gt_output("facet_summary_gt"),
          plotlyOutput("response_dist_plot", height = "300px"),
          gt_output("data_quality_gt_preview")
        ),
        tabPanel("Summary", gt_output("summary_gt")),
        tabPanel(
          "Results Summary",
          conditionalPanel(
            "input.run > 0",
            h4("Quick Overview"),
            gt_output("results_summary_gt"),
            h4("Reliability Snapshot"),
            gt_output("results_reliability_snapshot_gt"),
            h4("Key Visuals"),
            plotlyOutput("results_reliability_plot", height = "320px"),
            plotlyOutput("results_fit_overview_plot", height = "320px"),
            h4("Flags (Low Reliability / Fit)"),
            gt_output("results_flags_gt"),
            h4("Quick Downloads"),
            downloadButton("download_summary_quick", "Summary (CSV)", class = "btn btn-primary"),
            downloadButton("download_reliability_quick", "Reliability (CSV)", class = "btn btn-default"),
            downloadButton("download_fit_quick", "Fit Statistics (CSV)", class = "btn btn-default"),
            downloadButton("download_steps_quick", "Thresholds (CSV)", class = "btn btn-default"),
            downloadButton("download_zip", "Download Key Outputs (ZIP)", class = "btn btn-info"),
            h4("Caution / Diagnostic Notes"),
            verbatimTextOutput("cautionText"),
            h4("Data Quality Checks"),
            gt_output("data_quality_gt"),
            h4("Interpretation Hints"),
            uiOutput("results_interpretation"),
            h4("Reproducibility Info"),
            uiOutput("results_repro")
          )
        ),
        tabPanel("Person", gt_output("person_gt"), plotlyOutput("person_hist")),
        tabPanel(
          "Facets",
          uiOutput("facets_report_tabs"),
          gt_output("facets_report_gt")
        ),
        tabPanel("Fit", gt_output("fit_overall_gt"), gt_output("fit_facet_gt")),
        tabPanel(
          "Reliability",
          gt_output("reliability_gt"),
          gt_output("reliability_chisq_gt"),
          h4("Inter-rater agreement"),
          gt_output("agreement_summary_gt"),
          gt_output("agreement_pairs_gt")
        ),
        tabPanel(
          "Bias",
          conditionalPanel(
            "input.bias_run == true",
            h4("Bias/Interaction (Facets-style)"),
            gt_output("bias_facets_gt"),
            h4("Bias iteration report"),
            gt_output("bias_iteration_gt"),
            h4("Bias summary"),
            gt_output("bias_summary_gt"),
            h4("Fixed (all = 0) chi-square"),
            gt_output("bias_chi_gt"),
            h4("Pairwise bias report"),
            uiOutput("bias_pairwise_ui"),
            gt_output("bias_pairwise_gt")
          ),
          h4("Legacy residual bias (main model)"),
          gt_output("bias_gt")
        ),
        tabPanel("Interactions", gt_output("bias_interaction_gt")),
        tabPanel("Steps", gt_output("steps_gt")),
        tabPanel(
          "Categories",
          plotlyOutput("category_usage_plot", height = "320px"),
          plotlyOutput("category_percent_plot", height = "320px"),
          plotlyOutput("category_resid_plot", height = "320px"),
          plotlyOutput("category_fit_plot", height = "320px"),
          gt_output("category_stats_gt"),
          gt_output("category_thresholds_gt"),
          h4("Category warnings"),
          verbatimTextOutput("category_warnings_text")
        ),
        tabPanel(
          "Constraints",
          gt_output("constraint_summary_gt"),
          gt_output("constraint_anchors_gt"),
          gt_output("constraint_groups_gt"),
          h4("Subsets / connectivity"),
          gt_output("subset_gt")
        ),
        tabPanel("Settings", gt_output("settings_gt")),
        tabPanel(
          "Downloads",
          p("Export current results as CSV tables."),
          downloadButton("download_zip", "Download all tables (ZIP)"),
          downloadButton("download_measures", "Download measures (CSV)"),
          downloadButton("download_scorefile", "Download scorefile (CSV)"),
          downloadButton("download_residuals", "Download residuals (CSV)")
        ),
        tabPanel(
          "Visualizations",
          tabsetPanel(
            tabPanel("Wright Map", plotlyOutput("wright_plot", height = "320px")),
            tabPanel("Pathway Map", plotlyOutput("pathway_plot", height = "320px")),
            tabPanel("Facet Distribution", plotlyOutput("facet_plot", height = "320px")),
            tabPanel("Steps/Thresholds", plotlyOutput("step_plot", height = "320px")),
            tabPanel("Category Curves", plotlyOutput("category_plot", height = "360px")),
            tabPanel("Observed vs Expected", plotlyOutput("obs_exp_plot", height = "320px")),
            tabPanel("Fit Scatter", plotlyOutput("fit_scatter_plot", height = "320px")),
            tabPanel("Fit ZSTD", plotlyOutput("fit_zstd_plot", height = "320px")),
            tabPanel("Misfit Levels", plotlyOutput("misfit_bar_plot", height = "360px"))
          )
        ),
        tabPanel(
          "Dimensionality",
          uiOutput("pca_facet_selector"),
          plotlyOutput("pca_scree_plot", height = "320px"),
          gt_output("pca_eigen_gt"),
          plotlyOutput("pca_loadings_plot", height = "360px"),
          plotlyOutput("pca_biplot", height = "360px"),
          gt_output("dimensionality_gt")
        ),
        tabPanel(
          "Help",
          tags$h4("Overview"),
          tags$ul(
            tags$li("Select columns: choose a Person column + at least two facet columns + a Score column (Person is treated as \"Person\" internally)."),
            tags$li("Data source: Sample / Paste / Upload. For Paste, fill the template."),
            tags$li("Facet selection: use Select all / Clear / Suggest, or add by regex pattern."),
            tags$li("Facet summary: the Data tab shows levels, missing counts, and example values."),
            tags$li("Model: choose RSM or PCM. PCM requires a Step Facet."),
            tags$li("Fit: Infit/Outfit and ZSTD (Wilson-Hilferty approximation)."),
            tags$li("Bias: mean residuals and standardized residuals; Interactions include Person x Facet pairs."),
            tags$li("Agreement: pairwise exact/adjacent agreement for a selected facet (e.g., Rater) within the same context."),
            tags$li("PTMEA: correlation between Observed and PersonMeasure by non-Person facet level."),
            tags$li("Pathway Map: Measure vs Infit ZSTD with labels for |ZSTD| >= 2."),
            tags$li("Fit visuals: Fit Scatter / Fit ZSTD / Misfit Levels."),
            tags$li("Fit filters: use Fit facets to display, plus quick buttons."),
            tags$li("Misfit threshold: adjust |ZSTD| cutoff and >= / <= mode."),
            tags$li("Top bias interactions / top misfit levels limit large tables."),
            tags$li("Dimensionality: residual PCA diagnostics."),
            tags$li("Subsets: connectivity check for disconnected groups of facet levels."),
            tags$li("Downloads tab provides CSV/ZIP outputs.")
          ),
          tags$h4("Formulas (based on prior literature)"),
          tags$p("The app follows standard formulations from FACETS/Winsteps manuals and key MFRM literature; the implementation is a simplified version for teaching and exploration."),
          tags$h5("MFRM (RSM / PCM)"),
          tags$pre(
            "RSM (common steps):\n",
            "  log( P_nik / P_ni(k-1) ) = theta_n - sum_f delta_f(level) - tau_k\n",
            "\n",
            "PCM (Step Facet = s):\n",
            "  log( P_nik / P_ni(k-1) ) = theta_n - sum_f delta_f(level) - tau_{s(level),k}\n",
            "\n",
            "eta_n = theta_n - sum_f delta_f(level)\n",
            "P_nik is the probability of category k (k = 0..K).\n",
            "\n",
            "Identification constraints: sum delta = 0; sum tau = 0."
          ),
          tags$h5("Expected score and residuals"),
          tags$pre(
            "E[Score] = r_min + sum_k k * P_k\n",
            "Var = sum_k k^2 * P_k - (sum_k k * P_k)^2\n",
            "Residual = Observed - Expected\n",
            "StdResidual = Residual / sqrt(Var)"
          ),
          tags$h5("Fit statistics (Infit / Outfit / ZSTD)"),
          tags$pre(
            "Infit MNSQ = sum(StdResidual^2 * Var) / sum Var\n",
            "Outfit MNSQ = mean(StdResidual^2)\n",
            "ZSTD approx = (MNSQ^(1/3) - (1 - 2/(9*df))) / sqrt(2/(9*df))\n",
            "df: Infit uses sum Var; Outfit uses N (Wilson-Hilferty approximation)."
          ),
          tags$h5("Reliability and separation"),
          tags$pre(
            "RMSE = sqrt( mean(SE^2) )\n",
            "Separation = SD / RMSE\n",
            "Strata = (4 * Separation + 1) / 3\n",
            "Reliability = Separation^2 / (1 + Separation^2)"
          ),
          tags$h5("Confidence intervals"),
          tags$pre(
            "95% CI = Estimate +/- 1.96 * SE"
          ),
          tags$h5("PTMEA / bias / interactions"),
          tags$pre(
            "PTMEA (non-Person) = cor(Observed, PersonMeasure)\n",
            "Bias = ObservedAverage - ExpectedAverage\n",
            "MeanResidual = mean(Residual)\n",
            "MeanStdResidual = mean(StdResidual)\n",
            "SE(MeanResidual) = sqrt(sum(Var)) / N\n",
            "SE(MeanStdResidual) = 1 / sqrt(N)\n",
            "t = Mean / SE (df approx = N - 1)\n",
            "p = 2 * pt(-|t|, df)\n",
            "Interactions: compute mean residuals and standardized residuals by facet x facet"
          ),
          tags$h5("FACETS-style reporting (Table 7)"),
          tags$ul(
            tags$li("Totalscore = Yes: use all responses (includes extreme elements). Totalscore = No: use active responses only (extreme elements excluded)."),
            tags$li("Fair(M) Average: expected score computed at facet-mean baselines; Fair(Z) Average: expected score computed at facet zero baselines."),
            tags$li("Umean/Uscale: reported Measure = logit * Uscale + Umean; SEs are scaled by |Uscale|."),
            tags$li("Omit unobserved: removes elements with zero observations from the Table 7 report."),
            tags$li("Xtreme correction: when > 0, extreme measures are approximated by solving for a target expected score (rating_min + Xtreme or rating_max - Xtreme).")
          ),
          tags$h5("Dimensionality (Residual PCA)"),
          tags$pre(
            "Person x facet-combination residual matrix -> correlation matrix -> PCA\n",
            "Large PC1 variance or small PC1/PC2 ratio may indicate secondary dimensions"
          ),
          tags$h4("Assumptions and limitations"),
          tags$ul(
            tags$li("Unidimensionality: the model assumes a dominant single latent dimension. Use the Dimensionality tab to check residual PCA."),
            tags$li("Local independence: residuals should not show strong systematic patterns after conditioning on measures."),
            tags$li("Ordered categories: thresholds should progress monotonically; severe disordering suggests category issues."),
            tags$li("Sparse cells: small N per level or empty category usage can destabilize estimates and fit statistics."),
            tags$li("Extreme scores: persons or facet levels with all min/max scores can yield unstable measures."),
            tags$li("Xtreme correction is an approximation based on inverse expected scores; results can differ from FACETS."),
            tags$li("Model fit depends on design quality (balanced designs generally provide more stable estimates).")
          ),
          tags$h4("Data requirements checklist"),
          tags$ul(
            tags$li("Minimum structure: Person + at least two facet columns + a Score column."),
            tags$li("Levels: 3+ levels per facet are recommended; 1-2 levels behave like fixed effects."),
            tags$li("Coverage: ensure each facet level appears with multiple Persons and across tasks/raters."),
            tags$li("Category usage: all score categories should be observed; consider collapsing if categories are unused."),
            tags$li("Missing data: avoid systematic missingness; check the Data tab for missing counts."),
            tags$li("Scale order: verify that the score variable is ordinal and correctly coded.")
          ),
          tags$h4("Step/threshold diagnostics"),
          tags$ul(
            tags$li("Ordered thresholds: steps should increase with category; disordering signals category problems."),
            tags$li("Category curves: each category should peak at some region of the latent trait."),
            tags$li("Adjacent categories: heavy overlap with no distinct peaks suggests collapsing categories."),
            tags$li("Step spacing: extremely small spacing indicates redundant categories; large gaps may indicate unused midpoints."),
            tags$li("PCM: inspect step patterns by Step Facet level to detect rater- or task-specific category issues.")
          ),
          tags$h4("Anchoring and linking (Linacre)"),
          tags$ul(
            tags$li("Anchoring uses anchor values so different analyses are directly comparable."),
            tags$li("Anchor values are pre-set logit values assigned to objects, agents, or steps as reference points for calibration."),
            tags$li("Under JMLE (joint maximum likelihood), anchored values are held fixed during iterations; convergence is evaluated on unanchored estimates, and fit is computed as if anchored values had converged."),
            tags$li("Include the Umean (origin and scaling) from the anchoring analysis so the baseline for Fair Average and centering is consistent."),
            tags$li("One facet is typically designated as noncentered; other facets are centered unless anchored or group-anchored."),
            tags$li("Linking relates measures from one test with another so measures can be directly compared; equating puts measures from two tests in the same frame of reference.")
          ),
          tags$h4("Equating across datasets (Linacre)"),
          tags$ul(
            tags$li("Equate all facets except one (usually Persons, unless rater behavior is the focus)."),
            tags$li("Use at least five common elements per facet as a rule of thumb; fewer gives weak quality control."),
            tags$li("Run separate analyses, cross-plot common elements, and drop elements with clearly discrepant measures."),
            tags$li("If remaining common elements align near a line parallel to the identity line, anchor those elements or combine datasets."),
            tags$li("If common elements still form a cloud, use group-anchoring for that facet."),
            tags$li("Monitor criterion-level success rates and other indicators to confirm the equating makes sense.")
          ),
          tags$h4("Linkage and design requirements (Linacre)"),
          tags$ul(
            tags$li("The judging plan must provide enough linkage among all facets so parameters are estimable in one frame of reference."),
            tags$li("For robust estimation, a rule of thumb is about 30 observations per element and at least 10 observations per rating-scale category."),
            tags$li("Reduced designs can still be estimable if overlap is maintained, but precision decreases as observations drop."),
            tags$li("Disjoint or weakly connected subsets can fail to converge or converge to non-repeatable values.")
          ),
          tags$h4("Bias and interaction interpretation"),
          tags$ul(
            tags$li("Bias tables summarize mean residuals (raw and standardized) by facet level."),
            tags$li("Large absolute mean standardized residuals indicate systematic over- or under-scoring."),
            tags$li("Interaction tables highlight differential effects for specific facet pairs (e.g., Person x Rater)."),
            tags$li("Always check N: small cell counts can inflate apparent bias."),
            tags$li("p-values are approximate (two-tailed t with df = N - 1)."),
            tags$li("Use bias patterns to guide follow-up (rater training, rubric refinement, or task review).")
          ),
          tags$h4("Beginner Q&A (Linacre-based)"),
          tags$ul(
            tags$li("Q: What does anchoring do? A: Anchoring fixes selected measures or steps so different analyses are on a common scale."),
            tags$li("Q: Does anchoring change estimation? A: Under JMLE, anchored values are held fixed; convergence is evaluated on unanchored estimates."),
            tags$li("Q: Why include Umean from the anchoring run? A: It preserves the same origin and scaling so fair averages and centering are consistent."),
            tags$li("Q: What does linking/equating mean? A: Linking relates measures from different tests so they are comparable; equating puts them in the same frame of reference."),
            tags$li("Q: How many common elements are needed? A: Linacre gives a rule of thumb of at least five common elements per facet for quality control."),
            tags$li("Q: Why does convergence fail? A: Common reasons include over- or under-constraints, weak linkage, nested facets without anchoring, or very low category frequencies."),
            tags$li("Q: How should I interpret standardized residuals? A: When data fit the model, about 5% exceed +/-2 and about 1% exceed +/-3; concentrations indicate local misfit.")
          ),
          tags$h4("Worked examples (interpretation)"),
          tags$ul(
            tags$li("Example 1 (anchoring across administrations): Analyze Administration A, output anchor values, then analyze Administration B with those anchor values and the Umean from A so measures remain on the same scale."),
            tags$li("Example 2 (equating with common elements): Run separate analyses, cross-plot common elements, drop inconsistent elements, then anchor common elements or group-anchor to align datasets."),
            tags$li("Example 3 (minimal linkage): A reduced design can still estimate parameters if overlap connects all facets, but precision drops as observations are reduced."),
            tags$li("Example 4 (fit interpretation): If many elements show |ZSTD| > 3, investigate design, category use, or anchoring assumptions before interpreting measures."),
            tags$li("Example 5 (category issues): Disordered steps or very low category counts can destabilize estimation; consider collapsing categories or revising the rubric.")
          ),
          tags$h4("Story-based examples (beginner friendly)"),
          tags$ul(
            tags$li("Story 1: A writing program runs Session 1 and Session 2 with partly different raters. The analyst anchors a set of common tasks and includes the Umean from Session 1 so Session 2 measures sit on the same scale. The team can now compare person measures across sessions without shifting the origin."),
            tags$li("Story 2: Two schools share a few common tasks. The analyst runs separate analyses, cross-plots the common task measures, removes tasks that drifted, and anchors the stable tasks. The remaining network links the schools into one frame of reference."),
            tags$li("Story 3: A small pilot uses a sparse judging plan. Because each person is still linked to multiple tasks and raters, estimation is possible, but fit and SEs are less stable. The team treats results as preliminary and expands overlap in the next round."),
            tags$li("Story 4: A rater training session shows many |ZSTD| values above 3 for one rater. The team checks category usage and identifies a rubric misunderstanding, then retrains and re-estimates before reporting results.")
          ),
          tags$h4("Facet selection examples"),
          tags$ol(
            tags$li("Use Select all facets to review available columns."),
            tags$li("Remove unwanted columns or Clear facets to reset."),
            tags$li("Use Add by pattern (e.g., \"occasion|group\")."),
            tags$li("Check levels and missing values in the Data tab.")
          ),
          tags$h4("Key sources (APA style)"),
          tags$ul(
            tags$li("Linacre, J. M. (2024). A user's guide to FACETS (64-bit): Rasch-model computer programs (Program manual 4.2.3). Winsteps.com."),
            tags$li("Linacre, J. M. (2024). A user's guide to WINSTEPS/MINISTEP: Rasch-model computer programs (Program manual 5.8.1). Winsteps.com."),
            tags$li("Myford, C. M., & Wolfe, E. W. (2003). Detecting and measuring rater effects using many-facet Rasch measurement: Part I. Journal of Applied Measurement, 4(4), 386-422."),
            tags$li("Myford, C. M., & Wolfe, E. W. (2004). Detecting and measuring rater effects using many-facet Rasch measurement: Part II. Journal of Applied Measurement, 5(2), 189-227."),
            tags$li("Koizumi, R., Kaneko, E., Setoguchi, E., In'nami, Y., & Naganuma, N. (2019). Examination of CEFR-J spoken interaction tasks using many-facet Rasch measurement and generalizability theory. Papers in Language Testing and Assessment, 8(2), 1-33."),
            tags$li("Uto, M., & Ueno, M. (2020). A generalized many-facet Rasch model and its Bayesian estimation using Hamiltonian Monte Carlo. Behaviormetrika, 47, 469-496."),
            tags$li("Eckes, T., & Jin, K.-Y. (2021). Measuring rater centrality effects in writing assessment: A Bayesian facets modeling approach. Psychological Test and Assessment Modeling, 63(1), 65-94.")
          ),
          tags$h4("Quick workflow"),
          tags$pre(
            "1) Data source -> 2) Person/Facets/Score -> 3) Model/Method\n",
            "            |\n",
            "            v\n",
            "4) Run estimation -> 5) Fit/Bias/Plots -> 6) Download\n"
          )
        ),
        tabPanel(
          "References",
          tags$ul(
            tags$li("Bond, T. G., & Fox, C. M. (2007). Applying the Rasch model: Fundamental measurement in the human sciences (2nd ed.). Lawrence Erlbaum Associates."),
            tags$li("Wright, B. D., & Masters, G. N. (1982). Rating scale analysis. MESA Press."),
            tags$li("Linacre, J. M. (2024). A user's guide to FACETS (64-bit): Rasch-model computer programs (Program manual 4.2.3). Winsteps.com."),
            tags$li("Linacre, J. M. (2024). A user's guide to WINSTEPS/MINISTEP: Rasch-model computer programs (Program manual 5.8.1). Winsteps.com."),
            tags$li("Myford, C. M., & Wolfe, E. W. (2003). Detecting and measuring rater effects using many-facet Rasch measurement: Part I. Journal of Applied Measurement, 4(4), 386-422."),
            tags$li("Myford, C. M., & Wolfe, E. W. (2004). Detecting and measuring rater effects using many-facet Rasch measurement: Part II. Journal of Applied Measurement, 5(2), 189-227."),
            tags$li("Koizumi, R., Kaneko, E., Setoguchi, E., In'nami, Y., & Naganuma, N. (2019). Examination of CEFR-J spoken interaction tasks using many-facet Rasch measurement and generalizability theory. Papers in Language Testing and Assessment, 8(2), 1-33."),
            tags$li("Uto, M., & Ueno, M. (2020). A generalized many-facet Rasch model and its Bayesian estimation using Hamiltonian Monte Carlo. Behaviormetrika, 47, 469-496."),
            tags$li("Eckes, T., & Jin, K.-Y. (2021). Measuring rater centrality effects in writing assessment: A Bayesian facets modeling approach. Psychological Test and Assessment Modeling, 63(1), 65-94.")
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  raw_data <- reactive({
    if (identical(input$data_source, "sample")) {
      return(sample_mfrm_data())
    }
    if (identical(input$data_source, "paste")) {
      pasted <- input$pasted_data
      if (is.null(pasted)) pasted <- ""
      if (!nzchar(trimws(pasted))) {
        return(template_tab_source_demo[0, ])
      }
      df <- tryCatch({
        read.table(
          text = pasted,
          header = isTRUE(input$paste_header),
          sep = input$paste_sep,
          stringsAsFactors = FALSE,
          check.names = FALSE,
          comment.char = "",
          strip.white = TRUE,
          na.strings = c("", "NA")
        )
      }, error = function(e) {
        showNotification(paste("Paste parsing error:", e$message), type = "error", duration = 6)
        template_tab_source_demo[0, ]
      })
      if (ncol(df) == 0) {
        return(template_tab_source_demo[0, ])
      }
      blank_idx <- which(names(df) == "")
      if (length(blank_idx) > 0) {
        names(df)[blank_idx] <- paste0("Column_", blank_idx)
      }
      return(df)
    }
    req(input$file)
    read.csv(
      input$file$datapath,
      header = input$header,
      sep = input$sep,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })

  anchor_df <- reactive({
    df <- read_flexible_table(input$anchor_table, input$anchor_file)
    normalize_anchor_df(df)
  })

  group_anchor_df <- reactive({
    df <- read_flexible_table(input$group_anchor_table, input$group_anchor_file)
    normalize_group_anchor_df(df)
  })

  output$col_selectors <- renderUI({
    req(raw_data())
    cols <- names(raw_data())
    if (length(cols) == 0) {
      return(helpText("No columns detected. Paste or upload data first."))
    }
    person_sel <- guess_col(cols, c("person", "examinee", "student", "id"), 1)
    score_sel <- guess_col(cols, c("score", "rating", "result"), 5)
    facet_defaults <- unique(na.omit(c(
      intersect(cols, c("Rater", "rater"))[1],
      intersect(cols, c("Task", "task"))[1],
      intersect(cols, c("Criterion", "criteria", "criterion"))[1]
    )))
    if (length(facet_defaults) < 2) {
      facet_defaults <- setdiff(cols, c(person_sel, score_sel))[seq_len(min(2, length(cols)))]
    }
    tagList(
      selectInput("person_col", "Person column", choices = cols, selected = person_sel),
      selectizeInput(
        "facet_cols",
        "Facet columns (2+)",
        choices = cols,
        selected = facet_defaults,
        multiple = TRUE,
        options = list(plugins = list("remove_button"))
      ),
      selectInput("score_col", "Score column", choices = cols, selected = score_sel),
      selectInput("weight_col", "Weight column (optional)", choices = c("(None)", cols), selected = "(None)")
    )
  })

  observeEvent(input$paste_clear, {
    updateTextAreaInput(session, "pasted_data", value = template_header_text)
  })

  observeEvent(input$paste_demo, {
    updateTextAreaInput(session, "pasted_data", value = template_tab_text)
  })

  observeEvent(input$paste_toy, {
    updateTextAreaInput(session, "pasted_data", value = template_tab_text_toy)
  })

  output$facet_actions <- renderUI({
    req(raw_data())
    tagList(
      tags$div(
        actionButton("facet_select_all", "Select all facets"),
        actionButton("facet_clear", "Clear facets"),
        actionButton("facet_suggest", "Suggest facets")
      ),
      textInput(
        "facet_pattern",
        "Facet column pattern (regex)",
        value = "rater|task|criterion|occasion|group"
      ),
      actionButton("facet_apply_pattern", "Add by pattern")
    )
  })

  output$constraint_ui <- renderUI({
    facets <- c("Person", selected_facets())
    if (length(facets) == 0) {
      return(helpText("Select facets to enable constraints."))
    }
    tagList(
      tags$hr(),
      tags$h4("Constraints"),
      selectInput(
        "noncenter_facet",
        "Noncenter facet",
        choices = facets,
        selected = "Person"
      ),
      selectizeInput(
        "dummy_facets",
        "Dummy facets (anchor all levels at 0)",
        choices = facets,
        selected = character(0),
        multiple = TRUE,
        options = list(plugins = list("remove_button"))
      ),
      selectizeInput(
        "positive_facets",
        "Positive facets (higher values increase score)",
        choices = setdiff(facets, "Person"),
        selected = character(0),
        multiple = TRUE,
        options = list(plugins = list("remove_button"))
      ),
      helpText("Anchors: Facet, Level, Anchor (numeric)."),
      textAreaInput(
        "anchor_table",
        "Anchor table (CSV/TSV)",
        value = "Facet\tLevel\tAnchor\n",
        rows = 4
      ),
      fileInput("anchor_file", "Upload anchor CSV/TSV", accept = c(".csv", ".tsv", ".txt")),
      helpText("Group anchors: Facet, Level, Group, GroupValue."),
      textAreaInput(
        "group_anchor_table",
        "Group anchor table (CSV/TSV)",
        value = "Facet\tLevel\tGroup\tGroupValue\n",
        rows = 4
      ),
      fileInput("group_anchor_file", "Upload group anchor CSV/TSV", accept = c(".csv", ".tsv", ".txt"))
    )
  })

  output$fit_facet_filter_ui <- renderUI({
    facets <- c("Person", selected_facets())
    if (length(facets) == 0) {
      return(helpText("Select facets to enable fit filters."))
    }
    selectizeInput(
      "fit_facet_filter",
      "Fit facets to display",
      choices = c("All", facets),
      selected = "All",
      multiple = TRUE,
      options = list(plugins = list("remove_button"))
    )
  })

  output$fit_facet_quick_ui <- renderUI({
    facets <- selected_facets()
    if (length(facets) == 0) {
      return(helpText("Add facets to enable quick fit filters."))
    }
    tagList(
      tags$div(
        actionButton("fit_filter_all", "Fit: All"),
        actionButton("fit_filter_person", "Fit: Person only")
      ),
      selectInput(
        "fit_facet_single",
        "Fit: Single facet",
        choices = facets,
        selected = facets[1]
      ),
      actionButton("fit_filter_single", "Use single facet")
    )
  })

  output$step_facet_selector <- renderUI({
    req(input$facet_cols)
    facets <- selected_facets()
    facets <- facets[facets %in% names(raw_data())]
    if (length(facets) == 0) {
      return(helpText("Select facet columns to choose a PCM step facet."))
    }
    selectInput(
      "step_facet",
      "PCM step facet",
      choices = facets,
      selected = facets[1]
    )
  })

  output$pcm_curve_selector <- renderUI({
    if (!identical(input$model, "PCM")) return(NULL)
    req(raw_data(), input$step_facet)
    if (!input$step_facet %in% names(raw_data())) {
      return(helpText("Step facet not found in current data. Re-select facets."))
    }
    vals <- raw_data()[[input$step_facet]]
    choices <- unique(as.character(vals))
    if (length(choices) == 0) {
      return(helpText("No step facet levels detected."))
    }
    selectInput(
      "pcm_curve_criterion",
      "PCM curve level",
      choices = choices,
      selected = choices[1]
    )
  })

  observeEvent(raw_data(), {
    df <- raw_data()
    if (nrow(df) == 0) return(NULL)
    score_guess <- df[[names(df)[min(5, ncol(df))]]]
    if (is.numeric(score_guess)) {
      updateNumericInput(session, "rating_min", value = min(score_guess, na.rm = TRUE))
      updateNumericInput(session, "rating_max", value = max(score_guess, na.rm = TRUE))
    }
  })

  output$data_gt <- render_gt({
    req(raw_data())
    df <- raw_data()
    df |>
      head(20) |>
      gt() |>
      tab_header(
        title = "Data preview",
        subtitle = paste("Showing first 20 of", nrow(df), "rows")
      )
  })

  data_quality_tbl <- reactive({
    req(raw_data())
    df <- raw_data()
    tibble(
      Column = names(df),
      Type = vapply(df, function(x) paste(class(x), collapse = ","), character(1)),
      Missing = vapply(df, function(x) sum(is.na(x)), integer(1)),
      Distinct = vapply(df, function(x) length(unique(x)), integer(1)),
      Example = vapply(df, function(x) {
        vals <- unique(as.character(x))
        vals <- vals[!is.na(vals)]
        paste(head(vals, 4), collapse = ", ")
      }, character(1))
    )
  })

  output$data_quality_gt_preview <- render_gt({
    req(data_quality_tbl())
    gt(data_quality_tbl())
  })

  output$data_quality_gt <- render_gt({
    req(data_quality_tbl())
    gt(data_quality_tbl())
  })

  output$response_dist_plot <- renderPlotly({
    req(raw_data(), input$score_col)
    df <- raw_data()
    if (!input$score_col %in% names(df)) return(NULL)
    score_vals <- df[[input$score_col]]
    weight_col <- if (!is.null(input$weight_col) && input$weight_col != "(None)") input$weight_col else NULL
    if (!is.null(weight_col) && weight_col %in% names(df)) {
      plot_df <- df |>
        transmute(
          Score = as.character(.data[[input$score_col]]),
          Weight = suppressWarnings(as.numeric(.data[[weight_col]]))
        ) |>
        filter(!is.na(Score), is.finite(Weight)) |>
        group_by(Score) |>
        summarize(n = sum(Weight, na.rm = TRUE), .groups = "drop")
    } else {
      plot_df <- tibble(Score = as.character(score_vals)) |>
        filter(!is.na(Score)) |>
        count(Score)
    }
    p <- ggplot(plot_df, aes(x = Score, y = n)) +
      geom_col(fill = "#4e79a7", alpha = 0.85) +
      labs(title = "Response distribution", x = "Score", y = "Count") +
      theme_minimal(base_size = 12)
    plotly_basic(plotly::ggplotly(p))
  })

  selected_facets <- reactive({
    req(input$facet_cols, input$person_col, input$score_col)
    weight_col <- if (!is.null(input$weight_col) && input$weight_col != "(None)") input$weight_col else NULL
    setdiff(input$facet_cols, c(input$person_col, input$score_col, weight_col))
  })

  output$facet_summary_gt <- render_gt({
    req(raw_data())
    facets <- selected_facets()
    if (length(facets) == 0) {
      return(gt(tibble(Message = "Select facet columns to see summary.")))
    }
    df <- raw_data()
    summary_tbl <- purrr::map_dfr(facets, function(facet) {
      if (!facet %in% names(df)) {
        return(tibble(
          Facet = facet,
          Levels = NA_integer_,
          Missing = NA_integer_,
          N = nrow(df),
          ExampleLevels = "Facet not found"
        ))
      }
      vals <- df[[facet]]
      non_missing <- vals[!is.na(vals)]
      levels_ex <- unique(as.character(non_missing))
      tibble(
        Facet = facet,
        Levels = length(unique(non_missing)),
        Missing = sum(is.na(vals)),
        N = length(vals),
        ExampleLevels = paste(head(levels_ex, 6), collapse = ", ")
      )
    })
    gt(summary_tbl)
  })

  observeEvent(input$facet_select_all, {
    req(raw_data(), input$person_col, input$score_col)
    weight_col <- if (!is.null(input$weight_col) && input$weight_col != "(None)") input$weight_col else NULL
    cols <- setdiff(names(raw_data()), c(input$person_col, input$score_col, weight_col))
    updateSelectizeInput(session, "facet_cols", selected = cols, server = TRUE)
  })

  observeEvent(input$facet_clear, {
    updateSelectizeInput(session, "facet_cols", selected = character(0), server = TRUE)
  })

  observeEvent(input$facet_suggest, {
    req(raw_data(), input$person_col, input$score_col)
    cols <- names(raw_data())
    patterns <- c("rater", "task", "criterion", "criteria", "occasion", "group", "judge", "item")
    suggested <- cols[stringr::str_detect(tolower(cols), paste(patterns, collapse = "|"))]
    weight_col <- if (!is.null(input$weight_col) && input$weight_col != "(None)") input$weight_col else NULL
    suggested <- setdiff(suggested, c(input$person_col, input$score_col, weight_col))
    if (length(suggested) < 2) {
      suggested <- setdiff(cols, c(input$person_col, input$score_col, weight_col))
      suggested <- suggested[seq_len(min(2, length(suggested)))]
    }
    updateSelectizeInput(session, "facet_cols", selected = suggested, server = TRUE)
  })

  observeEvent(input$facet_apply_pattern, {
    req(raw_data(), input$person_col, input$score_col)
    pattern <- input$facet_pattern
    if (!nzchar(pattern)) return(NULL)
    cols <- names(raw_data())
    matches <- tryCatch({
      cols[stringr::str_detect(tolower(cols), tolower(pattern))]
    }, error = function(e) {
      character(0)
    })
    weight_col <- if (!is.null(input$weight_col) && input$weight_col != "(None)") input$weight_col else NULL
    matches <- setdiff(matches, c(input$person_col, input$score_col, weight_col))
    selected <- unique(c(input$facet_cols, matches))
    updateSelectizeInput(session, "facet_cols", selected = selected, server = TRUE)
  })

  observeEvent(input$fit_filter_all, {
    updateSelectizeInput(session, "fit_facet_filter", selected = "All", server = TRUE)
  })

  observeEvent(input$fit_filter_person, {
    updateSelectizeInput(session, "fit_facet_filter", selected = "Person", server = TRUE)
  })

  observeEvent(input$fit_filter_single, {
    req(input$fit_facet_single)
    updateSelectizeInput(session, "fit_facet_filter", selected = input$fit_facet_single, server = TRUE)
  })

  pair_options <- reactive({
    facets <- c("Person", selected_facets())
    if (length(facets) < 2) {
      return(list(labels = character(), combos = list()))
    }
    combos <- combn(facets, 2, simplify = FALSE)
    labels <- vapply(combos, function(x) paste(x, collapse = " x "), character(1))
    list(labels = labels, combos = combos)
  })

  output$interaction_pair_selector <- renderUI({
    opts <- pair_options()
    if (length(opts$labels) == 0) {
      return(helpText("Select 2+ facet columns to enable interaction bias checks."))
    }
    selectizeInput(
      "interaction_pairs",
      "Bias interaction pairs",
      choices = opts$labels,
      selected = opts$labels,
      multiple = TRUE,
      options = list(plugins = list("remove_button"))
    )
  })

  output$bias_pair_selector <- renderUI({
    opts <- pair_options()
    if (length(opts$labels) == 0) {
      return(helpText("Select 2+ facet columns to enable bias estimation."))
    }
    selectInput(
      "bias_pair",
      "Bias pair",
      choices = c("(None)", opts$labels),
      selected = "(None)"
    )
  })

  output$agreement_rater_selector <- renderUI({
    facets <- selected_facets()
    if (length(facets) == 0) {
      return(helpText("Select facet columns to enable agreement checks."))
    }
    suggested <- facets[stringr::str_detect(tolower(facets), "rater|judge|assessor|grader")]
    default <- if (length(suggested) > 0) suggested[1] else facets[1]
    selectInput(
      "agreement_facet",
      "Agreement facet",
      choices = facets,
      selected = default
    )
  })

  interaction_pairs <- reactive({
    opts <- pair_options()
    sel <- input$interaction_pairs
    if (is.null(sel)) return(NULL)
    if (length(sel) == 0) return(list())
    idx <- match(sel, opts$labels)
    opts$combos[idx[!is.na(idx)]]
  })

  results <- eventReactive(input$run, {
    req(raw_data())
    facets <- selected_facets()
    if (length(facets) < 2) {
      showModal(modalDialog(
        title = "Facet selection error",
        "Select at least 2 facet columns (total facets >= 3, including Person).",
        easyClose = TRUE
      ))
      return(NULL)
    }
    withProgress(message = "Estimating MFRM...", value = 0, {
      incProgress(0.2)
      weight_col <- if (!is.null(input$weight_col) && input$weight_col != "(None)") input$weight_col else NULL
      res <- tryCatch({
        mfrm_estimate(
          data = raw_data(),
          person_col = input$person_col,
          facet_cols = facets,
          score_col = input$score_col,
          rating_min = input$rating_min,
          rating_max = input$rating_max,
          weight_col = weight_col,
          keep_original = isTRUE(input$keep_original),
          model = input$model,
          method = input$method,
          step_facet = input$step_facet,
          anchor_df = anchor_df(),
          group_anchor_df = group_anchor_df(),
          noncenter_facet = input$noncenter_facet,
          dummy_facets = input$dummy_facets,
          positive_facets = input$positive_facets,
          quad_points = input$quad_points,
          maxit = input$maxit
        )
      }, error = function(e) {
        showModal(modalDialog(
          title = "Estimation error",
          paste("Details:", e$message),
          easyClose = TRUE
        ))
        NULL
      })
      incProgress(0.8)
      res
    })
  })

  params_data <- reactive({
    req(results())
    config <- results()$config
    sizes <- build_param_sizes(config)
    params <- expand_params(results()$opt$par, sizes, config)
    list(config = config, sizes = sizes, params = params)
  })

  diagnostics <- reactive({
    req(results())
    mfrm_diagnostics(
      results(),
      interaction_pairs = interaction_pairs(),
      top_n_interactions = input$top_n_interactions,
      whexact = isTRUE(input$whexact)
    )
  })

  bias_results <- reactive({
    req(results(), diagnostics())
    if (!isTRUE(input$bias_run)) return(NULL)
    if (is.null(input$bias_pair) || input$bias_pair == "(None)") return(NULL)
    parts <- strsplit(input$bias_pair, " x ")[[1]]
    if (length(parts) != 2) return(NULL)
    estimate_bias_interaction(
      results(),
      diagnostics(),
      parts[1],
      parts[2],
      max_abs = input$bias_max_abs,
      omit_extreme = isTRUE(input$bias_omit_extreme)
    )
  })

  agreement_results <- reactive({
    req(diagnostics(), results())
    calc_interrater_agreement(
      diagnostics()$obs,
      c("Person", results()$config$facet_names),
      input$agreement_facet,
      res = results()
    )
  })

  pca_results <- reactive({
    req(diagnostics(), results())
    obs_df <- diagnostics()$obs
    facet_names <- results()$config$facet_names
    list(
      overall = compute_pca_overall(obs_df, facet_names),
      by_facet = compute_pca_by_facet(obs_df, facet_names)
    )
  })

  output$pca_facet_selector <- renderUI({
    req(results())
    choices <- c("Overall", results()$config$facet_names)
    selectInput("pca_facet", "Facet for residual PCA", choices = choices, selected = "Overall")
  })

  pca_selected <- reactive({
    req(pca_results())
    if (is.null(input$pca_facet) || input$pca_facet == "Overall") {
      pca_results()$overall
    } else {
      pca_results()$by_facet[[input$pca_facet]]
    }
  })

  fit_plot_df <- reactive({
    req(diagnostics())
    df <- diagnostics()$measures |>
      filter(!is.na(Infit), !is.na(Outfit)) |>
      mutate(Level = as.character(Level))
    sel <- input$fit_facet_filter
    if (!is.null(sel) && length(sel) > 0 && !("All" %in% sel)) {
      df <- df |> filter(Facet %in% sel)
    }
    df
  })

  facets_report_tbls <- reactive({
    req(results(), diagnostics())
    calc_facets_report_tbls(
      results(),
      diagnostics(),
      totalscore = identical(input$report_totalscore, "yes"),
      umean = input$report_umean,
      uscale = input$report_uscale,
      udecimals = input$report_udecimals,
      omit_unobserved = identical(input$report_omit_unobserved, "yes"),
      xtreme = input$report_xtreme
    )
  })

  output$facets_report_tabs <- renderUI({
    req(facets_report_tbls())
    tbls <- facets_report_tbls()
    if (length(tbls) == 0) {
      return(tags$em("Facet reports will appear after estimation."))
    }
    tabs <- lapply(names(tbls), function(facet) {
      tabPanel(facet, value = facet)
    })
    do.call(tabsetPanel, c(list(id = "facet_report_tab", selected = names(tbls)[1]), tabs))
  })

  output$facets_report_gt <- render_gt({
    req(facets_report_tbls())
    tbls <- facets_report_tbls()
    if (length(tbls) == 0) {
      return(gt(tibble(Message = "Facet reports will appear after estimation.")))
    }
    facet <- input$facet_report_tab
    if (is.null(facet) || !facet %in% names(tbls)) {
      facet <- names(tbls)[1]
    }
    format_facets_report_gt(
      tbls[[facet]],
      facet,
      decimals = input$report_udecimals,
      totalscore = identical(input$report_totalscore, "yes")
    )
  })

  summary_tbl <- reactive({
    req(results(), diagnostics())
    results()$summary |>
      mutate(
        Infit = diagnostics()$overall_fit$Infit,
        Outfit = diagnostics()$overall_fit$Outfit,
        InfitZSTD = diagnostics()$overall_fit$InfitZSTD,
        OutfitZSTD = diagnostics()$overall_fit$OutfitZSTD
      )
  })

  settings_tbl <- reactive({
    req(results())
    prep <- results()$prep
    facets <- results()$config$facet_names
    step_facet <- ifelse(is.null(results()$config$step_facet), "Common", results()$config$step_facet)
    anchor_summary <- results()$config$anchor_summary
    anchor_info <- if (!is.null(anchor_summary) && nrow(anchor_summary) > 0) {
      anchor_summary |>
        mutate(
          Label = paste0(Facet, ": ", AnchoredLevels, " anchors, ", GroupAnchors, " group anchors")
        ) |>
        pull(Label) |>
        paste(collapse = " | ")
    } else {
      "None"
    }
    tibble(
      Field = c("Rating min", "Rating max", "Persons", "Facet names", "Step facet",
                "Noncenter facet", "Dummy facets", "Positive facets", "Weight column",
                "Anchor summary", "Keep original categories",
                "Totalscore", "Omit unobserved", "Xtreme", "Umean", "Uscale", "Udecimals", "WHEXACT"),
      Value = c(
        prep$rating_min,
        prep$rating_max,
        length(prep$levels$Person),
        paste(facets, collapse = ", "),
        step_facet,
        results()$config$noncenter_facet,
        ifelse(length(results()$config$dummy_facets) == 0, "None",
               paste(results()$config$dummy_facets, collapse = ", ")),
        ifelse(length(results()$config$positive_facets) == 0, "None",
               paste(results()$config$positive_facets, collapse = ", ")),
        ifelse(is.null(results()$config$weight_col), "None", results()$config$weight_col),
        anchor_info,
        ifelse(isTRUE(input$keep_original), "Yes", "No"),
        ifelse(identical(input$report_totalscore, "yes"), "Yes", "No"),
        ifelse(identical(input$report_omit_unobserved, "yes"), "Yes", "No"),
        input$report_xtreme,
        input$report_umean,
        input$report_uscale,
        input$report_udecimals,
        ifelse(isTRUE(input$whexact), "Yes", "No")
      )
    )
  })

  fit_overview_tbl <- reactive({
    req(diagnostics())
    diagnostics()$fit |>
      group_by(Facet) |>
      summarize(
        InfitMean = mean(Infit, na.rm = TRUE),
        OutfitMean = mean(Outfit, na.rm = TRUE),
        .groups = "drop"
      )
  })

  output$results_summary_gt <- render_gt({
    req(summary_tbl())
    summary_tbl() |>
      gt() |>
      fmt_number(columns = c(LogLik, AIC, BIC, Infit, Outfit, InfitZSTD, OutfitZSTD), decimals = 3)
  })

  output$results_reliability_snapshot_gt <- render_gt({
    req(diagnostics())
    diagnostics()$reliability |>
      gt() |>
      fmt_number(columns = c(SD, RMSE, Separation, Strata, Reliability, MeanInfit, MeanOutfit), decimals = 3)
  })

  output$results_reliability_plot <- renderPlotly({
    req(diagnostics())
    plot_df <- diagnostics()$reliability |>
      select(Facet, Separation, Reliability) |>
      pivot_longer(cols = c(Separation, Reliability), names_to = "Metric", values_to = "Value")
    p <- ggplot(plot_df, aes(x = Facet, y = Value, fill = Metric)) +
      geom_col(position = "dodge") +
      labs(title = "Reliability & separation by facet", x = NULL, y = "Value") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "bottom")
    plotly_compact_legend(plotly::ggplotly(p), bottom_margin = 90)
  })

  output$results_fit_overview_plot <- renderPlotly({
    req(fit_overview_tbl())
    plot_df <- fit_overview_tbl() |>
      pivot_longer(cols = c(InfitMean, OutfitMean), names_to = "Metric", values_to = "Value")
    p <- ggplot(plot_df, aes(x = Facet, y = Value, color = Metric, group = Metric)) +
      geom_line(linewidth = 0.9) +
      geom_point(size = 2.2) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "gray60") +
      labs(title = "Mean fit by facet", x = NULL, y = "Mean MNSQ") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "bottom")
    plotly_compact_legend(plotly::ggplotly(p), bottom_margin = 90)
  })

  output$results_flags_gt <- render_gt({
    req(results(), diagnostics(), fit_overview_tbl())
    flags <- list()
    if (!isTRUE(results()$summary$Converged)) {
      flags <- append(flags, list(tibble(Issue = "Convergence", Details = "Model did not converge")))
    }
    rel <- diagnostics()$reliability
    low_rel <- rel |>
      filter(!is.na(Reliability), Reliability < 0.8)
    if (nrow(low_rel) > 0) {
      flags <- append(flags, list(tibble(
        Issue = "Low reliability",
        Details = paste(low_rel$Facet, collapse = ", ")
      )))
    }
    low_sep <- rel |>
      filter(!is.na(Separation), Separation < 2)
    if (nrow(low_sep) > 0) {
      flags <- append(flags, list(tibble(
        Issue = "Low separation",
        Details = paste(low_sep$Facet, collapse = ", ")
      )))
    }
    fit_over <- fit_overview_tbl() |>
      filter(
        (!is.na(InfitMean) & (InfitMean < 0.5 | InfitMean > 1.5)) |
          (!is.na(OutfitMean) & (OutfitMean < 0.5 | OutfitMean > 1.5))
      )
    if (nrow(fit_over) > 0) {
      flags <- append(flags, list(tibble(
        Issue = "Mean fit outside 0.5-1.5",
        Details = paste(fit_over$Facet, collapse = ", ")
      )))
    }
    misfit_levels <- diagnostics()$measures |>
      mutate(AbsZ = pmax(abs(InfitZSTD), abs(OutfitZSTD), na.rm = TRUE)) |>
      filter(is.finite(AbsZ), AbsZ >= 2)
    if (nrow(misfit_levels) > 0) {
      flags <- append(flags, list(tibble(
        Issue = "Levels with |ZSTD| >= 2",
        Details = paste0(nrow(misfit_levels), " levels")
      )))
    }
    subset_tbl <- diagnostics()$subsets$summary
    if (!is.null(subset_tbl) && nrow(subset_tbl) > 1) {
      flags <- append(flags, list(tibble(
        Issue = "Disconnected subsets",
        Details = paste0(nrow(subset_tbl), " subsets detected")
      )))
    }
    flag_tbl <- if (length(flags) == 0) {
      tibble(Issue = "None", Details = "No major flags")
    } else {
      bind_rows(flags)
    }
    gt(flag_tbl)
  })

  output$cautionText <- renderText({
    req(results(), diagnostics())
    msgs <- c()
    if (!isTRUE(results()$summary$Converged)) {
      msgs <- c(msgs, "Model did not converge.")
    }
    if (results()$config$n_cat < 3) {
      msgs <- c(msgs, "Response categories are fewer than 3; estimates may be unstable.")
    }
    facet_levels <- vapply(results()$config$facet_levels, length, integer(1))
    small_facets <- names(facet_levels[facet_levels < 3])
    if (length(small_facets) > 0) {
      msgs <- c(msgs, paste0("Facets with <3 levels: ", paste(small_facets, collapse = ", ")))
    }
    misfit_levels <- diagnostics()$measures |>
      mutate(AbsZ = pmax(abs(InfitZSTD), abs(OutfitZSTD), na.rm = TRUE)) |>
      filter(is.finite(AbsZ), AbsZ >= 2)
    if (nrow(misfit_levels) > 0) {
      msgs <- c(msgs, paste0("Levels with |ZSTD| >= 2: ", nrow(misfit_levels)))
    }
    cat_warn <- category_warnings_text(
      calc_category_stats(diagnostics()$obs, results(), whexact = isTRUE(input$whexact)),
      calc_step_order(results()$steps)
    )
    if (!is.null(cat_warn) && !identical(cat_warn, "No major category warnings detected.")) {
      msgs <- c(msgs, cat_warn)
    }
    subset_tbl <- diagnostics()$subsets$summary
    if (!is.null(subset_tbl) && nrow(subset_tbl) > 1) {
      msgs <- c(msgs, paste0("Disconnected subsets detected: ", nrow(subset_tbl)))
    }
    if (length(msgs) == 0) {
      "No major cautions detected."
    } else {
      paste(msgs, collapse = "\n")
    }
  })

  output$results_interpretation <- renderUI({
    req(diagnostics(), fit_overview_tbl())
    rel <- diagnostics()$reliability
    low_rel <- rel |>
      filter(!is.na(Reliability), Reliability < 0.8)
    fit_over <- fit_overview_tbl() |>
      filter(
        (!is.na(InfitMean) & (InfitMean < 0.5 | InfitMean > 1.5)) |
          (!is.na(OutfitMean) & (OutfitMean < 0.5 | OutfitMean > 1.5))
      )
    tags$ul(
      tags$li(if (nrow(low_rel) > 0) {
        paste("Low reliability facets:", paste(low_rel$Facet, collapse = ", "))
      } else {
        "Reliability looks acceptable across facets."
      }),
      tags$li(if (nrow(fit_over) > 0) {
        paste("Check fit for facets:", paste(fit_over$Facet, collapse = ", "))
      } else {
        "Mean fit values are within the 0.5-1.5 band."
      }),
      tags$li("Inspect Steps/Thresholds and Category Curves to validate rating scale ordering."),
      tags$li("Use Bias/Interactions tabs to explore systematic residual patterns."),
      tags$li("Check Dimensionality tab for residual PCA diagnostics.")
    )
  })

  output$results_repro <- renderUI({
    req(results())
    tags$ul(
      tags$li(paste("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))),
      tags$li(paste("Model:", results()$summary$Model)),
      tags$li(paste("Method:", results()$summary$Method)),
      tags$li(paste("Facets:", paste(results()$config$facet_names, collapse = ", "))),
      tags$li(paste("Step facet:", ifelse(is.null(results()$config$step_facet), "Common", results()$config$step_facet))),
      tags$li(paste("Noncenter facet:", results()$config$noncenter_facet)),
      tags$li(paste("Dummy facets:", ifelse(length(results()$config$dummy_facets) == 0, "None",
                                            paste(results()$config$dummy_facets, collapse = ", ")))),
      tags$li(paste("Positive facets:", ifelse(length(results()$config$positive_facets) == 0, "None",
                                              paste(results()$config$positive_facets, collapse = ", ")))),
      tags$li(paste("Weight column:", ifelse(is.null(results()$config$weight_col), "None", results()$config$weight_col)))
    )
  })

  output$pca_scree_plot <- renderPlotly({
    req(pca_selected())
    pca_bundle <- pca_selected()
    if (is.null(pca_bundle) || is.null(pca_bundle$pca)) {
      return(
        plotly::plotly_empty(type = "scatter", mode = "lines") |>
          plotly::layout(title = list(text = "PCA results not available"))
      )
    }
    pca <- pca_bundle$pca
    eigenvalues <- pca$values[1:min(20, length(pca$values))]
    scree_data <- data.frame(Component = 1:length(eigenvalues), Eigenvalue = eigenvalues)
    gp <- ggplot(scree_data, aes(x = Component, y = Eigenvalue)) +
      geom_line(linewidth = 1.1) +
      geom_point(size = 2.4) +
      geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.8) +
      scale_x_continuous(breaks = 1:length(eigenvalues)) +
      labs(
        title = paste0("Scree plot (", ifelse(is.null(input$pca_facet), "Overall", input$pca_facet), ")"),
        subtitle = "Dashed line indicates eigenvalue = 1 (Kaiser criterion)",
        x = "Principal component",
        y = "Eigenvalue"
      ) +
      theme_minimal(base_size = 12)
    plotly_basic(plotly::ggplotly(gp, tooltip = c("x", "y")))
  })

  output$pca_eigen_gt <- render_gt({
    req(pca_selected())
    pca_bundle <- pca_selected()
    if (is.null(pca_bundle) || is.null(pca_bundle$pca)) {
      return(
        data.frame(Message = "PCA results not available. Check data variability.") |>
          gt() |>
          tab_header(title = "Residual PCA") |>
          tab_source_note(gt::md("Residual PCA needs at least two non-constant columns."))
      )
    }
    pca <- pca_bundle$pca
    n_components <- min(10, ncol(pca$loadings))
    var_pct <- if (!is.null(pca$Vaccounted)) {
      pca$Vaccounted[2, 1:n_components] * 100
    } else {
      pca$values[1:n_components] / sum(pca$values, na.rm = TRUE) * 100
    }
    cum_pct <- cumsum(var_pct)
    eigen_df <- data.frame(
      Component = paste0("PC", 1:n_components),
      Eigenvalue = pca$values[1:n_components],
      Variance_Pct = var_pct,
      Cumulative_Pct = cum_pct
    )
    eigen_df |>
      gt() |>
      tab_header(title = "Residual PCA", subtitle = "Eigenvalues and variance explained") |>
      fmt_number(columns = c(Eigenvalue), decimals = 3) |>
      fmt_number(columns = c(Variance_Pct, Cumulative_Pct), decimals = 1) |>
      cols_label(
        Component = "Component",
        Eigenvalue = "Eigenvalue",
        Variance_Pct = "% Variance",
        Cumulative_Pct = "Cumulative %"
      ) |>
      tab_style(
        style = cell_fill(color = "lightyellow"),
        locations = cells_body(columns = everything(), rows = Eigenvalue > 1)
      ) |>
      tab_source_note(gt::md("Highlighted rows have eigenvalues > 1 (Kaiser criterion)."))
  })

  output$pca_loadings_plot <- renderPlotly({
    req(pca_selected())
    pca_bundle <- pca_selected()
    if (is.null(pca_bundle) || is.null(pca_bundle$pca)) {
      return(
        plotly::plotly_empty(type = "bar") |>
          plotly::layout(title = list(text = "PCA loadings not available"))
      )
    }
    pca <- pca_bundle$pca
    loadings_mat <- as.data.frame(unclass(pca$loadings))
    if (ncol(loadings_mat) < 1) {
      return(plotly::plotly_empty(type = "bar") |>
               plotly::layout(title = list(text = "PCA loadings not available")))
    }
    loadings_pc1 <- data.frame(
      Item = rownames(loadings_mat),
      Loading = loadings_mat[, 1]
    ) |>
      arrange(desc(abs(Loading))) |>
      head(20)
    gp <- ggplot(loadings_pc1, aes(x = reorder(Item, Loading), y = Loading,
                                   text = sprintf("%s<br>Loading: %.3f", Item, Loading))) +
      geom_col(fill = ifelse(loadings_pc1$Loading > 0, "#4e79a7", "#e15759")) +
      coord_flip() +
      geom_hline(yintercept = 0, color = "gray40") +
      labs(
        title = "PC1 loadings (top 20)",
        x = NULL,
        y = "Loading"
      ) +
      theme_minimal(base_size = 12)
    plotly_basic(plotly::ggplotly(gp, tooltip = "text"))
  })

  output$pca_biplot <- renderPlotly({
    req(pca_selected())
    pca_bundle <- pca_selected()
    if (is.null(pca_bundle) || is.null(pca_bundle$pca)) {
      return(
        plotly::plotly_empty(type = "scatter", mode = "markers") |>
          plotly::layout(title = list(text = "PCA biplot not available"))
      )
    }
    pca <- pca_bundle$pca
    loadings_mat <- as.data.frame(unclass(pca$loadings))
    if (ncol(loadings_mat) < 2) {
      return(
        plotly::plotly_empty(type = "scatter", mode = "markers") |>
          plotly::layout(title = list(text = "Insufficient components for biplot"))
      )
    }
    loadings_df <- data.frame(
      Item = rownames(loadings_mat),
      PC1 = loadings_mat[, 1],
      PC2 = loadings_mat[, 2]
    )
    loadings_df$Distance <- sqrt(loadings_df$PC1^2 + loadings_df$PC2^2)
    top_loadings <- loadings_df |>
      arrange(desc(Distance)) |>
      head(if (requireNamespace("ggrepel", quietly = TRUE)) 30 else 20)
    var_pc1 <- if (!is.null(pca$Vaccounted)) pca$Vaccounted[2, 1] * 100 else NA_real_
    var_pc2 <- if (!is.null(pca$Vaccounted)) pca$Vaccounted[2, 2] * 100 else NA_real_
    use_repel <- requireNamespace("ggrepel", quietly = TRUE)
    gp <- ggplot(top_loadings, aes(x = PC1, y = PC2,
                                   label = Item,
                                   text = sprintf("%s<br>PC1: %.3f<br>PC2: %.3f", Item, PC1, PC2))) +
      geom_point(color = "#2f4b7c", size = 2) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
      labs(
        title = paste0("Residual PCA biplot (top ", nrow(top_loadings), " loadings)"),
        x = ifelse(is.na(var_pc1), "PC1", paste0("PC1 (", round(var_pc1, 1), "%)")),
        y = ifelse(is.na(var_pc2), "PC2", paste0("PC2 (", round(var_pc2, 1), "%)"))
      ) +
      theme_minimal(base_size = 12)
    if (nrow(top_loadings) > 0) {
      if (use_repel) {
        gp <- gp + ggrepel::geom_text_repel(
          size = 3,
          max.overlaps = Inf,
          min.segment.length = 0,
          box.padding = 0.35,
          point.padding = 0.2
        )
      } else {
        gp <- gp + geom_text(size = 3, hjust = 0, vjust = 0, check_overlap = TRUE)
      }
    }
    plotly_basic(plotly::ggplotly(gp, tooltip = "text"))
  })

  output$dimensionality_gt <- render_gt({
    req(pca_selected())
    pca_bundle <- pca_selected()
    if (is.null(pca_bundle) || is.null(pca_bundle$pca)) {
      return(
        data.frame(Message = "PCA analysis not completed. Unable to assess dimensionality.") |>
          gt() |>
          tab_header(title = "Unidimensionality Assessment") |>
          tab_source_note(gt::md("Run the analysis with sufficient observations to enable residual PCA."))
      )
    }
    pca <- pca_bundle$pca
    pc1_variance <- if (!is.null(pca$Vaccounted)) pca$Vaccounted[2, 1] * 100 else NA_real_
    eigenvalue_ratio <- if (length(pca$values) >= 2) pca$values[1] / pca$values[2] else NA_real_
    eigenvalues_above_1 <- sum(pca$values > 1, na.rm = TRUE)
    tibble(
      Metric = c("PC1 variance", "Eigenvalue ratio (PC1/PC2)", "Eigenvalues > 1", "Assessment"),
      Value = c(
        ifelse(is.na(pc1_variance), "N/A", paste0(round(pc1_variance, 1), "%")),
        ifelse(is.na(eigenvalue_ratio), "N/A", round(eigenvalue_ratio, 2)),
        eigenvalues_above_1,
        ifelse(!is.na(pc1_variance) && (pc1_variance > 20 || eigenvalue_ratio < 3),
               "Potential multidimensionality",
               "Acceptable unidimensionality")
      )
    ) |>
      gt() |>
      tab_header(
        title = "Unidimensionality Assessment",
        subtitle = "Based on PCA of standardized residuals"
      ) |>
      tab_source_note(gt::md("Heuristics: PC1 variance > 20% or PC1/PC2 < 3 may indicate secondary dimensions."))
  })

  output$summary_gt <- render_gt({
    req(summary_tbl())
    summary_tbl() |>
      gt() |>
      fmt_number(columns = c(LogLik, AIC, BIC, Infit, Outfit, InfitZSTD, OutfitZSTD), decimals = 3)
  })

  output$person_gt <- render_gt({
    req(diagnostics())
    tbl <- diagnostics()$measures |>
      filter(Facet == "Person") |>
      select(Level, Estimate, SE, CI_Lower, CI_Upper, Infit, Outfit, InfitZSTD, OutfitZSTD,
             MeanResidual, MeanStdResidual, N)
    gt(tbl) |>
      fmt_number(
        columns = c(Estimate, SE, CI_Lower, CI_Upper, Infit, Outfit, InfitZSTD, OutfitZSTD,
                    MeanResidual, MeanStdResidual),
        decimals = 3
      )
  })

  output$person_hist <- renderPlotly({
    req(results())
    tbl <- results()$facets$person
    p <- ggplot(tbl, aes(x = Estimate)) +
      geom_histogram(bins = 20, fill = "#4682b4", color = "white") +
      theme_minimal(base_size = 12) +
      labs(title = "Person estimates", x = "Theta", y = "Count")
    plotly_basic(plotly::ggplotly(p))
  })

  output$facets_gt <- render_gt({
    req(diagnostics())
    tbl <- diagnostics()$measures |>
      filter(Facet != "Person") |>
      select(Facet, Level, Estimate, SE, CI_Lower, CI_Upper, Infit, Outfit, InfitZSTD, OutfitZSTD,
             PTMEA, MeanResidual, MeanStdResidual, N)
    if (nrow(tbl) == 0) {
      return(gt(tibble(Message = "Facet rows are empty. Check facet selection and rerun.")))
    }
    gt(tbl) |>
      fmt_number(
        columns = c(Estimate, SE, CI_Lower, CI_Upper, Infit, Outfit, InfitZSTD, OutfitZSTD,
                    PTMEA, MeanResidual, MeanStdResidual),
        decimals = 3
      )
  })

  output$steps_gt <- render_gt({
    req(results())
    results()$steps |>
      gt() |>
      fmt_number(columns = Estimate, decimals = 3)
  })

  output$category_stats_gt <- render_gt({
    req(diagnostics(), results())
    tbl <- calc_category_stats(diagnostics()$obs, results(), whexact = isTRUE(input$whexact))
    if (nrow(tbl) == 0) {
      return(gt(tibble(Message = "No category statistics available.")))
    }
    gt(tbl) |>
      cols_label(
        ExpectedCount = "Expected Count",
        ExpectedPercent = "Expected %",
        DiffCount = "Obs - Exp",
        DiffPercent = "Obs - Exp %",
        AvgPersonMeasure = "Avg Measure",
        ExpectedAverage = "Expected Avg",
        InfitZSTD = "Infit ZSTD",
        OutfitZSTD = "Outfit ZSTD",
        MeanResidual = "Mean Residual"
      ) |>
      fmt_number(columns = c(Count, ExpectedCount, DiffCount), decimals = 0) |>
      fmt_number(columns = c(Percent, ExpectedPercent, DiffPercent,
                             AvgPersonMeasure, ExpectedAverage,
                             Infit, Outfit, InfitZSTD, OutfitZSTD,
                             MeanResidual),
                 decimals = 3) |>
      tab_style(
        style = list(cell_fill(color = "#fff4e5")),
        locations = cells_body(rows = LowCount)
      ) |>
      tab_style(
        style = list(cell_fill(color = "#ffe6e6")),
        locations = cells_body(rows = InfitFlag | OutfitFlag)
      ) |>
      tab_style(
        style = list(cell_fill(color = "#ffe6e6")),
        locations = cells_body(rows = ZSTDFlag)
      )
  })

  output$category_usage_plot <- renderPlotly({
    req(diagnostics(), results())
    tbl <- calc_category_stats(diagnostics()$obs, results(), whexact = isTRUE(input$whexact))
    if (nrow(tbl) == 0) return(NULL)
    plot_df <- tbl |>
      select(Category, Count, ExpectedCount) |>
      mutate(Category = factor(Category, levels = sort(unique(Category))))
    if (all(!is.finite(plot_df$ExpectedCount))) {
      p <- ggplot(plot_df, aes(x = Category, y = Count)) +
        geom_col(fill = "#4e79a7", alpha = 0.85) +
        labs(title = "Category usage (Observed)", x = "Category", y = "Count") +
        theme_minimal(base_size = 12)
      return(plotly_basic(plotly::ggplotly(p)))
    }
    plot_long <- plot_df |>
      pivot_longer(cols = c(Count, ExpectedCount), names_to = "Type", values_to = "Value") |>
      mutate(Type = recode(Type, Count = "Observed", ExpectedCount = "Expected"))
    p <- ggplot(plot_long, aes(x = Category, y = Value, fill = Type)) +
      geom_col(position = "dodge", alpha = 0.85) +
      labs(title = "Category usage: Observed vs Expected", x = "Category", y = "Count") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "bottom")
    plotly_compact_legend(plotly::ggplotly(p), bottom_margin = 90)
  })

  output$category_percent_plot <- renderPlotly({
    req(diagnostics(), results())
    tbl <- calc_category_stats(diagnostics()$obs, results(), whexact = isTRUE(input$whexact))
    if (nrow(tbl) == 0) return(NULL)
    plot_df <- tbl |>
      select(Category, Percent, ExpectedPercent) |>
      mutate(Category = factor(Category, levels = sort(unique(Category))))
    if (all(!is.finite(plot_df$ExpectedPercent))) {
      p <- ggplot(plot_df, aes(x = Category, y = Percent)) +
        geom_col(fill = "#59a14f", alpha = 0.85) +
        labs(title = "Category usage (%)", x = "Category", y = "Percent") +
        theme_minimal(base_size = 12)
      return(plotly_basic(plotly::ggplotly(p)))
    }
    plot_long <- plot_df |>
      pivot_longer(cols = c(Percent, ExpectedPercent), names_to = "Type", values_to = "Value") |>
      mutate(Type = recode(Type, Percent = "Observed %", ExpectedPercent = "Expected %"))
    p <- ggplot(plot_long, aes(x = Category, y = Value, fill = Type)) +
      geom_col(position = "dodge", alpha = 0.85) +
      labs(title = "Category usage: Observed vs Expected (%)", x = "Category", y = "Percent") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "bottom")
    plotly_compact_legend(plotly::ggplotly(p), bottom_margin = 90)
  })

  output$category_resid_plot <- renderPlotly({
    req(diagnostics())
    obs_df <- diagnostics()$obs
    if (nrow(obs_df) == 0) return(NULL)
    plot_df <- obs_df |>
      mutate(Category = factor(Observed, levels = sort(unique(Observed)))) |>
      filter(is.finite(StdResidual))
    if (nrow(plot_df) == 0) return(NULL)
    p <- ggplot(plot_df, aes(x = Category, y = StdResidual)) +
      geom_hline(yintercept = c(-2, 0, 2), linetype = c("dashed", "solid", "dashed"), color = "gray60") +
      geom_boxplot(outlier.shape = NA, alpha = 0.25, fill = "#76b7b2") +
      geom_jitter(width = 0.15, alpha = 0.4, size = 1, color = "#4e79a7") +
      labs(title = "Standardized residuals by category", x = "Category", y = "StdResidual") +
      theme_minimal(base_size = 12)
    plotly_basic(plotly::ggplotly(p))
  })

  output$category_fit_plot <- renderPlotly({
    req(diagnostics(), results())
    tbl <- calc_category_stats(diagnostics()$obs, results(), whexact = isTRUE(input$whexact))
    if (nrow(tbl) == 0) return(NULL)
    plot_df <- tbl |>
      select(Category, InfitZSTD, OutfitZSTD) |>
      pivot_longer(cols = c(InfitZSTD, OutfitZSTD), names_to = "Statistic", values_to = "ZSTD") |>
      filter(is.finite(ZSTD)) |>
      mutate(Category = factor(Category, levels = sort(unique(tbl$Category))))
    if (nrow(plot_df) == 0) return(NULL)
    p <- ggplot(plot_df, aes(x = Category, y = ZSTD, color = Statistic, group = Statistic)) +
      geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "gray60") +
      geom_point(size = 2.4) +
      geom_line(alpha = 0.4) +
      labs(title = "Category fit (ZSTD)", x = "Category", y = "ZSTD") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "bottom")
    plotly_compact_legend(plotly::ggplotly(p), bottom_margin = 90)
  })

  output$category_thresholds_gt <- render_gt({
    req(results())
    tbl <- calc_step_order(results()$steps)
    if (nrow(tbl) == 0) {
      return(gt(tibble(Message = "No threshold statistics available.")))
    }
    gt(tbl) |>
      fmt_number(columns = c(Estimate, Spacing), decimals = 3) |>
      tab_style(
        style = list(cell_fill(color = "#ffe6e6")),
        locations = cells_body(rows = !is.na(Ordered) & Ordered == FALSE)
      )
  })

  output$category_warnings_text <- renderText({
    req(diagnostics(), results())
    category_warnings_text(
      calc_category_stats(diagnostics()$obs, results(), whexact = isTRUE(input$whexact)),
      calc_step_order(results()$steps)
    )
  })

  output$constraint_summary_gt <- render_gt({
    req(results())
    summary_tbl <- results()$config$anchor_summary
    if (is.null(summary_tbl) || nrow(summary_tbl) == 0) {
      return(gt(tibble(Message = "No constraint summary available.")))
    }
    gt(summary_tbl)
  })

  output$constraint_anchors_gt <- render_gt({
    req(results())
    tbls <- extract_anchor_tables(results()$config)
    tbl <- tbls$anchors
    if (is.null(tbl) || nrow(tbl) == 0) {
      return(gt(tibble(Message = "No anchored levels.")))
    }
    gt(tbl) |>
      fmt_number(columns = Anchor, decimals = 3)
  })

  output$constraint_groups_gt <- render_gt({
    req(results())
    tbls <- extract_anchor_tables(results()$config)
    tbl <- tbls$groups
    if (is.null(tbl) || nrow(tbl) == 0) {
      return(gt(tibble(Message = "No group anchors.")))
    }
    gt(tbl) |>
      fmt_number(columns = GroupValue, decimals = 3)
  })

  output$subset_gt <- render_gt({
    req(diagnostics())
    tbl <- diagnostics()$subsets$summary
    if (is.null(tbl) || nrow(tbl) == 0) {
      return(gt(tibble(Message = "Subset analysis not available.")))
    }
    gt(tbl) |>
      fmt_number(columns = where(is.numeric), decimals = 0)
  })

  output$settings_gt <- render_gt({
    req(settings_tbl())
    settings_tbl() |>
      gt()
  })

  output$fit_overall_gt <- render_gt({
    req(diagnostics())
    diagnostics()$overall_fit |>
      gt() |>
      fmt_number(columns = c(Infit, Outfit, InfitZSTD, OutfitZSTD), decimals = 3)
  })

  output$fit_facet_gt <- render_gt({
    req(diagnostics())
    diagnostics()$fit |>
      gt() |>
      fmt_number(columns = c(Infit, Outfit, InfitZSTD, OutfitZSTD), decimals = 3)
  })

  output$reliability_gt <- render_gt({
    req(diagnostics())
    diagnostics()$reliability |>
      gt() |>
      fmt_number(columns = c(SD, RMSE, Separation, Strata, Reliability, MeanInfit, MeanOutfit), decimals = 3)
  })

  output$reliability_chisq_gt <- render_gt({
    req(diagnostics())
    chisq_tbl <- calc_facets_chisq(diagnostics()$measures)
    if (is.null(chisq_tbl) || nrow(chisq_tbl) == 0) {
      return(gt(tibble(Message = "Chi-square statistics not available.")))
    }
    chisq_tbl |>
      gt() |>
      cols_label(
        Levels = "Levels",
        MeanMeasure = "Mean",
        SD = "S.D.",
        FixedChiSq = "Fixed ChiSq",
        FixedDF = "Fixed df",
        FixedProb = "Fixed Prob.",
        RandomChiSq = "Random ChiSq",
        RandomDF = "Random df",
        RandomProb = "Random Prob."
      ) |>
      fmt_number(columns = c(MeanMeasure, SD, FixedChiSq, FixedProb, RandomChiSq, RandomProb), decimals = 3) |>
      fmt_number(columns = c(FixedDF, RandomDF, Levels), decimals = 0)
  })

  output$agreement_summary_gt <- render_gt({
    req(agreement_results())
    tbl <- agreement_results()$summary
    if (is.null(tbl) || nrow(tbl) == 0) {
      return(gt(tibble(Message = "Agreement summary not available.")))
    }
    gt(tbl) |>
      fmt_number(columns = c(Raters, Pairs, Contexts, TotalPairs, ExactAgreements, ExpectedAgreements), decimals = 0) |>
      fmt_percent(columns = c(ExactAgreement, ExpectedExactAgreement, AdjacentAgreement), decimals = 1) |>
      fmt_number(columns = c(MeanAbsDiff, MeanCorr), decimals = 3)
  })

  output$agreement_pairs_gt <- render_gt({
    req(agreement_results())
    tbl <- agreement_results()$pairs
    if (is.null(tbl) || nrow(tbl) == 0) {
      return(gt(tibble(Message = "Pairwise agreement not available.")))
    }
    gt(tbl) |>
      fmt_number(columns = c(N), decimals = 0) |>
      fmt_percent(columns = c(Exact, ExpectedExact, Adjacent), decimals = 1) |>
      fmt_number(columns = c(MeanDiff, MAD, Corr), decimals = 3)
  })

  output$bias_pairwise_ui <- renderUI({
    br <- bias_results()
    if (is.null(br) || is.null(br$table) || nrow(br$table) == 0) return(NULL)
    selectInput(
      "bias_target_facet",
      "Target facet (pairwise)",
      choices = c(br$facet_a, br$facet_b),
      selected = br$facet_a
    )
  })

  output$bias_facets_gt <- render_gt({
    br <- bias_results()
    if (is.null(br) || is.null(br$table) || nrow(br$table) == 0) {
      return(gt(tibble(Message = "Bias/interaction results not available.")))
    }
    gt(br$table) |>
      fmt_number(columns = where(is.numeric), decimals = 3)
  })

  output$bias_iteration_gt <- render_gt({
    br <- bias_results()
    if (is.null(br) || is.null(br$iteration) || nrow(br$iteration) == 0) {
      return(gt(tibble(Message = "Bias iteration report not available.")))
    }
    gt(br$iteration) |>
      fmt_number(columns = where(is.numeric), decimals = 3)
  })

  output$bias_summary_gt <- render_gt({
    br <- bias_results()
    if (is.null(br) || is.null(br$summary) || nrow(br$summary) == 0) {
      return(gt(tibble(Message = "Bias summary not available.")))
    }
    gt(br$summary) |>
      fmt_number(columns = where(is.numeric), decimals = 3)
  })

  output$bias_chi_gt <- render_gt({
    br <- bias_results()
    if (is.null(br) || is.null(br$chi_sq) || nrow(br$chi_sq) == 0) {
      return(gt(tibble(Message = "Bias chi-square not available.")))
    }
    gt(br$chi_sq) |>
      fmt_number(columns = where(is.numeric), decimals = 3)
  })

  output$bias_pairwise_gt <- render_gt({
    br <- bias_results()
    if (is.null(br) || is.null(br$table) || nrow(br$table) == 0) {
      return(gt(tibble(Message = "Pairwise bias report not available.")))
    }
    target <- input$bias_target_facet
    if (is.null(target) || !target %in% c(br$facet_a, br$facet_b)) {
      target <- br$facet_a
    }
    context <- ifelse(target == br$facet_a, br$facet_b, br$facet_a)
    pair_tbl <- calc_bias_pairwise(br$table, target, context)
    if (is.null(pair_tbl) || nrow(pair_tbl) == 0) {
      return(gt(tibble(Message = "Pairwise bias report not available.")))
    }
    gt(pair_tbl) |>
      fmt_number(columns = where(is.numeric), decimals = 3)
  })

  output$bias_gt <- render_gt({
    req(diagnostics())
    diagnostics()$bias |>
      select(Facet, Level, N, ObservedAverage, ExpectedAverage, Bias,
             MeanAbsStdResidual, ChiSq, ChiDf, ChiP,
             SE_Residual, t_Residual, p_Residual,
             MeanStdResidual, SE_StdResidual, t_StdResidual, p_StdResidual, DF) |>
      gt() |>
      cols_label(
        ObservedAverage = "Observed Avg",
        ExpectedAverage = "Expected Avg",
        Bias = "Bias",
        MeanAbsStdResidual = "Mean |ZSTD|",
        ChiSq = "ChiSq",
        ChiDf = "Chi df",
        ChiP = "Chi p",
        SE_Residual = "SE (Bias)",
        t_Residual = "t (Bias)",
        p_Residual = "p (Bias)",
        MeanStdResidual = "Mean ZSTD",
        SE_StdResidual = "SE (ZSTD)",
        t_StdResidual = "t (ZSTD)",
        p_StdResidual = "p (ZSTD)"
      ) |>
      fmt_number(columns = c(ObservedAverage, ExpectedAverage, Bias, SE_Residual, t_Residual, p_Residual,
                             MeanAbsStdResidual, ChiSq, ChiDf, ChiP,
                             MeanStdResidual, SE_StdResidual, t_StdResidual, p_StdResidual, DF),
                 decimals = 3)
  })

  output$bias_interaction_gt <- render_gt({
    req(diagnostics())
    diagnostics()$interactions |>
      select(Pair, Level, N, ObservedAverage, ExpectedAverage, Bias,
             MeanAbsStdResidual, ChiSq, ChiDf, ChiP,
             SE_Residual, t_Residual, p_Residual,
             MeanStdResidual, SE_StdResidual, t_StdResidual, p_StdResidual, DF) |>
      gt() |>
      cols_label(
        ObservedAverage = "Observed Avg",
        ExpectedAverage = "Expected Avg",
        Bias = "Bias",
        MeanAbsStdResidual = "Mean |ZSTD|",
        ChiSq = "ChiSq",
        ChiDf = "Chi df",
        ChiP = "Chi p",
        SE_Residual = "SE (Bias)",
        t_Residual = "t (Bias)",
        p_Residual = "p (Bias)",
        MeanStdResidual = "Mean ZSTD",
        SE_StdResidual = "SE (ZSTD)",
        t_StdResidual = "t (ZSTD)",
        p_StdResidual = "p (ZSTD)"
      ) |>
      fmt_number(columns = c(ObservedAverage, ExpectedAverage, Bias, SE_Residual, t_Residual, p_Residual,
                             MeanAbsStdResidual, ChiSq, ChiDf, ChiP,
                             MeanStdResidual, SE_StdResidual, t_StdResidual, p_StdResidual, DF),
                 decimals = 3)
  })

  output$download_sample_tsv <- downloadHandler(
    filename = function() paste0("mfrm_sample_", Sys.Date(), ".tsv"),
    content = function(file) {
      readr::write_tsv(download_sample_data, file)
    }
  )

  output$download_sample_csv <- downloadHandler(
    filename = function() paste0("mfrm_sample_", Sys.Date(), ".csv"),
    content = function(file) {
      readr::write_csv(download_sample_data, file)
    }
  )

  output$download_summary_quick <- downloadHandler(
    filename = function() paste0("mfrm_summary_", Sys.Date(), ".csv"),
    content = function(file) {
      req(summary_tbl())
      readr::write_csv(summary_tbl(), file)
    }
  )

  output$download_reliability_quick <- downloadHandler(
    filename = function() paste0("mfrm_reliability_", Sys.Date(), ".csv"),
    content = function(file) {
      req(diagnostics())
      readr::write_csv(diagnostics()$reliability, file)
    }
  )

  output$download_fit_quick <- downloadHandler(
    filename = function() paste0("mfrm_fit_", Sys.Date(), ".csv"),
    content = function(file) {
      req(diagnostics())
      readr::write_csv(diagnostics()$fit, file)
    }
  )

  output$download_steps_quick <- downloadHandler(
    filename = function() paste0("mfrm_thresholds_", Sys.Date(), ".csv"),
    content = function(file) {
      req(results())
      readr::write_csv(results()$steps, file)
    }
  )

  output$download_measures <- downloadHandler(
    filename = function() paste0("mfrm_measures_", Sys.Date(), ".csv"),
    content = function(file) {
      req(diagnostics())
      measures_download <- diagnostics()$measures |>
        mutate(PTMEA = ifelse(Facet == "Person", NA_real_, PTMEA))
      readr::write_csv(measures_download, file, na = "")
    }
  )

  output$download_scorefile <- downloadHandler(
    filename = function() paste0("mfrm_scorefile_", Sys.Date(), ".csv"),
    content = function(file) {
      req(results())
      scorefile <- compute_scorefile(results())
      readr::write_csv(scorefile, file, na = "")
    }
  )

  output$download_residuals <- downloadHandler(
    filename = function() paste0("mfrm_residuals_", Sys.Date(), ".csv"),
    content = function(file) {
      req(results())
      residuals <- compute_residual_file(results())
      readr::write_csv(residuals, file, na = "")
    }
  )

  output$download_zip <- downloadHandler(
    filename = function() paste0("mfrm_results_", Sys.Date(), ".zip"),
    content = function(file) {
      req(results(), diagnostics(), summary_tbl(), settings_tbl())
      tmpdir <- tempfile("mfrm_results_")
      dir.create(tmpdir)
      files <- character()
      write_tbl <- function(df, name, na_value = "NA") {
        if (is.null(df)) df <- tibble()
        path <- file.path(tmpdir, name)
        readr::write_csv(df, path, na = na_value)
        files <<- c(files, path)
      }
      write_tbl(summary_tbl(), "summary.csv")
      measures_download <- diagnostics()$measures |>
        mutate(PTMEA = ifelse(Facet == "Person", NA_real_, PTMEA))
      write_tbl(measures_download, "measures.csv", na_value = "")
      write_tbl(compute_scorefile(results()), "scorefile.csv", na_value = "")
      write_tbl(compute_residual_file(results()), "residuals.csv", na_value = "")
      write_tbl(diagnostics()$overall_fit, "fit_overall.csv")
      write_tbl(diagnostics()$fit, "fit_facet.csv")
      write_tbl(diagnostics()$reliability, "reliability.csv")
      write_tbl(calc_facets_chisq(diagnostics()$measures), "reliability_chisq.csv")
      write_tbl(diagnostics()$bias, "bias.csv")
      write_tbl(diagnostics()$interactions, "interactions.csv")
      br <- bias_results()
      if (!is.null(br) && !is.null(br$table)) {
        write_tbl(br$table, "bias_facets.csv")
        if (!is.null(br$summary)) write_tbl(br$summary, "bias_summary.csv")
        if (!is.null(br$chi_sq)) write_tbl(br$chi_sq, "bias_chi_sq.csv")
        if (!is.null(br$iteration)) write_tbl(br$iteration, "bias_iteration.csv")
        pair_tbl <- calc_bias_pairwise(br$table, br$facet_a, br$facet_b)
        if (!is.null(pair_tbl) && nrow(pair_tbl) > 0) {
          write_tbl(pair_tbl, "bias_pairwise.csv")
        }
      }
      agreement_tbls <- agreement_results()
      if (!is.null(agreement_tbls)) {
        write_tbl(agreement_tbls$summary, "agreement_summary.csv")
        write_tbl(agreement_tbls$pairs, "agreement_pairs.csv")
      }
      write_tbl(calc_category_stats(diagnostics()$obs, results(), whexact = isTRUE(input$whexact)), "category_stats.csv")
      write_tbl(calc_step_order(results()$steps), "category_thresholds.csv")
      write_tbl(diagnostics()$subsets$summary, "subset_summary.csv")
      constraint_tbls <- extract_anchor_tables(results()$config)
      write_tbl(results()$config$anchor_summary, "constraint_summary.csv")
      write_tbl(constraint_tbls$anchors, "constraint_anchors.csv")
      write_tbl(constraint_tbls$groups, "constraint_groups.csv")
      write_tbl(results()$steps, "steps.csv")
      write_tbl(settings_tbl(), "settings.csv")
      write_tbl(diagnostics()$obs, "obs_table.csv")
      facet_reports <- facets_report_tbls()
      if (!is.null(facet_reports) && length(facet_reports) > 0) {
        for (facet in names(facet_reports)) {
          safe_name <- stringr::str_replace_all(facet, "[^A-Za-z0-9]", "_")
          write_tbl(facet_reports[[facet]], paste0("facet_report_", safe_name, ".csv"))
        }
      }
      utils::zip(zipfile = file, files = files, flags = "-j")
    }
  )

  output$wright_plot <- renderPlotly({
    req(results(), params_data())
    res <- results()
    params <- params_data()$params
    theta_hat <- if (res$config$method == "JMLE") {
      params$theta
    } else {
      res$facets$person$Estimate
    }
    person_df <- tibble(panel = "Persons", Estimate = theta_hat)
    facet_df <- purrr::map_dfr(names(params$facets), function(facet) {
      tibble(panel = "Facets", Type = facet, Estimate = params$facets[[facet]])
    })
    p <- ggplot() +
      geom_density(
        data = person_df,
        aes(x = Estimate, y = ..density..),
        fill = "#6baed6",
        alpha = 0.5
      ) +
      geom_point(
        data = facet_df,
        aes(x = Estimate, y = Type, color = Type),
        position = position_jitter(height = 0.12),
        size = 2,
        alpha = 0.85
      ) +
      facet_grid(. ~ panel, scales = "free_y", space = "free") +
      scale_color_brewer(palette = "Dark2") +
      labs(
        title = "Wright map (side-by-side)",
        x = "Logit scale",
        y = NULL
      ) +
      theme_minimal(base_size = 12) +
      theme(legend.position = "bottom")
    plotly_legend_top(plotly::ggplotly(p), top_margin = 90)
  })

  output$pathway_plot <- renderPlotly({
    req(diagnostics())
    plot_df <- diagnostics()$measures |>
      filter(!is.na(InfitZSTD), !is.na(Estimate)) |>
      mutate(
        Facet = ifelse(Facet == "Person", "Person", Facet),
        LabelFull = Level,
        Label = truncate_label(Level, width = 16),
        Tooltip = paste0(
          "Facet: ", Facet,
          "<br>Level: ", LabelFull,
          "<br>Measure: ", round(Estimate, 3),
          "<br>Infit ZSTD: ", round(InfitZSTD, 2)
        )
      )
    label_df <- plot_df |>
      filter(abs(InfitZSTD) >= 2) |>
      arrange(desc(abs(InfitZSTD))) |>
      slice_head(n = 20)
    use_repel <- requireNamespace("ggrepel", quietly = TRUE)
    p <- ggplot(plot_df, aes(x = InfitZSTD, y = Estimate, color = Facet, text = Tooltip)) +
      geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "gray60") +
      geom_errorbar(aes(ymin = Estimate - SE, ymax = Estimate + SE), alpha = 0.3, width = 0) +
      geom_point(alpha = 0.8, size = 2) +
      scale_color_brewer(palette = "Dark2") +
      labs(
        title = "Pathway map (measure vs. infit t)",
        x = "Infit ZSTD",
        y = "Measure"
      ) +
      theme_minimal(base_size = 12) +
      theme(legend.position = "bottom")
    if (nrow(label_df) > 0) {
      if (use_repel) {
        p <- p + ggrepel::geom_text_repel(
          data = label_df,
          aes(label = Label),
          size = 3,
          max.overlaps = Inf,
          min.segment.length = 0,
          box.padding = 0.35,
          point.padding = 0.2,
          show.legend = FALSE
        )
      } else {
        p <- p + geom_text(
          data = label_df,
          aes(label = Label),
          hjust = -0.1,
          size = 3,
          check_overlap = TRUE,
          show.legend = FALSE
        )
      }
    }
    plotly_compact_legend(plotly::ggplotly(p, tooltip = "text"), bottom_margin = 100)
  })

  output$facet_plot <- renderPlotly({
    req(results(), params_data())
    params <- params_data()$params
    facet_df <- purrr::map_dfr(names(params$facets), function(facet) {
      tibble(Type = facet, Estimate = params$facets[[facet]])
    })
    p <- ggplot(facet_df, aes(x = Type, y = Estimate, color = Type)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.2) +
      geom_jitter(width = 0.15, size = 2, alpha = 0.85) +
      scale_color_brewer(palette = "Dark2") +
      labs(title = "Facet estimate distribution", x = NULL, y = "Estimate") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "none")
    plotly_basic(plotly::ggplotly(p))
  })

  output$step_plot <- renderPlotly({
    req(results(), params_data())
    res <- results()
    step_tbl <- res$steps |>
      mutate(StepIndex = as.integer(stringr::str_extract(Step, "\\d+")))
    if (!"StepFacet" %in% names(step_tbl)) {
      step_tbl <- step_tbl |> mutate(StepFacet = "Common")
    }
    show_legend <- length(unique(step_tbl$StepFacet)) > 1
    p <- ggplot(step_tbl, aes(x = StepIndex, y = Estimate, group = StepFacet, color = StepFacet)) +
      geom_line(linewidth = 0.9) +
      geom_point(size = 2.2) +
      labs(
        title = "Step / threshold estimates",
        x = "Step",
        y = "Estimate"
      ) +
      theme_minimal(base_size = 12) +
      theme(legend.position = if (show_legend) "bottom" else "none") +
      guides(color = guide_legend(title = NULL))
    p <- plotly::ggplotly(p)
    if (show_legend) {
      p <- plotly_compact_legend(p, bottom_margin = 100)
    } else {
      p <- plotly::style(p, showlegend = FALSE)
      p <- plotly_basic(p)
    }
    p
  })

  output$category_plot <- renderPlotly({
    req(results(), params_data())
    res <- results()
    params <- params_data()$params
    theta_grid <- seq(-4, 4, length.out = 81)
    score_vals <- seq(res$prep$rating_min, res$prep$rating_max)

    if (res$config$model == "RSM") {
      step_cum <- c(0, cumsum(params$steps))
      probs <- category_prob_rsm(theta_grid, step_cum)
      df <- as_tibble(probs)
      names(df) <- paste0("Score_", score_vals)
      plot_df <- df |>
        mutate(Theta = theta_grid) |>
        pivot_longer(cols = starts_with("Score_"), names_to = "Score", values_to = "Prob") |>
        mutate(Score = factor(Score, levels = paste0("Score_", score_vals)))
      p <- ggplot(plot_df, aes(x = Theta, y = Prob, color = Score)) +
        geom_line(linewidth = 0.9) +
        labs(title = "Category probability curves (RSM)", x = "Theta", y = "Probability") +
        theme_minimal(base_size = 12) +
        theme(legend.position = "bottom")
      plotly_compact_legend(plotly::ggplotly(p), bottom_margin = 100)
    } else {
      step_facet <- res$config$step_facet
      if (is.null(step_facet) || !step_facet %in% names(res$prep$levels)) {
        return(NULL)
      }
      step_levels <- res$prep$levels[[step_facet]]
      crit_sel <- input$pcm_curve_criterion
      if (is.null(crit_sel) || !crit_sel %in% step_levels) {
        crit_sel <- step_levels[1]
      }
      step_cum_mat <- t(apply(params$steps_mat, 1, function(x) c(0, cumsum(x))))
      plot_list <- lapply(seq_len(nrow(step_cum_mat)), function(i) {
        probs <- category_prob_rsm(theta_grid, step_cum_mat[i, ])
        df <- as_tibble(probs)
        names(df) <- paste0("Score_", score_vals)
        df |>
          mutate(Theta = theta_grid, StepFacet = step_levels[i]) |>
          pivot_longer(cols = starts_with("Score_"), names_to = "Score", values_to = "Prob")
      })
      plot_df <- bind_rows(plot_list) |>
        mutate(Score = factor(Score, levels = paste0("Score_", score_vals))) |>
        filter(StepFacet == crit_sel)
      p <- ggplot(plot_df, aes(x = Theta, y = Prob, color = Score)) +
        geom_line(linewidth = 0.85) +
        labs(
          title = glue("Category probability curves (PCM: {crit_sel})"),
          x = "Theta",
          y = "Probability"
        ) +
        theme_minimal(base_size = 12) +
        theme(legend.position = "bottom")
      plotly_compact_legend(plotly::ggplotly(p), bottom_margin = 100)
    }
  })

  output$obs_exp_plot <- renderPlotly({
    req(results())
    obs_exp <- expected_score_table(results())
    plot_df <- obs_exp |>
      mutate(bin = ntile(Expected, 10)) |>
      group_by(bin) |>
      summarise(
        Expected = mean(Expected),
        Observed = mean(Observed),
        n = n(),
        .groups = "drop"
      )
    p <- ggplot(plot_df, aes(x = Expected, y = Observed, size = n)) +
      geom_point(color = "#1f77b4", alpha = 0.85) +
      geom_line(color = "#1f77b4", alpha = 0.6) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
      scale_size_continuous(range = c(2, 6)) +
      labs(
        title = "Observed vs expected score",
        x = "Expected score",
        y = "Observed score",
        size = "N (bin)"
      ) +
      theme_minimal(base_size = 12) +
      theme(legend.position = "bottom")
    plotly_compact_legend(plotly::ggplotly(p), bottom_margin = 100)
  })

  output$fit_scatter_plot <- renderPlotly({
    req(fit_plot_df())
    df <- fit_plot_df()
    if (!"N" %in% names(df)) df$N <- 1
    p <- ggplot(df, aes(
      x = Infit,
      y = Outfit,
      color = Facet,
      text = paste0(
        "Facet: ", Facet,
        "<br>Level: ", Level,
        "<br>Infit: ", round(Infit, 3),
        "<br>Outfit: ", round(Outfit, 3),
        "<br>N: ", N
      )
    )) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "gray60") +
      geom_hline(yintercept = 1, linetype = "dashed", color = "gray60") +
      geom_point(aes(size = N), alpha = 0.8) +
      scale_size_continuous(range = c(1.5, 6)) +
      labs(
        title = "Fit scatter (Infit vs Outfit)",
        x = "Infit MNSQ",
        y = "Outfit MNSQ",
        size = "N"
      ) +
      theme_minimal(base_size = 12) +
      theme(legend.position = "bottom") +
      guides(size = "none")
    plotly_compact_legend(plotly::ggplotly(p, tooltip = "text"), bottom_margin = 100)
  })

  output$fit_zstd_plot <- renderPlotly({
    req(fit_plot_df())
    plot_df <- fit_plot_df() |>
      select(Facet, InfitZSTD, OutfitZSTD) |>
      pivot_longer(cols = c(InfitZSTD, OutfitZSTD), names_to = "Statistic", values_to = "ZSTD") |>
      filter(is.finite(ZSTD))
    p <- ggplot(plot_df, aes(x = ZSTD, fill = Statistic)) +
      geom_density(alpha = 0.4) +
      geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "gray60") +
      facet_wrap(~ Facet, scales = "free_y") +
      labs(title = "ZSTD distribution by facet", x = "ZSTD", y = "Density") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "bottom")
    plotly_compact_legend(plotly::ggplotly(p), bottom_margin = 100)
  })

  output$misfit_bar_plot <- renderPlotly({
    req(fit_plot_df(), input$top_n_misfit, input$misfit_threshold, input$misfit_compare)
    df <- fit_plot_df() |>
      mutate(
        AbsZSTD = pmax(abs(InfitZSTD), abs(OutfitZSTD), na.rm = TRUE),
        AbsZSTD = ifelse(is.infinite(AbsZSTD), NA_real_, AbsZSTD)
      )
    df <- df |>
      filter(is.finite(AbsZSTD)) |>
      (\(x) {
        if (input$misfit_compare == "gte") {
          filter(x, AbsZSTD >= input$misfit_threshold)
        } else {
          filter(x, AbsZSTD <= input$misfit_threshold)
        }
      })() |>
      arrange(desc(AbsZSTD)) |>
      slice_head(n = input$top_n_misfit) |>
      mutate(
        LabelFull = paste(Facet, Level, sep = ": "),
        Label = truncate_label(LabelFull, width = 32),
        Tooltip = paste0(
          LabelFull,
          "<br>|ZSTD|: ", round(AbsZSTD, 2),
          "<br>Infit ZSTD: ", round(InfitZSTD, 2),
          "<br>Outfit ZSTD: ", round(OutfitZSTD, 2)
        )
      )
    if (nrow(df) == 0) return(NULL)
    p <- ggplot(df, aes(x = reorder(Label, AbsZSTD), y = AbsZSTD, fill = Facet, text = Tooltip)) +
      geom_col(alpha = 0.85) +
      coord_flip() +
      geom_hline(yintercept = 2, linetype = "dashed", color = "gray60") +
      labs(title = "Top misfit levels (max |ZSTD|)", x = NULL, y = "|ZSTD|") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "bottom")
    plotly_compact_legend(plotly::ggplotly(p, tooltip = "text"), bottom_margin = 110)
  })
}

shinyApp(ui, server)
