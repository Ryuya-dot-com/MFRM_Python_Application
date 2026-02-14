# Core MFRM computation functions extracted from app.R (non-UI section)
# Source: ../../app.R
# NOTE: This file intentionally keeps function names aligned with the app implementation.
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
  if (n == 1) {
    return(list(nodes = 0, weights = 1))
  }
  i <- seq_len(n - 1)
  a <- rep(0, n)
  b <- sqrt(i / 2)
  jmat <- matrix(0, nrow = n, ncol = n)
  diag(jmat) <- a
  jmat[cbind(seq_len(n - 1), seq_len(n - 1) + 1)] <- b
  jmat[cbind(seq_len(n - 1) + 1, seq_len(n - 1))] <- b
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
  mnsq <- as.numeric(mnsq)
  df <- as.numeric(df)

  if (length(df) == 1L && length(mnsq) > 1L) {
    df <- rep(df, length(mnsq))
  } else if (length(mnsq) == 1L && length(df) > 1L) {
    mnsq <- rep(mnsq, length(df))
  }

  n <- min(length(mnsq), length(df))
  if (n == 0L) return(numeric(0))

  out <- rep(NA_real_, n)
  m <- mnsq[seq_len(n)]
  d <- df[seq_len(n)]
  ok <- is.finite(m) & is.finite(d) & (d > 0)

  if (isTRUE(whexact)) {
    out[ok] <- (m[ok] - 1) * sqrt(d[ok] / 2)
  } else {
    out[ok] <- (m[ok]^(1 / 3) - (1 - 2 / (9 * d[ok]))) / sqrt(2 / (9 * d[ok]))
  }
  out
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

resolve_pcm_step_facet <- function(model, step_facet, facet_names) {
  if (model != "PCM") return(NULL)
  resolved <- if (is.null(step_facet)) facet_names[1] else step_facet
  if (!resolved %in% facet_names) {
    stop("Selected step facet is not in the facet list.")
  }
  resolved
}

sanitize_noncenter_facet <- function(noncenter_facet, facet_names) {
  if (!is.null(noncenter_facet) && noncenter_facet %in% c("Person", facet_names)) {
    return(noncenter_facet)
  }
  "Person"
}

sanitize_dummy_facets <- function(dummy_facets, facet_names) {
  if (is.null(dummy_facets)) return(character(0))
  intersect(dummy_facets, c("Person", facet_names))
}

build_facet_signs <- function(facet_names, positive_facets = character(0)) {
  positives <- intersect(positive_facets, facet_names)
  signs <- setNames(ifelse(facet_names %in% positives, 1, -1), facet_names)
  list(signs = signs, positive_facets = positives)
}

build_estimation_config <- function(prep,
                                    model,
                                    method,
                                    step_facet,
                                    weight_col,
                                    facet_signs,
                                    positive_facets,
                                    noncenter_facet,
                                    dummy_facets,
                                    anchor_df,
                                    group_anchor_df) {
  config <- list(
    model = model,
    method = method,
    n_person = length(prep$levels$Person),
    n_cat = prep$rating_max - prep$rating_min + 1,
    facet_names = prep$facet_names,
    facet_levels = prep$levels[prep$facet_names],
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

  list(
    config = config,
    sizes = build_param_sizes(config)
  )
}

build_initial_param_vector <- function(config, sizes) {
  n_cat <- config$n_cat
  step_init <- if (n_cat > 1) {
    seq(-1, 1, length.out = n_cat - 1)
  } else {
    numeric(0)
  }

  facet_starts <- unlist(lapply(config$facet_names, function(f) rep(0, sizes[[f]])))
  c(
    rep(0, sizes$theta),
    facet_starts,
    if (config$model == "RSM") {
      step_init
    } else {
      rep(step_init, length(config$facet_levels[[config$step_facet]]))
    }
  )
}

run_mfrm_optimization <- function(start,
                                  method,
                                  idx,
                                  config,
                                  sizes,
                                  quad_points,
                                  maxit,
                                  reltol) {
  control <- list(maxit = maxit, reltol = reltol)
  if (method == "JMLE") {
    return(optim(
      start,
      fn = mfrm_loglik_jmle,
      idx = idx,
      config = config,
      sizes = sizes,
      method = "BFGS",
      control = control
    ))
  }

  quad <- gauss_hermite_normal(quad_points)
  optim(
    start,
    fn = mfrm_loglik_mml,
    idx = idx,
    config = config,
    sizes = sizes,
    quad = quad,
    method = "BFGS",
    control = control
  )
}

build_person_table <- function(method, idx, config, params, prep, quad_points) {
  if (method == "MML") {
    quad <- gauss_hermite_normal(quad_points)
    return(
      compute_person_eap(idx, config, params, quad) |>
        mutate(Person = prep$levels$Person) |>
        select(Person, Estimate, SD)
    )
  }

  tibble(
    Person = prep$levels$Person,
    Estimate = params$theta
  )
}

build_other_facet_table <- function(config, prep, params) {
  facet_tbls <- lapply(config$facet_names, function(facet) {
    tibble(Level = prep$levels[[facet]], Estimate = params$facets[[facet]]) |>
      mutate(Facet = facet, .before = 1)
  })
  bind_rows(facet_tbls)
}

build_step_table <- function(config, prep, params) {
  if (config$model == "RSM") {
    return(
      tibble(
        Step = paste0("Step_", seq_len(config$n_cat - 1)),
        Estimate = params$steps
      )
    )
  }

  expand_grid(
    StepFacet = prep$levels[[config$step_facet]],
    Step = paste0("Step_", seq_len(config$n_cat - 1))
  ) |>
    mutate(Estimate = as.vector(t(params$steps_mat)))
}

build_estimation_summary <- function(model, method, prep, config, sizes, opt) {
  k_params <- sum(unlist(sizes))
  loglik <- -opt$value
  n_obs <- if (!is.null(config$weight_col) && "Weight" %in% names(prep$data)) {
    sum(prep$data$Weight, na.rm = TRUE)
  } else {
    nrow(prep$data)
  }
  aic <- 2 * k_params - 2 * loglik
  bic <- log(n_obs) * k_params - 2 * loglik

  tibble(
    Model = model,
    Method = method,
    N = n_obs,
    Persons = config$n_person,
    Facets = length(config$facet_names),
    Categories = config$n_cat,
    LogLik = loglik,
    AIC = aic,
    BIC = bic,
    Converged = opt$convergence == 0,
    Iterations = opt$counts[["function"]]
  )
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
  # Stage 1: Normalize model options and input data.
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

  # Stage 2: Resolve facet-level modeling choices.
  step_facet <- resolve_pcm_step_facet(model, step_facet, prep$facet_names)
  noncenter_facet <- sanitize_noncenter_facet(noncenter_facet, prep$facet_names)
  dummy_facets <- sanitize_dummy_facets(dummy_facets, prep$facet_names)
  sign_info <- build_facet_signs(prep$facet_names, positive_facets)

  # Stage 3: Build reusable structures for optimization.
  idx <- build_indices(prep, step_facet = step_facet)
  cfg <- build_estimation_config(
    prep = prep,
    model = model,
    method = method,
    step_facet = step_facet,
    weight_col = weight_col,
    facet_signs = sign_info$signs,
    positive_facets = sign_info$positive_facets,
    noncenter_facet = noncenter_facet,
    dummy_facets = dummy_facets,
    anchor_df = anchor_df,
    group_anchor_df = group_anchor_df
  )
  config <- cfg$config
  sizes <- cfg$sizes
  start <- build_initial_param_vector(config, sizes)

  # Stage 4: Optimize model parameters.
  opt <- run_mfrm_optimization(
    start = start,
    method = method,
    idx = idx,
    config = config,
    sizes = sizes,
    quad_points = quad_points,
    maxit = maxit,
    reltol = reltol
  )

  # Stage 5: Build human-readable output tables.
  params <- expand_params(opt$par, sizes, config)
  person_tbl <- build_person_table(method, idx, config, params, prep, quad_points)
  facet_tbl <- build_other_facet_table(config, prep, params)
  step_tbl <- build_step_table(config, prep, params)
  summary_tbl <- build_estimation_summary(model, method, prep, config, sizes, opt)

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
  # Stage 1: Validate required inputs and extract core model objects.
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

  # Stage 2: Build step structures used for fair-average calculations.
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

  # Stage 3: Detect extreme-only levels and cache row-level flags.
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

  # Stage 4: Build one FACETS-style table per facet.
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
      # FACETS convention: remove rows with multiple extreme flags,
      # but keep rows where the focal facet is the only extreme flag.
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

    # Convert any xtreme-adjusted eta back to facet-specific measure units.
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

compute_pca_overall <- function(obs_df, facet_names, max_factors = 10L) {
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

  n_factors <- max(1, min(as.integer(max_factors), ncol(cor_matrix) - 1, nrow(cor_matrix) - 1))
  pca_result <- tryCatch(psych::principal(cor_matrix, nfactors = n_factors, rotate = "none"), error = function(e) NULL)
  list(pca = pca_result, residual_matrix = residual_matrix_wide, cor_matrix = cor_matrix)
}

compute_pca_by_facet <- function(obs_df, facet_names, max_factors = 10L) {
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

    n_factors <- max(1, min(as.integer(max_factors), ncol(cor_mat) - 1, nrow(cor_mat) - 1))
    pca_obj <- tryCatch(psych::principal(cor_mat, nfactors = n_factors, rotate = "none"), error = function(e) NULL)
    out[[facet]] <- list(pca = pca_obj, cor_matrix = cor_mat, residual_matrix = wide)
  }
  out
}

mfrm_diagnostics <- function(res,
                             interaction_pairs = NULL,
                             top_n_interactions = 20,
                             whexact = FALSE,
                             residual_pca = c("none", "overall", "facet", "both"),
                             pca_max_factors = 10L) {
  residual_pca <- match.arg(tolower(residual_pca), c("none", "overall", "facet", "both"))
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

  pca_overall <- NULL
  pca_by_facet <- NULL
  if (residual_pca %in% c("overall", "both")) {
    pca_overall <- compute_pca_overall(
      obs_df = obs_df,
      facet_names = res$config$facet_names,
      max_factors = pca_max_factors
    )
  }
  if (residual_pca %in% c("facet", "both")) {
    pca_by_facet <- compute_pca_by_facet(
      obs_df = obs_df,
      facet_names = res$config$facet_names,
      max_factors = pca_max_factors
    )
  }

  list(
    obs = obs_df,
    facet_names = res$config$facet_names,
    overall_fit = overall_fit,
    measures = measures,
    fit = fit_tbl,
    reliability = reliability_tbl,
    bias = bias_tbl,
    interactions = interaction_tbl,
    subsets = subset_tbls,
    residual_pca_mode = residual_pca,
    residual_pca_overall = pca_overall,
    residual_pca_by_facet = pca_by_facet
  )
}
