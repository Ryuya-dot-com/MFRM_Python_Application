args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("usage: run_immer_cmle_parity.R wide.csv mapping.csv output.csv RSM|PCM")
}
if (!requireNamespace("immer", quietly = TRUE)) {
  stop("R package immer is required")
}

wide <- utils::read.csv(args[[1]], stringsAsFactors = FALSE, check.names = FALSE)
mapping <- utils::read.csv(args[[2]], stringsAsFactors = FALSE, check.names = FALSE)
model <- toupper(args[[4]])
if (!model %in% c("RSM", "PCM")) stop("model must be RSM or PCM")
person_col <- "Person"
resp <- wide[, setdiff(names(wide), person_col), drop = FALSE]
mapping <- mapping[match(names(resp), mapping$VirtualUnit), , drop = FALSE]
if (anyNA(mapping$VirtualUnit)) {
  stop("virtual-unit mapping is incomplete")
}

prep <- immer:::lpcm_data_prep(resp, weights = NULL, a = NULL)
pars_info <- prep$pars_info
rater_levels <- sort(unique(mapping$Rater))
criterion_levels <- sort(unique(mapping$Criterion))
step_levels <- if (model == "RSM") "__shared__" else criterion_levels
n_steps <- max(prep$maxK)

sum_zero <- function(n) {
  if (n <= 1L) return(matrix(0, nrow = n, ncol = 0L))
  out <- matrix(0, nrow = n, ncol = n - 1L)
  out[seq_len(n - 1L), ] <- diag(n - 1L)
  out[n, ] <- -1
  out
}

rater_contrast <- sum_zero(length(rater_levels))
criterion_contrast <- sum_zero(length(criterion_levels))
step_contrast <- sum_zero(n_steps)
parameter_names <- c(
  paste0("facet:Rater:free:", rater_levels[seq_len(ncol(rater_contrast))]),
  paste0("facet:Criterion:free:", criterion_levels[seq_len(ncol(criterion_contrast))]),
  unlist(lapply(step_levels, function(level) {
    paste0("step:", level, ":free:", seq_len(ncol(step_contrast)))
  }), use.names = FALSE)
)
W <- matrix(0, nrow = nrow(pars_info), ncol = length(parameter_names))
colnames(W) <- parameter_names
rownames(W) <- rownames(pars_info)

rater_offset <- 0L
criterion_offset <- ncol(rater_contrast)
step_offset <- criterion_offset + ncol(criterion_contrast)
for (row_index in seq_len(nrow(pars_info))) {
  unit <- as.character(pars_info$item[row_index])
  category <- as.integer(pars_info$cat[row_index])
  meta <- mapping[mapping$VirtualUnit == unit, , drop = FALSE]
  rater_index <- match(meta$Rater[[1]], rater_levels)
  criterion_index <- match(meta$Criterion[[1]], criterion_levels)
  if (ncol(rater_contrast)) {
    cols <- rater_offset + seq_len(ncol(rater_contrast))
    W[row_index, cols] <- category * rater_contrast[rater_index, ]
  }
  if (ncol(criterion_contrast)) {
    cols <- criterion_offset + seq_len(ncol(criterion_contrast))
    W[row_index, cols] <- category * criterion_contrast[criterion_index, ]
  }
  if (ncol(step_contrast)) {
    step_level_index <- if (model == "RSM") 1L else criterion_index
    level_offset <- (step_level_index - 1L) * ncol(step_contrast)
    cols <- step_offset + level_offset + seq_len(ncol(step_contrast))
    W[row_index, cols] <- colSums(step_contrast[seq_len(category), , drop = FALSE])
  }
}

fit <- immer::immer_cml(
  resp,
  W = W,
  par_init = rep(0, ncol(W)),
  use_rcpp = FALSE,
  control = list(maxit = 2000L, reltol = 1e-12)
)
out <- data.frame(
  Parameter = colnames(W),
  Estimate = as.numeric(fit$coefficients),
  ConditionalLogLik = rep(as.numeric(fit$loglike), length(fit$coefficients)),
  Convergence = rep(as.integer(fit$result_optim$convergence), length(fit$coefficients)),
  stringsAsFactors = FALSE
)
utils::write.csv(out, args[[3]], row.names = FALSE)
