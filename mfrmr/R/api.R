#' Fit a many-facet Rasch model with a flexible number of facets
#'
#' This is the package entry point. It wraps `mfrm_estimate()` and defaults to
#' `method = "MML"`. Any number of facet columns can be supplied via `facets`.
#'
#' @param data A data.frame with one row per observation.
#' @param person Column name for the person (character scalar).
#' @param facets Character vector of facet column names.
#' @param score Column name for observed category score.
#' @param rating_min Optional minimum category value.
#' @param rating_max Optional maximum category value.
#' @param weight Optional weight column name.
#' @param keep_original Keep original category values.
#' @param model `"RSM"` or `"PCM"`.
#' @param method `"MML"` (default) or `"JML"` / `"JMLE"`.
#' @param step_facet Step facet for PCM.
#' @param anchors Optional anchor table.
#' @param group_anchors Optional group-anchor table.
#' @param noncenter_facet One facet to leave non-centered.
#' @param dummy_facets Facets to fix at zero.
#' @param positive_facets Facets with positive orientation.
#' @param anchor_policy How to handle anchor-audit issues: `"warn"` (default),
#'   `"error"`, or `"silent"`.
#' @param min_common_anchors Minimum anchored levels per linking facet used in
#'   anchor-audit recommendations.
#' @param min_obs_per_element Minimum weighted observations per facet level used
#'   in anchor-audit recommendations.
#' @param min_obs_per_category Minimum weighted observations per score category
#'   used in anchor-audit recommendations.
#' @param quad_points Quadrature points for MML.
#' @param maxit Maximum optimizer iterations.
#' @param reltol Optimization tolerance.
#'
#' @details
#' Data must be in **long format** (one row per observed rating event).
#'
#' @section Input requirements:
#' Minimum required columns are:
#' - person identifier (`person`)
#' - one or more facet identifiers (`facets`)
#' - observed score (`score`)
#'
#' Scores are treated as ordered categories.
#' If your observed categories do not start at 0, set `rating_min`/`rating_max`
#' explicitly to avoid unintended recoding assumptions.
#'
#' Supported model/estimation combinations:
#' - `model = "RSM"` with `method = "MML"` or `"JML"/"JMLE"`
#' - `model = "PCM"` with a designated `step_facet` (defaults to first facet)
#'
#' Anchor inputs are optional:
#' - `anchors` should contain facet/level/fixed-value information.
#' - `group_anchors` should contain facet/level/group/group-value information.
#' Both are normalized internally, so column names can be flexible
#' (`facet`, `level`, `anchor`, `group`, `groupvalue`, etc.).
#'
#' Anchor audit behavior:
#' - `fit_mfrm()` runs an internal anchor audit.
#' - invalid rows are removed before estimation.
#' - duplicate rows keep the last occurrence for each key.
#' - `anchor_policy` controls whether detected issues are warned, treated as
#'   errors, or kept silent.
#'
#' Facet sign orientation:
#' - facets listed in `positive_facets` are treated as `+1`
#' - all other facets are treated as `-1`
#' This affects interpretation of reported facet measures.
#'
#' @section Interpreting output:
#' A typical first-pass read is:
#' 1. `fit$summary` for convergence and global fit indicators.
#' 2. `summary(fit)` for human-readable overviews.
#' 3. `diagnose_mfrm(fit)` for element-level fit, separation, and warning tables.
#'
#' @section Typical workflow:
#' 1. Fit the model with `fit_mfrm(...)`.
#' 2. Validate convergence and scale structure with `summary(fit)`.
#' 3. Run [diagnose_mfrm()] and proceed to reporting with [build_apa_outputs()].
#'
#' @return
#' An object of class `mfrm_fit` (named list) with:
#' - `summary`: one-row model summary (`LogLik`, `AIC`, `BIC`, convergence)
#' - `facets$person`: person estimates (`Estimate`; plus `SD` for MML)
#' - `facets$others`: facet-level estimates for each facet
#' - `steps`: estimated threshold/step parameters
#' - `config`: resolved model configuration used for estimation
#'   (includes `config$anchor_audit`)
#' - `prep`: preprocessed data/level metadata
#' - `opt`: raw optimizer result from [stats::optim()]
#'
#' @seealso [diagnose_mfrm()], [estimate_bias()], [build_apa_outputs()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#'
#' fit <- fit_mfrm(
#'   data = toy,
#'   person = "Person",
#'   facets = c("Rater", "Criterion"),
#'   score = "Score",
#'   method = "JML",
#'   model = "RSM",
#'   maxit = 25
#' )
#' fit$summary
#' s_fit <- summary(fit)
#' s_fit$overview[, c("Model", "Method", "Converged")]
#' p_fit <- plot(fit, draw = FALSE)
#' class(p_fit)
#'
#' # MML is the default:
#' fit_mml <- fit_mfrm(
#'   data = toy,
#'   person = "Person",
#'   facets = c("Rater", "Criterion"),
#'   score = "Score",
#'   model = "RSM",
#'   quad_points = 7,
#'   maxit = 25
#' )
#' summary(fit_mml)
#' @export
fit_mfrm <- function(data,
                     person,
                     facets,
                     score,
                     rating_min = NULL,
                     rating_max = NULL,
                     weight = NULL,
                     keep_original = FALSE,
                     model = c("RSM", "PCM"),
                     method = c("MML", "JML", "JMLE"),
                     step_facet = NULL,
                     anchors = NULL,
                     group_anchors = NULL,
                     noncenter_facet = "Person",
                     dummy_facets = NULL,
                     positive_facets = NULL,
                     anchor_policy = c("warn", "error", "silent"),
                     min_common_anchors = 5L,
                     min_obs_per_element = 30,
                     min_obs_per_category = 10,
                     quad_points = 15,
                     maxit = 400,
                     reltol = 1e-6) {
  # -- input validation --
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame. Got: ", class(data)[1], ". ",
         "Convert with as.data.frame() if needed.", call. = FALSE)
  }
  if (nrow(data) == 0) {
    stop("`data` has zero rows. ",
         "Supply a data.frame with at least one observation.", call. = FALSE)
  }
  if (!is.character(person) || length(person) != 1 || !nzchar(person)) {
    stop("`person` must be a single non-empty character string ",
         "naming the person column.", call. = FALSE)
  }
  if (!is.character(facets) || length(facets) == 0) {
    stop("`facets` must be a character vector of one or more facet column names.",
         call. = FALSE)
  }
  if (!is.character(score) || length(score) != 1 || !nzchar(score)) {
    stop("`score` must be a single non-empty character string ",
         "naming the score column.", call. = FALSE)
  }
  if (!is.null(weight) && (!is.character(weight) || length(weight) != 1)) {
    stop("`weight` must be NULL or a single character string ",
         "naming the weight column.", call. = FALSE)
  }
  if (!is.numeric(maxit) || length(maxit) != 1 || maxit < 1) {
    stop("`maxit` must be a positive integer. Got: ", deparse(maxit), ".",
         call. = FALSE)
  }
  if (!is.numeric(reltol) || length(reltol) != 1 || reltol <= 0) {
    stop("`reltol` must be a positive number. Got: ", deparse(reltol), ".",
         call. = FALSE)
  }
  if (!is.numeric(quad_points) || length(quad_points) != 1 || quad_points < 1) {
    stop("`quad_points` must be a positive integer. Got: ", deparse(quad_points), ".",
         call. = FALSE)
  }

  model <- toupper(match.arg(model))
  method <- toupper(match.arg(method))
  method <- ifelse(method == "JML", "JMLE", method)
  anchor_policy <- tolower(match.arg(anchor_policy))

  anchor_audit <- audit_mfrm_anchors(
    data = data,
    person = person,
    facets = facets,
    score = score,
    rating_min = rating_min,
    rating_max = rating_max,
    weight = weight,
    keep_original = keep_original,
    anchors = anchors,
    group_anchors = group_anchors,
    min_common_anchors = min_common_anchors,
    min_obs_per_element = min_obs_per_element,
    min_obs_per_category = min_obs_per_category,
    noncenter_facet = noncenter_facet,
    dummy_facets = dummy_facets
  )

  anchors <- anchor_audit$anchors
  group_anchors <- anchor_audit$group_anchors

  issue_counts <- anchor_audit$issue_counts
  issue_total <- if (is.null(issue_counts) || nrow(issue_counts) == 0) 0L else sum(issue_counts$N, na.rm = TRUE)
  if (issue_total > 0) {
    msg <- format_anchor_audit_message(anchor_audit)
    if (anchor_policy == "error") {
      stop(msg, call. = FALSE)
    } else if (anchor_policy == "warn") {
      warning(msg, call. = FALSE)
    }
  }

  fit <- mfrm_estimate(
    data = data,
    person_col = person,
    facet_cols = facets,
    score_col = score,
    rating_min = rating_min,
    rating_max = rating_max,
    weight_col = weight,
    keep_original = keep_original,
    model = model,
    method = method,
    step_facet = step_facet,
    anchor_df = anchors,
    group_anchor_df = group_anchors,
    noncenter_facet = noncenter_facet,
    dummy_facets = dummy_facets,
    positive_facets = positive_facets,
    quad_points = quad_points,
    maxit = maxit,
    reltol = reltol
  )

  fit$config$anchor_audit <- anchor_audit

  class(fit) <- c("mfrm_fit", class(fit))
  fit
}

format_anchor_audit_message <- function(anchor_audit) {
  if (is.null(anchor_audit$issue_counts) || nrow(anchor_audit$issue_counts) == 0) {
    return("Anchor audit detected no issues.")
  }
  nonzero <- anchor_audit$issue_counts |>
    dplyr::filter(.data$N > 0)
  if (nrow(nonzero) == 0) {
    return("Anchor audit detected no issues.")
  }
  labels <- paste0(nonzero$Issue, "=", nonzero$N)
  paste0(
    "Anchor audit detected ", sum(nonzero$N), " issue row(s): ",
    paste(labels, collapse = "; "),
    ". Invalid rows were removed; duplicate keys keep the last row."
  )
}

summarize_linkage_by_facet <- function(df, facet) {
  by_level <- df |>
    dplyr::group_by(.data[[facet]]) |>
    dplyr::summarize(
      PersonsPerLevel = dplyr::n_distinct(.data$Person),
      Observations = dplyr::n(),
      WeightedN = sum(.data$Weight, na.rm = TRUE),
      .groups = "drop"
    )

  by_person <- df |>
    dplyr::group_by(.data$Person) |>
    dplyr::summarize(
      LevelsPerPerson = dplyr::n_distinct(.data[[facet]]),
      .groups = "drop"
    )

  tibble::tibble(
    Facet = facet,
    Levels = nrow(by_level),
    MinPersonsPerLevel = min(by_level$PersonsPerLevel, na.rm = TRUE),
    MedianPersonsPerLevel = stats::median(by_level$PersonsPerLevel, na.rm = TRUE),
    MinLevelsPerPerson = min(by_person$LevelsPerPerson, na.rm = TRUE),
    MedianLevelsPerPerson = stats::median(by_person$LevelsPerPerson, na.rm = TRUE)
  )
}

#' Summarize MFRM input data (TAM-style descriptive snapshot)
#'
#' @param data A data.frame in long format (one row per rating event).
#' @param person Column name for person IDs.
#' @param facets Character vector of facet column names.
#' @param score Column name for observed score.
#' @param weight Optional weight/frequency column name.
#' @param rating_min Optional minimum category value.
#' @param rating_max Optional maximum category value.
#' @param keep_original Keep original category values.
#' @param include_person_facet If `TRUE`, include person-level rows in
#'   `facet_level_summary`.
#' @param include_agreement If `TRUE`, include an observed-score inter-rater
#'   agreement bundle (summary/pairs/settings) in the output.
#' @param rater_facet Optional rater facet name used for agreement summaries.
#'   If `NULL`, inferred from facet names.
#' @param context_facets Optional facets used to define matched contexts for
#'   agreement. If `NULL`, all remaining facets (including `Person`) are used.
#' @param agreement_top_n Optional maximum number of agreement pair rows.
#'
#' @details
#' This function provides a compact descriptive bundle similar to the
#' pre-fit summaries commonly checked in TAM workflows:
#' sample size, score distribution, per-facet coverage, and linkage counts.
#' `psych::describe()` is used for numeric descriptives of score and weight.
#'
#' @section Interpreting output:
#' Recommended order:
#' - `overview`: confirms sample size, facet count, and category span.
#' - `missing_by_column`: identifies immediate data-quality risks.
#' - `score_distribution`: checks sparse/unused score categories.
#' - `facet_level_summary` and `linkage_summary`: checks per-level support and
#'   person-facet connectivity.
#' - `agreement`: optional observed inter-rater consistency summary.
#'
#' @section Typical workflow:
#' 1. Run `describe_mfrm_data()` on long-format input.
#' 2. Review `summary(ds)` and `plot(ds, ...)`.
#' 3. Resolve missingness/sparsity issues before [fit_mfrm()].
#'
#' @return A list of class `mfrm_data_description` with:
#' - `overview`: one-row run-level summary
#' - `missing_by_column`: missing counts in selected input columns
#' - `score_descriptives`: output from [psych::describe()] for score
#' - `weight_descriptives`: output from [psych::describe()] for weight
#' - `score_distribution`: weighted and raw score frequencies
#' - `facet_level_summary`: per-level usage and score summaries
#' - `linkage_summary`: person-facet connectivity diagnostics
#' - `agreement`: observed-score inter-rater agreement bundle
#'
#' @seealso [fit_mfrm()], [audit_mfrm_anchors()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' ds <- describe_mfrm_data(
#'   data = toy,
#'   person = "Person",
#'   facets = c("Rater", "Criterion"),
#'   score = "Score"
#' )
#' s_ds <- summary(ds)
#' s_ds$overview
#' p_ds <- plot(ds, draw = FALSE)
#' class(p_ds)
#' @export
describe_mfrm_data <- function(data,
                               person,
                               facets,
                               score,
                               weight = NULL,
                               rating_min = NULL,
                               rating_max = NULL,
                               keep_original = FALSE,
                               include_person_facet = FALSE,
                               include_agreement = TRUE,
                               rater_facet = NULL,
                               context_facets = NULL,
                               agreement_top_n = NULL) {
  prep <- prepare_mfrm_data(
    data = data,
    person_col = person,
    facet_cols = facets,
    score_col = score,
    rating_min = rating_min,
    rating_max = rating_max,
    weight_col = weight,
    keep_original = keep_original
  )

  df <- prep$data |>
    dplyr::mutate(
      Person = as.character(.data$Person),
      dplyr::across(dplyr::all_of(prep$facet_names), as.character)
    )

  selected_cols <- unique(c(person, facets, score, if (!is.null(weight)) weight))
  missing_by_column <- tibble::tibble(
    Column = selected_cols,
    Missing = vapply(selected_cols, function(col) sum(is.na(data[[col]])), integer(1))
  )

  score_desc <- psych::describe(df$Score, fast = TRUE)
  weight_desc <- psych::describe(df$Weight, fast = TRUE)

  total_weight <- sum(df$Weight, na.rm = TRUE)
  score_distribution <- df |>
    dplyr::group_by(.data$Score) |>
    dplyr::summarize(
      RawN = dplyr::n(),
      WeightedN = sum(.data$Weight, na.rm = TRUE),
      Percent = ifelse(total_weight > 0, 100 * .data$WeightedN / total_weight, NA_real_),
      .groups = "drop"
    ) |>
    dplyr::arrange(.data$Score)

  report_facets <- prep$facet_names
  if (isTRUE(include_person_facet)) {
    report_facets <- c("Person", report_facets)
  }

  facet_level_summary <- purrr::map_dfr(report_facets, function(facet) {
    df |>
      dplyr::group_by(.data[[facet]]) |>
      dplyr::summarize(
        RawN = dplyr::n(),
        WeightedN = sum(.data$Weight, na.rm = TRUE),
        MeanScore = weighted_mean(.data$Score, .data$Weight),
        SDScore = stats::sd(.data$Score, na.rm = TRUE),
        MinScore = min(.data$Score, na.rm = TRUE),
        MaxScore = max(.data$Score, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::rename(Level = dplyr::all_of(facet)) |>
      dplyr::mutate(
        Facet = facet,
        Level = as.character(.data$Level),
        .before = 1
      )
  })

  linkage_summary <- if (length(prep$facet_names) == 0) {
    tibble::tibble()
  } else {
    purrr::map_dfr(prep$facet_names, function(facet) summarize_linkage_by_facet(df, facet))
  }

  agreement_bundle <- list(
    summary = data.frame(),
    pairs = data.frame(),
    settings = list(
      included = FALSE,
      rater_facet = NA_character_,
      context_facets = character(0),
      expected_exact_from_model = FALSE,
      top_n = if (is.null(agreement_top_n)) NA_integer_ else max(1L, as.integer(agreement_top_n))
    )
  )

  if (isTRUE(include_agreement) && length(prep$facet_names) > 0) {
    known_facets <- c("Person", prep$facet_names)
    if (is.null(rater_facet) || !nzchar(as.character(rater_facet[1]))) {
      rater_facet <- infer_default_rater_facet(prep$facet_names)
    } else {
      rater_facet <- as.character(rater_facet[1])
    }
    if (is.null(rater_facet) || !rater_facet %in% known_facets) {
      stop("`rater_facet` must match one of: ", paste(known_facets, collapse = ", "))
    }
    if (identical(rater_facet, "Person")) {
      stop("`rater_facet = 'Person'` is not supported. Use a non-person facet.")
    }

    if (is.null(context_facets)) {
      facet_cols <- known_facets
      resolved_context <- setdiff(facet_cols, rater_facet)
    } else {
      context_facets <- unique(as.character(context_facets))
      unknown <- setdiff(context_facets, known_facets)
      if (length(unknown) > 0) {
        stop("Unknown `context_facets`: ", paste(unknown, collapse = ", "))
      }
      resolved_context <- setdiff(context_facets, rater_facet)
      if (length(resolved_context) == 0) {
        stop("`context_facets` must include at least one facet different from `rater_facet`.")
      }
      facet_cols <- c(rater_facet, resolved_context)
    }

    obs_agreement <- df |>
      dplyr::select(dplyr::all_of(unique(c("Person", prep$facet_names, "Score", "Weight")))) |>
      dplyr::rename(Observed = "Score")

    agreement <- calc_interrater_agreement(
      obs_df = obs_agreement,
      facet_cols = facet_cols,
      rater_facet = rater_facet,
      res = NULL
    )
    agreement_pairs <- as.data.frame(agreement$pairs, stringsAsFactors = FALSE)
    if (!is.null(agreement_top_n) && nrow(agreement_pairs) > 0) {
      agreement_pairs <- agreement_pairs |>
        dplyr::slice_head(n = max(1L, as.integer(agreement_top_n)))
    }

    agreement_bundle <- list(
      summary = as.data.frame(agreement$summary, stringsAsFactors = FALSE),
      pairs = agreement_pairs,
      settings = list(
        included = TRUE,
        rater_facet = rater_facet,
        context_facets = resolved_context,
        expected_exact_from_model = FALSE,
        top_n = if (is.null(agreement_top_n)) NA_integer_ else max(1L, as.integer(agreement_top_n))
      )
    )
  }

  overview <- tibble::tibble(
    Observations = nrow(df),
    TotalWeight = total_weight,
    Persons = length(prep$levels$Person),
    Facets = length(prep$facet_names),
    Categories = prep$rating_max - prep$rating_min + 1,
    RatingMin = prep$rating_min,
    RatingMax = prep$rating_max
  )

  out <- list(
    overview = overview,
    missing_by_column = missing_by_column,
    score_descriptives = score_desc,
    weight_descriptives = weight_desc,
    score_distribution = score_distribution,
    facet_level_summary = facet_level_summary,
    linkage_summary = linkage_summary,
    agreement = agreement_bundle
  )
  class(out) <- c("mfrm_data_description", class(out))
  out
}

#' @export
print.mfrm_data_description <- function(x, ...) {
  cat("mfrm data description\n")
  if (!is.null(x$overview) && nrow(x$overview) > 0) {
    print(x$overview, row.names = FALSE)
  }
  if (!is.null(x$score_distribution) && nrow(x$score_distribution) > 0) {
    cat("\nScore distribution\n")
    print(x$score_distribution, row.names = FALSE)
  }
  if (!is.null(x$agreement$summary) && nrow(x$agreement$summary) > 0) {
    cat("\nInter-rater agreement (observed)\n")
    print(x$agreement$summary, row.names = FALSE)
  }
  invisible(x)
}

#' Summarize a data-description object
#'
#' @param object Output from [describe_mfrm_data()].
#' @param digits Number of digits for numeric rounding.
#' @param top_n Maximum rows shown in preview blocks.
#' @param ... Reserved for generic compatibility.
#'
#' @details
#' This summary is intended as a compact pre-fit quality snapshot for
#' manuscripts and analysis logs.
#'
#' @section Interpreting output:
#' Recommended read order:
#' - `overview`: sample size, persons/facets/categories.
#' - `missing`: missingness hotspots by selected input columns.
#' - `score_distribution`: category usage balance.
#' - `facet_overview`: coverage per facet (minimum/maximum weighted counts).
#' - `agreement`: observed-score inter-rater agreement (when available).
#'
#' Very low `MinWeightedN` in `facet_overview` is a practical warning for
#' unstable downstream facet estimates.
#'
#' @section Typical workflow:
#' 1. Run [describe_mfrm_data()] on raw long-format data.
#' 2. Inspect `summary(ds)` before model fitting.
#' 3. Resolve sparse/missing issues, then run [fit_mfrm()].
#'
#' @return An object of class `summary.mfrm_data_description`.
#' @seealso [describe_mfrm_data()], [summary.mfrm_fit()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' ds <- describe_mfrm_data(toy, "Person", c("Rater", "Criterion"), "Score")
#' summary(ds)
#' @export
summary.mfrm_data_description <- function(object, digits = 3, top_n = 10, ...) {
  digits <- max(0L, as.integer(digits))
  top_n <- max(1L, as.integer(top_n))

  overview <- as.data.frame(object$overview %||% data.frame(), stringsAsFactors = FALSE)
  missing_tbl <- as.data.frame(object$missing_by_column %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(missing_tbl) > 0 && all(c("Column", "Missing") %in% names(missing_tbl))) {
    missing_tbl <- missing_tbl |>
      dplyr::arrange(dplyr::desc(.data$Missing), .data$Column) |>
      dplyr::slice_head(n = top_n)
  }

  score_dist <- as.data.frame(object$score_distribution %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(score_dist) > 0) {
    score_dist <- utils::head(score_dist, n = top_n)
  }

  facet_tbl <- as.data.frame(object$facet_level_summary %||% data.frame(), stringsAsFactors = FALSE)
  facet_overview <- data.frame()
  if (nrow(facet_tbl) > 0 && all(c("Facet", "Level", "WeightedN") %in% names(facet_tbl))) {
    facet_overview <- facet_tbl |>
      dplyr::group_by(.data$Facet) |>
      dplyr::summarise(
        Levels = dplyr::n_distinct(.data$Level),
        TotalWeightedN = sum(.data$WeightedN, na.rm = TRUE),
        MeanWeightedN = mean(.data$WeightedN, na.rm = TRUE),
        MinWeightedN = min(.data$WeightedN, na.rm = TRUE),
        MaxWeightedN = max(.data$WeightedN, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::arrange(.data$Facet) |>
      as.data.frame(stringsAsFactors = FALSE)
  }

  agreement_tbl <- as.data.frame(object$agreement$summary %||% data.frame(), stringsAsFactors = FALSE)
  notes <- if (nrow(missing_tbl) > 0 && any(suppressWarnings(as.numeric(missing_tbl$Missing)) > 0, na.rm = TRUE)) {
    "Missing values were detected in one or more input columns."
  } else {
    "No missing values were detected in selected input columns."
  }

  out <- list(
    overview = overview,
    missing = missing_tbl,
    score_distribution = score_dist,
    facet_overview = facet_overview,
    agreement = agreement_tbl,
    notes = notes,
    digits = digits
  )
  class(out) <- "summary.mfrm_data_description"
  out
}

#' @export
print.summary.mfrm_data_description <- function(x, ...) {
  digits <- as.integer(x$digits %||% 3L)
  if (!is.finite(digits)) digits <- 3L
  cat("mfrm Data Description Summary\n")
  if (!is.null(x$overview) && nrow(x$overview) > 0) {
    cat("\nOverview\n")
    print(round_numeric_df(as.data.frame(x$overview), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$missing) && nrow(x$missing) > 0) {
    cat("\nMissing by column\n")
    print(round_numeric_df(as.data.frame(x$missing), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$score_distribution) && nrow(x$score_distribution) > 0) {
    cat("\nScore distribution\n")
    print(round_numeric_df(as.data.frame(x$score_distribution), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$facet_overview) && nrow(x$facet_overview) > 0) {
    cat("\nFacet coverage\n")
    print(round_numeric_df(as.data.frame(x$facet_overview), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$agreement) && nrow(x$agreement) > 0) {
    cat("\nInter-rater agreement\n")
    print(round_numeric_df(as.data.frame(x$agreement), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$notes) && nzchar(x$notes)) {
    cat("\nNotes\n")
    cat(" - ", x$notes, "\n", sep = "")
  }
  invisible(x)
}

#' Plot a data-description object
#'
#' @param x Output from [describe_mfrm_data()].
#' @param y Reserved for generic compatibility.
#' @param type Plot type: `"score_distribution"`, `"facet_levels"`, or `"missing"`.
#' @param main Optional title override.
#' @param palette Optional named colors (`score`, `facet`, `missing`).
#' @param label_angle X-axis label angle for bar plots.
#' @param draw If `TRUE`, draw using base graphics.
#' @param ... Reserved for generic compatibility.
#'
#' @details
#' This method draws quick pre-fit quality views from [describe_mfrm_data()]:
#' - score distribution balance
#' - facet-level structure size
#' - missingness by selected columns
#'
#' @section Interpreting output:
#' - `"score_distribution"`: identifies sparse/unused categories.
#' - `"facet_levels"`: reveals highly imbalanced facet granularity.
#' - `"missing"`: pinpoints columns with potential data-quality bottlenecks.
#'
#' @section Typical workflow:
#' 1. Run [describe_mfrm_data()] before fitting.
#' 2. Inspect `summary(ds)` and `plot(ds, type = "missing")`.
#' 3. Check category/facet balance with other plot types.
#' 4. Fit model after resolving obvious data issues.
#'
#' @return A plotting-data object of class `mfrm_plot_data`.
#' @seealso [describe_mfrm_data()], `plot()`
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' ds <- describe_mfrm_data(toy, "Person", c("Rater", "Criterion"), "Score")
#' p <- plot(ds, draw = FALSE)
#' @export
plot.mfrm_data_description <- function(x,
                                       y = NULL,
                                       type = c("score_distribution", "facet_levels", "missing"),
                                       main = NULL,
                                       palette = NULL,
                                       label_angle = 45,
                                       draw = TRUE,
                                       ...) {
  type <- match.arg(tolower(as.character(type[1])), c("score_distribution", "facet_levels", "missing"))
  pal <- resolve_palette(
    palette = palette,
    defaults = c(
      score = "#2b8cbe",
      facet = "#31a354",
      missing = "#756bb1"
    )
  )

  if (type == "score_distribution") {
    tbl <- as.data.frame(x$score_distribution %||% data.frame(), stringsAsFactors = FALSE)
    if (nrow(tbl) == 0 || !all(c("Score", "WeightedN") %in% names(tbl))) {
      stop("Score distribution is not available. Ensure describe_mfrm_data() was run on valid data.", call. = FALSE)
    }
    if (isTRUE(draw)) {
      barplot_rot45(
        height = suppressWarnings(as.numeric(tbl$WeightedN)),
        labels = as.character(tbl$Score),
        col = pal["score"],
        main = if (is.null(main)) "Score distribution" else as.character(main[1]),
        ylab = "Weighted N",
        label_angle = label_angle,
        mar_bottom = 7.8
      )
    }
    return(invisible(new_mfrm_plot_data(
      "data_description",
      list(plot = "score_distribution", table = tbl)
    )))
  }

  if (type == "facet_levels") {
    tbl <- as.data.frame(x$facet_level_summary %||% data.frame(), stringsAsFactors = FALSE)
    if (nrow(tbl) == 0 || !all(c("Facet", "Level") %in% names(tbl))) {
      stop("Facet level summary is not available. Ensure describe_mfrm_data() was run on valid data.", call. = FALSE)
    }
    agg <- tbl |>
      dplyr::group_by(.data$Facet) |>
      dplyr::summarise(Levels = dplyr::n_distinct(.data$Level), .groups = "drop") |>
      dplyr::arrange(.data$Facet)
    if (isTRUE(draw)) {
      barplot_rot45(
        height = suppressWarnings(as.numeric(agg$Levels)),
        labels = as.character(agg$Facet),
        col = pal["facet"],
        main = if (is.null(main)) "Facet levels" else as.character(main[1]),
        ylab = "Levels",
        label_angle = label_angle,
        mar_bottom = 7.8
      )
    }
    return(invisible(new_mfrm_plot_data(
      "data_description",
      list(plot = "facet_levels", table = as.data.frame(agg, stringsAsFactors = FALSE))
    )))
  }

  tbl <- as.data.frame(x$missing_by_column %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(tbl) == 0 || !all(c("Column", "Missing") %in% names(tbl))) {
    stop("Missing-by-column table is not available. Ensure describe_mfrm_data() was run on valid data.", call. = FALSE)
  }
  if (isTRUE(draw)) {
    barplot_rot45(
      height = suppressWarnings(as.numeric(tbl$Missing)),
      labels = as.character(tbl$Column),
      col = pal["missing"],
      main = if (is.null(main)) "Missing values by column" else as.character(main[1]),
      ylab = "Missing",
      label_angle = label_angle,
      mar_bottom = 8.0
    )
  }
  invisible(new_mfrm_plot_data(
    "data_description",
    list(plot = "missing", table = tbl)
  ))
}

#' Audit and normalize anchor/group-anchor tables
#'
#' @param data A data.frame in long format (one row per rating event).
#' @param person Column name for person IDs.
#' @param facets Character vector of facet column names.
#' @param score Column name for observed score.
#' @param anchors Optional anchor table (Facet, Level, Anchor).
#' @param group_anchors Optional group-anchor table
#'   (Facet, Level, Group, GroupValue).
#' @param weight Optional weight/frequency column name.
#' @param rating_min Optional minimum category value.
#' @param rating_max Optional maximum category value.
#' @param keep_original Keep original category values.
#' @param min_common_anchors Minimum anchored levels per linking facet used in
#'   recommendations (default `5`).
#' @param min_obs_per_element Minimum weighted observations per facet level used
#'   in recommendations (default `30`).
#' @param min_obs_per_category Minimum weighted observations per score category
#'   used in recommendations (default `10`).
#' @param noncenter_facet One facet to leave non-centered.
#' @param dummy_facets Facets to fix at zero.
#'
#' @details
#' This function applies the same preprocessing and key-resolution rules as
#' `fit_mfrm()`, but returns an audit object so constraints can be checked
#' before estimation.
#'
#' FACETS-style behaviors used here:
#' - direct anchors fix level values
#' - grouped anchors constrain group means to `GroupValue`
#' - overlapping rows prefer direct anchors for that level
#' - missing `GroupValue` defaults to 0
#'
#' @section Interpreting output:
#' - `issue_counts`/`issues`: concrete data or specification problems.
#' - `facet_summary`: constraint coverage by facet.
#' - `design_checks`: whether anchor targets have enough observations.
#' - `recommendations`: action items before estimation.
#'
#' @section Typical workflow:
#' 1. Build candidate anchors (e.g., with [make_anchor_table()]).
#' 2. Run `audit_mfrm_anchors(...)`.
#' 3. Resolve issues, then fit with [fit_mfrm()].
#'
#' @return A list of class `mfrm_anchor_audit` with:
#' - `anchors`: cleaned anchor table used by estimation
#' - `group_anchors`: cleaned group-anchor table used by estimation
#' - `facet_summary`: counts of levels, constrained levels, and free levels
#' - `design_checks`: observation-count checks by level/category
#' - `thresholds`: active threshold settings used for recommendations
#' - `issue_counts`: issue-type counts
#' - `issues`: list of issue tables
#' - `recommendations`: FACETS-oriented guidance strings
#'
#' @seealso [fit_mfrm()], [describe_mfrm_data()], [make_anchor_table()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#'
#' anchors <- data.frame(
#'   Facet = c("Rater", "Rater"),
#'   Level = c("R1", "R1"),
#'   Anchor = c(0, 0.1),
#'   stringsAsFactors = FALSE
#' )
#' aud <- audit_mfrm_anchors(
#'   data = toy,
#'   person = "Person",
#'   facets = c("Rater", "Criterion"),
#'   score = "Score",
#'   anchors = anchors
#' )
#' aud$issue_counts
#' summary(aud)
#' p_aud <- plot(aud, draw = FALSE)
#' class(p_aud)
#' @export
audit_mfrm_anchors <- function(data,
                               person,
                               facets,
                               score,
                               anchors = NULL,
                               group_anchors = NULL,
                               weight = NULL,
                               rating_min = NULL,
                               rating_max = NULL,
                               keep_original = FALSE,
                               min_common_anchors = 5L,
                               min_obs_per_element = 30,
                               min_obs_per_category = 10,
                               noncenter_facet = "Person",
                               dummy_facets = NULL) {
  prep <- prepare_mfrm_data(
    data = data,
    person_col = person,
    facet_cols = facets,
    score_col = score,
    rating_min = rating_min,
    rating_max = rating_max,
    weight_col = weight,
    keep_original = keep_original
  )

  noncenter_facet <- sanitize_noncenter_facet(noncenter_facet, prep$facet_names)
  dummy_facets <- sanitize_dummy_facets(dummy_facets, prep$facet_names)

  audit <- audit_anchor_tables(
    prep = prep,
    anchor_df = anchors,
    group_anchor_df = group_anchors,
    min_common_anchors = min_common_anchors,
    min_obs_per_element = min_obs_per_element,
    min_obs_per_category = min_obs_per_category,
    noncenter_facet = noncenter_facet,
    dummy_facets = dummy_facets
  )
  class(audit) <- c("mfrm_anchor_audit", class(audit))
  audit
}

#' @export
print.mfrm_anchor_audit <- function(x, ...) {
  issue_total <- if (!is.null(x$issue_counts) && nrow(x$issue_counts) > 0) sum(x$issue_counts$N) else 0
  cat("mfrm anchor audit\n")
  cat("  issue rows: ", issue_total, "\n", sep = "")

  if (!is.null(x$issue_counts) && nrow(x$issue_counts) > 0) {
    nonzero <- x$issue_counts |>
      dplyr::filter(.data$N > 0)
    if (nrow(nonzero) > 0) {
      cat("\nIssue counts\n")
      print(nonzero, row.names = FALSE)
    }
  }

  if (!is.null(x$facet_summary) && nrow(x$facet_summary) > 0) {
    cat("\nFacet summary\n")
    print(x$facet_summary, row.names = FALSE)
  }

  if (!is.null(x$design_checks) &&
      !is.null(x$design_checks$level_observation_summary) &&
      nrow(x$design_checks$level_observation_summary) > 0) {
    cat("\nLevel observation summary\n")
    print(x$design_checks$level_observation_summary, row.names = FALSE)
  }

  invisible(x)
}

#' Summarize an anchor-audit object
#'
#' @param object Output from [audit_mfrm_anchors()].
#' @param digits Number of digits for numeric rounding.
#' @param top_n Maximum rows shown in issue previews.
#' @param ... Reserved for generic compatibility.
#'
#' @details
#' This summary provides a compact pre-estimation audit of anchor and
#' group-anchor specifications.
#'
#' @section Interpreting output:
#' Recommended order:
#' - `issue_counts`: primary triage table (non-zero issues first).
#' - `facet_summary`: anchored/grouped/free-level balance by facet.
#' - `level_observation_summary` and `category_counts`: sparse-cell diagnostics.
#' - `recommendations`: concrete remediation suggestions.
#'
#' If `issue_counts` is non-empty, treat anchor constraints as provisional and
#' resolve issues before final estimation.
#'
#' @section Typical workflow:
#' 1. Run [audit_mfrm_anchors()] with intended anchors/group anchors.
#' 2. Review `summary(aud)` and recommendations.
#' 3. Revise anchor tables, then call [fit_mfrm()].
#'
#' @return An object of class `summary.mfrm_anchor_audit`.
#' @seealso [audit_mfrm_anchors()], [fit_mfrm()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' aud <- audit_mfrm_anchors(toy, "Person", c("Rater", "Criterion"), "Score")
#' summary(aud)
#' @export
summary.mfrm_anchor_audit <- function(object, digits = 3, top_n = 10, ...) {
  digits <- max(0L, as.integer(digits))
  top_n <- max(1L, as.integer(top_n))

  issue_counts <- as.data.frame(object$issue_counts %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(issue_counts) > 0 && all(c("Issue", "N") %in% names(issue_counts))) {
    issue_counts <- issue_counts |>
      dplyr::filter(.data$N > 0) |>
      dplyr::arrange(dplyr::desc(.data$N), .data$Issue) |>
      dplyr::slice_head(n = top_n)
  }

  facet_summary <- as.data.frame(object$facet_summary %||% data.frame(), stringsAsFactors = FALSE)
  level_summary <- as.data.frame(object$design_checks$level_observation_summary %||% data.frame(), stringsAsFactors = FALSE)
  category_summary <- as.data.frame(object$design_checks$category_counts %||% data.frame(), stringsAsFactors = FALSE)

  recommendations <- as.character(object$recommendations %||% character(0))
  if (length(recommendations) > top_n) {
    recommendations <- recommendations[seq_len(top_n)]
  }

  notes <- if (nrow(issue_counts) > 0) {
    "Anchor-audit issues were detected. Review issue counts and recommendations."
  } else {
    "No anchor-table issue rows were detected."
  }

  out <- list(
    issue_counts = issue_counts,
    facet_summary = facet_summary,
    level_observation_summary = level_summary,
    category_counts = category_summary,
    recommendations = recommendations,
    notes = notes,
    digits = digits
  )
  class(out) <- "summary.mfrm_anchor_audit"
  out
}

#' @export
print.summary.mfrm_anchor_audit <- function(x, ...) {
  digits <- as.integer(x$digits %||% 3L)
  if (!is.finite(digits)) digits <- 3L

  cat("mfrm Anchor Audit Summary\n")
  if (!is.null(x$issue_counts) && nrow(x$issue_counts) > 0) {
    cat("\nIssue counts\n")
    print(round_numeric_df(as.data.frame(x$issue_counts), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$facet_summary) && nrow(x$facet_summary) > 0) {
    cat("\nFacet summary\n")
    print(round_numeric_df(as.data.frame(x$facet_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$level_observation_summary) && nrow(x$level_observation_summary) > 0) {
    cat("\nLevel observation summary\n")
    print(round_numeric_df(as.data.frame(x$level_observation_summary), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$category_counts) && nrow(x$category_counts) > 0) {
    cat("\nCategory counts\n")
    print(round_numeric_df(as.data.frame(x$category_counts), digits = digits), row.names = FALSE)
  }
  if (length(x$recommendations) > 0) {
    cat("\nRecommendations\n")
    for (line in x$recommendations) cat(" - ", line, "\n", sep = "")
  }
  if (!is.null(x$notes) && nzchar(x$notes)) {
    cat("\nNotes\n")
    cat(" - ", x$notes, "\n", sep = "")
  }
  invisible(x)
}

#' Plot an anchor-audit object
#'
#' @param x Output from [audit_mfrm_anchors()].
#' @param y Reserved for generic compatibility.
#' @param type Plot type: `"issue_counts"`, `"facet_constraints"`,
#'   or `"level_observations"`.
#' @param main Optional title override.
#' @param palette Optional named colors.
#' @param label_angle X-axis label angle for bar plots.
#' @param draw If `TRUE`, draw using base graphics.
#' @param ... Reserved for generic compatibility.
#'
#' @details
#' Base-R visualization helper for anchor audit outputs.
#'
#' @section Interpreting output:
#' - `"issue_counts"`: volume of each issue class.
#' - `"facet_constraints"`: anchored/grouped/free mix by facet.
#' - `"level_observations"`: observation support across levels.
#'
#' @section Typical workflow:
#' 1. Run [audit_mfrm_anchors()].
#' 2. Start with `plot(aud, type = "issue_counts")`.
#' 3. Inspect constraint and support plots before fitting.
#'
#' @return A plotting-data object of class `mfrm_plot_data`.
#' @seealso [audit_mfrm_anchors()], [make_anchor_table()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' aud <- audit_mfrm_anchors(toy, "Person", c("Rater", "Criterion"), "Score")
#' p <- plot(aud, draw = FALSE)
#' @export
plot.mfrm_anchor_audit <- function(x,
                                   y = NULL,
                                   type = c("issue_counts", "facet_constraints", "level_observations"),
                                   main = NULL,
                                   palette = NULL,
                                   label_angle = 45,
                                   draw = TRUE,
                                   ...) {
  type <- match.arg(tolower(as.character(type[1])), c("issue_counts", "facet_constraints", "level_observations"))
  pal <- resolve_palette(
    palette = palette,
    defaults = c(
      issues = "#cb181d",
      anchored = "#756bb1",
      grouped = "#9ecae1",
      levels = "#2b8cbe"
    )
  )

  if (type == "issue_counts") {
    tbl <- as.data.frame(x$issue_counts %||% data.frame(), stringsAsFactors = FALSE)
    if (nrow(tbl) == 0 || !all(c("Issue", "N") %in% names(tbl))) {
      stop("Issue-count table is not available. Ensure audit_mfrm_anchors() was run with valid anchor inputs.", call. = FALSE)
    }
    tbl <- tbl |>
      dplyr::arrange(dplyr::desc(.data$N), .data$Issue)
    if (isTRUE(draw)) {
      barplot_rot45(
        height = suppressWarnings(as.numeric(tbl$N)),
        labels = as.character(tbl$Issue),
        col = pal["issues"],
        main = if (is.null(main)) "Anchor-audit issue counts" else as.character(main[1]),
        ylab = "Rows",
        label_angle = label_angle,
        mar_bottom = 9.2
      )
    }
    return(invisible(new_mfrm_plot_data(
      "anchor_audit",
      list(plot = "issue_counts", table = tbl)
    )))
  }

  if (type == "facet_constraints") {
    tbl <- as.data.frame(x$facet_summary %||% data.frame(), stringsAsFactors = FALSE)
    if (nrow(tbl) == 0 || !all(c("Facet", "AnchoredLevels", "GroupedLevels", "FreeLevels") %in% names(tbl))) {
      stop("Facet summary with constraint columns is not available. Ensure audit_mfrm_anchors() was run with valid anchor inputs.", call. = FALSE)
    }
    if (isTRUE(draw)) {
      old_mar <- graphics::par("mar")
      on.exit(graphics::par(mar = old_mar), add = TRUE)
      mar <- old_mar
      mar[1] <- max(mar[1], 8.8)
      graphics::par(mar = mar)
      mat <- rbind(
        Anchored = suppressWarnings(as.numeric(tbl$AnchoredLevels)),
        Grouped = suppressWarnings(as.numeric(tbl$GroupedLevels)),
        Free = suppressWarnings(as.numeric(tbl$FreeLevels))
      )
      mids <- graphics::barplot(
        height = mat,
        beside = FALSE,
        names.arg = FALSE,
        col = c(pal["anchored"], pal["grouped"], "#d9d9d9"),
        border = "white",
        ylab = "Levels",
        main = if (is.null(main)) "Constraint profile by facet" else as.character(main[1])
      )
      draw_rotated_x_labels(
        at = mids,
        labels = as.character(tbl$Facet),
        srt = label_angle,
        cex = 0.82,
        line_offset = 0.085
      )
      graphics::legend(
        "topright",
        legend = c("Anchored", "Grouped", "Free"),
        fill = c(pal["anchored"], pal["grouped"], "#d9d9d9"),
        bty = "n",
        cex = 0.85
      )
    }
    return(invisible(new_mfrm_plot_data(
      "anchor_audit",
      list(plot = "facet_constraints", table = tbl)
    )))
  }

  tbl <- as.data.frame(x$design_checks$level_observation_summary %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(tbl) == 0 || !all(c("Facet", "MinObsPerLevel") %in% names(tbl))) {
    stop("Level observation summary is not available. Ensure audit_mfrm_anchors() was run with valid anchor inputs.", call. = FALSE)
  }
  if (isTRUE(draw)) {
    barplot_rot45(
      height = suppressWarnings(as.numeric(tbl$MinObsPerLevel)),
      labels = as.character(tbl$Facet),
      col = pal["levels"],
      main = if (is.null(main)) "Minimum observations per level" else as.character(main[1]),
      ylab = "Min observations",
      label_angle = label_angle,
      mar_bottom = 8.0
    )
    if ("RecommendedMinObs" %in% names(tbl)) {
      r <- suppressWarnings(as.numeric(tbl$RecommendedMinObs))
      r <- r[is.finite(r)]
      if (length(r) > 0) graphics::abline(h = unique(r)[1], lty = 2, col = "gray45")
    }
  }
  invisible(new_mfrm_plot_data(
    "anchor_audit",
    list(plot = "level_observations", table = tbl)
  ))
}

#' Build an anchor table from fitted estimates
#'
#' @param fit Output from [fit_mfrm()].
#' @param facets Optional subset of facets to include.
#' @param include_person Include person estimates as anchors.
#' @param digits Rounding digits for anchor values.
#'
#' @details
#' This helper supports FACETS-style linking workflows:
#' estimate one run, export stable reference levels, and reuse them as anchors
#' in subsequent calibrations.
#'
#' @section Interpreting output:
#' - `Facet`: facet name to be anchored in later runs.
#' - `Level`: specific element/level name inside that facet.
#' - `Anchor`: fixed logit value (rounded by `digits`).
#'
#' @section Typical workflow:
#' 1. Fit a reference run with [fit_mfrm()].
#' 2. Export anchors with `make_anchor_table(fit)`.
#' 3. Pass selected rows back into `fit_mfrm(..., anchors = ...)`.
#'
#' @return A data.frame with `Facet`, `Level`, and `Anchor`.
#' @seealso [fit_mfrm()], [audit_mfrm_anchors()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' anchors_tbl <- make_anchor_table(fit)
#' head(anchors_tbl)
#' summary(anchors_tbl$Anchor)
#' @export
make_anchor_table <- function(fit,
                              facets = NULL,
                              include_person = FALSE,
                              digits = 6) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().")
  }

  digits <- max(0L, as.integer(digits))
  out <- tibble::tibble()

  if (isTRUE(include_person) && !is.null(fit$facets$person) && nrow(fit$facets$person) > 0) {
    per <- tibble::as_tibble(fit$facets$person)
    if ("Person" %in% names(per) && "Estimate" %in% names(per)) {
      out <- dplyr::bind_rows(
        out,
        per |>
          dplyr::transmute(
            Facet = "Person",
            Level = as.character(.data$Person),
            Anchor = round(as.numeric(.data$Estimate), digits = digits)
          )
      )
    }
  }

  others <- tibble::as_tibble(fit$facets$others)
  if (nrow(others) > 0 && all(c("Facet", "Level", "Estimate") %in% names(others))) {
    out <- dplyr::bind_rows(
      out,
      others |>
        dplyr::transmute(
          Facet = as.character(.data$Facet),
          Level = as.character(.data$Level),
          Anchor = round(as.numeric(.data$Estimate), digits = digits)
        )
    )
  }

  if (!is.null(facets)) {
    keep <- as.character(facets)
    out <- out |>
      dplyr::filter(.data$Facet %in% keep)
  }

  out |>
    dplyr::arrange(.data$Facet, .data$Level)
}

#' Compute diagnostics for an `mfrm_fit` object
#'
#' @param fit Output from [fit_mfrm()].
#' @param interaction_pairs Optional list of facet pairs.
#' @param top_n_interactions Number of top interactions.
#' @param whexact Use exact ZSTD transformation.
#' @param residual_pca Residual PCA mode: `"none"`, `"overall"`, `"facet"`, or `"both"`.
#' @param pca_max_factors Maximum number of PCA factors to retain per matrix.
#'
#' @details
#' This function computes a diagnostic bundle used by downstream reporting.
#'
#' `interaction_pairs` controls which facet interactions are summarized.
#' Each element can be:
#' - a length-2 character vector such as `c("Rater", "Criterion")`, or
#' - omitted (`NULL`) to let the function select top interactions automatically.
#'
#' Residual PCA behavior:
#' - `"none"`: skip PCA
#' - `"overall"`: compute only overall residual PCA
#' - `"facet"`: compute only facet-specific residual PCA
#' - `"both"`: compute both sets
#'
#' @section Reading key components:
#' Practical interpretation often starts with:
#' - `overall_fit`: global infit/outfit and degrees of freedom.
#' - `reliability`: separation/reliability by facet (not a single pooled value).
#' - `fit`: element-level misfit scan (`Infit`, `Outfit`, `ZSTD`).
#' - `unexpected`, `fair_average`, `displacement`: targeted QC bundles.
#'
#' @section Interpreting output:
#' Start with `overall_fit` and `reliability`, then move to element-level
#' diagnostics (`fit`) and targeted bundles (`unexpected`, `fair_average`,
#' `displacement`, `interrater`, `facets_chisq`).
#'
#' Consistent signals across multiple components are typically more robust than
#' a single isolated warning.
#'
#' @section Typical workflow:
#' 1. Run `diagnose_mfrm(fit, residual_pca = "none")` for baseline diagnostics.
#' 2. Inspect `summary(diag)` and targeted tables/plots.
#' 3. If needed, rerun with residual PCA (`"overall"` or `"both"`).
#'
#' @return
#' An object of class `mfrm_diagnostics` including:
#' - `obs`: observed/expected/residual-level table
#' - `measures`: facet/person fit table (`Infit`, `Outfit`, `ZSTD`, `PTMEA`)
#' - `overall_fit`: overall fit summary
#' - `fit`: element-level fit diagnostics
#' - `reliability`: separation/reliability by facet
#' - `facets_chisq`: FACETS-style fixed/random facet chi-square summary
#' - `interactions`: top interaction diagnostics
#' - `interrater`: inter-rater agreement bundle (`summary`, `pairs`)
#' - `unexpected`: FACETS Table 4-style unexpected-response bundle
#' - `fair_average`: FACETS Table 12-style fair-average bundle
#' - `displacement`: displacement diagnostics bundle
#' - `residual_pca_overall`: optional overall PCA object
#' - `residual_pca_by_facet`: optional facet PCA objects
#'
#' @seealso [fit_mfrm()], [analyze_residual_pca()], [build_visual_summaries()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#' s_diag <- summary(diag)
#' s_diag$overview[, c("Observations", "Facets", "Categories")]
#' p_qc <- plot_qc_dashboard(fit, diagnostics = diag, draw = FALSE)
#' class(p_qc)
#'
#' # Optional: include residual PCA in the diagnostic bundle
#' diag_pca <- diagnose_mfrm(fit, residual_pca = "overall")
#' pca <- analyze_residual_pca(diag_pca, mode = "overall")
#' head(pca$overall_table)
#' @export
diagnose_mfrm <- function(fit,
                          interaction_pairs = NULL,
                          top_n_interactions = 20,
                          whexact = FALSE,
                          residual_pca = c("none", "overall", "facet", "both"),
                          pca_max_factors = 10L) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm(). ",
         "Got: ", paste(class(fit), collapse = "/"), ".", call. = FALSE)
  }
  residual_pca <- match.arg(tolower(residual_pca), c("none", "overall", "facet", "both"))

  out <- mfrm_diagnostics(
    fit,
    interaction_pairs = interaction_pairs,
    top_n_interactions = top_n_interactions,
    whexact = whexact,
    residual_pca = residual_pca,
    pca_max_factors = pca_max_factors
  )
  class(out) <- c("mfrm_diagnostics", class(out))
  out
}

#' Build a FACETS Table 10-style inter-rater agreement report
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()].
#' @param rater_facet Name of the rater facet. If `NULL`, inferred from facet names.
#' @param context_facets Optional context facets used to match observations for
#'   agreement. If `NULL`, all remaining facets (including `Person`) are used.
#' @param exact_warn Warning threshold for exact agreement.
#' @param corr_warn Warning threshold for pairwise correlation.
#' @param top_n Optional maximum number of pair rows to keep.
#'
#' @details
#' This helper computes pairwise rater agreement on matched contexts
#' and returns both a pair-level table and a one-row summary.
#'
#' @section Interpreting output:
#' - `summary`: overall agreement level, number/share of flagged pairs.
#' - `pairs`: pairwise exact agreement, correlation, and direction/size gaps.
#' - `settings`: applied facet matching and warning thresholds.
#'
#' Pairs flagged by both low exact agreement and low correlation generally
#' deserve highest calibration priority.
#'
#' @section Typical workflow:
#' 1. Run with explicit `rater_facet` (and `context_facets` if needed).
#' 2. Review `summary(ir)` and top flagged rows in `ir$pairs`.
#' 3. Visualize with [plot_interrater_agreement()].
#'
#' @section Output columns:
#' The `pairs` data.frame contains:
#' \describe{
#'   \item{Rater1, Rater2}{Rater pair identifiers.}
#'   \item{N}{Number of matched-context observations for this pair.}
#'   \item{Exact}{Proportion of exact score agreements.}
#'   \item{ExpectedExact}{Expected exact agreement under chance.}
#'   \item{Adjacent}{Proportion of adjacent (+/- 1 category) agreements.}
#'   \item{MeanDiff}{Signed mean score difference (Rater1 - Rater2).}
#'   \item{MAD}{Mean absolute score difference.}
#'   \item{Corr}{Pearson correlation between paired scores.}
#'   \item{Flag}{Logical; `TRUE` when Exact < `exact_warn` or Corr < `corr_warn`.}
#' }
#'
#' The `summary` data.frame contains:
#' \describe{
#'   \item{RaterFacet}{Name of the rater facet analyzed.}
#'   \item{TotalPairs}{Number of rater pairs evaluated.}
#'   \item{ExactAgreement}{Mean exact agreement across all pairs.}
#'   \item{MeanCorr}{Mean pairwise correlation.}
#'   \item{FlaggedPairs, FlaggedShare}{Count and proportion of flagged pairs.}
#' }
#'
#' @return A named list with:
#' - `summary`: one-row inter-rater summary
#' - `pairs`: pair-level agreement table
#' - `settings`: applied options and thresholds
#'
#' @seealso [diagnose_mfrm()], [facets_chisq_table()], [plot_interrater_agreement()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' ir <- interrater_agreement_table(fit, rater_facet = "Rater")
#' summary(ir)
#' p_ir <- plot(ir, draw = FALSE)
#' class(p_ir)
#' @export
interrater_agreement_table <- function(fit,
                                       diagnostics = NULL,
                                       rater_facet = NULL,
                                       context_facets = NULL,
                                       exact_warn = 0.50,
                                       corr_warn = 0.30,
                                       top_n = NULL) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().")
  }
  if (is.null(diagnostics)) {
    diagnostics <- diagnose_mfrm(fit, residual_pca = "none")
  }
  if (is.null(diagnostics$obs) || nrow(diagnostics$obs) == 0) {
    stop("`diagnostics$obs` is empty. Run diagnose_mfrm() first.")
  }

  known_facets <- c("Person", fit$config$facet_names)
  if (is.null(rater_facet) || !nzchar(as.character(rater_facet[1]))) {
    if (!is.null(diagnostics$interrater$summary) &&
        nrow(diagnostics$interrater$summary) > 0 &&
        "RaterFacet" %in% names(diagnostics$interrater$summary)) {
      rater_facet <- as.character(diagnostics$interrater$summary$RaterFacet[1])
    } else {
      rater_facet <- infer_default_rater_facet(fit$config$facet_names)
    }
  } else {
    rater_facet <- as.character(rater_facet[1])
  }
  if (is.null(rater_facet) || !rater_facet %in% known_facets) {
    stop("`rater_facet` must match one of: ", paste(known_facets, collapse = ", "))
  }
  if (identical(rater_facet, "Person")) {
    stop("`rater_facet = 'Person'` is not supported. Use a non-person facet.")
  }

  if (is.null(context_facets)) {
    facet_cols <- known_facets
  } else {
    context_facets <- unique(as.character(context_facets))
    unknown <- setdiff(context_facets, known_facets)
    if (length(unknown) > 0) {
      stop("Unknown `context_facets`: ", paste(unknown, collapse = ", "))
    }
    context_facets <- setdiff(context_facets, rater_facet)
    if (length(context_facets) == 0) {
      stop("`context_facets` must include at least one facet different from `rater_facet`.")
    }
    facet_cols <- c(rater_facet, context_facets)
  }

  agreement <- calc_interrater_agreement(
    obs_df = diagnostics$obs,
    facet_cols = facet_cols,
    rater_facet = rater_facet,
    res = fit
  )

  pairs <- as.data.frame(agreement$pairs, stringsAsFactors = FALSE)
  flagged_n <- 0L
  if (nrow(pairs) > 0) {
    pairs <- pairs |>
      mutate(
        Pair = paste(.data$Rater1, .data$Rater2, sep = " | "),
        ExactGap = ifelse(is.finite(.data$ExpectedExact), .data$Exact - .data$ExpectedExact, NA_real_),
        LowExactFlag = is.finite(.data$Exact) & .data$Exact < exact_warn,
        LowCorrFlag = is.finite(.data$Corr) & .data$Corr < corr_warn,
        Flag = .data$LowExactFlag | .data$LowCorrFlag
      ) |>
      arrange(desc(.data$Flag), .data$Exact, .data$Corr)
    flagged_n <- sum(pairs$Flag, na.rm = TRUE)
    if (!is.null(top_n)) {
      pairs <- pairs |>
        slice_head(n = max(1L, as.integer(top_n)))
    }
  }

  summary_tbl <- as.data.frame(agreement$summary, stringsAsFactors = FALSE)
  if (nrow(summary_tbl) > 0) {
    summary_tbl$FlaggedPairs <- flagged_n
    summary_tbl$FlaggedShare <- ifelse(summary_tbl$Pairs > 0, flagged_n / summary_tbl$Pairs, NA_real_)
  }

  out <- list(
    summary = summary_tbl,
    pairs = pairs,
    settings = list(
      rater_facet = rater_facet,
      context_facets = setdiff(facet_cols, rater_facet),
      exact_warn = exact_warn,
      corr_warn = corr_warn
    )
  )
  as_mfrm_bundle(out, "mfrm_interrater")
}

#' Build FACETS-style facet chi-square and random-effect diagnostics
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()].
#' @param fixed_p_max Warning cutoff for fixed-effect chi-square p-values.
#' @param random_p_max Warning cutoff for random-effect chi-square p-values.
#' @param top_n Optional maximum number of facet rows to keep.
#'
#' @details
#' This helper summarizes facet-level variability with fixed and random
#' chi-square indices, aligned to FACETS-style facet variance checks.
#'
#' @section Interpreting output:
#' - `table`: facet-level fixed/random chi-square and p-value flags.
#' - `summary`: number of significant facets and overall magnitude indicators.
#' - `thresholds`: p-value criteria used for flagging.
#'
#' Use this table together with inter-rater and displacement diagnostics to
#' distinguish global facet effects from local anomalies.
#'
#' @section Typical workflow:
#' 1. Run `facets_chisq_table(fit, ...)`.
#' 2. Inspect `summary(chi)` then facet rows in `chi$table`.
#' 3. Visualize with [plot_facets_chisq()].
#'
#' @section Output columns:
#' The `table` data.frame contains:
#' \describe{
#'   \item{Facet}{Facet name.}
#'   \item{Levels}{Number of estimated levels in this facet.}
#'   \item{MeanMeasure, SD}{Mean and standard deviation of level measures.}
#'   \item{FixedChiSq, FixedDF, FixedProb}{Fixed-effect chi-square test
#'     (null hypothesis: all levels equal). Significant result means the
#'     facet elements differ more than measurement error alone.}
#'   \item{RandomChiSq, RandomDF, RandomProb, RandomVar}{Random-effect test
#'     (null hypothesis: variation equals that of a random sample from a
#'     single population). Significant result suggests systematic
#'     heterogeneity beyond sampling variation.}
#'   \item{FixedFlag, RandomFlag}{Logical flags for significance.}
#' }
#'
#' @return A named list with:
#' - `table`: facet-level chi-square diagnostics
#' - `summary`: one-row summary
#' - `thresholds`: applied p-value thresholds
#'
#' @seealso [diagnose_mfrm()], [interrater_agreement_table()], [plot_facets_chisq()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' chi <- facets_chisq_table(fit)
#' summary(chi)
#' p_chi <- plot(chi, draw = FALSE)
#' class(p_chi)
#' @export
facets_chisq_table <- function(fit,
                               diagnostics = NULL,
                               fixed_p_max = 0.05,
                               random_p_max = 0.05,
                               top_n = NULL) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().")
  }
  if (is.null(diagnostics)) {
    diagnostics <- diagnose_mfrm(fit, residual_pca = "none")
  }
  if (is.null(diagnostics$measures) || nrow(diagnostics$measures) == 0) {
    stop("`diagnostics$measures` is empty. Run diagnose_mfrm() first.")
  }

  tbl <- calc_facets_chisq(diagnostics$measures)
  if (nrow(tbl) > 0) {
    tbl <- tbl |>
      mutate(
        FixedFlag = is.finite(.data$FixedProb) & .data$FixedProb < fixed_p_max,
        RandomFlag = is.finite(.data$RandomProb) & .data$RandomProb < random_p_max
      ) |>
      arrange(desc(.data$FixedChiSq))
    if (!is.null(top_n)) {
      tbl <- tbl |>
        slice_head(n = max(1L, as.integer(top_n)))
    }
  }

  summary_tbl <- if (nrow(tbl) == 0) {
    data.frame()
  } else {
    safe_max <- function(x) {
      x <- x[is.finite(x)]
      if (length(x) == 0) NA_real_ else max(x)
    }
    safe_mean <- function(x) {
      x <- x[is.finite(x)]
      if (length(x) == 0) NA_real_ else mean(x)
    }
    data.frame(
      Facets = nrow(tbl),
      FixedSignificant = sum(tbl$FixedFlag, na.rm = TRUE),
      RandomSignificant = sum(tbl$RandomFlag, na.rm = TRUE),
      MeanRandomVar = safe_mean(tbl$RandomVar),
      MaxFixedChiSq = safe_max(tbl$FixedChiSq),
      MaxRandomChiSq = safe_max(tbl$RandomChiSq),
      stringsAsFactors = FALSE
    )
  }

  out <- list(
    table = as.data.frame(tbl, stringsAsFactors = FALSE),
    summary = summary_tbl,
    thresholds = list(
      fixed_p_max = fixed_p_max,
      random_p_max = random_p_max
    )
  )
  as_mfrm_bundle(out, "mfrm_facets_chisq")
}

#' Build a FACETS Table 4-style unexpected-response report
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()].
#' @param abs_z_min Absolute standardized-residual cutoff.
#' @param prob_max Maximum observed-category probability cutoff.
#' @param top_n Maximum number of rows to return.
#' @param rule Flagging rule: `"either"` (default) or `"both"`.
#'
#' @details
#' A response is flagged as unexpected when:
#' - `rule = "either"`: `|StdResidual| >= abs_z_min` OR `ObsProb <= prob_max`
#' - `rule = "both"`: both conditions must be met.
#'
#' The table includes row-level observed/expected values, residuals,
#' observed-category probability, most-likely category, and a composite
#' severity score for sorting.
#'
#' @section Interpreting output:
#' - `summary`: prevalence of unexpected responses under current thresholds.
#' - `table`: ranked row-level diagnostics for case review.
#' - `thresholds`: active cutoffs and flagging rule.
#'
#' Compare results across `rule = "either"` and `rule = "both"` to assess how
#' conservative your screening should be.
#'
#' @section Typical workflow:
#' 1. Start with `rule = "either"` for broad screening.
#' 2. Re-run with `rule = "both"` for strict subset.
#' 3. Inspect top rows and visualize with [plot_unexpected()].
#'
#' @section Output columns:
#' The `table` data.frame contains:
#' \describe{
#'   \item{Row}{Original row index in the prepared data.}
#'   \item{Person}{Person identifier (plus one column per facet).}
#'   \item{Score}{Observed score category.}
#'   \item{Observed, Expected}{Observed and model-expected score values.}
#'   \item{Residual, StdResidual}{Raw and standardized residuals.}
#'   \item{ObsProb}{Probability of the observed category under the model.}
#'   \item{MostLikely, MostLikelyProb}{Most probable category and its
#'     probability.}
#'   \item{Severity}{Composite severity index (higher = more unexpected).}
#'   \item{Direction}{"Higher than expected" or "Lower than expected".}
#'   \item{FlagLowProbability, FlagLargeResidual}{Logical flags for each
#'     criterion.}
#' }
#'
#' The `summary` data.frame contains:
#' \describe{
#'   \item{TotalObservations}{Total observations analyzed.}
#'   \item{UnexpectedN, UnexpectedPercent}{Count and share of flagged rows.}
#'   \item{AbsZThreshold, ProbThreshold}{Applied cutoff values.}
#'   \item{Rule}{"either" or "both".}
#' }
#'
#' @return A named list with:
#' - `table`: flagged response rows
#' - `summary`: one-row overview
#' - `thresholds`: applied thresholds
#'
#' @seealso [diagnose_mfrm()], [displacement_table()], [fair_average_table()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' t4 <- unexpected_response_table(fit, abs_z_min = 1.5, prob_max = 0.4, top_n = 10)
#' summary(t4)
#' p_t4 <- plot(t4, draw = FALSE)
#' class(p_t4)
#' @export
unexpected_response_table <- function(fit,
                                      diagnostics = NULL,
                                      abs_z_min = 2,
                                      prob_max = 0.30,
                                      top_n = 100,
                                      rule = c("either", "both")) {
  rule <- match.arg(tolower(rule), c("either", "both"))
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().")
  }
  if (is.null(diagnostics)) {
    diagnostics <- diagnose_mfrm(fit, residual_pca = "none")
  }
  if (is.null(diagnostics$obs) || nrow(diagnostics$obs) == 0) {
    stop("`diagnostics$obs` is empty. Run diagnose_mfrm() first.")
  }

  tbl <- calc_unexpected_response_table(
    obs_df = diagnostics$obs,
    probs = compute_prob_matrix(fit),
    facet_names = fit$config$facet_names,
    rating_min = fit$prep$rating_min,
    abs_z_min = abs_z_min,
    prob_max = prob_max,
    top_n = top_n,
    rule = rule
  )
  summary_tbl <- summarize_unexpected_response_table(
    unexpected_tbl = tbl,
    total_observations = nrow(diagnostics$obs),
    abs_z_min = abs_z_min,
    prob_max = prob_max,
    rule = rule
  )

  out <- list(
    table = tbl,
    summary = summary_tbl,
    thresholds = list(
      abs_z_min = abs_z_min,
      prob_max = prob_max,
      rule = rule
    )
  )
  as_mfrm_bundle(out, "mfrm_unexpected")
}

#' Build a FACETS Table 12-style fair-average table bundle
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()].
#' @param facets Optional subset of facets.
#' @param totalscore Include all observations for score totals (`TRUE`) or apply
#'   FACETS-style extreme-row exclusion (`FALSE`).
#' @param umean Additive score-to-report origin shift.
#' @param uscale Multiplicative score-to-report scale.
#' @param udecimals Rounding digits used in formatted output.
#' @param omit_unobserved If `TRUE`, remove unobserved levels.
#' @param xtreme Extreme-score adjustment amount.
#'
#' @details
#' This function wraps internal FACETS-style fair-average calculations and
#' returns both facet-wise and stacked tables, including `Fair(M) Average`
#' and `Fair(Z) Average`.
#'
#' @section Interpreting output:
#' - `stacked`: cross-facet table for global comparison.
#' - `by_facet`: per-facet formatted tables for reporting.
#' - `raw_by_facet`: unformatted values for custom analyses/plots.
#' - `settings`: scoring-transformation and filtering options used.
#'
#' Larger observed-vs-fair gaps can indicate systematic scoring tendencies by
#' specific facet levels.
#'
#' @section Typical workflow:
#' 1. Run `fair_average_table(fit, ...)`.
#' 2. Inspect `summary(t12)` and `t12$stacked`.
#' 3. Visualize with [plot_fair_average()].
#'
#' @section Output columns:
#' The `stacked` data.frame contains:
#' \describe{
#'   \item{Facet}{Facet name for this row.}
#'   \item{Level}{Element label within the facet.}
#'   \item{Obsvd Average}{Observed raw-score average.}
#'   \item{Fair(M) Average}{Model-adjusted fair average (mean-based).}
#'   \item{Fair(Z) Average}{Model-adjusted fair average (z-score-based).}
#'   \item{Measure}{Estimated logit measure for this level.}
#'   \item{SE}{Standard error of the measure.}
#'   \item{Infit MnSq, Outfit MnSq}{Fit statistics for this level.}
#' }
#'
#' @return A named list with:
#' - `by_facet`: named list of formatted data.frames
#' - `stacked`: one stacked data.frame across facets
#' - `raw_by_facet`: unformatted internal tables
#' - `settings`: resolved options
#'
#' @seealso [diagnose_mfrm()], [unexpected_response_table()], [displacement_table()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' t12 <- fair_average_table(fit, udecimals = 2)
#' summary(t12)
#' p_t12 <- plot(t12, draw = FALSE)
#' class(p_t12)
#' @export
fair_average_table <- function(fit,
                               diagnostics = NULL,
                               facets = NULL,
                               totalscore = TRUE,
                               umean = 0,
                               uscale = 1,
                               udecimals = 2,
                               omit_unobserved = FALSE,
                               xtreme = 0) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().")
  }
  if (is.null(diagnostics)) {
    diagnostics <- diagnose_mfrm(fit, residual_pca = "none")
  }
  if (is.null(diagnostics$obs) || is.null(diagnostics$measures)) {
    stop("`diagnostics` must include both `obs` and `measures`.")
  }

  bundle <- calc_fair_average_bundle(
    res = fit,
    diagnostics = diagnostics,
    facets = facets,
    totalscore = totalscore,
    umean = umean,
    uscale = uscale,
    udecimals = udecimals,
    omit_unobserved = omit_unobserved,
    xtreme = xtreme
  )
  bundle$settings <- list(
    facets = if (is.null(facets)) NULL else as.character(facets),
    totalscore = totalscore,
    umean = umean,
    uscale = uscale,
    udecimals = udecimals,
    omit_unobserved = omit_unobserved,
    xtreme = xtreme
  )
  as_mfrm_bundle(bundle, "mfrm_fair_average")
}

#' Compute displacement diagnostics for facet levels
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()].
#' @param facets Optional subset of facets.
#' @param anchored_only If `TRUE`, keep only directly/group anchored levels.
#' @param abs_displacement_warn Absolute displacement warning threshold.
#' @param abs_t_warn Absolute displacement t-value warning threshold.
#' @param top_n Optional maximum number of rows to keep after sorting.
#'
#' @details
#' Displacement is computed as a one-step Newton update:
#' `sum(residual) / sum(information)` for each facet level.
#' This approximates how much a level would move if constraints were relaxed.
#'
#' @section Interpreting output:
#' - `table`: level-wise displacement and flag indicators.
#' - `summary`: count/share of flagged levels.
#' - `thresholds`: displacement and t-value cutoffs.
#'
#' Large absolute displacement in anchored levels suggests potential instability
#' in anchor assumptions.
#'
#' @section Typical workflow:
#' 1. Run `displacement_table(fit, anchored_only = TRUE)` for anchor checks.
#' 2. Inspect `summary(disp)` then detailed rows.
#' 3. Visualize with [plot_displacement()].
#'
#' @section Output columns:
#' The `table` data.frame contains:
#' \describe{
#'   \item{Facet, Level}{Facet name and element label.}
#'   \item{Displacement}{One-step Newton displacement estimate (logits).}
#'   \item{DisplacementSE}{Standard error of the displacement.}
#'   \item{DisplacementT}{Displacement / SE ratio.}
#'   \item{Estimate, SE}{Current measure estimate and its standard error.}
#'   \item{N}{Number of observations involving this level.}
#'   \item{AnchorValue, AnchorStatus, AnchorType}{Anchor metadata.}
#'   \item{Flag}{Logical; `TRUE` when displacement exceeds thresholds.}
#' }
#'
#' @return A named list with:
#' - `table`: displacement diagnostics by level
#' - `summary`: one-row summary
#' - `thresholds`: applied thresholds
#'
#' @seealso [diagnose_mfrm()], [unexpected_response_table()], [fair_average_table()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' disp <- displacement_table(fit, anchored_only = FALSE)
#' summary(disp)
#' p_disp <- plot(disp, draw = FALSE)
#' class(p_disp)
#' @export
displacement_table <- function(fit,
                               diagnostics = NULL,
                               facets = NULL,
                               anchored_only = FALSE,
                               abs_displacement_warn = 0.5,
                               abs_t_warn = 2,
                               top_n = NULL) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().")
  }
  if (is.null(diagnostics)) {
    diagnostics <- diagnose_mfrm(fit, residual_pca = "none")
  }
  if (is.null(diagnostics$obs) || nrow(diagnostics$obs) == 0) {
    stop("`diagnostics$obs` is empty. Run diagnose_mfrm() first.")
  }

  tbl <- calc_displacement_table(
    obs_df = diagnostics$obs,
    res = fit,
    measures = diagnostics$measures,
    abs_displacement_warn = abs_displacement_warn,
    abs_t_warn = abs_t_warn
  )
  if (!is.null(facets) && nrow(tbl) > 0) {
    tbl <- tbl |>
      filter(.data$Facet %in% as.character(facets))
  }
  if (isTRUE(anchored_only) && nrow(tbl) > 0) {
    tbl <- tbl |>
      filter(.data$AnchorType %in% c("Anchor", "Group"))
  }
  if (!is.null(top_n) && nrow(tbl) > 0) {
    top_n <- max(1L, as.integer(top_n))
    tbl <- tbl |>
      slice_head(n = top_n)
  }

  summary_tbl <- summarize_displacement_table(
    displacement_tbl = tbl,
    abs_displacement_warn = abs_displacement_warn,
    abs_t_warn = abs_t_warn
  )

  out <- list(
    table = tbl,
    summary = summary_tbl,
    thresholds = list(
      abs_displacement_warn = abs_displacement_warn,
      abs_t_warn = abs_t_warn
    )
  )
  as_mfrm_bundle(out, "mfrm_displacement")
}

#' Build a FACETS Table 5-style measurable data summary
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()].
#'
#' @details
#' This helper consolidates measurable-data diagnostics into a dedicated
#' report bundle: run-level summary, facet coverage, category usage, and
#' subset (connected-component) information.
#'
#' `summary(t5)` is supported through `summary()`.
#' `plot(t5)` is dispatched through `plot()` for class
#' `mfrm_measurable` (`type = "facet_coverage"`, `"category_counts"`,
#' `"subset_observations"`).
#'
#' @section Interpreting output:
#' - `summary`: overall measurable design status.
#' - `facet_coverage`: spread/precision by facet.
#' - `category_stats`: category usage and fit context.
#' - `subsets`: connectivity diagnostics (fragmented subsets reduce comparability).
#'
#' @section Typical workflow:
#' 1. Run `measurable_summary_table(fit)`.
#' 2. Check `summary(t5)` for subset/connectivity warnings.
#' 3. Use `plot(t5, ...)` to inspect facet/category/subset views.
#'
#' @section Output columns:
#' The `summary` data.frame (one row) contains:
#' \describe{
#'   \item{Observations, TotalWeight}{Total observations and summed weight.}
#'   \item{Persons, Facets, Categories}{Design dimensions.}
#'   \item{ConnectedSubsets}{Number of connected subsets.}
#'   \item{LargestSubsetObs, LargestSubsetPct}{Largest subset coverage.}
#' }
#'
#' The `facet_coverage` data.frame contains:
#' \describe{
#'   \item{Facet}{Facet name.}
#'   \item{Levels}{Number of estimated levels.}
#'   \item{MeanSE}{Mean standard error across levels.}
#'   \item{MeanInfit, MeanOutfit}{Mean fit statistics across levels.}
#'   \item{MinEstimate, MaxEstimate}{Measure range for this facet.}
#' }
#'
#' The `category_stats` data.frame contains:
#' \describe{
#'   \item{Category}{Score category value.}
#'   \item{Count, Percent}{Observed count and percentage.}
#'   \item{Infit, Outfit, InfitZSTD, OutfitZSTD}{Category-level fit.}
#'   \item{ExpectedCount, DiffCount, LowCount}{Expected-observed comparison
#'     and low-count flag.}
#' }
#'
#' @return A named list with:
#' - `summary`: one-row measurable-data summary
#' - `facet_coverage`: per-facet coverage summary
#' - `category_stats`: category-level usage/fit summary
#' - `subsets`: subset summary table (when available)
#'
#' @seealso [diagnose_mfrm()], [rating_scale_table()], [describe_mfrm_data()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' t5 <- measurable_summary_table(fit)
#' summary(t5)
#' p_t5 <- plot(t5, draw = FALSE)
#' class(p_t5)
#' @export
measurable_summary_table <- function(fit, diagnostics = NULL) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().")
  }
  if (is.null(diagnostics)) {
    diagnostics <- diagnose_mfrm(fit, residual_pca = "none")
  }
  if (is.null(diagnostics$obs) || nrow(diagnostics$obs) == 0) {
    stop("`diagnostics$obs` is empty. Run diagnose_mfrm() first.")
  }

  obs_df <- diagnostics$obs
  total_weight <- if ("Weight" %in% names(obs_df)) sum(obs_df$Weight, na.rm = TRUE) else nrow(obs_df)
  subset_tbl <- if (!is.null(diagnostics$subsets$summary)) as.data.frame(diagnostics$subsets$summary) else data.frame()
  subset_n <- if (nrow(subset_tbl) > 0 && "Subset" %in% names(subset_tbl)) nrow(subset_tbl) else NA_integer_
  largest_subset_obs <- if (nrow(subset_tbl) > 0 && "Observations" %in% names(subset_tbl)) max(subset_tbl$Observations, na.rm = TRUE) else NA_real_

  summary_tbl <- data.frame(
    Observations = nrow(obs_df),
    TotalWeight = total_weight,
    Persons = length(fit$prep$levels$Person),
    Facets = length(fit$config$facet_names),
    Categories = fit$prep$rating_max - fit$prep$rating_min + 1,
    ConnectedSubsets = subset_n,
    LargestSubsetObs = largest_subset_obs,
    LargestSubsetPct = ifelse(is.finite(largest_subset_obs) && nrow(obs_df) > 0, 100 * largest_subset_obs / nrow(obs_df), NA_real_),
    stringsAsFactors = FALSE
  )

  facet_coverage <- if (!is.null(diagnostics$measures) && nrow(diagnostics$measures) > 0) {
    as.data.frame(diagnostics$measures) |>
      group_by(.data$Facet) |>
      summarize(
        Levels = n(),
        MeanSE = mean(.data$SE, na.rm = TRUE),
        MeanInfit = mean(.data$Infit, na.rm = TRUE),
        MeanOutfit = mean(.data$Outfit, na.rm = TRUE),
        MinEstimate = min(.data$Estimate, na.rm = TRUE),
        MaxEstimate = max(.data$Estimate, na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    data.frame()
  }

  category_stats <- as.data.frame(calc_category_stats(obs_df, res = fit, whexact = FALSE), stringsAsFactors = FALSE)

  out <- list(
    summary = summary_tbl,
    facet_coverage = facet_coverage,
    category_stats = category_stats,
    subsets = subset_tbl
  )
  as_mfrm_bundle(out, "mfrm_measurable")
}

#' Build a FACETS Table 8.1-style rating-scale report
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()].
#' @param whexact Use exact ZSTD transformation for category fit.
#' @param drop_unused If `TRUE`, remove categories with zero count.
#'
#' @details
#' This helper provides category usage/fit statistics and threshold summaries
#' in a format aligned with FACETS Table 8.1 style checks.
#'
#' Typical checks:
#' - sparse category usage (`Count`, `ExpectedCount`)
#' - category fit (`Infit`, `Outfit`, `ZStd`)
#' - threshold ordering (`threshold_table$Estimate`, `GapFromPrev`)
#'
#' @section Interpreting output:
#' Start with `summary`:
#' - `UsedCategories` close to total `Categories` suggests stable category usage.
#' - very small `MinCategoryCount` indicates potential instability.
#' - `ThresholdMonotonic = FALSE` indicates disordered thresholds.
#'
#' Then inspect:
#' - `category_table` for category-level misfit/sparsity.
#' - `threshold_table` for adjacent-step gaps and ordering.
#'
#' @section Typical workflow:
#' 1. Fit model: [fit_mfrm()].
#' 2. Build diagnostics: [diagnose_mfrm()].
#' 3. Run `rating_scale_table()` and review `summary()`.
#' 4. Use `plot()` to visualize category profile quickly.
#'
#' @section Output columns:
#' The `category_table` data.frame contains:
#' \describe{
#'   \item{Category}{Score category value.}
#'   \item{Count, Percent}{Observed count and percentage of total.}
#'   \item{AvgPersonMeasure}{Mean person measure for respondents in this
#'     category.}
#'   \item{Infit, Outfit}{Category-level fit statistics.}
#'   \item{InfitZSTD, OutfitZSTD}{Standardized fit values.}
#'   \item{ExpectedCount, DiffCount}{Expected count and observed-expected
#'     difference.}
#'   \item{LowCount}{Logical; `TRUE` if count is below minimum threshold.}
#'   \item{InfitFlag, OutfitFlag, ZSTDFlag}{Fit-based warning flags.}
#' }
#'
#' The `threshold_table` data.frame contains:
#' \describe{
#'   \item{Step}{Step label (e.g., "1-2", "2-3").}
#'   \item{Estimate}{Estimated threshold/step difficulty (logits).}
#'   \item{GapFromPrev}{Difference from the previous threshold.  Gaps below
#'     1.4 logits may indicate category underuse; gaps above 5.0 may
#'     indicate wide unused regions (Linacre, 2002).}
#' }
#'
#' @return A named list with:
#' - `category_table`: category-level counts, expected counts, fit, and ZSTD
#' - `threshold_table`: model step/threshold estimates
#' - `summary`: one-row summary (usage and threshold monotonicity)
#'
#' @seealso [diagnose_mfrm()], [measurable_summary_table()], [plot.mfrm_fit()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' t8 <- rating_scale_table(fit)
#' summary(t8)
#' summary(t8)$summary
#' p_t8 <- plot(t8, draw = FALSE)
#' class(p_t8)
#' @export
rating_scale_table <- function(fit,
                               diagnostics = NULL,
                               whexact = FALSE,
                               drop_unused = FALSE) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().")
  }
  if (is.null(diagnostics)) {
    diagnostics <- diagnose_mfrm(fit, residual_pca = "none")
  }
  if (is.null(diagnostics$obs) || nrow(diagnostics$obs) == 0) {
    stop("`diagnostics$obs` is empty. Run diagnose_mfrm() first.")
  }

  cat_tbl <- as.data.frame(calc_category_stats(diagnostics$obs, res = fit, whexact = whexact), stringsAsFactors = FALSE)
  if (isTRUE(drop_unused) && nrow(cat_tbl) > 0 && "Count" %in% names(cat_tbl)) {
    cat_tbl <- cat_tbl[cat_tbl$Count > 0, , drop = FALSE]
  }

  step_tbl <- as.data.frame(fit$steps, stringsAsFactors = FALSE)
  if (nrow(step_tbl) > 0 && all(c("Step", "Estimate") %in% names(step_tbl))) {
    ord <- order(step_index_from_label(step_tbl$Step))
    step_tbl <- step_tbl[ord, , drop = FALSE]
    step_tbl$GapFromPrev <- c(NA_real_, diff(suppressWarnings(as.numeric(step_tbl$Estimate))))
  }

  threshold_monotonic <- if (nrow(step_tbl) > 1 && "Estimate" %in% names(step_tbl)) {
    est <- suppressWarnings(as.numeric(step_tbl$Estimate))
    all(diff(est) >= -sqrt(.Machine$double.eps), na.rm = TRUE)
  } else {
    NA
  }

  summary_tbl <- data.frame(
    Categories = nrow(cat_tbl),
    UsedCategories = if ("Count" %in% names(cat_tbl)) sum(cat_tbl$Count > 0, na.rm = TRUE) else NA_integer_,
    MinCategoryCount = if ("Count" %in% names(cat_tbl) && nrow(cat_tbl) > 0) min(cat_tbl$Count, na.rm = TRUE) else NA_real_,
    MaxCategoryCount = if ("Count" %in% names(cat_tbl) && nrow(cat_tbl) > 0) max(cat_tbl$Count, na.rm = TRUE) else NA_real_,
    MeanCategoryInfit = if ("Infit" %in% names(cat_tbl)) mean(cat_tbl$Infit, na.rm = TRUE) else NA_real_,
    MeanCategoryOutfit = if ("Outfit" %in% names(cat_tbl)) mean(cat_tbl$Outfit, na.rm = TRUE) else NA_real_,
    ThresholdMonotonic = threshold_monotonic,
    stringsAsFactors = FALSE
  )

  out <- list(
    category_table = cat_tbl,
    threshold_table = step_tbl,
    summary = summary_tbl
  )
  as_mfrm_bundle(out, "mfrm_rating_scale")
}

#' Build a FACETS Table 11-style bias-count report
#'
#' @param bias_results Output from [estimate_bias()].
#' @param min_count_warn Minimum count threshold for flagging sparse bias cells.
#' @param branch Output branch:
#'   `"facets"` keeps FACETS-style naming, `"original"` returns compact QC-oriented names.
#' @param fit Optional [fit_mfrm()] result used to attach run context metadata.
#'
#' @details
#' FACETS Table 11 emphasizes how many observations contribute to each
#' bias-cell estimate. This helper extracts those counts and summarizes
#' sparsity patterns.
#'
#' Branch behavior:
#' - `"facets"`: keeps FACETS-like column labels (`Sq`, `Observd Count`,
#'   `Obs-Exp Average`, `Model S.E.`) for side-by-side interpretation.
#' - `"original"`: keeps compact field names (`Count`, `BiasSize`, `SE`) for
#'   custom QC workflows and scripting.
#'
#' @section Interpreting output:
#' - `table`: cell-level contribution counts and low-count flags.
#' - `by_facet`: sparse-cell structure by each interaction facet.
#' - `summary`: overall low-count prevalence.
#' - `fit_overview`: optional run context (when `fit` is supplied).
#'
#' Low-count cells should be interpreted cautiously because bias-size estimates
#' can become unstable with sparse support.
#'
#' @section Typical workflow:
#' 1. Estimate bias with [estimate_bias()].
#' 2. Build `bias_count_table(...)` in desired branch.
#' 3. Review low-count flags before interpreting bias magnitudes.
#'
#' @section Output columns:
#' The `table` data.frame contains (FACETS-style branch):
#' \describe{
#'   \item{<FacetA>, <FacetB>}{Interaction facet level identifiers.}
#'   \item{Sq}{Sequential row number.}
#'   \item{Observd Count}{Number of observations for this cell.}
#'   \item{Obs-Exp Average}{Observed minus expected average for this cell.}
#'   \item{Model S.E.}{Standard error of the bias estimate.}
#'   \item{Infit, Outfit}{Fit statistics for this cell.}
#'   \item{LowCountFlag}{Logical; `TRUE` when count < `min_count_warn`.}
#' }
#'
#' The `summary` data.frame contains:
#' \describe{
#'   \item{InteractionFacets}{Names of the interaction facets.}
#'   \item{Cells, TotalCount}{Number of cells and total observations.}
#'   \item{LowCountCells, LowCountPercent}{Number and share of low-count
#'     cells.}
#' }
#'
#' @return A named list with:
#' - `table`: cell-level counts with low-count flags
#' - `by_facet`: named list of counts aggregated by each interaction facet
#' - `by_facet_a`, `by_facet_b`: first two facet summaries (legacy compatibility)
#' - `summary`: one-row summary
#' - `thresholds`: applied thresholds
#' - `branch`, `style`: output branch metadata
#' - `fit_overview`: optional one-row fit metadata when `fit` is supplied
#'
#' @seealso [estimate_bias()], [unexpected_after_bias_table()], [build_fixed_reports()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#' bias <- estimate_bias(fit, diag, facet_a = "Rater", facet_b = "Criterion", max_iter = 2)
#' t11 <- bias_count_table(bias)
#' t11_facets <- bias_count_table(bias, branch = "facets", fit = fit)
#' summary(t11)
#' p <- plot(t11, draw = FALSE)
#' p2 <- plot(t11, type = "lowcount_by_facet", draw = FALSE)
#' if (interactive()) {
#'   plot(
#'     t11,
#'     type = "cell_counts",
#'     draw = TRUE,
#'     main = "Bias Cell Counts (Customized)",
#'     palette = c(count = "#2b8cbe", low = "#cb181d"),
#'     label_angle = 45
#'   )
#' }
#' @export
bias_count_table <- function(bias_results,
                             min_count_warn = 10,
                             branch = c("original", "facets"),
                             fit = NULL) {
  branch <- match.arg(tolower(as.character(branch[1])), c("original", "facets"))

  if (is.null(bias_results) || is.null(bias_results$table) || nrow(bias_results$table) == 0) {
    stop("`bias_results` must be output from estimate_bias() with non-empty `table`.")
  }
  spec <- extract_bias_facet_spec(bias_results)
  tbl <- as.data.frame(bias_results$table, stringsAsFactors = FALSE)
  if (is.null(spec) || length(spec$facets) < 2) {
    stop("`bias_results$table` does not include recognizable interaction facet columns.")
  }
  req_cols <- c(spec$level_cols, "Observd Count", "Bias Size", "S.E.")
  if (!all(req_cols %in% names(tbl))) {
    stop("`bias_results$table` does not include required columns.")
  }

  min_count_warn <- max(0, as.numeric(min_count_warn))
  level_tbl <- tbl[, spec$level_cols, drop = FALSE]
  level_tbl[] <- lapply(level_tbl, as.character)
  names(level_tbl) <- spec$facets

  cell_tbl <- dplyr::bind_cols(
    level_tbl,
    tbl |>
      dplyr::transmute(
        Count = suppressWarnings(as.numeric(.data$`Observd Count`)),
        BiasSize = suppressWarnings(as.numeric(.data$`Bias Size`)),
        SE = suppressWarnings(as.numeric(.data$`S.E.`)),
        Infit = suppressWarnings(as.numeric(.data$Infit)),
        Outfit = suppressWarnings(as.numeric(.data$Outfit))
      )
  ) |>
    dplyr::mutate(
      LowCountFlag = is.finite(.data$Count) & .data$Count < min_count_warn
    ) |>
    dplyr::arrange(.data$Count)

  by_facet <- stats::setNames(vector("list", length(spec$facets)), spec$facets)
  for (facet in spec$facets) {
    by_facet[[facet]] <- cell_tbl |>
      dplyr::group_by(Level = .data[[facet]]) |>
      dplyr::summarize(
        Cells = dplyr::n(),
        TotalCount = sum(.data$Count, na.rm = TRUE),
        MeanCount = mean(.data$Count, na.rm = TRUE),
        MinCount = min(.data$Count, na.rm = TRUE),
        LowCountCells = sum(.data$LowCountFlag, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::mutate(Facet = facet, .before = 1)
  }

  summary_tbl <- data.frame(
    InteractionFacets = paste(spec$facets, collapse = " x "),
    InteractionOrder = spec$interaction_order,
    InteractionMode = spec$interaction_mode,
    Branch = branch,
    Style = ifelse(branch == "facets", "facets_manual", "original"),
    FacetA = spec$facets[1],
    FacetB = spec$facets[2],
    Cells = nrow(cell_tbl),
    TotalCount = sum(cell_tbl$Count, na.rm = TRUE),
    MeanCount = mean(cell_tbl$Count, na.rm = TRUE),
    MedianCount = stats::median(cell_tbl$Count, na.rm = TRUE),
    MinCount = min(cell_tbl$Count, na.rm = TRUE),
    MaxCount = max(cell_tbl$Count, na.rm = TRUE),
    LowCountCells = sum(cell_tbl$LowCountFlag, na.rm = TRUE),
    LowCountPercent = ifelse(nrow(cell_tbl) > 0, 100 * sum(cell_tbl$LowCountFlag, na.rm = TRUE) / nrow(cell_tbl), NA_real_),
    stringsAsFactors = FALSE
  )

  table_out <- cell_tbl
  if (branch == "facets") {
    table_out <- cell_tbl |>
      dplyr::mutate(
        Sq = dplyr::row_number(),
        `Observd Count` = .data$Count,
        `Obs-Exp Average` = .data$BiasSize,
        `Model S.E.` = .data$SE
      ) |>
      dplyr::select(
        dplyr::all_of(spec$facets),
        "Sq",
        "Observd Count",
        "Obs-Exp Average",
        "Model S.E.",
        "Infit",
        "Outfit",
        "LowCountFlag"
      )
  }

  fit_overview <- data.frame()
  if (!is.null(fit) && inherits(fit, "mfrm_fit") &&
      is.data.frame(fit$summary) && nrow(fit$summary) > 0) {
    fit_overview <- fit$summary[1, , drop = FALSE]
    fit_overview$InteractionFacets <- paste(spec$facets, collapse = " x ")
    fit_overview$Branch <- branch
    fit_overview <- fit_overview[, c("InteractionFacets", "Branch", setdiff(names(fit_overview), c("InteractionFacets", "Branch"))), drop = FALSE]
  }

  out <- list(
    table = table_out,
    by_facet = by_facet,
    by_facet_a = by_facet[[spec$facets[1]]],
    by_facet_b = by_facet[[spec$facets[2]]],
    summary = summary_tbl,
    thresholds = list(min_count_warn = min_count_warn),
    branch = branch,
    style = ifelse(branch == "facets", "facets_manual", "original"),
    fit_overview = fit_overview
  )
  out <- as_mfrm_bundle(out, "mfrm_bias_count")
  class(out) <- unique(c(paste0("mfrm_bias_count_", branch), class(out)))
  out
}

#' Build a FACETS Table 10-style unexpected-after-bias report
#'
#' @param fit Output from [fit_mfrm()].
#' @param bias_results Output from [estimate_bias()].
#' @param diagnostics Optional output from [diagnose_mfrm()] for baseline comparison.
#' @param abs_z_min Absolute standardized-residual cutoff.
#' @param prob_max Maximum observed-category probability cutoff.
#' @param top_n Maximum number of rows to return.
#' @param rule Flagging rule: `"either"` or `"both"`.
#'
#' @details
#' FACETS Table 10 reports responses that remain unexpected after bias terms
#' are introduced. This helper recomputes expected values and residuals using
#' the estimated bias contrasts from [estimate_bias()].
#'
#' `summary(t10)` is supported through `summary()`.
#' `plot(t10)` is dispatched through `plot()` for class
#' `mfrm_unexpected_after_bias` (`type = "scatter"`, `"severity"`,
#' `"comparison"`).
#'
#' @section Interpreting output:
#' - `summary`: before/after unexpected counts and reduction metrics.
#' - `table`: residual unexpected responses after bias adjustment.
#' - `thresholds`: screening settings used in this comparison.
#'
#' Large reductions indicate bias terms explain part of prior unexpectedness;
#' persistent unexpected rows indicate remaining model-data mismatch.
#'
#' @section Typical workflow:
#' 1. Run [unexpected_response_table()] as baseline.
#' 2. Estimate bias via [estimate_bias()].
#' 3. Run `unexpected_after_bias_table(...)` and compare reductions.
#'
#' @section Output columns:
#' The `table` data.frame has the same structure as
#' [unexpected_response_table()] output, with an additional
#' `BiasAdjustment` column showing the bias correction applied to each
#' observation's expected value.
#'
#' The `summary` data.frame contains:
#' \describe{
#'   \item{TotalObservations}{Total observations analyzed.}
#'   \item{BaselineUnexpectedN}{Unexpected count before bias adjustment.}
#'   \item{AfterBiasUnexpectedN}{Unexpected count after adjustment.}
#'   \item{ReducedBy, ReducedPercent}{Reduction in unexpected count.}
#' }
#'
#' @return A named list with:
#' - `table`: unexpected responses after bias adjustment
#' - `summary`: one-row summary (includes baseline-vs-after counts)
#' - `thresholds`: applied thresholds
#' - `facets`: analyzed bias facet pair
#'
#' @seealso [estimate_bias()], [unexpected_response_table()], [bias_count_table()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#' bias <- estimate_bias(fit, diag, facet_a = "Rater", facet_b = "Criterion", max_iter = 2)
#' t10 <- unexpected_after_bias_table(fit, bias, diagnostics = diag, top_n = 20)
#' summary(t10)
#' p_t10 <- plot(t10, draw = FALSE)
#' class(p_t10)
#' @export
unexpected_after_bias_table <- function(fit,
                                        bias_results,
                                        diagnostics = NULL,
                                        abs_z_min = 2,
                                        prob_max = 0.30,
                                        top_n = 100,
                                        rule = c("either", "both")) {
  rule <- match.arg(tolower(rule), c("either", "both"))
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().")
  }
  if (is.null(bias_results) || is.null(bias_results$table) || nrow(bias_results$table) == 0) {
    stop("`bias_results` must be output from estimate_bias() with non-empty `table`.")
  }
  if (is.null(diagnostics)) {
    diagnostics <- diagnose_mfrm(fit, residual_pca = "none")
  }

  obs_adj <- compute_obs_table_with_bias(fit, bias_results = bias_results)
  probs_adj <- compute_prob_matrix_with_bias(fit, bias_results = bias_results)

  tbl <- calc_unexpected_response_table(
    obs_df = obs_adj,
    probs = probs_adj,
    facet_names = fit$config$facet_names,
    rating_min = fit$prep$rating_min,
    abs_z_min = abs_z_min,
    prob_max = prob_max,
    top_n = top_n,
    rule = rule
  )
  if (nrow(tbl) > 0 && "Row" %in% names(tbl) && "BiasAdjustment" %in% names(obs_adj)) {
    adj_df <- data.frame(Row = seq_len(nrow(obs_adj)), BiasAdjustment = obs_adj$BiasAdjustment, stringsAsFactors = FALSE)
    tbl <- tbl |>
      left_join(adj_df, by = "Row")
  }

  summary_tbl <- summarize_unexpected_response_table(
    unexpected_tbl = tbl,
    total_observations = nrow(obs_adj),
    abs_z_min = abs_z_min,
    prob_max = prob_max,
    rule = rule
  )

  baseline <- unexpected_response_table(
    fit = fit,
    diagnostics = diagnostics,
    abs_z_min = abs_z_min,
    prob_max = prob_max,
    top_n = max(top_n, nrow(obs_adj)),
    rule = rule
  )
  baseline_n <- if (is.null(baseline$table)) NA_integer_ else nrow(baseline$table)
  after_n <- nrow(tbl)
  summary_tbl <- summary_tbl |>
    mutate(
      BaselineUnexpectedN = baseline_n,
      AfterBiasUnexpectedN = after_n,
      ReducedBy = ifelse(is.finite(baseline_n), baseline_n - after_n, NA_real_),
      ReducedPercent = ifelse(is.finite(baseline_n) && baseline_n > 0, 100 * (baseline_n - after_n) / baseline_n, NA_real_)
    )

  out <- list(
    table = tbl,
    summary = summary_tbl,
    thresholds = list(
      abs_z_min = abs_z_min,
      prob_max = prob_max,
      rule = rule
    ),
    facets = list(
      facet_a = bias_results$facet_a,
      facet_b = bias_results$facet_b,
      interaction_facets = bias_results$interaction_facets %||% c(bias_results$facet_a, bias_results$facet_b),
      interaction_order = bias_results$interaction_order %||% 2L,
      interaction_mode = bias_results$interaction_mode %||% "pairwise"
    )
  )
  as_mfrm_bundle(out, "mfrm_unexpected_after_bias")
}

resolve_table2_source_columns <- function(fit,
                                          person = NULL,
                                          facets = NULL,
                                          score = NULL,
                                          weight = NULL) {
  source <- fit$config$source_columns
  resolved_person <- if (!is.null(person)) as.character(person[1]) else as.character(source$person %||% "Person")
  resolved_facets <- if (!is.null(facets)) as.character(facets) else as.character(source$facets %||% fit$config$facet_names)
  resolved_score <- if (!is.null(score)) as.character(score[1]) else as.character(source$score %||% "Score")
  resolved_weight <- if (!is.null(weight)) as.character(weight[1]) else as.character(source$weight %||% NA_character_)
  if (length(resolved_weight) == 0 || is.na(resolved_weight) || !nzchar(resolved_weight)) {
    resolved_weight <- NA_character_
  }
  list(
    person = resolved_person,
    facets = resolved_facets,
    score = resolved_score,
    weight = resolved_weight
  )
}

compute_iteration_state <- function(par, idx, prep, config, sizes, quad_points = 15L) {
  params <- expand_params(par, sizes, config)
  if (config$method == "MML") {
    quad <- gauss_hermite_normal(quad_points)
    theta_tbl <- compute_person_eap(idx, config, params, quad)
    theta_diag <- suppressWarnings(as.numeric(theta_tbl$Estimate))
    if (length(theta_diag) != config$n_person || any(!is.finite(theta_diag))) {
      theta_diag <- rep(0, config$n_person)
    }
  } else {
    theta_diag <- suppressWarnings(as.numeric(params$theta))
  }

  eta <- compute_eta(idx, params, config, theta_override = theta_diag)
  if (config$model == "RSM") {
    step_cum <- c(0, cumsum(params$steps))
    probs <- category_prob_rsm(eta, step_cum)
  } else {
    step_cum_mat <- t(apply(params$steps_mat, 1, function(x) c(0, cumsum(x))))
    probs <- category_prob_pcm(eta, step_cum_mat, idx$step_idx)
  }
  k_vals <- 0:(ncol(probs) - 1)
  expected_k <- as.vector(probs %*% k_vals)
  expected <- prep$rating_min + expected_k
  observed <- as.numeric(prep$data$Score)
  weight <- if (!is.null(idx$weight)) as.numeric(idx$weight) else rep(1, length(observed))

  obs_df <- prep$data |>
    mutate(
      Observed = observed,
      Expected = expected,
      Residual = .data$Observed - .data$Expected,
      .Weight = weight
    )

  facet_cols <- c("Person", config$facet_names)
  max_elem_resid <- 0
  for (facet in facet_cols) {
    if (!facet %in% names(obs_df)) next
    sub <- obs_df |>
      group_by(.data[[facet]]) |>
      summarize(ResidualSum = sum(.data$Residual * .data$.Weight, na.rm = TRUE), .groups = "drop")
    if (nrow(sub) == 0) next
    idx_max <- which.max(abs(sub$ResidualSum))
    v <- sub$ResidualSum[idx_max]
    if (is.finite(v) && abs(v) > abs(max_elem_resid)) max_elem_resid <- v
  }
  score_range <- prep$rating_max - prep$rating_min
  max_elem_resid_pct <- ifelse(score_range > 0, 100 * max_elem_resid / score_range, NA_real_)

  obs_k <- observed - prep$rating_min
  n_cat <- ncol(probs)
  obs_count <- rep(0, n_cat)
  valid_k <- is.finite(obs_k) & obs_k >= 0 & obs_k < n_cat
  if (any(valid_k)) {
    obs_idx <- as.integer(obs_k[valid_k]) + 1L
    obs_count <- tapply(weight[valid_k], obs_idx, sum)
    obs_count <- {
      out <- rep(0, n_cat)
      out[as.integer(names(obs_count))] <- as.numeric(obs_count)
      out
    }
  }
  exp_count <- colSums(probs * weight, na.rm = TRUE)
  cat_diff <- obs_count - exp_count
  if (length(cat_diff) == 0 || all(!is.finite(cat_diff))) {
    max_cat_resid <- NA_real_
  } else {
    max_cat_resid <- cat_diff[which.max(abs(cat_diff))]
  }

  facet_params <- unlist(params$facets, use.names = FALSE)
  element_vec <- c(theta_diag, facet_params)
  step_vec <- if (config$model == "RSM") {
    as.numeric(params$steps)
  } else {
    as.numeric(as.vector(params$steps_mat))
  }

  list(
    max_score_resid_elements = as.numeric(max_elem_resid),
    max_score_resid_pct = as.numeric(max_elem_resid_pct),
    max_score_resid_categories = as.numeric(max_cat_resid),
    element_vec = element_vec,
    step_vec = step_vec
  )
}

signal_legacy_name_deprecation <- function(old_name,
                                           new_name,
                                           suppress_if_called_from = NULL,
                                           when = "0.1.0") {
  if (isTRUE(getOption("mfrmr.suppress_legacy_name_warning", FALSE))) {
    return(invisible(NULL))
  }

  caller_call <- tryCatch(sys.call(-2), error = function(e) NULL)
  caller <- ""
  if (!is.null(caller_call) && length(caller_call) > 0) {
    caller <- as.character(caller_call[[1]])
    if (length(caller) == 0 || !is.finite(nchar(caller[1])) || is.na(caller[1])) {
      caller <- ""
    } else {
      caller <- caller[1]
    }
  }
  if (!is.null(suppress_if_called_from) &&
      is.character(suppress_if_called_from) &&
      length(suppress_if_called_from) > 0 &&
      nzchar(caller) &&
      caller %in% suppress_if_called_from) {
    return(invisible(NULL))
  }

  dep_env <- tryCatch(rlang::caller_env(2), error = function(e) rlang::caller_env())
  user_env <- tryCatch(rlang::caller_env(3), error = function(e) dep_env)
  lifecycle::deprecate_soft(
    when = when,
    what = paste0(old_name, "()"),
    with = paste0(new_name, "()"),
    id = paste0("mfrmr_", old_name, "_legacy_name"),
    env = dep_env,
    user_env = user_env
  )
  invisible(NULL)
}

with_legacy_name_warning_suppressed <- function(expr) {
  old_opt <- options(mfrmr.suppress_legacy_name_warning = TRUE)
  on.exit(options(old_opt), add = TRUE)
  eval.parent(substitute(expr))
}

as_mfrm_bundle <- function(x, class_name) {
  if (!is.list(x)) return(x)
  class(x) <- unique(c(class_name, "mfrm_bundle", class(x)))
  x
}

#' Build a FACETS Table 1-style specification summary
#'
#' @param fit Output from [fit_mfrm()].
#' @param title Optional analysis title.
#' @param data_file Optional data-file label (for reporting only).
#' @param output_file Optional output-file label (for reporting only).
#' @param include_fixed If `TRUE`, include a FACETS-style fixed-width text block.
#'
#' @details
#' FACETS Table 1 groups model settings by function (data, output, convergence).
#' This helper assembles those settings from a fitted object.
#'
#' @section Lifecycle:
#' Soft-deprecated. Prefer [specifications_report()].
#'
#' @section Output columns:
#' The `header` data.frame contains:
#' \describe{
#'   \item{Engine}{Estimation engine identifier (always `"mfrmr"`).}
#'   \item{Model}{`"RSM"` or `"PCM"`.}
#'   \item{Method}{`"JMLE"` or `"MML"`.}
#' }
#'
#' The `facet_labels` data.frame contains:
#' \describe{
#'   \item{Facet}{Facet column name.}
#'   \item{Elements}{Number of levels in this facet.}
#' }
#'
#' @return A named list with:
#' - `header`: one-row run header
#' - `data_spec`: key-value table for data/model settings
#' - `facet_labels`: facet names with level counts
#' - `output_spec`: key-value table for reporting-related defaults
#' - `convergence_control`: key-value table for optimizer controls/results
#' - `anchor_summary`: anchor/group-anchor summary by facet
#' - `fixed`: fixed-width report text (when `include_fixed = TRUE`)
#'
#' @seealso [fit_mfrm()], [data_quality_report()], [estimation_iteration_report()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' t1 <- specifications_report(fit, title = "Toy run")
#' @keywords internal
#' @noRd
table1_specifications <- function(fit,
                                  title = NULL,
                                  data_file = NULL,
                                  output_file = NULL,
                                  include_fixed = FALSE) {
  signal_legacy_name_deprecation(
    old_name = "table1_specifications",
    new_name = "specifications_report",
    suppress_if_called_from = "specifications_report"
  )
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().")
  }

  cfg <- fit$config
  prep <- fit$prep
  ov <- as.data.frame(fit$summary, stringsAsFactors = FALSE)
  if (nrow(ov) == 0) stop("`fit$summary` is empty.")

  header <- data.frame(
    Engine = paste0("mfrmr ", as.character(utils::packageVersion("mfrmr"))),
    Title = ifelse(is.null(title), "", as.character(title[1])),
    DataFile = ifelse(is.null(data_file), "", as.character(data_file[1])),
    OutputFile = ifelse(is.null(output_file), "", as.character(output_file[1])),
    Model = as.character(ov$Model[1]),
    Method = as.character(ov$Method[1]),
    stringsAsFactors = FALSE
  )

  all_facets <- c("Person", cfg$facet_names)
  facet_labels <- data.frame(
    FacetIndex = seq_along(all_facets),
    Facet = all_facets,
    Elements = vapply(prep$levels[all_facets], length, integer(1)),
    stringsAsFactors = FALSE
  )

  est_ctl <- cfg$estimation_control %||% list()
  positive <- cfg$positive_facets %||% character(0)
  dummy <- cfg$dummy_facets %||% character(0)
  data_spec <- data.frame(
    Setting = c(
      "Facets",
      "Persons",
      "Categories",
      "RatingMin",
      "RatingMax",
      "NonCenteredFacet",
      "PositiveFacets",
      "DummyFacets",
      "StepFacet",
      "WeightColumn"
    ),
    Value = c(
      as.character(length(cfg$facet_names)),
      as.character(cfg$n_person),
      as.character(cfg$n_cat),
      as.character(prep$rating_min),
      as.character(prep$rating_max),
      as.character(cfg$noncenter_facet %||% ""),
      ifelse(length(positive) == 0, "", paste(positive, collapse = ", ")),
      ifelse(length(dummy) == 0, "", paste(dummy, collapse = ", ")),
      as.character(cfg$step_facet %||% ""),
      as.character(cfg$weight_col %||% "")
    ),
    stringsAsFactors = FALSE
  )

  output_spec <- data.frame(
    Setting = c(
      "UnexpectedAbsZThreshold",
      "UnexpectedProbThreshold",
      "DisplacementWarnLogit",
      "DisplacementWarnT",
      "FairScoreDefault"
    ),
    Value = c("2", "0.30", "0.5", "2", "Mean"),
    stringsAsFactors = FALSE
  )

  convergence_control <- data.frame(
    Setting = c(
      "MaxIterations",
      "RelativeTolerance",
      "QuadPoints",
      "Converged",
      "FunctionEvaluations",
      "OptimizerMessage"
    ),
    Value = c(
      as.character(est_ctl$maxit %||% NA_integer_),
      as.character(est_ctl$reltol %||% NA_real_),
      as.character(est_ctl$quad_points %||% NA_integer_),
      as.character(isTRUE(ov$Converged[1])),
      as.character(fit$opt$counts[["function"]] %||% NA_integer_),
      as.character(fit$opt$message %||% "")
    ),
    stringsAsFactors = FALSE
  )

  anchor_summary <- if (!is.null(cfg$anchor_summary)) {
    as.data.frame(cfg$anchor_summary, stringsAsFactors = FALSE)
  } else {
    data.frame()
  }

  out <- list(
    header = header,
    data_spec = data_spec,
    facet_labels = facet_labels,
    output_spec = output_spec,
    convergence_control = convergence_control,
    anchor_summary = anchor_summary
  )

  if (isTRUE(include_fixed)) {
    out$fixed <- build_sectioned_fixed_report(
      title = "FACETS-style Table 1 Specification Summary",
      sections = list(
        list(title = "Header", data = header),
        list(title = "Data specification", data = data_spec),
        list(title = "Facet labels", data = facet_labels),
        list(title = "Output specification", data = output_spec),
        list(title = "Convergence control", data = convergence_control),
        list(title = "Anchor summary", data = anchor_summary)
      ),
      max_col_width = 24
    )
  }
  out
}

#' Build a FACETS Table 2-style data summary report
#'
#' @param fit Output from [fit_mfrm()].
#' @param data Optional raw data frame used for additional row-level audit.
#' @param person Optional person column name in `data`.
#' @param facets Optional facet column names in `data`.
#' @param score Optional score column name in `data`.
#' @param weight Optional weight column name in `data`.
#' @param include_fixed If `TRUE`, include a FACETS-style fixed-width text block.
#'
#' @details
#' When `data` is supplied, this function performs row-level validity checks
#' (missing identifiers, missing score, non-positive weight, out-of-range score)
#' and reports dropped rows similarly to FACETS Table 2 diagnostics.
#'
#' @section Lifecycle:
#' Soft-deprecated. Prefer [data_quality_report()].
#'
#' @section Output columns:
#' The `summary` data.frame contains:
#' \describe{
#'   \item{TotalLinesInData}{Total rows in the input data.}
#'   \item{ValidResponsesUsedForEstimation}{Rows used after filtering.}
#'   \item{MissingScoreRows, MissingFacetRows, MissingPersonRows}{Counts
#'     of rows excluded due to missing values.}
#' }
#'
#' The `row_audit` data.frame contains:
#' \describe{
#'   \item{Status}{Row-status category (e.g., "Valid", "MissingScore").}
#'   \item{N}{Number of rows in this category.}
#' }
#'
#' @return A named list with:
#' - `summary`: one-row data summary
#' - `model_match`: model-match counts
#' - `row_audit`: counts by row-status category
#' - `unknown_elements`: elements in `data` not present in fitted levels
#' - `category_counts`: observed category counts in fitted data
#' - `fixed`: fixed-width report text (when `include_fixed = TRUE`)
#'
#' @seealso [fit_mfrm()], [specifications_report()], [describe_mfrm_data()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' t2 <- data_quality_report(
#'   fit, data = toy, person = "Person",
#'   facets = c("Rater", "Criterion"),
#'   score = "Score"
#' )
#' @keywords internal
#' @noRd
table2_data_summary <- function(fit,
                                data = NULL,
                                person = NULL,
                                facets = NULL,
                                score = NULL,
                                weight = NULL,
                                include_fixed = FALSE) {
  signal_legacy_name_deprecation(
    old_name = "table2_data_summary",
    new_name = "data_quality_report",
    suppress_if_called_from = "data_quality_report"
  )
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().")
  }

  prep <- fit$prep
  cfg <- fit$config
  src <- resolve_table2_source_columns(
    fit = fit,
    person = person,
    facets = facets,
    score = score,
    weight = weight
  )

  fitted_df <- prep$data
  valid_used <- nrow(fitted_df)
  category_counts <- fitted_df |>
    group_by(.data$Score) |>
    summarize(
      Count = n(),
      WeightedCount = sum(.data$Weight, na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(.data$Score)

  if (is.null(data)) {
    summary_tbl <- data.frame(
      TotalLinesInData = NA_integer_,
      TotalDataLines = NA_integer_,
      TotalNonBlankResponsesFound = NA_integer_,
      ValidResponsesUsedForEstimation = valid_used,
      stringsAsFactors = FALSE
    )
    model_match <- data.frame(
      Model = paste(cfg$model, cfg$method, sep = " / "),
      MatchedResponses = valid_used,
      stringsAsFactors = FALSE
    )
    out <- list(
      summary = summary_tbl,
      model_match = model_match,
      row_audit = data.frame(),
      unknown_elements = data.frame(),
      category_counts = as.data.frame(category_counts, stringsAsFactors = FALSE)
    )
    if (isTRUE(include_fixed)) {
      out$fixed <- build_sectioned_fixed_report(
        title = "FACETS-style Table 2 Data Summary",
        sections = list(
          list(title = "Summary", data = summary_tbl),
          list(title = "Model match", data = model_match),
          list(title = "Row audit", data = data.frame()),
          list(title = "Unknown elements", data = data.frame()),
          list(title = "Category counts", data = as.data.frame(category_counts, stringsAsFactors = FALSE))
        ),
        max_col_width = 24
      )
    }
    return(out)
  }

  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.")
  }
  required <- c(src$person, src$facets, src$score)
  if (!all(required %in% names(data))) {
    stop("`data` is missing required columns: ", paste(setdiff(required, names(data)), collapse = ", "))
  }
  if (!is.na(src$weight) && !src$weight %in% names(data)) {
    stop("`weight` column not found in `data`.")
  }

  raw <- data
  n_total <- nrow(raw)
  person_ok <- !is.na(raw[[src$person]])
  facet_ok <- if (length(src$facets) == 0) rep(TRUE, n_total) else {
    apply(!is.na(raw[, src$facets, drop = FALSE]), 1, all)
  }
  score_num <- suppressWarnings(as.numeric(raw[[src$score]]))
  score_nonblank <- !is.na(score_num)
  if (!is.na(src$weight)) {
    w_num <- suppressWarnings(as.numeric(raw[[src$weight]]))
    weight_ok <- is.finite(w_num) & w_num > 0
  } else {
    weight_ok <- rep(TRUE, n_total)
  }
  range_ok <- score_nonblank & score_num >= prep$rating_min & score_num <= prep$rating_max

  row_status <- rep("valid", n_total)
  row_status[!person_ok] <- "missing_person"
  row_status[person_ok & !facet_ok] <- "missing_facet"
  row_status[person_ok & facet_ok & !score_nonblank] <- "missing_score"
  row_status[person_ok & facet_ok & score_nonblank & !weight_ok] <- "invalid_weight"
  row_status[person_ok & facet_ok & score_nonblank & weight_ok & !range_ok] <- "score_out_of_range"

  row_audit <- data.frame(Status = row_status, stringsAsFactors = FALSE) |>
    count(.data$Status, name = "N") |>
    arrange(desc(.data$N), .data$Status)

  model_match <- data.frame(
    Model = paste(cfg$model, cfg$method, sep = " / "),
    MatchedResponses = sum(person_ok & facet_ok & score_nonblank, na.rm = TRUE),
    ValidResponses = sum(row_status == "valid", na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  unknown_rows <- list()
  check_facets <- c("Person", cfg$facet_names)
  raw_names <- c(src$person, src$facets)
  names(raw_names) <- check_facets
  for (facet in check_facets) {
    col <- raw_names[[facet]]
    if (is.null(col) || !col %in% names(raw)) next
    known <- as.character(prep$levels[[facet]])
    vals <- as.character(raw[[col]])
    bad <- vals[!is.na(vals) & !(vals %in% known)]
    if (length(bad) == 0) next
    unknown_rows[[facet]] <- data.frame(
      Facet = facet,
      Level = sort(unique(bad)),
      stringsAsFactors = FALSE
    )
  }
  unknown_tbl <- if (length(unknown_rows) == 0) data.frame() else bind_rows(unknown_rows)

  summary_tbl <- data.frame(
    TotalLinesInData = n_total,
    TotalDataLines = n_total,
    TotalNonBlankResponsesFound = sum(score_nonblank, na.rm = TRUE),
    MissingScoreRows = sum(row_status == "missing_score", na.rm = TRUE),
    MissingFacetRows = sum(row_status == "missing_facet", na.rm = TRUE),
    MissingPersonRows = sum(row_status == "missing_person", na.rm = TRUE),
    InvalidWeightRows = sum(row_status == "invalid_weight", na.rm = TRUE),
    OutOfRangeScoreRows = sum(row_status == "score_out_of_range", na.rm = TRUE),
    ValidResponsesUsedForEstimation = valid_used,
    stringsAsFactors = FALSE
  )

  out <- list(
    summary = summary_tbl,
    model_match = model_match,
    row_audit = row_audit,
    unknown_elements = unknown_tbl,
    category_counts = as.data.frame(category_counts, stringsAsFactors = FALSE)
  )
  if (isTRUE(include_fixed)) {
    out$fixed <- build_sectioned_fixed_report(
      title = "FACETS-style Table 2 Data Summary",
      sections = list(
        list(title = "Summary", data = summary_tbl),
        list(title = "Model match", data = model_match),
        list(title = "Row audit", data = row_audit),
        list(title = "Unknown elements", data = unknown_tbl),
        list(title = "Category counts", data = as.data.frame(category_counts, stringsAsFactors = FALSE))
      ),
      max_col_width = 24
    )
  }
  out
}

#' Build a FACETS Table 3-style iteration report
#'
#' @param fit Output from [fit_mfrm()].
#' @param max_iter Maximum replay iterations (excluding optional initial row).
#' @param reltol Stopping tolerance for replayed max-logit change.
#' @param include_prox If `TRUE`, include an initial pseudo-row labeled `PROX`.
#' @param include_fixed If `TRUE`, include a FACETS-style fixed-width text block.
#'
#' @details
#' FACETS prints per-iteration score residual and logit-change diagnostics.
#' The underlying optimizer in this package does not expose the exact internal
#' per-iteration path, so this function reconstructs an approximation by
#' repeatedly running one-iteration updates from the current parameter vector.
#'
#' @section Lifecycle:
#' Soft-deprecated. Prefer [estimation_iteration_report()].
#'
#' @section Output columns:
#' The `table` data.frame contains:
#' \describe{
#'   \item{Method}{Estimation method label.}
#'   \item{Iteration}{Iteration number.}
#'   \item{MaxLogitChangeElements}{Largest parameter change across
#'     elements in this iteration.}
#'   \item{MaxLogitChangeSteps}{Largest step-parameter change.}
#'   \item{Objective}{Objective function value (negative log-likelihood).}
#' }
#'
#' The `summary` data.frame contains:
#' \describe{
#'   \item{FinalConverged}{Logical; `TRUE` if the optimizer converged.}
#'   \item{FinalIterations}{Number of iterations used.}
#' }
#'
#' @return A named list with:
#' - `table`: iteration rows (residual and logit-change metrics)
#' - `summary`: one-row convergence summary
#' - `settings`: replay settings used
#' - `fixed`: fixed-width report text (when `include_fixed = TRUE`)
#'
#' @seealso [fit_mfrm()], [specifications_report()], [data_quality_report()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' t3 <- estimation_iteration_report(fit, max_iter = 5)
#' @keywords internal
#' @noRd
table3_iteration_report <- function(fit,
                                    max_iter = 20,
                                    reltol = NULL,
                                    include_prox = TRUE,
                                    include_fixed = FALSE) {
  signal_legacy_name_deprecation(
    old_name = "table3_iteration_report",
    new_name = "estimation_iteration_report",
    suppress_if_called_from = "estimation_iteration_report"
  )
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().")
  }
  cfg <- fit$config
  prep <- fit$prep
  sizes <- build_param_sizes(cfg)
  idx <- build_indices(prep, step_facet = cfg$step_facet)
  est_ctl <- cfg$estimation_control %||% list()
  if (is.null(reltol) || !is.finite(reltol)) reltol <- as.numeric(est_ctl$reltol %||% 1e-6)
  quad_points <- as.integer(est_ctl$quad_points %||% 15L)
  max_iter <- max(1L, as.integer(max_iter))

  current_par <- build_initial_param_vector(cfg, sizes)
  prev_state <- compute_iteration_state(
    par = current_par,
    idx = idx,
    prep = prep,
    config = cfg,
    sizes = sizes,
    quad_points = quad_points
  )

  rows <- list()
  row_id <- 1L
  if (isTRUE(include_prox)) {
    rows[[row_id]] <- tibble(
      Method = "PROX",
      Iteration = 1L,
      MaxScoreResidualElements = prev_state$max_score_resid_elements,
      MaxScoreResidualPercent = prev_state$max_score_resid_pct,
      MaxScoreResidualCategories = prev_state$max_score_resid_categories,
      MaxLogitChangeElements = NA_real_,
      MaxLogitChangeSteps = NA_real_,
      Objective = NA_real_
    )
    row_id <- row_id + 1L
  }

  for (it in seq_len(max_iter)) {
    opt_step <- run_mfrm_optimization(
      start = current_par,
      method = cfg$method,
      idx = idx,
      config = cfg,
      sizes = sizes,
      quad_points = quad_points,
      maxit = 1L,
      reltol = reltol
    )
    current_par <- opt_step$par
    state <- compute_iteration_state(
      par = current_par,
      idx = idx,
      prep = prep,
      config = cfg,
      sizes = sizes,
      quad_points = quad_points
    )

    elem_change <- if (length(state$element_vec) == length(prev_state$element_vec) && length(state$element_vec) > 0) {
      max(abs(state$element_vec - prev_state$element_vec), na.rm = TRUE)
    } else {
      NA_real_
    }
    step_change <- if (length(state$step_vec) == length(prev_state$step_vec) && length(state$step_vec) > 0) {
      max(abs(state$step_vec - prev_state$step_vec), na.rm = TRUE)
    } else {
      NA_real_
    }

    rows[[row_id]] <- tibble(
      Method = cfg$method,
      Iteration = if (isTRUE(include_prox)) it + 1L else it,
      MaxScoreResidualElements = state$max_score_resid_elements,
      MaxScoreResidualPercent = state$max_score_resid_pct,
      MaxScoreResidualCategories = state$max_score_resid_categories,
      MaxLogitChangeElements = elem_change,
      MaxLogitChangeSteps = step_change,
      Objective = -opt_step$value
    )
    row_id <- row_id + 1L
    prev_state <- state

    change_vec <- c(elem_change, step_change)
    change_vec <- change_vec[is.finite(change_vec)]
    max_change <- if (length(change_vec) == 0) NA_real_ else max(change_vec)
    if (is.finite(max_change) && max_change < reltol) break
  }

  tbl <- bind_rows(rows)
  subset_tbl <- calc_subsets(compute_obs_table(fit), c("Person", cfg$facet_names))$summary
  connected <- if (!is.null(subset_tbl) && nrow(subset_tbl) > 0) nrow(subset_tbl) == 1 else NA

  summary_tbl <- data.frame(
    FinalConverged = isTRUE(fit$summary$Converged[1]),
    FinalIterations = as.integer(fit$summary$Iterations[1]),
    ReplayRows = nrow(tbl),
    ConnectedSubset = connected,
    stringsAsFactors = FALSE
  )

  out <- list(
    table = as.data.frame(tbl, stringsAsFactors = FALSE),
    summary = summary_tbl,
    settings = list(
      max_iter = max_iter,
      reltol = reltol,
      include_prox = include_prox,
      quad_points = quad_points,
      include_fixed = isTRUE(include_fixed)
    )
  )
  if (isTRUE(include_fixed)) {
    out$fixed <- build_sectioned_fixed_report(
      title = "FACETS-style Table 3 Iteration Report",
      sections = list(
        list(title = "Iteration rows", data = as.data.frame(tbl, stringsAsFactors = FALSE), max_rows = 200L),
        list(title = "Summary", data = summary_tbl),
        list(title = "Settings", data = as.data.frame(out$settings, stringsAsFactors = FALSE))
      ),
      max_col_width = 24
    )
  }
  out
}

#' Build a FACETS Table 6.0.0-style subset/disjoint-element listing
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()].
#' @param top_n_subsets Optional maximum number of subset rows to keep.
#' @param min_observations Minimum observations required to keep a subset row.
#'
#' @details
#' FACETS Table 6.0.0 reports disjoint subsets when the design is not fully
#' connected. This helper exposes the same design-check idea with:
#' 1) subset-level counts and percentages, and
#' 2) per-subset/per-facet element listings.
#'
#' @section Lifecycle:
#' Soft-deprecated. Prefer [subset_connectivity_report()].
#'
#' @section Output columns:
#' The `summary` data.frame contains:
#' \describe{
#'   \item{Subset}{Subset index.}
#'   \item{Observations}{Number of observations in this subset.}
#'   \item{ObservationPercent}{Percentage of total observations.}
#' }
#'
#' The `listing` data.frame contains:
#' \describe{
#'   \item{Subset}{Subset index.}
#'   \item{Facet}{Facet name.}
#'   \item{LevelsN}{Number of levels in this facet within this subset.}
#'   \item{Levels}{Comma-separated level labels.}
#'   \item{Ruler}{ASCII ruler for visual comparison.}
#' }
#'
#' @return A named list with:
#' - `summary`: subset-level counts (including `ObservationPercent`)
#' - `listing`: facet-level element listing by subset
#' - `nodes`: node-level table from subset detection
#' - `settings`: applied filters and flags
#'
#' @seealso [diagnose_mfrm()], [measurable_summary_table()], [data_quality_report()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' t6 <- subset_connectivity_report(fit)
#' @keywords internal
#' @noRd
table6_subsets_listing <- function(fit,
                                   diagnostics = NULL,
                                   top_n_subsets = NULL,
                                   min_observations = 0) {
  signal_legacy_name_deprecation(
    old_name = "table6_subsets_listing",
    new_name = "subset_connectivity_report",
    suppress_if_called_from = "subset_connectivity_report"
  )
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().")
  }
  if (is.null(diagnostics)) {
    diagnostics <- diagnose_mfrm(fit, residual_pca = "none")
  }

  summary_tbl <- if (!is.null(diagnostics$subsets$summary)) {
    as.data.frame(diagnostics$subsets$summary, stringsAsFactors = FALSE)
  } else {
    data.frame()
  }
  nodes_tbl <- if (!is.null(diagnostics$subsets$nodes)) {
    as.data.frame(diagnostics$subsets$nodes, stringsAsFactors = FALSE)
  } else {
    data.frame()
  }

  if (nrow(summary_tbl) == 0 || nrow(nodes_tbl) == 0) {
    return(list(
      summary = summary_tbl,
      listing = data.frame(),
      nodes = nodes_tbl,
      settings = list(
        top_n_subsets = if (is.null(top_n_subsets)) NA_integer_ else max(1L, as.integer(top_n_subsets)),
        min_observations = as.numeric(min_observations),
        is_disjoint = FALSE
      )
    ))
  }

  summary_tbl$Subset <- suppressWarnings(as.integer(summary_tbl$Subset))
  summary_tbl$Observations <- suppressWarnings(as.numeric(summary_tbl$Observations))
  total_obs <- sum(summary_tbl$Observations, na.rm = TRUE)
  summary_tbl$ObservationPercent <- ifelse(
    is.finite(total_obs) && total_obs > 0,
    100 * summary_tbl$Observations / total_obs,
    NA_real_
  )
  summary_tbl <- summary_tbl |>
    dplyr::arrange(dplyr::desc(.data$Observations), .data$Subset)

  min_observations <- as.numeric(min_observations)
  if (is.finite(min_observations) && min_observations > 0) {
    summary_tbl <- summary_tbl |>
      dplyr::filter(.data$Observations >= min_observations)
  }
  if (!is.null(top_n_subsets)) {
    summary_tbl <- summary_tbl |>
      dplyr::slice_head(n = max(1L, as.integer(top_n_subsets)))
  }

  keep_subsets <- unique(summary_tbl$Subset)
  nodes_tbl <- nodes_tbl |>
    dplyr::mutate(Subset = suppressWarnings(as.integer(.data$Subset))) |>
    dplyr::filter(.data$Subset %in% keep_subsets)

  listing_tbl <- if (nrow(nodes_tbl) == 0) {
    data.frame()
  } else {
    nodes_tbl |>
      dplyr::group_by(.data$Subset, .data$Facet) |>
      dplyr::summarise(
        LevelsN = dplyr::n_distinct(.data$Level),
        Levels = paste(sort(unique(as.character(.data$Level))), collapse = ", "),
        .groups = "drop"
      ) |>
      dplyr::left_join(
        summary_tbl |>
          dplyr::select("Subset", "Observations", "ObservationPercent"),
        by = "Subset"
      ) |>
      dplyr::group_by(.data$Facet) |>
      dplyr::mutate(MaxLevelsN = max(.data$LevelsN, na.rm = TRUE)) |>
      dplyr::ungroup() |>
      dplyr::mutate(
        Ruler = purrr::map2_chr(
          .data$LevelsN,
          .data$MaxLevelsN,
          function(levels_n, max_levels_n) {
            width <- 20L
            if (!is.finite(max_levels_n) || max_levels_n <= 0 || !is.finite(levels_n)) {
              return("")
            }
            fill <- as.integer(round(width * levels_n / max_levels_n))
            fill <- min(width, max(1L, fill))
            paste0("[", strrep("=", fill), strrep(".", width - fill), "]")
          }
        )
      ) |>
      dplyr::select(-"MaxLevelsN") |>
      dplyr::arrange(.data$Subset, .data$Facet)
  }

  list(
    summary = as.data.frame(summary_tbl, stringsAsFactors = FALSE),
    listing = as.data.frame(listing_tbl, stringsAsFactors = FALSE),
    nodes = as.data.frame(nodes_tbl, stringsAsFactors = FALSE),
    settings = list(
      top_n_subsets = if (is.null(top_n_subsets)) NA_integer_ else max(1L, as.integer(top_n_subsets)),
      min_observations = min_observations,
      is_disjoint = nrow(summary_tbl) > 1
    )
  )
}

#' Build a FACETS Table 6.2-style facet-statistics graphic summary
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()].
#' @param metrics Numeric columns in `diagnostics$measures` to summarize.
#' @param ruler_width Width of the fixed-width ruler used for `M/S/Q/X` marks.
#'
#' @details
#' FACETS Table 6.2 describes each facet with a compact graphical ruler where:
#' - `M` marks the mean
#' - `S` marks one SD from the mean
#' - `Q` marks two SD from the mean
#' - `X` marks three SD from the mean
#'
#' Rulers for a given metric share the same min/max scale across facets.
#'
#' @section Lifecycle:
#' Soft-deprecated. Prefer [facet_statistics_report()].
#'
#' @section Output columns:
#' The `table` data.frame contains:
#' \describe{
#'   \item{Metric}{Statistic name (e.g., "Measure", "Infit", "Outfit").}
#'   \item{Facet}{Facet name.}
#'   \item{Levels}{Number of levels.}
#'   \item{Mean, SD, Min, Max}{Summary statistics for this metric-facet
#'     combination.}
#'   \item{Ruler}{ASCII ruler string showing the distribution of level
#'     values. Markers: M = mean, S = +/- 1 SD, Q = +/- 2 SD.}
#' }
#'
#' @return A named list with:
#' - `table`: facet-by-metric summary with ruler strings
#' - `ranges`: global metric ranges used to draw rulers
#' - `settings`: applied metrics and ruler settings
#'
#' @seealso [diagnose_mfrm()], [summary.mfrm_fit()], [plot_facets_chisq()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' t62 <- facet_statistics_report(fit)
#' @keywords internal
#' @noRd
table6_2_facet_statistics <- function(fit,
                                      diagnostics = NULL,
                                      metrics = c("Estimate", "SE", "Infit", "Outfit"),
                                      ruler_width = 41) {
  signal_legacy_name_deprecation(
    old_name = "table6_2_facet_statistics",
    new_name = "facet_statistics_report",
    suppress_if_called_from = "facet_statistics_report"
  )
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().")
  }
  if (is.null(diagnostics)) {
    diagnostics <- diagnose_mfrm(fit, residual_pca = "none")
  }
  if (is.null(diagnostics$measures) || nrow(diagnostics$measures) == 0) {
    stop("`diagnostics$measures` is empty. Run diagnose_mfrm() first.")
  }

  measure_tbl <- as.data.frame(diagnostics$measures, stringsAsFactors = FALSE)
  if (!"Facet" %in% names(measure_tbl)) {
    stop("`diagnostics$measures` must include a `Facet` column.")
  }

  metrics <- unique(as.character(metrics))
  available_metrics <- metrics[metrics %in% names(measure_tbl)]
  if (length(available_metrics) == 0) {
    stop("None of `metrics` were found in `diagnostics$measures`.")
  }
  numeric_metrics <- available_metrics[
    vapply(measure_tbl[available_metrics], is.numeric, logical(1))
  ]
  if (length(numeric_metrics) == 0) {
    stop("Selected `metrics` must be numeric columns in `diagnostics$measures`.")
  }

  ruler_width <- max(21L, as.integer(ruler_width))

  metric_ranges <- purrr::map_dfr(numeric_metrics, function(metric) {
    vals <- suppressWarnings(as.numeric(measure_tbl[[metric]]))
    tibble::tibble(
      Metric = metric,
      GlobalMin = if (any(is.finite(vals))) min(vals, na.rm = TRUE) else NA_real_,
      GlobalMax = if (any(is.finite(vals))) max(vals, na.rm = TRUE) else NA_real_
    )
  })

  make_ruler <- function(mean_val, sd_val, global_min, global_max, width) {
    chars <- rep(".", width)
    priority <- c("." = 0L, "X" = 1L, "Q" = 2L, "S" = 3L, "M" = 4L)

    scale_pos <- function(x) {
      if (!is.finite(x) || !is.finite(global_min) || !is.finite(global_max)) {
        return(NA_integer_)
      }
      if (global_max <= global_min) {
        return(as.integer((width + 1L) / 2L))
      }
      pos <- 1L + as.integer(round((width - 1L) * (x - global_min) / (global_max - global_min)))
      as.integer(min(width, max(1L, pos)))
    }

    place_marker <- function(x, marker) {
      pos <- scale_pos(x)
      if (!is.finite(pos)) return(invisible(NULL))
      current <- chars[pos]
      if (priority[[marker]] >= priority[[current]]) {
        chars[pos] <<- marker
      }
      invisible(NULL)
    }

    if (is.finite(mean_val) && is.finite(sd_val) && sd_val >= 0) {
      for (k in c(-3, 3)) place_marker(mean_val + k * sd_val, "X")
      for (k in c(-2, 2)) place_marker(mean_val + k * sd_val, "Q")
      for (k in c(-1, 1)) place_marker(mean_val + k * sd_val, "S")
      place_marker(mean_val, "M")
    } else if (is.finite(mean_val)) {
      place_marker(mean_val, "M")
    }

    paste0("[", paste(chars, collapse = ""), "]")
  }

  out_tbl <- purrr::map_dfr(numeric_metrics, function(metric) {
    gmin <- metric_ranges$GlobalMin[match(metric, metric_ranges$Metric)]
    gmax <- metric_ranges$GlobalMax[match(metric, metric_ranges$Metric)]

    measure_tbl |>
      dplyr::group_by(.data$Facet) |>
      dplyr::summarise(
        Levels = dplyr::n(),
        Mean = mean(.data[[metric]], na.rm = TRUE),
        SD = stats::sd(.data[[metric]], na.rm = TRUE),
        Min = min(.data[[metric]], na.rm = TRUE),
        Max = max(.data[[metric]], na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::mutate(
        Metric = metric,
        GlobalMin = gmin,
        GlobalMax = gmax,
        Ruler = purrr::pmap_chr(
          list(.data$Mean, .data$SD, .data$GlobalMin, .data$GlobalMax),
          function(mean_val, sd_val, global_min, global_max) {
            make_ruler(
              mean_val = mean_val,
              sd_val = sd_val,
              global_min = global_min,
              global_max = global_max,
              width = ruler_width
            )
          }
        ),
        .before = 1
      )
  }) |>
    dplyr::arrange(.data$Metric, .data$Facet)

  list(
    table = as.data.frame(out_tbl, stringsAsFactors = FALSE),
    ranges = as.data.frame(metric_ranges, stringsAsFactors = FALSE),
    settings = list(
      metrics = numeric_metrics,
      ruler_width = ruler_width,
      marker_legend = c(M = "mean", S = "+/-1 SD", Q = "+/-2 SD", X = "+/-3 SD")
    )
  )
}

closest_theta_for_target <- function(theta, y, target) {
  if (length(theta) == 0 || length(y) == 0 || !is.finite(target)) return(NA_real_)
  ok <- is.finite(theta) & is.finite(y)
  if (!any(ok)) return(NA_real_)
  theta <- theta[ok]
  y <- y[ok]
  theta[which.min(abs(y - target))]
}

#' Build a FACETS Table 8 bar-chart style scale-structure export
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()].
#' @param theta_range Theta/logit range used to derive mode/mean transition points.
#' @param theta_points Number of grid points used for transition-point search.
#' @param drop_unused If `TRUE`, remove zero-count categories from `category_table`.
#' @param include_fixed If `TRUE`, include a FACETS-style fixed-width text block.
#' @param fixed_max_rows Maximum rows per section in the fixed-width text output.
#'
#' @details
#' FACETS Table 8 bar-chart output describes category structure with
#' mode / median / mean landmarks along the latent scale.
#' This helper returns those landmarks as numeric tables:
#' - mode peaks and mode transition points
#' - median thresholds (step calibrations)
#' - mean half-score transition points
#'
#' @section Lifecycle:
#' Soft-deprecated. Prefer [category_structure_report()].
#'
#' @return A named list with:
#' - `category_table`: observed/expected category counts and fit
#' - `mode_peaks`: peak theta/probability by group and category
#' - `mode_boundaries`: theta points where modal category changes
#' - `median_thresholds`: threshold table (step-based)
#' - `mean_halfscore_points`: theta points where expected score crosses half-scores
#' - `fixed`: fixed-width report text (when `include_fixed = TRUE`)
#' - `settings`: applied options
#'
#' @seealso [rating_scale_table()], [category_curves_report()], [plot.mfrm_fit()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' t8b <- category_structure_report(fit)
#' @keywords internal
#' @noRd
table8_barchart_export <- function(fit,
                                   diagnostics = NULL,
                                   theta_range = c(-6, 6),
                                   theta_points = 241,
                                   drop_unused = FALSE,
                                   include_fixed = FALSE,
                                   fixed_max_rows = 200) {
  signal_legacy_name_deprecation(
    old_name = "table8_barchart_export",
    new_name = "category_structure_report",
    suppress_if_called_from = "category_structure_report"
  )
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().")
  }
  if (is.null(diagnostics)) {
    diagnostics <- diagnose_mfrm(fit, residual_pca = "none")
  }
  if (is.null(diagnostics$obs) || nrow(diagnostics$obs) == 0) {
    stop("`diagnostics$obs` is empty. Run diagnose_mfrm() first.")
  }
  theta_points <- max(51L, as.integer(theta_points))
  theta_range <- as.numeric(theta_range)
  if (length(theta_range) != 2 || !all(is.finite(theta_range)) || theta_range[1] >= theta_range[2]) {
    stop("`theta_range` must be a numeric length-2 vector with increasing values.")
  }

  category_table <- as.data.frame(calc_category_stats(diagnostics$obs, res = fit, whexact = FALSE), stringsAsFactors = FALSE)
  if (isTRUE(drop_unused) && nrow(category_table) > 0 && "Count" %in% names(category_table)) {
    category_table <- category_table[category_table$Count > 0, , drop = FALSE]
  }

  curve_spec <- build_step_curve_spec(fit)
  theta_grid <- seq(theta_range[1], theta_range[2], length.out = theta_points)
  curve_tbl <- build_curve_tables(curve_spec, theta_grid)
  prob_df <- as.data.frame(curve_tbl$probabilities, stringsAsFactors = FALSE)
  exp_df <- as.data.frame(curve_tbl$expected, stringsAsFactors = FALSE)

  mode_peaks <- prob_df |>
    dplyr::group_by(.data$CurveGroup, .data$Category) |>
    dplyr::arrange(dplyr::desc(.data$Probability), .by_group = TRUE) |>
    dplyr::slice_head(n = 1) |>
    dplyr::transmute(
      CurveGroup = .data$CurveGroup,
      Category = .data$Category,
      PeakTheta = .data$Theta,
      PeakProbability = .data$Probability
    ) |>
    dplyr::ungroup() |>
    as.data.frame(stringsAsFactors = FALSE)

  mode_boundaries <- prob_df |>
    dplyr::group_by(.data$CurveGroup, .data$Theta) |>
    dplyr::arrange(dplyr::desc(.data$Probability), .by_group = TRUE) |>
    dplyr::slice_head(n = 1) |>
    dplyr::transmute(
      CurveGroup = .data$CurveGroup,
      Theta = .data$Theta,
      ModalCategory = .data$Category
    ) |>
    dplyr::ungroup() |>
    dplyr::group_by(.data$CurveGroup) |>
    dplyr::arrange(.data$Theta, .by_group = TRUE) |>
    dplyr::mutate(
      PrevCategory = dplyr::lag(.data$ModalCategory),
      PrevTheta = dplyr::lag(.data$Theta)
    ) |>
    dplyr::filter(
      !is.na(.data$PrevCategory),
      .data$ModalCategory != .data$PrevCategory
    ) |>
    dplyr::transmute(
      CurveGroup = .data$CurveGroup,
      LowerCategory = .data$PrevCategory,
      UpperCategory = .data$ModalCategory,
      ModeBoundaryTheta = (.data$PrevTheta + .data$Theta) / 2
    ) |>
    dplyr::ungroup() |>
    as.data.frame(stringsAsFactors = FALSE)

  median_thresholds <- curve_spec$step_points |>
    dplyr::mutate(
      StepIndex = as.integer(.data$StepIndex),
      LowerCategory = as.character(curve_spec$categories[.data$StepIndex]),
      UpperCategory = as.character(curve_spec$categories[.data$StepIndex + 1]),
      MedianThreshold = .data$Threshold
    ) |>
    dplyr::select("CurveGroup", "Step", "StepIndex", "LowerCategory", "UpperCategory", "MedianThreshold") |>
    as.data.frame(stringsAsFactors = FALSE)

  cat_values <- suppressWarnings(as.numeric(curve_spec$categories))
  if (!all(is.finite(cat_values))) {
    cat_values <- seq_along(curve_spec$categories) - 1
  }
  mean_halfscore_points <- purrr::map_dfr(unique(exp_df$CurveGroup), function(g) {
    sub <- exp_df[exp_df$CurveGroup == g, , drop = FALSE]
    if (nrow(sub) == 0) return(tibble::tibble())
    mids <- (head(cat_values, -1) + tail(cat_values, -1)) / 2
    purrr::map_dfr(seq_along(mids), function(i) {
      tibble::tibble(
        CurveGroup = g,
        LowerCategory = as.character(curve_spec$categories[i]),
        UpperCategory = as.character(curve_spec$categories[i + 1]),
        MeanHalfScore = mids[i],
        MeanBoundaryTheta = closest_theta_for_target(sub$Theta, sub$ExpectedScore, mids[i])
      )
    })
  }) |>
    as.data.frame(stringsAsFactors = FALSE)

  out <- list(
    category_table = category_table,
    mode_peaks = mode_peaks,
    mode_boundaries = mode_boundaries,
    median_thresholds = median_thresholds,
    mean_halfscore_points = mean_halfscore_points,
    settings = list(
      theta_range = theta_range,
      theta_points = theta_points,
      drop_unused = isTRUE(drop_unused),
      include_fixed = isTRUE(include_fixed),
      fixed_max_rows = max(10L, as.integer(fixed_max_rows))
    )
  )
  if (isTRUE(include_fixed)) {
    fixed_max_rows <- max(10L, as.integer(fixed_max_rows))
    out$fixed <- build_sectioned_fixed_report(
      title = "FACETS-style Table 8 Bar-chart Export",
      sections = list(
        list(title = "Category table", data = category_table, max_rows = fixed_max_rows),
        list(title = "Mode peaks", data = mode_peaks, max_rows = fixed_max_rows),
        list(title = "Mode boundaries", data = mode_boundaries, max_rows = fixed_max_rows),
        list(title = "Median thresholds", data = median_thresholds, max_rows = fixed_max_rows),
        list(title = "Mean half-score transition points", data = mean_halfscore_points, max_rows = fixed_max_rows)
      ),
      max_col_width = 24
    )
  }
  out
}

#' Build a FACETS Graphfile-style Table 8 curves export
#'
#' @param fit Output from [fit_mfrm()].
#' @param theta_range Theta/logit range for curve coordinates.
#' @param theta_points Number of points on the theta grid.
#' @param digits Rounding digits for numeric graph output.
#' @param include_fixed If `TRUE`, include a FACETS-style fixed-width text block.
#' @param fixed_max_rows Maximum rows shown in the fixed-width graph table.
#'
#' @details
#' FACETS `Graphfile=` output for Table 8 contains curve coordinates:
#' scale number, measure, expected score, expected category, and
#' category probabilities by measure. This helper returns the same
#' information as data frames ready for CSV export.
#'
#' @section Lifecycle:
#' Soft-deprecated. Prefer [category_curves_report()].
#'
#' @return A named list with:
#' - `graphfile`: wide table with `Scale, Measure, Expected, ExpCat, Prob:*`
#' - `graphfile_syntactic`: same content with syntactic probability column names
#' - `probabilities`: long probability table (`Theta`, `Category`, `CurveGroup`)
#' - `expected_ogive`: long expected-score table
#' - `fixed`: fixed-width report text (when `include_fixed = TRUE`)
#' - `settings`: applied options
#'
#' @seealso [category_structure_report()], [rating_scale_table()], [plot.mfrm_fit()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' t8c <- category_curves_report(fit, theta_points = 101)
#' @keywords internal
#' @noRd
table8_curves_export <- function(fit,
                                 theta_range = c(-6, 6),
                                 theta_points = 241,
                                 digits = 4,
                                 include_fixed = FALSE,
                                 fixed_max_rows = 400) {
  signal_legacy_name_deprecation(
    old_name = "table8_curves_export",
    new_name = "category_curves_report",
    suppress_if_called_from = c("category_curves_report", "facets_output_file_bundle")
  )
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().")
  }
  theta_points <- max(51L, as.integer(theta_points))
  theta_range <- as.numeric(theta_range)
  if (length(theta_range) != 2 || !all(is.finite(theta_range)) || theta_range[1] >= theta_range[2]) {
    stop("`theta_range` must be a numeric length-2 vector with increasing values.")
  }
  digits <- max(0L, as.integer(digits))

  curve_spec <- build_step_curve_spec(fit)
  theta_grid <- seq(theta_range[1], theta_range[2], length.out = theta_points)
  curve_tbl <- build_curve_tables(curve_spec, theta_grid)
  prob_df <- as.data.frame(curve_tbl$probabilities, stringsAsFactors = FALSE)
  exp_df <- as.data.frame(curve_tbl$expected, stringsAsFactors = FALSE)

  groups <- unique(as.character(exp_df$CurveGroup))
  scale_map <- setNames(seq_along(groups), groups)

  prob_wide <- prob_df |>
    dplyr::mutate(
      Category = as.character(.data$Category),
      ProbCol = paste0("Prob:", .data$Category)
    ) |>
    dplyr::select("CurveGroup", "Theta", "ProbCol", "Probability") |>
    tidyr::pivot_wider(names_from = "ProbCol", values_from = "Probability")

  graph_tbl <- exp_df |>
    dplyr::left_join(prob_wide, by = c("CurveGroup", "Theta")) |>
    dplyr::mutate(
      Scale = as.integer(scale_map[as.character(.data$CurveGroup)]),
      Measure = .data$Theta
    )

  prob_cols <- grep("^Prob:", names(graph_tbl), value = TRUE)
  if (length(prob_cols) == 0) {
    stop("No probability columns were generated for Table 8 curves export.")
  }
  category_values <- suppressWarnings(as.numeric(sub("^Prob:", "", prob_cols)))
  if (!all(is.finite(category_values))) {
    category_values <- seq_along(prob_cols) - 1
  }
  expected_vec <- suppressWarnings(as.numeric(graph_tbl$ExpectedScore))
  nearest_idx <- vapply(expected_vec, function(ev) {
    if (!is.finite(ev)) return(NA_integer_)
    as.integer(which.min(abs(ev - category_values)))
  }, integer(1))
  exp_cat <- ifelse(is.na(nearest_idx), NA_real_, category_values[nearest_idx])

  graph_tbl <- graph_tbl |>
    dplyr::mutate(ExpCat = exp_cat) |>
    dplyr::select("Scale", "CurveGroup", "Measure", "ExpectedScore", "ExpCat", dplyr::all_of(prob_cols)) |>
    dplyr::arrange(.data$Scale, .data$Measure)

  graph_tbl_syntactic <- graph_tbl
  names(graph_tbl_syntactic) <- gsub("^Prob:", "Prob_", names(graph_tbl_syntactic))
  names(graph_tbl_syntactic) <- sub("^ExpectedScore$", "Expected", names(graph_tbl_syntactic))

  graph_tbl_facets <- graph_tbl
  names(graph_tbl_facets)[names(graph_tbl_facets) == "ExpectedScore"] <- "Expected"
  graph_tbl_facets <- as.data.frame(graph_tbl_facets, check.names = FALSE, stringsAsFactors = FALSE)

  round_numeric <- function(df) {
    if (!is.data.frame(df) || nrow(df) == 0) return(df)
    out <- df
    num_cols <- vapply(out, is.numeric, logical(1))
    out[num_cols] <- lapply(out[num_cols], round, digits = digits)
    out
  }

  out <- list(
    graphfile = round_numeric(graph_tbl_facets),
    graphfile_syntactic = round_numeric(graph_tbl_syntactic),
    probabilities = round_numeric(as.data.frame(prob_df, stringsAsFactors = FALSE)),
    expected_ogive = round_numeric(as.data.frame(exp_df, stringsAsFactors = FALSE)),
    settings = list(
      theta_range = theta_range,
      theta_points = theta_points,
      digits = digits,
      include_fixed = isTRUE(include_fixed),
      fixed_max_rows = max(25L, as.integer(fixed_max_rows)),
      scales = as.data.frame(
        data.frame(
          Scale = as.integer(scale_map),
          CurveGroup = names(scale_map),
          stringsAsFactors = FALSE
        )
      )
    )
  )
  if (isTRUE(include_fixed)) {
    fixed_max_rows <- max(25L, as.integer(fixed_max_rows))
    out$fixed <- build_sectioned_fixed_report(
      title = "FACETS-style Table 8 Curves Graphfile",
      sections = list(
        list(title = "Graphfile wide table", data = out$graphfile, max_rows = fixed_max_rows),
        list(title = "Scale map", data = out$settings$scales),
        list(title = "Expected ogive", data = out$expected_ogive, max_rows = fixed_max_rows)
      ),
      max_col_width = 24
    )
  }
  out
}

#' Build FACETS-style output-file bundle (`GRAPH=` / `SCORE=`)
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()] (used for score file).
#' @param include Output components to include: `"graph"` and/or `"score"`.
#' @param theta_range Theta/logit range for graph coordinates.
#' @param theta_points Number of points on the theta grid for graph coordinates.
#' @param digits Rounding digits for numeric fields.
#' @param include_fixed If `TRUE`, include fixed-width text mirrors of output tables.
#' @param fixed_max_rows Maximum rows shown in fixed-width text blocks.
#' @param write_files If `TRUE`, write selected outputs to files in `output_dir`.
#' @param output_dir Output directory used when `write_files = TRUE`.
#' @param file_prefix Prefix used for output file names.
#' @param overwrite If `FALSE`, existing output files are not overwritten.
#'
#' @details
#' FACETS output files often include:
#' - graph coordinates for Table 8 curves (`GRAPH=` / `Graphfile=`), and
#' - observation-level modeled score lines (`SCORE=`-style inspection).
#'
#' This helper returns both as data frames and can optionally write
#' CSV/fixed-width text files to disk.
#'
#' `summary(out)` is supported through `summary()`.
#' `plot(out)` is dispatched through `plot()` for class
#' `mfrm_output_bundle` (`type = "graph_expected"`, `"score_residuals"`,
#' `"obs_probability"`).
#'
#' @section Interpreting output:
#' - `graphfile`: FACETS-style wide curve coordinates (human-readable labels).
#' - `graphfile_syntactic`: same curves with syntactic column names for programmatic use.
#' - `scorefile`: observation-level observed/expected/residual diagnostics.
#' - `written_files`: audit trail of files produced when `write_files = TRUE`.
#'
#' For reproducible pipelines, prefer `graphfile_syntactic` and keep
#' `written_files` in run logs.
#'
#' @section Typical workflow:
#' 1. Fit and diagnose model.
#' 2. Generate bundle with `include = c("graph", "score")`.
#' 3. Validate with `summary(out)` / `plot(out)`.
#' 4. Export with `write_files = TRUE` for reporting handoff.
#'
#' @return A named list including:
#' - `graphfile` / `graphfile_syntactic` when `"graph"` is requested
#' - `scorefile` when `"score"` is requested
#' - `graphfile_fixed` / `scorefile_fixed` when `include_fixed = TRUE`
#' - `written_files` when `write_files = TRUE`
#' - `settings`: applied options
#'
#' @seealso [category_curves_report()], [diagnose_mfrm()], [unexpected_response_table()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' out <- facets_output_file_bundle(fit, diagnostics = diagnose_mfrm(fit, residual_pca = "none"))
#' summary(out)
#' p_out <- plot(out, draw = FALSE)
#' class(p_out)
#' @export
facets_output_file_bundle <- function(fit,
                                      diagnostics = NULL,
                                      include = c("graph", "score"),
                                      theta_range = c(-6, 6),
                                      theta_points = 241,
                                      digits = 4,
                                      include_fixed = FALSE,
                                      fixed_max_rows = 400,
                                      write_files = FALSE,
                                      output_dir = NULL,
                                      file_prefix = "mfrmr_output",
                                      overwrite = FALSE) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().")
  }
  include <- unique(tolower(as.character(include)))
  allowed <- c("graph", "score")
  bad <- setdiff(include, allowed)
  if (length(bad) > 0) {
    stop("Unsupported `include` values: ", paste(bad, collapse = ", "), ". Use: ", paste(allowed, collapse = ", "))
  }
  if (length(include) == 0) {
    stop("`include` must contain at least one of: graph, score.")
  }
  digits <- max(0L, as.integer(digits))
  include_fixed <- isTRUE(include_fixed)
  fixed_max_rows <- max(25L, as.integer(fixed_max_rows))
  write_files <- isTRUE(write_files)
  overwrite <- isTRUE(overwrite)
  file_prefix <- as.character(file_prefix[1] %||% "mfrmr_output")
  if (!nzchar(file_prefix)) file_prefix <- "mfrmr_output"
  output_dir <- if (is.null(output_dir)) NULL else as.character(output_dir[1])

  if (write_files) {
    if (is.null(output_dir) || !nzchar(output_dir)) {
      stop("`output_dir` must be supplied when `write_files = TRUE`.")
    }
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }
    if (!dir.exists(output_dir)) {
      stop("Failed to create `output_dir`: ", output_dir)
    }
  }

  written_files <- data.frame(
    Component = character(0),
    Format = character(0),
    Path = character(0),
    stringsAsFactors = FALSE
  )
  add_written <- function(component, format, path) {
    written_files <<- rbind(
      written_files,
      data.frame(
        Component = as.character(component),
        Format = as.character(format),
        Path = as.character(path),
        stringsAsFactors = FALSE
      )
    )
    invisible(NULL)
  }
  write_csv_if_needed <- function(df, filename, component) {
    if (!write_files) return(invisible(NULL))
    path <- file.path(output_dir, filename)
    if (file.exists(path) && !overwrite) {
      stop("Output file already exists: ", path, ". Set `overwrite = TRUE` to replace it.")
    }
    utils::write.csv(df, file = path, row.names = FALSE, na = "")
    add_written(component = component, format = "csv", path = path)
    invisible(NULL)
  }
  write_txt_if_needed <- function(text, filename, component) {
    if (!write_files) return(invisible(NULL))
    path <- file.path(output_dir, filename)
    if (file.exists(path) && !overwrite) {
      stop("Output file already exists: ", path, ". Set `overwrite = TRUE` to replace it.")
    }
    writeLines(as.character(text), con = path, useBytes = TRUE)
    add_written(component = component, format = "txt", path = path)
    invisible(NULL)
  }

  out <- list()
  if ("graph" %in% include) {
    graph <- with_legacy_name_warning_suppressed(
      table8_curves_export(
        fit = fit,
        theta_range = theta_range,
        theta_points = theta_points,
        digits = digits,
        include_fixed = include_fixed,
        fixed_max_rows = fixed_max_rows
      )
    )
    out$graphfile <- graph$graphfile
    out$graphfile_syntactic <- graph$graphfile_syntactic
    if (include_fixed) {
      out$graphfile_fixed <- graph$fixed
    }
    write_csv_if_needed(out$graphfile, paste0(file_prefix, "_graphfile.csv"), "graphfile")
    write_csv_if_needed(out$graphfile_syntactic, paste0(file_prefix, "_graphfile_syntactic.csv"), "graphfile_syntactic")
    if (include_fixed && !is.null(out$graphfile_fixed)) {
      write_txt_if_needed(out$graphfile_fixed, paste0(file_prefix, "_graphfile_fixed.txt"), "graphfile_fixed")
    }
  }

  if ("score" %in% include) {
    if (is.null(diagnostics)) {
      diagnostics <- diagnose_mfrm(fit, residual_pca = "none")
    }
    if (is.null(diagnostics$obs) || nrow(diagnostics$obs) == 0) {
      stop("`diagnostics$obs` is empty. Run diagnose_mfrm() first.")
    }

    obs <- as.data.frame(diagnostics$obs, stringsAsFactors = FALSE)
    keep_cols <- c(
      "Person", fit$config$facet_names, "Observed", "Expected", "Residual",
      "StdResidual", "Var", "Weight", "score_k", "PersonMeasure"
    )
    keep_cols <- keep_cols[keep_cols %in% names(obs)]
    scorefile <- obs[, keep_cols, drop = FALSE]

    probs <- compute_prob_matrix(fit)
    if (!is.null(probs) && nrow(probs) == nrow(obs) && "score_k" %in% names(obs)) {
      k_idx <- suppressWarnings(as.integer(obs$score_k)) + 1L
      obs_prob <- vapply(seq_len(nrow(obs)), function(i) {
        k <- k_idx[i]
        if (is.finite(k) && k >= 1L && k <= ncol(probs)) {
          return(as.numeric(probs[i, k]))
        }
        NA_real_
      }, numeric(1))
      scorefile$ObsProb <- obs_prob
    }
    out$scorefile <- round_numeric_df(scorefile, digits = digits)
    if (include_fixed) {
      out$scorefile_fixed <- build_sectioned_fixed_report(
        title = "FACETS-style SCORE Output",
        sections = list(
          list(title = "Score file", data = out$scorefile, max_rows = fixed_max_rows)
        ),
        max_col_width = 24
      )
    }
    write_csv_if_needed(out$scorefile, paste0(file_prefix, "_scorefile.csv"), "scorefile")
    if (include_fixed && !is.null(out$scorefile_fixed)) {
      write_txt_if_needed(out$scorefile_fixed, paste0(file_prefix, "_scorefile_fixed.txt"), "scorefile_fixed")
    }
  }

  if (write_files) {
    out$written_files <- written_files
  }

  out$settings <- list(
    include = include,
    theta_range = as.numeric(theta_range),
    theta_points = as.integer(theta_points),
    digits = digits,
    include_fixed = include_fixed,
    fixed_max_rows = fixed_max_rows,
    write_files = write_files,
    output_dir = output_dir,
    file_prefix = file_prefix,
    overwrite = overwrite
  )
  as_mfrm_bundle(out, "mfrm_output_bundle")
}

extract_pca_eigenvalues <- function(pca_bundle) {
  if (is.null(pca_bundle)) return(numeric(0))

  eig <- numeric(0)
  if (!is.null(pca_bundle$pca) && "values" %in% names(pca_bundle$pca)) {
    eig <- suppressWarnings(as.numeric(pca_bundle$pca$values))
  }
  if (length(eig) == 0 && !is.null(pca_bundle$cor_matrix)) {
    eig <- tryCatch(
      suppressWarnings(as.numeric(eigen(pca_bundle$cor_matrix, symmetric = TRUE, only.values = TRUE)$values)),
      error = function(e) numeric(0)
    )
  }

  eig[is.finite(eig)]
}

build_pca_variance_table <- function(pca_bundle, facet = NULL) {
  eig <- extract_pca_eigenvalues(pca_bundle)
  if (length(eig) == 0) return(data.frame())

  total <- sum(eig, na.rm = TRUE)
  prop <- if (is.finite(total) && total > 0) eig / total else rep(NA_real_, length(eig))
  out <- data.frame(
    Component = seq_along(eig),
    Eigenvalue = eig,
    Proportion = prop,
    Cumulative = cumsum(prop),
    stringsAsFactors = FALSE
  )
  if (!is.null(facet)) out$Facet <- facet
  out
}

infer_facet_names <- function(diagnostics) {
  if (!is.null(diagnostics$facet_names) && length(diagnostics$facet_names) > 0) {
    return(unique(as.character(diagnostics$facet_names)))
  }

  if (!is.null(diagnostics$measures) && "Facet" %in% names(diagnostics$measures)) {
    f <- unique(as.character(diagnostics$measures$Facet))
    f <- setdiff(f, "Person")
    if (length(f) > 0) return(f)
  }

  if (!is.null(diagnostics$obs)) {
    nm <- names(diagnostics$obs)
    skip <- c(
      "Person", "Score", "Weight", "score_k", "PersonMeasure",
      "Observed", "Expected", "Var", "Residual", "StdResidual", "StdSq"
    )
    f <- setdiff(nm, skip)
    if (length(f) > 0) return(f)
  }

  character(0)
}

#' Run residual PCA for unidimensionality checks
#'
#' FACETS-style residual diagnostics can be inspected in two ways:
#' 1) overall residual PCA on the person x combined-facet matrix
#' 2) facet-specific residual PCA on person x facet-level matrices
#'
#' @param diagnostics Output from [diagnose_mfrm()] or [fit_mfrm()].
#' @param mode `"overall"`, `"facet"`, or `"both"`.
#' @param facets Optional subset of facets for facet-specific PCA.
#' @param pca_max_factors Maximum number of retained components.
#'
#' @details
#' The function works on standardized residual structures derived from
#' [diagnose_mfrm()]. When a fitted object from [fit_mfrm()] is supplied,
#' diagnostics are computed internally.
#'
#' Output tables use:
#' - `Component`: principal-component index (1, 2, ...)
#' - `Eigenvalue`: eigenvalue for each component
#' - `Proportion`: component variance proportion
#' - `Cumulative`: cumulative variance proportion
#'
#' For `mode = "facet"` or `"both"`, `by_facet_table` additionally includes
#' a `Facet` column.
#'
#' `summary(pca)` is supported through `summary()`.
#' `plot(pca)` is dispatched through `plot()` for class
#' `mfrm_residual_pca`. Available types include `"overall_scree"`,
#' `"facet_scree"`, `"overall_loadings"`, and `"facet_loadings"`.
#'
#' @section Interpreting output:
#' Use `overall_table` first:
#' - high PC1 eigenvalue and high `Proportion` suggest stronger residual
#'   structure beyond the main Rasch dimension.
#'
#' Then inspect `by_facet_table`:
#' - helps localize which facet contributes most to residual structure.
#'
#' Finally, inspect loadings via [plot_residual_pca()] to identify which
#' variables/elements drive each component.
#'
#' @section Typical workflow:
#' 1. Fit model and run [diagnose_mfrm()] with `residual_pca = "none"` or `"both"`.
#' 2. Call `analyze_residual_pca(..., mode = "both")`.
#' 3. Review `summary(pca)`, then plot scree/loadings.
#' 4. Cross-check with fit/misfit diagnostics before conclusions.
#'
#' @return
#' A named list with:
#' - `mode`: resolved mode used for computation
#' - `facet_names`: facets analyzed
#' - `overall`: overall PCA bundle (or `NULL`)
#' - `by_facet`: named list of facet PCA bundles
#' - `overall_table`: variance table for overall PCA
#' - `by_facet_table`: stacked variance table across facets
#'
#' @seealso [diagnose_mfrm()], [plot_residual_pca()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' diag <- diagnose_mfrm(fit, residual_pca = "both")
#' pca <- analyze_residual_pca(diag, mode = "both")
#' pca2 <- analyze_residual_pca(fit, mode = "both")
#' summary(pca)
#' p <- plot_residual_pca(pca, mode = "overall", plot_type = "scree", draw = FALSE)
#' class(p)
#' head(p$data)
#' head(pca$overall_table)
#' @export
analyze_residual_pca <- function(diagnostics,
                                 mode = c("overall", "facet", "both"),
                                 facets = NULL,
                                 pca_max_factors = 10L) {
  mode <- match.arg(tolower(mode), c("overall", "facet", "both"))

  if (inherits(diagnostics, "mfrm_fit")) {
    diagnostics <- diagnose_mfrm(
      diagnostics,
      residual_pca = "none",
      pca_max_factors = pca_max_factors
    )
  }

  if (!is.list(diagnostics)) {
    stop("`diagnostics` must be output from diagnose_mfrm() or fit_mfrm().")
  }

  if (is.null(diagnostics$obs) || nrow(diagnostics$obs) == 0) {
    stop("diagnostics$obs is empty. Run diagnose_mfrm() first.")
  }

  facet_names <- infer_facet_names(diagnostics)
  if (!is.null(facets)) {
    facet_names <- intersect(facet_names, as.character(facets))
    if (length(facet_names) == 0) {
      stop("No matching facets found in diagnostics for `facets`.")
    }
  }

  out_overall <- NULL
  out_by_facet <- list()

  if (mode %in% c("overall", "both")) {
    can_reuse <- is.null(facets) && !is.null(diagnostics$residual_pca_overall)
    out_overall <- if (can_reuse) {
      diagnostics$residual_pca_overall
    } else {
      compute_pca_overall(
        obs_df = diagnostics$obs,
        facet_names = facet_names,
        max_factors = pca_max_factors
      )
    }
  }

  if (mode %in% c("facet", "both")) {
    can_reuse <- is.null(facets) && !is.null(diagnostics$residual_pca_by_facet)
    out_by_facet <- if (can_reuse) {
      diagnostics$residual_pca_by_facet
    } else {
      compute_pca_by_facet(
        obs_df = diagnostics$obs,
        facet_names = facet_names,
        max_factors = pca_max_factors
      )
    }

    if (!is.null(facets) && length(out_by_facet) > 0) {
      out_by_facet <- out_by_facet[intersect(names(out_by_facet), facet_names)]
    }
  }

  overall_table <- build_pca_variance_table(out_overall)
  by_facet_table <- if (length(out_by_facet) == 0) {
    data.frame()
  } else {
    tbls <- lapply(names(out_by_facet), function(f) {
      build_pca_variance_table(out_by_facet[[f]], facet = f)
    })
    tbls <- tbls[vapply(tbls, nrow, integer(1)) > 0]
    if (length(tbls) == 0) data.frame() else dplyr::bind_rows(tbls)
  }

  out <- list(
    mode = mode,
    facet_names = facet_names,
    overall = out_overall,
    by_facet = out_by_facet,
    overall_table = overall_table,
    by_facet_table = by_facet_table
  )
  as_mfrm_bundle(out, "mfrm_residual_pca")
}

resolve_pca_input <- function(x) {
  if (is.null(x)) stop("Input cannot be NULL.")
  if (!is.null(x$overall_table) || !is.null(x$by_facet_table)) return(x)
  if (inherits(x, "mfrm_fit")) return(analyze_residual_pca(x, mode = "both"))
  if (!is.null(x$obs)) return(analyze_residual_pca(x, mode = "both"))
  stop("Input must be fit from fit_mfrm(), diagnostics from diagnose_mfrm(), or result from analyze_residual_pca().")
}

extract_loading_table <- function(pca_bundle, component = 1L, top_n = 20L) {
  if (is.null(pca_bundle) || is.null(pca_bundle$pca) || is.null(pca_bundle$pca$loadings)) {
    return(data.frame())
  }

  loads <- tryCatch(as.matrix(unclass(pca_bundle$pca$loadings)), error = function(e) NULL)
  if (is.null(loads) || nrow(loads) == 0) return(data.frame())
  if (component > ncol(loads) || component < 1) return(data.frame())

  vals <- suppressWarnings(as.numeric(loads[, component]))
  vars <- rownames(loads)
  if (is.null(vars)) vars <- paste0("V", seq_along(vals))

  ok <- is.finite(vals)
  if (!any(ok)) return(data.frame())

  tbl <- data.frame(
    Variable = vars[ok],
    Loading = vals[ok],
    stringsAsFactors = FALSE
  )
  tbl <- tbl[order(abs(tbl$Loading), decreasing = TRUE), , drop = FALSE]
  top_n <- max(1L, as.integer(top_n))
  head(tbl, n = min(nrow(tbl), top_n))
}

#' Visualize residual PCA results
#'
#' @param x Output from [analyze_residual_pca()], [diagnose_mfrm()], or [fit_mfrm()].
#' @param mode `"overall"` or `"facet"`.
#' @param facet Facet name for `mode = "facet"`.
#' @param plot_type `"scree"` or `"loadings"`.
#' @param component Component index for loadings plot.
#' @param top_n Maximum number of variables shown in loadings plot.
#' @param draw If `TRUE`, draws the plot using base graphics.
#'
#' @details
#' `x` can be either:
#' - output of [analyze_residual_pca()], or
#' - a diagnostics object from [diagnose_mfrm()] (PCA is computed internally), or
#' - a fitted object from [fit_mfrm()] (diagnostics and PCA are computed internally).
#'
#' Plot types:
#' - `"scree"`: component vs eigenvalue line plot
#' - `"loadings"`: horizontal bar chart of top absolute loadings
#'
#' For `mode = "facet"` and `facet = NULL`, the first available facet is used.
#'
#' @section Interpreting output:
#' - `plot_type = "scree"`: look for dominant early components
#'   (eigenvalue notably above 1 as a practical warning sign).
#' - `plot_type = "loadings"`: identifies variables/elements driving each
#'   component; inspect both sign and absolute magnitude.
#'
#' Facet mode (`mode = "facet"`) helps localize residual structure to a
#' specific facet after global PCA review.
#'
#' @section Typical workflow:
#' 1. Run [diagnose_mfrm()] with `residual_pca = "overall"` or `"both"`.
#' 2. Build PCA object via [analyze_residual_pca()] (or pass diagnostics directly).
#' 3. Use scree plot first, then loadings plot for targeted interpretation.
#'
#' @return
#' A named list of plotting data (class `mfrm_plot_data`) with:
#' - `plot`: `"scree"` or `"loadings"`
#' - `mode`: `"overall"` or `"facet"`
#' - `facet`: facet name (or `NULL`)
#' - `title`: plot title text
#' - `data`: underlying table used for plotting
#'
#' @seealso [analyze_residual_pca()], [diagnose_mfrm()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' diag <- diagnose_mfrm(fit, residual_pca = "both")
#' pca <- analyze_residual_pca(diag, mode = "both")
#' plt <- plot_residual_pca(pca, mode = "overall", plot_type = "scree", draw = FALSE)
#' head(plt$data)
#' plt_load <- plot_residual_pca(
#'   pca, mode = "overall", plot_type = "loadings", component = 1, draw = FALSE
#' )
#' head(plt_load$data)
#' @export
plot_residual_pca <- function(x,
                              mode = c("overall", "facet"),
                              facet = NULL,
                              plot_type = c("scree", "loadings"),
                              component = 1L,
                              top_n = 20L,
                              draw = TRUE) {
  mode <- match.arg(tolower(mode), c("overall", "facet"))
  plot_type <- match.arg(tolower(plot_type), c("scree", "loadings"))
  pca_obj <- resolve_pca_input(x)

  pca_bundle <- NULL
  title_suffix <- ""

  if (mode == "overall") {
    pca_bundle <- pca_obj$overall
    title_suffix <- "Overall Residual PCA"
  } else {
    if (is.null(pca_obj$by_facet) || length(pca_obj$by_facet) == 0) {
      stop("No facet-level PCA results available.")
    }
    if (is.null(facet)) facet <- names(pca_obj$by_facet)[1]
    if (!facet %in% names(pca_obj$by_facet)) {
      stop("Requested facet not found in PCA results.")
    }
    pca_bundle <- pca_obj$by_facet[[facet]]
    title_suffix <- paste0("Residual PCA - ", facet)
  }

  if (plot_type == "scree") {
    tbl <- build_pca_variance_table(pca_bundle)
    if (nrow(tbl) == 0) stop("No eigenvalues available for scree plot.")
    title <- paste0(title_suffix, " (Scree)")
    if (isTRUE(draw)) {
      graphics::plot(
        x = tbl$Component,
        y = tbl$Eigenvalue,
        type = "b",
        pch = 16,
        col = "black",
        xlab = "Component",
        ylab = "Eigenvalue",
        main = title
      )
      graphics::abline(h = 1, lty = 2, col = "gray40")
    }

    out <- list(
      plot = "scree",
      mode = mode,
      facet = if (mode == "facet") facet else NULL,
      title = title,
      data = tbl
    )
    class(out) <- c("mfrm_plot_data", class(out))
    return(invisible(out))
  }

  load_tbl <- extract_loading_table(
    pca_bundle = pca_bundle,
    component = as.integer(component),
    top_n = as.integer(top_n)
  )
  if (nrow(load_tbl) == 0) stop("No loadings available for the requested component.")

  load_tbl <- load_tbl[order(load_tbl$Loading), , drop = FALSE]
  title <- paste0(title_suffix, " (Loadings: PC", as.integer(component), ")")

  if (isTRUE(draw)) {
    cols <- ifelse(load_tbl$Loading >= 0, "#1b9e77", "#d95f02")
    graphics::barplot(
      height = load_tbl$Loading,
      names.arg = load_tbl$Variable,
      horiz = TRUE,
      las = 1,
      col = cols,
      xlab = "Loading",
      main = title
    )
    graphics::abline(v = 0, lty = 2, col = "gray40")
  }

  out <- list(
    plot = "loadings",
    mode = mode,
    facet = if (mode == "facet") facet else NULL,
    title = title,
    component = as.integer(component),
    data = load_tbl
  )
  class(out) <- c("mfrm_plot_data", class(out))
  invisible(out)
}

#' Estimate FACETS-style bias/interaction terms iteratively
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Output from [diagnose_mfrm()].
#' @param facet_a First facet name.
#' @param facet_b Second facet name.
#' @param interaction_facets Character vector of two or more facets to model as
#'   one interaction effect. When supplied, this takes precedence over
#'   `facet_a`/`facet_b`.
#' @param max_abs Bound for absolute bias size.
#' @param omit_extreme Omit extreme-only elements.
#' @param max_iter Iteration cap.
#' @param tol Convergence tolerance.
#'
#' @details
#' The function estimates interaction contrasts with iterative recalibration in
#' a FACETS-like style.
#'
#' - For two-way mode, use `facet_a` and `facet_b` (or `interaction_facets`
#'   with length 2).
#' - For higher-order mode, provide `interaction_facets` with length >= 3.
#'
#' Typical usage:
#' 1. fit model via [fit_mfrm()]
#' 2. compute diagnostics via [diagnose_mfrm()]
#' 3. call `estimate_bias()` for one 2-way or higher-order interaction
#' 4. format output with [build_fixed_reports()]
#'
#' @section Interpreting output:
#' Use `summary` for global magnitude, then inspect `table` for cell-level
#' interaction effects.
#'
#' Prioritize rows with:
#' - larger `|Bias Size|`
#' - larger `|t|`
#' - smaller `Prob.`
#'
#' `iteration` helps verify whether iterative recalibration stabilized.
#'
#' @section Typical workflow:
#' 1. Fit and diagnose model.
#' 2. Run `estimate_bias(...)` for target interaction facets.
#' 3. Review `summary(bias)` and `bias$table`.
#' 4. Visualize/report via [plot_bias_interaction()] and [build_fixed_reports()].
#'
#' @section Interpreting key output columns:
#' In `bias$table`, the most-used columns are:
#' - `Bias Size`: estimated interaction effect (logit scale)
#' - `t` and `Prob.`: significance-style screening metrics
#' - `Obs-Exp Average`: direction and practical size of observed-vs-expected gap
#'
#' `bias$iteration` records recalibration trajectory and can be checked to assess
#' whether additional iterations may be needed.
#'
#' @return
#' An object of class `mfrm_bias` with:
#' - `table`: interaction rows with effect size, SE, t, p, fit columns
#' - `summary`: compact summary statistics
#' - `chi_sq`: fixed-effect chi-square style summary
#' - `facet_a`, `facet_b`: first two analyzed facet names (legacy compatibility)
#' - `interaction_facets`, `interaction_order`, `interaction_mode`: full
#'   interaction metadata
#' - `iteration`: iteration history/metadata
#'
#' @seealso [build_fixed_reports()], [build_apa_outputs()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#' bias <- estimate_bias(fit, diag, facet_a = "Rater", facet_b = "Criterion", max_iter = 2)
#' summary(bias)
#' p_bias <- plot_bias_interaction(bias, draw = FALSE)
#' class(p_bias)
#' @export
estimate_bias <- function(fit,
                          diagnostics,
                          facet_a = NULL,
                          facet_b = NULL,
                          interaction_facets = NULL,
                          max_abs = 10,
                          omit_extreme = TRUE,
                          max_iter = 4,
                          tol = 1e-3) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm(). ",
         "Got: ", paste(class(fit), collapse = "/"), ".", call. = FALSE)
  }
  if (!is.list(diagnostics) || is.null(diagnostics$obs)) {
    stop("`diagnostics` must be the output of diagnose_mfrm(). ",
         "Run: diagnostics <- diagnose_mfrm(fit)", call. = FALSE)
  }
  if (is.null(interaction_facets)) {
    if (is.null(facet_a) || is.null(facet_b)) {
      stop("Provide either `interaction_facets` (length >= 2) ",
           "or both `facet_a` and `facet_b`.", call. = FALSE)
    }
    interaction_facets <- c(as.character(facet_a[1]), as.character(facet_b[1]))
  }
  interaction_facets <- as.character(interaction_facets)
  interaction_facets <- interaction_facets[!is.na(interaction_facets) & nzchar(interaction_facets)]
  interaction_facets <- unique(interaction_facets)
  if (length(interaction_facets) < 2) {
    stop("`interaction_facets` must contain at least two facet names.")
  }
  if (is.null(facet_a) && length(interaction_facets) >= 1) facet_a <- interaction_facets[1]
  if (is.null(facet_b) && length(interaction_facets) >= 2) facet_b <- interaction_facets[2]

  out <- estimate_bias_interaction(
    res = fit,
    diagnostics = diagnostics,
    facet_a = facet_a,
    facet_b = facet_b,
    interaction_facets = interaction_facets,
    max_abs = max_abs,
    omit_extreme = omit_extreme,
    max_iter = max_iter,
    tol = tol
  )
  if (is.list(out) && length(out) > 0) {
    class(out) <- c("mfrm_bias", class(out))
  }
  out
}

#' Build FACETS-style fixed-width text reports
#'
#' @param bias_results Output from [estimate_bias()].
#' @param target_facet Optional target facet for pairwise contrast table.
#' @param branch Output branch:
#'   `"facets"` keeps FACETS-like fixed-width text;
#'   `"original"` returns compact sectioned fixed-width text for internal reporting.
#'
#' @details
#' This function generates plain-text, fixed-width output intended to be read in
#' console/log environments or exported into text reports.
#'
#' The pairwise section (Table 14 style) is only generated for 2-way bias runs.
#' For higher-order interactions (`interaction_facets` length >= 3), the function
#' returns the bias table text and a note explaining why pairwise contrasts were
#' skipped.
#'
#' @section Interpreting output:
#' - `bias_fixed`: fixed-width table of interaction effects.
#' - `pairwise_fixed`: pairwise contrast text (2-way only).
#' - `pairwise_table`: machine-readable contrast table.
#' - `interaction_label`: facets used for the bias run.
#'
#' @section Typical workflow:
#' 1. Run [estimate_bias()].
#' 2. Build text bundle with `build_fixed_reports(...)`.
#' 3. Use `summary()`/`plot()` for quick checks, then export text blocks.
#'
#' @return
#' A named list with:
#' - `bias_fixed`: fixed-width interaction table text
#' - `pairwise_fixed`: fixed-width pairwise contrast text
#' - `pairwise_table`: underlying pairwise data.frame
#'
#' @seealso [estimate_bias()], [build_apa_outputs()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#' bias <- estimate_bias(fit, diag, facet_a = "Rater", facet_b = "Criterion", max_iter = 2)
#' fixed <- build_fixed_reports(bias)
#' fixed_original <- build_fixed_reports(bias, branch = "original")
#' summary(fixed)
#' p <- plot(fixed, draw = FALSE)
#' p2 <- plot(fixed, type = "pvalue", draw = FALSE)
#' if (interactive()) {
#'   plot(
#'     fixed,
#'     type = "contrast",
#'     draw = TRUE,
#'     main = "Pairwise Contrasts (Customized)",
#'     palette = c(pos = "#1b9e77", neg = "#d95f02"),
#'     label_angle = 45
#'   )
#' }
#' @export
build_fixed_reports <- function(bias_results,
                                target_facet = NULL,
                                branch = c("facets", "original")) {
  branch <- match.arg(tolower(as.character(branch[1])), c("facets", "original"))
  style <- ifelse(branch == "facets", "facets_manual", "original")

  make_empty_bundle <- function(msg) {
    out <- list(
      bias_fixed = msg,
      pairwise_fixed = "No pairwise data",
      pairwise_table = tibble::tibble(),
      branch = branch,
      style = style,
      interaction_label = NA_character_,
      target_facet = as.character(target_facet %||% NA_character_)
    )
    out <- as_mfrm_bundle(out, "mfrm_fixed_reports")
    class(out) <- unique(c(paste0("mfrm_fixed_reports_", branch), class(out)))
    out
  }

  if (is.null(bias_results) || is.null(bias_results$table) || nrow(bias_results$table) == 0) {
    return(make_empty_bundle("No bias data"))
  }

  spec <- extract_bias_facet_spec(bias_results)
  if (is.null(spec) || length(spec$facets) < 2) {
    return(make_empty_bundle("No bias data"))
  }

  facets <- spec$facets
  interaction_label <- paste(facets, collapse = " x ")
  tbl <- as.data.frame(bias_results$table, stringsAsFactors = FALSE)

  core_cols <- c(
    "Sq", "Observd Score", "Expctd Score", "Observd Count", "Obs-Exp Average",
    "Bias Size", "S.E.", "t", "d.f.", "Prob.", "Infit", "Outfit"
  )
  detail_cols <- c(spec$index_cols, spec$level_cols, spec$measure_cols)
  keep_cols <- c(core_cols, detail_cols)
  keep_cols <- keep_cols[keep_cols %in% names(tbl)]
  tbl_display <- tbl |>
    dplyr::select(dplyr::all_of(keep_cols))

  if ("S.E." %in% names(tbl_display)) {
    tbl_display <- dplyr::rename(tbl_display, `Model S.E.` = `S.E.`)
  }
  if ("Infit" %in% names(tbl_display)) {
    tbl_display <- dplyr::rename(tbl_display, `Infit MnSq` = Infit)
  }
  if ("Outfit" %in% names(tbl_display)) {
    tbl_display <- dplyr::rename(tbl_display, `Outfit MnSq` = Outfit)
  }

  for (i in seq_along(facets)) {
    facet_i <- facets[i]
    idx_col <- spec$index_cols[i]
    lvl_col <- spec$level_cols[i]
    meas_col <- spec$measure_cols[i]
    if (idx_col %in% names(tbl_display)) {
      names(tbl_display)[names(tbl_display) == idx_col] <- paste0(facet_i, " N")
    }
    if (lvl_col %in% names(tbl_display)) {
      names(tbl_display)[names(tbl_display) == lvl_col] <- facet_i
    }
    if (meas_col %in% names(tbl_display)) {
      names(tbl_display)[names(tbl_display) == meas_col] <- paste0(facet_i, " measr")
    }
  }

  bias_cols <- names(tbl_display)
  bias_formats <- list(
    Sq = "{}",
    `Observd Score` = "{:.2f}",
    `Expctd Score` = "{:.2f}",
    `Observd Count` = "{:.0f}",
    `Obs-Exp Average` = "{:.2f}",
    `Bias Size` = "{:.2f}",
    `Model S.E.` = "{:.2f}",
    t = "{:.2f}",
    `d.f.` = "{:.0f}",
    `Prob.` = "{:.4f}",
    `Infit MnSq` = "{:.2f}",
    `Outfit MnSq` = "{:.2f}"
  )
  for (facet_i in facets) {
    bias_formats[[paste0(facet_i, " N")]] <- "{:.0f}"
    bias_formats[[paste0(facet_i, " measr")]] <- "{:.2f}"
  }

  if (branch == "facets") {
    bias_fixed <- build_bias_fixed_text(
      table_df = tbl_display,
      summary_df = bias_results$summary,
      chi_df = bias_results$chi_sq,
      facet_a = facets[1],
      facet_b = if (length(facets) >= 2) facets[2] else "",
      interaction_label = interaction_label,
      columns = bias_cols,
      formats = bias_formats
    )
  } else {
    bias_fixed <- build_sectioned_fixed_report(
      title = paste0("Bias interaction summary: ", interaction_label),
      sections = list(
        list(
          title = "Top interaction rows",
          data = tbl_display,
          columns = bias_cols,
          formats = bias_formats,
          max_rows = 40L
        ),
        list(
          title = "Summary",
          data = as.data.frame(bias_results$summary %||% data.frame(), stringsAsFactors = FALSE)
        ),
        list(
          title = "Chi-square",
          data = as.data.frame(bias_results$chi_sq %||% data.frame(), stringsAsFactors = FALSE)
        )
      ),
      max_col_width = 18,
      min_col_width = 6
    )
  }

  if (length(facets) != 2) {
    pairwise_tbl <- tibble::tibble()
    pairwise_fixed <- paste0(
      "Pairwise Table 14-style contrasts are available only for 2-way interactions.\n",
      "Current interaction: ", interaction_label, " (order ", length(facets), ")."
    )
  } else {
    facet_a <- facets[1]
    facet_b <- facets[2]
    if (is.null(target_facet)) target_facet <- facet_a
    context_facet <- ifelse(target_facet == facet_a, facet_b, facet_a)
    pairwise_tbl <- calc_bias_pairwise(bias_results$table, target_facet, context_facet)

    pairwise_fixed <- if (nrow(pairwise_tbl) == 0) {
      "No pairwise data"
    } else {
      pair_display <- pairwise_tbl |>
        dplyr::select(
          `Target N`,
          Target,
          `Target Measure`,
          `Target S.E.`,
          `Context1 N`,
          Context1,
          `Local Measure1`,
          SE1,
          `Obs-Exp Avg1`,
          Count1,
          `Context2 N`,
          Context2,
          `Local Measure2`,
          SE2,
          `Obs-Exp Avg2`,
          Count2,
          Contrast,
          SE,
          t,
          `d.f.`,
          `Prob.`
        ) |>
        dplyr::rename(
          `Target Measr` = `Target Measure`,
          `Context1 Measr` = `Local Measure1`,
          `Context1 S.E.` = SE1,
          `Context2 Measr` = `Local Measure2`,
          `Context2 S.E.` = SE2
        )

      pair_cols <- c(
        "Target N", "Target", "Target Measr", "Target S.E.",
        "Context1 N", "Context1", "Context1 Measr", "Context1 S.E.",
        "Obs-Exp Avg1", "Count1", "Context2 N", "Context2", "Context2 Measr", "Context2 S.E.",
        "Obs-Exp Avg2", "Count2", "Contrast", "SE", "t", "d.f.", "Prob."
      )

      pair_formats <- list(
        `Target N` = "{:.0f}",
        `Target Measr` = "{:.2f}",
        `Target S.E.` = "{:.2f}",
        `Context1 N` = "{:.0f}",
        `Context1 Measr` = "{:.2f}",
        `Context1 S.E.` = "{:.2f}",
        `Obs-Exp Avg1` = "{:.2f}",
        Count1 = "{:.0f}",
        `Context2 N` = "{:.0f}",
        `Context2 Measr` = "{:.2f}",
        `Context2 S.E.` = "{:.2f}",
        `Obs-Exp Avg2` = "{:.2f}",
        Count2 = "{:.0f}",
        Contrast = "{:.2f}",
        SE = "{:.2f}",
        t = "{:.2f}",
        `d.f.` = "{:.0f}",
        `Prob.` = "{:.4f}"
      )

      if (branch == "facets") {
        build_pairwise_fixed_text(
          pair_df = pair_display,
          target_facet = target_facet,
          context_facet = context_facet,
          columns = pair_cols,
          formats = pair_formats
        )
      } else {
        build_sectioned_fixed_report(
          title = paste0("Pairwise contrast summary: ", target_facet, " within ", context_facet),
          sections = list(
            list(
              title = "Pairwise rows",
              data = pair_display,
              columns = pair_cols,
              formats = pair_formats,
              max_rows = 40L
            )
          ),
          max_col_width = 18,
          min_col_width = 6
        )
      }
    }
  }

  out <- list(
    bias_fixed = bias_fixed,
    pairwise_fixed = pairwise_fixed,
    pairwise_table = pairwise_tbl,
    branch = branch,
    style = style,
    interaction_label = interaction_label,
    target_facet = as.character(target_facet %||% NA_character_)
  )
  out <- as_mfrm_bundle(out, "mfrm_fixed_reports")
  class(out) <- unique(c(paste0("mfrm_fixed_reports_", branch), class(out)))
  out
}

normalize_bias_plot_input <- function(x,
                                      diagnostics = NULL,
                                      facet_a = NULL,
                                      facet_b = NULL,
                                      interaction_facets = NULL,
                                      max_abs = 10,
                                      omit_extreme = TRUE,
                                      max_iter = 4,
                                      tol = 1e-3) {
  if (is.list(x) && !is.null(x$table) && !is.null(x$summary) && !is.null(x$chi_sq)) {
    return(x)
  }
  if (inherits(x, "mfrm_fit")) {
    if (is.null(interaction_facets)) {
      if (is.null(facet_a) || is.null(facet_b)) {
        stop("When `x` is mfrm_fit, provide `interaction_facets` or both `facet_a` and `facet_b`.")
      }
      interaction_facets <- c(as.character(facet_a[1]), as.character(facet_b[1]))
    }
    interaction_facets <- as.character(interaction_facets)
    interaction_facets <- interaction_facets[!is.na(interaction_facets) & nzchar(interaction_facets)]
    interaction_facets <- unique(interaction_facets)
    if (length(interaction_facets) < 2) {
      stop("`interaction_facets` must contain at least two facet names.")
    }
    if (is.null(facet_a)) facet_a <- interaction_facets[1]
    if (is.null(facet_b) && length(interaction_facets) >= 2) facet_b <- interaction_facets[2]
    if (is.null(facet_b)) {
      stop("`interaction_facets` must contain at least two facet names.")
    }
    if (is.null(diagnostics)) {
      diagnostics <- diagnose_mfrm(x, residual_pca = "none")
    }
    return(estimate_bias(
      fit = x,
      diagnostics = diagnostics,
      facet_a = facet_a,
      facet_b = facet_b,
      interaction_facets = interaction_facets,
      max_abs = max_abs,
      omit_extreme = omit_extreme,
      max_iter = max_iter,
      tol = tol
    ))
  }
  stop("`x` must be output from estimate_bias() or an mfrm_fit object.")
}

#' Build a FACETS Table 13-style bias-plot export bundle
#'
#' @param x Output from [estimate_bias()] or [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()] (used when `x` is fit).
#' @param facet_a First facet name (required when `x` is fit).
#' @param facet_b Second facet name (required when `x` is fit).
#' @param interaction_facets Character vector of two or more facets (required
#'   when `x` is fit and higher-order interaction output is needed).
#' @param max_abs Bound for absolute bias size when estimating from fit.
#' @param omit_extreme Omit extreme-only elements when estimating from fit.
#' @param max_iter Iteration cap for bias estimation when `x` is fit.
#' @param tol Convergence tolerance for bias estimation when `x` is fit.
#' @param top_n Maximum number of ranked rows to keep.
#' @param abs_t_warn Warning cutoff for absolute t statistics.
#' @param abs_bias_warn Warning cutoff for absolute bias size.
#' @param p_max Warning cutoff for p-values.
#' @param sort_by Ranking key: `"abs_t"`, `"abs_bias"`, or `"prob"`.
#'
#' @details
#' FACETS Table 13 is often inspected graphically (bias size, observed-minus-expected
#' average, significance). This helper prepares a plotting-ready bundle with:
#' - ranked table for lollipop/strip displays
#' - scatter table (`Obs-Exp Average` vs `Bias Size`)
#' - summary and threshold metadata
#'
#' @section Lifecycle:
#' Soft-deprecated. Prefer [bias_interaction_report()].
#'
#' @return A named list with:
#' - `ranked_table`: top-ranked bias rows with flags
#' - `scatter_data`: bias scatter data with flags
#' - `facet_profile`: per-facet level profile (`MeanAbsBias`, `FlagRate`)
#' - `summary`: one-row overview
#' - `thresholds`: applied cutoffs
#' - `facet_a`, `facet_b`: first two analyzed facet names
#' - `interaction_facets`, `interaction_order`, `interaction_mode`: full
#'   interaction metadata
#'
#' @seealso [estimate_bias()], [plot_bias_interaction()], [build_fixed_reports()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#' bias <- estimate_bias(fit, diag, facet_a = "Rater", facet_b = "Criterion", max_iter = 2)
#' t13 <- bias_interaction_report(bias, top_n = 10)
#' @keywords internal
#' @noRd
table13_bias_plot_export <- function(x,
                                     diagnostics = NULL,
                                     facet_a = NULL,
                                     facet_b = NULL,
                                     interaction_facets = NULL,
                                     max_abs = 10,
                                     omit_extreme = TRUE,
                                     max_iter = 4,
                                     tol = 1e-3,
                                     top_n = 50,
                                     abs_t_warn = 2,
                                     abs_bias_warn = 0.5,
                                     p_max = 0.05,
                                     sort_by = c("abs_t", "abs_bias", "prob")) {
  signal_legacy_name_deprecation(
    old_name = "table13_bias_plot_export",
    new_name = "bias_interaction_report",
    suppress_if_called_from = c("bias_interaction_report", "plot_table13_bias", "plot_bias_interaction")
  )
  sort_by <- match.arg(tolower(sort_by), c("abs_t", "abs_bias", "prob"))
  top_n <- max(1L, as.integer(top_n))
  abs_t_warn <- abs(as.numeric(abs_t_warn))
  abs_bias_warn <- abs(as.numeric(abs_bias_warn))
  p_max <- max(0, min(1, as.numeric(p_max)))

  bias_results <- normalize_bias_plot_input(
    x = x,
    diagnostics = diagnostics,
    facet_a = facet_a,
    facet_b = facet_b,
    interaction_facets = interaction_facets,
    max_abs = max_abs,
    omit_extreme = omit_extreme,
    max_iter = max_iter,
    tol = tol
  )
  spec <- extract_bias_facet_spec(bias_results)
  if (is.null(spec) || length(spec$facets) < 2) {
    stop("`bias_results$table` does not include recognizable interaction facet columns.")
  }
  interaction_facets <- spec$facets
  interaction_order <- spec$interaction_order
  interaction_mode <- spec$interaction_mode
  facet_a <- interaction_facets[1]
  facet_b <- interaction_facets[2]

  tbl <- as.data.frame(bias_results$table, stringsAsFactors = FALSE)
  req <- c(spec$level_cols, "Obs-Exp Average", "Bias Size", "t", "Prob.")
  if (!all(req %in% names(tbl))) {
    stop("`bias_results$table` does not include required Table 13 columns.")
  }

  level_df <- tbl[, spec$level_cols, drop = FALSE]
  level_df[] <- lapply(level_df, as.character)
  names(level_df) <- paste0("Level", seq_along(spec$level_cols))

  tbl2 <- dplyr::bind_cols(
    data.frame(
      InteractionFacets = paste(interaction_facets, collapse = " x "),
      InteractionOrder = interaction_order,
      InteractionMode = interaction_mode,
      FacetA = facet_a,
      FacetB = facet_b,
      stringsAsFactors = FALSE
    )[rep(1, nrow(tbl)), , drop = FALSE],
    level_df,
    tbl |>
      dplyr::transmute(
        ObsExpAverage = suppressWarnings(as.numeric(.data$`Obs-Exp Average`)),
        BiasSize = suppressWarnings(as.numeric(.data$`Bias Size`)),
        SE = if ("S.E." %in% names(tbl)) suppressWarnings(as.numeric(.data$`S.E.`)) else NA_real_,
        t = suppressWarnings(as.numeric(.data$t)),
        Prob = suppressWarnings(as.numeric(.data$`Prob.`)),
        ObservedCount = if ("Observd Count" %in% names(tbl)) suppressWarnings(as.numeric(.data$`Observd Count`)) else NA_real_
      )
  ) |>
    dplyr::mutate(
      Pair = do.call(paste, c(level_df, sep = " | ")),
      AbsT = abs(.data$t),
      AbsBias = abs(.data$BiasSize),
      TFlag = is.finite(.data$AbsT) & .data$AbsT >= abs_t_warn,
      BiasFlag = is.finite(.data$AbsBias) & .data$AbsBias >= abs_bias_warn,
      PFlag = is.finite(.data$Prob) & .data$Prob <= p_max,
      Flag = .data$TFlag | .data$BiasFlag | .data$PFlag
    )

  for (i in seq_along(interaction_facets)) {
    tbl2[[paste0("Facet", i)]] <- interaction_facets[i]
    tbl2[[paste0("Facet", i, "_Level")]] <- level_df[[i]]
  }
  # Keep legacy aliases for 2-way compatibility.
  tbl2$FacetA_Level <- tbl2$Facet1_Level
  tbl2$FacetB_Level <- tbl2$Facet2_Level

  ord <- switch(
    sort_by,
    abs_t = order(tbl2$AbsT, decreasing = TRUE, na.last = NA),
    abs_bias = order(tbl2$AbsBias, decreasing = TRUE, na.last = NA),
    prob = order(tbl2$Prob, decreasing = FALSE, na.last = NA)
  )
  ranked <- if (length(ord) == 0) tbl2[0, , drop = FALSE] else tbl2[ord, , drop = FALSE]
  if (nrow(ranked) > top_n) ranked <- ranked[seq_len(top_n), , drop = FALSE]

  scatter_cols <- c(
    "InteractionFacets", "InteractionOrder", "InteractionMode",
    "FacetA", "FacetA_Level", "FacetB", "FacetB_Level",
    paste0("Facet", seq_along(interaction_facets)),
    paste0("Facet", seq_along(interaction_facets), "_Level"),
    "Pair", "ObsExpAverage", "BiasSize", "SE", "t", "Prob",
    "ObservedCount", "Flag", "TFlag", "BiasFlag", "PFlag"
  )
  scatter_cols <- unique(scatter_cols)
  scatter <- tbl2 |>
    dplyr::select(dplyr::all_of(scatter_cols))

  profile_rows <- lapply(seq_along(interaction_facets), function(i) {
    facet_i <- interaction_facets[i]
    level_col <- paste0("Facet", i, "_Level")
    tbl2 |>
      dplyr::group_by(Level = .data[[level_col]]) |>
      dplyr::summarize(
        Cells = dplyr::n(),
        MeanAbsBias = mean(.data$AbsBias, na.rm = TRUE),
        MeanAbsT = mean(.data$AbsT, na.rm = TRUE),
        Flagged = sum(.data$Flag, na.rm = TRUE),
        FlagRate = ifelse(dplyr::n() > 0, 100 * sum(.data$Flag, na.rm = TRUE) / dplyr::n(), NA_real_),
        .groups = "drop"
      ) |>
      dplyr::mutate(Facet = facet_i, .before = 1)
  })
  facet_profile <- dplyr::bind_rows(profile_rows)

  summary_tbl <- data.frame(
    InteractionFacets = paste(interaction_facets, collapse = " x "),
    InteractionOrder = interaction_order,
    InteractionMode = interaction_mode,
    FacetA = facet_a,
    FacetB = facet_b,
    Cells = nrow(tbl2),
    Flagged = sum(tbl2$Flag, na.rm = TRUE),
    FlaggedPercent = ifelse(nrow(tbl2) > 0, 100 * sum(tbl2$Flag, na.rm = TRUE) / nrow(tbl2), NA_real_),
    MeanAbsT = mean(tbl2$AbsT, na.rm = TRUE),
    MeanAbsBias = mean(tbl2$AbsBias, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  list(
    ranked_table = as.data.frame(ranked, stringsAsFactors = FALSE),
    scatter_data = as.data.frame(scatter, stringsAsFactors = FALSE),
    facet_profile = as.data.frame(facet_profile, stringsAsFactors = FALSE),
    summary = summary_tbl,
    thresholds = list(
      abs_t_warn = abs_t_warn,
      abs_bias_warn = abs_bias_warn,
      p_max = p_max,
      sort_by = sort_by,
      top_n = top_n
    ),
    facet_a = facet_a,
    facet_b = facet_b,
    interaction_facets = interaction_facets,
    interaction_order = interaction_order,
    interaction_mode = interaction_mode
  )
}

#' Plot FACETS Table 13-style bias diagnostics using base R
#'
#' @param x Output from [bias_interaction_report()], [estimate_bias()], or [fit_mfrm()].
#' @param plot Plot type: `"scatter"`, `"ranked"`, `"abs_t_hist"`, or
#'   `"facet_profile"`.
#' @param diagnostics Optional output from [diagnose_mfrm()] when `x` is fit.
#' @param facet_a First facet name (required when `x` is fit).
#' @param facet_b Second facet name (required when `x` is fit).
#' @param interaction_facets Character vector of two or more facets (required
#'   when `x` is fit and higher-order interaction output is needed).
#' @param top_n Maximum number of ranked rows to show.
#' @param abs_t_warn Warning cutoff for absolute t.
#' @param abs_bias_warn Warning cutoff for absolute bias size.
#' @param p_max Warning cutoff for p-values.
#' @param sort_by Ranking key for `"ranked"` plot.
#' @param main Optional plot title override.
#' @param palette Optional named color overrides (`normal`, `flag`, `hist`,
#'   `profile`).
#' @param label_angle Label angle hint for ranked/profile labels.
#' @param draw If `TRUE`, draw with base graphics.
#'
#' @section Lifecycle:
#' Soft-deprecated. Prefer [plot_bias_interaction()].
#'
#' @return A plotting-data object of class `mfrm_plot_data`.
#' @seealso [bias_interaction_report()], [estimate_bias()], [plot_displacement()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' p13 <- plot_bias_interaction(
#'   fit,
#'   diagnostics = diagnose_mfrm(fit, residual_pca = "none"),
#'   facet_a = "Rater",
#'   facet_b = "Criterion",
#'   draw = FALSE
#' )
#' @keywords internal
#' @noRd
plot_table13_bias <- function(x,
                              plot = c("scatter", "ranked", "abs_t_hist", "facet_profile"),
                              diagnostics = NULL,
                              facet_a = NULL,
                              facet_b = NULL,
                              interaction_facets = NULL,
                              top_n = 40,
                              abs_t_warn = 2,
                              abs_bias_warn = 0.5,
                              p_max = 0.05,
                              sort_by = c("abs_t", "abs_bias", "prob"),
                              main = NULL,
                              palette = NULL,
                              label_angle = 45,
                              draw = TRUE) {
  signal_legacy_name_deprecation(
    old_name = "plot_table13_bias",
    new_name = "plot_bias_interaction",
    suppress_if_called_from = "plot_bias_interaction"
  )
  plot <- match.arg(tolower(plot), c("scatter", "ranked", "abs_t_hist", "facet_profile"))
  sort_by <- match.arg(tolower(sort_by), c("abs_t", "abs_bias", "prob"))
  top_n <- max(1L, as.integer(top_n))
  label_angle <- suppressWarnings(as.numeric(label_angle[1]))
  if (!is.finite(label_angle)) label_angle <- 45
  las_rank <- if (label_angle >= 45) 2 else 1
  pal <- resolve_palette(
    palette = palette,
    defaults = c(
      normal = "#2b8cbe",
      flag = "#cb181d",
      hist = "#c7e9c0",
      profile = "#756bb1"
    )
  )

  bundle <- if (is.list(x) && all(c("ranked_table", "scatter_data", "summary", "thresholds") %in% names(x))) {
    x
  } else {
    with_legacy_name_warning_suppressed(
      table13_bias_plot_export(
        x = x,
        diagnostics = diagnostics,
        facet_a = facet_a,
        facet_b = facet_b,
        interaction_facets = interaction_facets,
        top_n = top_n,
        abs_t_warn = abs_t_warn,
        abs_bias_warn = abs_bias_warn,
        p_max = p_max,
        sort_by = sort_by
      )
    )
  }

  ranked <- as.data.frame(bundle$ranked_table, stringsAsFactors = FALSE)
  scatter <- as.data.frame(bundle$scatter_data, stringsAsFactors = FALSE)
  profile <- as.data.frame(bundle$facet_profile %||% data.frame(), stringsAsFactors = FALSE)
  thr <- bundle$thresholds

  if (isTRUE(draw)) {
    if (plot == "scatter") {
      if (nrow(scatter) == 0) {
        graphics::plot.new()
        graphics::title(main = main %||% "Table 13 bias scatter")
        graphics::text(0.5, 0.5, "No data")
      } else {
        col <- ifelse(as.logical(scatter$Flag), pal["flag"], pal["normal"])
        graphics::plot(
          x = scatter$ObsExpAverage,
          y = scatter$BiasSize,
          pch = 16,
          col = col,
          xlab = "Obs-Exp Average",
          ylab = "Bias Size (logits)",
          main = main %||% "Table 13: Bias scatter"
        )
        graphics::abline(h = c(-thr$abs_bias_warn, 0, thr$abs_bias_warn), lty = c(2, 1, 2), col = c("gray45", "gray30", "gray45"))
        graphics::abline(v = 0, lty = 2, col = "gray45")
      }
    } else if (plot == "ranked") {
      if (nrow(ranked) == 0) {
        graphics::plot.new()
        graphics::title(main = main %||% "Table 13 ranked bias")
        graphics::text(0.5, 0.5, "No data")
      } else {
        sub <- ranked[seq_len(min(nrow(ranked), top_n)), , drop = FALSE]
        y <- rev(seq_len(nrow(sub)))
        vals <- rev(suppressWarnings(as.numeric(sub$BiasSize)))
        lbl <- truncate_axis_label(rev(as.character(sub$Pair)), width = 28L)
        col <- ifelse(rev(as.logical(sub$Flag)), pal["flag"], pal["normal"])
        graphics::plot(
          x = vals,
          y = y,
          type = "n",
          yaxt = "n",
          ylab = "",
          xlab = "Bias Size (logits)",
          main = main %||% "Table 13: Ranked bias size"
        )
        graphics::segments(0, y, vals, y, col = "gray60")
        graphics::points(vals, y, pch = 16, col = col)
        graphics::axis(side = 2, at = y, labels = lbl, las = las_rank, cex.axis = 0.75)
        graphics::abline(v = c(-thr$abs_bias_warn, 0, thr$abs_bias_warn), lty = c(2, 1, 2), col = c("gray45", "gray30", "gray45"))
      }
    } else if (plot == "abs_t_hist") {
      tvals <- abs(suppressWarnings(as.numeric(scatter$t)))
      tvals <- tvals[is.finite(tvals)]
      if (length(tvals) == 0) {
        graphics::plot.new()
        graphics::title(main = main %||% "Table 13 |t| distribution")
        graphics::text(0.5, 0.5, "No data")
      } else {
        graphics::hist(
          x = tvals,
          breaks = "FD",
          col = pal["hist"],
          border = "white",
          xlab = "|t|",
          main = main %||% "Table 13: |t| distribution"
        )
        graphics::abline(v = thr$abs_t_warn, lty = 2, col = "gray45")
      }
    } else {
      if (nrow(profile) == 0 || !all(c("Facet", "Level", "MeanAbsBias", "FlagRate") %in% names(profile))) {
        graphics::plot.new()
        graphics::title(main = main %||% "Table 13 facet profile")
        graphics::text(0.5, 0.5, "No data")
      } else {
        prof <- profile |>
          dplyr::mutate(
            MeanAbsBias = suppressWarnings(as.numeric(.data$MeanAbsBias)),
            FlagRate = suppressWarnings(as.numeric(.data$FlagRate))
          ) |>
          dplyr::filter(is.finite(.data$MeanAbsBias)) |>
          dplyr::arrange(.data$Facet, dplyr::desc(.data$MeanAbsBias))
        if (nrow(prof) == 0) {
          graphics::plot.new()
          graphics::title(main = main %||% "Table 13 facet profile")
          graphics::text(0.5, 0.5, "No data")
        } else {
          lbl <- truncate_axis_label(paste0(prof$Facet, ": ", prof$Level), width = 34L)
          cols <- ifelse(is.finite(prof$FlagRate) & prof$FlagRate > 0, pal["flag"], pal["profile"])
          graphics::dotchart(
            x = prof$MeanAbsBias,
            labels = lbl,
            pch = 16,
            col = cols,
            cex = 0.8,
            cex.axis = 0.75,
            xlab = "Mean |Bias Size| (logits)",
            main = main %||% "Table 13: Facet-level mean |Bias|"
          )
          graphics::abline(v = thr$abs_bias_warn, lty = 2, col = "gray45")
        }
      }
    }
  }

  out <- new_mfrm_plot_data(
    "table13_bias",
    list(
      plot = plot,
      ranked_table = ranked,
      scatter_data = scatter,
      facet_profile = profile,
      summary = bundle$summary,
      thresholds = thr
    )
  )
  invisible(out)
}

# Human-friendly API aliases (preferred names)

#' Build a specification summary report (preferred alias)
#'
#' @param fit Output from [fit_mfrm()].
#' @param title Optional analysis title.
#' @param data_file Optional data-file label (for reporting only).
#' @param output_file Optional output-file label (for reporting only).
#' @param include_fixed If `TRUE`, include a FACETS-style fixed-width text block.
#' @details
#' `summary(out)` is supported through `summary()`.
#' `plot(out)` is dispatched through `plot()` for class
#' `mfrm_specifications` (`type = "facet_elements"`,
#' `"anchor_constraints"`, `"convergence"`).
#'
#' @section Interpreting output:
#' - `header` / `data_spec`: run identity and model settings.
#' - `facet_labels`: facet sizes and labels.
#' - `convergence_control`: optimizer configuration and status.
#'
#' @section Typical workflow:
#' 1. Generate `specifications_report(fit)`.
#' 2. Verify model settings and convergence metadata.
#' 3. Use as Table 1-style documentation in reports.
#' @return A named list with specification-report components. Class:
#'   `mfrm_specifications`.
#' @seealso [fit_mfrm()], [data_quality_report()], [estimation_iteration_report()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' out <- specifications_report(fit, title = "Toy run")
#' summary(out)
#' p_spec <- plot(out, draw = FALSE)
#' class(p_spec)
#' @export
specifications_report <- function(fit,
                                  title = NULL,
                                  data_file = NULL,
                                  output_file = NULL,
                                  include_fixed = FALSE) {
  out <- with_legacy_name_warning_suppressed(
    table1_specifications(
      fit = fit,
      title = title,
      data_file = data_file,
      output_file = output_file,
      include_fixed = include_fixed
    )
  )
  as_mfrm_bundle(out, "mfrm_specifications")
}

#' Build a data quality summary report (preferred alias)
#'
#' @param fit Output from [fit_mfrm()].
#' @param data Optional raw data frame used for row-level audit.
#' @param person Optional person column name in `data`.
#' @param facets Optional facet column names in `data`.
#' @param score Optional score column name in `data`.
#' @param weight Optional weight column name in `data`.
#' @param include_fixed If `TRUE`, include a FACETS-style fixed-width text block.
#' @details
#' `summary(out)` is supported through `summary()`.
#' `plot(out)` is dispatched through `plot()` for class
#' `mfrm_data_quality` (`type = "row_audit"`, `"category_counts"`,
#' `"missing_rows"`).
#'
#' @section Interpreting output:
#' - `summary`: retained/dropped row overview.
#' - `row_audit`: reason-level breakdown for data issues.
#' - `category_counts`: post-filter category usage.
#' - `unknown_elements`: facet levels in raw data but not in fitted design.
#'
#' @section Typical workflow:
#' 1. Run `data_quality_report(...)` with raw data.
#' 2. Check row-audit and missing/unknown element sections.
#' 3. Resolve issues before final estimation/reporting.
#' @return A named list with data-quality report components. Class:
#'   `mfrm_data_quality`.
#' @seealso [fit_mfrm()], [describe_mfrm_data()], [specifications_report()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' out <- data_quality_report(
#'   fit, data = toy, person = "Person",
#'   facets = c("Rater", "Criterion"), score = "Score"
#' )
#' summary(out)
#' p_dq <- plot(out, draw = FALSE)
#' class(p_dq)
#' @export
data_quality_report <- function(fit,
                                data = NULL,
                                person = NULL,
                                facets = NULL,
                                score = NULL,
                                weight = NULL,
                                include_fixed = FALSE) {
  out <- with_legacy_name_warning_suppressed(
    table2_data_summary(
      fit = fit,
      data = data,
      person = person,
      facets = facets,
      score = score,
      weight = weight,
      include_fixed = include_fixed
    )
  )
  as_mfrm_bundle(out, "mfrm_data_quality")
}

#' Build an estimation-iteration report (preferred alias)
#'
#' @param fit Output from [fit_mfrm()].
#' @param max_iter Maximum replay iterations (excluding optional initial row).
#' @param reltol Stopping tolerance for replayed max-logit change.
#' @param include_prox If `TRUE`, include an initial pseudo-row labeled `PROX`.
#' @param include_fixed If `TRUE`, include a FACETS-style fixed-width text block.
#' @details
#' `summary(out)` is supported through `summary()`.
#' `plot(out)` is dispatched through `plot()` for class
#' `mfrm_iteration_report` (`type = "residual"`, `"logit_change"`,
#' `"objective"`).
#'
#' @section Interpreting output:
#' - `iterations`: trajectory of convergence indicators by iteration.
#' - `summary`: final status and stopping diagnostics.
#' - optional `PROX` row: pseudo-initial reference point when enabled.
#'
#' @section Typical workflow:
#' 1. Run `estimation_iteration_report(fit)`.
#' 2. Inspect plateau/stability patterns in summary/plot.
#' 3. Adjust optimization settings if convergence looks weak.
#' @return A named list with iteration-report components. Class:
#'   `mfrm_iteration_report`.
#' @seealso [fit_mfrm()], [specifications_report()], [data_quality_report()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' out <- estimation_iteration_report(fit, max_iter = 5)
#' summary(out)
#' p_iter <- plot(out, draw = FALSE)
#' class(p_iter)
#' @export
estimation_iteration_report <- function(fit,
                                        max_iter = 20,
                                        reltol = NULL,
                                        include_prox = TRUE,
                                        include_fixed = FALSE) {
  out <- with_legacy_name_warning_suppressed(
    table3_iteration_report(
      fit = fit,
      max_iter = max_iter,
      reltol = reltol,
      include_prox = include_prox,
      include_fixed = include_fixed
    )
  )
  as_mfrm_bundle(out, "mfrm_iteration_report")
}

#' Build a subset connectivity report (preferred alias)
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()].
#' @param top_n_subsets Optional maximum number of subset rows to keep.
#' @param min_observations Minimum observations required to keep a subset row.
#' @details
#' `summary(out)` is supported through `summary()`.
#' `plot(out)` is dispatched through `plot()` for class
#' `mfrm_subset_connectivity` (`type = "subset_observations"`,
#' `"facet_levels"`).
#'
#' @section Interpreting output:
#' - `summary`: number and size of connected subsets.
#' - subset table: whether data are fragmented into disconnected components.
#' - facet-level columns: where connectivity bottlenecks occur.
#'
#' @section Typical workflow:
#' 1. Run `subset_connectivity_report(fit)`.
#' 2. Confirm near-single-subset structure when possible.
#' 3. Use results to justify linking/anchoring strategy.
#' @return A named list with subset-connectivity components. Class:
#'   `mfrm_subset_connectivity`.
#' @seealso [diagnose_mfrm()], [measurable_summary_table()], [data_quality_report()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' out <- subset_connectivity_report(fit)
#' summary(out)
#' p_sub <- plot(out, draw = FALSE)
#' class(p_sub)
#' @export
subset_connectivity_report <- function(fit,
                                       diagnostics = NULL,
                                       top_n_subsets = NULL,
                                       min_observations = 0) {
  out <- with_legacy_name_warning_suppressed(
    table6_subsets_listing(
      fit = fit,
      diagnostics = diagnostics,
      top_n_subsets = top_n_subsets,
      min_observations = min_observations
    )
  )
  as_mfrm_bundle(out, "mfrm_subset_connectivity")
}

#' Build a facet statistics report (preferred alias)
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()].
#' @param metrics Numeric columns in `diagnostics$measures` to summarize.
#' @param ruler_width Width of the fixed-width ruler used for `M/S/Q/X` marks.
#' @details
#' `summary(out)` is supported through `summary()`.
#' `plot(out)` is dispatched through `plot()` for class
#' `mfrm_facet_statistics` (`type = "means"`, `"sds"`, `"ranges"`).
#'
#' @section Interpreting output:
#' - facet-level means/SD/ranges of selected metrics (`Estimate`, fit indices, `SE`).
#' - fixed-width ruler rows (`M/S/Q/X`) for FACETS-like profile scanning.
#'
#' @section Typical workflow:
#' 1. Run `facet_statistics_report(fit)`.
#' 2. Inspect summary/ranges for anomalous facets.
#' 3. Cross-check flagged facets with fit and chi-square diagnostics.
#' @return A named list with facet-statistics components. Class:
#'   `mfrm_facet_statistics`.
#' @seealso [diagnose_mfrm()], [summary.mfrm_fit()], [plot_facets_chisq()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' out <- facet_statistics_report(fit)
#' summary(out)
#' p_fs <- plot(out, draw = FALSE)
#' class(p_fs)
#' @export
facet_statistics_report <- function(fit,
                                    diagnostics = NULL,
                                    metrics = c("Estimate", "Infit", "Outfit", "SE"),
                                    ruler_width = 41) {
  out <- with_legacy_name_warning_suppressed(
    table6_2_facet_statistics(
      fit = fit,
      diagnostics = diagnostics,
      metrics = metrics,
      ruler_width = ruler_width
    )
  )
  as_mfrm_bundle(out, "mfrm_facet_statistics")
}

#' Build a category structure report (preferred alias)
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()].
#' @param theta_range Theta/logit range used to derive transition points.
#' @param theta_points Number of grid points used for transition-point search.
#' @param drop_unused If `TRUE`, remove zero-count categories from outputs.
#' @param include_fixed If `TRUE`, include a FACETS-style fixed-width text block.
#' @param fixed_max_rows Maximum rows per fixed-width section.
#'
#' @details
#' Preferred high-level API for category-structure diagnostics.
#' This wraps the legacy Table 8 bar/transition export and returns a stable
#' bundle interface for reporting and plotting.
#'
#' @section Interpreting output:
#' Key components include:
#' - category usage/fit table (count, expected, infit/outfit, ZSTD)
#' - threshold ordering and adjacent threshold gaps
#' - category transition-point table on the requested theta grid
#'
#' Practical read order:
#' 1. `summary(out)` for compact warnings and threshold ordering.
#' 2. `out$category_table` for sparse/misfitting categories.
#' 3. `out$transition_points` to inspect crossing structure.
#' 4. `plot(out)` for quick visual check.
#'
#' @section Typical workflow:
#' 1. [fit_mfrm()] -> model.
#' 2. [diagnose_mfrm()] -> residual/fit diagnostics (optional argument here).
#' 3. `category_structure_report()` -> category health snapshot.
#' 4. `summary()` and `plot()` for report-ready interpretation.
#' @return A named list with category-structure components. Class:
#'   `mfrm_category_structure`.
#' @seealso [rating_scale_table()], [category_curves_report()], [plot.mfrm_fit()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' out <- category_structure_report(fit)
#' summary(out)
#' names(out)
#' p_cs <- plot(out, draw = FALSE)
#' class(p_cs)
#' @export
category_structure_report <- function(fit,
                                      diagnostics = NULL,
                                      theta_range = c(-6, 6),
                                      theta_points = 241,
                                      drop_unused = FALSE,
                                      include_fixed = FALSE,
                                      fixed_max_rows = 200) {
  out <- with_legacy_name_warning_suppressed(
    table8_barchart_export(
      fit = fit,
      diagnostics = diagnostics,
      theta_range = theta_range,
      theta_points = theta_points,
      drop_unused = drop_unused,
      include_fixed = include_fixed,
      fixed_max_rows = fixed_max_rows
    )
  )
  as_mfrm_bundle(out, "mfrm_category_structure")
}

#' Build a category curve export bundle (preferred alias)
#'
#' @param fit Output from [fit_mfrm()].
#' @param theta_range Theta/logit range for curve coordinates.
#' @param theta_points Number of points on the theta grid.
#' @param digits Rounding digits for numeric graph output.
#' @param include_fixed If `TRUE`, include a FACETS-style fixed-width text block.
#' @param fixed_max_rows Maximum rows shown in fixed-width graph tables.
#'
#' @details
#' Preferred high-level API for category-probability curve exports.
#' Returns tidy curve coordinates and summary metadata for quick
#' plotting/report integration without calling low-level helpers directly.
#'
#' @section Interpreting output:
#' Use this report to inspect:
#' - where each category has highest probability across theta
#' - whether adjacent categories cross in expected order
#' - whether probability bands look compressed (often sparse categories)
#'
#' Recommended read order:
#' 1. `summary(out)` for compact diagnostics.
#' 2. `out$curve_points` (or equivalent curve table) for downstream graphics.
#' 3. `plot(out)` for a default visual check.
#'
#' @section Typical workflow:
#' 1. Fit model with [fit_mfrm()].
#' 2. Run `category_curves_report()` with suitable `theta_points`.
#' 3. Use `summary()` and `plot()`; export tables for manuscripts/dashboard use.
#' @return A named list with category-curve components. Class:
#'   `mfrm_category_curves`.
#' @seealso [category_structure_report()], [rating_scale_table()], [plot.mfrm_fit()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' out <- category_curves_report(fit, theta_points = 101)
#' summary(out)
#' names(out)
#' p_cc <- plot(out, draw = FALSE)
#' class(p_cc)
#' @export
category_curves_report <- function(fit,
                                   theta_range = c(-6, 6),
                                   theta_points = 241,
                                   digits = 4,
                                   include_fixed = FALSE,
                                   fixed_max_rows = 400) {
  out <- with_legacy_name_warning_suppressed(
    table8_curves_export(
      fit = fit,
      theta_range = theta_range,
      theta_points = theta_points,
      digits = digits,
      include_fixed = include_fixed,
      fixed_max_rows = fixed_max_rows
    )
  )
  as_mfrm_bundle(out, "mfrm_category_curves")
}

#' Build a bias-interaction plot-data bundle (preferred alias)
#'
#' @param x Output from [estimate_bias()] or [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()] (used when `x` is fit).
#' @param facet_a First facet name (required when `x` is fit and
#'   `interaction_facets` is not supplied).
#' @param facet_b Second facet name (required when `x` is fit and
#'   `interaction_facets` is not supplied).
#' @param interaction_facets Character vector of two or more facets.
#' @param max_abs Bound for absolute bias size when estimating from fit.
#' @param omit_extreme Omit extreme-only elements when estimating from fit.
#' @param max_iter Iteration cap for bias estimation when `x` is fit.
#' @param tol Convergence tolerance for bias estimation when `x` is fit.
#' @param top_n Maximum number of ranked rows to keep.
#' @param abs_t_warn Warning cutoff for absolute t statistics.
#' @param abs_bias_warn Warning cutoff for absolute bias size.
#' @param p_max Warning cutoff for p-values.
#' @param sort_by Ranking key: `"abs_t"`, `"abs_bias"`, or `"prob"`.
#'
#' @details
#' Preferred bundle API for interaction-bias diagnostics. The function can:
#' - use a precomputed bias object from [estimate_bias()], or
#' - estimate internally from `mfrm_fit` + facet specification.
#'
#' @section Interpreting output:
#' Focus on ranked rows where multiple criteria converge:
#' - large absolute t statistic
#' - large absolute bias size
#' - small p-value
#'
#' The bundle is optimized for downstream `summary()` and
#' [plot_bias_interaction()] views.
#'
#' @section Typical workflow:
#' 1. Run [estimate_bias()] (or provide `mfrm_fit` here).
#' 2. Build `bias_interaction_report(...)`.
#' 3. Review `summary(out)` and visualize with [plot_bias_interaction()].
#' @return A named list with bias-interaction plotting/report components. Class:
#'   `mfrm_bias_interaction`.
#' @seealso [estimate_bias()], [build_fixed_reports()], [plot_bias_interaction()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#' bias <- estimate_bias(fit, diag, facet_a = "Rater", facet_b = "Criterion", max_iter = 2)
#' out <- bias_interaction_report(bias, top_n = 10)
#' summary(out)
#' p_bi <- plot(out, draw = FALSE)
#' class(p_bi)
#' @export
bias_interaction_report <- function(x,
                                    diagnostics = NULL,
                                    facet_a = NULL,
                                    facet_b = NULL,
                                    interaction_facets = NULL,
                                    max_abs = 10,
                                    omit_extreme = TRUE,
                                    max_iter = 4,
                                    tol = 1e-3,
                                    top_n = 50,
                                    abs_t_warn = 2,
                                    abs_bias_warn = 0.5,
                                    p_max = 0.05,
                                    sort_by = c("abs_t", "abs_bias", "prob")) {
  out <- with_legacy_name_warning_suppressed(
    table13_bias_plot_export(
      x = x,
      diagnostics = diagnostics,
      facet_a = facet_a,
      facet_b = facet_b,
      interaction_facets = interaction_facets,
      max_abs = max_abs,
      omit_extreme = omit_extreme,
      max_iter = max_iter,
      tol = tol,
      top_n = top_n,
      abs_t_warn = abs_t_warn,
      abs_bias_warn = abs_bias_warn,
      p_max = p_max,
      sort_by = sort_by
    )
  )
  as_mfrm_bundle(out, "mfrm_bias_interaction")
}

#' Plot bias interaction diagnostics (preferred alias)
#'
#' @inheritParams bias_interaction_report
#' @param plot Plot type: `"scatter"`, `"ranked"`, `"abs_t_hist"`,
#'   or `"facet_profile"`.
#' @param main Optional plot title override.
#' @param palette Optional named color overrides (`normal`, `flag`, `hist`,
#'   `profile`).
#' @param label_angle Label angle hint for ranked/profile labels.
#' @param draw If `TRUE`, draw with base graphics.
#'
#' @details
#' Visualization front-end for [bias_interaction_report()] with multiple views.
#'
#' @section Plot types:
#' \describe{
#'   \item{`"scatter"` (default)}{Scatter plot of bias size (x) vs
#'     t-statistic (y).  Points colored by flag status.  Dashed reference
#'     lines at `abs_bias_warn` and `abs_t_warn`.  Use for overall triage
#'     of interaction effects.}
#'   \item{`"ranked"`}{Ranked bar chart of top `top_n` interactions sorted
#'     by `sort_by` criterion (absolute t, absolute bias, or probability).
#'     Bars colored red for flagged cells.}
#'   \item{`"abs_t_hist"`}{Histogram of absolute t-statistics across all
#'     interaction cells.  Dashed reference line at `abs_t_warn`.  Use for
#'     assessing the overall distribution of interaction effect sizes.}
#'   \item{`"facet_profile"`}{Per-facet-level aggregation showing mean
#'     absolute bias and flag rate.  Useful for identifying which
#'     individual facet levels drive systematic interaction patterns.}
#' }
#'
#' @section Interpreting output:
#' Start with `"scatter"` or `"ranked"` for triage, then confirm pattern shape
#' using `"abs_t_hist"` and `"facet_profile"`.
#'
#' Consistent flags across multiple views are stronger evidence of systematic
#' interaction bias than a single extreme row.
#'
#' @section Typical workflow:
#' 1. Estimate bias with [estimate_bias()] or pass `mfrm_fit` directly.
#' 2. Plot with `plot = "ranked"` for top interactions.
#' 3. Cross-check using `plot = "scatter"` and `plot = "facet_profile"`.
#' @return A plotting-data object of class `mfrm_plot_data`.
#' @seealso [bias_interaction_report()], [estimate_bias()], [plot_displacement()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' p <- plot_bias_interaction(
#'   fit,
#'   diagnostics = diagnose_mfrm(fit, residual_pca = "none"),
#'   facet_a = "Rater",
#'   facet_b = "Criterion",
#'   draw = FALSE
#' )
#' @export
plot_bias_interaction <- function(x,
                                  plot = c("scatter", "ranked", "abs_t_hist", "facet_profile"),
                                  diagnostics = NULL,
                                  facet_a = NULL,
                                  facet_b = NULL,
                                  interaction_facets = NULL,
                                  top_n = 40,
                                  abs_t_warn = 2,
                                  abs_bias_warn = 0.5,
                                  p_max = 0.05,
                                  sort_by = c("abs_t", "abs_bias", "prob"),
                                  main = NULL,
                                  palette = NULL,
                                  label_angle = 45,
                                  draw = TRUE) {
  with_legacy_name_warning_suppressed(
    plot_table13_bias(
      x = x,
      plot = plot,
      diagnostics = diagnostics,
      facet_a = facet_a,
      facet_b = facet_b,
      interaction_facets = interaction_facets,
      top_n = top_n,
      abs_t_warn = abs_t_warn,
      abs_bias_warn = abs_bias_warn,
      p_max = p_max,
      sort_by = sort_by,
      main = main,
      palette = palette,
      label_angle = label_angle,
      draw = draw
    )
  )
}

#' Build APA text outputs from model results
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Output from [diagnose_mfrm()].
#' @param bias_results Optional output from [estimate_bias()].
#' @param context Optional named list for report context.
#' @param whexact Use exact ZSTD transformation.
#'
#' @details
#' `context` is an optional named list for narrative customization.
#' Frequently used fields include:
#' - `assessment`, `setting`, `scale_desc`
#' - `rater_training`, `raters_per_response`
#' - `rater_facet` (used for targeted reliability note text)
#' - `line_width` (optional text wrapping width for `report_text`; default = 92)
#'
#' Output text includes residual PCA interpretation if PCA diagnostics are
#' available in `diagnostics`.
#'
#' By default, `report_text` includes:
#' - model/data design summary (N, facet counts, scale range)
#' - optimization/convergence metrics (`Converged`, `Iterations`, `LogLik`, `AIC`, `BIC`)
#' - anchor/constraint summary (`noncenter_facet`, anchored levels, group anchors, dummy facets)
#' - category/threshold diagnostics (including disordered-step details when present)
#' - overall fit, misfit count, and top misfit levels
#' - facet reliability/separation, residual PCA summary, and bias-screen counts
#'
#' @section Interpreting output:
#' - `report_text`: manuscript-ready narrative core.
#' - `table_figure_notes`: reusable note blocks for table/figure appendices.
#' - `table_figure_captions`: caption candidates aligned to generated outputs.
#'
#' @section Typical workflow:
#' 1. Build diagnostics (and optional bias results).
#' 2. Run `build_apa_outputs(...)`.
#' 3. Check `summary(apa)` for completeness.
#' 4. Insert `apa$report_text` and note/caption fields into manuscript drafts.
#'
#' @section Context template:
#' A minimal `context` list can include fields such as:
#' - `assessment`: name of the assessment task
#' - `setting`: administration context
#' - `scale_desc`: short description of the score scale
#' - `rater_facet`: rater facet label used in narrative reliability text
#'
#' @return
#' An object of class `mfrm_apa_outputs` with:
#' - `report_text`: APA-style Method/Results prose
#' - `table_figure_notes`: consolidated notes for tables/visuals
#' - `table_figure_captions`: caption-ready labels without figure numbering
#'
#' @seealso [build_visual_summaries()], [estimate_bias()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' diag <- diagnose_mfrm(fit, residual_pca = "both")
#' apa <- build_apa_outputs(
#'   fit,
#'   diag,
#'   context = list(
#'     assessment = "Toy writing task",
#'     setting = "Demonstration dataset",
#'     scale_desc = "0-2 rating scale",
#'     rater_facet = "Rater"
#'   )
#' )
#' names(apa)
#' class(apa)
#' s_apa <- summary(apa)
#' s_apa$overview
#' cat(apa$report_text)
#' @export
build_apa_outputs <- function(fit,
                              diagnostics,
                              bias_results = NULL,
                              context = list(),
                              whexact = FALSE) {
  report_text <- build_apa_report_text(
    res = fit,
    diagnostics = diagnostics,
    bias_results = bias_results,
    context = context,
    whexact = whexact
  )

  out <- list(
    report_text = structure(
      as.character(report_text),
      class = c("mfrm_apa_text", "character")
    ),
    table_figure_notes = build_apa_table_figure_notes(
      res = fit,
      diagnostics = diagnostics,
      bias_results = bias_results,
      context = context,
      whexact = whexact
    ),
    table_figure_captions = build_apa_table_figure_captions(
      res = fit,
      diagnostics = diagnostics,
      bias_results = bias_results,
      context = context
    )
  )
  class(out) <- c("mfrm_apa_outputs", "list")
  out
}

#' Print APA narrative text with preserved line breaks
#'
#' @param x Character text object from `build_apa_outputs()$report_text`.
#' @param ... Reserved for generic compatibility.
#'
#' @details
#' Prints APA narrative text with preserved paragraph breaks using `cat()`.
#' This is preferred over bare `print()` when you want readable multi-line
#' report output in the console.
#'
#' @section Interpreting output:
#' The printed text is the same content stored in
#' `build_apa_outputs(...)$report_text`, but with explicit paragraph breaks.
#'
#' @section Typical workflow:
#' 1. Generate `apa <- build_apa_outputs(...)`.
#' 2. Print readable narrative with `apa$report_text`.
#' 3. Use `summary(apa)` to check completeness before manuscript use.
#'
#' @return The input object (invisibly).
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' diag <- diagnose_mfrm(fit, residual_pca = "both")
#' apa <- build_apa_outputs(fit, diag)
#' apa$report_text
#' @export
print.mfrm_apa_text <- function(x, ...) {
  cat(as.character(x), "\n", sep = "")
  invisible(x)
}

#' Summarize APA report-output bundles
#'
#' @param object Output from [build_apa_outputs()].
#' @param top_n Maximum non-empty lines shown in each component preview.
#' @param preview_chars Maximum characters shown in each preview cell.
#' @param ... Reserved for generic compatibility.
#'
#' @details
#' This summary is a diagnostics layer for APA text products, not a replacement
#' for the full narrative.
#'
#' It reports component completeness, line/character volume, and a compact
#' preview for quick QA before manuscript insertion.
#'
#' @section Interpreting output:
#' - `overview`: total coverage across standard text components.
#' - `components`: per-component density and mention checks
#'   (including residual-PCA mentions).
#' - `preview`: first non-empty lines for fast visual review.
#'
#' @section Typical workflow:
#' 1. Build outputs via [build_apa_outputs()].
#' 2. Run `summary(apa)` to screen for empty/short components.
#' 3. Use `apa$report_text`, `apa$table_figure_notes`,
#'    and `apa$table_figure_captions` directly for final text.
#'
#' @return An object of class `summary.mfrm_apa_outputs`.
#' @seealso [build_apa_outputs()], [summary()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' diag <- diagnose_mfrm(fit, residual_pca = "both")
#' apa <- build_apa_outputs(fit, diag)
#' summary(apa)
#' @export
summary.mfrm_apa_outputs <- function(object, top_n = 3, preview_chars = 160, ...) {
  if (!inherits(object, "mfrm_apa_outputs")) {
    stop("`object` must be an mfrm_apa_outputs object from build_apa_outputs().", call. = FALSE)
  }

  top_n <- max(1L, as.integer(top_n))
  preview_chars <- max(40L, as.integer(preview_chars))

  text_line_count <- function(text) {
    if (!nzchar(text)) return(0L)
    length(strsplit(text, "\n", fixed = TRUE)[[1]])
  }
  nonempty_line_count <- function(text) {
    if (!nzchar(text)) return(0L)
    lines <- strsplit(text, "\n", fixed = TRUE)[[1]]
    sum(nzchar(trimws(lines)))
  }
  text_preview <- function(text, top_n, preview_chars) {
    if (!nzchar(text)) return("")
    lines <- strsplit(text, "\n", fixed = TRUE)[[1]]
    lines <- trimws(lines)
    lines <- lines[nzchar(lines)]
    if (length(lines) == 0) return("")
    pv <- paste(utils::head(lines, n = top_n), collapse = " | ")
    if (nchar(pv) > preview_chars) {
      pv <- paste0(substr(pv, 1, preview_chars - 3), "...")
    }
    pv
  }

  components <- c("report_text", "table_figure_notes", "table_figure_captions")
  stats_tbl <- do.call(
    rbind,
    lapply(components, function(comp) {
      text_vec <- as.character(object[[comp]] %||% character(0))
      text <- paste(text_vec, collapse = "\n")
      data.frame(
        Component = comp,
        NonEmpty = nzchar(trimws(text)),
        Characters = nchar(text),
        Lines = text_line_count(text),
        NonEmptyLines = nonempty_line_count(text),
        ResidualPCA_Mentions = stringr::str_count(
          text,
          stringr::regex("Residual\\s*PCA", ignore_case = TRUE)
        ),
        stringsAsFactors = FALSE
      )
    })
  )

  preview_tbl <- do.call(
    rbind,
    lapply(components, function(comp) {
      text_vec <- as.character(object[[comp]] %||% character(0))
      text <- paste(text_vec, collapse = "\n")
      data.frame(
        Component = comp,
        Preview = text_preview(text, top_n = top_n, preview_chars = preview_chars),
        stringsAsFactors = FALSE
      )
    })
  )

  overview <- data.frame(
    Components = nrow(stats_tbl),
    NonEmptyComponents = sum(stats_tbl$NonEmpty),
    TotalCharacters = sum(stats_tbl$Characters),
    TotalNonEmptyLines = sum(stats_tbl$NonEmptyLines),
    stringsAsFactors = FALSE
  )

  empty_components <- stats_tbl$Component[!stats_tbl$NonEmpty]
  notes <- if (length(empty_components) == 0) {
    c(
      "All standard APA text components are populated.",
      "Use object fields directly for full text; summary provides compact diagnostics."
    )
  } else {
    c(
      paste0("Empty components: ", paste(empty_components, collapse = ", "), "."),
      "Use object fields directly for full text; summary provides compact diagnostics."
    )
  }

  out <- list(
    overview = overview,
    components = stats_tbl,
    preview = preview_tbl,
    notes = notes,
    top_n = top_n,
    preview_chars = preview_chars
  )
  class(out) <- "summary.mfrm_apa_outputs"
  out
}

#' @export
print.summary.mfrm_apa_outputs <- function(x, ...) {
  cat("mfrmr APA Outputs Summary\n")

  if (!is.null(x$overview) && nrow(x$overview) > 0) {
    cat("\nOverview\n")
    print(round_numeric_df(as.data.frame(x$overview), digits = 0), row.names = FALSE)
  }
  if (!is.null(x$components) && nrow(x$components) > 0) {
    cat("\nComponent stats\n")
    print(round_numeric_df(as.data.frame(x$components), digits = 0), row.names = FALSE)
  }
  if (!is.null(x$preview) && nrow(x$preview) > 0) {
    cat("\nPreview\n")
    print(as.data.frame(x$preview), row.names = FALSE)
  }
  if (length(x$notes) > 0) {
    cat("\nNotes\n")
    cat(" - ", x$notes, "\n", sep = "")
  }
  invisible(x)
}

#' Build APA-style table output using base R structures
#'
#' @param x A data.frame, `mfrm_fit`, diagnostics list, or bias-result list.
#' @param which Optional table selector when `x` has multiple tables.
#' @param diagnostics Optional diagnostics from [diagnose_mfrm()] (used when
#'   `x` is `mfrm_fit` and `which` targets diagnostics tables).
#' @param digits Number of rounding digits for numeric columns.
#' @param caption Optional caption text.
#' @param note Optional note text.
#' @param branch Output branch:
#'   `"apa"` for manuscript-oriented labels, `"facets"` for FACETS-aligned labels.
#'
#' @details
#' This helper avoids styling dependencies and returns a reproducible base
#' `data.frame` plus metadata.
#'
#' Supported `which` values:
#' - For `mfrm_fit`: `"summary"`, `"person"`, `"facets"`, `"steps"`
#' - For diagnostics list: `"overall_fit"`, `"measures"`, `"fit"`,
#'   `"reliability"`, `"facets_chisq"`, `"bias"`, `"interactions"`,
#'   `"interrater_summary"`, `"interrater_pairs"`, `"obs"`
#' - For bias-result list: `"table"`, `"summary"`, `"chi_sq"`
#'
#' @section Interpreting output:
#' - `table`: plain data.frame ready for export or further formatting.
#' - `which`: source component that produced the table.
#' - `caption`/`note`: manuscript-oriented metadata stored with the table.
#'
#' @section Typical workflow:
#' 1. Build table object with `apa_table(...)`.
#' 2. Inspect quickly with `summary(tbl)`.
#' 3. Render base preview via `plot(tbl, ...)` or export `tbl$table`.
#'
#' @return A list of class `apa_table` with fields:
#' - `table` (`data.frame`)
#' - `which`
#' - `caption`
#' - `note`
#' - `digits`
#' - `branch`, `style`
#' @seealso [fit_mfrm()], [diagnose_mfrm()], [build_apa_outputs()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' tbl <- apa_table(fit, which = "summary", caption = "Model summary", note = "Toy example")
#' tbl_facets <- apa_table(fit, which = "summary", branch = "facets")
#' summary(tbl)
#' p <- plot(tbl, draw = FALSE)
#' p_facets <- plot(tbl_facets, type = "numeric_profile", draw = FALSE)
#' if (interactive()) {
#'   plot(
#'     tbl,
#'     type = "numeric_profile",
#'     main = "APA Table Numeric Profile (Customized)",
#'     palette = c(numeric_profile = "#2b8cbe", grid = "#d9d9d9"),
#'     label_angle = 45
#'   )
#' }
#' class(tbl)
#' @export
apa_table <- function(x,
                      which = NULL,
                      diagnostics = NULL,
                      digits = 2,
                      caption = NULL,
                      note = NULL,
                      branch = c("apa", "facets")) {
  branch <- match.arg(tolower(as.character(branch[1])), c("apa", "facets"))
  style <- ifelse(branch == "facets", "facets_manual", "apa")
  digits <- max(0L, as.integer(digits))
  table_out <- NULL
  source_type <- "data.frame"

  if (is.data.frame(x)) {
    table_out <- x
    source_type <- "data.frame"
  } else if (inherits(x, "mfrm_fit")) {
    source_type <- "mfrm_fit"
    opts <- c("summary", "person", "facets", "steps")
    diag_opts <- c(
      "overall_fit",
      "measures",
      "fit",
      "reliability",
      "facets_chisq",
      "bias",
      "interactions",
      "interrater_summary",
      "interrater_pairs",
      "obs"
    )
    if (is.null(which)) which <- "summary"
    which <- tolower(as.character(which[1]))

    if (which %in% opts) {
      table_out <- switch(
        which,
        summary = x$summary,
        person = x$facets$person,
        facets = x$facets$others,
        steps = x$steps
      )
    } else if (which %in% diag_opts) {
      if (is.null(diagnostics)) {
        diagnostics <- diagnose_mfrm(x, residual_pca = "none")
      }
      if (which == "interrater_summary") {
        table_out <- diagnostics$interrater$summary
      } else if (which == "interrater_pairs") {
        table_out <- diagnostics$interrater$pairs
      } else {
        table_out <- diagnostics[[which]]
      }
    } else {
      stop("Unsupported `which` for mfrm_fit. Use one of: ", paste(c(opts, diag_opts), collapse = ", "))
    }
  } else if (is.list(x) && !is.null(names(x))) {
    source_type <- "list"
    candidate <- names(x)
    if (is.null(which)) {
      pref <- c(
        "summary", "table", "overall_fit", "measures", "fit", "reliability", "facets_chisq",
        "bias", "interactions", "interrater_summary", "interrater_pairs", "obs", "chi_sq"
      )
      hit <- pref[pref %in% candidate]
      if (length(hit) == 0) {
        stop("Could not infer `which` from list input. Please specify `which`.")
      }
      which <- hit[1]
    }
    which <- as.character(which[1])
    if (!which %in% names(x)) {
      stop("Requested `which` not found in list input.")
    }
    table_out <- x[[which]]
  } else {
    stop("`x` must be a data.frame, mfrm_fit, or named list.")
  }

  if (is.null(table_out)) {
    table_out <- data.frame()
  }
  table_out <- as.data.frame(table_out, stringsAsFactors = FALSE)
  if (nrow(table_out) > 0) {
    num_cols <- vapply(table_out, is.numeric, logical(1))
    table_out[num_cols] <- lapply(table_out[num_cols], round, digits = digits)
  }

  out <- list(
    table = table_out,
    which = if (is.null(which)) source_type else as.character(which),
    caption = if (is.null(caption)) {
      if (branch == "facets") {
        paste0("FACETS-aligned table: ", if (is.null(which)) source_type else as.character(which))
      } else {
        ""
      }
    } else {
      as.character(caption)
    },
    note = if (is.null(note)) "" else as.character(note),
    digits = digits,
    branch = branch,
    style = style
  )
  class(out) <- c(paste0("apa_table_", branch), "apa_table", class(out))
  out
}

#' @export
print.apa_table <- function(x, ...) {
  if (!is.null(x$caption) && nzchar(x$caption)) {
    cat(x$caption, "\n", sep = "")
  }
  if (is.data.frame(x$table) && nrow(x$table) > 0) {
    print(x$table, row.names = FALSE)
  } else {
    cat("<empty table>\n")
  }
  if (!is.null(x$note) && nzchar(x$note)) {
    cat("Note. ", x$note, "\n", sep = "")
  }
  invisible(x)
}

#' Summarize an APA/FACETS table object
#'
#' @param object Output from [apa_table()].
#' @param digits Number of digits used for numeric summaries.
#' @param top_n Maximum numeric columns shown in `numeric_profile`.
#' @param ... Reserved for generic compatibility.
#'
#' @details
#' Compact summary helper for QA of table payloads before manuscript export.
#'
#' @section Interpreting output:
#' - `overview`: table size/composition and missingness.
#' - `numeric_profile`: quick distribution summary of numeric columns.
#' - `caption`/`note`: text metadata readiness.
#'
#' @section Typical workflow:
#' 1. Build table with [apa_table()].
#' 2. Run `summary(tbl)` and inspect `overview`.
#' 3. Use [plot.apa_table()] for quick numeric checks if needed.
#'
#' @return An object of class `summary.apa_table`.
#' @seealso [apa_table()], [plot()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' tbl <- apa_table(fit, which = "summary")
#' summary(tbl)
#' @export
summary.apa_table <- function(object, digits = 3, top_n = 8, ...) {
  digits <- max(0L, as.integer(digits))
  top_n <- max(1L, as.integer(top_n))
  tbl <- as.data.frame(object$table %||% data.frame(), stringsAsFactors = FALSE)

  num_cols <- names(tbl)[vapply(tbl, is.numeric, logical(1))]
  numeric_profile <- data.frame()
  if (length(num_cols) > 0) {
    numeric_profile <- do.call(
      rbind,
      lapply(num_cols, function(nm) {
        vals <- suppressWarnings(as.numeric(tbl[[nm]]))
        vals <- vals[is.finite(vals)]
        data.frame(
          Column = nm,
          N = length(vals),
          Mean = if (length(vals) > 0) mean(vals) else NA_real_,
          SD = if (length(vals) > 1) stats::sd(vals) else NA_real_,
          Min = if (length(vals) > 0) min(vals) else NA_real_,
          Max = if (length(vals) > 0) max(vals) else NA_real_,
          stringsAsFactors = FALSE
        )
      })
    )
    numeric_profile <- numeric_profile |>
      dplyr::arrange(dplyr::desc(.data$SD), .data$Column) |>
      dplyr::slice_head(n = top_n)
  }

  overview <- data.frame(
    Branch = as.character(object$branch %||% "apa"),
    Style = as.character(object$style %||% "apa"),
    Which = as.character(object$which %||% ""),
    Rows = nrow(tbl),
    Columns = ncol(tbl),
    NumericColumns = length(num_cols),
    MissingValues = sum(is.na(tbl)),
    stringsAsFactors = FALSE
  )

  out <- list(
    overview = overview,
    numeric_profile = numeric_profile,
    caption = as.character(object$caption %||% ""),
    note = as.character(object$note %||% ""),
    digits = digits
  )
  class(out) <- "summary.apa_table"
  out
}

#' @export
print.summary.apa_table <- function(x, ...) {
  digits <- as.integer(x$digits %||% 3L)
  if (!is.finite(digits)) digits <- 3L

  cat("APA Table Summary\n")
  if (!is.null(x$overview) && nrow(x$overview) > 0) {
    print(round_numeric_df(as.data.frame(x$overview), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$caption) && nzchar(x$caption)) {
    cat("\nCaption\n")
    cat(" - ", x$caption, "\n", sep = "")
  }
  if (!is.null(x$note) && nzchar(x$note)) {
    cat("\nNote\n")
    cat(" - ", x$note, "\n", sep = "")
  }
  if (!is.null(x$numeric_profile) && nrow(x$numeric_profile) > 0) {
    cat("\nNumeric profile\n")
    print(round_numeric_df(as.data.frame(x$numeric_profile), digits = digits), row.names = FALSE)
  }
  invisible(x)
}

#' Plot an APA/FACETS table object using base R
#'
#' @param x Output from [apa_table()].
#' @param y Reserved for generic compatibility.
#' @param type Plot type: `"numeric_profile"` (column means) or
#'   `"first_numeric"` (distribution of the first numeric column).
#' @param main Optional title override.
#' @param palette Optional named color overrides.
#' @param label_angle Axis-label rotation angle for bar-type plots.
#' @param draw If `TRUE`, draw using base graphics.
#' @param ... Reserved for generic compatibility.
#'
#' @details
#' Quick visualization helper for numeric columns in [apa_table()] output.
#' It is intended for table QA and exploratory checks, not final publication
#' graphics.
#'
#' @section Interpreting output:
#' - `"numeric_profile"`: compares column means to spot scale/centering mismatches.
#' - `"first_numeric"`: checks distribution shape of the first numeric column.
#'
#' @section Typical workflow:
#' 1. Build table with [apa_table()].
#' 2. Run `summary(tbl)` for metadata.
#' 3. Use `plot(tbl, type = "numeric_profile")` for quick numeric QC.
#'
#' @return A plotting-data object of class `mfrm_plot_data`.
#' @seealso [apa_table()], [summary()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' tbl <- apa_table(fit, which = "summary")
#' p <- plot(tbl, draw = FALSE)
#' p2 <- plot(tbl, type = "first_numeric", draw = FALSE)
#' if (interactive()) {
#'   plot(
#'     tbl,
#'     type = "numeric_profile",
#'     main = "APA Numeric Profile (Customized)",
#'     palette = c(numeric_profile = "#2b8cbe", grid = "#d9d9d9"),
#'     label_angle = 45
#'   )
#' }
#' @export
plot.apa_table <- function(x,
                           y = NULL,
                           type = c("numeric_profile", "first_numeric"),
                           main = NULL,
                           palette = NULL,
                           label_angle = 45,
                           draw = TRUE,
                           ...) {
  type <- match.arg(tolower(as.character(type[1])), c("numeric_profile", "first_numeric"))
  tbl <- as.data.frame(x$table %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(tbl) == 0) {
    stop("`x$table` is empty.")
  }
  num_cols <- names(tbl)[vapply(tbl, is.numeric, logical(1))]
  if (length(num_cols) == 0) {
    stop("`x$table` has no numeric columns to plot.")
  }
  pal <- resolve_palette(
    palette = palette,
    defaults = c(
      numeric_profile = "#1f78b4",
      first_numeric = "#33a02c",
      grid = "#ececec"
    )
  )

  if (type == "numeric_profile") {
    vals <- vapply(num_cols, function(nm) {
      v <- suppressWarnings(as.numeric(tbl[[nm]]))
      mean(v[is.finite(v)])
    }, numeric(1))
    ord <- order(abs(vals), decreasing = TRUE, na.last = NA)
    vals <- vals[ord]
    labels <- num_cols[ord]
    if (isTRUE(draw)) {
      barplot_rot45(
        height = vals,
        labels = labels,
        col = pal["numeric_profile"],
        main = if (is.null(main)) "APA table numeric profile (column means)" else as.character(main[1]),
        ylab = "Mean",
        label_angle = label_angle,
        mar_bottom = 8.8
      )
      graphics::abline(h = 0, col = pal["grid"], lty = 2)
    }
    out <- new_mfrm_plot_data(
      "apa_table",
      list(plot = "numeric_profile", column = labels, mean = vals)
    )
    return(invisible(out))
  }

  nm <- num_cols[1]
  vals <- suppressWarnings(as.numeric(tbl[[nm]]))
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0) {
    stop("First numeric column does not contain finite values.")
  }
  if (isTRUE(draw)) {
    graphics::hist(
      x = vals,
      breaks = "FD",
      col = pal["first_numeric"],
      border = "white",
      main = if (is.null(main)) paste0("Distribution of ", nm) else as.character(main[1]),
      xlab = nm,
      ylab = "Count"
    )
  }
  out <- new_mfrm_plot_data(
    "apa_table",
    list(plot = "first_numeric", column = nm, values = vals)
  )
  invisible(out)
}

#' List literature-based warning threshold profiles
#'
#' @return An object of class `mfrm_threshold_profiles` with
#'   `profiles` (`strict`, `standard`, `lenient`) and `pca_reference_bands`.
#' @details
#' Use this function to inspect available profile presets before calling
#' [build_visual_summaries()].
#'
#' `profiles` contains thresholds used by warning logic
#' (sample size, fit ratios, PCA cutoffs, etc.).
#' `pca_reference_bands` contains literature-oriented descriptive bands used in
#' summary text.
#'
#' @section Interpreting output:
#' - `profiles`: numeric threshold presets (`strict`, `standard`, `lenient`).
#' - `pca_reference_bands`: narrative reference bands for PCA interpretation.
#'
#' @section Typical workflow:
#' 1. Review presets with `mfrm_threshold_profiles()`.
#' 2. Pick a default profile for project policy.
#' 3. Override only selected fields in [build_visual_summaries()] when needed.
#'
#' @seealso [build_visual_summaries()]
#' @examples
#' profiles <- mfrm_threshold_profiles()
#' names(profiles)
#' names(profiles$profiles)
#' class(profiles)
#' s_profiles <- summary(profiles)
#' s_profiles$overview
#' @export
mfrm_threshold_profiles <- function() {
  out <- warning_threshold_profiles()
  class(out) <- c("mfrm_threshold_profiles", "list")
  out
}

#' Summarize threshold-profile presets for visual warning logic
#'
#' @param object Output from [mfrm_threshold_profiles()].
#' @param digits Number of digits used for numeric summaries.
#' @param ... Reserved for generic compatibility.
#'
#' @details
#' Summarizes available warning presets and their PCA reference bands used by
#' [build_visual_summaries()].
#'
#' @section Interpreting output:
#' - `thresholds`: raw preset values by profile (`strict`, `standard`, `lenient`).
#' - `threshold_ranges`: per-threshold span across profiles (sensitivity to profile choice).
#' - `pca_reference`: literature bands used for PCA narrative labeling.
#'
#' Larger `Span` in `threshold_ranges` indicates settings that most change
#' warning behavior between strict and lenient modes.
#'
#' @section Typical workflow:
#' 1. Inspect `summary(mfrm_threshold_profiles())`.
#' 2. Choose profile (`strict` / `standard` / `lenient`) for project policy.
#' 3. Override selected thresholds in [build_visual_summaries()] only when justified.
#'
#' @return An object of class `summary.mfrm_threshold_profiles`.
#' @seealso [mfrm_threshold_profiles()], [build_visual_summaries()]
#' @examples
#' profiles <- mfrm_threshold_profiles()
#' summary(profiles)
#' @export
summary.mfrm_threshold_profiles <- function(object, digits = 3, ...) {
  if (!inherits(object, "mfrm_threshold_profiles")) {
    stop("`object` must be an mfrm_threshold_profiles object from mfrm_threshold_profiles().", call. = FALSE)
  }
  digits <- max(0L, as.integer(digits))

  profiles <- object$profiles %||% list()
  profile_names <- names(profiles)
  if (is.null(profile_names)) profile_names <- character(0)

  threshold_names <- sort(unique(unlist(lapply(profiles, names), use.names = FALSE)))
  thresholds_tbl <- if (length(threshold_names) == 0) {
    data.frame()
  } else {
    tbl <- data.frame(Threshold = threshold_names, stringsAsFactors = FALSE)
    for (nm in profile_names) {
      vals <- vapply(
        threshold_names,
        function(key) {
          val <- profiles[[nm]][[key]]
          val <- suppressWarnings(as.numeric(val))
          ifelse(length(val) == 0, NA_real_, val[1])
        },
        numeric(1)
      )
      tbl[[nm]] <- vals
    }
    tbl
  }

  thresholds_range_tbl <- data.frame()
  if (nrow(thresholds_tbl) > 0 && length(profile_names) > 0) {
    mat <- as.matrix(thresholds_tbl[, profile_names, drop = FALSE])
    suppressWarnings(storage.mode(mat) <- "numeric")
    row_stats <- t(apply(mat, 1, function(v) {
      vv <- suppressWarnings(as.numeric(v))
      vv <- vv[is.finite(vv)]
      if (length(vv) == 0) return(c(Min = NA_real_, Median = NA_real_, Max = NA_real_, Span = NA_real_))
      c(
        Min = min(vv),
        Median = stats::median(vv),
        Max = max(vv),
        Span = max(vv) - min(vv)
      )
    }))
    thresholds_range_tbl <- data.frame(
      Threshold = thresholds_tbl$Threshold,
      row_stats,
      stringsAsFactors = FALSE
    )
  }

  band_tbl <- data.frame()
  bands <- object$pca_reference_bands %||% list()
  if (length(bands) > 0) {
    band_rows <- lapply(names(bands), function(band_name) {
      vals <- bands[[band_name]]
      if (is.null(vals) || length(vals) == 0) return(NULL)
      keys <- names(vals)
      if (is.null(keys) || length(keys) != length(vals)) {
        keys <- paste0("value_", seq_along(vals))
      }
      data.frame(
        Band = band_name,
        Key = as.character(keys),
        Value = suppressWarnings(as.numeric(vals)),
        stringsAsFactors = FALSE
      )
    })
    band_rows <- Filter(Negate(is.null), band_rows)
    if (length(band_rows) > 0) {
      band_tbl <- do.call(rbind, band_rows)
    }
  }

  overview <- data.frame(
    Profiles = length(profile_names),
    ThresholdCount = nrow(thresholds_tbl),
    PCAReferenceCount = nrow(band_tbl),
    DefaultProfile = if ("standard" %in% profile_names) "standard" else ifelse(length(profile_names) > 0, profile_names[1], ""),
    stringsAsFactors = FALSE
  )

  notes <- c(
    "Profiles tune warning strictness for build_visual_summaries().",
    "Use `thresholds` in build_visual_summaries() to override selected values."
  )
  required_profiles <- c("strict", "standard", "lenient")
  missing_profiles <- setdiff(required_profiles, profile_names)
  if (length(missing_profiles) > 0) {
    notes <- c(notes, paste0("Missing presets: ", paste(missing_profiles, collapse = ", "), "."))
  }

  out <- list(
    overview = overview,
    thresholds = thresholds_tbl,
    threshold_ranges = thresholds_range_tbl,
    pca_reference = band_tbl,
    notes = notes,
    digits = digits
  )
  class(out) <- "summary.mfrm_threshold_profiles"
  out
}

#' @export
print.summary.mfrm_threshold_profiles <- function(x, ...) {
  digits <- as.integer(x$digits %||% 3L)
  if (!is.finite(digits)) digits <- 3L

  cat("mfrmr Threshold Profile Summary\n")

  if (!is.null(x$overview) && nrow(x$overview) > 0) {
    cat("\nOverview\n")
    print(as.data.frame(x$overview), row.names = FALSE)
  }
  if (!is.null(x$thresholds) && nrow(x$thresholds) > 0) {
    cat("\nProfile thresholds\n")
    print(round_numeric_df(as.data.frame(x$thresholds), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$threshold_ranges) && nrow(x$threshold_ranges) > 0) {
    cat("\nThreshold ranges across profiles\n")
    print(round_numeric_df(as.data.frame(x$threshold_ranges), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$pca_reference) && nrow(x$pca_reference) > 0) {
    cat("\nPCA reference bands\n")
    print(round_numeric_df(as.data.frame(x$pca_reference), digits = digits), row.names = FALSE)
  }
  if (length(x$notes) > 0) {
    cat("\nNotes\n")
    cat(" - ", x$notes, "\n", sep = "")
  }
  invisible(x)
}

#' Build warning and narrative summaries for visual outputs
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Output from [diagnose_mfrm()].
#' @param threshold_profile Threshold profile name (`strict`, `standard`, `lenient`).
#' @param thresholds Optional named overrides for profile thresholds.
#' @param summary_options Summary options for `build_visual_summary_map()`.
#' @param whexact Use exact ZSTD transformation.
#' @param branch Output branch:
#'   `"facets"` adds FACETS crosswalk metadata for manual-aligned reporting;
#'   `"original"` keeps package-native summary output.
#'
#' @details
#' This function returns visual-keyed text maps
#' to support dashboard/report rendering without hard-coding narrative strings
#' in UI code.
#'
#' `thresholds` can override any profile field by name. Common overrides:
#' - `n_obs_min`, `n_person_min`
#' - `misfit_ratio_warn`, `zstd2_ratio_warn`, `zstd3_ratio_warn`
#' - `pca_first_eigen_warn`, `pca_first_prop_warn`
#'
#' `summary_options` supports:
#' - `detail`: `"standard"` or `"detailed"`
#' - `max_facet_ranges`: max facet-range snippets shown in visual summaries
#' - `top_misfit_n`: number of top misfit entries included
#'
#' @section Interpreting output:
#' - `warning_map`: rule-triggered warning text by visual key.
#' - `summary_map`: descriptive narrative text by visual key.
#' - `warning_counts` / `summary_counts`: message-count tables for QA checks.
#'
#' @section Typical workflow:
#' 1. inspect defaults with [mfrm_threshold_profiles()]
#' 2. choose `threshold_profile` (`strict` / `standard` / `lenient`)
#' 3. optionally override selected fields via `thresholds`
#' 4. pass result maps to report/dashboard rendering logic
#'
#' @return
#' An object of class `mfrm_visual_summaries` with:
#' - `warning_map`: visual-level warning text vectors
#' - `summary_map`: visual-level descriptive text vectors
#' - `warning_counts`, `summary_counts`: message counts by visual key
#' - `crosswalk`: FACETS-reference mapping for main visual keys
#' - `branch`, `style`, `threshold_profile`: branch metadata
#'
#' @seealso [mfrm_threshold_profiles()], [build_apa_outputs()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' diag <- diagnose_mfrm(fit, residual_pca = "both")
#' vis <- build_visual_summaries(fit, diag, threshold_profile = "strict")
#' vis2 <- build_visual_summaries(
#'   fit,
#'   diag,
#'   threshold_profile = "standard",
#'   thresholds = c(misfit_ratio_warn = 0.20, pca_first_eigen_warn = 2.0),
#'   summary_options = list(detail = "detailed", top_misfit_n = 5)
#' )
#' vis_facets <- build_visual_summaries(fit, diag, branch = "facets")
#' vis_facets$branch
#' summary(vis)
#' p <- plot(vis, type = "comparison", draw = FALSE)
#' p2 <- plot(vis, type = "warning_counts", draw = FALSE)
#' if (interactive()) {
#'   plot(
#'     vis,
#'     type = "comparison",
#'     draw = TRUE,
#'     main = "Warning vs Summary Counts (Customized)",
#'     palette = c(warning = "#cb181d", summary = "#3182bd"),
#'     label_angle = 45
#'   )
#' }
#' @export
build_visual_summaries <- function(fit,
                                   diagnostics,
                                   threshold_profile = "standard",
                                   thresholds = NULL,
                                   summary_options = NULL,
                                   whexact = FALSE,
                                   branch = c("original", "facets")) {
  branch <- match.arg(tolower(as.character(branch[1])), c("original", "facets"))
  style <- ifelse(branch == "facets", "facets_manual", "original")

  warning_map <- build_visual_warning_map(
    res = fit,
    diagnostics = diagnostics,
    whexact = whexact,
    thresholds = thresholds,
    threshold_profile = threshold_profile
  )
  summary_map <- build_visual_summary_map(
    res = fit,
    diagnostics = diagnostics,
    whexact = whexact,
    options = summary_options,
    thresholds = thresholds,
    threshold_profile = threshold_profile
  )

  count_map_messages <- function(x) {
    if (is.null(x) || length(x) == 0) return(0L)
    vals <- unlist(x, use.names = FALSE)
    vals <- trimws(as.character(vals))
    sum(nzchar(vals))
  }
  to_count_table <- function(x) {
    keys <- names(x)
    if (is.null(keys) || length(keys) == 0) {
      return(tibble::tibble(Visual = character(0), Messages = integer(0)))
    }
    tibble::tibble(
      Visual = keys,
      Messages = vapply(x, count_map_messages, integer(1))
    ) |>
      dplyr::arrange(dplyr::desc(.data$Messages), .data$Visual)
  }

  crosswalk <- tibble::tibble(
    Visual = c(
      "unexpected",
      "fair_average",
      "displacement",
      "interrater",
      "facets_chisq",
      "residual_pca_overall",
      "residual_pca_by_facet"
    ),
    FACETS = c(
      "Table 4 / Table 10",
      "Table 12",
      "Table 9",
      "Inter-rater outputs",
      "Facet fixed/random chi-square",
      "Residual PCA (overall)",
      "Residual PCA (by facet)"
    )
  )

  out <- list(
    warning_map = warning_map,
    summary_map = summary_map,
    warning_counts = to_count_table(warning_map),
    summary_counts = to_count_table(summary_map),
    crosswalk = crosswalk,
    branch = branch,
    style = style,
    threshold_profile = as.character(threshold_profile[1])
  )
  out <- as_mfrm_bundle(out, "mfrm_visual_summaries")
  class(out) <- unique(c(paste0("mfrm_visual_summaries_", branch), class(out)))
  out
}

resolve_facets_contract_path <- function(contract_file = NULL) {
  if (!is.null(contract_file)) {
    path <- as.character(contract_file[1])
    if (file.exists(path)) return(path)
    stop("`contract_file` does not exist: ", path)
  }

  installed <- system.file("references", "facets_column_contract.csv", package = "mfrmr")
  if (nzchar(installed) && file.exists(installed)) return(installed)

  source_path <- file.path("inst", "references", "facets_column_contract.csv")
  if (file.exists(source_path)) return(source_path)

  stop(
    "Could not locate `facets_column_contract.csv`.\n",
    "Set `contract_file` explicitly or ensure the package was installed with `inst/references`."
  )
}

read_facets_contract <- function(contract_file = NULL, branch = c("facets", "original")) {
  branch <- match.arg(tolower(as.character(branch[1])), c("facets", "original"))
  path <- resolve_facets_contract_path(contract_file)
  contract <- utils::read.csv(path, stringsAsFactors = FALSE)
  need <- c("table_id", "function_name", "object_id", "component", "required_columns")
  if (!all(need %in% names(contract))) {
    stop("FACETS contract file is missing required columns: ", paste(setdiff(need, names(contract)), collapse = ", "))
  }

  # Original branch uses compact Table 11 column names.
  if (identical(branch, "original")) {
    idx <- contract$object_id == "t11" & contract$component == "table"
    contract$required_columns[idx] <- "Count|BiasSize|SE|LowCountFlag"
  }

  list(path = path, contract = contract)
}

split_contract_tokens <- function(required_columns) {
  vals <- strsplit(as.character(required_columns[1]), "|", fixed = TRUE)[[1]]
  vals <- trimws(vals)
  vals[nzchar(vals)]
}

contract_token_present <- function(token, columns) {
  token <- as.character(token[1])
  if (!nzchar(token)) return(TRUE)
  if (endsWith(token, "*")) {
    prefix <- substr(token, 1L, nchar(token) - 1L)
    return(any(startsWith(columns, prefix)))
  }
  token %in% columns
}

make_metric_row <- function(table_id, check, pass, actual = NA_real_, expected = NA_real_, note = "") {
  data.frame(
    Table = as.character(table_id),
    Check = as.character(check),
    Pass = if (is.na(pass)) NA else as.logical(pass),
    Actual = as.character(actual),
    Expected = as.character(expected),
    Note = as.character(note),
    stringsAsFactors = FALSE
  )
}

safe_num <- function(x) suppressWarnings(as.numeric(x))

build_parity_metric_audit <- function(outputs, tol = 1e-8) {
  rows <- list()

  add_row <- function(table_id, check, pass, actual = NA_real_, expected = NA_real_, note = "") {
    rows[[length(rows) + 1L]] <<- make_metric_row(table_id, check, pass, actual, expected, note)
  }

  t4 <- outputs$t4
  if (!is.null(t4) && is.data.frame(t4$summary) && nrow(t4$summary) > 0) {
    s4 <- t4$summary[1, , drop = FALSE]
    total <- safe_num(s4$TotalObservations)
    unexpected_n <- safe_num(s4$UnexpectedN)
    pct <- safe_num(s4$UnexpectedPercent)
    calc <- if (is.finite(total) && total > 0) 100 * unexpected_n / total else NA_real_
    pass <- if (is.finite(calc) && is.finite(pct)) abs(calc - pct) <= 1e-6 else NA
    add_row("T4", "UnexpectedPercent consistency", pass, pct, calc)
  }

  t10 <- outputs$t10
  if (!is.null(t10) && is.data.frame(t10$summary) && nrow(t10$summary) > 0) {
    s10 <- t10$summary[1, , drop = FALSE]
    baseline <- safe_num(s10$BaselineUnexpectedN)
    after <- safe_num(s10$AfterBiasUnexpectedN)
    reduced <- safe_num(s10$ReducedBy)
    reduced_pct <- safe_num(s10$ReducedPercent)
    calc_reduced <- if (all(is.finite(c(baseline, after)))) baseline - after else NA_real_
    calc_pct <- if (is.finite(baseline) && baseline > 0 && is.finite(reduced)) 100 * reduced / baseline else NA_real_
    pass_reduced <- if (is.finite(calc_reduced) && is.finite(reduced)) abs(calc_reduced - reduced) <= tol else NA
    pass_pct <- if (is.finite(calc_pct) && is.finite(reduced_pct)) abs(calc_pct - reduced_pct) <= 1e-6 else NA
    add_row("T10", "ReducedBy consistency", pass_reduced, reduced, calc_reduced)
    add_row("T10", "ReducedPercent consistency", pass_pct, reduced_pct, calc_pct)
  }

  t11 <- outputs$t11
  if (!is.null(t11) && is.data.frame(t11$summary) && nrow(t11$summary) > 0) {
    s11 <- t11$summary[1, , drop = FALSE]
    cells <- safe_num(s11$Cells)
    low <- safe_num(s11$LowCountCells)
    low_pct <- safe_num(s11$LowCountPercent)
    calc <- if (is.finite(cells) && cells > 0 && is.finite(low)) 100 * low / cells else NA_real_
    pass <- if (is.finite(calc) && is.finite(low_pct)) abs(calc - low_pct) <= 1e-6 else NA
    add_row("T11", "LowCountPercent consistency", pass, low_pct, calc)
  }

  t7a <- outputs$t7agree
  if (!is.null(t7a) && is.data.frame(t7a$summary) && nrow(t7a$summary) > 0) {
    s <- t7a$summary[1, , drop = FALSE]
    exact <- safe_num(s$ExactAgreement)
    expected_exact <- safe_num(s$ExpectedExactAgreement)
    adjacent <- safe_num(s$AdjacentAgreement)
    in_range <- function(v) is.finite(v) && v >= -tol && v <= 1 + tol
    add_row("T7", "ExactAgreement range", in_range(exact), exact, "[0,1]")
    add_row("T7", "ExpectedExactAgreement range", in_range(expected_exact), expected_exact, "[0,1]")
    add_row("T7", "AdjacentAgreement range", in_range(adjacent), adjacent, "[0,1]")
  }

  t7c <- outputs$t7chisq
  if (!is.null(t7c) && is.data.frame(t7c$table) && nrow(t7c$table) > 0) {
    fp <- safe_num(t7c$table$FixedProb)
    rp <- safe_num(t7c$table$RandomProb)
    in_unit <- function(v) {
      vals <- v[is.finite(v)]
      if (length(vals) == 0) return(NA)
      all(vals >= -tol & vals <= 1 + tol)
    }
    add_row("T7", "FixedProb range", in_unit(fp), "all", "[0,1]")
    add_row("T7", "RandomProb range", in_unit(rp), "all", "[0,1]")
  }

  disp <- outputs$disp
  if (!is.null(disp) && is.data.frame(disp$summary) && nrow(disp$summary) > 0) {
    s <- disp$summary[1, , drop = FALSE]
    levels_n <- safe_num(s$Levels)
    anchored <- safe_num(s$AnchoredLevels)
    flagged <- safe_num(s$FlaggedLevels)
    flagged_anch <- safe_num(s$FlaggedAnchoredLevels)
    pass1 <- if (all(is.finite(c(levels_n, anchored)))) anchored <= levels_n + tol else NA
    pass2 <- if (all(is.finite(c(levels_n, flagged)))) flagged <= levels_n + tol else NA
    pass3 <- if (all(is.finite(c(anchored, flagged_anch)))) flagged_anch <= anchored + tol else NA
    add_row("T9", "AnchoredLevels <= Levels", pass1, anchored, levels_n)
    add_row("T9", "FlaggedLevels <= Levels", pass2, flagged, levels_n)
    add_row("T9", "FlaggedAnchoredLevels <= AnchoredLevels", pass3, flagged_anch, anchored)
  }

  t81 <- outputs$t81
  if (!is.null(t81) && is.data.frame(t81$summary) && nrow(t81$summary) > 0) {
    s <- t81$summary[1, , drop = FALSE]
    cats <- safe_num(s$Categories)
    used <- safe_num(s$UsedCategories)
    pass_used <- if (all(is.finite(c(cats, used)))) used <= cats + tol else NA
    add_row("T8.1", "UsedCategories <= Categories", pass_used, used, cats)

    tt <- t81$threshold_table
    if (is.data.frame(tt) && nrow(tt) > 1 && "GapFromPrev" %in% names(tt)) {
      gaps <- safe_num(tt$GapFromPrev)
      monotonic_calc <- !any(gaps[is.finite(gaps)] < -tol)
      monotonic_flag <- isTRUE(s$ThresholdMonotonic)
      add_row("T8.1", "ThresholdMonotonic consistency", monotonic_flag == monotonic_calc, monotonic_flag, monotonic_calc)
    }
  }

  if (length(rows) == 0) {
    return(data.frame(
      Table = character(0),
      Check = character(0),
      Pass = logical(0),
      Actual = character(0),
      Expected = character(0),
      Note = character(0),
      stringsAsFactors = FALSE
    ))
  }
  dplyr::bind_rows(rows)
}

#' Build a FACETS parity report (column + metric contracts)
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()]. If omitted,
#'   diagnostics are computed internally with `residual_pca = "none"`.
#' @param bias_results Optional output from [estimate_bias()]. If omitted and
#'   at least two facets exist, a 2-way bias run is computed internally.
#' @param branch Contract branch. `"facets"` checks FACETS-style columns.
#'   `"original"` adapts branch-sensitive contracts (currently Table 11) to the
#'   package's compact naming.
#' @param contract_file Optional path to a custom contract CSV.
#' @param include_metrics If `TRUE`, run additional numerical consistency checks.
#' @param top_n_missing Number of lowest-coverage contract rows to keep in
#'   `missing_preview`.
#'
#' @details
#' This function compares produced report components to a contract specification
#' (`inst/references/facets_column_contract.csv`) and returns:
#' - column-level coverage per contract row
#' - table-level coverage summaries
#' - optional metric-level consistency checks
#'
#' Coverage interpretation in `overall`:
#' - `MeanColumnCoverage` and `MinColumnCoverage` are computed across all
#'   contract rows (unavailable rows count as 0 coverage).
#' - `MeanColumnCoverageAvailable` and `MinColumnCoverageAvailable` summarize
#'   only rows whose source component is available.
#'
#' `summary(out)` is supported through `summary()`.
#' `plot(out)` is dispatched through `plot()` for class
#' `mfrm_parity_report` (`type = "column_coverage"`, `"table_coverage"`,
#' `"metric_status"`, `"metric_by_table"`).
#'
#' @section Interpreting output:
#' - `overall`: high-level contract coverage and metric-check pass rates.
#' - `column_summary` / `column_audit`: where schema mismatches occur.
#' - `metric_summary` / `metric_audit`: numerical consistency checks.
#' - `missing_preview`: quickest path to unresolved parity gaps.
#'
#' @section Typical workflow:
#' 1. Run `facets_parity_report(fit, branch = "facets")`.
#' 2. Inspect `summary(parity)` and `missing_preview`.
#' 3. Patch upstream table builders, then rerun parity report.
#'
#' @return
#' An object of class `mfrm_parity_report` with:
#' - `overall`: one-row overall parity summary
#' - `column_summary`: coverage summary by table ID
#' - `column_audit`: row-level contract audit
#' - `missing_preview`: lowest-coverage rows
#' - `metric_summary`: one-row metric-check summary
#' - `metric_by_table`: metric-check summary by table ID
#' - `metric_audit`: row-level metric checks
#' - `settings`: branch/contract metadata
#'
#' @seealso [fit_mfrm()], [diagnose_mfrm()], [build_fixed_reports()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#' parity <- facets_parity_report(fit, diagnostics = diag, branch = "facets")
#' summary(parity)
#' p <- plot(parity, draw = FALSE)
#' @export
facets_parity_report <- function(fit,
                                 diagnostics = NULL,
                                 bias_results = NULL,
                                 branch = c("facets", "original"),
                                 contract_file = NULL,
                                 include_metrics = TRUE,
                                 top_n_missing = 15L) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().")
  }
  branch <- match.arg(tolower(as.character(branch[1])), c("facets", "original"))
  include_metrics <- isTRUE(include_metrics)
  top_n_missing <- max(1L, as.integer(top_n_missing))

  if (is.null(diagnostics)) {
    diagnostics <- diagnose_mfrm(fit, residual_pca = "none")
  }

  facet_names <- as.character(fit$config$facet_names %||% character(0))
  if (is.null(bias_results) && length(facet_names) >= 2) {
    bias_results <- estimate_bias(
      fit = fit,
      diagnostics = diagnostics,
      facet_a = facet_names[1],
      facet_b = facet_names[2],
      max_iter = 2
    )
  }

  contract_info <- read_facets_contract(contract_file = contract_file, branch = branch)
  contract <- as.data.frame(contract_info$contract, stringsAsFactors = FALSE)

  outputs <- list(
    t1 = specifications_report(fit),
    t2 = data_quality_report(
      fit = fit,
      data = fit$prep$data,
      person = fit$config$person_col,
      facets = fit$config$facet_names,
      score = fit$config$score_col,
      weight = fit$config$weight_col
    ),
    t3 = estimation_iteration_report(fit, max_iter = 5),
    t4 = unexpected_response_table(fit, diagnostics = diagnostics, top_n = 50),
    t5 = measurable_summary_table(fit, diagnostics = diagnostics),
    t6 = subset_connectivity_report(fit, diagnostics = diagnostics),
    t62 = facet_statistics_report(fit, diagnostics = diagnostics),
    t7chisq = facets_chisq_table(fit, diagnostics = diagnostics),
    t7agree = interrater_agreement_table(fit, diagnostics = diagnostics),
    t81 = rating_scale_table(fit, diagnostics = diagnostics),
    t8bar = category_structure_report(fit, diagnostics = diagnostics),
    t8curves = category_curves_report(fit, theta_points = 101),
    out = facets_output_file_bundle(fit, diagnostics = diagnostics, include = c("graph", "score"), theta_points = 81),
    t12 = fair_average_table(fit, diagnostics = diagnostics),
    disp = displacement_table(fit, diagnostics = diagnostics)
  )
  if (!is.null(bias_results) && is.data.frame(bias_results$table) && nrow(bias_results$table) > 0) {
    outputs$t10 <- unexpected_after_bias_table(fit, bias_results, diagnostics = diagnostics, top_n = 50)
    outputs$t11 <- bias_count_table(bias_results, branch = branch)
    outputs$t13 <- bias_interaction_report(bias_results)
    outputs$t14 <- build_fixed_reports(bias_results, branch = branch)
  } else {
    outputs$t10 <- NULL
    outputs$t11 <- NULL
    outputs$t13 <- NULL
    outputs$t14 <- NULL
  }

  audit_rows <- lapply(seq_len(nrow(contract)), function(i) {
    row <- contract[i, , drop = FALSE]
    tokens <- split_contract_tokens(row$required_columns)
    obj <- outputs[[row$object_id]]
    if (is.null(obj)) {
      return(data.frame(
        table_id = row$table_id,
        function_name = row$function_name,
        object_id = row$object_id,
        component = row$component,
        required_n = length(tokens),
        present_n = NA_integer_,
        coverage = NA_real_,
        available = FALSE,
        full_match = FALSE,
        status = "missing_object",
        missing = paste(tokens, collapse = " | "),
        stringsAsFactors = FALSE
      ))
    }
    comp <- obj[[row$component]]
    if (!is.data.frame(comp)) {
      return(data.frame(
        table_id = row$table_id,
        function_name = row$function_name,
        object_id = row$object_id,
        component = row$component,
        required_n = length(tokens),
        present_n = NA_integer_,
        coverage = NA_real_,
        available = FALSE,
        full_match = FALSE,
        status = "missing_component",
        missing = paste(tokens, collapse = " | "),
        stringsAsFactors = FALSE
      ))
    }
    cols <- names(comp)
    present <- vapply(tokens, contract_token_present, logical(1), columns = cols)
    missing <- tokens[!present]
    cov <- if (length(tokens) == 0) 1 else sum(present) / length(tokens)
    data.frame(
      table_id = row$table_id,
      function_name = row$function_name,
      object_id = row$object_id,
      component = row$component,
      required_n = length(tokens),
      present_n = sum(present),
      coverage = cov,
      available = TRUE,
      full_match = isTRUE(all(present)),
      status = if (isTRUE(all(present))) "match" else "partial",
      missing = paste(missing, collapse = " | "),
      stringsAsFactors = FALSE
    )
  })
  column_audit <- dplyr::bind_rows(audit_rows)

  summarize_coverage <- function(v, fn) {
    vals <- suppressWarnings(as.numeric(v))
    vals <- vals[is.finite(vals)]
    if (length(vals) == 0) return(NA_real_)
    fn(vals)
  }

  # Contract-level coverage should treat unavailable rows as zero coverage.
  # This avoids reporting perfect mean/min coverage when some contract rows
  # are entirely missing from available outputs.
  contract_coverage_values <- ifelse(
    column_audit$available %in% TRUE,
    suppressWarnings(as.numeric(column_audit$coverage)),
    0
  )
  contract_coverage_values[!is.finite(contract_coverage_values)] <- 0

  column_summary <- column_audit |>
    dplyr::group_by(.data$table_id, .data$function_name) |>
    dplyr::summarize(
      Components = dplyr::n(),
      Available = sum(.data$available, na.rm = TRUE),
      FullMatch = sum(.data$full_match, na.rm = TRUE),
      MeanCoverage = summarize_coverage(.data$coverage, mean),
      MinCoverage = summarize_coverage(.data$coverage, min),
      .groups = "drop"
    ) |>
    dplyr::arrange(.data$table_id, .data$function_name)

  missing_preview <- column_audit |>
    dplyr::filter(!.data$full_match | !.data$available) |>
    dplyr::arrange(.data$coverage, .data$table_id, .data$component) |>
    dplyr::slice_head(n = top_n_missing)

  metric_audit <- if (isTRUE(include_metrics)) {
    build_parity_metric_audit(outputs = outputs)
  } else {
    data.frame(
      Table = character(0),
      Check = character(0),
      Pass = logical(0),
      Actual = character(0),
      Expected = character(0),
      Note = character(0),
      stringsAsFactors = FALSE
    )
  }

  metric_summary <- if (nrow(metric_audit) == 0) {
    data.frame(
      Checks = 0L,
      Evaluated = 0L,
      Passed = 0L,
      Failed = 0L,
      PassRate = NA_real_,
      stringsAsFactors = FALSE
    )
  } else {
    ev <- metric_audit$Pass[!is.na(metric_audit$Pass)]
    data.frame(
      Checks = nrow(metric_audit),
      Evaluated = length(ev),
      Passed = sum(ev %in% TRUE),
      Failed = sum(ev %in% FALSE),
      PassRate = if (length(ev) > 0) sum(ev %in% TRUE) / length(ev) else NA_real_,
      stringsAsFactors = FALSE
    )
  }

  metric_by_table <- if (nrow(metric_audit) == 0) {
    data.frame(
      Table = character(0),
      Checks = integer(0),
      Evaluated = integer(0),
      Passed = integer(0),
      Failed = integer(0),
      PassRate = numeric(0),
      stringsAsFactors = FALSE
    )
  } else {
    metric_audit |>
      dplyr::group_by(.data$Table) |>
      dplyr::summarize(
        Checks = dplyr::n(),
        Evaluated = sum(!is.na(.data$Pass)),
        Passed = sum(.data$Pass %in% TRUE, na.rm = TRUE),
        Failed = sum(.data$Pass %in% FALSE, na.rm = TRUE),
        PassRate = ifelse(sum(!is.na(.data$Pass)) > 0, sum(.data$Pass %in% TRUE, na.rm = TRUE) / sum(!is.na(.data$Pass)), NA_real_),
        .groups = "drop"
      ) |>
      dplyr::arrange(.data$Table)
  }

  mean_cov_all <- summarize_coverage(contract_coverage_values, mean)
  min_cov_all <- summarize_coverage(contract_coverage_values, min)
  mean_cov_available <- summarize_coverage(column_audit$coverage, mean)
  min_cov_available <- summarize_coverage(column_audit$coverage, min)
  contract_rows <- nrow(column_audit)
  mismatches <- sum(!column_audit$full_match, na.rm = TRUE)
  overall <- data.frame(
    Branch = branch,
    ContractRows = contract_rows,
    AvailableRows = sum(column_audit$available, na.rm = TRUE),
    FullMatchRows = sum(column_audit$full_match, na.rm = TRUE),
    ColumnMismatches = mismatches,
    ColumnMismatchRate = if (contract_rows > 0) mismatches / contract_rows else NA_real_,
    MeanColumnCoverage = mean_cov_all,
    MinColumnCoverage = min_cov_all,
    MeanColumnCoverageAvailable = mean_cov_available,
    MinColumnCoverageAvailable = min_cov_available,
    MetricChecks = metric_summary$Checks[1],
    MetricEvaluated = metric_summary$Evaluated[1],
    MetricFailed = metric_summary$Failed[1],
    MetricPassRate = metric_summary$PassRate[1],
    stringsAsFactors = FALSE
  )

  out <- list(
    overall = overall,
    column_summary = as.data.frame(column_summary, stringsAsFactors = FALSE),
    column_audit = as.data.frame(column_audit, stringsAsFactors = FALSE),
    missing_preview = as.data.frame(missing_preview, stringsAsFactors = FALSE),
    metric_summary = metric_summary,
    metric_by_table = as.data.frame(metric_by_table, stringsAsFactors = FALSE),
    metric_audit = as.data.frame(metric_audit, stringsAsFactors = FALSE),
    settings = list(
      branch = branch,
      contract_path = contract_info$path,
      include_metrics = include_metrics,
      top_n_missing = top_n_missing,
      bias_included = !is.null(outputs$t10)
    )
  )
  as_mfrm_bundle(out, "mfrm_parity_report")
}

bundle_settings_table <- function(settings) {
  if (is.null(settings) || !is.list(settings) || length(settings) == 0) return(data.frame())
  keys <- names(settings)
  if (is.null(keys) || any(!nzchar(keys))) {
    keys <- paste0("Setting", seq_along(settings))
  }
  vals <- vapply(settings, function(v) {
    if (is.null(v)) return("NULL")
    if (is.data.frame(v)) return(paste0("<table ", nrow(v), "x", ncol(v), ">"))
    if (is.list(v)) return(paste0("<list ", length(v), ">"))
    paste(as.character(v), collapse = ", ")
  }, character(1))
  data.frame(Setting = keys, Value = vals, stringsAsFactors = FALSE)
}

bundle_preview_table <- function(object, top_n = 10L) {
  keys <- c(
    "table", "pairs", "stacked", "ranked_table", "facet_profile", "graphfile",
    "category_table", "facet_coverage", "listing", "overall_table", "by_facet_table",
    "missing_preview", "column_audit", "metric_audit", "column_summary", "metric_summary"
  )
  nm <- names(object)
  if (is.null(nm) || length(nm) == 0) {
    return(list(name = NA_character_, table = data.frame()))
  }
  key <- keys[keys %in% nm][1]
  if (is.na(key) || length(key) == 0) {
    return(list(name = NA_character_, table = data.frame()))
  }
  tbl <- object[[key]]
  if (!is.data.frame(tbl) || nrow(tbl) == 0) {
    return(list(name = key, table = data.frame()))
  }
  top_n <- max(1L, as.integer(top_n))
  list(name = key, table = utils::head(as.data.frame(tbl, stringsAsFactors = FALSE), n = top_n))
}

summarize_bias_count_bundle <- function(object, digits = 3, top_n = 10) {
  tbl <- as.data.frame(object$table %||% data.frame(), stringsAsFactors = FALSE)
  if ("Observd Count" %in% names(tbl) && !"Count" %in% names(tbl)) {
    tbl$Count <- suppressWarnings(as.numeric(tbl$`Observd Count`))
  }
  if (!"LowCountFlag" %in% names(tbl)) {
    tbl$LowCountFlag <- FALSE
  }
  if (!is.logical(tbl$LowCountFlag)) {
    tbl$LowCountFlag <- as.logical(tbl$LowCountFlag)
  }
  if (!"Count" %in% names(tbl)) {
    tbl$Count <- suppressWarnings(as.numeric(tbl$Count))
  }

  cnt <- suppressWarnings(as.numeric(tbl$Count))
  cnt <- cnt[is.finite(cnt)]
  count_distribution <- if (length(cnt) == 0) {
    data.frame()
  } else {
    data.frame(
      Min = min(cnt),
      Q1 = stats::quantile(cnt, 0.25, names = FALSE),
      Median = stats::median(cnt),
      Mean = mean(cnt),
      Q3 = stats::quantile(cnt, 0.75, names = FALSE),
      Max = max(cnt),
      stringsAsFactors = FALSE
    )
  }

  low_tbl <- tbl[tbl$LowCountFlag %in% TRUE, , drop = FALSE]
  if (nrow(low_tbl) > 0 && "Count" %in% names(low_tbl)) {
    low_tbl <- low_tbl |>
      dplyr::arrange(.data$Count) |>
      dplyr::slice_head(n = top_n)
  }

  summary_tbl <- as.data.frame(object$summary %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(summary_tbl) == 0) {
    summary_tbl <- data.frame(
      Branch = as.character(object$branch %||% "original"),
      Cells = nrow(tbl),
      LowCountCells = sum(tbl$LowCountFlag %in% TRUE, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }

  out <- list(
    summary_kind = "bias_count",
    overview = summary_tbl,
    count_distribution = count_distribution,
    low_count_cells = low_tbl,
    thresholds = bundle_settings_table(object$thresholds),
    notes = if (identical(object$branch, "facets")) {
      "FACETS branch: table columns mirror FACETS Table 11 naming."
    } else {
      "Original branch: compact count/bias columns for QC screening."
    },
    digits = digits
  )
  class(out) <- "summary.mfrm_bundle"
  out
}

summarize_fixed_reports_bundle <- function(object, digits = 3, top_n = 10) {
  pair_tbl <- as.data.frame(object$pairwise_table %||% data.frame(), stringsAsFactors = FALSE)
  n_bias_lines <- length(strsplit(as.character(object$bias_fixed %||% ""), "\n", fixed = TRUE)[[1]])
  n_pair_lines <- length(strsplit(as.character(object$pairwise_fixed %||% ""), "\n", fixed = TRUE)[[1]])

  overview <- data.frame(
    Branch = as.character(object$branch %||% "facets"),
    Style = as.character(object$style %||% "facets_manual"),
    Interaction = as.character(object$interaction_label %||% ""),
    PairwiseRows = nrow(pair_tbl),
    BiasTextLines = n_bias_lines,
    PairwiseTextLines = n_pair_lines,
    stringsAsFactors = FALSE
  )

  out <- list(
    summary_kind = "fixed_reports",
    overview = overview,
    summary = data.frame(),
    preview_name = if (nrow(pair_tbl) > 0) "pairwise_table" else "",
    preview = utils::head(pair_tbl, n = top_n),
    settings = data.frame(),
    notes = if (nrow(pair_tbl) == 0) {
      "No pairwise contrasts available in this interaction mode."
    } else if (identical(object$branch, "facets")) {
      "FACETS branch: fixed-width text follows FACETS-like layout."
    } else {
      "Original branch: sectioned fixed-width text optimized for quick review."
    },
    digits = digits
  )
  class(out) <- "summary.mfrm_bundle"
  out
}

summarize_visual_summaries_bundle <- function(object, digits = 3, top_n = 10) {
  warning_counts <- as.data.frame(object$warning_counts %||% data.frame(), stringsAsFactors = FALSE)
  summary_counts <- as.data.frame(object$summary_counts %||% data.frame(), stringsAsFactors = FALSE)
  crosswalk <- as.data.frame(object$crosswalk %||% data.frame(), stringsAsFactors = FALSE)

  overview <- data.frame(
    Branch = as.character(object$branch %||% "original"),
    Style = as.character(object$style %||% "original"),
    ThresholdProfile = as.character(object$threshold_profile %||% ""),
    WarningVisuals = nrow(warning_counts),
    SummaryVisuals = nrow(summary_counts),
    stringsAsFactors = FALSE
  )

  preview_tbl <- warning_counts
  if (nrow(preview_tbl) == 0) preview_tbl <- summary_counts
  preview_tbl <- utils::head(preview_tbl, n = top_n)

  notes <- if (identical(object$branch, "facets")) {
    "FACETS branch includes crosswalk metadata to manual-oriented output names."
  } else {
    "Original branch keeps package-native warning/summary map organization."
  }

  out <- list(
    summary_kind = "visual_summaries",
    overview = overview,
    summary = warning_counts,
    preview_name = if (nrow(preview_tbl) > 0) "warning_counts" else "",
    preview = preview_tbl,
    settings = crosswalk,
    notes = notes,
    digits = digits,
    summary_counts = summary_counts
  )
  class(out) <- "summary.mfrm_bundle"
  out
}

bundle_component_table <- function(object, name) {
  if (!is.list(object) || is.null(name) || !nzchar(name) || !name %in% names(object)) {
    return(data.frame())
  }
  value <- object[[name]]
  if (!is.data.frame(value)) return(data.frame())
  as.data.frame(value, stringsAsFactors = FALSE)
}

bundle_first_table <- function(object, candidates, top_n = 10L) {
  candidates <- as.character(candidates %||% character(0))
  if (length(candidates) == 0) {
    return(list(name = NA_character_, table = data.frame()))
  }
  for (nm in candidates) {
    tbl <- bundle_component_table(object, nm)
    if (nrow(tbl) > 0) {
      return(list(name = nm, table = utils::head(tbl, n = top_n)))
    }
  }
  for (nm in candidates) {
    tbl <- bundle_component_table(object, nm)
    if (ncol(tbl) > 0) {
      return(list(name = nm, table = tbl))
    }
  }
  list(name = NA_character_, table = data.frame())
}

bundle_known_overview <- function(object, obj_class, preview_name, preview_rows) {
  comp_names <- names(object)
  if (is.null(comp_names)) comp_names <- character(0)
  data.frame(
    Class = obj_class,
    Components = length(comp_names),
    ComponentNames = if (length(comp_names) == 0) "" else paste(comp_names, collapse = ", "),
    PreviewComponent = ifelse(is.na(preview_name), "", preview_name),
    PreviewRows = as.integer(preview_rows),
    stringsAsFactors = FALSE
  )
}

summarize_known_bundle <- function(object,
                                   obj_class,
                                   summary_candidates = "summary",
                                   preview_candidates = NULL,
                                   settings_candidates = "settings",
                                   notes = NULL,
                                   digits = 3,
                                   top_n = 10,
                                   summary_override = NULL) {
  top_n <- max(1L, as.integer(top_n))

  summary_tbl <- if (!is.null(summary_override)) {
    as.data.frame(summary_override, stringsAsFactors = FALSE)
  } else {
    data.frame()
  }
  if (nrow(summary_tbl) == 0 && ncol(summary_tbl) == 0) {
    summary_pick <- bundle_first_table(object, summary_candidates, top_n = 1L)
    summary_tbl <- summary_pick$table
  }

  preview_pick <- bundle_first_table(object, preview_candidates, top_n = top_n)
  if (is.na(preview_pick$name) || nrow(preview_pick$table) == 0) {
    preview_pick <- bundle_preview_table(object, top_n = top_n)
  }

  settings_tbl <- data.frame()
  for (nm in as.character(settings_candidates %||% character(0))) {
    if (!nm %in% names(object)) next
    value <- object[[nm]]
    if (is.data.frame(value)) {
      settings_tbl <- as.data.frame(value, stringsAsFactors = FALSE)
      break
    }
    if (is.list(value)) {
      settings_tbl <- bundle_settings_table(value)
      break
    }
  }

  notes <- as.character(notes %||% "")
  notes <- notes[nzchar(notes)]
  if (length(notes) == 0) {
    if (nrow(summary_tbl) > 0 && nrow(preview_pick$table) > 0) {
      notes <- "Summary and preview tables were extracted for this bundle."
    } else if (nrow(preview_pick$table) > 0) {
      notes <- "Preview rows were extracted from the main table component."
    } else {
      notes <- "No populated table components were found in this bundle."
    }
  }

  out <- list(
    summary_kind = obj_class,
    overview = bundle_known_overview(
      object = object,
      obj_class = obj_class,
      preview_name = preview_pick$name,
      preview_rows = nrow(preview_pick$table)
    ),
    summary = summary_tbl,
    preview_name = preview_pick$name,
    preview = preview_pick$table,
    settings = settings_tbl,
    notes = notes,
    digits = digits
  )
  class(out) <- "summary.mfrm_bundle"
  out
}

summarize_measurable_bundle <- function(object, digits = 3, top_n = 10) {
  summarize_known_bundle(
    object = object,
    obj_class = "mfrm_measurable",
    summary_candidates = "summary",
    preview_candidates = c("facet_coverage", "category_stats", "subsets"),
    settings_candidates = character(0),
    notes = "FACETS Table 5-style measurable-data summary with facet coverage and category diagnostics.",
    digits = digits,
    top_n = top_n
  )
}

summarize_unexpected_after_bias_bundle <- function(object, digits = 3, top_n = 10) {
  facet_note <- if (!is.null(object$facets) && length(object$facets) > 0) {
    paste("Bias interaction:", paste(as.character(object$facets), collapse = " x "))
  } else {
    "Bias interaction facets are not attached in this object."
  }
  summarize_known_bundle(
    object = object,
    obj_class = "mfrm_unexpected_after_bias",
    summary_candidates = "summary",
    preview_candidates = "table",
    settings_candidates = "thresholds",
    notes = c(
      "FACETS Table 10-style unexpected-response summary after bias adjustment.",
      facet_note
    ),
    digits = digits,
    top_n = top_n
  )
}

summarize_output_bundle <- function(object, digits = 3, top_n = 10) {
  settings <- object$settings %||% list()
  summary_tbl <- data.frame(
    GraphRows = nrow(bundle_component_table(object, "graphfile")),
    ScoreRows = nrow(bundle_component_table(object, "scorefile")),
    WrittenFiles = nrow(bundle_component_table(object, "written_files")),
    IncludeFixed = as.logical(settings$include_fixed %||% FALSE),
    WriteFiles = as.logical(settings$write_files %||% FALSE),
    stringsAsFactors = FALSE
  )
  summarize_known_bundle(
    object = object,
    obj_class = "mfrm_output_bundle",
    summary_candidates = character(0),
    preview_candidates = c("scorefile", "graphfile", "graphfile_syntactic", "written_files"),
    settings_candidates = "settings",
    notes = "Graphfile/SCORE-style export bundle (table output and optional file-write metadata).",
    digits = digits,
    top_n = top_n,
    summary_override = summary_tbl
  )
}

summarize_residual_pca_bundle <- function(object, digits = 3, top_n = 10) {
  mode <- as.character(object$mode %||% "unknown")
  facet_names <- as.character(object$facet_names %||% character(0))
  summary_tbl <- data.frame(
    Mode = mode,
    Facets = length(facet_names),
    OverallComponents = nrow(bundle_component_table(object, "overall_table")),
    FacetComponentRows = nrow(bundle_component_table(object, "by_facet_table")),
    stringsAsFactors = FALSE
  )
  summarize_known_bundle(
    object = object,
    obj_class = "mfrm_residual_pca",
    summary_candidates = character(0),
    preview_candidates = c("overall_table", "by_facet_table"),
    settings_candidates = character(0),
    notes = "Residual PCA summary for unidimensionality checks (overall and/or by facet).",
    digits = digits,
    top_n = top_n,
    summary_override = summary_tbl
  )
}

summarize_specifications_bundle <- function(object, digits = 3, top_n = 10) {
  summarize_known_bundle(
    object = object,
    obj_class = "mfrm_specifications",
    summary_candidates = "header",
    preview_candidates = c("data_spec", "facet_labels", "output_spec", "convergence_control", "anchor_summary"),
    settings_candidates = character(0),
    notes = "FACETS Table 1-style model specification summary.",
    digits = digits,
    top_n = top_n
  )
}

summarize_data_quality_bundle <- function(object, digits = 3, top_n = 10) {
  summarize_known_bundle(
    object = object,
    obj_class = "mfrm_data_quality",
    summary_candidates = "summary",
    preview_candidates = c("row_audit", "category_counts", "model_match", "unknown_elements"),
    settings_candidates = character(0),
    notes = "FACETS Table 2-style data quality summary and row-level audit.",
    digits = digits,
    top_n = top_n
  )
}

summarize_iteration_report_bundle <- function(object, digits = 3, top_n = 10) {
  summarize_known_bundle(
    object = object,
    obj_class = "mfrm_iteration_report",
    summary_candidates = "summary",
    preview_candidates = "table",
    settings_candidates = "settings",
    notes = "FACETS Table 3-style replay of estimation iterations.",
    digits = digits,
    top_n = top_n
  )
}

summarize_subset_connectivity_bundle <- function(object, digits = 3, top_n = 10) {
  summarize_known_bundle(
    object = object,
    obj_class = "mfrm_subset_connectivity",
    summary_candidates = "summary",
    preview_candidates = c("listing", "nodes"),
    settings_candidates = "settings",
    notes = "FACETS Table 6 subset/connectivity report with subset and node listings.",
    digits = digits,
    top_n = top_n
  )
}

summarize_facet_statistics_bundle <- function(object, digits = 3, top_n = 10) {
  table_tbl <- bundle_component_table(object, "table")
  range_tbl <- bundle_component_table(object, "ranges")
  summary_tbl <- data.frame(
    Facets = if ("Facet" %in% names(table_tbl)) length(unique(table_tbl$Facet)) else NA_integer_,
    Rows = nrow(table_tbl),
    Metrics = if ("Metric" %in% names(table_tbl)) length(unique(table_tbl$Metric)) else NA_integer_,
    stringsAsFactors = FALSE
  )
  summarize_known_bundle(
    object = object,
    obj_class = "mfrm_facet_statistics",
    summary_candidates = character(0),
    preview_candidates = c("table", "ranges"),
    settings_candidates = "settings",
    notes = if (nrow(range_tbl) > 0) {
      "FACETS Table 6.2-style facet statistics including range summaries."
    } else {
      "FACETS Table 6.2-style facet statistics summary."
    },
    digits = digits,
    top_n = top_n,
    summary_override = summary_tbl
  )
}

summarize_parity_bundle <- function(object, digits = 3, top_n = 10) {
  overall_tbl <- as.data.frame(object$overall %||% data.frame(), stringsAsFactors = FALSE)
  missing_tbl <- as.data.frame(object$missing_preview %||% data.frame(), stringsAsFactors = FALSE)
  metric_summary <- as.data.frame(object$metric_summary %||% data.frame(), stringsAsFactors = FALSE)

  notes <- character(0)
  if (nrow(overall_tbl) > 0) {
    mismatch <- suppressWarnings(as.integer(overall_tbl$ColumnMismatches[1]))
    if (is.finite(mismatch) && mismatch == 0) {
      notes <- c(notes, "All contract rows reached full column coverage.")
    } else if (is.finite(mismatch)) {
      notes <- c(notes, paste0("Column mismatches detected: ", mismatch, "."))
    }
  }
  if (nrow(metric_summary) > 0) {
    failed <- suppressWarnings(as.integer(metric_summary$Failed[1]))
    if (is.finite(failed) && failed == 0) {
      notes <- c(notes, "All evaluated metric checks passed.")
    } else if (is.finite(failed)) {
      notes <- c(notes, paste0("Metric checks failed: ", failed, "."))
    }
  }
  if (length(notes) == 0) {
    notes <- "Parity checks completed."
  }

  summarize_known_bundle(
    object = object,
    obj_class = "mfrm_parity_report",
    summary_candidates = character(0),
    preview_candidates = c("missing_preview", "column_audit", "metric_audit"),
    settings_candidates = "settings",
    notes = notes,
    digits = digits,
    top_n = top_n,
    summary_override = overall_tbl
  )
}

summarize_unexpected_bundle <- function(object, digits = 3, top_n = 10) {
  summarize_known_bundle(
    object = object,
    obj_class = "mfrm_unexpected",
    summary_candidates = "summary",
    preview_candidates = "table",
    settings_candidates = "thresholds",
    notes = "FACETS Table 4-style unexpected-response summary.",
    digits = digits,
    top_n = top_n
  )
}

summarize_fair_average_bundle <- function(object, digits = 3, top_n = 10) {
  stacked <- bundle_component_table(object, "stacked")
  obs_avg <- if ("Obsvd Average" %in% names(stacked)) suppressWarnings(as.numeric(stacked[["Obsvd Average"]])) else numeric(0)
  fair_m <- if ("Fair(M) Average" %in% names(stacked)) suppressWarnings(as.numeric(stacked[["Fair(M) Average"]])) else numeric(0)
  mean_abs_gap <- NA_real_
  if (length(obs_avg) == length(fair_m) && length(obs_avg) > 0) {
    dif <- abs(obs_avg - fair_m)
    dif <- dif[is.finite(dif)]
    if (length(dif) > 0) mean_abs_gap <- mean(dif)
  }
  summary_tbl <- data.frame(
    Facets = if ("Facet" %in% names(stacked)) length(unique(as.character(stacked$Facet))) else length(object$by_facet %||% list()),
    Levels = nrow(stacked),
    MeanAbsObservedFairM = mean_abs_gap,
    stringsAsFactors = FALSE
  )
  summarize_known_bundle(
    object = object,
    obj_class = "mfrm_fair_average",
    summary_candidates = character(0),
    preview_candidates = c("stacked", "raw_by_facet"),
    settings_candidates = "settings",
    notes = "FACETS Table 12-style fair-average comparison by facet level.",
    digits = digits,
    top_n = top_n,
    summary_override = summary_tbl
  )
}

summarize_displacement_bundle <- function(object, digits = 3, top_n = 10) {
  summarize_known_bundle(
    object = object,
    obj_class = "mfrm_displacement",
    summary_candidates = "summary",
    preview_candidates = "table",
    settings_candidates = "thresholds",
    notes = "FACETS-style displacement diagnostics for anchor drift checks.",
    digits = digits,
    top_n = top_n
  )
}

summarize_interrater_bundle <- function(object, digits = 3, top_n = 10) {
  summarize_known_bundle(
    object = object,
    obj_class = "mfrm_interrater",
    summary_candidates = "summary",
    preview_candidates = "pairs",
    settings_candidates = "settings",
    notes = "FACETS Table 7 agreement-style inter-rater summary.",
    digits = digits,
    top_n = top_n
  )
}

summarize_facets_chisq_bundle <- function(object, digits = 3, top_n = 10) {
  summarize_known_bundle(
    object = object,
    obj_class = "mfrm_facets_chisq",
    summary_candidates = "summary",
    preview_candidates = "table",
    settings_candidates = "thresholds",
    notes = "FACETS Table 7 summary-statistics style facet chi-square report.",
    digits = digits,
    top_n = top_n
  )
}

summarize_bias_interaction_bundle <- function(object, digits = 3, top_n = 10) {
  summarize_known_bundle(
    object = object,
    obj_class = "mfrm_bias_interaction",
    summary_candidates = "summary",
    preview_candidates = c("ranked_table", "facet_profile"),
    settings_candidates = "thresholds",
    notes = "FACETS Table 13/14-style bias interaction report.",
    digits = digits,
    top_n = top_n
  )
}

summarize_rating_scale_bundle <- function(object, digits = 3, top_n = 10) {
  summarize_known_bundle(
    object = object,
    obj_class = "mfrm_rating_scale",
    summary_candidates = "summary",
    preview_candidates = c("category_table", "threshold_table"),
    settings_candidates = character(0),
    notes = "FACETS Table 8.1-style rating scale diagnostics.",
    digits = digits,
    top_n = top_n
  )
}

summarize_category_structure_bundle <- function(object, digits = 3, top_n = 10) {
  cat_tbl <- bundle_component_table(object, "category_table")
  flags <- integer(0)
  for (nm in c("LowCount", "InfitFlag", "OutfitFlag", "ZSTDFlag")) {
    if (nm %in% names(cat_tbl)) {
      v <- as.logical(cat_tbl[[nm]])
      flags <- c(flags, sum(v, na.rm = TRUE))
    }
  }
  summary_tbl <- data.frame(
    Categories = nrow(cat_tbl),
    UsedCategories = if ("Count" %in% names(cat_tbl)) sum(suppressWarnings(as.numeric(cat_tbl$Count)) > 0, na.rm = TRUE) else NA_integer_,
    FlaggedStats = if (length(flags) > 0) sum(flags, na.rm = TRUE) else NA_integer_,
    ModeBoundaries = nrow(bundle_component_table(object, "mode_boundaries")),
    MeanHalfscorePoints = nrow(bundle_component_table(object, "mean_halfscore_points")),
    stringsAsFactors = FALSE
  )
  summarize_known_bundle(
    object = object,
    obj_class = "mfrm_category_structure",
    summary_candidates = character(0),
    preview_candidates = c("category_table", "mode_boundaries", "mean_halfscore_points"),
    settings_candidates = "settings",
    notes = "FACETS Table 8 bar-chart style category structure diagnostics.",
    digits = digits,
    top_n = top_n,
    summary_override = summary_tbl
  )
}

summarize_category_curves_bundle <- function(object, digits = 3, top_n = 10) {
  graph_tbl <- bundle_component_table(object, "graphfile")
  prob_cols <- grep("^Prob:", names(graph_tbl), value = TRUE)
  summary_tbl <- data.frame(
    Rows = nrow(graph_tbl),
    CurveGroups = if ("CurveGroup" %in% names(graph_tbl)) length(unique(as.character(graph_tbl$CurveGroup))) else NA_integer_,
    ThetaPoints = if ("Scale" %in% names(graph_tbl)) length(unique(suppressWarnings(as.numeric(graph_tbl$Scale)))) else NA_integer_,
    ProbabilityColumns = length(prob_cols),
    stringsAsFactors = FALSE
  )
  summarize_known_bundle(
    object = object,
    obj_class = "mfrm_category_curves",
    summary_candidates = character(0),
    preview_candidates = c("expected_ogive", "graphfile", "probabilities"),
    settings_candidates = "settings",
    notes = "FACETS Table 8 curves style expected-score and category-probability bundle.",
    digits = digits,
    top_n = top_n,
    summary_override = summary_tbl
  )
}

#' Summarize report/table bundles in a user-friendly format
#'
#' @param object Any report bundle produced by `mfrmr` table/report helpers.
#' @param digits Number of digits for printed numeric values.
#' @param top_n Number of preview rows shown from the main table component.
#' @param ... Reserved for generic compatibility.
#'
#' @details
#' This method provides a compact summary for bundle-like outputs
#' (for example: unexpected-response, fair-average, chi-square, and
#' category report objects). It extracts:
#' - object class and available components
#' - one-row summary table when available
#' - preview rows from the main data component
#' - resolved settings/options
#'
#' Branch-aware summaries are provided for:
#' - `mfrm_bias_count` (`branch = "original"` / `"facets"`)
#' - `mfrm_fixed_reports` (`branch = "original"` / `"facets"`)
#' - `mfrm_visual_summaries` (`branch = "original"` / `"facets"`)
#'
#' Additional class-aware summaries are provided for:
#' - `mfrm_unexpected`, `mfrm_fair_average`, `mfrm_displacement`
#' - `mfrm_interrater`, `mfrm_facets_chisq`, `mfrm_bias_interaction`
#' - `mfrm_rating_scale`, `mfrm_category_structure`, `mfrm_category_curves`
#' - `mfrm_measurable`, `mfrm_unexpected_after_bias`, `mfrm_output_bundle`
#' - `mfrm_residual_pca`, `mfrm_specifications`, `mfrm_data_quality`
#' - `mfrm_iteration_report`, `mfrm_subset_connectivity`, `mfrm_facet_statistics`
#' - `mfrm_parity_report`
#'
#' @section Interpreting output:
#' - `overview`: class, component count, and selected preview component.
#' - `summary`: one-row aggregate block when supplied by the bundle.
#' - `preview`: first `top_n` rows from the main table-like component.
#' - `settings`: resolved option values if available.
#'
#' @section Typical workflow:
#' 1. Generate a bundle table/report helper output.
#' 2. Run `summary(bundle)` for compact QA.
#' 3. Drill into specific components via `$` and visualize with `plot(bundle, ...)`.
#'
#' @return An object of class `summary.mfrm_bundle`.
#' @seealso [unexpected_response_table()], [fair_average_table()], `plot()`
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' t4 <- unexpected_response_table(fit, abs_z_min = 1.5, prob_max = 0.4, top_n = 10)
#' summary(t4)
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#' bias <- estimate_bias(fit, diag, facet_a = "Rater", facet_b = "Criterion", max_iter = 2)
#' t11 <- bias_count_table(bias, branch = "facets")
#' summary(t11)
#' @export
summary.mfrm_bundle <- function(object, digits = 3, top_n = 10, ...) {
  if (!is.list(object)) {
    stop("`object` must be a bundle-like list output.")
  }
  digits <- max(0L, as.integer(digits))
  top_n <- max(1L, as.integer(top_n))

  if (inherits(object, "mfrm_bias_count")) {
    return(summarize_bias_count_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_fixed_reports")) {
    return(summarize_fixed_reports_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_visual_summaries")) {
    return(summarize_visual_summaries_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_unexpected")) {
    return(summarize_unexpected_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_fair_average")) {
    return(summarize_fair_average_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_displacement")) {
    return(summarize_displacement_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_interrater")) {
    return(summarize_interrater_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_facets_chisq")) {
    return(summarize_facets_chisq_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_bias_interaction")) {
    return(summarize_bias_interaction_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_rating_scale")) {
    return(summarize_rating_scale_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_category_structure")) {
    return(summarize_category_structure_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_category_curves")) {
    return(summarize_category_curves_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_measurable")) {
    return(summarize_measurable_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_unexpected_after_bias")) {
    return(summarize_unexpected_after_bias_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_output_bundle")) {
    return(summarize_output_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_residual_pca")) {
    return(summarize_residual_pca_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_specifications")) {
    return(summarize_specifications_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_data_quality")) {
    return(summarize_data_quality_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_iteration_report")) {
    return(summarize_iteration_report_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_subset_connectivity")) {
    return(summarize_subset_connectivity_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_facet_statistics")) {
    return(summarize_facet_statistics_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_parity_report")) {
    return(summarize_parity_bundle(object, digits = digits, top_n = top_n))
  }

  cls <- class(object)
  cls <- cls[!cls %in% c("list", "mfrm_bundle")]
  obj_class <- if (length(cls) == 0) "mfrm_bundle" else cls[1]

  comp_names <- names(object)
  if (is.null(comp_names)) comp_names <- character(0)

  summary_tbl <- if ("summary" %in% comp_names && is.data.frame(object$summary)) {
    as.data.frame(object$summary, stringsAsFactors = FALSE)
  } else {
    data.frame()
  }
  preview <- bundle_preview_table(object, top_n = top_n)
  settings_tbl <- if ("settings" %in% comp_names) bundle_settings_table(object$settings) else data.frame()

  overview <- data.frame(
    Class = obj_class,
    Components = length(comp_names),
    ComponentNames = if (length(comp_names) == 0) "" else paste(comp_names, collapse = ", "),
    PreviewComponent = ifelse(is.na(preview$name), "", preview$name),
    PreviewRows = nrow(preview$table),
    stringsAsFactors = FALSE
  )

  notes <- if (nrow(summary_tbl) > 0) {
    "Summary table and preview rows were extracted."
  } else if (nrow(preview$table) > 0) {
    "No `summary` component found; showing preview rows from the main table."
  } else {
    "No tabular components available for preview."
  }

  out <- list(
    overview = overview,
    summary = summary_tbl,
    preview_name = preview$name,
    preview = preview$table,
    settings = settings_tbl,
    notes = notes,
    digits = digits
  )
  class(out) <- "summary.mfrm_bundle"
  out
}

bundle_summary_labels <- function(summary_kind, overview = NULL) {
  class_name <- NA_character_
  if (!is.null(overview) && is.data.frame(overview) && nrow(overview) > 0 && "Class" %in% names(overview)) {
    class_name <- as.character(overview$Class[1])
  }
  key <- as.character(summary_kind %||% class_name %||% "")
  if (!nzchar(key) || identical(key, "NA")) key <- "mfrm_bundle"

  defaults <- list(
    title = "mfrmr Bundle Summary",
    summary = "Summary table",
    preview = "Preview",
    settings = "Settings"
  )

  maps <- list(
    mfrm_unexpected = list(title = "mfrmr Unexpected Response Summary", summary = "Threshold summary", preview = "Flagged responses"),
    mfrm_fair_average = list(title = "mfrmr Fair Average Summary", summary = "Overview", preview = "Facet-level fair averages"),
    mfrm_displacement = list(title = "mfrmr Displacement Summary", summary = "Displacement summary", preview = "Displacement rows"),
    mfrm_interrater = list(title = "mfrmr Inter-rater Agreement Summary", summary = "Agreement summary", preview = "Rater-pair rows"),
    mfrm_facets_chisq = list(title = "mfrmr Facet Chi-square Summary", summary = "Facet chi-square summary", preview = "Facet rows"),
    mfrm_bias_interaction = list(title = "mfrmr Bias Interaction Summary", summary = "Interaction summary", preview = "Ranked interaction rows"),
    mfrm_rating_scale = list(title = "mfrmr Rating Scale Summary", summary = "Category/threshold summary", preview = "Category rows"),
    mfrm_category_structure = list(title = "mfrmr Category Structure Summary", summary = "Category structure overview", preview = "Category structure rows"),
    mfrm_category_curves = list(title = "mfrmr Category Curves Summary", summary = "Curve grid summary", preview = "Expected-score / curve rows"),
    mfrm_measurable = list(title = "mfrmr Measurable Summary", summary = "Run overview", preview = "Facet/category rows"),
    mfrm_unexpected_after_bias = list(title = "mfrmr Unexpected-after-Bias Summary", summary = "After-bias threshold summary", preview = "After-bias flagged rows"),
    mfrm_output_bundle = list(title = "mfrmr Output File Bundle Summary", summary = "Output overview", preview = "Output preview rows"),
    mfrm_residual_pca = list(title = "mfrmr Residual PCA Summary", summary = "PCA overview", preview = "Eigenvalue / loading rows"),
    mfrm_specifications = list(title = "mfrmr Specifications Summary", summary = "Specification header", preview = "Specification rows"),
    mfrm_data_quality = list(title = "mfrmr Data Quality Summary", summary = "Data quality overview", preview = "Audit rows"),
    mfrm_iteration_report = list(title = "mfrmr Iteration Report Summary", summary = "Iteration overview", preview = "Iteration rows"),
    mfrm_subset_connectivity = list(title = "mfrmr Subset Connectivity Summary", summary = "Subset overview", preview = "Subset/node rows"),
    mfrm_facet_statistics = list(title = "mfrmr Facet Statistics Summary", summary = "Facet-statistics overview", preview = "Facet-statistics rows"),
    mfrm_parity_report = list(title = "mfrmr FACETS Parity Summary", summary = "Overall parity", preview = "Lowest-coverage components")
  )

  if (key %in% names(maps)) {
    out <- utils::modifyList(defaults, maps[[key]])
  } else {
    out <- defaults
  }
  out
}

print_bundle_section <- function(title, table, digits = 3, round_numeric = TRUE) {
  if (is.null(table) || !is.data.frame(table) || nrow(table) == 0) return(invisible(NULL))
  cat("\n", title, "\n", sep = "")
  if (isTRUE(round_numeric)) {
    print(round_numeric_df(as.data.frame(table), digits = digits), row.names = FALSE)
  } else {
    print(as.data.frame(table), row.names = FALSE)
  }
  invisible(NULL)
}

#' @export
print.summary.mfrm_bundle <- function(x, ...) {
  digits <- as.integer(x$digits %||% 3L)
  if (!is.finite(digits)) digits <- 3L

  if (identical(x$summary_kind, "bias_count")) {
    cat("mfrmr Bias Count Summary\n")
    if (!is.null(x$overview) && nrow(x$overview) > 0) {
      cat("\nOverview\n")
      print(round_numeric_df(as.data.frame(x$overview), digits = digits), row.names = FALSE)
    }
    if (!is.null(x$count_distribution) && nrow(x$count_distribution) > 0) {
      cat("\nCount distribution\n")
      print(round_numeric_df(as.data.frame(x$count_distribution), digits = digits), row.names = FALSE)
    }
    if (!is.null(x$low_count_cells) && nrow(x$low_count_cells) > 0) {
      cat("\nLow-count cells (preview)\n")
      print(round_numeric_df(as.data.frame(x$low_count_cells), digits = digits), row.names = FALSE)
    }
    if (!is.null(x$thresholds) && nrow(x$thresholds) > 0) {
      cat("\nThresholds\n")
      print(as.data.frame(x$thresholds), row.names = FALSE)
    }
    if (length(x$notes) > 0) {
      cat("\nNotes\n")
      cat(" - ", x$notes, "\n", sep = "")
    }
    return(invisible(x))
  }

  if (identical(x$summary_kind, "visual_summaries")) {
    cat("mfrmr Visual Summary Bundle\n")
    if (!is.null(x$overview) && nrow(x$overview) > 0) {
      cat("\nOverview\n")
      print(round_numeric_df(as.data.frame(x$overview), digits = digits), row.names = FALSE)
    }
    if (!is.null(x$summary) && nrow(x$summary) > 0) {
      cat("\nWarning counts\n")
      print(round_numeric_df(as.data.frame(x$summary), digits = digits), row.names = FALSE)
    }
    if (!is.null(x$summary_counts) && nrow(x$summary_counts) > 0) {
      cat("\nSummary counts\n")
      print(round_numeric_df(as.data.frame(x$summary_counts), digits = digits), row.names = FALSE)
    }
    if (!is.null(x$settings) && nrow(x$settings) > 0) {
      cat("\nFACETS crosswalk\n")
      print(as.data.frame(x$settings), row.names = FALSE)
    }
    if (length(x$notes) > 0) {
      cat("\nNotes\n")
      cat(" - ", x$notes, "\n", sep = "")
    }
    return(invisible(x))
  }

  if (identical(x$summary_kind, "fixed_reports")) {
    cat("mfrmr Fixed-Report Bundle\n")
    if (!is.null(x$overview) && nrow(x$overview) > 0) {
      cat("\nOverview\n")
      print(round_numeric_df(as.data.frame(x$overview), digits = digits), row.names = FALSE)
    }
    if (!is.null(x$preview) && nrow(x$preview) > 0) {
      cat("\nPairwise preview\n")
      print(round_numeric_df(as.data.frame(x$preview), digits = digits), row.names = FALSE)
    }
    if (length(x$notes) > 0) {
      cat("\nNotes\n")
      cat(" - ", x$notes, "\n", sep = "")
    }
    return(invisible(x))
  }

  labels <- bundle_summary_labels(summary_kind = x$summary_kind, overview = x$overview)
  cat(labels$title, "\n")
  if (!is.null(x$overview) && nrow(x$overview) > 0) {
    ov <- x$overview[1, , drop = FALSE]
    cat(sprintf("  Class: %s\n", ov$Class))
    cat(sprintf("  Components (%s): %s\n", ov$Components, ov$ComponentNames))
  }
  print_bundle_section(labels$summary, x$summary, digits = digits, round_numeric = TRUE)
  if (!is.null(x$preview) && nrow(x$preview) > 0) {
    preview_title <- labels$preview
    if (!is.null(x$preview_name) && !is.na(x$preview_name) && nzchar(x$preview_name)) {
      preview_title <- paste0(preview_title, ": ", x$preview_name)
    }
    print_bundle_section(preview_title, x$preview, digits = digits, round_numeric = TRUE)
  }
  if (!is.null(x$settings) && nrow(x$settings) > 0) {
    print_bundle_section(labels$settings, x$settings, digits = digits, round_numeric = FALSE)
  }
  if (length(x$notes) > 0) {
    cat("\nNotes\n")
    cat(" - ", x$notes, "\n", sep = "")
  }
  invisible(x)
}

draw_category_structure_bundle <- function(x,
                                           type = c("counts", "mode_boundaries", "mean_halfscore"),
                                           draw = TRUE,
                                           main = NULL,
                                           palette = NULL,
                                           label_angle = 45) {
  type <- match.arg(tolower(type), c("counts", "mode_boundaries", "mean_halfscore"))
  pal <- resolve_palette(
    palette = palette,
    defaults = c(
      counts = "#9ecae1",
      expected = "#08519c",
      mode = "#2b8cbe",
      mean = "#238b45"
    )
  )
  cat_tbl <- as.data.frame(x$category_table %||% data.frame(), stringsAsFactors = FALSE)
  mode_tbl <- as.data.frame(x$mode_boundaries %||% data.frame(), stringsAsFactors = FALSE)
  half_tbl <- as.data.frame(x$mean_halfscore_points %||% data.frame(), stringsAsFactors = FALSE)

  if (isTRUE(draw)) {
    if (type == "counts") {
      if (nrow(cat_tbl) == 0 || !all(c("Category", "Count") %in% names(cat_tbl))) stop("No category count data available.")
      bp <- barplot_rot45(
        height = suppressWarnings(as.numeric(cat_tbl$Count)),
        labels = as.character(cat_tbl$Category),
        col = pal["counts"],
        main = if (is.null(main)) "Category counts" else as.character(main[1]),
        ylab = "Count",
        label_angle = label_angle,
        mar_bottom = 8.2
      )
      if ("ExpectedCount" %in% names(cat_tbl)) {
        exp_ct <- suppressWarnings(as.numeric(cat_tbl$ExpectedCount))
        if (any(is.finite(exp_ct))) {
          graphics::points(bp, exp_ct, pch = 21, bg = "white", col = pal["expected"])
          graphics::lines(bp, exp_ct, col = pal["expected"], lwd = 1.4)
        }
      }
    } else if (type == "mode_boundaries") {
      if (nrow(mode_tbl) == 0 || !all(c("CurveGroup", "ModeBoundaryTheta") %in% names(mode_tbl))) {
        stop("No mode-boundary data available.")
      }
      grp <- as.factor(mode_tbl$CurveGroup)
      y <- as.numeric(grp)
      graphics::plot(
        x = suppressWarnings(as.numeric(mode_tbl$ModeBoundaryTheta)),
        y = y,
        pch = 16,
        col = pal["mode"],
        xlab = "Theta / Logit",
        ylab = "",
        yaxt = "n",
        main = if (is.null(main)) "Mode boundaries" else as.character(main[1])
      )
      graphics::axis(side = 2, at = seq_along(levels(grp)), labels = levels(grp), las = 2)
    } else {
      if (nrow(half_tbl) == 0 || !all(c("CurveGroup", "MeanBoundaryTheta") %in% names(half_tbl))) {
        stop("No mean half-score data available.")
      }
      grp <- as.factor(half_tbl$CurveGroup)
      y <- as.numeric(grp)
      graphics::plot(
        x = suppressWarnings(as.numeric(half_tbl$MeanBoundaryTheta)),
        y = y,
        pch = 16,
        col = pal["mean"],
        xlab = "Theta / Logit",
        ylab = "",
        yaxt = "n",
        main = if (is.null(main)) "Mean half-score boundaries" else as.character(main[1])
      )
      graphics::axis(side = 2, at = seq_along(levels(grp)), labels = levels(grp), las = 2)
    }
  }

  new_mfrm_plot_data(
    "category_structure",
    list(
      plot = type,
      category_table = cat_tbl,
      mode_boundaries = mode_tbl,
      mean_halfscore_points = half_tbl
    )
  )
}

draw_category_curves_bundle <- function(x,
                                        type = c("ogive", "ccc"),
                                        draw = TRUE,
                                        main = NULL,
                                        palette = NULL) {
  type <- match.arg(tolower(type), c("ogive", "ccc"))
  ogive <- as.data.frame(x$expected_ogive %||% data.frame(), stringsAsFactors = FALSE)
  probs <- as.data.frame(x$probabilities %||% data.frame(), stringsAsFactors = FALSE)

  if (isTRUE(draw)) {
    if (type == "ogive") {
      if (nrow(ogive) == 0 || !all(c("Theta", "ExpectedScore", "CurveGroup") %in% names(ogive))) {
        stop("No expected-ogive data available.")
      }
      groups <- unique(as.character(ogive$CurveGroup))
      defaults <- stats::setNames(grDevices::hcl.colors(max(3L, length(groups)), "Dark 3")[seq_along(groups)], groups)
      cols <- resolve_palette(palette = palette, defaults = defaults)
      graphics::plot(
        x = range(ogive$Theta, finite = TRUE),
        y = range(ogive$ExpectedScore, finite = TRUE),
        type = "n",
        xlab = "Theta / Logit",
        ylab = "Expected score",
        main = if (is.null(main)) "Expected-score ogive" else as.character(main[1])
      )
      for (i in seq_along(groups)) {
        sub <- ogive[ogive$CurveGroup == groups[i], , drop = FALSE]
        graphics::lines(sub$Theta, sub$ExpectedScore, col = cols[groups[i]], lwd = 2)
      }
      graphics::legend("topleft", legend = groups, col = cols[groups], lty = 1, lwd = 2, bty = "n")
    } else {
      if (nrow(probs) == 0 || !all(c("Theta", "Probability", "Category", "CurveGroup") %in% names(probs))) {
        stop("No category-curve data available.")
      }
      traces <- unique(paste(probs$CurveGroup, probs$Category, sep = " | Cat "))
      defaults <- stats::setNames(grDevices::hcl.colors(max(3L, length(traces)), "Dark 3")[seq_along(traces)], traces)
      cols <- resolve_palette(palette = palette, defaults = defaults)
      graphics::plot(
        x = range(probs$Theta, finite = TRUE),
        y = c(0, 1),
        type = "n",
        xlab = "Theta / Logit",
        ylab = "Probability",
        main = if (is.null(main)) "Category characteristic curves" else as.character(main[1])
      )
      for (i in seq_along(traces)) {
        parts <- strsplit(traces[i], " \\| Cat ")[[1]]
        sub <- probs[probs$CurveGroup == parts[1] & probs$Category == parts[2], , drop = FALSE]
        graphics::lines(sub$Theta, sub$Probability, col = cols[traces[i]], lwd = 1.4)
      }
    }
  }

  new_mfrm_plot_data(
    "category_curves",
    list(plot = type, expected_ogive = ogive, probabilities = probs)
  )
}

draw_rating_scale_bundle <- function(x,
                                     type = c("counts", "thresholds"),
                                     draw = TRUE,
                                     main = NULL,
                                     palette = NULL,
                                     label_angle = 45) {
  type <- match.arg(tolower(type), c("counts", "thresholds"))
  pal <- resolve_palette(
    palette = palette,
    defaults = c(
      counts = "#c7e9c0",
      expected = "#08519c",
      step_line = "#1B9E77"
    )
  )
  cat_tbl <- as.data.frame(x$category_table %||% data.frame(), stringsAsFactors = FALSE)
  thr_tbl <- as.data.frame(x$threshold_table %||% data.frame(), stringsAsFactors = FALSE)

  if (isTRUE(draw)) {
    if (type == "counts") {
      if (nrow(cat_tbl) == 0 || !all(c("Category", "Count") %in% names(cat_tbl))) {
        stop("No category count data available.")
      }
      bp <- barplot_rot45(
        height = suppressWarnings(as.numeric(cat_tbl$Count)),
        labels = as.character(cat_tbl$Category),
        col = pal["counts"],
        main = if (is.null(main)) "Rating-scale category counts" else as.character(main[1]),
        ylab = "Count",
        label_angle = label_angle,
        mar_bottom = 8.2
      )
      if ("ExpectedCount" %in% names(cat_tbl)) {
        exp_ct <- suppressWarnings(as.numeric(cat_tbl$ExpectedCount))
        if (any(is.finite(exp_ct))) {
          graphics::points(bp, exp_ct, pch = 21, bg = "white", col = pal["expected"])
          graphics::lines(bp, exp_ct, col = pal["expected"], lwd = 1.3)
        }
      }
    } else {
      if (nrow(thr_tbl) == 0 || !all(c("Step", "Estimate") %in% names(thr_tbl))) {
        stop("No threshold data available.")
      }
      draw_step_plot(
        thr_tbl,
        title = if (is.null(main)) "Rating-scale thresholds" else as.character(main[1]),
        palette = c(step_line = pal["step_line"]),
        label_angle = label_angle
      )
    }
  }

  new_mfrm_plot_data(
    "rating_scale",
    list(plot = type, category_table = cat_tbl, threshold_table = thr_tbl)
  )
}

draw_measurable_bundle <- function(x,
                                   type = c("facet_coverage", "category_counts", "subset_observations"),
                                   draw = TRUE,
                                   main = NULL,
                                   palette = NULL,
                                   label_angle = 45) {
  type <- match.arg(tolower(as.character(type[1])), c("facet_coverage", "category_counts", "subset_observations"))
  pal <- resolve_palette(
    palette = palette,
    defaults = c(
      facet = "#2b8cbe",
      category = "#31a354",
      subset = "#756bb1"
    )
  )
  facet_tbl <- as.data.frame(x$facet_coverage %||% data.frame(), stringsAsFactors = FALSE)
  cat_tbl <- as.data.frame(x$category_stats %||% data.frame(), stringsAsFactors = FALSE)
  sub_tbl <- as.data.frame(x$subsets %||% data.frame(), stringsAsFactors = FALSE)

  if (type == "facet_coverage") {
    if (nrow(facet_tbl) == 0 || !all(c("Facet", "Levels") %in% names(facet_tbl))) {
      stop("No facet-coverage table available.")
    }
    vals <- suppressWarnings(as.numeric(facet_tbl$Levels))
    labels <- as.character(facet_tbl$Facet)
    if (isTRUE(draw)) {
      barplot_rot45(
        height = vals,
        labels = labels,
        col = pal["facet"],
        main = if (is.null(main)) "Facet coverage (levels per facet)" else as.character(main[1]),
        ylab = "Levels",
        label_angle = label_angle,
        mar_bottom = 7.8
      )
    }
    return(invisible(new_mfrm_plot_data(
      "measurable",
      list(plot = "facet_coverage", table = facet_tbl)
    )))
  }

  if (type == "category_counts") {
    if (nrow(cat_tbl) == 0 || !all(c("Category", "Count") %in% names(cat_tbl))) {
      stop("No category-statistics table available.")
    }
    vals <- suppressWarnings(as.numeric(cat_tbl$Count))
    labels <- as.character(cat_tbl$Category)
    if (isTRUE(draw)) {
      barplot_rot45(
        height = vals,
        labels = labels,
        col = pal["category"],
        main = if (is.null(main)) "Category counts (measurable data)" else as.character(main[1]),
        ylab = "Count",
        label_angle = label_angle,
        mar_bottom = 7.8
      )
    }
    return(invisible(new_mfrm_plot_data(
      "measurable",
      list(plot = "category_counts", table = cat_tbl)
    )))
  }

  if (nrow(sub_tbl) == 0 || !all(c("Subset", "Observations") %in% names(sub_tbl))) {
    stop("No subset summary available.")
  }
  vals <- suppressWarnings(as.numeric(sub_tbl$Observations))
  labels <- paste0("Subset ", as.character(sub_tbl$Subset))
  if (isTRUE(draw)) {
    barplot_rot45(
      height = vals,
      labels = labels,
      col = pal["subset"],
      main = if (is.null(main)) "Observations by subset" else as.character(main[1]),
      ylab = "Observations",
      label_angle = label_angle,
      mar_bottom = 7.8
    )
  }
  invisible(new_mfrm_plot_data(
    "measurable",
    list(plot = "subset_observations", table = sub_tbl)
  ))
}

draw_unexpected_after_bias_bundle <- function(x,
                                              type = c("scatter", "severity", "comparison"),
                                              top_n = 40,
                                              draw = TRUE,
                                              main = NULL,
                                              palette = NULL,
                                              label_angle = 45) {
  type <- match.arg(tolower(as.character(type[1])), c("scatter", "severity", "comparison"))
  top_n <- max(1L, as.integer(top_n))
  pal <- resolve_palette(
    palette = palette,
    defaults = c(
      higher = "#d95f02",
      lower = "#1b9e77",
      severity = "#2b8cbe",
      baseline = "#9ecae1",
      after = "#3182bd"
    )
  )
  tbl <- as.data.frame(x$table %||% data.frame(), stringsAsFactors = FALSE)
  summary_tbl <- as.data.frame(x$summary %||% data.frame(), stringsAsFactors = FALSE)
  thr <- x$thresholds %||% list(abs_z_min = 2, prob_max = 0.30)

  if (type == "comparison") {
    if (nrow(summary_tbl) == 0) stop("No summary table available.")
    base_n <- suppressWarnings(as.numeric(summary_tbl$BaselineUnexpectedN[1] %||% NA_real_))
    after_n <- suppressWarnings(as.numeric(summary_tbl$AfterBiasUnexpectedN[1] %||% NA_real_))
    vals <- c(base_n, after_n)
    if (!all(is.finite(vals))) stop("Baseline/after-bias counts are not available.")
    labels <- c("Baseline", "After bias")
    if (isTRUE(draw)) {
      mids <- graphics::barplot(
        height = vals,
        col = c(pal["baseline"], pal["after"]),
        names.arg = labels,
        ylab = "Unexpected responses",
        main = if (is.null(main)) "Unexpected responses: baseline vs after bias" else as.character(main[1]),
        border = "white"
      )
      graphics::text(mids, vals, labels = as.integer(vals), pos = 3, cex = 0.85)
    }
    return(invisible(new_mfrm_plot_data(
      "unexpected_after_bias",
      list(plot = "comparison", baseline = base_n, after = after_n)
    )))
  }

  if (nrow(tbl) == 0) stop("No unexpected-after-bias rows available.")

  if (type == "scatter") {
    x_vals <- suppressWarnings(as.numeric(tbl$StdResidual))
    y_vals <- -log10(pmax(suppressWarnings(as.numeric(tbl$ObsProb)), .Machine$double.xmin))
    dirs <- as.character(tbl$Direction %||% rep(NA_character_, nrow(tbl)))
    cols <- ifelse(dirs == "Higher than expected", pal["higher"], pal["lower"])
    cols[!is.finite(x_vals) | !is.finite(y_vals)] <- "gray60"
    if (isTRUE(draw)) {
      graphics::plot(
        x = x_vals,
        y = y_vals,
        xlab = "Standardized residual",
        ylab = expression(-log[10](P[obs])),
        main = if (is.null(main)) "Unexpected responses after bias adjustment" else as.character(main[1]),
        pch = 16,
        col = cols
      )
      z_thr <- as.numeric(thr$abs_z_min %||% 2)
      p_thr <- as.numeric(thr$prob_max %||% 0.30)
      graphics::abline(v = c(-z_thr, z_thr), lty = 2, col = "gray45")
      graphics::abline(h = -log10(p_thr), lty = 2, col = "gray45")
      graphics::legend(
        "topleft",
        legend = c("Higher than expected", "Lower than expected"),
        col = c(pal["higher"], pal["lower"]),
        pch = 16,
        bty = "n",
        cex = 0.85
      )
    }
    return(invisible(new_mfrm_plot_data(
      "unexpected_after_bias",
      list(plot = "scatter", table = tbl, thresholds = thr)
    )))
  }

  sev <- suppressWarnings(as.numeric(tbl$Severity))
  sev <- sev[is.finite(sev)]
  if (length(sev) == 0) stop("No finite severity values available.")
  ord <- order(suppressWarnings(as.numeric(tbl$Severity)), decreasing = TRUE, na.last = NA)
  use <- ord[seq_len(min(length(ord), top_n))]
  sub <- tbl[use, , drop = FALSE]
  labels <- if ("Row" %in% names(sub)) paste0("Row ", sub$Row) else paste0("Case ", seq_len(nrow(sub)))
  vals <- suppressWarnings(as.numeric(sub$Severity))
  if (isTRUE(draw)) {
    barplot_rot45(
      height = vals,
      labels = labels,
      col = pal["severity"],
      main = if (is.null(main)) "Unexpected-response severity after bias" else as.character(main[1]),
      ylab = "Severity",
      label_angle = label_angle,
      mar_bottom = 8.2
    )
  }
  invisible(new_mfrm_plot_data(
    "unexpected_after_bias",
    list(plot = "severity", table = sub)
  ))
}

draw_output_bundle <- function(x,
                               type = c("graph_expected", "score_residuals", "obs_probability"),
                               draw = TRUE,
                               main = NULL,
                               palette = NULL) {
  type <- match.arg(tolower(as.character(type[1])), c("graph_expected", "score_residuals", "obs_probability"))
  graph_tbl <- as.data.frame(x$graphfile %||% data.frame(), stringsAsFactors = FALSE)
  score_tbl <- as.data.frame(x$scorefile %||% data.frame(), stringsAsFactors = FALSE)

  if (type == "graph_expected") {
    if (nrow(graph_tbl) == 0 || !all(c("Measure", "Expected") %in% names(graph_tbl))) {
      stop("Graphfile table with `Measure` and `Expected` is not available. Run the full workflow (fit_mfrm -> diagnose_mfrm) first.", call. = FALSE)
    }
    groups <- if ("CurveGroup" %in% names(graph_tbl)) unique(as.character(graph_tbl$CurveGroup)) else "All"
    if (!"CurveGroup" %in% names(graph_tbl)) graph_tbl$CurveGroup <- "All"
    defaults <- stats::setNames(grDevices::hcl.colors(max(3L, length(groups)), "Dark 3")[seq_along(groups)], groups)
    cols <- resolve_palette(palette = palette, defaults = defaults)
    if (isTRUE(draw)) {
      graphics::plot(
        x = range(graph_tbl$Measure, finite = TRUE),
        y = range(graph_tbl$Expected, finite = TRUE),
        type = "n",
        xlab = "Theta / Logit",
        ylab = "Expected score",
        main = if (is.null(main)) "Graphfile expected-score curves" else as.character(main[1])
      )
      for (g in groups) {
        sub <- graph_tbl[as.character(graph_tbl$CurveGroup) == g, , drop = FALSE]
        sub <- sub[order(sub$Measure), , drop = FALSE]
        graphics::lines(sub$Measure, sub$Expected, col = cols[g], lwd = 1.8)
      }
      if (length(groups) > 1) {
        graphics::legend("topleft", legend = groups, col = cols[groups], lty = 1, lwd = 2, bty = "n", cex = 0.85)
      }
    }
    return(invisible(new_mfrm_plot_data(
      "output_bundle",
      list(plot = "graph_expected", table = graph_tbl)
    )))
  }

  if (nrow(score_tbl) == 0) stop("Scorefile table is not available. Run the full workflow (fit_mfrm -> diagnose_mfrm) first.", call. = FALSE)

  if (type == "score_residuals") {
    if (!"Residual" %in% names(score_tbl)) stop("`Residual` column is not available in scorefile. Run the full workflow (fit_mfrm -> diagnose_mfrm) first.", call. = FALSE)
    vals <- suppressWarnings(as.numeric(score_tbl$Residual))
    vals <- vals[is.finite(vals)]
    if (length(vals) == 0) stop("No finite residual values available.")
    if (isTRUE(draw)) {
      graphics::hist(
        x = vals,
        breaks = "FD",
        col = "#9ecae1",
        border = "white",
        main = if (is.null(main)) "Scorefile residual distribution" else as.character(main[1]),
        xlab = "Residual",
        ylab = "Count"
      )
      graphics::abline(v = 0, lty = 2, col = "gray45")
    }
    return(invisible(new_mfrm_plot_data(
      "output_bundle",
      list(plot = "score_residuals", values = vals)
    )))
  }

  if (!"ObsProb" %in% names(score_tbl)) stop("`ObsProb` column is not available in scorefile. Run the full workflow (fit_mfrm -> diagnose_mfrm) first.", call. = FALSE)
  vals <- suppressWarnings(as.numeric(score_tbl$ObsProb))
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0) stop("No finite observed-probability values available.")
  if (isTRUE(draw)) {
    graphics::hist(
      x = vals,
      breaks = "FD",
      col = "#c7e9c0",
      border = "white",
      main = if (is.null(main)) "Observed probability distribution" else as.character(main[1]),
      xlab = "Observed probability",
      ylab = "Count"
    )
  }
  invisible(new_mfrm_plot_data(
    "output_bundle",
    list(plot = "obs_probability", values = vals)
  ))
}

draw_specifications_bundle <- function(x,
                                       type = c("facet_elements", "anchor_constraints", "convergence"),
                                       draw = TRUE,
                                       main = NULL,
                                       palette = NULL,
                                       label_angle = 45) {
  type <- match.arg(tolower(as.character(type[1])), c("facet_elements", "anchor_constraints", "convergence"))
  pal <- resolve_palette(
    palette = palette,
    defaults = c(
      facet = "#2b8cbe",
      anchor = "#756bb1",
      group = "#9ecae1",
      free = "#d9d9d9",
      convergence = "#31a354"
    )
  )
  facet_tbl <- as.data.frame(x$facet_labels %||% data.frame(), stringsAsFactors = FALSE)
  anchor_tbl <- as.data.frame(x$anchor_summary %||% data.frame(), stringsAsFactors = FALSE)
  conv_tbl <- as.data.frame(x$convergence_control %||% data.frame(), stringsAsFactors = FALSE)

  if (type == "facet_elements") {
    if (nrow(facet_tbl) == 0 || !all(c("Facet", "Elements") %in% names(facet_tbl))) {
      stop("Facet-label table is not available. Run the full workflow (fit_mfrm -> diagnose_mfrm) first.", call. = FALSE)
    }
    vals <- suppressWarnings(as.numeric(facet_tbl$Elements))
    labels <- as.character(facet_tbl$Facet)
    if (isTRUE(draw)) {
      barplot_rot45(
        height = vals,
        labels = labels,
        col = pal["facet"],
        main = if (is.null(main)) "Facet elements in model specification" else as.character(main[1]),
        ylab = "Elements",
        label_angle = label_angle,
        mar_bottom = 7.8
      )
    }
    return(invisible(new_mfrm_plot_data(
      "specifications",
      list(plot = "facet_elements", table = facet_tbl)
    )))
  }

  if (type == "anchor_constraints") {
    if (nrow(anchor_tbl) == 0 || !all(c("Facet", "AnchoredLevels", "GroupAnchors") %in% names(anchor_tbl))) {
      stop("Anchor summary table is not available. Run the full workflow (fit_mfrm -> diagnose_mfrm) first.", call. = FALSE)
    }
    base_tbl <- anchor_tbl |>
      dplyr::transmute(
        Facet = as.character(.data$Facet),
        Anchored = suppressWarnings(as.numeric(.data$AnchoredLevels)),
        Grouped = suppressWarnings(as.numeric(.data$GroupAnchors))
      )
    if (nrow(facet_tbl) > 0 && all(c("Facet", "Elements") %in% names(facet_tbl))) {
      base_tbl <- base_tbl |>
        dplyr::left_join(
          facet_tbl |>
            dplyr::transmute(Facet = as.character(.data$Facet), Elements = suppressWarnings(as.numeric(.data$Elements))),
          by = "Facet"
        )
      base_tbl$Free <- pmax(0, base_tbl$Elements - base_tbl$Anchored - base_tbl$Grouped)
    } else {
      base_tbl$Elements <- NA_real_
      base_tbl$Free <- NA_real_
    }
    base_tbl <- base_tbl[order(base_tbl$Facet), , drop = FALSE]
    if (isTRUE(draw)) {
      old_mar <- graphics::par("mar")
      on.exit(graphics::par(mar = old_mar), add = TRUE)
      mar <- old_mar
      mar[1] <- max(mar[1], 8.8)
      graphics::par(mar = mar)
      mat <- rbind(
        Anchored = base_tbl$Anchored,
        Grouped = base_tbl$Grouped,
        Free = ifelse(is.finite(base_tbl$Free), base_tbl$Free, 0)
      )
      mids <- graphics::barplot(
        height = mat,
        beside = FALSE,
        names.arg = FALSE,
        col = c(pal["anchor"], pal["group"], pal["free"]),
        border = "white",
        ylab = "Levels",
        main = if (is.null(main)) "Anchor constraints by facet" else as.character(main[1])
      )
      draw_rotated_x_labels(
        at = mids,
        labels = base_tbl$Facet,
        srt = label_angle,
        cex = 0.82,
        line_offset = 0.085
      )
      graphics::legend(
        "topright",
        legend = c("Anchored", "Grouped", "Free"),
        fill = c(pal["anchor"], pal["group"], pal["free"]),
        bty = "n",
        cex = 0.85
      )
    }
    return(invisible(new_mfrm_plot_data(
      "specifications",
      list(plot = "anchor_constraints", table = base_tbl)
    )))
  }

  if (nrow(conv_tbl) == 0 || !all(c("Setting", "Value") %in% names(conv_tbl))) {
    stop("Convergence-control table is not available. Run the full workflow (fit_mfrm -> diagnose_mfrm) first.", call. = FALSE)
  }
  keep <- c("MaxIterations", "RelativeTolerance", "QuadPoints", "FunctionEvaluations")
  sub <- conv_tbl[as.character(conv_tbl$Setting) %in% keep, , drop = FALSE]
  if (nrow(sub) == 0) stop("No numeric convergence settings found.")
  vals <- suppressWarnings(as.numeric(sub$Value))
  ok <- is.finite(vals)
  if (!any(ok)) stop("No finite numeric values in convergence settings.")
  sub <- sub[ok, , drop = FALSE]
  vals <- vals[ok]
  labels <- as.character(sub$Setting)
  if (isTRUE(draw)) {
    barplot_rot45(
      height = vals,
      labels = labels,
      col = pal["convergence"],
      main = if (is.null(main)) "Convergence controls and counts" else as.character(main[1]),
      ylab = "Value",
      label_angle = label_angle,
      mar_bottom = 8.2
    )
  }
  invisible(new_mfrm_plot_data(
    "specifications",
    list(plot = "convergence", table = sub)
  ))
}

draw_data_quality_bundle <- function(x,
                                     type = c("row_audit", "category_counts", "missing_rows"),
                                     draw = TRUE,
                                     main = NULL,
                                     palette = NULL,
                                     label_angle = 45) {
  type <- match.arg(tolower(as.character(type[1])), c("row_audit", "category_counts", "missing_rows"))
  pal <- resolve_palette(
    palette = palette,
    defaults = c(
      row_audit = "#2b8cbe",
      category = "#31a354",
      missing = "#756bb1"
    )
  )
  row_tbl <- as.data.frame(x$row_audit %||% data.frame(), stringsAsFactors = FALSE)
  cat_tbl <- as.data.frame(x$category_counts %||% data.frame(), stringsAsFactors = FALSE)
  sum_tbl <- as.data.frame(x$summary %||% data.frame(), stringsAsFactors = FALSE)

  if (type == "row_audit") {
    if (nrow(row_tbl) == 0 || !all(c("Status", "N") %in% names(row_tbl))) {
      stop("Row-audit table is not available. Run the full workflow (fit_mfrm -> diagnose_mfrm) first.", call. = FALSE)
    }
    vals <- suppressWarnings(as.numeric(row_tbl$N))
    labels <- as.character(row_tbl$Status)
    if (isTRUE(draw)) {
      barplot_rot45(
        height = vals,
        labels = labels,
        col = pal["row_audit"],
        main = if (is.null(main)) "Row-audit status counts" else as.character(main[1]),
        ylab = "Rows",
        label_angle = label_angle,
        mar_bottom = 8.2
      )
    }
    return(invisible(new_mfrm_plot_data(
      "data_quality",
      list(plot = "row_audit", table = row_tbl)
    )))
  }

  if (type == "category_counts") {
    if (nrow(cat_tbl) == 0 || !all(c("Score", "Count") %in% names(cat_tbl))) {
      stop("Category-count table is not available. Run the full workflow (fit_mfrm -> diagnose_mfrm) first.", call. = FALSE)
    }
    vals <- suppressWarnings(as.numeric(cat_tbl$Count))
    labels <- as.character(cat_tbl$Score)
    if (isTRUE(draw)) {
      barplot_rot45(
        height = vals,
        labels = labels,
        col = pal["category"],
        main = if (is.null(main)) "Observed category counts" else as.character(main[1]),
        ylab = "Count",
        label_angle = label_angle,
        mar_bottom = 7.8
      )
    }
    return(invisible(new_mfrm_plot_data(
      "data_quality",
      list(plot = "category_counts", table = cat_tbl)
    )))
  }

  if (nrow(sum_tbl) == 0) stop("Summary table is not available. Run the full workflow (fit_mfrm -> diagnose_mfrm) first.", call. = FALSE)
  row_cols <- grep("Rows$", names(sum_tbl), value = TRUE)
  if (length(row_cols) == 0) stop("No row-count columns found in summary table.")
  vals <- suppressWarnings(as.numeric(sum_tbl[1, row_cols, drop = TRUE]))
  labels <- row_cols
  if (isTRUE(draw)) {
    barplot_rot45(
      height = vals,
      labels = labels,
      col = pal["missing"],
      main = if (is.null(main)) "Missing/invalid row counts" else as.character(main[1]),
      ylab = "Rows",
      label_angle = label_angle,
      mar_bottom = 9.0
    )
  }
  invisible(new_mfrm_plot_data(
    "data_quality",
    list(plot = "missing_rows", table = data.frame(Field = labels, Count = vals, stringsAsFactors = FALSE))
  ))
}

draw_iteration_report_bundle <- function(x,
                                         type = c("residual", "logit_change", "objective"),
                                         draw = TRUE,
                                         main = NULL,
                                         palette = NULL) {
  type <- match.arg(tolower(as.character(type[1])), c("residual", "logit_change", "objective"))
  pal <- resolve_palette(
    palette = palette,
    defaults = c(
      residual_element = "#2b8cbe",
      residual_category = "#31a354",
      change_element = "#756bb1",
      change_step = "#d95f02",
      objective = "#1b9e77"
    )
  )
  tbl <- as.data.frame(x$table %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(tbl) == 0 || !"Iteration" %in% names(tbl)) {
    stop("Iteration table is not available. Run the full workflow (fit_mfrm -> diagnose_mfrm) first.", call. = FALSE)
  }
  it <- suppressWarnings(as.numeric(tbl$Iteration))
  if (!all(is.finite(it))) it <- seq_len(nrow(tbl))

  if (type == "residual") {
    y1 <- suppressWarnings(as.numeric(tbl$MaxScoreResidualElements))
    y2 <- suppressWarnings(as.numeric(tbl$MaxScoreResidualCategories))
    if (!any(is.finite(y1)) && !any(is.finite(y2))) {
      stop("No residual metrics available.")
    }
    if (isTRUE(draw)) {
      yr <- range(c(y1, y2), finite = TRUE)
      graphics::plot(
        x = it,
        y = y1,
        type = "b",
        pch = 16,
        col = pal["residual_element"],
        ylim = yr,
        xlab = "Iteration",
        ylab = "Residual metric",
        main = if (is.null(main)) "Iteration residual trajectory" else as.character(main[1])
      )
      graphics::lines(it, y2, type = "b", pch = 17, col = pal["residual_category"])
      graphics::legend(
        "topright",
        legend = c("Elements", "Categories"),
        col = c(pal["residual_element"], pal["residual_category"]),
        pch = c(16, 17),
        lty = 1,
        bty = "n",
        cex = 0.85
      )
    }
    return(invisible(new_mfrm_plot_data(
      "iteration_report",
      list(plot = "residual", iteration = it, element = y1, category = y2)
    )))
  }

  if (type == "logit_change") {
    y1 <- suppressWarnings(as.numeric(tbl$MaxLogitChangeElements))
    y2 <- suppressWarnings(as.numeric(tbl$MaxLogitChangeSteps))
    if (!any(is.finite(y1)) && !any(is.finite(y2))) {
      stop("No logit-change metrics available.")
    }
    if (isTRUE(draw)) {
      yr <- range(c(y1, y2), finite = TRUE)
      graphics::plot(
        x = it,
        y = y1,
        type = "b",
        pch = 16,
        col = pal["change_element"],
        ylim = yr,
        xlab = "Iteration",
        ylab = "Max absolute change",
        main = if (is.null(main)) "Iteration logit-change trajectory" else as.character(main[1])
      )
      graphics::lines(it, y2, type = "b", pch = 17, col = pal["change_step"])
      graphics::legend(
        "topright",
        legend = c("Elements", "Steps"),
        col = c(pal["change_element"], pal["change_step"]),
        pch = c(16, 17),
        lty = 1,
        bty = "n",
        cex = 0.85
      )
    }
    return(invisible(new_mfrm_plot_data(
      "iteration_report",
      list(plot = "logit_change", iteration = it, element = y1, step = y2)
    )))
  }

  vals <- suppressWarnings(as.numeric(tbl$Objective))
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0) stop("No objective values available.")
  it2 <- it[is.finite(suppressWarnings(as.numeric(tbl$Objective)))]
  if (isTRUE(draw)) {
    graphics::plot(
      x = it2,
      y = vals,
      type = "b",
      pch = 16,
      col = pal["objective"],
      xlab = "Iteration",
      ylab = "Objective (log-likelihood proxy)",
      main = if (is.null(main)) "Iteration objective trajectory" else as.character(main[1])
    )
  }
  invisible(new_mfrm_plot_data(
    "iteration_report",
    list(plot = "objective", iteration = it2, objective = vals)
  ))
}

draw_subset_connectivity_bundle <- function(x,
                                            type = c("subset_observations", "facet_levels"),
                                            draw = TRUE,
                                            main = NULL,
                                            palette = NULL,
                                            label_angle = 45) {
  type <- match.arg(tolower(as.character(type[1])), c("subset_observations", "facet_levels"))
  pal <- resolve_palette(
    palette = palette,
    defaults = c(
      subset = "#756bb1",
      facet = "#2b8cbe"
    )
  )
  summary_tbl <- as.data.frame(x$summary %||% data.frame(), stringsAsFactors = FALSE)
  listing_tbl <- as.data.frame(x$listing %||% data.frame(), stringsAsFactors = FALSE)

  if (type == "subset_observations") {
    if (nrow(summary_tbl) == 0 || !all(c("Subset", "Observations") %in% names(summary_tbl))) {
      stop("Subset summary table is not available. Run the full workflow (fit_mfrm -> diagnose_mfrm) first.", call. = FALSE)
    }
    vals <- suppressWarnings(as.numeric(summary_tbl$Observations))
    labels <- paste0("Subset ", as.character(summary_tbl$Subset))
    if (isTRUE(draw)) {
      barplot_rot45(
        height = vals,
        labels = labels,
        col = pal["subset"],
        main = if (is.null(main)) "Observations by subset" else as.character(main[1]),
        ylab = "Observations",
        label_angle = label_angle,
        mar_bottom = 8.0
      )
    }
    return(invisible(new_mfrm_plot_data(
      "subset_connectivity",
      list(plot = "subset_observations", table = summary_tbl)
    )))
  }

  if (nrow(listing_tbl) == 0 || !all(c("Subset", "Facet", "LevelsN") %in% names(listing_tbl))) {
    stop("Subset facet-listing table is not available. Run the full workflow (fit_mfrm -> diagnose_mfrm) first.", call. = FALSE)
  }
  vals <- suppressWarnings(as.numeric(listing_tbl$LevelsN))
  labels <- paste0("S", listing_tbl$Subset, ":", listing_tbl$Facet)
  if (isTRUE(draw)) {
    barplot_rot45(
      height = vals,
      labels = labels,
      col = pal["facet"],
      main = if (is.null(main)) "Facet levels by subset" else as.character(main[1]),
      ylab = "Levels",
      label_angle = label_angle,
      mar_bottom = 8.8
    )
  }
  invisible(new_mfrm_plot_data(
    "subset_connectivity",
    list(plot = "facet_levels", table = listing_tbl)
  ))
}

draw_facet_statistics_bundle <- function(x,
                                         type = c("means", "sds", "ranges"),
                                         metric = NULL,
                                         draw = TRUE,
                                         main = NULL,
                                         palette = NULL,
                                         label_angle = 45) {
  type <- match.arg(tolower(as.character(type[1])), c("means", "sds", "ranges"))
  pal <- resolve_palette(
    palette = palette,
    defaults = c(
      mean = "#2b8cbe",
      sd = "#756bb1",
      range = "#9ecae1"
    )
  )
  tbl <- as.data.frame(x$table %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(tbl) == 0 || !all(c("Metric", "Facet", "Mean", "SD", "Min", "Max") %in% names(tbl))) {
    stop("Facet-statistics table is not available. Run the full workflow (fit_mfrm -> diagnose_mfrm) first.", call. = FALSE)
  }
  metrics <- unique(as.character(tbl$Metric))
  if (is.null(metric)) metric <- metrics[1]
  metric <- as.character(metric[1])
  if (!metric %in% metrics) {
    stop("Requested `metric` not found. Available: ", paste(metrics, collapse = ", "))
  }
  sub <- tbl[as.character(tbl$Metric) == metric, , drop = FALSE]
  sub <- sub[order(as.character(sub$Facet)), , drop = FALSE]

  if (type == "means") {
    vals <- suppressWarnings(as.numeric(sub$Mean))
    labels <- as.character(sub$Facet)
    if (isTRUE(draw)) {
      barplot_rot45(
        height = vals,
        labels = labels,
        col = pal["mean"],
        main = if (is.null(main)) paste0("Facet means (", metric, ")") else as.character(main[1]),
        ylab = "Mean",
        label_angle = label_angle,
        mar_bottom = 7.8
      )
    }
    return(invisible(new_mfrm_plot_data(
      "facet_statistics",
      list(plot = "means", metric = metric, table = sub)
    )))
  }

  if (type == "sds") {
    vals <- suppressWarnings(as.numeric(sub$SD))
    labels <- as.character(sub$Facet)
    if (isTRUE(draw)) {
      barplot_rot45(
        height = vals,
        labels = labels,
        col = pal["sd"],
        main = if (is.null(main)) paste0("Facet SDs (", metric, ")") else as.character(main[1]),
        ylab = "SD",
        label_angle = label_angle,
        mar_bottom = 7.8
      )
    }
    return(invisible(new_mfrm_plot_data(
      "facet_statistics",
      list(plot = "sds", metric = metric, table = sub)
    )))
  }

  y <- seq_len(nrow(sub))
  mn <- suppressWarnings(as.numeric(sub$Min))
  mx <- suppressWarnings(as.numeric(sub$Max))
  md <- suppressWarnings(as.numeric(sub$Mean))
  if (isTRUE(draw)) {
    xr <- range(c(mn, mx), finite = TRUE)
    graphics::plot(
      x = xr,
      y = c(1, nrow(sub)),
      type = "n",
      yaxt = "n",
      xlab = metric,
      ylab = "",
      main = if (is.null(main)) paste0("Facet ranges (", metric, ")") else as.character(main[1])
    )
    graphics::segments(x0 = mn, y0 = y, x1 = mx, y1 = y, col = pal["range"], lwd = 2)
    graphics::points(md, y, pch = 16, col = pal["mean"])
    graphics::axis(side = 2, at = y, labels = as.character(sub$Facet), las = 2, cex.axis = 0.8)
  }
  invisible(new_mfrm_plot_data(
    "facet_statistics",
    list(plot = "ranges", metric = metric, table = sub)
  ))
}

draw_residual_pca_bundle <- function(x,
                                     type = c("overall_scree", "facet_scree", "overall_loadings", "facet_loadings"),
                                     facet = NULL,
                                     component = 1L,
                                     top_n = 20L,
                                     draw = TRUE) {
  type <- match.arg(tolower(as.character(type[1])), c("overall_scree", "facet_scree", "overall_loadings", "facet_loadings"))
  if (type == "overall_scree") {
    return(invisible(plot_residual_pca(
      x,
      mode = "overall",
      plot_type = "scree",
      component = component,
      top_n = top_n,
      draw = draw
    )))
  }
  if (type == "facet_scree") {
    return(invisible(plot_residual_pca(
      x,
      mode = "facet",
      facet = facet,
      plot_type = "scree",
      component = component,
      top_n = top_n,
      draw = draw
    )))
  }
  if (type == "overall_loadings") {
    return(invisible(plot_residual_pca(
      x,
      mode = "overall",
      plot_type = "loadings",
      component = component,
      top_n = top_n,
      draw = draw
    )))
  }
  invisible(plot_residual_pca(
    x,
    mode = "facet",
    facet = facet,
    plot_type = "loadings",
    component = component,
    top_n = top_n,
    draw = draw
  ))
}

draw_parity_bundle <- function(x,
                               type = c("column_coverage", "table_coverage", "metric_status", "metric_by_table"),
                               top_n = 40,
                               draw = TRUE,
                               main = NULL,
                               palette = NULL,
                               label_angle = 45) {
  type <- match.arg(tolower(as.character(type[1])), c("column_coverage", "table_coverage", "metric_status", "metric_by_table"))
  top_n <- max(1L, as.integer(top_n))
  pal <- resolve_palette(
    palette = palette,
    defaults = c(
      pass = "#31a354",
      fail = "#cb181d",
      missing = "#969696",
      coverage = "#3182bd",
      metric = "#756bb1"
    )
  )

  column_audit <- as.data.frame(x$column_audit %||% data.frame(), stringsAsFactors = FALSE)
  column_summary <- as.data.frame(x$column_summary %||% data.frame(), stringsAsFactors = FALSE)
  metric_audit <- as.data.frame(x$metric_audit %||% data.frame(), stringsAsFactors = FALSE)
  metric_by_table <- as.data.frame(x$metric_by_table %||% data.frame(), stringsAsFactors = FALSE)

  if (type == "column_coverage") {
    if (nrow(column_audit) == 0 || !all(c("table_id", "component", "coverage", "available", "full_match") %in% names(column_audit))) {
      stop("Column-audit table is not available. Run the full workflow (fit_mfrm -> diagnose_mfrm) first.", call. = FALSE)
    }
    tbl <- column_audit
    tbl$coverage <- suppressWarnings(as.numeric(tbl$coverage))
    tbl <- tbl |>
      dplyr::arrange(.data$coverage, .data$table_id, .data$component)
    if (nrow(tbl) > top_n) tbl <- tbl |> dplyr::slice_head(n = top_n)
    vals <- ifelse(is.finite(tbl$coverage), tbl$coverage, 0)
    labels <- paste0(tbl$table_id, ":", tbl$component)
    cols <- ifelse(!tbl$available, pal["missing"], ifelse(tbl$full_match, pal["pass"], pal["fail"]))
    if (isTRUE(draw)) {
      barplot_rot45(
        height = vals,
        labels = labels,
        col = cols,
        main = if (is.null(main)) "Column contract coverage (lowest first)" else as.character(main[1]),
        ylab = "Coverage",
        label_angle = label_angle,
        mar_bottom = 9.2
      )
      graphics::abline(h = 1, lty = 3, col = "#999999")
    }
    return(invisible(new_mfrm_plot_data(
      "parity_report",
      list(plot = "column_coverage", table = tbl, labels = labels)
    )))
  }

  if (type == "table_coverage") {
    if (nrow(column_summary) == 0 || !all(c("table_id", "MeanCoverage") %in% names(column_summary))) {
      stop("Column-summary table is not available. Run the full workflow (fit_mfrm -> diagnose_mfrm) first.", call. = FALSE)
    }
    tbl <- column_summary |>
      dplyr::arrange(.data$table_id)
    vals <- suppressWarnings(as.numeric(tbl$MeanCoverage))
    vals[!is.finite(vals)] <- 0
    labels <- as.character(tbl$table_id)
    if (isTRUE(draw)) {
      barplot_rot45(
        height = vals,
        labels = labels,
        col = pal["coverage"],
        main = if (is.null(main)) "Mean column coverage by table" else as.character(main[1]),
        ylab = "Mean coverage",
        label_angle = label_angle,
        mar_bottom = 7.8
      )
      graphics::abline(h = 1, lty = 3, col = "#999999")
    }
    return(invisible(new_mfrm_plot_data(
      "parity_report",
      list(plot = "table_coverage", table = tbl)
    )))
  }

  if (type == "metric_status") {
    if (nrow(metric_audit) == 0 || !"Pass" %in% names(metric_audit)) {
      stop("Metric-audit table is not available. Run the full workflow (fit_mfrm -> diagnose_mfrm) first.", call. = FALSE)
    }
    status <- ifelse(is.na(metric_audit$Pass), "Not evaluated", ifelse(metric_audit$Pass %in% TRUE, "Pass", "Fail"))
    cnt <- table(factor(status, levels = c("Pass", "Fail", "Not evaluated")))
    vals <- as.numeric(cnt)
    labels <- names(cnt)
    cols <- c(pal["pass"], pal["fail"], pal["missing"])
    if (isTRUE(draw)) {
      graphics::barplot(
        height = vals,
        names.arg = labels,
        col = cols,
        las = 2,
        ylab = "Checks",
        main = if (is.null(main)) "Metric-check status counts" else as.character(main[1])
      )
    }
    return(invisible(new_mfrm_plot_data(
      "parity_report",
      list(plot = "metric_status", table = data.frame(Status = labels, Checks = vals, stringsAsFactors = FALSE))
    )))
  }

  if (nrow(metric_by_table) == 0 || !all(c("Table", "PassRate") %in% names(metric_by_table))) {
    stop("Metric-by-table summary is not available. Run the full workflow (fit_mfrm -> diagnose_mfrm) first.", call. = FALSE)
  }
  tbl <- metric_by_table |>
    dplyr::arrange(.data$Table)
  vals <- suppressWarnings(as.numeric(tbl$PassRate))
  vals[!is.finite(vals)] <- 0
  labels <- as.character(tbl$Table)
  if (isTRUE(draw)) {
    barplot_rot45(
      height = vals,
      labels = labels,
      col = pal["metric"],
      main = if (is.null(main)) "Metric pass rate by table" else as.character(main[1]),
      ylab = "Pass rate",
      label_angle = label_angle,
      mar_bottom = 7.8
    )
    graphics::abline(h = 1, lty = 3, col = "#999999")
  }
  invisible(new_mfrm_plot_data(
    "parity_report",
    list(plot = "metric_by_table", table = tbl)
  ))
}

plot_bias_count_bundle <- function(x,
                                   plot_type = c("cell_counts", "lowcount_by_facet"),
                                   top_n = 40,
                                   draw = TRUE,
                                   main = NULL,
                                   palette = NULL,
                                   label_angle = 45) {
  plot_type <- match.arg(tolower(as.character(plot_type[1])), c("cell_counts", "lowcount_by_facet"))
  top_n <- max(1L, as.integer(top_n))
  pal <- resolve_palette(
    palette = palette,
    defaults = c(
      count = "#2b8cbe",
      low = "#cb181d",
      rate = "#756bb1",
      grid = "#ececec"
    )
  )

  tbl <- as.data.frame(x$table %||% data.frame(), stringsAsFactors = FALSE)
  if ("Observd Count" %in% names(tbl) && !"Count" %in% names(tbl)) {
    tbl$Count <- suppressWarnings(as.numeric(tbl$`Observd Count`))
  }
  if (!"Count" %in% names(tbl)) {
    stop("Bias-count table does not include a count column.")
  }
  if (!"LowCountFlag" %in% names(tbl)) {
    tbl$LowCountFlag <- FALSE
  }
  tbl$LowCountFlag <- as.logical(tbl$LowCountFlag)

  if (plot_type == "cell_counts") {
    tbl <- tbl[is.finite(suppressWarnings(as.numeric(tbl$Count))), , drop = FALSE]
    if (nrow(tbl) == 0) stop("No finite count rows available.")
    tbl$Count <- suppressWarnings(as.numeric(tbl$Count))
    ord <- order(tbl$Count, decreasing = TRUE, na.last = NA)
    use <- ord[seq_len(min(length(ord), top_n))]
    tbl <- tbl[use, , drop = FALSE]

    facet_cols <- names(x$by_facet %||% list())
    facet_cols <- facet_cols[facet_cols %in% names(tbl)]
    if (length(facet_cols) == 0) {
      facet_cols <- names(tbl)[vapply(tbl, is.character, logical(1))]
      facet_cols <- setdiff(facet_cols, c("Count", "LowCountFlag"))
      facet_cols <- facet_cols[seq_len(min(2L, length(facet_cols)))]
    }
    labels <- if (length(facet_cols) > 0) {
      apply(tbl[, facet_cols, drop = FALSE], 1, paste, collapse = " | ")
    } else {
      paste0("Cell ", seq_len(nrow(tbl)))
    }

    if (isTRUE(draw)) {
      barplot_rot45(
        height = tbl$Count,
        labels = labels,
        col = ifelse(tbl$LowCountFlag %in% TRUE, pal["low"], pal["count"]),
        main = if (is.null(main)) "Bias cell counts" else as.character(main[1]),
        ylab = "Observed count",
        label_angle = label_angle,
        mar_bottom = 9.0
      )
    }
    return(invisible(new_mfrm_plot_data(
      "bias_count",
      list(plot = "cell_counts", table = tbl, labels = labels)
    )))
  }

  by_facet <- x$by_facet %||% list()
  rate_tbl <- lapply(names(by_facet), function(facet_nm) {
    df <- as.data.frame(by_facet[[facet_nm]], stringsAsFactors = FALSE)
    if (!all(c("Level", "Cells", "LowCountCells") %in% names(df))) return(NULL)
    data.frame(
      Facet = facet_nm,
      Level = as.character(df$Level),
      Cells = suppressWarnings(as.numeric(df$Cells)),
      LowCountCells = suppressWarnings(as.numeric(df$LowCountCells)),
      LowCountRate = ifelse(
        suppressWarnings(as.numeric(df$Cells)) > 0,
        suppressWarnings(as.numeric(df$LowCountCells)) / suppressWarnings(as.numeric(df$Cells)),
        NA_real_
      ),
      stringsAsFactors = FALSE
    )
  })
  rate_tbl <- rate_tbl[!vapply(rate_tbl, is.null, logical(1))]
  if (length(rate_tbl) == 0) {
    stop("No by-facet low-count summary available.")
  }
  rate_tbl <- dplyr::bind_rows(rate_tbl)
  rate_tbl <- rate_tbl[is.finite(rate_tbl$LowCountRate), , drop = FALSE]
  if (nrow(rate_tbl) == 0) {
    stop("No finite low-count rates available.")
  }
  rate_tbl <- rate_tbl |>
    dplyr::arrange(dplyr::desc(.data$LowCountRate), dplyr::desc(.data$LowCountCells), .data$Facet, .data$Level) |>
    dplyr::slice_head(n = top_n)
  labels <- paste0(rate_tbl$Facet, ":", rate_tbl$Level)

  if (isTRUE(draw)) {
    barplot_rot45(
      height = rate_tbl$LowCountRate,
      labels = labels,
      col = pal["rate"],
      main = if (is.null(main)) "Low-count rate by facet level" else as.character(main[1]),
      ylab = "Low-count rate",
      label_angle = label_angle,
      mar_bottom = 9.0
    )
    graphics::abline(h = 0, col = pal["grid"], lty = 1)
  }
  invisible(new_mfrm_plot_data(
    "bias_count",
    list(plot = "lowcount_by_facet", table = rate_tbl, labels = labels)
  ))
}

plot_fixed_reports_bundle <- function(x,
                                      plot_type = c("contrast", "pvalue"),
                                      top_n = 30,
                                      draw = TRUE,
                                      main = NULL,
                                      palette = NULL,
                                      label_angle = 45) {
  plot_type <- match.arg(tolower(as.character(plot_type[1])), c("contrast", "pvalue"))
  top_n <- max(1L, as.integer(top_n))
  pair_tbl <- as.data.frame(x$pairwise_table %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(pair_tbl) == 0) {
    stop("Pairwise table is empty; no plot is available.")
  }
  pal <- resolve_palette(
    palette = palette,
    defaults = c(
      pos = "#1b9e77",
      neg = "#d95f02",
      hist = "#756bb1"
    )
  )

  if (plot_type == "contrast") {
    if (!"Contrast" %in% names(pair_tbl)) {
      stop("Pairwise table does not include `Contrast`.")
    }
    pair_tbl$Contrast <- suppressWarnings(as.numeric(pair_tbl$Contrast))
    pair_tbl <- pair_tbl[is.finite(pair_tbl$Contrast), , drop = FALSE]
    if (nrow(pair_tbl) == 0) stop("No finite contrast values available.")
    pair_tbl <- pair_tbl |>
      dplyr::mutate(.abs = abs(.data$Contrast)) |>
      dplyr::arrange(dplyr::desc(.data$.abs)) |>
      dplyr::slice_head(n = top_n)
    labels <- if (all(c("Target", "Context1", "Context2") %in% names(pair_tbl))) {
      paste0(pair_tbl$Target, ": ", pair_tbl$Context1, " vs ", pair_tbl$Context2)
    } else {
      paste0("Pair ", seq_len(nrow(pair_tbl)))
    }
    if (isTRUE(draw)) {
      barplot_rot45(
        height = pair_tbl$Contrast,
        labels = labels,
        col = ifelse(pair_tbl$Contrast >= 0, pal["pos"], pal["neg"]),
        main = if (is.null(main)) "Pairwise contrasts" else as.character(main[1]),
        ylab = "Contrast (logit)",
        label_angle = label_angle,
        mar_bottom = 9.2
      )
      graphics::abline(h = 0, lty = 2, col = "gray50")
    }
    return(invisible(new_mfrm_plot_data(
      "fixed_reports",
      list(plot = "contrast", table = pair_tbl, labels = labels)
    )))
  }

  p_col <- if ("Prob." %in% names(pair_tbl)) "Prob." else if ("p.value" %in% names(pair_tbl)) "p.value" else NA_character_
  if (is.na(p_col)) {
    stop("Pairwise table does not include p-value column (`Prob.` or `p.value`).")
  }
  p_vals <- suppressWarnings(as.numeric(pair_tbl[[p_col]]))
  p_vals <- p_vals[is.finite(p_vals)]
  if (length(p_vals) == 0) stop("No finite p-values available.")
  if (isTRUE(draw)) {
    graphics::hist(
      x = p_vals,
      breaks = "FD",
      col = pal["hist"],
      border = "white",
      main = if (is.null(main)) "Pairwise p-value distribution" else as.character(main[1]),
      xlab = "p-value",
      ylab = "Count"
    )
    graphics::abline(v = 0.05, lty = 2, col = "gray45")
  }
  invisible(new_mfrm_plot_data(
    "fixed_reports",
    list(plot = "pvalue", p_values = p_vals)
  ))
}

plot_visual_summaries_bundle <- function(x,
                                         plot_type = c("comparison", "warning_counts", "summary_counts"),
                                         draw = TRUE,
                                         main = NULL,
                                         palette = NULL,
                                         label_angle = 45) {
  plot_type <- match.arg(tolower(as.character(plot_type[1])), c("comparison", "warning_counts", "summary_counts"))
  warning_counts <- as.data.frame(x$warning_counts %||% data.frame(), stringsAsFactors = FALSE)
  summary_counts <- as.data.frame(x$summary_counts %||% data.frame(), stringsAsFactors = FALSE)
  pal <- resolve_palette(
    palette = palette,
    defaults = c(
      warning = "#cb181d",
      summary = "#2b8cbe",
      single = "#756bb1"
    )
  )

  if (plot_type == "warning_counts" || plot_type == "summary_counts") {
    tbl <- if (plot_type == "warning_counts") warning_counts else summary_counts
    if (nrow(tbl) == 0 || !all(c("Visual", "Messages") %in% names(tbl))) {
      stop("Requested count table is empty.")
    }
    if (isTRUE(draw)) {
      barplot_rot45(
        height = suppressWarnings(as.numeric(tbl$Messages)),
        labels = as.character(tbl$Visual),
        col = pal["single"],
        main = if (is.null(main)) {
          if (plot_type == "warning_counts") "Warning message counts by visual" else "Summary message counts by visual"
        } else {
          as.character(main[1])
        },
        ylab = "Messages",
        label_angle = label_angle,
        mar_bottom = 8.8
      )
    }
    return(invisible(new_mfrm_plot_data(
      "visual_summaries",
      list(plot = plot_type, table = tbl)
    )))
  }

  vis <- sort(unique(c(as.character(warning_counts$Visual), as.character(summary_counts$Visual))))
  if (length(vis) == 0) {
    stop("No warning/summary counts available.")
  }
  warn <- stats::setNames(rep(0, length(vis)), vis)
  summ <- stats::setNames(rep(0, length(vis)), vis)
  if (nrow(warning_counts) > 0) {
    warn[as.character(warning_counts$Visual)] <- suppressWarnings(as.numeric(warning_counts$Messages))
  }
  if (nrow(summary_counts) > 0) {
    summ[as.character(summary_counts$Visual)] <- suppressWarnings(as.numeric(summary_counts$Messages))
  }
  mat <- rbind(warn, summ)
  rownames(mat) <- c("Warning", "Summary")
  if (isTRUE(draw)) {
    old_mar <- graphics::par("mar")
    on.exit(graphics::par(mar = old_mar), add = TRUE)
    mar <- old_mar
    mar[1] <- max(mar[1], 9.0)
    graphics::par(mar = mar)
    mids <- graphics::barplot(
      height = mat,
      beside = TRUE,
      names.arg = FALSE,
      col = c(pal["warning"], pal["summary"]),
      ylab = "Messages",
      main = if (is.null(main)) "Warning vs summary counts by visual" else as.character(main[1]),
      border = "white"
    )
    centers <- vapply(split(as.numeric(mids), rep(seq_along(vis), each = 2L)), mean, numeric(1))
    draw_rotated_x_labels(
      at = centers,
      labels = vis,
      srt = label_angle,
      cex = 0.82,
      line_offset = 0.085
    )
    graphics::legend(
      "topright",
      legend = c("Warning", "Summary"),
      fill = c(pal["warning"], pal["summary"]),
      bty = "n",
      cex = 0.85
    )
  }
  invisible(new_mfrm_plot_data(
    "visual_summaries",
    list(plot = "comparison", matrix = mat, visuals = vis)
  ))
}

#' Plot report/table bundles with base R defaults
#'
#' @param x A bundle object returned by mfrmr table/report helpers.
#' @param y Reserved for generic compatibility.
#' @param type Optional plot type. Available values depend on bundle class.
#' @param ... Additional arguments forwarded to class-specific plotters.
#'
#' @details
#' `plot()` dispatches by bundle class:
#' - `mfrm_unexpected` -> [plot_unexpected()]
#' - `mfrm_fair_average` -> [plot_fair_average()]
#' - `mfrm_displacement` -> [plot_displacement()]
#' - `mfrm_interrater` -> [plot_interrater_agreement()]
#' - `mfrm_facets_chisq` -> [plot_facets_chisq()]
#' - `mfrm_bias_interaction` -> [plot_bias_interaction()]
#' - `mfrm_bias_count` -> bias-count plots (cell counts / low-count rates)
#' - `mfrm_fixed_reports` -> pairwise-contrast diagnostics
#' - `mfrm_visual_summaries` -> warning/summary message count plots
#' - `mfrm_category_structure` -> default base-R category plots
#' - `mfrm_category_curves` -> default ogive/CCC plots
#' - `mfrm_rating_scale` -> category-counts/threshold plots
#' - `mfrm_measurable` -> measurable-data coverage/count plots
#' - `mfrm_unexpected_after_bias` -> post-bias unexpected-response plots
#' - `mfrm_output_bundle` -> graph/score output-file diagnostics
#' - `mfrm_residual_pca` -> residual PCA scree/loadings via [plot_residual_pca()]
#' - `mfrm_specifications` -> facet/anchor/convergence plots
#' - `mfrm_data_quality` -> row-audit/category/missing-row plots
#' - `mfrm_iteration_report` -> replayed-iteration trajectories
#' - `mfrm_subset_connectivity` -> subset-observation/connectivity plots
#' - `mfrm_facet_statistics` -> facet statistic profile plots
#'
#' If a class is outside these families, use dedicated plotting helpers
#' or custom base R graphics on component tables.
#'
#' @section Interpreting output:
#' The returned object is plotting data (`mfrm_plot_data`) that captures
#' the selected route and payload; set `draw = TRUE` for immediate base graphics.
#'
#' @section Typical workflow:
#' 1. Create bundle output (e.g., `unexpected_response_table()`).
#' 2. Inspect routing with `summary(bundle)` if needed.
#' 3. Call `plot(bundle, type = ..., draw = FALSE)` to obtain reusable plot data.
#'
#' @return A plotting-data object of class `mfrm_plot_data`.
#' @seealso `summary()`, [plot_unexpected()], [plot_fair_average()], [plot_displacement()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' t4 <- unexpected_response_table(fit, abs_z_min = 1.5, prob_max = 0.4, top_n = 10)
#' p <- plot(t4, draw = FALSE)
#' vis <- build_visual_summaries(fit, diagnose_mfrm(fit, residual_pca = "none"))
#' p_vis <- plot(vis, type = "comparison", draw = FALSE)
#' spec <- specifications_report(fit)
#' p_spec <- plot(spec, type = "facet_elements", draw = FALSE)
#' if (interactive()) {
#'   plot(
#'     t4,
#'     type = "severity",
#'     draw = TRUE,
#'     main = "Unexpected Response Severity (Customized)",
#'     palette = c(higher = "#d95f02", lower = "#1b9e77", bar = "#2b8cbe"),
#'     label_angle = 45
#'   )
#'   plot(
#'     vis,
#'     type = "comparison",
#'     draw = TRUE,
#'     main = "Warning vs Summary Counts (Customized)",
#'     palette = c(warning = "#cb181d", summary = "#3182bd"),
#'     label_angle = 45
#'   )
#' }
#' @export
plot.mfrm_bundle <- function(x, y = NULL, type = NULL, ...) {
  dots <- list(...)

  if (inherits(x, "mfrm_unexpected")) {
    args <- c(list(x = x), dots)
    if (!is.null(type)) args$plot_type <- type
    return(do.call(plot_unexpected, args))
  }
  if (inherits(x, "mfrm_fair_average")) {
    args <- c(list(x = x), dots)
    if (!is.null(type)) args$plot_type <- type
    return(do.call(plot_fair_average, args))
  }
  if (inherits(x, "mfrm_displacement")) {
    args <- c(list(x = x), dots)
    if (!is.null(type)) args$plot_type <- type
    return(do.call(plot_displacement, args))
  }
  if (inherits(x, "mfrm_interrater")) {
    args <- c(list(x = x), dots)
    if (!is.null(type)) args$plot_type <- type
    return(do.call(plot_interrater_agreement, args))
  }
  if (inherits(x, "mfrm_facets_chisq")) {
    args <- c(list(x = x), dots)
    if (!is.null(type)) args$plot_type <- type
    return(do.call(plot_facets_chisq, args))
  }
  if (inherits(x, "mfrm_bias_interaction")) {
    args <- c(list(x = x), dots)
    if (!is.null(type)) args$plot <- type
    return(do.call(plot_bias_interaction, args))
  }
  if (inherits(x, "mfrm_bias_count")) {
    args <- c(list(x = x), dots)
    if (!is.null(type)) args$plot_type <- type
    return(do.call(plot_bias_count_bundle, args))
  }
  if (inherits(x, "mfrm_fixed_reports")) {
    args <- c(list(x = x), dots)
    if (!is.null(type)) args$plot_type <- type
    return(do.call(plot_fixed_reports_bundle, args))
  }
  if (inherits(x, "mfrm_visual_summaries")) {
    args <- c(list(x = x), dots)
    if (!is.null(type)) args$plot_type <- type
    return(do.call(plot_visual_summaries_bundle, args))
  }
  if (inherits(x, "mfrm_parity_report")) {
    draw <- if ("draw" %in% names(dots)) isTRUE(dots$draw) else TRUE
    ptype <- if (is.null(type)) "column_coverage" else as.character(type[1])
    top_n <- if ("top_n" %in% names(dots)) dots$top_n else 40L
    main <- dots$main %||% NULL
    palette <- dots$palette %||% NULL
    label_angle <- as.numeric(dots$label_angle %||% 45)
    return(invisible(draw_parity_bundle(
      x,
      type = ptype,
      top_n = top_n,
      draw = draw,
      main = main,
      palette = palette,
      label_angle = label_angle
    )))
  }
  if (inherits(x, "mfrm_category_structure")) {
    draw <- if ("draw" %in% names(dots)) isTRUE(dots$draw) else TRUE
    ptype <- if (is.null(type)) "counts" else as.character(type[1])
    main <- dots$main %||% NULL
    palette <- dots$palette %||% NULL
    label_angle <- as.numeric(dots$label_angle %||% 45)
    return(invisible(draw_category_structure_bundle(
      x,
      type = ptype,
      draw = draw,
      main = main,
      palette = palette,
      label_angle = label_angle
    )))
  }
  if (inherits(x, "mfrm_category_curves")) {
    draw <- if ("draw" %in% names(dots)) isTRUE(dots$draw) else TRUE
    ptype <- if (is.null(type)) "ogive" else as.character(type[1])
    main <- dots$main %||% NULL
    palette <- dots$palette %||% NULL
    return(invisible(draw_category_curves_bundle(
      x,
      type = ptype,
      draw = draw,
      main = main,
      palette = palette
    )))
  }
  if (inherits(x, "mfrm_rating_scale")) {
    draw <- if ("draw" %in% names(dots)) isTRUE(dots$draw) else TRUE
    ptype <- if (is.null(type)) "counts" else as.character(type[1])
    main <- dots$main %||% NULL
    palette <- dots$palette %||% NULL
    label_angle <- as.numeric(dots$label_angle %||% 45)
    return(invisible(draw_rating_scale_bundle(
      x,
      type = ptype,
      draw = draw,
      main = main,
      palette = palette,
      label_angle = label_angle
    )))
  }
  if (inherits(x, "mfrm_measurable")) {
    draw <- if ("draw" %in% names(dots)) isTRUE(dots$draw) else TRUE
    ptype <- if (is.null(type)) "facet_coverage" else as.character(type[1])
    main <- dots$main %||% NULL
    palette <- dots$palette %||% NULL
    label_angle <- as.numeric(dots$label_angle %||% 45)
    return(invisible(draw_measurable_bundle(
      x,
      type = ptype,
      draw = draw,
      main = main,
      palette = palette,
      label_angle = label_angle
    )))
  }
  if (inherits(x, "mfrm_unexpected_after_bias")) {
    draw <- if ("draw" %in% names(dots)) isTRUE(dots$draw) else TRUE
    ptype <- if (is.null(type)) "scatter" else as.character(type[1])
    main <- dots$main %||% NULL
    palette <- dots$palette %||% NULL
    label_angle <- as.numeric(dots$label_angle %||% 45)
    top_n <- if ("top_n" %in% names(dots)) dots$top_n else 40L
    return(invisible(draw_unexpected_after_bias_bundle(
      x,
      type = ptype,
      top_n = top_n,
      draw = draw,
      main = main,
      palette = palette,
      label_angle = label_angle
    )))
  }
  if (inherits(x, "mfrm_output_bundle")) {
    draw <- if ("draw" %in% names(dots)) isTRUE(dots$draw) else TRUE
    ptype <- if (is.null(type)) "graph_expected" else as.character(type[1])
    main <- dots$main %||% NULL
    palette <- dots$palette %||% NULL
    return(invisible(draw_output_bundle(
      x,
      type = ptype,
      draw = draw,
      main = main,
      palette = palette
    )))
  }
  if (inherits(x, "mfrm_residual_pca")) {
    draw <- if ("draw" %in% names(dots)) isTRUE(dots$draw) else TRUE
    ptype <- if (is.null(type)) "overall_scree" else as.character(type[1])
    facet <- dots$facet %||% NULL
    component <- if ("component" %in% names(dots)) dots$component else 1L
    top_n <- if ("top_n" %in% names(dots)) dots$top_n else 20L
    return(invisible(draw_residual_pca_bundle(
      x,
      type = ptype,
      facet = facet,
      component = component,
      top_n = top_n,
      draw = draw
    )))
  }
  if (inherits(x, "mfrm_specifications")) {
    draw <- if ("draw" %in% names(dots)) isTRUE(dots$draw) else TRUE
    ptype <- if (is.null(type)) "facet_elements" else as.character(type[1])
    main <- dots$main %||% NULL
    palette <- dots$palette %||% NULL
    label_angle <- as.numeric(dots$label_angle %||% 45)
    return(invisible(draw_specifications_bundle(
      x,
      type = ptype,
      draw = draw,
      main = main,
      palette = palette,
      label_angle = label_angle
    )))
  }
  if (inherits(x, "mfrm_data_quality")) {
    draw <- if ("draw" %in% names(dots)) isTRUE(dots$draw) else TRUE
    ptype <- if (is.null(type)) "row_audit" else as.character(type[1])
    main <- dots$main %||% NULL
    palette <- dots$palette %||% NULL
    label_angle <- as.numeric(dots$label_angle %||% 45)
    return(invisible(draw_data_quality_bundle(
      x,
      type = ptype,
      draw = draw,
      main = main,
      palette = palette,
      label_angle = label_angle
    )))
  }
  if (inherits(x, "mfrm_iteration_report")) {
    draw <- if ("draw" %in% names(dots)) isTRUE(dots$draw) else TRUE
    ptype <- if (is.null(type)) "residual" else as.character(type[1])
    main <- dots$main %||% NULL
    palette <- dots$palette %||% NULL
    return(invisible(draw_iteration_report_bundle(
      x,
      type = ptype,
      draw = draw,
      main = main,
      palette = palette
    )))
  }
  if (inherits(x, "mfrm_subset_connectivity")) {
    draw <- if ("draw" %in% names(dots)) isTRUE(dots$draw) else TRUE
    ptype <- if (is.null(type)) "subset_observations" else as.character(type[1])
    main <- dots$main %||% NULL
    palette <- dots$palette %||% NULL
    label_angle <- as.numeric(dots$label_angle %||% 45)
    return(invisible(draw_subset_connectivity_bundle(
      x,
      type = ptype,
      draw = draw,
      main = main,
      palette = palette,
      label_angle = label_angle
    )))
  }
  if (inherits(x, "mfrm_facet_statistics")) {
    draw <- if ("draw" %in% names(dots)) isTRUE(dots$draw) else TRUE
    ptype <- if (is.null(type)) "means" else as.character(type[1])
    metric <- dots$metric %||% NULL
    main <- dots$main %||% NULL
    palette <- dots$palette %||% NULL
    label_angle <- as.numeric(dots$label_angle %||% 45)
    return(invisible(draw_facet_statistics_bundle(
      x,
      type = ptype,
      metric = metric,
      draw = draw,
      main = main,
      palette = palette,
      label_angle = label_angle
    )))
  }

  stop(
    "No default plot method for class `", class(x)[1], "`.\n",
    "Use a dedicated plot helper (for example, `plot_unexpected()`, `plot_fair_average()`, or `plot_bias_interaction()`)."
  )
}

#' Summarize an `mfrm_diagnostics` object in a user-friendly format
#'
#' @param object Output from [diagnose_mfrm()].
#' @param digits Number of digits for printed numeric values.
#' @param top_n Number of highest-absolute-Z fit rows to keep.
#' @param ... Reserved for generic compatibility.
#'
#' @details
#' This method returns a compact diagnostics summary designed for quick review:
#' - design overview (observations, persons, facets, categories, subsets)
#' - global fit statistics
#' - reliability/separation by facet
#' - top facet/person fit rows by absolute ZSTD
#' - counts of flagged diagnostics (unexpected, displacement, interactions)
#'
#' @section Interpreting output:
#' - `overview`: analysis scale, subset count, and residual-PCA mode.
#' - `overall_fit`: global fit indices.
#' - `reliability`: facet separation/reliability block.
#' - `top_fit`: highest `|ZSTD|` elements for immediate inspection.
#' - `flags`: compact counts for key warning domains.
#'
#' @section Typical workflow:
#' 1. Run diagnostics with [diagnose_mfrm()].
#' 2. Review `summary(diag)` for major warnings.
#' 3. Follow up with dedicated tables/plots for flagged domains.
#'
#' @return An object of class `summary.mfrm_diagnostics` with:
#' - `overview`: design-level counts and residual-PCA mode
#' - `overall_fit`: global fit block
#' - `reliability`: facet-level separation/reliability summary
#' - `top_fit`: top `|ZSTD|` rows
#' - `flags`: compact flag counts for major diagnostics
#' - `notes`: short interpretation notes
#' @seealso [diagnose_mfrm()], [summary.mfrm_fit()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#' summary(diag)
#' @export
summary.mfrm_diagnostics <- function(object, digits = 3, top_n = 10, ...) {
  if (!is.list(object) || is.null(object$obs)) {
    stop("`object` must be output from diagnose_mfrm().")
  }

  digits <- max(0L, as.integer(digits))
  top_n <- max(1L, as.integer(top_n))

  obs_tbl <- tibble::as_tibble(object$obs)
  fit_tbl <- tibble::as_tibble(object$fit %||% tibble::tibble())
  reliability_tbl <- tibble::as_tibble(object$reliability %||% tibble::tibble())
  overall_fit <- tibble::as_tibble(object$overall_fit %||% tibble::tibble())
  subset_summary <- tibble::as_tibble(object$subsets$summary %||% tibble::tibble())

  n_obs <- nrow(obs_tbl)
  n_person <- if ("Person" %in% names(obs_tbl)) dplyr::n_distinct(obs_tbl$Person) else NA_integer_
  n_cat <- if ("Observed" %in% names(obs_tbl)) dplyr::n_distinct(obs_tbl$Observed) else NA_integer_
  n_subsets <- if ("Subset" %in% names(subset_summary)) dplyr::n_distinct(subset_summary$Subset) else 0L

  overview <- tibble::tibble(
    Observations = n_obs,
    Persons = n_person,
    Facets = length(object$facet_names %||% character(0)),
    Categories = n_cat,
    Subsets = n_subsets,
    ResidualPCA = as.character(object$residual_pca_mode %||% "none")
  )

  reliability_overview <- tibble::tibble()
  keep_rel <- c("Facet", "Levels", "Separation", "Strata", "Reliability", "MeanInfit", "MeanOutfit")
  if (nrow(reliability_tbl) > 0) {
    keep <- intersect(keep_rel, names(reliability_tbl))
    reliability_overview <- reliability_tbl |>
      dplyr::select(dplyr::all_of(keep)) |>
      dplyr::arrange(.data$Facet)
  }

  top_fit <- tibble::tibble()
  fit_need <- c("Facet", "Level", "Infit", "Outfit", "InfitZSTD", "OutfitZSTD")
  if (nrow(fit_tbl) > 0 && all(fit_need %in% names(fit_tbl))) {
    top_fit <- fit_tbl |>
      dplyr::mutate(
        AbsZ = pmax(abs(.data$InfitZSTD), abs(.data$OutfitZSTD), na.rm = TRUE)
      ) |>
      dplyr::arrange(dplyr::desc(.data$AbsZ)) |>
      dplyr::slice_head(n = top_n) |>
      dplyr::select("Facet", "Level", "Infit", "Outfit", "InfitZSTD", "OutfitZSTD", "AbsZ")
  }

  unexpected_n <- suppressWarnings(as.integer(object$unexpected$summary$UnexpectedN[1] %||% NA_integer_))
  displacement_flagged <- suppressWarnings(as.integer(object$displacement$summary$FlaggedLevels[1] %||% NA_integer_))
  interaction_n <- if (!is.null(object$interactions)) nrow(object$interactions) else NA_integer_
  interrater_pairs <- suppressWarnings(as.integer(object$interrater$summary$Pairs[1] %||% NA_integer_))

  flags <- tibble::tibble(
    Metric = c(
      "Unexpected responses",
      "Flagged displacement levels",
      "Interaction rows",
      "Inter-rater pairs"
    ),
    Count = c(unexpected_n, displacement_flagged, interaction_n, interrater_pairs)
  )

  notes <- character(0)
  if (isTRUE(n_subsets > 1L)) {
    notes <- c(notes, "Multiple disconnected subsets were detected.")
  }
  if (isTRUE(!is.na(unexpected_n) && unexpected_n > 0L)) {
    notes <- c(notes, "Unexpected responses were flagged under current thresholds.")
  }
  if (length(notes) == 0) {
    notes <- "No immediate warnings from diagnostics summary."
  }

  out <- list(
    overview = overview,
    overall_fit = overall_fit,
    reliability = reliability_overview,
    top_fit = top_fit,
    flags = flags,
    notes = notes,
    digits = digits
  )
  class(out) <- "summary.mfrm_diagnostics"
  out
}

#' @export
print.summary.mfrm_diagnostics <- function(x, ...) {
  digits <- as.integer(x$digits %||% 3L)
  if (!is.finite(digits)) digits <- 3L

  cat("Many-Facet Rasch Diagnostics Summary\n")
  if (!is.null(x$overview) && nrow(x$overview) > 0) {
    ov <- round_numeric_df(as.data.frame(x$overview), digits = digits)[1, , drop = FALSE]
    cat(sprintf(
      "  Observations: %s | Persons: %s | Facets: %s | Categories: %s | Subsets: %s\n",
      ov$Observations, ov$Persons, ov$Facets, ov$Categories, ov$Subsets
    ))
    cat(sprintf("  Residual PCA mode: %s\n", ov$ResidualPCA))
  }

  if (!is.null(x$overall_fit) && nrow(x$overall_fit) > 0) {
    cat("\nOverall fit\n")
    print(round_numeric_df(as.data.frame(x$overall_fit), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$reliability) && nrow(x$reliability) > 0) {
    cat("\nReliability by facet\n")
    print(round_numeric_df(as.data.frame(x$reliability), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$top_fit) && nrow(x$top_fit) > 0) {
    cat("\nLargest |ZSTD| rows\n")
    print(round_numeric_df(as.data.frame(x$top_fit), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$flags) && nrow(x$flags) > 0) {
    cat("\nFlag counts\n")
    print(as.data.frame(x$flags), row.names = FALSE)
  }

  if (length(x$notes) > 0) {
    cat("\nNotes\n")
    for (line in x$notes) cat(" - ", line, "\n", sep = "")
  }
  invisible(x)
}

#' Summarize an `mfrm_bias` object in a user-friendly format
#'
#' @param object Output from [estimate_bias()].
#' @param digits Number of digits for printed numeric values.
#' @param top_n Number of strongest bias rows to keep.
#' @param p_cut Significance cutoff used for counting flagged rows.
#' @param ... Reserved for generic compatibility.
#'
#' @details
#' This method returns a compact interaction-bias summary:
#' - interaction facets/order and analyzed cell counts
#' - effect-size profile (`|bias|` mean/max, significant cell count)
#' - fixed-effect chi-square block
#' - iteration-end convergence indicators
#' - top rows ranked by absolute t
#'
#' @section Interpreting output:
#' - `overview`: interaction order, analyzed cells, and effect-size profile.
#' - `chi_sq`: fixed-effect test block.
#' - `final_iteration`: end-of-loop status from the bias routine.
#' - `top_rows`: strongest bias contrasts by `|t|`.
#'
#' @section Typical workflow:
#' 1. Estimate interactions with [estimate_bias()].
#' 2. Check `summary(bias)` for significant and unstable cells.
#' 3. Use [bias_interaction_report()] or [plot_bias_interaction()] for details.
#'
#' @return An object of class `summary.mfrm_bias` with:
#' - `overview`: interaction facets/order, cell counts, and effect-size profile
#' - `chi_sq`: fixed-effect chi-square block
#' - `final_iteration`: end-of-iteration status row
#' - `top_rows`: highest-`|t|` interaction rows
#' - `notes`: short interpretation notes
#' @seealso [estimate_bias()], [bias_interaction_report()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#' bias <- estimate_bias(fit, diag, facet_a = "Rater", facet_b = "Criterion", max_iter = 2)
#' summary(bias)
#' @export
summary.mfrm_bias <- function(object, digits = 3, top_n = 10, p_cut = 0.05, ...) {
  if (!is.list(object) || is.null(object$table) || nrow(object$table) == 0) {
    stop("`object` must be non-empty output from estimate_bias().")
  }

  digits <- max(0L, as.integer(digits))
  top_n <- max(1L, as.integer(top_n))
  p_cut <- max(0, min(1, as.numeric(p_cut[1])))

  bias_tbl <- tibble::as_tibble(object$table)
  chi_tbl <- tibble::as_tibble(object$chi_sq %||% tibble::tibble())
  iter_tbl <- tibble::as_tibble(object$iteration %||% tibble::tibble())
  spec <- extract_bias_facet_spec(object)
  interaction_facets <- if (!is.null(spec)) spec$facets else unique(c(
    as.character(object$facet_a[1] %||% NA_character_),
    as.character(object$facet_b[1] %||% NA_character_)
  ))
  interaction_facets <- interaction_facets[!is.na(interaction_facets) & nzchar(interaction_facets)]
  if (length(interaction_facets) == 0) interaction_facets <- c("Unknown")
  interaction_label <- paste(interaction_facets, collapse = " x ")
  interaction_order <- length(interaction_facets)
  interaction_mode <- ifelse(interaction_order > 2, "higher_order", "pairwise")

  abs_bias <- abs(suppressWarnings(as.numeric(bias_tbl$`Bias Size`)))
  p_vals <- suppressWarnings(as.numeric(bias_tbl$`Prob.`))
  sig_n <- sum(is.finite(p_vals) & p_vals <= p_cut, na.rm = TRUE)

  overview <- tibble::tibble(
    FacetPair = interaction_label,
    InteractionOrder = interaction_order,
    InteractionMode = interaction_mode,
    Cells = nrow(bias_tbl),
    MeanAbsBias = mean(abs_bias, na.rm = TRUE),
    MaxAbsBias = max(abs_bias, na.rm = TRUE),
    Significant = sig_n,
    SignificantCut = p_cut
  )

  final_iteration <- tibble::tibble()
  if (nrow(iter_tbl) > 0) {
    final_iteration <- iter_tbl |>
      dplyr::slice_tail(n = 1)
  }

  top_rows <- tibble::tibble()
  level_cols <- if (!is.null(spec)) {
    spec$level_cols
  } else if (all(c("FacetA_Level", "FacetB_Level") %in% names(bias_tbl))) {
    c("FacetA_Level", "FacetB_Level")
  } else {
    character(0)
  }
  keep <- c(level_cols, "Bias Size", "S.E.", "t", "Prob.", "Obs-Exp Average")
  if (all(keep %in% names(bias_tbl))) {
    top_rows <- bias_tbl |>
      dplyr::mutate(AbsT = abs(.data$t)) |>
      dplyr::arrange(dplyr::desc(.data$AbsT)) |>
      dplyr::slice_head(n = top_n) |>
      dplyr::select(dplyr::all_of(c(keep, "AbsT")))
    if (length(level_cols) == length(interaction_facets)) {
      names(top_rows)[seq_along(level_cols)] <- interaction_facets
      top_rows <- dplyr::mutate(
        top_rows,
        Pair = do.call(paste, c(top_rows[interaction_facets], sep = " | ")),
        .before = 1
      )
    }
  }

  notes <- character(0)
  if (nrow(iter_tbl) > 0) {
    tail_row <- iter_tbl[nrow(iter_tbl), , drop = FALSE]
    tail_cells <- suppressWarnings(as.numeric(tail_row$BiasCells[1]))
    if (is.finite(tail_cells) && tail_cells > 0) {
      notes <- c(notes, "Bias iteration may not have fully stabilized (BiasCells > 0 at final step).")
    }
  }
  if (interaction_order > 2) {
    notes <- c(notes, "Higher-order interaction mode is active; pairwise contrasts should be interpreted from dedicated 2-way runs.")
  }
  if (length(notes) == 0) {
    notes <- "No immediate warnings from bias summary."
  }

  out <- list(
    overview = overview,
    chi_sq = chi_tbl,
    final_iteration = final_iteration,
    top_rows = top_rows,
    notes = notes,
    digits = digits
  )
  class(out) <- "summary.mfrm_bias"
  out
}

#' @export
print.summary.mfrm_bias <- function(x, ...) {
  digits <- as.integer(x$digits %||% 3L)
  if (!is.finite(digits)) digits <- 3L

  cat("Many-Facet Rasch Bias Summary\n")
  if (!is.null(x$overview) && nrow(x$overview) > 0) {
    ov <- round_numeric_df(as.data.frame(x$overview), digits = digits)[1, , drop = FALSE]
    cat(sprintf("  Interaction facets: %s | Cells: %s\n", ov$FacetPair, ov$Cells))
    if ("InteractionOrder" %in% names(ov) && "InteractionMode" %in% names(ov)) {
      cat(sprintf("  Order: %s | Mode: %s\n", ov$InteractionOrder, ov$InteractionMode))
    }
    cat(sprintf(
      "  Mean |Bias|: %s | Max |Bias|: %s | Significant (p <= %.3f): %s\n",
      ov$MeanAbsBias, ov$MaxAbsBias, as.numeric(ov$SignificantCut), ov$Significant
    ))
  }

  if (!is.null(x$chi_sq) && nrow(x$chi_sq) > 0) {
    cat("\nFixed-effect chi-square\n")
    print(round_numeric_df(as.data.frame(x$chi_sq), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$final_iteration) && nrow(x$final_iteration) > 0) {
    cat("\nFinal iteration status\n")
    print(round_numeric_df(as.data.frame(x$final_iteration), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$top_rows) && nrow(x$top_rows) > 0) {
    cat("\nTop |t| bias rows\n")
    print(round_numeric_df(as.data.frame(x$top_rows), digits = digits), row.names = FALSE)
  }

  if (length(x$notes) > 0) {
    cat("\nNotes\n")
    for (line in x$notes) cat(" - ", line, "\n", sep = "")
  }
  invisible(x)
}

#' Summarize an `mfrm_fit` object in a user-friendly format
#'
#' @param object Output from [fit_mfrm()].
#' @param digits Number of digits for printed numeric values.
#' @param top_n Number of extreme facet/person rows shown in summaries.
#' @param ... Reserved for generic compatibility.
#'
#' @details
#' This method provides a compact, human-readable summary oriented to reporting.
#' It returns a structured object and prints:
#' - model fit overview (N, LogLik, AIC/BIC, convergence)
#' - facet-level estimate distribution (mean/SD/range)
#' - person measure distribution
#' - step/threshold checks
#' - high/low person measures and extreme facet levels
#'
#' @section Interpreting output:
#' - `overview`: convergence and information criteria.
#' - `facet_overview`: per-facet spread and range of estimates.
#' - `person_overview`: distribution of person measures.
#' - `step_overview`: threshold spread and monotonicity checks.
#' - `top_person` / `top_facet`: extreme estimates for quick triage.
#'
#' @section Typical workflow:
#' 1. Fit model with [fit_mfrm()].
#' 2. Run `summary(fit)` for first-pass diagnostics.
#' 3. Continue with [diagnose_mfrm()] for element-level fit checks.
#'
#' @return An object of class `summary.mfrm_fit` with:
#' - `overview`: global model/fit indicators
#' - `facet_overview`: per-facet estimate distribution summary
#' - `person_overview`: person-measure distribution summary
#' - `step_overview`: threshold/step diagnostics
#' - `top_person`: highest/lowest person measures
#' - `top_facet`: extreme facet-level estimates
#' - `notes`: short interpretation notes
#' @seealso [fit_mfrm()], [diagnose_mfrm()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' summary(fit)
#' @export
summary.mfrm_fit <- function(object, digits = 3, top_n = 5, ...) {
  if (is.null(object$summary) || nrow(object$summary) == 0) {
    stop("`object` does not contain fit summary information.")
  }

  digits <- max(0L, as.integer(digits))
  top_n <- max(1L, as.integer(top_n))

  overview <- tibble::as_tibble(object$summary)
  person_raw <- object$facets$person
  if (is.null(person_raw)) person_raw <- tibble::tibble()
  other_raw <- object$facets$others
  if (is.null(other_raw)) other_raw <- tibble::tibble()
  step_raw <- object$steps
  if (is.null(step_raw)) step_raw <- tibble::tibble()

  person_tbl <- tibble::as_tibble(person_raw)
  other_tbl <- tibble::as_tibble(other_raw)
  step_tbl <- tibble::as_tibble(step_raw)

  facet_overview <- tibble::tibble()
  if (nrow(other_tbl) > 0 && all(c("Facet", "Estimate") %in% names(other_tbl))) {
    facet_overview <- other_tbl |>
      dplyr::group_by(.data$Facet) |>
      dplyr::summarise(
        Levels = dplyr::n(),
        MeanEstimate = mean(.data$Estimate, na.rm = TRUE),
        SDEstimate = stats::sd(.data$Estimate, na.rm = TRUE),
        MinEstimate = min(.data$Estimate, na.rm = TRUE),
        MaxEstimate = max(.data$Estimate, na.rm = TRUE),
        Span = .data$MaxEstimate - .data$MinEstimate,
        .groups = "drop"
      ) |>
      dplyr::arrange(.data$Facet)
  }

  person_overview <- tibble::tibble()
  if (nrow(person_tbl) > 0 && "Estimate" %in% names(person_tbl)) {
    person_overview <- tibble::tibble(
      Persons = nrow(person_tbl),
      Mean = mean(person_tbl$Estimate, na.rm = TRUE),
      SD = stats::sd(person_tbl$Estimate, na.rm = TRUE),
      Median = stats::median(person_tbl$Estimate, na.rm = TRUE),
      Min = min(person_tbl$Estimate, na.rm = TRUE),
      Max = max(person_tbl$Estimate, na.rm = TRUE),
      Span = max(person_tbl$Estimate, na.rm = TRUE) - min(person_tbl$Estimate, na.rm = TRUE)
    )

    if ("SD" %in% names(person_tbl)) {
      person_overview$MeanPosteriorSD <- mean(person_tbl$SD, na.rm = TRUE)
    }
  }

  step_overview <- tibble::tibble()
  if (nrow(step_tbl) > 0 && all(c("Step", "Estimate") %in% names(step_tbl))) {
    ord <- order(step_tbl$Step)
    step_vals <- as.numeric(step_tbl$Estimate[ord])
    monotonic <- if (length(step_vals) <= 1) TRUE else all(diff(step_vals) >= -sqrt(.Machine$double.eps))
    step_overview <- tibble::tibble(
      Steps = nrow(step_tbl),
      Min = min(step_tbl$Estimate, na.rm = TRUE),
      Max = max(step_tbl$Estimate, na.rm = TRUE),
      Span = max(step_tbl$Estimate, na.rm = TRUE) - min(step_tbl$Estimate, na.rm = TRUE),
      Monotonic = monotonic
    )
  }

  facet_extremes <- tibble::tibble()
  if (nrow(other_tbl) > 0 && all(c("Facet", "Level", "Estimate") %in% names(other_tbl))) {
    facet_extremes <- other_tbl |>
      dplyr::mutate(AbsEstimate = abs(.data$Estimate)) |>
      dplyr::arrange(dplyr::desc(.data$AbsEstimate)) |>
      dplyr::slice_head(n = top_n) |>
      dplyr::select("Facet", "Level", "Estimate")
  }

  person_high <- tibble::tibble()
  person_low <- tibble::tibble()
  if (nrow(person_tbl) > 0 && all(c("Person", "Estimate") %in% names(person_tbl))) {
    person_high <- person_tbl |>
      dplyr::arrange(dplyr::desc(.data$Estimate)) |>
      dplyr::slice_head(n = top_n)
    person_low <- person_tbl |>
      dplyr::arrange(.data$Estimate) |>
      dplyr::slice_head(n = top_n)
  }

  notes <- character(0)
  if ("Converged" %in% names(overview) && !isTRUE(overview$Converged[1])) {
    notes <- c(notes, "Optimization did not converge; interpret parameter estimates cautiously.")
  }
  if (nrow(step_overview) > 0 && !isTRUE(step_overview$Monotonic[1])) {
    notes <- c(notes, "Step estimates are not monotonic; verify category functioning.")
  }
  if (length(notes) == 0) {
    notes <- "No immediate warnings from fit-level summary checks."
  }

  out <- list(
    overview = overview,
    facet_overview = facet_overview,
    person_overview = person_overview,
    step_overview = step_overview,
    facet_extremes = facet_extremes,
    person_high = person_high,
    person_low = person_low,
    notes = notes,
    digits = digits
  )
  class(out) <- "summary.mfrm_fit"
  out
}

round_numeric_df <- function(df, digits = 3L) {
  if (!is.data.frame(df) || nrow(df) == 0) return(df)
  out <- df
  numeric_cols <- vapply(out, is.numeric, logical(1))
  out[numeric_cols] <- lapply(out[numeric_cols], round, digits = digits)
  out
}

#' @export
print.summary.mfrm_fit <- function(x, ...) {
  digits <- x$digits
  if (is.null(digits) || !is.finite(digits)) digits <- 3L
  overview <- round_numeric_df(as.data.frame(x$overview), digits = digits)

  cat("Many-Facet Rasch Model Summary\n")
  if (nrow(overview) > 0) {
    ov <- overview[1, , drop = FALSE]
    cat(sprintf("  Model: %s | Method: %s\n", ov$Model, ov$Method))
    cat(sprintf("  N: %s | Persons: %s | Facets: %s | Categories: %s\n", ov$N, ov$Persons, ov$Facets, ov$Categories))
    cat(sprintf("  LogLik: %s | AIC: %s | BIC: %s\n", ov$LogLik, ov$AIC, ov$BIC))
    cat(sprintf("  Converged: %s | Iterations: %s\n", ifelse(isTRUE(ov$Converged), "Yes", "No"), ov$Iterations))
  }

  if (nrow(x$facet_overview) > 0) {
    cat("\nFacet overview\n")
    print(round_numeric_df(as.data.frame(x$facet_overview), digits = digits), row.names = FALSE)
  }

  if (nrow(x$person_overview) > 0) {
    cat("\nPerson measure distribution\n")
    print(round_numeric_df(as.data.frame(x$person_overview), digits = digits), row.names = FALSE)
  }

  if (nrow(x$step_overview) > 0) {
    cat("\nStep parameter summary\n")
    print(round_numeric_df(as.data.frame(x$step_overview), digits = digits), row.names = FALSE)
  }

  if (nrow(x$facet_extremes) > 0) {
    cat("\nMost extreme facet levels (|estimate|)\n")
    print(round_numeric_df(as.data.frame(x$facet_extremes), digits = digits), row.names = FALSE)
  }

  if (nrow(x$person_high) > 0) {
    cat("\nHighest person measures\n")
    print(round_numeric_df(as.data.frame(x$person_high), digits = digits), row.names = FALSE)
  }

  if (nrow(x$person_low) > 0) {
    cat("\nLowest person measures\n")
    print(round_numeric_df(as.data.frame(x$person_low), digits = digits), row.names = FALSE)
  }

  if (length(x$notes) > 0) {
    cat("\nNotes\n")
    for (line in x$notes) cat(" - ", line, "\n", sep = "")
  }

  invisible(x)
}

# Plot helpers for `plot.mfrm_fit()`.
step_index_from_label <- function(step_labels) {
  step_labels <- as.character(step_labels)
  idx <- suppressWarnings(as.integer(gsub("[^0-9]+", "", step_labels)))
  if (all(is.na(idx))) return(seq_along(step_labels))
  if (any(is.na(idx))) {
    fill_start <- max(idx, na.rm = TRUE)
    idx[is.na(idx)] <- fill_start + seq_len(sum(is.na(idx)))
  }
  idx
}

build_step_curve_spec <- function(x) {
  step_tbl <- x$steps
  if (is.null(step_tbl) || nrow(step_tbl) == 0 || !"Estimate" %in% names(step_tbl)) {
    stop("Step estimates are required for pathway/CCC plots.")
  }
  step_tbl <- tibble::as_tibble(step_tbl)

  model <- toupper(as.character(x$config$model[1]))
  rating_min <- suppressWarnings(min(as.numeric(x$prep$data$Score), na.rm = TRUE))
  if (!is.finite(rating_min)) rating_min <- 0

  n_cat <- suppressWarnings(as.integer(x$config$n_cat[1]))
  if (!is.finite(n_cat) || n_cat < 2) {
    n_cat <- max(2L, nrow(step_tbl) + 1L)
  }
  categories <- rating_min + 0:(n_cat - 1L)

  groups <- list()
  step_points <- tibble::tibble()

  if (model == "RSM") {
    if (!"Step" %in% names(step_tbl)) {
      step_tbl$Step <- paste0("Step_", seq_len(nrow(step_tbl)))
    }
    ord <- order(step_index_from_label(step_tbl$Step))
    tau <- as.numeric(step_tbl$Estimate[ord])
    groups[["Common"]] <- list(name = "Common", step_cum = c(0, cumsum(tau)), tau = tau)
    step_points <- tibble::tibble(
      CurveGroup = "Common",
      Step = as.character(step_tbl$Step[ord]),
      StepIndex = seq_along(tau),
      Threshold = tau
    )
  } else {
    if (!all(c("StepFacet", "Step", "Estimate") %in% names(step_tbl))) {
      stop("PCM step table must include StepFacet, Step, and Estimate columns.")
    }
    step_facet <- as.character(x$config$step_facet[1])
    ordered_levels <- x$prep$levels[[step_facet]]
    if (is.null(ordered_levels)) {
      ordered_levels <- unique(as.character(step_tbl$StepFacet))
    }
    for (lvl in ordered_levels) {
      sub <- step_tbl[as.character(step_tbl$StepFacet) == as.character(lvl), , drop = FALSE]
      if (nrow(sub) == 0) next
      ord <- order(step_index_from_label(sub$Step))
      tau <- as.numeric(sub$Estimate[ord])
      groups[[as.character(lvl)]] <- list(
        name = as.character(lvl),
        step_cum = c(0, cumsum(tau)),
        tau = tau
      )
      step_points <- dplyr::bind_rows(
        step_points,
        tibble::tibble(
          CurveGroup = as.character(lvl),
          Step = as.character(sub$Step[ord]),
          StepIndex = seq_along(tau),
          Threshold = tau
        )
      )
    }
  }

  if (length(groups) == 0) stop("No step groups available for pathway/CCC plots.")

  list(
    model = model,
    categories = categories,
    rating_min = rating_min,
    n_cat = n_cat,
    groups = groups,
    step_points = step_points
  )
}

build_curve_tables <- function(curve_spec, theta_grid) {
  prob_tables <- list()
  exp_tables <- list()
  idx_prob <- 1L
  idx_exp <- 1L
  for (g in names(curve_spec$groups)) {
    grp <- curve_spec$groups[[g]]
    probs <- category_prob_rsm(theta_grid, grp$step_cum)
    k_vals <- as.numeric(curve_spec$categories)
    expected <- as.numeric(probs %*% matrix(k_vals, ncol = 1))
    exp_tables[[idx_exp]] <- tibble::tibble(
      Theta = theta_grid,
      ExpectedScore = expected,
      CurveGroup = grp$name
    )
    idx_exp <- idx_exp + 1L
    for (k in seq_len(ncol(probs))) {
      prob_tables[[idx_prob]] <- tibble::tibble(
        Theta = theta_grid,
        Probability = probs[, k],
        Category = as.character(curve_spec$categories[k]),
        CurveGroup = grp$name
      )
      idx_prob <- idx_prob + 1L
    }
  }
  list(
    probabilities = dplyr::bind_rows(prob_tables),
    expected = dplyr::bind_rows(exp_tables)
  )
}

compute_se_for_plot <- function(x, ci_level = 0.95) {
  tryCatch({
    obs_df <- compute_obs_table(x)
    facet_cols <- x$config$source_columns$facets
    if (is.null(facet_cols) || length(facet_cols) == 0) {
      facet_cols <- x$config$facet_names
    }
    se_tbl <- calc_facet_se(obs_df, facet_cols)
    z <- stats::qnorm(1 - (1 - ci_level) / 2)
    se_tbl$CI_Lower <- -z * se_tbl$SE
    se_tbl$CI_Upper <- z * se_tbl$SE
    se_tbl$CI_Level <- ci_level
    as.data.frame(se_tbl, stringsAsFactors = FALSE)
  }, error = function(e) NULL)
}

build_wright_map_data <- function(x, top_n = 30L, se_tbl = NULL) {
  person_tbl <- tibble::as_tibble(x$facets$person)
  person_tbl <- person_tbl[is.finite(person_tbl$Estimate), , drop = FALSE]
  if (nrow(person_tbl) == 0) stop("Person estimates are not available for Wright map.")

  facet_tbl <- tibble::as_tibble(x$facets$others)
  if (!all(c("Facet", "Level", "Estimate") %in% names(facet_tbl))) {
    stop("Facet-level estimates are not available for Wright map.")
  }
  facet_tbl <- facet_tbl[is.finite(facet_tbl$Estimate), , drop = FALSE]

  step_points <- tryCatch(build_step_curve_spec(x)$step_points, error = function(e) tibble::tibble())
  step_tbl <- if (nrow(step_points) > 0) {
    tibble::tibble(
      PlotType = "Step threshold",
      Group = if ("CurveGroup" %in% names(step_points)) paste0("Step:", step_points$CurveGroup) else "Step",
      Label = step_points$Step,
      Estimate = step_points$Threshold
    )
  } else {
    tibble::tibble()
  }

  facet_pts <- facet_tbl |>
    dplyr::transmute(
      PlotType = "Facet level",
      Group = .data$Facet,
      Label = .data$Level,
      Estimate = .data$Estimate
    )
  if (!is.null(se_tbl) && is.data.frame(se_tbl) && nrow(se_tbl) > 0 &&
      all(c("Facet", "Level", "SE") %in% names(se_tbl))) {
    se_join <- se_tbl[, intersect(c("Facet", "Level", "SE", "CI_Lower", "CI_Upper"), names(se_tbl)), drop = FALSE]
    se_join$Level <- as.character(se_join$Level)
    facet_pts$Group <- as.character(facet_pts$Group)
    facet_pts$Label <- as.character(facet_pts$Label)
    facet_pts <- merge(facet_pts, se_join,
                       by.x = c("Group", "Label"), by.y = c("Facet", "Level"),
                       all.x = TRUE, sort = FALSE)
  }
  point_tbl <- dplyr::bind_rows(facet_pts, step_tbl)
  if (nrow(point_tbl) == 0) stop("No facet/step locations available for Wright map.")

  if (nrow(point_tbl) > top_n) {
    point_tbl <- point_tbl |>
      dplyr::mutate(.Abs = abs(.data$Estimate)) |>
      dplyr::arrange(dplyr::desc(.data$.Abs)) |>
      dplyr::slice_head(n = top_n) |>
      dplyr::select(-.data$.Abs)
  }

  group_levels <- unique(point_tbl$Group)
  point_tbl <- point_tbl |>
    dplyr::mutate(Group = factor(.data$Group, levels = group_levels)) |>
    dplyr::group_by(.data$Group) |>
    dplyr::arrange(.data$Estimate, .by_group = TRUE) |>
    dplyr::mutate(
      XBase = as.numeric(.data$Group),
      X = if (dplyr::n() == 1) .data$XBase else .data$XBase + seq(-0.18, 0.18, length.out = dplyr::n())
    ) |>
    dplyr::ungroup()

  list(
    title = "Wright Map",
    person = person_tbl,
    locations = point_tbl,
    group_levels = group_levels
  )
}

draw_wright_map <- function(plot_data,
                            title = NULL,
                            palette = NULL,
                            label_angle = 45,
                            show_ci = FALSE,
                            ci_level = 0.95) {
  pal <- resolve_palette(
    palette = palette,
    defaults = c(
      facet_level = "#1b9e77",
      step_threshold = "#d95f02",
      person_hist = "gray80",
      grid = "#ececec"
    )
  )
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::layout(matrix(c(1, 2), nrow = 1), widths = c(1.05, 2.35))
  graphics::par(mgp = c(2.6, 0.85, 0))

  graphics::par(mar = c(5.2, 5, 3.2, 1))
  graphics::hist(
    x = plot_data$person$Estimate,
    breaks = "FD",
    main = "Person distribution",
    xlab = "Logit scale",
    ylab = "Count",
    col = pal["person_hist"],
    border = "white"
  )

  loc <- plot_data$locations
  xr <- c(0.5, length(plot_data$group_levels) + 0.5)
  yr <- range(c(loc$Estimate, plot_data$person$Estimate), finite = TRUE)
  if (!all(is.finite(yr)) || diff(yr) <= sqrt(.Machine$double.eps)) {
    center <- mean(c(loc$Estimate, plot_data$person$Estimate), na.rm = TRUE)
    yr <- center + c(-0.5, 0.5)
  }
  graphics::par(mar = c(8.8, 5.2, 3.2, 1.2))
  graphics::plot(
    x = loc$X,
    y = loc$Estimate,
    type = "n",
    xlim = xr,
    ylim = yr,
    xaxt = "n",
    xlab = "",
    ylab = "Logit scale",
    main = if (is.null(title)) "Facet and step locations" else as.character(title[1])
  )
  pretty_y <- pretty(yr, n = 6)
  graphics::abline(h = pretty_y, col = pal["grid"], lty = 1)
  x_at <- seq_along(plot_data$group_levels)
  x_lab <- truncate_axis_label(plot_data$group_levels, width = 16L)
  draw_rotated_x_labels(
    at = x_at,
    labels = x_lab,
    srt = label_angle,
    cex = 0.82,
    line_offset = 0.085
  )
  cols <- ifelse(loc$PlotType == "Step threshold", pal["step_threshold"], pal["facet_level"])
  pch <- ifelse(loc$PlotType == "Step threshold", 17, 16)
  if (isTRUE(show_ci) && "SE" %in% names(loc)) {
    z <- stats::qnorm(1 - (1 - ci_level) / 2)
    ci_ok <- is.finite(loc$SE) & loc$SE > 0
    if (any(ci_ok)) {
      ci_lo <- loc$Estimate[ci_ok] - z * loc$SE[ci_ok]
      ci_hi <- loc$Estimate[ci_ok] + z * loc$SE[ci_ok]
      graphics::arrows(
        x0 = loc$X[ci_ok], y0 = ci_lo,
        x1 = loc$X[ci_ok], y1 = ci_hi,
        angle = 90, code = 3, length = 0.04,
        col = grDevices::adjustcolor(cols[ci_ok], alpha.f = 0.6), lwd = 1.2
      )
    }
  }
  graphics::points(loc$X, loc$Estimate, pch = pch, col = cols)
  graphics::legend(
    "topleft",
    legend = c("Facet level", "Step threshold"),
    pch = c(16, 17),
    col = c(pal["facet_level"], pal["step_threshold"]),
    bty = "n",
    cex = 0.9
  )
}

build_pathway_map_data <- function(x, theta_range = c(-6, 6), theta_points = 241L) {
  curve_spec <- build_step_curve_spec(x)
  theta_grid <- seq(theta_range[1], theta_range[2], length.out = theta_points)
  curve_tbl <- build_curve_tables(curve_spec, theta_grid)
  step_df <- curve_spec$step_points |>
    dplyr::mutate(PathY = curve_spec$rating_min + .data$StepIndex - 0.5)
  list(
    title = "Pathway Map (Expected Score by Theta)",
    expected = curve_tbl$expected,
    steps = step_df
  )
}

draw_pathway_map <- function(plot_data, title = NULL, palette = NULL) {
  exp_df <- plot_data$expected
  groups <- unique(exp_df$CurveGroup)
  defaults <- stats::setNames(grDevices::hcl.colors(max(3L, length(groups)), "Dark 3")[seq_along(groups)], groups)
  cols <- resolve_palette(palette = palette, defaults = defaults)
  graphics::plot(
    x = range(exp_df$Theta, finite = TRUE),
    y = range(exp_df$ExpectedScore, finite = TRUE),
    type = "n",
    xlab = "Theta / Logit",
    ylab = "Expected score",
    main = if (is.null(title)) plot_data$title else as.character(title[1])
  )
  for (i in seq_along(groups)) {
    sub <- exp_df[exp_df$CurveGroup == groups[i], , drop = FALSE]
    graphics::lines(sub$Theta, sub$ExpectedScore, col = cols[groups[i]], lwd = 2)
  }
  if (nrow(plot_data$steps) > 0) {
    step_groups <- as.character(plot_data$steps$CurveGroup)
    step_cols <- cols[step_groups]
    graphics::points(plot_data$steps$Threshold, plot_data$steps$PathY, pch = 18, col = step_cols)
  }
  graphics::legend("topleft", legend = groups, col = cols[groups], lty = 1, lwd = 2, bty = "n")
}

build_ccc_data <- function(x, theta_range = c(-6, 6), theta_points = 241L) {
  curve_spec <- build_step_curve_spec(x)
  theta_grid <- seq(theta_range[1], theta_range[2], length.out = theta_points)
  prob_df <- build_curve_tables(curve_spec, theta_grid)$probabilities
  list(
    title = "Category Characteristic Curves",
    probabilities = prob_df
  )
}

draw_ccc <- function(plot_data, title = NULL, palette = NULL) {
  prob_df <- plot_data$probabilities
  traces <- unique(paste(prob_df$CurveGroup, prob_df$Category, sep = " | Cat "))
  defaults <- stats::setNames(grDevices::hcl.colors(max(3L, length(traces)), "Dark 3")[seq_along(traces)], traces)
  cols <- resolve_palette(palette = palette, defaults = defaults)
  graphics::plot(
    x = range(prob_df$Theta, finite = TRUE),
    y = c(0, 1),
    type = "n",
    xlab = "Theta / Logit",
    ylab = "Probability",
    main = if (is.null(title)) plot_data$title else as.character(title[1])
  )
  for (i in seq_along(traces)) {
    parts <- strsplit(traces[i], " \\| Cat ", fixed = FALSE)[[1]]
    sub <- prob_df[prob_df$CurveGroup == parts[1] & prob_df$Category == parts[2], , drop = FALSE]
    graphics::lines(sub$Theta, sub$Probability, col = cols[traces[i]], lwd = 1.5)
  }
}

draw_person_plot <- function(person_tbl, bins, title = "Person measure distribution", palette = NULL) {
  pal <- resolve_palette(palette = palette, defaults = c(person_hist = "#2C7FB8"))
  graphics::hist(
    x = person_tbl$Estimate,
    breaks = bins,
    main = title,
    xlab = "Person estimate",
    ylab = "Count",
    col = pal["person_hist"],
    border = "white"
  )
}

draw_step_plot <- function(step_tbl,
                           title = "Step parameter estimates",
                           palette = NULL,
                           label_angle = 45) {
  pal <- resolve_palette(
    palette = palette,
    defaults = c(
      step_line = "#1B9E77",
      grid = "#ececec"
    )
  )
  x_idx <- step_index_from_label(step_tbl$Step)
  ord <- order(x_idx)
  step_tbl <- step_tbl[ord, , drop = FALSE]
  old_mar <- graphics::par("mar")
  on.exit(graphics::par(mar = old_mar), add = TRUE)
  mar <- old_mar
  mar[1] <- max(mar[1], 8.6)
  graphics::par(mar = mar)
  graphics::plot(
    x = seq_len(nrow(step_tbl)),
    y = step_tbl$Estimate,
    type = "b",
    pch = 16,
    xaxt = "n",
    xlab = "",
    ylab = "Estimate",
    main = title,
    col = pal["step_line"]
  )
  graphics::abline(h = pretty(step_tbl$Estimate, n = 5), col = pal["grid"], lty = 1)
  draw_rotated_x_labels(
    at = seq_len(nrow(step_tbl)),
    labels = truncate_axis_label(step_tbl$Step, width = 16L),
    srt = label_angle,
    cex = 0.84,
    line_offset = 0.085
  )
}

draw_facet_plot <- function(facet_tbl, title, palette = NULL,
                            show_ci = FALSE, ci_level = 0.95) {
  pal <- resolve_palette(palette = palette, defaults = c(facet_bar = "gray70"))
  labels <- paste0(facet_tbl$Facet, ": ", facet_tbl$Level)
  mids <- graphics::barplot(
    height = facet_tbl$Estimate,
    names.arg = labels,
    horiz = TRUE,
    las = 1,
    main = title,
    xlab = "Estimate",
    col = pal["facet_bar"],
    border = "white"
  )
  graphics::abline(v = 0, lty = 2, col = "gray40")
  if (isTRUE(show_ci) && "SE" %in% names(facet_tbl)) {
    z <- stats::qnorm(1 - (1 - ci_level) / 2)
    ci_ok <- is.finite(facet_tbl$SE) & facet_tbl$SE > 0
    if (any(ci_ok)) {
      ci_lo <- facet_tbl$Estimate[ci_ok] - z * facet_tbl$SE[ci_ok]
      ci_hi <- facet_tbl$Estimate[ci_ok] + z * facet_tbl$SE[ci_ok]
      graphics::arrows(
        x0 = ci_lo, y0 = mids[ci_ok],
        x1 = ci_hi, y1 = mids[ci_ok],
        angle = 90, code = 3, length = 0.04,
        col = "gray30", lwd = 1.2
      )
    }
  }
}

new_mfrm_plot_data <- function(name, data) {
  out <- list(name = name, data = data)
  class(out) <- c("mfrm_plot_data", class(out))
  out
}

truncate_axis_label <- function(x, width = 28L) {
  x <- as.character(x)
  width <- max(8L, as.integer(width))
  ifelse(nchar(x) > width, paste0(substr(x, 1, width - 3L), "..."), x)
}

draw_rotated_x_labels <- function(at,
                                  labels,
                                  srt = 45,
                                  cex = 0.85,
                                  line_offset = 0.08) {
  at <- as.numeric(at)
  labels <- as.character(labels)
  ok <- is.finite(at) & nzchar(labels)
  if (!any(ok)) return(invisible(NULL))

  at <- at[ok]
  labels <- labels[ok]
  graphics::axis(side = 1, at = at, labels = FALSE, tck = -0.02)

  usr <- graphics::par("usr")
  y <- usr[3] - line_offset * diff(usr[3:4])
  graphics::text(
    x = at,
    y = y,
    labels = labels,
    srt = srt,
    adj = 1,
    xpd = NA,
    cex = cex
  )
  invisible(NULL)
}

resolve_palette <- function(palette = NULL, defaults = character(0)) {
  defaults <- stats::setNames(as.character(defaults), names(defaults))
  if (length(defaults) == 0) return(defaults)
  if (is.null(palette) || length(palette) == 0) return(defaults)

  palette <- stats::setNames(as.character(palette), names(palette))
  nm <- names(palette)
  if (is.null(nm) || any(!nzchar(nm))) {
    take <- seq_len(min(length(defaults), length(palette)))
    defaults[take] <- palette[take]
    return(defaults)
  }
  hit <- intersect(names(defaults), nm)
  if (length(hit) > 0) defaults[hit] <- palette[hit]
  defaults
}

barplot_rot45 <- function(height,
                          labels,
                          col,
                          border = "white",
                          main = NULL,
                          ylab = NULL,
                          label_angle = 45,
                          label_cex = 0.84,
                          mar_bottom = 8.2,
                          label_width = 22L,
                          add_grid = FALSE,
                          ...) {
  old_mar <- graphics::par("mar")
  on.exit(graphics::par(mar = old_mar), add = TRUE)
  mar <- old_mar
  mar[1] <- max(mar[1], mar_bottom)
  graphics::par(mar = mar)

  mids <- graphics::barplot(
    height = height,
    names.arg = FALSE,
    col = col,
    border = border,
    main = main,
    ylab = ylab,
    ...
  )
  if (isTRUE(add_grid)) {
    ylim <- graphics::par("usr")[3:4]
    graphics::abline(h = pretty(ylim, n = 5), col = "#ececec", lty = 1)
  }
  draw_rotated_x_labels(
    at = mids,
    labels = truncate_axis_label(labels, width = label_width),
    srt = label_angle,
    cex = label_cex,
    line_offset = 0.085
  )
  invisible(mids)
}

stack_fair_raw_tables <- function(raw_by_facet) {
  if (is.null(raw_by_facet) || length(raw_by_facet) == 0) return(data.frame())
  out <- lapply(names(raw_by_facet), function(facet) {
    df <- raw_by_facet[[facet]]
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df <- as.data.frame(df, stringsAsFactors = FALSE)
    df$Facet <- facet
    df
  })
  out <- out[!vapply(out, is.null, logical(1))]
  if (length(out) == 0) data.frame() else dplyr::bind_rows(out)
}

resolve_unexpected_bundle <- function(x,
                                      diagnostics = NULL,
                                      abs_z_min = 2,
                                      prob_max = 0.30,
                                      top_n = 100,
                                      rule = "either") {
  if (inherits(x, "mfrm_fit")) {
    return(unexpected_response_table(
      fit = x,
      diagnostics = diagnostics,
      abs_z_min = abs_z_min,
      prob_max = prob_max,
      top_n = top_n,
      rule = rule
    ))
  }
  if (is.list(x) && all(c("table", "summary", "thresholds") %in% names(x))) {
    return(x)
  }
  stop("`x` must be an mfrm_fit object or output from unexpected_response_table().")
}

resolve_fair_bundle <- function(x,
                                diagnostics = NULL,
                                facets = NULL,
                                totalscore = TRUE,
                                umean = 0,
                                uscale = 1,
                                udecimals = 2,
                                omit_unobserved = FALSE,
                                xtreme = 0) {
  if (inherits(x, "mfrm_fit")) {
    return(fair_average_table(
      fit = x,
      diagnostics = diagnostics,
      facets = facets,
      totalscore = totalscore,
      umean = umean,
      uscale = uscale,
      udecimals = udecimals,
      omit_unobserved = omit_unobserved,
      xtreme = xtreme
    ))
  }
  if (is.list(x) && all(c("raw_by_facet", "by_facet", "stacked") %in% names(x))) {
    return(x)
  }
  stop("`x` must be an mfrm_fit object or output from fair_average_table().")
}

resolve_displacement_bundle <- function(x,
                                        diagnostics = NULL,
                                        facets = NULL,
                                        anchored_only = FALSE,
                                        abs_displacement_warn = 0.5,
                                        abs_t_warn = 2,
                                        top_n = NULL) {
  if (inherits(x, "mfrm_fit")) {
    return(displacement_table(
      fit = x,
      diagnostics = diagnostics,
      facets = facets,
      anchored_only = anchored_only,
      abs_displacement_warn = abs_displacement_warn,
      abs_t_warn = abs_t_warn,
      top_n = top_n
    ))
  }
  if (is.list(x) && all(c("table", "summary", "thresholds") %in% names(x))) {
    return(x)
  }
  stop("`x` must be an mfrm_fit object or output from displacement_table().")
}

resolve_interrater_bundle <- function(x,
                                      diagnostics = NULL,
                                      rater_facet = NULL,
                                      context_facets = NULL,
                                      exact_warn = 0.50,
                                      corr_warn = 0.30,
                                      top_n = NULL) {
  if (inherits(x, "mfrm_fit")) {
    return(interrater_agreement_table(
      fit = x,
      diagnostics = diagnostics,
      rater_facet = rater_facet,
      context_facets = context_facets,
      exact_warn = exact_warn,
      corr_warn = corr_warn,
      top_n = top_n
    ))
  }
  if (is.list(x) && all(c("summary", "pairs", "settings") %in% names(x))) {
    return(x)
  }
  stop("`x` must be an mfrm_fit object or output from interrater_agreement_table().")
}

resolve_facets_chisq_bundle <- function(x,
                                        diagnostics = NULL,
                                        fixed_p_max = 0.05,
                                        random_p_max = 0.05,
                                        top_n = NULL) {
  if (inherits(x, "mfrm_fit")) {
    return(facets_chisq_table(
      fit = x,
      diagnostics = diagnostics,
      fixed_p_max = fixed_p_max,
      random_p_max = random_p_max,
      top_n = top_n
    ))
  }
  if (is.list(x) && all(c("table", "summary", "thresholds") %in% names(x))) {
    return(x)
  }
  stop("`x` must be an mfrm_fit object or output from facets_chisq_table().")
}

#' Plot unexpected responses using base R
#'
#' @param x Output from [fit_mfrm()] or [unexpected_response_table()].
#' @param diagnostics Optional output from [diagnose_mfrm()] when `x` is `mfrm_fit`.
#' @param abs_z_min Absolute standardized-residual cutoff.
#' @param prob_max Maximum observed-category probability cutoff.
#' @param top_n Maximum rows used from the unexpected table.
#' @param rule Flagging rule (`"either"` or `"both"`).
#' @param plot_type `"scatter"` or `"severity"`.
#' @param main Optional custom plot title.
#' @param palette Optional named color overrides (`higher`, `lower`, `bar`).
#' @param label_angle X-axis label angle for `"severity"` bar plot.
#' @param draw If `TRUE`, draw with base graphics.
#'
#' @details
#' This helper visualizes flagged observations from [unexpected_response_table()]:
#' - `"scatter"`: standardized residual vs `-log10(P_obs)`
#' - `"severity"`: ranked severity index bar chart
#'
#' It accepts either a fitted object (builds the unexpected table internally)
#' or an existing unexpected-response bundle.
#'
#' @section Plot types:
#' \describe{
#'   \item{`"scatter"` (default)}{X-axis: standardized residual.
#'     Y-axis: `-log10(P_obs)` (negative log of observed-category
#'     probability).  Points colored by direction (higher/lower than
#'     expected).  Dashed lines mark `abs_z_min` and `prob_max`
#'     thresholds.  Use for global pattern detection.}
#'   \item{`"severity"`}{Ranked bar chart of the composite severity index
#'     for the `top_n` most unexpected responses.  Use for QC triage
#'     and case-level prioritization.}
#' }
#'
#' @section Interpreting output:
#' Scatter plot: farther from zero on x-axis = larger residual mismatch;
#' higher y-axis = lower observed-category probability.
#'
#' Severity plot: focuses on the most extreme observations for targeted
#' case review.
#'
#' @section Typical workflow:
#' 1. Fit model and run [diagnose_mfrm()].
#' 2. Start with `"scatter"` to assess global unexpected pattern.
#' 3. Switch to `"severity"` for case prioritization.
#'
#' @return A plotting-data object of class `mfrm_plot_data`.
#' @seealso [unexpected_response_table()], [plot_fair_average()], [plot_displacement()], [plot_qc_dashboard()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' p <- plot_unexpected(fit, abs_z_min = 1.5, prob_max = 0.4, top_n = 10, draw = FALSE)
#' if (interactive()) {
#'   plot_unexpected(
#'     fit,
#'     abs_z_min = 1.5,
#'     prob_max = 0.4,
#'     top_n = 10,
#'     plot_type = "severity",
#'     main = "Unexpected Response Severity (Customized)",
#'     palette = c(higher = "#d95f02", lower = "#1b9e77", bar = "#2b8cbe"),
#'     label_angle = 45
#'   )
#' }
#' @export
plot_unexpected <- function(x,
                            diagnostics = NULL,
                            abs_z_min = 2,
                            prob_max = 0.30,
                            top_n = 100,
                            rule = c("either", "both"),
                            plot_type = c("scatter", "severity"),
                            main = NULL,
                            palette = NULL,
                            label_angle = 45,
                            draw = TRUE) {
  rule <- match.arg(tolower(rule), c("either", "both"))
  plot_type <- match.arg(tolower(plot_type), c("scatter", "severity"))
  top_n <- max(1L, as.integer(top_n))
  pal <- resolve_palette(
    palette = palette,
    defaults = c(
      higher = "#d95f02",
      lower = "#1b9e77",
      bar = "#1f78b4"
    )
  )

  bundle <- resolve_unexpected_bundle(
    x = x,
    diagnostics = diagnostics,
    abs_z_min = abs_z_min,
    prob_max = prob_max,
    top_n = top_n,
    rule = rule
  )
  tbl <- as.data.frame(bundle$table, stringsAsFactors = FALSE)
  if (nrow(tbl) == 0) {
    stop("No unexpected responses were flagged under the current thresholds.")
  }
  tbl <- tbl[seq_len(min(nrow(tbl), top_n)), , drop = FALSE]

  if (isTRUE(draw)) {
    if (plot_type == "scatter") {
      x_vals <- suppressWarnings(as.numeric(tbl$StdResidual))
      y_vals <- -log10(pmax(suppressWarnings(as.numeric(tbl$ObsProb)), .Machine$double.xmin))
      dirs <- as.character(tbl$Direction)
      cols <- ifelse(dirs == "Higher than expected", pal["higher"], pal["lower"])
      cols[!is.finite(x_vals) | !is.finite(y_vals)] <- "gray60"
      graphics::plot(
        x = x_vals,
        y = y_vals,
        xlab = "Standardized residual",
        ylab = expression(-log[10](P[obs])),
        main = if (is.null(main)) "Unexpected responses" else as.character(main[1]),
        pch = 16,
        col = cols
      )
      graphics::abline(v = c(-abs_z_min, abs_z_min), lty = 2, col = "gray45")
      graphics::abline(h = -log10(prob_max), lty = 2, col = "gray45")
      graphics::legend(
        "topleft",
        legend = c("Higher than expected", "Lower than expected"),
        col = c(pal["higher"], pal["lower"]),
        pch = 16,
        bty = "n",
        cex = 0.85
      )
    } else {
      sev <- suppressWarnings(as.numeric(tbl$Severity))
      ord <- order(sev, decreasing = TRUE, na.last = NA)
      use <- ord[seq_len(min(length(ord), top_n))]
      sev <- sev[use]
      labels <- if ("Row" %in% names(tbl)) {
        paste0("Row ", tbl$Row[use])
      } else {
        paste0("Case ", seq_along(use))
      }
      barplot_rot45(
        height = sev,
        labels = labels,
        col = pal["bar"],
        main = if (is.null(main)) "Unexpected response severity" else as.character(main[1]),
        ylab = "Severity index",
        label_angle = label_angle,
        mar_bottom = 8.2
      )
    }
  }

  out <- new_mfrm_plot_data(
    "unexpected",
    list(
      plot = plot_type,
      table = tbl,
      summary = bundle$summary,
      thresholds = bundle$thresholds
    )
  )
  invisible(out)
}

#' Plot fair-average diagnostics using base R
#'
#' @param x Output from [fit_mfrm()] or [fair_average_table()].
#' @param diagnostics Optional output from [diagnose_mfrm()] when `x` is `mfrm_fit`.
#' @param facet Optional facet name for level-wise lollipop plots.
#' @param metric Fair-average metric (`"FairM"` or `"FairZ"`).
#' @param plot_type `"difference"` or `"scatter"`.
#' @param top_n Maximum levels shown for `"difference"` plot.
#' @param draw If `TRUE`, draw with base graphics.
#' @param ... Additional arguments passed to [fair_average_table()] when `x` is `mfrm_fit`.
#'
#' @details
#' Fair-average plots compare observed scoring tendency against model-based
#' fair metrics.
#'
#' @section Plot types:
#' \describe{
#'   \item{`"difference"` (default)}{Lollipop chart showing the gap between
#'     observed and fair-average score for each element.  X-axis:
#'     Observed - Fair metric.  Y-axis: element labels.  Points colored
#'     teal (lenient, gap >= 0) or orange (severe, gap < 0).  Ordered by
#'     absolute gap.}
#'   \item{`"scatter"`}{Scatter plot of fair metric (x) vs observed average
#'     (y) with an identity line.  Points colored by facet.  Useful for
#'     checking overall alignment between observed and model-adjusted
#'     scores.}
#' }
#'
#' @section Interpreting output:
#' Difference plot: ranked element-level gaps (`Observed - Fair`), useful
#' for triage of potentially lenient/severe levels.
#'
#' Scatter plot: global agreement pattern relative to the identity line.
#'
#' Larger absolute gaps suggest stronger divergence between observed and
#' model-adjusted scoring.
#'
#' @section Typical workflow:
#' 1. Start with `plot_type = "difference"` to find largest discrepancies.
#' 2. Use `plot_type = "scatter"` to check overall alignment pattern.
#' 3. Follow up with facet-level diagnostics for flagged levels.
#'
#' @return A plotting-data object of class `mfrm_plot_data`.
#' @seealso [fair_average_table()], [plot_unexpected()], [plot_displacement()], [plot_qc_dashboard()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' p <- plot_fair_average(fit, metric = "FairM", draw = FALSE)
#' if (interactive()) {
#'   plot_fair_average(fit, metric = "FairM", plot_type = "difference")
#' }
#' @export
plot_fair_average <- function(x,
                              diagnostics = NULL,
                              facet = NULL,
                              metric = c("FairM", "FairZ"),
                              plot_type = c("difference", "scatter"),
                              top_n = 40,
                              draw = TRUE,
                              ...) {
  metric <- match.arg(metric, c("FairM", "FairZ"))
  plot_type <- match.arg(tolower(plot_type), c("difference", "scatter"))
  top_n <- max(1L, as.integer(top_n))

  bundle <- if (inherits(x, "mfrm_fit")) {
    fair_average_table(x, diagnostics = diagnostics, ...)
  } else {
    resolve_fair_bundle(x)
  }

  fair_df <- stack_fair_raw_tables(bundle$raw_by_facet)
  if (nrow(fair_df) == 0) stop("No fair-average data available.")
  needed <- c("Facet", "Level", "ObservedAverage", "FairM", "FairZ")
  if (!all(needed %in% names(fair_df))) {
    stop("Fair-average table does not include required columns.")
  }
  fair_df <- fair_df[is.finite(fair_df$ObservedAverage) & is.finite(fair_df[[metric]]), , drop = FALSE]
  if (nrow(fair_df) == 0) stop("No finite fair-average rows available.")

  if (!is.null(facet)) {
    fair_df <- fair_df[as.character(fair_df$Facet) == as.character(facet[1]), , drop = FALSE]
    if (nrow(fair_df) == 0) stop("Requested `facet` was not found in fair-average output.")
  }
  fair_df$Gap <- fair_df$ObservedAverage - fair_df[[metric]]

  if (isTRUE(draw)) {
    if (plot_type == "difference") {
      ord <- order(abs(fair_df$Gap), decreasing = TRUE, na.last = NA)
      use <- ord[seq_len(min(length(ord), top_n))]
      sub <- fair_df[use, , drop = FALSE]
      y <- seq_len(nrow(sub))
      lbl <- paste0(sub$Facet, ":", sub$Level)
      lbl <- truncate_axis_label(lbl, width = 26L)
      graphics::plot(
        x = sub$Gap,
        y = y,
        type = "n",
        xlab = paste0("Observed - ", metric),
        ylab = "",
        yaxt = "n",
        main = paste0("Fair-average gaps (", metric, ")")
      )
      graphics::segments(x0 = 0, y0 = y, x1 = sub$Gap, y1 = y, col = "gray55")
      cols <- ifelse(sub$Gap >= 0, "#1b9e77", "#d95f02")
      graphics::points(sub$Gap, y, pch = 16, col = cols)
      graphics::axis(side = 2, at = y, labels = lbl, las = 2, cex.axis = 0.75)
      graphics::abline(v = 0, lty = 2, col = "gray40")
    } else {
      fac <- as.character(fair_df$Facet)
      fac_levels <- unique(fac)
      col_idx <- match(fac, fac_levels)
      cols <- grDevices::hcl.colors(length(fac_levels), "Dark 3")[col_idx]
      graphics::plot(
        x = fair_df[[metric]],
        y = fair_df$ObservedAverage,
        xlab = metric,
        ylab = "Observed average",
        main = paste0("Observed vs ", metric),
        pch = 16,
        col = cols
      )
      lims <- range(c(fair_df[[metric]], fair_df$ObservedAverage), finite = TRUE)
      graphics::abline(a = 0, b = 1, lty = 2, col = "gray45")
      graphics::legend("topleft", legend = fac_levels, col = grDevices::hcl.colors(length(fac_levels), "Dark 3"), pch = 16, bty = "n", cex = 0.85)
      graphics::segments(x0 = lims[1], y0 = lims[1], x1 = lims[2], y1 = lims[2], col = "gray70", lty = 3)
    }
  }

  out <- new_mfrm_plot_data(
    "fair_average",
    list(
      plot = plot_type,
      metric = metric,
      data = fair_df,
      settings = bundle$settings
    )
  )
  invisible(out)
}

#' Plot displacement diagnostics using base R
#'
#' @param x Output from [fit_mfrm()] or [displacement_table()].
#' @param diagnostics Optional output from [diagnose_mfrm()] when `x` is `mfrm_fit`.
#' @param anchored_only Keep only anchored/group-anchored levels.
#' @param facets Optional subset of facets.
#' @param plot_type `"lollipop"` or `"hist"`.
#' @param top_n Maximum levels shown in `"lollipop"` mode.
#' @param draw If `TRUE`, draw with base graphics.
#' @param ... Additional arguments passed to [displacement_table()] when `x` is `mfrm_fit`.
#'
#' @details
#' Displacement plots focus on anchor stability and estimate movement.
#'
#' @section Plot types:
#' \describe{
#'   \item{`"lollipop"` (default)}{Dot-and-line chart of displacement values.
#'     X-axis: displacement (logits).  Y-axis: element labels.  Points
#'     colored red when flagged.  Dashed lines at +/- threshold.  Ordered
#'     by absolute displacement.}
#'   \item{`"hist"`}{Histogram of displacement values with Freedman-Diaconis
#'     breaks.  Dashed reference lines at +/- threshold.  Use for
#'     inspecting the overall distribution shape.}
#' }
#'
#' @section Interpreting output:
#' Lollipop: top absolute displacement levels; flagged points indicate
#' larger movement from anchor expectations.
#'
#' Histogram: overall displacement distribution and threshold lines.
#'
#' Use `anchored_only = TRUE` when your main question is anchor robustness.
#'
#' @section Typical workflow:
#' 1. Run with `plot_type = "lollipop"` and `anchored_only = TRUE`.
#' 2. Inspect distribution with `plot_type = "hist"`.
#' 3. Drill into flagged rows via [displacement_table()].
#'
#' @return A plotting-data object of class `mfrm_plot_data`.
#' @seealso [displacement_table()], [plot_unexpected()], [plot_fair_average()], [plot_qc_dashboard()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' p <- plot_displacement(fit, anchored_only = FALSE, draw = FALSE)
#' if (interactive()) {
#'   plot_displacement(fit, anchored_only = FALSE, plot_type = "lollipop")
#' }
#' @export
plot_displacement <- function(x,
                              diagnostics = NULL,
                              anchored_only = FALSE,
                              facets = NULL,
                              plot_type = c("lollipop", "hist"),
                              top_n = 40,
                              draw = TRUE,
                              ...) {
  plot_type <- match.arg(tolower(plot_type), c("lollipop", "hist"))
  top_n <- max(1L, as.integer(top_n))

  bundle <- if (inherits(x, "mfrm_fit")) {
    displacement_table(
      fit = x,
      diagnostics = diagnostics,
      facets = facets,
      anchored_only = anchored_only,
      ...
    )
  } else {
    resolve_displacement_bundle(x)
  }

  tbl <- as.data.frame(bundle$table, stringsAsFactors = FALSE)
  if (nrow(tbl) == 0) stop("No displacement rows available.")
  tbl <- tbl[is.finite(tbl$Displacement), , drop = FALSE]
  if (nrow(tbl) == 0) stop("No finite displacement values available.")

  if (!is.null(facets)) {
    tbl <- tbl[as.character(tbl$Facet) %in% as.character(facets), , drop = FALSE]
  }
  if (isTRUE(anchored_only) && "AnchorType" %in% names(tbl)) {
    tbl <- tbl[tbl$AnchorType %in% c("Anchor", "Group"), , drop = FALSE]
  }
  if (nrow(tbl) == 0) stop("No rows left after filtering.")

  if (isTRUE(draw)) {
    if (plot_type == "lollipop") {
      ord <- order(abs(tbl$Displacement), decreasing = TRUE, na.last = NA)
      use <- ord[seq_len(min(length(ord), top_n))]
      sub <- tbl[use, , drop = FALSE]
      y <- seq_len(nrow(sub))
      lbl <- truncate_axis_label(paste0(sub$Facet, ":", sub$Level), width = 26L)
      cols <- ifelse(isTRUE(sub$Flag), "#d73027", "#1b9e77")
      graphics::plot(
        x = sub$Displacement,
        y = y,
        type = "n",
        xlab = "Displacement (logit)",
        ylab = "",
        yaxt = "n",
        main = "Displacement diagnostics"
      )
      graphics::segments(0, y, sub$Displacement, y, col = "gray60")
      graphics::points(sub$Displacement, y, pch = 16, col = cols)
      graphics::axis(side = 2, at = y, labels = lbl, las = 2, cex.axis = 0.75)
      d_thr <- as.numeric(bundle$thresholds$abs_displacement_warn %||% 0.5)
      graphics::abline(v = c(-d_thr, 0, d_thr), lty = c(2, 1, 2), col = c("gray45", "gray30", "gray45"))
    } else {
      vals <- suppressWarnings(as.numeric(tbl$Displacement))
      graphics::hist(
        x = vals,
        breaks = "FD",
        col = "#9ecae1",
        border = "white",
        main = "Displacement distribution",
        xlab = "Displacement (logit)"
      )
      d_thr <- as.numeric(bundle$thresholds$abs_displacement_warn %||% 0.5)
      graphics::abline(v = c(-d_thr, d_thr), lty = 2, col = "gray45")
    }
  }

  out <- new_mfrm_plot_data(
    "displacement",
    list(
      plot = plot_type,
      table = tbl,
      summary = bundle$summary,
      thresholds = bundle$thresholds
    )
  )
  invisible(out)
}

#' Plot inter-rater agreement diagnostics using base R
#'
#' @param x Output from [fit_mfrm()] or [interrater_agreement_table()].
#' @param diagnostics Optional output from [diagnose_mfrm()] when `x` is `mfrm_fit`.
#' @param rater_facet Name of the rater facet when `x` is `mfrm_fit`.
#' @param context_facets Optional context facets when `x` is `mfrm_fit`.
#' @param exact_warn Warning threshold for exact agreement.
#' @param corr_warn Warning threshold for pairwise correlation.
#' @param plot_type `"exact"`, `"corr"`, or `"difference"`.
#' @param top_n Maximum pairs displayed for bar-style plots.
#' @param main Optional custom plot title.
#' @param palette Optional named color overrides (`ok`, `flag`, `expected`).
#' @param label_angle X-axis label angle for bar-style plots.
#' @param draw If `TRUE`, draw with base graphics.
#'
#' @details
#' Inter-rater agreement plots summarize pairwise consistency for a chosen
#' rater facet under matched contexts.
#'
#' @section Plot types:
#' \describe{
#'   \item{`"exact"` (default)}{Bar chart of exact agreement proportion by
#'     rater pair.  Expected agreement overlaid as connected circles.
#'     Horizontal reference line at `exact_warn`.  Bars colored red when
#'     flagged.}
#'   \item{`"corr"`}{Bar chart of pairwise Pearson correlation by rater
#'     pair.  Reference line at `corr_warn`.  Ordered by correlation
#'     (lowest first).}
#'   \item{`"difference"`}{Scatter plot.  X-axis: mean signed score
#'     difference (Rater1 - Rater2).  Y-axis: mean absolute difference.
#'     Points colored by flag status.  Vertical reference at 0.}
#' }
#'
#' @section Interpreting output:
#' Pairs below `exact_warn` and/or `corr_warn` should be prioritized for
#' calibration review.
#'
#' @section Typical workflow:
#' 1. Select rater facet and run `"exact"` view.
#' 2. Confirm with `"corr"` view.
#' 3. Use `"difference"` to inspect directional disagreement.
#'
#' @return A plotting-data object of class `mfrm_plot_data`.
#' @seealso [interrater_agreement_table()], [plot_facets_chisq()], [plot_qc_dashboard()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' p <- plot_interrater_agreement(fit, rater_facet = "Rater", draw = FALSE)
#' if (interactive()) {
#'   plot_interrater_agreement(
#'     fit,
#'     rater_facet = "Rater",
#'     draw = TRUE,
#'     plot_type = "exact",
#'     main = "Inter-rater Agreement (Customized)",
#'     palette = c(ok = "#2b8cbe", flag = "#cb181d"),
#'     label_angle = 45
#'   )
#' }
#' @export
plot_interrater_agreement <- function(x,
                                      diagnostics = NULL,
                                      rater_facet = NULL,
                                      context_facets = NULL,
                                      exact_warn = 0.50,
                                      corr_warn = 0.30,
                                      plot_type = c("exact", "corr", "difference"),
                                      top_n = 20,
                                      main = NULL,
                                      palette = NULL,
                                      label_angle = 45,
                                      draw = TRUE) {
  plot_type <- match.arg(tolower(plot_type), c("exact", "corr", "difference"))
  top_n <- max(1L, as.integer(top_n))
  pal <- resolve_palette(
    palette = palette,
    defaults = c(
      ok = "#2b8cbe",
      flag = "#cb181d",
      expected = "#08519c"
    )
  )

  bundle <- resolve_interrater_bundle(
    x = x,
    diagnostics = diagnostics,
    rater_facet = rater_facet,
    context_facets = context_facets,
    exact_warn = exact_warn,
    corr_warn = corr_warn,
    top_n = NULL
  )

  tbl <- as.data.frame(bundle$pairs, stringsAsFactors = FALSE)
  if (nrow(tbl) == 0) stop("No inter-rater pair rows are available.")
  if (!all(c("Rater1", "Rater2", "Exact", "Corr", "MeanDiff", "MAD") %in% names(tbl))) {
    stop("Inter-rater table does not include required columns.")
  }

  ord_exact <- order(tbl$Exact, na.last = NA)
  use <- ord_exact[seq_len(min(length(ord_exact), top_n))]
  sub <- tbl[use, , drop = FALSE]
  labels <- truncate_axis_label(paste0(sub$Rater1, " | ", sub$Rater2), width = 28L)
  cols <- if ("Flag" %in% names(sub)) ifelse(sub$Flag, pal["flag"], pal["ok"]) else pal["ok"]

  if (isTRUE(draw)) {
    if (plot_type == "exact") {
      bp <- barplot_rot45(
        height = suppressWarnings(as.numeric(sub$Exact)),
        labels = labels,
        col = cols,
        main = if (is.null(main)) "Inter-rater exact agreement" else as.character(main[1]),
        ylab = "Exact agreement",
        label_angle = label_angle,
        mar_bottom = 8.2
      )
      exp_vals <- suppressWarnings(as.numeric(sub$ExpectedExact))
      if (any(is.finite(exp_vals))) {
        graphics::points(bp, exp_vals, pch = 21, bg = "white", col = pal["expected"])
        graphics::lines(bp, exp_vals, col = pal["expected"], lwd = 1.3)
      }
      graphics::abline(h = exact_warn, lty = 2, col = "gray45")
    } else if (plot_type == "corr") {
      corr_ord <- order(tbl$Corr, na.last = NA)
      use_corr <- corr_ord[seq_len(min(length(corr_ord), top_n))]
      sub_corr <- tbl[use_corr, , drop = FALSE]
      lbl_corr <- truncate_axis_label(paste0(sub_corr$Rater1, " | ", sub_corr$Rater2), width = 28L)
      col_corr <- if ("Flag" %in% names(sub_corr)) ifelse(sub_corr$Flag, pal["flag"], pal["ok"]) else pal["ok"]
      barplot_rot45(
        height = suppressWarnings(as.numeric(sub_corr$Corr)),
        labels = lbl_corr,
        col = col_corr,
        main = if (is.null(main)) "Inter-rater correlation" else as.character(main[1]),
        ylab = "Correlation",
        label_angle = label_angle,
        mar_bottom = 8.2
      )
      graphics::abline(h = corr_warn, lty = 2, col = "gray45")
    } else {
      graphics::plot(
        x = suppressWarnings(as.numeric(tbl$MeanDiff)),
        y = suppressWarnings(as.numeric(tbl$MAD)),
        pch = 16,
        col = if ("Flag" %in% names(tbl)) ifelse(tbl$Flag, pal["flag"], pal["ok"]) else pal["ok"],
        xlab = "Mean score difference (Rater1 - Rater2)",
        ylab = "Mean absolute difference",
        main = if (is.null(main)) "Inter-rater difference profile" else as.character(main[1])
      )
      graphics::abline(v = 0, lty = 2, col = "gray45")
    }
  }

  out <- new_mfrm_plot_data(
    "interrater",
    list(
      plot = plot_type,
      pairs = tbl,
      summary = bundle$summary,
      settings = bundle$settings
    )
  )
  invisible(out)
}

#' Plot FACETS-style facet chi-square diagnostics using base R
#'
#' @param x Output from [fit_mfrm()] or [facets_chisq_table()].
#' @param diagnostics Optional output from [diagnose_mfrm()] when `x` is `mfrm_fit`.
#' @param fixed_p_max Warning cutoff for fixed-effect chi-square p-values.
#' @param random_p_max Warning cutoff for random-effect chi-square p-values.
#' @param plot_type `"fixed"`, `"random"`, or `"variance"`.
#' @param main Optional custom plot title.
#' @param palette Optional named color overrides (`fixed_ok`, `fixed_flag`,
#' `random_ok`, `random_flag`, `variance`).
#' @param label_angle X-axis label angle for bar-style plots.
#' @param draw If `TRUE`, draw with base graphics.
#'
#' @details
#' Facet chi-square plots provide facet-level global fit diagnostics using
#' fixed/random chi-square summaries and random-variance magnitudes.
#'
#' @section Plot types:
#' \describe{
#'   \item{`"fixed"` (default)}{Bar chart of fixed-effect chi-square by
#'     facet.  Bars colored red when the null hypothesis (all elements
#'     equal) is rejected at `fixed_p_max`.}
#'   \item{`"random"`}{Bar chart of random-effect chi-square by facet.
#'     Bars colored red when rejected at `random_p_max`.}
#'   \item{`"variance"`}{Bar chart of estimated random variance by facet.
#'     Reference line at 0.  Larger variance indicates greater
#'     heterogeneity among elements.}
#' }
#'
#' @section Interpreting output:
#' Colored flags reflect configured p-value thresholds (`fixed_p_max`,
#' `random_p_max`) when available.
#'
#' @section Typical workflow:
#' 1. Review `"fixed"` and `"random"` panels for flagged facets.
#' 2. Check `"variance"` to contextualize heterogeneity.
#' 3. Cross-check with inter-rater and element-level fit diagnostics.
#'
#' @return A plotting-data object of class `mfrm_plot_data`.
#' @seealso [facets_chisq_table()], [plot_interrater_agreement()], [plot_qc_dashboard()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' p <- plot_facets_chisq(fit, draw = FALSE)
#' if (interactive()) {
#'   plot_facets_chisq(
#'     fit,
#'     draw = TRUE,
#'     plot_type = "fixed",
#'     main = "Facet Chi-square (Customized)",
#'     palette = c(fixed_ok = "#2b8cbe", fixed_flag = "#cb181d"),
#'     label_angle = 45
#'   )
#' }
#' @export
plot_facets_chisq <- function(x,
                              diagnostics = NULL,
                              fixed_p_max = 0.05,
                              random_p_max = 0.05,
                              plot_type = c("fixed", "random", "variance"),
                              main = NULL,
                              palette = NULL,
                              label_angle = 45,
                              draw = TRUE) {
  plot_type <- match.arg(tolower(plot_type), c("fixed", "random", "variance"))
  pal <- resolve_palette(
    palette = palette,
    defaults = c(
      fixed_ok = "#2b8cbe",
      fixed_flag = "#cb181d",
      random_ok = "#31a354",
      random_flag = "#cb181d",
      variance = "#9ecae1"
    )
  )
  bundle <- resolve_facets_chisq_bundle(
    x = x,
    diagnostics = diagnostics,
    fixed_p_max = fixed_p_max,
    random_p_max = random_p_max
  )

  tbl <- as.data.frame(bundle$table, stringsAsFactors = FALSE)
  if (nrow(tbl) == 0) stop("No facet chi-square rows are available.")
  if (!all(c("Facet", "FixedChiSq", "RandomChiSq", "RandomVar") %in% names(tbl))) {
    stop("Facet chi-square table does not include required columns.")
  }

  if (isTRUE(draw)) {
    facet_labels <- truncate_axis_label(as.character(tbl$Facet), width = 20L)
    if (plot_type == "fixed") {
      ord <- order(tbl$FixedChiSq, decreasing = TRUE, na.last = NA)
      sub <- tbl[ord, , drop = FALSE]
      col_fixed <- if ("FixedFlag" %in% names(sub)) ifelse(sub$FixedFlag, pal["fixed_flag"], pal["fixed_ok"]) else pal["fixed_ok"]
      barplot_rot45(
        height = suppressWarnings(as.numeric(sub$FixedChiSq)),
        labels = truncate_axis_label(as.character(sub$Facet), width = 20L),
        col = col_fixed,
        main = if (is.null(main)) "Facet fixed-effect chi-square" else as.character(main[1]),
        ylab = expression(chi^2),
        label_angle = label_angle,
        mar_bottom = 8.2
      )
    } else if (plot_type == "random") {
      ord <- order(tbl$RandomChiSq, decreasing = TRUE, na.last = NA)
      sub <- tbl[ord, , drop = FALSE]
      col_random <- if ("RandomFlag" %in% names(sub)) ifelse(sub$RandomFlag, pal["random_flag"], pal["random_ok"]) else pal["random_ok"]
      barplot_rot45(
        height = suppressWarnings(as.numeric(sub$RandomChiSq)),
        labels = truncate_axis_label(as.character(sub$Facet), width = 20L),
        col = col_random,
        main = if (is.null(main)) "Facet random-effect chi-square" else as.character(main[1]),
        ylab = expression(chi^2),
        label_angle = label_angle,
        mar_bottom = 8.2
      )
    } else {
      vals <- suppressWarnings(as.numeric(tbl$RandomVar))
      barplot_rot45(
        height = vals,
        labels = facet_labels,
        col = pal["variance"],
        main = if (is.null(main)) "Facet random variance" else as.character(main[1]),
        ylab = "Variance",
        label_angle = label_angle,
        mar_bottom = 8.2
      )
      graphics::abline(h = 0, lty = 2, col = "gray45")
    }
  }

  out <- new_mfrm_plot_data(
    "facets_chisq",
    list(
      plot = plot_type,
      table = tbl,
      summary = bundle$summary,
      thresholds = bundle$thresholds
    )
  )
  invisible(out)
}

#' Plot a base-R QC dashboard
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()].
#' @param threshold_profile Threshold profile name (`strict`, `standard`, `lenient`).
#' @param thresholds Optional named threshold overrides.
#' @param abs_z_min Absolute standardized-residual cutoff for unexpected panel.
#' @param prob_max Maximum observed-category probability cutoff for unexpected panel.
#' @param rater_facet Optional rater facet used in inter-rater panel.
#' @param interrater_exact_warn Warning threshold for inter-rater exact agreement.
#' @param interrater_corr_warn Warning threshold for inter-rater correlation.
#' @param fixed_p_max Warning cutoff for fixed-effect facet chi-square p-values.
#' @param random_p_max Warning cutoff for random-effect facet chi-square p-values.
#' @param top_n Maximum elements displayed in displacement panel.
#' @param draw If `TRUE`, draw with base graphics.
#'
#' @details
#' The dashboard draws nine QC panels:
#' - observed vs expected category counts
#' - infit vs outfit scatter
#' - |ZSTD| histogram
#' - unexpected-response scatter
#' - fair-average gap by facet
#' - displacement lollipop (largest absolute values)
#' - inter-rater exact agreement by pair
#' - facet fixed-effect chi-square by facet
#' - separation reliability by facet
#'
#' `threshold_profile` controls warning overlays used in applicable panels.
#' Use `thresholds` to override any profile value with named entries.
#'
#' @section Plot types:
#' This function draws a fixed 3x3 panel grid (no `plot_type` argument).
#' For individual panel control, use the dedicated helpers:
#' [plot_unexpected()], [plot_fair_average()], [plot_displacement()],
#' [plot_interrater_agreement()], [plot_facets_chisq()].
#'
#' @section Interpreting output:
#' Recommended panel order for fast review:
#' 1. category counts and infit/outfit scatter (global health)
#' 2. unexpected responses and displacement (element-level outliers)
#' 3. inter-rater and facet chi-square panels (facet-level comparability)
#' 4. reliability panel (precision/separation context)
#'
#' Treat this dashboard as a screening layer; follow up with dedicated helpers
#' (`plot_unexpected()`, `plot_displacement()`, `plot_interrater_agreement()`,
#' `plot_facets_chisq()`) for detailed diagnosis.
#'
#' @section Typical workflow:
#' 1. Fit and diagnose model.
#' 2. Run `plot_qc_dashboard()` for one-page triage.
#' 3. Drill into flagged panels using dedicated functions.
#'
#' @return A plotting-data object of class `mfrm_plot_data`.
#' @seealso [plot_unexpected()], [plot_fair_average()], [plot_displacement()], [plot_interrater_agreement()], [plot_facets_chisq()], [build_visual_summaries()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' qc <- plot_qc_dashboard(fit, draw = FALSE)
#' if (interactive()) {
#'   plot_qc_dashboard(fit, rater_facet = "Rater")
#' }
#' @export
plot_qc_dashboard <- function(fit,
                              diagnostics = NULL,
                              threshold_profile = "standard",
                              thresholds = NULL,
                              abs_z_min = 2,
                              prob_max = 0.30,
                              rater_facet = NULL,
                              interrater_exact_warn = 0.50,
                              interrater_corr_warn = 0.30,
                              fixed_p_max = 0.05,
                              random_p_max = 0.05,
                              top_n = 20,
                              draw = TRUE) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().")
  }
  top_n <- max(5L, as.integer(top_n))
  if (is.null(diagnostics)) {
    diagnostics <- diagnose_mfrm(fit, residual_pca = "none")
  }
  if (is.null(diagnostics$obs) || nrow(diagnostics$obs) == 0) {
    stop("`diagnostics$obs` is empty. Run diagnose_mfrm() first.")
  }

  resolved <- resolve_warning_thresholds(thresholds = thresholds, threshold_profile = threshold_profile)
  cat_tbl <- calc_category_stats(diagnostics$obs, res = fit, whexact = FALSE)
  fit_tbl <- as.data.frame(diagnostics$fit, stringsAsFactors = FALSE)
  zstd <- if (nrow(fit_tbl) > 0) {
    pmax(abs(suppressWarnings(as.numeric(fit_tbl$InfitZSTD))), abs(suppressWarnings(as.numeric(fit_tbl$OutfitZSTD))), na.rm = TRUE)
  } else {
    numeric(0)
  }
  zstd <- zstd[is.finite(zstd)]

  unexpected <- unexpected_response_table(
    fit = fit,
    diagnostics = diagnostics,
    abs_z_min = abs_z_min,
    prob_max = prob_max,
    top_n = max(top_n, 20),
    rule = "either"
  )
  fair <- fair_average_table(fit = fit, diagnostics = diagnostics)
  fair_df <- stack_fair_raw_tables(fair$raw_by_facet)
  fair_gap <- if (nrow(fair_df) > 0 && all(c("ObservedAverage", "FairM") %in% names(fair_df))) {
    fair_df$ObservedAverage - fair_df$FairM
  } else {
    numeric(0)
  }
  disp <- displacement_table(
    fit = fit,
    diagnostics = diagnostics,
    anchored_only = FALSE
  )
  disp_tbl <- as.data.frame(disp$table, stringsAsFactors = FALSE)
  interrater <- interrater_agreement_table(
    fit = fit,
    diagnostics = diagnostics,
    rater_facet = rater_facet,
    exact_warn = interrater_exact_warn,
    corr_warn = interrater_corr_warn,
    top_n = max(top_n, 20)
  )
  inter_tbl <- as.data.frame(interrater$pairs, stringsAsFactors = FALSE)
  fchi <- facets_chisq_table(
    fit = fit,
    diagnostics = diagnostics,
    fixed_p_max = fixed_p_max,
    random_p_max = random_p_max
  )
  fchi_tbl <- as.data.frame(fchi$table, stringsAsFactors = FALSE)
  rel_tbl <- as.data.frame(diagnostics$reliability, stringsAsFactors = FALSE)

  if (isTRUE(draw)) {
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    graphics::par(mfrow = c(3, 3), mar = c(4, 4, 3, 1))

    # 1) Category counts
    if (nrow(cat_tbl) > 0) {
      cat_lbl <- as.character(cat_tbl$Category)
      obs_ct <- suppressWarnings(as.numeric(cat_tbl$Count))
      exp_ct <- suppressWarnings(as.numeric(cat_tbl$ExpectedCount))
      bp <- barplot_rot45(
        height = obs_ct,
        labels = cat_lbl,
        col = "#9ecae1",
        main = "QC: Category counts",
        ylab = "Count",
        label_angle = 45,
        mar_bottom = 6.4,
        label_cex = 0.72,
        label_width = 14L
      )
      if (all(is.finite(exp_ct))) {
        graphics::points(bp, exp_ct, pch = 21, bg = "white", col = "#08519c")
        graphics::lines(bp, exp_ct, col = "#08519c", lwd = 1.5)
      }
    } else {
      graphics::plot.new()
      graphics::title(main = "QC: Category counts")
      graphics::text(0.5, 0.5, "No data")
    }

    # 2) Infit/Outfit scatter
    if (nrow(fit_tbl) > 0) {
      infit <- suppressWarnings(as.numeric(fit_tbl$Infit))
      outfit <- suppressWarnings(as.numeric(fit_tbl$Outfit))
      ok <- is.finite(infit) & is.finite(outfit)
      graphics::plot(
        x = infit[ok],
        y = outfit[ok],
        pch = 16,
        col = "#2c7fb8",
        xlab = "Infit MnSq",
        ylab = "Outfit MnSq",
        main = "QC: Infit vs Outfit"
      )
      graphics::abline(v = c(0.5, 1, 1.5), h = c(0.5, 1, 1.5), lty = c(2, 1, 2), col = "gray50")
    } else {
      graphics::plot.new()
      graphics::title(main = "QC: Infit vs Outfit")
      graphics::text(0.5, 0.5, "No data")
    }

    # 3) |ZSTD| histogram
    if (length(zstd) > 0) {
      graphics::hist(
        x = zstd,
        breaks = "FD",
        col = "#c7e9c0",
        border = "white",
        main = "QC: |ZSTD| distribution",
        xlab = "|ZSTD|"
      )
      graphics::abline(v = c(2, 3), lty = 2, col = "gray45")
    } else {
      graphics::plot.new()
      graphics::title(main = "QC: |ZSTD| distribution")
      graphics::text(0.5, 0.5, "No data")
    }

    # 4) Unexpected response scatter
    if (nrow(unexpected$table) > 0) {
      ut <- unexpected$table
      x_u <- suppressWarnings(as.numeric(ut$StdResidual))
      y_u <- -log10(pmax(suppressWarnings(as.numeric(ut$ObsProb)), .Machine$double.xmin))
      graphics::plot(
        x = x_u,
        y = y_u,
        pch = 16,
        col = "#756bb1",
        xlab = "Std residual",
        ylab = expression(-log[10](P[obs])),
        main = "QC: Unexpected responses"
      )
      graphics::abline(v = c(-abs_z_min, abs_z_min), h = -log10(prob_max), lty = 2, col = "gray45")
    } else {
      graphics::plot.new()
      graphics::title(main = "QC: Unexpected responses")
      graphics::text(0.5, 0.5, "No flagged rows")
    }

    # 5) Fair-average gap
    if (length(fair_gap) > 0) {
      fac <- as.character(fair_df$Facet)
      split_gap <- split(fair_gap, fac)
      old_mar <- graphics::par("mar")
      mar <- old_mar
      mar[1] <- max(mar[1], 6.4)
      graphics::par(mar = mar)
      graphics::boxplot(
        split_gap,
        xaxt = "n",
        col = "#fdd0a2",
        main = "QC: Observed - Fair(M)",
        ylab = "Gap"
      )
      draw_rotated_x_labels(
        at = seq_along(split_gap),
        labels = truncate_axis_label(names(split_gap), width = 14L),
        srt = 45,
        cex = 0.72,
        line_offset = 0.085
      )
      graphics::mtext("Facet", side = 1, line = 4.8, cex = 0.82)
      graphics::par(mar = old_mar)
      graphics::abline(h = 0, lty = 2, col = "gray45")
    } else {
      graphics::plot.new()
      graphics::title(main = "QC: Observed - Fair(M)")
      graphics::text(0.5, 0.5, "No data")
    }

    # 6) Displacement lollipop
    if (nrow(disp_tbl) > 0 && all(c("Facet", "Level", "Displacement") %in% names(disp_tbl))) {
      ord <- order(abs(suppressWarnings(as.numeric(disp_tbl$Displacement))), decreasing = TRUE, na.last = NA)
      use <- ord[seq_len(min(length(ord), top_n))]
      sub <- disp_tbl[use, , drop = FALSE]
      y <- seq_len(nrow(sub))
      lbl <- truncate_axis_label(paste0(sub$Facet, ":", sub$Level), width = 24L)
      disp_vals <- suppressWarnings(as.numeric(sub$Displacement))
      cols <- if ("Flag" %in% names(sub)) ifelse(as.logical(sub$Flag), "#cb181d", "#238b45") else "#238b45"
      graphics::plot(
        x = disp_vals,
        y = y,
        type = "n",
        xlab = "Displacement",
        ylab = "",
        yaxt = "n",
        main = "QC: Displacement"
      )
      graphics::segments(0, y, disp_vals, y, col = "gray60")
      graphics::points(disp_vals, y, pch = 16, col = cols)
      graphics::axis(side = 2, at = y, labels = lbl, las = 2, cex.axis = 0.7)
      d_thr <- as.numeric(disp$thresholds$abs_displacement_warn %||% 0.5)
      graphics::abline(v = c(-d_thr, 0, d_thr), lty = c(2, 1, 2), col = c("gray45", "gray30", "gray45"))
    } else {
      graphics::plot.new()
      graphics::title(main = "QC: Displacement")
      graphics::text(0.5, 0.5, "No data")
    }

    # 7) Inter-rater exact agreement
    if (nrow(inter_tbl) > 0 && all(c("Rater1", "Rater2", "Exact") %in% names(inter_tbl))) {
      ord <- order(suppressWarnings(as.numeric(inter_tbl$Exact)), na.last = NA)
      use <- ord[seq_len(min(length(ord), top_n))]
      sub <- inter_tbl[use, , drop = FALSE]
      pair_lbl <- truncate_axis_label(paste0(sub$Rater1, " | ", sub$Rater2), width = 20L)
      cols <- if ("Flag" %in% names(sub)) ifelse(as.logical(sub$Flag), "#cb181d", "#2b8cbe") else "#2b8cbe"
      bp <- barplot_rot45(
        height = suppressWarnings(as.numeric(sub$Exact)),
        labels = pair_lbl,
        col = cols,
        main = "QC: Inter-rater exact",
        ylab = "Exact agreement",
        label_angle = 45,
        mar_bottom = 6.4,
        label_cex = 0.72,
        label_width = 14L
      )
      exp_vals <- suppressWarnings(as.numeric(sub$ExpectedExact))
      if (any(is.finite(exp_vals))) {
        graphics::points(bp, exp_vals, pch = 21, bg = "white", col = "#08519c")
        graphics::lines(bp, exp_vals, col = "#08519c", lwd = 1.3)
      }
      graphics::abline(h = interrater_exact_warn, lty = 2, col = "gray45")
    } else {
      graphics::plot.new()
      graphics::title(main = "QC: Inter-rater exact")
      graphics::text(0.5, 0.5, "No data")
    }

    # 8) Facet fixed-effect chi-square
    if (nrow(fchi_tbl) > 0 && all(c("Facet", "FixedChiSq") %in% names(fchi_tbl))) {
      ord <- order(suppressWarnings(as.numeric(fchi_tbl$FixedChiSq)), decreasing = TRUE, na.last = NA)
      sub <- fchi_tbl[ord, , drop = FALSE]
      labels <- truncate_axis_label(as.character(sub$Facet), width = 20L)
      cols <- if ("FixedFlag" %in% names(sub)) ifelse(as.logical(sub$FixedFlag), "#cb181d", "#31a354") else "#31a354"
      barplot_rot45(
        height = suppressWarnings(as.numeric(sub$FixedChiSq)),
        labels = labels,
        col = cols,
        main = "QC: Facet fixed chi-square",
        ylab = expression(chi^2),
        label_angle = 45,
        mar_bottom = 6.4,
        label_cex = 0.72,
        label_width = 14L
      )
    } else {
      graphics::plot.new()
      graphics::title(main = "QC: Facet fixed chi-square")
      graphics::text(0.5, 0.5, "No data")
    }

    # 9) Separation reliability
    if (nrow(rel_tbl) > 0 && all(c("Facet", "Separation") %in% names(rel_tbl))) {
      ord <- order(suppressWarnings(as.numeric(rel_tbl$Separation)), decreasing = TRUE, na.last = NA)
      sub <- rel_tbl[ord, , drop = FALSE]
      barplot_rot45(
        height = suppressWarnings(as.numeric(sub$Separation)),
        labels = truncate_axis_label(as.character(sub$Facet), width = 20L),
        col = "#9e9ac8",
        main = "QC: Separation by facet",
        ylab = "Separation",
        label_angle = 45,
        mar_bottom = 6.4,
        label_cex = 0.72,
        label_width = 14L
      )
      graphics::abline(h = 1, lty = 2, col = "gray45")
    } else {
      graphics::plot.new()
      graphics::title(main = "QC: Separation by facet")
      graphics::text(0.5, 0.5, "No data")
    }
  }

  out <- new_mfrm_plot_data(
    "qc_dashboard",
    list(
      threshold_profile = resolved$profile_name,
      thresholds = resolved$thresholds,
      category_stats = cat_tbl,
      fit = fit_tbl,
      zstd = zstd,
      unexpected = unexpected,
      fair_average = fair,
      displacement = disp,
      interrater = interrater,
      facets_chisq = fchi,
      reliability = rel_tbl
    )
  )
  invisible(out)
}

#' Plot an `mfrm_fit` object
#'
#' @param x Output from [fit_mfrm()].
#' @param type Plot type.
#' If omitted (`NULL`), returns a bundle with Wright map, Pathway map,
#' and category characteristic curves.
#' Single-plot options are `"facet"`, `"person"`, `"step"`, `"wright"`,
#' `"pathway"`, and `"ccc"`.
#' @param facet Optional facet name when `type = "facet"`.
#' @param top_n Maximum number of facet levels shown for `type = "facet"`.
#' @param theta_range Theta/logit plotting range for pathway and CCC plots.
#' @param theta_points Number of grid points for pathway and CCC curves.
#' @param title Optional custom title for the selected plot.
#' @param palette Optional named color overrides.
#' Supported names include `facet_level`, `step_threshold`, `person_hist`,
#' `step_line`, `facet_bar`, and `grid`.
#' @param label_angle X-axis label angle for categorical axes.
#' @param show_ci If `TRUE`, draws confidence-interval whiskers on
#'   Wright map facet/step locations and facet bar plots.
#'   Requires SE information computed from the observation table.
#' @param ci_level Confidence level for whiskers (default 0.95).
#' @param draw If `TRUE`, draws with base graphics.
#' @param ... Reserved for generic compatibility.
#'
#' @details
#' `plot(fit)` returns a named list with three plotting-data objects:
#' - `wright_map`
#' - `pathway_map`
#' - `category_characteristic_curves`
#'
#' Palette behavior:
#' - For `"wright"`/bundle Wright panel, use named keys such as
#'   `facet_level`, `step_threshold`, `person_hist`, `grid`.
#' - For `"step"`, use `step_line` (and optionally `grid`).
#' - For `"person"`, use `person_hist`.
#' - For `"pathway"` and `"ccc"`, unnamed or named vectors are matched to
#'   detected curve groups in plotting order.
#'
#' If `draw = TRUE`, base R graphics are produced.
#'
#' @section Interpreting output:
#' - bundle mode (`type = NULL`) returns all three core plotting payloads.
#' - single mode (`type = "facet"`, `"person"`, `"step"`, `"wright"`,
#'   `"pathway"`, `"ccc"`) returns one plotting payload.
#'
#' @section Typical workflow:
#' 1. Run `plot(fit, draw = FALSE)` to collect reusable plot data.
#' 2. Switch to one panel with `type = ...` for focused review.
#' 3. Set `draw = TRUE` with palette options for direct base-R rendering.
#'
#' @return
#' If `type` is `NULL`, a named list (`mfrm_plot_bundle`) containing plotting data.
#' Otherwise, a single plotting-data object (`mfrm_plot_data`).
#' @seealso [fit_mfrm()], [summary.mfrm_fit()], [plot_residual_pca()]
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' maps <- plot(fit, draw = FALSE)
#' p1 <- plot(fit, type = "wright", draw = FALSE)
#' p2 <- plot(fit, type = "pathway", draw = FALSE)
#' p3 <- plot(fit, type = "ccc", draw = FALSE)
#' if (interactive()) {
#'   plot(
#'     fit,
#'     type = "wright",
#'     title = "Customized Wright Map",
#'     palette = c(facet_level = "#1f78b4", step_threshold = "#d95f02"),
#'     label_angle = 45
#'   )
#'   plot(
#'     fit,
#'     type = "step",
#'     title = "Customized Step Parameters",
#'     palette = c(step_line = "#008b8b"),
#'     label_angle = 45
#'   )
#'   plot(
#'     fit,
#'     type = "pathway",
#'     title = "Customized Pathway Map",
#'     palette = c("#1f78b4")
#'   )
#'   plot(
#'     fit,
#'     type = "ccc",
#'     title = "Customized Category Characteristic Curves",
#'     palette = c("#1b9e77", "#d95f02", "#7570b3")
#'   )
#' }
#' @export
plot.mfrm_fit <- function(x,
                          type = NULL,
                          facet = NULL,
                          top_n = 30,
                          theta_range = c(-6, 6),
                          theta_points = 241,
                          title = NULL,
                          palette = NULL,
                          label_angle = 45,
                          show_ci = FALSE,
                          ci_level = 0.95,
                          draw = TRUE,
                          ...) {
  if (!inherits(x, "mfrm_fit")) {
    stop("`x` must be an mfrm_fit object from fit_mfrm().")
  }
  top_n <- max(1L, as.integer(top_n))
  theta_points <- max(51L, as.integer(theta_points))
  theta_range <- as.numeric(theta_range)
  if (length(theta_range) != 2 || !all(is.finite(theta_range)) || theta_range[1] >= theta_range[2]) {
    stop("`theta_range` must be a numeric length-2 vector with increasing values.")
  }

  as_plot_data <- function(name, data) {
    out <- list(name = name, data = data)
    class(out) <- c("mfrm_plot_data", class(out))
    out
  }

  se_tbl_ci <- if (isTRUE(show_ci)) compute_se_for_plot(x, ci_level = ci_level) else NULL

  default_bundle <- missing(type) || is.null(type)
  if (default_bundle || tolower(as.character(type[1])) %in% c("bundle", "all", "default")) {
    out <- list(
      wright_map = as_plot_data("wright_map", build_wright_map_data(x, top_n = top_n, se_tbl = se_tbl_ci)),
      pathway_map = as_plot_data("pathway_map", build_pathway_map_data(x, theta_range = theta_range, theta_points = theta_points)),
      category_characteristic_curves = as_plot_data("category_characteristic_curves", build_ccc_data(x, theta_range = theta_range, theta_points = theta_points))
    )
    class(out) <- c("mfrm_plot_bundle", class(out))
    if (isTRUE(draw)) {
      draw_wright_map(out$wright_map$data, title = title, palette = palette, label_angle = label_angle,
                      show_ci = show_ci, ci_level = ci_level)
      draw_pathway_map(out$pathway_map$data, title = title, palette = palette)
      draw_ccc(out$category_characteristic_curves$data, title = title, palette = palette)
    }
    return(invisible(out))
  }

  type <- match.arg(tolower(as.character(type[1])), c("facet", "person", "step", "wright", "pathway", "ccc"))

  if (type == "wright") {
    out <- as_plot_data("wright_map", build_wright_map_data(x, top_n = top_n, se_tbl = se_tbl_ci))
    if (isTRUE(draw)) draw_wright_map(out$data, title = title, palette = palette, label_angle = label_angle,
                                       show_ci = show_ci, ci_level = ci_level)
    return(invisible(out))
  }
  if (type == "pathway") {
    out <- as_plot_data("pathway_map", build_pathway_map_data(x, theta_range = theta_range, theta_points = theta_points))
    if (isTRUE(draw)) draw_pathway_map(out$data, title = title, palette = palette)
    return(invisible(out))
  }
  if (type == "ccc") {
    out <- as_plot_data("category_characteristic_curves", build_ccc_data(x, theta_range = theta_range, theta_points = theta_points))
    if (isTRUE(draw)) draw_ccc(out$data, title = title, palette = palette)
    return(invisible(out))
  }
  if (type == "person") {
    person_tbl <- tibble::as_tibble(x$facets$person)
    person_tbl <- person_tbl[is.finite(person_tbl$Estimate), , drop = FALSE]
    if (nrow(person_tbl) == 0) stop("No finite person estimates available for plotting.")
    bins <- max(10L, min(35L, as.integer(round(sqrt(nrow(person_tbl))))))
    this_title <- if (is.null(title)) "Person measure distribution" else as.character(title[1])
    out <- as_plot_data("person", list(person = person_tbl, bins = bins, title = this_title))
    if (isTRUE(draw)) draw_person_plot(person_tbl, bins = bins, title = this_title, palette = palette)
    return(invisible(out))
  }
  if (type == "step") {
    step_tbl <- tibble::as_tibble(x$steps)
    if (nrow(step_tbl) == 0 || !all(c("Step", "Estimate") %in% names(step_tbl))) {
      stop("Step estimates are not available in this fit object.")
    }
    this_title <- if (is.null(title)) "Step parameter estimates" else as.character(title[1])
    out <- as_plot_data("step", list(steps = step_tbl, title = this_title))
    if (isTRUE(draw)) draw_step_plot(step_tbl, title = this_title, palette = palette, label_angle = label_angle)
    return(invisible(out))
  }

  facet_tbl <- tibble::as_tibble(x$facets$others)
  if (nrow(facet_tbl) == 0 || !all(c("Facet", "Level", "Estimate") %in% names(facet_tbl))) {
    stop("Facet-level estimates are not available in this fit object.")
  }
  if (!is.null(facet)) {
    facet_tbl <- dplyr::filter(facet_tbl, .data$Facet == as.character(facet[1]))
    if (nrow(facet_tbl) == 0) stop("Requested `facet` not found in facet-level estimates.")
  }
  if (isTRUE(show_ci) && !is.null(se_tbl_ci) && is.data.frame(se_tbl_ci) &&
      nrow(se_tbl_ci) > 0 && all(c("Facet", "Level", "SE") %in% names(se_tbl_ci))) {
    se_join <- se_tbl_ci[, intersect(c("Facet", "Level", "SE"), names(se_tbl_ci)), drop = FALSE]
    se_join$Level <- as.character(se_join$Level)
    facet_tbl$Level <- as.character(facet_tbl$Level)
    facet_tbl <- merge(facet_tbl, se_join, by = c("Facet", "Level"), all.x = TRUE, sort = FALSE)
    facet_tbl <- tibble::as_tibble(facet_tbl)
  }
  facet_tbl <- dplyr::arrange(facet_tbl, .data$Estimate)
  if (nrow(facet_tbl) > top_n) {
    facet_tbl <- facet_tbl |>
      dplyr::mutate(AbsEstimate = abs(.data$Estimate)) |>
      dplyr::arrange(dplyr::desc(.data$AbsEstimate)) |>
      dplyr::slice_head(n = top_n) |>
      dplyr::arrange(.data$Estimate) |>
      dplyr::select(-.data$AbsEstimate)
  }
  facet_title <- if (is.null(facet)) "Facet-level estimates" else paste0("Facet-level estimates: ", as.character(facet[1]))
  if (!is.null(title)) facet_title <- as.character(title[1])
  out <- as_plot_data("facet", list(facets = facet_tbl, title = facet_title))
  if (isTRUE(draw)) draw_facet_plot(facet_tbl, title = facet_title, palette = palette,
                                     show_ci = show_ci, ci_level = ci_level)
  invisible(out)
}

# ---- Bubble Chart ----

resolve_bubble_measures <- function(x, diagnostics = NULL) {
  if (inherits(x, "mfrm_diagnostics") ||
      (is.list(x) && "measures" %in% names(x) && is.data.frame(x$measures))) {
    return(as.data.frame(x$measures, stringsAsFactors = FALSE))
  }
  if (inherits(x, "mfrm_fit")) {
    if (is.null(diagnostics)) {
      diagnostics <- diagnose_mfrm(x, residual_pca = "none")
    }
    return(as.data.frame(diagnostics$measures, stringsAsFactors = FALSE))
  }
  stop("`x` must be an mfrm_fit object or output from diagnose_mfrm().")
}

#' Bubble chart of measure estimates and fit statistics
#'
#' Produces a Rasch-convention bubble chart where each element is a circle
#' positioned at its measure estimate (x) and fit mean-square (y).
#' Bubble radius reflects measurement precision or sample size.
#'
#' @param x Output from \code{\link{fit_mfrm}} or \code{\link{diagnose_mfrm}}.
#' @param diagnostics Optional output from \code{\link{diagnose_mfrm}} when
#'   \code{x} is an \code{mfrm_fit} object. If omitted, diagnostics are
#'   computed automatically.
#' @param fit_stat Fit statistic for the y-axis: \code{"Infit"} (default) or
#'   \code{"Outfit"}.
#' @param bubble_size Variable controlling bubble radius: \code{"SE"} (default),
#'   \code{"N"} (observation count), or \code{"equal"} (uniform size).
#' @param facets Character vector of facets to include. \code{NULL} (default)
#'   includes all non-person facets.
#' @param fit_range Numeric length-2 vector defining the acceptable fit range
#'   shown as a shaded band (default \code{c(0.5, 1.5)}).
#' @param top_n Maximum number of elements to plot (default 60).
#' @param main Optional custom plot title.
#' @param palette Optional named colour vector keyed by facet name.
#' @param draw If \code{TRUE} (default), render the plot using base graphics.
#'
#' @section Interpreting the plot:
#' Points near the horizontal reference line at 1.0 fit the model well.
#' Points above 1.5 show underfit (unpredictable response patterns).
#' Points below 0.5 show overfit (Guttman-like patterns).
#'
#' @return Invisibly, an object of class \code{mfrm_plot_data}.
#' @seealso \code{\link{diagnose_mfrm}}, \code{\link{plot_unexpected}},
#'   \code{\link{plot_fair_average}}
#' @export
plot_bubble <- function(x,
                        diagnostics = NULL,
                        fit_stat = c("Infit", "Outfit"),
                        bubble_size = c("SE", "N", "equal"),
                        facets = NULL,
                        fit_range = c(0.5, 1.5),
                        top_n = 60,
                        main = NULL,
                        palette = NULL,
                        draw = TRUE) {
  fit_stat <- match.arg(fit_stat)
  bubble_size <- match.arg(bubble_size)
  top_n <- max(1L, as.integer(top_n))

  measures <- resolve_bubble_measures(x, diagnostics)
  measures <- measures[measures$Facet != "Person", , drop = FALSE]
  if (!is.null(facets)) {
    measures <- measures[measures$Facet %in% as.character(facets), , drop = FALSE]
  }
  if (nrow(measures) == 0) stop("No measures available for bubble chart.")

  needed <- c("Facet", "Level", "Estimate", fit_stat)
  missing_cols <- setdiff(needed, names(measures))
  if (length(missing_cols) > 0) {
    stop("Missing columns in measures: ", paste(missing_cols, collapse = ", "))
  }

  ok <- is.finite(measures$Estimate) & is.finite(measures[[fit_stat]])
  measures <- measures[ok, , drop = FALSE]
  if (nrow(measures) == 0) stop("No finite measure/fit values for bubble chart.")

  if (nrow(measures) > top_n) {
    measures <- measures[order(abs(measures[[fit_stat]] - 1), decreasing = TRUE), ]
    measures <- measures[seq_len(top_n), , drop = FALSE]
  }

  radius <- switch(bubble_size,
    SE = {
      se_vals <- if ("SE" %in% names(measures)) measures$SE else rep(0.1, nrow(measures))
      se_vals[!is.finite(se_vals)] <- stats::median(se_vals[is.finite(se_vals)], na.rm = TRUE)
      se_vals / max(se_vals, na.rm = TRUE) * 0.15
    },
    N = {
      n_vals <- if ("N" %in% names(measures)) measures$N else rep(1, nrow(measures))
      n_vals[!is.finite(n_vals)] <- 1
      sqrt(n_vals) / max(sqrt(n_vals), na.rm = TRUE) * 0.15
    },
    equal = rep(0.08, nrow(measures))
  )

  unique_facets <- unique(measures$Facet)
  default_cols <- stats::setNames(
    grDevices::hcl.colors(max(3L, length(unique_facets)), "Dark 3")[seq_along(unique_facets)],
    unique_facets
  )
  cols <- resolve_palette(palette = palette, defaults = default_cols)
  point_cols <- cols[as.character(measures$Facet)]

  if (isTRUE(draw)) {
    xr <- range(measures$Estimate, na.rm = TRUE)
    xr <- xr + diff(xr) * c(-0.15, 0.15)
    yr <- range(c(measures[[fit_stat]], fit_range), na.rm = TRUE)
    yr <- yr + diff(yr) * c(-0.1, 0.1)

    graphics::plot(
      x = measures$Estimate, y = measures[[fit_stat]], type = "n",
      xlim = xr, ylim = yr,
      xlab = "Measure (logits)",
      ylab = paste0(fit_stat, " Mean Square"),
      main = if (is.null(main)) paste0("Bubble Chart: ", fit_stat) else as.character(main[1])
    )
    graphics::rect(
      xleft = xr[1] - 1, ybottom = fit_range[1],
      xright = xr[2] + 1, ytop = fit_range[2],
      col = grDevices::adjustcolor("#d9f0d3", alpha.f = 0.4), border = NA
    )
    graphics::abline(h = 1, lty = 2, col = "gray40", lwd = 1.5)
    graphics::abline(h = fit_range, lty = 3, col = "gray60")
    graphics::symbols(
      x = measures$Estimate, y = measures[[fit_stat]],
      circles = radius, inches = FALSE, add = TRUE,
      fg = point_cols,
      bg = grDevices::adjustcolor(point_cols, alpha.f = 0.45)
    )
    graphics::legend(
      "topleft", legend = unique_facets,
      col = cols[unique_facets], pch = 16, bty = "n", cex = 0.85
    )
  }

  out <- new_mfrm_plot_data(
    "bubble",
    list(fit_stat = fit_stat, bubble_size = bubble_size,
         fit_range = fit_range, table = measures, radius = radius)
  )
  invisible(out)
}

# ---- CSV Export ----

#' Export MFRM results to CSV files
#'
#' Writes tidy CSV files suitable for import into spreadsheet software or
#' further analysis in other tools.
#'
#' @param fit Output from \code{\link{fit_mfrm}}.
#' @param diagnostics Optional output from \code{\link{diagnose_mfrm}}.
#'   When provided, enriches facet estimates with SE, fit statistics, and
#'   writes the full measures table.
#' @param output_dir Directory for CSV files. Created if it does not exist.
#' @param prefix Filename prefix (default \code{"mfrm"}).
#' @param tables Character vector of tables to export. Any subset of
#'   \code{"person"}, \code{"facets"}, \code{"summary"}, \code{"steps"},
#'   \code{"measures"}. Default exports all available tables.
#' @param overwrite If \code{FALSE} (default), refuse to overwrite existing
#'   files.
#'
#' @section Exported files:
#' \describe{
#'   \item{\code{{prefix}_person_estimates.csv}}{Person ID, Estimate, SD.}
#'   \item{\code{{prefix}_facet_estimates.csv}}{Facet, Level, Estimate,
#'     and optionally SE, Infit, Outfit, PTMEA when diagnostics supplied.}
#'   \item{\code{{prefix}_fit_summary.csv}}{One-row model summary.}
#'   \item{\code{{prefix}_step_parameters.csv}}{Step/threshold parameters.}
#'   \item{\code{{prefix}_measures.csv}}{Full measures table (requires
#'     diagnostics).}
#' }
#'
#' @return Invisibly, a data.frame listing written files with columns
#'   \code{Table} and \code{Path}.
#' @seealso \code{\link{fit_mfrm}}, \code{\link{diagnose_mfrm}},
#'   \code{\link{as.data.frame.mfrm_fit}}
#' @export
export_mfrm <- function(fit,
                        diagnostics = NULL,
                        output_dir = ".",
                        prefix = "mfrm",
                        tables = c("person", "facets", "summary", "steps", "measures"),
                        overwrite = FALSE) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().")
  }
  tables <- unique(tolower(as.character(tables)))
  allowed <- c("person", "facets", "summary", "steps", "measures")
  bad <- setdiff(tables, allowed)
  if (length(bad) > 0) {
    stop("Unknown table names: ", paste(bad, collapse = ", "),
         ". Allowed: ", paste(allowed, collapse = ", "))
  }
  prefix <- as.character(prefix[1])
  if (!nzchar(prefix)) prefix <- "mfrm"
  overwrite <- isTRUE(overwrite)
  output_dir <- as.character(output_dir[1])

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(output_dir)) {
    stop("Could not create output directory: ", output_dir)
  }

  written <- data.frame(Table = character(0), Path = character(0),
                        stringsAsFactors = FALSE)

  write_one <- function(df, filename, table_name) {
    path <- file.path(output_dir, filename)
    if (file.exists(path) && !overwrite) {
      stop("File already exists: ", path, ". Set overwrite = TRUE to replace.")
    }
    utils::write.csv(df, file = path, row.names = FALSE, na = "")
    written <<- rbind(written, data.frame(Table = table_name, Path = path,
                                          stringsAsFactors = FALSE))
  }

  if ("person" %in% tables) {
    person_df <- as.data.frame(fit$facets$person, stringsAsFactors = FALSE)
    write_one(person_df, paste0(prefix, "_person_estimates.csv"), "person")
  }

  if ("facets" %in% tables) {
    facet_df <- as.data.frame(fit$facets$others, stringsAsFactors = FALSE)
    if (!is.null(diagnostics) && !is.null(diagnostics$measures)) {
      enrich_cols <- intersect(c("SE", "Infit", "Outfit", "PTMEA", "N"),
                               names(diagnostics$measures))
      if (length(enrich_cols) > 0) {
        enrich <- diagnostics$measures[diagnostics$measures$Facet != "Person",
                                       c("Facet", "Level", enrich_cols), drop = FALSE]
        enrich <- as.data.frame(enrich, stringsAsFactors = FALSE)
        enrich$Level <- as.character(enrich$Level)
        facet_df$Level <- as.character(facet_df$Level)
        facet_df <- merge(facet_df, enrich, by = c("Facet", "Level"), all.x = TRUE)
      }
    }
    write_one(facet_df, paste0(prefix, "_facet_estimates.csv"), "facets")
  }

  if ("summary" %in% tables) {
    summary_df <- as.data.frame(fit$summary, stringsAsFactors = FALSE)
    write_one(summary_df, paste0(prefix, "_fit_summary.csv"), "summary")
  }

  if ("steps" %in% tables) {
    step_df <- as.data.frame(fit$steps, stringsAsFactors = FALSE)
    if (nrow(step_df) > 0) {
      write_one(step_df, paste0(prefix, "_step_parameters.csv"), "steps")
    }
  }

  if ("measures" %in% tables && !is.null(diagnostics) && !is.null(diagnostics$measures)) {
    measures_df <- as.data.frame(diagnostics$measures, stringsAsFactors = FALSE)
    write_one(measures_df, paste0(prefix, "_measures.csv"), "measures")
  }

  invisible(written)
}

#' Convert mfrm_fit to a tidy data.frame
#'
#' Returns all facet-level estimates (person and others) in a single
#' tidy data.frame. Useful for quick interactive export:
#' \code{write.csv(as.data.frame(fit), "results.csv")}.
#'
#' @param x An \code{mfrm_fit} object from \code{\link{fit_mfrm}}.
#' @param row.names Ignored (included for S3 generic compatibility).
#' @param optional Ignored (included for S3 generic compatibility).
#' @param ... Additional arguments (ignored).
#'
#' @return A data.frame with columns \code{Facet}, \code{Level},
#'   \code{Estimate}.
#' @seealso \code{\link{fit_mfrm}}, \code{\link{export_mfrm}}
#' @export
as.data.frame.mfrm_fit <- function(x, row.names = NULL, optional = FALSE, ...) {
  person_df <- data.frame(
    Facet = "Person",
    Level = as.character(x$facets$person$Person),
    Estimate = x$facets$person$Estimate,
    stringsAsFactors = FALSE
  )
  facet_df <- as.data.frame(
    x$facets$others[, c("Facet", "Level", "Estimate")],
    stringsAsFactors = FALSE
  )
  facet_df$Level <- as.character(facet_df$Level)
  rbind(person_df, facet_df)
}

#' @export
print.mfrm_plot_bundle <- function(x, ...) {
  cat("mfrm plot bundle\n")
  cat("  - wright_map\n")
  cat("  - pathway_map\n")
  cat("  - category_characteristic_curves\n")
  cat("Use `$` to access each plotting-data object.\n")
  invisible(x)
}

#' @export
print.mfrm_fit <- function(x, ...) {
  if (is.list(x) && !is.null(x$summary) && nrow(x$summary) > 0) {
    cat("mfrm_fit object\n")
    print(x$summary)
  } else {
    cat("mfrm_fit object (empty summary)\n")
  }
  invisible(x)
}
