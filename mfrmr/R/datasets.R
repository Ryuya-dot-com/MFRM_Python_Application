#' Simulated MFRM datasets based on Eckes and Jin (2021)
#'
#' Synthetic many-facet rating datasets in long format.
#' All datasets include one row per observed rating.
#'
#' Available data objects:
#' - `ej2021_study1`
#' - `ej2021_study2`
#' - `ej2021_combined`
#' - `ej2021_study1_itercal`
#' - `ej2021_study2_itercal`
#' - `ej2021_combined_itercal`
#'
#' @format A data.frame with 5 columns:
#' \describe{
#'   \item{Study}{Study label (`"Study1"` or `"Study2"`).}
#'   \item{Person}{Person/respondent identifier.}
#'   \item{Rater}{Rater identifier.}
#'   \item{Criterion}{Criterion facet label.}
#'   \item{Score}{Observed category score.}
#' }
#' @source Simulated for this package with design settings informed by Eckes and Jin (2021).
#' @details
#' Naming convention:
#' - `study1` / `study2`: separate simulation studies
#' - `combined`: row-bind of study1 and study2
#' - `_itercal`: iterative-calibration variant
#'
#' Use [load_mfrmr_data()] for programmatic selection by key.
#'
#' @section Data dimensions:
#' \tabular{lrrrr}{
#'   \strong{Dataset} \tab \strong{Rows} \tab \strong{Persons} \tab \strong{Raters} \tab \strong{Criteria} \cr
#'   study1 \tab 2160 \tab 60 \tab 6 \tab 6 \cr
#'   study2 \tab 2160 \tab 60 \tab 6 \tab 6 \cr
#'   combined \tab 4320 \tab 120 \tab 6 \tab 6 \cr
#'   study1_itercal \tab 2160 \tab 60 \tab 6 \tab 6 \cr
#'   study2_itercal \tab 2160 \tab 60 \tab 6 \tab 6 \cr
#'   combined_itercal \tab 4320 \tab 120 \tab 6 \tab 6
#' }
#' Score range: 1--5 (five-category rating scale).
#'
#' @section Simulation design:
#' Person ability is drawn from N(0, 1).  Rater severity effects span
#' approximately -0.5 to +0.5 logits.  Criterion difficulty effects span
#' approximately -0.3 to +0.3 logits.  Scores are generated from the
#' resulting linear predictor plus Gaussian noise, then discretized into
#' five categories.  The `_itercal` variants use a second iteration of
#' calibrated rater severity parameters.
#'
#' @section Interpreting output:
#' Each dataset is already in long format and can be passed directly to
#' [fit_mfrm()] after confirming column-role mapping.
#'
#' @section Typical workflow:
#' 1. Inspect available datasets with [list_mfrmr_data()].
#' 2. Load one dataset using [load_mfrmr_data()].
#' 3. Fit and diagnose with [fit_mfrm()] and [diagnose_mfrm()].
#'
#' @examples
#' data("ej2021_study1", package = "mfrmr")
#' head(ej2021_study1)
#' table(ej2021_study1$Study)
#' @name ej2021_data
#' @aliases ej2021_study1 ej2021_study2 ej2021_combined
#' @aliases ej2021_study1_itercal ej2021_study2_itercal ej2021_combined_itercal
NULL

#' List packaged simulation datasets
#'
#' @return Character vector of dataset keys accepted by [load_mfrmr_data()].
#' @details
#' Use this helper when you want to select packaged data programmatically
#' (e.g., inside scripts, loops, or shiny/streamlit wrappers).
#'
#' Typical pattern:
#' 1. call `list_mfrmr_data()` to see available keys.
#' 2. pass one key to [load_mfrmr_data()].
#'
#' @section Interpreting output:
#' Returned values are canonical dataset keys accepted by [load_mfrmr_data()].
#'
#' @section Typical workflow:
#' 1. Capture keys in a script (`keys <- list_mfrmr_data()`).
#' 2. Select one key by index or name.
#' 3. Load data via [load_mfrmr_data()] and continue analysis.
#'
#' @seealso [load_mfrmr_data()], [ej2021_data]
#' @examples
#' keys <- list_mfrmr_data()
#' keys
#' d <- load_mfrmr_data(keys[1])
#' head(d)
#' @export
list_mfrmr_data <- function() {
  c(
    "study1",
    "study2",
    "combined",
    "study1_itercal",
    "study2_itercal",
    "combined_itercal"
  )
}

#' Load a packaged simulation dataset
#'
#' @param name Dataset key. One of values from [list_mfrmr_data()].
#'
#' @return A data.frame in long format.
#' @details
#' This helper is useful in scripts/functions where you want to choose a dataset
#' by string key instead of calling `data()` manually.
#'
#' Returned columns are always:
#' `Study`, `Person`, `Rater`, `Criterion`, `Score`.
#'
#' @section Interpreting output:
#' The return value is a plain long-format `data.frame`, ready for direct use
#' in [fit_mfrm()] without additional reshaping.
#'
#' @section Typical workflow:
#' 1. list valid names with [list_mfrmr_data()].
#' 2. load one dataset key with `load_mfrmr_data(name)`.
#' 3. fit a model with [fit_mfrm()] and inspect with `summary()` / `plot()`.
#'
#' @seealso [list_mfrmr_data()], [ej2021_data]
#' @examples
#' d <- load_mfrmr_data("study1")
#' head(d)
#' fit <- fit_mfrm(
#'   data = d,
#'   person = "Person",
#'   facets = c("Rater", "Criterion"),
#'   score = "Score",
#'   method = "JML",
#'   maxit = 25
#' )
#' summary(fit)
#' @export
load_mfrmr_data <- function(name = c(
                            "study1",
                            "study2",
                            "combined",
                            "study1_itercal",
                            "study2_itercal",
                            "combined_itercal"
                          )) {
  key <- match.arg(tolower(name), choices = list_mfrmr_data())

  obj_name <- switch(
    key,
    study1 = "ej2021_study1",
    study2 = "ej2021_study2",
    combined = "ej2021_combined",
    study1_itercal = "ej2021_study1_itercal",
    study2_itercal = "ej2021_study2_itercal",
    combined_itercal = "ej2021_combined_itercal"
  )

  utils::data(list = obj_name, package = "mfrmr", envir = environment())
  get(obj_name, envir = environment(), inherits = FALSE)
}
