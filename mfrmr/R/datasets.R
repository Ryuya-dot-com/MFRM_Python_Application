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
#' @name ej2021_data
#' @aliases ej2021_study1 ej2021_study2 ej2021_combined
#' @aliases ej2021_study1_itercal ej2021_study2_itercal ej2021_combined_itercal
NULL

#' List packaged simulation datasets
#'
#' @return Character vector of dataset keys accepted by [load_mfrmr_data()].
#' @seealso [load_mfrmr_data()], [ej2021_data]
#' @examples
#' list_mfrmr_data()
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
#' @seealso [list_mfrmr_data()], [ej2021_data]
#' @examples
#' d <- load_mfrmr_data("study1")
#' head(d)
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
