#!/usr/bin/env Rscript

# Build package data/ objects from CSV files in inst/extdata.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- sub("^--file=", "", file_arg[1])
pkg_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
ext_dir <- file.path(pkg_root, "inst", "extdata")
data_dir <- file.path(pkg_root, "data")

dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

csv_map <- list(
  ej2021_study1 = "eckes_jin_2021_study1_sim.csv",
  ej2021_study2 = "eckes_jin_2021_study2_sim.csv",
  ej2021_combined = "eckes_jin_2021_combined_sim.csv",
  ej2021_study1_itercal = "eckes_jin_2021_study1_itercal_sim.csv",
  ej2021_study2_itercal = "eckes_jin_2021_study2_itercal_sim.csv",
  ej2021_combined_itercal = "eckes_jin_2021_combined_itercal_sim.csv"
)

for (obj_name in names(csv_map)) {
  csv_file <- file.path(ext_dir, csv_map[[obj_name]])
  if (!file.exists(csv_file)) {
    stop("Missing CSV: ", csv_file)
  }

  obj <- utils::read.csv(csv_file, stringsAsFactors = FALSE)
  assign(obj_name, obj, envir = environment())

  save(
    list = obj_name,
    file = file.path(data_dir, paste0(obj_name, ".rda")),
    envir = environment(),
    compress = "bzip2"
  )
}

message("Saved data objects to: ", data_dir)
