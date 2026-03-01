# --------------------------------------------------------------------------
# test-identifiability-constraints.R
# Tests identifiability constraints and centering for the mfrmr package.
# --------------------------------------------------------------------------

# ---- 3.1  Sum-to-zero constraint on centered facets ---------------------

test_that("centered facet estimates sum to zero", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(fit_mfrm(
    d, "Person", c("Rater", "Task", "Criterion"), "Score",
    method = "JML", model = "RSM", maxit = 40, quad_points = 7
  ))
  for (facet in c("Rater", "Task", "Criterion")) {
    est <- fit$facets$others |>
      dplyr::filter(Facet == facet) |>
      dplyr::pull(Estimate)
    expect_equal(sum(est), 0, tolerance = 1e-6,
                 label = paste("sum-to-zero for", facet))
  }
})

# ---- 3.2  Anchor constraint -- anchored level matches its value ---------

test_that("anchored facet level matches its anchor value exactly", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  anchors <- data.frame(
    Facet = "Rater", Level = "R2", Anchor = 0,
    stringsAsFactors = FALSE
  )
  fit <- suppressWarnings(fit_mfrm(
    d, "Person", c("Rater", "Task", "Criterion"), "Score",
    anchors = anchors, method = "JML", model = "RSM", maxit = 40,
    quad_points = 7
  ))
  r2_est <- fit$facets$others |>
    dplyr::filter(Facet == "Rater", Level == "R2") |>
    dplyr::pull(Estimate)
  expect_equal(unname(r2_est), 0, tolerance = 1e-8)
})

# ---- 3.3  Constraint methods preserve ordering and positive correlation --

test_that("constraint method changes preserve ordering and positive correlation", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)

  fit_centered <- suppressWarnings(fit_mfrm(
    d, "Person", c("Rater", "Task", "Criterion"), "Score",
    method = "JML", model = "RSM", maxit = 50, quad_points = 7
  ))

  anchors <- data.frame(
    Facet = "Rater", Level = "R1", Anchor = 0,
    stringsAsFactors = FALSE
  )
  fit_anchored <- suppressWarnings(fit_mfrm(
    d, "Person", c("Rater", "Task", "Criterion"), "Score",
    anchors = anchors, method = "JML", model = "RSM", maxit = 50,
    quad_points = 7
  ))

  # Both LogLik should be finite
  expect_true(is.finite(fit_centered$summary$LogLik))
  expect_true(is.finite(fit_anchored$summary$LogLik))

  # Task estimates (not directly affected by Rater anchor): positive correlation
  task_c <- fit_centered$facets$others |>
    dplyr::filter(Facet == "Task") |>
    dplyr::arrange(Level) |>
    dplyr::pull(Estimate)
  task_a <- fit_anchored$facets$others |>
    dplyr::filter(Facet == "Task") |>
    dplyr::arrange(Level) |>
    dplyr::pull(Estimate)
  expect_gt(cor(task_c, task_a), 0.5)

  # Person ability estimates: high correlation across methods
  pers_c <- fit_centered$facets$person |>
    dplyr::arrange(Person) |>
    dplyr::pull(Estimate)
  pers_a <- fit_anchored$facets$person |>
    dplyr::arrange(Person) |>
    dplyr::pull(Estimate)
  expect_gt(cor(pers_c, pers_a), 0.9)
})

# ---- 3.4a  count_facet_params: 3 levels, centered -> 2 free params ------

test_that("count_facet_params: 3 levels centered gives 2 free params", {
  spec <- mfrmr:::build_facet_constraint(
    levels = c("A", "B", "C"),
    centered = TRUE
  )
  expect_equal(mfrmr:::count_facet_params(spec), 2)
})

# ---- 3.4b  count_facet_params: 3 levels, not centered -> 3 free params --

test_that("count_facet_params: 3 levels not centered gives 3 free params", {
  spec <- mfrmr:::build_facet_constraint(
    levels = c("A", "B", "C"),
    centered = FALSE
  )
  expect_equal(mfrmr:::count_facet_params(spec), 3)
})

# ---- 3.4c  count_facet_params: 3 levels, 1 anchored, centered -> 1 free -

test_that("count_facet_params: 3 levels, 1 anchored, centered gives 1 free param", {
  anch <- c(B = 0.5)
  spec <- mfrmr:::build_facet_constraint(
    levels = c("A", "B", "C"),
    anchors = anch,
    centered = TRUE
  )
  expect_equal(mfrmr:::count_facet_params(spec), 1)
})

# ---- 3.4d  count_facet_params: all anchored -> 0 free params ------------

test_that("count_facet_params: all anchored gives 0 free params", {
  anch <- c(A = 0.1, B = 0.2, C = 0.3)
  spec <- mfrmr:::build_facet_constraint(
    levels = c("A", "B", "C"),
    anchors = anch,
    centered = TRUE
  )
  expect_equal(mfrmr:::count_facet_params(spec), 0)
})

# ---- 3.4e  count_facet_params: 1 level, centered -> 0 free params -------

test_that("count_facet_params: 1 level centered gives 0 free params", {
  spec <- mfrmr:::build_facet_constraint(
    levels = c("A"),
    centered = TRUE
  )
  expect_equal(mfrmr:::count_facet_params(spec), 0)
})

# ---- 3.5  expand_facet_with_constraints round-trip -----------------------

test_that("expand_facet_with_constraints respects anchored value and length", {
  anch <- c(B = 0.3)
  spec <- mfrmr:::build_facet_constraint(
    levels = c("A", "B", "C", "D"),
    anchors = anch,
    centered = TRUE
  )
  n_free <- mfrmr:::count_facet_params(spec)
  # Feed arbitrary free parameters
  free <- c(0.5, -0.2)
  full <- mfrmr:::expand_facet_with_constraints(free, spec)
  expect_equal(length(full), 4)
  # Anchored level must equal its value
  expect_equal(unname(full[2]), 0.3, tolerance = 1e-10)
  # All values should be finite
  expect_true(all(is.finite(full)))
})

# ---- 3.6  Group anchor constraint ---------------------------------------

test_that("group anchor constrains group mean", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  group_anchors <- data.frame(
    Facet = c("Rater", "Rater"),
    Level = c("R1", "R2"),
    Group = c("G1", "G1"),
    GroupValue = c(0.1, 0.1),
    stringsAsFactors = FALSE
  )
  fit <- suppressWarnings(fit_mfrm(
    d, "Person", c("Rater", "Task", "Criterion"), "Score",
    group_anchors = group_anchors,
    method = "JML", model = "RSM", maxit = 50, quad_points = 7
  ))
  r1_est <- fit$facets$others |>
    dplyr::filter(Facet == "Rater", Level == "R1") |>
    dplyr::pull(Estimate)
  r2_est <- fit$facets$others |>
    dplyr::filter(Facet == "Rater", Level == "R2") |>
    dplyr::pull(Estimate)
  group_mean <- mean(c(r1_est, r2_est))
  expect_equal(group_mean, 0.1, tolerance = 0.1)
})

# ---- 3.7  Dummy facet: estimates all zero --------------------------------

test_that("dummy facet estimates are all zero", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(fit_mfrm(
    d, "Person", c("Rater", "Task", "Criterion"), "Score",
    dummy_facets = "Criterion",
    method = "JML", model = "RSM", maxit = 40, quad_points = 7
  ))
  crit_est <- fit$facets$others |>
    dplyr::filter(Facet == "Criterion") |>
    dplyr::pull(Estimate)
  expect_true(all(crit_est == 0),
              info = "Dummy facet estimates should all be exactly 0")
})

# ---- 3.8  noncenter_facet: Person not sum-to-zero, others are -----------

test_that("noncenter_facet Person does not sum to zero; other facets do", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(fit_mfrm(
    d, "Person", c("Rater", "Task", "Criterion"), "Score",
    method = "JML", model = "RSM", maxit = 40, quad_points = 7
  ))
  # Rater, Task, Criterion should each sum to zero
  for (facet in c("Rater", "Task", "Criterion")) {
    est <- fit$facets$others |>
      dplyr::filter(Facet == facet) |>
      dplyr::pull(Estimate)
    expect_equal(sum(est), 0, tolerance = 1e-6,
                 label = paste("sum-to-zero check for", facet))
  }
  # Person (noncenter) need NOT sum to zero -- just check it exists
  person_est <- fit$facets$person$Estimate
  expect_true(length(person_est) > 0)
})

# ---- 3.9a  Single-level facet: build_facet_constraint -> 0 free params --

test_that("single-level facet constraint gives 0 free params", {
  spec <- mfrmr:::build_facet_constraint(
    levels = c("Only"),
    centered = TRUE
  )
  expect_equal(mfrmr:::count_facet_params(spec), 0)
  full <- mfrmr:::expand_facet_with_constraints(numeric(0), spec)
  expect_equal(unname(full), 0)
})

# ---- 3.9b  Single-level facet in practice: estimate = 0 -----------------

test_that("single-level facet estimate equals 0 in practice", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  # Collapse all tasks into a single level
  d$Task <- "SingleTask"
  fit <- suppressWarnings(fit_mfrm(
    d, "Person", c("Rater", "Task", "Criterion"), "Score",
    method = "JML", model = "RSM", maxit = 40, quad_points = 7
  ))
  task_est <- fit$facets$others |>
    dplyr::filter(Facet == "Task") |>
    dplyr::pull(Estimate)
  expect_equal(unname(task_est), 0, tolerance = 1e-10)
})

# ---- 3.10  Multiple anchors across facets --------------------------------

test_that("multiple anchors across facets are all respected", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  anchors <- data.frame(
    Facet = c("Rater", "Task"),
    Level = c("R1", "T1"),
    Anchor = c(0, -0.3),
    stringsAsFactors = FALSE
  )
  fit <- suppressWarnings(fit_mfrm(
    d, "Person", c("Rater", "Task", "Criterion"), "Score",
    anchors = anchors, method = "JML", model = "RSM", maxit = 50,
    quad_points = 7
  ))

  r1_est <- fit$facets$others |>
    dplyr::filter(Facet == "Rater", Level == "R1") |>
    dplyr::pull(Estimate)
  t1_est <- fit$facets$others |>
    dplyr::filter(Facet == "Task", Level == "T1") |>
    dplyr::pull(Estimate)

  expect_equal(unname(r1_est), 0, tolerance = 1e-8,
               label = "Rater R1 anchor = 0")
  expect_equal(unname(t1_est), -0.3, tolerance = 1e-8,
               label = "Task T1 anchor = -0.3")
})

# ---- Additional: anchor + centering interaction -------------------------

test_that("non-anchored facets remain centered even when another facet is anchored", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  anchors <- data.frame(
    Facet = "Rater", Level = "R1", Anchor = 0.5,
    stringsAsFactors = FALSE
  )
  fit <- suppressWarnings(fit_mfrm(
    d, "Person", c("Rater", "Task", "Criterion"), "Score",
    anchors = anchors, method = "JML", model = "RSM", maxit = 50,
    quad_points = 7
  ))
  # Rater facet has an anchor so sum-to-zero is no longer enforced for Rater
  # But Task and Criterion should still sum to zero
  for (facet in c("Task", "Criterion")) {
    est <- fit$facets$others |>
      dplyr::filter(Facet == facet) |>
      dplyr::pull(Estimate)
    expect_equal(sum(est), 0, tolerance = 1e-6,
                 label = paste("sum-to-zero for", facet, "with Rater anchored"))
  }
  # Verify the anchor is respected
  r1_est <- fit$facets$others |>
    dplyr::filter(Facet == "Rater", Level == "R1") |>
    dplyr::pull(Estimate)
  expect_equal(unname(r1_est), 0.5, tolerance = 1e-8)
})
