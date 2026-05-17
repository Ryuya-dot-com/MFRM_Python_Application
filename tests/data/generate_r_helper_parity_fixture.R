# Regenerate tests/data/r_helper_parity_output.json from the deterministic
# CSV at tests/data/r_bias_parity_input.csv via mfrmr 0.2.0.
#
# Usage from the repository root:
#   Rscript tests/data/generate_r_helper_parity_fixture.R
#
# The output JSON is the reference for four R parity tests in
# tests/test_r_helper_parity.py covering the FACETS d.f. / ZSTD
# alignment, the rater severity directed network, the rater halo
# correlation network, and the design connectivity diagnostics.
# Run this after either side of the math contract changes (the
# helper code, the mfrmr R package version, or the deterministic
# CSV input) to refresh the reference.

suppressMessages({
  library(mfrmr)
  library(jsonlite)
  library(lme4)
})

input_path <- "tests/data/r_bias_parity_input.csv"
output_path <- "tests/data/r_helper_parity_output.json"

stopifnot(file.exists(input_path))

df <- read.csv(input_path, stringsAsFactors = FALSE)
cat("read", nrow(df), "rows from", input_path, "\n")

# A RSM JMLE fit is the simplest path that produces per-observation
# diagnostics. We use RSM (not GPCM) so the FACETS d.f. comparison is
# independent of slope-estimation noise; both Python and R compute the
# same per-observation Var / FourthCentralMoment from the same
# converged step thresholds.
fit <- fit_mfrm(
  data = df,
  person = "Person",
  facets = c("Rater", "Criterion"),
  score = "Score",
  rating_min = 0L,
  rating_max = 2L,
  model = "RSM",
  method = "JML",
  maxit = 200L,
  reltol = 1e-6
)
cat("fit converged:", fit$summary$Converged[1], "; log_lik:",
    sprintf("%.6f", fit$summary$LogLik[1]), "\n")

dx <- diagnose_mfrm(
  fit, residual_pca = "none",
  diagnostic_mode = "legacy",
  fit_df_method = "both"
)
obs <- as.data.frame(dx$obs)
cat("obs rows:", nrow(obs), "\n")

# ---------------------------------------------------------------------------
# FACETS d.f. / ZSTD alignment
# ---------------------------------------------------------------------------
# diagnose_mfrm() with fit_df_method = "both" populates dx$fit with both
# engine + FACETS d.f. and ZSTD columns; Python's _apply_fit_df_method
# follows the same dispatch.
fit_tbl <- as.data.frame(dx$fit)
facets_alignment_rows <- fit_tbl[, c(
  "Facet", "Level", "N", "Infit", "Outfit",
  "DF_Infit", "DF_Outfit",
  "DF_Infit_FACETS", "DF_Outfit_FACETS",
  "InfitZSTD_ENGINE", "OutfitZSTD_ENGINE",
  "InfitZSTD_FACETS", "OutfitZSTD_FACETS"
)]

# ---------------------------------------------------------------------------
# Rater severity directed network
# ---------------------------------------------------------------------------
# Use severity_direction mode so each pair returns Rater1HigherProp /
# Rater1HigherCount / DirectionN. The Python helper computes the same
# closed form on raw scores.
rn <- rater_network_analysis(
  fit, diagnostics = dx,
  rater_facet = "Rater",
  mode = "severity_direction",
  min_pair_n = 1L,
  score_diff_tolerance = 0
)
rn_pairs <- as.data.frame(rn$pair_metrics)
severity_rows <- rn_pairs[, c(
  "Rater1", "Rater2", "N", "MeanDiff", "MAD",
  "Rater1HigherCount", "Rater2HigherCount",
  "Rater1HigherProp", "Rater2HigherProp", "DirectionN"
)]

# ---------------------------------------------------------------------------
# Rater halo network: Spearman correlations across (Rater, Criterion)
# nodes over shared Person contexts.
# ---------------------------------------------------------------------------
halo <- rater_halo_network_analysis(
  fit, diagnostics = dx,
  rater_facet = "Rater",
  criterion_facet = "Criterion",
  method = "spearman",
  min_pair_n = 5L,
  alpha = 0.05,
  p_adjust = "bonferroni",
  positive_only = TRUE
)
halo_pairs <- as.data.frame(halo$pair_metrics)
halo_rows <- halo_pairs[, c(
  "From", "To", "Rater1", "Criterion1", "Rater2", "Criterion2",
  "EdgeType", "Estimate", "N", "PValue", "PAdjusted"
)]

# ---------------------------------------------------------------------------
# Design network connectivity (undirected co-observation graph)
# ---------------------------------------------------------------------------
dn <- mfrm_network_analysis(fit, diagnostics = dx)
dn_summary <- as.data.frame(dn$summary)
dn_nodes <- as.data.frame(dn$node_metrics)
dn_summary_keep <- intersect(
  c("Nodes", "Edges", "Components", "LargestComponentNodes",
    "LargestComponentShare", "Density", "MeanDegree", "MeanStrength",
    "ArticulationPoints", "Bridges", "Connected"),
  colnames(dn_summary)
)
dn_summary <- dn_summary[, dn_summary_keep]
dn_node_keep <- intersect(
  c("Node", "Facet", "Level", "Degree", "Strength", "IsArticulationPoint"),
  colnames(dn_nodes)
)
dn_nodes <- dn_nodes[, dn_node_keep]

# ---------------------------------------------------------------------------
# G-study variance components via lme4 REML (Brennan 2001 Eq. 3.14)
# ---------------------------------------------------------------------------
# Fit the full 3-way random-effects ANOVA on the same CSV via lme4 with
# explicit two-way interaction terms; this is the canonical reference
# for the Python helper's method-of-moments decomposition. Under the
# balanced one-observation-per-cell design the MoM and REML estimators
# are algebraically equivalent (Searle, Casella, & McCulloch 1992
# Ch. 3.4 / 4.6) so machine-precision agreement is expected.
gstudy_df <- read.csv(input_path, stringsAsFactors = FALSE)
gstudy_df$Person <- as.factor(gstudy_df$Person)
gstudy_df$Rater <- as.factor(gstudy_df$Rater)
gstudy_df$Criterion <- as.factor(gstudy_df$Criterion)

gstudy_formula <- Score ~ 1 +
  (1 | Person) + (1 | Rater) + (1 | Criterion) +
  (1 | Person:Rater) + (1 | Person:Criterion) +
  (1 | Rater:Criterion)
lmer_fit <- lmer(gstudy_formula, data = gstudy_df, REML = TRUE)
vc <- as.data.frame(VarCorr(lmer_fit))
vc <- vc[is.na(vc$var2), c("grp", "vcov")]

# Map lme4's group labels (e.g. "Person:Rater") to the Python "Source"
# labels (same separator); Residual maps to itself.
gstudy_rows <- data.frame(
  Source = as.character(vc$grp),
  Variance = as.numeric(vc$vcov),
  stringsAsFactors = FALSE
)
gstudy_rows <- gstudy_rows[order(gstudy_rows$Source), ]
cat("lme4 G-study variance components:\n")
print(gstudy_rows)

# G / Phi coefficients from the lme4 REML variance components
# (Brennan 2001 Eq. 3.18). The parity test compares these against
# the Python helper's MoM-derived G and Phi; small differences are
# expected when MoM clamping is active at the variance-component
# level, but the coefficients themselves should agree at a
# manuscript-readable precision.
vc_lookup <- setNames(gstudy_rows$Variance, gstudy_rows$Source)
n_p_obs <- length(unique(as.character(gstudy_df$Person)))
n_r_obs <- length(unique(as.character(gstudy_df$Rater)))
n_c_obs <- length(unique(as.character(gstudy_df$Criterion)))
sigma2_p_R <- vc_lookup["Person"]
sigma2_r_R <- vc_lookup["Rater"]
sigma2_c_R <- vc_lookup["Criterion"]
sigma2_pr_R <- vc_lookup["Person:Rater"]
sigma2_pc_R <- vc_lookup["Person:Criterion"]
sigma2_rc_R <- vc_lookup["Rater:Criterion"]
sigma2_residual_R <- vc_lookup["Residual"]
sigma2_delta_R <- (
  sigma2_pr_R / n_r_obs
  + sigma2_pc_R / n_c_obs
  + sigma2_residual_R / (n_r_obs * n_c_obs)
)
sigma2_big_delta_R <- (
  sigma2_delta_R
  + sigma2_r_R / n_r_obs
  + sigma2_c_R / n_c_obs
  + sigma2_rc_R / (n_r_obs * n_c_obs)
)
g_lme4 <- as.numeric(sigma2_p_R / (sigma2_p_R + sigma2_delta_R))
phi_lme4 <- as.numeric(sigma2_p_R / (sigma2_p_R + sigma2_big_delta_R))
cat(sprintf("lme4 REML observed-design G = %.6f, Phi = %.6f\n", g_lme4, phi_lme4))
gstudy_coefficients <- list(
  G = g_lme4,
  Phi = phi_lme4,
  n_p = n_p_obs, n_r = n_r_obs, n_c = n_c_obs,
  decomposition = "lme4_REML_with_two_way_interactions"
)

fixture <- list(
  metadata = list(
    r_version = R.version.string,
    mfrmr_version = as.character(packageVersion("mfrmr")),
    fit_summary = list(
      model = "RSM",
      method = "JML",
      maxit = 200L,
      reltol = 1e-6,
      converged = as.logical(fit$summary$Converged[1]),
      log_lik = as.numeric(fit$summary$LogLik[1])
    )
  ),
  facets_alignment = facets_alignment_rows,
  rater_severity = severity_rows,
  rater_halo = halo_rows,
  design_network_summary = dn_summary,
  design_network_nodes = dn_nodes,
  g_study_variance_components = gstudy_rows,
  g_study_coefficients = gstudy_coefficients
)
write_json(fixture, output_path,
           auto_unbox = TRUE, na = "null", digits = 17, pretty = TRUE)
cat("wrote 5-helper fixture to", output_path, "\n")
cat("  facets_alignment rows :", nrow(facets_alignment_rows), "\n")
cat("  rater_severity pairs  :", nrow(severity_rows), "\n")
cat("  rater_halo pairs      :", nrow(halo_rows), "\n")
cat("  design_network nodes  :", nrow(dn_nodes), "\n")
cat("  g_study var. components rows :", nrow(gstudy_rows), "\n")
cat("R version :", R.version.string, "\n")
cat("mfrmr     :", as.character(packageVersion("mfrmr")), "\n")
