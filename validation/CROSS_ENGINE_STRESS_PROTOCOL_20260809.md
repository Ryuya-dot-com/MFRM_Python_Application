# Python–R cross-engine MFRM stress protocol (2026-08-09)

## Status and purpose

This is a deterministic feasibility and failure-mode pilot. It checks whether
the standalone Python application, the in-development `mfrmr` 0.2.3 snapshot,
TAM, immer, and sirt recover the same RSM/PCM response surface when their
estimands can be matched. It is not a high-replication coverage study, a package
ranking, or a sample-size recommendation.

The Python target is local commit `e729f19`. The R target is an isolated install
of the dirty `mfrmr` 0.2.3 development snapshot at source commit `7ee1fd5`;
the working source repository was not modified. Exact loaded-function hashes
and package versions are recorded in `cross_engine_stress_20260809/`.

## Common estimand

For rater \(r\), criterion \(c\), and response category \(k > 0\), every
supported fit is converted to the cumulative-difficulty surface

\[
D_{rck}=k(\rho_r+\delta_c)+\sum_{h=1}^{k}\tau_{ch}.
\]

RSM uses one shared step vector; PCM uses criterion-specific steps. Because
software packages impose different location constraints, the fitted surface is
aligned to truth by the least-squares location shift

\[
\hat a=\arg\min_a\sum_{r,c,k}\{\hat D_{rck}-D_{rck}-ka\}^2.
\]

Bias, RMSE, maximum absolute error, and correlation are then computed on
`EstimateAligned`. JML person comparisons exclude observed all-minimum and
all-maximum persons. MML person comparisons are descriptive EAP recovery, not
recovery of freely estimated person parameters.

## Design cells

Twenty generated datasets cover:

- RSM and PCM;
- crossed balanced panels;
- connected two-rater cycle assignments;
- deliberately disconnected two-component assignments;
- 20% forced all-minimum/all-maximum persons;
- true latent SD 2.2;
- rare boundary categories;
- person-rater local dependence;
- score-dependent MNAR deletion;
- a 28-person small-sample cell.

Balanced, sparse-connected, and high-variance cells use two deterministic
replicates. Other adversarial cells use one. This replication is sufficient for
software and failure-mode comparisons only; it is not sufficient for interval
coverage or rare-failure-rate estimation.

The long generated ratings, truth surfaces, person truths, realized category
support, seeds, and realized graph components are retained with the results.

## Estimator modes

### JML

- Python application: unadjusted JML finite optimizer trace.
- mfrmr 0.2.3: unadjusted JML with typed unbounded-person status.
- TAM: `adj=0, bias=FALSE` and separately `adj=0.3, bias=FALSE`.
- immer: `est_method="jml"`, which still applies epsilon handling to extreme
  persons, and separately `est_method="eps_adj"`.

These modes are never pooled as if they used the same extreme-score or
finite-item-bias convention.

### MML

- Python: free population SD at q=31; q=61 in high-variance/rare-category
  cells; fixed-SD/person-noncentered q=15 and q=31 balanced controls.
- mfrmr: intercept-only latent regression with free residual SD at q=31/q=61;
  fixed-standard-normal q=15/q=31 balanced controls.
- TAM: estimated variance using 21 or 61 equally spaced nodes on `[-8, 8]`.
- sirt `rm.facets`: PCM/rater-facet MML using 30 or 61 equally spaced nodes on
  `[-8, 8]`.

Point counts across Gauss–Hermite and equally spaced rules are sensitivity
labels, not claims that the integration rules are identical.

## Negative-control rules

- A disconnected JML design is expected to be rejected or explicitly labelled
  unidentified. A returned optimizer success is a false-ready signal.
- A disconnected MML design may be numerically estimable through the shared
  population distribution, but it must remain review-only because the observed
  rater graph is disconnected.
- A finite JML trace for an extreme person is not called a finite MLE.
- A TAM fit that changes the pseudoitem category map is not counted as a strict
  common-surface parity result even if it returns a likelihood.
- Optimizer completion, terminal-gradient quality, and inference readiness are
  stored separately.

## Reproduction

From the Python application repository:

```text
python3 validation/cross_engine_stress.py generate \
  --output validation/cross_engine_stress_20260809

python3 validation/cross_engine_stress.py fit-python \
  --input validation/cross_engine_stress_20260809 \
  --output validation/cross_engine_stress_20260809

Rscript validation/cross_engine_stress.R \
  --input validation/cross_engine_stress_20260809 \
  --output validation/cross_engine_stress_20260809 \
  --mfrmr-lib <isolated-mfrmr-0.2.3-library> \
  --mfrmr-git-head <source-commit> \
  --mfrmr-source-state <snapshot-state>

python3 validation/cross_engine_stress.py summarize \
  --input validation/cross_engine_stress_20260809 \
  --output validation/cross_engine_stress_20260809
```

The mfrmr development tree should first be copied to temporary storage and
installed into an isolated temporary R library. Do not compile or install in
the dirty development working tree.

## Scope boundary

No common GPCM comparison is made. The Python/mfrmr bounded owner-slope model,
TAM's separate GPCM routes, and sirt's multiplicative item/rater slope model do
not impose the same response kernel. Existing mfrmr GPCM boundary experiments
therefore remain separate evidence.

