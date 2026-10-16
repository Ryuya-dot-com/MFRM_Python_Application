# Read the quadrature comparison and recheck locally

Start with both the EAP and posterior-SD changes under “Compare quadrature settings.”
**This holds the original estimates and population SD fixed. It does not refit the model.**

## What do you need to check?

| Question | Comparison needed |
|---|---|
| Has estimation progressed sufficiently under the current quadrature rule? | Inspect the stopping reason and gradient; if needed, optimize further under the same rule and examine the change in estimates |
| Do scores change when integration changes at the current estimates? | The same-point comparison on screen; examine EAP and posterior SD separately |
| Would refitting with different node counts change the conclusions? | Refit with the same data, model and population SD, then compare estimates and scores |
| How accurate is the integration approximation? | Evaluate the same estimates against an accuracy-checked reference, such as continuous integration |

The baseline row is zero because it is compared with itself; this does not establish accuracy.
EAP is the posterior mean person score. Posterior SD describes score uncertainty with calibration fixed;
it is not the standard error of population SD or structural parameters. The largest changes in different
columns may belong to different people. NLL changes refer to the summed negative log-likelihood.
A lower NLL under a different quadrature rule does not by itself mean more accurate estimation.

A small gradient describes the calculation under the chosen quadrature rule. It does not establish
that the rule accurately approximates the continuous population distribution. EAP differences can
also be nearly zero while posterior-SD differences remain. After refitting, score changes include
both integration differences and the effect of changed calibration estimates.

If differences could affect your conclusions, save the original results before further comparisons.
Small differences are not an automatic pass either. Define the accuracy needed for the study before
examining comparison results; increasing the node count alone does not complete that check.

## Agreement across programs and integration accuracy

When results differ between programs such as R and Python, match the data, model, constraints,
calibration, population SD, nodes, weights and scoring method before attributing the difference
to the language. Agreement under the same conditions supports reproduction of that calculation;
the accuracy of its approximation needs a separate check. Also distinguish evaluation at the
same estimates from comparison after each program refits the model.

## Save the records and inputs

1. Use “Save check results (JSON)” for the same-point comparison.
2. Open “Report & Export → Review checks → Supporting material → Compare with another program” and save the ZIP.
3. Extract it, keeping `data.csv` and `analysis_snapshot.json` together. Preserve their data, identifiers and
   settings unchanged, privately alongside the original analysis.

The JSON alone cannot run the check below. The ZIP's `comparison_report.html` is an unassessed
worksheet for external comparisons. Downloading it does not run a comparison.

## Python refitting and continuous-integration check (development use)

You need the application source tree containing `validation/mml_fixed_sd_accuracy_check.py` and
the Python dependencies specified by that source. The comparison ZIP does not contain this Python
check engine. Use source and an environment corresponding to the original analysis; retain their
versions and recorded source hashes.

Supported inputs: MML PCM with fixed normal mean zero and fixed population SD, three two-level facets,
scores 0–3, unit weights, and no anchors, regression or penalties. Only the step facet is uncentered;
the other facets are centered. Every facet combination must contain all four observed categories.
Missing person/cell observations are retained; repeated observations in the same person/cell are unsupported.
The fitted GH order must be 31, 61 or 121, SD must be 0.05–10, and responses are capped at 1,000 rows.
These are input limits, not an accuracy guarantee across that scope. For unsupported inputs, use a
model-appropriate validation rather than changing the data or model simply to pass these restrictions.

Replace the paths below with actual locations and run in Terminal. Choose an output folder that **does not yet exist**.

```sh
cd "/path/to/MFRM_Python_Application"
python3 validation/mml_fixed_sd_accuracy_check.py \
  --input "/path/to/extracted/MFRM_Cross_Engine_Validation_Bundle" \
  --output "/path/to/new/pcm_check_01"
```

This refits GH31, GH61 and GH121 from both the original estimates and zero. It evaluates the original
and refitted points by continuous integration, without optimizing the continuous-integral likelihood.
It preserves the original files and saves every trajectory, failures, inputs, source and runtime details.
Before running, the analyst should review the source settings that will be saved in `protocol.json`.
Its numerical targets are development checks, not study-specific equivalence margins.

## Read the output

- `summary.json`: check whether all six planned trajectories were evaluated, their final gradients and
  cross-start differences. Optimizer termination does not establish integration accuracy or qualify SEs/CIs.
- `q*_same_point.json`: evaluations at the original estimates under each rule.
- `q*_original.json` and `q*_zero.json`: refitted results from each start. The `decomposition` separates
  integration differences at the same refitted point from recalibration movement on a continuous-integral
  scale. Total change is relative to **continuous evaluation of the original estimates**, not the original
  finite-GH value on screen. Signed person-level differences add; their maximum absolute values do not.
- If `reference_checks_pass` is false, inspect refinement differences and `allowances`. This development
  script uses conservative error estimates aggregated over people. A failed criterion is not a measurement
  of the true integration error. Investigate the reason, including person-level assessment where needed,
  while preserving the original record. Good results on another development example do not qualify this fit.

For TAM/ConQuest, see `README_cross_engine.md` in the ZIP. Its fixed grid and ConQuest scoring use
different calculations from the GH check above. Importing external results or refitted candidates into
the app as the official analysis is not available.
