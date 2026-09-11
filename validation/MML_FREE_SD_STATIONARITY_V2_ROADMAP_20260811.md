# Free-SD MML stationarity v2 roadmap

Status: **engineering development; not scientifically qualified**.

## Why v2 exists

The frozen `kac200` audit showed that the legacy free-SD EM stopping rule and
finite-Q31 stationarity are different claims. The legacy rule stops on relative
marginal-likelihood movement. It does not gate the joint structural and
`log(sigma)` score. The existing study remains frozen and valid under its
registered v1 gates; v2 is a new estimator implementation and cannot replace or
retroactively upgrade v1.

FACETS 4.5 remains the JMLE calibration comparator. FACETS does not validate
the MML likelihood. Free-SD MML requires a separate independent-R likelihood,
gradient, and Q31/Q61 evidence layer. FACETS display-rounded fit statistics are
never used to calculate raw parity metrics.

## Current engineering state

- Estimator label reserved for future work:
  `PYTHON_MML_FREE_SD_Q31_STATIONARITY_V2`.
- Legacy EM is retained only as a warm start.
- The finite-Q objective is polished jointly in structural coordinates and
  `log(sigma)`, then restarted once from the polished solution.
- Optimizer termination, finite solution, objective non-worsening, score,
  finite-difference agreement, information curvature, Newton correction,
  constraint residual, and restart movement are recorded separately.
- Phase 1 is restricted to unregularized RSM/PCM with unbounded structural
  coordinates and an interior sigma. GPCM and penalized models fail closed.
- A stationarity assessment can report a numerical gate, but it always reports
  `InferenceReady=False`. No bare Q31/Q61 boolean can promote the result.
- A typed Q-fit/Q-check layer now recomputes both stationarity assessments,
  evaluates each solution under both objectives, and retains structural,
  `log(sigma)`, optimization-gain, and back-evaluation differences. It reports
  only `NumericalSensitivityPass`; `ScientificInferenceReady` remains false.
- The Q-independent problem digest binds exact index arrays and shapes,
  observations, effective facet signs, category support, identification
  constraints, parameter and sigma bounds, population design, implementation
  source, and package runtime. Unmodelled nonfinite values fail closed. The
  downstream quadrature adapter is a separate module, so adding it did not
  invalidate the frozen stationarity-v2 source identity.
- An exact-denominator batch layer now binds the canonical `RunId -> problem
  digest` map, unique problem evidence, one common optimizer configuration,
  and one common stationarity and Q sensitivity contract. It reconstructs both
  stationarity assessments from raw runs, replays all four Q cross-objective
  evaluations through one implementation-hash-bound evaluator, and then
  reconstructs every sensitivity summary.
  Missing, unexpected, duplicate, failed, replaced, mixed-contract, or
  relabelled records cannot yield `NumericalBatchPass`. Even a complete
  numerical batch still reports `ScientificInferenceReady=False` until
  prospective registration.
- The current batch API accepts only `DEVELOPMENT_ONLY`. Objective replay is
  resolved outside each record by the auditor. A future registered scope must
  provide a strict registrar-owned resolver that rebuilds each problem from
  frozen raw input and verifies the exact source/runtime manifest; an arbitrary
  caller-supplied resolver is not registration evidence.
- The known KAC case and all current synthetic/unit fixtures are registered as
  observed development evidence and excluded from qualification.

## Three evidence layers

1. **Per-dataset Q31 numerical evidence**

   Every planned dataset retains the warm start, both polish results, full
   joint score at two finite-difference steps, exact constraint residual,
   objective/value-gradient consistency, and failure reason. Optimizer success
   alone is never a readiness criterion.

2. **Study-level numerical qualification**

   A separately registered bundle must bind the exact estimator source,
   category mapping, PCM threshold parameterization, normal population model,
   Q31 primary basis, Q61 sensitivity basis, independent-R evaluator, runtime,
   input hashes, and planned denominator. Q61 is a sensitivity basis, not a raw
   likelihood parity target. A subset qualification cannot be copied into 600
   per-dataset boolean fields.

3. **Scientific inference**

   Scientific endpoints may run only on fresh data generated after layers 1
   and 2, their thresholds, and their failure policy are hash-frozen. The v1
   and v2 endpoint estimates are not pooled. FACETS operational qualification
   remains a separate status from MML scientific inference.

The future study-level registrar, not the numerical kernel, is the only layer
allowed to create `ScientificInferenceReady`.

## Prospective registrar schema kernel

A development-only registrar schema now separates four identities:

1. `T0`: the ordered RunId/seed/DGP/cohort plan, numerical contracts, source
   and runtime manifests, no-replacement policy, independent-R plan, and
   FACETS operational plan, sealed before qualification-data generation;
2. `T1`: the realized input instance for every planned RunId, including
   preregistered input-invalid and generator-failure cases that remain in the
   denominator;
3. `T1.5`: a single-attempt fit authorization that binds each generated input
   instance to a likelihood-problem digest; and
4. `T2`: the complete attempt ledger, including typed no-fit and execution
   failures, positive-cohort outcomes, an explicit versioned negative-control
   gate registry, and separate Q31, Q61, independent-R, and FACETS evidence.

The schema distinguishes a dataset-instance digest from a likelihood-content
digest. Independently generated datasets may coincide in likelihood content
without becoming the same registered case. Positive qualification cases and
intentional negative controls are frozen as different cohorts; a failed
positive case cannot be relabelled as a successful negative control.

FACETS evidence has its own operational status and cannot change the MML
numerical status. Conversely, a FACETS failure remains visible and cannot be
silently removed from its registered comparison subset. FACETS display values
remain prohibited as raw calculation inputs.

This is only a schema and publication kernel. All four stages remain
development-only, all scientific-readiness flags are forced false, and the
`T1.5` status explicitly says that the pure pre-fit problem factory is not yet
qualified. Two blockers therefore remain before any real prospective
registration:

- construct the Q-independent likelihood problem from retained raw input and
  frozen configuration without reading an optimizer result, estimated sigma,
  summary table, or other post-fit state; and
- add strict artifact transport and registrar-owned reconstruction from the
  retained bytes, live source/runtime manifest, and batch record, rather than
  accepting hashes or evaluator callables supplied by the worker.

The current fault tests cover denominator preservation, RunId/seed and cohort
relabeling, invalid-input authorization, numerical retry, unregistered and
wrong-gate negative controls, missing independent-R evidence, FACETS/MML status
separation, duplicate/nonfinite/overflowing JSON, bool-as-integer input,
forged publication parents and publisher identity, validation-time file
mutation, identity-last publication, and cleanup after a simulated second-file
replace failure. They are development fixtures and are excluded from any
future qualification pool.

An external receipt or append-only timestamp is still required before the word
"prospective" can describe an actual run. A self-consistent local hash chain
detects later alteration; it cannot prove that its T0 plan existed before a
researcher saw generated inputs or fitted results. Likewise, the current T1
hash fields are not evidence until a registrar reads retained raw bytes,
verifies their exact file identity, and constructs the fit-free problem from
that same byte snapshot. These are promotion blockers, not optional hardening.

## Prospective sequence

### Development

- Complete negative controls for optimizer failure, nonfinite callbacks,
  inconsistent objectives, indefinite/asymmetric information, fractional Q,
  sigma evidence mismatch, and artifact tampering.
- Add a small PCM development fixture; it is automatically excluded from
  qualification once observed.
- Profile fast first-order diagnostics separately from full Hessian diagnostics.
  The full N=960, Q31 known replay took about 195 seconds in the current serial
  implementation, so this path belongs in batch validation, not synchronous UI.
- Implement a typed, artifact-bound Q31/Q61 study assessment. A boolean is not
  evidence.
- Implement append-only publication with identity last, read-back validation,
  fixed expected input hashes, and failure cleanup.

### Freeze before qualification data exist

- Freeze estimator source and runtime/package identities.
- Freeze optimizer options and all finite-difference steps.
- Freeze objective scale, coordinate order, constraints, sigma bounds, and
  category/threshold map.
- Freeze numerical thresholds and their rationale without using KAC endpoint
  values or current development fixtures.
- Freeze whether the full Hessian is required for every qualification dataset.
- Freeze Q31/Q61 parameter, sigma, score, and denominator rules.
- Freeze independent-R cross-evaluation and gradient rules.
- Freeze no replacement, no optional extension, and all-planned-cases policy.

### Fresh qualification

- Generate fresh PCM qualification data only after the freeze.
- Include every registered gamma condition and fresh Person-vector triplets.
- Include separate sparse, boundary, poor-start, non-PD, and early-stop negative
  controls.
- Do not calculate KA1, effect signs, p-values, or scientific endpoints.
- Require all planned cases. Exceptions remain in the denominator as NOT_READY.
- A failed protocol may be revised only as a new version with new fresh data.

### Fresh confirmatory study

- Generate a new 200-triplet/600-dataset study after qualification passes.
- Keep Q31 per-dataset stationarity, study-level Q31/Q61 qualification,
  independent R, FACETS operational evidence, and scientific endpoints as
  separately reported statuses.
- Preserve the original v1 `kac200` result and its strict audit without edits.

## Current claim boundary

The code now demonstrates that the observed legacy plateau can be polished to
a much smaller joint score on one known PCM case, with an answer matching the
retained independent-R Q31 solution. This is regression evidence only. It is
not a threshold justification, prospective qualification, endpoint
re-analysis, or a reason to relabel the frozen v1 result.
