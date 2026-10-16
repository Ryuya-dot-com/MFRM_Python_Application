# Anchor User Guidelines (MFRM Streamlit App)

## 1) Purpose
Use anchors when you need stable scale linkage across runs, cohorts, prompts, or administrations.

Typical use cases:
- Reusing previously calibrated rater/task measures.
- Linking two administrations with only partial overlap.
- Stabilizing known reference elements.

## 2) Which option to use
- Anchor table: fix specific element measures to exact values.
- Group anchor table: fix a group's mean measure; its individual levels may differ. A stable group mean must be justified by the linking design.
- Dummy facet: fix all levels of a facet to 0. Differences between these levels are not estimated or adjusted for.
- Non-centered facet: omit sum-to-zero centering for the selected facet. The full model still needs identifying constraints; inspect the identification audit.

## 3) Input formats
### 3.1 Anchor table
Required columns:
- Facet
- Level
- Anchor

Example:
```csv
Facet,Level,Anchor
Rater,R01,0.50
Rater,R02,-0.20
Task,T01,0.30
```

### 3.2 Group anchor table
Required columns:
- Facet
- Level
- Group
- GroupValue

Interpretation:
- For each group, the app enforces:
  - sum(level measures in group) = GroupValue * group size

Example:
```csv
Facet,Level,Group,GroupValue
Rater,R01,Core,0.00
Rater,R02,Core,0.00
Rater,R03,Core,0.00
Rater,R04,New,0.20
Rater,R05,New,0.20
```

## 4) Matching rules (critical)
- Facet and Level must exactly match the values in your analysis data.
- Matching is string-based.
- Leading/trailing spaces can break matching.
- Case mismatches can break matching.

Rows with an unknown facet or level are excluded from the constraints and
recorded in the anchor audit. Missing or non-finite anchor values are also
excluded. The default `warn` policy retains these issues for review; `error`
stops estimation when the audit contains warnings or errors. Always inspect
the audit and applied-anchor counts before interpreting a linked scale.

## 5) Practical workflow
1. Run once without anchors and inspect fit and counts.
2. Pick stable reference elements (reasonable fit, enough observations).
3. Build anchor table or group anchor table.
4. Rerun with anchors.
5. Confirm results:
- The anchor/linking audit: issues and applied-anchor counts by facet.
- Exported `anchor_audit_summary.csv` and the report's Anchor status column.

## 6) Anchor status in report
- A: directly anchored.
- G: group-anchor constrained.
- X: group-anchor constrained and extreme-only element.

## 7) Recommended guardrails
- Start with a small number of high-quality anchors.
- Avoid anchoring many weak or misfitting elements.
- Use a group anchor only when the shared scale and group mean are defensible; grouping does not establish stability or equivalence by itself.
- Recheck fit after applying anchors.
- Keep a versioned record of anchor files used in each run.

## 8) Common mistakes
- Wrong Facet name (for example, `Raters` vs `Rater`).
- Level label mismatch (`R1` vs `R01`).
- Hidden spaces in pasted text.
- Mixing scales from incompatible runs without common elements.

## 9) Notes for estimation mode
- In JMLE, person parameters are explicitly estimated with constraints.
- In MML, person estimates are EAP scores under a population model. Do not treat these scores as independently fixed person parameters; verify the supported constraints and their audit.

## 10) Minimal checklist before publishing results
- Anchor files archived with the analysis run.
- Anchor summary non-zero where expected.
- Report shows intended A/G/X flags.
- Key facet ranges and ordering are stable across reruns.
- Shared content, connectedness, drift, and uncertainty support the intended comparison; counts and A/G/X flags alone do not establish valid linking.
