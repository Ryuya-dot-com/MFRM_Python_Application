# MFRM Estimation (Streamlit)

A Streamlit application for estimating **Many-Facet Rasch Models (MFRM)** without FACETS. The app supports RSM/PCM, JMLE/MML, fit and bias diagnostics, rich visualizations, and APA-style reporting outputs.

## Table of Contents
- Overview
- Key Features
- Requirements & Installation
- Run the App
- Data Requirements
- Workflow (UI Flow)
- Models and Estimation
- Core Equations (Quick Reference)
- Constraints & Anchoring
- Reporting Options
- Diagnostics (Quick Guide)
- Tabs and Outputs
- Downloaded Files
- Performance Tips
- Troubleshooting
- Files

## Overview
- 1 row = 1 observation (one rating instance).
- Estimate measures for **Person** and multiple **facets** (e.g., Rater, Task, Criterion) on a common logit scale.
- Export tables, diagnostic files, and APA-ready text.

## Key Features
- MFRM estimation: **RSM / PCM**, **JMLE / MML**
- Constraints: anchors, group anchors, non-centering, dummy facets, positive-orientation facets
- Fit statistics: Infit/Outfit, ZSTD, WH exact option
- Bias/interaction (FACETS-style re-estimation)
- Category diagnostics, residual PCA (dimensionality), subset connectivity diagnostics
- Visuals: Wright map, pathway map, category curves, step plots
- Exports: CSV/TXT including Facets-style report and APA drafts

## Requirements & Installation
- Python 3.x
- Packages: `streamlit`, `numpy`, `pandas`, `scipy`, `plotly`

### Recommended (virtual environment)
```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install streamlit numpy pandas scipy plotly
```

> Note: On macOS/Homebrew Python, system-wide `pip install` is blocked (PEP 668). Use a venv as above.

## Run the App
```bash
cd MFRM_App/MFRM_R_Estimation_NoFACETS
streamlit run streamlit_app.py
```

## Data Requirements
### Required Columns
- **Person**: unit being rated (e.g., student ID)
- **Score**: ordered rating category (integer)
- **Facet columns (2+)**: e.g., Rater, Task, Criterion

### Optional Column
- **Weight**: positive numeric weight/frequency; rows with weight ≤ 0 are dropped

### Important Notes
- **1 row = 1 observation**.
- Column names must be unique.
- Score must be numeric; it is cast to `int`.
- Missing values are dropped.

### Non-contiguous categories
If your categories are non-contiguous (e.g., 1, 3, 5), enable **Keep original category values (K)**. Otherwise, the app remaps to contiguous integers.

### Minimal Example (CSV)
```csv
Person,Rater,Task,Criterion,Score
P01,R1,T1,C1,3
P01,R2,T1,C1,4
P02,R1,T2,C2,2
```

### With Weight
```csv
Person,Rater,Task,Criterion,Score,Weight
P01,R1,T1,C1,3,2
P01,R2,T1,C1,4,1
```

### Data Input Modes
- **Sample data**: synthetic demo
- **Paste table**: tab-delimited recommended
- **Upload file**: CSV/TSV/TXT

## Workflow (UI Flow)
1. Choose **Data source**
2. Select **Person**, **Score**, and **Facet columns (2+)**
3. Select **Model** (RSM/PCM) and **Estimation** (JMLE/MML)
4. Set constraints / anchors if needed
5. Click **Run estimation**
6. Review results per tab and export files

## Models and Estimation
### Model
- **RSM**: common thresholds across elements
- **PCM**: thresholds vary by a specified **step facet**

### Estimation
- **JMLE**: traditional many-facet estimation
- **MML**: EAP/SD for persons using Gaussian–Hermite quadrature

### Behavior Notes
- PCM requires a **Step facet**.
- Step parameters are centered (sum to zero within a set).

## Core Equations (Quick Reference)
### Linear predictor (MFRM)
For observation *i* (person *p* and facet levels *l_f*):
```
eta_i = theta_p + sum_f (s_f * beta_{f, l_f})
```
- `s_f = +1` for **Positive facets**, otherwise `-1` (default).

### Category probability
Let `k = 0..K-1` and `step_cum[0] = 0`, `step_cum[k] = sum_{j<=k} tau_j`.

**RSM**
```
P(Y_i = k) = exp(k*eta_i - step_cum[k]) / sum_m exp(m*eta_i - step_cum[m])
```

**PCM** (step facet = c)
```
P(Y_i = k | c_i) = exp(k*eta_i - step_cum_c[k]) / sum_m exp(m*eta_i - step_cum_c[m])
```

### Weighted log-likelihood
```
log L = sum_i w_i * log P(Y_i = k_i)
```

### Separation / Strata / Reliability (facet summaries)
Let `mv = Var(Estimate)` (ddof=1), `ev = mean(SE^2)`, `TV = max(mv - ev, 0)`.
```
Separation = sqrt(TV / ev)
Strata     = (4*Separation + 1) / 3
Reliability = TV / mv
RMSE       = sqrt(ev)
```

## Constraints & Anchoring
### Non-centering
- Exactly one facet can be non-centered; others are centered to sum 0.

### Dummy facet
- All elements fixed to 0 (useful for adjustment or interaction-only facets).

### Positive orientation
- By default, higher facet values are *more severe* (decrease scores).
- **Positive facets** reverse the sign (higher = increases scores).

### Anchors
- **Anchor**: fixes specific element values
- **Group anchor**: fixes a group mean

**Anchor CSV**
```csv
Facet,Level,Anchor
Rater,R1,0.2
Rater,R2,-0.1
```

**Group Anchor CSV**
```csv
Facet,Level,Group,GroupValue
Rater,R1,GroupA,0.0
Rater,R2,GroupA,0.0
Rater,R3,GroupB,0.2
```

## Reporting Options
- **Total scores reported**: include extreme elements (Totalscore=Yes)
- **Omit unobserved elements**: hide elements with zero observations
- **Xtreme correction**: fractional adjustment for extreme-only elements
- **Umean / Uscale / Udecimals**: transform logits into user units

## Diagnostics (Quick Guide)
- **Infit/Outfit**: near 1.0 is expected; >1 suggests underfit; <1 suggests overfit
- **ZSTD**: exact-fit z score; can inflate with large samples
- **Reliability / Separation**: how well elements are distinguishable
- **Bias/Interaction**: systematic effects between facet pairs
- **Category diagnostics**: low counts or disordered thresholds signal scale issues

## Tabs and Outputs
- **Data**: preview, facet summary, category distribution
- **Summary**: estimation summary, anchor summary
- **Results**: quick overview + reliability snapshot
- **Person / Facets**: measures & distributions
- **Report**: Facets-style measurement report
- **Fit / Reliability**: fit statistics, chi-square, agreement
- **Bias / Interactions**: FACETS-style bias tables and pairwise reports
- **Steps / Categories**: step estimates and category diagnostics
- **Dimensionality**: residual PCA
- **Subsets**: connectivity diagnostics
- **Visuals**: APA notes, warnings, summaries
- **Downloads**: export tables and text
- **Help**: glossary, reporting guidance, examples

## Downloaded Files
- `mfrm_summary.csv`: log-likelihood, AIC/BIC, N, categories
- `mfrm_measures.csv`: measures, SE, fit, bias, reliability
- `mfrm_fit.csv`: facet-level fit statistics
- `mfrm_steps.csv`: step/threshold estimates
- `mfrm_scorefile.csv`: expected scores, residuals, category probabilities
- `mfrm_residuals.csv`: residuals + standardized residuals
- `mfrm_facets_report_all.csv`: Facets-style report
- Bias/interaction CSVs and fixed-width TXT
- APA method/results and table/figure notes (TXT)

## Performance Tips
- Bias estimation, PCM, MML, and WHEXACT can be slow.
- Start with **RSM + JMLE** for large datasets.
- Collapse sparse categories or simplify facets if convergence is unstable.

## Troubleshooting
- **Need 2+ facets**: Person/Score 제외 2つ以上必須
- **Duplicate columns**: rename columns to be unique
- **Score not numeric**: convert to integers before upload
- **Only extreme scores**: consider Xtreme correction; interpret with caution
- **No convergence**: simplify model, reduce facets, or add anchors

## Files
- App: `streamlit_app.py`
- R version (reference): `app.R`
