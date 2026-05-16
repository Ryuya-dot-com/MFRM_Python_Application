# Advanced Models Guide — v0.2.0-beta

This document describes the seven advanced response models that the
Streamlit app can export as Stan code, and the workflow for running
them locally and bringing the draws back into the app's Posterior
Viewer.

The core estimator (Python-side EM / JMLE) still covers **RSM, PCM,
and GPCM** and runs entirely in the browser. The models below require
Stan and run on the user's own machine.

## Model catalogue

| Key | Family | Response | Needs Q-matrix | Needs testlet map | Reference |
|---|---|---|---|---|---|
| `DINA` | CDM | Binary | **Yes** | No | de la Torre (2009) |
| `HRM` | Rater | Ordinal | No | No | Patz et al. (2002) |
| `TESTLET_RI` | Testlet | Ordinal | No | **Yes** | Bradlow, Wainer & Wang (1999) |
| `TESTLET_BIFACTOR` | Testlet | Ordinal | No | **Yes** | DeMars (2006) |
| `MIXTURE_RASCH` | Mixture | Binary | No | No | Rost (1990) |
| `IRT_2PL_BINARY` | IRT | Binary | No | No | Birnbaum (1968) |
| `PAIRWISE_BTL` | IRT | Binary | No | No | Bradley & Terry (1952) |

## Workflow

```
 ┌────────────────────────────────────────┐
 │ Streamlit app                          │
 │                                        │
 │ 1. Sidebar →                           │
 │    Advanced models (Stan, DL only)     │
 │ 2. Enable + pick family                │
 │ 3. Download the model data template     │
 │ 4. Upload Q-matrix (DINA) /             │
 │    set class count (Mixture)           │
 │ 5. Generate + Download .stan           │
 └────────────┬───────────────────────────┘
              │  file: mfrm_<model>.stan
              ▼
 ┌────────────────────────────────────────┐
 │ User's machine                         │
 │                                        │
 │ cmdstanpy / cmdstanr / rstan:          │
 │   model.compile()                      │
 │   model.sample(data=...)               │
 │                                        │
 │ Output: CmdStan CSVs / InferenceData   │
 └────────────┬───────────────────────────┘
              │  posterior draws
              ▼
 ┌────────────────────────────────────────┐
 │ Streamlit app                          │
 │                                        │
 │ 1. Sidebar → App mode =               │
 │    Posterior Viewer                   │
 │ 2. Upload CmdStan CSV / parquet / .nc  │
 │ 3. Inspect:                            │
 │    summary, HMC diagnostics,            │
 │    trace / ridge / pair / forest       │
 └────────────────────────────────────────┘
```

## Q-matrix format (DINA)

Rows = items, columns = attributes, cells = 0 or 1.

```csv
item,attr_A,attr_B,attr_C
I1,1,0,0
I2,0,1,0
I3,1,1,0
I4,0,0,1
I5,1,0,1
```

The sidebar uploader accepts either a leading ID column (recommended)
or plain 0/1 matrices. The app validates:

- every cell is 0 or 1
- every item (row) uses at least one attribute
- every attribute (column) is used by at least one item

Failing Q-matrices still render a warning but do not block the
download — the generator just uses whatever shape arrived.

## Data templates

The sidebar now provides a CSV template for each Stan-first model. These
templates are deliberately small and readable:

- DINA uses a wide binary response matrix plus a separate Q-matrix.
- HRM uses long ordinal rows: `person`, `item`, `rater`, `score`.
- Testlet RI / Bifactor use long ordinal rows with a required `testlet`
  column.
- Mixture Rasch and 2PL use binary long rows: `person`, `item`, `score`.
- Pairwise BTL uses `object_a`, `object_b`, `wins_a`.

Before passing data to Stan, convert text labels to 1-based integer
arrays that match the generated Stan `data` block.

## Local dependence and mixture scope

For the generated testlet models, **local dependence** means residual
association among items sharing a passage, station, prompt bundle, or
other testlet after the main trait is accounted for. The generated Stan
program handles this through testlet random effects. This is different
from FACETS-mode fixed bias/interaction tables between non-person
facets.

For the generated Mixture Rasch model, **mixture** means unobserved
latent response-pattern classes with class-specific item difficulties.
It does not mean response time, speededness, or inattentive rating. If
response time is part of the study, extend the Stan data block and
likelihood explicitly rather than treating the current mixture template
as a response-time model.

## Running Stan locally

With `cmdstanpy`:

```python
from cmdstanpy import CmdStanModel
import pandas as pd, json

model = CmdStanModel(stan_file="mfrm_dina.stan")
data = {
    "J": n_persons,
    "I": n_items,
    "K": n_attributes,
    "Y": y_matrix.tolist(),
    "Q": q_matrix.values.tolist(),
}
fit = model.sample(data=data, chains=4, iter_sampling=1000, iter_warmup=500)
fit.save_csvfiles(dir="draws/")
```

or `cmdstanr`:

```R
library(cmdstanr)
mod <- cmdstan_model("mfrm_dina.stan")
fit <- mod$sample(data = list(J=..., I=..., K=..., Y=..., Q=...),
                  chains = 4, iter_sampling = 1000, iter_warmup = 500)
fit$save_output_files(dir = "draws/")
```

Then upload the per-chain CSVs into the app's Posterior Viewer mode.

## Priors used in generated Stan

| Model | Core priors |
|---|---|
| DINA | `slip, guess ~ Beta(2, 10)`, `class_probs ~ Dirichlet(1)` |
| HRM | `theta ~ N(0, 1)`, `beta_item ~ N(0, 1)`, `kappa ~ N(0, 2)`, `phi ~ LogNormal(0, sigma_phi)`, `eta ~ N(0, sigma_eta)` |
| Testlet RI | standard + `gamma_testlet ~ N(0, sigma_testlet)`, `sigma_testlet ~ N+(0, 1)` |
| Testlet bifactor | same + `theta_testlet_general ~ N(0, 1)` |
| Mixture Rasch | `class_probs ~ Dirichlet(1)`, `theta ~ N(0, 1)`, `beta_class ~ N(0, 2)` |
| 2PL | `theta ~ N(0, 1)`, `beta_item ~ N(0, 1)`, `alpha_item ~ LogNormal(0, 0.5)` |
| BTL | `ability ~ N(0, 1)` |

All generated programs emit `log_lik` in `generated quantities` for
LOO / WAIC via `arviz.loo`, and (where applicable) `y_rep` for
posterior predictive checks.

## Limitations

- The app does **not** emit per-model runner scripts yet — the Stan
  Code sub-tab in **Report → Exports** provides a generic
  cmdstanpy / cmdstanr runner that you adapt to each model's data
  shape.
- Hyperparameters are intentionally visible in generated Stan code.
  Results can change when priors, class counts, or sampler settings
  change; set these values based on study design and sensitivity checks.
- Mixture Rasch and DINA both enumerate over classes/profiles; for
  large data you may want to factor the likelihood yourself.
- HRM assumes a global set of step thresholds; for rater-specific
  thresholds you'd extend `kappa` to be rater-indexed.

## References

- Birnbaum, A. (1968). *Some latent trait models and their use in inferring an examinee's ability.*
- Bradley, R. A., & Terry, M. E. (1952). Rank analysis of incomplete block designs. *Biometrika*, 39(3/4), 324–345.
- Bradlow, E. T., Wainer, H., & Wang, X. (1999). A Bayesian random effects model for testlets. *Psychometrika*, 64, 153–168.
- de la Torre, J. (2009). DINA model and parameter estimation: A didactic. *Journal of Educational and Behavioral Statistics*, 34(1), 115–130.
- DeMars, C. E. (2006). Application of the bi-factor model to testlet-based tests. *JECM*, 43(2), 145–168.
- Patz, R. J., Junker, B. W., Johnson, M. S., & Mariano, L. T. (2002). The hierarchical rater model. *JEBS*, 27(4), 341–384.
- Rost, J. (1990). Rasch models in latent classes. *Applied Psychological Measurement*, 14(3), 271–282.
