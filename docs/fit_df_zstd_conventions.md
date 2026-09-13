# Fit d.f. and ZSTD conventions

The app uses a FACETS-primary reporting profile for new RSM and PCM analyses:

- Wright-Masters fourth-moment degrees of freedom;
- the Wilson-Hilferty cube-root transformation; and
- an absolute ZSTD cap of 9.

Changing the convention changes d.f., ZSTD, and ZSTD-based flags. It does not
refit the model and does not change measures, standard errors, Infit MnSq, or
Outfit MnSq.

## FACETS-primary formulas

For observation-level conditional variance \(V_i\), fourth central moment
\(C_i\), and analysis weight \(w_i\), the implementation uses

\[
df_{\mathrm{Infit}} =
\frac{2(\sum_i V_iw_i)^2}
     {\sum_i w_i(C_i - V_i^2)}
\]

and

\[
df_{\mathrm{Outfit}} =
\frac{2(\sum_i w_i)^2}
     {\sum_i w_i(C_i/V_i^2 - 1)}.
\]

Primary `DF_Infit`, `DF_Outfit`, `InfitZSTD`, and `OutfitZSTD` columns use
these d.f. values. The corresponding `*_FACETS` columns make the convention
explicit. The previous homogeneous-variance engine convention remains in
`DF_Infit_ENGINE`, `DF_Outfit_ENGINE`, `InfitZSTD_ENGINE`, and
`OutfitZSTD_ENGINE` for sensitivity analyses.

`FitDfMethod`, `FitDfFormula`, `FitZSTDTransform`, `FitZSTDCap`, and
`FitDfApplicability` travel with exported fit tables. Category diagnostics and
the generated Python/R reproduction scripts use the same profile.

## GPCM boundary

Bounded GPCM has free discrimination parameters and is not the strict
Rasch-family model implemented by FACETS. The app therefore labels its
fourth-moment result `facets_style_approximation_for_gpcm`, displays a warning,
and preserves the engine sidecars. Do not report a GPCM ZSTD as exact FACETS
equivalence.

## Reproducing the previous convention

In the Fit Details panel, choose **Engine** to make the historical engine d.f.
and ZSTD primary. Choose **Engine-primary comparison + FACETS sidecars** to
inspect both conventions while keeping engine values in the unsuffixed columns.
The default **FACETS-primary + engine sidecars** option is recommended for RSM
and PCM reporting.

References: Wright and Masters (1982); Wilson and Hilferty (1931); Linacre
(2002).
