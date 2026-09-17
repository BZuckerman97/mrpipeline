# Cochran's Q from the correlated fits

The generalised heterogeneity statistic for correlated instruments,
`Q = r' O^-1 r` on the GLS residuals, as
`mr_ivw(correl = TRUE)@Heter.Stat` (`n - 1` df; identical for fixed and
random effects, since Q depends only on the estimate and the weight
matrix) and `mr_egger(correl = TRUE)@Heter.Stat` (`n - 2` df, needs
`n >= 3`). Same columns and row order as
[`TwoSampleMR::mr_heterogeneity()`](https://mrcieu.github.io/TwoSampleMR/reference/mr_heterogeneity.html)
(Egger first), plus `ld_corrected`; the `method` labels are
TwoSampleMR's so the column is the only thing that differs between the
two arms.

## Usage

``` r
heterogeneity_correlated(ld_input, harmonised)
```

## Arguments

- ld_input:

  An `MRInput` carrying the correlation matrix.

- harmonised:

  Harmonised instrument data frame (for the id columns).

## Value

A data frame with columns `id.exposure`, `id.outcome`, `outcome`,
`exposure`, `method`, `Q`, `Q_df`, `Q_pval`, `ld_corrected`.
