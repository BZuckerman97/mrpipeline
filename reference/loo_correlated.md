# Leave-one-out IVW with correlated instruments

The correlated counterpart of
[`TwoSampleMR::mr_leaveoneout()`](https://mrcieu.github.io/TwoSampleMR/reference/mr_leaveoneout.html):
the random-effects IVW estimate with each instrument dropped in turn,
plus the pooled `"All"` row, fitted with the GLS weight matrix. There is
no upstream implementation –
[`MendelianRandomization::mr_loo()`](https://rdrr.io/pkg/MendelianRandomization/man/mr_loo.html)
takes no `correl` argument and never reads the correlation slot.

## Usage

``` r
loo_correlated(harmonised, ld_matrix, rcond_min = 1e-10)
```

## Arguments

- harmonised:

  Harmonised instrument data frame, in the matrix's row order.

- ld_matrix:

  Signed, aligned LD correlation matrix.

- rcond_min:

  Reciprocal-condition-number threshold below which direct per-SNP
  solves replace the block-inverse update.

## Value

A data frame in the columns of
[`TwoSampleMR::mr_leaveoneout()`](https://mrcieu.github.io/TwoSampleMR/reference/mr_leaveoneout.html)
(`exposure`, `outcome`, `id.exposure`, `id.outcome`, `samplesize`,
`SNP`, `b`, `se`, `p`) plus `ld_corrected`; one row per dropped SNP and
a final `"All"` row.

## Details

Not n refits of `mr_ivw()`, which would be n solves of an (n-1)x(n-1)
system (O(n^4) overall: 15 s at n = 300, minutes at n = 600). The full
inverse `W = O^-1` is computed once and the inverse with SNP `i` removed
is obtained by the Schur-complement identity
`W[-i, -i] - W[-i, i] W[i, -i] / W[i, i]` – O(n^2) per SNP, O(n^3) in
total, exact (agrees with refits to machine precision). With that in
hand the estimate is exactly what
`mr_ivw(correl = TRUE, model = "random")` computes:
`b = (bx' W bx)^-1 bx' W by`, residual SE `sqrt(r' W r / (m - 1))` with
`m` instruments in the fit, standard error
`sqrt(1 / (bx' W bx)) * max(rse, 1)`, normal p-value. `model = "random"`
matches `mr_leaveoneout()`'s default, so the two arms differ only in the
weight matrix.

The update divides by `W[i, i]`, so on a near-singular `O` it inherits
the full inverse's rounding error; below `rcond_min` the function falls
back to a direct solve per SNP instead. (An `O` that bad already makes
the headline `mr_ivw()` fit unstable;
[`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md)
warns about it at the source.)
