# Summarise what happened during harmonisation

Counts, from the unfiltered
[`TwoSampleMR::harmonise_data()`](https://mrcieu.github.io/TwoSampleMR/reference/harmonise_data.html)
output, how many variants were carried forward and why the rest were
not. Used by
[`summary.mr_result()`](https://github.com/BZuckerman97/mrpipeline/reference/summary.mr_result.md)
and
[`summary.coloc_result()`](https://github.com/BZuckerman97/mrpipeline/reference/summary.coloc_result.md).

## Usage

``` r
harmonisation_summary(raw)
```

## Arguments

- raw:

  The `raw` element of
  [`harmonise_and_filter()`](https://github.com/BZuckerman97/mrpipeline/reference/harmonise_and_filter.md)'s
  return value.

## Value

A named list of integers: `n_candidates`, `n_kept`, `n_dropped`,
`n_duplicate`, `n_palindromic`, `n_ambiguous`, `n_incompatible` and
`n_incomplete`. All zero when `raw` is empty or lacks the flag columns.

## Details

A variant can carry more than one flag – an ambiguous variant is by
definition palindromic – so the reason counts are *not* a partition of
`n_dropped` and must not be presented as one. Which flags actually cost
a variant its place depends on the harmonisation action: `remove` at
every level, `ambiguous` from level 2, `palindromic` only at level 3
(see
[`validate_harmonise_action()`](https://github.com/BZuckerman97/mrpipeline/reference/validate_harmonise_action.md)).
`n_incomplete` covers the separate case of a variant dropped by
`harmonise_data()` for missing beta/se rather than for any allele
problem.
