# Build one `$results` row from a registry entry

Every `$results` row
[`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md)
produces goes through here, so the schema – including the `ld_corrected`
and `model` columns that say which estimator ran – cannot differ between
dispatch branches (which would make the final `rbind` fail).

## Usage

``` r
mr_result_row(
  entry,
  exposure_id,
  outcome_id,
  nsnp,
  b,
  se,
  pval,
  ld_corrected,
  label = entry$label
)
```

## Arguments

- entry:

  One-row data frame from
  [`mr_method_entry()`](https://github.com/BZuckerman97/mrpipeline/reference/mr_method_entry.md).

- exposure_id, outcome_id:

  Labels for the `exposure`/`outcome` columns.

- nsnp, b, se, pval:

  The estimate.

- ld_corrected:

  Logical. Whether the LD matrix was used for this fit.

- label:

  The `method` label. Defaults to the registry label; the raw
  TwoSampleMR passthrough supplies TwoSampleMR's own.

## Value

A one-row data frame.
