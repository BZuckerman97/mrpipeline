# Forest plot of MR methods for one or more results

One row per method (e.g. IVW fixed/random effects, MR Egger, weighted
median), for a single `mr_result` or several sectioned together. Wraps
[`TwoSampleMR::forest_plot_1_to_many()`](https://mrcieu.github.io/TwoSampleMR/reference/forest_plot_1_to_many.html).
Requires the `ggplot2` package.

## Usage

``` r
forest_plot(
  results,
  methods = c("IVW (fixed effects)", "Inverse variance weighted", "MR Egger",
    "Weighted median"),
  relabel = c(`Inverse variance weighted` = "IVW (random effects)"),
  exponentiate = TRUE,
  xlab = NULL,
  trans = NULL,
  xlim = NULL,
  shape_points = 15,
  colour_scheme = "black",
  col_text_size = 5,
  subheading_size = 6,
  col1_width = 1,
  ...
)
```

## Arguments

- results:

  An `mr_result`, or a list of `mr_result` objects. If the list has more
  than one element, every element must be named – names become section
  labels (e.g. `list(Primary = res1, "Positive control" = res2)`). A
  single result (bare, or a length-1 list) is never sectioned,
  regardless of any name given.

- methods:

  Character vector of `method` values to include, in display order.
  Default puts fixed effects above random effects:
  `c("IVW (fixed effects)", "Inverse variance weighted", "MR Egger", "Weighted median")`.

- relabel:

  Named character vector, `c(old = new)`, applied to the `method` column
  for display only, after row filtering/ordering – so matching against
  `results$results$method` is unaffected. Default relabels
  `"Inverse variance weighted"` to `"IVW (random effects)"`.

- exponentiate:

  Logical. Plot on the OR scale. Default `TRUE`.

- xlab:

  X-axis label. Defaults to `"OR (95% CI)"` if `exponentiate`, else
  `"Effect (95% CI)"`.

- trans:

  X-axis scale transform, e.g. `"identity"`, `"log2"`. Defaults to
  `"log2"` if `exponentiate`, else `"identity"`.

- xlim:

  Numeric vector `c(lower, upper)` giving x-axis limits, or `NULL` for
  automatic.

- shape_points:

  Point shape (see
  [`ggplot2::geom_point()`](https://ggplot2.tidyverse.org/reference/geom_point.html)'s
  `shape`).

- colour_scheme:

  Text/point colour. Note: font *family* is not configurable –
  [`TwoSampleMR::forest_plot_1_to_many()`](https://mrcieu.github.io/TwoSampleMR/reference/forest_plot_1_to_many.html)
  has no such argument.

- col_text_size:

  Row-label font size.

- subheading_size:

  Section-heading font size.

- col1_width:

  Width of the row-label column.

- ...:

  Forwarded to
  [`TwoSampleMR::forest_plot_1_to_many()`](https://mrcieu.github.io/TwoSampleMR/reference/forest_plot_1_to_many.html)
  (e.g. `addcols`, `addcol_widths`, `addcol_titles`, `weight`).

## Value

A grid plot object, or `NULL` invisibly if no result has successful,
matching rows.

## Examples

``` r
if (FALSE) { # \dontrun{
result <- run_mr(
  exposure = cd40_exposure, exposure_id = "CD40",
  outcome = sjogren_outcome, outcome_id = "SjD",
  instrument_region = list(chromosome = "20", start = 44746911, end = 44758502),
  bfile = system.file("extdata", "ld_ref", package = "mrpipeline"),
  methods = c("ivw", "ivw_fe", "egger", "weighted_median")
)

# Single result: one row per method
forest_plot(result)

# Several results sectioned together, e.g. primary analysis plus controls
forest_plot(list(
  Primary = result,
  "Positive control" = positive_control_result,
  "Negative control" = negative_control_result
))
} # }
```
