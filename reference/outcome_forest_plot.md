# Forest plot of MR outcomes, optionally coloured/shaped by a grouping column

One row per outcome per requested method. `method` is a plain inclusion
filter, not a priority/fallback selector: request one model to get one
row per outcome, or both to get two rows dodged apart on the same
outcome wherever both were computed. `colour_by`/`shape_by` map any
column in `mr_res` (already computed by the caller) to either aesthetic
independently, so e.g. instrument and model can be shown together
without colliding. Requires the `ggplot2` package.

## Usage

``` r
outcome_forest_plot(
  mr_res,
  xlab,
  method = c("Inverse variance weighted", "IVW (fixed effects)"),
  relabel = c(`Inverse variance weighted` = "Random effects", `IVW (fixed effects)` =
    "Fixed effects"),
  colour_by = NULL,
  shape_by = NULL,
  colour_values = NULL,
  shape_values = NULL,
  section_order = NULL,
  dodge_width = 0.5,
  scale = c("or", "beta"),
  or_col = "or",
  or_lci_col = "or_lci95",
  or_uci_col = "or_uci95",
  b_col = "b",
  se_col = "se",
  group_by = NULL,
  group_order = NULL
)
```

## Arguments

- mr_res:

  Data frame with columns `exposure`, `outcome`, `method`,
  `subcategory`, plus whichever value columns `scale` requires (see
  below) – typically built by the caller by binding together the
  `$results` of several `mr_result` objects and attaching a
  `subcategory` (and any grouping columns) with
  [`dplyr::mutate()`](https://dplyr.tidyverse.org/reference/mutate.html).
  See examples.

- xlab:

  X-axis label.

- method:

  Character vector of `method` values to include. Default
  `c("Inverse variance weighted", "IVW (fixed effects)")` includes both
  wherever computed. Outcomes lacking a requested method simply
  contribute no row for it – there is no silent substitution.

- relabel:

  Named character vector, `c(old = new)`, applied to the `method` column
  for display, so it reads sensibly if `colour_by` or `shape_by` is
  `"method"`. Default relabels the two IVW variants to
  `"Random effects"` / `"Fixed effects"`.

- colour_by, shape_by:

  Column name in `mr_res` to map to colour/shape respectively, or `NULL`
  for none.

- colour_values, shape_values:

  Named vectors passed to
  [`ggplot2::scale_colour_manual()`](https://ggplot2.tidyverse.org/reference/scale_manual.html)/[`ggplot2::scale_shape_manual()`](https://ggplot2.tidyverse.org/reference/scale_manual.html),
  or `NULL` to use ggplot2's defaults.

- section_order:

  Character vector giving the `subcategory` display order, or `NULL` to
  use first-appearance order.

- dodge_width:

  Vertical dodge width applied when two rows share an outcome (e.g.
  fixed and random effects for the same outcome).

- scale:

  `"or"` (default) plots an odds-ratio scale: `or_col`/
  `or_lci_col`/`or_uci_col` on a log2 x-axis with a reference line at 1.
  `"beta"` plots a linear effect-size scale: `b_col`/`se_col` (CI
  derived internally as `b \eqn{\pm} 1.96*se`) on a linear x-axis with a
  reference line at 0. Required columns depend on this choice – only the
  3 OR columns are required for `scale="or"`, only the 2 beta/SE columns
  for `scale="beta"`.

- or_col, or_lci_col, or_uci_col:

  Column names to use for the OR point estimate and 95% CI when
  `scale="or"`. Defaults (`"or"`/`"or_lci95"`/ `"or_uci95"`) match every
  existing call site – override to plot a rescaled column (e.g.
  `or_col="or_scaled"`) without needing to rename columns in the caller.

- b_col, se_col:

  Column names to use for the point estimate and standard error when
  `scale="beta"`. Defaults `"b"`/`"se"` – override to plot a rescaled
  column (e.g. `b_col="b_scaled", se_col="se_scaled"`).

- group_by:

  Optional column name in `mr_res` giving an outer grouping above
  `subcategory` (e.g. "Primary" vs "Sensitivity Analysis"). When
  supplied, facets on `vars(.data[[group_by]], subcategory)` instead of
  just `vars(subcategory)` –
  [`ggplot2::facet_grid()`](https://ggplot2.tidyverse.org/reference/facet_grid.html)'s
  native handling of multiple row-facetting variables spans the outer
  variable's strip across its contiguous inner rows, giving a nested
  header with no extra dependency. `NULL` (default) reproduces today's
  single-facet behaviour exactly.

- group_order:

  Character vector giving `group_by`'s display order (parallel to
  `section_order` for `subcategory`), or `NULL` for first-appearance
  order. Ignored if `group_by` is `NULL`.

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{
primary <- run_mr(
  exposure = cd40_exposure, exposure_id = "CD40",
  outcome = sjogren_outcome, outcome_id = "SjD",
  instrument_region = list(chromosome = "20", start = 44746911, end = 44758502),
  bfile = system.file("extdata", "ld_ref", package = "mrpipeline"),
  methods = c("ivw", "ivw_fe")
)

mr_res <- dplyr::bind_rows(
  primary$results |> dplyr::mutate(subcategory = "Primary"),
  positive_control_result$results |> dplyr::mutate(subcategory = "Positive control"),
  negative_control_result$results |> dplyr::mutate(subcategory = "Negative control")
) |>
  dplyr::mutate(
    instrument = dplyr::if_else(exposure == "CD40", "GWS instrument", "Functional instrument")
  )

# Both models, coloured by instrument, shaped by model
outcome_forest_plot(
  mr_res,
  xlab = "OR (95% CI)",
  colour_by = "instrument",
  shape_by = "method",
  colour_values = c("GWS instrument" = "#3c6e5e", "Functional instrument" = "#b5541f"),
  shape_values = c("Random effects" = 16, "Fixed effects" = 17),
  section_order = c("Primary", "Positive control", "Negative control")
)

# Linear beta/SE scale on a rescaled column, nested Primary/Sensitivity header
outcome_forest_plot(
  mr_res,
  xlab = "Beta (95% CI)",
  scale = "beta",
  b_col = "b_scaled", se_col = "se_scaled",
  group_by = "tier_group", group_order = c("Primary", "Sensitivity Analysis"),
  section_order = c("Cohort A", "Cohort B")
)
} # }
```
