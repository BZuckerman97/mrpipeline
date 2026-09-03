# Table-style forest plot (forestplot package backend)

One row per outcome, with an inline result/CI/p-value table alongside
the plotted estimate and confidence interval. Requires the `forestplot`
package.

## Usage

``` r
table_forest_plot(dat, null_value, xlab, box_colour = "#35568a")
```

## Arguments

- dat:

  Data frame with columns `exposure`, `outcome`, `n_snps`, `estimate`,
  `lower`, `upper`, `p_value`.

- null_value:

  Reference line value (1 for OR scale, 0 for beta scale).

- xlab:

  X-axis label.

- box_colour:

  Point/line colour.

## Value

A forestplot object (render with
[`plot()`](https://rdrr.io/r/graphics/plot.default.html)/[`print()`](https://rdrr.io/r/base/print.html)).

## Examples

``` r
if (FALSE) { # \dontrun{
dat <- data.frame(
  exposure = "Genetically-proxied NLRP3 inhibition",
  outcome = c("CHD", "Stroke", "Type 2 diabetes"),
  n_snps = c(8, 8, 8),
  estimate = c(0.92, 0.95, 1.03),
  lower = c(0.85, 0.88, 0.94),
  upper = c(0.99, 1.02, 1.13),
  p_value = c(0.02, 0.15, 0.5)
)
table_forest_plot(dat, null_value = 1, xlab = "OR (95% CI)")
} # }
```
