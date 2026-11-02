# Stack a gene track under a regional plot

[`locuszoomr::genetracks_grob()`](https://rdrr.io/pkg/locuszoomr/man/genetracks_grob.html)
fills its panel with an x scale of the locus range plus 4% each side, so
both plots use exactly that axis (see
[`gene_track_x_scale()`](https://github.com/BZuckerman97/mrpipeline/reference/gene_track_x_scale.md)).
With ggplot2's default 5% expansion the genes would sit slightly off the
points.

## Usage

``` r
add_gene_track(p, dat, ens_db, title, subtitle)
```

## Arguments

- p:

  The regional plot, with positions in Mb.

- dat:

  `harmonised_data` of the `coloc_result`.

- ens_db:

  As for
  [`plot_coloc_regional()`](https://github.com/BZuckerman97/mrpipeline/reference/plot_coloc_regional.md).

- title, subtitle:

  Moved from `p` to the composed plot, where a later `labs(title = ...)`
  replaces them.

## Value

A composed ggplot object.
