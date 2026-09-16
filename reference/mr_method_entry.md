# Look up one registry row

Look up one registry row

## Usage

``` r
mr_method_entry(shortcut = NULL, label = NULL, engine = NULL)
```

## Arguments

- shortcut, label, engine:

  Exactly one of these identifies the row: `shortcut` for the ten named
  methods, `label = "Wald ratio"` for the automatic single-instrument
  path, `engine = "TwoSampleMR::mr"` for the raw passthrough.

## Value

A one-row data frame from
[`mr_method_registry()`](https://github.com/BZuckerman97/mrpipeline/reference/mr_method_registry.md).
