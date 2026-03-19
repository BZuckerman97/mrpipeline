# Resolve sample size from multiple sources

Attempts to determine sample size from (in order of priority):

1.  An explicitly provided value

2.  The median of a data column

3.  `NULL` (caller decides whether to error or warn)

## Usage

``` r
resolve_sample_size(explicit_n = NULL, data_column = NULL, label = "dataset")
```

## Arguments

- explicit_n:

  Explicit sample size (numeric scalar or `NULL`).

- data_column:

  Numeric vector (e.g. `samplesize.exposure` column), or `NULL`.

- label:

  Character label for messages (e.g. `"exposure"`).

## Value

Integer sample size, or `NULL` if unavailable.
