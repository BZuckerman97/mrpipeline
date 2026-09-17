# mrpipeline – Package Developer Guide for Claude Code

## Package Overview

`mrpipeline` is an R package for Mendelian randomisation (MR) and
colocalization analysis, currently focused on proteomic GWAS data.
However, the plan is to integrate eQTL, scQTL and other GWAS data.

**Core exported functions:**

- [`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md)
  – Cis-MR, genome-wide MR, or manual-instrument MR with sensitivity
  methods.
- [`run_coloc()`](https://github.com/BZuckerman97/mrpipeline/reference/run_coloc.md)
  – Colocalization (coloc.abf, SuSiE, coloc.signals, colocPropTest)
- [`format_pqtl_decode()`](https://github.com/BZuckerman97/mrpipeline/reference/format_pqtl_decode.md)
  – Format deCODE proteomics GWAS to TwoSampleMR exposure format
- [`format_pqtl_ukbppp()`](https://github.com/BZuckerman97/mrpipeline/reference/format_pqtl_ukbppp.md)
  – Format UKB-PPP pQTL data to TwoSampleMR exposure format
- [`format_single_cell_onek1k()`](https://github.com/BZuckerman97/mrpipeline/reference/format_single_cell_onek1k.md)
  – Format OneK1K single-cell eQTL data
- [`last_allele_check()`](https://github.com/BZuckerman97/mrpipeline/reference/last_allele_check.md)
  – Record of the most recent harmonisation allele orientation check
  (detects A1/A2 = REF/ALT files, issue \#18)
- [`mr_methods()`](https://github.com/BZuckerman97/mrpipeline/reference/mr_methods.md)
  – Table of every method
  [`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md)
  accepts: estimator (fixed/random), LD-correctability, minimum
  instruments, output field. Rendered from the internal registry that
  drives dispatch, so it cannot drift (issue \#27)

**Internal helpers** (not exported, in `R/helpers.R`):
[`harmonise_and_filter()`](https://github.com/BZuckerman97/mrpipeline/reference/harmonise_and_filter.md),
[`check_allele_orientation()`](https://github.com/BZuckerman97/mrpipeline/reference/check_allele_orientation.md),
[`check_allele_orientation_gwas()`](https://github.com/BZuckerman97/mrpipeline/reference/check_allele_orientation_gwas.md),
[`compute_ld_matrix()`](https://github.com/BZuckerman97/mrpipeline/reference/compute_ld_matrix.md),
[`clump_instruments()`](https://github.com/BZuckerman97/mrpipeline/reference/clump_instruments.md),
[`align_to_ld_matrix()`](https://github.com/BZuckerman97/mrpipeline/reference/align_to_ld_matrix.md),
[`eaf_to_maf()`](https://github.com/BZuckerman97/mrpipeline/reference/eaf_to_maf.md),
[`resolve_sample_size()`](https://github.com/BZuckerman97/mrpipeline/reference/resolve_sample_size.md),
[`warn_no_ld_correction()`](https://github.com/BZuckerman97/mrpipeline/reference/warn_no_ld_correction.md),
[`gls_weight_matrix()`](https://github.com/BZuckerman97/mrpipeline/reference/gls_weight_matrix.md),
[`diag_ids()`](https://github.com/BZuckerman97/mrpipeline/reference/diag_ids.md),
[`heterogeneity_correlated()`](https://github.com/BZuckerman97/mrpipeline/reference/heterogeneity_correlated.md),
[`pleiotropy_correlated()`](https://github.com/BZuckerman97/mrpipeline/reference/pleiotropy_correlated.md),
[`loo_correlated()`](https://github.com/BZuckerman97/mrpipeline/reference/loo_correlated.md)
(the last three compute the diagnostics from the correlated fits under
`ld_correct = TRUE`;
[`loo_correlated()`](https://github.com/BZuckerman97/mrpipeline/reference/loo_correlated.md)
is a block-inverse GLS refit, O(n^3) – do not replace it with a per-SNP
`mr_ivw()` loop, issue \#31)

**Method registry** (in `R/mr_methods.R`):
[`mr_method_registry()`](https://github.com/BZuckerman97/mrpipeline/reference/mr_method_registry.md)
is the single source of truth for what
[`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md)
can run – shortcut, label, fixed/random model, `ld_correctable`,
`min_instruments`, output field, engine per LD path.
[`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md)
validates, dispatches, computes skip reasons and emits the “no
LD-corrected form” warnings from it;
[`mr_methods()`](https://github.com/BZuckerman97/mrpipeline/reference/mr_methods.md)
exports it;
[`rd_method_table()`](https://github.com/BZuckerman97/mrpipeline/reference/rd_method_table.md)
renders it into
[`?run_mr`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md)
via roxygen `@eval`;
[`mr_method_entry()`](https://github.com/BZuckerman97/mrpipeline/reference/mr_method_entry.md)
/
[`mr_result_row()`](https://github.com/BZuckerman97/mrpipeline/reference/mr_result_row.md)
build every `$results` row from it. **Adding a method means adding a
registry row** – never a new hard-coded name list.
[`empty_mr_results()`](https://github.com/BZuckerman97/mrpipeline/reference/empty_mr_results.md)
(in `R/mr_result.R`) is the zero-row `$results` schema.

**S3 classes:** `mr_result` (from
[`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md)),
`coloc_result` (from
[`run_coloc()`](https://github.com/BZuckerman97/mrpipeline/reference/run_coloc.md))

## `run_mr()` has no `flip_beta` argument – do not reintroduce one

[`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md)
previously had a `flip_beta` argument that negated `beta.exposure` on
the harmonised data before any method/sensitivity computation, so a
downstream project could redefine “increase in exposure” (e.g. modelling
a drug’s inhibition mechanism rather than a GWAS trait’s raw increasing
direction). It has been **removed** (2026-08-01) and must not be
reintroduced in any form – as an argument to
[`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md),
a wrapper around it, or an option threaded through from a caller.

The design was judged unsafe: it only ever touched `beta.exposure`, and
nothing about the function signature made that scope visible at the call
site, which invites misuse by anyone who reaches for `flip_beta`
expecting it to behave like a general allele/effect-direction
reorientation (it never flipped `effect_allele.exposure`,
`other_allele.exposure`, or `eaf.exposure`, nor `beta.outcome`). A
shared library function is the wrong place for a transformation whose
correctness depends entirely on the caller understanding exactly which
columns it does and doesn’t touch.

**If a project needs to redefine exposure direction, write a small
`flipped_beta()`-style helper in that project’s own script** (not in
`mrpipeline`) that makes explicit, at the call site, exactly what is
being transformed – typically negating `beta.exposure` on the
harmonised/exposure data before it reaches
[`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md),
with a comment stating which columns are deliberately left untouched and
why (`effect_allele`/`other_allele`/`eaf` should almost never change for
a phenotype-direction flip; see the distinction between a phenotype
transformation and a genuine allele-reorientation, which also requires
negating `beta` and replacing `eaf` with `1 - eaf`). Keeping this in
project-level code, next to the specific analysis that needs it, keeps
the transformation visible and auditable instead of hidden behind a
generic library flag.

This must be documented consistently across every project’s `CLAUDE.md`
in `MR_projects/` – see the root `CLAUDE.md`’s cross-project recipe for
the canonical wording, and each project’s own `CLAUDE.md` for how (or
whether) that project is affected.

## Build and Check Commands

``` r

devtools::load_all()      # load package for interactive development
devtools::document()      # regenerate NAMESPACE and man/ from roxygen2
devtools::check()         # R CMD Check (must pass: 0 errors, 0 warnings)
devtools::test()          # run testthat tests
pkgdown::build_site()     # build documentation site (must succeed)
```

``` bash
air format .              # auto-format R code (run before every commit)
```

``` r

lintr::lint_package()     # style/lint checks (run before every commit)
```

## Code Conventions

**ASCII only:** Never use non-ASCII characters anywhere in `.R` or `.Rd`
files – R CMD check warns on them and CI will fail. Common offenders: -
Em dash (Unicode U+2014) – write `--` instead - Box-drawing horizontal
(Unicode U+2500, used in RStudio section headers) – write `-` instead -
Right arrow (Unicode U+2192) – write `->` instead - Any other Unicode
punctuation or symbols

Section dividers in code comments must use plain hyphens:
`# -- Section name --`.

**Messages / warnings / errors:** Use `cli` package exclusively.

``` r

cli::cli_inform("Loading {protein} data...")   # informational
cli::cli_warn("Only {n} SNPs retained.")       # warning
cli::cli_abort("bfile is required for LD-corrected MR.")  # error
```

**Pipe:** Native `|>` only. Never `%>%`.

**String operations:** Use `stringr` functions.

``` r

stringr::str_remove(x, "_.*")      # not gsub("_.*", "", x)
stringr::str_detect(x, "^chr")     # not grepl("^chr", x)
stringr::str_starts(x, "rs")       # not startsWith(x, "rs")
```

**Namespace:** Prefer `pkg::fn()` over `@importFrom`.

**Internal helpers:** Tag with `@keywords internal`; do NOT use
`@export`.

**Documentation:** Use roxygen2 for all functions. Vignettes
(`mrpipeline-user-guide.Rmd`, `mrpipeline-developer-guide.Rmd`) are
**living documents** – update them in the same commit as any API
change: - When adding/changing/removing a function parameter: update
roxygen, user guide (usage examples), developer guide
(architecture/internals) - When adding/removing internal helpers: update
the developer guide’s helpers section - When modifying S3 class
structure: update the developer guide’s S3 classes section

**`verbose` argument:** Planned future feature. When added, gate all
[`cli::cli_inform()`](https://cli.r-lib.org/reference/cli_abort.html)
calls behind it.

## Documentation Rule – ALWAYS regenerate man/ after changing function signatures

Any time you add, remove, or rename a parameter in ANY function
(exported or internal), you MUST regenerate the man/ files before
committing. Failure to do this causes R CMD check to error on all CI
platforms.

**After every function signature change, remind the user to run:**

``` bash
air format .           # must be run before devtools::document()
```

``` r

devtools::document()   # must be run in an R session, not via Rscript
```

Then stage and commit the updated `man/*.Rd` files in the same commit as
the code change. Never let a signature change and its .Rd update go in
separate commits.

If `devtools::document()` fails (e.g. missing dependency), patch the
`.Rd` file manually: update the `\usage{}` block and add a `\item{}` in
`\arguments{}`.

## PR Checklist

Before opening any pull request:

1.  `air format .` – auto-format R files (run before
    devtools::document())
2.  `devtools::document()` – regenerate man/ files (run in R, not
    terminal)
3.  `lintr::lint_package()` – fix any lint warnings
4.  [`pkgdown::build_site()`](https://pkgdown.r-lib.org/reference/build_site.html)
    – confirm site builds without errors
5.  `devtools::check()` – must produce 0 errors, 0 warnings
6.  `devtools::test()` – all tests must pass

## Workflow

Feature branches are created off `main` and merged back to `main` via PR
(the former `dev` branch was retired once it fell behind `main`). Each
feature should have a corresponding GitHub issue.

Branch naming: `phase-N/short-description`
(e.g. `phase-1/shared-helpers`).

## Session Logging (Obsidian)

At the end of each session, when the user requests it, write a session
log for their Obsidian vault. Output it as markdown text in the chat for
the user to copy – do not attempt to write to a file path.

**Format – use these four headings every time:**

    ## YYYY-MM-DD -- mrpipeline session

    ### What we did
    - Bullet list of concrete actions taken (code written, targets run, results reviewed)

    ### Why this matters
    - How the session's work connects to the broader analysis goals

    ### What we noticed / questions
    - Interesting findings, unexpected results, open questions, things to investigate

    ### Next steps
    - Specific actionable tasks for the next session, in priority order

## Architecture Pointers

- Full architectural context: `vignettes/mrpipeline-developer-guide.Rmd`
- End-to-end usage examples: `vignettes/mrpipeline-user-guide.Rmd`
- Test data: `R/data.R` (lazy-loaded datasets), `inst/extdata/` (LD
  reference panel)
- Tests requiring plink/LD reference: use
  `testthat::skip_if_not(file.exists(bfile))` for CI safety
