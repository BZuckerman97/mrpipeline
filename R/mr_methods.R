#' The MR method registry
#'
#' The single source of truth for every method [run_mr()] can run: its
#' shortcut name, the label written to `$results$method`, whether it is a
#' fixed- or random-effects estimator, whether `ld_correct = TRUE` applies to
#' it, the minimum number of instruments it needs, where its output lands on
#' the `mr_result`, and the function that actually runs on each LD path.
#' `run_mr()` reads validation, dispatch, skip reasons and the LD warnings
#' from this table, and the method tables in `?run_mr` and the vignettes are
#' rendered from it, so the documentation cannot drift from the code. Adding
#' a method means adding a row here.
#'
#' Two rows have no `shortcut`: the Wald ratio, which `run_mr()` uses
#' automatically when exactly one instrument survives, and the raw
#' TwoSampleMR passthrough, one row standing for every
#' `TwoSampleMR::mr_method_list()$obj` name that has no shortcut.
#'
#' `ld_correctable` is declared, not inferred: nothing inspects upstream
#' function signatures at run time. The row is kept honest by
#' `engine_ld` -- a method claims LD support only by naming the function
#' that provides it, and the tests assert `ld_correctable == !is.na(engine_ld)`.
#' There is deliberately no free-text notes column: notes about upstream
#' internals go stale silently when a dependency changes, whereas
#' `engine`/`engine_ld` name what mrpipeline calls and so stay true as long
#' as the dispatch does.
#'
#' @return A data frame with one row per method path and columns `shortcut`,
#'   `description`, `label`, `output`, `model`, `ld_correctable`,
#'   `min_instruments`, `engine`, `engine_ld`.
#'
#' @keywords internal
mr_method_registry <- function() {
  data.frame(
    shortcut = c(
      "ivw_random",
      "ivw_fixed",
      "egger",
      "weighted_median",
      "presso",
      "conmix",
      NA_character_,
      NA_character_,
      "steiger",
      "pleiotropy",
      "heterogeneity",
      "loo"
    ),
    description = c(
      "IVW, multiplicative random effects",
      "IVW, fixed effects",
      "MR Egger regression",
      "Weighted median",
      "MR-PRESSO outlier test",
      "Contamination mixture",
      "Wald ratio, single instrument",
      "Any other TwoSampleMR method",
      "Steiger directionality test",
      "Egger intercept (pleiotropy) test",
      "Cochran's Q heterogeneity test",
      "Leave-one-out IVW"
    ),
    label = c(
      "IVW (random effects)",
      "IVW (fixed effects)",
      "MR Egger",
      "Weighted median",
      "MR-PRESSO",
      "ConMix",
      "Wald ratio",
      NA_character_,
      NA_character_,
      NA_character_,
      NA_character_,
      NA_character_
    ),
    output = c(
      rep("$results", 8),
      "$steiger",
      "$pleiotropy",
      "$heterogeneity",
      "$loo"
    ),
    model = c(
      "random",
      "fixed",
      "random",
      rep(NA_character_, 9)
    ),
    ld_correctable = c(TRUE, TRUE, TRUE, rep(FALSE, 6), TRUE, TRUE, TRUE),
    min_instruments = c(2L, 2L, 3L, 3L, 3L, 2L, 1L, 2L, 1L, 3L, 2L, 3L),
    engine = c(
      "TwoSampleMR::mr_ivw",
      "TwoSampleMR::mr_ivw_fe",
      "TwoSampleMR::mr_egger_regression",
      "TwoSampleMR::mr_weighted_median",
      "TwoSampleMR::run_mr_presso",
      "MendelianRandomization::mr_conmix",
      "TwoSampleMR::mr_wald_ratio",
      "TwoSampleMR::mr",
      "TwoSampleMR::steiger_filtering",
      "TwoSampleMR::mr_pleiotropy_test",
      "TwoSampleMR::mr_heterogeneity",
      "TwoSampleMR::mr_leaveoneout"
    ),
    engine_ld = c(
      "MendelianRandomization::mr_ivw(model = \"random\")",
      "MendelianRandomization::mr_ivw(model = \"fixed\")",
      "MendelianRandomization::mr_egger",
      rep(NA_character_, 6),
      "MendelianRandomization::mr_egger()@Intercept",
      "MendelianRandomization::mr_ivw()@Heter.Stat / mr_egger()@Heter.Stat",
      "mrpipeline:::loo_correlated (block-inverse GLS per dropped SNP)"
    ),
    stringsAsFactors = FALSE
  )
}

#' Methods available to run_mr()
#'
#' Returns the table of methods [run_mr()] understands: what each shortcut
#' runs, whether it is a fixed- or random-effects estimator, whether
#' `ld_correct = TRUE` applies to it, how many instruments it needs, and
#' where its result lands on the returned `mr_result`. The table is the same
#' object `run_mr()` dispatches from (see [mr_method_registry()]), so it is
#' always current.
#'
#' @param detail `"concise"` (default) lists the ten shortcuts that can be
#'   passed to `methods=`. `"full"` adds the two paths that have no single
#'   name -- the automatic Wald ratio used at exactly one instrument, and the
#'   raw TwoSampleMR passthrough -- plus the `engine` and `engine_ld` columns
#'   naming the function that runs on each LD path.
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{`shortcut`}{Name accepted by `run_mr(methods = )`. `NA` for the
#'       two rows that cannot be requested by a single name.}
#'     \item{`description`}{What the method does.}
#'     \item{`label`}{The `method` value written to `$results`; `NA` for
#'       diagnostics, which produce no results row.}
#'     \item{`output`}{Where the result lands: `$results`, `$steiger`,
#'       `$pleiotropy`, `$heterogeneity` or `$loo`.}
#'     \item{`model`}{`"random"`, `"fixed"`, or `NA` where the distinction
#'       does not apply. Random effects are multiplicative: the standard
#'       error is inflated by `max(RSE, 1)` and never deflated.}
#'     \item{`ld_correctable`}{Whether `ld_correct = TRUE` changes this
#'       method. Methods with `FALSE` run on uncorrected data and warn.}
#'     \item{`min_instruments`}{Minimum number of instruments; below it the
#'       method is skipped and `$methods_skipped` records why.}
#'     \item{`engine`, `engine_ld`}{(`detail = "full"` only) The function
#'       mrpipeline calls when `ld_correct` is `FALSE` / `TRUE`. Documentary
#'       -- users never pass an engine; `run_mr()` picks it from the
#'       shortcut and `ld_correct`.}
#'   }
#'
#' @examples
#' mr_methods()
#' mr_methods(detail = "full")
#'
#' @seealso [run_mr()]
#' @export
mr_methods <- function(detail = c("concise", "full")) {
  detail <- rlang::arg_match(detail)
  reg <- mr_method_registry()
  if (detail == "full") {
    return(reg)
  }
  concise <- reg[
    !is.na(reg$shortcut),
    c(
      "shortcut",
      "description",
      "label",
      "output",
      "model",
      "ld_correctable",
      "min_instruments"
    )
  ]
  rownames(concise) <- NULL
  concise
}

#' Look up one registry row
#'
#' @param shortcut,label,engine Exactly one of these identifies the row:
#'   `shortcut` for the ten named methods, `label = "Wald ratio"` for the
#'   automatic single-instrument path, `engine = "TwoSampleMR::mr"` for the
#'   raw passthrough.
#'
#' @return A one-row data frame from [mr_method_registry()].
#'
#' @keywords internal
mr_method_entry <- function(shortcut = NULL, label = NULL, engine = NULL) {
  reg <- mr_method_registry()
  hit <- if (!is.null(shortcut)) {
    reg$shortcut %in% shortcut
  } else if (!is.null(label)) {
    reg$label %in% label
  } else {
    reg$engine %in% engine
  }
  reg[which(hit)[1], , drop = FALSE]
}

#' Build one `$results` row from a registry entry
#'
#' Every `$results` row `run_mr()` produces goes through here, so the schema
#' -- including the `ld_corrected` and `model` columns that say which
#' estimator ran -- cannot differ between dispatch branches (which would make
#' the final `rbind` fail).
#'
#' @param entry One-row data frame from [mr_method_entry()].
#' @param exposure_id,outcome_id Labels for the `exposure`/`outcome` columns.
#' @param nsnp,b,se,pval The estimate.
#' @param ld_corrected Logical. Whether the LD matrix was used for this fit.
#' @param label The `method` label. Defaults to the registry label; the raw
#'   TwoSampleMR passthrough supplies TwoSampleMR's own.
#'
#' @return A one-row data frame.
#'
#' @keywords internal
mr_result_row <- function(
  entry,
  exposure_id,
  outcome_id,
  nsnp,
  b,
  se,
  pval,
  ld_corrected,
  label = entry$label
) {
  data.frame(
    exposure = exposure_id,
    outcome = outcome_id,
    method = label,
    nsnp = nsnp,
    b = b,
    se = se,
    pval = pval,
    ld_corrected = ld_corrected,
    model = entry$model,
    stringsAsFactors = FALSE
  )
}

#' Render the concise method table as roxygen lines
#'
#' Used through `@eval` in [run_mr()]'s documentation, so `?run_mr` is
#' regenerated from the registry on every `devtools::document()` and cannot
#' drift from what `run_mr()` accepts.
#'
#' @return Character vector of roxygen lines (a markdown table inside a
#'   section).
#'
#' @keywords internal
rd_method_table <- function() {
  tab <- mr_methods()
  cell <- function(x) ifelse(is.na(x), "--", as.character(x))
  rows <- sprintf(
    "| `%s` | %s | %s | `%s` | %s | %s | %s |",
    tab$shortcut,
    tab$description,
    cell(tab$label),
    tab$output,
    cell(tab$model),
    ifelse(tab$ld_correctable, "yes", "no"),
    tab$min_instruments
  )
  c(
    "@section Available methods:",
    "Rendered from [mr_methods()]. `label` is the `method` value written to",
    "`$results`; `model` the effects model; `LD` whether `ld_correct = TRUE`",
    "applies; `min n` the minimum number of instruments.",
    "",
    "| shortcut | description | label | output | model | LD | min n |",
    "|---|---|---|---|---|---|---|",
    rows
  )
}
