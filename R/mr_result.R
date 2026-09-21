#' Create an mr_result object
#'
#' @param results Data frame with columns: exposure, outcome, method, nsnp,
#'   b, se, pval, ld_corrected (logical: whether the LD matrix was used for
#'   that fit), model (`"random"`, `"fixed"` or `NA`), plus or, or_lci95,
#'   or_uci95 (and lo_ci, up_ci) from [TwoSampleMR::generate_odds_ratios()]
#'   on a successful run. Defaults to [empty_mr_results()].
#' @param instruments Data frame of harmonised (and clumped) instrument data:
#'   the kept variants the MR estimates are computed from.
#' @param harmonisation Data frame. The complete, unfiltered
#'   [TwoSampleMR::harmonise_data()] output -- every candidate variant, with
#'   the `mr_keep`, `palindromic`, `ambiguous` and `remove` flags that explain
#'   why each one was or was not carried forward. `instruments` is the subset
#'   of this that survived; see [harmonisation_summary()].
#' @param f_stats List with elements `per_snp` (numeric vector),
#'   `mean` (numeric scalar), `min` (numeric scalar).
#' @param steiger Output of [TwoSampleMR::steiger_filtering()], or `NULL`.
#' @param pleiotropy Output of [TwoSampleMR::mr_pleiotropy_test()] plus an
#'   `ld_corrected` column, or `NULL`. When `ld_corrected` is `TRUE` the
#'   intercept is the correlated Egger fit's ([pleiotropy_correlated()]).
#' @param heterogeneity Output of [TwoSampleMR::mr_heterogeneity()] (Cochran's
#'   Q per method) plus an `ld_corrected` column, or `NULL`. When
#'   `ld_corrected` is `TRUE` the Q values are the generalised statistics
#'   from the correlated fits ([heterogeneity_correlated()]).
#' @param loo Output of [TwoSampleMR::mr_leaveoneout()] plus an
#'   `ld_corrected` column, or `NULL`. When `ld_corrected` is `TRUE` each row
#'   is a correlated random-effects refit ([loo_correlated()]).
#' @param methods_skipped Named character vector: names are method names,
#'   values are reasons for skipping.
#' @param ld_matrix LD correlation matrix if `ld_correct = TRUE`, or `NULL`.
#' @param params List of all input parameters to `run_mr()`.
#' @param status Character. One of `"success"`, `"no_instruments"`,
#'   `"no_harmonised_variants"`, `"singular_ld_matrix"`. Default
#'   `"success"`.
#' @param status_reason Character or `NULL`. Human-readable explanation when
#'   `status != "success"`.
#' @param timing Named numeric vector of elapsed times (seconds) for each
#'   major step inside `run_mr()`. Empty by default.
#'
#' @return An object of class `mr_result`.
#'
#' @keywords internal
new_mr_result <- function(
  results = empty_mr_results(),
  instruments = data.frame(),
  harmonisation = data.frame(),
  f_stats = list(per_snp = numeric(), mean = NA_real_, min = NA_real_),
  steiger = NULL,
  pleiotropy = NULL,
  heterogeneity = NULL,
  loo = NULL,
  methods_skipped = character(),
  ld_matrix = NULL,
  params = list(),
  status = "success",
  status_reason = NULL,
  timing = numeric(0)
) {
  structure(
    list(
      results = results,
      instruments = instruments,
      harmonisation = harmonisation,
      f_stats = f_stats,
      steiger = steiger,
      pleiotropy = pleiotropy,
      heterogeneity = heterogeneity,
      loo = loo,
      methods_skipped = methods_skipped,
      ld_matrix = ld_matrix,
      params = params,
      status = status,
      status_reason = status_reason,
      timing = timing
    ),
    class = "mr_result"
  )
}

#' Empty `$results` frame carrying the full column schema
#'
#' Used as the `results` default of [new_mr_result()] (so every early return
#' has the same shape as a successful run) and by `run_mr()` when no method
#' produced an estimate. The OR columns are appended by
#' `TwoSampleMR::generate_odds_ratios()` in `run_mr()`, not here.
#'
#' @return A zero-row data frame with columns `exposure`, `outcome`,
#'   `method`, `nsnp`, `b`, `se`, `pval`, `ld_corrected`, `model`.
#'
#' @keywords internal
empty_mr_results <- function() {
  data.frame(
    exposure = character(),
    outcome = character(),
    method = character(),
    nsnp = integer(),
    b = numeric(),
    se = numeric(),
    pval = numeric(),
    ld_corrected = logical(),
    model = character(),
    stringsAsFactors = FALSE
  )
}

#' Print an mr_result object
#'
#' Displays a one-line summary: exposure -> outcome, primary estimate,
#' number of SNPs, and mean F-statistic.
#'
#' @param x An `mr_result` object.
#' @param ... Ignored.
#'
#' @return Invisibly returns `x`.
#'
#' @examples
#' \dontrun{
#' bfile <- sub("\\.bed$", "", system.file("extdata", "ld_ref.bed", package = "mrpipeline"))
#' result <- run_mr(
#'   exposure = cd40_exposure, exposure_id = "CD40",
#'   outcome = sjogren_outcome, outcome_id = "SjD",
#'   instrument_region = list(chromosome = "20", start = 44746911, end = 44758502),
#'   rsq_thresh = 0.3,
#'   bfile = bfile
#' )
#' print(result)
#' }
#'
#' @export
print.mr_result <- function(x, ...) {
  if (x$status != "success") {
    reason <- x$status_reason %||% "unknown reason" # nolint: object_usage_linter.
    cli::cli_inform("mr_result: {x$status} \u2014 {reason}")
    return(invisible(x))
  }

  if (nrow(x$results) == 0) {
    cli::cli_inform("mr_result: no results (all methods failed or skipped)")
    return(invisible(x))
  }

  # Use first result row as primary estimate

  # object_usage_linter false positives below: every variable here is used
  # via cli's glue-style "{var}" string interpolation, which lintr's static
  # analysis cannot trace.
  primary <- x$results[1, ]
  nsnp <- primary$nsnp # nolint: object_usage_linter.
  mean_f <- x$f_stats$mean # nolint: object_usage_linter.

  # nolint next: object_usage_linter.
  or_str <- if ("or" %in% names(primary) && !is.na(primary$or)) {
    paste0(
      ", OR = ",
      round(primary$or, 3),
      " [",
      round(primary$or_lci95, 3),
      "-",
      round(primary$or_uci95, 3),
      "]"
    )
  } else {
    ""
  }

  # The shortest view still says which estimator ran.
  ld_tag <- if (isTRUE(primary$ld_corrected)) " [LD-corrected]" else "" # nolint: object_usage_linter.

  cli::cli_inform(c(
    "{primary$exposure} -> {primary$outcome}",
    "i" = paste0(
      "{primary$method}{ld_tag}: b = {round(primary$b, 4)}, ",
      "se = {round(primary$se, 4)}, p = {signif(primary$pval, 3)}{or_str}"
    ),
    "i" = "{nsnp} SNP{?s}, mean F = {round(mean_f, 1)}"
  ))

  invisible(x)
}

#' Summarise an mr_result object
#'
#' Displays the full results table, F-statistics, Steiger summary,
#' Egger intercept (if available), and skipped methods.
#'
#' @param object An `mr_result` object.
#' @param ... Ignored.
#'
#' @return Invisibly returns `object`.
#'
#' @examples
#' \dontrun{
#' bfile <- sub("\\.bed$", "", system.file("extdata", "ld_ref.bed", package = "mrpipeline"))
#' result <- run_mr(
#'   exposure = cd40_exposure, exposure_id = "CD40",
#'   outcome = sjogren_outcome, outcome_id = "SjD",
#'   instrument_region = list(chromosome = "20", start = 44746911, end = 44758502),
#'   rsq_thresh = 0.3,
#'   bfile = bfile
#' )
#' summary(result)
#' }
#'
#' @export
summary.mr_result <- function(object, ...) {
  cli::cli_h1(
    "MR Results: {object$params$exposure_id} -> {object$params$outcome_id}"
  )

  if (object$status != "success") {
    reason <- object$status_reason %||% "unknown reason" # nolint: object_usage_linter.
    cli::cli_alert_warning("Status: {object$status} \u2014 {reason}")
    return(invisible(object))
  }

  if (nrow(object$results) == 0) {
    cli::cli_alert_warning("No results available.")
    return(invisible(object))
  }

  # Results table
  cli::cli_h2("Method estimates")
  res <- object$results
  ld_run <- !is.null(object$ld_matrix)
  for (i in seq_len(nrow(res))) {
    # nolint next: object_usage_linter.
    or_str <- if ("or" %in% names(res) && !is.na(res$or[i])) {
      paste0(
        ", OR = ",
        round(res$or[i], 3),
        " [",
        round(res$or_lci95[i], 3),
        "-",
        round(res$or_uci95[i], 3),
        "]"
      )
    } else {
      ""
    }
    # Say which estimator ran: the effects model where the label does not
    # already carry it (MR Egger), and the LD status on an LD-corrected run.
    tags <- character()
    model_i <- (res$model %||% NA_character_)[i]
    if (
      !is.na(model_i) &&
        !stringr::str_detect(res$method[i], stringr::fixed(model_i))
    ) {
      tags <- c(tags, paste(model_i, "effects"))
    }
    if (ld_run) {
      tags <- c(
        tags,
        if (isTRUE(res$ld_corrected[i])) "LD-corrected" else "not LD-corrected"
      )
    }
    # nolint next: object_usage_linter.
    tag_str <- if (length(tags) > 0) {
      paste0(" [", paste(tags, collapse = ", "), "]")
    } else {
      ""
    }
    cli::cli_bullets(c(
      "*" = paste0(
        "{res$method[i]}{tag_str}: b = {round(res$b[i], 4)}, ",
        "se = {round(res$se[i], 4)}, p = {signif(res$pval[i], 3)}{or_str} ",
        "({res$nsnp[i]} SNPs)"
      )
    ))
  }

  print_harmonisation_summary(object$harmonisation)

  # F-statistics
  cli::cli_h2("Instrument strength")
  cli::cli_bullets(c(
    "*" = "Mean F-statistic: {round(object$f_stats$mean, 1)}",
    "*" = "Min F-statistic: {round(object$f_stats$min, 1)}",
    "*" = "N instruments: {length(object$f_stats$per_snp)}"
  ))

  # Steiger: $steiger is TwoSampleMR::steiger_filtering() output, one row
  # per SNP with a logical steiger_dir, so summarise it per SNP.
  if (!is.null(object$steiger)) {
    st <- object$steiger
    # nolint start: object_usage_linter.
    n_dir <- sum(st$steiger_dir, na.rm = TRUE)
    wrong <- st$SNP[!st$steiger_dir %in% TRUE]
    max_p <- signif(max(st$steiger_pval, na.rm = TRUE), 3)
    # nolint end
    cli::cli_h2("Steiger filtering")
    cli::cli_bullets(c(
      "*" = paste0(
        "{n_dir}/{nrow(st)} SNP{?s} explain more variance in the exposure ",
        "than in the outcome"
      ),
      "*" = "Largest Steiger p-value: {max_p}"
    ))
    if (length(wrong) > 0) {
      cli::cli_bullets(c("!" = "Not in the expected direction: {.val {wrong}}"))
    }
  }

  # Each diagnostics frame says whether it came from the correlated fits;
  # the heading repeats that so it cannot be read as the uncorrected one.
  # lintr cannot see uses inside cli's {} strings; a block rather than a
  # trailing # nolint, which air format moves off the flagged line.
  # nolint start: object_usage_linter.
  diag_tag <- function(frame) {
    if (isTRUE(frame$ld_corrected[1])) " [LD-corrected]" else ""
  }
  # nolint end

  # Pleiotropy test
  if (!is.null(object$pleiotropy)) {
    pt <- object$pleiotropy # nolint: object_usage_linter.
    cli::cli_h2("Pleiotropy test (Egger intercept){diag_tag(pt)}")
    cli::cli_bullets(c(
      "*" = "Intercept: {round(pt$egger_intercept, 4)}",
      "*" = "SE: {round(pt$se, 4)}",
      "*" = "p-value: {signif(pt$pval, 3)}"
    ))
  }

  # Heterogeneity test (Cochran's Q)
  if (!is.null(object$heterogeneity)) {
    ht <- object$heterogeneity
    cli::cli_h2("Heterogeneity test (Cochran's Q){diag_tag(ht)}")
    for (i in seq_len(nrow(ht))) {
      cli::cli_bullets(c(
        "*" = "{ht$method[i]}: Q = {round(ht$Q[i], 3)}, df = {ht$Q_df[i]}, p = {signif(ht$Q_pval[i], 3)}"
      ))
    }
  }

  # Leave-one-out analysis
  if (!is.null(object$loo)) {
    cli::cli_h2("Leave-one-out analysis{diag_tag(object$loo)}")
    cli::cli_bullets(c(
      "*" = paste0(
        "{nrow(object$loo)} row{?s} (per-SNP estimates plus the pooled ",
        "'All' row); see {.code $loo} for the full table."
      )
    ))
  }

  # Skipped methods
  if (length(object$methods_skipped) > 0) {
    cli::cli_h2("Skipped methods")
    for (nm in names(object$methods_skipped)) {
      cli::cli_bullets(c("!" = "{nm}: {object$methods_skipped[[nm]]}"))
    }
  }

  # LD correction -- keyed off what each row records, not off the presence
  # of the matrix, so a run whose methods all lacked a correlated form is
  # never announced as "LD-corrected".
  if (ld_run) {
    cli::cli_h2("LD correction")
    n_ld <- nrow(object$ld_matrix) # nolint: object_usage_linter.
    if (length(object$f_stats$per_snp) == 1) {
      cli::cli_bullets(c(
        "i" = "Not applicable: 1 instrument (Wald ratio)"
      ))
    } else {
      applied <- res$method[res$ld_corrected %in% TRUE]
      not_applied <- res$method[!(res$ld_corrected %in% TRUE)]
      applied_str <- paste(applied, collapse = ", ") # nolint: object_usage_linter.
      if (length(applied) == 0) {
        applied_str <- "none"
      }
      lines <- c(
        "i" = "{n_ld}-SNP LD matrix from the reference panel",
        "v" = "Applied to: {applied_str}"
      )
      if (length(not_applied) > 0) {
        not_applied_str <- paste(not_applied, collapse = ", ") # nolint: object_usage_linter.
        lines <- c(lines, "!" = "Not applied to: {not_applied_str}")
      }
      diag_names <- c(
        pleiotropy = "Egger intercept",
        heterogeneity = "Cochran's Q",
        loo = "leave-one-out"
      )
      diag_on <- vapply(
        names(diag_names),
        function(f) isTRUE(object[[f]]$ld_corrected[1]),
        logical(1)
      )
      if (any(diag_on)) {
        diag_str <- paste(diag_names[diag_on], collapse = ", ") # nolint: object_usage_linter.
        lines <- c(
          lines,
          "v" = "Diagnostics from the correlated fits: {diag_str}"
        )
      }
      cli::cli_bullets(lines)
    }
  }

  invisible(object)
}
