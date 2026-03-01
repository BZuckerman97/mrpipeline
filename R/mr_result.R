#' Create an mr_result object
#'
#' @param results Data frame with columns: exposure, outcome, method, nsnp,
#'   b, se, pval.
#' @param instruments Data frame of harmonised (and clumped) instrument data.
#' @param f_stats List with elements `per_snp` (numeric vector),
#'   `mean` (numeric scalar), `min` (numeric scalar).
#' @param steiger Output of [TwoSampleMR::steiger_filtering()], or `NULL`.
#' @param methods_skipped Named character vector: names are method names,
#'   values are reasons for skipping.
#' @param ld_matrix LD correlation matrix if `ld_correct = TRUE`, or `NULL`.
#' @param params List of all input parameters to `run_mr()`.
#' @param status Character. One of `"success"`, `"no_instruments"`,
#'   `"no_harmonised_variants"`. Default `"success"`.
#' @param status_reason Character or `NULL`. Human-readable explanation when
#'   `status != "success"`.
#' @param timing Named numeric vector of elapsed times (seconds) for each
#'   major step inside `run_mr()`. Empty by default.
#'
#' @return An object of class `mr_result`.
#'
#' @keywords internal
new_mr_result <- function(
  results = data.frame(),
  instruments = data.frame(),
  f_stats = list(per_snp = numeric(), mean = NA_real_, min = NA_real_),
  steiger = NULL,
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
      f_stats = f_stats,
      steiger = steiger,
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
#' @export
print.mr_result <- function(x, ...) {
  if (x$status != "success") {
    reason <- x$status_reason %||% "unknown reason"
    cli::cli_inform("mr_result: {x$status} \u2014 {reason}")
    return(invisible(x))
  }

  if (nrow(x$results) == 0) {
    cli::cli_inform("mr_result: no results (all methods failed or skipped)")
    return(invisible(x))
  }

  # Use first result row as primary estimate

  primary <- x$results[1, ]
  nsnp <- primary$nsnp
  mean_f <- x$f_stats$mean

  cli::cli_inform(c(
    "{primary$exposure} -> {primary$outcome}",
    "i" = "{primary$method}: b = {round(primary$b, 4)}, se = {round(primary$se, 4)}, p = {signif(primary$pval, 3)}",
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
#' @export
summary.mr_result <- function(object, ...) {
  cli::cli_h1(
    "MR Results: {object$params$exposure_id} -> {object$params$outcome_id}"
  )

  if (object$status != "success") {
    reason <- object$status_reason %||% "unknown reason"
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
  for (i in seq_len(nrow(res))) {
    cli::cli_bullets(c(
      "*" = "{res$method[i]}: b = {round(res$b[i], 4)}, se = {round(res$se[i], 4)}, p = {signif(res$pval[i], 3)} ({res$nsnp[i]} SNPs)"
    ))
  }

  # F-statistics
  cli::cli_h2("Instrument strength")
  cli::cli_bullets(c(
    "*" = "Mean F-statistic: {round(object$f_stats$mean, 1)}",
    "*" = "Min F-statistic: {round(object$f_stats$min, 1)}",
    "*" = "N instruments: {length(object$f_stats$per_snp)}"
  ))

  # Steiger
  if (!is.null(object$steiger)) {
    cli::cli_h2("Steiger filtering")
    cli::cli_bullets(c(
      "*" = "Correct direction: {object$steiger$correct_causal_direction}",
      "*" = "Steiger p-value: {signif(object$steiger$steiger_pval, 3)}"
    ))
  }

  # Skipped methods
  if (length(object$methods_skipped) > 0) {
    cli::cli_h2("Skipped methods")
    for (nm in names(object$methods_skipped)) {
      cli::cli_bullets(c("!" = "{nm}: {object$methods_skipped[[nm]]}"))
    }
  }

  # LD correction
  if (!is.null(object$ld_matrix)) {
    cli::cli_alert_info(
      "LD-corrected analysis ({nrow(object$ld_matrix)} SNPs in LD matrix)"
    )
  }

  invisible(object)
}
