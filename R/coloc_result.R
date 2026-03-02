#' Create a coloc_result object
#'
#' @param coloc_abf Output of [coloc::coloc.abf()], or `NULL`.
#' @param coloc_susie Output of [coloc::coloc.susie()], or `NULL`.
#' @param coloc_signals Output of [coloc::coloc.signals()], or `NULL`.
#' @param coloc_prop_test Output of `colocPropTest::coloc.prop.test()`, or
#'   `NULL`.
#' @param n_snps Integer. Number of SNPs used in the analysis.
#' @param harmonised_data Data frame of harmonised data used for analysis.
#' @param methods_skipped Named character vector: names are method names,
#'   values are reasons for skipping.
#' @param params List of all input parameters to `run_coloc()`.
#' @param status Character. One of `"success"`, `"no_snps_in_region"`,
#'   `"too_few_snps"`, `"no_harmonised_variants"`. Default `"success"`.
#' @param status_reason Character or `NULL`. Human-readable explanation when
#'   `status != "success"`.
#' @param timing Named numeric vector of elapsed times (seconds) for each
#'   major step inside `run_coloc()`. Empty by default.
#'
#' @return An object of class `coloc_result`.
#'
#' @keywords internal
new_coloc_result <- function(
  coloc_abf = NULL,
  coloc_susie = NULL,
  coloc_signals = NULL,
  coloc_prop_test = NULL,
  n_snps = 0L,
  harmonised_data = data.frame(),
  methods_skipped = character(),
  params = list(),
  status = "success",
  status_reason = NULL,
  timing = numeric(0)
) {
  structure(
    list(
      coloc_abf = coloc_abf,
      coloc_susie = coloc_susie,
      coloc_signals = coloc_signals,
      coloc_prop_test = coloc_prop_test,
      n_snps = as.integer(n_snps),
      harmonised_data = harmonised_data,
      methods_skipped = methods_skipped,
      params = params,
      status = status,
      status_reason = status_reason,
      timing = timing
    ),
    class = "coloc_result"
  )
}

#' Print a coloc_result object
#'
#' Displays a one-line summary: number of SNPs, ABF PP.H4 if available,
#' SuSiE max PP.H4 and credible set count if available. Shows status if
#' not success.
#'
#' @param x A `coloc_result` object.
#' @param ... Ignored.
#'
#' @return Invisibly returns `x`.
#'
#' @examples
#' \dontrun{
#' result <- run_coloc(
#'   exposure = cd40_exposure,
#'   outcome = sjogren_outcome,
#'   gene_chr = "20", gene_start = 44746911, gene_end = 44758502,
#'   bfile = system.file("extdata", "ld_ref", package = "mrpipeline")
#' )
#' print(result)
#' }
#'
#' @export
print.coloc_result <- function(x, ...) {
  if (x$status != "success") {
    reason <- x$status_reason %||% "unknown reason"
    cli::cli_inform("coloc_result: {x$status} \u2014 {reason}")
    return(invisible(x))
  }

  parts <- character()

  # ABF summary
  if (!is.null(x$coloc_abf)) {
    pp_h4 <- x$coloc_abf$summary["PP.H4.abf"]
    parts <- c(parts, "ABF PP.H4 = {round(pp_h4, 3)}")
  }

  # SuSiE summary
  if (!is.null(x$coloc_susie) && !is.null(x$coloc_susie$summary)) {
    susie_summary <- x$coloc_susie$summary
    if (is.data.frame(susie_summary) && nrow(susie_summary) > 0) {
      max_h4 <- max(susie_summary$PP.H4.abf, na.rm = TRUE)
      n_pairs <- nrow(susie_summary)
      parts <- c(
        parts,
        "SuSiE max PP.H4 = {round(max_h4, 3)} ({n_pairs} credible set pair{?s})"
      )
    }
  }

  if (length(parts) == 0) {
    cli::cli_inform("coloc_result: {x$n_snps} SNPs, no method results")
  } else {
    header <- "coloc_result: {x$n_snps} SNPs"
    cli::cli_inform(c(header, stats::setNames(parts, rep("i", length(parts)))))
  }

  invisible(x)
}

#' Summarise a coloc_result object
#'
#' Displays full details: region info, all posterior probabilities (H0-H4),
#' SuSiE credible set pairs, signals hits, and skipped methods.
#'
#' @param object A `coloc_result` object.
#' @param ... Ignored.
#'
#' @return Invisibly returns `object`.
#'
#' @examples
#' \dontrun{
#' result <- run_coloc(
#'   exposure = cd40_exposure,
#'   outcome = sjogren_outcome,
#'   gene_chr = "20", gene_start = 44746911, gene_end = 44758502,
#'   bfile = system.file("extdata", "ld_ref", package = "mrpipeline")
#' )
#' summary(result)
#' }
#'
#' @export
summary.coloc_result <- function(object, ...) {
  cli::cli_h1("Colocalization Results")

  if (object$status != "success") {
    reason <- object$status_reason %||% "unknown reason"
    cli::cli_alert_warning("Status: {object$status} \u2014 {reason}")
    return(invisible(object))
  }

  cli::cli_alert_info("{object$n_snps} SNPs in analysis")

  # ABF results
  if (!is.null(object$coloc_abf)) {
    cli::cli_h2("coloc.abf")
    s <- object$coloc_abf$summary
    cli::cli_bullets(c(
      "*" = "PP.H0 = {round(s['PP.H0.abf'], 4)}",
      "*" = "PP.H1 = {round(s['PP.H1.abf'], 4)}",
      "*" = "PP.H2 = {round(s['PP.H2.abf'], 4)}",
      "*" = "PP.H3 = {round(s['PP.H3.abf'], 4)}",
      "*" = "PP.H4 = {round(s['PP.H4.abf'], 4)}",
      "*" = "N SNPs = {s['nsnps']}"
    ))
  }

  # SuSiE results
  if (!is.null(object$coloc_susie) && !is.null(object$coloc_susie$summary)) {
    cli::cli_h2("coloc.susie")
    susie_summary <- object$coloc_susie$summary
    if (is.data.frame(susie_summary) && nrow(susie_summary) > 0) {
      sort_cols <- intersect(
        c("PP.H4.abf", "PP.H3.abf", "PP.H2.abf", "PP.H1.abf"),
        names(susie_summary)
      )
      ord <- do.call(order, lapply(susie_summary[sort_cols], `-`))
      susie_sorted <- susie_summary[ord, , drop = FALSE]
      cap <- 20L
      n_total <- nrow(susie_sorted)
      shown <- utils::head(susie_sorted, cap)
      bullets <- stats::setNames(
        paste0(
          "Pair ",
          shown$idx1,
          "-",
          shown$idx2,
          ": PP.H4 = ",
          round(shown$PP.H4.abf, 4)
        ),
        rep("*", nrow(shown))
      )
      if (n_total > cap) {
        n_more <- n_total - cap
        rest <- utils::tail(susie_sorted, n_more)
        suffix <- if (all(rest$PP.H4.abf == 0)) " (all PP.H4 = 0)" else ""
        bullets <- c(
          bullets,
          c("i" = paste0("... and ", n_more, " more pair(s)", suffix))
        )
      }
      cli::cli_bullets(bullets)
    } else {
      cli::cli_alert_info("No credible set pairs found.")
    }
  }

  # Signals results
  if (!is.null(object$coloc_signals)) {
    cli::cli_h2("coloc.signals")
    signals_summary <- object$coloc_signals$summary
    if (is.data.frame(signals_summary) && nrow(signals_summary) > 0) {
      cli::cli_alert_info("{nrow(signals_summary)} signal pair{?s} tested")
      sort_cols <- intersect(
        c("PP.H4.abf", "PP.H3.abf", "PP.H2.abf", "PP.H1.abf"),
        names(signals_summary)
      )
      ord <- do.call(order, lapply(signals_summary[sort_cols], `-`))
      signals_sorted <- signals_summary[ord, , drop = FALSE]
      cap <- 20L
      n_total <- nrow(signals_sorted)
      shown <- utils::head(signals_sorted, cap)
      bullets <- stats::setNames(
        paste0(
          "Hit ",
          shown$hit1,
          "-",
          shown$hit2,
          ": PP.H4 = ",
          round(shown$PP.H4.abf, 4)
        ),
        rep("*", nrow(shown))
      )
      if (n_total > cap) {
        n_more <- n_total - cap
        rest <- utils::tail(signals_sorted, n_more)
        suffix <- if (all(rest$PP.H4.abf == 0)) " (all PP.H4 = 0)" else ""
        bullets <- c(
          bullets,
          c("i" = paste0("... and ", n_more, " more pair(s)", suffix))
        )
      }
      cli::cli_bullets(bullets)
    } else {
      cli::cli_alert_info("No signal pairs found.")
    }
  }

  # Prop test results
  if (!is.null(object$coloc_prop_test)) {
    cli::cli_h2("colocPropTest")
    cli::cli_alert_info("Proportionality test completed.")
  }

  # Skipped methods
  if (length(object$methods_skipped) > 0) {
    cli::cli_h2("Skipped methods")
    bullets <- stats::setNames(
      paste0(
        names(object$methods_skipped),
        ": ",
        unname(object$methods_skipped)
      ),
      rep("!", length(object$methods_skipped))
    )
    cli::cli_bullets(bullets)
  }

  invisible(object)
}
