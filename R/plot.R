#' Plot an mr_result object
#'
#' Creates diagnostic plots for Mendelian randomisation results using
#' TwoSampleMR plotting functions. Requires the `ggplot2` package.
#'
#' @param x An `mr_result` object.
#' @param type Character. Plot type: `"scatter"` (default), `"forest"`, or
#'   `"funnel"`.
#' @param ... Ignored.
#'
#' @return A ggplot object (or list of ggplot objects for `"scatter"`).
#'   Returns `NULL` invisibly if the result status is not `"success"` or if
#'   no results are available.
#'
#' @export
plot.mr_result <- function(x, type = c("scatter", "forest", "funnel"), ...) {
  rlang::check_installed("ggplot2", reason = "to plot MR results.")
  type <- match.arg(type)

  if (x$status != "success" || nrow(x$results) == 0) {
    cli::cli_inform("Cannot plot: no successful MR results available.")
    return(invisible(NULL))
  }

  harmonised <- x$instruments

  # TwoSampleMR plot functions expect id.exposure and id.outcome in the

  # results data frame; run_mr() stores only exposure/outcome names, so
  # we pull the IDs from the harmonised instruments.
  mr_res <- x$results
  mr_res$id.exposure <- harmonised$id.exposure[1]
  mr_res$id.outcome <- harmonised$id.outcome[1]

  if (type == "scatter") {
    TwoSampleMR::mr_scatter_plot(mr_res, harmonised)
  } else {
    singlesnp <- TwoSampleMR::mr_singlesnp(harmonised)
    if (type == "forest") {
      TwoSampleMR::mr_forest_plot(singlesnp)
    } else {
      TwoSampleMR::mr_funnel_plot(singlesnp)
    }
  }
}

#' Plot a coloc_result object
#'
#' Creates visualisations for colocalization results. Requires the `ggplot2`
#' package.
#'
#' @param x A `coloc_result` object.
#' @param type Character. Plot type: `"pp_bar"` (default) for a bar chart of
#'   ABF posterior probabilities (H0--H4), or `"regional"` for side-by-side
#'   regional association plots.
#' @param ... Ignored.
#'
#' @return A ggplot object. Returns `NULL` invisibly if the result status is
#'   not `"success"` or if the required data is unavailable.
#'
#' @export
plot.coloc_result <- function(x, type = c("pp_bar", "regional"), ...) {
  rlang::check_installed("ggplot2", reason = "to plot coloc results.")
  type <- match.arg(type)

  if (x$status != "success") {
    cli::cli_inform("Cannot plot: no successful coloc results available.")
    return(invisible(NULL))
  }

  switch(
    type,
    pp_bar = plot_coloc_pp_bar(x),
    regional = plot_coloc_regional(x)
  )
}

#' Bar chart of ABF posterior probabilities
#' @param x A `coloc_result` object.
#' @return A ggplot object or `NULL` invisibly.
#' @keywords internal
plot_coloc_pp_bar <- function(x) {
  if (is.null(x$coloc_abf)) {
    cli::cli_inform("No ABF results available for plotting.")
    return(invisible(NULL))
  }

  s <- x$coloc_abf$summary
  pp_names <- c("PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")
  pp_labels <- c("H0", "H1", "H2", "H3", "H4")

  plot_data <- data.frame(
    hypothesis = factor(pp_labels, levels = pp_labels),
    probability = as.numeric(s[pp_names])
  )

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data$hypothesis, y = .data$probability)
  ) +
    ggplot2::geom_col(fill = "steelblue") +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      title = "Colocalization posterior probabilities (ABF)",
      x = "Hypothesis",
      y = "Posterior probability"
    ) +
    ggplot2::theme_minimal()
}

#' Regional association plot for coloc results
#' @param x A `coloc_result` object.
#' @return A ggplot object or `NULL` invisibly.
#' @keywords internal
plot_coloc_regional <- function(x) {
  dat <- x$harmonised_data
  if (is.null(dat) || nrow(dat) == 0) {
    cli::cli_inform("No harmonised data available for regional plot.")
    return(invisible(NULL))
  }

  # Build long-format data for faceted plot
  exposure_dat <- data.frame(
    pos = dat$pos.exposure,
    log10p = -log10(dat$pval.exposure),
    panel = "Exposure"
  )
  outcome_dat <- data.frame(
    pos = dat$pos.outcome,
    log10p = -log10(dat$pval.outcome),
    panel = "Outcome"
  )
  plot_data <- rbind(exposure_dat, outcome_dat)

  ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$pos, y = .data$log10p)) +
    ggplot2::geom_point(alpha = 0.6) +
    ggplot2::facet_wrap(~panel, ncol = 1, scales = "free_y") +
    ggplot2::labs(
      title = "Regional association",
      x = "Position (bp)",
      y = expression(-log[10](p))
    ) +
    ggplot2::theme_minimal()
}
