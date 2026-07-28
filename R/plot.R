#' Plot an mr_result object
#'
#' Creates diagnostic plots for Mendelian randomisation results using
#' TwoSampleMR plotting functions. Requires the `ggplot2` package.
#'
#' @param x An `mr_result` object.
#' @param type Character. Plot type: `"scatter"` (default), `"forest"`,
#'   `"funnel"`, or `"loo"` (leave-one-out, via
#'   [TwoSampleMR::mr_leaveoneout_plot()] on `x$loo`; requires `"loo"` to
#'   have been in `methods` when `run_mr()` was called).
#' @param ... Ignored.
#'
#' @return A ggplot object (or list of ggplot objects for `"scatter"` and
#'   `"loo"`). Returns `NULL` invisibly if the result status is not
#'   `"success"`, if no results are available, or (for `type = "loo"`) if
#'   `x$loo` was not computed.
#'
#' @examples
#' \dontrun{
#' result <- run_mr(
#'   exposure = cd40_exposure, exposure_id = "CD40",
#'   outcome = sjogren_outcome, outcome_id = "SjD",
#'   instrument_region = list(chromosome = "20", start = 44746911, end = 44758502),
#'   bfile = system.file("extdata", "ld_ref", package = "mrpipeline"),
#'   methods = c("ivw", "egger", "weighted_median", "loo")
#' )
#' plot(result, type = "scatter")
#' plot(result, type = "forest")
#' plot(result, type = "funnel")
#' plot(result, type = "loo")
#' }
#'
#' @export
plot.mr_result <- function(x, type = c("scatter", "forest", "funnel", "loo"), ...) {
  rlang::check_installed("ggplot2", reason = "to plot MR results.")
  type <- match.arg(type)

  if (x$status != "success" || nrow(x$results) == 0) {
    cli::cli_inform("Cannot plot: no successful MR results available.")
    return(invisible(NULL))
  }

  if (type == "loo") {
    if (is.null(x$loo) || nrow(x$loo) == 0) {
      cli::cli_inform("Cannot plot: leave-one-out results not available (add \"loo\" to methods).")
      return(invisible(NULL))
    }
    return(TwoSampleMR::mr_leaveoneout_plot(x$loo))
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
#' @examples
#' \dontrun{
#' result <- run_coloc(
#'   exposure = cd40_exposure,
#'   outcome = sjogren_outcome,
#'   gene_chr = "20", gene_start = 44746911, gene_end = 44758502,
#'   bfile = system.file("extdata", "ld_ref", package = "mrpipeline")
#' )
#' plot(result, type = "pp_bar")
#' plot(result, type = "regional")
#' }
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

#' Forest plot of MR methods for one or more results
#'
#' One row per method (e.g. IVW fixed/random effects, MR Egger, weighted
#' median), for a single `mr_result` or several sectioned together. Wraps
#' [TwoSampleMR::forest_plot_1_to_many()]. Requires the `ggplot2` package.
#'
#' @param results An `mr_result`, or a list of `mr_result` objects. If the
#'   list has more than one element, every element must be named -- names
#'   become section labels (e.g. `list(Primary = res1, "Positive control" =
#'   res2)`). A single result (bare, or a length-1 list) is never sectioned,
#'   regardless of any name given.
#' @param methods Character vector of `method` values to include, in display
#'   order. Default puts fixed effects above random effects: `c("IVW (fixed
#'   effects)", "Inverse variance weighted", "MR Egger", "Weighted
#'   median")`.
#' @param relabel Named character vector, `c(old = new)`, applied to the
#'   `method` column for display only, after row filtering/ordering -- so
#'   matching against `results$results$method` is unaffected. Default
#'   relabels `"Inverse variance weighted"` to `"IVW (random effects)"`.
#' @param exponentiate Logical. Plot on the OR scale. Default `TRUE`.
#' @param xlab X-axis label. Defaults to `"OR (95% CI)"` if `exponentiate`,
#'   else `"Effect (95% CI)"`.
#' @param trans X-axis scale transform, e.g. `"identity"`, `"log2"`.
#'   Defaults to `"log2"` if `exponentiate`, else `"identity"`.
#' @param xlim Numeric vector `c(lower, upper)` giving x-axis limits, or
#'   `NULL` for automatic.
#' @param shape_points Point shape (see `ggplot2::geom_point()`'s `shape`).
#' @param colour_scheme Text/point colour. Note: font *family* is not
#'   configurable -- `TwoSampleMR::forest_plot_1_to_many()` has no such
#'   argument.
#' @param col_text_size Row-label font size.
#' @param subheading_size Section-heading font size.
#' @param col1_width Width of the row-label column.
#' @param ... Forwarded to `TwoSampleMR::forest_plot_1_to_many()` (e.g.
#'   `addcols`, `addcol_widths`, `addcol_titles`, `weight`).
#'
#' @return A grid plot object, or `NULL` invisibly if no result has
#'   successful, matching rows.
#'
#' @examples
#' \dontrun{
#' result <- run_mr(
#'   exposure = cd40_exposure, exposure_id = "CD40",
#'   outcome = sjogren_outcome, outcome_id = "SjD",
#'   instrument_region = list(chromosome = "20", start = 44746911, end = 44758502),
#'   bfile = system.file("extdata", "ld_ref", package = "mrpipeline"),
#'   methods = c("ivw", "ivw_fe", "egger", "weighted_median")
#' )
#'
#' # Single result: one row per method
#' forest_plot(result)
#'
#' # Several results sectioned together, e.g. primary analysis plus controls
#' forest_plot(list(
#'   Primary = result,
#'   "Positive control" = positive_control_result,
#'   "Negative control" = negative_control_result
#' ))
#' }
#'
#' @export
forest_plot <- function(
  results,
  methods = c(
    "IVW (fixed effects)",
    "Inverse variance weighted",
    "MR Egger",
    "Weighted median"
  ),
  relabel = c("Inverse variance weighted" = "IVW (random effects)"),
  exponentiate = TRUE,
  xlab = NULL,
  trans = NULL,
  xlim = NULL,
  shape_points = 15,
  colour_scheme = "black",
  col_text_size = 5,
  subheading_size = 6,
  col1_width = 1,
  ...
) {
  rlang::check_installed("ggplot2", reason = "to plot MR results.")

  if (inherits(results, "mr_result")) {
    results <- list(results)
  }
  if (length(results) > 1) {
    nm <- names(results)
    if (is.null(nm) || any(!nzchar(nm))) {
      cli::cli_abort(c(
        "All elements of {.arg results} must be named when plotting more than one result.",
        "i" = "e.g. {.code list(Primary = res1, \"Positive control\" = res2)}"
      ))
    }
  }
  sectioned <- length(results) > 1
  labels <- names(results) %||% rep("", length(results))

  pieces <- vector("list", length(results))
  for (i in seq_along(results)) {
    res <- results[[i]]
    label <- labels[i]
    if (res$status != "success" || nrow(res$results) == 0) {
      cli::cli_warn("Skipping {label}: no successful MR results available.")
      next
    }
    # Filtering/ordering uses the RAW method labels (as they appear in
    # res$results); relabelling happens afterwards, display-only.
    dat <- res$results |> dplyr::filter(method %in% methods)
    dat <- dat[match(methods, dat$method), , drop = FALSE]
    dat <- dat[!is.na(dat$method), , drop = FALSE]
    dat$subcategory <- label
    pieces[[i]] <- dat
  }
  pieces <- purrr::compact(pieces)

  if (length(pieces) == 0) {
    cli::cli_inform("Cannot plot: no successful MR results available.")
    return(invisible(NULL))
  }

  combined <- dplyr::bind_rows(pieces)
  if (length(relabel) > 0) {
    combined$method <- dplyr::recode(combined$method, !!!relabel)
  }

  TwoSampleMR::forest_plot_1_to_many(
    mr_res = combined,
    b = "b",
    se = "se",
    TraitM = "method",
    by = if (sectioned) "subcategory" else NULL,
    ao_slc = FALSE,
    exponentiate = exponentiate,
    trans = trans %||% if (exponentiate) "log2" else "identity",
    xlab = xlab %||% if (exponentiate) "OR (95% CI)" else "Effect (95% CI)",
    lo = if (!is.null(xlim)) xlim[1] else NULL,
    up = if (!is.null(xlim)) xlim[2] else NULL,
    shape_points = shape_points,
    colour_scheme = colour_scheme,
    col_text_size = col_text_size,
    subheading_size = subheading_size,
    col1_width = col1_width,
    col1_title = "Method",
    ...
  )
}

#' Forest plot of MR outcomes, optionally coloured/shaped by a grouping column
#'
#' One row per outcome per requested method. `method` is a plain inclusion
#' filter, not a priority/fallback selector: request one model to get one
#' row per outcome, or both to get two rows dodged apart on the same outcome
#' wherever both were computed. `colour_by`/`shape_by` map any column in
#' `mr_res` (already computed by the caller) to either aesthetic
#' independently, so e.g. instrument and model can be shown together without
#' colliding. Requires the `ggplot2` package.
#'
#' @param mr_res Data frame with columns `exposure`, `outcome`, `method`,
#'   `or`, `or_lci95`, `or_uci95`, `subcategory` -- typically built by the
#'   caller by binding together the `$results` of several `mr_result`
#'   objects and attaching a `subcategory` (and any grouping columns) with
#'   `dplyr::mutate()`. See examples.
#' @param xlab X-axis label.
#' @param method Character vector of `method` values to include. Default
#'   `c("Inverse variance weighted", "IVW (fixed effects)")` includes both
#'   wherever computed. Outcomes lacking a requested method simply
#'   contribute no row for it -- there is no silent substitution.
#' @param relabel Named character vector, `c(old = new)`, applied to the
#'   `method` column for display, so it reads sensibly if `colour_by` or
#'   `shape_by` is `"method"`. Default relabels the two IVW variants to
#'   `"Random effects"` / `"Fixed effects"`.
#' @param colour_by,shape_by Column name in `mr_res` to map to colour/shape
#'   respectively, or `NULL` for none.
#' @param colour_values,shape_values Named vectors passed to
#'   `ggplot2::scale_colour_manual()`/`ggplot2::scale_shape_manual()`, or
#'   `NULL` to use ggplot2's defaults.
#' @param section_order Character vector giving the `subcategory` display
#'   order, or `NULL` to use first-appearance order.
#' @param dodge_width Vertical dodge width applied when two rows share an
#'   outcome (e.g. fixed and random effects for the same outcome).
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' primary <- run_mr(
#'   exposure = cd40_exposure, exposure_id = "CD40",
#'   outcome = sjogren_outcome, outcome_id = "SjD",
#'   instrument_region = list(chromosome = "20", start = 44746911, end = 44758502),
#'   bfile = system.file("extdata", "ld_ref", package = "mrpipeline"),
#'   methods = c("ivw", "ivw_fe")
#' )
#'
#' mr_res <- dplyr::bind_rows(
#'   primary$results |> dplyr::mutate(subcategory = "Primary"),
#'   positive_control_result$results |> dplyr::mutate(subcategory = "Positive control"),
#'   negative_control_result$results |> dplyr::mutate(subcategory = "Negative control")
#' ) |>
#'   dplyr::mutate(
#'     instrument = dplyr::if_else(exposure == "CD40", "GWS instrument", "Functional instrument")
#'   )
#'
#' # Both models, coloured by instrument, shaped by model
#' outcome_forest_plot(
#'   mr_res,
#'   xlab = "OR (95% CI)",
#'   colour_by = "instrument",
#'   shape_by = "method",
#'   colour_values = c("GWS instrument" = "#3c6e5e", "Functional instrument" = "#b5541f"),
#'   shape_values = c("Random effects" = 16, "Fixed effects" = 17),
#'   section_order = c("Primary", "Positive control", "Negative control")
#' )
#' }
#'
#' @export
outcome_forest_plot <- function(
  mr_res,
  xlab,
  method = c("Inverse variance weighted", "IVW (fixed effects)"),
  relabel = c(
    "Inverse variance weighted" = "Random effects",
    "IVW (fixed effects)" = "Fixed effects"
  ),
  colour_by = NULL,
  shape_by = NULL,
  colour_values = NULL,
  shape_values = NULL,
  section_order = NULL,
  dodge_width = 0.5
) {
  rlang::check_installed("ggplot2", reason = "to plot MR results.")

  required_cols <- c(
    "exposure",
    "outcome",
    "method",
    "or",
    "or_lci95",
    "or_uci95",
    "subcategory"
  )
  missing_cols <- setdiff(required_cols, names(mr_res))
  if (length(missing_cols) > 0) {
    cli::cli_abort("{.arg mr_res} is missing column{?s}: {missing_cols}.")
  }

  # `method` is a plain inclusion filter, not a priority/fallback selector:
  # request one model to get one row per outcome, or both to get two rows
  # (dodged apart) wherever both were computed. Outcomes lacking a requested
  # model simply don't contribute a row for it -- no silent substitution.
  dat <- mr_res |>
    dplyr::filter(.data$method %in% .env$method)

  if (nrow(dat) == 0) {
    cli::cli_abort("No rows match any of {.val {method}} in {.arg mr_res}.")
  }

  if (length(relabel) > 0) {
    dat$method <- dplyr::recode(dat$method, !!!relabel)
  }

  dat$subcategory <- if (!is.null(section_order)) {
    factor(dat$subcategory, levels = section_order)
  } else {
    factor(dat$subcategory, levels = unique(dat$subcategory))
  }
  dat$outcome <- factor(dat$outcome, levels = rev(unique(dat$outcome)))

  # Always group by method so that two rows sharing the same outcome (one
  # random-effects, one fixed-effects) dodge apart instead of overplotting,
  # regardless of whether `method` is also mapped to colour/shape.
  aes_list <- list(
    x = rlang::expr(.data$or),
    y = rlang::expr(.data$outcome),
    group = rlang::expr(.data$method)
  )
  if (!is.null(colour_by)) {
    aes_list$colour <- rlang::expr(.data[[!!colour_by]])
  }
  if (!is.null(shape_by)) {
    aes_list$shape <- rlang::expr(.data[[!!shape_by]])
  }
  mapping <- rlang::inject(ggplot2::aes(!!!aes_list))

  dodge <- ggplot2::position_dodge(width = dodge_width)

  p <- ggplot2::ggplot(dat, mapping) +
    ggplot2::geom_vline(
      xintercept = 1,
      linetype = "dashed",
      colour = "grey45",
      linewidth = 0.4
    ) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = .data$or_lci95, xmax = .data$or_uci95),
      height = 0,
      linewidth = 0.7,
      position = dodge
    ) +
    ggplot2::geom_point(size = 3, position = dodge) +
    ggplot2::scale_x_continuous(trans = "log2") +
    ggplot2::facet_grid(
      rows = ggplot2::vars(subcategory),
      scales = "free_y",
      space = "free_y",
      switch = "y"
    ) +
    ggplot2::labs(x = xlab, y = NULL, colour = colour_by, shape = shape_by) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      strip.placement = "outside",
      strip.background = ggplot2::element_rect(fill = "grey93", colour = NA),
      strip.text.y.left = ggplot2::element_text(
        angle = 0,
        hjust = 0,
        face = "bold",
        size = 11
      ),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.spacing.y = grid::unit(0.4, "lines"),
      legend.position = "top",
      legend.title = ggplot2::element_text(face = "bold")
    )

  if (!is.null(colour_values)) {
    p <- p + ggplot2::scale_colour_manual(values = colour_values)
  }
  if (!is.null(shape_values)) {
    p <- p + ggplot2::scale_shape_manual(values = shape_values)
  }

  p
}
