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
#' @return A list of ggplot objects, one per exposure-outcome pair, as
#'   returned by TwoSampleMR's plotting functions -- for a single analysis,
#'   take `[[1]]`. Returns `NULL` invisibly if the result status is not
#'   `"success"`, if no results are available, or (for `type = "loo"`) if
#'   `x$loo` was not computed.
#'
#' @examples
#' \dontrun{
#' bfile <- sub("\\.bed$", "", system.file("extdata", "ld_ref.bed", package = "mrpipeline"))
#' result <- run_mr(
#'   exposure = cd40_exposure, exposure_id = "CD40",
#'   outcome = sjogren_outcome, outcome_id = "SjD",
#'   instrument_region = list(chromosome = "20", start = 44746911, end = 44758502),
#'   rsq_thresh = 0.3,
#'   bfile = bfile,
#'   methods = c("ivw_random", "egger", "weighted_median", "loo")
#' )
#' plot(result, type = "scatter")
#' plot(result, type = "forest")
#' plot(result, type = "funnel")
#' plot(result, type = "loo")
#' }
#'
#' @export
plot.mr_result <- function(
  x,
  type = c("scatter", "forest", "funnel", "loo"),
  ...
) {
  rlang::check_installed("ggplot2", reason = "to plot MR results.")
  type <- match.arg(type)

  if (x$status != "success" || nrow(x$results) == 0) {
    cli::cli_inform("Cannot plot: no successful MR results available.")
    return(invisible(NULL))
  }

  if (type == "loo") {
    if (is.null(x$loo) || nrow(x$loo) == 0) {
      cli::cli_inform(
        "Cannot plot: leave-one-out results not available (add \"loo\" to methods)."
      )
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
#'   ABF posterior probabilities (H0--H4), or `"regional"` for exposure and
#'   outcome regional association plots, coloured by LD when r2 is available
#'   and with an optional gene track (see [plot_coloc_regional()]).
#'   `"locuszoom"` is deprecated in favour of `"regional"` with `ens_db`.
#' @param ... For `type = "regional"`, forwarded to [plot_coloc_regional()]
#'   (`ens_db`, `bfile`, `plink_bin` and `index_snp`). For `type = "locuszoom"`,
#'   forwarded to [plot_coloc_locuszoom()]. Ignored for `"pp_bar"`.
#'
#' @return For `type = "pp_bar"`/`"regional"`, a ggplot object (or `NULL`
#'   invisibly if the result status is not `"success"` or the required data
#'   is unavailable). For `type = "locuszoom"`, `NULL` invisibly always --
#'   see [plot_coloc_locuszoom()], which draws directly to the current
#'   graphics device rather than returning a re-renderable object.
#'
#' @examples
#' \dontrun{
#' bfile <- sub("\\.bed$", "", system.file("extdata", "ld_ref.bed", package = "mrpipeline"))
#' result <- run_coloc(
#'   exposure = cd40_exposure,
#'   outcome = sjogren_outcome,
#'   gene_chr = "20", gene_start = 44746911, gene_end = 44758502,
#'   bfile = bfile
#' )
#' plot(result, type = "pp_bar")
#' plot(result, type = "regional")
#'
#' # Colour by LD with another SNP (recomputed, so it needs the panel)
#' plot(result, type = "regional", index_snp = "rs1883832", bfile = bfile)
#'
#' # Add a gene track. A character ens_db must be attached, not just installed.
#' library(EnsDb.Hsapiens.v75)
#' plot(result, type = "regional", ens_db = "EnsDb.Hsapiens.v75")
#' }
#'
#' @export
plot.coloc_result <- function(
  x,
  type = c("pp_bar", "regional", "locuszoom"),
  ...
) {
  rlang::check_installed("ggplot2", reason = "to plot coloc results.")
  type <- match.arg(type)

  if (x$status != "success") {
    cli::cli_inform("Cannot plot: no successful coloc results available.")
    return(invisible(NULL))
  }

  switch(
    type,
    pp_bar = plot_coloc_pp_bar(x),
    regional = plot_coloc_regional(x, ...),
    locuszoom = plot_coloc_locuszoom(x, ...)
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
#'
#' Plots -log10(p) against position, with exposure and outcome panels. Points
#' are coloured by LD (r2) with an index SNP, using the r2 stored by
#' [run_coloc()], so no PLINK or panel is needed at plot time. With no stored
#' r2 (a [run_coloc()] call with no `bfile`, or an object saved before
#' mrpipeline 0.10.0), passing `bfile` computes it for this plot, and otherwise
#' the plot is uncoloured, with any recorded index SNP still marked.
#'
#' The same index SNP and r2 are used in both panels, so the outcome panel shows
#' LD with the exposure lead by default, not with the outcome lead.
#'
#' With `ens_db`, a gene track from `locuszoomr` is added under the panels, the
#' x axis is in Mb rather than bp, and the plots are composed with
#' `cowplot::plot_grid()`. Adding `ggplot2::labs(title = ...)` still sets the
#' title, but other `ggplot2` additions, such as `theme()` changes to the
#' panels, do not reach a composed plot.
#'
#' @param x A `coloc_result` object.
#' @param ens_db `NULL` (the default) for no gene track, or an
#'   `ensembldb::EnsDb` object, or the name of an Ensembl annotation package as
#'   a string (e.g. `"EnsDb.Hsapiens.v75"`), matching the genome build of `x`'s
#'   data. If given as a string, the package must be **attached** first
#'   (`library(EnsDb.Hsapiens.v75)`), not merely installed.
#' @param bfile Path to PLINK bfile prefix (without .bed/.bim/.fam), the local
#'   LD reference panel. Only used when r2 has to be computed: `x` has no stored
#'   r2, or `index_snp` differs from the stored index SNP. SNPs absent from the
#'   panel stay on the plot in grey.
#' @param plink_bin Path to PLINK binary. If `NULL`, auto-detected via
#'   [genetics.binaRies::get_plink_binary()].
#' @param index_snp rsID to colour LD against, overriding any index SNP stored
#'   by [run_coloc()] for this plot. Defaults to the stored index SNP, or the
#'   exposure lead when there is none: the lowest `pval.exposure`, with ties
#'   going to the largest |beta/se|. Must be one of the SNPs in `x`.
#' @param ... Forwarded to [compute_ld_to_index()] (`plink_threads`,
#'   `plink_memory`).
#'
#' @return A ggplot object, or `NULL` invisibly if there is no harmonised data.
#' @keywords internal
plot_coloc_regional <- function(
  x,
  ens_db = NULL,
  bfile = NULL,
  plink_bin = NULL,
  index_snp = NULL,
  ...
) {
  dat <- x$harmonised_data
  if (is.null(dat) || nrow(dat) == 0) {
    cli::cli_inform("No harmonised data available for regional plot.")
    return(invisible(NULL))
  }

  if (!is.null(ens_db)) {
    rlang::check_installed(
      c("locuszoomr", "ensembldb", "cowplot"),
      reason = "to add a gene track to regional plots."
    )
  }

  ld <- regional_ld(x, index_snp, bfile, plink_bin, ...)

  # Build long-format data for faceted plot
  pos_unit <- if (is.null(ens_db)) 1 else 1e6
  exposure_dat <- data.frame(
    pos = dat$pos.exposure / pos_unit,
    log10p = -log10(dat$pval.exposure),
    panel = "Exposure"
  )
  outcome_dat <- data.frame(
    pos = dat$pos.outcome / pos_unit,
    log10p = -log10(dat$pval.outcome),
    panel = "Outcome"
  )
  plot_data <- rbind(exposure_dat, outcome_dat)

  title <- "Regional association"
  subtitle <- NULL
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data$pos, y = .data$log10p)
  )

  if (is.null(ld)) {
    p <- p + ggplot2::geom_point(alpha = 0.6)
  } else {
    is_index <- rep(dat$SNP == ld$index_snp, 2)
    points <- plot_data[!is_index, ]
    subtitle <- paste("Purple diamond: index SNP", ld$index_snp)
    if (is.null(ld$r2)) {
      p <- p + ggplot2::geom_point(data = points, alpha = 0.6)
    } else {
      points$r2_bin <- rep(bin_ld_r2(ld$r2), 2)[!is_index]
      p <- p +
        # Sorted by bin so that strong LD is drawn on top
        ggplot2::geom_point(
          data = points[order(points$r2_bin), ],
          ggplot2::aes(colour = .data$r2_bin),
          show.legend = TRUE
        ) +
        ggplot2::scale_colour_manual(
          name = expression(r^2),
          values = ld_r2_colours(),
          drop = FALSE,
          guide = ggplot2::guide_legend(reverse = TRUE)
        )
    }
    p <- p +
      ggplot2::geom_point(
        data = plot_data[is_index, ],
        shape = 23,
        size = 3,
        colour = "black",
        fill = "#9632B8"
      )
  }

  p <- p +
    ggplot2::facet_wrap(~panel, ncol = 1, scales = "free_y") +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Position (bp)",
      y = expression(-log[10](p))
    ) +
    ggplot2::theme_minimal()

  if (is.null(ens_db)) {
    return(p)
  }
  add_gene_track(p, dat, ens_db, title, subtitle)
}

#' Stack a gene track under a regional plot
#'
#' `locuszoomr::genetracks_grob()` fills its panel with an x scale of the locus
#' range plus 4% each side, so both plots use exactly that axis (see
#' [gene_track_x_scale()]). With ggplot2's default 5% expansion the genes would
#' sit slightly off the points.
#'
#' @param p The regional plot, with positions in Mb.
#' @param dat `harmonised_data` of the `coloc_result`.
#' @param ens_db As for [plot_coloc_regional()].
#' @param title,subtitle Moved from `p` to the composed plot, where a later
#'   `labs(title = ...)` replaces them.
#'
#' @return A composed ggplot object.
#' @keywords internal
add_gene_track <- function(p, dat, ens_db, title, subtitle) {
  # Only the genes are taken from the locus: no data goes in, so the plotted
  # points and r2 are exactly those of the plot without a track
  loc <- suppressMessages(locuszoomr::locus(
    seqname = dat$chr.exposure[1],
    xrange = range(dat$pos.exposure, dat$pos.outcome),
    ens_db = ens_db
  ))
  x_scale <- gene_track_x_scale(loc$xrange)

  panels <- p +
    x_scale +
    ggplot2::labs(title = NULL, subtitle = NULL, x = NULL) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )
  track <- locuszoomr::gg_genetracks(loc) + x_scale

  # plot_grid() returns a ggplot whose theme blanks the title, so labs() would
  # silently draw nothing unless the outer plot has title elements of its own
  cowplot::plot_grid(
    panels,
    track,
    ncol = 1,
    rel_heights = c(3, 1),
    align = "v",
    axis = "lr"
  ) +
    ggplot2::labs(title = title, subtitle = subtitle) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0, size = 13),
      plot.subtitle = ggplot2::element_text(hjust = 0, size = 10),
      plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 5.5)
    )
}

#' Shared x axis for a regional plot and its gene track
#'
#' @param xrange Plotted region in bp, as `c(start, end)`.
#'
#' @return A ggplot2 scale, in Mb.
#' @keywords internal
gene_track_x_scale <- function(xrange) {
  ggplot2::scale_x_continuous(
    limits = xrange / 1e6,
    expand = ggplot2::expansion(mult = 0.04)
  )
}

#' r2 and index SNP for a regional plot
#'
#' The index SNP is `index_snp`, or else the one stored by [run_coloc()], or
#' else the exposure lead. Stored r2 is used when it is relative to that SNP,
#' and otherwise r2 is computed from `bfile`.
#'
#' @inheritParams plot_coloc_regional
#'
#' @return `NULL` for no LD colouring, or a list with `index_snp` and `r2`,
#'   aligned to the rows of `harmonised_data`. `r2` is `NULL` when an index SNP
#'   is recorded but there is no r2 for it and no `bfile`.
#' @keywords internal
regional_ld <- function(
  x,
  index_snp = NULL,
  bfile = NULL,
  plink_bin = NULL,
  ...
) {
  dat <- x$harmonised_data
  stored_index <- x$ld_index_snp

  if (!is.null(index_snp) && !index_snp %in% dat$SNP) {
    cli::cli_abort(
      "Index SNP {.val {index_snp}} is not among the SNPs in this result."
    )
  }

  index <- index_snp %||% stored_index
  is_stored <- !is.null(index) && identical(index, stored_index)

  if (is_stored && "ld_r2_index" %in% names(dat)) {
    return(list(r2 = dat$ld_r2_index, index_snp = index))
  }

  if (is.null(bfile)) {
    if (!is.null(index) && !is_stored) {
      cli::cli_abort(c(
        "{.arg bfile} is required to colour by LD with {.val {index}}.",
        "i" = "Any other SNP than a stored index SNP needs r2 computed from an LD reference panel."
      ))
    }
    cli::cli_inform(
      "No stored r2, so the regional plot has no LD colouring. Pass {.arg bfile} to compute it."
    )
    if (is_stored) {
      return(list(r2 = NULL, index_snp = index))
    }
    return(NULL)
  }

  index <- index %||% exposure_lead(dat)
  ld <- compute_ld_to_index(dat$SNP, index, bfile, plink_bin = plink_bin, ...)
  list(r2 = ld$r2[match(dat$SNP, ld$SNP)], index_snp = index)
}

#' Bin r2 for LD colouring
#'
#' @param r2 Numeric vector of r2 values, `NA` where unknown.
#'
#' @return A factor with the LocusZoom bins and a level for missing r2, ordered
#'   from no data up to strongest LD.
#' @keywords internal
bin_ld_r2 <- function(r2) {
  levels <- names(ld_r2_colours())
  bins <- cut(
    pmin(r2, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = levels[-1],
    include.lowest = TRUE
  )
  bins <- as.character(bins)
  bins[is.na(bins)] <- levels[1]
  factor(bins, levels = levels)
}

#' LocusZoom colours for r2 bins
#' @return Named character vector, from no data up to strongest LD.
#' @keywords internal
ld_r2_colours <- function() {
  c(
    "No LD data" = "#B8B8B8",
    "0-0.2" = "#357EBD",
    "0.2-0.4" = "#46B8DA",
    "0.4-0.6" = "#5CB85C",
    "0.6-0.8" = "#EEA236",
    "0.8-1" = "#D43F3A"
  )
}

#' LD-coloured regional plots for coloc results, via locuszoomr
#'
#' Deprecated: use `plot(x, type = "regional", ens_db = ...)`, which returns a
#' ggplot. `plot_coloc_locuszoom()` keeps working for now, and warns.
#'
#' Builds exposure and outcome regional plots using `locuszoomr::locus()`,
#' coloured by LD (r2) to a chosen index SNP, and stacks them with
#' `locuszoomr::multi_layout()`. LD is computed from a local reference panel
#' via [compute_ld_to_index()] (reusing [compute_ld_matrix()] -- no network
#' lookups). Without `bfile`, a warning says how to add LD colouring and the
#' plots are drawn uncoloured.
#'
#' SNPs absent from the panel are kept, with no LD value. That happens when
#' `x` came from [run_coloc()] with no `bfile`, or was aligned to another panel.
#'
#' There is no genome-build field stored on a `coloc_result` -- `ens_db` must
#' match whatever build the exposure/outcome data (and `bfile`) actually use
#' (e.g. `EnsDb.Hsapiens.v75` for GRCh37, `EnsDb.Hsapiens.v86` for GRCh38).
#'
#' @param x A `coloc_result` object.
#' @param ens_db Either an `ensembldb::EnsDb` object, or the name of an
#'   Ensembl annotation package as a string (e.g. `"EnsDb.Hsapiens.v75"`),
#'   matching the genome build of `x`'s data. If given as a string, the
#'   package must be **attached** first (`library(EnsDb.Hsapiens.v75)`), not
#'   merely installed -- `locuszoomr::locus()` resolves a character `ens_db`
#'   via `get(ens_db)` on the search path, and errors
#'   ("Ensembl database not loaded") if it isn't attached.
#' @param bfile Path to PLINK bfile prefix (without .bed/.bim/.fam), used as
#'   the local LD reference panel, or `NULL` (the default) for no LD colouring.
#' @param plink_bin Path to PLINK binary. If `NULL`, auto-detected via
#'   [genetics.binaRies::get_plink_binary()].
#' @param index_snp rsID to colour LD against and centre the plotted region
#'   on. Defaults to the SNP with the lowest `pval.exposure`, matching
#'   `locuszoomr::locus()`'s own default when no index SNP is given.
#' @param ... Forwarded to [compute_ld_to_index()] (`plink_threads`,
#'   `plink_memory`).
#'
#' @return `NULL`, invisibly. Like base/grid plotting functions (and unlike
#'   this package's other `plot_coloc_*()` functions, which return a ggplot
#'   object), this draws directly to the current graphics device as a side
#'   effect -- open a device (`pdf()`, `png()`, ...) before calling, and
#'   `dev.off()` after, to save the result.
#'
#' @keywords internal
plot_coloc_locuszoom <- function(
  x,
  ens_db,
  bfile = NULL,
  plink_bin = NULL,
  index_snp = NULL,
  ...
) {
  rlang::check_installed(
    c("locuszoomr", "ensembldb"),
    reason = "to build locuszoom regional plots."
  )

  cli::cli_warn(c(
    "{.code type = \"locuszoom\"} is deprecated.",
    "i" = "Use {.code plot(x, type = \"regional\", ens_db = ens_db)} instead, which returns a ggplot."
  ))

  dat <- x$harmonised_data
  if (is.null(dat) || nrow(dat) == 0) {
    cli::cli_inform("No harmonised data available for locuszoom plot.")
    return(invisible(NULL))
  }

  if (is.null(index_snp)) {
    index_snp <- dat$SNP[which.min(dat$pval.exposure)]
  }

  ld <- NULL
  if (is.null(bfile)) {
    cli::cli_warn(c(
      "No {.arg bfile} supplied, so the locuszoom plots have no LD colouring.",
      "i" = paste0(
        "To colour by LD, pass an LD reference panel: ",
        "{.code plot(x, type = \"locuszoom\", ens_db = ens_db, bfile = \"/path/to/ld_ref\")}."
      ),
      "i" = "See \"Setting Up an LD Reference Panel\" in {.vignette mrpipeline::mrpipeline-user-guide}."
    ))
  } else {
    ld <- compute_ld_to_index(
      dat$SNP,
      index_snp,
      bfile,
      plink_bin = plink_bin,
      ...
    )
  }

  exposure_dat <- dat[, c(
    "SNP",
    "chr.exposure",
    "pos.exposure",
    "pval.exposure"
  )]
  names(exposure_dat) <- c("SNP", "chrom", "pos", "pval")

  outcome_dat <- dat[, c("SNP", "chr.outcome", "pos.outcome", "pval.outcome")]
  names(outcome_dat) <- c("SNP", "chrom", "pos", "pval")

  ld_col <- NULL
  if (!is.null(ld)) {
    exposure_dat <- merge(exposure_dat, ld, by = "SNP", all.x = TRUE)
    outcome_dat <- merge(outcome_dat, ld, by = "SNP", all.x = TRUE)
    ld_col <- "r2"
  }

  loc_exp <- locuszoomr::locus(
    data = exposure_dat,
    ens_db = ens_db,
    chrom = "chrom",
    pos = "pos",
    p = "pval",
    labs = "SNP",
    LD = ld_col,
    index_snp = index_snp
  )
  loc_out <- locuszoomr::locus(
    data = outcome_dat,
    ens_db = ens_db,
    chrom = "chrom",
    pos = "pos",
    p = "pval",
    labs = "SNP",
    LD = ld_col,
    index_snp = index_snp
  )

  # multi_layout()'s nrow/ncol both default to 1 -- a 1x1 grid silently
  # drops the second plot rather than erroring, so nrow must be set
  # explicitly to actually stack both panels. Like base/grid plotting
  # functions, multi_layout() draws directly to the current graphics device
  # as a side effect and returns NULL -- it is not a re-renderable object
  # like a ggplot; wrap in pdf()/png()/dev.off() (or similar) to save it.
  invisible(locuszoomr::multi_layout(
    list(loc_exp, loc_out),
    nrow = 2,
    ncol = 1
  ))
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
#'   effects)", "IVW (random effects)", "MR Egger", "Weighted median")`.
#'   Labels no longer encode LD correction (that lives in the `ld_corrected`
#'   column), so LD-corrected results match the same defaults.
#' @param relabel Named character vector, `c(old = new)`, applied to the
#'   `method` column for display only, after row filtering/ordering -- so
#'   matching against `results$results$method` is unaffected. Empty by
#'   default.
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
#' bfile <- sub("\\.bed$", "", system.file("extdata", "ld_ref.bed", package = "mrpipeline"))
#' result <- run_mr(
#'   exposure = cd40_exposure, exposure_id = "CD40",
#'   outcome = sjogren_outcome, outcome_id = "SjD",
#'   instrument_region = list(chromosome = "20", start = 44746911, end = 44758502),
#'   rsq_thresh = 0.3,
#'   bfile = bfile,
#'   methods = c("ivw_random", "ivw_fixed", "egger", "weighted_median")
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
#'
#' # Uncorrected vs LD-corrected arms of the same analysis
#' corrected <- run_mr(
#'   exposure = cd40_exposure, exposure_id = "CD40",
#'   outcome = sjogren_outcome, outcome_id = "SjD",
#'   instrument_region = list(chromosome = "20", start = 44746911, end = 44758502),
#'   rsq_thresh = 0.3,
#'   bfile = bfile, ld_correct = TRUE,
#'   methods = c("ivw_random", "ivw_fixed", "egger", "weighted_median")
#' )
#' forest_plot(list("Uncorrected" = result, "LD-corrected" = corrected))
#' }
#'
#' @export
forest_plot <- function(
  results,
  methods = c(
    "IVW (fixed effects)",
    "IVW (random effects)",
    "MR Egger",
    "Weighted median"
  ),
  relabel = character(0),
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
    dat <- res$results |> dplyr::filter(.data$method %in% methods)
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
#'   `subcategory`, plus whichever value columns `scale` requires (see
#'   below) -- typically built by the caller by binding together the
#'   `$results` of several `mr_result` objects and attaching a `subcategory`
#'   (and any grouping columns) with `dplyr::mutate()`. See examples.
#' @param xlab X-axis label.
#' @param method Character vector of `method` values to include. Default
#'   `c("IVW (random effects)", "IVW (fixed effects)")` includes both
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
#' @param scale `"or"` (default) plots an odds-ratio scale: `or_col`/
#'   `or_lci_col`/`or_uci_col` on a log2 x-axis with a reference line at 1.
#'   `"beta"` plots a linear effect-size scale: `b_col`/`se_col` (CI derived
#'   internally as `b \eqn{\pm} 1.96*se`) on a linear x-axis with a
#'   reference line at 0. Required columns depend on this choice -- only the
#'   3 OR columns are required for `scale="or"`, only the 2 beta/SE columns
#'   for `scale="beta"`.
#' @param or_col,or_lci_col,or_uci_col Column names to use for the OR point
#'   estimate and 95% CI when `scale="or"`. Defaults (`"or"`/`"or_lci95"`/
#'   `"or_uci95"`) match every existing call site -- override to plot a
#'   rescaled column (e.g. `or_col="or_scaled"`) without needing to rename
#'   columns in the caller.
#' @param b_col,se_col Column names to use for the point estimate and
#'   standard error when `scale="beta"`. Defaults `"b"`/`"se"` -- override to
#'   plot a rescaled column (e.g. `b_col="b_scaled", se_col="se_scaled"`).
#' @param group_by Optional column name in `mr_res` giving an outer grouping
#'   above `subcategory` (e.g. "Primary" vs "Sensitivity Analysis"). When
#'   supplied, facets on `vars(.data[[group_by]], subcategory)` instead of
#'   just `vars(subcategory)` -- `ggplot2::facet_grid()`'s native handling of
#'   multiple row-facetting variables spans the outer variable's strip across
#'   its contiguous inner rows, giving a nested header with no extra
#'   dependency. `NULL` (default) reproduces today's single-facet behaviour
#'   exactly.
#' @param group_order Character vector giving `group_by`'s display order
#'   (parallel to `section_order` for `subcategory`), or `NULL` for
#'   first-appearance order. Ignored if `group_by` is `NULL`.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' bfile <- sub("\\.bed$", "", system.file("extdata", "ld_ref.bed", package = "mrpipeline"))
#' primary <- run_mr(
#'   exposure = cd40_exposure, exposure_id = "CD40",
#'   outcome = sjogren_outcome, outcome_id = "SjD",
#'   instrument_region = list(chromosome = "20", start = 44746911, end = 44758502),
#'   rsq_thresh = 0.3,
#'   bfile = bfile,
#'   methods = c("ivw_random", "ivw_fixed")
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
#'
#' # Linear beta/SE scale on a rescaled column, nested Primary/Sensitivity header
#' outcome_forest_plot(
#'   mr_res,
#'   xlab = "Beta (95% CI)",
#'   scale = "beta",
#'   b_col = "b_scaled", se_col = "se_scaled",
#'   group_by = "tier_group", group_order = c("Primary", "Sensitivity Analysis"),
#'   section_order = c("Cohort A", "Cohort B")
#' )
#' }
#'
#' @export
outcome_forest_plot <- function(
  mr_res,
  xlab,
  method = c("IVW (random effects)", "IVW (fixed effects)"),
  relabel = c(
    "IVW (random effects)" = "Random effects",
    "IVW (fixed effects)" = "Fixed effects"
  ),
  colour_by = NULL,
  shape_by = NULL,
  colour_values = NULL,
  shape_values = NULL,
  section_order = NULL,
  dodge_width = 0.5,
  scale = c("or", "beta"),
  or_col = "or",
  or_lci_col = "or_lci95",
  or_uci_col = "or_uci95",
  b_col = "b",
  se_col = "se",
  group_by = NULL,
  group_order = NULL
) {
  rlang::check_installed("ggplot2", reason = "to plot MR results.")
  scale <- match.arg(scale)

  value_cols <- if (scale == "or") {
    c(or_col, or_lci_col, or_uci_col)
  } else {
    c(b_col, se_col)
  }
  required_cols <- c("exposure", "outcome", "method", "subcategory", value_cols)
  if (!is.null(group_by)) {
    required_cols <- c(required_cols, group_by)
  }
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

  if (!is.null(group_by)) {
    dat[[group_by]] <- if (!is.null(group_order)) {
      factor(dat[[group_by]], levels = group_order)
    } else {
      factor(dat[[group_by]], levels = unique(dat[[group_by]]))
    }
  }

  # scale-dependent x/xmin/xmax, trans, and reference line -- everything
  # downstream (aes, geoms, facets, theme) is otherwise identical between
  # the two scales.
  if (scale == "or") {
    dat$.x <- dat[[or_col]]
    dat$.xmin <- dat[[or_lci_col]]
    dat$.xmax <- dat[[or_uci_col]]
    x_trans <- "log2"
    vline_at <- 1
  } else {
    dat$.x <- dat[[b_col]]
    dat$.xmin <- dat[[b_col]] - 1.96 * dat[[se_col]]
    dat$.xmax <- dat[[b_col]] + 1.96 * dat[[se_col]]
    x_trans <- "identity"
    vline_at <- 0
  }

  # Always group by method so that two rows sharing the same outcome (one
  # random-effects, one fixed-effects) dodge apart instead of overplotting,
  # regardless of whether `method` is also mapped to colour/shape.
  aes_list <- list(
    x = rlang::expr(.data$.x),
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

  # facet_grid()'s native handling of multiple row-facetting variables spans
  # the outer (group_by) variable's strip across its contiguous inner
  # (subcategory) rows -- a nested header with no extra dependency (e.g.
  # ggh4x) needed.
  facet_rows <- if (!is.null(group_by)) {
    ggplot2::vars(.data[[group_by]], .data$subcategory)
  } else {
    ggplot2::vars(.data$subcategory)
  }

  p <- ggplot2::ggplot(dat, mapping) +
    ggplot2::geom_vline(
      xintercept = vline_at,
      linetype = "dashed",
      colour = "grey45",
      linewidth = 0.4
    ) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = .data$.xmin, xmax = .data$.xmax),
      height = 0,
      linewidth = 0.7,
      position = dodge
    ) +
    ggplot2::geom_point(size = 3, position = dodge) +
    ggplot2::scale_x_continuous(trans = x_trans) +
    ggplot2::facet_grid(
      rows = facet_rows,
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

#' Table-style forest plot (forestplot package backend)
#'
#' One row per outcome, with an inline result/CI/p-value table alongside the
#' plotted estimate and confidence interval. Requires the `forestplot`
#' package.
#'
#' @param dat Data frame with columns `exposure`, `outcome`, `n_snps`,
#'   `estimate`, `lower`, `upper`, `p_value`.
#' @param null_value Reference line value (1 for OR scale, 0 for beta scale).
#' @param xlab X-axis label.
#' @param box_colour Point/line colour.
#'
#' @return A forestplot object (render with `plot()`/`print()`).
#'
#' @examples
#' \dontrun{
#' dat <- data.frame(
#'   exposure = "Genetically-proxied NLRP3 inhibition",
#'   outcome = c("CHD", "Stroke", "Type 2 diabetes"),
#'   n_snps = c(8, 8, 8),
#'   estimate = c(0.92, 0.95, 1.03),
#'   lower = c(0.85, 0.88, 0.94),
#'   upper = c(0.99, 1.02, 1.13),
#'   p_value = c(0.02, 0.15, 0.5)
#' )
#' table_forest_plot(dat, null_value = 1, xlab = "OR (95% CI)")
#' }
#'
#' @export
table_forest_plot <- function(dat, null_value, xlab, box_colour = "#35568a") {
  rlang::check_installed(
    "forestplot",
    reason = "to build table-style forest plots."
  )

  required_cols <- c(
    "exposure",
    "outcome",
    "n_snps",
    "estimate",
    "lower",
    "upper",
    "p_value"
  )
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0) {
    cli::cli_abort("{.arg dat} is missing column{?s}: {missing_cols}.")
  }

  result <- sprintf("%.3f [%.3f; %.3f]", dat$estimate, dat$lower, dat$upper)
  p_text <- ifelse(
    dat$p_value < 0.001,
    format(dat$p_value, scientific = TRUE, digits = 2),
    sprintf("%.3f", dat$p_value)
  )

  table_text <- cbind(
    c("Exposure", dat$exposure),
    c("Outcome", as.character(dat$outcome)),
    c("# SNPs", dat$n_snps),
    c("Estimate [95% CI]", result),
    c("P value", p_text)
  )

  forestplot::forestplot(
    labeltext = table_text,
    mean = c(NA, dat$estimate),
    lower = c(NA, dat$lower),
    upper = c(NA, dat$upper),
    new_page = FALSE,
    zero = null_value,
    xlog = FALSE,
    xlab = xlab,
    graph.pos = 3,
    graphwidth = grid::unit(35, "mm"),
    boxsize = 0.14,
    lwd.ci = 1,
    ci.vertices = TRUE,
    fn.ci_norm = forestplot::fpDrawCircleCI,
    col = forestplot::fpColors(
      box = box_colour,
      lines = box_colour,
      zero = "grey50"
    ),
    txt_gp = forestplot::fpTxtGp(
      label = grid::gpar(cex = 0.8),
      ticks = grid::gpar(cex = 0.8),
      xlab = grid::gpar(cex = 0.8)
    )
  )
}
