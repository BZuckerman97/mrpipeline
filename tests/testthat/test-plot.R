# Tests for plot methods

skip_if_not_installed("ggplot2")

# --- mr_result plots ---------------------------------------------------------

test_that("plot.mr_result returns NULL for non-success status", {
  res <- new_mr_result(
    status = "no_instruments",
    status_reason = "No instruments found"
  )

  expect_message(p <- plot(res), "Cannot plot")
  expect_null(p)
})

test_that("plot.mr_result returns NULL for empty results", {
  res <- new_mr_result(
    results = data.frame(
      exposure = character(),
      outcome = character(),
      method = character(),
      nsnp = integer(),
      b = numeric(),
      se = numeric(),
      pval = numeric()
    ),
    status = "success"
  )

  expect_message(p <- plot(res), "Cannot plot")
  expect_null(p)
})

# --- coloc_result plots ------------------------------------------------------

test_that("plot.coloc_result returns NULL for non-success status", {
  res <- new_coloc_result(
    status = "no_snps_in_region",
    status_reason = "No SNPs in the region"
  )

  expect_message(p <- plot(res), "Cannot plot")
  expect_null(p)
})

test_that("plot.coloc_result pp_bar returns ggplot with ABF results", {
  abf_summary <- c(
    nsnps = 50,
    PP.H0.abf = 0.01,
    PP.H1.abf = 0.05,
    PP.H2.abf = 0.02,
    PP.H3.abf = 0.12,
    PP.H4.abf = 0.80
  )

  res <- new_coloc_result(
    coloc_abf = list(summary = abf_summary),
    n_snps = 50L,
    status = "success"
  )

  p <- plot(res, type = "pp_bar")
  expect_s3_class(p, "ggplot")
})

test_that("plot.coloc_result pp_bar returns NULL without ABF", {
  res <- new_coloc_result(
    coloc_abf = NULL,
    n_snps = 50L,
    status = "success"
  )

  expect_message(p <- plot(res, type = "pp_bar"), "No ABF results")
  expect_null(p)
})

test_that("plot.coloc_result regional returns ggplot with harmonised data", {
  dat <- data.frame(
    pos.exposure = 1:10 * 1000,
    pval.exposure = runif(10, 1e-8, 0.1),
    pos.outcome = 1:10 * 1000,
    pval.outcome = runif(10, 1e-8, 0.1)
  )

  res <- new_coloc_result(
    harmonised_data = dat,
    n_snps = 10L,
    status = "success"
  )

  p <- plot(res, type = "regional")
  expect_s3_class(p, "ggplot")
})

test_that("plot.coloc_result regional returns NULL without data", {
  res <- new_coloc_result(
    harmonised_data = data.frame(),
    n_snps = 0L,
    status = "success"
  )

  expect_message(p <- plot(res, type = "regional"), "No harmonised data")
  expect_null(p)
})

# --- plot_coloc_regional: LD colouring ----------------------------------------

# Five SNPs with a stored r2 unless `r2 = NULL`. `index` is the stored index SNP.
regional_result <- function(r2 = c(1, 0.9, 0.5, 0.1, NA), index = "rs1") {
  dat <- data.frame(
    SNP = paste0("rs", 1:5),
    chr.exposure = "20",
    pos.exposure = 1:5 * 1000,
    beta.exposure = c(0.5, 0.4, 0.3, 0.2, 0.1),
    se.exposure = 0.05,
    pval.exposure = c(1e-20, 1e-10, 1e-5, 0.01, 0.1),
    pos.outcome = 1:5 * 1000,
    pval.outcome = c(1e-3, 1e-2, 0.1, 0.2, 0.3)
  )
  dat$ld_r2_index <- r2
  new_coloc_result(
    harmonised_data = dat,
    n_snps = 5L,
    status = "success",
    ld_index_snp = index
  )
}

# Records the index SNP compute_ld_to_index() is called with
mock_ld_to_index <- function(env = parent.frame()) {
  seen <- new.env()
  testthat::local_mocked_bindings(
    compute_ld_to_index = function(snps, index_snp, bfile, ...) {
      seen$index <- index_snp
      data.frame(
        SNP = rev(snps),
        r2 = rev(seq(1, 0, length.out = length(snps)))
      )
    },
    .env = env
  )
  seen
}

test_that("bin_ld_r2 uses the LocusZoom bins, with a level for missing r2", {
  bins <- bin_ld_r2(c(0, 0.2, 0.21, 0.4, 0.6, 0.8, 1, 1 + 1e-8, NA))

  expect_equal(
    as.character(bins),
    c(
      "0-0.2",
      "0-0.2",
      "0.2-0.4",
      "0.2-0.4",
      "0.4-0.6",
      "0.6-0.8",
      "0.8-1",
      "0.8-1",
      "No LD data"
    )
  )
  expect_equal(levels(bins), names(ld_r2_colours()))
})

test_that("regional plot colours by stored r2 and marks the index SNP", {
  p <- plot(regional_result(), type = "regional")

  expect_s3_class(p, "ggplot")
  points <- ggplot2::layer_data(p, 1)
  index <- ggplot2::layer_data(p, 2)

  # Both panels, without the index SNP in the coloured layer
  expect_equal(nrow(points), 8)
  expect_equal(nrow(index), 2)
  expect_equal(unique(index$x), 1000)
  expect_equal(unique(index$shape), 23)
  expect_equal(
    sort(unique(points$colour)),
    sort(unname(ld_r2_colours()[c("No LD data", "0-0.2", "0.4-0.6", "0.8-1")]))
  )
  # Same r2 at the same position in each panel
  by_pos <- split(points$colour, points$x)
  expect_true(all(lengths(lapply(by_pos, unique)) == 1))
  # Strong LD is drawn last
  expect_equal(
    points$colour[1:2],
    rep(unname(ld_r2_colours()["No LD data"]), 2)
  )
  expect_equal(p$labels$subtitle, "Purple diamond: index SNP rs1")
})

test_that("regional plot keeps every r2 bin in the legend", {
  p <- plot(
    regional_result(r2 = c(1, 0.9, 0.9, 0.9, 0.9)),
    type = "regional"
  )
  colour_scale <- ggplot2::ggplot_build(p)$plot$scales$get_scales("colour")

  expect_equal(colour_scale$get_limits(), names(ld_r2_colours()))
})

test_that("regional plot takes labs(), as rps_amd's figures do", {
  p <- plot(regional_result(), type = "regional") +
    ggplot2::labs(title = "CD40") +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "CD40")
})

test_that("regional plot does not need locuszoomr without ens_db", {
  local_mocked_bindings(
    check_installed = function(...) stop("check_installed called"),
    .package = "rlang"
  )

  expect_s3_class(plot_coloc_regional(regional_result()), "ggplot")
})

test_that("regional plot is uncoloured with no stored r2 and no bfile", {
  res <- regional_result(r2 = NULL, index = NULL)

  expect_message(
    p <- plot(res, type = "regional"),
    "No stored r2.*bfile"
  )

  expect_s3_class(p, "ggplot")
  expect_length(p$layers, 1)
  expect_null(p$labels$subtitle)
})

test_that("regional plot accepts objects saved before ld_index_snp existed", {
  res <- regional_result(r2 = NULL, index = NULL)
  res$ld_index_snp <- NULL
  expect_false("ld_index_snp" %in% names(res))

  expect_no_error(suppressMessages(p <- plot(res, type = "regional")))
  expect_s3_class(p, "ggplot")
})

test_that("regional plot uses stored r2 without PLINK", {
  seen <- mock_ld_to_index()

  plot(regional_result(), type = "regional")
  plot(regional_result(), type = "regional", bfile = "unused")
  plot(regional_result(), type = "regional", index_snp = "rs1")

  expect_null(seen$index)
})

test_that("regional plot recomputes r2 for another index_snp", {
  seen <- mock_ld_to_index()

  p <- plot(
    regional_result(),
    type = "regional",
    index_snp = "rs2",
    bfile = "panel"
  )

  expect_equal(seen$index, "rs2")
  expect_equal(p$labels$subtitle, "Purple diamond: index SNP rs2")
  expect_equal(unique(ggplot2::layer_data(p, 2)$x), 2000)
})

test_that("regional plot computes r2 at plot time when none is stored", {
  seen <- mock_ld_to_index()
  res <- regional_result(r2 = NULL, index = NULL)

  # To the exposure lead by default, or to index_snp
  p <- plot(res, type = "regional", bfile = "panel")
  expect_equal(seen$index, "rs1")
  expect_length(p$layers, 2)

  p <- plot(res, type = "regional", bfile = "panel", index_snp = "rs3")
  expect_equal(seen$index, "rs3")
})

test_that("regional plot marks a recorded index SNP with no r2", {
  seen <- mock_ld_to_index()
  res <- regional_result(r2 = NULL, index = "rs3")

  expect_message(
    p <- plot(res, type = "regional"),
    "No stored r2.*bfile"
  )

  expect_length(p$layers, 2)
  expect_equal(unique(ggplot2::layer_data(p, 2)$x), 3000)
  expect_equal(nrow(ggplot2::layer_data(p, 1)), 8)
  expect_equal(p$labels$subtitle, "Purple diamond: index SNP rs3")
  expect_null(seen$index)
})

test_that("regional plot prefers index_snp, then the recorded SNP, then the lead", {
  seen <- mock_ld_to_index()
  res <- regional_result(r2 = NULL, index = "rs3")

  plot(res, type = "regional", bfile = "panel")
  expect_equal(seen$index, "rs3")

  plot(res, type = "regional", bfile = "panel", index_snp = "rs2")
  expect_equal(seen$index, "rs2")

  plot(regional_result(r2 = NULL, index = NULL), type = "regional", bfile = "x")
  expect_equal(seen$index, "rs1")
})

test_that("regional plot's index_snp overrides a stored one for that plot only", {
  seen <- mock_ld_to_index()
  res <- regional_result(index = "rs2")

  p <- plot(res, type = "regional", index_snp = "rs4", bfile = "panel")

  expect_equal(seen$index, "rs4")
  expect_equal(unique(ggplot2::layer_data(p, 2)$x), 4000)
  expect_equal(res$ld_index_snp, "rs2")
})

test_that("regional plot matches computed r2 to SNPs by name", {
  mock_ld_to_index()
  res <- regional_result(r2 = NULL, index = NULL)

  # The mock returns rs5 first, with r2 rising from 0 at rs5 to 1 at rs1
  ld <- regional_ld(res, bfile = "panel")

  expect_equal(ld$r2, c(1, 0.75, 0.5, 0.25, 0))
})

test_that("regional plot requires bfile for an index_snp with no stored r2", {
  seen <- mock_ld_to_index()

  expect_error(
    plot(regional_result(), type = "regional", index_snp = "rs2"),
    "bfile.*is required"
  )
  expect_error(
    plot(
      regional_result(r2 = NULL, index = NULL),
      type = "regional",
      index_snp = "rs2"
    ),
    "bfile.*is required"
  )
  expect_null(seen$index)
})

test_that("regional plot rejects an index_snp not in the result", {
  seen <- mock_ld_to_index()

  expect_error(
    plot(
      regional_result(),
      type = "regional",
      index_snp = "rsMissing",
      bfile = "panel"
    ),
    "rsMissing.*not among the SNPs"
  )
  expect_null(seen$index)
})

test_that("regional plot computes real r2 and keeps SNPs absent from the panel", {
  skip_if_not_installed("TwoSampleMR")
  bfile <- sub(
    "\\.bed$",
    "",
    system.file("extdata", "ld_ref.bed", package = "mrpipeline")
  )
  skip_if_not(
    file.exists(paste0(bfile, ".bed")),
    "LD reference panel not available"
  )
  bim <- utils::read.table(
    paste0(bfile, ".bim"),
    header = FALSE,
    stringsAsFactors = FALSE
  )
  res <- regional_result(r2 = NULL, index = NULL)
  res$harmonised_data$SNP <- c(head(bim$V2, 4), "rs999999901")

  p <- plot(res, type = "regional", bfile = bfile)

  grey <- unname(ld_r2_colours()["No LD data"])
  points <- ggplot2::layer_data(p, 1)
  expect_equal(sum(points$colour == grey), 2)
  expect_equal(unique(points$x[points$colour == grey]), 5000)
})

test_that("coloured and uncoloured regional plots combine with cowplot", {
  skip_if_not_installed("cowplot")

  panels <- list(
    plot(regional_result(), type = "regional") +
      ggplot2::labs(title = "A"),
    suppressMessages(
      plot(regional_result(r2 = NULL, index = NULL), type = "regional")
    ) +
      ggplot2::labs(title = "B")
  )

  combined <- cowplot::plot_grid(
    plotlist = panels,
    ncol = 2,
    align = "hv",
    axis = "tblr"
  )

  expect_s3_class(combined, "ggplot")
  expect_no_error(ggplot2::ggplotGrob(combined))
})

# --- plot_coloc_regional: gene track ------------------------------------------

test_that("regional plot with a gene track keeps the same points, in Mb", {
  skip_if_not_installed("locuszoomr")
  skip_if_not_installed("ensembldb")
  skip_if_not_installed("cowplot")
  plain <- plot(regional_result(), type = "regional")
  local_mocked_bindings(
    add_gene_track = function(p, ...) p
  )

  p <- plot(regional_result(), type = "regional", ens_db = "unused")

  expect_equal(
    ggplot2::layer_data(p, 1)$x,
    ggplot2::layer_data(plain, 1)$x / 1e6
  )
  expect_equal(ggplot2::layer_data(p, 1)$y, ggplot2::layer_data(plain, 1)$y)
  expect_equal(
    ggplot2::layer_data(p, 1)$colour,
    ggplot2::layer_data(plain, 1)$colour
  )
})

test_that("regional plot's axis matches the gene track's own scale", {
  skip_if_not_installed("locuszoomr")
  skip_if_not_installed("ensembldb")
  skip_if_not_installed("EnsDb.Hsapiens.v75")

  loc <- suppressMessages(locuszoomr::locus(
    seqname = "20",
    xrange = c(44746911, 44758502),
    ens_db = EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75
  ))
  expect_true("CD40" %in% loc$TX$gene_name)

  p <- ggplot2::ggplot(data.frame(x = loc$xrange / 1e6, y = 0:1)) +
    ggplot2::aes(.data$x, .data$y) +
    ggplot2::geom_point() +
    gene_track_x_scale(loc$xrange)
  axis_range <- ggplot2::ggplot_build(p)$layout$panel_params[[1]]$x.range

  # Fails if locuszoomr changes its padding, which would misplace the genes
  expect_equal(
    axis_range,
    locuszoomr::genetracks_grob(loc)$childrenvp$xscale
  )
})

test_that("regional plot with a gene track renders and takes a title", {
  skip_if_not_installed("locuszoomr")
  skip_if_not_installed("ensembldb")
  skip_if_not_installed("cowplot")
  skip_if_not_installed("EnsDb.Hsapiens.v75")
  res <- regional_result()
  res$harmonised_data$pos.exposure <- 44746911 + 0:4 * 2000
  res$harmonised_data$pos.outcome <- res$harmonised_data$pos.exposure

  p <- plot(
    res,
    type = "regional",
    ens_db = EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75
  ) +
    ggplot2::labs(title = "CD40")

  expect_s3_class(p, "ggplot")
  # On a bare plot_grid() labs() draws nothing, without any error
  grob <- ggplot2::ggplotGrob(p)
  title <- grob$grobs[[which(grob$layout$name == "title")]]
  expect_false(inherits(title, "zeroGrob"))

  out_file <- tempfile(fileext = ".png")
  ggplot2::ggsave(out_file, p, width = 8, height = 7, dpi = 72)
  expect_gt(file.size(out_file), 1000)
})

test_that("plot.coloc_result regional ens_db requires the gene track packages", {
  skip_if(
    requireNamespace("locuszoomr", quietly = TRUE) &&
      requireNamespace("ensembldb", quietly = TRUE) &&
      requireNamespace("cowplot", quietly = TRUE),
    "gene track packages are installed -- guard cannot be exercised"
  )

  expect_error(
    plot(regional_result(), type = "regional", ens_db = "EnsDb.Hsapiens.v75"),
    "locuszoomr|ensembldb|cowplot"
  )
})

test_that("plot.coloc_result locuszoom requires locuszoomr and ensembldb", {
  skip_if(
    requireNamespace("locuszoomr", quietly = TRUE) &&
      requireNamespace("ensembldb", quietly = TRUE),
    "locuszoomr and ensembldb are installed -- guard cannot be exercised"
  )

  dat <- data.frame(
    SNP = c("rs1", "rs2", "rs3"),
    chr.exposure = 1,
    pos.exposure = c(1000, 2000, 3000),
    pval.exposure = c(0.01, 1e-6, 0.2),
    chr.outcome = 1,
    pos.outcome = c(1000, 2000, 3000),
    pval.outcome = c(0.05, 1e-4, 0.3)
  )

  res <- new_coloc_result(
    harmonised_data = dat,
    n_snps = 3L,
    status = "success"
  )

  expect_error(
    plot(
      res,
      type = "locuszoom",
      ens_db = "EnsDb.Hsapiens.v75",
      bfile = "unused"
    ),
    "locuszoomr.*ensembldb|ensembldb.*locuszoomr"
  )
})

test_that("plot.coloc_result locuszoom renders with locuszoomr/ensembldb installed", {
  skip_if_not_installed("locuszoomr")
  skip_if_not_installed("ensembldb")
  skip_if_not_installed("EnsDb.Hsapiens.v75")
  # locuszoomr::locus() resolves a character `ens_db` via `get(ens_db)` on
  # the search path, so the annotation package must be attached, not just
  # installed -- matching what a real caller must also do.
  library(EnsDb.Hsapiens.v75)

  bfile <- sub(
    "\\.bed$",
    "",
    system.file("extdata", "ld_ref.bed", package = "mrpipeline")
  )
  skip_if_not(
    file.exists(paste0(bfile, ".bed")),
    "LD reference panel not available"
  )

  bim <- utils::read.table(
    paste0(bfile, ".bim"),
    header = FALSE,
    stringsAsFactors = FALSE
  )
  test_snps <- head(bim$V2, 5)

  dat <- data.frame(
    SNP = test_snps,
    chr.exposure = bim$V1[seq_along(test_snps)],
    pos.exposure = bim$V4[seq_along(test_snps)],
    pval.exposure = c(1e-8, 0.01, 0.02, 0.03, 0.04),
    chr.outcome = bim$V1[seq_along(test_snps)],
    pos.outcome = bim$V4[seq_along(test_snps)],
    pval.outcome = c(0.05, 1e-6, 0.06, 0.07, 0.08)
  )

  res <- new_coloc_result(
    harmonised_data = dat,
    n_snps = length(test_snps),
    status = "success"
  )

  # Like base/grid plotting functions, this draws directly to the current
  # graphics device and returns NULL -- confirm it actually drew something
  # (a non-trivial file), not just that the call completed without error.
  out_file <- tempfile(fileext = ".png")
  grDevices::png(out_file, width = 800, height = 800)
  expect_warning(
    p <- plot(
      res,
      type = "locuszoom",
      ens_db = "EnsDb.Hsapiens.v75",
      bfile = bfile
    ),
    "deprecated.*regional"
  )
  grDevices::dev.off()

  expect_null(p)
  expect_gt(file.size(out_file), 1000)
})

test_that("plot.coloc_result locuszoom without bfile warns and draws uncoloured", {
  skip_if_not_installed("locuszoomr")
  skip_if_not_installed("ensembldb")
  skip_if_not_installed("EnsDb.Hsapiens.v75")
  library(EnsDb.Hsapiens.v75)

  bfile <- sub(
    "\\.bed$",
    "",
    system.file("extdata", "ld_ref.bed", package = "mrpipeline")
  )
  skip_if_not(
    file.exists(paste0(bfile, ".bed")),
    "LD reference panel not available"
  )

  bim <- utils::read.table(
    paste0(bfile, ".bim"),
    header = FALSE,
    stringsAsFactors = FALSE
  )
  test_snps <- head(bim$V2, 5)

  dat <- data.frame(
    SNP = test_snps,
    chr.exposure = bim$V1[seq_along(test_snps)],
    pos.exposure = bim$V4[seq_along(test_snps)],
    pval.exposure = c(1e-8, 0.01, 0.02, 0.03, 0.04),
    chr.outcome = bim$V1[seq_along(test_snps)],
    pos.outcome = bim$V4[seq_along(test_snps)],
    pval.outcome = c(0.05, 1e-6, 0.06, 0.07, 0.08)
  )
  res <- new_coloc_result(
    harmonised_data = dat,
    n_snps = length(test_snps),
    status = "success"
  )

  out_file <- tempfile(fileext = ".png")
  grDevices::png(out_file, width = 800, height = 800)
  expect_warning(
    expect_warning(
      plot(res, type = "locuszoom", ens_db = "EnsDb.Hsapiens.v75"),
      "deprecated.*regional"
    ),
    "no LD colouring.*bfile.*mrpipeline-user-guide"
  )
  grDevices::dev.off()
  expect_gt(file.size(out_file), 1000)

  # A result from run_coloc() with no bfile can hold SNPs absent from a panel
  # used for plotting: they stay on the plot, without an LD value
  dat$SNP[5] <- "rsNotInPanel"
  res <- new_coloc_result(
    harmonised_data = dat,
    n_snps = length(test_snps),
    status = "success"
  )
  out_file <- tempfile(fileext = ".png")
  grDevices::png(out_file, width = 800, height = 800)
  expect_warning(
    p <- plot(
      res,
      type = "locuszoom",
      ens_db = "EnsDb.Hsapiens.v75",
      bfile = bfile
    ),
    "deprecated"
  )
  grDevices::dev.off()
  expect_null(p)
  expect_gt(file.size(out_file), 1000)
})

# --- forest_plot ---------------------------------------------------------

test_that("forest_plot on a single result returns a plot with default methods", {
  res <- new_mr_result(
    results = data.frame(
      exposure = "CD40",
      outcome = "SjD",
      method = c(
        "IVW (fixed effects)",
        "IVW (random effects)",
        "MR Egger",
        "Weighted median"
      ),
      nsnp = 10,
      b = c(0.10, 0.12, 0.15, 0.11),
      se = c(0.03, 0.04, 0.05, 0.04),
      pval = 0.01
    ),
    status = "success"
  )

  p <- forest_plot(res)
  expect_s3_class(p, "ggplot")
})

test_that("forest_plot sections a named list of results", {
  res1 <- new_mr_result(
    results = data.frame(
      exposure = "CD40",
      outcome = "SjD",
      method = "IVW (random effects)",
      nsnp = 10,
      b = 0.1,
      se = 0.03,
      pval = 0.01
    ),
    status = "success"
  )
  res2 <- new_mr_result(
    results = data.frame(
      exposure = "CD40",
      outcome = "Control",
      method = "IVW (random effects)",
      nsnp = 10,
      b = 0.02,
      se = 0.05,
      pval = 0.6
    ),
    status = "success"
  )

  p <- forest_plot(list(Primary = res1, "Negative control" = res2))
  expect_s3_class(p, "ggplot")
})

test_that("forest_plot errors on an unnamed list of length > 1", {
  res <- new_mr_result(
    results = data.frame(
      exposure = "CD40",
      outcome = "SjD",
      method = "IVW (random effects)",
      nsnp = 10,
      b = 0.1,
      se = 0.03,
      pval = 0.01
    ),
    status = "success"
  )

  expect_error(forest_plot(list(res, res)), "must be named")
})

test_that("forest_plot returns NULL when no result has successful rows", {
  res <- new_mr_result(
    status = "no_instruments",
    status_reason = "No instruments found"
  )

  expect_message(p <- forest_plot(res), "Cannot plot")
  expect_null(p)
})

test_that("forest_plot's relabel argument does not break method matching", {
  res <- new_mr_result(
    results = data.frame(
      exposure = "CD40",
      outcome = "SjD",
      method = c("IVW (fixed effects)", "IVW (random effects)"),
      nsnp = 10,
      b = c(0.1, 0.12),
      se = c(0.03, 0.04),
      pval = 0.01
    ),
    status = "success"
  )

  # relabel is applied after matching, so both rows still match.
  p <- forest_plot(
    res,
    methods = c("IVW (fixed effects)", "IVW (random effects)"),
    relabel = c("IVW (random effects)" = "Random")
  )
  expect_s3_class(p, "ggplot")

  # Relabelling disabled (the default) should also still work.
  p_raw <- forest_plot(res, relabel = character(0))
  expect_s3_class(p_raw, "ggplot")
})

# --- outcome_forest_plot ---------------------------------------------------

outcome_forest_plot_fixture <- data.frame(
  exposure = c(
    "RPS",
    "RPS",
    "RPS_functional",
    "RPS_functional",
    "RPS_functional"
  ),
  outcome = c("Melanoma", "Cataract", "Melanoma", "Cataract", "Cataract"),
  method = c(
    "IVW (random effects)",
    "IVW (random effects)",
    "IVW (random effects)",
    "IVW (random effects)",
    "IVW (fixed effects)"
  ),
  or = c(1.4, 0.9, 1.3, 0.95, 0.95),
  or_lci95 = c(1.1, 0.8, 0.9, 0.85, 0.88),
  or_uci95 = c(1.8, 1.0, 1.9, 1.05, 1.02),
  subcategory = c(
    "Positive control",
    "Negative control",
    "Positive control",
    "Negative control",
    "Negative control"
  )
)

test_that("outcome_forest_plot errors when required columns are missing", {
  expect_error(
    outcome_forest_plot(
      dplyr::select(outcome_forest_plot_fixture, -"or"),
      xlab = "OR (95% CI)"
    ),
    "missing column"
  )
})

test_that("outcome_forest_plot with one method includes only matching rows", {
  p <- outcome_forest_plot(
    outcome_forest_plot_fixture,
    xlab = "OR (95% CI)",
    method = "IVW (random effects)"
  )
  expect_s3_class(p, "ggplot")
  # The extra "IVW (fixed effects)" row for Cataract (RPS_functional) is
  # excluded; the other 4 outcome/method combinations remain.
  expect_equal(nrow(p$data), 4)
})

test_that("outcome_forest_plot with both methods keeps rows unique to one model", {
  p <- outcome_forest_plot(
    outcome_forest_plot_fixture,
    xlab = "OR (95% CI)",
    method = c("IVW (random effects)", "IVW (fixed effects)")
  )
  # All 5 rows survive -- Cataract (RPS_functional) now has two.
  expect_equal(nrow(p$data), 5)
  # Relabelled for display.
  expect_setequal(unique(p$data$method), c("Random effects", "Fixed effects"))
})

test_that("outcome_forest_plot maps colour_by/shape_by to the requested columns", {
  dat <- outcome_forest_plot_fixture
  dat$instrument <- ifelse(
    dat$exposure == "RPS",
    "GWS instrument",
    "Functional instrument"
  )

  p <- outcome_forest_plot(
    dat,
    xlab = "OR (95% CI)",
    colour_by = "instrument",
    shape_by = "method"
  )
  expect_identical(rlang::as_label(p$mapping$colour), "instrument")
  expect_identical(rlang::as_label(p$mapping$shape), "method")
})

test_that("outcome_forest_plot's section_order controls subcategory levels", {
  p <- outcome_forest_plot(
    outcome_forest_plot_fixture,
    xlab = "OR (95% CI)",
    section_order = c("Negative control", "Positive control")
  )
  expect_equal(
    levels(p$data$subcategory),
    c("Negative control", "Positive control")
  )
})

test_that("outcome_forest_plot scale='beta' requires b/se, not the OR columns", {
  dat <- outcome_forest_plot_fixture
  dat$b <- log(dat$or)
  dat$se <- (log(dat$or_uci95) - log(dat$or_lci95)) / (2 * 1.96)

  # OR columns absent entirely -- should still work under scale="beta".
  dat_no_or <- dplyr::select(dat, -"or", -"or_lci95", -"or_uci95")
  p <- outcome_forest_plot(dat_no_or, xlab = "Beta (95% CI)", scale = "beta")
  expect_s3_class(p, "ggplot")
  expect_equal(p$data$.x, dat_no_or$b)

  expect_error(
    outcome_forest_plot(
      dplyr::select(dat, -"se"),
      xlab = "Beta (95% CI)",
      scale = "beta"
    ),
    "missing column"
  )
})

test_that("outcome_forest_plot scale='beta' uses a linear x-axis with vline at 0", {
  dat <- outcome_forest_plot_fixture
  dat$b <- log(dat$or)
  dat$se <- (log(dat$or_uci95) - log(dat$or_lci95)) / (2 * 1.96)

  p <- outcome_forest_plot(dat, xlab = "Beta (95% CI)", scale = "beta")
  expect_equal(p$scales$scales[[1]]$trans$name, "identity")
  # ggplot2 4.0 moved geom_vline()'s xintercept from geom_params into a
  # one-row `data` frame on the layer -- not a functional change, just where
  # the value now lives internally.
  expect_equal(p$layers[[1]]$data$xintercept, 0)
})

test_that("outcome_forest_plot's or_col/or_lci_col/or_uci_col let a rescaled column be plotted", {
  dat <- outcome_forest_plot_fixture
  dat$or_scaled <- dat$or * 2
  dat$or_scaled_lci95 <- dat$or_lci95 * 2
  dat$or_scaled_uci95 <- dat$or_uci95 * 2

  p <- outcome_forest_plot(
    dat,
    xlab = "OR (95% CI)",
    or_col = "or_scaled",
    or_lci_col = "or_scaled_lci95",
    or_uci_col = "or_scaled_uci95"
  )
  expect_equal(p$data$.x, dat$or_scaled)
})

test_that("outcome_forest_plot's b_col/se_col let a rescaled column be plotted", {
  dat <- outcome_forest_plot_fixture
  dat$b_scaled <- log(dat$or) * 0.5
  dat$se_scaled <- ((log(dat$or_uci95) - log(dat$or_lci95)) / (2 * 1.96)) * 0.5

  p <- outcome_forest_plot(
    dat,
    xlab = "Beta (95% CI)",
    scale = "beta",
    b_col = "b_scaled",
    se_col = "se_scaled"
  )
  expect_equal(p$data$.x, dat$b_scaled)
})

test_that("outcome_forest_plot's group_by nests an outer nested facet above subcategory", {
  dat <- outcome_forest_plot_fixture
  dat$tier_group <- ifelse(
    dat$subcategory == "Positive control",
    "Primary",
    "Sensitivity Analysis"
  )

  p <- outcome_forest_plot(
    dat,
    xlab = "OR (95% CI)",
    group_by = "tier_group",
    group_order = c("Primary", "Sensitivity Analysis")
  )
  expect_equal(
    levels(p$data$tier_group),
    c("Primary", "Sensitivity Analysis")
  )
  # facet_grid's rows should include both grouping variables.
  facet_vars <- names(p$facet$params$rows)
  expect_setequal(facet_vars, c("tier_group", "subcategory"))
})

test_that("outcome_forest_plot without group_by behaves exactly as before (single facet)", {
  p <- outcome_forest_plot(outcome_forest_plot_fixture, xlab = "OR (95% CI)")
  facet_vars <- names(p$facet$params$rows)
  expect_setequal(facet_vars, "subcategory")
})

# --- table_forest_plot ------------------------------------------------------

table_forest_plot_fixture <- data.frame(
  exposure = "Genetically-proxied NLRP3 inhibition",
  outcome = c("CHD", "Stroke", "Type 2 diabetes"),
  n_snps = c(8, 8, 8),
  estimate = c(0.92, 0.95, 1.03),
  lower = c(0.85, 0.88, 0.94),
  upper = c(0.99, 1.02, 1.13),
  p_value = c(0.02, 0.15, 0.0001)
)

test_that("table_forest_plot errors when required columns are missing", {
  skip_if_not_installed("forestplot")

  expect_error(
    table_forest_plot(
      dplyr::select(table_forest_plot_fixture, -"estimate"),
      null_value = 1,
      xlab = "OR (95% CI)"
    ),
    "missing column"
  )
})

test_that("table_forest_plot returns a forestplot object", {
  skip_if_not_installed("forestplot")

  fp <- table_forest_plot(
    table_forest_plot_fixture,
    null_value = 1,
    xlab = "OR (95% CI)"
  )
  expect_s3_class(fp, "gforge_forestplot")
})

test_that("table_forest_plot does not mutate the caller's data frame", {
  skip_if_not_installed("forestplot")

  dat_before <- table_forest_plot_fixture
  table_forest_plot(table_forest_plot_fixture, null_value = 1, xlab = "x")
  expect_identical(table_forest_plot_fixture, dat_before)
})

test_that("table_forest_plot formats small p-values in scientific notation", {
  skip_if_not_installed("forestplot")

  fp <- table_forest_plot(
    table_forest_plot_fixture,
    null_value = 1,
    xlab = "OR (95% CI)"
  )
  # fp$labels is a forestplot_labeltext list of columns; the last column is
  # "P value", with the header as its first element.
  p_text_col <- unlist(fp$labels[[length(fp$labels)]])[-1]
  expect_true(grepl("e-04", p_text_col[3]))
  expect_false(any(grepl("e", p_text_col[1:2])))
})
