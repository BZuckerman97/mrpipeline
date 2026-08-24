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
    plot(res, type = "locuszoom", ens_db = "EnsDb.Hsapiens.v75", bfile = "unused"),
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
  p <- plot(res, type = "locuszoom", ens_db = "EnsDb.Hsapiens.v75", bfile = bfile)
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
        "Inverse variance weighted",
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
      method = "Inverse variance weighted",
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
      method = "Inverse variance weighted",
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
      method = "Inverse variance weighted",
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
      method = c("IVW (fixed effects)", "Inverse variance weighted"),
      nsnp = 10,
      b = c(0.1, 0.12),
      se = c(0.03, 0.04),
      pval = 0.01
    ),
    status = "success"
  )

  # Default relabel (IVW -> "IVW (random effects)") still matches both rows.
  p <- forest_plot(
    res,
    methods = c("IVW (fixed effects)", "Inverse variance weighted")
  )
  expect_s3_class(p, "ggplot")

  # Relabelling disabled should also still work (no-op).
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
    "Inverse variance weighted",
    "Inverse variance weighted",
    "Inverse variance weighted",
    "Inverse variance weighted",
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
    method = "Inverse variance weighted"
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
    method = c("Inverse variance weighted", "IVW (fixed effects)")
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
