# LD-corrected diagnostics: $heterogeneity, $pleiotropy, $loo (issue #31)

# --- loo_correlated() against the naive refit loop --------------------------

# A random dense correlation structure, as an MRInput-shaped frame.
random_correlated_set <- function(n, seed = 11) {
  set.seed(seed)
  a <- matrix(stats::rnorm(n * 3), n)
  ld <- stats::cov2cor(tcrossprod(a) + diag(n) * 3)
  dimnames(ld) <- list(paste0("rs", seq_len(n)), paste0("rs", seq_len(n)))
  harmonised <- data.frame(
    SNP = paste0("rs", seq_len(n)),
    beta.exposure = stats::rnorm(n, 0, 0.3),
    se.exposure = stats::runif(n, 0.02, 0.05),
    beta.outcome = NA_real_,
    se.outcome = stats::runif(n, 0.005, 0.01),
    exposure = "X",
    outcome = "Y",
    id.exposure = "x",
    id.outcome = "y",
    samplesize.outcome = 5000,
    stringsAsFactors = FALSE
  )
  harmonised$beta.outcome <- 0.1 *
    harmonised$beta.exposure +
    stats::rnorm(n, 0, 0.01)
  list(harmonised = harmonised, ld = ld)
}

# The reference implementation: one mr_ivw(correl = TRUE) refit per SNP.
naive_loo <- function(h, ld) {
  n <- nrow(h)
  rows <- t(vapply(
    seq_len(n),
    function(i) {
      s <- seq_len(n)[-i]
      fit <- MendelianRandomization::mr_ivw(
        MendelianRandomization::mr_input(
          bx = h$beta.exposure[s],
          bxse = h$se.exposure[s],
          by = h$beta.outcome[s],
          byse = h$se.outcome[s],
          correlation = ld[s, s]
        ),
        correl = TRUE,
        model = "random"
      )
      c(fit@Estimate, fit@StdError, fit@Pvalue)
    },
    numeric(3)
  ))
  all <- MendelianRandomization::mr_ivw(
    MendelianRandomization::mr_input(
      bx = h$beta.exposure,
      bxse = h$se.exposure,
      by = h$beta.outcome,
      byse = h$se.outcome,
      correlation = ld
    ),
    correl = TRUE,
    model = "random"
  )
  rbind(rows, c(all@Estimate, all@StdError, all@Pvalue))
}

test_that("loo_correlated() reproduces per-SNP mr_ivw(correl = TRUE) refits", {
  d <- random_correlated_set(30)
  got <- loo_correlated(d$harmonised, d$ld)
  ref <- naive_loo(d$harmonised, d$ld)

  expect_equal(nrow(got), 31)
  expect_equal(got$SNP, c(d$harmonised$SNP, "All"))
  expect_equal(got$b, ref[, 1], tolerance = 1e-10)
  expect_equal(got$se, ref[, 2], tolerance = 1e-10)
  expect_equal(got$p, ref[, 3], tolerance = 1e-10)
  expect_true(all(got$ld_corrected))
  expect_equal(
    names(got),
    c(
      "exposure",
      "outcome",
      "id.exposure",
      "id.outcome",
      "samplesize",
      "SNP",
      "b",
      "se",
      "p",
      "ld_corrected"
    )
  )
})

test_that("loo_correlated()'s direct-solve fallback gives the same answer", {
  d <- random_correlated_set(12, seed = 3)
  fast <- loo_correlated(d$harmonised, d$ld)
  # rcond_min = 1 forces the per-SNP solve path on a well-conditioned matrix
  slow <- loo_correlated(d$harmonised, d$ld, rcond_min = 1)
  expect_equal(fast, slow, tolerance = 1e-10)
})

# --- End to end, on a panel with real LD ------------------------------------

# Correlated quantities hand-computed from the result's own instruments and
# LD matrix: the GLS estimate, the generalised Q, and the correlated fits.
correlated_reference <- function(res) {
  h <- res$instruments
  ld <- res$ld_matrix
  omega <- diag(h$se.outcome) %*% ld %*% diag(h$se.outcome)
  w <- solve(omega)
  bx <- h$beta.exposure
  by <- h$beta.outcome
  b <- as.numeric(solve(t(bx) %*% w %*% bx) %*% (t(bx) %*% w %*% by))
  r <- by - b * bx
  input <- MendelianRandomization::mr_input(
    bx = bx,
    bxse = h$se.exposure,
    by = by,
    byse = h$se.outcome,
    correlation = ld
  )
  list(
    n = nrow(h),
    b = b,
    q = as.numeric(t(r) %*% w %*% r),
    egger = MendelianRandomization::mr_egger(input, correl = TRUE),
    ivw_random = MendelianRandomization::mr_ivw(
      input,
      correl = TRUE,
      model = "random"
    )
  )
}

test_that("with instruments in LD, the diagnostics differ from the uncorrected run and match the correlated fits", {
  skip_if_not_installed("TwoSampleMR")
  f <- synthetic_five_snps()
  bfile <- make_ld_panel(f$exposure, ld_pair = c(1, 2), copy_prob = 0.85)
  methods <- c("ivw_fixed", "ivw_random", "egger", "heterogeneity", "loo")

  corr <- run_quiet(
    f$exposure,
    f$outcome,
    instruments = f$exposure$SNP,
    bfile = bfile,
    ld_correct = TRUE,
    methods = methods
  )$result
  unc <- run_quiet(
    f$exposure,
    f$outcome,
    instruments = f$exposure$SNP,
    bfile = bfile,
    methods = methods
  )$result
  expect_equal(corr$status, "success")
  expect_equal(unc$status, "success")

  # The panel really does carry LD
  expect_gt(abs(corr$ld_matrix["rs1", "rs2"]), 0.5)
  ref <- correlated_reference(corr)
  expect_equal(ref$n, 5)

  # The reporter's case: run_mr(ld_correct = TRUE, methods =
  # c("ivw_fixed", "heterogeneity")) must not return the uncorrected Q.
  q_corr <- corr$heterogeneity$Q[
    corr$heterogeneity$method == "Inverse variance weighted"
  ]
  q_unc <- unc$heterogeneity$Q[
    unc$heterogeneity$method == "Inverse variance weighted"
  ]
  expect_gt(abs(q_corr - q_unc) / q_unc, 0.05)
  expect_equal(q_corr, ref$q, tolerance = 1e-8)
  expect_true(all(corr$heterogeneity$ld_corrected))
  expect_false(any(unc$heterogeneity$ld_corrected))

  # ... and ivw_fixed's own SE moves too -- the assertion the bundled panel
  # could never support.
  se_corr <- corr$results$se[corr$results$method == "IVW (fixed effects)"]
  se_unc <- unc$results$se[unc$results$method == "IVW (fixed effects)"]
  expect_gt(abs(se_corr - se_unc) / se_unc, 0.05)

  # Pleiotropy: the correlated intercept, and a different one
  expect_equal(
    corr$pleiotropy$egger_intercept,
    ref$egger@Intercept,
    tolerance = 1e-8
  )
  expect_equal(corr$pleiotropy$se, ref$egger@StdError.Int, tolerance = 1e-8)
  expect_true(corr$pleiotropy$ld_corrected)
  expect_false(unc$pleiotropy$ld_corrected)
  expect_false(isTRUE(all.equal(
    corr$pleiotropy$egger_intercept,
    unc$pleiotropy$egger_intercept
  )))

  # Leave-one-out: correlated refits, differing from the uncorrected ones
  expect_true(all(corr$loo$ld_corrected))
  expect_false(any(unc$loo$ld_corrected))
  expect_equal(nrow(corr$loo), 6)
  expect_false(isTRUE(all.equal(corr$loo$b[1:5], unc$loo$b[1:5])))
  all_row <- corr$loo[corr$loo$SNP == "All", ]
  expect_equal(all_row$b, ref$ivw_random@Estimate, tolerance = 1e-8)
  expect_equal(all_row$se, ref$ivw_random@StdError, tolerance = 1e-8)
})

test_that("a singular LD matrix is reported at the source", {
  skip_if_not_installed("TwoSampleMR")
  f <- synthetic_five_snps()
  # copy_prob = 1: rs1 and rs2 are identical, so R has a zero eigenvalue
  bfile <- make_ld_panel(f$exposure, ld_pair = c(1, 2), copy_prob = 1)

  # Only a method with no GLS fit, so the run survives the singular matrix
  # and the warning about it is what we see.
  out <- run_quiet(
    f$exposure,
    f$outcome,
    instruments = f$exposure$SNP,
    bfile = bfile,
    ld_correct = TRUE,
    methods = "weighted_median"
  )
  expect_true(any(grepl("near-singular", out$warnings)))
  expect_true(any(grepl("identical or near-identical", out$warnings)))
})
