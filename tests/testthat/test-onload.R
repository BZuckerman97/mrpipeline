# Loading mrpipeline loads ieugwasr without moving the caller's RNG
# (issue #43). Only a fresh R process can show this -- in the test session
# ieugwasr is long since loaded -- and only for an installed package:
# pkgload::load_all() loads every Imports package before .onLoad() runs, so
# under devtools::test() ieugwasr would already have drawn from the RNG.
# R CMD check (and covr) test the installed package, so this runs there.

load_mrpipeline_fresh <- function(seed, via = c("library", "colon")) {
  via <- match.arg(via)
  callr::r(
    function(seed, via) {
      if (!is.null(seed)) {
        set.seed(seed)
      }
      before <- get0(".Random.seed", globalenv())
      # .onLoad() runs when the namespace loads, so `::` alone (no attach)
      # must behave the same as library().
      if (via == "library") {
        suppressPackageStartupMessages(library(mrpipeline))
      } else {
        invisible(mrpipeline::mr_methods)
      }
      list(
        ieugwasr_loaded = "ieugwasr" %in% loadedNamespaces(),
        seed_unchanged = identical(before, get0(".Random.seed", globalenv()))
      )
    },
    args = list(seed = seed, via = via)
  )
}

test_that("loading mrpipeline, by library() or ::, leaves the RNG state alone", {
  skip_if_not_installed("callr")
  skip_if(
    file.exists(testthat::test_path("..", "..", "DESCRIPTION")),
    "needs the installed package; runs under R CMD check"
  )

  seeded <- load_mrpipeline_fresh(seed = 1)
  expect_true(seeded$ieugwasr_loaded)
  expect_true(seeded$seed_unchanged)

  # With no seed set yet, loading must not create one either.
  unseeded <- load_mrpipeline_fresh(seed = NULL)
  expect_true(unseeded$ieugwasr_loaded)
  expect_true(unseeded$seed_unchanged)

  # Calling through mrpipeline:: without library() loads the namespace too.
  colon <- load_mrpipeline_fresh(seed = 1, via = "colon")
  expect_true(colon$ieugwasr_loaded)
  expect_true(colon$seed_unchanged)
})
