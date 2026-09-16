# Tests for the method registry and mr_methods()

# A run_mr() call that gets past method validation and then fails
# deterministically on the empty exposure -- so an error mentioning "not
# found" proves the method name was accepted, and "Unknown method" (or the
# shortcut error) proves it was not.
run_validation_only <- function(methods) {
  run_mr(
    exposure = data.frame(),
    exposure_id = "e",
    outcome = data.frame(),
    outcome_id = "o",
    methods = methods,
    instruments = "rs1",
    instruments_strict = TRUE,
    verbose = FALSE
  )
}

# --- Registry integrity ------------------------------------------------------

test_that("registry is internally consistent", {
  reg <- mr_method_registry()
  expect_s3_class(reg, "data.frame")
  expect_setequal(
    names(reg),
    c(
      "shortcut",
      "description",
      "label",
      "output",
      "model",
      "ld_correctable",
      "min_instruments",
      "engine",
      "engine_ld"
    )
  )

  shortcuts <- reg$shortcut[!is.na(reg$shortcut)]
  expect_length(shortcuts, 10L)
  expect_equal(anyDuplicated(shortcuts), 0L)
  expect_true(all(nzchar(reg$description)))

  # Every $results-producing row has a label, except the raw passthrough,
  # whose label is whatever TwoSampleMR supplies; diagnostics have none.
  results_rows <- reg$output == "$results"
  passthrough <- reg$engine == "TwoSampleMR::mr"
  expect_false(anyNA(reg$label[results_rows & !passthrough]))
  expect_true(all(is.na(reg$label[!results_rows])))

  # LD support is claimed only by naming the function that provides it,
  # and never claimed without one.
  expect_equal(reg$ld_correctable, !is.na(reg$engine_ld))

  # An effects model is stated exactly where the distinction applies.
  expect_equal(
    !is.na(reg$model),
    reg$shortcut %in% c("ivw_random", "ivw_fixed", "egger")
  )
  expect_equal(reg$model[reg$shortcut %in% "ivw_random"], "random")
  expect_equal(reg$model[reg$shortcut %in% "ivw_fixed"], "fixed")

  expect_true(all(reg$min_instruments >= 1L))
  expect_true(all(
    reg$output %in%
      c("$results", "$steiger", "$pleiotropy", "$heterogeneity", "$loo")
  ))
})

test_that("mr_methods() concise is the requestable subset of full", {
  concise <- mr_methods()
  full <- mr_methods(detail = "full")

  expect_equal(nrow(concise), 10L)
  expect_equal(nrow(full), 12L)
  expect_false(anyNA(concise$shortcut))
  expect_setequal(
    names(concise),
    c(
      "shortcut",
      "description",
      "label",
      "output",
      "model",
      "ld_correctable",
      "min_instruments"
    )
  )
  expect_true(all(names(concise) %in% names(full)))
  expect_equal(
    concise,
    full[!is.na(full$shortcut), names(concise)],
    ignore_attr = TRUE
  )
  expect_error(mr_methods("verbose"), "detail")
})

# --- The registry is what run_mr() accepts -----------------------------------

test_that("every registry shortcut passes run_mr() validation", {
  for (m in mr_methods()$shortcut) {
    expect_error(run_validation_only(m), "not found")
  }
  expect_error(run_validation_only("banana"), "Unknown method")
})

test_that("raw names that are a shortcut's engine are refused", {
  for (raw in c(
    "mr_ivw",
    "mr_ivw_fe",
    "mr_egger_regression",
    "mr_weighted_median"
  )) {
    expect_error(run_validation_only(raw), "shortcut")
  }
  # A raw name with no shortcut still passes through
  expect_error(run_validation_only("mr_ivw_mre"), "not found")
})

test_that("deprecated ivw / ivw_fe map to the new shortcuts with a warning", {
  expect_warning(
    try(run_validation_only(c("ivw", "ivw_fe")), silent = TRUE),
    "deprecated"
  )
})
