# Unit tests for get_gene_coords()
# All tests require biomaRt namespace (even mocked tests use local_mocked_bindings)

skip_if_not_installed("biomaRt")

# Helper: fake getBM output
fake_getBM <- function(attributes, filters, values, mart) {
  genes <- data.frame(
    hgnc_symbol = c("CD40", "CD40", "APOE", "FAKEGENE"),
    chromosome_name = c("20", "20", "19", "CHR_HSCHR1_1_CTG3"),
    start_position = c(44746911L, 44746500L, 44905791L, 1000L),
    end_position = c(44758502L, 44760000L, 44909393L, 2000L),
    stringsAsFactors = FALSE
  )
  genes[genes$hgnc_symbol %in% values, ]
}

test_that("get_gene_coords returns correct structure", {
  local_mocked_bindings(
    useEnsembl = function(...) "mock_mart",
    getBM = fake_getBM,
    .package = "biomaRt"
  )

  result <- get_gene_coords(c("CD40", "APOE"))

  expect_s3_class(result, "tbl_df")
  expect_named(result, c("hgnc_symbol", "chromosome", "start", "end"))
  expect_equal(nrow(result), 2)

  cd40 <- result[result$hgnc_symbol == "CD40", ]
  expect_equal(cd40$start, 44746500L)
  expect_equal(cd40$end, 44760000L)
  expect_equal(cd40$chromosome, "20")
})

test_that("non-standard chromosomes are filtered out", {
  local_mocked_bindings(
    useEnsembl = function(...) "mock_mart",
    getBM = fake_getBM,
    .package = "biomaRt"
  )

  result <- get_gene_coords("FAKEGENE")

  # FAKEGENE is on CHR_HSCHR1_1_CTG3, should be filtered out
  expect_equal(nrow(result), 0)
  expect_warning(
    get_gene_coords("FAKEGENE"),
    "not found"
  )
})

test_that("unknown gene produces warning", {
  local_mocked_bindings(
    useEnsembl = function(...) "mock_mart",
    getBM = function(...) {
      data.frame(
        hgnc_symbol = character(),
        chromosome_name = character(),
        start_position = integer(),
        end_position = integer()
      )
    },
    .package = "biomaRt"
  )

  expect_warning(
    result <- get_gene_coords("NOTAREALGENE"),
    "not found"
  )
  expect_equal(nrow(result), 0)
})

test_that("build = 'grch37' passes correct host", {
  captured_host <- NULL
  local_mocked_bindings(
    useEnsembl = function(...) {
      args <- list(...)
      captured_host <<- args$host
      "mock_mart"
    },
    getBM = fake_getBM,
    .package = "biomaRt"
  )

  get_gene_coords("CD40", build = "grch37")
  expect_equal(captured_host, "grch37.ensembl.org")
})

test_that("build = 'grch38' does not pass grch37 host", {
  captured_host <- "sentinel"
  local_mocked_bindings(
    useEnsembl = function(...) {
      args <- list(...)
      captured_host <<- args$host
      "mock_mart"
    },
    getBM = fake_getBM,
    .package = "biomaRt"
  )

  get_gene_coords("CD40", build = "grch38")
  expect_null(captured_host)
})

test_that("autosomes preferred over sex chromosomes", {
  local_mocked_bindings(
    useEnsembl = function(...) "mock_mart",
    getBM = function(attributes, filters, values, mart) {
      data.frame(
        hgnc_symbol = c("GENEX", "GENEX"),
        chromosome_name = c("X", "1"),
        start_position = c(1000L, 2000L),
        end_position = c(3000L, 4000L),
        stringsAsFactors = FALSE
      )
    },
    .package = "biomaRt"
  )

  result <- get_gene_coords("GENEX")
  expect_equal(result$chromosome, "1")
})

# Integration test — also needs network
test_that("get_gene_coords works with real biomaRt (integration)", {
  skip_on_cran()
  skip_if_offline()

  result <- tryCatch(
    get_gene_coords("CD40", build = "grch38"),
    error = function(e) {
      skip(paste("Ensembl unavailable:", conditionMessage(e)))
    }
  )

  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 1)
  expect_equal(result$hgnc_symbol, "CD40")
  expect_equal(result$chromosome, "20")
  expect_true(result$start > 0)
  expect_true(result$end > result$start)
})
