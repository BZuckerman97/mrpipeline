#' Look up gene coordinates from Ensembl via biomaRt
#'
#' Queries Ensembl for genomic coordinates of one or more HGNC gene symbols.
#' Results are filtered to standard chromosomes (1--22, X, Y), deduplicated
#' (widest range per gene, preferring autosomes), and returned as a tibble.
#'
#' @param genes Character vector of HGNC gene symbols.
#' @param build Genome build: `"grch38"` (default) or `"grch37"`.
#'
#' @return A tibble with columns `hgnc_symbol`, `chromosome`, `start`, `end`.
#'   Genes not found in Ensembl are dropped with a warning.
#'
#' @examples
#' \dontrun{
#' get_gene_coords("CD40")
#' get_gene_coords(c("CD40", "APOE"), build = "grch37")
#' }
#'
#' @export
get_gene_coords <- function(genes, build = c("grch38", "grch37")) {
  rlang::check_installed("biomaRt", reason = "to look up gene coordinates.")
  build <- match.arg(build)

  host <- if (build == "grch37") "grch37.ensembl.org" else NULL

  mart_args <- list(
    biomart = "genes",
    dataset = "hsapiens_gene_ensembl"
  )
  if (!is.null(host)) {
    mart_args$host <- host
  }
  mart <- do.call(biomaRt::useEnsembl, mart_args)

  raw <- biomaRt::getBM(
    attributes = c(
      "hgnc_symbol",
      "chromosome_name",
      "start_position",
      "end_position"
    ),
    filters = "hgnc_symbol",
    values = unique(genes),
    mart = mart
  )

  # Filter to standard chromosomes
  standard_chr <- c(as.character(1:22), "X", "Y")
  result <- raw |>
    dplyr::filter(.data$chromosome_name %in% standard_chr) |>
    # Collapse duplicates to widest range per gene + chromosome
    dplyr::group_by(.data$hgnc_symbol, .data$chromosome_name) |>
    dplyr::summarise(
      start = min(.data$start_position),
      end = max(.data$end_position),
      .groups = "drop"
    ) |>
    # Prefer autosomes over sex chromosomes, then keep one row per gene
    dplyr::arrange(
      .data$hgnc_symbol,
      !.data$chromosome_name %in% as.character(1:22)
    ) |>
    dplyr::distinct(.data$hgnc_symbol, .keep_all = TRUE) |>
    dplyr::rename(chromosome = "chromosome_name")

  # Warn about genes not found
  missing <- setdiff(unique(genes), result$hgnc_symbol)
  if (length(missing) > 0) {
    cli::cli_warn(
      "{length(missing)} gene{?s} not found in Ensembl ({build}): {.val {missing}}."
    )
  }

  result
}
