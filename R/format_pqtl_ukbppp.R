#' Formats pQTL data from UKB-PPP for MR analysis.
#'
#' This function processes pQTL summary statistics from UKB-PPP,
#' matches them with rsID information, handles non-Mendelian
#' chromosomes, standardizes column names, and formats the data
#' for use with the TwoSampleMR package.
#'
#' @param ukbppp Dataframe, file path to the file containing the ukbppp GWAS data
#' @param ukbppp_rsid Dataframe, file path to the file containing the ukbppp GWAS data rsids
#' @param pqtl_assay String, of the ukbppp protein assay
#' @param x_y_chr_file Data frame or string file path containing rsIDs for X
#'   and Y chromosomes, or NULL. When a data frame is supplied, it is used
#'   directly (skipping `fread()`). When a string path is supplied, the file is
#'   read via `data.table::fread()`.
#' @param pos_build Character, genome build to use for position matching against
#'   the rsid file. Either `"b38"` (GRCh38, uses `POS38` column; default) or
#'   `"b37"` (GRCh37, uses `POS19` column). Must match the build used in the
#'   UKB-PPP GWAS `GENPOS` column.
#' @param type Character, either `"exposure"` (default) or `"outcome"`.
#'   `"exposure"` runs the data through [TwoSampleMR::format_data()].
#'   `"outcome"` skips `format_data()` and returns a plain data frame
#'   normalised to the same schema as `format_gwas(type = "outcome")`, so
#'   [run_mr()] can pre-filter by rsids before calling `format_data()`
#'   internally.
#'
#' @return A plain data frame (not a list) in both cases: for
#'   `type = "exposure"`, the output of TwoSampleMR::format_data(); for
#'   `type = "outcome"`, data normalised to `format_gwas(type = "outcome")`'s
#'   schema.
#'
#' @export
#'
#' @examples
#' # See the test script for example usage.
format_pqtl_ukbppp <- function(
  ukbppp,
  ukbppp_rsid,
  pqtl_assay,
  x_y_chr_file = NULL,
  pos_build = c("b38", "b37"),
  type = c("exposure", "outcome")
) {
  type <- match.arg(type)
  pos_build <- match.arg(pos_build)
  # read from filepath
  if (is.character(ukbppp)) {
    # Read in files using data.table::fread()
    stopifnot(file.exists(ukbppp))

    ukbppp <- ukbppp |>
      purrr::map(\(x) {
        data.table::fread(x, nThread = parallel::detectCores())
      }) |>
      dplyr::bind_rows()
  } else {
    stopifnot(is.data.frame(ukbppp))
  }

  if (is.character(ukbppp_rsid)) {
    # Read in files using data.table::fread()
    ukbppp_rsid <- ukbppp_rsid |>
      purrr::map(\(x) {
        data.table::fread(x, nThread = parallel::detectCores())
      }) |>
      dplyr::bind_rows()
  } else {
    stopifnot(is.data.frame(ukbppp_rsid))
  }

  # Standardise column names
  # UKB-PPP
  ukbppp <- ukbppp |>
    dplyr::mutate(phenotype = pqtl_assay) |>
    dplyr::rename(
      # all_of() (rather than a bare name) avoids an error if a column is
      # missing.
      beta = dplyr::all_of("BETA"),
      sebeta = dplyr::all_of("SE"),
      af_alt = dplyr::all_of("A1FREQ"),
      effect_allele = dplyr::all_of("ALLELE1"),
      other_allele = dplyr::all_of("ALLELE0"),
      pval = dplyr::all_of("LOG10P"),
      chr = dplyr::all_of("CHROM"),
      # GENPOS is UKB-PPP's raw column name regardless of genome build --
      # pos_build only affects which rsid-file column (POS38/POS19) it's
      # matched against below, not this rename.
      pos = dplyr::all_of("GENPOS")
    ) |>
    # N deliberately left un-renamed here (unlike the other columns above) --
    # it must stay optional for type = "outcome" (see below), so renaming it
    # unconditionally via dplyr::all_of() would hard-error for outcome files
    # that lack a sample-size column, exactly the case the outcome branch's
    # NA_integer_ fallback exists to handle.
    dplyr::mutate(
      chr = dplyr::if_else(.data$chr == "23", "X", as.character(.data$chr))
    ) |> #change 23 to X if needed
    dplyr::mutate(pval = 10^-.data$pval) # Convert LOG10P into P

  # UKB-PPP RSID
  pos_col_rsid <- if (pos_build == "b38") "POS38" else "POS19"
  ukbppp_rsid <- ukbppp_rsid |>
    dplyr::rename(
      effect_allele_rsid_file = dplyr::all_of("ALT"),
      other_allele_rsid_file = dplyr::all_of("REF"),
      pos_rsid_file = dplyr::all_of(pos_col_rsid)
    ) |>
    dplyr::mutate(chr_rsid_file = sub(":.*", "", .data$ID))

  # Match by ID or create it
  if ("ID" %in% colnames(ukbppp) && "ID" %in% colnames(ukbppp_rsid)) {
    ukbppp <- dplyr::inner_join(ukbppp, ukbppp_rsid, by = "ID")
  } else {
    ukbppp <- ukbppp |>
      dplyr::mutate(
        ID = paste(
          .data$chr,
          .data$pos,
          .data$effect_allele,
          .data$other_allele,
          sep = ":"
        )
      )
    ukbppp_rsid <- ukbppp_rsid |>
      dplyr::mutate(
        ID = paste(
          .data$chr_rsid_file,
          .data$pos_rsid_file,
          .data$effect_allele_rsid_file,
          .data$other_allele_rsid_file,
          sep = ":"
        )
      )

    ukbppp <- dplyr::inner_join(
      ukbppp,
      ukbppp_rsid[, c("ID", "rsid")],
      by = "ID"
    )
  }

  # Handle non-Mendelian chromosomes rsIDs
  # Check if the chromosome is X
  if (!is.null(x_y_chr_file) && "X" %in% unique(ukbppp$chr)) {
    if (is.character(x_y_chr_file)) {
      stopifnot(file.exists(x_y_chr_file))
      x_y_rsid <- data.table::fread(
        x_y_chr_file,
        nThread = parallel::detectCores()
      )
    } else {
      stopifnot(is.data.frame(x_y_chr_file))
      x_y_rsid <- x_y_chr_file
    }
    # Rename columns to match
    x_y_rsid <- x_y_rsid |>
      dplyr::rename(
        pos = dplyr::all_of("V2"),
        rsids = dplyr::all_of("V3"),
        chr = dplyr::all_of("V1")
      )
    # Merge with ukbppp by position
    ukbppp <- dplyr::left_join(
      ukbppp,
      x_y_rsid[, c("pos", "rsids")],
      by = "pos"
    )

    # Update rsid column with rsids from x_y_rsid
    ukbppp <- ukbppp |>
      dplyr::mutate(rsid = dplyr::coalesce(.data$rsids, .data$rsid)) |>
      dplyr::select(-dplyr::all_of("rsids"))
  }
  # Convert data table to data frame for format_data()
  ukbppp <- as.data.frame(ukbppp)

  # Outcome: return normalised data frame (same schema as format_gwas(type="outcome"))
  # so that run_mr() can pre-filter by rsids and then call format_data() internally.
  # n is optional here -- rename N -> n if present, else default to NA_integer_
  # (mirroring format_pqtl_decode()'s outcome branch).
  if (type == "outcome") {
    ukbppp <- ukbppp |>
      dplyr::rename(rsids = "rsid", se = "sebeta", eaf = "af_alt")
    if ("N" %in% names(ukbppp) && !"n" %in% names(ukbppp)) {
      ukbppp <- dplyr::rename(ukbppp, n = "N")
    }
    if (!"n" %in% names(ukbppp)) {
      ukbppp$n <- NA_integer_
    }
    return(ukbppp)
  }

  # Format data using TwoSampleMR::format_data()
  # samplesize_col references the raw "N" column directly -- unlike the other
  # columns, N is not renamed above (see comment there). format_data() does
  # not hard-require the named samplesize_col to exist (it silently leaves
  # samplesize.exposure unset if absent), matching format_pqtl_decode()'s
  # exposure path, which points samplesize_col at its own un-renamed "N"
  # the same way.
  result <- TwoSampleMR::format_data(
    ukbppp,
    type = "exposure",
    header = TRUE,
    phenotype_col = "phenotype",
    snp_col = "rsid",
    beta_col = "beta",
    se_col = "sebeta",
    eaf_col = "af_alt",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "pval",
    chr_col = "chr",
    samplesize_col = "N",
    pos_col = "pos",
    log_pval = FALSE
  )

  result
}

#' UKBPPP_PQTL_FILE_NAME
#'
#' @param synapse_id String, synapse id to access from olink linker file
#' @param olink_linker_file String or Dataframe, file containing the olink linker file,
#' or the dataframe of the linker file. This function is designed to handle the
#' full linker file and will perform the necessary filtering and validation internally.
#' It can handle cases where a single `synapse_id` maps to multiple rows (e.g., for
#' protein complexes), as long as the underlying file paths are consistent.
#' @param olink_dir String, directory of the olink files
#' @param olink_rsid_dir String, directory of the olink rsid files
#' @param region_chr String or NULL. If supplied, overrides the protein's own
#'   gene chromosome (from the linker file) when building both file paths.
#'   Use this when reading the protein as an *outcome* at a locus other than
#'   its own gene (e.g. reversed-role drug-target MR, where a fixed exposure
#'   region is being looked up across many outcome proteins). Default `NULL`
#'   preserves the original behaviour exactly (uses the protein's own gene
#'   chromosome).
#'
#' @returns a list of the 2 filepaths, one for the ukbppp_pqtl data and the other is
#' the corresponding rsID metadata file, as well as the name of the assay
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming you have a linker file and the data directory
#' synapse_id <- "syn12345678"
#' olink_linker_file <- "path/to/your/olink_linker_file.csv"
#' olink_dir <- "path/to/olink_dir"
#' olink_rsid_dir <- "path/to/olink_rsid_dir"
#' file_paths <- ukbppp_pqtl_file_name(synapse_id = synapse_id,
#'                                   olink_linker_file = olink_linker_file,
#'                                   olink_dir = olink_dir,
#'                                   olink_rsid_dir = olink_rsid_dir)
#' print(file_paths)
#' }
ukbppp_pqtl_file_name <- function(
  synapse_id,
  olink_linker_file,
  olink_dir,
  olink_rsid_dir,
  region_chr = NULL
) {
  if (rlang::is_string(olink_linker_file)) {
    stopifnot(file.exists(olink_linker_file))
    olink_linker_file <- data.table::fread(olink_linker_file)
  } else {
    stopifnot(is.data.frame(olink_linker_file))
  }

  # get relevant metadata
  metadata <- olink_linker_file |>
    dplyr::filter(.data$Code == .env$synapse_id)

  if (nrow(metadata) == 0) {
    stop(paste0(
      "No entry found for synapse_id '",
      synapse_id,
      "' in the linker file."
    ))
  }

  # Handle cases where one synapse ID maps to multiple proteins (e.g. AMY1A/B/C).
  # As long as the file paths and assay name are the same, we can proceed.
  if (nrow(metadata) > 1) {
    # Check for consistency in columns required for file paths and phenotype ID.
    cols_to_check <- c("Docname", "chr", "UKBPPP_ProteinID", "Panel", "Assay")

    # Use vapply for a safer and more explicit check for consistency.
    is_consistent <- vapply(
      metadata[, cols_to_check],
      function(x) data.table::uniqueN(x) == 1,
      logical(1)
    )

    if (!all(is_consistent)) {
      inconsistent_cols <- names(is_consistent)[!is_consistent]
      stop(paste0(
        "Multiple entries for synapse_id '",
        synapse_id,
        "' have conflicting metadata. Inconsistent columns: ",
        paste(inconsistent_cols, collapse = ", "),
        ". Cannot determine a unique file path or assay name."
      ))
    }
    # All good, just take the first row.
    metadata <- metadata[1, ]
  }
  metadata <- as.list(metadata)

  chr_to_use <- if (!is.null(region_chr)) region_chr else metadata$chr

  # result
  list(
    ukbppp = file.path(
      olink_dir,
      paste0(
        gsub(".tar", "", metadata$Docname),
        "/",
        "discovery_chr",
        chr_to_use,
        "_",
        metadata$UKBPPP_ProteinID,
        ":",
        metadata$Panel,
        ".gz"
      )
    ),
    ukbppp_rsid = file.path(
      olink_rsid_dir,
      paste0(
        "olink_rsid_map_mac5_info03_b0_7_chr",
        chr_to_use,
        "_patched_v2.tsv.gz"
      )
    ),
    id = metadata$Assay
  )
}
