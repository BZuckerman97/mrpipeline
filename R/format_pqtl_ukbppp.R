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
#' @param x_y_chr_file String, file path to the file containing rsids for X and Y chromosomes
#'
#' @return A list with two elements:
#'   - `exposure`: Formatted exposure data frame (output of TwoSampleMR::format_data).
#'
#' @export
#'
#' @examples
#' # See the test script for example usage.
format_pqtl_ukbppp <- function(ukbppp,
                               ukbppp_rsid,
                               pqtl_assay,
                               x_y_chr_file = NULL) {

  # read from filepath
  if (is.character(ukbppp)) {
    # Read in files using data.table::fread()
    stopifnot(file.exists(ukbppp))

    ukbppp <- ukbppp |>
      purrr::map(\(x) data.table::fread(x, nThread = parallel::detectCores())) |>
      dplyr::bind_rows()
  } else {
    stopifnot(is.data.frame(ukbppp))
  }

  if (is.character(ukbppp_rsid)) {
    # Read in files using data.table::fread()
    ukbppp_rsid <- ukbppp_rsid |>
      purrr::map(\(x) data.table::fread(x, nThread = parallel::detectCores())) |>
      dplyr::bind_rows()
  } else {
    stopifnot(is.data.frame(ukbppp_rsid))
  }

  # Standardise column names
  # UKB-PPP
  ukbppp <- ukbppp |>
    dplyr::mutate(phenotype = pqtl_assay) |>
    dplyr::rename(
      beta = dplyr::all_of("BETA"), # using all_of to prevent errors if column is missing
      sebeta = dplyr::all_of("SE"),
      af_alt = dplyr::all_of("A1FREQ"),
      effect_allele = dplyr::all_of("ALLELE1"),
      other_allele = dplyr::all_of("ALLELE0"),
      pval = dplyr::all_of("LOG10P"),
      chr = dplyr::all_of("CHROM"),
      n = dplyr::all_of("N"),
      pos = dplyr::all_of("GENPOS") #' Is there a way of altering this dependent on whether you use a build37 or build38 data
    ) |>
    dplyr::mutate(chr = dplyr::if_else(chr == "23", "X", as.character(chr))) |>  #change 23 to X if needed
    dplyr::mutate(pval = 10 ^ -pval) # Convert LOG10P into P

  # UKB-PPP RSID
  ukbppp_rsid <- ukbppp_rsid |>
    dplyr::rename(
      effect_allele_rsid_file = dplyr::all_of("ALT"),
      other_allele_rsid_file = dplyr::all_of("REF"),
      pos_rsid_file = dplyr::all_of("POS38")) |>
    dplyr::mutate(chr_rsid_file = sub(":.*", "", ID))

  # Match by ID or create it
  if ("ID" %in% colnames(ukbppp) & "ID" %in% colnames(ukbppp_rsid)) {
    ukbppp <- dplyr::inner_join(ukbppp, ukbppp_rsid, by = "ID")
  } else {
    ukbppp <- ukbppp |>
      dplyr::mutate(ID = paste(chr, pos, effect_allele, other_allele, sep = ":"))
    ukbppp_rsid <- ukbppp_rsid |>
      dplyr::mutate(ID = paste(chr_rsid_file, pos_rsid_file, effect_allele_rsid_file, other_allele_rsid_file, sep = ":"))

    ukbppp <- dplyr::inner_join(ukbppp, ukbppp_rsid, by = "ID")
  }

  # Handle non-Mendelian chromosomes rsIDs
  # Check if the chromosome is X
  if (!is.null(x_y_chr_file)) {
    stopifnot(file.exists(x_y_chr_file))
    if ("X" %in% unique(ukbppp$chr)) {
      # Load x_y_rsid
      x_y_rsid <- data.table::fread(x_y_chr_file)
      # Rename columns to match
      x_y_rsid <- x_y_rsid |>
        dplyr::rename(
          pos = dplyr::all_of("V2"),
          rsids = dplyr::all_of("V3"),
          chr = dplyr::all_of("V1")
        )
      # Merge with ukbppp by position
      ukbppp <- dplyr::left_join(ukbppp, x_y_rsid[, c("pos", "rsids")], by = "pos")

      # Update rsid column with rsids from x_y_rsid
      ukbppp <- ukbppp |>
        dplyr::mutate(rsid = dplyr::coalesce(rsids, rsid)) |>
        dplyr::select(-rsids)
    }
  }
  # Convert data table to data frame for format_data()
  ukbppp <- as.data.frame(ukbppp)

  # Format data using TwoSampleMR::format_data()
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
    samplesize_col = "n",
    pos_col = "pos",
    log_pval = FALSE
  )

  return(result)
}

#' UKBPPP_PQTL_FILE_NAME
#'
#' @param synapse_id String, synapse id to access from olink linker file
#' @param olink_linker_file String or Dataframe, file containing the olink linker file, or the dataframe of the linker file
#' @param olink_dir String, directory of the olink files
#' @param olink_rsid_dir String, directory of the olink rsid files
#'
#' @returns a list of the 2 filepaths, one for the ukbppp_pqtl data and the other is
#' the corresponding rsID metadata file, as well as the name of the assay
#' @export
#'
#' @examples
#' \dontrun{
#' Assuming you have a linker file and the deCODE data directory
#' synapse_id
#' olink_linker_file <- "path/to/your/olink_linker_file.csv"
#' olink_dir <- "path/to/olink_dir"
#' olink_rsid_dir <- "path/to/olink_rsid_dir"
#' file_paths <- ukbppp_pqtl_file_name(synapse_id = synapse_id,
#'                                   olink_linker_file = olink_linker_file,
#'                                   olink_dir = olink_dir,
#'                                   olink_rsid_dir = olink_rsid_dir)
#' print(file_paths)
#' }
ukbppp_pqtl_file_name <- function(synapse_id,
                                  olink_linker_file,
                                  olink_dir,
                                  olink_rsid_dir) {

  if (rlang::is_string(olink_linker_file)) {
    stopifnot(file.exists(olink_linker_file))
    olink_linker_file <- data.table::fread(olink_linker_file)
  } else {
    stopifnot(is.data.frame(olink_linker_file))
  }

  # get relevant metadata
  metadata <- olink_linker_file |>
    dplyr::filter(Code == synapse_id)

  if (nrow(metadata) == 0) {
    stop(paste0("No entry found for synapse_id '", synapse_id, "' in the linker file."))
  }

  # Handle cases where one synapse ID maps to multiple proteins (e.g. AMY1A/B/C).
  # As long as the file paths and assay name are the same, we can proceed.
  if (nrow(metadata) > 1) {
    # Check for consistency in columns required for file paths and phenotype ID.
    cols_to_check <- c("Docname", "chr", "UKBPPP_ProteinID", "Panel", "Assay")

    # Use vapply for a safer and more explicit check for consistency.
    is_consistent <- vapply(metadata[, cols_to_check], function(x) data.table::uniqueN(x) == 1, logical(1))

    if (!all(is_consistent)) {
      inconsistent_cols <- names(is_consistent)[!is_consistent]
      stop(paste0("Multiple entries for synapse_id '", synapse_id,
                  "' have conflicting metadata. Inconsistent columns: ",
                  paste(inconsistent_cols, collapse = ", "),
                  ". Cannot determine a unique file path or assay name."))
    }
    # All good, just take the first row.
    metadata <- metadata[1, ]
  }
  metadata <- as.list(metadata)

  # result
  list(
    ukbppp = file.path(olink_dir,
                       paste0(gsub(".tar", "", metadata$Docname),
                              "/",
                              "discovery_chr",
                              metadata$chr,
                              "_",
                              metadata$UKBPPP_ProteinID,
                              ":",
                              metadata$Panel,
                              ".gz"
                       )
    ),
    ukbppp_rsid = file.path(olink_rsid_dir,
                            paste0("olink_rsid_map_mac5_info03_b0_7_chr", metadata$chr, "_patched_v2.tsv.gz")
    ),
    id = metadata$Assay
  )
}
