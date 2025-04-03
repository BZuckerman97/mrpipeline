#' Formats pQTL data from UKB-PPP for MR analysis.
#'
#' This function processes pQTL summary statistics from UKB-PPP,
#' matches them with rsID information, handles non-Mendelian
#' chromosomes, standardizes column names, and formats the data
#' for use with the TwoSampleMR package.
#'
#' @param ukbppp Dataframe, the ukbppp data
#' @param ukbppp_rsid Dataframe, of ukbppp rsids
#' @param pqtl_assay String, of the ukbppp protein assayed
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
                               pqtl_assay) {

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
    dplyr::rename(
      phenotype = dplyr::all_of(pqtl_assay), # using all_of to prevent errors if column is missing
      beta = dplyr::all_of("BETA"),
      sebeta = dplyr::all_of("SE"),
      af_alt = dplyr::all_of("A1FREQ"),
      effect_allele = dplyr::all_of("ALLELE1"),
      other_allele = dplyr::all_of("ALLELE0"),
      pval = dplyr::all_of("P"),
      chr = dplyr::all_of("CHROM"),
      pos = dplyr::all_of("GENPOS") #' Is there a way of altering this dependent on whether you use a build37 or build38 data
    ) %>%
    dplyr::select(phenotype, rsid, beta, sebeta, af_alt, effect_allele, other_allele, pval, chr, pos) %>%
    dplyr::mutate(chr = dplyr::if_else(chr == "23", "X", chr)) |>  #change 23 to X if needed
    dplyr::mutate(pval = 10 ^ -pval) # Convert LOG10P into P

  # UKB-PPP RSID
  ukbppp_rsid <- ukbppp_rsid |>
    dplyr::rename(
      effect_allele = dplyr::all_of("ALT"),
      other_allele = dplyr::all_of("REF")
      chr = dplyr::mutate(chr = as.numeric(gsub("chr", "", chr))),
      pos = dplyr::all_of("GENPOS")
    )

  # Handle non-Mendelian chromosomes
  # Handle X chromosome rsIDs
  # Check if the chromosome is X
  if (!is.null(ref_rsid_file)) {
    stopifnot(file.exists(ref_rsid_file))
    if ("X" %in% unique(ukbppp$chr)) {
      # Load ref_rsid
      ref_rsid <- data.table::fread(ref_rsid_file)
      # Rename columns to match
      ref_rsid <- ref_rsid |>
        dplyr::rename(
          pos = dplyr::all_of("V2"),
          rsids = dplyr::all_of("V3"),
          chr = dplyr::all_of("V1")
        )
      # Merge with ukbppp by position
      ukbppp <- dplyr::left_join(ukbppp, ref_rsid[, c("pos", "rsids")], by = "pos")

      # Update rsid column with rsids from ref_rsid
      ukbppp <- ukbppp |>
        dplyr::mutate(rsid = dplyr::coalesce(rsids, rsid)) |>
        dplyr::select(-rsids)
    }
  }

  # Match by ID or create it
  if ("ID" %in% colnames(ukbppp) & "ID" %in% colnames(ukbppp_rsid)) {
    ukbppp <- dplyr::inner_join(ukbppp, ukbppp_rsid, by = "ID")
  } else {
    ukbppp <- ukbppp %>%
      dplyr::mutate(ID = paste(chr, pos, effect_allele, other_allele, sep = ":"))
    ukbppp_rsid <- ukbppp_rsid %>%
      dplyr::mutate(ID = paste(chr, pos, alt, ref, sep = ":"))

    ukbppp <- dplyr::inner_join(ukbppp, ukbppp_rsid, by = "ID")
  }

  # Format data using TwoSampleMR::format_data()
  result <- TwoSampleMR::format_data(
    ukbppp,
    type = "exposure",
    phenotype_col = "phenotype",
    snp_col = "rsid",
    beta_col = "beta",
    se_col = "sebeta",
    eaf_col = "af_alt",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "pval",
    chr_col = "chr",
    pos_col = "pos"
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
ukbppp_pqtl_file_name <- function(synapse_id, olink_linker_file, olink_dir, olink_rsid_dir) {

  if (rlang::is_string(olink_linker_file)) {
    stopifnot(file.exists(olink_linker_file))
    olink_linker_file <- data.table::fread(olink_linker_file)
  } else {
    stopifnot(is.data.frame(olink_linker_file))
  }

  # get relevant metadata
  metadata <- olink_linker_file[olink_linker_file$Code == synapse_id, ]
  stopifnot(identical(nrow(metadata), 1L))
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
                            paste0("olink_rsid_map_mac4_info03_b0_7_chr", metadata$chr, "_patched_v2.tsv.gz")
    ),
    id = metadata$Assay
  )
}

# A---------
# This code was not part of any function and has now been removed:
# # TEMP workflow
#  synapse_id <- c("A", "B")

#  synapse_id |>
#    purrr::set_names() |>
#    purrr::map(\(x) ukbppp_pqtl_file_name(synapse_id = x, olink_linker_file = olink_linker_file, olink_dir = "TODO", olink_rsid_dir = "TODO") |>
#    dplyr::bind_rows() |>
#    purrr::pmap(format_pqtl_ukbppp) # this bit to be amended
# A---------
