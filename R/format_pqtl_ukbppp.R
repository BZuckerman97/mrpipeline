#' Formats pQTL data from UKB-PPP for MR analysis.
#'
#' This function processes pQTL summary statistics from UKB-PPP,
#' matches them with rsID information, handles non-Mendelian
#' chromosomes, standardizes column names, and formats the data
#' for use with the TwoSampleMR package.
#'
#' @param olink_linker_file Data frame containing linker information
#'   for Olink proteins. Must contain columns: `Code`, `Docname`, `chr`,
#'   `UKBPPP_ProteinID`, `Panel`, `Assay`.
#' @param exposure Character, the UKBPPP_ProteinID of the exposure protein.
#' @param outcome_file_path Character, file path to the outcome data.
#' @param outcome_id Character, ID of the outcome trait.
#' @param ref_rsid Data frame of reference rsIDs for non-Mendelian
#'   chromosomes. Must contain columns: `rsid`, `chr`, `pos`, `ref`, `alt`.
#'
#' @return A list with two elements:
#'   - `exposure`: Formatted exposure data frame (output of TwoSampleMR::format_data).
#'   - `outcome`: Formatted outcome data frame (output of TwoSampleMR::format_data).
#'
#' @export
#'
#' @examples
#' # See the test script for example usage.
#'
format_pqtl_ukbppp <- function(olink_linker_file,
                               exposure,
                               outcome_file_path,
                               outcome_id,
                               ref_rsid) {

    ukbppp_pqtl_file_name <- function(olink_linker_file, synapse_id, olink_dir = "/tmp/olink") {
      file.path(olink_dir,
                paste0(gsub(".tar", "", olink_linker_file[olink_linker_file$Code == synapse_id, ]$Docname[1]),
                       "/",
                       "discovery_chr",
                       olink_linker_file[olink_linker_file$Code == synapse_id, ]$chr[1],
                       "_",
                       olink_linker_file[olink_linker_file$Code == synapse_id, ]$UKBPPP_ProteinID[1],
                       ":",
                       olink_linker_file[olink_linker_file$Code == synapse_id, ]$Panel[1],
                       ".gz"
                )
      )
    }

  ukbppp_rsid_file_name <- function(olink_linker_file, metadata_dir = "/tmp/metadata") {
    file.path(metadata_dir,
              paste0(olink_linker_file$chr,
                     "_patched_v2.tsv.gz"))
  }

  # 1. & 2. Create file paths using the ukbppp functions
  exposure_synapse_id <- olink_linker_file[olink_linker_file$UKBPPP_ProteinID == exposure, ]$Code[1]
  if (length(exposure_synapse_id) == 0) {
    stop("Exposure ID not found in linker file. Ensure the UKBPPP_ProteinID is present in the linker file")
  }

  exposure_file_path <- ukbppp_pqtl_file_name(olink_linker_file, exposure_synapse_id)
  exposure_chr <- olink_linker_file[olink_linker_file$UKBPPP_ProteinID == exposure, ]$chr[1]
  rsid_file_path <- ukbppp_rsid_file_name(olink_linker_file[olink_linker_file$UKBPPP_ProteinID == exposure, ])

  # 3. Read in files using data.table::fread()
  exposure_data <- data.table::fread(exposure_file_path, nThread = parallel::detectCores())
  rsid_data <- data.table::fread(rsid_file_path, nThread = parallel::detectCores())

  if (!"Assay" %in% colnames(olink_linker_file)) {
    stop("Assay variable not present in the linker file")
  }

  #create the exposure_id variable
  exposure_id <- olink_linker_file[olink_linker_file$UKBPPP_ProteinID == exposure,]$Assay[1]

  # 4. Standardize column names
  exposure_data <- exposure_data %>%
    dplyr::rename(
      phenotype = dplyr::all_of("Assay"), # using all_of to prevent errors if column is missing
      rsid = dplyr::all_of("rsid"),
      beta = dplyr::all_of("BETA"),
      sebeta = dplyr::coalesce(dplyr::all_of(c("SE","sebeta","se"))), #coalesce is used to select the first non-NA value from the list of columns
      af_alt = dplyr::all_of("A1FREQ"),
      effect_allele = dplyr::coalesce(dplyr::all_of(c("ALLELE1","alt"))),
      other_allele = dplyr::coalesce(dplyr::all_of(c("ALLELE0","ref"))),
      pval = dplyr::all_of("P"),
      chr = dplyr::all_of("CHROM"),
      pos = dplyr::coalesce(dplyr::all_of(c("POS19","GENPOS","pos")))
    ) %>%
    dplyr::mutate(phenotype = exposure_id) %>%
    dplyr::select(phenotype, rsid, beta, sebeta, af_alt, effect_allele, other_allele, pval, chr, pos) %>%
    dplyr::mutate(chr = dplyr::if_else(chr == "23", "X", chr)) #change 23 to X if needed

   # 5. Handle non-Mendelian chromosomes
  if (exposure_chr %in% c("X", "Y")) {
    exposure_data <- dplyr::inner_join(exposure_data, ref_rsid, by = c("chr", "pos"))
  }
    # Match by ID or create it
  if ("ID" %in% colnames(exposure_data) & "ID" %in% colnames(rsid_data)) {
    exposure_data <- dplyr::inner_join(exposure_data, rsid_data, by = "ID")
  } else {
    exposure_data <- exposure_data %>%
      dplyr::mutate(ID = paste(chr, pos, effect_allele, other_allele, sep = ":"))
    rsid_data <- rsid_data %>%
      dplyr::mutate(ID = paste(chr, pos, alt, ref, sep = ":"))

    exposure_data <- dplyr::inner_join(exposure_data, rsid_data, by = "ID")
  }

  # Read in outcome data using data.table::fread()
  outcome_data <- data.table::fread(outcome_file_path, nThread = parallel::detectCores())
  outcome_data$outcome <- outcome_id

  # 4. Standardize column names
  outcome_data <- outcome_data %>%
    dplyr::rename(
      rsid = dplyr::all_of("rsid"),
      beta = dplyr::all_of("beta"),
      sebeta = dplyr::coalesce(dplyr::all_of(c("sebeta", "SE", "se"))),
      af_alt = dplyr::all_of("A1FREQ"),
      effect_allele = dplyr::coalesce(dplyr::all_of(c("ALLELE1","alt"))),
      other_allele = dplyr::coalesce(dplyr::all_of(c("ALLELE0","ref"))),
      pval = dplyr::all_of("P"),
      chr = dplyr::all_of("CHROM"),
      pos = dplyr::coalesce(dplyr::all_of(c("GENPOS", "POS19","pos")))
    ) %>%
    dplyr::select(outcome, rsid, beta, sebeta, af_alt, effect_allele, other_allele, pval, chr, pos) %>%
    dplyr::rename(phenotype = outcome) %>%
    dplyr::mutate(chr = dplyr::if_else(chr == "23", "X", chr)) #change 23 to X if needed

  # 5. Format data using TwoSampleMR::format_data()
  exposure_formatted <- TwoSampleMR::format_data(
    exposure_data,
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

  outcome_formatted <- TwoSampleMR::format_data(
    outcome_data,
    type = "outcome",
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

  return(list(exposure = exposure_formatted, outcome = outcome_formatted))
}
