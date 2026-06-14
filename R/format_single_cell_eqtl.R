# -- Internal helper ----------------------------------------------------------

#' Skip VCF meta-information lines and return data as a data.frame
#'
#' Counts lines beginning with `##`, reads the `#CHROM` header, then uses
#' [data.table::fread()] to read the remaining rows with the correct column
#' names.
#'
#' @param file Path to a VCF file (plain or gzip-compressed).
#' @return A data frame of VCF records.
#' @keywords internal
read_vcf_data <- function(file) {
  con <- file(file, "r")
  on.exit(close(con), add = TRUE)

  n_skip <- 0L
  header_line <- NULL

  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (length(line) == 0L) {
      break
    }
    if (stringr::str_starts(line, "##")) {
      n_skip <- n_skip + 1L
    } else {
      header_line <- line
      break
    }
  }

  if (is.null(header_line)) {
    cli::cli_abort(
      "No VCF column header found (expected a line starting with '#CHROM')."
    )
  }

  col_names <- stringr::str_split(
    stringr::str_remove(header_line, "^#"),
    "\t"
  )[[1]]

  data.table::fread(
    file,
    skip = n_skip + 1L,
    header = FALSE,
    col.names = col_names,
    nThread = parallel::detectCores()
  ) |>
    as.data.frame()
}


# -- Exported format functions ------------------------------------------------

#' Format single-cell RNA eQTL data from OneK1K for MR analysis
#'
#' Reads and formats eQTL summary statistics from the OneK1K project for a
#' specific cell type, producing a data frame ready for [run_mr()] or
#' [run_coloc()].
#'
#' OneK1K data is available from <https://onek1k.org/>. The expected input
#' format is the per-cell-type or combined eQTL table with columns:
#' `CHR`, `POS`, `RSID`, `A1`, `A2`, `A2_FREQ_ONEK1K`,
#' `SPEARMANS_RHO`, `P_VALUE`, `GENE`, and optionally `CELL_ID`.
#'
#' `A2` is the effect allele (the allele whose frequency is reported in
#' `A2_FREQ_ONEK1K`) and `A1` is the other allele. Spearman's rho
#' (`SPEARMANS_RHO`) is used as the effect size (beta). Standard errors are
#' derived from the p-value:
#' `se = |rho| / qnorm(P_VALUE / 2, lower.tail = FALSE)`.
#' Rows with non-finite or zero SE are dropped. The exposure phenotype label
#' is `GENE___<cell_type>`.
#'
#' If the file contains a `CELL_ID` column (combined multi-cell-type file),
#' rows are filtered to those matching `onek1k_cell_type`.
#'
#' @param onek1k_mapping A data frame with columns `cell_type` and
#'   `path_to_eqtl_file`, mapping each cell type to its eQTL summary
#'   statistics file.
#' @param onek1k_cell_type Character. The cell type to format (e.g.,
#'   `"cd4nc"`). Must match an entry in the `cell_type` column of
#'   `onek1k_mapping`. When the file contains a `CELL_ID` column, rows are
#'   filtered to this value.
#' @param sample_size Integer or `NULL`. Sample size for this cell type,
#'   stored as `samplesize.exposure`. Required by [run_coloc()]; if `NULL`
#'   (default), no sample size column is added.
#'
#' @return A data frame formatted as a TwoSampleMR exposure dataset, including
#'   `chr.exposure`, `pos.exposure`, and `eaf.exposure` columns required by
#'   [run_coloc()].
#'
#' @seealso [run_mr()], [run_coloc()]
#' @family sceqtl-format
#' @export
#'
#' @examples
#' \dontrun{
#' mapping_df <- data.frame(
#'   cell_type         = "cd4nc",
#'   path_to_eqtl_file = "esnp_table.tsv.gz"
#' )
#' exposure <- format_single_cell_onek1k(
#'   onek1k_mapping   = mapping_df,
#'   onek1k_cell_type = "cd4nc",
#'   sample_size      = 463528L
#' )
#' }
format_single_cell_onek1k <- function(
  onek1k_mapping,
  onek1k_cell_type,
  sample_size = NULL
) {
  if (!is.data.frame(onek1k_mapping)) {
    cli::cli_abort("{.arg onek1k_mapping} must be a data frame.")
  }
  if (!rlang::is_string(onek1k_cell_type)) {
    cli::cli_abort("{.arg onek1k_cell_type} must be a single character string.")
  }

  missing_cols <- setdiff(
    c("cell_type", "path_to_eqtl_file"),
    names(onek1k_mapping)
  )
  if (length(missing_cols) > 0L) {
    cli::cli_abort(
      "{.arg onek1k_mapping} is missing column{?s}: {.val {missing_cols}}."
    )
  }

  cell_type_info <- onek1k_mapping |>
    dplyr::filter(.data$cell_type == .env$onek1k_cell_type)

  if (nrow(cell_type_info) == 0L) {
    cli::cli_abort(
      "Cell type {.val {onek1k_cell_type}} not found in {.arg onek1k_mapping}."
    )
  }
  if (nrow(cell_type_info) > 1L) {
    cli::cli_warn(
      "Multiple entries found for {.val {onek1k_cell_type}}; using the first."
    )
    cell_type_info <- cell_type_info[1L, ]
  }

  eqtl_file_path <- cell_type_info$path_to_eqtl_file
  cli::cli_inform("Reading OneK1K eQTL data for {.val {onek1k_cell_type}}...")

  eqtl_data <- data.table::fread(
    eqtl_file_path,
    nThread = parallel::detectCores()
  ) |>
    as.data.frame()

  # Filter by CELL_ID when the file spans multiple cell types
  if ("CELL_ID" %in% names(eqtl_data)) {
    eqtl_data <- eqtl_data |>
      dplyr::filter(.data$CELL_ID == .env$onek1k_cell_type)
    if (nrow(eqtl_data) == 0L) {
      cli::cli_abort(
        "No rows for {.val {onek1k_cell_type}} after filtering {.field CELL_ID}."
      )
    }
  }

  # Derive SE; keep only rows with finite, positive SE
  eqtl_data <- eqtl_data |>
    dplyr::mutate(
      se_derived = abs(.data$SPEARMANS_RHO) /
        stats::qnorm(.data$P_VALUE / 2, lower.tail = FALSE),
      phenotype_label = paste(.data$GENE, .env$onek1k_cell_type, sep = "___")
    ) |>
    dplyr::filter(is.finite(.data$se_derived) & .data$se_derived > 0)

  if (nrow(eqtl_data) == 0L) {
    cli::cli_abort("No variants with valid SE remain after filtering.")
  }

  if (!is.null(sample_size)) {
    eqtl_data <- eqtl_data |>
      dplyr::mutate(samplesize_col = as.integer(.env$sample_size))
    ss_col <- "samplesize_col"
  } else {
    ss_col <- NULL
  }

  TwoSampleMR::format_data(
    dat = eqtl_data,
    type = "exposure",
    phenotype_col = "phenotype_label",
    snp_col = "RSID",
    beta_col = "SPEARMANS_RHO",
    se_col = "se_derived",
    eaf_col = "A2_FREQ_ONEK1K",
    effect_allele_col = "A2",
    other_allele_col = "A1",
    pval_col = "P_VALUE",
    chr_col = "CHR",
    pos_col = "POS",
    samplesize_col = ss_col
  )
}


#' Format single-cell RNA eQTL data from 1M-scBloodNL for MR analysis
#'
#' Reads and formats eQTL summary statistics from the 1 Million Immune Cells
#' single-cell blood cohort (1M-scBloodNL; Oelen et al. 2022) for a single
#' cell-type file.
#'
#' Beta and SE are derived from the Z-score and total sample size:
#' `se = 1 / sqrt(N_total)`, `beta = OverallZScore * se`, where `N_total` is
#' the sum of per-dataset sample counts from the `DatasetsNrSamples` column
#' (semicolon-separated integers). The other allele is derived by splitting
#' `SNPType` on `/` and taking the allele not equal to `AlleleAssessed`. Rows
#' where the other allele cannot be unambiguously determined are dropped. The
#' exposure phenotype label is `HGNCName___<cell_type>`.
#'
#' @param file Character. Path to a 1M-scBloodNL per-cell-type eQTL file
#'   (e.g., `"CD4T_expression_eQTLsFDR-ProbeLevel.txt.gz"`).
#' @param cell_type Character or `NULL`. Cell type label used in the phenotype
#'   string (`GENE___cell_type`). If `NULL`, derived from the filename by
#'   stripping `_expression_eQTLsFDR-ProbeLevel.txt.gz`.
#' @param cis_only Logical. If `TRUE` (default), retain only rows where
#'   `CisTrans == "cis"`.
#'
#' @return A data frame formatted as a TwoSampleMR exposure dataset, including
#'   `chr.exposure` and `pos.exposure` columns required by [run_coloc()].
#'   Note: EAF is not available in 1M-scBloodNL files; `eaf.exposure` will be
#'   absent and `exposure_n` must be supplied explicitly to [run_coloc()].
#'
#' @seealso [run_mr()], [run_coloc()]
#' @family sceqtl-format
#' @export
#'
#' @examples
#' \dontrun{
#' exposure <- format_sceqtl_1m_scbloodnl(
#'   file = "CD4T_expression_eQTLsFDR-ProbeLevel.txt.gz"
#' )
#' }
format_sceqtl_1m_scbloodnl <- function(
  file,
  cell_type = NULL,
  cis_only = TRUE
) {
  if (!file.exists(file)) {
    cli::cli_abort("File not found: {.path {file}}")
  }

  if (is.null(cell_type)) {
    cell_type <- stringr::str_remove(
      basename(file),
      "_expression_eQTLsFDR-ProbeLevel\\.txt\\.gz$"
    )
    if (identical(cell_type, basename(file))) {
      cli::cli_warn(
        "Could not infer cell type from filename; set {.arg cell_type} explicitly."
      )
    }
    cli::cli_inform("Inferred cell type: {.val {cell_type}}")
  }

  cli::cli_inform("Reading 1M-scBloodNL data from {.path {basename(file)}}...")

  dt <- data.table::fread(file, nThread = parallel::detectCores()) |>
    as.data.frame()

  if (nrow(dt) == 0L) {
    cli::cli_abort("File is empty: {.path {file}}")
  }

  if (cis_only && "CisTrans" %in% names(dt)) {
    dt <- dt[dt$CisTrans == "cis", , drop = FALSE]
    if (nrow(dt) == 0L) {
      cli::cli_abort(
        "No cis eQTLs remain after filtering on {.field CisTrans}."
      )
    }
  }

  # Derive total N from semicolon-separated per-dataset sample counts
  sample_counts <- stringr::str_split(dt$DatasetsNrSamples, ";")
  dt$n_total <- vapply(
    sample_counts,
    function(x) sum(suppressWarnings(as.numeric(x)), na.rm = TRUE),
    numeric(1)
  )

  # Derive beta and se from Z-score and N
  dt$se <- 1 / sqrt(dt$n_total)
  dt$beta_derived <- as.numeric(dt$OverallZScore) * dt$se

  # Derive alleles from SNPType ("A/G" -> allele1 = A, allele2 = G)
  allele_split <- stringr::str_split_fixed(dt$SNPType, "/", n = 2L)
  dt$allele1 <- allele_split[, 1]
  dt$allele2 <- allele_split[, 2]

  dt <- dt |>
    dplyr::mutate(
      effect_allele = .data$AlleleAssessed,
      other_allele = ifelse(
        .data$AlleleAssessed == .data$allele1,
        .data$allele2,
        ifelse(
          .data$AlleleAssessed == .data$allele2,
          .data$allele1,
          NA_character_
        )
      )
    )

  n_dropped <- sum(is.na(dt$other_allele))
  if (n_dropped > 0L) {
    cli::cli_inform(
      "Dropping {n_dropped} row{?s} where other allele could not be determined."
    )
    dt <- dt[!is.na(dt$other_allele), , drop = FALSE]
  }

  if (nrow(dt) == 0L) {
    cli::cli_abort("No variants remain after allele derivation.")
  }

  dt$phenotype_label <- paste(dt$HGNCName, cell_type, sep = "___")

  TwoSampleMR::format_data(
    dat = dt,
    type = "exposure",
    phenotype_col = "phenotype_label",
    snp_col = "SNPName",
    beta_col = "beta_derived",
    se_col = "se",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "PValue",
    chr_col = "SNPChr",
    pos_col = "SNPChrPos"
  )
}


#' Format single-cell RNA eQTL data from DICE for MR analysis
#'
#' Reads and formats eQTL summary statistics from the Database of Immune Cell
#' Expression, eQTLs and Epigenomics (DICE; Schmiedel et al. 2018) for a
#' single cell-type VCF file.
#'
#' The VCF INFO field is parsed with regex to extract `Gene`, `GeneSymbol`,
#' `Pvalue`, and `Beta`. SE is derived as
#' `se = |Beta| / sqrt(qchisq(Pvalue, df = 1, lower.tail = FALSE))`. Only
#' biallelic SNPs (single-nucleotide, A/C/G/T) are retained. Palindromic SNPs
#' (A/T and C/G pairs) are dropped because DICE provides no allele frequency
#' data to resolve strand ambiguity. `ALT` is treated as the effect allele and
#' `REF` as the other allele. The exposure phenotype label is
#' `GeneSymbol___<cell_type>`, falling back to `Gene` when `GeneSymbol` is
#' absent.
#'
#' @param file Character. Path to a DICE eQTL VCF file
#'   (e.g., `"t_cell_cd4_naive.vcf"`). Plain and gzip-compressed files are
#'   both supported.
#' @param cell_type Character or `NULL`. Cell type label used in the phenotype
#'   string (`GENE___cell_type`). If `NULL`, derived from the filename by
#'   stripping the `.vcf` (or `.vcf.gz`) extension.
#'
#' @return A data frame formatted as a TwoSampleMR exposure dataset, including
#'   `chr.exposure` and `pos.exposure` columns required by [run_coloc()].
#'   Note: EAF is not available in DICE VCF files; `eaf.exposure` will be
#'   absent and `exposure_n` must be supplied explicitly to [run_coloc()].
#'
#' @seealso [run_mr()], [run_coloc()]
#' @family sceqtl-format
#' @export
#'
#' @examples
#' \dontrun{
#' exposure <- format_sceqtl_dice(file = "t_cell_cd4_naive.vcf")
#' }
format_sceqtl_dice <- function(file, cell_type = NULL) {
  if (!file.exists(file)) {
    cli::cli_abort("File not found: {.path {file}}")
  }

  if (is.null(cell_type)) {
    cell_type <- stringr::str_remove(basename(file), "\\.vcf(\\.gz)?$")
    cli::cli_inform("Inferred cell type: {.val {cell_type}}")
  }

  cli::cli_inform("Reading DICE VCF from {.path {basename(file)}}...")

  dt <- read_vcf_data(file)

  if (nrow(dt) == 0L) {
    cli::cli_abort("No data rows found in {.path {file}}.")
  }

  # Parse INFO field key=value pairs via regex
  # Use anchored pattern for Gene= to avoid matching GeneSymbol=
  dt$Gene <- stringr::str_match(dt$INFO, "(?:^|;)Gene=([^;]+)")[, 2]
  dt$GeneSymbol <- stringr::str_match(dt$INFO, "GeneSymbol=([^;]+)")[, 2]
  dt$Pvalue <- as.numeric(
    stringr::str_match(dt$INFO, "Pvalue=([^;]+)")[, 2]
  )
  dt$Beta <- as.numeric(
    stringr::str_match(dt$INFO, "Beta=([^;]+)")[, 2]
  )

  dt$POS <- as.integer(dt$POS)
  dt$CHROM <- stringr::str_remove(dt$CHROM, "^[Cc][Hh][Rr]")

  # Retain biallelic single-nucleotide variants only
  valid_bases <- c("A", "C", "G", "T")
  dt <- dt |>
    dplyr::filter(
      nchar(.data$REF) == 1L &
        nchar(.data$ALT) == 1L &
        .data$REF %in% valid_bases &
        .data$ALT %in% valid_bases
    )

  if (nrow(dt) == 0L) {
    cli::cli_abort("No biallelic SNP-only variants remain.")
  }

  # Drop palindromic SNPs (no EAF available to resolve strand ambiguity)
  dt <- dt |>
    dplyr::filter(
      !((.data$REF == "A" & .data$ALT == "T") |
        (.data$REF == "T" & .data$ALT == "A") |
        (.data$REF == "C" & .data$ALT == "G") |
        (.data$REF == "G" & .data$ALT == "C"))
    )

  if (nrow(dt) == 0L) {
    cli::cli_abort("No non-palindromic variants remain.")
  }

  # Derive SE from beta and chi-squared p-value; drop rows with invalid SE
  dt <- dt |>
    dplyr::mutate(
      se_derived = abs(.data$Beta) /
        sqrt(stats::qchisq(.data$Pvalue, df = 1L, lower.tail = FALSE))
    ) |>
    dplyr::filter(is.finite(.data$se_derived) & .data$se_derived > 0)

  if (nrow(dt) == 0L) {
    cli::cli_abort("No variants with valid SE remain.")
  }

  # Phenotype: prefer GeneSymbol, fall back to Gene
  dt <- dt |>
    dplyr::mutate(
      gene_label = dplyr::if_else(
        !is.na(.data$GeneSymbol) & nchar(.data$GeneSymbol) > 0L,
        .data$GeneSymbol,
        .data$Gene
      ),
      phenotype_label = paste(.data$gene_label, cell_type, sep = "___")
    )

  TwoSampleMR::format_data(
    dat = dt,
    type = "exposure",
    phenotype_col = "phenotype_label",
    snp_col = "ID",
    beta_col = "Beta",
    se_col = "se_derived",
    effect_allele_col = "ALT",
    other_allele_col = "REF",
    pval_col = "Pvalue",
    chr_col = "CHROM",
    pos_col = "POS"
  )
}


#' Format single-cell RNA eQTL data from dynamic_cseqtl for MR analysis
#'
#' Reads and formats eQTL summary statistics from the dynamic cis single-cell
#' eQTL (dynamic_cseqtl) dataset for a single cell-type file.
#'
#' Input files are expected to be pre-formatted tab-separated files with
#' columns: `gene`, `SNP`, `beta`, `se`, `effect_allele`, `other_allele`,
#' `pval`, `chr`, `pos`. The naming convention for files is
#' `<cell_type>_500kb_combined.MR.tsv.gz`. The exposure phenotype label is
#' `gene___<cell_type>`.
#'
#' @param file Character. Path to a dynamic_cseqtl MR file
#'   (e.g., `"CD4T_500kb_combined.MR.tsv.gz"`).
#' @param cell_type Character or `NULL`. Cell type label used in the phenotype
#'   string (`gene___cell_type`). If `NULL`, derived from the filename by
#'   stripping `_500kb_combined.MR.tsv.gz`.
#'
#' @return A data frame formatted as a TwoSampleMR exposure dataset, including
#'   `chr.exposure` and `pos.exposure` columns required by [run_coloc()].
#'   Note: EAF is not available in dynamic_cseqtl files; `eaf.exposure` will
#'   be absent and `exposure_n` must be supplied explicitly to [run_coloc()].
#'
#' @seealso [run_mr()], [run_coloc()]
#' @family sceqtl-format
#' @export
#'
#' @examples
#' \dontrun{
#' exposure <- format_sceqtl_dynamic_cseqtl(
#'   file = "CD4T_500kb_combined.MR.tsv.gz"
#' )
#' }
format_sceqtl_dynamic_cseqtl <- function(file, cell_type = NULL) {
  if (!file.exists(file)) {
    cli::cli_abort("File not found: {.path {file}}")
  }

  if (is.null(cell_type)) {
    cell_type <- stringr::str_remove(
      basename(file),
      "_500kb_combined\\.MR\\.tsv\\.gz$"
    )
    if (identical(cell_type, basename(file))) {
      cli::cli_warn(
        "Could not infer cell type from filename; set {.arg cell_type} explicitly."
      )
    }
    cli::cli_inform("Inferred cell type: {.val {cell_type}}")
  }

  cli::cli_inform(
    "Reading dynamic_cseqtl data from {.path {basename(file)}}..."
  )

  dt <- data.table::fread(file, nThread = parallel::detectCores()) |>
    as.data.frame()

  if (nrow(dt) == 0L) {
    cli::cli_abort("File is empty: {.path {file}}")
  }

  dt <- dt |>
    dplyr::mutate(
      phenotype_label = paste(.data$gene, .env$cell_type, sep = "___")
    )

  TwoSampleMR::format_data(
    dat = dt,
    type = "exposure",
    phenotype_col = "phenotype_label",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "pval",
    chr_col = "chr",
    pos_col = "pos"
  )
}
