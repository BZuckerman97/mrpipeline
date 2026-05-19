#' Format GWAS summary statistics to a standard schema
#'
#' Normalises column names from any GWAS summary statistics file to the
#' canonical schema expected by [run_mr()] and [run_coloc()]. Handles
#' datasets where rsIDs are absent by looking them up from a PLINK bim file,
#' and supports extracting chromosome and position from compound marker ID
#' columns (e.g. SCALLOP `MarkerName` format: `"CHR:POS:A1_A2"`).
#'
#' @section Column normalisation:
#' The function renames source columns to a fixed canonical schema by checking
#' a built-in table of known aliases for each target column:
#'
#' | Canonical | Built-in aliases recognised automatically |
#' |-----------|------------------------------------------|
#' | `rsids`   | `rsid`, `rs_id`, `rsID`, `SNP` |
#' | `chr`     | `chromosome`, `Chr`, `CHROM`, `#CHROM`, `CHR` |
#' | `pos`     | `base_pair_location`, `PosB37`, `PosB38`, `BP`, `POS`, `position`, `GENPOS` |
#' | `beta`    | `Beta`, `Effect`, `BETA` |
#' | `se`      | `standard_error`, `StdErr`, `SE`, `sebeta` |
#' | `or`      | `odds_ratio`, `OR` |
#' | `eaf`     | `effect_allele_frequency`, `Freq1`, `EAFrq`, `A1FREQ`, `af_alt`, `EAF` |
#' | `pval`    | `p_value`, `P-value`, `P`, `Pval`, `p.value` |
#' | `n`       | `N`, `TotalSampleSize`, `n_total` |
#' | `effect_allele` | `Allele1`, `EA`, `A1`, `ALLELE1`, `effectAllele`, `ALT` |
#' | `other_allele`  | `Allele2`, `OA`, `A2`, `ALLELE0`, `otherAllele`, `REF` |
#'
#' Supply `col_map` **only** when your dataset uses a column name that does not
#' appear in the table above -- for example, if your p-value column is called
#' `"PVALUE"`, add `col_map = list(pval = "PVALUE")`. Inspect `names()` of your
#' loaded data to check. User-supplied aliases are checked before the built-in
#' list, so they take precedence in the event of ambiguity.
#'
#' @section Automatic odds-ratio to log-odds conversion:
#' Some GWAS files (particularly older EBI deposits) report effect sizes as odds
#' ratios rather than log-odds. When `beta` is absent after column normalisation
#' but an `or` column is present, the function automatically derives:
#'
#' - `beta = log(or)` (natural log; converts OR to the log-odds scale MR requires)
#' - `se = |beta| / qnorm(pval / 2)` (Z-score back-calculation from p-value)
#'
#' The Z-score method requires `pval` to be present and is accurate when effect
#' sizes are small (OR close to 1), which is typical for common-variant GWAS.
#' An informative message is emitted whenever the conversion is applied.
#' If your file has an OR column under a non-standard name, add it via
#' `col_map = list(or = "MY_OR_COLUMN")`.
#'
#' @section rsID lookup from bim file:
#' When `rsids` is absent (or all `NA`) after column normalisation, and
#' `bim_path` is supplied, the function inner-joins the data to the PLINK bim
#' file by chromosome and position to recover rsIDs. Rows without a bim match
#' are dropped -- they are absent from the reference panel and cannot be used
#' in LD-based analyses. A message reports how many SNPs were retained.
#'
#' @section Marker column parsing:
#' Set `marker_col` to the name of a compound marker ID column whose values
#' have the form `"CHR<sep>POS<sep>..."` (e.g. SCALLOP `"MarkerName"`).
#' `chr` and `pos` are extracted from the first two fields. This step runs
#' before the rsID lookup so that the extracted coordinates are available for
#' the bim join.
#'
#' @param path Character file path (`.tsv`, `.tsv.gz`, `.txt.gz`, etc. --
#'   `data.table::fread()` auto-detects compression) or a pre-loaded data frame.
#' @param phenotype_id Character. Trait / phenotype identifier (e.g. `"IL-18"`,
#'   `"CAD"`). Stored in the `phenotype` output column.
#' @param type `"outcome"` (default) or `"exposure"`. `"outcome"` returns the
#'   normalised data frame ready for [run_mr()]'s `outcome` argument.
#'   `"exposure"` additionally calls [TwoSampleMR::format_data()] and returns
#'   TwoSampleMR-formatted data with `exposure.` column suffixes.
#' @param col_map Named list of extra column-name aliases, e.g.
#'   `list(pval = "PVALUE", n = "SampleSize")`. Only needed for column names
#'   not already covered by the built-in alias table (see *Column normalisation*
#'   section). User entries are checked before the built-in list.
#' @param bim_path Character. Path to a PLINK bfile prefix (without `.bim`)
#'   used to recover rsIDs when the data lacks them. Required whenever the
#'   rsids column is absent.
#' @param marker_col Character. Name of a compound marker ID column in
#'   `"CHR<sep>POS<sep>..."` format (e.g. `"MarkerName"` for SCALLOP files).
#'   When supplied, `chr` and `pos` are parsed from this column.
#' @param marker_sep Character. Field separator used in `marker_col`. Default
#'   `":"`.
#' @param log10_pval Logical. If `TRUE`, the p-value column is in `-log10`
#'   scale and is back-transformed via `10^-x`. Default `FALSE`.
#' @param flip_beta Logical. If `TRUE`, multiplies `beta` by `-1` - use when
#'   the source file encodes the inverse direction of the intended exposure
#'   (e.g. modelling NLRP3 inhibition rather than activation). Default `FALSE`.
#' @param n Integer. Explicit sample size. Added as the `n` column only when no
#'   sample-size column is already present in the data.
#'
#' @return
#' - `type = "outcome"`: a data frame with columns `rsids`, `chr`, `pos`,
#'   `beta`, `se`, `eaf`, `pval`, `n`, `effect_allele`, `other_allele`,
#'   `phenotype` (plus any extra columns from the source file).
#' - `type = "exposure"`: a TwoSampleMR-formatted data frame (output of
#'   [TwoSampleMR::format_data()]) with `exposure.`-suffixed columns, suitable
#'   for [run_mr()]'s `exposure` argument.
#'
#' @examples
#' \dontrun{
#' # Outcome GWAS whose columns are already in the built-in alias table
#' cad <- format_gwas(
#'   path         = "genomics_data/outcome_GWAS/CAD/cad_gwas.tsv.gz",
#'   phenotype_id = "CAD"
#' )
#'
#' # SCALLOP outcome: no rsIDs, chr+pos embedded in MarkerName column
#' scallop_il6 <- format_gwas(
#'   path         = "genomics_data/outcome_GWAS/SCALLOP/CVD1_IL6.tsv.gz",
#'   phenotype_id = "IL6",
#'   marker_col   = "MarkerName",
#'   bim_path     = "LD_ref/g1000_eur"
#' )
#'
#' # Dataset with a non-standard p-value column not in the alias table
#' ebi_il18 <- format_gwas(
#'   path         = "genomics_data/outcome_GWAS/EBI/GCST90428399.tsv.gz",
#'   phenotype_id = "IL-18",
#'   col_map      = list(pval = "PVALUE")
#' )
#'
#' # Exposure GWAS -- flip beta to model NLRP3 activation not suppression
#' exposure <- format_gwas(
#'   path         = "NLRP3/Output/NLRP3_CRP_IVs_300kb.tsv",
#'   phenotype_id = "NLRP3",
#'   type         = "exposure",
#'   flip_beta    = TRUE
#' )
#' }
#'
#' @importFrom rlang .data .env
#' @export
format_gwas <- function(
  path,
  phenotype_id,
  type = c("outcome", "exposure"),
  col_map = NULL,
  bim_path = NULL,
  marker_col = NULL,
  marker_sep = ":",
  log10_pval = FALSE,
  flip_beta = FALSE,
  n = NULL
) {
  type <- match.arg(type)
  path_label <- if (is.character(path)) path else "<data frame>"

  # -- Read data ----------------------------------------------------------------
  if (is.character(path)) {
    if (!file.exists(path)) {
      cli::cli_abort("File not found: {.path {path}}")
    }
    dat <- as.data.frame(data.table::fread(path))
  } else if (is.data.frame(path)) {
    dat <- as.data.frame(path)
  } else {
    cli::cli_abort("{.arg path} must be a file path or data frame.")
  }

  # -- Parse chr + pos from compound marker ID column ---------------------------
  if (!is.null(marker_col)) {
    if (!marker_col %in% names(dat)) {
      cli::cli_abort(
        c(
          "{.arg marker_col} {.val {marker_col}} not found in data.",
          "i" = "Available columns: {.val {names(dat)}}"
        )
      )
    }
    parts <- stringr::str_split_fixed(
      dat[[marker_col]],
      stringr::fixed(marker_sep),
      n = 3L
    )
    dat <- dplyr::mutate(
      dat,
      chr = parts[, 1L],
      pos = suppressWarnings(as.integer(parts[, 2L]))
    )
  }

  # -- Column normalisation -----------------------------------------------------
  col_aliases <- list(
    rsids = c("rsids", "rsid", "rs_id", "rsID", "SNP"),
    chr = c("chr", "chromosome", "Chr", "CHROM", "#CHROM", "CHR"),
    pos = c(
      "pos",
      "base_pair_location",
      "PosB37",
      "PosB38",
      "BP",
      "POS",
      "position",
      "GENPOS"
    ),
    beta = c("beta", "Beta", "Effect", "BETA"),
    se = c("se", "standard_error", "StdErr", "SE", "sebeta"),
    or = c("or", "odds_ratio", "OR"),
    eaf = c(
      "eaf",
      "effect_allele_frequency",
      "Freq1",
      "EAFrq",
      "A1FREQ",
      "af_alt",
      "EAF"
    ),
    pval = c("pval", "p_value", "P-value", "P", "Pval", "p.value"),
    n = c("n", "N", "TotalSampleSize", "n_total"),
    effect_allele = c(
      "effect_allele",
      "Allele1",
      "EA",
      "A1",
      "ALLELE1",
      "effectAllele",
      "ALT"
    ),
    other_allele = c(
      "other_allele",
      "Allele2",
      "OA",
      "A2",
      "ALLELE0",
      "otherAllele",
      "REF"
    )
  )

  # User-supplied aliases are prepended so they are found before built-in ones
  if (!is.null(col_map)) {
    for (nm in names(col_map)) {
      col_aliases[[nm]] <- c(col_map[[nm]], col_aliases[[nm]])
    }
  }

  for (canonical in names(col_aliases)) {
    if (canonical %in% names(dat)) {
      next
    }
    found <- intersect(col_aliases[[canonical]], names(dat))
    if (length(found) > 0L) {
      dat <- dplyr::rename(
        dat,
        dplyr::all_of(stats::setNames(found[[1L]], canonical))
      )
    }
  }

  # -- Uppercase allele columns -------------------------------------------------
  dat <- dplyr::mutate(
    dat,
    dplyr::across(dplyr::any_of(c("effect_allele", "other_allele")), toupper)
  )

  # -- rsID lookup from bim file ------------------------------------------------
  has_rsids <- "rsids" %in% names(dat) && !all(is.na(dat[["rsids"]]))

  if (!has_rsids) {
    if (is.null(bim_path)) {
      cli::cli_abort(
        c(
          "No rsID column found for {.val {phenotype_id}} after column normalisation.",
          "i" = "Supply {.arg bim_path} to look up rsIDs by chromosome and position.",
          "i" = "Or add the source rsID column name to {.arg col_map}."
        )
      )
    }

    if (!all(c("chr", "pos") %in% names(dat))) {
      cli::cli_abort(
        c(
          "rsID lookup requires {.val chr} and {.val pos} columns, which are absent.",
          "i" = "Use {.arg marker_col} to parse them from a compound marker ID column."
        )
      )
    }

    bim_file <- paste0(bim_path, ".bim")
    if (!file.exists(bim_file)) {
      cli::cli_abort("bim file not found: {.path {bim_file}}")
    }

    bim <- data.table::fread(
      bim_file,
      header = FALSE,
      select = c(1L, 2L, 4L),
      col.names = c("chr", "rsids", "pos")
    ) |>
      as.data.frame() |>
      # One rsID per chr:pos -- multi-allelic sites share a position in the bim
      dplyr::distinct(.data$chr, .data$pos, .keep_all = TRUE) |>
      dplyr::mutate(
        chr = as.character(.data$chr),
        pos = as.integer(.data$pos)
      )

    # Drop any existing all-NA rsids column before joining to avoid name clash
    if ("rsids" %in% names(dat)) {
      dat <- dplyr::select(dat, -"rsids")
    }

    dat <- dplyr::mutate(
      dat,
      chr = as.character(.data$chr),
      pos = as.integer(.data$pos)
    )

    n_before <- nrow(dat)
    dat <- dplyr::inner_join(dat, bim, by = c("chr", "pos"))
    n_after <- nrow(dat)

    cli::cli_inform(
      "{.val {phenotype_id}}: {n_after}/{n_before} SNPs matched in {.path {bim_file}}."
    )
    if (n_after == 0L) {
      cli::cli_abort(
        "No SNPs matched in bim file. Check chromosome format and genome build."
      )
    }
  }

  # -- Normalise chromosome values ----------------------------------------------
  # Strip "chr" prefix if present (e.g. UCSC-format "chr1" -> "1", "chrX" -> "X").
  # Keeps chr as character so sex chromosomes (X, Y, MT) are preserved correctly.
  if ("chr" %in% names(dat)) {
    dat[["chr"]] <- sub("^chr", "", dat[["chr"]], ignore.case = FALSE)
  }

  # -- Coerce canonical numeric columns ----------------------------------------
  # fread reads entirely-NA or mixed-character columns as logical or character
  # (e.g. an eaf column that is wholly NA in an EBI harmonised file, or a pval
  # column that contains "NA" strings alongside numeric values). Coercing here
  # means TwoSampleMR::format_data() always receives the expected types and does
  # not emit "column is not numeric. Coercing..." warnings.
  for (col in intersect(c("beta", "se", "or", "eaf", "pval"), names(dat))) {
    dat[[col]] <- suppressWarnings(as.numeric(dat[[col]]))
  }
  for (col in intersect(c("pos", "n"), names(dat))) {
    dat[[col]] <- suppressWarnings(as.integer(dat[[col]]))
  }

  # -- Transformations ----------------------------------------------------------
  if (log10_pval && "pval" %in% names(dat)) {
    dat <- dplyr::mutate(dat, pval = 10^-.data$pval)
  }

  if (flip_beta && "beta" %in% names(dat)) {
    dat <- dplyr::mutate(dat, beta = -.data$beta)
  }

  # -- Derive beta + se from odds ratio when beta is absent ---------------------
  # Triggered when an or column exists but beta does not -- e.g. Rashkin 2020
  # cancer GWASes (NHL, melanoma) which publish odds_ratio + p_value only.
  # Formula: beta = log(OR);  se = |beta| / qnorm(pval / 2)  (Z-score method).
  if ("or" %in% names(dat) && !"beta" %in% names(dat)) {
    if (!"pval" %in% names(dat)) {
      cli::cli_abort(
        c(
          "{.val {phenotype_id}}: cannot derive {.val se} from OR without a p-value column.",
          "i" = "The Z-score method requires: se = |log(OR)| / qnorm(pval / 2).",
          "i" = "Supply a p-value column via {.arg col_map}."
        )
      )
    }
    if (!"se" %in% names(dat)) {
      dat <- dplyr::mutate(
        dat,
        beta = log(.data$or),
        se   = abs(.data$beta) / stats::qnorm(.data$pval / 2, lower.tail = FALSE)
      )
      cli::cli_inform(
        "{.val {phenotype_id}}: no beta/se columns found -- derived beta = log(OR) and se via Z-score method."
      )
    } else {
      dat <- dplyr::mutate(dat, beta = log(.data$or))
      cli::cli_inform(
        "{.val {phenotype_id}}: derived beta = log(OR); using existing se column."
      )
    }
  }

  if (!is.null(n) && !"n" %in% names(dat)) {
    dat <- dplyr::mutate(dat, n = as.integer(.env$n))
  }
  # TwoSampleMR::format_data() does `if (samplesize_col %in% names(dat))` which
  # errors when samplesize_col is NULL (NULL %in% x returns logical(0), not FALSE).
  # Guarantee an n column so we always pass a string, never NULL.
  if (!"n" %in% names(dat)) {
    dat$n <- NA_integer_
  }

  dat <- dplyr::mutate(dat, phenotype = .env$phenotype_id)

  # -- Validate required columns ------------------------------------------------
  required <- c("rsids", "beta", "se", "pval", "effect_allele", "other_allele")
  missing_cols <- setdiff(required, names(dat))
  if (length(missing_cols) > 0L) {
    cli::cli_abort(
      c(
        "{.val {phenotype_id}}: {length(missing_cols)} required column{?s} missing after normalisation.",
        "x" = "Missing: {.val {missing_cols}}",
        "i" = "Add the source column name(s) to {.arg col_map}.",
        "i" = "File: {.path {path_label}}"
      )
    )
  }

  # -- Exposure: return TwoSampleMR-formatted data ------------------------------
  if (type == "exposure") {
    return(TwoSampleMR::format_data(
      dat,
      type = "exposure",
      header = TRUE,
      phenotype_col = "phenotype",
      snp_col = "rsids",
      beta_col = "beta",
      se_col = "se",
      eaf_col = if ("eaf" %in% names(dat)) "eaf" else NULL,
      effect_allele_col = "effect_allele",
      other_allele_col = "other_allele",
      pval_col = "pval",
      chr_col = if ("chr" %in% names(dat)) "chr" else NULL,
      pos_col = if ("pos" %in% names(dat)) "pos" else NULL,
      samplesize_col = "n",
      log_pval = FALSE
    ))
  }

  # -- Outcome: return normalised data frame ------------------------------------
  dat
}
