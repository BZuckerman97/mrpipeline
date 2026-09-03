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
#' | `zscore`  | `zscore`, `Zscore`, `Z-score`, `Z_score` |
#' | `eaf`     | `effect_allele_frequency`, `Freq1`, `EAFrq`, `A1FREQ`, `af_alt`, `EAF`, `FRQ` |
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
#' ratios rather than log-odds. When `beta` is absent -- or present but entirely
#' `NA` (some GWAS Catalog templates carry an empty `beta` column alongside a
#' populated `odds_ratio` column, e.g. Ji 2016 PSC) -- after column
#' normalisation, and an `or` column is present, the function automatically
#' derives:
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
#' @section Automatic Z-score to beta/se conversion:
#' Some GWAS files (typically METAL output) report only a Z-score rather than
#' beta/se -- e.g. Hysi 2020 refractive error, He strabismus, Cuellar-Partida
#' 2021 handedness. When `beta` is absent after column normalisation but a
#' `zscore` column is present, the function derives it automatically:
#'
#' - If `se` is already present: `beta = zscore * se`.
#' - Otherwise, `se` is first derived from allele frequency and sample size
#'   via the standard GWAS Z-score formula, then `beta = zscore * se`:
#'   `se = 1 / sqrt(2 * eaf * (1 - eaf) * (n + zscore^2))`. This path requires
#'   both `eaf` and `n` columns.
#'
#' An informative message is emitted whenever the conversion is applied. If
#' your file has a Z-score column under a non-standard name, add it via
#' `col_map = list(zscore = "MY_Z_COLUMN")`; likewise use `col_map` for a
#' non-standard sample-size column (e.g. METAL's `Weight`), via
#' `col_map = list(n = "Weight")`.
#'
#' @section Automatic se derivation from beta + zscore:
#' Some GWAS files report a real `beta` and a `zscore` but no `se` column --
#' e.g. the D-dimer Suhre GCST90100910 file, which publishes `Beta` + `Z`
#' (its `z` column needs `col_map = list(zscore = "z")`, since a bare `"z"`
#' is not in the built-in `zscore` alias table). When this happens, the
#' function derives `se` as an exact algebraic inversion of `z = beta / se`:
#'
#' - `se = |beta / zscore|`
#'
#' This is preferred over the beta+p-value back-calculation below whenever a
#' `zscore` column is available, since it does not depend on a (possibly
#' rounded) p-value -- the p-value method is undefined when `pval == 1`, a
#' real risk for a true-null SNP. An informative message is emitted whenever
#' this derivation is applied.
#'
#' @section Automatic se derivation from beta + p-value:
#' Some GWAS files report a real `beta` and `pval` but no `se`, `or`, or
#' `zscore` column at all -- e.g. the deCODE haematinics files (Iron_TSAT,
#' Iron_Ferritin). When this happens (and the OR/Z-score derivations above
#' don't apply, since neither an `or` nor a `zscore` column is present), the
#' function derives `se` directly via the same Z-score back-calculation used
#' in the OR conversion above:
#'
#' - `se = |beta| / qnorm(pval / 2)`
#'
#' An informative message is emitted whenever this derivation is applied. If
#' your file has `se` under a non-standard name, add it via
#' `col_map = list(se = "MY_SE_COLUMN")` instead of relying on this fallback.
#'
#' @section rsID lookup from bim file:
#' When `rsids` is absent (or all `NA`) after column normalisation, and
#' `bim_path` is supplied, the function inner-joins the data to the PLINK bim
#' file by chromosome and position to recover rsIDs. Rows without a bim match
#' are dropped -- they are absent from the reference panel and cannot be used
#' in LD-based analyses. A message reports how many SNPs were retained.
#'
#' @section chr/pos lookup from bim file (reverse direction):
#' Some files report a real rsID but no positional columns at all (e.g. the
#' Corbin & Timpson 2020 cytokine GWASes, whose `MarkerName` column already
#' contains rsIDs). When `rsids` is present but `chr`/`pos` are absent, and
#' `bim_path` is supplied, the function instead inner-joins to the bim file
#' **by rsID** to recover chromosome and position. This is the mirror image of
#' the lookup above and is needed because [run_coloc()] requires `chr`/`pos`
#' to window data around a gene region, even though [run_mr()] does not need
#' them. If `bim_path` is not supplied in this situation, `chr`/`pos` are
#' simply left absent -- `format_gwas()` does not require them for `run_mr()`,
#' but a coloc analysis on that outcome will fail downstream until they are
#' added (e.g. by re-running with `bim_path` set).
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
#'   used to recover rsIDs when the data lacks them (joins by chr + pos), or
#'   to recover chr/pos when the data has rsIDs but no positional columns
#'   (joins by rsID) -- see *rsID lookup from bim file* and *chr/pos lookup
#'   from bim file (reverse direction)*. Required whenever the rsids column
#'   is absent; optional but recommended when rsids are present without
#'   chr/pos and a coloc analysis is planned.
#' @param marker_col Character. Name of a compound marker ID column in
#'   `"CHR<sep>POS<sep>..."` format (e.g. `"MarkerName"` for SCALLOP files).
#'   When supplied, `chr` and `pos` are parsed from this column.
#' @param marker_sep Character. Field separator used in `marker_col`. Default
#'   `":"`.
#' @param log10_pval Logical. If `TRUE`, the p-value column is in `-log10`
#'   scale and is back-transformed via `10^-x`. Default `FALSE`.
#' @param n Integer. Explicit sample size. Added as the `n` column only when no
#'   sample-size column is already present in the data.
#' @param ref_frq Character. Path to a PLINK `.frq` file (produced by
#'   `plink --freq`, columns `CHR SNP A1 A2 MAF NCHROBS`). When supplied,
#'   patches `eaf` for any row where it is missing (creating the `eaf` column
#'   first if the data has none at all) by matching `rsids` against `SNP` and
#'   aligning `MAF`/`1-MAF` to `effect_allele`. Existing non-missing `eaf`
#'   values are never overwritten. Useful for pre-clumped instrument files
#'   that carry no allele-frequency column, where palindromic SNPs would
#'   otherwise be unresolvable (and typically dropped) during harmonisation.
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
#' # Exposure GWAS, ready for run_mr()
#' exposure <- format_gwas(
#'   path         = "NLRP3/Output/NLRP3_CRP_IVs_300kb.tsv",
#'   phenotype_id = "NLRP3",
#'   type         = "exposure"
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
  n = NULL,
  ref_frq = NULL
) {
  type <- match.arg(type)
  path_label <- if (is.character(path)) path else "<data frame>" # nolint: object_usage_linter.

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
    zscore = c("zscore", "Zscore", "Z-score", "Z_score"),
    eaf = c(
      "eaf",
      "effect_allele_frequency",
      "Freq1",
      "EAFrq",
      "A1FREQ",
      "af_alt",
      "EAF",
      "FRQ"
    ),
    pval = c(
      "pval",
      "p_value",
      "P-value",
      "P",
      "Pval",
      "p.value",
      "neg_log_10_p_value"
    ),
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

  # Detect neg_log10 p-value before any renaming so we can back-transform later
  auto_neg_log_pval <- "neg_log_10_p_value" %in%
    names(dat) &&
    !"pval" %in% names(dat)

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

  # -- Clean compound "rsID:allele:allele" values in rsids column ---------------
  # Some consortium files (e.g. GIANT_UKBB_2018 BMI/WHR/WHRadjBMI) name a
  # compound "rsID:allele:allele" column "SNP", which the built-in alias
  # table above auto-maps to `rsids`. Left as-is, these compound strings
  # never match a clean instrument rsID during harmonisation -- discovered
  # when GIANT_UKBB_2018 branches fell through to "no overlapping SNPs"
  # despite the raw file having full coverage in the queried region. Strip
  # the trailing ":allele:allele" suffix wherever present so rsids is always
  # a bare rsID; rows that already carry a bare rsID (a minority even within
  # GIANT's own file, e.g. some indel sites) pass through unchanged. A
  # simple strip is used rather than routing through the bim lookup below,
  # since that would needlessly discard any SNP absent from the 1000G
  # reference panel even though its native rsID is perfectly usable.
  if ("rsids" %in% names(dat)) {
    compound <- grepl("^rs[0-9]+:", dat[["rsids"]])
    if (any(compound, na.rm = TRUE)) {
      n_compound <- sum(compound, na.rm = TRUE) # nolint: object_usage_linter.
      dat[["rsids"]][compound] <- sub(
        "^(rs[0-9]+):.*$",
        "\\1",
        dat[["rsids"]][compound]
      )
      cli::cli_inform(paste0(
        "{.val {phenotype_id}}: cleaned {n_compound} compound ",
        "\"rsID:allele:allele\" value{?s} in the rsids column to bare rsIDs."
      ))
    }
  }

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

    n_before <- nrow(dat) # nolint: object_usage_linter.
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
  } else if (!all(c("chr", "pos") %in% names(dat)) && !is.null(bim_path)) {
    # -- chr/pos lookup from bim file (via rsID) --------------------------------
    # Triggered when rsids are already present but chr/pos are absent -- e.g.
    # Corbin & Timpson 2020 cytokine files, whose MarkerName column is a real
    # rsID but the file carries no positional columns at all. run_mr() does not
    # need chr/pos (it works by rsID alone), but run_coloc() requires them to
    # window the data around the target gene region, so this recovers them by
    # joining on rsID instead of chromosome + position.
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
      # One row per rsID -- guards against any duplicate rsIDs in the bim
      dplyr::distinct(.data$rsids, .keep_all = TRUE) |>
      dplyr::mutate(
        chr = as.character(.data$chr),
        pos = as.integer(.data$pos)
      )

    n_before <- nrow(dat)
    dat <- dplyr::inner_join(dat, bim, by = "rsids")
    n_after <- nrow(dat)

    cli::cli_inform(
      "{.val {phenotype_id}}: {n_after}/{n_before} SNPs matched chr/pos via rsID in {.path {bim_file}}."
    )
    if (n_after == 0L) {
      cli::cli_abort(
        "No SNPs matched in bim file by rsID. Check rsID format and genome build."
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
  for (col in intersect(
    c("beta", "se", "or", "zscore", "eaf", "pval"),
    names(dat)
  )) {
    dat[[col]] <- suppressWarnings(as.numeric(dat[[col]]))
  }
  for (col in intersect(c("pos", "n"), names(dat))) {
    dat[[col]] <- suppressWarnings(as.integer(dat[[col]]))
  }

  # -- Patch missing EAF from a PLINK .frq reference file -----------------------
  # Same lookup run_coloc() already does inline for its own harmonised data
  # (see ref_frq there) -- offered here too because some sources have no EAF
  # at all (e.g. a pre-clumped instrument file that only carries beta/se/pval),
  # which otherwise leaves palindromic SNPs unresolvable during harmonisation.
  # .frq columns: CHR, SNP, A1, A2, MAF, NCHROBS. A1 is the frq file's coded
  # allele: EAF = MAF when effect_allele == A1, else EAF = 1 - MAF. Only fills
  # rows where eaf is missing; never overwrites an EAF already in the data.
  if (!is.null(ref_frq)) {
    if (!file.exists(ref_frq)) {
      cli::cli_abort("{.arg ref_frq} file not found: {.path {ref_frq}}")
    }
    if (!"eaf" %in% names(dat)) {
      dat$eaf <- NA_real_
    }
    if (all(c("rsids", "effect_allele") %in% names(dat))) {
      needs_eaf <- is.na(dat$eaf)
      if (any(needs_eaf)) {
        frq <- data.table::fread(ref_frq, data.table = FALSE)
        frq <- frq[!duplicated(frq$SNP), c("SNP", "A1", "MAF")]
        idx <- match(dat$rsids[needs_eaf], frq$SNP)
        found <- !is.na(idx)
        if (any(found)) {
          rows_need <- which(needs_eaf)[found]
          frq_matched <- frq[idx[found], ]
          eaf_ref <- ifelse(
            toupper(frq_matched$A1) == toupper(dat$effect_allele[rows_need]),
            frq_matched$MAF,
            1 - frq_matched$MAF
          )
          dat$eaf[rows_need] <- eaf_ref
        }
        cli::cli_inform(paste0(
          "{.val {phenotype_id}}: patched EAF for {sum(found)} of ",
          "{sum(needs_eaf)} SNP{?s} missing it, from {.path {ref_frq}}."
        ))
      }
    }
  }

  # -- Transformations ----------------------------------------------------------
  if (log10_pval && "pval" %in% names(dat)) {
    dat <- dplyr::mutate(dat, pval = 10^-.data$pval)
  }
  # Auto back-transform when the source column was neg_log_10_p_value
  if (auto_neg_log_pval && "pval" %in% names(dat)) {
    dat <- dplyr::mutate(dat, pval = 10^-.data$pval)
  }

  # -- Derive beta + se from odds ratio when beta is absent ---------------------
  # Triggered when an or column exists but beta does not -- e.g. Rashkin 2020
  # cancer GWASes (NHL, melanoma) which publish odds_ratio + p_value only.
  # Also triggered when a beta column exists but is entirely NA -- e.g. Ji 2016
  # PSC GWAS, whose template carries an empty beta column alongside a populated
  # odds_ratio column. Without this, the all-NA column would satisfy
  # `"beta" %in% names(dat)` and silently skip derivation, leaving every row's
  # effect estimate as NA rather than erroring.
  # Formula: beta = log(OR);  se = |beta| / qnorm(pval / 2)  (Z-score method).
  beta_absent <- !"beta" %in% names(dat) || all(is.na(dat[["beta"]]))
  if ("or" %in% names(dat) && beta_absent) {
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
        se = abs(.data$beta) / stats::qnorm(.data$pval / 2, lower.tail = FALSE)
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

  # -- Derive beta + se from Z-score when beta is absent ------------------------
  # Triggered when a zscore column exists but beta does not -- e.g. Cuellar-
  # Partida 2021 handedness (zscore + se already given) or Hysi 2020 refractive
  # error / He strabismus METAL output (zscore + eaf + n, no se at all).
  if ("zscore" %in% names(dat) && !"beta" %in% names(dat)) {
    if ("se" %in% names(dat)) {
      dat <- dplyr::mutate(dat, beta = .data$zscore * .data$se)
      cli::cli_inform(
        "{.val {phenotype_id}}: no beta column found -- derived beta = zscore * se."
      )
    } else if (all(c("eaf", "n") %in% names(dat))) {
      dat <- dplyr::mutate(
        dat,
        se = 1 /
          sqrt(2 * .data$eaf * (1 - .data$eaf) * (.data$n + .data$zscore^2)),
        beta = .data$zscore * .data$se
      )
      cli::cli_inform(
        "{.val {phenotype_id}}: no beta/se columns found -- derived se from eaf + n, and beta = zscore * se."
      )
    } else {
      cli::cli_abort(
        c(
          "{.val {phenotype_id}}: cannot derive {.val beta} from zscore.",
          "i" = "Need either an {.val se} column, or both {.val eaf} and {.val n} columns.",
          "i" = "Supply the missing column name(s) via {.arg col_map}."
        )
      )
    }
  }

  # -- Derive se from beta + zscore when se is absent ---------------------------
  # Triggered when beta and zscore are both present but se is not -- e.g. the
  # D-dimer Suhre GCST90100910 file, which publishes Beta + Z with no SE
  # column at all (col_map = list(zscore = "z") needed, since a bare "z"
  # column is not in the built-in zscore alias table). Preferred over the
  # beta+pval back-calculation below whenever a zscore column is available,
  # since it is an exact algebraic inversion of z = beta / se rather than a
  # back-calculation from a (possibly rounded) p-value -- which is undefined
  # when pval == 1, a real risk for any true-null SNP in the file.
  # Formula: se = |beta / zscore|.
  if (
    "beta" %in%
      names(dat) &&
      !all(is.na(dat[["beta"]])) &&
      "zscore" %in% names(dat) &&
      !"se" %in% names(dat)
  ) {
    dat <- dplyr::mutate(dat, se = abs(.data$beta / .data$zscore))
    cli::cli_inform(
      "{.val {phenotype_id}}: no se column found -- derived se = |beta / zscore|."
    )
  }

  # -- Derive se from beta + pval when se is absent -----------------------------
  # Triggered when beta and pval are both present but se, or, and zscore are
  # not -- e.g. the deCODE haematinics files (Iron_TSAT, Iron_Ferritin), which
  # publish Beta + P with no SE/OR/Z-score column at all. Same back-calculation
  # formula as the OR-derivation block above.
  if (
    "beta" %in%
      names(dat) &&
      !all(is.na(dat[["beta"]])) &&
      !"zscore" %in% names(dat) &&
      !"se" %in% names(dat)
  ) {
    if (!"pval" %in% names(dat)) {
      cli::cli_abort(
        c(
          paste0(
            "{.val {phenotype_id}}: cannot derive {.val se} -- no se/or/",
            "zscore column, and no pval to back-calculate from."
          ),
          "i" = paste0(
            "Supply an se column via {.arg col_map}, or a p-value column ",
            "so se can be derived: se = |beta| / qnorm(pval / 2)."
          )
        )
      )
    }
    dat <- dplyr::mutate(
      dat,
      se = abs(.data$beta) / stats::qnorm(.data$pval / 2, lower.tail = FALSE)
    )
    cli::cli_inform(
      "{.val {phenotype_id}}: no se column found -- derived se from beta + pval via Z-score method."
    )
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
