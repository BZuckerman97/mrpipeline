# -- mrpipeline single-cell eQTL analysis: configuration ---------------------
#
# Edit all paths and parameters in this file before running the pipeline.
# Source this file at the top of _targets.R; do not run it directly.

cfg <- list(
  # -- OneK1K ----------------------------------------------------------------
  # esnp_table.tsv    : significant eQTLs only (small, used for MR instruments)
  # eqtl_table.tsv.gz : full cis summary stats (21 GB, streamed per gene for coloc)
  onek1k_esnp_file = "/Users/Benjamin/Documents/Research/MR_projects/genomics_data/single cell/onek1k/esnp_table.tsv",
  onek1k_eqtl_file = "/Users/Benjamin/Documents/Research/MR_projects/genomics_data/single cell/onek1k/eqtl_table.tsv.gz",

  # Cell types present in the OneK1K data (all 14)
  onek1k_cell_types = c(
    "bin",
    "bmem",
    "cd4et",
    "cd4nc",
    "cd4sox4",
    "cd8et",
    "cd8nc",
    "cd8s100b",
    "dc",
    "monoc",
    "mononc",
    "nk",
    "nkr",
    "plasma"
  ),

  # Number of cells per cell type -- used as samplesize.exposure for coloc.
  # Fill in from Yazar et al. 2022 Science supplementary tables.
  # Only cd4nc is confirmed from the test run; others are NA until verified.
  onek1k_n = c(
    bin = NA_integer_,
    bmem = NA_integer_,
    cd4et = NA_integer_,
    cd4nc = 463528L,
    cd4sox4 = NA_integer_,
    cd8et = NA_integer_,
    cd8nc = NA_integer_,
    cd8s100b = NA_integer_,
    dc = NA_integer_,
    monoc = NA_integer_,
    mononc = NA_integer_,
    nk = NA_integer_,
    nkr = NA_integer_,
    plasma = NA_integer_
  ),

  # -- 1M-scBloodNL ----------------------------------------------------------
  # Per-cell-type FDR-level eQTL files (one file per cell type per condition).
  # File naming: <cell_type>_expression_eQTLsFDR-ProbeLevel.txt.gz
  #
  # NOTE: these files contain only FDR-significant cis associations, not the
  # full cis summary statistics. MR is valid, but coloc requires full summary
  # stats (not available in the public FDR-level files) -- coloc is therefore
  # skipped for 1M-scBloodNL in this pipeline.
  scbloodnl_dir = "/Users/Benjamin/Documents/Research/MR_projects/genomics_data/single cell/1M-scBloodNL/sc_eqtls_20201106_genome_wide",
  scbloodnl_condition = "UT", # subdirectory: UT, 3hCA, 3hMTB, 3hPA, 24hCA, etc.

  # Cell types available in the UT condition (file stems before
  # "_expression_eQTLsFDR-ProbeLevel.txt.gz").
  # Fill in sample sizes from Oelen et al. 2022 supplementary tables.
  scbloodnl_cell_types = c(
    "B",
    "CD4T",
    "CD8T",
    "DC",
    "NK",
    "monocyte",
    "megakaryocyte"
  ),
  scbloodnl_n = c(
    B = NA_integer_,
    CD4T = NA_integer_,
    CD8T = NA_integer_,
    DC = NA_integer_,
    NK = NA_integer_,
    monocyte = NA_integer_,
    megakaryocyte = NA_integer_
  ),

  # -- DICE (template -- no data locally yet) --------------------------------
  # dice_dir = "/path/to/dice/vcf_files",
  # dice_cell_types = c("t_cell_cd4_naive", "t_cell_cd8_naive", "monocyte"),
  # dice_n = c(t_cell_cd4_naive = 91L, t_cell_cd8_naive = 91L, monocyte = 91L),

  # -- dynamic_cseqtl (template -- no data locally yet) ---------------------
  # cseqtl_dir = "/path/to/dynamic_cseqtl",
  # cseqtl_cell_types = c("CD4T", "CD8T", "NK"),
  # cseqtl_n = c(CD4T = NA_integer_, CD8T = NA_integer_, NK = NA_integer_),

  # -- LD reference ----------------------------------------------------------
  bfile = "/Users/Benjamin/Documents/Research/MR_projects/LD_ref/g1000_eur",

  # -- Outcomes --------------------------------------------------------------
  # Add further outcomes by copying the ra_eur block and adjusting values.
  outcomes = list(
    ra_eur = list(
      file = "/Users/Benjamin/Documents/Research/MR_projects/genomics_data/outcome_GWAS/RA_Ishigaki2022/GCST90132223_RA_EUR_seropos_seroneg/GCST90132223_buildGRCh37.tsv.gz",
      id = "RA_EUR",
      n = 36458L,
      n_cases = 14551L,
      type = "cc",
      s = 14551L / 36458L # proportion of cases
    )
  ),

  # -- MR parameters ---------------------------------------------------------
  # pval_thresh: instrument selection p-value. Note: esnp_table contains
  # FDR-significant eQTLs (not all genome-wide significant), so you may want
  # to set this to 1e-5 or even 1 and rely on clumping alone.
  mr_pval_thresh = 5e-8,
  mr_rsq_thresh = 0.001,
  mr_window = 100000L, # bp either side of cis region for instrument search
  mr_methods = c("ivw", "egger", "weighted_median"),

  # Exclude MHC from instruments (chr6:26-34 Mb)
  mhc_region = data.frame(chr = "6", start = 26e6, end = 34e6),

  # -- Coloc parameters ------------------------------------------------------
  # Run coloc only for gene x cell_type pairs where any MR method gives p < this
  coloc_mr_pval = 0.05,
  coloc_window = 10000L,
  coloc_methods = c("abf", "susie", "signals"),

  # -- Output ----------------------------------------------------------------
  output_dir = "output"
)
