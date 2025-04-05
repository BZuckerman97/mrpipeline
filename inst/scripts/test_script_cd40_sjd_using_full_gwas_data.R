sjd <- data.table::fread("/Users/benzuckerman/Library/Mobile Documents/com~apple~CloudDocs/Research/mr_work/data/sjogrens_6098_34928.txt")

sjd <- as.data.frame(sjd)
sjd <- sjd |>
  dplyr::mutate(phenotype == "SjD")

sjd <- TwoSampleMR::format_data(
  sjd,
  type = "outcome",
  phenotype_col = "phenotype"
  header = "TRUE",
  snp_col = "rsids",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  beta_col = "beta",
  se_col = "se",
  pval_col = "pval"
)

cd40 <- format_pqtl_ukbppp(ukbppp = "data/CD40_data_50.csv",
                   ukbppp_rsid = "data/CHR20_data_50.csv",
                   pqtl_assay = "CD40")

run_mr(exposure = cd40,
       exposure_id = "CD40",
       outcome = sjd,
       outcome_id = "SjD",
       instrument_region = list(chromosome = 20, start = 44746911, end = 44758502),
       pval_thresh = 5e-06,
       rsq_thresh = 0.1,
       bfile = "/Users/benzuckerman/LD_folder/g1000_eur")
