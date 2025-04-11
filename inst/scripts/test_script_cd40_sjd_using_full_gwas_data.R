sjd <- data.table::fread("/Users/benzuckerman/Library/Mobile Documents/com~apple~CloudDocs/Research/mr_work/data/sjogrens_6098_34928.txt")

cd40 <- format_pqtl_ukbppp(ukbppp = "data/CD40_data_50.csv",
                   ukbppp_rsid = "data/CHR20_data_50.csv",
                   pqtl_assay = "CD40")

#' Need to find 50 cd40 SNPs which are genome wide significant

run_mr(exposure = cd40,
       exposure_id = "CD40",
       outcome = sjd,
       outcome_id = "SjD",
       instrument_region = list(chromosome = 20, start = 44746911, end = 44758502),
       pval_thresh = 5e-6,
       rsq_thresh = 0.1,
       bfile = "/Users/benzuckerman/LD_folder/g1000_eur")
