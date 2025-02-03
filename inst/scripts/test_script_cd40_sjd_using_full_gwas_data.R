#sjd <- data.table::fread("/Users/benzuckerman/Library/Mobile Documents/com~apple~CloudDocs/Research/mr_work/data/sjogrens_6098_34928.txt")

#cd40 <- data.table::fread("/Users/benzuckerman/Downloads/cd40_rsids_olink_h37.csv")

#sjd_1 <- sjd |>
#  dplyr::mutate(chrom = chr,
#                rsids = SNP,
#                alt = effect_allele,
#                ref = other_allele,
#                sebeta = se,
#                af_alt = eaf) |>
#  dplyr::mutate(phenotype = "SjD")

#sjd_1 <- as.data.frame(sjd_1)
#cd40 <- as.data.frame(cd40)

#run_mr(exposure = cd40,
#       exposure_id = "CD40",
#       outcome = sjd_1,
#       outcome_id = "SjD",
#       instrument_region = list(chromosome = 20, start = 44746911, end = 44758502),
#       pval_thresh = 5e-06,
#       rsq_thresh = 0.1,
#       bfile = "/Users/benzuckerman/LD_folder/g1000_eur")
