#' This test script needs to be altered to incorporate format_pqtl_ukbppp() and run_mr()
#' The main thing I need to do is get the cd40_sumstats set up appropriately for format_pqtl_ukbppp


run_mr(exposure = cd40_sumstats,
       exposure_id = "CD40",
       outcome = sjogren_sumstats,
       outcome_id = "SjD",
       instrument_region = list(chromosome = 20, start = 44746911, end = 44758502),
       pval_thresh = 5e-06,
       rsq_thresh = 0.1,
       bfile = "~/LD_folder/g1000_eur")
