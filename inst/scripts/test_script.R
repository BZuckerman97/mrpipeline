run_mr(exposure = cd40_sumstats,
       exposure_id = "CD40",
       outcome = sjogren_sumstats,
       outcome_id = "SjD",
       instrument_region = list(chromosome = 20, start = 44746911, end = 44758502),
       pval_thresh = 5e-06,
       rsq_thresh = 0.1,
       bfile = "/Users/benzuckerman/Library/Mobile Documents/com~apple~CloudDocs/Research/mr_work/LD_ref/")
