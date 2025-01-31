run_mr(exposure = cd40_sumstats,
       exposure_id = "CD40",
       outcome = sjogren_sumstats,
       outcome_id = "SjD",
       instrument_region = list(chromosome = 1L, start = 1L, end = 200L),
       pval_thresh = 5e-06,
       rsq_thresh = 0.1)
