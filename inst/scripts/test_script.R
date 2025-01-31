run_mr(expoure = cd40_sumstats,
       exposure_id = "exposure",
       outcome = sjogren_sumstats,
       outcome_id = "outcome",
       instrument_region = list(chromosome = 1L, start = 1L, end = 200L),
       pval_thresh = 5e-06,
       rsq_thresh = 0.1)
