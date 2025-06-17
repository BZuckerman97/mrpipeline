# Colocalization ----------------------------------------------------------
library(coloc)
library(data.table)
library(dplyr)
library(TwoSampleMR)

# Assuming these files exist and are accessible
listformr <- fread("output/decode_sjd_proteins_robust_from_egger_and_bhp_0103225.csv")
decode_sumstats <- fread("data/decode_sumstats_weij_bz_comb.csv")
incl_decode_variants <- fread("~/Shared2/rmgpbzu/deCODE/assocvariants.annotated.txt.gz") #' Add in effectAllelefreq
outcome <- fread("data/h38_sjogrens_6098_34928.tsv.gz") # Assuming this file exists
ref_rsid <- fread("hg38_common_chrpos_X.txt") # Assuming this file exists

# Filter for IL11RA
decode_sumstats_filter <- decode_sumstats %>%
  filter(Symbol %in% listformr$exp) %>%
  filter(Symbol == "IL11RA")

if (nrow(decode_sumstats_filter) == 0) {
  stop("No data found for IL11RA in decode_sumstats_filter.")
}

i <- decode_sumstats_filter$seqID

print(i)
timestamp()
print("******")

exposure <- fread(paste0("~/Shared2/rmgpbzu/deCODE/", decode_sumstats$identifier[decode_sumstats$seqID == i]))
exposure <- exposure |>
  inner_join(incl_decode_variants %>% select(Name, effectAlleleFreq), by = "Name")

exposure <- exposure %>%
  mutate(phen = decode_sumstats$proteinID[decode_sumstats$seqID == i])

#' Cleaning exposure

exposure <- exposure %>%
  rename(P = Pval)

exposure <- exposure %>%
  rename(
    chr = Chrom,
    pos = Pos,
    effect_allele = effectAllele,
    other_allele = otherAllele,
    beta = Beta,
    se = SE,
    eaf = effectAlleleFreq
  ) %>%
  mutate(chr = gsub("chr", "", chr))

exposure <- exposure %>%
  dplyr::filter(chr == decode_sumstats[decode_sumstats$seqID == i,]$chr)

exposure <- exposure[exposure$pos > (decode_sumstats[decode_sumstats$seqID == i,]$gene_start[1] - 200000) &
                       exposure$pos < (decode_sumstats[decode_sumstats$seqID ==
                                                         i,]$gene_end[1] + 200000), ]

outcome <- outcome |>
  mutate(phen = "SjD") |>
  rename(
    rsids = SNP,
    chr = CHR,
    sebeta = SE,
    af_alt = FRQ,
    alt = A2,
    ref = A1,
    beta = BETA,
    pval = P,
    pos = BP
  )

outcome_overlap <- outcome %>%
  dplyr::filter(chr %in% exposure$chr)
outcome_overlap <- outcome_overlap %>%
  dplyr::filter(pos %in% exposure$pos)

outcome_overlap <- as.data.frame(outcome_overlap)

outcome_overlap <- format_data(outcome_overlap,
                               type = "outcome",
                               phenotype_col = "phen",
                               snp_col = "rsids",
                               beta_col = "beta",
                               se_col = "sebeta",
                               eaf_col = "af_alt",
                               effect_allele_col = "alt",
                               other_allele_col = "ref",
                               pval_col = "pval",
                               chr_col = "chr",
                               pos_col = "pos")

exposure <- as.data.frame(exposure)
exposure <-
  format_data(
    exposure,
    type = "exposure",
    phenotype_col = "phen",
    snp_col = "rsids",
    beta_col = "beta",
    se_col = "se",
    eaf_col = "eaf",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "P",
    chr_col = "chr",
    samplesize_col = "N",
    pos_col = "pos"
  )

dat <- harmonise_data(exposure_dat = exposure, outcome_dat = outcome_overlap)

if (decode_sumstats_filter[decode_sumstats_filter$seqID == i,]$chr[1] == "X") {
  dat <- merge(dat, ref_rsid[, c("V1", "V2", "V3")], by.x = "pos.exposure", by.y = "V2", all.x = T)
  colnames(dat)[colnames(dat) %in% c("V1", "V3")] <- c("chr", "SNP")
}

dat <- dat[order(dat$pval.exposure),]
dat <- dat[!duplicated(dat$SNP),]
ld <- ld_matrix(dat$SNP, bfile = "LD_ref/g1000_eur", plink_bin = genetics.binaRies::get_plink_binary())
rownames(ld) <- gsub("_.*", "", rownames(ld))
colnames(ld) <- gsub("_.*", "", colnames(ld))
dat <- dat[dat$SNP %in% rownames(ld),]
dat <- dat[match(rownames(ld), dat$SNP),]

dat_exp <- list(beta = dat$beta.exposure, varbeta = dat$se.exposure^2, snp = dat$SNP, position = dat$pos.exposure, type = "quant", sdY = 1, LD = ld, N = 35350)
dat_outc <- list(beta = dat$beta.outcome, varbeta = dat$se.outcome^2, snp = dat$SNP, position = dat$pos.outcome, type = "cc", LD = ld, N = 41420)

s1 <- runsusie(dat_exp)
s2 <- runsusie(dat_outc)

susie_res <- coloc.susie(s1, s2)

coloc_outc <- coloc::coloc.abf(dataset1 = dat_exp, dataset2 = dat_outc)

results <- data.frame(
  exp = decode_sumstats_filter[decode_sumstats_filter$seqID == i,]$Symbol[1],
  outc = "SjD",
  nsnps = coloc_outc$summary["nsnps"],
  pp_h0 = coloc_outc$summary["PP.H0.abf"],
  pp_h1 = coloc_outc$summary["PP.H1.abf"],
  pp_h2 = coloc_outc$summary["PP.H2.abf"],
  pp_h3 = coloc_outc$summary["PP.H3.abf"],
  pp_h4 = coloc_outc$summary["PP.H4.abf"]
)

# Initialize Data Frames Correctly
if (!exists("df_sum")) {
  df_sum <- data.frame(exp = character(), outc = character(), nsnps = numeric(),
                       pp_h0 = numeric(), pp_h1 = numeric(), pp_h2 = numeric(),
                       pp_h3 = numeric(), pp_h4 = numeric())
}

df_sum <- rbind(df_sum, results)

write.csv(df_sum, paste0("output/coloc_h38_pqtl_sjd_decode_", Sys.Date(), ".csv"), row.names = FALSE)
rm(results,  dat, dat_exp, dat_outc, coloc_outc, ld)
print(paste(decode_sumstats[decode_sumstats$seqID == i,]$Symbol[1], "done"))
