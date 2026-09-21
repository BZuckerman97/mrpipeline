# Generate the bundled LD reference panel: inst/extdata/ld_ref.{bed,bim,fam}
#
# Real 1000 Genomes phase 3 (release 20130502, GRCh37) genotypes for the 503
# EUR individuals at the 50 CD40-region SNPs in `cd40_exposure` (which are
# also the SNPs in `sjogren_outcome`). 1000 Genomes data are open access.
#
#   VCF:     <base>/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.
#            20130502.genotypes.vcf.gz (queried remotely by tabix; only the
#            ~120 kb CD40 slice is transferred)
#   Samples: <base>/integrated_call_samples_v3.20130502.ALL.panel
#   <base> = https://s3.amazonaws.com/1000genomes/release/20130502
#
# Allele handling: for each SNP the VCF ALT whose {REF, ALT} set equals the
# exposure's {effect, other} allele set is kept. rs13045469 is multi-allelic
# in 1000 Genomes (REF G, ALT A,C; EUR frequencies ~0.66/0.22/0.12) and the
# exposure describes it as C/G, so it is coded C vs not-C: the A allele is
# written as G. That gives the exact LD of the C-allele dosage. (Setting A
# carriers to missing instead drops 40% of individuals and leaves
# rs35922919 monomorphic among the rest, i.e. an NaN r.) No genotype is
# missing. A1 is PLINK's default (the minor allele in this sample).
#
# Requirements (deliberately not in DESCRIPTION -- data-raw/ is
# build-ignored): VariantAnnotation, Rsamtools, GenomicRanges, IRanges,
# genetics.binaRies, devtools.
#
# Usage, from the package root:
#   Rscript data-raw/ld_ref.R
#
# The fetched genotypes are cached under tools::R_user_dir("mrpipeline",
# "cache"), so a re-run works offline; delete the cache file to re-download.

base_url <- "https://s3.amazonaws.com/1000genomes/release/20130502"
vcf_url <- file.path(
  base_url,
  "ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
)
panel_url <- file.path(
  base_url,
  "integrated_call_samples_v3.20130502.ALL.panel"
)
out_prefix <- file.path("inst", "extdata", "ld_ref")

stopifnot(file.exists("DESCRIPTION"))
devtools::load_all(quiet = TRUE)

snps <- cd40_exposure[
  order(cd40_exposure$pos.exposure),
  c(
    "SNP",
    "chr.exposure",
    "pos.exposure",
    "effect_allele.exposure",
    "other_allele.exposure"
  )
]
stopifnot(
  !anyDuplicated(snps$SNP),
  all(snps$chr.exposure == "20"),
  setequal(snps$SNP, sjogren_outcome$rsids)
)

# -- Fetch (with retries and an on-disk cache) --

with_retries <- function(f, what, attempts = 4) {
  for (i in seq_len(attempts)) {
    result <- tryCatch(f(), error = function(e) e)
    if (!inherits(result, "error")) {
      return(result)
    }
    message(sprintf(
      "Fetching %s failed (attempt %d/%d): %s",
      what,
      i,
      attempts,
      conditionMessage(result)
    ))
    if (i < attempts) Sys.sleep(5 * i)
  }
  stop("Could not fetch ", what, "; check the connection and re-run.")
}

cache_dir <- tools::R_user_dir("mrpipeline", "cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
cache_file <- file.path(cache_dir, "ld_ref_1kg_eur_cd40.rds")

if (file.exists(cache_file)) {
  message("Using cached genotypes: ", cache_file)
  fetched <- readRDS(cache_file)
} else {
  panel <- with_retries(
    function() utils::read.table(panel_url, header = TRUE),
    "the 1000 Genomes sample panel"
  )
  eur <- panel[panel$super_pop == "EUR", ]
  stopifnot(nrow(eur) == 503)

  # Rsamtools downloads the remote .tbi into the working directory; keep it
  # out of the package tree
  index_dir <- tempfile("ld_ref_tbi")
  dir.create(index_dir)
  vcf <- with_retries(
    function() {
      old <- setwd(index_dir)
      on.exit(setwd(old))
      VariantAnnotation::readVcf(
        Rsamtools::TabixFile(vcf_url),
        "hg19",
        param = VariantAnnotation::ScanVcfParam(
          which = GenomicRanges::GRanges(
            "20",
            IRanges::IRanges(min(snps$pos.exposure), max(snps$pos.exposure))
          ),
          info = NA,
          geno = "GT",
          samples = eur$sample
        )
      )
    },
    "the chr20 CD40 slice of the 1000 Genomes VCF"
  )
  vcf <- vcf[names(vcf) %in% snps$SNP]
  ranges <- SummarizedExperiment::rowRanges(vcf)
  fetched <- list(
    samples = eur[match(colnames(vcf), eur$sample), ],
    variants = data.frame(
      SNP = names(vcf),
      pos = GenomicRanges::start(ranges),
      ref = as.character(ranges$REF),
      alt = vapply(
        ranges$ALT,
        function(a) paste(as.character(a), collapse = ","),
        character(1)
      )
    ),
    gt = VariantAnnotation::geno(vcf)$GT
  )
  saveRDS(fetched, cache_file)
}

# -- Check and recode --

variants <- fetched$variants
gt <- fetched$gt
samples <- fetched$samples
stopifnot(
  identical(samples$sample, colnames(gt)),
  !anyDuplicated(variants$SNP),
  setequal(variants$SNP, snps$SNP)
)
variants <- variants[match(snps$SNP, variants$SNP), ]
gt <- gt[snps$SNP, , drop = FALSE]
stopifnot(all(variants$pos == snps$pos.exposure))

geno <- vapply(
  seq_len(nrow(snps)),
  function(j) {
    alts <- strsplit(variants$alt[j], ",", fixed = TRUE)[[1]]
    wanted <- c(snps$effect_allele.exposure[j], snps$other_allele.exposure[j])
    k <- which(vapply(
      alts,
      function(a) setequal(c(variants$ref[j], a), wanted),
      logical(1)
    ))
    if (length(k) != 1) {
      stop(
        snps$SNP[j],
        ": no single 1000 Genomes ALT matches the exposure alleles"
      )
    }
    codes <- strsplit(gt[j, ], "[|/]")
    stopifnot(all(lengths(codes) == 2))
    codes <- do.call(rbind, codes)
    stopifnot(all(codes %in% as.character(0:length(alts))))
    other <- codes != "0" & codes != as.character(k)
    if (any(other)) {
      message(sprintf(
        "%s: %d allele(s) other than %s/%s coded as %s",
        snps$SNP[j],
        sum(other),
        variants$ref[j],
        alts[k],
        variants$ref[j]
      ))
    }
    allele <- ifelse(codes == as.character(k), alts[k], variants$ref[j])
    paste(allele[, 1], allele[, 2])
  },
  character(ncol(gt))
)

# -- Write PLINK files --

plink <- genetics.binaRies::get_plink_binary()
work <- tempfile("ld_ref")
dir.create(work)
text_prefix <- file.path(work, "ld_ref")

ped <- data.frame(
  fid = samples$sample,
  iid = samples$sample,
  pat = 0,
  mat = 0,
  sex = ifelse(samples$gender == "male", 1, 2),
  pheno = -9,
  geno
)
utils::write.table(
  ped,
  paste0(text_prefix, ".ped"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
utils::write.table(
  data.frame(snps$chr.exposure, snps$SNP, 0, snps$pos.exposure),
  paste0(text_prefix, ".map"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

status <- system2(
  plink,
  c("--file", text_prefix, "--make-bed", "--out", text_prefix),
  stdout = FALSE,
  stderr = FALSE
)
stopifnot(status == 0)
for (ext in c("bed", "bim", "fam")) {
  file.copy(
    paste0(text_prefix, ".", ext),
    paste0(out_prefix, ".", ext),
    overwrite = TRUE
  )
}

# -- Summary --

bim <- utils::read.table(paste0(out_prefix, ".bim"), colClasses = "character")
fam <- utils::read.table(paste0(out_prefix, ".fam"))
stopifnot(
  nrow(fam) == 503,
  !anyNA(geno) && !any(geno == "0 0"),
  identical(bim$V2, snps$SNP),
  # test-run_mr.R's single-SNP test relies on this A1 (ieugwasr reads a lone
  # "T" allele column as logical TRUE)
  bim$V5[bim$V2 == "rs1883832"] == "T"
)

ld <- mrpipeline:::compute_ld_matrix(snps$SNP, out_prefix, plink)$ld
stopifnot(!anyNA(ld))
diag(ld) <- 0
message(sprintf(
  "Wrote %s.{bed,bim,fam}: %d individuals x %d SNPs; max |r| %.3f; %d pairs with |r| > 0.8; .bed %d bytes",
  out_prefix,
  nrow(fam),
  nrow(bim),
  max(abs(ld)),
  sum(abs(ld) > 0.8) / 2,
  file.size(paste0(out_prefix, ".bed"))
))
