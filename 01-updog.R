####################
##### Using updog to convert to polyploid format
####################
# Givanildo Rodrigues da Silva
#February, 18th, 2026

library(updog)
library(vcfR)

vcf <- read.vcfR("USP_191101_RAW_SNPs.vcf")
vcf
# Variant-level fields (CHROM POS ID REF ALT ...)
fix <- getFIX(vcf)

### removing no-bialelic SNPs
bial <- !grepl(",", fix[, "ALT"])
vcf <- vcf[bial, ]
vcf


alleles <- data.frame(
  CHROM = fix[, "CHROM"],
  POS   = as.integer(fix[, "POS"]),
  ID    = fix[, "ID"],
  REF   = fix[, "REF"],
  ALT   = fix[, "ALT"],
  stringsAsFactors = FALSE
)


ROmat <- extract.gt(vcf, element = "RO", as.numeric = TRUE)
dim(ROmat)
ROmat[1:10,1:5]

DPmat <- extract.gt(vcf, element = "DP", as.numeric = TRUE)
dim(DPmat)
DPmat[1:10,1:5]

#rm(vcf)
gc()
save(ROmat, DPmat, file = "mat.rda")

sum(rowMeans(DPmat, na.rm = TRUE) > 40) # number of SNPs average read depth > 40
ROmat <- ROmat[which(rowMeans(DPmat, na.rm = TRUE) > 40),]; dim(ROmat)
DPmat <- DPmat[which(rowMeans(DPmat, na.rm = TRUE) > 40),]; dim(DPmat)

save(ROmat, DPmat, file = "mat_mean_filtered.rda")

parallel::detectCores()
parallelly::availableCores()

mout <- multidog(refmat = ROmat,
                 sizemat = DPmat,
                 ploidy = 4,
                 model = "norm",
                 nc = 15)
save(mout, file="mout.rda")
rm(ROmat, DPmat)
gc()

snp <- t(format_multidog(mout, varname = "geno")); dim(snp)

key_snp_base <- sub("_[^_]+$", "", colnames(snp))   # chr01_434399

key_vcf_base <- paste0(alleles$CHROM, "_", alleles$POS)

idx <- match(key_snp_base, key_vcf_base)

colnames(snp) <- make.unique(
  paste0(colnames(snp), "_", alleles$REF[idx], "_", alleles$ALT[idx])
)

snp[1:10,1:10]
save(snp, file = "snp.rda")

write.csv(snp, "snp_LMH.csv")
