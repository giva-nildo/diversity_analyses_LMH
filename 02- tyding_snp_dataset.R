##### 2 - Tidying the geno dataset
### GIvanildo Rodrigues da Silva
# February 18th, 2026

# =========================
# Libraries
# =========================
library(dplyr)
library(tibble)
library(tidyr)
library(readxl)
library(writexl)

# =========================
# Parameters
# =========================
ploidy <- 4
maf_min_marker <- 0.01       # marker MAF threshold
sample_missing_max <- 0.30   # max missing rate allowed per sample (0-1)
locus_missing_max  <- 0.50    # max missing rate allowed per locus/marker (0-1)

# =========================
# 0) Basic info + keep a copy for marker parsing (REF/ALT)
# =========================
stopifnot(!is.null(rownames(snp)))
snp_copy <- snp

n_samples_start <- nrow(snp)
n_markers_start <- ncol(snp)

# =========================
# 1) Rename samples using key file
# =========================
snp_tbl <- snp %>%
  as.data.frame() %>%
  rownames_to_column("sample_code")

sample_key <- read_excel("sample_names_genotyping_process.xlsx") %>%
  as.data.frame() %>%
  rename(sample_code = `RAPiD Genomics Sample Code`) %>%
  select(sample_code, new_name)

new_snp_mat <- snp_tbl %>%
  left_join(sample_key, by = "sample_code") %>%
  mutate(sample_code = coalesce(new_name, sample_code)) %>%
  select(-new_name) %>%
  column_to_rownames("sample_code") %>%
  as.matrix()

stopifnot(!any(duplicated(rownames(new_snp_mat))))
stopifnot(!any(duplicated(colnames(new_snp_mat))))

# =========================
# 2) Clean marker names -> keep CHROM_POS, transpose -> remove duplicated markers
#    chr01_535192_12_C_T -> chr01_535192
# =========================
colnames(new_snp_mat) <- sub("_[^_]+_[^_]+_[^_]+$", "", colnames(new_snp_mat))

G_mk <- t(new_snp_mat)  # markers x samples
dup_mk <- duplicated(rownames(G_mk))
G_nodup <- G_mk[!dup_mk, , drop = FALSE]

# =========================
# 3) Build marker metadata (REF/ALT) from ORIGINAL marker names and join
#    Remove markers with ALT == NA / "" / "."
# =========================
marker_info <- tibble(marker_full = colnames(snp_copy)) %>%
  separate(
    marker_full,
    into = c("CHROM", "POS", "ID", "REF", "ALT"),
    sep = "_",
    remove = FALSE,
    extra = "merge",
    fill  = "right"
  ) %>%
  mutate(
    POS = suppressWarnings(as.integer(POS)),
    REF = na_if(REF, "NA"),
    ALT = na_if(ALT, "NA"),
    ALT = na_if(ALT, "."),
    ALT = na_if(ALT, "")
  ) %>%
  filter(!is.na(CHROM), !is.na(POS)) %>%
  mutate(marker_key = paste0(CHROM, "_", POS)) %>%
  group_by(marker_key) %>%
  slice(1) %>%
  ungroup() %>%
  filter(!is.na(ALT))

snp_duplicated_rm_df <- as.data.frame(G_nodup) %>%
  rownames_to_column("marker_key") %>%
  left_join(marker_info %>% select(marker_key, CHROM, POS, ID, REF, ALT), by = "marker_key") %>%
  filter(!is.na(ALT)) %>%
  relocate(CHROM, POS, ID, REF, ALT, .after = marker_key) %>%
  column_to_rownames("marker_key")

# =========================
# 3.5) Filter loci by missingness (max NA per locus)
# =========================
meta_cols <- intersect(colnames(snp_duplicated_rm_df), c("CHROM","POS","ID","REF","ALT"))
geno_cols <- setdiff(colnames(snp_duplicated_rm_df), meta_cols)

G0 <- as.matrix(snp_duplicated_rm_df[, geno_cols, drop = FALSE])
storage.mode(G0) <- "numeric"

locus_missing_rate <- rowMeans(is.na(G0))
keep_locus_missing <- locus_missing_rate <= locus_missing_max

df_locusmiss <- snp_duplicated_rm_df[keep_locus_missing, , drop = FALSE]

# =========================
# 4) Filter markers by MAF (tetraploid dosage)
# =========================
G <- as.matrix(df_locusmiss[, geno_cols, drop = FALSE])
storage.mode(G) <- "numeric"

p_alt <- rowMeans(G, na.rm = TRUE) / ploidy
maf_marker <- pmin(p_alt, 1 - p_alt)

keep_marker <- !is.na(maf_marker) & maf_marker >= maf_min_marker
df_maf <- df_locusmiss[keep_marker, , drop = FALSE]
maf_kept <- maf_marker[keep_marker]

# =========================
# 5) Filter samples by missingness (NAs per sample)
# =========================
G2 <- as.matrix(df_maf[, geno_cols, drop = FALSE])
storage.mode(G2) <- "numeric"

sample_missing_rate <- colMeans(is.na(G2))
keep_sample <- sample_missing_rate <= sample_missing_max

df_final <- df_maf[, c(meta_cols, geno_cols[keep_sample]), drop = FALSE]

# =========================
# 6) QC report
# =========================
qc_report <- list(
  n_samples_start = n_samples_start,
  n_markers_start = n_markers_start,
  n_markers_after_chrpos_and_nodup = nrow(G_nodup),
  n_markers_after_refalt_nonNA = nrow(snp_duplicated_rm_df),
  n_markers_after_locus_missing = nrow(df_locusmiss),
  n_markers_after_maf = nrow(df_maf),
  n_samples_after_missing = sum(keep_sample),
  maf_summary = summary(maf_kept),
  locus_missing_summary = summary(locus_missing_rate),
  sample_missing_summary = summary(sample_missing_rate)
)

print(qc_report)

# =========================
# 7) Save outputs
# =========================
saveRDS(df_final, "final_snp_qc_maf_missing_refalt.rds")
write_xlsx(df_final, "final_snp_qc_maf_missing_refalt.xlsx")








hardy_Fz_eq8 <- function(G, ploidy = 4) {
  G <- as.matrix(G)
  storage.mode(G) <- "numeric"
  k <- ploidy
  
  # per-locus sample size (non-missing individuals)
  Nl <- colSums(!is.na(G))
  
  # hi matrix (same dims as G); NA stays NA
  hi <- 2 * G * (k - G) / (k * (k - 1))
  
  # sums over individuals (per locus)
  sum_hi <- colSums(hi, na.rm = TRUE)
  
  # allele frequencies from dosages (biallelic)
  p_alt <- colSums(G, na.rm = TRUE) / (k * Nl)
  p_ref <- 1 - p_alt
  sum_p2 <- p_alt^2 + p_ref^2
  
  # eqn (8) in “second form”: Fz = 1 - (Σ_i hi) / [ N^2/(N-1) * (1 - Σ_a p_a^2 - (k-1)/(kN^2) Σ_i hi ) ]
  denom <- (Nl^2 / (Nl - 1)) * (1 - sum_p2 - ((k - 1) / (k * Nl^2)) * sum_hi)
  
  Fz_locus <- 1 - (sum_hi / denom)
  
  # guard against loci with too few samples or zero/negative denom
  Fz_locus[Nl < 2 | !is.finite(Fz_locus) | denom <= 0] <- NA_real_
  
  # multilocus estimate by summing numerator and denominator across loci (Hardy’s multilocus idea; eqn 9)
  keep <- which(!is.na(Fz_locus) & denom > 0)
  Fz_multi <- 1 - sum(sum_hi[keep]) / sum(denom[keep])
  
  list(
    Fz_multilocus = Fz_multi,
    per_locus = data.frame(
      locus = colnames(G),
      N = Nl,
      p_alt = p_alt,
      sum_p2 = sum_p2,
      sum_hi = sum_hi,
      denom = denom,
      Fz = Fz_locus,
      stringsAsFactors = FALSE
    )
  )
}

# Example:
test <- df_final[,-c(1:5)]

res <- hardy_Fz_eq8(test, ploidy = 4)
res$Fz_multilocus
head(res$per_locus[order(res$per_locus$Fz, decreasing = TRUE), ], 10)
write_xlsx(res$per_locus, "testing.xlsx")
