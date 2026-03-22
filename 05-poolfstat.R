###############################################################################
# poolfstat genome scan (tetraploid dosage 0–4) — two populations
# Author: Givanildo Rodrigues da Silva
# Date:   2025-03-01
#
# Goal
# ----
# Build a poolfstat `countdata` object from a dosage matrix (0–4) and run:
#   (i) SNP-wise FST and genome-wide FST
#  (ii) sliding-window multi-locus FST scan (windows defined by SNP count)
#
# Key assumptions
# ---------------
# 1) Input genotype table contains biallelic SNPs with allele dosage values:
#    0..4 = number of ALT alleles per individual at each SNP (autotetraploid).
# 2) Exactly two populations are defined in `meta` (pop labels).
# 3) Marker IDs follow the pattern "chr##_POS" (e.g., "chr01_432462") so we can
#    align REF/ALT and map information robustly.
#
# Notes
# -----
# - poolfstat uses allele counts. We convert dosage to allele counts per pop:
#       ALT_count = sum(dosage) across individuals in the population
#       N_alleles = ploidy * (# non-missing individuals)
#       REF_count = N_alleles - ALT_count
# - Sliding windows in computeFST(sliding.window.size=...) are defined by
#   *number of SNPs*, not bp. Physical sizes vary with SNP density.
###############################################################################

library(poolfstat)
library(data.table)
library(dplyr)
library(readxl)
library(tibble)


# =========================
# USER SETTINGS
# =========================

# Working directory
wd <- "C:/Users/givan/OneDrive/Área de Trabalho/Repositorios/diversity_LMH/"
setwd(wd)

# Input files
final_xlsx   <- "final_snp_qc_maf_missing_refalt.xlsx"
origin_xlsx  <- "origin_genotypes.xlsx"

# Output directory
outdir <- file.path(
  wd,
  "selective_sweep_analyses/poolfstat/poolfstat_out"
)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Ploidy
ploidy <- 4L

# Sliding-window scan (number of SNPs per window)
window_snps <- 50L

# =========================
# 1) LOAD DATA
# =========================

df_final <- read_excel(final_xlsx)

# -------------------------
# Identify columns
# -------------------------
# Expect the genotype table to include these metadata columns.
meta_cols <- intersect(colnames(df_final), c("marker_ID", "CHROM", "POS", "ID", "REF", "ALT"))

# If marker_ID is not in meta_cols but exists as "ID" or similar, adapt here.
# This script assumes there is a "marker_ID" column.
stopifnot("marker_ID" %in% colnames(df_final))
stopifnot(all(c("CHROM", "POS", "REF", "ALT") %in% colnames(df_final)))

geno_cols <- setdiff(colnames(df_final), meta_cols)

# -------------------------
# Dosage table (marker_ID + individuals)
# -------------------------
dos <- as.data.table(df_final[, c("marker_ID", geno_cols)])
# Ensure numeric dosage
for (cc in geno_cols) set(dos, j = cc, value = as.integer(dos[[cc]]))

# -------------------------
# Map table (marker_ID, CHROM, POS)
# -------------------------
map <- as.data.table(df_final[, c("marker_ID", "CHROM", "POS")])
map[, CHROM := as.character(CHROM)]
map[, POS   := as.integer(POS)]
map <- unique(map)

# Sort map by physical position (needed for meaningful windows across genome)
setorder(map, CHROM, POS)

# -------------------------
# Meta table (sample_ID, pop)
# -------------------------
# Build sample list from genotype columns.
genos <- data.table(sample_ID = geno_cols)

# Read origin info and join
origins <- as.data.table(read_xlsx(origin_xlsx))
setnames(origins, old = "Clone Name", new = "sample_ID")
setnames(origins, old = "Cluster by origin", new = "pop")

meta <- merge(genos, origins[, .(sample_ID, pop)], by = "sample_ID", all.x = TRUE)

# Remove samples without pop assignment
meta <- meta[!is.na(pop)]
meta[, pop := as.character(pop)]
meta <- unique(meta)

# Check exactly two populations
pops <- unique(meta$pop)
if (length(pops) != 2) {
  stop("This script expects exactly TWO populations in meta$pop. Found: ",
       paste(pops, collapse = ", "))
}
setorder(meta, pop, sample_ID)

# =========================
# 2) ALIGN SAMPLES AND MARKERS (CRITICAL)
# =========================

# Keep only sample columns present in `dos`
sample_cols <- setdiff(names(dos), "marker_ID")
meta <- meta[sample_ID %in% sample_cols]
if (nrow(meta) == 0) stop("No sample_ID in meta matches dosage columns.")

# Reorder dosage columns to match meta order
dos <- dos[, c("marker_ID", meta$sample_ID), with = FALSE]

# Join map onto dos without losing dosage order; keep only markers present in map.
# Using map[dos, on=...] preserves the order of `dos` (i.e., marker order in dosage table).
dt0 <- map[dos, on = "marker_ID", nomatch = 0]

# Track original order in dosage table
dt0[, orig_idx := .I]

# Now reorder to physical order (recommended for genome scan visualization)
dt <- copy(dt0)
setorder(dt, CHROM, POS)
dt[, phys_idx := .I]

# Build dosage matrix X in the *same* order as dt (SNPs x individuals)
X <- as.matrix(dt[, meta$sample_ID, with = FALSE])
storage.mode(X) <- "integer"

# Sanity checks
stopifnot(nrow(X) == nrow(dt))
stopifnot(ncol(X) == nrow(meta))
stopifnot(all(X[!is.na(X)] >= 0 & X[!is.na(X)] <= ploidy))

# =========================
# 3) BUILD ALLELE COUNT MATRICES (ref_mat / tot_mat)
# =========================

idxA <- which(meta$pop == pops[1])
idxB <- which(meta$pop == pops[2])

# Number of called individuals per SNP per pop
calledA <- rowSums(!is.na(X[, idxA, drop = FALSE]))
calledB <- rowSums(!is.na(X[, idxB, drop = FALSE]))

# ALT allele counts
altA <- rowSums(X[, idxA, drop = FALSE], na.rm = TRUE)
altB <- rowSums(X[, idxB, drop = FALSE], na.rm = TRUE)

# Total allele copies sampled
nA <- ploidy * calledA
nB <- ploidy * calledB

# REF allele counts
refA <- nA - altA
refB <- nB - altB

ref_mat <- cbind(as.integer(refA), as.integer(refB))
tot_mat <- cbind(as.integer(nA),   as.integer(nB))
colnames(ref_mat) <- colnames(tot_mat) <- pops

# Basic validity checks
stopifnot(all(ref_mat >= 0, na.rm = TRUE))
stopifnot(all(tot_mat >= 0, na.rm = TRUE))

# =========================
# 4) BUILD snp.info (aligned to dt order)
# =========================
# Your df_final already has CHROM, POS, REF, ALT, but its order may differ from dt.
# Align by key = "CHROM_POS", assuming marker_ID follows "CHROM_POS".

snp_info2 <- as.data.table(df_final[, c("CHROM", "POS", "REF", "ALT")])
snp_info2[, CHROM := as.character(CHROM)]
snp_info2[, POS   := as.integer(POS)]
snp_info2[, REF   := as.character(REF)]
snp_info2[, ALT   := as.character(ALT)]
snp_info2[, key := paste0(CHROM, "_", POS)]

dt[, key := as.character(marker_ID)]

m <- match(dt$key, snp_info2$key)
if (anyNA(m)) {
  message("Some markers could not be matched to REF/ALT using CHROM_POS key.")
  print(dt[is.na(m), .(marker_ID, CHROM, POS)][1:10])
  stop("Fix marker_ID / CHROM / POS consistency before proceeding.")
}

snp_aligned <- snp_info2[m]

snp_info_pf <- data.frame(
  Contig   = snp_aligned$CHROM,
  Position = as.integer(snp_aligned$POS),
  Ref      = snp_aligned$REF,
  Alt      = snp_aligned$ALT,
  stringsAsFactors = FALSE
)

# Alignment checks (must pass)
stopifnot(nrow(snp_info_pf) == nrow(dt))
stopifnot(identical(as.character(snp_info_pf$Contig), dt$CHROM))
stopifnot(identical(as.integer(snp_info_pf$Position), as.integer(dt$POS)))

# =========================
# 5) CREATE countdata OBJECT
# =========================
cd <- new("countdata",
          npops = as.integer(2),
          nsnp  = as.integer(nrow(ref_mat)),
          refallele.count = ref_mat,
          total.count     = tot_mat,
          snp.info        = snp_info_pf,
          popnames        = as.character(pops))

# Final validation
stopifnot(cd@nsnp == nrow(cd@refallele.count),
          cd@nsnp == nrow(cd@total.count),
          cd@nsnp == nrow(cd@snp.info))

# =========================
# 6) RUN FST (SNP-wise and window scan)
# =========================

# Genome-wide + SNP-specific FST (no jackknife)
fst_global <- computeFST(cd, method = "Anova", nsnp.per.bjack.block = 0, verbose = TRUE)
saveRDS(fst_global, file = file.path(outdir, "fst_global.rds"))
fwrite(as.data.table(fst_global$snp.Fstats), file.path(outdir, "fst_snp.tsv"), sep = "\t")

# Sliding-window scan (windows defined by SNP count)
scan <- computeFST(cd, method = "Anova", sliding.window.size = window_snps,
                   nsnp.per.bjack.block = 0, verbose = TRUE)
saveRDS(scan, file = file.path(outdir, paste0("fst_scan_window", window_snps, ".rds")))

# Export window table if available
if (!is.null(scan$sliding.windows.fvalues)) {
  fwrite(as.data.table(scan$sliding.windows.fvalues),
         file.path(outdir, paste0("fst_windows_nsnp_", window_snps, ".tsv")),
         sep = "\t")
}

# =========================
# 7) PLOT WINDOW SCAN
# =========================

w <- scan$sliding.windows.fvalues
chr_col <- as.integer(factor(w$Chr))

png(file.path(outdir, paste0("fst_windows_nsnp_", window_snps, ".png")),
    width = 1100, height = 550, res = 150)

plot(w$CumMidPos / 1e6, w$MultiLocusFst,
     xlab = "Cumulated position (Mb)",
     ylab = "Multi-locus FST",
     pch  = 16, col = chr_col, cex = 0.6)

abline(h = fst_global$Fst, lty = 2)

dev.off()

message("Done. Outputs written to: ", outdir)

## ====================
### Extracting the SNPs
## ====================

library(data.table)

# =========================
# 0) Preconditions
# =========================
stopifnot(exists("cd"))
stopifnot(exists("w") || exists("scan"))

# If you did not create `w` explicitly earlier:
if (!exists("w")) w <- scan$sliding.windows.fvalues

# Convert window table to data.table
w_dt <- as.data.table(w)

# poolfstat window outputs sometimes use column names:
#   Chr, Start, End, MidPos, CumMidPos, nsnp, MultiLocusFst
# Confirm required columns are present:
need_cols <- c("Chr", "MultiLocusFst")
stopifnot(all(need_cols %in% names(w_dt)))

# If Start/End don't exist (rare), you cannot do interval joins.
# (computeFST sliding windows should have Start/End in bp.)
if (!all(c("Start", "End") %in% names(w_dt))) {
  stop("Window table w_dt must contain Start and End columns to map SNPs to windows.")
}

# =========================
# 1) Build SNP table from cd@snp.info
# =========================
# cd@snp.info MUST be aligned with allele count matrices.
# We create a SNP table with Chr/Pos/Ref/Alt + marker_ID.

snp_dt <- as.data.table(cd@snp.info)

# Normalize column names (poolfstat uses Contig/Position; we rename to Chr/Pos)
if (all(c("Contig","Position") %in% names(snp_dt))) {
  setnames(snp_dt, c("Contig","Position"), c("Chr","Pos"))
} else if (all(c("Chr","Pos") %in% names(snp_dt))) {
  # already ok
} else {
  stop("cd@snp.info must contain Contig/Position (or Chr/Pos).")
}

# Ensure correct types
snp_dt[, Chr := as.character(Chr)]
snp_dt[, Pos := as.integer(Pos)]

# If Ref/Alt columns exist, keep them; otherwise create placeholders
if (!("Ref" %in% names(snp_dt))) snp_dt[, Ref := NA_character_]
if (!("Alt" %in% names(snp_dt))) snp_dt[, Alt := NA_character_]

# Create marker_ID (must match your convention: "chr##_POS")
snp_dt[, marker_ID := paste0(Chr, "_", Pos)]

# Keep only useful columns (and remove missing positions)
snp_dt <- snp_dt[!is.na(Pos), .(marker_ID, Chr, Pos, Ref, Alt)]

# =========================
# 2) Define "significant" windows (empirical threshold)
# =========================
# Example: top 5% windows (q=0.95).
q <- 0.95
thr <- as.numeric(quantile(w_dt$MultiLocusFst, probs = q, na.rm = TRUE))

sig_win <- w_dt[MultiLocusFst >= thr]
setorder(sig_win, Chr, Start, End)

# Add a region_id (one id per window, before merging)
sig_win[, region_id := .I]

cat("Threshold (quantile):", q, "\n")
cat("MultiLocusFst threshold:", thr, "\n")
cat("Significant windows:", nrow(sig_win), "\n")

plot(w$CumMidPos / 1e6, w$MultiLocusFst,
     xlab = "Cumulated position (Mb)",
     ylab = "Multi-locus FST",
     pch  = 16, col = chr_col, cex = 0.6)

abline(h = thr, lty = 2)



# =========================
# 3) Interval join: SNP within [Start, End] on same chromosome
# =========================
# Set keys to enable fast non-equi join.
setkey(sig_win, Chr, Start, End)
setkey(snp_dt,  Chr, Pos)

snps_by_region <- snp_dt[
  sig_win,
  on = .(Chr, Pos >= Start, Pos <= End),
  nomatch = 0L,
  allow.cartesian = TRUE
]

# Add window-level stats to each SNP hit (i.* columns come from sig_win)
snps_by_region <- snps_by_region[, .(
  region_id = region_id,
  Chr       = Chr,
  MultiLocusFst = MultiLocusFst,
  marker_ID = marker_ID,
  Pos       = Pos,
  Ref       = Ref,
  Alt       = Alt
)]

cat("SNP×window hits:", nrow(snps_by_region), "\n")

# =========================
# 4) OPTIONAL: Merge adjacent significant windows into broader regions
# =========================
# This is often nicer for reporting (Chr:Start–End).
# Here we merge windows that overlap or are contiguous (gap <= 1 bp).
sig_win2 <- copy(sig_win)
setorder(sig_win2, Chr, Start, End)

sig_win2[, merged_id := {
  # new region whenever:
  #  (a) chromosome changes OR
  #  (b) next window starts after previous ends + 1
  new_region <- c(TRUE, Chr[-1] != Chr[-.N] | Start[-1] > (End[-.N] + 1L))
  cumsum(new_region)
}]

regions_merged <- sig_win2[, .(
  Chr = first(Chr),
  Start = min(Start),
  End   = max(End),
  n_windows = .N,
  maxFst = max(MultiLocusFst, na.rm = TRUE),
  meanFst = mean(MultiLocusFst, na.rm = TRUE)
), by = merged_id][order(-maxFst)]


### for plotting in the circos
write_csv(regions_merged, "stack_line_regions_merged_poolfstat_95percent.csv")

# Map SNPs to merged regions (optional)
# (Join SNPs to merged regions instead of per-window)
regions_for_join <- regions_merged[, .(merged_id, Chr, Start, End)]
setkey(regions_for_join, Chr, Start, End)



# =========================
# 5) Export (optional)
# =========================
# fwrite(snps_by_region, file.path(outdir, "snps_in_sig_windows.tsv"), sep="\t")
# fwrite(regions_merged, file.path(outdir, "sig_regions_merged.tsv"), sep="\t")
# fwrite(snps_by_merged_region, file.path(outdir, "snps_in_merged_regions.tsv"), sep="\t")

# Final objects:
# - snps_by_region
# - regions_merged
# - snps_by_merged_region
###############################################################################





library(data.table)

# 1) SNP table vinda do cd@snp.info
snp_dt <- as.data.table(cd@snp.info)
setnames(snp_dt, c("Contig","Position"), c("Chr","Pos"))
snp_dt[, Pos := as.integer(Pos)]

# 2) Recriar marker_ID
snp_dt[, marker_ID := paste0(Chr, "_", Pos)]

# 3) Preparar para overlap join (SNP como intervalo 1-bp)
snp_dt2 <- snp_dt[, .(marker_ID, Chr, Pos, Ref, Alt)]
snp_dt2[, `:=`(Start = Pos, End = Pos)]

# 4) Janelas significativas
w_dt <- as.data.table(w)
q <- 0.95
thr <- quantile(w_dt$MultiLocusFst, probs = q, na.rm = TRUE)

sig_win <- w_dt[MultiLocusFst >= thr][order(Chr, Start)]

library(data.table)

sig_win <- as.data.table(sig_win)
snp_dt  <- as.data.table(snp_dt)

# 1) enumerar regiões (janelas significativas)
sig_win[, region_id := .I]

# 2) (rápido) join por intervalo: SNP dentro de [Start, End] no mesmo Chr
setkey(sig_win, Chr, Start, End)
setkey(snp_dt,  Chr, Pos)

snps_by_region <- snp_dt[sig_win, on = .(Chr, Pos >= Start, Pos <= End),
                         nomatch = 0L, allow.cartesian = TRUE]

# 3) selecionar colunas úteis
snps_by_region <- snps_by_region[, .(
  region_id,
  Chr,
  MultiLocusFst,
  marker_ID,
  Pos,
  Ref,
  Alt
)]

# pronto
snps_by_region
