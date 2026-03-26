####################################
# pcadapt scan on tetraploid dosage
####################################
# Givanildo Rodrigues da Silva
# 03/26/2026

install.packages("pcadapt")
suppressPackageStartupMessages({
  library(data.table)
  library(pcadapt)
})

# -----------------------
# USER SETTINGS
# -----------------------
ploidy <- 4
maf_min <- 0.01          # dosage-aware MAF filter
missing_max <- 0.30      # max SNP missing rate
K <- 4                # number of PCs used by pcadapt (choose via scree plot)
alpha_snp <- 0.05        # SNP-level BH-FDR threshold for outliers

# window settings (regions)
window_bp <- 1e6         # 1 Mb
overlap <- 0.5           # 50% overlap -> step = 0.5 Mb
alpha_win <- 0.05        # window-level BH-FDR threshold

# -----------------------
# 0) CHECK INPUTS
# -----------------------
map <- read_table("C:/Users/givan/OneDrive/Área de Trabalho/Repositorios/diversity_LMH/selective_sweep_analyses/poolfstat/map_file.tsv")
stopifnot(exists("G_mx"), exists("map"))
G_mx <- as.matrix(G_mx)
storage.mode(G_mx) <- "numeric"

map <- as.data.table(map)
stopifnot(all(c("marker_ID","CHROM","POS") %in% names(map)))

# Align map to G_mx row order
if (is.null(rownames(G_mx))) stop("G_mx must have rownames = marker_ID.")
map <- unique(map[, .(marker_ID, CHROM, POS)])
map[, `:=`(marker_ID = as.character(marker_ID),
           CHROM = as.character(CHROM),
           POS = as.integer(POS))]

m <- match(rownames(G_mx), map$marker_ID)
if (anyNA(m)) {
  cat("Some G_mx markers not found in map. Examples:\n")
  print(head(rownames(G_mx)[is.na(m)], 10))
  stop("Fix map vs G_mx marker_ID mismatch.")
}
map <- map[m]
stopifnot(identical(map$marker_ID, rownames(G_mx)))

# -----------------------
# 1) QC FILTERS (dosage-aware)
# -----------------------
# Missingness per SNP
miss_rate <- rowMeans(is.na(G_mx))
keep_miss <- miss_rate <= missing_max



keep <- keep_miss & keep_maf

G <- G_mx[drop = FALSE]
map2 <- map

cat("SNPs retained after QC:", nrow(G), "\n")

cat("SNPs retained after removing zero variance:", nrow(G), "\n")

# -----------------------
# 2) Run pcadapt
# -----------------------
# pcadapt expects: SNPs in rows, individuals in columns (this matches G_scaled)
library(pcadapt)

# G_mx: SNP x IND dosage matrix (0..4), NA allowed
G_in <- as.matrix(G)
storage.mode(G_in) <- "integer"

# Convert to pcadapt input object (type="pcadapt" = individuals in columns, SNPs in lines)
X <- read.pcadapt(G_in, type = "pcadapt")

# Run pcadapt (set ploidy=4 for tetraploid dosage)
K <- 2
pc <- pcadapt(X, K = K, method = "mahalanobis", ploidy = 4, min.maf = 0.01)

# choose K by inspecting scree plot 
#plot(pc, option = "screeplot")

pop_name <- as.list(y)
plot(pc, option = "scores", pop = y)

# Results
p <- pc$pvalues
q <- p.adjust(p, "BH")


res <- data.table(
  marker_ID = rownames(G_in),
  CHROM = map2$CHROM,
  POS = map2$POS,
  p = p,
  q = q
)
setorder(res, CHROM, POS)

# Identify outlier SNPs
hits <- res[q <= 0.01]

cat("Outlier SNPs (BH q<= ", alpha_snp, "): ", nrow(hits), "\n", sep="")
write_csv(hits, "pcadapt_significative.csv")



# Export SNP results
# fwrite(res, "pcadapt_snp_results.tsv", sep="\t")

# -----------------------
# 3) Manhattan-like plot (base R)
# -----------------------
# Build cumulative positions for plotting
chr_levels <- unique(res$CHROM)
# order chromosomes naturally if like chr01..chr12
chr_levels <- chr_levels[order(chr_levels)]
res[, CHROM := factor(CHROM, levels = chr_levels)]
setorder(res, CHROM, POS)

chr_max <- res[, .(maxPOS = max(POS, na.rm = TRUE)), by = CHROM]
chr_max[, offset := c(0, cumsum(head(maxPOS, -1)))]
res <- merge(res, chr_max[, .(CHROM, offset)], by = "CHROM")
res[, cumPOS := POS + offset]

res[, neglogq := -log10(pmax(q, .Machine$double.xmin))]

plot(res$cumPOS/1e6, res$neglogq,
     pch = 16, cex = 0.5,
     col = as.integer(res$CHROM),
     xlab = "Genome position (Mb, cumulative)",
     ylab = expression(-log[10](q)))

abline(h = -log10(0.01), lty = 2)

plot(pc, option = "manhattan")
plot(pc, option = "qqplot")
plot(pc, option = "stat.distribution")

qval <- qvalue(pc$pvalues)$qvalues
alpha <- 0.05
outliers <- which(qval < alpha)
length(outliers)

# BH - just verif
padj <- p.adjust(pc$pvalues,method="BH")
alpha <- 0.05
outliers <- which(padj < alpha)
length(outliers)

#Bonferroni
padj <- p.adjust(pc$pvalues,method="bonferroni")
alpha <- 0.05
outliers <- which(padj < alpha)
length(outliers)


