library(dplyr)
library(readxl)
library(GenomicRanges)
library(rtracklayer)
library(tibble)

# =========================
# USER SETTINGS
# =========================
window_size <- 3730000 #considered as the upstream and downstream regions

setwd("C:/Users/givan/OneDrive/Área de Trabalho/Repositorios/diversity_LMH")
input_xlsx <- "C:/Users/givan/OneDrive/Área de Trabalho/Repositorios/diversity_LMH/circos_plot_data/all_points_data.xlsx"
input_sheet <- "all_for_annotation"

gff_file <- "C:/Users/givan/OneDrive/Área de Trabalho/Repositorios/diversity_LMH/reference_genome/DM_1-3_516_R44_potato.v6.1.hc_gene_models.gff3.gz"

output_csv_overlap <- "circos_plot_data/annotated_regions_overlap.csv"
output_csv_nearest <- "circos_plot_data/annotated_regions_nearest.csv"

# =========================
# 1) LOAD SIGNIFICANT REGIONS
# =========================
sig_regions <- read_excel(input_xlsx, sheet = input_sheet) %>%
  as.data.frame()

# Check available columns
print(colnames(sig_regions))

# -------------------------
# Standardize coordinates
# -------------------------
# This block accepts either:
#   CHROM + POS
# or
#   CHROM + start + end
# and creates start/end consistently.

if (!"CHROM" %in% names(sig_regions)) {
  stop("Column 'CHROM' was not found in the annotation sheet.")
}

if ("POS" %in% names(sig_regions)) {
  sig_regions <- sig_regions %>%
    mutate(
      CHROM = as.character(CHROM),
      start = as.integer(POS),
      end   = as.integer(POS)
    )
} else if (all(c("start", "end") %in% names(sig_regions))) {
  sig_regions <- sig_regions %>%
    mutate(
      CHROM = as.character(CHROM),
      start = as.integer(start),
      end   = as.integer(end)
    )
} else {
  stop("Your file must contain either 'POS' or both 'start' and 'end'.")
}

# Create marker/region ID if absent
if (!"marker_ID" %in% names(sig_regions)) {
  sig_regions <- sig_regions %>%
    mutate(marker_ID = paste0(CHROM, "_", start))
}

# Add region ID for traceability
sig_regions <- sig_regions %>%
  mutate(region_ID = paste0(marker_ID, "_", row_number()))

# Expand region upstream/downstream
sig_regions <- sig_regions %>%
  mutate(
    query_start = pmax(start - window_size, 1L),
    query_end   = end + window_size
  )

# Quick check
print(head(sig_regions[, c("region_ID", "CHROM", "start", "end", "query_start", "query_end")]))

# =========================
# 2) LOAD GFF
# =========================
gff <- import(gff_file)

cat("\nGFF loaded.\n")
print(gff)
cat("\nFeature types in GFF:\n")
print(table(gff$type))

# -------------------------
# Choose annotation feature
# -------------------------
# Prefer "gene" if available; otherwise fall back to "mRNA"
if ("gene" %in% unique(gff$type)) {
  features <- gff[gff$type == "gene"]
  feature_type_used <- "gene"
} else if ("mRNA" %in% unique(gff$type)) {
  features <- gff[gff$type == "mRNA"]
  feature_type_used <- "mRNA"
} else {
  stop("Neither 'gene' nor 'mRNA' features were found in the GFF.")
}

cat("\nUsing feature type for annotation:", feature_type_used, "\n")

if (length(features) == 0) {
  stop("No annotation features found in the GFF.")
}

# =========================
# 3) HARMONIZE CHROMOSOME NAMES
# =========================
cat("\nFirst seqlevels in sig_regions:\n")
print(unique(sig_regions$CHROM)[1:min(10, length(unique(sig_regions$CHROM)))])

cat("\nFirst seqlevels in GFF:\n")
print(head(seqlevels(features), 10))

# Simple harmonization helper
normalize_chr <- function(x) {
  x <- as.character(x)
  x <- sub("^0+", "", sub("^chr", "", x, ignore.case = TRUE))
  paste0("chr", sprintf("%02d", as.integer(x)))
}

# Try to standardize both sides if possible
sig_regions$CHROM_std <- normalize_chr(sig_regions$CHROM)

feature_seq <- as.character(seqnames(features))
feature_seq_std <- suppressWarnings(normalize_chr(feature_seq))

# If normalization produced usable chromosome names, replace them
if (sum(!is.na(feature_seq_std)) > 0) {
  seqlevels(features) <- unique(feature_seq)
  seqnames(features) <- feature_seq_std
}

sig_regions$CHROM_use <- sig_regions$CHROM_std

cat("\nStandardized seqlevels in sig_regions:\n")
print(unique(sig_regions$CHROM_use)[1:min(10, length(unique(sig_regions$CHROM_use)))])

cat("\nStandardized seqlevels in GFF features:\n")
print(unique(as.character(seqnames(features)))[1:min(10, length(unique(as.character(seqnames(features)))))])

# Keep only chromosomes present in both
common_chr <- intersect(unique(sig_regions$CHROM_use), unique(as.character(seqnames(features))))
if (length(common_chr) == 0) {
  stop("No matching chromosome names between sig_regions and GFF after harmonization.")
}

sig_regions <- sig_regions %>% filter(CHROM_use %in% common_chr)
features <- features[as.character(seqnames(features)) %in% common_chr]

# =========================
# 4) CONVERT TO GRanges
# =========================
sig_gr <- GRanges(
  seqnames = sig_regions$CHROM_use,
  ranges = IRanges(start = sig_regions$query_start, end = sig_regions$query_end)
)

mcols(sig_gr) <- sig_regions

if (length(sig_gr) == 0) {
  stop("No GRanges objects were created from sig_regions.")
}

# =========================
# 5) FIND OVERLAPS
# =========================
overlaps <- findOverlaps(sig_gr, features)

cat("\nNumber of overlaps found:", length(overlaps), "\n")

# =========================
# 6) EXTRACT OVERLAP ANNOTATIONS
# =========================
feature_df <- as.data.frame(features)
query_df   <- as.data.frame(sig_gr)

# Safe extraction of gene/transcript identifiers
pick_col <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) return(rep(NA_character_, nrow(df)))
  as.character(df[[hit[1]]])
}

feature_df$feature_ID <- pick_col(feature_df, c("ID", "gene_id", "transcript_id", "Name"))
feature_df$feature_Name <- pick_col(feature_df, c("Name", "gene_name", "gene", "ID"))
feature_df$feature_biotype <- pick_col(feature_df, c("biotype", "gene_biotype"))

if (length(overlaps) > 0) {
  annotated_overlap <- cbind(
    query_df[queryHits(overlaps), , drop = FALSE],
    feature_df[subjectHits(overlaps),
               c("seqnames", "start", "end", "width", "strand", "type",
                 "feature_ID", "feature_Name", "feature_biotype"),
               drop = FALSE]
  )
  
  colnames(annotated_overlap)[colnames(annotated_overlap) == "seqnames"] <- "feature_CHROM"
  colnames(annotated_overlap)[colnames(annotated_overlap) == "start"] <- "feature_start"
  colnames(annotated_overlap)[colnames(annotated_overlap) == "end"] <- "feature_end"
  colnames(annotated_overlap)[colnames(annotated_overlap) == "width"] <- "feature_width"
  colnames(annotated_overlap)[colnames(annotated_overlap) == "strand"] <- "feature_strand"
  colnames(annotated_overlap)[colnames(annotated_overlap) == "type"] <- "feature_type"
  
  write.csv(annotated_overlap, output_csv_overlap, row.names = FALSE)
  cat("\nOverlap annotation saved to:\n", output_csv_overlap, "\n")
} else {
  annotated_overlap <- data.frame()
  cat("\nNo direct overlaps found.\n")
}

# =========================
# 7) FIND NEAREST FEATURE FOR EACH REGION
# =========================
nearest_idx <- nearest(sig_gr, features, ignore.strand = TRUE)

annotated_nearest <- cbind(
  query_df,
  feature_df[nearest_idx,
             c("seqnames", "start", "end", "width", "strand", "type",
               "feature_ID", "feature_Name", "feature_biotype"),
             drop = FALSE]
)

colnames(annotated_nearest)[colnames(annotated_nearest) == "seqnames"] <- "nearest_CHROM"
colnames(annotated_nearest)[colnames(annotated_nearest) == "start"] <- "nearest_start"
colnames(annotated_nearest)[colnames(annotated_nearest) == "end"] <- "nearest_end"
colnames(annotated_nearest)[colnames(annotated_nearest) == "width"] <- "nearest_width"
colnames(annotated_nearest)[colnames(annotated_nearest) == "strand"] <- "nearest_strand"
colnames(annotated_nearest)[colnames(annotated_nearest) == "type"] <- "nearest_type"

# Distance from query interval to nearest feature
dist_to_nearest <- distance(sig_gr, features[nearest_idx], ignore.strand = TRUE)
annotated_nearest$distance_to_nearest <- dist_to_nearest

write.csv(annotated_nearest, output_csv_nearest, row.names = FALSE)
cat("\nNearest-feature annotation saved to:\n", output_csv_nearest, "\n")
