############################################################
# Genetic diversity parameters for biallelic polyploid SNPs
# Givanildo Rodrigues da Silva
#February, 25th, 2026
############################################################

# =========================
# Libraries
# =========================
library(dplyr)
library(tibble)
library(tidyr)
library(readxl)
library(writexl)

# =========================
# INPUTS
# =========================

# final dataset:
df_final <- read_excel("final_snp_qc_maf_missing_refalt.xlsx")

# =========================
# Extracting genotype matrix (markers x samples)
# =========================
meta_cols <- intersect(colnames(df_final), c("CHROM","POS","ID","REF","ALT"))
geno_cols <- setdiff(colnames(df_final), meta_cols)

G_mx <- as.matrix(df_final[, geno_cols, drop = FALSE])  # markers x samples
G_mx <- column_to_rownames(as.data.frame(G_mx), var = "marker_ID")
G_mx <- G_mx[,-1]
G_mx <- as.matrix(G_mx)
storage.mode(G_mx) <- "numeric"

# After imputation, all loci have N individuals (if there were at least 1 call)
N_j <- rowSums(!is.na(G_mx))

# Ploidy for the dosage matrix (2 or 4)
P <- 4

# =========================
# Population allele frequencies per locus: p_hat_j and q_hat_j
# =========================
# p_hat_j = sum_i g_ij / (P * N_j)
p_hat_tetra <- rowSums(G_mx, na.rm = TRUE) / (P * N_j) #tetra
q_hat_tetra <- 1 - p_hat_tetra

# =========================
# Observed heterozygosity per locus (diploidized): H_O,j = n_H,j / N_j
# =========================
# Diploidize dosages into 0/1/2 coding
G_dip <- G_mx
G_dip[!is.na(G_dip) & G_dip == 0] <- 0 #if 0 -> 0
G_dip[!is.na(G_dip) & (G_dip >= 1 & G_dip <= (P-1))] <- 1 #if 1,2,3 -> 1
G_dip[!is.na(G_dip) & G_dip == P] <- 2 #if 4 -> 2 

n_Hj <- rowSums(G_dip ==1, na.rm = T)
H_Oj <- n_Hj / N_j

# =========================
# Moody/“gametic heterozygosity” per locus: H_Moody,j
#    h_ij = 2*g_ij*(P-g_ij) / (P*(P-1))
#    H_Moody,j = mean_i(h_ij) over non-missing individuals at locus j
# =========================
h_ij <- (2 * G_mx * (P - G_mx)) / (P * (P - 1))
H_Moody_j <- rowSums(h_ij, na.rm = TRUE)/ N_j

# =========================
# Expected heterozygosity per locus (Nei gene diversity): H_S,j = 2*p_hat*q_hat
# =========================
#Tetraploid - He (Nei genetic diversity)
H_Sj_tetra <- 2 * p_hat_tetra * q_hat_tetra # <- Tetraploid

#Diploid - allelic frequency
p_hat_diplo <- rowSums(G_dip, na.rm = TRUE) / (2 * N_j) #<- diploid
q_hat_diplo <- 1 - p_hat_diplo
#Diploid - He (Nei genetic diversity)
H_Sj_diplo <- 2 * p_hat_diplo * q_hat_diplo # <- Diploid

# =====================================
# PIC - polymorphic information content --> PIC = 1-(p^2+q^2)-2*(p^2)*(q^2)
# =====================================

#Formula
pic_biallelic <- function(p) {
  q <- 1 - p
  1 - (p^2 + q^2) - 2*(p^2)*(q^2)
}

#tetraploid
N_j_tetra <- rowSums(!is.na(G_mx))
p_tetra <- rowSums(G_mx, na.rm = TRUE) / (P * N_j_tetra) #tetra
PIC_tetra <- pic_biallelic(p_tetra)

#diploid
N_j_dip <- rowSums(!is.na(G_dip))
p_dip <- rowSums(G_dip, na.rm = TRUE) / (2 * N_j_dip) #diplo
PIC_dip <- pic_biallelic(p_dip)

# ===============
# MAF - pmin(p,q)
# ===============

#Tetraploid
MAF_tetra = pmin(p_hat_tetra, 1 - p_hat_tetra)

#Diploid
MAF_diplo = pmin(p_hat_diplo, 1 - p_hat_diplo)

# =====================================================
# RESULTS -- Per-locus results table (population-level)
# =====================================================
per_locus <- df_final %>%
  select(all_of(meta_cols)) %>%
  mutate(
    p_hat_tetra      = p_hat_tetra,
    q_hat_tetra      = q_hat_tetra,
    p_hat_diplo      = p_hat_diplo,
    q_hat_diplo      = q_hat_diplo,
    H_O              = H_Oj,
    H_Moody          = H_Moody_j,
    H_S_tetra        = H_Sj_tetra,
    H_S_diplo        = H_Sj_diplo,
    MAF_tetra        = MAF_tetra,
    MAF_diplo        = MAF_diplo,
    PIC_dip          = PIC_dip,
    PIC_tetra        = PIC_tetra
  )

#saving the table
write_xlsx(per_locus, "per_locus_diversity.xlsx")

#population parameters summary
library(dplyr)
library(tidyr)

per_locus_summary_long <- per_locus %>%
  summarise(across(
    where(is.numeric),
    list(
      mean = ~mean(.x, na.rm = TRUE),
      min  = ~min(.x,  na.rm = TRUE),
      max  = ~max(.x,  na.rm = TRUE)
    ),
    .names = "{.col}__{.fn}"   # <-- safe separator
  )) %>%
  pivot_longer(
    cols = everything(),
    names_to = c("variable", "stat"),
    names_sep = "__"
  ) %>%
  pivot_wider(
    names_from = stat,
    values_from = value
  ) %>%
  arrange(variable)

per_locus_summary_long
write_xlsx(per_locus_summary_long, "population_summary.xlsx")

# =========================
# Per-individual summaries
# =========================
# per-individual counts
N_i  <- colSums(!is.na(G_dip))            # called loci per individual
n_Hi <- colSums(G_dip == 1, na.rm = TRUE) # heterozygous loci per individual

# per-individual observed heterozygosity (diploidized)
H_Oi <- n_Hi / N_i

# H_Moody,i (dosage-aware): mean of h_ij across loci 
H_Moody_i <- colSums(h_ij, na.rm = T) / N_i   

# =========================
# GRM (VanRaden) ### genomic inbreeding F_G = diag(G) - 1
# =========================
library(AGHmatrix)

#Tetraploid
G_VanRaden_tetra <- Gmatrix(t(G_mx), method="VanRaden", ploidy=4)
F_G_Tetra <- diag(G_VanRaden_tetra) - 1

#Diploid
G_VanRaden_diplo <- Gmatrix(t(G_dip), method="VanRaden", ploidy=2)
F_G_diplo <- diag(G_VanRaden_diplo) - 1


per_individual <- data.frame(
  Clone    = colnames(G_mx),
  H_O       = as.numeric(H_Oi),
  H_Moody   = as.numeric(H_Moody_i),
  F_G_diplo = as.numeric(F_G_diplo),
  F_G_Tetra = as.numeric(F_G_Tetra),
  stringsAsFactors = FALSE
)

#saving
write_xlsx(per_individual, "per_individual.xlsx")

# =========================
# Tajima's D in PopGenome (diploidized calling)
# =========================

library(readxl)
library(dplyr)
library(vcfR)

# ---- tyding the things ----
df_final <- as.data.frame(df_final)

#names to keep
new_names <- colnames(G_mx)  
map_new_to_vcf <- setNames(sample_key$vcf_id, sample_key$new_name) # mapa new_name -> vcf_id
vcf_ids_to_keep <- unname(map_new_to_vcf[new_names])

write.table(vcf_ids_to_keep,
  file = "vcf_ids_to_keep.txt",
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)

# positions_to_keep.txt (CHROM POS)
pos <- do.call(rbind, strsplit(rownames(G_mx), "_", fixed = TRUE))
positions_df <- data.frame(CHROM = pos[,1], POS = as.integer(pos[,2]))

write.table(
  positions_df,
  file = "positions_to_keep.txt",
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)


# VCF sample names do NOT match colnames(G_mx), load mapping:
sample_key <- read_excel("sample_names_genotyping_process.xlsx") %>%
  as.data.frame() %>%
  rename(vcf_id = `RAPiD Genomics Sample Code`) %>%
  select(vcf_id, new_name)

# ---- Inputs ----
#install.packages("devtools")
#library(devtools)
#devtools::install_github("pievos101/PopGenome")

library(PopGenome)

vcf_directory <- "C:/Users/givan/OneDrive/Área de Trabalho/Repositorios/diversity_LMH/popgenome_run/"
genome <- readData(vcf_directory, format = "VCF")

# Check the summary of the genome data
print(summary(genome))

# Calculate Tajima's D
genome1 <- neutrality.stats(genome)

# Extract Tajima's D values
tajimas_d <- get.neutrality(genome1, theta = TRUE)$Tajima.D

# Print the Tajima's D value
print(genome1@Tajima.D)

