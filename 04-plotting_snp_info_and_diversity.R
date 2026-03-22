##########################################
### Plotting the SNP details and diversity
##########################################
# Givanildo Rodrigues da Silva
#February, 26th, 2026

##

###################
#### Calculating the mean and distance among the SNPs
###################

library(dplyr)

df_final <- read_excel("final_snp_qc_maf_missing_refalt.xlsx")

# SNPs per chromosome
snps_per_chrom <- df_final %>%
  group_by(CHROM) %>%
  summarise(total_snp_chr = n())

snps_per_chrom
snps_per_chrom$total_snp_chr
mean(snps_per_chrom$total_snp_chr)


# df_final has CHROM and POS columns :contentReference[oaicite:1]{index=1}
snp_dist <- df_final %>%
  mutate(POS = as.integer(POS)) %>%
  arrange(CHROM, POS) %>%
  group_by(CHROM) %>%
  mutate(dist_bp = POS - lag(POS)) %>%      # distance to previous SNP on same chromosome
  filter(!is.na(dist_bp) & dist_bp >= 0) %>%
  ungroup()

# Mean distance overall (bp)
mean_dist_overall <- mean(snp_dist$dist_bp)

# Mean distance per chromosome
mean_dist_by_chr <- snp_dist %>%
  group_by(CHROM) %>%
  summarise(
    n_snps = n() + 1,
    mean_dist_bp = mean(dist_bp),
    .groups = "drop"
  )

mean_dist_overall
mean_dist_by_chr







# Calcular a média de distância entre os SNPs por cromossomo
mean_distance_per_chrom <- df_final %>%
  arrange(CHROM, POS) %>%
  group_by(CHROM) %>%
  summarise(mean_distance = mean(diff(POS), na.rm = TRUE))

mean(mean_dist_by_chr$mean_dist_bp)

# Exibir os resultados
print(snps_per_chrom) #SNPs per chromosomes
mean(snps_per_chrom$mean_snps)
min(snps_per_chrom$mean_snps)
max(snps_per_chrom$mean_snps)

print(mean_distance_per_chrom) #mean distance per chromosome
mean(mean_distance_per_chrom$mean_distance)

###############
####### Plotting the percentage by snps distance
###############
library(dplyr)
library(ggplot2)

# Calcular a distância entre os SNPs
results_merged <- df_final %>%
  arrange(CHROM, POS) %>%
  group_by(CHROM) %>%
  mutate(distance = c(NA, diff(POS))) %>%
  ungroup()

# Remover NAs
results_merged <- results_merged %>%
  filter(!is.na(distance))

# Definir os intervalos de distância
breaks <- c(0, 1000, 10000, 100000, Inf)
labels <- c("<1 kb", "1-10 kb", "10-100 kb", ">100 kb")

# Agrupar as distâncias
results_merged <- results_merged %>%
  mutate(distance_group = cut(distance, breaks = breaks, labels = labels, right = FALSE))

# Calcular a porcentagem de SNPs em cada intervalo
distance_summary <- results_merged %>%
  group_by(distance_group) %>%
  summarise(count = n()) %>%
  mutate(percentage = (count / sum(count)) * 100)

# Plotar as informações
ggplot(distance_summary, aes(x = distance_group, y = percentage)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = " ",
       x = "Marker distance",
       y = "SNPs percentage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Ploting the same information with labels
# Plotar as informações
ggplot(distance_summary, aes(x = distance_group, y = percentage)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%.2f%%", percentage)), vjust = -0.3, size = 4) +
  theme_minimal() +
  labs(title = " ",
       x = "Marker distance",
       y = "Percentage of SNPs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

########

# Plotting a SNP density 

#############
#### plotting a SNP density
#############

install.packages("CMplot")
library(CMplot)
library(metan)

#genodf <- (genodf[!duplicated(genodf$Marker), ]) #Removing duplicates
genodf <- column_to_rownames(genodf, var = "Marker")
genodf <- rownames_to_column(genodf, var = "Marker")
dim(df_final)
table(df_final$CHROM)


#genodf <- rownames_to_column(genodf, var = "Marker") #rownames to column
genodf <- df_final %>% 
  mutate_at(vars(marker_ID, CHROM, POS),
            as.factor)

genodf <- genodf %>% 
  relocate(marker_ID, everything())

#Creating the SNP density plot
genodf$CHROM <- gsub("^chr", "", genodf$CHROM) #Removing the prefix "chr"
CMplot(genodf[,1:3], plot.type="d",bin.size=1e6,chr.den.col=c("blue", "red"),
       file="jpg", #or tiff
       dpi=300)


########################################
### Plotting the populational parameters
########################################

library(dplyr)
library(tidyr)
library(ggplot2)

pyramid_plot <- function(df, left_col, right_col, breaks = NULL,
                         binwidth = NULL, n_bins = 20,
                         proportion = FALSE, title = NULL,
                         digits = 2) {
  
  stopifnot(left_col %in% names(df), right_col %in% names(df))
  
  x <- df %>%
    dplyr::select(dplyr::all_of(c(left_col, right_col))) %>%
    tidyr::pivot_longer(everything(), names_to = "group", values_to = "value") %>%
    dplyr::filter(!is.na(value)) 
  
  x <- x %>% mutate(group = factor(group, levels = c(left_col, right_col)))
  cols <- setNames(c("#00BFC4", "#F8766D"), c(left_col, right_col))
  
  rng <- range(x$value, na.rm = TRUE)
  
  if (is.null(breaks)) {
    if (!is.null(binwidth)) {
      breaks <- seq(floor(rng[1] / binwidth) * binwidth,
                    ceiling(rng[2] / binwidth) * binwidth,
                    by = binwidth)
    } else {
      breaks <- seq(rng[1], rng[2], length.out = n_bins + 1)
    }
  }
  
  # ---- rounded labels for bins ----
  fmt <- function(z) sprintf(paste0("%.", digits, "f"), z)
  bin_labels <- paste0(fmt(breaks[-length(breaks)]), "-", fmt(breaks[-1]))
  
  x <- x %>%
    dplyr::mutate(bin = cut(value,
                            breaks = breaks,
                            include.lowest = TRUE,
                            right = FALSE,
                            labels = bin_labels)) %>%
    dplyr::count(group, bin, name = "n")
  
  if (proportion) {
    x <- x %>% dplyr::group_by(group) %>% dplyr::mutate(n = n / sum(n)) %>% dplyr::ungroup()
    ylab <- "Proportion"
  } else {
    ylab <- "Count"
  }
  
  x <- x %>% dplyr::mutate(n_signed = ifelse(group == left_col, -n, n))
  
  ggplot2::ggplot(x, ggplot2::aes(x = bin, y = n_signed, fill = group)) +
    ggplot2::geom_col(width = 0.9) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(labels = function(z) abs(z)) +
    ggplot2::labs(
      title = title,
      x = NULL,
      y = ylab,
      fill = NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::scale_fill_manual(values = cols) +
    theme(legend.position = "none")
}

colnames(per_locus)

p1 <- pyramid_plot(per_locus, "H_S_tetra", "H_S_diplo",
             title = expression(H[S-Tetra]~vs~H[S-Diplo]), n_bins = 20, digits = 2) 
p1 <- p1 + theme(axis.title.x = element_blank())

p2 <- pyramid_plot(per_locus, "H_Moody", "H_O",
                   title = expression(H[Moody]~vs~H[O]), n_bins = 20, digits = 2)
p2 <- p2 + theme(axis.title.x = element_blank())

p3 <- pyramid_plot(per_locus, "MAF_tetra", "MAF_diplo",
                         title = expression(MAF[Tetra]~vs~MAF[Diplo]), n_bins = 20, digits = 2)
p3 <- p3 + theme(axis.title.x = element_blank())

p4 <- pyramid_plot(per_locus, "PIC_tetra", "PIC_dip",
                   title = expression(PIC[Tetra]~vs~PIC[Diplo]), n_bins = 20, digits = 2)
p4 <- p4 + theme(axis.title.x = element_blank())

#library(patchwork)
#(p1 | p2) / (p3 | p4) + plot_layout(guides = "collect")

library(cowplot)
final <- plot_grid(p1, p2, p3, p4, ncol = 2, labels = c("a", "b", "c", "d"))
final <- ggdraw(final) +
  draw_label("", x = 0.5, y = 0.01, vjust = 0) +
  theme(plot.margin = margin(5.5, 5.5, 20, 5.5))
final <- ggdraw(final) +
  draw_label("Count", x = 0.5, y = 0.01, vjust = 0) +
  theme(plot.margin = margin(5.5, 5.5, 20, 5.5))
final

ggsave("pyramids.png", final, width = 7.0, height = 10, dpi = 600)

#End

