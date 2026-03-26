## =============================
## Creating the manhattan plots with only significant markers
## =============================
# Givanildo Rodrigues da Silva
# 03/25/2026

library(readxl)

all_sig_markers <- read_xlsx("circos_plot_data/all_points_data.xlsx")
all_sig_markers <- all_sig_markers[,-c(3,6)]

library(dplyr)
library(ggplot2)
library(tibble)
library(patchwork)

# =========================================================
# 1) Chromosome lengths
# =========================================================
chr_info <- tribble(
  ~CHROM,   ~start, ~end_chr,
  "chr01",      1, 88495810,
  "chr02",      1, 46064988,
  "chr03",      1, 60649747,
  "chr04",      1, 69075021,
  "chr05",      1, 55203317,
  "chr06",      1, 59059153,
  "chr07",      1, 57620787,
  "chr08",      1, 59174247,
  "chr09",      1, 67528160,
  "chr10",      1, 60840980,
  "chr11",      1, 46727314,
  "chr12",      1, 59501733
) %>%
  mutate(
    CHROM   = factor(CHROM, levels = CHROM),
    chr_len = end_chr - start + 1,
    offset  = lag(cumsum(chr_len), default = 0),
    center  = offset + chr_len / 2
  )

# =========================================================
# 2) Colors
# =========================================================
base_cols <- c(
  "pcadapt"   = "gray40",  # verde-água
  "poolfstat" = "gray40",  # mostarda
  "xp-ehh"    = "gray40"   # azul claro
)

dual_col   <- "#F28E2B"     # compartilhado por 2 métodos
triple_col <- "red"     # compartilhado por 3 métodos

approach_levels <- c("pcadapt", "poolfstat", "xp-ehh")

# =========================================================
# 3) Build plotting table using ONLY CHROM + end
# =========================================================
df <- all_sig_markers %>%
  transmute(
    approach  = factor(approach, levels = approach_levels),
    CHROM     = factor(CHROM, levels = levels(chr_info$CHROM)),
    pos       = as.numeric(end),
    score     = as.numeric(neglogq),
    marker_ID = marker_ID
  ) %>%
  filter(!is.na(CHROM), !is.na(pos)) %>%
  distinct(approach, CHROM, pos, .keep_all = TRUE) %>%
  left_join(
    chr_info %>% select(CHROM, offset, center, chr_len),
    by = "CHROM"
  ) %>%
  mutate(
    bp_cum = pos + offset
  ) %>%
  arrange(CHROM, pos, approach)

# =========================================================
# 4) Count exact support by SNP position ONLY
#    No intervals, no windows
# =========================================================
support_tbl <- df %>%
  distinct(approach, CHROM, pos) %>%
  count(CHROM, pos, name = "n_models")

df_plot <- df %>%
  left_join(support_tbl, by = c("CHROM", "pos")) %>%
  mutate(
    hit_class = case_when(
      n_models == 3 ~ "3 models",
      n_models == 2 ~ "2 models",
      TRUE          ~ "Exclusive"
    ),
    hit_class = factor(hit_class, levels = c("Exclusive", "2 models", "3 models"))
  )

# =========================================================
# 5) Diagnostics
#    Shared SNPs must have one single cumulative position
# =========================================================
diag_shared <- df_plot %>%
  filter(hit_class != "Exclusive") %>%
  group_by(CHROM, pos) %>%
  summarise(
    n_methods = n_distinct(approach),
    n_bp_cum  = n_distinct(bp_cum),
    methods   = paste(sort(unique(as.character(approach))), collapse = ", "),
    .groups   = "drop"
  ) %>%
  arrange(CHROM, pos)

diag_problem <- diag_shared %>% filter(n_bp_cum != 1)

cat("\n================ Diagnostic ================\n")
cat("Shared SNPs with inconsistent cumulative x:\n")
print(diag_problem)
cat("===========================================\n\n")

# Optional: inspect all shared SNPs
shared_snps <- df_plot %>%
  filter(hit_class != "Exclusive") %>%
  arrange(CHROM, pos, approach) %>%
  select(approach, CHROM, pos, bp_cum, marker_ID, hit_class)

print(shared_snps, n = 200)

# =========================================================
# 6) Global x-axis settings
#    SAME limits for all Manhattan plots
# =========================================================
genome_max <- max(chr_info$offset + chr_info$chr_len)
x_breaks   <- chr_info$center
x_labels   <- gsub("chr", "", as.character(chr_info$CHROM))

# =========================================================
# 7) Plot function
# =========================================================
make_panel <- function(appr, ylab, show_x = FALSE) {
  
  dsub <- df_plot %>% filter(approach == appr)
  
  ggplot(dsub, aes(x = bp_cum, y = score)) +
    
    # Exclusive SNPs
    geom_point(
      data = dsub %>% filter(hit_class == "Exclusive"),
      color = base_cols[appr],
      size = 2.5,
      alpha = 0.85
    ) +
    
    # Shared by 2 methods
    geom_point(
      data = dsub %>% filter(hit_class == "2 models"),
      color = dual_col,
      size = 2.8,
      alpha = 0.95
    ) +
    
    # Shared by 3 methods
    geom_point(
      data = dsub %>% filter(hit_class == "3 models"),
      shape = 21,
      fill = triple_col,
      color = "black",
      stroke = 0.40,
      size = 2.8,
      alpha = 1
    ) +
    
    scale_x_continuous(
      breaks = x_breaks,
      labels = x_labels,
      limits = c(1, genome_max),
      expand = c(0, 0)
    ) +
    
    labs(
      title = appr,
      x = if (show_x) "Chromosome" else NULL,
      y = ylab
    ) +
    
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0),
      axis.title.y = element_text(size = 11),
      axis.title.x = element_text(size = 12),
      axis.text.x  = if (show_x) element_text(size = 9) else element_blank(),
      axis.ticks.x = if (show_x) element_line() else element_blank(),
      axis.line.x  = if (show_x) element_line() else element_blank(),
      plot.margin  = margin(5.5, 8, 5.5, 5.5)
    )
}

# =========================================================
# 8) Build final panels
# =========================================================
p1 <- make_panel(
  appr   = "pcadapt",
  ylab   = expression(-log[10](q-value)),
  show_x = FALSE
)

p2 <- make_panel(
  appr   = "poolfstat",
  ylab   = "FST(95th percentile)",
  show_x = FALSE
)

p3 <- make_panel(
  appr   = "xp-ehh",
  ylab   = "XP-EHH score",
  show_x = TRUE
)

final_manhattan <- p1 / p2 / p3 +
  plot_layout(heights = c(1, 1, 1))

final_manhattan

# =========================================================
# 9) Save
# =========================================================

ggsave(
  filename = "three_manhattan_exact_overlap.png",
  plot = final_manhattan,
  width = 13,
  height = 8,
  dpi = 600,
  bg = "white"
)
