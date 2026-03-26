## =================
## PCA 
## =================
# Givanildo Rodrigues da Silva
# 03/25/2026


library(AGHmatrix)
library(ggplot2)

# Genomic relationship matrix
G_VanRaden <- Gmatrix(
  as.matrix(t(Gmx)),
  method = "VanRaden",
  ploidy = 4,
  impute.method = "mean"
)

# PCA
res.pca <- prcomp(G_VanRaden, scale. = FALSE)

# Variance explained
var_exp <- (res.pca$sdev^2 / sum(res.pca$sdev^2)) * 100

# PCA table
pca_scores <- data.frame(
  ID  = rownames(res.pca$x),
  PC1 = res.pca$x[, 1],
  PC2 = res.pca$x[, 2],
  stringsAsFactors = FALSE
)

# Align metadata by genotype name
correct_classes2 <- correct_classes[match(pca_scores$ID, correct_classes$newest_name), ]

# Check alignment
stopifnot(all(pca_scores$ID == correct_classes2$newest_name))

# Add metadata
pca_scores$Cluster <- factor(correct_classes2$Cluster_STRUCTURE, levels = c("1", "2", "3", "4"))
pca_scores$Origin  <- factor(correct_classes2$Origin, levels = c("Brazilian", "Foreign"))

# Colors
group_colors <- c(
  "1" = "green3",
  "2" = "blue",
  "3" = "red",
  "4" = "orange"
)

# Shapes
shape_values <- c(
  "Brazilian" = 16,
  "Foreign"   = 17
)

# Plot
ggplot(pca_scores, aes(x = PC1, y = PC2, color = Cluster, shape = Origin)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(values = group_colors) +
  scale_shape_manual(values = shape_values) +
  labs(
    x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
    y = paste0("PC2 (", round(var_exp[2], 1), "%)"),
    color = "Cluster",
    shape = "Origin"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )



