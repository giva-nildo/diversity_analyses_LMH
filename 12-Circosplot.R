# =========================================================
# Circular dendrogram + STRUCTURE membership using circlize
# Better layout:
#   outer labels
#   outer STRUCTURE ring
#   thin assigned-cluster ring
#   inner dendrogram
# =========================================================
# Givanildo Rodrigues
# 03/25/2026

library(circlize)
library(dendextend)
library(readxl)
library(dplyr)

# -------------------------
# 1) Build tree
# -------------------------

gi4 #The genotype

nei_dis <- nei.dist(gi4@tab, warning = TRUE)
x <- as.matrix(nei_dis)
heatmap(x) # just visualizing


hc <- hclust(as.dist(x), method = "ward.D2")
dend <- as.dendrogram(hc)

lab_order <- labels(dend)
max_h <- attr(dend, "height")


K_dend <- 4 # choose dendrogram K

# colors for dendrogram groups
dend_cols <- c("forestgreen", "blue3", "orange3", "red3")[1:K_dend]
names(dend_cols) <- as.character(1:K_dend)

# cut tree
cl_dend <- cutree(hc, k = K_dend)
cl_dend <- cl_dend[lab_order]

# color dendrogram branches
dend_k <- color_branches(dend, k = K_dend, col = dend_cols)

# -------------------------
# 2) Read STRUCTURE table
# -------------------------
membership_struc <- read_excel(
  "C:/Users/givan/OneDrive/LMH/Painel batatas - LMH.xls",
  sheet = "IDs",
  range = "M1:T451"
)

membership_struc <- membership_struc[, -7]
membership_struc <- membership_struc[-140, ]

Q <- membership_struc %>%
  transmute(
    label    = newest_name,
    K1       = `1`,
    K2       = `2`,
    K3       = `3`,
    K4       = `4`,
    assigned = factor(Cluster_STRUCTURE)
  ) %>%
  filter(label %in% lab_order)

# reorder to match dendrogram order
Q <- Q[match(lab_order, Q$label), , drop = FALSE]
stopifnot(all(Q$label == lab_order))

Qm <- as.matrix(Q[, c("K1", "K2", "K3", "K4")])
mode(Qm) <- "numeric"
rownames(Qm) <- Q$label

assigned <- as.character(Q$assigned)

# -------------------------
# 3) Colors
# -------------------------
q_cols <- c(
  K1 = "green3",
  K2 = "blue2",
  K3 = "red2",
  K4 = "orange2"
)

assign_cols <- c(
  "1" = "green3",
  "2" = "blue2",
  "3" = "red2",
  "4" = "orange2"
)

# -------------------------
# 4) Options
# -------------------------
# label mode: "none", "selected", or "all"
label_mode <- "selected"

# manually chosen labels
selected_samples <- c("Atlantic", "Michune_Negra", "Antarctita", "Grao_Mogol", "Maria_Bonita", "Barna", "Slaney", "Baraka", "Monalisa", "Caesar", "Sifra", "Prince", "Agata", "Mondial", "Aracy_Ruiva", "Aracy", "Aziza", "Markies", "Asterix", "Ana", "Cupido", "Russet_Burbank")
selected_samples <- intersect(selected_samples, Q$label)

# automatic example:
# selected_samples <- Q$label[apply(Qm, 1, max) < 0.70]

show_assigned_ring   <- TRUE
show_membership_edge <- TRUE

#Starting to save
pdf("circos_plot.pdf", width = 10, height = 10)

# label size
cex_selected <- 0.55
cex_all      <- 0.18   # for many samples this will still be tiny

# track heights
h_labels    <- if (label_mode == "none") 0.001 else if (label_mode == "all") 0.10 else 0.14
h_structure <- 0.16
h_assigned  <- 0.035
h_dend      <- 0.50

# -------------------------
# 5) Plot
# -------------------------
circos.clear()
par(mar = c(1, 1, 1, 8))

circos.par(
  start.degree = 90,
  gap.after = 2,
  cell.padding = c(0, 0, 0, 0),
  track.margin = c(0.002, 0.002),
  points.overflow.warning = FALSE,
  canvas.xlim = c(-1.25, 1.75),
  canvas.ylim = c(-1.25, 1.25)
)

# ONE single sector spanning all individuals
circos.initialize(
  factors = "all",
  xlim = c(0, nrow(Q))
)

# -------------------------
# 5A) Outer label track
# -------------------------
if (label_mode != "none") {
  circos.trackPlotRegion(
    ylim = c(0, 1),
    bg.border = NA,
    track.height = h_labels,
    panel.fun = function(x, y) {
      
      for (i in seq_len(nrow(Q))) {
        lab <- Q$label[i]
        
        draw_this <- FALSE
        lab_cex <- cex_selected
        
        if (label_mode == "all") {
          draw_this <- TRUE
          lab_cex <- cex_all
        }
        
        if (label_mode == "selected" && lab %in% selected_samples) {
          draw_this <- TRUE
          lab_cex <- cex_selected
        }
        
        if (draw_this) {
          circos.text(
            x = i - 0.5,
            y = 0,
            labels = lab,
            facing = "clockwise",
            niceFacing = TRUE,
            adj = c(0, 0.5),
            cex = lab_cex,
            font = ifelse(label_mode == "selected", 2, 1)
          )
        }
      }
    }
  )
}

# -------------------------
# 5B) STRUCTURE membership ring
# -------------------------
circos.trackPlotRegion(
  ylim = c(0, 1),
  bg.border = NA,
  track.height = h_structure,
  panel.fun = function(x, y) {
    
    for (i in seq_len(nrow(Qm))) {
      bottom <- 0
      
      for (k in seq_len(ncol(Qm))) {
        top <- bottom + Qm[i, k]
        
        circos.rect(
          xleft   = i - 1,
          ybottom = bottom,
          xright  = i,
          ytop    = top,
          col     = q_cols[colnames(Qm)[k]],
          border  = NA
        )
        
        bottom <- top
      }
      
      # outline of the WHOLE membership bar
      if (show_membership_edge) {
        circos.rect(
          xleft   = i - 1,
          ybottom = 0,
          xright  = i,
          ytop    = 1,
          col     = NA,
          border  = "gray3",
          lwd     = 0.55
        )
      }
    }
  }
)

# -------------------------
# 5C) Thin assigned-cluster ring
# -------------------------
if (show_assigned_ring) {
  circos.trackPlotRegion(
    ylim = c(0, 1),
    bg.border = NA,
    track.height = h_assigned,
    panel.fun = function(x, y) {
      
      for (i in seq_len(nrow(Q))) {
        circos.rect(
          xleft   = i - 1,
          ybottom = 0,
          xright  = i,
          ytop    = 1,
          col     = assign_cols[assigned[i]],
          border  = "white",
          lwd     = 0.25
        )
      }
    }
  )
}

# -------------------------
# 5D) Inner dendrogram track
# -------------------------
circos.trackPlotRegion(
  ylim = c(0, max_h),
  bg.border = NA,
  track.height = h_dend,
  panel.fun = function(x, y) {
    circos.dendrogram(
      dend_k,
      facing = "outside",
      max_height = max_h,
      use_x_attr = FALSE
    )
  }
)

# -------------------------
# 6) Legends
# -------------------------
legend(
  "right",
  inset = c(-0.02, 0),
  legend = colnames(Qm),
  fill = q_cols,
  border = NA,
  bty = "n",
  cex = 1.0,
  title = "STRUCTURE"
)

legend(
  "right",
  inset = c(-0.02, -0.28),
  legend = names(assign_cols),
  fill = assign_cols,
  border = NA,
  bty = "n",
  cex = 1.0,
  title = "Assigned cluster"
)

dev.off() # Finishing the saving at the directory

# -------------------------
# 7) Clean up
# -------------------------
circos.clear()

# ------------
# Saving
# ------------

ggsave(
  filename = "circos_plot.png",
  plot = p_final,
  width = 12,
  height = 12,
  units = "in",
  dpi = 600,
  bg = "white"
)





