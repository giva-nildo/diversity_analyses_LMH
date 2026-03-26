# ===================================
# Estimating FST betweeen groups 
# ===================================
# Givanildo Rodrigues 
# 03/25/2026

install.packages("StAMPP")
library(StAMPP)

G_matrix <- t(Gmx) #Original matrix


recode_tetra <- function(x) {
  ifelse(is.na(x), NA,
         ifelse(x == 0, "A/A/A/A",
                ifelse(x == 1, "A/A/A/B",
                       ifelse(x == 2, "A/A/B/B",
                              ifelse(x == 3, "A/B/B/B",
                                     ifelse(x == 4, "B/B/B/B", NA))))))
}

# obtaining the Nei distance


##
geno4 <- as.data.frame(apply(G_matrix, 2, recode_tetra))
rownames(geno4) <- rownames(G_matrix)

gi4 <- df2genind(
  X = geno4,
  sep = "/",
  ncode = 1,
  ind.names = rownames(geno4),
  ploidy = 4,
  type = "codom"
)

labels <- read_excel("C:/Users/givan/OneDrive/Área de Trabalho/Repositorios/diversity_LMH/sample_names_genotyping_process.xlsx")
labels <- labels[, c(6, 8, 9)]
labels$structure <- as.factor(labels$structure)
labels$origin <- as.factor(labels$origin)

labels$new_name <- as.character(labels$new_name)

idx <- match(indNames(gi4), labels$new_name)

if (any(is.na(idx))) {
  stop("Some individuals in gi4 were not found in labels$new_name: ",
       paste(indNames(gi4)[is.na(idx)], collapse = ", "))
}

pop(gi4) <- droplevels(labels$structure[idx]) # altern the origin and structure!

metadata <- data.frame(
  pop = pop(gi4),
  row.names = indNames(gi4)
)
metadata

# Step 5: Convert to genlight object
genlight_obj <- new("genlight", gi4@tab)
pop(genlight_obj) <- pop(gi4)

# Check the structure of the genlight object
print(genlight_obj)

# Step 6: Prepare the genlight object for StAMPP
genlight_obj_stampp <- stamppConvert(genlight_obj, type = "genlight")

# Step 7: Compute pairwise Fst
fst_results <- stamppFst(genlight_obj_stampp, nboots = 1000, percent = 95, nclusters = 3)

# Print Fst results
print(fst_results)

# Extracting pairwise Fst values
pairwise_fst <- fst_results$Fsts
print(pairwise_fst)
