# =================
# AMOVA 
# =================
# Givanildo Rodrigues 
#03/25/2026

######################
###### Converting the matrix to genind object
#####################

library(adegenet)
library(dplyr)
library(readxl)
library(poppr)
library(pegas)

G_matrix <- t(Gmx) #Original matrix

G_matrix_imp <- G_matrix #THat it will be used for imputation by mean

#imputation by mean
for (j in 1:ncol(G_matrix_imp)) {
  m <- round(mean(G_matrix_imp[, j], na.rm = TRUE))
  G_matrix_imp[is.na(G_matrix_imp[, j]), j] <- m
}


recode_tetra <- function(x) {
  ifelse(is.na(x), NA,
         ifelse(x == 0, "A/A/A/A",
                ifelse(x == 1, "A/A/A/B",
                       ifelse(x == 2, "A/A/B/B",
                              ifelse(x == 3, "A/B/B/B",
                                     ifelse(x == 4, "B/B/B/B", NA))))))
}

# obtaining the Nei distance

geno4 <- as.data.frame(apply(G_matrix_imp, 2, recode_tetra))
rownames(geno4) <- rownames(G_matrix_imp)

gi4 <- df2genind(
  X = G_matrix_imp,
  sep = "/",
  ncode = 1,
  ind.names = rownames(G_matrix_imp),
  ploidy = 4,
  type = "codom"
)

dist_matrix <- nei.dist(G_matrix_imp, warning = TRUE) #The Nei Distance
dist_matrix #Nei distance matrix

### Now will be performed the AMOVA
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

pop(gi4) <- droplevels(labels$structure[idx])

metadata <- data.frame(
  pop = pop(gi4),
  row.names = indNames(gi4)
)

amova_result <- pegas::amova(
  dist_matrix ~ pop, #AMOVA formula
  data = metadata, #MEtadata with pop info
  nperm = 1000 # number of permutations for hypotesis testing
)

amova_result