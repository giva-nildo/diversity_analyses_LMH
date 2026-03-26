## ====================================================================
## Normalized Mutual Information (NMI) index between sctructure and HC
## ====================================================================
# Givanildo Rodrigues
# 03/25/2026

cl_dend1 <- as.data.frame(cl_dend)
cl_dend1 <- rownames_to_column(cl_dend1, var = "newest_name")
cl_dend1$newest_name <- as.factor(cl_dend1$newest_name)

membership_struc <- as.data.frame(membership_struc)
membership_struc$newest_name <- as.factor(membership_struc$newest_name)

classes_struc_hc <- inner_join(cl_dend1, membership_struc, by = "newest_name")
#The clusters should be verified
writexl::write_xlsx(classes_struc_hc, "classes_struc_hc.xlsx") # it will be fixed manually in the Excel program

## calling the corrected file
correct_classes <- read_excel("classes_struc_hc.xlsx")
library(aricode)
nmi_score <- NMI(correct_classes$cl_dend_corrected, correct_classes$Cluster_STRUCTURE) 
nmi_score