# The aim of this script is to extract the raw data from the RDS object 
# into a csv, to be later loaded in the notebook

setwd("/media/seth/SETH_DATA/SETH_Alex/kranocito/data/leinroth")


# Install Seurat (for R==3.6.3 the last supported Seurat version is v3.2.3)
install.packages('remotes')  # to install seurat later
remotes::install_version("Seurat")

library(Seurat)

# Load the RDS
rds <- readRDS("GSE200234_RAW/GSM6025297_leinrothtomato.rds")

# extract counts
rds_data <- rds[["RNA"]]@counts

rds_data_dense <- as.data.frame(as.matrix(rds_data))

write.csv(rds_data_dense, "leinroth.csv")
write.csv(rds@meta.data, "leinroth_metadata.csv")
