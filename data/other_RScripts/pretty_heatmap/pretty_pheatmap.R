library(pheatmap)
library("dplyr")

# Load DESeq2 VST normalized data
data2 <- read.delim("heatmap_Evolution_deseq2data_subset.txt", header=T, row.names="miRNA")
data_subset <- data2[2:ncol(data2)]

# miRNA and locus information
my_gene_col <- read.delim("heatmapMetadataEvolution.txt", header=T, row.names="miRNA")

# Sample and cell class information
my_sample_col <- read.delim("class_evolution_subset.txt", header=T, row.names="miRNA")

# Select required cell class
required_data <- c("Stem", "Immune","Platelet","Endothelial","RBC","Epithelial","Brain","Plasma","Sperm","Fat","Fibroblast","Muscle")

# Sub select the data 
my_sample_col_subset <- my_sample_col %>% filter((sample %in% required_data))
my_sample_col_subset <-my_sample_col_subset[order(my_sample_col_subset$sample),,drop=FALSE]
data_subset <- data_subset[, rownames(my_sample_col_subset)]
data_subset2 = data_subset[,(colnames(data_subset) %in% rownames(my_sample_col_subset))]

# Plot pretty heatmap
pheatmap(data_subset2, annotation_row = my_gene_col, annotation_col = my_sample_col_subset, cluster_rows  = FALSE, cluster_cols  = F, show_rownames=F, show_colnames=F)
