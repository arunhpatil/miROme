# This R script runs DESeq2 VST, ComBat-Seq and other normalizing methods or combinations of the methods on the counts data and plot interactive UMAP plot 
# Install all the required R packages and call the libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# If your computer has restrictions over admin previliages then choose a different library path to install your packages
.libPaths(c("C:/Program Files/R/R-4.0.2/library"))

BiocManager::install("RUVSeq")
library(RUVSeq)
install.packages("DescTools")
library(DescTools)
library(ellipsis)
library(edgeR)
library(plotly)
library(tidyverse)
library(htmlwidgets)
# Attach the `DESeq2` library
library(DESeq2)
# Attach the `umap` library
library(umap)
# Attach the `ggplot2` library for plotting
library(ggplot2)
# We will need this so we can use the pipe: %>%
library(magrittr)
library(sva)
library(pamr)
library(limma)


# Please choose your working directory accordingly. This is a path based on my working directory on my computer
setwd("D:/Halushka_lab/Arun/datasets/microRNAOme2/UMAP/masterDirectory_072921/FINAL_PART2/RUVr_DESEQ2_122221/")

### IMPORTANT - miRNAs to retain (max(RPM) >= 100), list. Basically, we are retainig and plotting UMAP on only ligit miRNAs based on overlap between MirGeneDB and miRBase
required_miRNAs <- read.table("D:/Halushka_lab/Arun/datasets/microRNAOme2/UMAP/masterDirectory_072921/FINAL_PART2/RUVr_DESEQ2_122221/req_miRNA_MirGeneDB_mkh_ge100_YESFilter.csv")

# Reading the miRNA names
mireq <- as.vector(required_miRNAs$V1)

# This part of code is to run RUVg, if you want to test RUVg
zfGenes <- read.table("Complete_primaryCells_FINAL_Counts_021122.csv", row.names=1, header=TRUE, sep = ",")
expression_df2 <- zfGenes
# miRNAs which are ubiquituously expressed
spikes <- c("hsa-let-7a-5p/7c-5p","hsa-let-7f-5p","hsa-miR-103a-3p/107","hsa-miR-125a-5p","hsa-miR-181a-5p",
            "hsa-miR-186-5p","hsa-miR-191-5p","hsa-miR-22-3p","hsa-miR-27a-3p/27b-3p","hsa-miR-30d-5p")
genes <- rownames(expression_df2)[!grepl(paste(spikes,collapse="|"), rownames(expression_df2))]
set.seed(123)
idx <- c(genes, spikes)
expression_df2 <- expression_df2 + 1
seq <- newSeqExpressionSet(as.matrix(expression_df2[idx,]))
seqRUVg <- RUVg(seq, spikes, k=1)
expression_df <- as.data.frame(normCounts(seqRUVg))
write.table(expression_df,"ruvg_normalized.csv", sep=",")


# this is a again reading the counts file, which can be used to test DESeq2 or other normalizing methods
expression_df <- readr::read_csv("Complete_primaryCells_FINAL_Counts_021122.csv") %>%
  # Tuck away the gene ID  column as row names, leaving only numeric values
  tibble::column_to_rownames("miRNA")

# metadata <- readr::read_csv("PrimaryCells_metadata_090921_without_outliers2.csv")
metadata <- readr::read_csv("PrimaryCells_metadata_021122.csv")

expression_df <- expression_df %>%
  dplyr::select(metadata$Run)

# Check if this is in the same order
all.equal(colnames(expression_df), metadata$Run)

# convert the columns we will be using for annotation into factors
metadata <- metadata %>%
  dplyr::select( # select only the columns that we will need for plotting
    Run,
    SRA,
    Assay,
    LibrarySelection,
    CellType,
    Class,
    All_miRNA_Reads,
    batch,
    group,
    Unique_miRNAs
  ) %>%
  dplyr::mutate( # Now let's convert the annotation variables into factors
    refinebio_treatment = factor(
      Class,
      # specify the possible levels in the order we want them to appear
      levels = c("Muscle","Epithelial","Brain","Other","Fibroblast","Immune","RBC","Stem","Endothelial","Fat","Plasma","Platelet","Sperm")
    ),
    refinebio_disease = as.factor(Class) #### Arun what do you want here? replace with in braces!!
  )

# filter and retain only ligit miRNAs 
filtered_expression_df <- expression_df[rownames(expression_df) %in% mireq, ]
filtered_expression_df <- filtered_expression_df + 1
dim(filtered_expression_df)

# Define parameters from metadata
batch = metadata$batch
group = metadata$group

design <- model.matrix(~batch)

#y <- DGEList(counts=counts(as.matrix(filtered_expression_df)), group=metadata$batch)
y <- DGEList(counts=as.matrix(filtered_expression_df), group=metadata$batch)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)

fit <- glmFit(y, design)
res <- residuals(fit, type="deviance")

# RUVr normalization (after UQ)
seq <- as.matrix(filtered_expression_df)
seqUQ <- betweenLaneNormalization(seq, which="upper")
controls <- rownames(seq)
seqRUVr <- RUVr(seqUQ, controls, k=1, res)

filtered_expression_df <- seqRUVr$normalizedCounts 
filtered_expression_df <- filtered_expression_df %>% replace(is.na(.), 0)
write.table(tibble::rownames_to_column(as.data.frame(filtered_expression_df), "miRNA"), "fineNAs_expressionDF_020622.txt", sep = "\t")
# adjusted_counts <- ComBat_seq(as.matrix(filtered_expression_df), batch=batch, group=group, full_mod=FALSE)
# #filtered_expression_df <- adjusted_counts
# #write.table(adjusted_counts,"ruvg_combatSeq_normalized.csv", sep=",")
# #write.table(adjusted_counts,"ruvr_combatSeq_normalized.csv", sep=",")
# write.table(adjusted_counts,"combatSeq_normalized.csv", sep=",")

# Create a `DESeqDataSet` object
dds <- DESeqDataSetFromMatrix(
  #countData = round(filtered_expression_df/3), # the counts values for all samples in our dataset
  countData = filtered_expression_df, # the counts values for all samples in our dataset
  colData = metadata, # annotation data for the samples in the counts data frame
  design = ~ group # Here we are not specifying a model
  # Replace with an appropriate design variable for your analysis
)

# Error in validObject(.Object) :
#   invalid class "DESeqDataSet" object: NA values are not allowed in the count matrix
# In addition: Warning message:
#   In mde(x) : NAs introduced by coercion to integer range
# > max(filtered_expression_df)
# [1] 4812424899
# > as.integer(max(filtered_expression_df))
# [1] NA
## https://stackoverflow.com/questions/64553944/deseq2-na-values-are-not-allowed-error-when-anyis-nacounts-false
## DESeqDataSetFromMatrix(round(counts/2), colData = data.frame(colnames(counts)), design = ~1)
# > max(round(filtered_expression_df/2))
# [1] 2406212450
# > .Machine$integer.max
# [1] 2147483647
# # > 


# Normalize and transform the data in the `DESeqDataSet` object
# using the `vst()` function from the `DESeq2` R package
#
#dds_norm <- vst(dds, nsub=nrow(dds))
dds_norm <- varianceStabilizingTransformation(dds, blind=FALSE)

# First we are going to retrieve the normalized data
# from the `DESeqDataSet` object using the `assay()` function
# 
# adjusted_counts <- ComBat_seq(as.matrix(filtered_expression_df), batch=batch, group=group, full_mod=FALSE)
# filtered_expression_df <- adjusted_counts

# normalized_counts <- assay(dds_norm)
# write.table(normalized_counts, "RUVr_DESeq_batchenr_enh_grp_exall.csv", sep = ",")

normalized_counts <- assay(dds_norm) %>%
  t() # We need to transpose this data so each row is a sample
 
# normalized_counts <- filtered_expression_df %>%
#   t() # We need to transpose this data so each row is a sample

write.table(tibble::rownames_to_column(as.data.frame(normalized_counts), "miRNA"), "normalizedCounts_DESeq_MirGeneDB_Select_miRNAs_021122.txt", sep = "\t", row.names = FALSE)

##################### ALL THE PLOTTING PART HERE BELOW #####################################

#umap_results <- umap::umap(log2(abs(normalized_counts+1)))
# Now perform UMAP on the normalized data
umap_results <- umap::umap(normalized_counts)
# Make into data frame for plotting with `ggplot2`
# The UMAP values we need for plotting are stored in the `layout` element
umap_plot_df <- data.frame(umap_results$layout) %>%
  # Turn sample IDs stored as row names into a column
  tibble::rownames_to_column("Run") %>%
  # Add the metadata into this data frame; match by sample IDs
  dplyr::inner_join(metadata, by = "Run")


# Plot using `ggplot()` function
ggplot(
  umap_plot_df,
  aes(
    x = X1,
    y = X2,
    color = refinebio_treatment # label points with different colors for each `subgroup`
  )
) +
  geom_point(size = 1.5)+labs(col="Class") # This tells R that we want a scatterplot


# Plot using `ggplot()` function and save to an object
final_annotated_umap_plot <- ggplot(
  umap_plot_df,
  aes(
    x = X1,
    y = X2,
    # plot points with different colors for each `refinebio_treatment` group
    color = refinebio_treatment,
    # plot points with different shapes for each `refinebio_disease` groupd
    shape = refinebio_disease,
    text = paste(
      "SRA: ", SRA, "\n",
      "Run: ", Run, "\n",
      "Cell: ", CellType, "\n",
      "Assay: ", Assay, "\n",
      "miRNA reads: ", All_miRNA_Reads, "\n",
      "miRNAs: ", Unique_miRNAs, "\n",
      sep = ""
    )
  )
) +
  scale_shape_manual(values=1:nlevels(metadata$refinebio_disease))+ labs(col="Class", shape="Class")+
  geom_point(size=2) +
  labs(color  = "Guide name", linetype = "Guide name", shape = "Guide name") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


# Display the plot that we saved above
final_annotated_umap_plot
p = ggplotly(final_annotated_umap_plot, tooltip = "text")
p
htmlwidgets::saveWidget(as_widget(p), "DESeq_MirGeneDB_Select_miRNAs_021122.html")
#########################################################################################################
