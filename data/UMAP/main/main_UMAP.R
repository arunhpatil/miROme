# This R script runs DESeq2 VST on select class and cell types 
# Install all the required R packages and call the libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("RUVSeq")
library(RUVSeq)
#install.packages("DescTools")
library(DescTools)
library(ellipsis)
library(edgeR)
library(plotly)
library(tidyverse)
library(htmlwidgets)
library(DESeq2)
library(umap)
library(ggplot2)
library(magrittr)
library(sva)
library(pamr)
library(limma)
library(dplyr)

source('main_umap_ggplotly.R')

## Step 1: To show all individual cell types on the legends ##
load("data.rda")
#all cell class
req_levels = c("Muscle","Epithelial","Brain","Other","Fibroblast","Immune","RBC","Stem","Endothelial","Fat","Plasma","Platelet","Sperm")
# Check if this is in the same order
all.equal(colnames(expression_df), metadata$Sample)
# Extract rows whos sum is >= 5000 across samples
filtered_expression_df <- expression_df %>%
  dplyr::filter(rowSums(.) >= 5000)
p <- main_umap_ggplotly(metadata, expression_df, filtered_expression_df)
htmlwidgets::saveWidget(as_widget(p), "DESEq2_allClass_CellType.html")



source('main_umap_ggplotly_class.R')
## Step 2: To show only cell class as the legends ##
load("data.rda")
#all cell class
req_levels = c("Muscle","Epithelial","Brain","Other","Fibroblast","Immune","RBC","Stem","Endothelial","Fat","Plasma","Platelet","Sperm")
# Check if this is in the same order
all.equal(colnames(expression_df), metadata$Sample)
# Extract rows whos sum is >= 5000 across samples
filtered_expression_df <- expression_df %>%
  dplyr::filter(rowSums(.) >= 5000)
p <- main_umap_ggplotly_class(metadata, expression_df, filtered_expression_df)
htmlwidgets::saveWidget(as_widget(p), "DESEq2_allClass_Class.html")


