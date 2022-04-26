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

source('umap_ggplotly.R')

## Step 1: Brain Cells ##
load("data.rda")
metadata <-metadata[metadata$Class %in% c('Brain', 'Fat') | metadata$CellType == 'Red_blood_cell', ]
req_levels = c("Brain","Fat","RBC")
#EXTRACT DATA FOR SAMPLES BELONG TO PARTICULAR CLASS OF CELLS ALONG WITH FAT AND RBC's 
expression_df <-expression_df[, c(metadata$Sample)]
# Check if this is in the same order
all.equal(colnames(expression_df), metadata$Sample)
# Extract rows whos sum is >= 5000 across samples
filtered_expression_df <- expression_df %>%
  dplyr::filter(rowSums(.) >= 5000)
#Call the function to plot umap and return ggplotly html file.
p <- umap_ggplotly(metadata, expression_df, filtered_expression_df)
#Save the file
htmlwidgets::saveWidget(as_widget(p), "SuppFig_1_Brain.html")


## Step 2: Endothelial Cells ##
load("data.rda")
metadata <-metadata[metadata$Class %in% c('Endothelial', 'Fat') | metadata$CellType == 'Red_blood_cell', ]
req_levels = c("Endothelial","Fat","RBC")
#EXTRACT DATA FOR SAMPLES BELONG TO PARTICULAR CLASS OF CELLS ALONG WITH FAT AND RBC's 
expression_df <-expression_df[, c(metadata$Sample)]
# Check if this is in the same order
all.equal(colnames(expression_df), metadata$Sample)
# Extract rows whos sum is >= 5000 across samples
filtered_expression_df <- expression_df %>%
  dplyr::filter(rowSums(.) >= 5000)
#Call the function to plot umap and return ggplotly html file.
p <- umap_ggplotly(metadata, expression_df, filtered_expression_df)
#Save the file
htmlwidgets::saveWidget(as_widget(p), "SuppFig_2_Endothelial.html")


## Step 3: Epithelial Cells ##
load("data.rda")
metadata <-metadata[metadata$Class %in% c('Epithelial', 'Fat') | metadata$CellType == 'Red_blood_cell', ]
req_levels = c("Epithelial","Fat","RBC")
#EXTRACT DATA FOR SAMPLES BELONG TO PARTICULAR CLASS OF CELLS ALONG WITH FAT AND RBC's 
expression_df <-expression_df[, c(metadata$Sample)]
# Check if this is in the same order
all.equal(colnames(expression_df), metadata$Sample)
# Extract rows whos sum is >= 5000 across samples
filtered_expression_df <- expression_df %>%
  dplyr::filter(rowSums(.) >= 5000)
#Call the function to plot umap and return ggplotly html file.
p <- umap_ggplotly(metadata, expression_df, filtered_expression_df)
#Save the file
htmlwidgets::saveWidget(as_widget(p), "SuppFig_3_Epithelial.html")


## Step 4: Fibroblast ##
load("data.rda")
metadata <-metadata[metadata$Class %in% c('Fibroblast', 'Fat') | metadata$CellType == 'Red_blood_cell', ]
req_levels = c("Fibroblast","Fat","RBC")
#EXTRACT DATA FOR SAMPLES BELONG TO PARTICULAR CLASS OF CELLS ALONG WITH FAT AND RBC's 
expression_df <-expression_df[, c(metadata$Sample)]
# Check if this is in the same order
all.equal(colnames(expression_df), metadata$Sample)
# Extract rows whos sum is >= 5000 across samples
filtered_expression_df <- expression_df %>%
  dplyr::filter(rowSums(.) >= 5000)
#Call the function to plot umap and return ggplotly html file.
p <- umap_ggplotly(metadata, expression_df, filtered_expression_df)
#Save the file
htmlwidgets::saveWidget(as_widget(p), "SuppFig_4_Fibroblasts.html")


## Step 5: IMMUNE CELLS ##
load("data.rda")
metadata <-metadata[metadata$Class %in% c('Immune', 'Fat') | metadata$CellType == 'Red_blood_cell', ]
req_levels = c("Immune","Fat","RBC")
#EXTRACT DATA FOR SAMPLES BELONG TO PARTICULAR CLASS OF CELLS ALONG WITH FAT AND RBC's 
expression_df <-expression_df[, c(metadata$Sample)]
# Check if this is in the same order
all.equal(colnames(expression_df), metadata$Sample)
# Extract rows whos sum is >= 5000 across samples
filtered_expression_df <- expression_df %>%
  dplyr::filter(rowSums(.) >= 5000)
#Call the function to plot umap and return ggplotly html file.
p <- umap_ggplotly(metadata, expression_df, filtered_expression_df)
#Save the file
htmlwidgets::saveWidget(as_widget(p), "SuppFig_5_Immune.html")


## Step 6: Muscle cells ##
load("data.rda")
metadata <-metadata[metadata$Class %in% c('Muscle', 'Fat') | metadata$CellType == 'Red_blood_cell', ]
req_levels = c("Muscle","Fat","RBC")
#EXTRACT DATA FOR SAMPLES BELONG TO PARTICULAR CLASS OF CELLS ALONG WITH FAT AND RBC's 
expression_df <-expression_df[, c(metadata$Sample)]
# Check if this is in the same order
all.equal(colnames(expression_df), metadata$Sample)
# Extract rows whos sum is >= 5000 across samples
filtered_expression_df <- expression_df %>%
  dplyr::filter(rowSums(.) >= 5000)
#Call the function to plot umap and return ggplotly html file.
p <- umap_ggplotly(metadata, expression_df, filtered_expression_df)
#Save the file
htmlwidgets::saveWidget(as_widget(p), "SuppFig_6_Muscle.html")


## Step 7: Other cells ##
load("data.rda")
metadata <-metadata[metadata$Class %in% c('Other', 'Fat') | metadata$CellType == 'Red_blood_cell', ]
req_levels = c("Other","Fat","RBC")
#EXTRACT DATA FOR SAMPLES BELONG TO PARTICULAR CLASS OF CELLS ALONG WITH FAT AND RBC's 
expression_df <-expression_df[, c(metadata$Sample)]
# Check if this is in the same order
all.equal(colnames(expression_df), metadata$Sample)
# Extract rows whos sum is >= 5000 across samples
filtered_expression_df <- expression_df %>%
  dplyr::filter(rowSums(.) >= 5000)
#Call the function to plot umap and return ggplotly html file.
p <- umap_ggplotly(metadata, expression_df, filtered_expression_df)
#Save the file
htmlwidgets::saveWidget(as_widget(p), "SuppFig_7_Other.html")


## Step 8: Plasma ##
load("data.rda")
metadata <-metadata[metadata$Class %in% c('Plasma', 'Fat') | metadata$CellType == 'Red_blood_cell', ]
req_levels = c("Plasma","Fat","RBC")
#EXTRACT DATA FOR SAMPLES BELONG TO PARTICULAR CLASS OF CELLS ALONG WITH FAT AND RBC's 
expression_df <-expression_df[, c(metadata$Sample)]
# Check if this is in the same order
all.equal(colnames(expression_df), metadata$Sample)
# Extract rows whos sum is >= 5000 across samples
filtered_expression_df <- expression_df %>%
  dplyr::filter(rowSums(.) >= 5000)
#Call the function to plot umap and return ggplotly html file.
p <- umap_ggplotly(metadata, expression_df, filtered_expression_df)
#Save the file
htmlwidgets::saveWidget(as_widget(p), "SuppFig_8_Plasma.html")


## Step 9: Platelet ##
load("data.rda")
metadata <-metadata[metadata$Class %in% c('Platelet', 'Fat') | metadata$CellType == 'Red_blood_cell', ]
req_levels = c("Platelet","Fat","RBC")
#EXTRACT DATA FOR SAMPLES BELONG TO PARTICULAR CLASS OF CELLS ALONG WITH FAT AND RBC's 
expression_df <-expression_df[, c(metadata$Sample)]
# Check if this is in the same order
all.equal(colnames(expression_df), metadata$Sample)
# Extract rows whos sum is >= 5000 across samples
filtered_expression_df <- expression_df %>%
  dplyr::filter(rowSums(.) >= 5000)
#Call the function to plot umap and return ggplotly html file.
p <- umap_ggplotly(metadata, expression_df, filtered_expression_df)
#Save the file
htmlwidgets::saveWidget(as_widget(p), "SuppFig_9_Platelet.html")


## Step 10: Sperm ##
load("data.rda")
metadata <-metadata[metadata$Class %in% c('Sperm', 'Fat') | metadata$CellType == 'Red_blood_cell', ]
req_levels = c("Sperm","Fat","RBC")
#EXTRACT DATA FOR SAMPLES BELONG TO PARTICULAR CLASS OF CELLS ALONG WITH FAT AND RBC's 
expression_df <-expression_df[, c(metadata$Sample)]
# Check if this is in the same order
all.equal(colnames(expression_df), metadata$Sample)
# Extract rows whos sum is >= 5000 across samples
filtered_expression_df <- expression_df %>%
  dplyr::filter(rowSums(.) >= 5000)
#Call the function to plot umap and return ggplotly html file.
p <- umap_ggplotly(metadata, expression_df, filtered_expression_df)
#Save the file
htmlwidgets::saveWidget(as_widget(p), "SuppFig_10_Sperm.html")


## Step 11: Stem cells ##
load("data.rda")
metadata <-metadata[metadata$Class %in% c('Stem', 'Fat') | metadata$CellType == 'Red_blood_cell', ]
req_levels = c("Stem","Fat","RBC")
#EXTRACT DATA FOR SAMPLES BELONG TO PARTICULAR CLASS OF CELLS ALONG WITH FAT AND RBC's 
expression_df <-expression_df[, c(metadata$Sample)]
# Check if this is in the same order
all.equal(colnames(expression_df), metadata$Sample)
# Extract rows whos sum is >= 5000 across samples
filtered_expression_df <- expression_df %>%
  dplyr::filter(rowSums(.) >= 5000)
#Call the function to plot umap and return ggplotly html file.
p <- umap_ggplotly(metadata, expression_df, filtered_expression_df)
#Save the file
htmlwidgets::saveWidget(as_widget(p), "SuppFig_11_Stem.html")
