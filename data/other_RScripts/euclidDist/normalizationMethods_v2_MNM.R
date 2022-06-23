if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#Call all the required libraries
library(RUVSeq)
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

# load expression counts matrix, tissue_type, CV_fold, metadata, spikes 
# NOTE: Start here
load("mainFiles_v2.rda")

#Load function to run DESeq2 VST that returns normalized counts
deseq_call <- function(input_df, metadata){
  # Create a `DESeqDataSet` object
  dds <- DESeqDataSetFromMatrix(
    countData = input_df, 
    colData = metadata, # annotation data for the samples in the counts data frame
    design = ~ group # Here we are not specifying a model
  )
  dds_norm <- varianceStabilizingTransformation(dds, blind=FALSE)
  normalized_counts <- assay(dds_norm)
  return(normalized_counts)
}

#Function to evaluate best normalization methods 
euclidean_measure <- function(gene_counts, CV_fold, tissue_type){
  n_fold=max(CV_fold)
  gene_counts_all=gene_counts
  tissue_type_all=tissue_type
  
  tab_mat_tissue=as.data.frame(table(tissue_type_all))
  names(tab_mat_tissue)[1] <- "tissue_type_all" # required for R version 4.2.0
  tissue_types=tab_mat_tissue$tissue_type_all
  n_tissues=dim(tab_mat_tissue)[1]
  n_samples=dim(gene_counts_all)[2]
  n_genes=dim(gene_counts_all)[1]
  
  n_tiss_truth=matrix(, nrow=n_samples, ncol=1)
  n_tiss_predict=matrix(, nrow=n_samples, ncol=1)
  for (ind_fold in 1:n_fold){
    
    #create training data
    idx=which(CV_fold!=ind_fold);
    #idx_test_arun=which(CV_fold!=ind_fold);
    #print(idx_test_arun)
    gene_counts_train=gene_counts_all[, idx];
    tissue_types_all_train=tissue_type_all[idx,1];
    
    #create test data
    idx_test=which(CV_fold==ind_fold)
    
    gene_counts_test=as.data.frame(gene_counts_all[, idx_test])
    tissue_types_all_test=tissue_type_all[idx_test,1]
    
    avg_tiss_count=matrix(, nrow=n_genes, ncol=n_tissues)
    for (ind in 1:n_tissues) {
      tissue_in=tissue_types[ind]
      idx=which(tissue_types_all_train==tissue_in);
      if (length(idx)==0){
        avg_tiss_count[, ind]=NaN
      } else {
        avg_tiss_count[,ind]=rowMeans(gene_counts_train[,idx])
      }
    }
    
    
    for (ind in 1:length(tissue_types_all_test)){
      str_in=tissue_types_all_test[ind]
      idx=which(tab_mat_tissue$tissue_type_all==str_in)
      n_tiss_truth[idx_test[ind],1]=idx
    }
    
    #print(idx_min)
    for (ind in 1:dim(gene_counts_test)[2]){
      sample_in=gene_counts_test[,ind]
      V  = sample_in - avg_tiss_count
      dist_all=matrix(, nrow=1, ncol=dim(V)[2])
      for (ind2 in 1:dim(V)[2]){
        dist_all[ind2] = sqrt(t(V[,ind2]) %*% V[,ind2]);
      }
      
      idx_min=which(dist_all==min(dist_all))
      n_tiss_predict[idx_test[ind],1]=idx_min
    }
    #print(idx_min)
  }
  
  sensitivity_by_tiss=matrix(, nrow=max(n_tiss_truth), ncol=1)
  for (ind in 1:max(n_tiss_truth)){
    idx_true=which(n_tiss_truth==ind);
    n_true=length(idx_true);
    predicted_vals=n_tiss_predict[idx_true];
    idx=which(predicted_vals==ind);
    sensitivity_by_tiss[ind,1]=length(idx)/n_true
  }
  #the rows in this matrix correspond to the rows in "tab_mat_tissue"
  
  idx=which(n_tiss_predict==n_tiss_truth);
  accuracy=length(idx)/length(n_tiss_predict)
  print(paste("Accuracy:", round(accuracy*100,1)))
  cat("\n")
  print("Predictability:")
  #table(n_tiss_truth,n_tiss_predict)
  temp <- table(n_tiss_truth,n_tiss_predict)
  colnames(temp) <- tab_mat_tissue$tissue_type_all
  rownames(temp)<- tab_mat_tissue$tissue_type_all
  print(temp)
  cat("\n")
  sens <- sensitivity_by_tiss
  rownames(sens) <- tab_mat_tissue$tissue_type_all
  sens <- sens * 100
  colnames(sens) <- "Sensitivity of prediction"
  print(sens)
}

# Step 1: Evaluate raw counts 
euclidean_measure(filtered_expression_df, CV_fold, tissue_type)

# Step 2: Evaluate log2 transformed raw counts 
log2Counts <- log2(filtered_expression_df + 1)
euclidean_measure(log2Counts, CV_fold, tissue_type)

# Step 3: Evaluate method - RUVg
expression_df2<- filtered_expression_df 
genes <- rownames(expression_df2)[!grepl(paste(spikes,collapse="|"), rownames(expression_df2))]
set.seed(123)
idx <- c(genes, spikes)
seq <- newSeqExpressionSet(as.matrix(expression_df2[idx,]))
seqRUVg <- RUVg(seq, spikes, k=1, isLog = FALSE)
expression_df <- as.data.frame(normCounts(seqRUVg))
log2expression_df <- log2(expression_df+1) ## MNM: need to log these before using Andrea's code
euclidean_measure(log2expression_df, CV_fold, tissue_type) ## updated to use logged data

# Step 4: Evaluate method - RUVg + DESeq2 and VST
normalized_counts <- deseq_call(expression_df, metadata)
euclidean_measure(normalized_counts, CV_fold, tissue_type)


# Defining batch, groups and design for RUVr and CombatSeq
batch = as.factor(metadata$batch)
group = metadata$group
design <- model.matrix(~batch)

# Step 4: Evaluate method - DESeq2 and VST
rm(filtered_expression_df)
load("mainFiles_v2.rda")
normalized_counts <- deseq_call(filtered_expression_df, metadata)
euclidean_measure(normalized_counts, CV_fold, tissue_type)

# Step 5: Evaluate method - CombatSeq
rm(filtered_expression_df)
load("mainFiles_v2.rda")
adjusted_counts <- ComBat_seq(as.matrix(filtered_expression_df), batch=batch, group=group, full_mod=FALSE)
euclidean_measure(adjusted_counts, CV_fold, tissue_type)
adjusted_counts2 <- ComBat_seq(as.matrix(filtered_expression_df), batch=batch, full_mod=FALSE)
identical(adjusted_counts, adjusted_counts2) 
log2adjusted_counts2 <- log2(adjusted_counts2+1) ## MNM: need to log these before using Andrea's code
euclidean_measure(log2adjusted_counts2, CV_fold, tissue_type) ## updated to use logged counts

# Step 6: Evaluate method - CombatSeq + DESeq2 and VST
rm(filtered_expression_df)
load("mainFiles_v2.rda")
normalized_counts <- deseq_call(adjusted_counts2, metadata)
euclidean_measure(normalized_counts, CV_fold, tissue_type)

# Step 7: Evaluate method - RUVr
rm(filtered_expression_df)
load("mainFiles_v2.rda")
filtered_expression_df <- filtered_expression_df + 1
y <- DGEList(counts=(as.matrix(filtered_expression_df)), group=metadata$batch)
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
log2filtered_expression_df <- log2(filtered_expression_df + 1) ## MNM: need to log these before using Andrea's code
euclidean_measure(log2filtered_expression_df, CV_fold, tissue_type) ## updated to use log counts

# Step 8: Evaluate method - RUVr + DESeq2 and VST
normalized_counts <- deseq_call(round(filtered_expression_df/4), metadata)
euclidean_measure(normalized_counts, CV_fold, tissue_type)
# IF you get an error of NA, countData = round(filtered_expression_df/4), You should use 
# normalized_counts <- deseq_call(round(filtered_expression_df/4), metadata)


## Few key notes to consider:
# First, the group variable in the metadata file is coded as numeric:
#   > load("mainFiles.rda")
# > class(metadata$group)
# [1] "numeric"
# This should be a factor, which can easily be done via:
#   metadata$group <- as.factor(metadata$group)


# Second, in several places you add 1 to all the counts. For example:
#   expression_df2<- filtered_expression_df + 1
# This is going to mess with the modeling for all the methods that expect counts (with some zeros as well). The only time you need to add a pseudo count is when taking the log, e.g. when calculating log(counts + 1).
