setwd("D:/Halushka_lab/Arun/datasets/microRNAOme2/UMAP/masterDirectory_072921/andrea/selectCells_fibroblasts_021322_n2406/")
### Euclidean distance classification with crossvalidation

# _100RPM_qualitymiRs, Run one input file at a time. 
gene_counts=read.csv('raw_MirGeneDB_overlapOnly.csv', header=FALSE)
gene_counts=read.csv('combatSeq_normalized_100RPM_qualitymiRs.csv', header=FALSE)
gene_counts=read.csv('combatSeq_DESeq_100RPM_qualitymiRs.csv', header=FALSE)
gene_counts=read.csv('DESeq_100RPM_qualitymiRs.csv', header=FALSE)
gene_counts=read.csv('ruvg_normalized_100RPM_qualitymiRs.csv', header=FALSE)
gene_counts=read.csv('ruvg_combatSeq_100RPM_qualitymiRs.csv', header=FALSE)
gene_counts=read.csv('RUVg_combatSeq_DESeq_100RPM_qualitymiRs.csv', header=FALSE)
gene_counts=read.csv('RUVg_DESeq_100RPM_qualitymiRs.csv', header=FALSE)
gene_counts=read.csv('ruvr_normalized_100RPM_qualitymiRs.csv', header=FALSE)
gene_counts=read.csv('ruvr_combatSeq_normalized_100RPM_qualitymiRs.csv', header=FALSE)
gene_counts=read.csv('RUVr_combatSeq_DESeq_100RPM_qualitymiRs.csv', header=FALSE)
gene_counts=read.csv('RUVr_DESeq_100RPM_qualitymiRs.csv', header=FALSE)


gene_counts <- gene_counts+1
#define tissue type for each sample in gene_counts
tissue_type=read.csv('tissue_type_miome.csv', header=FALSE)

#define crossvalidation folds (could be leave one study out, or stratified by tissue type, etc)
##for leave one study out, unique values of CV_fold are 1:n_studies, with all samples from a given study assigned the same value
CV_fold=read.csv('CV_microRNAome.csv', header=FALSE)

#CV_fold[42:363,1]=8
#CV_fold[364:374,1]=9
#CV_fold[1242:1393,1]=58


n_fold=max(CV_fold)

gene_counts_all=gene_counts
tissue_type_all=tissue_type


tab_mat_tissue=as.data.frame(table(tissue_type_all))
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
accuracy


#table(n_tiss_truth,n_tiss_predict)
temp <- table(n_tiss_truth,n_tiss_predict)
colnames(temp) <- tab_mat_tissue$tissue_type_all
#rowSums(temp)
rownames(temp)<- tab_mat_tissue$tissue_type_all
temp

 sens <- sensitivity_by_tiss
 rownames(sens) <- tab_mat_tissue$tissue_type_all
 sens <- sens * 100
 sens


