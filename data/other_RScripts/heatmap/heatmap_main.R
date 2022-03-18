setwd("D:/Halushka_lab/UCSC/CatGenEx-master/CatGenEx-master/")
library(ggplot2) # ggplot() for plotting
library(dplyr) # data reformatting
library(tidyr) # data reformatting
library(stringr) # string manipulation

#FINAL CODE USED FOR CLASS HEATMAP
####################################################################################
####################################################################################
####################################################################################
# CODE USED FOR CLASS FIGURE 4
setwd("D:/Halushka_lab/UCSC/CatGenEx-master/CatGenEx-master/")
library(pheatmap)

m <- read.table("forHeatmap_organ_012722.txt",header=T,stringsAsFactors=F, sep="\t")
m <- read.table("forHeatmap_tissue_012722.txt",header=T,stringsAsFactors=F, sep="\t")
m <- read.table("forHeatmap_tissue_012722_2.txt",header=T,stringsAsFactors=F, sep="\t")

m <- read.table("forHeatmap_marcSelected2_4paperFig_v1.txt",header=T, sep="\t")
samp2 <- m[,-1]
rownames(samp2) <- m[,1]
d <-samp2[,1:ncol(samp2)] + 1
dzt <- as.data.frame(t(d))

tryScale <- (apply(dzt, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X))))*100

pheatmap((as.matrix(tryScale)),cellwidth = 12,cellheight = 18,scale = "none", cluster_row = FALSE, cluster_cols = FALSE,border_color = "black",treeheight_col =0, treeheight_row = 0, fontsize_col = 10 )


m <- read.table("novel_CT_011422i_Q75_HM_tr.txt",header=T, sep="\t")
samp2 <- m[,-1]
rownames(samp2) <- m[,1]
d <-samp2[,1:ncol(samp2)] + 1
p = pheatmap(t(as.matrix(log(d,2))),cellwidth = 18,cellheight = 15,scale = "none", cluster_row = FALSE, cluster_cols = FALSE,border_color = "black")


m <- read.table("epithelial_panelB.txt",header=T, sep="\t")
m <- read.table("epithelial_panelB_v1.txt",header=T, sep="\t")
samp2 <- m[,-1]
rownames(samp2) <- m[,1]
d <-samp2[,1:ncol(samp2)] + 1
p = pheatmap(t(as.matrix(log(d,2))),cellwidth = 23,cellheight = 15,scale = "none", cluster_row = FALSE, cluster_cols = FALSE,border_color = "black")


pheatmap((as.matrix(tryScale)),cellwidth = 12,cellheight = 18,scale = "none", cluster_row = FALSE, cluster_cols = FALSE,border_color = "black",treeheight_col =0, treeheight_row = 0, fontsize_col = 10 )
m <- read.table("Blood_panelC.txt",header=T, sep="\t")
m <- read.table("Blood_panelC_v1.txt",header=T, sep="\t")
samp2 <- m[,-1]
rownames(samp2) <- m[,1]
d <-samp2[,1:ncol(samp2)] + 1
p = pheatmap(t(as.matrix(log(d,2))),cellwidth = 23,cellheight = 15,scale = "none", cluster_row = FALSE, cluster_cols = FALSE,border_color = "black")
p = pheatmap(t(as.matrix(tryScale)),cellwidth = 23,cellheight = 15,scale = "none", cluster_row = FALSE, cluster_cols = FALSE,border_color = "black")
p = pheatmap(t(as.matrix(tryScale)),cellwidth = 20,cellheight = 10,scale = "none", cluster_row = FALSE, cluster_cols = FALSE,border_color = "black")

p = pheatmap(t(as.matrix(tryScale)),cellwidth = 20,cellheight = 10,scale = "none", cluster_row = TRUE, cluster_cols = TRUE,border_color = "black")
p = pheatmap(t(as.matrix(tryScale)),cellwidth = 18,cellheight = 10,scale = "none", cluster_row = FALSE, cluster_cols = FALSE,border_color = "black")
########################################################################################33

library(pheatmap)
m <- read.table("Colon_top10miRs_Tissue.txt",header=T,stringsAsFactors=F, sep="\t")
samp2 <- m[,-1]
rownames(samp2) <- m[,1]
d <-samp2[,1:ncol(samp2)] + 1
#p = pheatmap(t(as.matrix(log(d,2))),cellwidth = 23,cellheight = 13,scale = "none", cluster_row = FALSE, cluster_cols = TRUE,border_color = "black")
p = pheatmap(t(as.matrix(log(d,2))),cellwidth = 23,cellheight = 15,scale = "none", cluster_row = FALSE, cluster_cols = FALSE,border_color = "black")
p = pheatmap(t(as.matrix(log(d,10))),cellwidth = 23,cellheight = 15,scale = "none", cluster_row = FALSE, cluster_cols = FALSE,border_color = "black")

# 
dzt <- as.data.frame(t(d))
# Scale between 0 to 100
tryScale <- (apply(dzt, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X))))*100
p = pheatmap(t(as.matrix(log(tryScale+1,2))),cellwidth = 23,cellheight = 15,scale = "none", cluster_row = FALSE, cluster_cols = FALSE,border_color = "black")
p = pheatmap((as.matrix(log(tryScale+1,2))),cellwidth = 23,cellheight = 15,scale = "none", cluster_row = FALSE, cluster_cols = FALSE,border_color = "black")

########################################################################################33

library(pheatmap)
library(ggplot2)
library(gplots)

m <- read.table("liver_mkhDesignedAvgRPMScaled_values.txt",header=T,stringsAsFactors=F, sep="\t")
m <- read.table("colon_mkhDesignedAvgRPMScaled_values.txt",header=T,stringsAsFactors=F, sep="\t")
m <- read.table("spleen_mkhDesignedAvgRPMScaled_values.txt",header=T,stringsAsFactors=F, sep="\t")
m <- read.table("lymphnode_mkhDesignedAvgRPMScaled_values2.txt",header=T,stringsAsFactors=F, sep="\t")

samp2 <- m[,-1]
rownames(samp2) <- m[,1]

maxv =  max(as.matrix(samp2))
colors = c(seq(0,0.5,length=100),seq(0.6,1.2,length=100),seq(1.21,1.25,length=100))

my_palette <- colorRampPalette(c("white", "red"))(n = 299)
p = pheatmap(t(as.matrix(samp2)),col=my_palette,breaks=colors,cellwidth = 23,cellheight = 15,scale = "none", cluster_row = FALSE, cluster_cols = FALSE,border_color = "black")

