.libPaths(c("C:/Program Files/R/R-4.0.2/library"))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# https://alexslemonade.github.io/refinebio-examples/03-rnaseq/dimension-reduction_rnaseq_02_umap.html
setwd("D:/Halushka_lab/Arun/datasets/microRNAOme2/UMAP/analysis_071721/Endothelial_fat_RBCs/")
# https://www.littlemissdata.com/blog/interactiveplots 
if (!("DESeq2" %in% installed.packages())) {
  # Install DESeq2
  BiocManager::install("DESeq2", update = FALSE)
}

if (!("umap" %in% installed.packages())) {
  # Install umap package
  BiocManager::install("umap", update = FALSE)
}

# Attach the `DESeq2` library
library(DESeq2)

# Attach the `umap` library
library(umap)

# Attach the `ggplot2` library for plotting
library(ggplot2)

# We will need this so we can use the pipe: %>%
library(magrittr)
# Set the seed so our results are reproducible:
set.seed(12345)

tumap_me <- function(input_Key){
  print(input_Key)
  str(input_Key)
}


install.packages("plotly")
install.packages("tidyverse")
install.packages("htmlwidgets")


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


metadata <- readr::read_csv("endothelial_fat_RBCs_metadata.csv")

expression_df <- readr::read_csv("endothelial_fat_RBCs.csv") %>%
  # Tuck away the gene ID  column as row names, leaving only numeric values
  tibble::column_to_rownames("miRNA")

expression_df <- expression_df %>%
  dplyr::select(metadata$Run)

# Check if this is in the same order
all.equal(colnames(expression_df), metadata$Run)

# change CellType to Class and CellDetail to CellType
###############################################################################################################
# convert the columns we will be using for annotation into factors
metadata <- metadata %>%
  dplyr::select( # select only the columns that we will need for plotting
    Run,
    SRA,
    Assay,
    CellType,
    Class,
    All_miRNA_Reads,
    Unique_miRNAs
  ) %>%
  dplyr::mutate( # Now let's convert the annotation variables into factors
    refinebio_treatment = factor(
      Class,
      # specify the possible levels in the order we want them to appear
      #levels = c("FALSE","CD4_lymphocyte_TRUE","Macrophage_TRUE","Dendritic_cell_TRUE","CD15_cell_TRUE","CD235a_cell_TRUE","CD14_cell_TRUE","CD8_lymphocyte_TRUE","CD56_cell_TRUE","Mononuclear_immune_cell_TRUE","CD19_lymphocyte_TRUE")
      #levels = c("Muscle","Epithelial","Brain","Other","Fibroblast","Immune","RBC","Stem","Endothelial","Fat","Plasma","Platelet","Sperm")
      levels = c("Endothelial","Fat","RBC")
      #levels = c("Skeletal_muscle_cell","Keratinocyte","Melanocyte","Keratinocyte_neonatal","Astrocyte","Smooth_muscle_cell_bladder","Smooth_muscle_cell_aorta","Myoblast","Mesangial_cell","Stromal_cell_prostate","Chondrocyte","Osteoblast","CD14_cell","CD56_cell","CD34_cell","T_lymphocyte","CD19_lymphocyte","CD10_cell","Dendritic_cell","CD15_cell","Monocyte","Neuron","Meningeal","Endothelial_cell_umbilical_vein","Endothelial_cell_aortic","Endothelial_cell_retinal_microvascular","Endothelial_cell_brain_microvascular","Endothelial_cell_brain","Endothelial_cell_microvascular","Endothelial_cell","Endothelial_cell_capillary","Endothelial_cell_internal_thoracic","Endothelial_cell_sinusoidal","Endothelial_cell_arterial","Endothelial_cell_kidney","Endothelial_cell_lymphatic","Endothelial_progenitor_cell","Nasal_epithelial_cell","Nasal_polyp_epithelial_cell","Bronchial_epithelial_cell","Small_airway_epithelial_cell","Cervical_epithelial_cell","Retinal_pigment_epithelial_cell","Beta_cell","Hepatocyte","Thymic_epithelial_cell","Conjunctival_epithelial_cell","Beta_cell_derived","Beta_cell_like_derived","Sebocyte","Urothelial_cell","Prostate_epithelial_cell","Retinal_pigment_epithelial_cell_fetal","Pancreas_epithelial_like","Breast_epithelial_cell","Retinal_pigment_epithelial_cell_derived","Colonic_epithelial","Intestinal_epithelial_cell","Renal_epithelial_cell","Renal_cortical_epithelial_cell","Renal_proximal_tubule_epithelial_cell","Tracheal_epithelial_cell","Placenta_epithelial_cell","Esophagus_epithelial_cell","Small_intestinal_epithelial_cell","Gum_epithelial_cell","Amniotic_epithelial_cell","Eye_corneal_epithelial_cell","Eye_ciliated_epithelial_cell","Preadipocyte","Mesenchymal_stromal_cell","Fibroblast_endometrial_stroma","Fibroblast_skin","Fibroblast_aortic_adventitial","Fibroblast_lymph_node","Fibroblast_breast","Fibroblast_gum_fetal","Fibroblast_periodontal_fetal","Fibroblast_lung_fetal","Stromal_cell_pancreas","Fibroblast_foreskin","Stellate_cell","Fibroblast_dermal","Myofibroblast","Fibroblast_dermal_fetal","Fibroblast_periodontal_ligament","Fibroblast","Fibroblast_ventricular_cardiac","Fibroblast_foreskin_neonatal","Fibroblast_dermal_derived","Fibroblast_heart_fetal","Fibroblast_choroid_plexus","Fibroblast_eye","Fibroblast_lung","Fibroblast_periodontal_ligament_fetal","Fibroblast_pulmonary_artery","Fibroblast_trophoblast","Macrophage","Mononuclear_immune_cell","CD4_lymphocyte","B_cell_germinal_center","CD8_lymphocyte","Red_blood_cell","CD27_IgG_cell","CD27_IgA_cell","CD27-_IgD_cell","CD27_IgD_cell","Thymocyte_CD4_CD8","Thymocyte_CD34","Neutrophil","B_cell_naive","Mast_cell","Hematopoietic_stem_cell","CD5_CD19_cell","Natural_killer_cell","Lymphocyte","Megakaryocyte_derived","Plasma_cell","Centroblast","Centrocyte","B_cell_pre_germinal_center","B_lymphycyte","Cardiomyocyte_derived","Smooth_muscle_cell_umbilical_artery","Smooth_muscle_subclavian_artery","Smooth_muscle_cell_brain","Smooth_muscle_cell_pulmonary_artery","Satellite_cell","Smooth_muscle_cell_brachiocephalic_artery","Smooth_muscle_umbilical_artery","Smooth_muscle_cell_vascular","Myotube","Smooth_muscle_cell","Muscle_cell","Smooth_muscle_cell_colon","Smooth_muscle_cell_coronary_artery","Smooth_muscle_cell_prostate","Smooth_muscle_cell_carotid_artery","Smooth_muscle_cell_esophagus","Smooth_muscle_cell_uterus","Smooth_muscle_cell_internal_thoracic","Smooth_muscle_cell_subclavian_artery","Smooth_muscle_cell_lung","Theca_cell","Schwann_cell","Osteocyte","Pericyte_brain","Nucleus_pulposus_cell","Sertoli_cell","Mural_granulosa_cell","Trabecular_meshwork_cell_eye","Chondroblast","Synovial_cell","Mesangioblast_derived","Granulosa_cell","Germ_cell","Mesothelial_cell","Plasma","Platelet","Erythroblast_derived","Sperm","Embryonic_stem_cell_H9","Bone_marrow_skeletal_stem_cell","Embryonic_stem_cell","iPSC_fibroblast","iPSC","iPSC_umbilical_cord_blood","Cardiac_progenitor_derived","Neural_progenitor_cell_derived","Mesoderm_progenitor_derived","Neural_stem_cell","Ventral_midbrain_derived","Embryonic_stem_cell_H7","H9_differentiated","Ectoderm_precursor_derived","Mesoderm_precursor_derived","Endoderm_precursor_derived","iPSC_bone_marrow","iPSC_skin","Blastocyst_derived","iPSC_foreskin","Embryonic_stem_cell_H1","iPSC_amniotic_fluid","Neuronal_stem_cell","Neural_stem_cell_derived","Neuroepithelial_stem_cell","Stem_cell_adipose_derived","Bone_marrow_mesenchymal_stem_cell","Mesenchymal_stem_cell_derived","Dental_pulp_stem_cell","Mesenchymal_stem_cell","Mesenchymal_stem_cell_derived_umbilical_cord","Mesenchymal_stem_cell_derived_bone_marrow","Mesenchymal_stem_cell_derived_adipose_tissue","Mesenchymal_stem_cell_derived_liver","Mesenchymal_stem_cell_derived_amnion","Mesenchymal_stem_cell_derived_spinal_cord","Unrestricted_somatic_stem_cell","HS578T_CD24+_cell","HS578T_CD44+_cell")
    ),
    refinebio_disease = as.factor(CellType) #### Arun what do you want here? replace with in braces!!
  )

## only for using normalized values for UMAP of signel class vs fat and macrophages
normalized_counts <- expression_df %>%
  t()
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
    color = refinebio_disease,
    # plot points with different shapes for each `refinebio_disease` groupd
    shape = refinebio_disease,
    text = paste(
      "SRA: ", SRA, "\n",
      "Run: ", Run, "\n",
      "Cell: ", CellType, "\n",
      "miRNA reads: ", All_miRNA_Reads, "\n",
      "miRNAs: ", Unique_miRNAs, "\n",
      sep = ""
    )
  )
) +
  scale_shape_manual(values=1:nlevels(metadata$refinebio_disease))+ 
  geom_point(size=2) +
  labs(color  = "Guide name", linetype = "Guide name", shape = "Guide name")+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



# Display the plot that we saved above
#final_annotated_umap_plot


p = ggplotly(final_annotated_umap_plot, tooltip = "text")


lvls  <- levels(umap_plot_df$refinebio_disease)
cols  <- setNames(hcl(h=seq(15, 375, length=length(lvls)+1), l=65, c=100),lvls)
# cols[c("Adipocyte","Preadipocyte","Sperm")] <- c("#808080","#808080", "#808080")
# cols[c("Adipocyte","Preadipocyte","Sperm")] <- c("#36454F","#36454F", "#808080")
cols[c("Adipocyte","Preadipocyte","Lipocyte","Red_blood_cell")] <- c("#000000","#000000","#000000", "#000000")
z <- final_annotated_umap_plot + scale_color_manual(values=cols)+
  labs(color  = "Guide name", linetype = "Guide name", shape = "Guide name")
p = ggplotly(z, tooltip = "text")
p
htmlwidgets::saveWidget(as_widget(p), "DESeq_endothelialFatRBCs_020822.html")

