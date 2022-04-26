
umap_ggplotly <- function(metadata, expression_df, filtered_expression_df){
  
  # convert the columns we will be using for annotation into factors
  metadata <- metadata %>%
    dplyr::select( # select only the columns that we will need for plotting
      Sample,
      SRA,
      Assay,
      LibrarySelection,
      CellType,
      Class,
      All_miRNA_Reads,
      batch,
      group,
      group2,
      Unique_miRNAs
    ) %>%
    dplyr::mutate( # Now let's convert the annotation variables into factors
      refinebio_treatment = factor(
        Class,
        levels =  req_levels
      ),
      refinebio_disease = as.factor(CellType) #### Arun what do you want here? replace with in braces!!
    )
  
  metadata$group2 <- as.factor(metadata$group2)
  metadata$group <- as.factor(metadata$group)
  #print(class(metadata$group))
  # Create a `DESeqDataSet` object
  dds <- DESeqDataSetFromMatrix(
    countData = filtered_expression_df, # the counts values for all samples in our dataset
    colData = metadata, # annotation data for the samples in the counts data frame
    design = ~ group # Here we are not specifying a model
  )
  
  dds_norm <- varianceStabilizingTransformation(dds, blind=FALSE)
  
  normalized_counts <- assay(dds_norm) %>%
    t() # We need to transpose this data so each row is a sample
  
  umap_results <- umap::umap(normalized_counts)
  # Make into data frame for plotting with `ggplot2`
  # The UMAP values we need for plotting are stored in the `layout` element
  umap_plot_df <- data.frame(umap_results$layout) %>%
    # Turn sample IDs stored as row names into a column
    tibble::rownames_to_column("Sample") %>%
    # Add the metadata into this data frame; match by sample IDs
    dplyr::inner_join(metadata, by = "Sample")
  
  
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
        "Sample: ", Sample, "\n",
        "Cell: ", CellType, "\n",
        "Assay: ", Assay, "\n",
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
  final_annotated_umap_plot
  p = ggplotly(final_annotated_umap_plot, tooltip = "text")
  p
    
  lvls  <- levels(umap_plot_df$refinebio_disease)
  cols  <- setNames(hcl(h=seq(15, 375, length=length(lvls)+1), l=65, c=100),lvls)
  cols[c("Adipocyte","Preadipocyte","Lipocyte","Red_blood_cell")] <- c("#000000","#000000","#000000", "#000000")
  z <- final_annotated_umap_plot + scale_color_manual(values=cols)+
    labs(color  = "Guide name", linetype = "Guide name", shape = "Guide name")
  p = ggplotly(z, tooltip = "text")
  p
  return(p)
}
