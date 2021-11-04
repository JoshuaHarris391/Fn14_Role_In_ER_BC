
BC2116_subset <- function(INPUT_DF, CLIN_DF, OUTPUT_CLIN, INPUT_GENES, RETURN_PROBES = F){

  # INPUT_DF <- EXP
  # CLIN_DF <- CLIN
  # INPUT_GENES = c("SCNN1A", "MKI67", "CDH1", "VIM")
  # OUTPUT_CLIN = T

  ########################################
  # Subsetting GE for two genes and clinical
  ########################################
  
  # Find INPUT_DF genes that are in INPUT_GENES
  m <- (INPUT_DF$GENE %in% INPUT_GENES)
  
  
  # Reporting genes not found
  gene_row <- match(INPUT_GENES, INPUT_DF$GENE)
  gene_row <- na.omit(gene_row)
  found_genes <- INPUT_DF[gene_row, "GENE"]
  if (length(found_genes) == length(INPUT_GENES)) {
    print("All Input Genes Found")  
  } else if (length(found_genes) < length(INPUT_GENES)){
    INPUT_GENES[!(INPUT_GENES %in% found_genes)] %>% paste(., "not found", sep = " ") %>% print()
  }
  
  
  # # Saving found gene names as vector
  # df_genes <- INPUT_DF[m, "GENE"]
  
  # Adding probe names to INPUT_DF
  INPUT_DF$PROBE <- rownames(INPUT_DF)
  
  # For each found gene find matching probe names
  probe_list <- list()
  for (gene_name in found_genes) {
    tmp_df <- INPUT_DF[INPUT_DF$GENE == gene_name, ]
    probe_list[[paste0(gene_name)]] <- rownames(tmp_df)
  }
  
  # Subsetting DF for all found probes
  PLOT_DF <- INPUT_DF[m, ] %>% 
    select(., !(GENE) & !(PROBE)) %>% 
    t() %>% 
    as.data.frame()

  
  # Subset PLOT_DF by gene, and find probe with highest avg expression
  optimal_probe <- list()
  for (gene_name in found_genes) {
    # extract probes for gene iteration
    probes <- probe_list[[gene_name]]
    # select probes
    tmp_df <- PLOT_DF %>% select(all_of(probes))
    # Calculate colmeans
    val_col_means <- apply(tmp_df, 2, mean)
    # Find col with highest average, save probe name to list
    optimal_probe[[paste0(gene_name)]] <- val_col_means[grep(max(val_col_means), val_col_means)] %>% names()
  }
  
  # Finding location of optimal probes
  m <- match(paste0(optimal_probe), INPUT_DF$PROBE)
  
  # creating df
  if (OUTPUT_CLIN == T) {
    PLOT_DF <- INPUT_DF[m, ] %>% 
    select(., !(GENE) & !(PROBE)) %>% 
      t() %>% 
      as.data.frame()
    colnames(PLOT_DF) <- names(optimal_probe)
  } else if (OUTPUT_CLIN == F){
    PLOT_DF <- INPUT_DF[m, ]
    rownames(PLOT_DF) <- PLOT_DF$GENE
    PLOT_DF <- select(PLOT_DF, !(GENE) & !(PROBE))
  }
  
  
  # Adding in clinical data if specified
  if (OUTPUT_CLIN == T) {
    # Adding pateint ID
    PLOT_DF$TITLE_ID <- rownames(PLOT_DF) %>% factor(., levels(CLIN_DF$TITLE_ID))
    # Combining Clinical DF with Plot DF
    PLOT_DF <- inner_join(PLOT_DF, CLIN_DF, by = "TITLE_ID")
  }
  
  # Returning DFs
  if(RETURN_PROBES == T){
    return(as.data.frame(optimal_probe))
  } else {
    return(PLOT_DF)
  }
  
  
}
  



