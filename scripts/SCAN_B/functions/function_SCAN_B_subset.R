SCAN_B_subset <- function(INPUT_DF, CLIN_DF, OUTPUT_CLIN = T, INPUT_GENES, CALCULATE_QUANTILE = F){
  
  # library(tidyverse)
  # # Input_DF
  # INPUT_DF <- GSE96058_GE
  # # Input Clinical DF
  # CLIN_DF <- GSE96058_Clinical_DF
  # # Input gene
  # INPUT_GENES = c("SCNN1A", "MKI67", "CDH1", "VIM")
  # INPUT_GENES = prolif_meta_gene_list
  # CALCULATE_QUANTILE = "SCNN1A"
  # OUTPUT_CLIN = T
  
  ########################################
  # Subsetting GE for two genes and clinical
  ########################################
  
  # Find rows of genes
  gene_row <- match(INPUT_GENES, INPUT_DF$GENE)
  gene_row <- na.omit(gene_row)
  # Reporting genes not found
  found_genes <- INPUT_DF[gene_row, "GENE"][[1]]
  if (length(found_genes) == length(INPUT_GENES)) {
    print("All Input Genes Found")  
  } else if (length(found_genes) < length(INPUT_GENES)){
    INPUT_GENES[!(INPUT_GENES %in% found_genes)] %>% paste(., "not found", sep = " ") %>% print()
  }
  
  # creating df
  if (OUTPUT_CLIN == T) {
    PLOT_DF <- INPUT_DF[gene_row, ] %>% 
      dplyr::select(., !(GENE)) %>% 
      t() %>% 
      as.data.frame()
    colnames(PLOT_DF) <-found_genes
  }else if (OUTPUT_CLIN == F){
    PLOT_DF <- INPUT_DF[gene_row, ]
  }
  
  
  # Adding in clinical data if specified
  if (OUTPUT_CLIN == T) {
    # Adding pateint ID
    PLOT_DF$TITLE_ID <- rownames(PLOT_DF) %>% factor()
    # Combining Clinical DF with Plot DF
    PLOT_DF <- inner_join(PLOT_DF, CLIN_DF, by = "TITLE_ID")
  }
  
  # Calculating tertiles for gene of intrest 
  if(CALCULATE_QUANTILE != F){
    
    # Quantile function
    n.quant <- function(x){
      frac <- 1/x 
      seq.int(0, by = frac, length.out = x )
    }
    
    # Calculating Tertile
    INPUT_QUANT = 3
    PLOT_DF$QUERY_TERTILE <- ntile(PLOT_DF[, CALCULATE_QUANTILE], INPUT_QUANT) %>% factor(., levels = seq(1:INPUT_QUANT))
    
    # Calculating UL quartile
    INPUT_QUANT = 4
    PLOT_DF$QUERY_UL25 <- ntile(PLOT_DF[, CALCULATE_QUANTILE], INPUT_QUANT) %>% factor(., levels = seq(1:INPUT_QUANT))
    PLOT_DF$QUERY_UL25 <- ifelse(PLOT_DF$QUERY_UL25 == 2 | PLOT_DF$QUERY_UL25 == 3, NA, ifelse(PLOT_DF$QUERY_UL25 == 1, 1, ifelse(PLOT_DF$QUERY_UL25 == 4, 4, NA)))
    PLOT_DF$QUERY_UL25 <- factor(PLOT_DF$QUERY_UL25, levels = c(1, 4))
    
    # Calculating median
    median_val <- median(PLOT_DF[, CALCULATE_QUANTILE])
    PLOT_DF$QUERY_MEDIAN <- ifelse(PLOT_DF[, CALCULATE_QUANTILE] >= median_val, 2, 1) %>%  factor(., levels = c(1, 2))
    
    
  }
  
  return(PLOT_DF)
}