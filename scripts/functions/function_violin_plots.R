library(tidyverse)
library(EnvStats)

jhf_clin_violin_plot <- function(INPUT_DF, XVAR, YVAR, XVAR_LAB, YVAR_LAB, LOG2_TRANSFROM = F){
  
  # # Variables
  # XVAR = VARS[1]
  # YVAR = "TNFRSF12A"
  # INPUT_DF = SCAN_B_QUERY_DF
  # XVAR_LAB = VARS[1]
  # YVAR_LAB = "Fn14 Log(2) FPKM"
  # LOG2_TRANSFROM = T
  
  # Log transforming
  if(LOG2_TRANSFROM == T){
    INPUT_DF[ ,YVAR] <- log2(INPUT_DF[ ,YVAR])
  }
  
  # Subsetting df
  INPUT_DF <- INPUT_DF %>% dplyr::select(., all_of(XVAR), all_of(YVAR))
  
  # Renaming columns
  colnames(INPUT_DF) <- c("XVAR", "YVAR")
  
  # Removing NA values
  INPUT_DF <- na.omit(INPUT_DF)
  
  # Calculating non-parametric stats
  if(length(levels(INPUT_DF$XVAR)) > 2){
    global_stat <- ggpubr::compare_means(YVAR ~ XVAR, method = "kruskal.test", data =INPUT_DF, p.adjust.method = "bonferroni")
    if (global_stat$p.signif == "ns"){
      compar_stat <- NULL
    } else {
      # conduct Wilcox multiple testing
      wilcox <- ggpubr::compare_means(YVAR ~ XVAR, method = "wilcox.test", data = INPUT_DF, p.adjust.method = "bonferroni")
      Sig_Compar <- wilcox[wilcox$p.signif != "ns", ]
      my_comparisons <- 0
      # Creating comparisons loop
      for (i in 1:nrow(Sig_Compar)) {
        my_comparisons[i] <- list(c(Sig_Compar$group1[i], Sig_Compar$group2[i]))
      }
      compar_stat <- ggpubr::stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", p.adjust.method = "bonferroni")
      
    }
    
  } else if (length(levels(INPUT_DF$XVAR)) <= 2){
      global_stat <- ggpubr::compare_means(YVAR ~ XVAR, method = "wilcox.test", data = INPUT_DF, p.adjust.method = "bonferroni")
      compar_stat <- NULL
  }
  
  
  # Creating boxplots
  g <- ggplot(INPUT_DF, aes(x = XVAR, y= YVAR)) +
    # geom_boxplot() +
    geom_violin() +
    labs(title= paste(XVAR_LAB," (", "p<", format(global_stat$p.adj, digits = 2), ")",  sep = "" ),
         x = NULL,
         y = paste(YVAR_LAB, sep = "") ) +
    # y = NULL) +
    compar_stat +
    stat_n_text(size = 3) +
    stat_summary(fun=median, geom="point", size=3, color="pink") +
    theme(title = element_text(size = 8))
  
  
  # Printing table
  print(table(INPUT_DF$XVAR))
  
  # Plotting Graph
  return(g)
  
}