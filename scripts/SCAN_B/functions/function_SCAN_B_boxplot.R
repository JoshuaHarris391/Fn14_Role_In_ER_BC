

# This function creates boxplots
SCAN_B_boxplot <- function(XVAR, YVAR, INPUT_DF, Y_LABEL = "FPKM Log2"){
  
  # # For troubleshooting
  # XVAR = Clin_Vars[2]
  # YVAR = INPUT_GENE
  # INPUT_DF = SCAN_B_QUERY_DF
  # Y_LABEL = "FPKM Log2"
  
  ################################################################
  ### Preprocessing and creating query DF
  ################################################################
  # INPUT_DF = GENE_EXP
  # XVAR = ""
  # # Loading environment
  # load(file = "Data/BC2116_DATA.RData")
  # Packages
  library(ggpubr)
  library(tidyverse)
  library(EnvStats)
  
  # Depends to Metagene Creation 
  GENE_EXP_Boxplot <- INPUT_DF
  
  ################################################################
  ### Creating Association plots
  ################################################################
  
  ### Subsetting Clinical factors
  # Viewing Clinical Factors
  # dplyr::glimpse(GENE_EXP_Boxplot)
  
  # choosing clinical factor to subset
  clin_query <- subset(GENE_EXP_Boxplot, select = c(XVAR, YVAR))
  table(clin_query[, XVAR])
  bfsub <- nrow(clin_query)
  # removing na values
  clin_query <- na.omit(clin_query)
  clin_query <- subset(clin_query, !(clin_query[, XVAR] == ""))
  aftsub <- nrow(clin_query)
  # Printing percentage of missing values excluded
  print(paste("Percentage of samples remaining after filtering   ", round((aftsub/bfsub)*100, digits = 1), "%", sep=" "))
  # Renaming column
  colnames(clin_query) <- c("XVAR", "YVAR")
  # Checking Class
  class(clin_query$XVAR)
  
  
  
  ################################################################
  ### Stats
  ################################################################
  
  ## Note: based on Wilcoxon rank sum test, which is a non parametric test equivalent to the two sampled t-test. This test is appropriate many of the sample groups do not have a normal distribution
  
  if(length(levels(clin_query$XVAR)) > 2){
    global_stat <- compare_means(YVAR ~ XVAR, method = "kruskal.test", data =clin_query, p.adjust.method = "BH")
    if (global_stat$p.signif == "ns"){
      compar_stat <- NULL
    } else {
      # conduct Wilcox multiple testing
      wilcox <- compare_means(YVAR ~ XVAR, method = "wilcox.test", data = clin_query, p.adjust.method = "BH")
      Sig_Compar <- wilcox[wilcox$p.signif != "ns", ]
      my_comparisons <- 0
      # Creating comparisons loop
      for (i in 1:nrow(Sig_Compar)) {
        my_comparisons[i] <- list(c(Sig_Compar$group1[i], Sig_Compar$group2[i]))
      }
      compar_stat <- stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", p.adjust.method = "BH")
      
    }
    
  } else if(length(levels(clin_query$XVAR)) <= 2){
      global_stat <- compare_means(YVAR ~ XVAR, method = "wilcox.test", data = clin_query, p.adjust.method = "BH")
      # compar_stat <- NULL
      
      # conduct Wilcox testing
      wilcox <- compare_means(YVAR ~ XVAR, method = "wilcox.test", data = clin_query, p.adjust.method = "BH")
      Sig_Compar <- wilcox[1, ]
      my_comparisons <- list()
      # Creating comparisons loop
      my_comparisons["Wilcox"] <- list(c(Sig_Compar$group1[1], Sig_Compar$group2[1]))
      compar_stat <- stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", p.adjust.method = "BH")
      
    }
  
  
  
  ################################################################
  ### Boxplots
  ################################################################
  
  # Creating boxplots
  g <- ggplot(clin_query, aes(x = XVAR, y= YVAR)) +
    geom_violin() +
    # geom_violin() +
    labs(title= paste(CLIN_VAR_NAMES[XVAR], "\n", "p<", format(global_stat$p.adj, digits = 2),  sep = "" ),
         x = NULL,
         y= paste(YVAR, " (", Y_LABEL, ")", sep = "") ) +
    # y = NULL) +
    compar_stat +
    stat_n_text(size = 4) + 
    stat_summary(fun=median, geom="point", size=4, color="grey") +
    theme(title = element_text(size = 12), 
          axis.text.x = element_text(size=12), 
          axis.text.y = element_text(size=12), 
          axis.title.y = element_text(size=12))
  
  
  # Printing table
  print(table(clin_query$XVAR))
  
  # Plotting Graph
  return(g)
  
}
