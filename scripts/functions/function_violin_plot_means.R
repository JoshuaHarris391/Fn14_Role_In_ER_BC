library(tidyverse)
library(EnvStats)

jhf_clin_violin_plot_means <- function(INPUT_DF, XVAR, YVAR, XVAR_LAB, YVAR_LAB, PATH, LOG2_TRANSFROM = T){
  
  # # Variables
  # XVAR = VARS[2]
  # YVAR = "TNFRSF12A_CPM"
  # INPUT_DF = TCGA_DF
  # XVAR_LAB = VARS[2]
  # YVAR_LAB = "Fn14 Log(2) CPM"
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
  
  # Calculating means
  mean_df <- aggregate(YVAR~XVAR, data = INPUT_DF, mean)
  
  # Adding unlog value
  mean_df$unlog <- 2^(mean_df$YVAR)
  
  # Adding back col names
  colnames(mean_df) <- c(XVAR, YVAR, "unlog")
  
  # Saving table 
  write.csv(mean_df, paste(PATH, "/", XVAR, "_violin_means.csv", sep = ""))
  # print(getwd())
  

}
