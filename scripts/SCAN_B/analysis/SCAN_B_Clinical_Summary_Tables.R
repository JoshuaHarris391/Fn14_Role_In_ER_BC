# Summary table of TCGA clinical factors
glimpse(SCAN_B_QUERY_DF)

# Loading function
source(file = "scripts/functions/function_clin_summary_table.R")

# Listing clinical factors for summary 
summary_vars <- c("AGE_GROUP", "NHG","LYMPH_NODE_STATUS", "SIZE_GROUP", "IHC_SUBTYPE", "PAM50_SUBTYPE", "ENDOCRINE_TREATED", "CHEMO_TREATED")

# Initialising Clinical summary table
SCAN_B_Clinical_Summary_DF <- data.frame()

# Creating summaries
for (i in summary_vars) {
  tmp_df <- jhf_clin_summary_table(i, SCAN_B_QUERY_DF)
  
  if(nrow(SCAN_B_Clinical_Summary_DF) < 1){
    SCAN_B_Clinical_Summary_DF <- tmp_df
  } else{
    SCAN_B_Clinical_Summary_DF <- rbind(SCAN_B_Clinical_Summary_DF, tmp_df)
  }
}

# Saving summary table
dir.create(path = "outputs/clinical_summary_tables/", recursive = T, showWarnings = F)
write.csv(SCAN_B_Clinical_Summary_DF, file = "outputs/clinical_summary_tables/SCAN_B_Clinical_Summary_DF.csv")
