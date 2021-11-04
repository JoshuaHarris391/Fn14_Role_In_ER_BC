# Function takes in a vector, generates table of values, calculates percentage for each factor, outputs a df with Clinvar, factor, number, percentage, total

jhf_clin_summary_table <- function(INPUT_VEC_NAME, DF){
  
  # INPUT_VEC_NAME = "ER_Status"
  # DF = TCGA_DF
  
  # Subsetting 
  INPUT_VEC = DF[, INPUT_VEC_NAME]
  
  # getting summary table
  sum_table <- INPUT_VEC %>% table()
  # Getting total cases
  total_cases <- sum(sum_table)
  # Getting factors
  factor_names <- names(sum_table)
  # Getting n_cases
  n_cases <- sum_table %>% as.numeric()
  # Getting percentages
  n_percent <- (n_cases/total_cases)*100 
  # Building results df
  Res_DF <- data.frame(VARIABLE = rep(INPUT_VEC_NAME, length(factor_names)), FACTOR = factor_names, CASES = n_cases, PERCENTAGE = n_percent, TOTAL_CASES = total_cases)
  # Returning df
  return(Res_DF)
}