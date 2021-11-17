jhf_calculate_tertiles <- function(INPUT_VEC, RETURN = "TERT_GROUP"){
  # Creating Tertiles for INPUT_VEC
  GENE_TERT <- quantile(INPUT_VEC, prob = c(1/3, 2/3, 3/3)) 
  
  if(RETURN == "TERT"){
    output <- ifelse(INPUT_VEC <= GENE_TERT[3] & INPUT_VEC > GENE_TERT[2], 3, 
                     ifelse(INPUT_VEC <= GENE_TERT[2] & INPUT_VEC > GENE_TERT[1], 2, 
                            ifelse(INPUT_VEC <= GENE_TERT[1], 1, NA)))
    table(output) %>% print()
    return(output)
  }
  
  if(RETURN == "TERT_GROUP"){
    # Creating tertile groups
    output <- ifelse(INPUT_VEC <= GENE_TERT[3] & INPUT_VEC > GENE_TERT[2], "High", 
                     ifelse(INPUT_VEC <= GENE_TERT[2] & INPUT_VEC > GENE_TERT[1], "Medium", 
                            ifelse(INPUT_VEC <= GENE_TERT[1], "Low", NA)))
    output <- factor(output, c("Low", "Medium", "High"))
    table(output) %>% print()
    return(output)
  } else {
    print("RETURN Parameter unrecognised, use TERT or TERT_GROUP")
  }
  
}