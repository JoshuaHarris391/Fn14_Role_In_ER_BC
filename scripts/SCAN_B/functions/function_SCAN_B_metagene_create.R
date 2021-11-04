######################################################
# Calculating Metagenes by SVD
######################################################
# Creating Subset Function 
SCAN_B_metagene_create <- function(METAGENE_DF, DIRECTION, NAME_PREFIX){
  
  # DIRECTION = "+"
  # METAGENE_DF = prolif_ge_df[, 2:ncol(prolif_ge_df)]
  # NAME_PREFIX = "PROLIF"
  
  ######################################################
  # Creating metagene by SVD
  ######################################################
  
  # Converting DF to matrix 
  METAGENE_DF <- as.matrix(METAGENE_DF)
  
  # Scaling Matrix
  METAGENE_DF <- scale(t(METAGENE_DF)) # transposing genes to columns
  
  # # Checking means and SD
  # colMeans(METAGENE_DF)[1:5]
  # apply(METAGENE_DF, 2, sd)[1:5]
  
  # Re transposing matrix 
  METAGENE_DF <- t(METAGENE_DF)
  
  # Creating metagene 
  METAGENE_VEC <- svd(METAGENE_DF)$v[, 1] # do SVD, take 1st column of v matrix
  
  # Scaling Metagene
  METAGENE_VEC <- (METAGENE_VEC-min(METAGENE_VEC))/(max(METAGENE_VEC)-min(METAGENE_VEC))
  # Defining direction 
  ifelse(DIRECTION == "+", METAGENE_VEC <- 1-METAGENE_VEC, 
         ifelse(DIRECTION == "-", METAGENE_VEC <- METAGENE_VEC, print("Direction invalid")))
  
  # Creating tertile groups
  Meta_Tert <- quantile(METAGENE_VEC, prob = c(1/3, 2/3, 3/3)) 
  METAGENE_VEC_TERT <- ifelse(METAGENE_VEC <= Meta_Tert[3] & METAGENE_VEC > Meta_Tert[2], "High", 
                              ifelse(METAGENE_VEC <= Meta_Tert[2] & METAGENE_VEC > Meta_Tert[1], "Medium", 
                                     ifelse(METAGENE_VEC <= Meta_Tert[1], "Low", NA)))
  
  # # Use if you are correlating high gene expression with low metagene values
  # METAGENE_VEC_TERT <- ifelse(METAGENE_VEC <= Meta_Tert[1] , "High", 
  #                             ifelse(METAGENE_VEC <= Meta_Tert[2] & METAGENE_VEC > Meta_Tert[1], "Medium", 
  #                                    ifelse(METAGENE_VEC <= Meta_Tert[3] & METAGENE_VEC > Meta_Tert[2], "Low", NA)))
  
  
  METAGENE_VEC_TERT <- factor(METAGENE_VEC_TERT, c("Low", "Medium", "High"))
  # Creating Median Groups 
  METAGENE_VEC_MED <- ifelse(METAGENE_VEC <= median(METAGENE_VEC), "Low", "High")
  METAGENE_VEC_MED <- factor(METAGENE_VEC_MED, c("Low", "High"))
  
  
  # Creating Rank DF 
  METAGENE_SVD_DF <- data.frame(colnames(METAGENE_DF), METAGENE_VEC, rank(METAGENE_VEC), METAGENE_VEC_TERT, METAGENE_VEC_MED)
  # Creating colnames
  metagene_colnames <- c("TITLE_ID", "METAGENE_SCORE", "METAGENE_RANK", "METAGENE_TERTILE", "METAGENE_MEDIAN")
  metagene_colnames <- gsub("METAGENE", paste(NAME_PREFIX, "METAGENE", sep = "_"), metagene_colnames)
  colnames(METAGENE_SVD_DF) <- metagene_colnames
  
  # Final Dataframe 
  METAGENE_SVD_DF
}