jhf_extract_cox <- function(COX_INPUT, SAVE_PATH = NULL, FILENAME = NULL, OUTPUT = "data.frame"){
  # Defining input
  Cox_Multi <- COX_INPUT
  # Prepare the columns
  beta <- coef(Cox_Multi)
  HR  <- exp(beta)
  se   <- sqrt(diag(Cox_Multi$var))
  p   <- 1 - pchisq((beta/se)^2, 1)
  p <- round(p, digits = 6)
  p_round <- ifelse(p < 0.001, "< 0.001", round(p, digits = 3))
  CI   <- confint(Cox_Multi, level = 0.95)
  CI <- exp(CI)
  CI <- round(CI, digits = 2)
  CI_ALL <- paste("(", round(CI[,1], digits = 2), " to ", round(CI[,2], digits = 2), ")", sep = "")
  HR_CI <-  paste(round(HR, digits = 2), "   (", round(CI[,1], digits = 2), " to ", round(CI[,2], digits = 2), ")", sep = "")
  # Bind columns together, and select desired rows
  Cox_Multi_DF <- data.frame(format(round(HR, 3), nsmall = 3), CI_ALL, CI, HR_CI, p, as.character(p_round))
  # Labeling Columns 
  colnames(Cox_Multi_DF) <- c("HR", "95% CI","Lower 95% CI", "Upper 95% CI", "Hazard Ratio (95% CI)", "pValue", "Rounded P-value")
  # Adding in significance Stars
  sig.func <- function(x){
    if(x < 0.05 & x >= 0.01){
      "*"
    } else if (x < 0.01 & x >= 0.001){
      "**"
    } else if (x < 0.001 & x >= 0.0001){
      "***"
    } else if (x < 0.0001){
      "****"
    } else if (x >= 0.05){
      "N.S"
    } else {
      "Non Numeric"
    }
  }
  Cox_Multi_DF$Sig <- lapply(Cox_Multi_DF$pValue, sig.func)
  Cox_Multi_DF$Sig <- as.character(Cox_Multi_DF$Sig)
  

  
  # Output or saving
  if(OUTPUT == "data.frame"){
    return(Cox_Multi_DF)
  } 
  
  if(OUTPUT == "save.file"){
    # saving Multivariate DFs
    write.table(Cox_Multi_DF, paste(SAVE_PATH , FILENAME, ".txt"), row.names = TRUE, sep = "\t")
    write.csv(Cox_Multi_DF, paste(SAVE_PATH , FILENAME, ".csv", sep = ""), row.names = TRUE)
    # Printing Model
    print(Cox_Multi_DF)
  } else {
    print("unknown output, use either data.frame or save.file")
  }
}









