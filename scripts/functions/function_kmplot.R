jhf_kmplot <- function(INPUT_DF,GROUP_VEC,COX_SURV,SURV, YLAB = "Proportion", XLAB = "Time", BINARY = F){
  
  # Variables
  # INPUT_DF = TAM_GENE_CLIN_DF
  # GROUP_VEC = TAM_GENE_CLIN_DF$TNFRSF12A_TERT_GROUP
  # COX_SURV <- Surv(t.rfs, e.rfs) ~ TNFRSF12A_TERT
  # SURV <- "Surv(t.rfs, e.rfs) ~ TNFRSF12A_TERT_GROUP"
  # KM_YLAB <- "DFS"
  # KM_XLAB <- "Time (months)"

  library(survival)
  library(survminer)
  
  ################################################################
  # Defining Color palette
  ################################################################
  # ## Orange, Dark Blue, Pink
  # Col_Pal <- c("#ff622d", "#002749", "#ff0081")
  ## Blue, Black, Red
  Col_Pal <- c("#3399FF", "#333333", "#FF4040")
  
  
  ################################################################
  ## Stratifying KMplots by IHC IHC_SUBTYPE for Fn14 Tertiles
  ################################################################
  
  if(BINARY == T){
    KM_TITLE <- NULL
    KM_LEDG <- c(paste("Negative", " (" , "N=", table(GROUP_VEC == 0)[2], ")", sep = ""), 
                 paste("Positive", " (" , "N=", table(GROUP_VEC == 1)[2], ")", sep = ""))
    
  }else {
    # Defining KM plot parameters
    KM_TITLE <- NULL
    KM_LEDG <- c(paste("Low", " (" , "N=", table(GROUP_VEC == "Low")[2], ")", sep = ""), 
                 paste("Medium", " (" , "N=", table(GROUP_VEC == "Medium")[2], ")", sep = ""), 
                 paste("High", " (" , "N=", table(GROUP_VEC == "High")[2], ")", sep = ""))
  }
  
  
  # Calculating univar cox
  Cox_Uni <- coxph(COX_SURV, data = INPUT_DF)
  # Prepare the columns
  beta <- coef(Cox_Uni)
  HR  <- exp(beta)
  se   <- sqrt(diag(Cox_Uni$var))
  p   <- 1 - pchisq((beta/se)^2, 1)
  p <- round(p, digits = 6)
  p_round <- ifelse(p < 0.001, "< 0.001", round(p, digits = 3))
  CI   <- exp(confint(Cox_Uni, level = 0.95))
  CI <-  round(CI, digits = 2)
  CI_ALL <- paste("(", round(CI[1], digits = 2), " to ", round(CI[2], digits = 2), ")", sep = "")
  HR_CI <-  paste(round(HR, digits = 2), " (", round(CI[1], digits = 2), " to ", round(CI[2], digits = 2), ")", sep = "")
  n <- Cox_Uni$n
  
  # Defining cox lab
  HR_Pval_lab <- paste("Univariate Cox", "\n","HR = ", HR_CI, "\n", "P = ", p_round, sep = "")
  
  SURV <<- SURV
  
  FIT = survfit(formula = as.formula(SURV), data = INPUT_DF)
  
  # Plot
  KMplot <- ggsurvplot(
    FIT,
    data = INPUT_DF,
    font.main = 20,
    title = KM_TITLE,
    pval = HR_Pval_lab,
    pval.size = 4.5,
    pval.coord = c(10 ,0.15),
    legend.title = "",
    legend.labs = KM_LEDG,
    legend = c(0.8,0.2),
    risk.table = F,
    risk.table.height = 0.3,
    tables.y.text = F, 
    risk.table.pos = "out",
    risk.table.col = "strata",
    break.time.by = 12*5,
    censor = T,
    palette = Col_Pal,
    # xlab = KM_XLAB,
    xlab = XLAB,
    ylab = YLAB
    
  ) 
  
  KMplot <- ggpar(KMplot, 
                  font.x = 15, 
                  font.xtickslab = 17,
                  font.y = 15, 
                  font.ytickslab = 17
  )
  
  KMplot$plot <- KMplot$plot + theme(legend.text = element_text(size = rel(0.8)))
  
  return(KMplot)
}


  

