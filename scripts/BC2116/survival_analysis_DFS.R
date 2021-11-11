################################################################
## Subsetting Data
################################################################

# Creating subset DF lists
KM_Plot_List <- list()
KM_Plot_List[["ER_Pos"]] <- GENE_EXP_ER_P %>% select(DFS_TIME, DFS_EVENT, Fn14_TERT_GROUP, Fn14_TERT)
KM_Plot_List[["PR_Pos"]] <- GENE_EXP_PR_P %>% select(DFS_TIME, DFS_EVENT, Fn14_TERT_GROUP, Fn14_TERT)
KM_Plot_List[["HR_Pos"]]<- GENE_EXP_HR_P %>% select(DFS_TIME, DFS_EVENT, Fn14_TERT_GROUP, Fn14_TERT)
KM_Plot_List[["HR_Pos_Endo_Mono_T"]]<- GENE_EXP_HR_P_Endo_Mono_T %>% select(DFS_TIME, DFS_EVENT, Fn14_TERT_GROUP, Fn14_TERT)
KM_Plot_List[["HR_Pos_Chemo_Endo"]]<- GENE_EXP_HR_P_Chemo_Endo %>% select(DFS_TIME, DFS_EVENT, Fn14_TERT_GROUP, Fn14_TERT)

# Inspecting List
summary(KM_Plot_List)

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

# Initialising KMplot list
IHC_KMplots_List <- list()

# Running loop
for (i in names(KM_Plot_List)) {
  # Loading Packages
  library(survival)
  library(survminer)
  
  # Storing Iteration in INPUT_DF
  INPUT_DF = KM_Plot_List[[i]]
  
  # Removing NA values 
  INPUT_DF <- na.omit(INPUT_DF)
  
  # Defining KM plot parameters
  KM_TITLE <- NULL
  KM_LEDG <- c(paste("Low", " (" , "N=", table(INPUT_DF$Fn14_TERT_GROUP == "Low")[2], ")", sep = ""), 
               paste("Medium", " (" , "N=", table(INPUT_DF$Fn14_TERT_GROUP == "Medium")[2], ")", sep = ""), 
               paste("High", " (" , "N=", table(INPUT_DF$Fn14_TERT_GROUP == "High")[2], ")", sep = "")
  )
  KM_YLAB <- "DFS"
  KM_XLAB <- "Time (Years)"
  
  attach(INPUT_DF)
  
  # Creating Survival Objects
  SURV <- Surv(DFS_TIME, DFS_EVENT) ~ Fn14_TERT_GROUP
  FIT <- survival::survfit(Surv(DFS_TIME, DFS_EVENT) ~ Fn14_TERT_GROUP, data = INPUT_DF)
  
  # # Saving 10 year survival 
  # year_time <- 10
  # tenyearsurv <- summary(FIT, times = year_time) %>% capture.output()
  # dir.create(path = "outputs/BC2116/KMplots/DFS/ten_year_surv_summary", recursive = T, showWarnings = F)
  # survfilepath <- paste("outputs/BC2116/KMplots/DFS/ten_year_surv_summary/", i, "_ten_year_surv.txt", sep = "")
  # write_lines(tenyearsurv, file = survfilepath)
  
  # Calculating univar cox
  Cox_Uni <- coxph(Surv(INPUT_DF$DFS_TIME, INPUT_DF$DFS_EVENT) ~ INPUT_DF[, "Fn14_TERT"], data = INPUT_DF)
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
  HR_Pval_lab <- paste("HR = ", HR_CI, "\n", "                       P = ", p_round, sep = "")
  
  # Plot
  KMplot <- ggsurvplot(
    FIT,
    data = INPUT_DF,
    font.main = 20,
    title = KM_TITLE,
    pval = HR_Pval_lab,
    pval.size = 3.9,
    pval.coord = c(max(DFS_TIME)*0.48 ,0.95),
    legend.title = "",
    legend.labs = KM_LEDG,
    legend = c(0.23,0.2),
    risk.table = F,
    risk.table.height = 0.3,
    tables.y.text = F, 
    risk.table.pos = "out",
    risk.table.col = "strata",
    break.time.by = 5,
    censor = T,
    palette = Col_Pal,
    # xlab = KM_XLAB,
    xlab = NULL,
    ylab = KM_YLAB
    
  ) 
  KMplot
  
  KMplot <- ggpar(KMplot, 
                  font.x = 15, 
                  font.xtickslab = 17,
                  font.y = 15, 
                  font.ytickslab = 17
  )
  
  KMplot$plot <- KMplot$plot + theme(legend.text = element_text(size = rel(0.8)))
  
  # Log Rank Test 
  Deg_Free <- 2
  STATS <- survdiff(SURV)
  Chi <- as.numeric(STATS[5])
  Chistat <- pchisq(Chi,Deg_Free,lower.tail=FALSE) 
  KMplot_LogRank <- paste("Log Rank P Value:", "  p < ", format(Chistat, digits = 3), sep = "")
  # print(KMplot_LogRank)
  
  # Creating list of plots 
  IHC_KMplots_List[[paste0(i)]] <- KMplot
  
  detach(INPUT_DF)
}


# Checking KMplot list 
summary(IHC_KMplots_List)


##########################################
# Saving IHC plots as square grid
##########################################


##########################################
# Saving IHC plots as individual plots
##########################################
# Variables
dir.create("outputs/BC2116/KMplots/DFS/", recursive = T, showWarnings = F)
File.Path <- "outputs/BC2116/KMplots/DFS/"
# Loop.Variable <- "Vector of Variables for the loop"
temp_width <- 12
temp_height <- 10
# Saving kmplots
for (i in names(IHC_KMplots_List)) {
  Plot.object <-  IHC_KMplots_List[[i]]
  F_Name <- paste(File.Path, "Fn14_tert_kmplot_by_", i, ".jpeg", sep = "")
  jpeg(F_Name, width = temp_width, height = temp_height, units = "cm", res = 500)
  print(Plot.object)
  dev.off()
}

