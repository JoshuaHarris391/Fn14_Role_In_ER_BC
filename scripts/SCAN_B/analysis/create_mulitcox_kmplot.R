# Multivariate KMplot cox by median for all
########################################
# Defining Variables
########################################

library('survminer')
library('survival')
# Input_DF
INPUT_DF <- SCAN_B_QUERY_DF
# Input Clinical DF
# CLIN_DF <- GSE96058_Clinical_DF
# Input gene
INPUT_GENE = "SCNN1A"
# Quantile Split
QUANTILE_TYPE = "UL25"
# Clinical Subset
CLIN_SUBSET_VAR = "ALL"
CLIN_SUBSET_FACTOR = NULL
# Legend position
LEG_POS <- c(0.24, 0.3)

# Hard defining quantiles
if (QUANTILE_TYPE == "TERTILE") {
  INPUT_QUANT = 3
} else if (QUANTILE_TYPE == "UL25"){
  INPUT_QUANT = 4
}



########################################
# Generating plot data frame
########################################
PLOT_DF = INPUT_DF

## Subsetting all variables
if (CLIN_SUBSET_VAR == "ALL") {
  PLOT_DF <- PLOT_DF
} else {
  # Defining essential variables
  essential_col <- c(paste(INPUT_GENE), "TITLE_ID", "OVERALL_SURVIVAL_MONTHS", "OVERALL_SURVIVAL_EVENT", "AGE_AT_DIAGNOSIS", "TUMOR_SIZE", "LYMPH_NODE_STATUS", "NHG")
  # Selecting essential variables and clinical variable
  PLOT_DF <- select(PLOT_DF, c(all_of(essential_col), all_of(CLIN_SUBSET_VAR)))
  # Subsetting clinical variable by factor
  PLOT_DF <- filter(PLOT_DF, PLOT_DF[, CLIN_SUBSET_VAR] == !!CLIN_SUBSET_FACTOR)
}

# Calculating Quantiles
n.quant <- function(x){
  frac <- 1/x 
  seq.int(0, by = frac, length.out = x )
}

# Classifying Quantiles
if (QUANTILE_TYPE == "MEDIAN") {
  median_val <- median(PLOT_DF[, INPUT_GENE])
  PLOT_DF$QUANTILE <- ifelse(PLOT_DF[, INPUT_GENE] >= median_val, 2, 1) %>%  factor(., levels = c(1, 2))
} else if (QUANTILE_TYPE == "UL25") {
  PLOT_DF$QUANTILE <- ntile(PLOT_DF[, INPUT_GENE], INPUT_QUANT) %>% factor(., levels = seq(1:INPUT_QUANT))
  # Subsetting upper and lower 25th percentile 
  PLOT_DF <- PLOT_DF %>% filter(QUANTILE == 1 | QUANTILE == 4)
  PLOT_DF$QUANTILE <- factor(PLOT_DF$QUANTILE, levels = c(1, 4))
} else if (QUANTILE_TYPE == "TERTILE"){
  PLOT_DF$QUANTILE <- ntile(PLOT_DF[, INPUT_GENE], INPUT_QUANT) %>% factor(., levels = seq(1:INPUT_QUANT))
}



########################################
# Creating KMplot
########################################
# Creating survival models

SURV <- Surv(PLOT_DF$OVERALL_SURVIVAL_MONTHS, PLOT_DF$OVERALL_SURVIVAL_EVENT) ~ PLOT_DF$QUANTILE
FIT <- survival::survfit(Surv(OVERALL_SURVIVAL_MONTHS, OVERALL_SURVIVAL_EVENT) ~ QUANTILE, data = PLOT_DF)

# Creating Univariate Cox Model 
SURV.int.uni <- Surv(PLOT_DF$OVERALL_SURVIVAL_MONTHS, PLOT_DF$OVERALL_SURVIVAL_EVENT) ~ 
  PLOT_DF[, INPUT_GENE]
cox.int.uni <- coxph(SURV.int.uni, data = PLOT_DF) ~ AGE_AT_DIAGNOSIS + TUMOR_SIZE + LYMPH_NODE_STATUS + NHG
# pulling cox stats
HR.int.uni <-summary(cox.int.uni)$coefficients[1, 2]
p.value.int.uni <- summary(cox.int.uni)$coefficients[1, 5]
Cox.graph.int.uni <- paste("Univariate Cox:  ", "HR Ratio", "=", format(HR.int.uni, digits = 4), "   ", "P", "<", format(p.value.int.uni, digits = 4)) 


# Conditional statements to set legend
if (QUANTILE_TYPE == "TERTILE") {
  # Defining Color palette
  ## Blue, Black, Red
  Col_Pal <- c("#3399FF", "#333333", "#FF4040")
  # Setting legend
  KM_LEDG <- c(paste("Low", " (" , "N=", table(PLOT_DF$QUANTILE == 1)[2], ")", sep = ""), 
               paste("Medium", " (" , "N=", table(PLOT_DF$QUANTILE == 2)[2], ")", sep = ""), 
               paste("High", " (" , "N=", table(PLOT_DF$QUANTILE == 3)[2], ")", sep = ""))
  
} else if (QUANTILE_TYPE == "MEDIAN" | QUANTILE_TYPE == "UL25"){
  # Defining Color palette
  ## Blue, Black, Red
  Col_Pal <- c("#3399FF", "#FF4040")
  # Setting legend
  KM_LEDG <- c(paste("Low", " (" , "N=", table(PLOT_DF$QUANTILE == 1)[2], ")", sep = ""), 
               paste("High", " (" , "N=", table(as.character(PLOT_DF$QUANTILE) > 1)[2], ")", sep = ""))
}


# Title 
if (CLIN_SUBSET_VAR == "ALL") {
  title_var <- paste(INPUT_GENE, sep = "")
} else {
  title_var <- paste(INPUT_GENE, " (", CLIN_SUBSET_VAR, " = ", CLIN_SUBSET_FACTOR, ")", sep = "")
}

# KM Plot
KMplot <- ggsurvplot(
  FIT,
  data = PLOT_DF,
  title = title_var,
  pval = F,
  pval.size = 5,
  pval.method = T,
  pval.method.coord = c(0.1,0.05),
  pval.coord = c(2,0.05, hjust = 0),
  legend.title = "",
  legend.labs = KM_LEDG,
  legend = c(0.22,0.22),
  risk.table = F,
  risk.table.height = 0.30,
  risk.table.col = "strata",
  censor = T,
  xlab = "Time (Months)",
  ylab = "Overall Survival", 
  break.x.by = 6,
  conf.int = F,
  palette = Col_Pal
) + theme_survminer(font.legend = c(14, "plain", "black"), legend = LEG_POS, font.main = 18)


# Adding Cox Stats
KMplot <- KMplot$plot+
  ggplot2::annotate("text", x = 0, y = 0.05, size = 5, label = paste(Cox.graph.int.uni), hjust = 0)