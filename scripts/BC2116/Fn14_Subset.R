################################################################
## Subsetting
################################################################

# setwd("/Users/joshua_harris/Dropbox/Research/PhD/Projects/Clinical_Signficance_of_Fn14/Bioinformatics/TNFRSF12A_Clin_Sig_Manuscript_Analysis")

# Loading environment
# load(file = "Data/BC2116_DATA.RData")

# Variables:
GENE_INPUT <- "TNFRSF12A"
LOG_BASE <- 2

## Extracting Gene/genes 
# Finding Probes corresponding to Genes 
match_probes <- match(GENE_INPUT, ANOT_DF$GENES)
GENE_PROBES <- ANOT_DF[match_probes, 2]
# Extracting Probes from Expression DF 
GENE_EXP <- EXP[GENE_PROBES, ]
# Calculating Row Means 
GENE_EXP$ROW_MEANS <- rowMeans(GENE_EXP)
# Annotating Probenames 
match_probes <- match(GENE_PROBES, ANOT_DF$PROBES)
GENE_EXP$GENE_NAMES <- ANOT_DF[match_probes, 1]
# Reordering Cols 
GENE_EXP <- GENE_EXP[, c(2118,2117,1:2116)]
# Removing Duplicate Probes by median 
GENE_EXP <- aggregate(. ~ GENE_NAMES, GENE_EXP, FUN=function(x) if(length(x)==2) max(x) else median(x))
# Transposing DF 
GENE_EXP <- as.data.frame(t(GENE_EXP[, 3:2118]))
# Renaming Col 
colnames(GENE_EXP) <- "GENE_INPUT"
# Creating ID_REF 
GENE_EXP$ID_REF <- rownames(GENE_EXP)
# Reordering Cols 
GENE_EXP <- GENE_EXP[, c("ID_REF", "GENE_INPUT")]

## Adding Clinical Data 
# Unfactoring Clin ID_REF 
CLIN$ID_REF <- as.character(CLIN$ID_REF)
# Joining Expression DF and Clinical DF 
GENE_EXP <- dplyr::full_join(GENE_EXP, CLIN, by = "ID_REF")
# refactoring Subtype so that "Normal" is the baseline
GENE_EXP$SUBTYPE <- factor(GENE_EXP$SUBTYPE, c("Normal-like", "Basal-like", "HER2-enriched", "Luminal A", "Luminal B"))
# Unlogging Gene expression values to get multiplicative scale
GENE_EXP[, "GENE_INPUT"] <- 2^(GENE_EXP[, "GENE_INPUT"])

# Re transforming data
log_base <- LOG_BASE
GENE_EXP[, "GENE_INPUT"] <- log(GENE_EXP[, "GENE_INPUT"], base = log_base)

# Renaming "Gene Input"
colnames(GENE_EXP)[2] <- "Fn14"

################################################################
## Subsetting by ER, PR, or HR
################################################################
GENE_EXP_ER_P <- GENE_EXP %>% filter(., ER_STATUS == 1)
GENE_EXP_PR_P <- GENE_EXP %>% filter(., PR_STATUS == 1)
GENE_EXP_HR_P <- GENE_EXP %>% filter(., ER_STATUS == 1 | PR_STATUS == 1)
GENE_EXP_HR_P_Endo_Mono_T <- GENE_EXP_HR_P %>% filter(., ENDOCRINE_TREATMENT == 1 & CHEMOTHERAPY_TREATMENT == 0)
GENE_EXP_HR_P_Chemo_Endo <- GENE_EXP_HR_P %>% filter(., ENDOCRINE_TREATMENT == 1 & CHEMOTHERAPY_TREATMENT == 1)

################################################################
## Creating Fn14 Semicontinuous Variable 
################################################################
# Defining function to calculate tertiles
jhf_calculate_tertiles <- function(GENE_EXP){
  # Creating Fn14 Tertiles for GENE_EXP
  Fn14_Tert <- quantile(GENE_EXP$Fn14, prob = c(1/3, 2/3, 3/3)) 
  GENE_EXP$Fn14_TERT <- ifelse(GENE_EXP$Fn14 <= Fn14_Tert[3] & GENE_EXP$Fn14 > Fn14_Tert[2], 3, 
                               ifelse(GENE_EXP$Fn14 <= Fn14_Tert[2] & GENE_EXP$Fn14 > Fn14_Tert[1], 2, 
                                      ifelse(GENE_EXP$Fn14 <= Fn14_Tert[1], 1, NA)))
  table(GENE_EXP$Fn14_TERT)
  # Creating Fn14 tertile groups
  GENE_EXP$Fn14_TERT_GROUP <- ifelse(GENE_EXP$Fn14 <= Fn14_Tert[3] & GENE_EXP$Fn14 > Fn14_Tert[2], "High", 
                                     ifelse(GENE_EXP$Fn14 <= Fn14_Tert[2] & GENE_EXP$Fn14 > Fn14_Tert[1], "Medium", 
                                            ifelse(GENE_EXP$Fn14 <= Fn14_Tert[1], "Low", NA)))
  GENE_EXP$Fn14_TERT_GROUP <- factor(GENE_EXP$Fn14_TERT_GROUP, c("Low", "Medium", "High"))
  table(GENE_EXP$Fn14_TERT_GROUP) %>% print()
  return(GENE_EXP)
}

# Adding Fn14 Tertile info
GENE_EXP <- jhf_calculate_tertiles(GENE_EXP)
GENE_EXP_ER_P <- jhf_calculate_tertiles(GENE_EXP_ER_P)
GENE_EXP_PR_P <- jhf_calculate_tertiles(GENE_EXP_PR_P)
GENE_EXP_HR_P <- jhf_calculate_tertiles(GENE_EXP_HR_P)
GENE_EXP_HR_P_Endo_Mono_T <- jhf_calculate_tertiles(GENE_EXP_HR_P_Endo_Mono_T)
GENE_EXP_HR_P_Chemo_Endo <- jhf_calculate_tertiles(GENE_EXP_HR_P_Chemo_Endo)



