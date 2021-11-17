################################################################
## Final Fn14 Model Construction
################################################################

# Converting GRADE to semicont
GENE_EXP$GRADE <- ifelse(GENE_EXP$GRADE == "G1", 1, ifelse(GENE_EXP$GRADE == "G2", 2, ifelse(GENE_EXP$GRADE == "G3", 3, NA))) %>% as.integer()

# Changing levels for ihc subtype
GENE_EXP$IHC_SUBTYPE <- factor(GENE_EXP$IHC_SUBTYPE, levels = c("HR-Pos_HER2-Neg", "HR-Pos_HER2-Pos", "HR-Neg_HER2-Pos", "TNBC"))

# Refactoring pam50 subtype so that Luminal A is the reference 
GENE_EXP$SUBTYPE <- factor(GENE_EXP$SUBTYPE, levels = c("Luminal A", "Basal-like", "HER2-enriched", "Normal-like", "Luminal B"))

# Defineing Formula
attach(GENE_EXP)
Cox_f <- as.formula(Surv(DMFS_TIME, DMFS_EVENT) ~ 
                      Fn14_TERT + 
                      IHC_SUBTYPE +
                      LN_STATUS +
                      GRADE + 
                      SIZE_GROUP_SEMICONT +
                      AGE_GROUP_SEMICONT +
                      ENDOCRINE_TREATMENT +
                      CHEMOTHERAPY_TREATMENT
)
# Constructing final model
Fn14_Cox_Model <- coxph(Cox_f, data = GENE_EXP)
detach(GENE_EXP)

print(paste("Fn14_Cox_Model", "n =", Fn14_Cox_Model$n))

################################################################
## Final Fn14 + PAM50 Model Construction
################################################################

# Defineing Formula
attach(GENE_EXP)
Cox_f <- as.formula(Surv(DMFS_TIME, DMFS_EVENT) ~ 
                      Fn14_TERT + 
                      IHC_SUBTYPE +
                      LN_STATUS +
                      GRADE + 
                      SIZE_GROUP_SEMICONT +
                      AGE_GROUP_SEMICONT +
                      SUBTYPE +
                      ENDOCRINE_TREATMENT +
                      CHEMOTHERAPY_TREATMENT
)
# Constructing final model
Fn14_PAM50_Cox_Model <- coxph(Cox_f, data = GENE_EXP)
detach(GENE_EXP)

print(paste("Fn14_PAM50_Cox_Model", "n =", Fn14_PAM50_Cox_Model$n))









################################################################
## Printing Final Multivariate Model Summaries
################################################################

# print("Fn14_Cox_Model")
# summary(Fn14_Cox_Model)
# 
# print("Fn14_PAM50_Cox_Model")
# summary(Fn14_PAM50_Cox_Model)



################################################################
## Storing Models in List
################################################################
COX_MODELS_LIST <- list()
COX_MODELS_LIST[["Fn14_Cox_Model"]] <- Fn14_Cox_Model
COX_MODELS_LIST[["Fn14_PAM50_Cox_Model"]] <- Fn14_PAM50_Cox_Model
# Viewing List
summary(COX_MODELS_LIST)












