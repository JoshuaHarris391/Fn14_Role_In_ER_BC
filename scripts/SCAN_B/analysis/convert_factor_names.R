##################################################################
# Converting some variables to publication ready factor names
##################################################################
# Converting Node status
GSE96058_Clinical_DF$LYMPH_NODE_STATUS <- gsub("NodeNegative", "Node Negative", GSE96058_Clinical_DF$LYMPH_NODE_STATUS)
GSE96058_Clinical_DF$LYMPH_NODE_STATUS <- gsub("NodePositive", "Node Positive", GSE96058_Clinical_DF$LYMPH_NODE_STATUS)
GSE96058_Clinical_DF$LYMPH_NODE_STATUS <- factor(GSE96058_Clinical_DF$LYMPH_NODE_STATUS, levels = c("Node Negative", "Node Positive"))

# Converting TNBC status
GSE96058_Clinical_DF$TNBC_STATUS <- gsub("Non_TNBC", "Non-TNBC", GSE96058_Clinical_DF$TNBC_STATUS)
GSE96058_Clinical_DF$TNBC_STATUS <- factor(GSE96058_Clinical_DF$TNBC_STATUS, levels = c("TNBC", "Non-TNBC"))


# Converting PAM50 subtype
GSE96058_Clinical_DF$PAM50_SUBTYPE <- gsub("LumA", "Luminal A", GSE96058_Clinical_DF$PAM50_SUBTYPE)
GSE96058_Clinical_DF$PAM50_SUBTYPE <- gsub("LumB", "Luminal B", GSE96058_Clinical_DF$PAM50_SUBTYPE)
GSE96058_Clinical_DF$PAM50_SUBTYPE <- gsub("Basal", "Basal-like", GSE96058_Clinical_DF$PAM50_SUBTYPE)
GSE96058_Clinical_DF$PAM50_SUBTYPE <- gsub("Her2", "HER2-enriched", GSE96058_Clinical_DF$PAM50_SUBTYPE)
GSE96058_Clinical_DF$PAM50_SUBTYPE <- gsub("Normal", "Normal-like", GSE96058_Clinical_DF$PAM50_SUBTYPE)
GSE96058_Clinical_DF$PAM50_SUBTYPE <- factor(GSE96058_Clinical_DF$PAM50_SUBTYPE, levels = c("Normal-like", "Basal-like", "HER2-enriched", "Luminal A", "Luminal B"))

# Creating Discrete Age Variable 
GSE96058_Clinical_DF$AGE_GROUP <- ifelse(GSE96058_Clinical_DF$AGE_AT_DIAGNOSIS > 50, ">50", 
                                  ifelse(GSE96058_Clinical_DF$AGE_AT_DIAGNOSIS >= 41 & GSE96058_Clinical_DF$AGE_AT_DIAGNOSIS <50, "41-50", 
                                         ifelse(GSE96058_Clinical_DF$AGE_AT_DIAGNOSIS <= 40, "<=40", NA)))
GSE96058_Clinical_DF$AGE_GROUP <- factor(GSE96058_Clinical_DF$AGE_GROUP, levels = c("<=40", "41-50", ">50"))
  
  
# Creating Discreate Tumour Sizes 
GSE96058_Clinical_DF$SIZE_GROUP <- ifelse(GSE96058_Clinical_DF$TUMOR_SIZE > 5, "T3", 
                                   ifelse(GSE96058_Clinical_DF$TUMOR_SIZE > 2 & GSE96058_Clinical_DF$TUMOR_SIZE <= 5, "T2", 
                                          ifelse(GSE96058_Clinical_DF$TUMOR_SIZE <= 2, "T1", NA)))
GSE96058_Clinical_DF$SIZE_GROUP <- factor(GSE96058_Clinical_DF$SIZE_GROUP, c("T1", "T2", "T3"))
  
  