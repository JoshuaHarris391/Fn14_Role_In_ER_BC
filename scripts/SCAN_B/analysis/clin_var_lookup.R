# Lookup table for variable names
CLIN_VAR_NAMES <- as.list(colnames(SCAN_B_QUERY_DF))
names(CLIN_VAR_NAMES) <- colnames(SCAN_B_QUERY_DF)

# Replacing variable names with final formatting
 CLIN_VAR_NAMES["ER_STATUS"] <- "ER"
 CLIN_VAR_NAMES["PGR_STATUS"] <- "PGR"
 CLIN_VAR_NAMES["HER2_STATUS"] <- "HER2"
 CLIN_VAR_NAMES["TNBC_STATUS"] <- "TNBC"
 CLIN_VAR_NAMES["NHG"] <- "Grade"
 CLIN_VAR_NAMES["KI67_STATUS"] <- "KI67"
 CLIN_VAR_NAMES["PAM50_SUBTYPE"] <- "PAM50 Subtype"
 CLIN_VAR_NAMES["PROLIF_METAGENE_SCORE"] <- "Proliferation Metagene Score"
 CLIN_VAR_NAMES["EMT_METAGENE_SCORE"] <- "EMT Metagene Score"
 CLIN_VAR_NAMES["LYMPH_NODE_STATUS"] <- "Lymph Node Status"