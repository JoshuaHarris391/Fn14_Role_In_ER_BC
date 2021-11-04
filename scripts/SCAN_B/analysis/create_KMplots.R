########################################
# Environment Prep
########################################
# load(file = "data/GSE96058/clean/GSE96058_Cleaned.RData")
library('tidyverse')
library('survminer')
library('survival')

# # Query Gene (Now defined in MASTER in local env)
# INPUT_GENE <- "SCNN1A"

# Defining plot save resolution
RES_VAR = 500



################################################################################
# Making kmplots split by TERTILEs
################################################################################
########################################
# Creating kmplots and storing in lists
########################################
# initialising lists
list_TERTILE <- list()
# Creating km plots for PAM50 subtypes
subtypes_ref <- c("Normal-like", "Basal-like", "HER2-enriched", "Luminal A", "Luminal B")

# Defining legend size
LEG_SIZE <- c(0.24, 0.3)

# Creating km plots and storing in lists
list_TERTILE[["ER_NEG"]] <- SCAN_B_kmplot(SCAN_B_QUERY_DF, INPUT_GENE, "ER_STATUS", "Negative", LEG_SIZE, "TERTILE")
list_TERTILE[["ER_POS"]] <- SCAN_B_kmplot(SCAN_B_QUERY_DF, INPUT_GENE, "ER_STATUS", "Positive", LEG_SIZE, "TERTILE")

list_TERTILE[["PGR_NEG"]] <- SCAN_B_kmplot(SCAN_B_QUERY_DF, INPUT_GENE, "PGR_STATUS", "Negative", LEG_SIZE, "TERTILE")
list_TERTILE[["PGR_POS"]] <- SCAN_B_kmplot(SCAN_B_QUERY_DF, INPUT_GENE, "PGR_STATUS", "Positive", LEG_SIZE, "TERTILE")

list_TERTILE[["HER2_NEG"]] <- SCAN_B_kmplot(SCAN_B_QUERY_DF, INPUT_GENE, "HER2_STATUS", "Negative", LEG_SIZE, "TERTILE")
list_TERTILE[["HER2_POS"]] <- SCAN_B_kmplot(SCAN_B_QUERY_DF, INPUT_GENE, "HER2_STATUS", "Positive", LEG_SIZE, "TERTILE")

list_TERTILE[["TNBC"]] <- SCAN_B_kmplot(SCAN_B_QUERY_DF, INPUT_GENE, "TNBC_STATUS", "TNBC", LEG_SIZE, "TERTILE")
list_TERTILE[["NON_TNBC"]] <- SCAN_B_kmplot(SCAN_B_QUERY_DF, INPUT_GENE, "TNBC_STATUS", "Non-TNBC", LEG_SIZE, "TERTILE")

list_TERTILE[["ALL"]] <- SCAN_B_kmplot(SCAN_B_QUERY_DF, INPUT_GENE, "ALL", NULL, LEG_SIZE, "TERTILE")

list_TERTILE[["CHEMO_TREATED"]] <- SCAN_B_kmplot(SCAN_B_QUERY_DF, INPUT_GENE, "CHEMO_TREATED", "1", LEG_SIZE, "TERTILE")
list_TERTILE[["NON_CHEMO_TREATED"]] <- SCAN_B_kmplot(SCAN_B_QUERY_DF, INPUT_GENE, "CHEMO_TREATED", "0", LEG_SIZE, "TERTILE")

list_TERTILE[["ENDOCRINE_TREATED"]] <- SCAN_B_kmplot(SCAN_B_QUERY_DF, INPUT_GENE, "ENDOCRINE_TREATED", "1", LEG_SIZE, "TERTILE")
list_TERTILE[["NON_ENDOCRINE_TREATED"]] <- SCAN_B_kmplot(SCAN_B_QUERY_DF, INPUT_GENE, "ENDOCRINE_TREATED", "0", LEG_SIZE, "TERTILE")

# PAM50 subtypes
for (subtype_input in subtypes_ref) {
  list_TERTILE[[paste("PAM50", subtype_input, sep = "_")]] <- SCAN_B_kmplot(SCAN_B_QUERY_DF, INPUT_GENE, "PAM50_SUBTYPE", subtype_input, LEG_SIZE, "TERTILE")
}


###################################
# Saving KM plots
###################################
# Creating DIR for kmplots
dir.create(path = paste("./outputs/km_plots", INPUT_GENE, "TERTILE", sep = "/"), recursive = T, showWarnings = F)

# Running Loop to save plots
for (i in 1:length(list_TERTILE)) {
  # Get Name of Figure
  subtype_name <- list_TERTILE[i] %>% names()
  # Defining file path
  file_path <-  paste("./outputs/km_plots", INPUT_GENE, "TERTILE", sep = "/")
  # Creating file path name
  F_Name <- paste(file_path, "/", INPUT_GENE, "_", subtype_name, ".jpeg", sep = "")
  # Defining plot object and saving
  Plot.object <- list_TERTILE[[i]]
  jpeg(F_Name, width = 15, height = 15, units = "cm", res = RES_VAR)
  plot(Plot.object)
  dev.off()
}










################################################################################
# Making kmplots split by UL25s
################################################################################
########################################
# Creating kmplots and storing in lists
########################################
# initialising lists
list_UL25 <- list()
# Defining legend size
LEG_SIZE <- c(0.24, 0.3)

# Creating km plots and storing in lists
list_UL25[["ER_NEG"]] <- SCAN_B_kmplot(SCAN_B_QUERY_DF, INPUT_GENE, "ER_STATUS", "Negative", LEG_SIZE, "UL25")
list_UL25[["ER_POS"]] <- SCAN_B_kmplot(SCAN_B_QUERY_DF, INPUT_GENE, "ER_STATUS", "Positive", LEG_SIZE, "UL25")

list_UL25[["PGR_NEG"]] <- SCAN_B_kmplot(SCAN_B_QUERY_DF, INPUT_GENE, "PGR_STATUS", "Negative", LEG_SIZE, "UL25")
list_UL25[["PGR_POS"]] <- SCAN_B_kmplot(SCAN_B_QUERY_DF, INPUT_GENE, "PGR_STATUS", "Positive", LEG_SIZE, "UL25")

list_UL25[["HER2_NEG"]] <- SCAN_B_kmplot(SCAN_B_QUERY_DF, INPUT_GENE, "HER2_STATUS", "Negative", LEG_SIZE, "UL25")
list_UL25[["HER2_POS"]] <- SCAN_B_kmplot(SCAN_B_QUERY_DF, INPUT_GENE, "HER2_STATUS", "Positive", LEG_SIZE, "UL25")

list_UL25[["TNBC"]] <- SCAN_B_kmplot(SCAN_B_QUERY_DF, INPUT_GENE, "TNBC_STATUS", "TNBC", LEG_SIZE, "UL25")
list_UL25[["NON_TNBC"]] <- SCAN_B_kmplot(SCAN_B_QUERY_DF, INPUT_GENE, "TNBC_STATUS", "Non-TNBC", LEG_SIZE, "UL25")

list_UL25[["ALL"]] <- SCAN_B_kmplot(SCAN_B_QUERY_DF, INPUT_GENE, "ALL", NULL, LEG_SIZE, "UL25")

list_UL25[["CHEMO_TREATED"]] <- SCAN_B_kmplot(SCAN_B_QUERY_DF, INPUT_GENE, "CHEMO_TREATED", "1", LEG_SIZE, "UL25")
list_UL25[["NON_CHEMO_TREATED"]] <- SCAN_B_kmplot(SCAN_B_QUERY_DF, INPUT_GENE, "CHEMO_TREATED", "0", LEG_SIZE, "UL25")

list_UL25[["ENDOCRINE_TREATED"]] <- SCAN_B_kmplot(SCAN_B_QUERY_DF, INPUT_GENE, "ENDOCRINE_TREATED", "1", LEG_SIZE, "UL25")
list_UL25[["NON_ENDOCRINE_TREATED"]] <- SCAN_B_kmplot(SCAN_B_QUERY_DF, INPUT_GENE, "ENDOCRINE_TREATED", "0", LEG_SIZE, "UL25")

# PAM50 subtypes
for (subtype_input in subtypes_ref) {
  list_UL25[[paste("PAM50", subtype_input, sep = "_")]] <- SCAN_B_kmplot(SCAN_B_QUERY_DF, INPUT_GENE, "PAM50_SUBTYPE", subtype_input, LEG_SIZE, "UL25")
}


###################################
# Saving KM plots
###################################
# Creating DIR for kmplots
dir.create(path = paste("./outputs/km_plots", INPUT_GENE, "UL25", sep = "/"), recursive = T, showWarnings = F)

# Running Loop to save plots
for (i in 1:length(list_UL25)) {
  # Get Name of Figure
  subtype_name <- list_UL25[i] %>% names()
  # Defining file path
  file_path <-  paste("./outputs/km_plots", INPUT_GENE, "UL25", sep = "/")
  # Creating file path name
  F_Name <- paste(file_path, "/", INPUT_GENE, "_", subtype_name, ".jpeg", sep = "")
  # Defining plot object and saving
  Plot.object <- list_UL25[[i]]
  jpeg(F_Name, width = 15, height = 15, units = "cm", res = RES_VAR)
  plot(Plot.object)
  dev.off()
}










################################################################################
# Making kmplots split by MEDIANs
################################################################################
########################################
# Creating kmplots and storing in lists
########################################
# initialising lists
list_MEDIAN <- list()
# Defining legend size
LEG_SIZE <- c(0.24, 0.3)

# Creating km plots and storing in lists
list_MEDIAN[["ER_NEG"]] <- SCAN_B_kmplot(SCAN_B_QUERY_DF, INPUT_GENE, "ER_STATUS", "Negative", LEG_SIZE, "MEDIAN")
list_MEDIAN[["ER_POS"]] <- SCAN_B_kmplot(SCAN_B_QUERY_DF, INPUT_GENE, "ER_STATUS", "Positive", LEG_SIZE, "MEDIAN")

list_MEDIAN[["PGR_NEG"]] <- SCAN_B_kmplot(SCAN_B_QUERY_DF, INPUT_GENE, "PGR_STATUS", "Negative", LEG_SIZE, "MEDIAN")
list_MEDIAN[["PGR_POS"]] <- SCAN_B_kmplot(SCAN_B_QUERY_DF, INPUT_GENE, "PGR_STATUS", "Positive", LEG_SIZE, "MEDIAN")

list_MEDIAN[["HER2_NEG"]] <- SCAN_B_kmplot(SCAN_B_QUERY_DF, INPUT_GENE, "HER2_STATUS", "Negative", LEG_SIZE, "MEDIAN")
list_MEDIAN[["HER2_POS"]] <- SCAN_B_kmplot(SCAN_B_QUERY_DF, INPUT_GENE, "HER2_STATUS", "Positive", LEG_SIZE, "MEDIAN")

list_MEDIAN[["TNBC"]] <- SCAN_B_kmplot(SCAN_B_QUERY_DF, INPUT_GENE, "TNBC_STATUS", "TNBC", LEG_SIZE, "MEDIAN")
list_MEDIAN[["NON_TNBC"]] <- SCAN_B_kmplot(SCAN_B_QUERY_DF, INPUT_GENE, "TNBC_STATUS", "Non-TNBC", LEG_SIZE, "MEDIAN")

list_MEDIAN[["ALL"]] <- SCAN_B_kmplot(SCAN_B_QUERY_DF, INPUT_GENE, "ALL", NULL, LEG_SIZE, "MEDIAN")

list_MEDIAN[["CHEMO_TREATED"]] <- SCAN_B_kmplot(SCAN_B_QUERY_DF, INPUT_GENE, "CHEMO_TREATED", "1", LEG_SIZE, "MEDIAN")
list_MEDIAN[["NON_CHEMO_TREATED"]] <- SCAN_B_kmplot(SCAN_B_QUERY_DF, INPUT_GENE, "CHEMO_TREATED", "0", LEG_SIZE, "MEDIAN")

list_MEDIAN[["ENDOCRINE_TREATED"]] <- SCAN_B_kmplot(SCAN_B_QUERY_DF, INPUT_GENE, "ENDOCRINE_TREATED", "1", LEG_SIZE, "MEDIAN")
list_MEDIAN[["NON_ENDOCRINE_TREATED"]] <- SCAN_B_kmplot(SCAN_B_QUERY_DF, INPUT_GENE, "ENDOCRINE_TREATED", "0", LEG_SIZE, "MEDIAN")

# PAM50 subtypes
for (subtype_input in subtypes_ref) {
  list_MEDIAN[[paste("PAM50", subtype_input, sep = "_")]] <- SCAN_B_kmplot(SCAN_B_QUERY_DF, INPUT_GENE, "PAM50_SUBTYPE", subtype_input, LEG_SIZE, "MEDIAN")
}


###################################
# Saving KM plots
###################################
# Creating DIR for kmplots
dir.create(path = paste("./outputs/km_plots", INPUT_GENE, "MEDIAN", sep = "/"), recursive = T, showWarnings = F)

# Running Loop to save plots
for (i in 1:length(list_MEDIAN)) {
  # Get Name of Figure
  subtype_name <- list_MEDIAN[i] %>% names()
  # Defining file path
  file_path <-  paste("./outputs/km_plots", INPUT_GENE, "MEDIAN", sep = "/")
  # Creating file path name
  F_Name <- paste(file_path, "/", INPUT_GENE, "_", subtype_name, ".jpeg", sep = "")
  # Defining plot object and saving
  Plot.object <- list_MEDIAN[[i]]
  jpeg(F_Name, width = 15, height = 15, units = "cm", res = RES_VAR)
  plot(Plot.object)
  dev.off()
}

