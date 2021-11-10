# Making dir
dir.create("outputs/SCAN_B/ER_Pos",  recursive = T, showWarnings = F)
# loading function
source("scripts/functions/function_violin_plots.R")
# Subsetting query DF to ER+
SCAN_B_QUERY_DF_ER_POS <- SCAN_B_QUERY_DF %>% filter(., ER_STATUS == "Positive")


# Initialising list 
TCGA_VIOLIN_PLOTS_LIST <- list()

# Variables
VARS <- c("LYMPH_NODE_STATUS",
          "AGE_GROUP",
          "SIZE_GROUP",
          # "ER_STATUS",
          "PGR_STATUS",
          "HER2_STATUS",
          "KI67_STATUS",
          "NHG", 
          # "TNBC_STATUS", 
          "PAM50_SUBTYPE",
          "IHC_SUBTYPE",
          "ENDOCRINE_TREATED",
          "CHEMO_TREATED")

# Creating plots
for (i in VARS) {
  print(paste0(i))
  TCGA_VIOLIN_PLOTS_LIST[[paste0(i)]] <- jhf_clin_violin_plot(SCAN_B_QUERY_DF_ER_POS, i, "TNFRSF12A", i, "Fn14 Log(2) FPKM", LOG2_TRANSFROM = F)
}

# Saving small plots
source(file = "scripts/functions/function_save_jpeg.R")
VARS <- c("LYMPH_NODE_STATUS",
          # "ER_STATUS",
          "PGR_STATUS",
          "SIZE_GROUP",
          "HER2_STATUS",
          "KI67_STATUS",
          "AGE_GROUP",
          "NHG",
          "ENDOCRINE_TREATED",
          "CHEMO_TREATED")

for (i in VARS) {
  jhf_save_jpeg(TCGA_VIOLIN_PLOTS_LIST[[paste0(i)]], 8, 10, paste("SCAN_B", i, sep = "_"), "./outputs/SCAN_B/ER_Pos/")
}

# Saving large plots
VARS <- c("PAM50_SUBTYPE", "IHC_SUBTYPE")
for (i in VARS) {
  jhf_save_jpeg(TCGA_VIOLIN_PLOTS_LIST[[paste0(i)]], 16, 12, paste("SCAN_B", i, sep = "_"), "./outputs/SCAN_B/ER_Pos/")
}

# Saving violin plot means
VARS <- c("LYMPH_NODE_STATUS",
          "AGE_GROUP",
          "SIZE_GROUP",
          # "ER_STATUS",
          "PGR_STATUS",
          "HER2_STATUS",
          "KI67_STATUS",
          "NHG", 
          # "TNBC_STATUS", 
          "PAM50_SUBTYPE",
          "IHC_SUBTYPE",
          "ENDOCRINE_TREATED",
          "CHEMO_TREATED")
source("scripts/functions/function_violin_plot_means.R")
dir.create(path = "outputs/SCAN_B/violin_means/ER_Pos/", recursive = T, showWarnings = F)
for (i in VARS) {
  jhf_clin_violin_plot_means(SCAN_B_QUERY_DF_ER_POS, i, "TNFRSF12A", i, "TNFRSF12A (Log2 CPM)", PATH = "./outputs/SCAN_B/violin_means/ER_Pos/", LOG2_TRANSFROM = F)
}


