#################################################################
# Transfer SCAN-B Data from source repository
#################################################################
# # Download from dropbox here: https://www.dropbox.com/sh/2uwqu2f3q2wwvu8/AABIiK9XxPnvxrDYUPxzpwrYa?dl=0
# # Copying data
# dir.create("data", recursive = T)
# system("cp -R /Users/joshua_harris/Dropbox/Research/PhD/Bioinformatics/Datasets/Dataset_SCAN_B/SCAN_B_Dataset_Cleaned/data/GSE96058/clean/GSE96058_Cleaned.RData ./data/")


#################################################################
# Clear Environment and Load Cleaned Data
#################################################################
rm(list = ls())
library(tidyverse)
load(file = "data/GSE96058_Cleaned.RData")

#################################################################
# Converting variable factors to publication ready labels
#################################################################
source("./scripts/SCAN_B/analysis/convert_factor_names.R")


#################################################################
# Defining input gene
#################################################################
INPUT_GENE = "TNFRSF12A"
QUERY_GENES = c(INPUT_GENE, "ESR1", "PGR", "MKI67", "ZEB1", "PCNA", "MCM2")



#################################################################
# Loading functions
#################################################################
functions_list <- c("./scripts/SCAN_B/functions/function_SCAN_B_boxplot.R",
                    "./scripts/SCAN_B/functions/function_SCAN_B_subset.R")
for (i in functions_list) {
  source(i)
}

#################################################################
# Creating a summary DF to create figures from (all functions query this df)
#################################################################
SCAN_B_QUERY_DF <- SCAN_B_subset(GSE96058_GE, GSE96058_Clinical_DF, T, QUERY_GENES, CALCULATE_QUANTILE = "TNFRSF12A")
# Converting endocrine and chemo to integer
SCAN_B_QUERY_DF$ENDOCRINE_TREATED <- SCAN_B_QUERY_DF$ENDOCRINE_TREATED %>% as.integer() %>% factor()
SCAN_B_QUERY_DF$CHEMO_TREATED <- SCAN_B_QUERY_DF$CHEMO_TREATED %>% as.integer() %>% factor()

#################################################################
# Adding HR/HER2 subtype
#################################################################
source("./scripts/SCAN_B/analysis/SCAN_B_add_HR_HER2_Subtype.R")


#################################################################
# Creating Clinical var label lookup
#################################################################
# Creating proliferation metagene
source("./scripts/SCAN_B/analysis/clin_var_lookup.R")


#################################################################
# creating clinical summary tables
#################################################################
# # Creating proliferation metagene
# source("./scripts/SCAN_B/analysis/SCAN_B_Clinical_Summary_Tables.R")


#################################################################
# Create violin plots of clinical factors
#################################################################
source("./scripts/SCAN_B/analysis/SCAN_B_ER_Pos_create_violinplots.R")
source("./scripts/SCAN_B/analysis/SCAN_B_PGR_Pos_create_violinplots.R")
source("./scripts/SCAN_B/analysis/SCAN_B_HR_Pos_create_violinplots.R")


#################################################################
# Creating scatter plot of Fn14 expression by ESR1 expression in HR+/HER2- who received endocrine therapy
#################################################################
source("./scripts/SCAN_B/functions/function_SCAN_B_scatter_plot.R")
# Subsetting DF
HR_Pos_HER2_Neg_EndoMonoT <- SCAN_B_QUERY_DF %>% filter(., IHC_SUBTYPE == "HR-Pos_HER2-Neg", ENDOCRINE_TREATED == 1, CHEMO_TREATED == 0)
# Creating scatter plot
HR_Pos_HER2_Neg_EndoMonoT_Scatter <- SCAN_B_scatter_plot(HR_Pos_HER2_Neg_EndoMonoT, "ESR1", "TNFRSF12A")
# Saving plot
source(file = "scripts/functions/function_save_jpeg.R")
dir.create("outputs/SCAN_B/Scatter_plots", recursive = T, showWarnings = F)
jhf_save_jpeg(HR_Pos_HER2_Neg_EndoMonoT_Scatter, 10, 10, "HR_Pos_HER2_Neg_EndoMonoT_Scatter", "./outputs/SCAN_B/Scatter_plots/")


#################################################################
# Graphing Fn14 tertiles by ER IHC status
#################################################################
Fn14_Tert_by_ER <- SCAN_B_QUERY_DF %>% group_by(ER_STATUS, QUERY_TERTILE) %>% summarise(n = n())
Fn14_Tert_by_ER <- Fn14_Tert_by_ER %>% na.omit()
Fn14_Tert_by_ER_Plot <- ggplot(Fn14_Tert_by_ER, aes(fill = QUERY_TERTILE, y=n, x=ER_STATUS))+
                              geom_bar(position="dodge", stat="identity") +
                              scale_fill_discrete(name = "Fn14 mRNA Tertile") +
                              labs(x = "ER IHC Status", y = "Cases")
# Saving plot
dir.create("outputs/SCAN_B/Barplots", recursive = T, showWarnings = F)
jhf_save_jpeg(Fn14_Tert_by_ER_Plot, 10, 10, "Fn14_Tert_by_ER_Plot", "./outputs/SCAN_B/Barplots/")
# Making tables
table(SCAN_B_QUERY_DF$ER_STATUS, SCAN_B_QUERY_DF$QUERY_TERTILE) 
table(SCAN_B_QUERY_DF$ER_STATUS, SCAN_B_QUERY_DF$QUERY_TERTILE) %>% chisq.test()
table(SCAN_B_QUERY_DF$ER_STATUS, SCAN_B_QUERY_DF$QUERY_TERTILE) %>% prop.table() %>% round(., 3)


# Making tables for HR+/HER2-
Fn14_Tert_by_HR <- SCAN_B_QUERY_DF %>% group_by(IHC_SUBTYPE, QUERY_TERTILE) %>% summarise(n = n())
Fn14_Tert_by_HR <- Fn14_Tert_by_HR %>% na.omit()
Fn14_Tert_by_HR_Plot <- ggplot(Fn14_Tert_by_HR, aes(fill = QUERY_TERTILE, y=n, x=IHC_SUBTYPE))+
  geom_bar(position="dodge", stat="identity") +
  scale_fill_discrete(name = "Fn14 mRNA Tertile") +
  labs(x = "IHC Subtype", y = "Cases")

# Saving plot
jhf_save_jpeg(Fn14_Tert_by_HR_Plot, 20, 15, "Fn14_Tert_by_HR_Plot", "./outputs/SCAN_B/Barplots/")

table(SCAN_B_QUERY_DF$IHC_SUBTYPE, SCAN_B_QUERY_DF$QUERY_TERTILE) 
table(SCAN_B_QUERY_DF$IHC_SUBTYPE, SCAN_B_QUERY_DF$QUERY_TERTILE) %>% chisq.test()
table(SCAN_B_QUERY_DF$IHC_SUBTYPE, SCAN_B_QUERY_DF$QUERY_TERTILE) %>% prop.table() %>% round(., 3)

# Saving session
writeLines(capture.output(sessionInfo()), "SCAN_B_sessionInfo.txt")


