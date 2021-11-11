##########################################
## Clinical Significance of Fn14 in Breast Cancer
##########################################

##########################################
## Clinical Significance of Fn14 in Breast Cancer
##########################################
# Deleting environment
rm(list=ls())
# Loading tidyverse
library(tidyverse)

##########################################
# Copying and loading data
##########################################
# system("rsync -av ~/Dropbox/Research/PhD/Bioinformatics/Datasets/Dataset_BC2116/Data/BC2116_DATA.RData ~/Dropbox/Research/PhD/Projects/Clinical_Signficance_of_Fn14/Bioinformatics/TNFRSF12A_Clin_Sig_Manuscript_Analysis/Data/BC2116/")

# Loading environment
load(file = "data/BC2116_DATA.RData")

##########################################
# Loading Packages
##########################################
# Loading Packages
library(survminer)
library(survival)
library(tidyverse)
library(ggpubr)
library(vtreat)
library(pROC)
library(ggstance)
library(RColorBrewer)
library(gplots)
library(EnvStats)


# Running script to subset data and calculate tertiles
source(file = "scripts/BC2116/Fn14_Subset.R")

# Running survival analysis
source(file = "scripts/BC2116/survival_analysis_DMFS.R")
source(file = "scripts/BC2116/survival_analysis_DFS.R")
source(file = "scripts/BC2116/survival_analysis_DSS.R")


# # Deleting Environment
# rm(list=ls())