#################################################################
# Clear Environment and Load Cleaned Data
#################################################################
rm(list = ls())
library(tidyverse)
load(file = "data/GSE6532_LUMINAL.RData")

# Checking data
glimpse(annot.tam) # Probe annotations
# View(data.tam[1:10, 1:10]) # Gene expression values as log2, Patients rows, probe cols
glimpse(demo.tam) # clinical info

# Converting annotation matrix to dataframe
annot.tam <- annot.tam %>% as.data.frame()

####################################################
# Getting basic stats about data
####################################################
data.tam %>% na.omit() %>% median()
data.tam %>% na.omit() %>% sd()
data.tam %>% na.omit() %>% mean()
data.tam %>% na.omit() %>% range()
plot(density((na.omit(data.tam))))

####################################################
# Extracting genes of interest
####################################################
INPUT_GENES <- c("TNFRSF12A", "ESR1")
# Finding probe names
m <- match(INPUT_GENES, annot.tam$HUGO.gene.symbol)
INPUT_GENES_PROBES <- annot.tam[m, "id"]
# Constructing DF with gene expression and clinical data
TAM_GENE_CLIN_DF <- data.tam[ , INPUT_GENES_PROBES] %>% as.data.frame()
# Renaming probes
colnames(TAM_GENE_CLIN_DF) <- INPUT_GENES
# Adding samplename column
TAM_GENE_CLIN_DF$samplename <- rownames(TAM_GENE_CLIN_DF)
# Adding clinical info to GENE_CLIN_DF
TAM_GENE_CLIN_DF <- TAM_GENE_CLIN_DF %>% left_join(demo.tam, by = "samplename")
# Converting survival time from days to months
TAM_GENE_CLIN_DF$t.rfs <- round(TAM_GENE_CLIN_DF$t.rfs/30.417, digit=0)


####################################################
# Calculating PAM50 Subtype 
####################################################
library(genefu)
# Editing annotation dfs
colnames(annot.tam)[1] <- "probe"
# Removing spaces and converting to integer
annot.tam$EntrezGene.ID <- gsub(" ", "", annot.tam$EntrezGene.ID) 
annot.tam$EntrezGene.ID <- gsub(" ", "", annot.tam$EntrezGene.ID) 
annot.tam$EntrezGene.ID <- annot.tam$EntrezGene.ID %>% as.integer()

# Calculating PAM50 using genefu
data(pam50.robust)
TAM_PAM50 <- molecular.subtyping(sbt.model = "pam50",data = data.tam, annot = annot.tam,do.mapping = TRUE)
# Adding PAM50 to GENE_CLIN_DF
TAM_GENE_CLIN_DF$PAM50 <- TAM_PAM50$subtype

####################################################
# Calculating expression tertiles
####################################################
source(file = "scripts/functions/function_calculate_tertiles.R")

TAM_GENE_CLIN_DF$TNFRSF12A_TERT <- jhf_calculate_tertiles(TAM_GENE_CLIN_DF$TNFRSF12A, RETURN = "TERT")
TAM_GENE_CLIN_DF$ESR1_TERT <- jhf_calculate_tertiles(TAM_GENE_CLIN_DF$ESR1, RETURN = "TERT")

TAM_GENE_CLIN_DF$TNFRSF12A_TERT_GROUP <- jhf_calculate_tertiles(TAM_GENE_CLIN_DF$TNFRSF12A, RETURN = "TERT_GROUP")
TAM_GENE_CLIN_DF$ESR1_TERT_GROUP <- jhf_calculate_tertiles(TAM_GENE_CLIN_DF$ESR1, RETURN = "TERT_GROUP")


####################################################
# Creating violin plots
####################################################
Fn14_Violin_Plot <- ggplot(TAM_GENE_CLIN_DF, aes(PAM50, TNFRSF12A)) +
  geom_violin() +
  stat_summary(fun=median, geom="point", size=2, color="grey") +
  labs(y = "TNFRSF12A (log2 Signal)")

prop_df <- TAM_GENE_CLIN_DF %>% 
  group_by(PAM50, TNFRSF12A_TERT_GROUP) %>% 
  summarise(n = n()) %>% 
  mutate(prop = n/sum(n))

Fn14_Tert_Prop_Subtype_Plot <- ggplot(prop_df, aes(PAM50, prop, fill = TNFRSF12A_TERT_GROUP)) + geom_bar(stat = "identity") 

# Saving plots
dir.create(path = "outputs/GSE6532/Expression_Plots/Tamoxifen/", recursive = T, showWarnings = F)
source(file = "scripts/functions/function_save_jpeg.R")
jhf_save_jpeg(Fn14_Violin_Plot, 15, 13, "Fn14_Log2_Violin_Plot_PAM50", "outputs/GSE6532/Expression_Plots/Tamoxifen/", RES = 500)
jhf_save_jpeg(Fn14_Tert_Prop_Subtype_Plot, 15, 13, "Fn14_Tert_Proportion_Plot_PAM50", "outputs/GSE6532/Expression_Plots/Tamoxifen/", RES = 500)

####################################################
# Conducting Cox PH modelling
####################################################
# Initialsing list
COX_MODELS_LIST <- list()
# Attaching DF
attach(TAM_GENE_CLIN_DF)

# Making cox models
Cox_f <- as.formula(Surv(t.rfs, e.rfs) ~ ESR1_TERT)
COX_MODELS_LIST[["Uni_Cox_ESR1_TERT"]] <- coxph(Cox_f, data = TAM_GENE_CLIN_DF)

Cox_f <- as.formula(Surv(t.rfs, e.rfs) ~ er)
COX_MODELS_LIST[["Uni_Cox_ER_IHC"]] <- coxph(Cox_f, data = TAM_GENE_CLIN_DF)

Cox_f <- as.formula(Surv(t.rfs, e.rfs) ~ TNFRSF12A_TERT)
COX_MODELS_LIST[["Uni_Cox_TNFRSF12A_TERT"]] <- coxph(Cox_f, data = TAM_GENE_CLIN_DF)

Cox_f <- as.formula(Surv(t.rfs, e.rfs) ~ er + age + node + size + grade)
COX_MODELS_LIST[["Multi_Cox_ER_IHC"]] <- coxph(Cox_f, data = TAM_GENE_CLIN_DF)

Cox_f <- as.formula(Surv(t.rfs, e.rfs) ~ ESR1_TERT + age + node + size + grade)
COX_MODELS_LIST[["Multi_Cox_ESR1_TERT"]] <- coxph(Cox_f, data = TAM_GENE_CLIN_DF)

Cox_f <- as.formula(Surv(t.rfs, e.rfs) ~ TNFRSF12A_TERT + age + node + size + grade)
COX_MODELS_LIST[["Multi_Cox_TNFRSF12A_TERT"]] <- coxph(Cox_f, data = TAM_GENE_CLIN_DF)

Cox_f <- as.formula(Surv(t.rfs, e.rfs) ~ ESR1_TERT + age + node + size + grade + PAM50)
COX_MODELS_LIST[["Multi_Cox_ESR1_TERT_PAM50"]] <- coxph(Cox_f, data = TAM_GENE_CLIN_DF)

Cox_f <- as.formula(Surv(t.rfs, e.rfs) ~ TNFRSF12A_TERT + age + node + size + grade +PAM50)
COX_MODELS_LIST[["Multi_Cox_TNFRSF12A_TERT_PAM50"]] <- coxph(Cox_f, data = TAM_GENE_CLIN_DF)

# Detaching DF
detach(TAM_GENE_CLIN_DF)

# Saving cox model results
source(file = "scripts/functions/function_extract_cox.R")
dir.create(path = "outputs/GSE6532/Cox/Tamoxifen/", showWarnings = F, recursive = T)
for (i in names(COX_MODELS_LIST)) {
  jhf_extract_cox(COX_MODELS_LIST[[i]], SAVE_PATH = "outputs/GSE6532/Cox/Tamoxifen/", FILENAME = i, OUTPUT = "save.file")
}

####################################################
# Creating KM plots
####################################################
source(file = "scripts/functions/function_kmplot.R")
args(jhf_kmplot)
kmplot_list <- list()
kmplot_list[["Fn14_RFS_kmplot"]] <- jhf_kmplot(INPUT_DF = TAM_GENE_CLIN_DF,
                              GROUP_VEC = TAM_GENE_CLIN_DF$TNFRSF12A_TERT_GROUP,
                              COX_SURV = survival::Surv(t.rfs, e.rfs) ~ TNFRSF12A_TERT,
                              SURV = survival::Surv(t.rfs, e.rfs) ~ TNFRSF12A_TERT_GROUP,
                              YLAB = "RFS",
                              XLAB = "Time (months)")

kmplot_list[["ESR1_RFS_kmplot"]] <- jhf_kmplot(INPUT_DF = TAM_GENE_CLIN_DF,
                              GROUP_VEC = TAM_GENE_CLIN_DF$ESR1_TERT_GROUP,
                              COX_SURV = survival::Surv(t.rfs, e.rfs) ~ ESR1_TERT,
                              SURV = survival::Surv(t.rfs, e.rfs) ~ ESR1_TERT_GROUP,
                              YLAB = "RFS",
                              XLAB = "Time (months)")

kmplot_list[["ER_IHC_RFS_kmplot"]] <- jhf_kmplot(INPUT_DF = TAM_GENE_CLIN_DF,
                              GROUP_VEC = TAM_GENE_CLIN_DF$er,
                              COX_SURV = survival::Surv(t.rfs, e.rfs) ~ er,
                              SURV = survival::Surv(t.rfs, e.rfs) ~ er,
                              YLAB = "RFS",
                              XLAB = "Time (months)", 
                              BINARY = T)

# Saving KMplots
dir.create(path = "outputs/GSE6532/kmplots/Tamoxifen/", recursive = T, showWarnings = F)
for (i in names(kmplot_list)) {
  jhf_save_jpeg(print(kmplot_list[[i]]), 15, 15, i, "outputs/GSE6532/kmplots/Tamoxifen/", RES = 500, MANUAL_PLOT = T )
}
