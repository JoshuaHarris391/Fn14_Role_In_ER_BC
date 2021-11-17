# Loading in data
library(tidyverse)
load(file = "data/GSE6532_LUMINAL.RData")

# Checking data
glimpse(annot.untreated) # Probe annotations
# View(data.untreated[1:10, 1:10]) # Gene expression values as log2, Patients rows, probe cols
glimpse(demo.untreated) # clinical info

# Converting annotation matrix to dataframe
annot.untreated <- annot.untreated %>% as.data.frame()

####################################################
# Getting basic stats about data
####################################################
data.untreated %>% na.omit() %>% median()
data.untreated %>% na.omit() %>% sd()
data.untreated %>% na.omit() %>% mean()
data.untreated %>% na.omit() %>% range()
plot(density((na.omit(data.untreated))))

####################################################
# Extracting genes of interest
####################################################
INPUT_GENES <- c("TNFRSF12A", "ESR1")
# Finding probe names
m <- match(INPUT_GENES, annot.untreated$HUGO.gene.symbol)
INPUT_GENES_PROBES <- annot.untreated[m, "id"]
# Constructing DF with gene expression and clinical data
UNTREAT_GENE_CLIN_DF <- data.untreated[ , INPUT_GENES_PROBES] %>% as.data.frame()
# Renaming probes
colnames(UNTREAT_GENE_CLIN_DF) <- INPUT_GENES
# Adding samplename column
UNTREAT_GENE_CLIN_DF$samplename <- rownames(UNTREAT_GENE_CLIN_DF)
# Adding clinical info to GENE_CLIN_DF
UNTREAT_GENE_CLIN_DF <- UNTREAT_GENE_CLIN_DF %>% left_join(demo.untreated, by = "samplename")
# Converting survival time from days to months
UNTREAT_GENE_CLIN_DF$t.rfs <- round(UNTREAT_GENE_CLIN_DF$t.rfs/30.417, digit=0)


####################################################
# Calculating PAM50 Subtype 
####################################################
library(genefu)
# Editing annotation dfs
colnames(annot.untreated)[1] <- "probe"
# Removing spaces and converting to integer
annot.untreated$EntrezGene.ID <- gsub(" ", "", annot.untreated$EntrezGene.ID) 
annot.untreated$EntrezGene.ID <- gsub(" ", "", annot.untreated$EntrezGene.ID) 
annot.untreated$EntrezGene.ID <- annot.untreated$EntrezGene.ID %>% as.integer()

# Calculating PAM50 using genefu
data(pam50.robust)
TAM_PAM50 <- molecular.subtyping(sbt.model = "pam50",data = data.untreated, annot = annot.untreated,do.mapping = TRUE)
# Adding PAM50 to GENE_CLIN_DF
UNTREAT_GENE_CLIN_DF$PAM50 <- TAM_PAM50$subtype

####################################################
# Calculating expression tertiles
####################################################
source(file = "scripts/functions/function_calculate_tertiles.R")

UNTREAT_GENE_CLIN_DF$TNFRSF12A_TERT <- jhf_calculate_tertiles(UNTREAT_GENE_CLIN_DF$TNFRSF12A, RETURN = "TERT")
UNTREAT_GENE_CLIN_DF$ESR1_TERT <- jhf_calculate_tertiles(UNTREAT_GENE_CLIN_DF$ESR1, RETURN = "TERT")

UNTREAT_GENE_CLIN_DF$TNFRSF12A_TERT_GROUP <- jhf_calculate_tertiles(UNTREAT_GENE_CLIN_DF$TNFRSF12A, RETURN = "TERT_GROUP")
UNTREAT_GENE_CLIN_DF$ESR1_TERT_GROUP <- jhf_calculate_tertiles(UNTREAT_GENE_CLIN_DF$ESR1, RETURN = "TERT_GROUP")


####################################################
# Creating violin plots
####################################################
Fn14_Violin_Plot <- ggplot(UNTREAT_GENE_CLIN_DF, aes(PAM50, TNFRSF12A)) +
  geom_violin() +
  stat_summary(fun=median, geom="point", size=2, color="grey") +
  labs(y = "TNFRSF12A (log2 Signal)")

prop_df <- UNTREAT_GENE_CLIN_DF %>% 
  group_by(PAM50, TNFRSF12A_TERT_GROUP) %>% 
  summarise(n = n()) %>% 
  mutate(prop = n/sum(n))

Fn14_Tert_Prop_Subtype_Plot <- ggplot(prop_df, aes(PAM50, prop, fill = TNFRSF12A_TERT_GROUP)) + geom_bar(stat = "identity") 

# Saving plots
dir.create(path = "outputs/GSE6532/Expression_Plots/Untreated/", recursive = T, showWarnings = F)
source(file = "scripts/functions/function_save_jpeg.R")
jhf_save_jpeg(Fn14_Violin_Plot, 15, 13, "Fn14_Log2_Violin_Plot_PAM50", "outputs/GSE6532/Expression_Plots/Untreated/", RES = 500)
jhf_save_jpeg(Fn14_Tert_Prop_Subtype_Plot, 15, 13, "Fn14_Tert_Proportion_Plot_PAM50", "outputs/GSE6532/Expression_Plots/Untreated/", RES = 500)

####################################################
# Conducting Cox PH modelling
####################################################
# Initialsing list
COX_MODELS_LIST <- list()
# Attaching DF
attach(UNTREAT_GENE_CLIN_DF)

# Making cox models
Cox_f <- as.formula(Surv(t.rfs, e.rfs) ~ ESR1_TERT)
COX_MODELS_LIST[["Uni_Cox_ESR1_TERT"]] <- coxph(Cox_f, data = UNTREAT_GENE_CLIN_DF)

Cox_f <- as.formula(Surv(t.rfs, e.rfs) ~ er)
COX_MODELS_LIST[["Uni_Cox_ER_IHC"]] <- coxph(Cox_f, data = UNTREAT_GENE_CLIN_DF)

Cox_f <- as.formula(Surv(t.rfs, e.rfs) ~ TNFRSF12A_TERT)
COX_MODELS_LIST[["Uni_Cox_TNFRSF12A_TERT"]] <- coxph(Cox_f, data = UNTREAT_GENE_CLIN_DF)

Cox_f <- as.formula(Surv(t.rfs, e.rfs) ~ er + age + size + grade)
COX_MODELS_LIST[["Multi_Cox_ER_IHC"]] <- coxph(Cox_f, data = UNTREAT_GENE_CLIN_DF)

Cox_f <- as.formula(Surv(t.rfs, e.rfs) ~ ESR1_TERT + age + size + grade)
COX_MODELS_LIST[["Multi_Cox_ESR1_TERT"]] <- coxph(Cox_f, data = UNTREAT_GENE_CLIN_DF)

Cox_f <- as.formula(Surv(t.rfs, e.rfs) ~ TNFRSF12A_TERT + age + size + grade)
COX_MODELS_LIST[["Multi_Cox_TNFRSF12A_TERT"]] <- coxph(Cox_f, data = UNTREAT_GENE_CLIN_DF)

Cox_f <- as.formula(Surv(t.rfs, e.rfs) ~ ESR1_TERT + age + size + grade + PAM50)
COX_MODELS_LIST[["Multi_Cox_ESR1_TERT_PAM50"]] <- coxph(Cox_f, data = UNTREAT_GENE_CLIN_DF)

Cox_f <- as.formula(Surv(t.rfs, e.rfs) ~ TNFRSF12A_TERT + age + size + grade +PAM50)
COX_MODELS_LIST[["Multi_Cox_TNFRSF12A_TERT_PAM50"]] <- coxph(Cox_f, data = UNTREAT_GENE_CLIN_DF)

# Detaching DF
detach(UNTREAT_GENE_CLIN_DF)

# Saving cox model results
source(file = "scripts/functions/function_extract_cox.R")
dir.create(path = "outputs/GSE6532/Cox/Untreated/", showWarnings = F, recursive = T)
for (i in names(COX_MODELS_LIST)) {
  jhf_extract_cox(COX_MODELS_LIST[[i]], SAVE_PATH = "outputs/GSE6532/Cox/Untreated/", FILENAME = i, OUTPUT = "save.file")
}

####################################################
# Creating KM plots
####################################################
source(file = "scripts/functions/function_kmplot.R")
args(jhf_kmplot)
kmplot_list <- list()
kmplot_list[["Fn14_RFS_kmplot"]] <- jhf_kmplot(INPUT_DF = UNTREAT_GENE_CLIN_DF,
                              GROUP_VEC = UNTREAT_GENE_CLIN_DF$TNFRSF12A_TERT_GROUP,
                              COX_SURV = survival::Surv(t.rfs, e.rfs) ~ TNFRSF12A_TERT,
                              SURV = survival::Surv(t.rfs, e.rfs) ~ TNFRSF12A_TERT_GROUP,
                              YLAB = "RFS",
                              XLAB = "Time (months)")

kmplot_list[["ESR1_RFS_kmplot"]] <- jhf_kmplot(INPUT_DF = UNTREAT_GENE_CLIN_DF,
                              GROUP_VEC = UNTREAT_GENE_CLIN_DF$ESR1_TERT_GROUP,
                              COX_SURV = survival::Surv(t.rfs, e.rfs) ~ ESR1_TERT,
                              SURV = survival::Surv(t.rfs, e.rfs) ~ ESR1_TERT_GROUP,
                              YLAB = "RFS",
                              XLAB = "Time (months)")

kmplot_list[["ER_IHC_RFS_kmplot"]] <- jhf_kmplot(INPUT_DF = UNTREAT_GENE_CLIN_DF,
                              GROUP_VEC = UNTREAT_GENE_CLIN_DF$er,
                              COX_SURV = survival::Surv(t.rfs, e.rfs) ~ er,
                              SURV = survival::Surv(t.rfs, e.rfs) ~ er,
                              YLAB = "RFS",
                              XLAB = "Time (months)", 
                              BINARY = T)

# Saving KMplots
dir.create(path = "outputs/GSE6532/kmplots/Untreated/", recursive = T, showWarnings = F)
for (i in names(kmplot_list)) {
  jhf_save_jpeg(print(kmplot_list[[i]]), 15, 15, i, "outputs/GSE6532/kmplots/Untreated/", RES = 500, MANUAL_PLOT = T )
}
