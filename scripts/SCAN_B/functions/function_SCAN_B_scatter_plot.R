###############################################
# Function to produce scatter plot
###############################################


SCAN_B_scatter_plot <- function(DF, X_INPUT, Y_INPUT, METRIC_LAB = "(FPKM Log2)"){

  # DF = SCAN_B_QUERY_DF
  # X_INPUT = "SCNN1A"
  # Y_INPUT = "EMT_METAGENE_SCORE"
  
  Comb_DF = DF
  
  # Calculating Stats
  Meta_Scat_Cor <- cor.test(Comb_DF[[X_INPUT]], Comb_DF[[Y_INPUT]], method = "spearman", alternative = 'two.sided', exact = F)
  p.val_n <- psych::pairwiseCount(Comb_DF[[X_INPUT]], Comb_DF[[Y_INPUT]])
  p_adj_bh <- p.adjust(Meta_Scat_Cor$p.value, method = "BH", n = p.val_n[1])
  
  # Setting labels
  XLAB <- ifelse(grepl("METAGENE", X_INPUT), paste(CLIN_VAR_NAMES[X_INPUT]), paste(CLIN_VAR_NAMES[X_INPUT], METRIC_LAB)) 
  YLAB <- ifelse(grepl("METAGENE", Y_INPUT), paste(CLIN_VAR_NAMES[Y_INPUT]), paste(CLIN_VAR_NAMES[Y_INPUT], METRIC_LAB))
  
  # Creating Scatterplot
  Scat_plot <- ggplot(Comb_DF, aes(Comb_DF[,X_INPUT], Comb_DF[,Y_INPUT], na.rm = TRUE)) +
    geom_point() +
    xlab(XLAB) +
    ylab(YLAB) +
    geom_smooth(method=lm , color="orange", se=FALSE, na.rm = TRUE) +
    ggtitle(label = paste("Spearman Correlation = ", round(Meta_Scat_Cor$estimate, digits = 3), "\n",
                          "P-value = ", format(p_adj_bh, digits = 3), sep = ""))
  
  Scat_plot
}