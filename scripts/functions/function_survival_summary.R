# Function to calculate survival metric
jhf_survival_clin_summary <- function(EVENT, TIME){
  
  # EVENT = GENE_EXP$DMFS_EVENT
  # TIME = GENE_EXP$DMFS_TIME
  
  survival_summary <- ifelse(is.na(EVENT), "Lost to follow up",     
                             ifelse(EVENT == 0 & TIME < 5, "No recurrence < 5", 
                                 ifelse(EVENT == 1 & TIME < 5, "Recurrance < 5",
                                        ifelse(EVENT == 0 & TIME >= 5, "No recurrance >= 5", 
                                               ifelse(EVENT == 1 & TIME >= 5, "Recurrance >= 5", NA)))))
  
  survival_summary <- factor(survival_summary, levels = c("No recurrence < 5", "Recurrance < 5", "No recurrance >= 5", "Recurrance >= 5", "Lost to follow up"))
  
  return(survival_summary)
}

