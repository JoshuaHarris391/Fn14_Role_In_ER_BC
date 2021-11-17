jhf_save_jpeg <- function(PLOT, W, H, NAME, DEST, RES = 1000, MANUAL_PLOT = F){
  
  # Saving Figure
  jpeg(paste(DEST, NAME,".jpeg", sep = ""), width = W, height = H, units = "cm", res = RES)
  
  if(MANUAL_PLOT == F){
    plot(PLOT)
    dev.off()
  }else if (MANUAL_PLOT == T){
    PLOT
    dev.off()
  }

  
}
