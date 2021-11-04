jhf_save_jpeg <- function(PLOT, W, H, NAME, DEST, RES = 1000){
  
  # Saving Figure
  jpeg(paste(DEST, NAME,".jpeg", sep = ""), width = W, height = H, units = "cm", res = RES)
  plot(PLOT)
  dev.off()
  
}
