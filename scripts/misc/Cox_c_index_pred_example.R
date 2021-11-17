

library(rms)
surv.obj=with(veteran,Surv(time,status))   ####This will be used for rcorr.cens
cox.mod=cph(surv.obj~celltype+karno,data=veteran,x=T,y=T,surv=TRUE, time.inc = 365) # calculate at single timepoint

##Here is the test data set that is the external "independent" data.
test_dat=data.frame(trt=replicate(137,NA), celltype=replicate(137,NA), time=replicate(137,NA), status=replicate(137,NA), karno=replicate(137, NA), diagtime=replicate(137,NA), age=replicate(137,NA), prior=replicate(137,NA))
for(i in seq(8)){
  test_dat[,i]=sample(veteran[,i],137,replace=T)
}

###Create your survival estimates
survest(cox.mod,newdata=test_dat, times = 365)
estimates=survest(cox.mod,newdata=test_dat, times = 365)$surv


###Determine concordance
rcorr.cens(x=estimates,S=surv.obj)
