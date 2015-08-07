######
# helper functions

######
#r.squared function
#created for data from weighted regression, but works for all models
#residuals are observed - predicted, but taken from model objects (see 'epc_mods.R')
rsq.fun<-function(resid,obs){
  
  require(Metrics)
  
  # get complete cases
  toeval <- data.frame(resid, obs)
  toeval <- na.omit(toeval)
  
  ssr<-sum(toeval$resid^2)
  sst<-sum(se(toeval$obs,mean(toeval$obs)))
  
  return(1 - (ssr/sst))

}

######
#variant of rmse fun in Metrics package but handles na values
#resid is obs - predicted
rmse.fun<-function(resid){
  
  out<-sqrt(mean(resid^2,na.rm=T))
    
  return(out) 
  
}

######
# average difference
ave.fun <- function(ts1, ts2){
  
  ts1 <- sum(ts1, na.rm = TRUE)
  ts2 <- sum(ts2, na.rm = TRUE)
  
  out <- 100 * (ts1 - ts2)/ts2 
  
  return(out)
  
}

######
# formatting of values in S expressions
form_fun <- function(x, rnd_val = 2, dig_val = 2, nsm_val = 2) {
  format(round(x, rnd_val), digits = dig_val, nsmall = nsm_val)
  }