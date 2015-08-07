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
# linear model function used for data prep in summary table
lm.fun <- function(ts1, ts2, origin = FALSE, alph = 0.05){
  
  tomod <- na.omit(data.frame(ts1, ts2))
  
  # through origin
  if(origin){ 
    
    # model coefs
    mod <- lm(ts1 ~ 0 + ts2, data = tomod)
    mod <- coef(summary(mod))
    
    # ests and se
    slo <- mod[, 'Estimate']
    slose <- mod[, 'Std. Error']
    
    # margin of error for slope
    qtval <- qt(1 - alph/2, df = nrow(tomod) - 1) 
    slom <- qtval * slose
    slorng <- c(slo - slom, slo + slom)
    slorng <- slorng[1] <= 1 & 1 <= slorng[2] # checks if confint for slope covers one
    
    # make bold, italic to indicate estimate is different from one
    if(!slorng) slo <- paste0('{\\bf \\textit{', form_fun(slo), '}}')
    else slo <- form_fun(slo)
    
    out <- slo
    
  # not through origin
  } else {      
   
    # model coefs
    mod <- lm(ts1 ~ ts2, data = tomod)
    mod <- coef(summary(mod))
        
    # ests and se
    int <- mod[1, 'Estimate']
    intse <- mod[1, 'Std. Error'] 
    slo <- mod[2, 'Estimate']
    slose <- mod[2, 'Std. Error'] 
    
    # margin of errors for each
    qtval <- qt(1 - alph/2, df = nrow(tomod) - 2) 
    intm <- qtval * intse
    slom <- qtval * slose
    intrng <- c(int - intm, int + intm)
    intrng <- intrng[1] <= 0 & 0 <= intrng[2] # checks if confint for intercept covers zero
    slorng <- c(slo - slom, slo + slom)
    slorng <- slorng[1] <= 1 & 1 <= slorng[2] # checks if confint for slope covers one
    
    # make bold, italic to indicate estimate is different from zero
    if(!intrng) int <- paste0('{\\bf \\textit{', form_fun(int), '}}')
    else int <- form_fun(int)
    
    # make bold, italic to indicate estimate is different from one
    if(!slorng) slo <- paste0('{\\bf \\textit{', form_fun(slo), '}}')
    else slo <- form_fun(slo)

    out <- c(int, slo)
    
  }
  
  return(out)
  
}

######
# formatting of values in S expressions
form_fun <- function(x, rnd_val = 2, dig_val = 2, nsm_val = 2) {
  format(round(x, rnd_val), digits = dig_val, nsmall = nsm_val)
  }