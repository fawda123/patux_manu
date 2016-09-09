library(tidyverse)
library(WRTDStidal)

seed_val <- 200

mytheme <- theme_minimal() + 
  theme(
    axis.line.x = element_line(), 
    axis.line.y = element_line(), 
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(),
    axis.ticks.x = element_line(),
    axis.ticks.y = element_line(),
    axis.ticks.length = unit(.1, "cm")
  )

##
# create data to plot

qsim <- lnQ_sim(daydat, comps = T, seed = seed_val)
errsim <- lnres_err(daydat, comps = T, seed = seed_val)
all <- lnQ_sim(daydat, seed = seed_val) %>% 
  lnres_err(seed = seed_val) %>% 
  lnres_sim

##
# q sim plots

toplo <- select(qsim[[2]], date, lnQ, seas_fit, seas_res, errs, sim_out)

# lnQ
p <- ggplot(toplo, aes(x = date, y = lnQ)) + 
  geom_line() + 
  mytheme
flnm <- paste0('figs/chlsim/qsim_lnq.tif')
tiff(flnm, height = 2, width = 7, units = 'in', compression = 'lzw', res = 300, family = 'serif')
print(p)
dev.off()
    
# seasfit
p <- ggplot(toplo, aes(x = date, y = lnQ)) + 
  geom_line() + 
  geom_line(aes(y = seas_fit), col = 'blue', size = 2) + 
  mytheme
flnm <- paste0('figs/chlsim/qsim_seasfit.tif')
tiff(flnm, height = 2, width = 7, units = 'in', compression = 'lzw', res = 300, family = 'serif')
print(p)
dev.off()
    
# seasonal residuals
p <- ggplot(toplo, aes(x = date, y = seas_res)) + 
  geom_point(size = 1) + 
  geom_hline(yintercept = 0) +
  mytheme
flnm <- paste0('figs/chlsim/qsim_seas_res.tif')
tiff(flnm, height = 2, width = 7, units = 'in', compression = 'lzw', res = 300, family = 'serif')
print(p)
dev.off()

# seasonal ARIMA sims
p <- ggplot(toplo, aes(x = date, y = as.numeric(errs))) + 
  geom_point(size = 1) + 
  geom_hline(yintercept = 0) +
  mytheme
flnm <- paste0('figs/chlsim/qsim_errs.tif')
tiff(flnm, height = 2, width = 7, units = 'in', compression = 'lzw', res = 300, family = 'serif')
print(p)
dev.off()

# simulated q
p <- ggplot(toplo, aes(x = date, y = sim_out)) + 
  geom_line() + 
  mytheme
flnm <- paste0('figs/chlsim/qsim_sim_out.tif')
tiff(flnm, height = 2, width = 7, units = 'in', compression = 'lzw', res = 300, family = 'serif')
print(p)
dev.off()
    
##
# err simulation plots
wrtdsmod <- errsim[[1]]

toplo <- select(wrtdsmod, date, res, fits, scls) %>% 
  mutate(resids = fits - res)

# res
p <- ggplot(toplo, aes(x = date, y = res)) + 
  geom_line() + 
  mytheme
flnm <- paste0('figs/chlsim/errmod_res.tif')
tiff(flnm, height = 2, width = 7, units = 'in', compression = 'lzw', res = 300, family = 'serif')
print(p)
dev.off()
   
# fits
p <- ggplot(toplo, aes(x = date, y = res)) + 
  geom_line() + 
  geom_line(aes(y = fits), col = 'blue', size = 2) + 
  mytheme
flnm <- paste0('figs/chlsim/errmod_fits.tif')
tiff(flnm, height = 2, width = 7, units = 'in', compression = 'lzw', res = 300, family = 'serif')
print(p)
dev.off() 

# chl resids
p <- ggplot(toplo, aes(x = date, y = resids)) + 
  geom_point(size = 1) + 
  geom_hline(yintercept = 0) +
  mytheme
flnm <- paste0('figs/chlsim/errmod_resids.tif')
tiff(flnm, height = 2, width = 7, units = 'in', compression = 'lzw', res = 300, family = 'serif')
print(p)
dev.off()

toplo <- select(all, date, errs, lnres, scls, lnres_seas, lnres_noQ, lnres_Q)

# chl errs whole time series
p <- ggplot(toplo, aes(x = date, y = errs)) + 
  geom_point(size = 1) + 
  geom_hline(yintercept = 0) +
  mytheme
flnm <- paste0('figs/chlsim/errmod_errs.tif')
tiff(flnm, height = 2, width = 7, units = 'in', compression = 'lzw', res = 300, family = 'serif')
print(p)
dev.off()

# chl scale whole time series
p <- ggplot(toplo, aes(x = date, y = scls)) + 
  geom_line() + 
  mytheme
flnm <- paste0('figs/chlsim/errmod_scls.tif')
tiff(flnm, height = 2, width = 7, units = 'in', compression = 'lzw', res = 300, family = 'serif')
print(p)
dev.off()

##
# all (complete) time series

# chl seas mod, whole time series
p <- ggplot(toplo, aes(x = date, y = lnres)) + 
  geom_line() + 
  geom_line(aes(y = lnres_seas), col = 'blue', size = 2) + 
  mytheme
flnm <- paste0('figs/chlsim/all_fits.tif')
tiff(flnm, height = 2, width = 7, units = 'in', compression = 'lzw', res = 300, family = 'serif')
print(p)
dev.off()

# chl bio
p <- ggplot(toplo, aes(x = date, y = lnres_noQ)) + 
  geom_line() +
  mytheme
flnm <- paste0('figs/chlsim/all_chlbio.tif')
tiff(flnm, height = 2, width = 7, units = 'in', compression = 'lzw', res = 300, family = 'serif')
print(p)
dev.off()

# chl obs
p <- ggplot(toplo, aes(x = date, y = lnres_Q)) + 
  geom_line() +
  mytheme
flnm <- paste0('figs/chlsim/all_chlobs.tif')
tiff(flnm, height = 2, width = 7, units = 'in', compression = 'lzw', res = 300, family = 'serif')
print(p)