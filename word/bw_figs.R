
library(httr)
library(tidyverse)
library(scales)
library(wesanderson)
library(grid)
library(gridExtra)
library(Hmisc)
library(mgcv)
library(WRTDStidal)

source('R/funcs.R')

# color palette
cols <- grey_pal()(100) %>% 
  as.character %>% 
  .[1:60]
  
# change default ggplot theme
theme_mine <- function (base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
  theme(
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 9), 
    plot.background = element_rect(fill='transparent', 
      colour = NA),
    panel.background = element_rect(fill='transparent', 
      colour = NA),
    legend.background = element_rect(fill='transparent', 
      colour = NA),
    strip.background = element_rect(fill = 
        alpha(cols[length(cols)], 0.5)),
    legend.key = element_rect(fill = 'transparent', 
      colour = NA)
    )   
}
theme_set(theme_mine())

######
# Appendix B, simulated data figure

# load data, from sim_dat.R in tidal_comp proj
load('data/sims_day.RData')
load('data/sims_mos.RData')

# organize daily and monthly ts for plots
toplo <- select(sims_day, date, lnchla_noQ, sim2, sim1, sim3) %>% 
  tidyr::gather('variable', 'value', -date)
toplo2 <- select(sims_mos, date, lnchla_noQ, sim2, sim1, sim3) %>% 
  tidyr::gather('variable', 'value', -date) 

# factor labels as expressions
labs <- c('Biological', 'No flow', 'Constant flow', 'Increasing flow')
levs <- c('lnchla_noQ', 'sim2', 'sim1', 'sim3')
toplo$variable <- factor(toplo$variable, levels = levs, labels = labs)
toplo2$variable <- factor(toplo2$variable, levels = levs, labels = labs)

# y axis labels
ylabs <- expression(paste('ln-Chl-',italic(a),' (',italic('\u03bc'),'g ',L^-1,')'))

# merge with sampled
p <- ggplot(toplo, aes(x = date, y = value, group = variable)) + 
  geom_point(pch = 1, size = 1, col = cols[40]) + 
  geom_line(data = toplo2, aes(x = date, y = value, group = variable), col = cols[1]) +
  scale_y_continuous(ylabs) + 
  theme(axis.title.x = element_blank()) +
  facet_wrap(~variable, ncol = 1)

tiff('word/FIGUREB1_bw.tif', width = 6, height = 7, units = 'in', compression = 'lzw', res = 500, family = 'serif')
print(p)
dev.off()

######
# Appendix C, wrtds v gam predictions by periods

data(bestLE12)
data(bestTF16)

bestLE12 <- select(bestLE12, 
    fits_wrtds, fits_gams, norm_wrtds, norm_gams, flcat, mocat, yrcat, dec_time
  ) %>% 
  gather(., 'var', 'val', fits_wrtds:norm_gams) %>% 
  separate(., var, c('res', 'mod'), by = '_') %>% 
  spread(mod, val) %>% 
  mutate(stat = 'LE1.2')

bestTF16 <- select(bestTF16, 
    fits_wrtds, fits_gams, norm_wrtds, norm_gams, flcat, mocat, yrcat, dec_time
  ) %>% 
  gather(., 'var', 'val', fits_wrtds:norm_gams) %>% 
  separate(., var, c('res', 'mod'), by = '_') %>% 
  spread(mod, val) %>% 
  mutate(stat = 'TF1.6')

toplo <- rbind(bestTF16, bestLE12)
toplo$res <- factor(toplo$res, levels = c('fits', 'norm'), labels = c('Pred', 'Norm'))

# margins
mars <- grid::unit(c(0, 0, 0.1, 0.1), 'cm')

p <- ggplot(toplo, aes(x = gams, y = wrtds, colour = res, group = res)) + 
  geom_abline(intercept = 0, slope = 1, colour = 'grey') + 
  geom_point(aes(shape = res), alpha = 0.6, size = 1) + 
  scale_colour_manual(values = cols[c(1, 50)]) +
  scale_shape_manual(values = c(1, 1)) + 
  scale_x_continuous(limits = c(0, 4)) + 
  scale_y_continuous(limits = c(0, 4)) + 
  geom_smooth(method = 'lm', se = F, size = 0.7, fullrange = TRUE, alpha = 0.7)  + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = 'top', 
    legend.title = element_blank()) 

pleg <- g_legend(p)
p <- p + 
  theme(legend.position = 'none')

p1 <- p + facet_grid(stat ~ yrcat) + theme(plot.margin = mars)
p2 <- p + facet_grid(stat ~ mocat) + theme(plot.margin = mars)
p3 <- p + facet_grid(stat ~ flcat) + theme(plot.margin = mars)

ylab <- expression(paste('WRTDS, ln-Chl-',italic(a),' (',italic('\u03bc'),'g ',L^-1,')'))
xlab <- expression(paste('GAM, ln-Chl-',italic(a),' (',italic('\u03bc'),'g ',L^-1,')'))

tiff('word/FIGUREC1_bw.tif', width = 6, height = 6.5, units = 'in', compression = 'lzw', res = 500, family = 'serif')
grid.arrange(pleg, p1, p2, p3, left = textGrob(ylab, rot = 90), bottom = textGrob(xlab),
  heights = c(0.12, 1, 1, 1))
dev.off()

######
# study site map

# load data
data(cb_poly)
data(usa_poly)
data(patux_poly)
data(pax_meta)

# lims <- bbox(patux_poly)

# color palette
mapcols <- grey_pal()(100) %>% 
  as.character %>% 
  .[1:60]
mapcols <- mapcols[round(seq(1, length(mapcols), length = 3))]

# prep for plots
cb_poly <- fortify(cb_poly)
usa_poly <- fortify(usa_poly)
patux_poly <- fortify(patux_poly) %>% 
  mutate(id = factor(id, levels = c('60', '74', '88'), labels = c('TF', 'OH', 'MH')))

pax_meta$lab <- with(pax_meta, paste0('"', STATION, ' (', round(dist_km, 1), ')"'))
pax_meta$lab[c(4, 8)] <- paste0('bold(', pax_meta$lab[c(4, 8)], ')') 

# base map
p1 <- ggplot(patux_poly, aes(x = long, y = lat)) + 
  geom_polygon(data = cb_poly, aes(x = long, y = lat, group = group), 
    fill = 'grey95', colour = 'grey50') +
  geom_polygon(aes(fill = factor(id), group = group), colour = NA) +
  scale_fill_manual(values = rev(mapcols)) +
  geom_point(data = pax_meta, aes(x = LONG, y = LAT), size = 5, alpha = 0.8) +
  geom_label(data = pax_meta, aes(x = LONG * 1.001, y = LAT, label = lab), alpha = 0.5, size = 5, 
    parse = TRUE, label.padding = grid::unit(0.15, "lines")) +
  theme_classic() + 
  theme(axis.line=element_blank(), axis.text.x=element_blank(),
          axis.text.y=element_blank(), axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(), 
          legend.position = 'top', 
          legend.title = element_blank()
    ) +
  coord_fixed(ratio = 1, xlim = c(-76.84, -76.38), ylim = c(38.27, 38.88))


# inset
p2 <- ggplot(patux_poly, aes(x = long, y = lat)) + 
  geom_polygon(data = usa_poly, aes(x = long, y = lat, group = group), 
    fill = 'white', colour = 'grey50') +
  geom_polygon(aes(group = group), colour = 'black', fill = 'black') +
  theme(axis.title =element_blank(), 
          axis.text.y = element_text(colour = 'black'), 
          axis.text.x = element_text(colour = 'black', angle = 90, vjust = 0.5), 
          panel.background=element_rect(fill = 'grey95'),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, size = 1)
    ) +
  coord_fixed(ratio = 1, xlim = c(-77.5, -75.5), ylim = c(36.7, 39.7))

tiff('word/FIGURE1_bw.tif', width = 5, height = 7, units = 'in', compression = 'lzw', res = 500, family = 'serif')
grid.newpage()
v1 <- viewport(width = 1, height = 1, x = 0.5, y = 0.5) 
v2 <- viewport(width = 0.55, height = 0.55, x = 0.7275, y = 0.64)
print(p1, vp = v1) 
print(p2, vp = v2)
dev.off()

######
# chlorophyll trends by year, month, flow, combined

load(file = 'data/pax_data.RData')
load(file = 'data/pax_meta.RData')
load(file = 'data/pax_clip.RData')

# color for water for continuity with last fig
wt_col <- grey_pal()(100)[10] %>% 
  alpha(0.3)

# format pax_meta for mrg with pax_plo
pax_meta <- select(pax_meta, STATION, LONG, LAT) %>% 
  mutate(STATION = as.character(STATION))

# add month, yr columns, make categories by yr, mo, fl
pax_data <- mutate(pax_data, 
  mo = as.numeric(strftime(date, '%m')),
  yr = as.numeric(strftime(date, '%Y'))
  ) %>% 
  mutate(STATION = as.character(STATION)) %>% 
  mutate(
    yrcat = cut(yr, breaks = c(-Inf, 1993, 2000, 2007, Inf), 
      labels = c('1986-1993', '1994-2000', '2001-2007', '2008-2014')),
    mocat = cut(mo, breaks = c(-Inf, 3, 6, 9, Inf), labels = c('JFM', 'AMJ', 'JAS', 'OND')), 
    flcat = cut(lnQ, breaks = c(-Inf, quantile(lnQ, c(0.25, 0.5, 0.75)), Inf), labels = c('Flow 1 (Low)', 'Flow 2', 'Flow 3', 'Flow 4 (High)'))
  )
  
# get yrcat medians by station, add coords
pax_yr <- group_by(pax_data, STATION, yrcat) %>% 
  summarise(lnchla = median(lnchla, na.rm = T)) %>% 
  left_join(., pax_meta, by = 'STATION') %>% 
  rename(cat = yrcat)

# get mocat medians by station, add coords
pax_mo <- group_by(pax_data, STATION, mocat) %>% 
  summarise(lnchla = median(lnchla, na.rm = T)) %>% 
  left_join(., pax_meta, by = 'STATION') %>% 
  rename(cat = mocat)

# get flocat medians by station, add coords
pax_fl <- group_by(pax_data, STATION, flcat) %>% 
  summarise(lnchla = median(lnchla, na.rm = T)) %>% 
  left_join(., pax_meta, by = 'STATION') %>% 
  rename(cat = flcat)

pax_plo <- rbind(pax_yr, pax_mo, pax_fl)

# plot label
ylabs <- expression(paste('ln-Chl-',italic(a),' (',italic('\u03bc'),'g ',L^-1,')'))

# station labels, first panel only
text_lab <- filter(data.frame(pax_plo), STATION %in% c('TF1.6', 'LE1.2') & cat %in% '1986-1993')

# plot
p1 <- ggplot(pax_meta, aes(x = LONG, y = LAT)) + 
  geom_polygon(data = pax_clip, aes(x = long, y = lat, group = group),
    fill = wt_col) +
  coord_map(
    xlim = c(-76.78, -76.36),
    ylim = c(38.27, 38.85)
  ) +
  geom_text(data = text_lab, aes(x = LONG + 0.09, y = LAT, label = STATION), size = 2.5) + 
  geom_point(data = pax_plo, aes(group = cat, size = lnchla, fill = lnchla), 
    pch = 21) +
  facet_wrap(~cat, ncol = 4) +
  theme_mine() + 
  scale_size(range = c(0.5, 7.5)) + 
  scale_fill_gradientn(colours = cols) +
  theme(
    legend.position = 'top', 
    axis.title = element_blank(), 
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 7, angle = 90, hjust = 1)
    ) + 
  guides(
    fill = guide_legend(title = ylabs), 
    size = guide_legend(title = ylabs)
    )

tiff('word/FIGURE2_bw.tif', width = 5, height = 7.75, units = 'in', compression = 'lzw', res = 500, family = 'serif')
print(p1)
dev.off()
  
######
# predicted trends for each mod - monthly time series of observed with predictions, RMSE in fig

library(dplyr)
library(ggplot2)
library(tidyr)

data(bestLE12)
data(bestTF16)

# some plots
bestTF16$site <- 'TF1.6'
bestLE12$site <- 'LE1.2'
names(bestLE12) <- gsub('gams$', 'GAM', names(bestLE12))
names(bestLE12) <- gsub('wrtds$', 'WRTDS', names(bestLE12))
names(bestTF16) <- gsub('gams$', 'GAM', names(bestTF16))
names(bestTF16) <- gsub('wrtds$', 'WRTDS', names(bestTF16))

ylabs <- expression(paste('ln-Chl-',italic(a),' (',italic('\u03bc'),'g ',L^-1,')'))

bestLE12 <- select(bestLE12, -flo, -lnQ, -dec_time)
bestTF16 <- select(bestTF16, -flo, -sal, -dec_time)
toplo <- rbind(bestTF16, bestLE12) %>% 
  gather(variable, value, fits_WRTDS:res_GAM) %>% 
  separate(variable, c('output', 'model'), sep = '_') %>% 
  mutate(output = factor(output, levels = c('fits', 'norm', 'res'), labels = c('Pred', 'Norm', 'Res')))

p1 <- ggplot(toplo[toplo$output != 'Res', ], aes(x = date, y = value, colour = output)) + 
  geom_point(aes(y = res, shape = 'Obs'), colour = 'black', size = 1.5, alpha = 0.15) +
  geom_line(size = 0.8, alpha = 0.85) + 
  facet_grid(site ~ model) + 
  ylab(ylabs) + 
  theme(
    axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5),
    axis.title.x = element_blank(), 
    legend.title = element_blank()
    ) +
  scale_y_continuous(limits = c(0, 4.5), breaks = c(0.5, 1.5, 2.5, 3.5, 4.5)) + 
  scale_colour_manual(values = cols[c(1, length(cols))]) + 
  guides(color=guide_legend(override.aes=list(shape=c(NA,NA),linetype=c(1,1)))) + 
  ggtitle('(a) Monthly')

toplo2 <- mutate(toplo, year = as.numeric(strftime(date, '%Y'))) %>% 
  group_by(year, site, output, model) %>% 
  summarise(value = mean(value, na.rm = TRUE)) %>% 
  filter(output != 'Res') %>% 
  spread(output, value) 

p2 <- ggplot(toplo2, aes(x = year, y = Pred, colour = 'Pred')) + 
  geom_line(aes(x = year, y = Norm, colour = 'Norm'), size = 1) + 
  geom_point() + 
  facet_grid(site ~ model) + 
  theme(
    axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5),
    axis.title.x = element_blank(), 
    legend.title = element_blank()
    ) +
  scale_y_continuous(limits = c(1.5, 3), breaks = seq(1.5, 3, 0.5)) + 
  scale_colour_manual(values = cols[c(length(cols), 1)]) + 
  ylab(ylabs) +
  guides(color=guide_legend(override.aes=list(shape=c(NA,16),linetype=c(1,0)))) + 
  ggtitle('(b) Annual')

tiff('word/FIGURE3_bw.tif', width = 8, height = 9, units = 'in', compression = 'lzw', res = 500, family = 'serif')
grid.arrange(p1, p2, ncol = 1)
dev.off()

######
# plot of seasonal variation from model predictions, by model and site, 

data(bestLE12)
data(bestTF16)

bestLE12 <- mutate(bestLE12, stat = 'LE1.2')
bestTF16 <- mutate(bestTF16, stat = 'TF1.6')

toplo <- rbind(bestLE12[, !names(bestLE12) %in% 'lnQ'], bestTF16[, !names(bestTF16) %in% 'sal']) %>% 
  select(date, fits_wrtds, fits_gams, res, stat) %>% 
  mutate(
    Year = '2014', 
    mo = as.numeric(strftime(date, '%m')),
    day = as.numeric(strftime(date, '%d')),
    fake_date = as.Date(paste(Year, mo, day, sep = '-'))
  ) %>% 
  gather('mod', 'val', fits_wrtds:res) %>% 
  mutate(
    mod = factor(mod, levels = c('res', 'fits_gams', 'fits_wrtds'), 
      labels = c('observed', 'GAM', 'WRTDS'))
  )
  
ylabs <- expression(paste('ln-Chl-',italic(a),' (',italic('\u03bc'),'g ',L^-1,')'))
xlab <- 'Day of year'

p <- ggplot(toplo, aes(x = fake_date, y = val, group = mod)) + 
  geom_point(pch = 1, col = cols[40]) + 
  stat_smooth(method = 'loess', span = 0.4, se = F, col = cols[1], size = 1) +
  scale_y_continuous(ylabs) + 
  scale_x_date(xlab, labels = date_format("%m/%d")) + 
  facet_grid(stat ~ mod, scales = 'free_y') + 
  theme_mine() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8))

tiff('word/FIGURE4_bw.tif', width = 9, height = 5, units = 'in', compression = 'lzw', res = 500, family = 'serif')
print(p)
dev.off()

######
# dynaplots for each mod (as in Fig. 8 Beck and Hagy 2015)

## LE12 models

data(bestLE12)

LE12_tomod <- select(bestLE12, date, dec_time, res, flo) %>% 
  mutate(
    doy = as.numeric(strftime(date, '%j'))
  )

bestLE12_gam <- gam(res~te(dec_time, doy, flo, bs=c("tp","cc","tp")), k = c(5, 8, 5), data = LE12_tomod, knots=list(doy=c(1,366)))

data(bestLE12_wrtds)

## TF16 models

data(bestTF16)

TF16_tomod <- select(bestTF16, date, dec_time, res, flo) %>% 
  mutate(
    doy = as.numeric(strftime(date, '%j'))
  )

bestTF16_gam <- gam(res~te(dec_time, doy, flo, bs=c("tp","cc","tp")), k = c(5, 10, 5), 
  data = TF16_tomod, 
  knots=list(doy=c(1,366)), na.action = na.pass)

data(bestTF16_wrtds)

## plots
margs <- grid::unit(c(0, 0, 0, 0), "mm") # margins
fac_txt <- 10 # facet text size
ylims <- c(0.6, 4)
col_vec <- grey_pal()(100) %>% 
  as.character %>% 
  .[1:100] %>% 
  rev
mos <- c(1, 4, 7, 10)
alph <- 0.7

p1 <- dynagam(bestLE12_gam, LE12_tomod, ncol = 1, use_bw = F, month = mos, col_vec = col_vec, alpha = alph) +
  theme(legend.position = 'none', axis.title.y = element_blank(), 
    strip.text.x = element_text(size = fac_txt), plot.margin = margs) +
  scale_y_continuous(limits = ylims) +
  scale_x_reverse('Salinity')

p2 <- dynaplot(bestLE12_wrtds, ncol = 1, use_bw = F, month = mos, col_vec = col_vec, alpha = alph) +
  theme(legend.position = 'none', axis.title.y = element_blank(), 
    strip.text.x = element_text(size = fac_txt), axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(), plot.margin = margs) +
  scale_y_continuous(limits = ylims) +
  scale_x_reverse('Salinity')

p3 <- dynagam(bestTF16_gam, TF16_tomod, ncol = 1, use_bw = F, month = mos, col_vec = col_vec, alpha = alph) +
  theme(legend.position = 'none', axis.title.y = element_blank(), 
    strip.text.x = element_text(size = fac_txt), axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(), plot.margin = margs) + 
  scale_x_continuous('Flow') +
  scale_y_continuous(limits = ylims)

p4 <- dynaplot(bestTF16_wrtds, ncol = 1, use_bw = F, month = mos, col_vec = col_vec, alpha = alph) +
  theme(axis.title.y = element_blank(), legend.position = 'right', 
    strip.text.x = element_text(size = fac_txt), axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(), plot.margin = margs) + 
  scale_x_continuous('Flow') +
  scale_y_continuous(limits = ylims) + 
  guides(colour = guide_colourbar(barwidth = 1, barheight = 10)) 
 
pleg <- g_legend(p4)
p4 <- p4 + theme(legend.position = 'none')

ylab <- expression(paste('ln-Chl-',italic(a),' (',italic('\u03bc'),'g ',L^-1,')'))

# arrange all as grob
grobwidths <- c(1, 1, 1, 1)
library(grid)

tiff('word/FIGURE5_bw.tif', width = 7.5, height = 6.5, units = 'in', compression = 'lzw', res = 500, family = 'serif')
grid.arrange(
  arrangeGrob(textGrob(ylab, rot = 90)), 
  arrangeGrob(
    arrangeGrob(textGrob('LE1.2'), textGrob('TF1.6'), ncol = 2), 
    arrangeGrob(textGrob('GAM'), textGrob('WRTDS'), textGrob('GAM'), textGrob('WRTDS'), ncol = 4), 
    arrangeGrob(p1, p2, p3, p4, ncol = 4), 
    heights = c(0.025, 0.05, 1)
  ),
  arrangeGrob(pleg), 
  ncol = 3, 
  widths = c(0.2, 2, 0.3)
  )
dev.off()

######
# dynaplot examples for mods with simulated data

# color palette
cols <- grey_pal()(100) %>% 
  as.character %>% 
  .[1:75]
modcols <- cols[c(1, 50)]
  
data(bestsim_wrtds)
data(bestsim_gam)

source('R/funcs.r')

ylims <- c(0.7, 3)
xlims <- c(2, 3.5)
mars <- grid::unit(c(0.1, 0.1, 0.1, 0.1), 'cm')

##
# wrtds figs

# no effect
p1_wrtds <- dynaplot(bestsim_wrtds[[2]], month = 8, col_vec = cols, fac_nms = 'No influence, WRTDS', floscl = F) + 
  scale_y_continuous(limits = ylims) +
  scale_x_continuous(limits = xlims) + 
  theme_mine() +
  theme(legend.position = 'none', axis.title = element_blank(), plot.margin = mars, 
    strip.text.x = element_text(size = NULL))

# constant
p2_wrtds <- dynaplot(bestsim_wrtds[[1]], month = 8, col_vec = cols, fac_nms = 'Constant flow, WRTDS', floscl = F) + 
  scale_y_continuous(limits = ylims) +
  scale_x_continuous(limits = xlims) + 
  theme_mine() +
  theme(legend.position = 'none', axis.title = element_blank(), plot.margin = mars)

# increasing
p3_wrtds <- dynaplot(bestsim_wrtds[[3]], month = 8, col_vec = cols, fac_nms = 'Increasing flow, WRTDS', floscl = F) + 
  scale_y_continuous(limits = ylims) +
  scale_x_continuous(limits = xlims) + 
  theme_mine() +
  theme(legend.position = 'top', axis.title = element_blank(), plot.margin = mars) +
  guides(colour = guide_colourbar(barwidth = 10, barheight = 1)) 

pleg <- g_legend(p3_wrtds)
p3_wrtds <- p3_wrtds + theme(legend.position = 'none')

##
# gam figs

# no effect
p1_gam <- dynagam(bestsim_gam[[2]]$mod, bestsim_gam[[2]]$dat, month = 8, col_vec = cols, fac_nms = 'No influence, GAM') + 
  scale_y_continuous(limits = ylims) +
  scale_x_continuous(limits = xlims) + 
  theme_mine() +
  theme(legend.position = 'none', axis.title = element_blank(), plot.margin = mars, 
    strip.text.x = element_text(size = NULL))

# constant
p2_gam <- dynagam(bestsim_gam[[1]]$mod, bestsim_gam[[1]]$dat, month = 8, col_vec = cols, fac_nms = 'Constant flow, GAM') + 
  scale_y_continuous(limits = ylims) +
  scale_x_continuous(limits = xlims) + 
  theme_mine() +
  theme(legend.position = 'none', axis.title = element_blank(), plot.margin = mars, 
    strip.text.x = element_text(size = NULL))

# increasing
p3_gam <- dynagam(bestsim_gam[[3]]$mod, bestsim_gam[[3]]$dat, month = 8, col_vec = cols, fac_nms = 'Increasing flow, GAM') + 
  scale_y_continuous(limits = ylims) +
  scale_x_continuous(limits = xlims) + 
  theme_mine() +
  theme(legend.position = 'none', axis.title = element_blank(), plot.margin = mars, 
    strip.text.x = element_text(size = NULL))

## 
# final plots 

ylabs <- expression(paste('ln-Chl-',italic(a),' (',italic('\u03bc'),'g ',L^-1,')'))

tiff('word/FIGURE6_bw.tif', width = 7, height = 4.5, units = 'in', compression = 'lzw', res = 500, family = 'serif')
grid.arrange(
  pleg, ncol = 1, heights = c(0.15, 1, 0.1),
  arrangeGrob(
    textGrob(ylabs, rot = 90), ncol = 2,
    arrangeGrob(p1_gam, p2_gam, p3_gam, p1_wrtds, p2_wrtds, p3_wrtds, ncol = 3), 
    widths = c(0.05, 1)
  ), 
  textGrob('Flow')
)
dev.off()
