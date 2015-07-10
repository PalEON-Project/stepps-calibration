library(rstan)
library(ggplot2)
library(fields, quietly=TRUE)
library(RColorBrewer)
library(grid)
library(reshape2)
library(plyr)

wd = getwd()

###############################################################################################################
# user pars
###############################################################################################################

path_utils = 'r/utils'
path_data  = 'r/dump'
path_out   = 'output'
path_figs  = 'figures'

suff_dat = '12taxa_mid_comp_ALL_v0.3'

cal_run  = 'cal_pl_Ka_Kgamma_EPs_ALL_v0.3'
pred_run = '12taxa_6341cells_120knots_cal_KW_KGAMMA_PL_umw_1by'

save_plots = TRUE
rescale    = 1e6

###############################################################################################################
# read in data and source utils
###############################################################################################################

source(file.path(path_utils, 'process_funs.r'))
source(file.path(path_utils, 'plot_funs.r'))
source(file.path(path_utils, 'paper_plot_funs.r'))
source(file.path(wd, 'r', 'runs.r'))

# load composition data; r and centers_veg
load(sprintf('%s/cal_data_%s.rdata', path_data, suff_dat))
colnames(centers_veg) = c('x', 'y')

# load potential pollen predictions; pp and centers_pp
load(paste0('figures/', cal_run, '/pp_all.rdata'))
colnames(centers_pp) = c('x', 'y')

# load composition predictions and uncertainty; r_mean, r_sd and centers_r_pred
load(paste0('../stepps-prediction/figures/', pred_run ,'/r_pred.rdata'))
centers_r_pred = centers_r_pred*rescale
colnames(centers_r_pred) = c('x', 'y')

# abb longer taxon name for plotting
taxa_abb = taxa
taxa_abb[which(taxa_abb == 'OTHER.CONIFER')] = 'OTHER CON'
taxa_abb[which(taxa_abb == 'OTHER.HARDWOOD')] = 'OTHER HW'

###############################################################################################################
# plot the modelled composition
###############################################################################################################
limits <- get_limits(centers_veg)

breaks = c(0, 0.01, 0.05, 0.10, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 1)
# plot_data_maps_binned(r_mean, centers=centers_pls, taxa=taxa, K, T, breaks, suff=suff4, save_plots=save_plots)

## plot_data_maps_binned(r, centers_veg, taxa_abb, K, breaks, limits, suff='veg', save_plots, fpath=path_figs)
# plot_data_maps_col(r, centers_comp, taxa_abb, K, limits, suff='veg', save_plots, fpath=path_figs)


###############################################################################################################
# plot the predictions
###############################################################################################################

limits <- get_limits(centers_veg)

# # plot_data_maps_binned(r, centers_veg, taxa_abb, K, breaks, limits, suff='veg', save_plots, fpath=path_figs)
# plot_data_maps_col(r_mean, centers_veg, taxa_abb, K, limits, suff='veg', save_plots, fpath=path_figs)

###############################################################################################################
# plot composition estimates (data), composition predictions, and potential pollen predictions
###############################################################################################################

# composition estimates: remove MI:LP and check ordering
dist_mat = rdist(as.matrix(centers_r_pred), as.matrix(centers_veg))
dist_mat[dist_mat < 1e-6] = 0
idx_in = unlist(apply(dist_mat, 1, function(x) if (any(x == 0)){ which(x==0)} else {0}))
r      = r[idx_in, ]
centers_veg = centers_veg[idx_in,]
centers_veg = centers_r_pred

# predicted pollen: remove MI:LP and check ordering
dist_mat = rdist(as.matrix(centers_r_pred), as.matrix(centers_pp))
dist_mat[dist_mat < 1e-6] = 0
idx_in = unlist(apply(dist_mat, 1, function(x) if (any(x == 0)){ which(x==0)} else {0}))
pp     = pp[idx_in, ]
centers_pp = centers_pp[idx_in,]
centers_pp = centers_r_pred

# names columns for merging
colnames(r_mean) = taxa_abb
colnames(r_sd)   = taxa_abb
colnames(r)      = taxa_abb
colnames(pp)     = taxa_abb

# melt composition estimates; r, centers_veg
dat_veg = data.frame(r, centers_veg)
dat_veg_melt = melt(dat_veg, id.vars=c('x', 'y'))
dat_veg_melt = cbind(dat_veg_melt, type=rep('comp', nrow(dat_veg_melt)))

# melt composition predictions; r_mean, centers_r_mean
dat_pred = data.frame(r_mean, centers_r_pred)
dat_pred_melt = melt(dat_pred, id.vars=c('x', 'y'))
dat_pred_melt = cbind(dat_pred_melt, type=rep('pred', nrow(dat_pred_melt)))

# melt composition predictions; r_mean, centers_r_mean
dat_pred_sd = data.frame(r_sd, centers_r_pred)
dat_pred_sd_melt = melt(dat_pred_sd, id.vars=c('x', 'y'))
dat_pred_sd_melt = cbind(dat_pred_sd_melt, type=rep('sd', nrow(dat_pred_melt)))


# melt potential pollen predictions; pp, centers_pp
dat_pp = data.frame(pp, centers_pp)
dat_pp_melt = melt(dat_pp, id.vars=c('x', 'y'))
dat_pp_melt = cbind(dat_pp_melt, type=rep('pollen', nrow(dat_pp_melt)))

# bind all three data types together
dat = rbind(dat_veg_melt, dat_pred_melt, dat_pp_melt)

# rename and reorder some factors
levels(dat$variable)[levels(dat$variable) == 'OTHER.CON'] = 'OTHER CON'
levels(dat$variable)[levels(dat$variable) == 'OTHER.HW'] = 'OTHER HW'
levels(dat$type) <- c('Data veg', 'Predicted veg', 'Predicted pollen')
dat$type <- factor(dat$type, levels=c('Data veg', 'Predicted veg', 'Predicted pollen'))
dat$type <- factor(dat$type, levels(dat$type)[c(2,1,3)])

# continuous
p <- ggplot() + geom_tile(data=dat, aes(x=x, y=y, fill=value)) + 
  scale_fill_gradientn(colours=tim.colors(), name='Proportions', limits=c(0,1)) + 
  coord_fixed() 
p <- add_map_albers(p, us.shp, limits, rescale)
p <- p + facet_grid(variable~type)
p <- theme_clean(p) 
p <- p + theme(strip.background = element_blank(), strip.text.x = element_blank())#+ theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.5)))
# p <- p + theme(strip.background = element_blank())
print(p)

# discrete binned
breaks = c(0, 0.01, 0.05, 0.10, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 1)
breaklabels = apply(cbind(breaks[1:(length(breaks)-1)], breaks[2:length(breaks)]), 1, 
                    function(r) { sprintf("%0.2f - %0.2f", r[1], r[2]) })

dat$value_binned = cut(dat$value, breaks, include.lowest=TRUE, labels=FALSE)

cols = tim.colors(length(breaks))

p <- ggplot() + geom_tile(data=dat, aes(x=x, y=y, fill=factor(value_binned))) + 
#   geom_tile(data=dat, aes(x=x, y=y, colour=factor(value_binned)), fill=NA) +
  scale_fill_manual(values = cols, labels=breaklabels, name='Proportions') + 
     coord_fixed() 
p <- add_map_albers(p, us.shp, limits, rescale)
p <- p + facet_grid(variable~type)
p <- theme_clean(p) 
# p <- p + theme(strip.background = element_blank(), strip.text.x = element_blank())
p <- p + theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.5)))
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p)

ggsave(p, file='figures/map_grid_plot.pdf', scale=3)#, width=14, units='in')
ggsave(p, file='figures/map_grid_plot.eps', scale=3)#, width=14, units='in')
ggsave(p, file='figures/map_grid_plot.png', scale=3)

####################################################################################################
# comp, preds, sd
####################################################################################################
limits <- get_limits(centers_veg)

# bind all three data types together
dat = rbind(dat_veg_melt, dat_pred_melt, dat_pred_sd_melt)

# rename and reorder some factors
levels(dat$variable)[levels(dat$variable) == 'OTHER.CON'] = 'OTHER CON'
levels(dat$variable)[levels(dat$variable) == 'OTHER.HW'] = 'OTHER HW'
levels(dat$type) <- c('PLS estimates', 'STEPPS veg', 'STEPPS SD')
dat$type <- factor(dat$type, levels=c('PLS estimates', 'STEPPS veg', 'STEPPS SD'))

# # continuous
# p <- ggplot() + geom_raster(data=dat, aes(x=x, y=y, fill=value)) + 
#   scale_fill_gradientn(colours=tim.colors(), name='Proportions', limits=c(0,1)) + 
#   coord_fixed() 
# p <- add_map_albers(p, us.shp, limits, rescale)
# p <- p + facet_grid(variable~type)
# p <- theme_clean(p) 
# p <- p + theme(strip.background = element_blank(), strip.text.x = element_blank())#+ theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.5)))
# # p <- p + theme(strip.background = element_blank())
# print(p)

# discrete binned
breaks = c(0, 0.01, 0.05, 0.10, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 1)
breaklabels = apply(cbind(breaks[1:(length(breaks)-1)], breaks[2:length(breaks)]), 1, 
                    function(r) { sprintf("%0.2f - %0.2f", r[1], r[2]) })

dat$value_binned = cut(dat$value, breaks, include.lowest=TRUE, labels=FALSE)

cols = tim.colors(length(breaks))

p <- ggplot() + geom_tile(data=dat, aes(x=x, y=y, fill=factor(value_binned))) + 
  #   geom_tile(data=dat, aes(x=x, y=y, colour=factor(value_binned)), fill=NA) +
  scale_fill_manual(values = cols, labels=breaklabels, name='Proportion or SD') + 
  coord_fixed() 
p <- add_map_albers(p, us.shp, limits)
p <- p + facet_grid(variable~type)
p <- theme_clean(p) 
# p <- p + theme(strip.background = element_blank(), strip.text.x = element_blank())
p <- p + theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.5)))
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p)

fname = paste0(getwd(), '/figures/map_assess_comp.pdf')
ggsave(p, file=fname, scale=3)#, width=14, units='in')
# ggsave(p, file='figures/map_assess_comp.eps', scale=3)#, width=14, units='in')
# ggsave(p, file='figures/map_assess_comp.png', scale=3)
sys_str = paste("pdfcrop", fname, fname, sep=' ')
system(sys_str)
