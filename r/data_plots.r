library(rstan)
library(ggplot2)
library(fields, quietly=TRUE)
library(RColorBrewer)

wd = getwd()

#####################################################################################
# user pars
#####################################################################################

path_utils = 'r/utils'
path_data  = 'r/dump'
path_out   = 'output'
path_figs  = 'figures'

suff_dat = '12taxa_mid_comp_v0.1'

save_plots = TRUE
rescale    = 1e6

#####################################################################################
# read in data and source utils
#####################################################################################

source(file.path(path_utils, 'processFuns.r'))
source(file.path(path_utils, 'plotFuns.r'))
# 
# path_figs1 = sprintf('%s/%s', path_figs, suff_fit)
# if (!file.exists(path_figs1)){
#   dir.create(file.path(wd, path_figs1))
# }



load(sprintf('%s/cal_data_%s.rdata', path_data, suff_dat))

limits <- get_limits(centers_veg)

# plot maps of the veg and pollen
breaks = c(0, 0.01, 0.05, 0.10, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 1)
plot_data_maps_binned(r, centers_veg, taxa, K, breaks, limits, suff='veg', save_plots, fpath=path_figs)
plot_pollen_maps_binned(y, centers_polA, taxa, K, breaks, limits, suff='', save_plots, fpath=path_figs)

# pine
breaks = c(0, 0.01, 0.05, 0.10, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 1)
# breaks = c(0, 0.1, 0.3, 0.4, 1)
plot_both_maps_binned(y,  r, centers_polA, centers_veg, taxa, taxa_list=10, K, breaks, limits, suff, save_plots, fpath=path_figs)
  
# birch
breaks = c(0, 0.01, 0.05, 0.10, 0.15, 0.2, 1)
plot_both_maps_binned(y,  r, centers_polA, centers_veg, taxa, taxa_list=3, K, breaks, limits, suff, save_plots, fpath=path_figs)

# birch
breaks = c(0, 0.01, 0.05, 0.10, 0.15, 0.2, 1)
plot_both_maps_binned(y,  r, centers_polA, centers_veg, taxa, taxa_list=4, K, breaks, limits, suff, save_plots, fpath=path_figs)
