library(rstan)
library(ggplot2)
library(fields, quietly=TRUE)

wd = getwd()

#####################################################################################
# user pars
#####################################################################################

path_utils = 'r/utils'
path_data  = 'r/dump'
path_out   = 'output'
path_figs  = 'figures'

suff=''

rerun = TRUE

# suff_dat = '12taxa_upper_comp_v0.1'
# suff_dat = '12taxa_upper_comp_v0.1'
# suff_dat = '12taxa_mid_comp_v0.1'
# suff_fit = '12taxa_mid_comp_v0.1_vary_psi_bigC_c1'
# suff_fit = '12taxa_mid_comp_v0.1_pl_bigC'
# suff_fit = '12taxa_mid_comp_v0.1_bigC_c1'

# run = list(suff_fit = '12taxa_mid_comp_g_v0.3', kernel = 'gaussian', one_psi = TRUE, one_gamma = TRUE, EPs = FALSE)
# run = list(suff_fit = '12taxa_mid_comp_vary_psi_EPs_v0.3', kernel = 'gaussian', one_psi = FALSE, one_gamma = TRUE, EPs = TRUE)
# run = list(suff_fit = '12taxa_mid_comp_vary_psi_gamma_EPs_v0.3', kernel = 'gaussian', one_psi = FALSE, one_gamma = FALSE, EPs = TRUE)
# run = list(suff_fit = '12taxa_mid_comp_vary_gamma_EPs_v0.3', kernel = 'gaussian', one_psi = TRUE, one_gamma = FALSE, EPs = TRUE)

# run = list(suff_fit = '12taxa_mid_comp_vary_psi_v0.3', kernel = 'gaussian', one_psi = FALSE, one_gamma = TRUE, EPs = FALSE)

# kernel     =  'pl'
# suff_fit = '12taxa_mid_comp_pl_v0.3'
# suff_fit = '12taxa_mid_comp_pl_vary_kernel_gamma_EPs_v0.3'

kernel   = run$kernel
suff_fit = run$suff_fit
suff_dat = run$suff_dat

if (kernel == 'gaussian'){
  one_psi    = run$one_psi
  one_gamma  = run$one_gamma
  EPs        = run$EPs
} else if (kernel == 'pl'){
  one_a      = run$one_a
  one_b      = run$one_b
  one_gamma  = run$one_gamma
  EPs        = run$EPs
}

save_plots = TRUE
rescale    = 1e6

#####################################################################################
# read in data and source utils
#####################################################################################

source(file.path(path_utils, 'process_funs.r'))
source(file.path(path_utils, 'plot_funs.r'))
  
path_figs1 = sprintf('%s/%s', path_figs, suff_fit)
if (!file.exists(path_figs1)){
  dir.create(file.path(wd, path_figs1))
}

load(sprintf('%s/cal_data_%s.rdata', path_data, suff_dat))

fname = sprintf('%s/%s.csv', path_out, suff_fit)

# fix fixup.pl
system(sprintf('r/fixup.pl %s', fname))
fit <- read_stan_csv(fname)
post = rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)

#####################################################################################
# read in data and source utils
#####################################################################################
# npars   = K # always have K phis
# if (kernel=='gaussian'){
#   if (one_psi){ npars = npars + 1 } else { npars = npars + K}
#   if (one_gamma){ npars = npars + 1 } else { npars = npars + K}
#   if (EPs & !one_psi){ npars = npars + 2} # mu and sigma
#   if (EPs & !one_gamma){ npars = npars + 2 + K} # mu and sigma, plus log_gamma
# } else if (kernel=='pl'){
#   if (one_gamma){npars = npars + 1} else {npars = npars + K}
#   if (one_a){npars = npars + 1} else {npars = npars + K}
#   if (one_b){npars = npars + 1} else {npars = npars + K}
#   if (EPs & !one_gamma){ npars = npars + 2} # mu and sigma, plus log_gamma
#   if (EPs & !one_a){ npars = npars + 2 + K} # mu and sigma, plus log_a
#   if (EPs & !one_b){ npars = npars + 2 + K} # mu and sigma, plus log_b
# }

npars   = K # always have K phis
if (kernel=='gaussian'){
  if (one_psi){ npars = npars + 1 } else { npars = npars + K}
  if (one_gamma){ npars = npars + 1 } else { npars = npars + K}
  if (EPs & !one_psi){ npars = npars + 2} # mu and sigma
  if (EPs & !one_gamma){ npars = npars + 2} # mu and sigma, plus log_gamma
} else if (kernel=='pl'){
  if (one_gamma){npars = npars + 1} else {npars = npars + K}
  if (one_a){npars = npars + 1} else {npars = npars + K}
  if (one_b){npars = npars + 1} else {npars = npars + K}
  if (EPs & !one_gamma){ npars = npars + 2} # mu and sigma, plus log_gamma
  if (EPs & !one_a){ npars = npars + 2 } # mu and sigma, plus log_a
  if (EPs & !one_b){ npars = npars + 2 } # mu and sigma, plus log_b
}


# phi <- extract(fit, "phi")$phi
# phi <- extract(fit, "psi")$psi

par_idx = c(seq(1,npars), ncol(post[,1,]))

print(fit)
summary(fit)$summary[,'mean'][par_idx]
ess(fit)
trace_plots(fit, npars, N_cores, suff, save_plots=save_plots, fpath=path_figs1)

print(waic(fit))
print(aic(fit, npars))

log_lik(fit)

# source('r/test_log_lik.r')
# log_lik_iter(fit, N_cores, d, idx_cores, r, N_pot, d_pot, run)

# sink(sprintf('%s/%s/summary.txt', wd, path_figs1), type='output')
sink(sprintf('%s/summary.txt', path_figs1), type='output')
suff_fit
cat('\n')
print('The taxa modelled are:')
print(taxa)
cat('\n')
print('Summary of posterior parameter vals:')
print(get_quants(fit, npars))
cat('\n')
print('WAIC:')
print(waic(fit))
print('AIC:')
print(aic(fit, npars))
cat('\n')
print('Log likelihood')
print(log_lik(fit))
# unlink(sprintf('%s/%s/summary.txt', wd, path_figs1))
# unlink(sprintf('%s/summary.txt', path_figs1))
sink()

#####################################################################################
# compute preds and plot results
#####################################################################################

plot_par_vals(post, parname='phi', taxa, wd, path_figs1)

if (kernel=='gaussian'){ 
  if (!one_psi){
    plot_par_vals(post, parname='psi', taxa, wd, path_figs1)
  }
}
if (!one_gamma){
  plot_par_vals(post, parname='gamma', taxa, wd, path_figs1)
}
if (kernel=='pl'){
  if (!one_a){
    plot_par_vals(post, parname='a', taxa, wd, path_figs1)
  }
  if (!one_b){
  plot_par_vals(post, parname='b', taxa, wd, path_figs1)
  }
}

pollen_props = compute_props(y, taxa)

# scale the veg by phi
local_preds  = phi_scale_veg(fit, N_cores, r, idx_cores)

local_pollen_veg_plot2(r, idx_cores, pollen_props, local_preds, taxa, suff, save_plots, fpath=path_figs1)

sum_w <- build_sumw_pot(post, K, N_pot, d_pot, run)

preds_out = pollen_preds(post, N_cores, d, idx_cores, r, sum_w, run)

alpha = preds_out$alpha # DM precision pars
preds = preds_out$preds
# sum_w = preds_out$sum_w

pollen_preds_plot(preds, pollen_props, N_cores, r, idx_cores, taxa, suff=suff, save_plots=save_plots, fpath=path_figs1)

# # THIS IS WRONG!
# vn_hood_props = sum_hood_props(post, C, N_pot, d_pot, kernel=kernel)
# vn_hood_props
# 
# #####################################################################################
# # potential pollen maps
# #####################################################################################

d_all = t(rdist(as.matrix(centers_veg), as.matrix(centers_veg))/rescale)

N_locs = nrow(d_all)

idx_locs = seq(1, N_locs)

if (rerun){
  preds_out = pollen_preds(post, N_locs, d_all, idx_locs, r, sum_w, run)
  preds     = preds_out$preds
  save(preds, file=paste0(path_figs1, '/preds_all.rdata'))
} else {
  load(file=paste0(path_figs1, '/preds_all.rdata'))
}
  
limits <- get_limits(centers_veg)

breaks = c(0, 0.01, 0.05, 0.10, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 1)
# plot_data_maps_binned(r_mean, centers=centers_pls, taxa=taxa, K, T, breaks, suff=suff4, save_plots=save_plots)

plot_data_maps_binned(preds, centers_veg, taxa, K, breaks, limits, suff='predicted_pollen', save_plots, fpath=path_figs1)

plot_data_maps(preds, centers_veg, taxa, K, limits, suff='predicted_pollen', save_plots, fpath=path_figs1)

# #####################################################################################
# # sum_w map
# #####################################################################################
# 
# dat = data.frame(centers_polA, sum_w=sum_w)
# plot_sumw(dat, fpath=path_figs1)
# 
# #####################################################################################
# # alpha map
# #####################################################################################
# 
# dat = data.frame(centers_polA, alpha=alpha)
# plot_alpha(dat, fpath=path_figs1)
# 
# #####################################################################################
# # proportion of pollen falling within a boundary versus radius
# #####################################################################################
# 
# radius = seq(8000,1000000, by=4000)
# 
# x_pot = seq(-528000, 528000, by=8000)
# y_pot = seq(-416000, 416000, by=8000)
# coord_pot = expand.grid(x_pot, y_pot)
# 
# dmat = t(rdist(matrix(c(0,0), ncol=2), as.matrix(coord_pot, ncol=2))/rescale)
# 
# r_int = dispersal_decay(post, dmat, C, radius, kernel=kernel)
# 
# fifty  = which.min(abs(r_int - 0.5))
# ninety = which.min(abs(r_int - 0.9))
# 
# segments = data.frame(x=c(0, 0, radius[fifty]/1e3, radius[ninety]/1e3), 
#                       xend=c(radius[fifty]/1e3, radius[ninety]/1e3, radius[fifty]/1e3, radius[ninety]/1e3),
#                       y=c(r_int[fifty], r_int[ninety], 0.2, 0.2),
#                       yend=c(r_int[fifty], r_int[ninety], r_int[fifty], r_int[ninety]))
# 
# dat = data.frame(radius=radius/1e3, pollen=r_int)
# 
# p <- ggplot(dat) + geom_line(aes(x=radius, y=pollen))
# p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#         xlab('Radius') + ylab('Proportion of pollen')
# p <- p + theme(axis.text= element_text(size=rel(1)), axis.title=element_text(size=14))
# p <- p + geom_segment(data=segments, aes(x=x, y=y, xend=xend, yend=yend), linetype=2, colour='royalblue4')
# p <- p + xlim(0, segments$xend[4] + segments$xend[4]/16) + ylim(0,1.0)
# p
# 
# ggsave(file=paste(path_figs1, '/dispersal_vs_distance.pdf', sep=''), scale=1)
# ggsave(file=paste(path_figs1, '/dispersal_vs_distance.eps', sep=''), scale=1)
# 
# 
# dvec = seq(0, 1, by=0.0001)
# 
# if (kernel=='gaussian'){
#   colsubstr = substr(colnames(post[,1,]),1,3)
#   psi = mean(post[,1,which(colsubstr == 'psi')])
#   px = gaussian(dvec, psi)
#   plot(dvec*1e3, px, type='l', ylab='Density', xlab='Distance')
# } else if (kernel=='pl'){
#   colsubstr = substr(colnames(post[,1,]),1,3)
#   a = mean(post[,1,which(colsubstr == 'a')])
#   b = mean(post[,1,which(colsubstr == 'b')])
#   px=power_law(dvec, a, b)
#   plot(dvec*1e3, px, type='l', ylab='Density', xlab='Distance')
# }
# 
# 
# #####################################################################################
# # core locations
# #####################################################################################
# # see http://stackoverflow.com/questions/23488022/ggmap-stamen-watercolor-png-error
# # library(ggmap)
# # 
# # centers = data.frame(x=limits$x, y=limits$y)
# # 
# # coordinates(centers) <- ~x + y
# # proj4string(centers) <- CRS('+init=epsg:3175')
# # 
# # centers_ll <- spTransform(centers, CRS('+proj=longlat +ellps=WGS84'))
# # bbox <- as.matrix(data.frame(centers_ll))
# # bbox <- c(bbox[1,1], bbox[1,2], bbox[2,1], bbox[2,2])
# # names(bbox) <- c('left','bottom','right','top')
# # stamen <- get_stamenmap(bbox, zoom = 18)
# # ggmap(stamen) +
# #   geom_point(aes(x = lon, y = lat), data = gc, colour = 'red', size = 2)