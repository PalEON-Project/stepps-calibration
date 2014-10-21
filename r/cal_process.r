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

suff_dat = '12taxa_lower_comp_v0.1'
# suff_fit = '12taxa_mid_comp_v0.1_chain2'
# suff_fit = '12taxa_mid_comp_v0.1_bigC_c3'
suff_fit = '12taxa_lower_comp_v0.1_bigC_c1'
# suff_fit = '12taxa_mid_comp_vary_psi_v0.1_chain1'

one_psi    = TRUE
save_plots = TRUE
rescale    = 1e6

#####################################################################################
# read in data and source utils
#####################################################################################

source(file.path(path_utils, 'processFuns.r'))
source(file.path(path_utils, 'plotFuns.r'))
  
path_figs1 = sprintf('%s/%s', path_figs, suff_fit)
if (!file.exists(path_figs1)){
  dir.create(file.path(wd, path_figs1))
}

load(sprintf('%s/cal_data_%s.rdata', path_data, suff_dat))
fit <- read_stan_csv(sprintf('%s/%s.csv', path_out, suff_fit))
post = rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)

#####################################################################################
# read in data and source utils
#####################################################################################

if (one_psi){
  npars = K+2
} else {
  npars = 2*K+2
}

print(fit)
summary(fit)$summary[,'mean'][1:npars]
ess(fit)
trace_plots(post, suff, save_plots=save_plots, fpath=path_figs1)

# sink(sprintf('%s/%s/summary.txt', wd, path_figs1), type='output')
sink(sprintf('%s/summary.txt', path_figs1), type='output')
print('The taxa modelled are:')
taxa
cat('\n')
print('Summary of posterior parameter vals:')
get_quants(fit)
# unlink(sprintf('%s/%s/summary.txt', wd, path_figs1))
# unlink(sprintf('%s/summary.txt', path_figs1))
sink()

#####################################################################################
# compute preds and plot results
#####################################################################################

colsubstr = substr(colnames(post[,1,]),1,3)
phi       = post[,1,which(colsubstr == 'phi')]
phi_mean  = colMeans(phi)
phi_lb = apply(phi, 2, function(x) quantile(x, probs=0.025))
phi_ub = apply(phi, 2, function(x) quantile(x, probs=0.975))

phi_stats = data.frame(name=taxa, mu=phi_mean, lb=phi_lb, ub=phi_ub)

p <- ggplot(data=phi_stats, aes(x=reorder(name, mu), y=mu)) + geom_point(size=4) + geom_errorbar(aes(ymin=lb, ymax=ub), width=.2) + 
       xlab("Taxa") + ylab(expression(phi)) +
       coord_flip() + theme_bw() + theme(axis.title.x=element_text(size=20), 
                                    axis.title.y=element_text(size=20), 
                                    axis.text.x=element_text(size=rel(1.3)),
                                    axis.text.y=element_text(size=rel(1.3)))

print(p)
#   theme(axis.text.x = element_text(angle=60, hjust=1)) + 
  #   theme(panel.background = element_blank())
ggsave(p, filename=sprintf('%s/%s/phi.pdf', wd, path_figs1), width=8, height=6)


if (!one_psi){
  colsubstr = substr(colnames(post[,1,]),1,3)
  psi       = post[,1,which(colsubstr == 'psi')]
  psi_mean  = colMeans(psi)
  psi_lb = apply(psi, 2, function(x) quantile(x, probs=0.025))
  psi_ub = apply(psi, 2, function(x) quantile(x, probs=0.975))
  
  psi_stats = data.frame(name=taxa, mu=psi_mean, lb=psi_lb, ub=psi_ub)
  
  p <- ggplot(data=psi_stats, aes(x=reorder(name, mu), y=mu)) + geom_point(size=4) + geom_errorbar(aes(ymin=lb, ymax=ub), width=.2) + 
    xlab("Taxa") + ylab(expression(psi)) +
    coord_flip() + theme_bw() + theme(axis.title.x=element_text(size=20), 
                                      axis.title.y=element_text(size=20), 
                                      axis.text.x=element_text(size=rel(1.3)),
                                      axis.text.y=element_text(size=rel(1.3)))
  
  print(p)
  #   theme(axis.text.x = element_text(angle=60, hjust=1)) + 
  #   theme(panel.background = element_blank())
  ggsave(p, filename=sprintf('%s/%s/psi.pdf', wd, path_figs1), width=8, height=6)
}


pollen_props = compute_props(y, taxa)

# scale the veg by phi
local_preds  = phi_scale_veg(fit, N_cores, r, idx_cores)

local_pollen_veg_plot2(r, idx_cores, pollen_props, local_preds, taxa, suff, save_plots, fpath=path_figs1)

C <- compute_C(post, N_pot, d_pot)

preds_out = pollen_preds(post, N_cores, d, idx_cores, r, C, one_psi)

alpha = preds_out$alpha # DM precision pars
preds = preds_out$preds
sum_w = preds_out$sum_w

#postscript(paste('calibration/figures/pollen_fit_plot_', suff_fit, '.eps', sep=''), width=10, height=10)
# pdf(paste('r/pollen/figures/pollen_fit_plot_', suff_fit, '.pdf', sep=''), width=10, height=10)
pollen_preds_plot(preds, pollen_props, N_cores, r, idx_cores, taxa, suff=suff, save_plots=save_plots, fpath=path_figs1)
# dev.off()

#####################################################################################
# potential pollen maps
#####################################################################################

d_all = t(rdist(as.matrix(centers_veg), as.matrix(centers_veg))/rescale)

N_locs = nrow(d_all)

idx_locs = seq(1, N_locs)

preds_out = pollen_preds(post, N_locs, d_all, idx_locs, r, C, one_psi)
preds = preds_out$preds


limits <- get_limits(centers_veg)

breaks = c(0, 0.01, 0.05, 0.10, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 1)
# plot_data_maps_binned(r_mean, centers=centers_pls, taxa=taxa, K, T, breaks, suff=suff4, save_plots=save_plots)

plot_data_maps_binned(preds, centers_veg, taxa, K, breaks, limits, suff='pollen', save_plots, fpath=path_figs1)

plot_data_maps_binned(r, centers_veg, taxa, K, breaks, limits, suff='veg', save_plots, fpath=path_figs1)

#####################################################################################
# sum_w map
#####################################################################################

foo = data.frame(centers_polA, sum_w=sum_w)
p <- ggplot() + geom_point(data=foo, aes(x=x, y=y, size=sum_w), shape=21, colour="black", fill="chartreuse4")
p <- add_map_albers(p, us.shp, limits, rescale) + theme(panel.grid.major = element_blank(), 
                                                        panel.grid.minor = element_blank())
p <- p + theme(axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank())
p 

fpath=path_figs1
ggsave(file=paste(fpath, '/sum_w_maps.pdf', sep=''), scale=1)
ggsave(file=paste(fpath, '/sum_w_maps.eps', sep=''), scale=1)

#####################################################################################
# alpha map
#####################################################################################

foo = data.frame(centers_polA, alpha=alpha)
p <- ggplot() + geom_point(data=foo, aes(x=x, y=y, size=alpha), shape=21, colour="black", fill="chartreuse4")
p <- add_map_albers(p, us.shp, limits, rescale) + theme(panel.grid.major = element_blank(), 
                                                        panel.grid.minor = element_blank())
p <- p + theme(axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank())
p 

fpath=path_figs1
ggsave(file=paste(fpath, '/alpha_maps.pdf', sep=''), scale=1)
ggsave(file=paste(fpath, '/alpha_maps.eps', sep=''), scale=1)

#####################################################################################
# core locations
#####################################################################################

radius = seq(8000,1000000, by=4000)
# radius = radius/1e6

# preds_out = pollen_preds_distance(post, N_cores, d, idx_cores, r, C, radius)
# 
# preds_tot  = preds_out$preds_tot
# preds_int  = preds_out$preds_int
# preds_dist = preds_out$preds_dist


x_pot = seq(-528000, 528000, by=8000)
y_pot = seq(-416000, 416000, by=8000)
coord_pot = expand.grid(x_pot, y_pot)

d_pot = t(rdist(matrix(c(0,0), ncol=2), as.matrix(coord_pot, ncol=2))/rescale)

r_int = dispersal_decay(post, d_pot, C, radius)

fifty  = which.min(abs(r_int - 0.5))
ninety = which.min(abs(r_int - 0.9))

segments = data.frame(x=c(0, 0, radius[fifty]/1e3, radius[ninety]/1e3), 
                      xend=c(radius[fifty]/1e3, radius[ninety]/1e3, radius[fifty]/1e3, radius[ninety]/1e3),
                      y=c(r_int[fifty], r_int[ninety], 0.2, 0.2),
                      yend=c(r_int[fifty], r_int[ninety], r_int[fifty], r_int[ninety]))

foo = data.frame(radius=radius/1e3, pollen=r_int)

p <- ggplot(foo) + geom_line(aes(x=radius, y=pollen))
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        xlab('Radius') + ylab('Proportion of pollen')
p <- p + theme(axis.text= element_text(size=rel(1)), axis.title=element_text(size=14))
p <- p + geom_segment(data=segments, aes(x=x, y=y, xend=xend, yend=yend), linetype=2, colour='royalblue4')
p <- p + xlim(0, 500) + ylim(0.2,1.0)
p

ggsave(file=paste(fpath, '/dispersal_vs_distance.pdf', sep=''), scale=1)
ggsave(file=paste(fpath, '/dispersal_vs_distance.eps', sep=''), scale=1)

# plot(radius/1e3, r_int, type='l')
#####################################################################################
# core locations
#####################################################################################
# see http://stackoverflow.com/questions/23488022/ggmap-stamen-watercolor-png-error
# library(ggmap)
# 
# centers = data.frame(x=limits$x, y=limits$y)
# 
# coordinates(centers) <- ~x + y
# proj4string(centers) <- CRS('+init=epsg:3175')
# 
# centers_ll <- spTransform(centers, CRS('+proj=longlat +ellps=WGS84'))
# bbox <- as.matrix(data.frame(centers_ll))
# bbox <- c(bbox[1,1], bbox[1,2], bbox[2,1], bbox[2,2])
# names(bbox) <- c('left','bottom','right','top')
# stamen <- get_stamenmap(bbox, zoom = 18)
# ggmap(stamen) +
#   geom_point(aes(x = lon, y = lat), data = gc, colour = 'red', size = 2)