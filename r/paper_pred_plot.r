library(rstan)
library(ggplot2)
library(fields, quietly=TRUE)
library(RColorBrewer)
library(grid)
library(reshape2)

wd = getwd()

###############################################################################################################
# user pars
###############################################################################################################

path_utils = 'r/utils'
path_data  = 'r/dump'
path_out   = 'output'
path_figs  = 'figures'

suff_dat = '12taxa_mid_comp_ALL_v0.3'

save_plots = TRUE
rescale    = 1e6

###############################################################################################################
# read in data and source utils
###############################################################################################################

source(file.path(path_utils, 'process_funs.r'))
source(file.path(path_utils, 'plot_funs.r'))
source(file.path(path_utils, 'paper_plot_funs.r'))
source(file.path(wd, 'r', 'runs.r'))

load(sprintf('%s/cal_data_%s.rdata', path_data, suff_dat))

###############################################################################################################
# plot the modelled composition
###############################################################################################################
limits <- get_limits(centers_veg)

breaks = c(0, 0.01, 0.05, 0.10, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 1)
# plot_data_maps_binned(r_mean, centers=centers_pls, taxa=taxa, K, T, breaks, suff=suff4, save_plots=save_plots)

plot_data_maps_binned(r, centers_veg, taxa, K, breaks, limits, suff='veg', save_plots, fpath=path_figs)

###############################################################################################################
###############################################################################################################

taxa[which(taxa == 'OTHER.CONIFER')] = 'OTHER CON'
taxa[which(taxa == 'OTHER.HARDWOOD')] = 'OTHER HW'

pollen_props      = compute_props(y, taxa)
pollen_props_melt = melt(compute_props(y, taxa))

colnames(r) = taxa
r_melt      = melt(r[idx_cores, ])

phi_scaled_dat = data.frame(matrix(0, nrow=0, ncol=6))
pred_dat       = data.frame(matrix(0, nrow=0, ncol=6))
for (run in runs){
  fname = sprintf('%s/%s.csv', path_out, run$suff_fit)
  system(sprintf('r/fixup.pl %s', fname))
  fit   = read_stan_csv(fname)
  post  = rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)
  
  sum_w <- build_sumw_pot(post, K, N_pot, d_pot, run)
  preds_out = pollen_preds(post, N_cores, d, idx_cores, r, sum_w, run)
  preds     = preds_out$preds
  colnames(preds) = taxa
  
  preds_melt = melt(preds)
  preds_melt = cbind(preds_melt, pollen_props_melt[,3], r_melt[,3], rep(run$handle, nrow(preds_melt)))
  colnames(preds_melt) = c('core', 'taxon', 'pred', 'data', 'veg', 'run')
  
  pred_dat = rbind(pred_dat, preds_melt)

  phi_out = phi_scale_veg(post, N_cores, r, idx_cores)
  colnames(phi_out) = taxa

  phi_out_melt = melt(phi_out)
  phi_out_melt = cbind(phi_out_melt, pollen_props_melt[,3], r_melt[,3], rep(run$handle, nrow(phi_out_melt)))
  colnames(phi_out_melt) = c('core', 'taxon', 'phi_pred', 'data', 'veg', 'run')
  phi_scaled_dat = rbind(phi_scaled_dat, phi_out_melt)
}

handles = as.vector(unique(pred_dat$run))
for (run in handles) {

  dat        = pred_dat[which(pred_dat$run %in% run),]
  phi_scaled = phi_scaled_dat[which(phi_scaled_dat$run %in% run),]
  
  plot_pollen_preds_paper(dat, path_figs, run)
  plot_pollen_phi_scaled_paper(phi_scaled, path_figs, run) 
}

############################################################################################################
# dispersal cdf
############################################################################################################

radius = seq(8000,1000000, by=4000)

# x_pot = seq(-528000, 528000, by=8000)
# y_pot = seq(-416000, 416000, by=8000)
# coord_pot = expand.grid(x_pot, y_pot)
# 
# dmat = t(rdist(matrix(c(0,0), ncol=2), as.matrix(coord_pot, ncol=2))/rescale)

coord_pot = seq(-700000, 700000, by=8000)
coord_pot = expand.grid(coord_pot, coord_pot)

d_pot = t(rdist(matrix(c(0,0), ncol=2), as.matrix(coord_pot, ncol=2))/dist.scale)

# want a circular region
idx_circ  = which(d_pot[,1] > 0.7)
coord_pot = coord_pot[idx_circ, ]
d_pot     = d_pot[idx_circ, ]

dat = data.frame(matrix(0, nrow=0, ncol=4))
for (run in runs){
  fname = sprintf('%s/%s.csv', path_out, run$suff_fit)
  fit   = read_stan_csv(fname)
  post  = rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)
  
  sum_w = build_sumw_pot(post, K, N_pot, d_pot, run)
  r_int = dispersal_decay(post, d_pot, sum_w, radius/rescale, run, taxa)
  
  run_dat = data.frame(radius = radius/1e3, r_int, handle=rep(run$handle, length(radius)))
  dat = rbind(dat, melt(run_dat, id=c('radius', 'handle')))
}

limits <- get_limits(centers_veg)

p <- ggplot(dat) + geom_line(data=dat, aes(x=radius, y=value, colour=factor(handle))) 
# p <- p + scale_colour_manual("Model", values=c('red', 'blue', 'black', 'green')) 
p <- p + facet_wrap(~variable, ncol=3) + theme_bw() + theme(panel.grid.major = element_blank(), 
                                                            panel.grid.minor = element_blank())

# p <- p + scale_fill_discrete(name="Model",
#                   breaks=c("G", "G_Kpsi_Kgamma", "PL", "PL_Ka_Kgamma"),
#                   labels=c("Gaussian", "Gaussian Kpsi Kgamma", "Power-law", "Power-law Ka Kgamma"))

print(p)

ggsave(p, filename=sprintf('%s/%s/%s.pdf', wd, path_figs, 'kernel_discrete_cdfs'), width=12, height=12)


# phi plot
dat = data.frame(matrix(0, nrow=0, ncol=3))

for (run in runs){
  fname = sprintf('%s/%s.csv', path_out, run$suff_fit)
  system(sprintf('r/fixup.pl %s', fname))
  fit   = read_stan_csv(fname)
  post  = rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)
  
  col_names = colnames(post[,1,])
  par_names  = unlist(lapply(col_names, function(x) strsplit(x, "\\[")[[1]][1]))
  par_vals  = post[,1,which(par_names == 'phi')]
  
  par_mean  = colMeans(par_vals)
  par_lb    = apply(par_vals, 2, function(x) quantile(x, probs=0.025))
  par_ub    = apply(par_vals, 2, function(x) quantile(x, probs=0.975))
  
  par_stats = data.frame(name=taxa, mu=par_mean, lb=par_lb, ub=par_ub, Model=rep(run$handle, length(taxa)))
  
  dat = rbind(dat, par_stats)
}


greys = brewer.pal(9, "Greys")[c(4,5,7,9)]

p <- ggplot(data=dat, aes(x=reorder(name, mu), y=mu, group=Model, colour=Model)) + 
  geom_point(size=4, position=position_dodge(width=0.5)) + 
  geom_errorbar(aes(ymin=lb, ymax=ub), width=.2, position=position_dodge(width=0.5)) + 
  scale_colour_manual(values=greys) +
  xlab("Taxon") + ylab(parse(text='phi')) +
  coord_flip() + theme_bw() + theme(axis.title.x=element_text(size=20), 
                                    axis.title.y=element_text(size=20), 
                                    axis.text.x=element_text(size=rel(1.3)),
                                    axis.text.y=element_text(size=rel(1.3)))

print(p)

ggsave(p, filename=sprintf('%s/%s/%s.pdf', wd, path_figs, 'phi'), width=12, height=10)

# gamma plot
dat = data.frame(matrix(0, nrow=0, ncol=3))

for (run in runs){
  fname = sprintf('%s/%s.csv', path_out, run$suff_fit)
  system(sprintf('r/fixup.pl %s', fname))
  fit   = read_stan_csv(fname)
  post  = rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)
  
  col_names = colnames(post[,1,])
  par_names = unlist(lapply(col_names, function(x) strsplit(x, "\\[")[[1]][1]))
  par_idx   = which(par_names == 'gamma')
  par_vals  = post[,1,par_idx]
  
  if (length(par_idx)>1){
    par_mean  = colMeans(par_vals)
    par_lb    = apply(par_vals, 2, function(x) quantile(x, probs=0.025))
    par_ub    = apply(par_vals, 2, function(x) quantile(x, probs=0.975))
  } else {
    par_mean  = rep(mean(par_vals), K)
    par_lb    = rep(quantile(par_vals, probs=0.025), K)
    par_ub    = rep(quantile(par_vals, probs=0.975), K)
  }
  
  
  par_stats = data.frame(name=taxa, mu=par_mean, lb=par_lb, ub=par_ub, handle=rep(run$handle, length(taxa)))
  
  dat = rbind(dat, par_stats)
}

p <- ggplot(data=dat, aes(x=reorder(name, mu), y=mu, group=handle, colour=handle)) + 
  geom_point(size=4, position=position_dodge(width=0.5)) + 
  geom_errorbar(aes(ymin=lb, ymax=ub), width=.2, position=position_dodge(width=0.5)) + 
  xlab("Taxon") + ylab(parse(text='gamma')) +
  coord_flip() + theme_bw() + theme(axis.title.x=element_text(size=20), 
                                    axis.title.y=element_text(size=20), 
                                    axis.text.x=element_text(size=rel(1.3)),
                                    axis.text.y=element_text(size=rel(1.3)))

ggsave(p, filename=sprintf('%s/%s/%s.pdf', wd, path_figs, 'gamma'), width=12, height=12)