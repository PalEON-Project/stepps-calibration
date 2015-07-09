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

taxa_abb = taxa
taxa_abb[which(taxa_abb == 'OTHER.CONIFER')] = 'OTHER CON'
taxa_abb[which(taxa_abb == 'OTHER.HARDWOOD')] = 'OTHER HW'

###############################################################################################################
# plot the modelled composition
###############################################################################################################
limits <- get_limits(centers_veg)

breaks = c(0, 0.01, 0.05, 0.10, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 1)
# plot_data_maps_binned(r_mean, centers=centers_pls, taxa=taxa, K, T, breaks, suff=suff4, save_plots=save_plots)

plot_data_maps_binned(r, centers_veg, taxa_abb, K, breaks, limits, suff='veg', save_plots, fpath=path_figs)

###############################################################################################################
###############################################################################################################

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

########################################################################################################################################
# dispersal cdf
########################################################################################################################################
radius = seq(8000,1000000, by=4000)

coord_pot = seq(-700000, 700000, by=8000)
coord_pot = expand.grid(coord_pot, coord_pot)

d_pot = t(rdist(matrix(c(0,0), ncol=2), as.matrix(coord_pot, ncol=2))/rescale)

# want a circular region
idx_circ  = which(d_pot[,1] > 0.7)
coord_pot = coord_pot[idx_circ, ]
d_pot     = d_pot[idx_circ, ]

d_pot = unname(as.matrix(count(d_pot)))
N_pot = nrow(d_pot)

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

handles = as.vector(unique(dat$handle))

capture_50 = data.frame(handle=character(0), taxon=character(0), r=numeric(0))
for (handle in handles){
  for (taxon in taxa){
    dat_sub = dat[which((dat$handle == handle) & (dat$variable == taxon)), ]
    rad = dat_sub[which.min(abs(dat_sub[, 'value'] - 0.9)), 'radius']
    capture_50 = rbind(capture_50, data.frame(handle=handle, taxon=taxon, r=rad))
  }
}

capture_50_sub = capture_50[capture_50$handle == 'pl_Ka_Kgamma',]
# capture_50_sub = capture_50[capture_50$handle == 'g_Kpsi_Kgamma',]

# svs       = read.table("data/svs.csv", sep=',', header=TRUE)
# # ppe_compiled = data.frame(taxon = taxa, ppe_stepps = ppe_stepps$mu, ppe_m = ppes[,2], ppe_s = ppes[,3])
# svs_compiled = data.frame(svs, capture_50=capture_50_sub$r)
# svs_compiled = svs_compiled[!(svs_compiled$taxon %in% c('OTHER.HARDWOOD')), ]
# 
# svs_melt     = melt(svs_compiled, id.vars=c('taxon', 'capture_50'))
# 
# p <- ggplot() + geom_point(data=svs_melt, aes(x=value, y=capture_50, colour=taxon), size=5) + 
#   scale_colour_manual(values=brewer.pal(nrow(svs_compiled), 'Paired')) + xlab('SV from lit') + ylab('Capture radius STEPPS') + 
#   #   geom_text(data = ppe_melt, aes(x=value,y=ppe_stepps, label = taxon), position=position_jitter(h=1,w=1)) + #hjust = 0.8, vjust=-1.5) +
#   facet_wrap(~variable, nrow=2)
# print(p)
# 
# ggsave(p, filename=sprintf('%s/%s/%s.pdf', wd, path_figs, 'SVs'), width=12, height=10)
# 
# 

svs       = read.table("data/svs_meta.csv", sep=',', header=TRUE)
# svs = rbind(svs[,c(1,2,6) ], data.frame(taxon=capture_50_sub[,2], sv=capture_50_sub[,3], code=rep('STEPPS', nrow(capture_50_sub))))

svs = cbind(svs, data.frame(stepps=capture_50_sub[match(svs$taxon, capture_50_sub$taxon),3]))

svs_melt = svs[!(svs$taxon %in% c('OTHER.HARDWOOD', 'OTHER.CONIFER')), ]
# svs_melt     = melt(svs_compiled, id.vars=c('taxon', 'capture_50'))

# svs_melt[svs_melt$taxon == 'MAPLE', 'stepps'] = svs_melt[svs_melt$taxon == 'MAPLE', 'stepps'] + 2
# svs_melt[svs_melt$taxon == 'SPRUCE', 'stepps'] = svs_melt[svs_melt$taxon == 'SPRUCE', 'stepps'] -2

levels(svs_melt$code)= c('Matthias (2014)', 'Jackson (1999)' )

p <- ggplot() + geom_point(data=svs_melt, aes(x=sv, y=stepps, colour=taxon), alpha=0.9, size=5) + 
  scale_colour_manual(values=brewer.pal(length(unique(svs_melt$taxon)), 'Paired')) + 
  xlab('SV from lit') + ylab('Capture radius STEPPS') +
  #   geom_text(data = ppe_melt, aes(x=value,y=ppe_stepps, label = taxon), position=position_jitter(h=1,w=1)) + #hjust = 0.8, vjust=-1.5) +
  facet_wrap(~code, nrow=2)
print(p)

ggsave(p, filename=sprintf('%s/%s/%s.pdf', wd, path_figs, 'SVs_v2'), width=10, height=14)

########################################################################################################################################
# phi plot
########################################################################################################################################

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

# compare with PPEs from lit
ppe_stepps = dat[which(dat$Model == 'pl_Ka_Kgamma'),]
ppes       = read.table("data/ppes.csv", sep=',', header=TRUE, stringsAsFactors=FALSE)
if (!all(as.character(ppe_stepps[,1]) == ppes[,1])){
  print("Taxa mismatch or ordering problem!")
}

# ppe_compiled = data.frame(taxon = taxa, ppe_stepps = ppe_stepps$mu, ppe_m = ppes[,2], ppe_s = ppes[,3])
ppe_compiled = data.frame(ppes, ppe_stepps=ppe_stepps$mu)
ppe_compiled = ppe_compiled[!(ppe_compiled$taxon %in% c('OTHER.HARDWOOD', 'HEMLOCK', 'TAMARACK', 'OTHER.CONIFER')), ]

ppe_compiled[,2:ncol(ppe_compiled)] = apply(ppe_compiled[,2:ncol(ppe_compiled)], 2, function(x) x/max(x, na.rm=TRUE))
ppe_melt     = melt(ppe_compiled, id.vars=c('taxon', 'ppe_stepps'))

p <- ggplot() + geom_point(data=ppe_melt, aes(x=value, y=ppe_stepps, colour=taxon), size=5) + 
     scale_colour_manual(values=brewer.pal(nrow(ppe_compiled), 'Paired')) + xlab('PPE from lit') + ylab('PPE STEPPS') + 
#   geom_text(data = ppe_melt, aes(x=value,y=ppe_stepps, label = taxon), position=position_jitter(h=1,w=1)) + #hjust = 0.8, vjust=-1.5) +
  facet_wrap(~variable, nrow=2)
# print(p)

ggsave(p, filename=sprintf('%s/%s/%s.pdf', wd, path_figs, 'PPEs_panels'), width=12, height=10)

p <- ggplot() + geom_point(data=ppe_melt, aes(x=value, y=ppe_stepps, colour=taxon, shape=variable), size=5) + 
  scale_colour_manual(values=brewer.pal(nrow(ppe_compiled), 'Paired')) + xlab('PPE from lit') + ylab('PPE STEPPS') #+ 
  #   geom_text(data = ppe_melt, aes(x=value,y=ppe_stepps, label = taxon), position=position_jitter(h=1,w=1)) + #hjust = 0.8, vjust=-1.5) +
#   facet_wrap(~variable, nrow=2)
# print(p)
ggsave(p, filename=sprintf('%s/%s/%s.pdf', wd, path_figs, 'PPEs_single'), width=12, height=10)

ppe_melt2 = melt(ppe_compiled, id.vars=c('taxon'))
q <- ggplot() + geom_point(data=ppe_melt2, aes(x=taxon, y=value))
print(q)



# p <- ggplot() + geom_point(data=ppe_compiled, aes(x=ppe_s, y=ppe_stepps), size=3) + 
#      geom_text(data = ppe_compiled, aes(x=ppe_s,y=ppe_stepps, label = taxon), hjust = 0.8, vjust=-1.5)
# print(p)
# direct.label(p, first.qp)
# 
# plot(ppes$marquer, ppe_stepps$mu)

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