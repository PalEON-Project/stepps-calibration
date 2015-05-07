library(rstan)
library(ggplot2)
library(fields, quietly=TRUE)
library(RColorBrewer)
library(reshape2)

wd = getwd()

#####################################################################################
# user pars
#####################################################################################

path_utils = 'r/utils'
path_data  = 'r/dump'
path_out   = 'output'
path_figs  = 'figures'

suff_dat = '12taxa_mid_comp_ALL_v0.2'

save_plots = TRUE
rescale    = 1e6

#####################################################################################
# read in data and source utils
#####################################################################################

source(file.path(path_utils, 'process_funs.r'))
source(file.path(path_utils, 'plot_funs.r'))

load(sprintf('%s/cal_data_%s.rdata', path_data, suff_dat))

path_out    = 'output'

run1 = list(suff_fit  = 'cal_g_ALL_v0.3', 
            kernel    = 'gaussian', 
            one_psi   = TRUE, 
            one_gamma = TRUE, 
            EPs       = FALSE,
            handle    = 'G')
run2 = list(suff_fit  = 'cal_g_Kpsi_Kgamma_EPs_ALL_v0.3', 
            kernel    = 'gaussian', 
            one_psi   = FALSE, 
            one_gamma = FALSE, 
            EPs       = TRUE,
            handle    = 'G_Kpsi_Kgamma')
run3 = list(suff_fit  = 'cal_pl_ALL_v0.3', 
            kernel    = 'pl', 
            one_a     = TRUE,
            one_b     = TRUE,
            one_gamma = TRUE, 
            EPs       = FALSE,
            handle    = 'PL')
run4 = list(suff_fit  = 'cal_pl_Ka_Kgamma_EPs_ALL_v0.3', 
            kernel    = 'pl', 
            one_a     = FALSE,
            one_b     = TRUE,
            one_gamma = FALSE, 
            EPs       = TRUE,
            handle    = 'PL_Ka_Kgamma')

runs = list(run1, run2, run3, run4)
# runs = list(run4)
# 
# fname = sprintf('%s/%s.csv', path_out, run5$suff_fit)
# fit_pl  = read_stan_csv(fname)
# post_pl = rstan::extract(fit_pl, permuted=FALSE, inc_warmup=FALSE)
# 
# fname = sprintf('%s/%s.csv', path_out, run1$suff_fit)
# fit_g  = read_stan_csv(fname)
# post_g = rstan::extract(fit_g, permuted=FALSE, inc_warmup=FALSE)

############################################################################################################
# dispersal cdf
############################################################################################################

radius = seq(8000,1000000, by=4000)

x_pot = seq(-528000, 528000, by=8000)
y_pot = seq(-416000, 416000, by=8000)
coord_pot = expand.grid(x_pot, y_pot)

dmat = t(rdist(matrix(c(0,0), ncol=2), as.matrix(coord_pot, ncol=2))/rescale)

dat = data.frame(matrix(0, nrow=0, ncol=4))
for (run in runs){
  fname = sprintf('%s/%s.csv', path_out, run$suff_fit)
  system(sprintf('r/fixup.pl %s', fname))
  fit   = read_stan_csv(fname)
  post  = rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)
  
  sum_w = build_sumw_pot(post, K, N_pot, d_pot, run)
  r_int = dispersal_decay(post, dmat, sum_w, radius/rescale, run, taxa)
  
  run_dat = data.frame(radius = radius/1e3, r_int, handle=rep(run$handle, length(radius)))
  dat = rbind(dat, melt(run_dat, id=c('radius', 'handle')))
}

pol_fifty  = data.frame(radius=numeric(0), handle=character(0), variable=character(0), value=numeric(0))
pol_ninety = data.frame(radius=numeric(0), handle=character(0), variable=character(0), value=numeric(0))

models = as.vector(unique(dat$handle))
for (model in models){
  for (taxon in taxa){
    dat_sub    = dat[(dat$handle == model) & (dat$variable == taxon), ]
    pol_fifty  = rbind(pol_fifty, dat_sub[which.min(abs(dat_sub$value - 0.5)), ])
    pol_ninety = rbind(pol_ninety, dat_sub[which.min(abs(dat_sub$value - 0.9)), ])
  }
}


limits <- get_limits(centers_veg)

p <- ggplot(dat) + 
  geom_line(data=dat, aes(x=radius, y=value, linetype=factor(handle), colour=factor(handle)), size=0.9) 
p <- p + scale_colour_manual("Model", values=c('darkorange', 'darkorange', 'royalblue4', 'royalblue4')) 
p <- p + scale_linetype_manual("Model", values=c("solid", "dashed", "solid", "dashed")) + xlab("Radius") +
  ylab("Cumulative density") 
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
  
  par_stats = data.frame(name=taxa, mu=par_mean, lb=par_lb, ub=par_ub, handle=rep(run$handle, length(taxa)))
  
  dat = rbind(dat, par_stats)
}

p <- ggplot(data=dat, aes(x=reorder(name, mu), y=mu, group=handle, colour=handle)) + 
  geom_point(size=4, position=position_dodge(width=0.5)) + 
  geom_errorbar(aes(ymin=lb, ymax=ub), width=.2, position=position_dodge(width=0.5)) + 
  xlab("Taxon") + ylab(parse(text='phi')) +
  coord_flip() + theme_bw() + theme(axis.title.x=element_text(size=20), 
                                    axis.title.y=element_text(size=20), 
                                    axis.text.x=element_text(size=rel(1.3)),
                                    axis.text.y=element_text(size=rel(1.3)))

print(p)

ggsave(p, filename=sprintf('%s/%s/%s.pdf', wd, path_figs, 'phi'), width=12, height=12)

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

print(p)

ggsave(p, filename=sprintf('%s/%s/%s.pdf', wd, path_figs, 'gamma'), width=12, height=12)

# 
# p <- ggplot(data=dat, aes(x=handle, y=mu, group=handle, colour=handle)) + 
#   geom_point(size=4, position=position_dodge(width=1)) + 
#   geom_errorbar(aes(ymin=lb, ymax=ub), width=.2, position=position_dodge(width=1)) + 
#   xlab("Taxon") + ylab(parse(text='phi')) +
#   coord_flip() + facet_wrap(~name) + theme_bw() + theme(axis.title.x=element_text(size=20), 
#                                     axis.title.y=element_text(size=20), 
#                                     axis.text.x=element_text(size=rel(1.3)),
#                                     axis.text.y=element_text(size=rel(1.3)))
# 
# print(p)


# # power law
# sum_w_pl = build_sumw_pot(post_pl, K, N_pot, d_pot, run5)
# r_int_pl = dispersal_decay(post_pl, dmat, sum_w_pl[1], radius/rescale, kernel='pl')
# # gaussian
# sum_w_g = build_sumw_pot(post_g, K, N_pot, d_pot, run1)
# r_int_g = dispersal_decay(post_g, dmat, sum_w_g[1], radius/rescale, kernel='gaussian')
# 
# fifty  = which.min(abs(r_int - 0.5))
# ninety = which.min(abs(r_int - 0.9))
# 
# # segments = data.frame(x=c(0, 0, radius[fifty]/1e3, radius[ninety]/1e3), 
# #                       xend=c(radius[fifty]/1e3, radius[ninety]/1e3, radius[fifty]/1e3, radius[ninety]/1e3),
# #                       y=c(r_int[fifty], r_int[ninety], 0.2, 0.2),
# #                       yend=c(r_int[fifty], r_int[ninety], r_int[fifty], r_int[ninety]))
# library(reshape)
# dat = data.frame(radius=radius/1e3, Gaussian=r_int_g, InversePowerLaw=r_int_pl)
# dat = melt(dat, id='radius')
# 
# p <- ggplot(dat) + geom_line(aes(x=radius, y=value, color=factor(variable), linetype=variable))
# p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   xlab('Radius') + ylab('Proportion of pollen')
# p <- p + theme(axis.text= element_text(size=rel(1)), axis.title=element_text(size=14))
# # p <- p + geom_segment(data=segments, aes(x=x, y=y, xend=xend, yend=yend), linetype=2, colour='royalblue4')
# p <- p + xlim(0, 600) + ylim(0,1.0)
# p <- p + labs(color='Kernel', linetype='Kernel')
# p
# 
# ggsave(file='figures/dispersal_vs_distance.pdf', scale=1)
# ggsave(file='figures/dispersal_vs_distance.eps', scale=1)
# 


# #####################################################################################
# # read in data and source utils
# #####################################################################################
# 
# # power law
# sum_w_pl = build_sumw_pot(post_pl, K, N_pot, d_pot, run5)
# sum_w_pl
# 
# # gaussian
# sum_w_g = build_sumw_pot(post_g, K, N_pot, d_pot, run1)
# sum_w_g
# 
# col_substr_g = substr(colnames(post_g[,1,]), 1, 3)
# col_substr_pl = substr(colnames(post_pl[,1,]), 1, 3)
# 
# one_psi = run1$one_psi
# if (one_psi){
#   psi   = rep(mean(post_g[,1,which(col_substr_g == 'psi')]), K)
# } else {
#   psi   = colMeans(post_g[,1,which(col_substr_g == 'psi')])
# }
# 
# 
# one_a = run5$one_a
# if (one_a){
#   a = rep(mean(post_pl[,1,which(col_substr_pl == 'a')]), K)
# } else {
#   a = colMeans(post_pl[,1,which(col_substr_pl == 'a')])
# }
# one_b = run5$one_b
# if (one_b){
#   b = rep(mean(post_pl[,1,which(col_substr_pl == 'b')]), K)
# } else {
#   b = colMeans(post_pl[,1,which(col_substr_pl == 'b')])
# }
# 
# 
# dx = 0.008
# dr = 0.001
# 
# dvec = seq(0, 1, dr)
# 
# # power law
# 
# C_pl=build_sumw_pot(post_pl, K, N_pot, d_pot, run5)*dx^2
# C_pl
# 
# # gaussian
# C_g=build_sumw_pot(post_g, K, N_pot, d_pot, run1)*dx^2
# C_g
# 
# p_g  = gaussian(dvec, psi[1])
# p_pl = power_law(dvec, a[1], b[1])
# 
# pdf('figures/kernel_pdfs.pdf')
# plot(dvec*1e3, p_pl*dvec/C_pl[1], type='l', xlab='Radius', ylab='Density')
# lines(dvec*1e3, p_g*dvec/C_g[1], lty=2, col='blue')
# legend('topright', c('Power law', 'Gaussian'), col=c('black', 'blue'), lty=c(1,2))
# dev.off()
# 
# # cdf
# c_pl = cumsum(p_pl*dvec/C_pl[1]*2*pi)*dr
# c_g = cumsum(p_g*dvec/C_g[1]*2*pi)*dr
# 
# pdf('figures/kernel_cdfs.pdf')
# plot(dvec*1e3, c_pl, type='l', xlab='Radius', ylab='Estimated cumulative density')
# lines(dvec*1e3, c_g, col='blue', lty=2)
# legend('bottomright', c('Power law', 'Gaussian'), col=c('black', 'blue'), lty=c(1,2))
# dev.off()
# 
# 
# sum(p_pl*dvec*2*pi)*dr/C_pl
# sum(p_g*dvec*2*pi)*dr/C_g