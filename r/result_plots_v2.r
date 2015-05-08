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

source(file.path(path_utils, 'process_funs.r'))
source(file.path(path_utils, 'plot_funs.r'))

load(sprintf('%s/cal_data_%s.rdata', path_data, suff_dat))

path_out    = 'output'

run1 = list(suff_fit  = 'cal_g_v0.3', 
            kernel    = 'gaussian', 
            one_psi   = TRUE, 
            one_gamma = TRUE, 
            EPs       = FALSE)
run2 = list(suff_fit  = 'cal_g_Kpsi_EPs_v0.3', 
            kernel    = 'gaussian', 
            one_psi   = FALSE, 
            one_gamma = TRUE, 
            EPs       = TRUE)
run3 = list(suff_fit  = 'cal_g_Kpsi_Kgamma_EPs_v0.3', 
            kernel    = 'gaussian', 
            one_psi   = FALSE, 
            one_gamma = FALSE, 
            EPs       = TRUE)
run4 = list(suff_fit  = 'cal_g_Kgamma_EPs_v0.3', 
            kernel    = 'gaussian', 
            one_psi   = TRUE, 
            one_gamma = FALSE, 
            EPs       = TRUE)
run5 = list(suff_fit  = 'cal_pl_v0.3', 
            kernel    = 'pl', 
            one_a     = TRUE,
            one_b     = TRUE,
            one_gamma = TRUE, 
            EPs       = FALSE)
run6 = list(suff_fit  = 'cal_pl_Kgamma_EPs_v0.3', 
            kernel    = 'pl', 
            one_a     = TRUE,
            one_b     = TRUE,
            one_gamma = FALSE, 
            EPs       = TRUE)

# runs = list(run1, run2, run3, run4, run5, run6)
runs = list(run1, run5)

fname = sprintf('%s/%s.csv', path_out, run5$suff_fit)
fit_pl  = read_stan_csv(fname)
post_pl = rstan::extract(fit_pl, permuted=FALSE, inc_warmup=FALSE)

fname = sprintf('%s/%s.csv', path_out, run1$suff_fit)
fit_g  = read_stan_csv(fname)
post_g = rstan::extract(fit_g, permuted=FALSE, inc_warmup=FALSE)

limits <- get_limits(centers_veg)

# #####################################################################################
# # trace plots
# #####################################################################################
# trace_plots(post_g, suff, save_plots=save_plots, fpath='figures')

#####################################################################################
# read in data and source utils
#####################################################################################

# power law
sum_w_pl = build_sumw_pot(post_pl, K, N_pot, d_pot, run5)
sum_w_pl

# gaussian
sum_w_g = build_sumw_pot(post_g, K, N_pot, d_pot, run1)
sum_w_g

col_substr_g = substr(colnames(post_g[,1,]), 1, 3)
col_substr_pl = substr(colnames(post_pl[,1,]), 1, 3)

one_psi = run1$one_psi
if (one_psi){
  psi   = rep(mean(post_g[,1,which(col_substr_g == 'psi')]), K)
} else {
  psi   = colMeans(post_g[,1,which(col_substr_g == 'psi')])
}


one_a = run5$one_a
if (one_a){
  a = rep(mean(post_pl[,1,which(col_substr_pl == 'a')]), K)
} else {
  a = colMeans(post_pl[,1,which(col_substr_pl == 'a')])
}
one_b = run5$one_b
if (one_b){
  b = rep(mean(post_pl[,1,which(col_substr_pl == 'b')]), K)
} else {
  b = colMeans(post_pl[,1,which(col_substr_pl == 'b')])
}


dx = 0.008
dr = 0.001

dvec = seq(0, 1, dr)

# power law

C_pl=build_sumw_pot(post_pl, K, N_pot, d_pot, run5)*dx^2
C_pl

# gaussian
C_g=build_sumw_pot(post_g, K, N_pot, d_pot, run1)*dx^2
C_g

p_g  = gaussian(dvec, psi[1])
p_pl = power_law(dvec, a[1], b[1])

pdf('figures/kernel_pdfs.pdf')
plot(dvec*1e3, p_pl*dvec/C_pl[1], type='l', xlab='Radius', ylab='Density')
lines(dvec*1e3, p_g*dvec/C_g[1], lty=2, col='blue')
legend('topright', c('Power law', 'Gaussian'), col=c('black', 'blue'), lty=c(1,2))
dev.off()

# cdf
c_pl = cumsum(p_pl*dvec/C_pl[1]*2*pi)*dr
c_g = cumsum(p_g*dvec/C_g[1]*2*pi)*dr

pdf('figures/kernel_cdfs.pdf')
plot(dvec*1e3, c_pl, type='l', xlab='Radius', ylab='Estimated cumulative density')
lines(dvec*1e3, c_g, col='blue', lty=2)
legend('bottomright', c('Power law', 'Gaussian'), col=c('black', 'blue'), lty=c(1,2))
dev.off()


sum(p_pl*dvec*2*pi)*dr/C_pl
sum(p_g*dvec*2*pi)*dr/C_g


############################################################################################################
# dispersal cdf
############################################################################################################

radius = seq(8000,1000000, by=4000)

x_pot = seq(-528000, 528000, by=8000)
y_pot = seq(-416000, 416000, by=8000)
coord_pot = expand.grid(x_pot, y_pot)

dmat = t(rdist(matrix(c(0,0), ncol=2), as.matrix(coord_pot, ncol=2))/rescale)

# power law
sum_w_pl = build_sumw_pot(post_pl, K, N_pot, d_pot, run5)
r_int_pl = dispersal_decay(post_pl, dmat, sum_w_pl[1], radius/rescale, kernel='pl')
# gaussian
sum_w_g = build_sumw_pot(post_g, K, N_pot, d_pot, run1)
r_int_g = dispersal_decay(post_g, dmat, sum_w_g[1], radius/rescale, kernel='gaussian')

fifty  = which.min(abs(r_int - 0.5))
ninety = which.min(abs(r_int - 0.9))

# segments = data.frame(x=c(0, 0, radius[fifty]/1e3, radius[ninety]/1e3), 
#                       xend=c(radius[fifty]/1e3, radius[ninety]/1e3, radius[fifty]/1e3, radius[ninety]/1e3),
#                       y=c(r_int[fifty], r_int[ninety], 0.2, 0.2),
#                       yend=c(r_int[fifty], r_int[ninety], r_int[fifty], r_int[ninety]))
library(reshape)
dat = data.frame(radius=radius/1e3, Gaussian=r_int_g, InversePowerLaw=r_int_pl)
dat = melt(dat, id='radius')

p <- ggplot(dat) + geom_line(aes(x=radius, y=value, color=factor(variable), linetype=variable))
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab('Radius') + ylab('Proportion of pollen')
p <- p + theme(axis.text= element_text(size=rel(1)), axis.title=element_text(size=14))
# p <- p + geom_segment(data=segments, aes(x=x, y=y, xend=xend, yend=yend), linetype=2, colour='royalblue4')
p <- p + xlim(0, 600) + ylim(0,1.0)
p <- p + labs(color='Kernel', linetype='Kernel')
p

ggsave(file='figures/dispersal_vs_distance.pdf', scale=1)
ggsave(file='figures/dispersal_vs_distance.eps', scale=1)

