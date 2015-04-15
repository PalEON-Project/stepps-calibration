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
suff_fit_pl  = '12taxa_mid_comp_v0.1_pl_bigC'
suff_fit_g = paste0('12taxa_mid_comp_v0.1_bigC_c', seq(1,3))

post_pl = read_stan_output(path_out, suff_fit_pl)
post_g = read_stan_output(path_out, suff_fit_g)

limits <- get_limits(centers_veg)

#####################################################################################
# trace plots
#####################################################################################
trace_plots(post_g, suff, save_plots=save_plots, fpath='figures')

#####################################################################################
# read in data and source utils
#####################################################################################

psi = 0.207299
a   = 0.01446661
b   = 2.020447

dx = 0.008
dr = 0.001

dvec = seq(0, 1, dr)

# power law

C_pl=sum(d_pot[,2]*power_law(d_pot[,1], a, b))*dx^2
C_pl

# gaussian
C_g=sum(d_pot[,2]*gaussian(d_pot[,1], psi))*dx^2
C_g

p_g = gaussian(dvec, psi)
p_pl = power_law(dvec, a, b)

pdf('figures/kernel_pdfs.pdf')
plot(dvec*1e3, p_pl*dvec/C_pl, type='l', xlab='Radius', ylab='Density')
lines(dvec*1e3, p_g*dvec/C_g, lty=2, col='blue')
legend('topright', c('Power law', 'Gaussian'), col=c('black', 'blue'), lty=c(1,2))
dev.off()

# cdf
c_pl = cumsum(p_pl*dvec/C_pl*2*pi)*dr
c_g = cumsum(p_g*dvec/C_g*2*pi)*dr

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

C_g <- compute_C(post_g, N_pot, d_pot, kernel='gaussian')
r_int_g  = dispersal_decay(post_g, dmat, C_g, radius/rescale, kernel='gaussian')
C_pl <- compute_C(post_pl, N_pot, d_pot, kernel='pl')
r_int_pl = dispersal_decay(post_pl, dmat, C_pl, radius/rescale, kernel='pl')

fifty  = which.min(abs(r_int - 0.5))
ninety = which.min(abs(r_int - 0.9))

# segments = data.frame(x=c(0, 0, radius[fifty]/1e3, radius[ninety]/1e3), 
#                       xend=c(radius[fifty]/1e3, radius[ninety]/1e3, radius[fifty]/1e3, radius[ninety]/1e3),
#                       y=c(r_int[fifty], r_int[ninety], 0.2, 0.2),
#                       yend=c(r_int[fifty], r_int[ninety], r_int[fifty], r_int[ninety]))

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

