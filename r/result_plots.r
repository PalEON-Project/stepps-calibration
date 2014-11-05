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

limits <- get_limits(centers_veg)

# let's do this!

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
plot(dvec, p_pl*dvec/C_pl, type='l', xlab='Radius', ylab='Density')
lines(dvec, p_g*dvec/C_g, col='blue')
legend('topright', c('Power law', 'Gaussian'), col=c('black', 'blue'), lty=c(1,1))
dev.off()

# cdf
c_pl = cumsum(p_pl*dvec/C_pl*2*pi)*dr
c_g = cumsum(p_g*dvec/C_g*2*pi)*dr

pdf('figures/kernel_cdfs.pdf')
plot(dvec, c_pl, type='l', xlab='Radius', ylab='Estimated cumulative density')
lines(dvec, c_g, col='blue')
legend('bottomright', c('Power law', 'Gaussian'), col=c('black', 'blue'), lty=c(1,1))
dev.off()


sum(p_pl*dvec*2*pi)*dr/C_pl
sum(p_g*dvec*2*pi)*dr/C_g
