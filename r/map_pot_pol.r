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

cal_runs  = list(list(run = 'cal_pl_Ka_Kgamma_EPs_ALL_v0.3', 
                      tag = 'Variable PL'),
                 list(run = 'cal_g_Kpsi_Kgamma_EPs_ALL_v0.3',
                      tag = 'Variable G'))

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
colnames(r)           = taxa

dat = data.frame(r, centers_veg, tag=rep('PLS data', nrow(r)))
for (cal_run in cal_runs){
  # load potential pollen predictions; pp and centers_pp
  load(paste0('figures/', cal_run$run, '/pp_all.rdata'))
  colnames(centers_pp) = c('x', 'y')
  colnames(pp) = taxa
  
  if (!(all(centers_pp == centers_veg))){
    # check ordering
    dist_mat = rdist(as.matrix(centers_veg), as.matrix(centers_pp))
    dist_mat[dist_mat < 1e-6] = 0
    idx_in = unlist(apply(dist_mat, 1, function(x) if (any(x == 0)){ which(x==0)} else {0}))
    pp     = pp[idx_in, ]
    centers_pp = centers_pp[idx_in,]
  }
  
  dat = rbind(dat, data.frame(pp, centers_pp, tag=rep(cal_run$tag, nrow(pp))))
  
}

dat = melt(dat, id.vars=c('x', 'y', 'tag'))

###############################################################################################################
# plot the modelled composition
###############################################################################################################
limits <- get_limits(centers_veg)

levels(dat$tag) = c("PLS composition \n estimates", "Predicted pollen \n Variable PL", "Predicted pollen \n Variable G")   

# rename and reorder some factors
levels(dat$variable)[levels(dat$variable) == 'OTHER.CONIFER'] = 'OTHER CON'
levels(dat$variable)[levels(dat$variable) == 'OTHER.HARDWOOD'] = 'OTHER HW'
dat$tag <- factor(dat$tag, levels=levels(dat$tag)[c(1,3,2)])

# discrete binned
breaks = c(0, 0.01, 0.05, 0.10, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 1)
breaklabels = apply(cbind(breaks[1:(length(breaks)-1)], breaks[2:length(breaks)]), 1, 
                    function(r) { sprintf("%0.2f - %0.2f", r[1], r[2]) })

dat$value_binned = cut(dat$value, breaks, include.lowest=TRUE, labels=FALSE)

pp_taxa = c('BIRCH', 'BEECH', 'MAPLE', 'OAK', 'PINE')
dat_sub = dat[dat$variable %in% pp_taxa,]

# continuous
p <- ggplot() + geom_tile(data=dat, aes(x=x, y=y, fill=value)) + 
  scale_fill_gradientn(colours=tim.colors(), name='Proportions', limits=c(0,1)) + 
  coord_fixed() 
p <- add_map_albers(p, us.shp, limits, rescale)
p <- p + facet_grid(variable~tag)
p <- theme_clean(p) 
p <- p + theme(strip.background = element_blank(), strip.text.x = element_blank())#+ theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.5)))
# p <- p + theme(strip.background = element_blank())
print(p)



cols = tim.colors(length(breaks))

p <- ggplot() + geom_tile(data=dat_sub, aes(x=x, y=y, fill=factor(value_binned))) + 
#   geom_tile(data=dat, aes(x=x, y=y, colour=factor(value_binned)), fill=NA) +
  scale_fill_manual(values = cols, labels=breaklabels, name='Proportions') + 
     coord_fixed() 
p <- add_map_albers(p, us.shp, limits, rescale)
p <- p + facet_grid(variable~tag)
p <- theme_clean(p) 
# p <- p + theme(strip.background = element_blank(), strip.text.x = element_blank())
p <- p + theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.5)))
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p)

fname = 'figures/map_pot_pol.pdf'
ggsave(p, file=fname, width=14)
# ggsave(p, file='figures/map_grid_plot.eps', scale=3)#, width=14, units='in')
# ggsave(p, file='figures/map_grid_plot.png', scale=3)
sys_str=paste('pdfcrop', fname, fname, sep=' ')
system(sys_str)

################################################################################################################

N_pol = nrow(centers_polA)

if (is.null(taxa)){taxa=seq(1,K)}

# discrete binned
breaks = c(0, 0.01, 0.05, 0.10, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 1)
breaklabels = apply(cbind(breaks[1:(length(breaks)-1)], breaks[2:length(breaks)]), 1, 
                    function(r) { sprintf("%0.2f - %0.2f", r[1], r[2]) })

colnames(y) = taxa
pol_dat = cbind(centers_polA, compute_props(y, taxa))
pol_dat = melt(pol_dat, id.vars=c('x', 'y'))
pol_dat = cbind(pol_dat, tag=rep('Pollen data', nrow(pol_dat)))[,c(1,2,5,3,4)]

all_dat = rbind(dat, pol_dat)

all_dat$tag = factor(all_dat$tag, levels = unique(all_dat$tag)[c(1,4,3,2)]) 

levels(all_dat$tag)


# discrete binned
breaks = c(0, 0.01, 0.05, 0.10, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 1)
breaklabels = apply(cbind(breaks[1:(length(breaks)-1)], breaks[2:length(breaks)]), 1, 
                    function(r) { sprintf("%0.2f - %0.2f", r[1], r[2]) })

all_dat$value_binned = cut(all_dat$value, breaks, include.lowest=TRUE, labels=FALSE)

pp_taxa = c('BIRCH', 'BEECH', 'MAPLE', 'OAK', 'PINE')
dat_sub = all_dat[all_dat$variable %in% pp_taxa,]

cols = tim.colors(length(breaks))

p <- ggplot() + geom_tile(data=subset(dat_sub, tag %in% c('PLS data', 'Variable PL', 'Variable G')), aes(x=x, y=y, fill=factor(value_binned))) + 
  #   geom_tile(data=dat, aes(x=x, y=y, colour=factor(value_binned)), fill=NA) +
  geom_point(data=subset(dat_sub, tag %in% c('Pollen data')), aes(x=x, y=y, colour=factor(value_binned)), size=2) +
  scale_fill_manual(values = cols, labels=breaklabels, name='Proportion') + 
  scale_colour_manual(values = cols, labels=breaklabels, name='Proportion') + 
  coord_fixed() 
p <- add_map_albers(p, us.shp, limits, rescale)
p <- p + facet_grid(variable~tag)
p <- theme_clean(p) 
# p <- p + theme(strip.background = element_blank(), strip.text.x = element_blank())
p <- p + theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.5)))
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p)

fname = 'figures/map_pls_pollen.pdf'
ggsave(p, file=fname, width=14)
# ggsave(p, file='figures/map_grid_plot.eps', scale=3)#, width=14, units='in')
# ggsave(p, file='figures/map_grid_plot.png', scale=3)
sys_str=paste('pdfcrop', fname, fname, sep=' ')
system(sys_str)