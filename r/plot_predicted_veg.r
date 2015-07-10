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

cal_run  = 'cal_pl_Ka_Kgamma_EPs_ALL_v0.3'
pred_runs = list('12taxa_6341cells_120knots_cal_KW_KGAMMA_PL_umw_1by',
                 '12taxa_6341cells_120knots_cal_KW_KGAMMA_G_umw_1by',
                 '12taxa_6341cells_120knots_cal_PL_umw_1by',
                 '12taxa_6341cells_120knots_cal_G_umw_1by')

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

###############################################################################################################
# read in the data, organize in a melted DF
###############################################################################################################
preds = list()
for (pred_run in pred_runs){

# load composition prediction; r_mean and centers_r_mean
load(paste0('../stepps-prediction/figures/', pred_run ,'/r_mean.rdata'))
centers_r_mean = centers_r_mean*rescale
colnames(centers_r_mean) = c('x', 'y')
colnames(r_mean) = taxa

preds = rbind(preds, cbind(r_mean, centers_r_mean, type=rep(pred_run, nrow(r_mean))))

}

preds=melt(preds, id.vars=c('x', 'y', 'type'))

###############################################################################################################
# fix names for plotting
###############################################################################################################
levels(preds$variable)[levels(preds$variable) == 'OTHER.CONIFER'] = 'OTHER CON'
levels(preds$variable)[levels(preds$variable) == 'OTHER.HARDWOOD'] = 'OTHER HW'
levels(preds$type) <- c('Variable PL', 'Variable G', 'Base PL', 'Base G')
preds$type <- factor(preds$type, levels=c('Base G', 'Base PL', 'Variable G', 'Variable PL'))

limits <- get_limits(cbind(preds$x, preds$y))


# 
# # continuous
# p <- ggplot() + geom_tile(data=preds_melt, aes(x=x, y=y, fill=value)) + 
#   scale_fill_gradientn(colours=tim.colors(), name='Proportions', limits=c(0,1)) + 
#   coord_fixed() 
# p <- add_map_albers(p, us.shp, limits, rescale)
# p <- p + facet_grid(variable~type)
# p <- theme_clean(p) 
# p <- p + theme(strip.background = element_blank(), strip.text.x = element_blank())#+ theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.5)))
# # p <- p + theme(strip.background = element_blank())
# print(p)

# discrete binned
breaks = c(0, 0.01, 0.05, 0.10, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 1)
breaklabels = apply(cbind(breaks[1:(length(breaks)-1)], breaks[2:length(breaks)]), 1, 
                    function(r) { sprintf("%0.2f - %0.2f", r[1], r[2]) })

preds$value_binned = cut(preds$value, breaks, include.lowest=TRUE, labels=FALSE)

cols = tim.colors(length(breaks))

p <- ggplot() + geom_tile(data=preds, aes(x=x, y=y, fill=factor(value_binned))) + 
  scale_fill_manual(values = cols, labels=breaklabels, name='Proportions') + 
  coord_fixed() 
p <- add_map_albers(p, us.shp, limits, rescale)
p <- p + facet_grid(variable~type)
p <- theme_clean(p) 
# p <- p + theme(strip.background = element_blank(), strip.text.x = element_blank())
p <- p + theme(strip.text.y = element_text(size = rel(1.5)), 
               strip.text.x = element_text(size = rel(1.5)))
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p)

fname = 'figures/maps_predicted_veg.pdf'
ggsave(p, file=fname, scale=3)
sys_str = paste('pdfcrop', fname, fname, sep=' ')
system(sys_str)

# ggsave(p, file='figures/comp_grid_plot.eps', scale=3)#, width=14, units='in')
# ggsave(p, file='figures/comp_grid_plot.png', scale=3)