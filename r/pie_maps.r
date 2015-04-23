# library(sp)
# library(rgdal)
library(maptools)

source('r/utils/dataPlotFuns.r')
source('r/utils/simDataFuns.r')

# load pls and pollen data
# pls      = read.table(file='data/pls_UMW.csv', sep=",", row.names=NULL, header=TRUE)
# pollen   = read.table('data/pollen_se_sum.csv', header=TRUE)
load('r/dump/cal_data_12taxa_mid_comp_UMW_v0.2.rdata')

# load('r/dump/cal_data_12taxa_mid_mi_sp.rdata')

v = '0_3'
# v = 'mi_sp'

########################################################################################################
## pls pie map
########################################################################################################
# reorder data by pollen abundance
new.order = order(colSums(y), decreasing=TRUE)

taxa = taxa[new.order]
y = y[,new.order]

r = r[,new.order]
colnames(r) = taxa

########################################################################################################
## pls pie map
########################################################################################################

ntaxa = length(taxa)

centers   = centers_veg
colnames(centers) = c('x', 'y')

xlo = min(centers[,1])
xhi = max(centers[,1])
ylo = min(centers[,2])
yhi = max(centers[,2])

subgrid = regular_subgrid(centers, dx=34000,dy=34000)
knots_in = knots_in_domain4(subgrid, centers, cell_width = 8000)

plot(centers[,1], centers[,2], asp=1)
points(knots_in[,1], knots_in[,2], col='red')

d = rdist(knots_in, centers)

pls_coarse = matrix(0, nrow=nrow(d), ncol=ntaxa)
colnames(pls_coarse) = taxa

for (i in 1: nrow(centers)){
  #print(min(d[,i]))
  close_knot = which.min(d[,i])
  
  #   print(as.double(pls_cut[i,]))
  #   print(pls_cut[i,])
  pls_coarse[close_knot,] = pls_coarse[close_knot,] + r[i,]
}

pls_props  = t(apply(pls_coarse, 1, function(x) if (sum(x) != 0){x/sum(x)} else {x}))

shift=30000

# postscript('r/data/figs/pie_plot_pls_UMW_v0.2.eps', width=8, height=6)
pdf(paste0('figures/pie_plot_pls_UMW_v', v, '.pdf'), width=12, height=10)
par(mfrow=c(1,1))
pieMap(proportions = pls_props, 
       centers  = knots_in,
       restrict = FALSE,
       inputRestricted = FALSE,
       xlim   = c(xlo+shift, xhi-shift),
       ylim   = c(ylo+shift, yhi-shift),
       radius = 14000,
       scale  = 1,
       xlab   = 'x',
       ylab   = 'y',
       add_legend = FALSE, 
       main_title='')
dev.off()


########################################################################################################
## pollen pie map
########################################################################################################

pollen = y
pollen_props  = t(apply(pollen, 1, function(x) x/sum(x)))
colnames(pollen_props) = taxa

centers   = centers_polA
colnames(centers) = c('x', 'y')

# centers = data.frame(centers)
# coordinates(centers) <- ~ x + y
# proj4string(centers) <- CRS('+proj=longlat +ellps=WGS84')

# centersA <- spTransform(centers, CRS('+init=epsg:3175'))
# centersA <- as.matrix(data.frame(centersA))

# props_cut = cbind(pollen_props[,1:11], rowSums(pollen_props[,12:ncol(pollen_props)]))
# colnames(props_cut)[ncol(props_cut)] = 'other'

# xlo = min(centersA[,1])
# xhi = max(centersA[,1])
# ylo = min(centersA[,2])
# yhi = max(centersA[,2])

#postscript('r/data/figs/pie_plot_pollen_UMW.eps', width=8, height=6)
pdf(paste0('figures/pie_plot_pollen_UMW_v', v, '.pdf'), width=12, height=10)
par(mfrow=c(1,1))
pieMap(proportions = pollen_props, 
       centers  = centers,
       restrict = FALSE,
       inputRestricted = FALSE,
       xlim   = c(xlo+shift, xhi-shift),
       ylim   = c(ylo+shift, yhi-shift),
       radius = 18000,
       scale  = 1,
       xlab   = 'x',
       ylab   = 'y', 
       add_legend=TRUE,
       main_title='')
dev.off()




# #compute proportions
# pls_props    = colSums(pls[,4:ncol(pls)], na.rm=TRUE)/sum(colSums(pls[,4:ncol(pls)], na.rm=TRUE))
# pollen_props = colSums(pollen[,7:ncol(pollen)], na.rm=TRUE)/sum(colSums(pollen[,7:ncol(pollen)], na.rm=TRUE))
# 
# #
# # make the plot!
# #
# 
# # condense to 12 taxa
# n = 12
# pls_props_cond = c(pls_props[1:(n-1)], sum(pls_props[n:length(pls_props)]))
# names(pls_props_cond)[n] = 'other'
# 
# pollen_props_cond = c(pollen_props[1:(n-1)], sum(pollen_props[n:length(pollen_props)]))
# names(pollen_props_cond)[n] = 'other'
# 
# col_list = c("#1F78B4", "#33A02C", "#FF7F00", "#E31A1C", "#6A3D9A", "#B15928", 
#              "#FFFF99", "#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F", "#CAB2D6")




