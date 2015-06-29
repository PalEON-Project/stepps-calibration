library(sp)
library(rgdal)
library(fields)
library(maptools)
gpclibPermit()

source('r/utils/build_data_funs.r')

#####################################################################################
# read in data
#####################################################################################

# set your working directory!
# AD: wd = '~/Documents/projects/stepps-calibration'

path_veg_old = 'data/composition_v0.2.csv'                                        # composition model results
path_veg     = 'data/composition_v0.3.csv'                                        # composition model results
path_convert = 'data/dict-comp2stepps.csv'                                        # conversion table

# append to output filenames
suff       = paste0(depth_type, '_ALL_v0.3')                                

states = c('michigan:north', 'minnesota', 'wisconsin', 'michigan:south')
# states = c('michigan:north', 'minnesota', 'wisconsin')#, 'michigan:south')
# states = c('michigan:south')

wd = getwd() 

us.shp  = readShapeLines(path_map, proj4string=CRS('+init=epsg:3175'))
# us.shp = readOGR('data/map_data/us_alb.shp', 'us_alb')

dist.scale = 1000000 # megamaters

#####################################################################################
# read in data
#####################################################################################
veg    = read.table(file=file.path(wd, path_veg),    sep=',', row.names=NULL, header=TRUE)

# conversion table
convert = read.table(file.path(wd, path_convert), sep=',', row.names=1, header=TRUE)

#####################################################################################
# clean up veg (taxa names and proportional values)
#####################################################################################

taxa.start.col = min(match(rownames(convert), colnames(veg)), na.rm=TRUE)
veg_props      = veg[,taxa.start.col:ncol(veg)]
veg.names      = colnames(veg)[taxa.start.col:ncol(veg)] 
veg_meta       = veg[,1:(taxa.start.col-1)]

veg_meta_tmp = split_mi(veg_meta, longlat=FALSE)

plot(veg_meta$x, veg_meta$y)

veg_meta  = veg_meta[veg_meta_tmp$state2 %in% states, ]
veg_props = veg_props[veg_meta_tmp$state2 %in% states, ]

plot(veg_meta$x, veg_meta$y)

colnames(veg_props) = as.vector(convert[match(veg.names, rownames(convert)),1])

veg_props.collapse = sapply(unique(colnames(veg_props)), 
                            function(x) rowSums( veg_props[ , grep(x, names(veg_props)), drop=FALSE]) )
veg_props          = veg_props.collapse[,sort(colnames(veg_props.collapse))]

# check for proportions
if (sum(rowSums(veg_props)) != nrow(veg_props)) {
  veg_props = as.data.frame(t(apply(veg_props, 1, compute_props)))
  print('Warning: provided data file not give in proportional abundance! 
        Rescaling to data to proportional values.')
}

#####################################################################################
# check for cells that can cause problems in the model
# ie: cells with no trees, or cells with missing data
# for raw PLS: expect some of each
# for comp model output: should not have any of either!
#####################################################################################

idx_none = which(rowSums(veg_props) == 0)
if (length(idx_none)>0) print('Warning: Some cells contain no trees!')

idx_na   = which(is.na(veg_props[,'OAK'])) # if one taxon is NA all are NA
if (length(idx_na)>0) print('Warning: Some cells have no data!')

idx_bad  = c(idx_none, idx_na)

if (length(idx_bad) > 0){
  plot(veg_props$y[,1], veg_props$y[,2])
  points(veg_props$x[idx_bad,1], veg_props$y[idx_bad,2], col='blue', pch=19)
  plot(us.shp, add=T, lwd=2)
}

# # remove problem cells?
# veg_props = veg_props[-idx_bad,]
# veg_meta  = veg_meta[-idx_bad,]

#####################################################################################
## define variables
#####################################################################################

N_cells = nrow(veg_props)
r_comp  = unname(as.matrix(veg_props))
taxa    = colnames(veg_props)
K       = length(taxa)

centers_comp = cbind(veg_meta$x, veg_meta$y)

#####################################################################################
# save the data; stan needs dump file, but use rdata for processing
#####################################################################################

save(K, N_cells,
     r_comp,  
     taxa, centers_comp,
     file=paste(wd, '/comp_data_', K, 'taxa_', suff, '.rdata', sep=""))
