library(sp)
library(rgdal)
library(fields)
library(maptools)
gpclibPermit()

# source('calibration/r/cal_build_config.r')
source('calibration/r/utils/build_data_funs.r')

#####################################################################################
# read in data
#####################################################################################

# set your working directory!
# AD: wd = '~/Documents/paleon/stepps2'

depth_type = 'lower'

data_date = '2014-10-20'

## relative to working directory!
path_pollen  = paste0('sh_analysis/data/cal_data_', depth_type, '_depth_', data_date, '.csv')    # sh elicitation results
path_veg     = 'data/pls/composition_v0.2.csv'              # composition model results
path_map     = 'data/map_data/us_alb.shp'                   # albers projected state map!
path_convert = 'calibration/r/meta/dict-comp2stepps.csv'    # conversion table

path_out = 'calibration/r/dump'                             # dump data file stored here

# append to output filenames
suff       = paste0(depth_type, '_comp_v0.1')                                

states = c('michigan:north', 'minnesota', 'wisconsin')

wd = getwd() 

us.shp     = readShapeLines(path_map, proj4string=CRS('+init=epsg:3175'))
# us.shp = readOGR('data/map_data/us_alb.shp', 'us_alb')

dist.scale = 1000000 # megamaters

#####################################################################################
# read in data
#####################################################################################

pollen = read.table(file=file.path(wd, path_pollen), sep=',', header=TRUE)
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
r       = unname(as.matrix(veg_props))
taxa    = colnames(veg_props)
K       = length(taxa)

centers_veg = cbind(veg_meta$x, veg_meta$y)

#####################################################################################
# reproject pollen coords from lat long to Albers
#####################################################################################

#XXX: add pollen to albers function
centers_pol = data.frame(cbind(pollen$long, pollen$lat))
colnames(centers_pol) = c('x', 'y')

coordinates(centers_pol) <- ~ x + y
proj4string(centers_pol) <- CRS('+proj=longlat +ellps=WGS84')

centers_polA <- spTransform(centers_pol, CRS('+init=epsg:3175'))
centers_polA <- as.matrix(data.frame(centers_polA))
# 
pollen_meta = data.frame(centers_polA, state=as.vector(pollen$state), stringsAsFactors=FALSE)
# pollen_meta$state[pollen_meta$state == 'michigan:north'] = 'michigan_north'
# pollen_meta$state[pollen_meta$state == 'michigan:south'] = 'michigan_south'
pollen_meta_tmp = split_mi(pollen_meta, longlat=FALSE)

pollen_meta  = pollen_meta_tmp[pollen_meta_tmp$state2 %in% states, ]
pollen       = pollen[pollen_meta_tmp$state2 %in% states, ]
centers_polA = pollen_meta[,c('x', 'y')]

pollen_counts = pollen[,min(match(taxa, colnames(pollen))):ncol(pollen)]

N_cores = nrow(pollen_counts)

# make sure the columns are in the right order!
y = pollen_counts[,taxa]
y = unname(round(as.matrix(y))) # stan needs integers and no names

#####################################################################################
# build idx_core and idx_hood
#####################################################################################

idx_cores = vector(length=N_cores)

for (i in 1:nrow(centers_polA)){
  core_site = centers_polA[i,]
  d = rdist(matrix(core_site, ncol=2), as.matrix(centers_veg))
  idx_cores[i] = which.min(d)
}

# visual check: cores and nearest cell should be close to each other!
plot(centers_veg[idx_cores,1], centers_veg[idx_cores,2])
points(centers_polA[,1], centers_polA[,2], col='blue', pch=8)
plot(us.shp, add=T, lwd=2)

# XXX: do we still need this? can we just use idx_cores
# whole domain as neighborhood
idx_hood   = matrix(0, nrow=N_cores, ncol=N_cells)

for (i in 1:nrow(centers_polA)){ 
  hood_cells = seq(1, N_cells)
  hood_cells = hood_cells[-idx_cores[i]]
  idx_hood[i,1:length(hood_cells)] = hood_cells
}

N_hood = ncol(idx_hood)

# distance matrix
d = t(rdist(as.matrix(centers_veg[idx_cores,], ncol=2), as.matrix(centers_veg))/dist.scale) # ponds at centroids - area based
# d = t(rdist(as.matrix(centers_polA, ncol=2), as.matrix(centers_veg))/dist.scale) # ponds as points in cell
d2 = d * d # XXX: is it more efficient to read in d2?
# d2t = t(d2)

#####################################################################################
# calculate potential d
# used to determine C normalizing constant in the non-local contribution term
#####################################################################################
library(plyr)

x_pot = seq(-528000, 528000, by=8000)
y_pot = seq(-416000, 416000, by=8000)
coord_pot = expand.grid(x_pot, y_pot)

d_pot = t(rdist(matrix(c(0,0), ncol=2), as.matrix(coord_pot, ncol=2))/dist.scale)
d_pot = unname(as.matrix(count(d_pot)))

N_pot = nrow(d_pot)

#####################################################################################
# save the data; stan needs dump file, but use rdata for processing
#####################################################################################

dump(c('K', 'N_cores', 'N_cells', 'N_hood', 
       'y', 'r', 
       'idx_cores', 'idx_hood', 
       'd','d2',
       'N_pot', 'd_pot',
       'taxa', 'centers_veg', 'centers_polA'), 
     file=paste(wd, '/', path_out, '/cal_data_', K, 'taxa_', suff,'.dump',sep=""))

save(K, N_cores, N_cells, N_hood, y, r, idx_cores, idx_hood, d, d2, N_pot, d_pot, taxa, centers_veg, centers_polA,
file=paste(wd, '/', path_out, '/cal_data_', K, 'taxa_', suff,'.rdata',sep=""))

#####################################################################################
# find cores for which the entire neighborhood is treeless
#####################################################################################

# idx_bad = vector(mode='numeric', length=0) 
# 
# for (i in 1:N_cores){
#   hood_cells = idx_hood[i,]
#   r_hood     = r[hood_cells,]
#   r_core     = r[idx_cores[i],]
#   
#   r_tot = rbind(r_hood, r_core)
#   
#   if (any(colSums(r_tot) == 0)){
#     print(paste('WARNING: no taxon near core ', i, sep=''))
#     idx_bad=c(idx_bad, i)
#   }
# }
# 
# bad_cores = centers_veg[idx_cores[idx_bad],]
# plot(centers_veg[,1], centers_veg[,2], col='black')
# points(bad_cores[,1], bad_cores[,2], pch=19, col='blue')
# 
