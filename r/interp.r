library(ggplot2)
library(fields, quietly=TRUE)
require(raster)
require(gstat)
require(automap)

wd = getwd()

#####################################################################################
# user pars
#####################################################################################

path_utils = 'calibration/r/utils'
path_figs  = 'calibration/figures'
path_data  = 'r/dump'
path_out   = 'calibration/output'

suff=''

suff_dat = '12taxa_mid_comp_v0.1'

rescale = 1e6

#####################################################################################
# read in data and source utils
#####################################################################################

load(sprintf('%s/cal_data_%s.rdata', path_data, suff_dat))
props = t(apply(y, 1, function(x) if(sum(x)==0){ rep(0, length(x)) } else { x / sum(x) } ))
# props = t(apply(y, 1, function(x) if(sum(x)==0){ rep(0, length(x)) } else { round(x / sum(x) * 500)} ))
# props=y
colnames(props)=taxa

#####################################################################################
# read in data and source utils
#####################################################################################

set.seed(234)
pol = data.frame(centers_polA, pine=props[,'PINE'])
coordinates(pol) =~x+y

centers_veg = data.frame(centers_veg)
colnames(centers_veg) = c('x', 'y')
coordinates(centers_veg) = ~x+y
gridded(centers_veg)<-TRUE

## create a grid onto which we will interpolate:
## first get the range in data
x.range <- as.integer(range(pol@coords[,1]))
y.range <- as.integer(range(pol@coords[,2]))

## now expand to a grid with 500 meter spacing:
grd <- expand.grid(x=seq(from=x.range[1], to=x.range[2], by=24000), y=seq(from=y.range[1], to=y.range[2], by=24000) )

## convert to SpatialPixel class
coordinates(grd) <- ~ x+y
gridded(grd) <- TRUE

## test it out:
plot(grd, cex=0.5)
points(pol, pch=1, col='red', cex=0.7)
title("Interpolation Grid and Sample Points")

## make gstat object:
g <- gstat(id="pine", formula=pine ~ 1, data=pol)

# variogram = autofitVariogram(pine~1,pol) 
variogram = autofitVariogram(log(pine)~1,pol) 
plot(variogram) 

## update the gstat object:
g <- gstat(g, id="pine", model=variogram )


# p = autoKrige(pine ~ 1, pol, centers_veg) 
### perform ordinary kriging prediction:
p = autoKrige(log(pine) ~ 1, pol, centers_veg) 
### perform ordinary kriging prediction:
# p <- predict(g, model=variogram, newdata=grd)

## visualize it:
plot(p)

krige_out = exp(p$krige_output$var1.pred)

out = data.frame(centers_veg, pred=krige_out)
# out = data.frame(centers_veg, pred=krige_out)
ggplot(out) + geom_tile(data=out, aes(x=x, y=y, fill=pred)) + scale_fill_gradientn(colours=tim.colors(), name="props") + coord_fixed() #+


breaks = c(0, 10, 50, 150, 250, 350)
breaklabels = apply(cbind(breaks[1:(length(breaks)-1)], breaks[2:length(breaks)]), 1, 
                    function(r) { sprintf("%0.2f - %0.2f", r[1], r[2]) })

out$pred = cut(out$pred, breaks, include.lowest=TRUE, labels=FALSE)

ggplot(out) + geom_tile(data=out, aes(x=x, y=y, fill=factor(pred))) +  scale_fill_manual(values = tim.colors(length(breaks)), labels=breaklabels, name='Proportions') + 
  coord_fixed() #+ scale_x_continuous(limits$xlims) + scale_y_continuous(limits$ylims)

#####################################################################################
# fit with gams
#####################################################################################

##############################################################################
library(reshape2)
library(mcgv)

load(sprintf('%s/cal_data_%s.rdata', path_data, suff_dat))
# props = t(apply(y, 1, function(x) if(sum(x)==0){ rep(0, length(x)) } else { x / sum(x) } ))
props = t(apply(y, 1, function(x) if(sum(x)==0){ rep(0, length(x)) } else { round(x / sum(x) * 500)} ))
colnames(props) = taxa
colnames(centers_veg) = c('x','y')



for (i in 1:length(taxa)){
  df = data.frame(n=rowSums(y), s=props[,i], f=500-props[,i], centers_polA)
  model = gam(cbind(s,f)~s(x,y), data=df, family=binomial)
  out   = predict(model, newdata=data.frame(centers_veg), type='response')
  #   preds = rbinom(rep(1, nrow(centers_veg)),rep(500, nrow(centers_veg)),out)
  # 
  #   foo = data.frame(centers_veg, pred=preds)
  foo = data.frame(centers_veg, pred=out)
  
  p <- ggplot(foo) + geom_tile(data=foo, aes(x=x, y=y, fill=pred)) + scale_fill_gradientn(colours=tim.colors(), name="props") + 
    coord_fixed() 
  print(p)
  
  breaks = c(0, 0.01, 0.05, 0.10, 0.15, 0.2, 1)
  
  breaklabels = apply(cbind(breaks[1:(length(breaks)-1)], breaks[2:length(breaks)]), 1, 
                      function(r) { sprintf("%0.2f - %0.2f", r[1], r[2]) })
  
  

  foo$discrete = cut(foo$pred, breaks, include.lowest=TRUE, labels=FALSE)

  
  cols = tim.colors(length(breaks))
  pol_facs = sort(unique(foo$discrete))
  
  p <- ggplot() + geom_tile(data=foo, aes(x=x, y=y, fill=factor(discrete)), alpha=1) + 
    scale_fill_manual(values = cols[pol_facs], labels=breaklabels, name='Proportions') + 
    coord_fixed() 
  p <- add_map_albers(p, us.shp, limits, rescale)
  p <- theme_clean(p) 
  p <- p + theme(strip.background = element_blank(), 
                 panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank())
  
  p
}

# 
# dat = data.frame(centers_polA, props)
# dat_long=list()
# for (i in 1:length(taxa)){
#   dat_long = rbind(dat_long, cbind(centers_polA, props[,i], rep(taxa[i], length=nrow(props))))
# }
# colnames(dat_long) = c('x', 'y', 'props', 'taxon')
# 
# dat_sub = dat_long[which(dat_long$taxon == 'PINE'),]
# model = gam(props~s(x,y), data=dat_sub, family=binomial)
# out = predict(model, newdata=data.frame(centers_veg))