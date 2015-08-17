# compute proportions
compute_props <- function(x){
  if (sum(x) == 0){
    props = x
  } else {
    props = x/sum(x)
  }
  return(props)
}

# split the data to separate michigan upper and lower peninsula
split_mi <- function(meta, longlat){

  if (any(colnames(meta)=='region')){
    meta$state = meta$region
  } 
  
  if (longlat){
    centers_ll = data.frame(x=meta$long, y=meta$lat)
  } else {
    centers = data.frame(x=meta$x, y=meta$y)
    
    coordinates(centers) <- ~x + y
    proj4string(centers) <- CRS('+init=epsg:3175')
    
    centers_ll <- spTransform(centers, CRS('+proj=longlat +ellps=WGS84'))
    centers_ll <- as.matrix(data.frame(centers_ll))
  }
  
  idx.mi = which(meta$state=='michigan_north')
  meta$state2 = as.vector(meta$state)
  meta$state2[idx.mi] = map.where(database="state", centers_ll[idx.mi,1], centers_ll[idx.mi,2])
  
  if (any(is.na(meta$state2))){
    idx.na = which(is.na(meta$state2))
  }
  idx.not.na = which(!is.na(meta$state2))
  
  idx.mi.s = which(meta$state=='michigan_south')
  meta$state2[idx.mi.s] = 'michigan:south'
  
  if (length(idx.na)>0){
    for (i in 1:length(idx.na)){
      idx = idx.na[i]
      
      centers = centers_ll[idx.not.na,]
      dmat    = rdist(matrix(centers_ll[idx,], nrow=1) , matrix(centers, ncol=2))
      min.val = dmat[1,which.min(dmat[which(dmat>1e-10)])]
      
      idx_close = which(dmat == min.val)
      state     = map.where(database="state", centers[idx_close,1], centers[idx_close,2])
      
      meta$state2[idx] = state
    }
  }
  
  meta$state2[which(meta$state2[idx.mi]=='minnesota')] = 'michigan:north'
  
  idx.bad = which((meta$state2=='michigan:north') & (meta$y<8e5))
  meta$state2[idx.bad] = 'michigan:south'
  
  return(meta)
  
}


# split the data to separate michigan upper and lower peninsula
get_state <- function(meta, longlat){
  
  if (longlat){
    centers_ll = data.frame(x=meta$long, y=meta$lat)
  } else {
    centers = data.frame(x=meta$x, y=meta$y)
    
    coordinates(centers) <- ~x + y
    proj4string(centers) <- CRS('+init=epsg:3175')
    
    centers_ll <- spTransform(centers, CRS('+proj=longlat +ellps=WGS84'))
    centers_ll <- as.matrix(data.frame(centers_ll))
  }
  
#   idx.mi = which(meta$state=='michigan_north')
#   meta$state2 = as.vector(meta$state)
  meta$state = map.where(database="state", centers_ll[,1], centers_ll[,2])
  
  meta$state[which(meta$state =='indiana')]  = 'michigan:south'
  meta$state[which(meta$state =='ohio')] = 'michigan:south'
  meta$state[which(meta$state =='north dakota')] = 'minnesota'
  meta$state[which(meta$state =='south dakota')] = 'minnesota'
  meta$state[which(meta$state =='illinois')] = 'wisconsin'

  if (any(is.na(meta$state))){
    idx.na = c(which(is.na(meta$state)), which(meta$state=='iowa'))
  }
  idx.not.na = which(!is.na(meta$state))
#   
#   idx.mi.s = which(meta$state=='michigan_south')
#   meta$state2[idx.mi.s] = 'michigan:south'
#   
  if (length(idx.na)>0){
    for (i in 1:length(idx.na)){
      idx = idx.na[i]
      
      centers = centers_ll[idx.not.na,]
      dmat    = rdist(matrix(centers_ll[idx,], nrow=1) , matrix(centers, ncol=2))
      min.val = dmat[1,which.min(dmat[which(dmat>1e-10)])]
      
      idx_close = which(dmat == min.val)
      state     = map.where(database="state", centers[idx_close,1], centers[idx_close,2])
      
      meta$state[idx] = state
    }
  }

  meta$state[which(meta$state =='iowa')] = 'wisconsin'

  return(meta)
  
}

