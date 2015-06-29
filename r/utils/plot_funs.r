library(ggplot2)
require(gridExtra)
library(maptools)
library(RColorBrewer)

us.shp <- readShapeLines('data/map_data/us_alb.shp',
                         proj4string=CRS('+init=epsg:3175'))
us.shp@data$id <- rownames(us.shp@data)
us.fort <- fortify(us.shp, region='id') 

# trace plots
# fit is a stanfit object
# can pass true values as a list if known
trace_plots <- function(fit, npars, N_cores, suff="", save_plots=TRUE, fpath){
  
  post = rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)
  
  n      = dim(post)[3]
  idx = c(seq(1,n-N_cores-1), n)
  labels = colnames(post[,1,idx])
  
  ntrace = length(idx)
  
  avg = summary(fit)$summary[,"mean"]
  
  if (save_plots){
    if (nchar(suff)>1){
      suff = paste0('_', suff)
    }
    pdf(paste(fpath, "/cal_trace", suff, ".pdf", sep=""), width=8, height=6)
  }
  
  par(mfrow=c(1,1))
  for (i in 1:ntrace){
    
    plot(post[,1,idx[i]], type="l", ylab=labels[i], xlab="iter")
    abline(h=avg[idx[i]], col="blue")
    abline(h=summary(fit)$summary[,"2.5%"][idx[i]], col='blue', lty=2)
    abline(h=summary(fit)$summary[,"97.5%"][idx[i]], col='blue', lty=2)
  }
  
  if (save_plots){
    dev.off()
  }
}


#plot of phi by taxa, and normalized log phi where mean of 
#phi is set to 1 at each iter
phi_plot <- function(fit, taxa, save_plots){
  
  post  = rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)
  iters = dim(post)[1]
  K     = dim(post)[3] - 3
  
  labels = colnames(post[,1,1:K])
  
  phi_post     = post[,1,1:K]
  phi_post_std = t(apply(log(phi_post), 1, function(x){x - mean(x) + 1}) )
  
  phi_stats <<- get_phi_stats(phi_post, taxa)
  phi_std_stats <<- get_phi_stats(phi_post_std, taxa)
  
  p1 <- ggplot(phi_stats, aes(y=mean, x=taxa)) + geom_point(size=4) + #theme_bw() + 
    geom_errorbar(aes(ymax=phi_stats[,'U'], ymin=phi_stats[,'L']), width=.05) + ylab("Phi") + xlab("Taxa") + 
    scale_x_discrete(limits=taxa)
  p2 <- ggplot(phi_std_stats, aes(y=mean, x=taxa)) + geom_point(size=4) + #theme_bw() + 
    geom_errorbar(aes(ymax=phi_std_stats[,'U'], ymin=phi_std_stats[,'L']), width=.05) + ylab("log shifted Phi") + xlab("Taxa") + 
    scale_x_discrete(limits=taxa)
  
  if (save_plots){
    pdf(paste("calibration/figures/pollen_phi_plot_", suff, ".pdf", sep=""), width=8, height=6)
  }
  grid.arrange(p1, p2, ncol=2)
  if (save_plots){
    dev.off()
  }
}

#plot raw versus phi scaled pond cell veg
local_pollen_veg_plot <- function(phi, preds, N_cores, r, idx_cores, taxa, suff, save_plots){  
  
  if (save_plots){
    pdf(paste("calibration/figures/pollen_focal_scaled_", suff, ".pdf", sep=""), width=12, height=6)
    #     postscript(paste('calibration/figures/pollen_focal_scaled_', suff, '.eps', sep=''), width=10, height=10)
  }
  
  par(pty="s")
  par(mfrow=c(2,6))
  for (k in 1:K){
    
    plot(0,0, type='n', xlim=c(0,1), ylim=c(0,1), xlab='veg or predicted pollen', ylab='raw pollen', 
         main = taxa[k], cex.main=2, cex.lab=2, cex.axis=1.5)
    
    for (i in 1:N_cores){
      
      #raw pollen data versus focal cell veg
      points(r[idx_cores[i],k], pollen_props[i,k], col='red', pch=3, cex=1.6) 
      
      #raw pollen data versus focal cell veg weighted by phi  
      points(preds[i,k], pollen_props[i,k], col='black', pch=20, cex=1.6)    
      
      abline(a=0, b=1, col='grey', lty=2, lwd=1.8)
    }
    
    #legend("topleft", c("veg", "pollen preds"), col=c('red', 'black'), pch=c(3,20))
  }
  
  if (save_plots){
    dev.off()
  }
}

#plot raw versus phi scaled pond cell veg
local_pollen_veg_plot2 <- function(r, idx_cores, pollen_props, local_preds, taxa, suff, save_plots, fpath){  
  
  N_cores = length(idx_cores)
  
  if (!is.null(taxa)){
    taxa[which(taxa == 'OTHER.HARDWOOD')] = 'OTHER HW'
    taxa[which(taxa == 'OTHER.CONIFER')] = 'OTHER CON'
    #     taxa[which(taxa == 'TAMARACK')] = 'LARCH'
    
    #     taxa = as.vector(sapply(taxa, simpleCap))
  }
  
  if (save_plots){
    if (nchar(suff)>1){
      suff = paste0('_', suff)
    }
    pdf(paste(fpath, "/pollen_focal_scaled", suff, ".pdf", sep=""), width=12, height=10)
    #     postscript(paste('calibration/figures/pollen_focal_scaled_', suff, '.eps', sep=''), width=10, height=10)
  }
  
  par(pty="s")
  par(mfrow=c(3,4))
  for (k in 1:K){
    
    plot(0,0, type='n', xlim=c(0,1), ylim=c(0,1), xlab='veg or predicted pollen', ylab='raw pollen', 
         main = taxa[k], cex.main=1.5, cex.lab=1.5, cex.axis=1.5)
    
    for (i in 1:N_cores){
      
      #raw pollen data versus focal cell veg
      points(r[idx_cores[i],k], pollen_props[i,k], col='red', pch=3, cex=1.6) 
      
      #raw pollen data versus focal cell veg weighted by phi  
      points(local_preds[i,k], pollen_props[i,k], col='black', pch=20, cex=1.6)    
      
      abline(a=0, b=1, col='grey', lty=2, lwd=1.8)
    }
    
    #legend("topleft", c("veg", "pollen preds"), col=c('red', 'black'), pch=c(3,20))
  }
  
  if (save_plots){
    dev.off()
  }
}


#plot raw versus weighted nhood preds
pollen_preds_plot <- function(preds, pollen_props, N_cores, r, idx_cores, taxa, suff, save_plots, fpath){  
  
  if (!is.null(taxa)){
    taxa[which(taxa == 'OTHER.HARDWOOD')] = 'OTHER HW'
    taxa[which(taxa == 'OTHER.CONIFER')] = 'OTHER CON'
    #     taxa[which(taxa == 'TAMARACK')] = 'LARCH'
    
    #     taxa = as.vector(sapply(taxa, simpleCap))
  }
  
  if (save_plots){
    if (nchar(suff)>1){
      suff = paste0('_', suff)
    }
    pdf(paste(fpath, "/pollen_preds", suff, ".pdf", sep=""), width=12, height=10)
    #     postscript(paste(fpath, '/pollen_preds_', suff, '.epse', sep=''), width=12, height=10)
  }
  
  par(pty="s")
  par(mfrow=c(3,4))
  #   par(mar=c(8,4.1,4.1,2.1))
  for (k in 1:K){
    
    plot(0,0, type='n', xlim=c(0,1), ylim=c(0,1), xlab='veg or predicted pollen', ylab='raw pollen', 
         main = taxa[k], cex.main=1.5, cex.lab=1.5, cex.axis=1.5)
    
    for (i in 1:N_cores){
      
      #raw pollen data
      points(r[idx_cores[i],k], pollen_props[i,k], col='red', pch=3, cex=1.6) 
      #     #weighted by phi for focal grid cell
      #     points(phi[k]*r[idx_cores[i],k], pollen_props[i,k], col='blue', pch=18)
      #weighted nhood preds
      
      points(preds[i,k], pollen_props[i,k], col='black', pch=20, cex=1.6)    
      
      abline(a=0, b=1, col='grey', lty=2, lwd=1.8)
    }
    
    #legend("topleft", c("veg", "pollen preds"), col=c('red', 'black'), pch=c(3,20))
  }
  
  if (save_plots){
    dev.off()
  }
}

plot_data_maps <- function(y, centers, taxa, K, limits, suff, save_plots, fpath=subDir){
  
  rescale=1000000
  N = nrow(centers)
  
  if (is.null(taxa)){taxa=seq(1,K)}
  
  #   y = t(y)
  
  props_data = t(apply(y, 1, function(x) if (sum(x) > 0) {x/sum(x)} else {x}))
  #colnames(props_data) = taxa
  
  prop_dat = data.frame(props=numeric(0), x=integer(0), y=integer(0), taxon=character())
  for (k in 1:K){
    prop_dat = rbind(prop_dat, data.frame(props = props_data[,k], 
                                          x     = centers[,1],#*rescale, 
                                          y     = centers[,2],#*rescale, 
                                          taxon = rep(taxa[k], N)))
  }
  
  prop_dat$type = rep('PLS', nrow(prop_dat))
  
  p <- ggplot() + geom_tile(data=prop_dat, aes(x=x, y=y, fill=props)) + 
    scale_fill_gradientn(colours=tim.colors(), name='Proportions') + 
    coord_fixed() #+ scale_x_continuous(limits$xlims) + scale_y_continuous(limits$ylims)
  p <- add_map_albers(p, us.shp, limits, rescale)
  p <- p + facet_wrap(~taxon, ncol=6)
  p <- theme_clean(p) #+ theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.5)))
  
  #   p <- p + theme(strip.text.x = element_blank(),
  #                  strip.text.y = element_blank())
  p <- p + theme(strip.background = element_blank())
  
  print(p)
  Sys.sleep(2)
  if (save_plots){
    ggsave(file=paste(fpath, '/maps_', suff, '_cont.pdf', sep=''), scale=1)
    ggsave(file=paste(fpath, '/maps_', suff, '_cont.eps', sep=''), scale=1)
    #     dev.off()
  }
  return(p)
}

plot_data_maps_binned <- function(y, centers, taxa, K, breaks, limits, suff, save_plots, fpath=subDir){
  
  rescale=1000000
  N = nrow(centers)
  
  if (is.null(taxa)){taxa=seq(1,K)}
  
  #   y = t(y)
  
  props_data = t(apply(y, 1, function(x) if (sum(x) > 0) {x/sum(x)} else {x}))
  #colnames(props_data) = taxa
  
  props_data_binned = matrix(0, nrow=nrow(props_data), ncol=ncol(props_data))
  colnames(props_data_binned) <- colnames(props_data)
  
  for (i in 1:ncol(props_data)){
    props_data_binned[,i] = cut(props_data[,i], breaks, include.lowest=TRUE, labels=FALSE)
  }
  
  breaklabels = apply(cbind(breaks[1:(length(breaks)-1)], breaks[2:length(breaks)]), 1, 
                      function(r) { sprintf("%0.2f - %0.2f", r[1], r[2]) })
  
  prop_dat = data.frame(props=numeric(0), x=integer(0), y=integer(0), taxon=character())
  for (k in 1:K){
    prop_dat = rbind(prop_dat, data.frame(props = props_data_binned[,k], 
                                          x     = centers[,1],#*rescale, 
                                          y     = centers[,2],#*rescale, 
                                          taxon = rep(taxa[k], N)))
  }
  
  prop_dat$type = rep('PLS', nrow(prop_dat))
  
  cols = rev(brewer.pal(length(breaks), name="BrBG"))
#   cols = tim.colors(length(breaks))
#   cols = gray.colors(length(breaks), start = 0.0, end = 1.0, gamma = 2.2, alpha = NULL)
  
  p <- ggplot() + geom_tile(data=prop_dat, aes(x=x, y=y, fill=factor(props))) + 
    scale_fill_manual(values = cols, labels=breaklabels, name='Proportions') + 
    coord_fixed() #+ scale_x_continuous(limits$xlims) + scale_y_continuous(limits$ylims)
  p <- add_map_albers(p, us.shp, limits, rescale)
  p <- p + facet_wrap(~taxon, ncol=6)
  p <- theme_clean(p) #+ theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.5)))
  
  #   p <- p + theme(strip.text.x = element_blank(),
  #                  strip.text.y = element_blank())
  p <- p + theme(strip.background = element_blank())
  
  print(p)
  Sys.sleep(2)
  
  if (save_plots){
    ggsave(file=paste(fpath, '/maps_', suff, '.pdf', sep=''), scale=1)
    ggsave(file=paste(fpath, '/maps_', suff, '.eps', sep=''), scale=1)
  }
  return(p)
}

plot_pollen_maps_binned <- function(y, centers, taxa, K, breaks, limits, suff, save_plots, fpath){
  
  N = nrow(centers)
  
  if (is.null(taxa)){taxa=seq(1,K)}
  
  props_data = t(apply(y, 1, function(x) if (sum(x) > 0) {x/sum(x)} else {x}))
  #colnames(props_data) = taxa
  
  props_data_binned = matrix(0, nrow=nrow(props_data), ncol=ncol(props_data))
  colnames(props_data_binned) <- colnames(props_data)
  
  for (i in 1:ncol(props_data)){
    props_data_binned[,i] = cut(props_data[,i], breaks, include.lowest=TRUE, labels=FALSE)
  }
  
  breaklabels = apply(cbind(breaks[1:(length(breaks)-1)], breaks[2:length(breaks)]), 1, 
                      function(r) { sprintf("%0.2f - %0.2f", r[1], r[2]) })
  
  prop_dat = data.frame(props=numeric(0), x=integer(0), y=integer(0), taxon=character())
  for (k in 1:K){
    prop_dat = rbind(prop_dat, data.frame(props = props_data_binned[,k], 
                                          x     = centers[,1],#*rescale, 
                                          y     = centers[,2],#*rescale, 
                                          taxon = rep(taxa[k], N)))
  }
  
#   cols = tim.colors(length(breaks))
  cols = rev(brewer.pal(length(breaks), 'BrBG'))
  
  p <- ggplot() + geom_point(data=prop_dat, aes(x=x, y=y, colour=factor(props)), shape=19) + 
#     scale_colour_manual(values = tim.colors(length(breaks)), labels=breaklabels, name='Proportions') + 
    scale_colour_manual(values = cols, labels=breaklabels, name='Proportions') + 
    coord_fixed() #+ scale_x_continuous(limits$xlims) + scale_y_continuous(limits$ylims)
  p <- add_map_albers(p, us.shp, limits, rescale)
  p <- p + facet_wrap(~taxon, ncol=6)
  p <- theme_clean(p) #+ theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.5)))
  
  #   p <- p + theme(strip.text.x = element_blank(),
  #                  strip.text.y = element_blank())
  p <- p + theme(strip.background = element_blank())
  
  print(p)
  Sys.sleep(2)
  if (save_plots){
    ggsave(file=paste(fpath, '/maps_pollen.pdf', sep=''), scale=1)
    ggsave(file=paste(fpath, '/maps_pollen.eps', sep=''), scale=1)
    #     dev.off()
  }
  return(p)
}

plot_smoothed_pollen_maps_binned <- function(y, centers_pol, centers_veg, taxa, K, breaks, limits, suff, save_plots, fpath){
  
  N = nrow(centers_veg)
  
  if (is.null(taxa)){taxa=seq(1,K)}
  
  colnames(y) = taxa
  colnames(centers_veg) = c('x', 'y')
  
  pred_dat = data.frame(pred=numeric(0), x=integer(0), y=integer(0), taxon=character())
  for (k in 1:K){
    
    dat   = data.frame(s=y[,k], f=500-y[,k], centers_pol)
    model = gam(cbind(s,f)~s(x,y), data=dat, family=binomial)
    out   = predict(model, newdata=data.frame(centers_veg), type='response')
    
    pred_dat = rbind(pred_dat, data.frame(pred = out, 
                                          x     = centers_veg[,1],#*rescale, 
                                          y     = centers_veg[,2],#*rescale, 
                                          taxon = rep(taxa[k], N)))
  }
  
  pred_dat$pred_bin = cut(pred_dat$pred, breaks, include.lowest=TRUE, labels=FALSE)
  
  breaklabels = apply(cbind(breaks[1:(length(breaks)-1)], breaks[2:length(breaks)]), 1, 
                      function(r) { sprintf("%0.2f - %0.2f", r[1], r[2]) })
  
  p <- ggplot() + geom_point(data=pred_dat, aes(x=x, y=y, colour=factor(pred_bin)), shape=19) + 
    scale_colour_manual(values = tim.colors(length(breaks)), labels=breaklabels, name='Proportions') + 
    coord_fixed() #+ scale_x_continuous(limits$xlims) + scale_y_continuous(limits$ylims)
  p <- add_map_albers(p, us.shp, limits, rescale)
  p <- p + facet_wrap(~taxon, ncol=6)
  p <- theme_clean(p) #+ theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.5)))
  
  #   p <- p + theme(strip.text.x = element_blank(),
  #                  strip.text.y = element_blank())
  p <- p + theme(strip.background = element_blank())
  
  print(p)
  Sys.sleep(2)
  if (save_plots){
    ggsave(file=paste(fpath, '/maps_smoothed_pollen.pdf', sep=''), scale=1)
    ggsave(file=paste(fpath, '/maps_smoothed_pollen.eps', sep=''), scale=1)
    #     dev.off()
  }
  return(p)
}


melt_dat <- function(y, centers, breaks, taxa){
  N=nrow(centers)
  props_data = t(apply(y, 1, function(x) if (sum(x) > 0) {x/sum(x)} else {x}))
  colnames(props_data) = taxa
  
  props_data_binned = matrix(0, nrow=nrow(y), ncol=ncol(y))
  colnames(props_data_binned) <- taxa
  
  for (i in 1:ncol(props_data)){
    props_data_binned[,i] = cut(props_data[,i], breaks, include.lowest=TRUE, labels=FALSE)
  }
  
  prop_dat = data.frame(props=numeric(0), x=integer(0), y=integer(0), taxon=character())
  for (k in 1:K){
    prop_dat = rbind(prop_dat, data.frame(props = props_data_binned[,k], 
                                          x     = centers[,1],#*rescale, 
                                          y     = centers[,2],#*rescale, 
                                          taxon = rep(taxa[k], N)))
    #     prop_dat = rbind(prop_dat, data.frame(props = seq(1,10),
    #                                           x = 0,
    #                                           y = 0, 
    #                                           taxon = rep(taxa[k], 10)))
  }
  
  return(prop_dat)
  
}

plot_both_maps_binned <- function(y_pol,  y_veg, centers_pol, centers_veg, taxa, taxa_list, K, breaks, limits, suff, save_plots, fpath){
  
  N_pol = nrow(centers_pol)
  
  if (is.null(taxa)){taxa=seq(1,K)}
  
  breaklabels = apply(cbind(breaks[1:(length(breaks)-1)], breaks[2:length(breaks)]), 1, 
                      function(r) { sprintf("%0.2f - %0.2f", r[1], r[2]) })
  
  prop_dat_veg = melt_dat(y_veg, centers_veg, breaks, taxa)
  prop_dat_pol = melt_dat(y_pol, centers_pol, breaks, taxa)
  
  for (k in taxa_list){#1:length(taxa)){
    
    veg = prop_dat_veg[prop_dat_veg$taxon == taxa[k], ]
    pol = prop_dat_pol[prop_dat_pol$taxon == taxa[k], ]
    
    cols = tim.colors(length(breaks))
#     cols = rev(brewer.pal(length(breaks), 'BrBG'))
#     cols = terrain.colors(length(breaks), alpha = 1)
    veg_facs = sort(unique(veg$props))
    pol_facs = sort(unique(pol$props))
    
    p <- ggplot() + geom_tile(data=veg, aes(x=x, y=y, fill=factor(props)), alpha=1) + 
      scale_fill_manual(values = cols[veg_facs], labels=breaklabels, name='Proportions') + 
      coord_fixed() 
    p <- add_map_albers(p, us.shp, limits, rescale)
    p <- theme_clean(p) 
    p <- p + theme(strip.background = element_blank(), 
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank())
    
    q <- ggplot() + geom_point(data=pol, aes(x=x, y=y, colour=factor(props)), shape=19, size=3) + 
      #       scale_colour_manual(values = brewer.pal(length(breaks),"Spectral")) + 
      scale_colour_manual(values = cols[pol_facs], labels=breaklabels, name='Proportions') + 
      coord_fixed() #+ scale_x_continuous(limits$xlims) + scale_y_continuous(limits$ylims)
    q <- add_map_albers(q, us.shp, limits, rescale)
    q <- theme_clean(q) 
    q <- q + theme(strip.background = element_blank(), 
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank())
    #q <- q + theme(legend.position="none")
    
    g <- arrangeGrob(p, q, nrow=2)
    print(g)
    
    Sys.sleep(2)
    if (save_plots){
      ggsave(file=paste(fpath, '/maps_compare_', taxa[k], '.pdf', sep=''), scale=1, plot=g)
      ggsave(file=paste(fpath, '/maps_compare_', taxa[k], '.eps', sep=''), scale=1, plot=g)
      #     dev.off()
    }
    
  }
  
}

plot_sumw <- function(dat, fpath){
  rescale = 1e6
  p <- ggplot() + geom_point(data=dat, aes(x=x, y=y, size=sum_w), shape=21, colour="black", fill="chartreuse4")
  p <- add_map_albers(p, us.shp, limits, rescale) 
  p <- p + theme(panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank())
  p <- p + theme(axis.ticks = element_blank(), 
                 axis.text = element_blank(), 
                 axis.title = element_blank())
  print(p) 
  
  Sys.sleep(2)
  if (save_plots){
    ggsave(file=paste(fpath, '/sum_w_maps.pdf', sep=''), scale=1)
    ggsave(file=paste(fpath, '/sum_w_maps.eps', sep=''), scale=1)
  }
}

plot_alpha <- function(dat, fpath){
  rescale = 1e6
  p <- ggplot() + geom_point(data=dat, aes(x=x, y=y, size=alpha), shape=21, colour="black", fill="chartreuse4")
  p <- add_map_albers(p, us.shp, limits, rescale) + theme(panel.grid.major = element_blank(), 
                                                          panel.grid.minor = element_blank())
  p <- p + theme(axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank())
  print(p) 
  
  Sys.sleep(2)
  if (save_plots){
    ggsave(file=paste(fpath, '/alpha_maps.pdf', sep=''), scale=1)
    ggsave(file=paste(fpath, '/alpha_maps.eps', sep=''), scale=1)
  }
}

plot_par_vals <- function(post, parname, taxa, wd, path_figs){
  
  col_names = colnames(post[,1,])
  par_names  = unlist(lapply(col_names, function(x) strsplit(x, "\\[")[[1]][1]))
  par_vals  = post[,1,which(par_names == parname)]
  
  par_mean  = colMeans(par_vals)
  par_lb    = apply(par_vals, 2, function(x) quantile(x, probs=0.025))
  par_ub    = apply(par_vals, 2, function(x) quantile(x, probs=0.975))
  
  par_stats = data.frame(name=taxa, mu=par_mean, lb=par_lb, ub=par_ub)
  
  p <- ggplot(data=par_stats, aes(x=reorder(name, mu), y=mu)) + 
    geom_point(size=4) + geom_errorbar(aes(ymin=lb, ymax=ub), width=.2) + 
    xlab("Taxon") + ylab(parse(text=parname)) +
    coord_flip() + theme_bw() + theme(axis.title.x=element_text(size=20), 
                                      axis.title.y=element_text(size=20), 
                                      axis.text.x=element_text(size=rel(1.3)),
                                      axis.text.y=element_text(size=rel(1.3)))
  
  print(p)
  #   theme(axis.text.x = element_text(angle=60, hjust=1)) + 
  #   theme(panel.background = element_blank())
  ggsave(p, filename=sprintf('%s/%s/%s.pdf', wd, path_figs, parname), width=8, height=6)
}

theme_clean <- function(plot_obj){
  plot_obj <- plot_obj + theme(axis.ticks = element_blank(), 
                               axis.text.y = element_blank(), 
                               axis.text.x = element_blank(),
                               axis.title.x = element_blank(),
                               axis.title.y = element_blank())
  
  return(plot_obj)
}

# add_map_albers <- function(plot_obj, map_data=us.fort, limits){
#   p <- plot_obj + geom_path(data=map_data, aes(x=long, y=lat, group=group),  colour='grey55') + 
#     #     scale_x_continuous(limits = c(min(umw.coord$x, na.rm=TRUE), max(umw.coord$x, na.rm=TRUE))) +
#     #     scale_y_continuous(limits = c(min(umw.coord$y, na.rm=TRUE), max(umw.coord$y, na.rm=TRUE)))#, colour = "black", size = 1, fill = "white", aes(x=long, y=lat, group = group))
#     # #   
#     #     scale_x_continuous(limits = c(min(dat[,1], na.rm=TRUE), max(dat[,1], na.rm=TRUE))) +
#     #     scale_y_continuous(limits = c(min(dat[,2], na.rm=TRUE), max(dat[,2], na.rm=TRUE)))
#     scale_x_continuous(limits = limits$xlims*1000000) +
#     scale_y_continuous(limits = limits$ylims*1000000) #+ coord_map("albers")
#   return(p)
#   
# }

get_limits <- function(centers){
  xlo = min(centers[,1])
  xhi = max(centers[,1])
  
  ylo = min(centers[,2])
  yhi = max(centers[,2])
  
  return(list(xlims=c(xlo,xhi),ylims=c(ylo, yhi)))
}  


add_map_albers <- function(plot_obj, map_data=us.fort, limits, rescale){
  p <- plot_obj + geom_path(data=map_data, aes(x=long, y=lat, group=group),  colour='grey55') + 
    #     scale_x_continuous(limits = c(min(umw.coord$x, na.rm=TRUE), max(umw.coord$x, na.rm=TRUE))) +
    #     scale_y_continuous(limits = c(min(umw.coord$y, na.rm=TRUE), max(umw.coord$y, na.rm=TRUE)))#, colour = "black", size = 1, fill = "white", aes(x=long, y=lat, group = group))
    # #   
    #     scale_x_continuous(limits = c(min(dat[,1], na.rm=TRUE), max(dat[,1], na.rm=TRUE))) +
    #     scale_y_continuous(limits = c(min(dat[,2], na.rm=TRUE), max(dat[,2], na.rm=TRUE)))
    scale_x_continuous(limits = limits$xlims) +
    scale_y_continuous(limits = limits$ylims) #+ coord_map("albers")
  return(p) 
}
