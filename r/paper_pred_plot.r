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

save_plots = TRUE
rescale    = 1e6

###############################################################################################################
# read in data and source utils
###############################################################################################################

source(file.path(path_utils, 'process_funs.r'))
source(file.path(path_utils, 'plot_funs.r'))
source(file.path(path_utils, 'paper_plot_funs.r'))
source(file.path(wd, 'r', 'runs.r'))
gpclibPermit()

load(sprintf('%s/cal_data_%s.rdata', path_data, suff_dat))

taxa_abb = taxa
taxa_abb[which(taxa_abb == 'OTHER.CONIFER')] = 'OTHER CON'
taxa_abb[which(taxa_abb == 'OTHER.HARDWOOD')] = 'OTHER HW'

# colors from pie maps; be consistent
load(file='r/stepps_cols.rdata')
stepps_cols = data.frame(lapply(stepps_cols, as.character), stringsAsFactors=FALSE)
stepps_cols = rbind(stepps_cols, c('FIR', as.vector(stepps_cols$cols[stepps_cols$taxa=='OTHER.CONIFER'])))


# col_list = c("#1F78B4", "#33A02C", "#E31A1C", "#6A3D9A", "#B15928", "#FF7F00", 
#              "#FFFF99", "#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F", "#CAB2D6")
###############################################################################################################
# plot the modelled composition
###############################################################################################################
limits <- get_limits(centers_veg)

breaks = c(0, 0.01, 0.05, 0.10, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 1)
# plot_data_maps_binned(r_mean, centers=centers_pls, taxa=taxa, K, T, breaks, suff=suff4, save_plots=save_plots)

plot_data_maps_binned(r, centers_veg, taxa_abb, K, breaks, limits, suff='veg', save_plots, fpath=path_figs)

###############################################################################################################
###############################################################################################################

pollen_props      = compute_props(y, taxa)
pollen_props_melt = melt(compute_props(y, taxa))

colnames(r) = taxa
r_melt      = melt(r[idx_cores, ])

phi_scaled_dat = data.frame(matrix(0, nrow=0, ncol=6))
pred_dat       = data.frame(matrix(0, nrow=0, ncol=6))
for (run in runs){
  fname = sprintf('%s/%s.csv', path_out, run$suff_fit)
  system(sprintf('r/fixup.pl %s', fname))
  fit   = read_stan_csv(fname)
  post  = rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)
  
  preds_fname = sprintf('%s/%s/%s.rdata', path_figs, run$suff_fit, 'pollen_preds')
  if (!file.exists(preds_fname)){
    sum_w <- build_sumw_pot(post, K, N_pot, d_pot, run)
    preds_out = pollen_preds(post, N_cores, d, idx_cores, r, sum_w, run)
    preds     = preds_out$preds
    colnames(preds) = taxa
    save(preds, file=preds_fname)
  } else {
    load(preds_fname)
  }
  
  preds_melt = melt(preds)
  preds_melt = cbind(preds_melt, pollen_props_melt[,3], r_melt[,3], rep(run$handle, nrow(preds_melt)))
  colnames(preds_melt) = c('core', 'taxon', 'pred', 'data', 'veg', 'run')
  
  pred_dat = rbind(pred_dat, preds_melt)

  phi_out = phi_scale_veg(post, N_cores, r, idx_cores)
  colnames(phi_out) = taxa

  phi_out_melt = melt(phi_out)
  phi_out_melt = cbind(phi_out_melt, pollen_props_melt[,3], r_melt[,3], rep(run$handle, nrow(phi_out_melt)))
  colnames(phi_out_melt) = c('core', 'taxon', 'phi_pred', 'data', 'veg', 'run')
  phi_scaled_dat = rbind(phi_scaled_dat, phi_out_melt)
}

handles = as.vector(unique(pred_dat$run))
for (run in handles) {

  dat        = pred_dat[which(pred_dat$run %in% run),]
  phi_scaled = phi_scaled_dat[which(phi_scaled_dat$run %in% run),]
  
  plot_pollen_preds_paper(dat, path_figs, run)
  plot_pollen_phi_scaled_paper(phi_scaled, path_figs, run) 
}

########################################################################################################################################
# dispersal cdf
########################################################################################################################################

limits <- get_limits(centers_veg)

radius = seq(8000,1000000, by=4000)

dispersal_cdf <- function(radius, runs, taxa, path_out){
  
  K = length(taxa)  
  
  coord_pot = seq(-700000, 700000, by=8000)
  coord_pot = expand.grid(coord_pot, coord_pot)
  
  d_pot = t(rdist(matrix(c(0,0), ncol=2), as.matrix(coord_pot, ncol=2))/rescale)
  
  # want a circular region
  idx_circ  = which(d_pot[,1] < 0.7)
  coord_pot = coord_pot[idx_circ, ]
  d_pot     = d_pot[idx_circ, ]
  
  d_pot = unname(as.matrix(count(d_pot)))
  N_pot = nrow(d_pot)
  
  dat = data.frame(matrix(0, nrow=0, ncol=4))
  r_all = list(length=length(runs))
  for (i in 1:length(runs)){
    run   = runs[[i]]
    fname = sprintf('%s/%s.csv', path_out, run$suff_fit)
    fit   = read_stan_csv(fname)
    post  = rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)
    
#     sum_w = build_sumw_pot(post, K, N_pot, d_pot, run)
#     r_int = dispersal_decay(post, d_pot, sum_w, radius/rescale, run, taxa)
#     colnames(r_int) = taxa
    
    sum_w = build_sumw_pot_ci(post, K, N_pot, d_pot, run)
    r_niter = dispersal_decay_ci(post, d_pot, sum_w, radius/rescale, run, taxa)
    
    r_all[[i]] = list(r_niter, run$handle)
#     r_all[[i]]$handle = run$handle

    r_int = t(apply(r_niter, 1, rowMeans))
    colnames(r_int) = taxa
    
    run_dat = data.frame(radius = radius/1e3, r_int, handle=rep(run$handle, length(radius)))
    dat = rbind(dat, melt(run_dat, id=c('radius', 'handle')))
  }
  
  levels(dat$handle) <- c("Base G", "Base PL", "Variable G", "Variable PL")
  
  return(list(dat=dat, r_all=r_all))
}

plot_dispersal_cdfs <- function(dat){
  
  levels(dat$variable)[levels(dat$variable) == "OTHER.CONIFER"] = "FIR"
  levels(dat$variable)[levels(dat$variable) == "OTHER.HARDWOOD"] = "OTHER HARDWOOD"
  
  dat$variable = factor(dat$variable, levels = sort(levels(dat$variable)))
  
  p <- ggplot(dat) + geom_line(data=dat, 
                               aes(x=radius, y=value, colour=factor(handle), linetype=factor(handle)), 
                               lwd=1.0) 
  p <- p + scale_colour_manual("Model", values=c('gray5', 'grey57', 'gray5', 'grey57')) 
  p <- p + scale_linetype_manual("Model", values = c("solid", "solid", "dashed", "dashed")) +
    xlab('Radius (km)') + ylab('Proportion of deposited pollen')
  p <- p + xlim(0, 800) + ylim(0,1.0)
  p <- p + facet_wrap(~variable, ncol=3) + theme_bw() + theme(panel.grid.major = element_blank(), 
                                                              panel.grid.minor = element_blank())
  print(p)
  
  ggsave(p, filename=sprintf('%s/%s/%s.pdf', wd, path_figs, 'kernel_discrete_cdfs'), width=12, height=12)
  
}

# do it!
disp_dat = dispersal_cdf(radius, runs, taxa, path_out)

disp_mean = disp_dat$dat
disp_all  = disp_dat$r_all
plot_dispersal_cdfs(disp_mean)

capture_radius <- function(dat, prop){
  
  handles = as.vector(unique(dat$handle))
  
  cr = data.frame(handle=character(0), taxon=character(0), r=numeric(0))
  for (handle in handles){
    for (taxon in unique(dat$variable)){
      dat_sub = dat[which((dat$handle == handle) & (dat$variable == taxon)), ]
      rad = dat_sub[which.min(abs(dat_sub[, 'value'] - prop)), 'radius']
      cr = rbind(cr, data.frame(handle=handle, taxon=taxon, r=rad))
    }
  }
  
  return(cr)
}

capture_25 = capture_radius(disp_mean, 0.25)
capture_50 = capture_radius(disp_mean, 0.5)
capture_70 = capture_radius(disp_mean, 0.7)
capture_90 = capture_radius(disp_mean, 0.9)

capture_radius_ci <- function(dat, prop, radius){
  
#   handles   = as.vector(unique(dat$handle))
  taxa_list = colnames(dat[[1]][[1]][,,1])
  ntaxa = length(colnames(dat[[1]][[1]][,,1]))
  
  cr = data.frame(handle=character(0), taxon=character(0), r_med=numeric(0), r_lb=numeric(0), r_ub=numeric(0))
  for (i in 1:length(dat)){
    for (j in 1:ntaxa){
      dat_sub = dat[[i]][[1]][,j,]
      cr_all    = apply(dat_sub, 2, function(x) radius[which.min(abs(x - prop))])
      cr_quants = quantile(cr_all, probs=c(0.025, 0.5, 0.975)) * 1e-3
      cr = rbind(cr, data.frame(handle = dat[[i]][[2]], 
                                taxon  = taxa_list[j], 
                                r_med  = cr_quants[2], 
                                r_lb   = cr_quants[1], 
                                r_ub   = cr_quants[3]))
    }
  }
  
  return(cr)
}

capture_50_ci = capture_radius_ci(disp_all, 0.5, radius)
capture_70_ci = capture_radius_ci(disp_all, 0.7, radius)
capture_90_ci = capture_radius_ci(disp_all, 0.9, radius)

# make latex-ready table
capture_50_gvar = capture_50_ci[capture_50_ci$handle == 'g_Kpsi_Kgamma',]
taxa_list = as.vector(capture_50_gvar[order(capture_50_gvar[,'r_med']),'taxon'])

handles = c('g_Kpsi_Kgamma', 'pl_Ka_Kgamma')
for (taxon in taxa_list){
  row = paste0(substr(taxon,1,1), tolower(substr(taxon,2, nchar(taxon))))
  for (handle in handles){
    r50 = capture_50_ci[((capture_50_ci$taxon==taxon) & (capture_50_ci$handle==handle)), 
                        c('r_med', 'r_lb', 'r_ub')]  
    r70 = capture_70_ci[((capture_70_ci$taxon==taxon) & (capture_70_ci$handle==handle)), 
                        c('r_med', 'r_lb', 'r_ub')]  
    r90 = capture_90_ci[((capture_90_ci$taxon==taxon) & (capture_90_ci$handle==handle)), 
                        c('r_med', 'r_lb', 'r_ub')] 
    #row = paste(row, r50, r90, collapse='&')
    row = paste(row, sprintf("%.0f & %3.0f & %3.0f) & %.0f & %3.0f & %3.0f) & %.0f & %3.0f & %3.0f)", 
                             r50[1], r50[2], r50[3], r70[1], r70[2], r70[3], r90[1], r90[2], r90[3]),
                sep='&')
  }
  row = paste0(row, ' \\')
  print(row)
}

################################################################################################################
# capture radius vs published sedimentation velocity
################################################################################################################

# use the 'best' model 
capture_25_sub = capture_25[capture_25$handle == 'Variable PL',]
capture_50_sub = capture_50[capture_50$handle == 'Variable PL',]
capture_90_sub = capture_90[capture_90$handle == 'Variable PL',]
capture_70_sub = capture_70[capture_70$handle == 'Variable PL',]
# capture_50_sub = capture_50[capture_50$handle == 'Variable G',]
# capture_90_sub = capture_90[capture_90$handle == 'Variable G',]

# # svs       = read.table("data/svs_meta.csv", sep=',', header=TRUE)
# svs = read.table("data/svs_meta2.csv", sep=',', header=TRUE)
# 
# svs_50 = cbind(svs, data.frame(stepps=capture_50_sub[match(svs$taxon, capture_50_sub$taxon),3]))
# svs_90 = cbind(svs, data.frame(stepps=capture_90_sub[match(svs$taxon, capture_90_sub$taxon),3]))
# 
# cr = 90
# 
# if (cr == 50){
#   svs = svs_50
# } else if (cr == 90){
#   svs = svs_90
# }
# 
# svs_melt = svs[!(svs$taxon %in% c('OTHER.HARDWOOD', 'OTHER.CONIFER')), ]
# svs_melt = droplevels(svs_melt)
# 
# levels(svs_melt$code)= c('Matthias (2014)', 'Jackson (1999)' )
# 
# # cols    = brewer.pal(12, 'Paired')
# taxa_sv = levels(svs_melt$taxon)
# # taxa_idx = match(taxa_sv, taxa)
# taxa_idx = match(taxa_sv, stepps_cols$taxa)
#                   
# # p <- ggplot() + geom_point(data=svs_melt, aes(x=sv, y=stepps, colour=taxon), alpha=0.9, size=5) + 
# #   scale_colour_manual(values=cols[taxa_idx]) + 
# #   xlab('SV from lit') + ylab('Capture radius STEPPS') +
# #   facet_wrap(~code, nrow=2)
# # print(p)
# # 
# # ggsave(p, filename=sprintf('%s/%s/%s.pdf', wd, path_figs, 'SVs_v2'), width=10, height=14)
# 
# p <- ggplot() + geom_point(data=svs_melt, aes(x=sv, y=stepps, colour=taxon, shape=code), alpha=0.9, size=5) + 
#   scale_colour_manual(values=as.vector(stepps_cols$cols[taxa_idx]), name='Taxon') + 
#   scale_shape_manual(values=c(16, 17), name='Reference') + 
#   xlab('Sedimentation velocity (m/s)') + ylab(paste0('STEPPS ', cr, '% capture radius (km)')) 
# print(p)
# ggsave(p, filename=sprintf('%s/%s/%s.pdf', wd, path_figs, paste0('SVs_', cr, '_single')), width=10)#, height=14)


# svs       = read.table("data/svs_meta.csv", sep=',', header=TRUE)
svs = read.table("data/svs_meta3.csv", sep=',', header=TRUE)

svs_25 = cbind(svs, data.frame(stepps=capture_25_sub[match(svs$taxon, capture_25_sub$taxon),3]))
svs_50 = cbind(svs, data.frame(stepps=capture_50_sub[match(svs$taxon, capture_50_sub$taxon),3]))
svs_90 = cbind(svs, data.frame(stepps=capture_90_sub[match(svs$taxon, capture_90_sub$taxon),3]))
svs_70 = cbind(svs, data.frame(stepps=capture_70_sub[match(svs$taxon, capture_70_sub$taxon),3]))

cr = 90

if (cr == 50){
  svs = svs_50
} else if (cr == 90){
  svs = svs_90
} else if (cr == 70){
  svs = svs_70
} else if (cr == 25){
  svs = svs_25
}


svs_melt = svs[!(svs$taxon %in% c('OTHER.HARDWOOD')), ]
svs_melt = droplevels(svs_melt)

colnames(svs_melt)[6] = 'code'
levels(svs_melt$code)= c('other', 'umw' )

levels(svs_melt$taxon)[levels(svs_melt$taxon) == "OTHER.CONIFER"] = "FIR"
svs_melt$taxon = factor(svs_melt$taxon,sort(levels(svs_melt$taxon)))
# cols    = brewer.pal(12, 'Paired')
taxa_sv = levels(svs_melt$taxon)
# taxa_idx = match(taxa_sv, taxa)
taxa_idx = match(taxa_sv, stepps_cols$taxa)

# svs_melt = svs_melt[svs_melt$code %in% c('umw'),]
# p <- ggplot() + geom_point(data=svs_melt, aes(x=sv, y=stepps, colour=taxon), alpha=0.9, size=5) + 
#   scale_colour_manual(values=as.vector(stepps_cols$cols[taxa_idx]), name='Taxon') + 
#   scale_shape_manual(values=c(16, 17), name='Reference') + 
#   xlab('Sedimentation velocity (m/s)') + ylab(paste0('STEPPS ', cr, '% capture radius (km)')) 
# print(p)



# p <- ggplot() + geom_point(data=svs_melt, aes(x=sv, y=stepps, colour=taxon), alpha=0.9, size=5) + 
#   scale_colour_manual(values=as.vector(stepps_cols$cols[taxa_idx]), name='Taxon') + 
#   #scale_shape_manual(values=c(16, 17), name='Reference') + 
#   xlab('Falling speed (m/s)') + ylab(paste0('STEPPS ', cr, '% capture radius (km)')) 
# print(p)
# ggsave(p, filename=sprintf('%s/%s/%s.pdf', wd, path_figs, paste0('SVs_', cr, '_all_single')), width=10)#, height=14)

# 
# # svs_melt = svs_melt[svs_melt$author %in% c('Eisenhut', 'Durham', 'Dyakowska', 'Bodmer'),]
# p <- ggplot() + geom_point(data=svs_melt, aes(x=sv, y=stepps, colour=taxon), alpha=0.9, size=5) + 
#   scale_colour_manual(values=as.vector(stepps_cols$cols[taxa_idx]), name='Taxon') + 
# #   scale_shape_manual(values=c(16, 17, 15, 18), name='Reference') +
#   xlab('Falling speed (m/s)') + ylab(paste0('STEPPS ', cr, '% capture radius (km)')) 
# print(p)

table(svs_melt$author)
# svs_melt = svs_melt[svs_melt$author %in% c('Eisenhut', 'Durham', 'Dyakowska', 'Bodmer'),]
# 
# svs_melt = svs_melt[!(svs_melt$author %in% c('Niklas')),]

svs_melt$shape[!(svs_melt$author %in% c('Niklas'))] = 'no'
svs_melt$shape[(svs_melt$author %in% c('Niklas'))] = 'yes'

p <- ggplot() + geom_point(data=svs_melt, aes(x=sv, y=stepps, colour=taxon), alpha=0.9, size=5) + 
  scale_colour_manual(values=as.vector(stepps_cols$cols[taxa_idx]), name='Taxon') + 
    scale_shape_manual(values=c(16, 17), name='Reference') +
  xlab('Falling speed (m/s)') + ylab(paste0('STEPPS ', cr, '% capture radius (km)')) 
print(p)
ggsave(p, filename=sprintf('%s/%s/%s.pdf', wd, path_figs, paste0('SVs_', cr, '_single')), width=10)#, height=14)


p <- ggplot() + geom_point(data=svs_melt, aes(x=sv, y=stepps, colour=taxon), alpha=0.9, size=5) + 
  scale_colour_manual(values=as.vector(stepps_cols$cols[taxa_idx]), name='Taxon') + 
  scale_shape_manual(values=c(16, 17, 15, 18), name='Reference') + 
  xlab('Falling speed (m/s)') + ylab(paste0('STEPPS ', cr, '% capture radius (km)')) +
  facet_wrap(~author)
print(p)

########################################################################################################################################
# phi plot
########################################################################################################################################

dat = data.frame(matrix(0, nrow=0, ncol=3))

for (run in runs){
  fname = sprintf('%s/%s.csv', path_out, run$suff_fit)
  system(sprintf('r/fixup.pl %s', fname))
  fit   = read_stan_csv(fname)
  post  = rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)
  
  col_names = colnames(post[,1,])
  par_names  = unlist(lapply(col_names, function(x) strsplit(x, "\\[")[[1]][1]))
  par_vals  = post[,1,which(par_names == 'phi')]
  
  par_mean  = colMeans(par_vals)
  par_lb    = apply(par_vals, 2, function(x) quantile(x, probs=0.025))
  par_ub    = apply(par_vals, 2, function(x) quantile(x, probs=0.975))
  
  par_stats = data.frame(name=taxa, mu=par_mean, lb=par_lb, ub=par_ub, Model=rep(run$handle, length(taxa)))
  
  dat = rbind(dat, par_stats)
}

levels(dat$Model) <- c("Base G", "Base PL", "Variable G", "Variable PL")
levels(dat$name)[levels(dat$name) == "OTHER.CONIFER"] = "FIR"
dat$name = factor(dat$name,sort(levels(dat$name)))
levels(dat$name)[levels(dat$name) == "OTHER.HARDWOOD"] = "OTHER HARDWOOD"


greys = brewer.pal(9, "Greys")[c(4,5,7,9)]

p <- ggplot(data=dat, aes(x=reorder(name, mu), y=mu, group=Model, colour=Model)) + 
  geom_point(size=4, position=position_dodge(width=0.5)) + 
  geom_errorbar(aes(ymin=lb, ymax=ub), width=.2, position=position_dodge(width=0.5)) + 
  scale_colour_manual(values=greys) +
  xlab("Taxon") + ylab(parse(text='phi')) +
  coord_flip() + theme_bw() + theme(axis.title.x=element_text(size=20), 
                                    axis.title.y=element_text(size=20), 
                                    axis.text.x=element_text(size=rel(1.3)),
                                    axis.text.y=element_text(size=rel(1.3)))

print(p)

ggsave(p, filename=sprintf('%s/%s/%s.pdf', wd, path_figs, 'phi'), width=12, height=8)

###############################################################################################################
# compare with PPEs from lit
ppe_stepps = dat[which(dat$Model == 'Variable PL'),]
ppe_stepps[,2:4] = ppe_stepps[,2:4]/max(ppe_stepps[,2], na.rm=TRUE)

# ppes     = read.table("data/ppes_meta2.csv", sep=',', header=TRUE, stringsAsFactors=FALSE)[,c(1,2,6)]
ppes     = read.table("data/ppes_meta2.csv", sep=',', header=TRUE, stringsAsFactors=FALSE)
ppes$taxon[which(ppes$taxon == "OTHER.CONIFER")]  = "FIR"
ppes$taxon[which(ppes$taxon == "OTHER.HARDWOOD")] = "OTHER HARDWOOD"

ppe_melt = data.frame(ppes, stepps=ppe_stepps[match(ppes$taxon, ppe_stepps$name), 'mu'])
ppe_melt = data.frame(ppe_melt[!(ppe_melt$taxon %in% c('OTHER.HARDWOOD')), ])

tags = unique(ppe_melt$tag)
for (tag in tags){
  idx = which(ppe_melt$tag == tag)
  max_ppe = max(ppe_melt[idx,'ppe'], na.rm=TRUE)
  ppe_melt[idx, 'ppe'] = ppe_melt[idx, 'ppe']/max_ppe
}

ppe_melt$taxon = factor(ppe_melt$taxon, levels=unique(ppe_melt$taxon))
ppe_melt$model = factor(ppe_melt$model, levels=unlist(unique(ppe_melt$model)[c(3,1,2)]))
ppe_melt$location_grouped = factor(ppe_melt$location_grouped, levels=unlist(unique(ppe_melt$location_grouped)[c(2,3,1)]))

# cols_ppe       = cols[match(levels(ppe_melt$taxon), taxa)]
levels(stepps_cols$taxa)[levels(stepps_cols$taxa) == "OTHER.CONIFER"]  = "FIR"
levels(stepps_cols$taxa)[levels(stepps_cols$taxa) == "OTHER.HARDWOOD"]  = "OTHER HARDWOOD"
ppe_melt$taxon = factor(ppe_melt$taxon,sort(levels(ppe_melt$taxon)))
taxa_idx = match(levels(ppe_melt$taxon), stepps_cols$taxa)

p <- ggplot() + geom_point(data=ppe_melt, aes(x=ppe, y=stepps, colour=taxon, shape=model), size=3) + 
  scale_colour_manual(values=as.vector(stepps_cols$cols[taxa_idx]), name="Taxon") + xlab('PPE lit') + ylab('PPE STEPPS') + 
  scale_shape_manual(name="Model", values=c(15,16,17)) + coord_fixed() +
  #   geom_text(data = ppe_melt, aes(x=value,y=ppe_stepps, label = taxon), position=position_jitter(h=1,w=1)) + #hjust = 0.8, vjust=-1.5) +
#   facet_wrap(~tag, nrow=2)
  facet_wrap(~location_grouped, nrow=3)
print(p)

ggsave(p, filename=sprintf('%s/%s/%s.pdf', wd, path_figs, 'PPEs_panels'), width=12)#, height=10)

################################################################################################################
# 
# # p <- ggplot() + geom_point(data=ppe_compiled, aes(x=ppe_s, y=ppe_stepps), size=3) + 
# #      geom_text(data = ppe_compiled, aes(x=ppe_s,y=ppe_stepps, label = taxon), hjust = 0.8, vjust=-1.5)
# # print(p)
# # direct.label(p, first.qp)
# # 
# # plot(ppes$marquer, ppe_stepps$mu)
# 
# # gamma plot
# dat = data.frame(matrix(0, nrow=0, ncol=3))
# 
# for (run in runs){
#   fname = sprintf('%s/%s.csv', path_out, run$suff_fit)
#   system(sprintf('r/fixup.pl %s', fname))
#   fit   = read_stan_csv(fname)
#   post  = rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)
#   
#   col_names = colnames(post[,1,])
#   par_names = unlist(lapply(col_names, function(x) strsplit(x, "\\[")[[1]][1]))
#   par_idx   = which(par_names == 'gamma')
#   par_vals  = post[,1,par_idx]
#   
#   if (length(par_idx)>1){
#     par_mean  = colMeans(par_vals)
#     par_lb    = apply(par_vals, 2, function(x) quantile(x, probs=0.025))
#     par_ub    = apply(par_vals, 2, function(x) quantile(x, probs=0.975))
#   } else {
#     par_mean  = rep(mean(par_vals), K)
#     par_lb    = rep(quantile(par_vals, probs=0.025), K)
#     par_ub    = rep(quantile(par_vals, probs=0.975), K)
#   }
#   
#   
#   par_stats = data.frame(name=taxa, mu=par_mean, lb=par_lb, ub=par_ub, handle=rep(run$handle, length(taxa)))
#   
#   dat = rbind(dat, par_stats)
# }
# 
# p <- ggplot(data=dat, aes(x=reorder(name, mu), y=mu, group=handle, colour=handle)) + 
#   geom_point(size=4, position=position_dodge(width=0.5)) + 
#   geom_errorbar(aes(ymin=lb, ymax=ub), width=.2, position=position_dodge(width=0.5)) + 
#   xlab("Taxon") + ylab(parse(text='gamma')) +
#   coord_flip() + theme_bw() + theme(axis.title.x=element_text(size=20), 
#                                     axis.title.y=element_text(size=20), 
#                                     axis.text.x=element_text(size=rel(1.3)),
#                                     axis.text.y=element_text(size=rel(1.3)))
# 
# ggsave(p, filename=sprintf('%s/%s/%s.pdf', wd, path_figs, 'gamma'), width=12, height=12)