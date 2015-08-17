plot_pollen_preds_paper <- function(dat, path_figs, type) {
  
#   levels(dat$taxon)[levels(dat$taxon) == "OTHER.CONIFER"] = "OTHER CONIFER"
  levels(dat$taxon)[levels(dat$taxon) == "OTHER.CONIFER"] = "FIR"
  levels(dat$taxon)[levels(dat$taxon) == "OTHER.HARDWOOD"] = "OTHER HARDWOOD"
  
  dat$taxon = factor(dat$taxon, levels = sort(levels(dat$taxon)))
  
  p <- ggplot(dat) + geom_point(data=dat, mapping=aes(x=veg, y=data, shape=3), colour='red') 
  p <- p + scale_shape_identity()
  p <- p + geom_point(data=dat, mapping=aes(x=pred, y=data, shape=16), colour="gray21") + 
    scale_colour_manual(values=c("lightgrey"), guide=FALSE)
  # p <- p + scale_alpha_manual(values=c(0,0.5))
  p <- p + geom_abline(intercept=0, slope=1, colour="grey56", linetype="dashed")
  p <- p + xlab("veg or predicted pollen") + ylab("raw pollen")
  p <- p + xlim(0,1) + ylim(0,1) + coord_fixed()
  p <- p + facet_wrap(~taxon)
  p <- p + theme_bw()
  p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.margin = unit(0.4, "lines"))
  print(p)
  
  ggsave(p, filename=sprintf('%s/%s/%s%s.pdf', getwd(), path_figs, 'pollen_preds_', type), width=12, scale=1)
  
}

plot_pollen_phi_scaled_paper <- function(dat, path_figs, type) {
  
#   levels(dat$taxon)[levels(dat$taxon) == "OTHER.CONIFER"] = "OTHER CONIFER"
  levels(dat$taxon)[levels(dat$taxon) == "OTHER.CONIFER"] = "FIR"
  levels(dat$taxon)[levels(dat$taxon) == "OTHER.HARDWOOD"] = "OTHER HARDWOOD"
  
  dat$taxon = factor(dat$taxon, levels = sort(levels(dat$taxon)))
  
  # plot the phi scaled veg versus 
  p <- ggplot(dat) + geom_point(data=dat, mapping=aes(x=veg, y=data, shape=3), colour='red') 
  p <- p + scale_shape_identity()
  p <- p + geom_point(data=dat, mapping=aes(x=phi_pred, y=data, shape=16), colour="gray21") + 
    scale_colour_manual(values=c("black"), guide=FALSE)
  p <- p + geom_abline(intercept=0, slope=1, colour="grey56", linetype="dashed")
  p <- p + xlab("veg or predicted pollen") + ylab("raw pollen")
  p <- p + xlim(0,1) + ylim(0,1) + coord_fixed()
  p <- p + facet_wrap(~taxon)
  p <- p + theme_bw()
  p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.margin = unit(0.4, "lines"))
  print(p)
  
  ggsave(p, filename=sprintf('%s/%s/%s%s.pdf', getwd(), path_figs, 'pollen_phi_scaled_', type), width=12, scale=1)
}


# plot_pollen_preds_paper <- function(dat, path_figs, type) {
#   
#   p <- ggplot(dat) + geom_point(data=dat, mapping=aes(x=veg, y=data, shape=3), colour='red') 
#   p <- p + scale_shape_identity()
#   p <- p + geom_point(data=dat, mapping=aes(x=pred, y=data, colour=run, shape=16), alpha=0.6) + 
#     scale_colour_manual(values=c("black", "deepskyblue2"), guide=FALSE)
#   # p <- p + scale_alpha_manual(values=c(0,0.5))
#   p <- p + geom_abline(intercept=0, slope=1, colour="grey56", linetype="dashed")
#   p <- p + xlab("veg or predicted pollen") + ylab("raw pollen")
#   p <- p + xlim(0,1) + ylim(0,1) + coord_fixed()
#   p <- p + facet_wrap(~taxon)
#   p <- p + theme_bw()
#   p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.margin = unit(0.4, "lines"))
#   print(p)
#   
#   ggsave(p, filename=sprintf('%s/%s/%s%s.pdf', getwd(), path_figs, 'pollen_preds_', type), width=12, scale=1)
#   
# }
# 
# 
# plot_pollen_phi_scaled_paper <- function(dat, path_figs, type) {
#   # plot the phi scaled veg versus 
#   p <- ggplot(dat) + geom_point(data=dat, mapping=aes(x=veg, y=data, shape=3), colour='red') 
#   p <- p + scale_shape_identity()
#   p <- p + geom_point(data=dat, mapping=aes(x=phi_pred, y=data, colour=run, shape=16), alpha=0.6) + 
#     scale_colour_manual(values=c("black", "deepskyblue2"), guide=FALSE)
#   p <- p + geom_abline(intercept=0, slope=1, colour="grey56", linetype="dashed")
#   p <- p + xlab("veg or predicted pollen") + ylab("raw pollen")
#   p <- p + xlim(0,1) + ylim(0,1) + coord_fixed()
#   p <- p + facet_wrap(~taxon)
#   p <- p + theme_bw()
#   p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.margin = unit(0.4, "lines"))
#   print(p)
#   
#   ggsave(p, filename=sprintf('%s/%s/%s.pdf', getwd(), path_figs, 'pollen_phi_scaled_', type), width=12, scale=1)
# }
