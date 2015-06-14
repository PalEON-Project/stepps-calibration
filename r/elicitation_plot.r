#####################################################################################
# user pars
#####################################################################################
library(neotoma)
# library(RColorBrewer)
library(analogue)

elicit = read.csv('../stepps-data/data/pollen_meta_2014-05-01_compiled.csv', header=TRUE, stringsAsFactor=FALSE)

# example of clear settlement horizon: datasetid 314, Brownsby
# example of unclear settlement horizon: datasetid 1884, Pennington
ids = c(314, 1884)

sh = elicit[(elicit$datasetID %in% ids), ]

dat <- get_download(ids)

for (i in 1:length(ids)){
  
  x = dat[[i]]
  
  age = x$sample.meta$age
  
  good.samps <- which(age < 2000)
  
  pol <- x$counts[good.samps,]
  pol.counts <- compile_taxa(pol, 'WhitmoreSmall')
  pol <- pol.counts / rowSums(pol.counts)
  pol <- pol[ ,!colnames(pol) %in% 'Other']
  
  markers = c('AMBROSIA', 'RUMEOXYR', 'POACEAE')
  
  most.common <- colnames(pol)[rank(colMeans(pol)) > (ncol(pol) - 10)]
  most.common <- unique(c(most.common, markers[markers  %in% colnames(pol)]))
  
  core.dat <- data.frame(Depth = x$sample.meta$depth[good.samps], pol.counts[,most.common])
  
#   Zones <- c()
#   Stratiplot(Age ~ . - Depth, data = chooseTaxa(abernethy, n.occ = 5, max.abun = 10),
#              type = c("h","l","g"), sort = "wa", zones = Zones,
#              zoneNames = c(LETTERS[1:6]))

  col1 ="steelblue3"
  col2 = "orange" 
#   col2 = "darkgrey"
  Zones <- core.dat$Depth[unique(unlist(sh[i, 9:12]))]

  fname=paste0('figures/', x$dataset$dataset.meta$collection.handle, '_', ids[i] )
  pdf(file=paste0(fname, ".pdf"), width=10)
  print(Stratiplot(Depth ~ . , data = chooseTaxa(core.dat, n.occ=5),
             type = c("h","l","g"), sort = "var", drawLegend = FALSE, 
             zones=Zones, zoneNames="", col.zones=col2, 
             xlab="Pollen percentages", lwd.zones = 2.5, lty.zones = 2,
             col.line = col1,
             col.symbol = 'white',
             col.refline = 'white',
             col.smooth = 'grey',
             col.poly = 'grey'))
  dev.off()

  sys_str=paste0("pdftk ", fname, ".pdf cat 2 output ", fname, "2.pdf") 
  system(sys_str)
  sys_str=paste0("mv ", fname, "2.pdf ", fname, ".pdf")
  system(sys_str)
  sys_str=paste0("gs -sDEVICE=png16m -r100 -o ", fname, ".png ", fname, ".pdf")
  system(sys_str)

}


# data(abernethy)
# 
# pdf(file='test.pdf', width=10)
# print(Stratiplot(Age ~ . - Depth, data = chooseTaxa(abernethy, n.occ = 5, max.abun = 10),
#                  type = c("h","l","g")))
# dev.off()
