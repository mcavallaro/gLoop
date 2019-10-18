#!/usr/bin/env r

# This is Free Software - You can use and distribute it under
# the terms of the GNU General Public License, version 3 or later
# (c) Massimo Cavallaro (m.cavallaro@warwick.ac.uk)

lim<-function(f, f1) {
  Min<-min(exprs(f)[,f1])
  Max<-max(exprs(f)[,f1])
  return(c(Min, Max))
}

DNA.threshold<-500

library("flowCore")
library("flowClust")

## decomment to use with littler:
#file.name<-argv
#file.name<-'./HBB_0742016/Experiment_pA- 80_007.fcs'

writeLines('########### 2')
writeLines(file.name)
writeLines('########### 2')

ff<-read.FCS(file.name)
file.base.name<-strsplit(base.name, '.fcs')


## 1st-step: remove reads out of scale
varNames<-colnames(ff)
rectGate<-rectangleGate(filterId = "Fluorescence Region", "FSC-A" = c(0, 200000))
ff1<-Subset(ff, rectGate)


## 2nd step: remove reads corresponding to dead cells or duplets/triplets
res2<-flowClust(ff1, varNames = varNames[c(2, 3, 5, 6)], K = 2, B = 1000)
level<-0.4
zc<-0.9
ruleOutliers(res2)<-list(level = level)
ruleOutliers(res2)<-list(z.cutoff = zc)

## 2.1: split the flowFrame and select the chunk
## with intermediate size (FSC-A) and granularity (SSC-A)
ff1.splitted<-split(ff1, res2, level = level, z.cutoff = zc, population = list(sc1 = 1, sc2 = 2))

f1<-summary(ff1)[c(4,22)]
f1.sc1<-summary(ff1.splitted$sc1)[c(4,22)]
f1.sc2<-summary(ff1.splitted$sc2)[c(4,22)]

dists<-c(dist(rbind(f1, f1.sc1)), dist(rbind(f1, f1.sc2)))

switch(which.min(dists),
  {ff2<-ff1.splitted$sc1; ff2_<-ff1.splitted$sc2},
  {ff2<-ff1.splitted$sc2; ff2_<-ff1.splitted$sc1}
  )


## 2.2 plot the clusters
sub.dir<-'G1G2__'
dir.create(file.path(dir.name, sub.dir), showWarnings = FALSE)
png(filename = file.path(dir.name, sub.dir, paste(file.base.name, "_k_2.png", sep = '')),
  width = 480*1.5, height = 680*1.5)
par(mfrow = c(4,2), mar = c(4,4,3,2) + 0.1)

plot(res2, data=ff1, subset = c(1, 3), level = level, z.cutoff = zc,
     ylim = c(0, 9000), xlim = c(0, 240000))
plot(res2, data=ff1, subset = c(2, 4), level = level, z.cutoff = zc,
     ylim = c(0, 140000), xlim = c(60000, 160000))
# plot(ff2, c("FSC-A", "SSC-A"), smooth=FALSE, ylim=c(0,9000), xlim=c(0,240000))
# plot(ff2, c("FSC-W", "SSC-W"), smooth=FALSE, ylim=c(50000,140000), xlim=c(60000,160000))

xlim<-lim(ff2_, 'R640-670/14-A')
ylim<-lim(ff2_, 'UV355-450/50-A')
plot(ff2_, c('R640-670/14-A', 'UV355-450/50-A'), smooth = FALSE, main = 'discarted', ylim = ylim, xlim = xlim)
xlim<-lim(ff2, 'R640-670/14-A')
ylim<-lim(ff2, 'UV355-450/50-A')
plot(ff2, c('R640-670/14-A', 'UV355-450/50-A'), smooth = FALSE, main = 'chosen', ylim = ylim, xlim = xlim)
abline(h = DNA.threshold)

xlim<-lim(ff2_, 'FSC-A')
ylim<-lim(ff2_, 'R640-670/14-A')
plot(ff2_, c('FSC-A', 'R640-670/14-A'), smooth = FALSE, main = 'discarted', ylim = ylim, xlim = xlim)
xlim<-lim(ff2, 'FSC-A')
ylim<-lim(ff2, 'R640-670/14-A')
plot(ff2, c('FSC-A', 'R640-670/14-A'), smooth = FALSE, main = 'chosen', ylim = ylim, xlim = xlim)

xlim<-lim(ff2_, 'FSC-A')
ylim<-lim(ff2_, 'UV355-450/50-A')
plot(ff2_, c('FSC-A', 'UV355-450/50-A'), smooth = FALSE, main = 'discarted', ylim = ylim, xlim = xlim)
xlim<-lim(ff2, 'FSC-A')
ylim<-lim(ff2, 'UV355-450/50-A')
plot(ff2, c('FSC-A', 'UV355-450/50-A'), smooth = FALSE, main = 'chosen', ylim = ylim, xlim = xlim)
abline(h = DNA.threshold)

dev.off()


## 3rd step: remove reads with no DNA UV355-450/50-A
rectGate<-rectangleGate(filterId="DNA.threshold",
  "UV355-450/50-A"=c(DNA.threshold, max(exprs(ff2)[,'UV355-450/50-A'])))
ff3<-Subset(ff2, rectGate)


## 4th step: export
write.csv(ff3.data.frame,
  file = file.path(dir.name, sub.dir, paste(file.base.name, "_k_2.csv", sep = ""))
  )

write.table(ff3.data.frame$"R640-670/14-A",
  file = file.path(dir.name, sub.dir, paste(file.base.name, "_k_2_mRNA.csv", sep = "")),
  row.names = FALSE, col.names = FALSE)

writeLines(length(ff3.data.frame$"R640-670/14-A"))
