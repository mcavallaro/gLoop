#!/usr/bin/env r

library(flowCore)
library(flowViz)

## decomment to use with littler:
#path <- argv

sampleFiles<-grep('.fcs', list.files(path), value = TRUE)
sampleNames<-sapply(strsplit(sampleFiles, "_"), "[", 2)

fs<-read.flowSet(path = path, pattern = ".fcs")

# set new sample names
fs@phenoData@data$name<-sampleNames 

pathCompensation<-paste(path, 'Compensation', sep='/')
# comp. controls had to be re-named so we can generate the correct matrix
fs.compensation<-read.flowSet(path = pathCompensation, pattern = ".fcs") 

comp.data<-fs.compensation[c(1,2,3), c("FSC-A", "SSC-A", "R640-670/14-A", "UV355-450/50-A")]

# generate spillover matrix
spill.mat<-spillover(comp.data, unstained = sampleNames(comp.data)[2], fsc = "FSC-A", ssc = "SSC-A")

# generate a new compensated flowSet
fs.compensated<-compensate(fs, spill.mat)

write.flowSet(fs.compensated, outdir = paste(path, 'Compensated', sep = '/'))
