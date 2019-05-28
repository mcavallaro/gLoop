#!/usr/bin/env r
# launch recursively the script to compensate flow cytometer reads.


dirs<-list.dirs(recursive = FALSE)

for (path in dirs){ 
  writeLines(path)
  source("CompensateFlowSet.R")
}
