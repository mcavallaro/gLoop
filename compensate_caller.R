#!/usr/bin/env r

# launch recursively the script to compensate flow cytometer reads.
# This is Free Software - You can use and distribute it under
# the terms of the GNU General Public License, version 3 or later
# (c) Massimo Cavallaro (m.cavallaro@warwick.ac.uk)

dirs<-list.dirs(recursive = FALSE)

for (path in dirs){ 
  writeLines(path)
  source("CompensateFlowSet.R")
}
