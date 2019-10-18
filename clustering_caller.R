#!/usr/bin/env r

# launch recursively the scripts to gate the flow cytometer with flowClust
# and export single-cell mRNA reads.
# This is Free Software - You can use and distribute it under
# the terms of the GNU General Public License, version 3 or later
# (c) Massimo Cavallaro (m.cavallaro@warwick.ac.uk)


# % This code can be redistributed and/or modified under the terms of the 
# % GNU General Public License as published by the Free Software Foundation, 
# % either version 3 of the License, or (at your option) any later version.
# % This program is distributed by the author in the hope that it will be 
# % useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
# % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# % Please cite:
# % (c) MC

dirs<-list.dirs(recursive = TRUE)
for (i in 1:length(dirs)){
  if(grepl("Compensated", dirs[i])){
        files<-list.files(path = dirs[i], pattern = '.fcs')
        for (j in 1:length(files)){
                file.name<-paste(dirs[i], files[j], sep = '/')
                base.name<-files[j]
                dir.name<-dirs[i]
                # writeLines(file.name)
                # writeLines(base.name)
                # writeLines(dir.name)
                source("clustering_2.R")
                source("clustering_3.R")
              }
            }
          }
