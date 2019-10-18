#loopscore functions
#written by Matt jones, Contact m.jones.18@warwick.ac.uk
filter_gene_list2=function(gene_list){
  require(dplyr)
gene_list = gene_list %>% group_by(txStart) %>% filter(n()==1)
gene_list = gene_list %>% group_by(txEnd) %>% filter(n()==1)


 #filter to remove chr in the list of chromosomes - just keep numebr
  gene_list$chrom=sub('chr', '',gene_list$chrom)

  #sort by chromosome and only inlcude 1-22 and x (not m ect) 
  gene_list2=gene_list[which(gene_list$chrom==1),]
  for (i in c(seq(from = 2, to = 22, by = 1),'X')) {
    gene_list2 <- rbind(gene_list2, gene_list[which(gene_list$chrom==i),])
  }

  gene_list<- gene_list2
  return(gene_list)
}

filter_gene_list=function(gene_list){
  require(dplyr)
  gene_list <- gene_list %>% group_by(name2) %>% filter(n()==1)
  #filter to remove chr in the list of chromosomes - just keep numebr
  gene_list$chrom=sub('chr', '',gene_list$chrom)
  
  #sort by chromosome and only inlcude 1-22 and x (not m ect)  gene_list2=gene_list[which(gene_list$chrom==1),]
  for (i in c(seq(from = 2, to = 22, by = 1),'X')) {
    gene_list2 <- rbind(gene_list2, gene_list[which(gene_list$chrom==i),])
  }
  
  gene_list<- gene_list2
  return(gene_list)
  
}
  initialize_bedgraph = function(bedgraph) {
    chrbins=matrix(c(seq(1:22),'X',rep(NA,46)), nrow=23)
    count=0
    colnames(chrbins) <- c('chrom','startbin','endbin')
    for (i in c(seq(1:22),'X')){
      count=count+1
      bedgraph_one_chr = bedgraph[which(bedgraph$V1==i),]
      
      first = min(bedgraph_one_chr$V4)
      last = max(bedgraph_one_chr$V4)
      chrbins[count,2]=bedgraph$V4[first]
      chrbins[count,3]=bedgraph$V4[last]
      
    }
    return(chrbins)
  }
  
 calculate_average_interactions = function(chia_pet_matrix,chrbins,resolution) {
   averageinteractions_all_chrom=c(NA, NA,NA)
  for( chromosome in  c(seq(1:22),'X')) {
    startbin= as.numeric(chrbins[which(chrbins[,1]==chromosome),2])
    endbin= as.numeric(chrbins[which(chrbins[,1]==chromosome),3])
    chr_pet <- chia_pet_matrix[which(chia_pet_matrix$V1 < endbin),]
    chr_pet <- chr_pet[which(chr_pet$V1> startbin),]
    chr_pet <- chr_pet[which(chr_pet$V2<endbin),]
    chr_pet <- chr_pet[which(chr_pet$V2> startbin),]
    
    genomic_distance <- (chr_pet$V2-chr_pet$V1)*resolution
    chr_distance <- cbind(chr_pet, genomic_distance)
    number_of_bins = (as.numeric(endbin)-as.numeric(startbin))
    distances = seq(from=0, to=number_of_bins*resolution, by = resolution)
    averageinteractions=c()
    for (i in 1:2100) {
      interactions=c()
      interactions = chr_distance[['V3']][which(chr_distance[['genomic_distance']]==distances[i])]
      number_of_zeros = number_of_bins-(length(interactions)+i) #the i is important, draw a matrix to realize
      interactions = c(interactions, rep(0, number_of_zeros))
      meaninteractions=mean(interactions)
      averageinteractions=c(averageinteractions, meaninteractions) 
    }
    averageinteractions=cbind(rep(chromosome, length(averageinteractions)), averageinteractions, distances[1:2100])
    averageinteractions_all_chrom <- rbind(averageinteractions_all_chrom, averageinteractions )
    
  }
  
  averageinteractions_all_chrom=na.omit(averageinteractions_all_chrom)
  
  return(averageinteractions_all_chrom)}
 
 calculate_loopscores = function (chia_pet_matrix, chrbins, gene_list) {
   loopscoreslist=matrix(c(NA,NA,NA,NA,NA),nrow=1)
   colnames(loopscoreslist) <- c('gene_name','chrom','txStart','txEnd','score')
   
   for (gene in 1:nrow(gene_list)){
     
     gene_name = as.character(gene_list$name2[gene])
     chrom=gene_list$chrom[gene]
     startbin = as.numeric(chrbins[which(as.character(chrbins[,1])==chrom),2]) + floor(gene_list$txStart[gene]/resolution) #note floor here, because of the start bin
     endbin = as.numeric(chrbins[which(as.character(chrbins[,1])==chrom),2]) + floor(gene_list$txEnd[gene]/resolution) #note floor here, because of the start bin
     instart = chia_pet_matrix[which(chia_pet_matrix$V1==startbin),]
     genescore=instart$V3[which(instart$V2==endbin)]
     
     if(length(genescore) > 0) {
       
       loopscoreslist=rbind(loopscoreslist,c(gene_name,chrom,as.numeric(gene_list$txStart[gene]),as.numeric(gene_list$txEnd[gene]),genescore))
     }
     if (length(genescore)==0) {
       genescore = 0
       loopscoreslist=rbind(loopscoreslist,c(gene_name,chrom,as.numeric(gene_list$txStart[gene]),as.numeric(gene_list$txEnd[gene]),genescore))
       
     }
   }
   return(loopscoreslist)}
  
 
 size_normalise_loopscore = function(loopscoreslist, averageinteractions_all_chrom) {
   loopscoreslist <- na.omit(loopscoreslist)
   adjusted_loopscore=c()
   #now to calculate the  loopscore size normalized
   for( chromosome in  c(seq(1:22),'X')) {
     averageinteractions <- averageinteractions_all_chrom[which(averageinteractions_all_chrom[,1] == chromosome),2]
     chr_loopscore_list = loopscoreslist[which(loopscoreslist[['chrom']]==chromosome),]
     
     startbin = as.numeric(chrbins[which(as.character(chrbins[,1])==chromosome),2]) + floor(as.numeric(chr_loopscore_list$txStart)/resolution) #note floor here, because of the start bin
     endbin = as.numeric(chrbins[which(as.character(chrbins[,1])==chromosome),2]) + floor(as.numeric(chr_loopscore_list$txEnd)/resolution) #
     binsize = endbin-startbin
     expectedloopscores = averageinteractions[binsize+1] #note the plus 1 is because if binsize = 0 we want the first row of average interactions
     chr_loopscore_list$score = as.numeric(chr_loopscore_list$score)/as.numeric(expectedloopscores)
     chr_loopscore_list <- cbind(chr_loopscore_list, expectedloopscores)
     adjusted_loopscore = rbind(adjusted_loopscore,chr_loopscore_list)
   }
   return(adjusted_loopscore)
 }
 
 compare_loop_to_noise = function(adjusted_loopscore, readcounts) {
   loopscoreslist<- adjusted_loopscore
   loopscoreslist=na.omit(loopscoreslist)
   noise_compared_loop=na.omit(loopscoreslist)
   means = rep(NA, nrow(noise_compared_loop))
   cv = rep(NA, nrow(noise_compared_loop))
   burst_freq = rep(NA, nrow(noise_compared_loop))
   burst_size = rep(NA, nrow(noise_compared_loop))
   
   noise_compared_loop = cbind(noise_compared_loop, means)
   noise_compared_loop = cbind(noise_compared_loop, cv)
   noise_compared_loop = cbind(noise_compared_loop, burst_freq)
   noise_compared_loop = cbind(noise_compared_loop, burst_size)
   
   for (i in 1:nrow(loopscoreslist)) {
     gene_name=as.character(loopscoreslist[i,1])
     whatgene = which(readcounts[,1]==gene_name)
     if (any(whatgene)==T) {
       variance = var(as.numeric(na.omit(as.numeric(readcounts[whatgene,2:ncol(readcounts)]))))
       meanreads=mean(as.numeric(na.omit(as.numeric(readcounts[whatgene,2:ncol(readcounts)]))))
       if (is.na(meanreads)==F){
         if (meanreads==0) {CVs=NA}
         if (abs(meanreads)>0) { CVs=sqrt(variance)/meanreads }
         noise_compared_loop[i,7]=meanreads
         noise_compared_loop[i,8]=CVs
         noise_compared_loop[i,9]=meanreads^2/(variance-meanreads)
         noise_compared_loop[i,10]=(variance-meanreads)/meanreads}
       if (any(whatgene)==F) {
         noise_compared_loop[i,7]=NA
         noise_compared_loop[i,8]=NA
         noise_compared_loop[i,9]=NA
         noise_compared_loop[i,10]=NA}
       
       #print(i/nrow(chr_loopscores))
       
     }}
   res <- na.omit(noise_compared_loop )
   return(res)
 }
 
 #new function (massimo's choice) to size adjust the loopscore  have not done this yet!
 #davids idea = get v3 as vector
 #make the rownames = paste the two numbers together we are looking for -> go straight to that
 calculate_and_adjust_loopscores = function (chia_pet_matrix, chrbins, gene_list,no_genes, resolution) {
   loopscoreslist=matrix(c(NA,NA,NA,NA,NA),nrow=1)
   colnames(loopscoreslist) <- c('gene_name','chrom','txStart','txEnd','score')
  list_of_interactions = chia_pet_matrix$V3
  my_rownames = paste(chia_pet_matrix$V1, chia_pet_matrix$V2, sep='_')
  rownames(chia_pet_matrix) <- my_rownames
  for (gene in 1:100) { #nrow(gene_list)){
     
     gene_name = as.character(gene_list$name2[gene])
     chrom=gene_list$chrom[gene]
     gene_size = gene_list$txEnd[gene]-gene_list$txStart[gene]
     startbin = as.numeric(chrbins[which(as.character(chrbins[,1])==chrom),2]) + floor(gene_list$txStart[gene]/resolution) #note floor here, because of the start bin
     endbin = as.numeric(chrbins[which(as.character(chrbins[,1])==chrom),2]) + floor(gene_list$txEnd[gene]/resolution) #note floor here, because of the start bin
    loopscore_row_name=paste(as.character(startbin), '_', as.character(endbin), sep='')
    genescore = chia_pet_matrix[loopscore_row_name,3]
     if(is.na(genescore)==F) {

      random_start_locations = sample(seq(from=as.numeric(chrbins[which(as.character(chrbins[,1])==chrom),2]),to=as.numeric(chrbins[which(as.character(chrbins[,1])==chrom),3]) -  ceiling(gene_size/resolution)),no_genes, replace =T)
    random_end_locations = random_start_locations + ceiling(gene_size/resolution)
    random_loopscore_rownames = paste(as.character(random_start_locations), '_', as.character(random_end_locations), sep='')
    random_scores = chia_pet_matrix[random_loopscore_rownames,3]
    random_scores[is.na(random_scores)] <- 0
    adjusted_genescore=genescore/mean(random_scores)
    loopscoreslist=rbind(loopscoreslist,c(gene_name,chrom,as.numeric(gene_list$txStart[gene]),as.numeric(gene_list$txEnd[gene]),adjusted_genescore))
    
     }
    if(is.na(genescore)==T) {
      adjusted_genescore=0
      loopscoreslist=rbind(loopscoreslist,c(gene_name,chrom,as.numeric(gene_list$txStart[gene]),as.numeric(gene_list$txEnd[gene]),adjusted_genescore))
      
    }}
   return(loopscoreslist)
    
 }

