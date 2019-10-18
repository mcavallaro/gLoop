 rm(list=ls())
# written by Matt Jones, contact m.jones.18@warwick.ac.uk
library(dplyr)
library(ggplot2)
source(file = '/home/u1762230/chia_pet/K562/rep_2_SRR372748_hg19/R_analysis/loopscore_functions.r')
chia_pet_matrix=read.table('/home/u1762230/chia_pet/K562/rep_2_SRR372748_hg19/rep_2_SRR372748_hg19_hic_matrix_res_2000.matrix')
bedgraph = read.table('/home/u1762230/chia_pet/K562/rep_2_SRR372748_hg19/rep_2_SRR372748_hg19_hic_matrix_res_2000_abs.bed')
bedgraph$V1=substr(bedgraph$V1,4,10) #onlydo if bedgraph has letters in
gene_list=read.table('/home/u1762230/chia_pet/K562/human_gene_list_hg19', header = T)
resolution=2000
gene_list = filter_gene_list(gene_list)
chrbins=initialize_bedgraph(bedgraph)


#massimos idea for no binning
no_genes=10000
 loopscoreslist=matrix(c(NA,NA,NA,NA,NA),nrow=1)
   colnames(loopscoreslist) <- c('gene_name','chrom','txStart','txEnd','score')
  list_of_interactions = chia_pet_matrix$V3
  my_rownames = paste(chia_pet_matrix$V1, chia_pet_matrix$V2, sep='_')
  rownames(chia_pet_matrix) <- my_rownames
   chrom_sizes=read.table('/home/u1762230/chia_pet/K562/hg19.chrom.sizes')
  chrom_sizes$V1=substr(chrom_sizes$V1,4,10) #onlydo if bedgraph has letters in
  chrom_sizes=chrom_sizes[1:23,]
for (gene in 1:nrow(gene_list)){

     gene_name = as.character(gene_list$name2[gene])
     chrom=gene_list$chrom[gene]
     gene_size = gene_list$txEnd[gene]-gene_list$txStart[gene]
     startbin = as.numeric(chrbins[which(as.character(chrbins[,1])==chrom),2]) + floor(gene_list$txStart[gene]/resolution) #note floor here, because of the start bin
     endbin = as.numeric(chrbins[which(as.character(chrbins[,1])==chrom),2]) + floor(gene_list$txEnd[gene]/resolution) #note floor here, because of the start bin
    loopscore_row_name=paste(as.character(startbin), '_', as.character(endbin), sep='')
    genescore = chia_pet_matrix[loopscore_row_name,3]
     if(is.na(genescore)==F) {

             chrom_size=chrom_sizes[which(as.character(chrom_sizes[,1])==chrom),2]
       #this location gives a start location - from 1 to chromosome size-gene_size

       random_tx_start = sample(seq(from =1 , to = chrom_size -gene_size), no_genes, replace = T )
       random_tx_end  = random_tx_start +gene_size
       random_start_locations = as.numeric(chrbins[which(as.character(chrbins[,1])==chrom),2]) + floor(random_tx_start/resolution) # floor here, because of the start bin
    random_end_locations = as.numeric(chrbins[which(as.character(chrbins[,1])==chrom),2]) + floor(random_tx_end/resolution) #note floor here, because of the start bin

    random_loopscore_rownames = paste(as.character(random_start_locations), '_', as.character(random_end_locations), sep='')
        ind <- rownames(chia_pet_matrix) %in% random_loopscore_rownames
    random_scores = chia_pet_matrix[ind,3]
    random_scores = c(random_scores, rep(0,no_genes-length(random_scores)))

    adjusted_genescore=genescore/mean(random_scores)
    loopscoreslist=rbind(loopscoreslist,c(gene_name,chrom,as.numeric(gene_list$txStart[gene]),as.numeric(gene_list$txEnd[gene]),adjusted_genescore))

     }
    if(is.na(genescore)==T) {
      adjusted_genescore=0
      loopscoreslist=rbind(loopscoreslist,c(gene_name,chrom,as.numeric(gene_list$txStart[gene]),as.numeric(gene_list$txEnd[gene]),adjusted_genescore))

    }
  }
write.csv(loopscoreslist, file = '/home/u1762230/chia_pet/K562/rep_2_SRR372748_hg19/R_analysis/rep2_hg19_adjusted_loopscores_all_genes_10000_randomlocs_res_2000_no_binning.csv')




