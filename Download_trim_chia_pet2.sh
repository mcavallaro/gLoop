#!/bin/bash

#written by Matt Jones, contact m.jones.18@warwick.ac.uk

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR372/SRR372747/SRR372747_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR372/SRR372747/SRR372747_2.fastq.gz
trim_galore --paired SRR372747_1.fastq.gz SRR372747_2.fastq.gz

ChIA-PET2 -s 1 -g /home/u1762230/reference/hs-bwa/hg19bwaidx -b /home/u1762230/bedtoolsgenomefile/hg19.chrom.sizes -f SRR372747_1_val_1.fq.gz -r SRR372747_2_val_2.fq.gz -t 20 --short 1 -Q 20 -o rep_1_SRR372747_hg19

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR372/SRR372748/SRR372748_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR372/SRR372748/SRR372748_2.fastq.gz
trim_galore --paired SRR372748_1.fastq.gz SRR372748_2.fastq.gz
ChIA-PET2 -s 1 -g /home/u1762230/reference/hs-bwa/hg19bwaidx -b /home/u1762230/bedtoolsgenomefile/hg19.chrom.sizes -f SRR372748_1_val_1.fq.gz -r SRR372748_2_val_2.fq.gz -t 20 --short 1 -Q 20 -o rep_2_SRR372748_hg19    





