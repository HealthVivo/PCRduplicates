## Computational method to estimate the PCR duplication rate in high-throughput DNA sequencing experiments 

PCR amplification is an important step in the preparation of DNA sequencing libraries prior to high-throughput sequencing. Existing computational methods for analysis of read duplicates assume that all read duplicates arise due to PCR amplification. However, a high rate of read duplicates is observed in deep sequencing experiments or experiments such as RNA-seq. We present a computational method that exploits the heterozygosity in diploid genomes to estimate the PCR duplication rate accounting for read duplicates that are not due to PCR amplification. 

A paper describing this method has been published in BMC Bioinformatics, March 2017: http://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1471-9 

# INPUT files for running program 
1. coordinate sorted BAM file with aligned reads 
2. VCF file with heterozygous variants called from BAM file using a variant calling tool such as GATK UnifiedGenotyper or samtools

# compile the CODE 

run 'make all'

#  Two step process to extract reads overlapping variant sites and then analyze read clusters to estimate PCR duplication rate 

1. ./extract\_duplicates --bam sample.bam --VCF variants.VCF > sample.hetreads 
2. python estimate\_PCRduprate.py -i sample.hetreads -f exome > sample.PCRdups 


# obtain FINAL estimate of the PCR duplication rate 

grep FINAL\_PCR\_RATE sample.PCRdups > sample.PCRduprate.estimate

# Sample Data 

see DATA folder


# FAQ

1. The program uses samtools (v 0.1.18) to parse BAM files. The source code for samtools is included in the github repository (directory parsebam/samtools-0.1.18). 
2. The program has been tested on exome-seq, targeted DNA seq and RNA-seq datasets. For RNA-seq, an independent set of heterozygous variants is needed. 
