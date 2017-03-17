## Computational method to estimate the PCR duplication rate in high-throughput DNA sequencing experiments 

PCR amplification is an important step in the preparation of DNA sequencing libraries prior to high-throughput sequencing. Existing computational methods for analysis of read duplicates assume that all read duplicates arise due to PCR amplification. However, a high rate of read duplicates is observed in deep sequencing experiments or experiments such as RNA-seq. We present a computational method that exploits the heterozygosity in diploid genomes to estimate the PCR duplication rate accounting for read duplicates that are not due to PCR amplification. 

A paper describing this method has been published in BMC Bioinformatics, March 2017: http://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1471-9 


### INPUT files required 
1. coordinate sorted BAM file with aligned reads 
2. VCF file with variants called from BAM file using a variant calling tool such as GATK UnifiedGenotyper or samtools

###  Two step process to extract reads overlapping variant sites and then analyze read clusters to estimate PCR duplication rate 

1. make all 
2. ./extract\_duplicates --bam sample.bam --VCF variants.VCF > sample.hetreads 
3. python estimate\_PCRduprate.py sample.hetreads rna > sample.PCRdups 
4. grep FINAL\_PCR\_RATE sample.PCRdups > sample.PCRduprate.estimate

### Sample Data 

see DATA folder
