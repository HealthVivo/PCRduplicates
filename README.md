# Computational method to estimate the PCR duplication rate in high-throughput DNA sequencing experiments 

PCR amplification is an important step in the preparation of DNA sequencing libraries prior to high-throughput sequencing. Existing computational methods for analysis of read duplicates assume that all read duplicates arise due to PCR amplification. However, a high rate of read duplicates is observed in deep sequencing experiments or experiments such as RNA-seq. We present a computational method that exploits the heterozygosity in diploid genomes to estimate the PCR duplication rate accounting for read duplicates that are not due to PCR amplification. 


# INPUT files required 
1. coordinate sorted BAM file with aligned reads 
2. VCF file with variants called from BAM file using a variant calling tool such as GATK UnifiedGenotyper or samtools

# TWO STEP PROCESS 

1. make all 

2. ./extract\_duplicates --bam sample.bam --VCF variants.VCF > sample.hetreads 
3. python estimate\_PCRduprate.py sample.hetreads rna > sample.PCRdups 
4. grep FINAL\_PCR\_RATE sample.PCRdups > sample.PCRduprate.estimate 
