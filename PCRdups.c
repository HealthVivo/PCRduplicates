
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>

#include "hashtable.h"
#include "readfasta.h"
#include "bamread.h"
#include "sam.h"
#include "readvariant.h"
#include "hapfragments.h"

#define OUTPUT_THIRD_ALLELE 

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


int MINQ = 13; // minimum base quality
int MIN_MQ = 20; // minimum read mapping quality
int MAX_IS =  1000; // maximum insert size
int MIN_IS =  0; // maximum insert size
int PEONLY = 0; // if this is set to 1, reads for which only one end is mapped are not considered for hairs 
int BSIZE = 500;
int IFLAG = 0;
int MAXFRAG = 500000;
int VARIANTS = 0;
int VCFformat = 0;
int PARSEINDELS =0;
int SINGLEREADS =1;
int FOSMIDS = 0;
//int QVoffset = 33; declared in samread.h
FILE* logfile;
int PFLAG = 1;
int PRINT_FRAGMENTS = 0;
int POOL_SIZE = 2;// 
int TRI_ALLELIC = 0;
FILE* fragment_file;
char* GROUPNAME;
int ADD_CHR = 0; // add 'chr' to chromosome names in VCF file, 1-> chr1
int FILTER_SE = 1; 

// last modified 02/18/2016, works for both PE and SE reads 

//int get_chrom_name(struct alignedread* read,HASHTABLE* ht,REFLIST* reflist);

#include "parsebamread.c"

void print_options();
int parse_bamfile_sorted(char* bamfile,HASHTABLE* ht,CHROMVARS* chromvars,VARIANT* varlist,REFLIST* reflist);

void print_options()
{
	fprintf(stderr,"\n PROGRAM TO analyze PCR duplicates in coordinate sorted BAM files \n\n");
	fprintf(stderr,"./output_het_reads [options] --bam reads.sorted.bam --variants variants.file   > output.fragments \n\n");
	fprintf(stderr,"=============== PROGRAM OPTIONS ======================================== \n\n");
	fprintf(stderr,"--qvoffset : quality value offset, 33/64 depending on how quality values were encoded, default is 33 \n");
	fprintf(stderr,"--mbq  : minimum base quality to consider a base, default 13\n");
	fprintf(stderr,"--mmq : minimum read mapping quality to consider a read, default 20\n");
	fprintf(stderr,"--variants : variant file with genotypes for a single individual in genotype format\n");
	fprintf(stderr,"--PEonly 0/1 : do not use single end reads, default is 1 (use all reads)\n");
	fprintf(stderr,"\n ============================================================================ \n\n");
	//fprintf(stderr,"--out : output file for haplotype informative fragments (hairs)\n\n");
}


// algorithm : 1. read through bam file, identify reads that overlap heterozygous SNPs and output their start position, IS, allele to text file... 
// sort the list of reads by "start position + IS", identify duplicate clusters and analyze...
// note that current bam file is only sorted by start position, not by IS (:
// we should consider a paired-end read as a single entity for correct analysis 

typedef struct
{
        int tid; int mtid; int pos; int mpos; int IS; int flag;

} PEread;  // data structure to store basic info about PE_read 

int compare_PEread(const void *a,const void *b)
{
        PEread* p1 = (PEread*)a; PEread* p2 = (PEread*)b;
	if (p1->IS == p2 ->IS) return 1; 
        if (p1->pos == p2->pos) return (p1->IS - p2->IS);
        else return (p1->pos -p2->pos);
}


// PE-reads with mate mapped to different chromosome, PE-reads with low mapping quality, not handled well right now... 
// low mapping quality cannot be handled correctly since the duplicate reads may be mapped to different locations even if they are PCR dups 

int parse_bamfile_sorted(char* bamfile,HASHTABLE* ht,CHROMVARS* chromvars,VARIANT* varlist,REFLIST* reflist)
{
	fprintf(stderr,"reading sorted bamfile to identify PCR duplicates %s \n",bamfile);
	int reads=0;
	struct alignedread* read = (struct alignedread*)malloc(sizeof(struct alignedread));
	
	int i=0,j=0; int sl=0; int chrom=0;
	int v1,v2; int absIS;
	int prevchrom=-1; int prevtid = -1; int delta_start =0;

	FRAGMENT* flist = (FRAGMENT*)malloc(sizeof(FRAGMENT)*MAXFRAG); int fragments =0; int prevfragments =0;
	FRAGMENT fragment; fragment.variants =0; fragment.alist = (allele*)malloc(sizeof(allele)*1000);
	int useful_reads = 0;

	int MAX_BUFF_SIZE = 4096; 
	PEread* readlist = (PEread*)malloc(sizeof(PEread)*MAX_BUFF_SIZE); int pereads = 0; int csize = 0; int prevpos = -1; 
	int PCRdups_stats[3] = {0,0,0}; 
	int total_reads= 0; int unique_reads=0;
        int MAX_CSIZE = 1024;
        int dup_clusters[MAX_CSIZE]; for (i=0;i<MAX_CSIZE;i++) dup_clusters[i] =0; i=0; // # of clusters  


	samfile_t *fp;
	if ((fp = samopen(bamfile, "rb", 0)) == 0) { fprintf(stderr, "Fail to open BAM file %s\n", bamfile); return -1; }
	bam1_t *b = bam_init1();
	int flag = 0; int reads1=0;

	while (samread(fp, b) >= 0)
	{
		reads1 +=1;
		flag = 0; 
		fetch_func(b, fp,read); 
		if ((read->flag & (BAM_FUNMAP)) || read->mquality < MIN_MQ || (read->flag & 256)) flag = 1; 
		// also ignore non-primary alignments, BWA MEM and output cigar for read
		if ((read->flag & 1) && (read->IS ==0 || (read->flag & 8) || read->tid != read->mtid) )  flag = 1; 
		// we should filter single-end reads in a file with both single-end and paired-end reads ?? unlikely 
		if (FILTER_SE ==1 && (read->flag & 1) ==0) flag = 1; 
		reads+=1; if (reads%2000000 ==0) fprintf(stderr,"processed %d reads, useful fragments %d, filtered-reads %d current read->tid %d\n",reads,useful_reads,PCRdups_stats[0],read->tid);

		if (flag ==1) 
		{
			PCRdups_stats[0] +=1; 	free_readmemory(read); continue;
		}
		if ((read->flag & 1024)) PCRdups_stats[1] +=1; PCRdups_stats[2] +=1; 


		// find the chromosome in reflist that matches read->chrom if the previous chromosome is different from current chromosome
		if (read->tid != prevtid)
		{
			chrom = getindex(ht,read->chrom); // doing this for every read, can replace this by string comparison ..april 4 2012
			i = read->tid;
			if (reflist->ns > 0)
			{
				reflist->current = i;
				if (i >= reflist->ns || i < 0 || strcmp(reflist->names[i],read->chrom) !=0)
				{
					reflist->current = -1;
					for (i=0;i<reflist->ns;i++)
					{
						if (strcmp(reflist->names[i],read->chrom) ==0) { reflist->current = i; break; }
					}
				}
			}
		}
		else chrom = prevchrom;

		absIS = (read->IS < 0) ? -1*read->IS: read->IS; 
		fragment.variants =0; v1 =0; 
		if (chrom >=0) 
		{
			fragment.id = read->readid;
			v1 = extract_variants_read(read,ht,chromvars,varlist,0,&fragment,chrom,reflist);
			if (fragment.variants > 0) 
			{
				//fprintf(stderr,"reads %d flag %d\n",reads1,flag);
				// consider overlapping paired-end reads, simple solution, if variant in read starts after mateposition, ignore it to avoid double counting
				//if (read->IS != 0 && (( read->flag & 8) == 0)) 
				{
					// print mate position if insert size is negative  
					delta_start=0;
					if (read->IS < 0) // second read in pair and soft clip at outer edge 
					{
						if ((read->cigarlist[read->cigs-1]&0xf) == BAM_CSOFT_CLIP) delta_start = read->cigarlist[read->cigs-1]>>4; 
						if (read->flag & 1) read->IS += delta_start; // modify only for PE reads 
						fprintf(stdout,"%s %s %d %d %d %d,",read->readid,read->chrom,read->mateposition,read->position,read->IS,read->flag);
					}
					else // first read in pair and soft clip at start
					{
						if ((read->cigarlist[0]&0xf) == BAM_CSOFT_CLIP) delta_start = read->cigarlist[0]>>4; 
						if (read->flag & 1) { read->IS += delta_start; read->position -= delta_start; } // modify the IS,read-pos only for PE reads 
						fprintf(stdout,"%s %s %d %d %d %d,",read->readid,read->chrom,read->position,read->mateposition,read->IS,read->flag);
					}
					for (j=0;j<read->cigs;++j) fprintf(stdout,"%d%c:",read->cigarlist[j]>>4,INT_CIGAROP[read->cigarlist[j]&0xf]); 
					for (i=0;i<fragment.variants;i++) fprintf(stdout,"\t%d:%c:%c",fragment.alist[i].varid,fragment.alist[i].allele,fragment.alist[i].qv);
					fprintf(stdout,"\n");
				}
				useful_reads++;
				//fprintf(stdout," useful frag \n");
			}
		}
	
		// CODE for PCR dup-clustering, we should only consider reads with IS > 0 to avoid double-counting clusters 
		if ((read->tid != prevtid || read->position != prevpos) && pereads > 0)
                {
                        qsort(readlist,pereads,sizeof(PEread),compare_PEread); csize = 1; 
                        for (i=0;i<pereads-1;i++)
                        {
				//if (readlist[i].IS < 0) continue; 
                                if (readlist[i].pos == readlist[i+1].pos && readlist[i].IS == readlist[i+1].IS)
				{
					//fprintf(stderr,"match %d %d %d i %d %d \n",readlist[i].pos,readlist[i].IS,readlist[i].mpos,i,pereads);
					 csize +=1; 
				}
                                else
                                {
                                        if (csize < MAX_CSIZE) dup_clusters[csize] +=1; 
                                        csize = 1; 
                                }
                        }
                        if (csize < MAX_CSIZE) dup_clusters[csize] +=1; 
                        pereads = 0; 
                }
                if (read->IS >0 || ((read->flag & 1) == 0 && FILTER_SE ==0))  // only reads with IS > 0 OR single-end reads 
		{ 
			readlist[pereads].pos = read->position; readlist[pereads++].IS = read->IS;  
			if (pereads >= MAX_BUFF_SIZE) { MAX_BUFF_SIZE *=2;  readlist = realloc(readlist,sizeof(PEread)*MAX_BUFF_SIZE);} 
		} 
		prevpos = read->position; 
		//////////////////////////////////////////////////////////////////////////////

		prevchrom = chrom; prevtid = read->tid;
		free_readmemory(read);
	}
	bam_destroy1(b); free(readlist);

	fprintf(stderr,"PCR duplicates marked %d total-reads %d frac %0.4f discarded %d\n",PCRdups_stats[1],PCRdups_stats[2],(float)PCRdups_stats[1]/PCRdups_stats[2],PCRdups_stats[0]);
	unique_reads = 0; total_reads =0; 
	fprintf(stdout,"#clusters "); fprintf(stderr,"#clusters ");
	for (i=1;i<MAX_CSIZE;i++)
        {	
		unique_reads += dup_clusters[i]; total_reads += dup_clusters[i]*i;
                if (dup_clusters[i] > 0) fprintf(stderr,"%d:%d ",i,dup_clusters[i]); 
                if (dup_clusters[i] > 0) fprintf(stdout,"%d:%d ",i,dup_clusters[i]); else break;
        }
	fprintf(stdout,"\n"); fprintf(stderr,"\n");
	fprintf(stderr,"total reads (PE=1) %d unique-reads %d duplicates:%d, duplication rate %0.5f \n",total_reads,unique_reads,total_reads-unique_reads,1.0-(float)unique_reads/total_reads);
}


int main (int argc, char** argv)
{
	char samfile[1024]; char bamfile[1024]; char variantfile[1024]; char fastafile[1024];
	strcpy(samfile,"None"); strcpy(bamfile,"None"); strcpy(variantfile,"None"); strcpy(fastafile,"None");
	GROUPNAME = NULL;
	int readsorted = 0;
	char* sampleid = (char*)malloc(1024); sampleid[0] = '-'; sampleid[1] = '\0';
	int samplecol=10; // default if there is a single sample in the VCF file
	int i=0,j=0,variants=0,hetvariants=0;
	char** bamfilelist = NULL; int bamfiles =0; 

	logfile = NULL;
	for (i=1;i<argc;i+=2)
	{
		if (strcmp(argv[i],"--bam") ==0 || strcmp(argv[i],"--bamfile") ==0)        bamfiles++; 
		else if (strcmp(argv[i],"--variants") ==0)        strcpy(variantfile,argv[i+1]);
		else if (strcmp(argv[i],"--reffile") ==0 || strcmp(argv[i],"--ref") ==0)        strcpy(fastafile,argv[i+1]);
		else if (strcmp(argv[i],"--VCF") ==0 || strcmp(argv[i],"--vcf") ==0)    {     strcpy(variantfile,argv[i+1]); VCFformat =1; }
		else if (strcmp(argv[i],"--sorted") ==0)       readsorted = atoi(argv[i+1]);
		else if (strcmp(argv[i],"--mbq") ==0)       MINQ = atoi(argv[i+1]);
		else if (strcmp(argv[i],"--mmq") ==0 || strcmp(argv[i],"--minmq") ==0 )       MIN_MQ = atoi(argv[i+1]);
		else if (strcmp(argv[i],"--maxIS") ==0)       MAX_IS = atoi(argv[i+1]);
		else if (strcmp(argv[i],"--minIS") ==0)       MIN_IS = atoi(argv[i+1]);
		else if (strcmp(argv[i],"--PEonly") ==0)       PEONLY = 1;  // discard single end mapped reads 
		else if (strcmp(argv[i],"--filterSE") ==0)       FILTER_SE = atoi(argv[i+1]);  // discard single end mapped reads 
		else if (strcmp(argv[i],"--indels") ==0)       PARSEINDELS = atoi(argv[i+1]);  // allow indels in hairs
		else if (strcmp(argv[i],"--pflag") ==0)      IFLAG  = atoi(argv[i+1]);  // allow indels in hairs
		else if (strcmp(argv[i],"--qvoffset") ==0)       QVoffset = atoi(argv[i+1]);
		else if (strcmp(argv[i],"--logfile")==0 || strcmp(argv[i],"--out") ==0) logfile = fopen(argv[i+1],"w");  
		else if (strcmp(argv[i],"--singlereads")==0) SINGLEREADS = atoi(argv[i+1]);  
		else if (strcmp(argv[i],"--maxfragments")==0) MAXFRAG = atoi(argv[i+1]);  
	}
	if (bamfiles > 0 && strcmp(variantfile,"None") !=0)
	{
		bamfilelist = (char**)malloc(sizeof(char*)*bamfiles); 
		for (i=0;i<bamfiles;i++) bamfilelist[i] = (char*)malloc(1024);
		bamfiles=0;
		for (i=1;i<argc;i+=2)
		{
			if (strcmp(argv[i],"--bam") ==0 || strcmp(argv[i],"--bamfile") ==0)     strcpy(bamfilelist[bamfiles++],argv[i+1]);
		}
		fprintf(stderr,"\n extracting reads covering SNPs from bamfile %s minQV %d minMQ %d maxIS %d \n\n",bamfilelist[0],MINQ,MIN_MQ,MAX_IS);
	}
	else
	{
		print_options(); return -1;
	}

	HASHTABLE ht; ht.htsize = 7919;  init_hashtable(&ht);
	VARIANT* varlist;
	int chromosomes=0;

	if (VCFformat ==1)// added 02/22/2016
        {
                variants = count_variants(variantfile,sampleid,&samplecol);
                if (variants < 0) return -1;
                varlist = (VARIANT*)malloc(sizeof(VARIANT)*variants);
                chromosomes = read_variantfile(variantfile,varlist,&ht,&hetvariants,samplecol);
        }
        else
	{	// variant format is old one !! not VCF 
		variants = count_variants_oldformat(variantfile);
		if (variants < 0) return -1; else varlist = (VARIANT*)malloc(sizeof(VARIANT)*variants);
		chromosomes = read_variantfile_oldformat(variantfile,varlist,&ht,variants);
	}

	for (i=0;i<variants;i++)
	{
		varlist[i].GLL = calloc(POOL_SIZE,sizeof(double));
		for (j=0;j<POOL_SIZE;j++) varlist[i].GLL[j] = 0.0;
	}

	// variants is set to hetvariants only, but this is not correct since 
	VARIANTS = variants;  
	// there are two options, we include all variants in the chromvars datastructure but only use heterozygous variants for outputting HAIRS 
	// variant-id should correspond to line-number in VCF file since that will be used for printing out variants in Hapcut 

	//	fprintf(stderr,"read %d variants from file %s chromosomes %d\n",snps,argv[1],chromosomes);
	CHROMVARS* chromvars  = (CHROMVARS*)malloc(sizeof(CHROMVARS)*chromosomes);
	build_intervalmap(chromvars,chromosomes,varlist,VARIANTS);

	// read reference fasta file for INDELS//
	REFLIST* reflist = (REFLIST*)malloc(sizeof(REFLIST)); 
	reflist->ns = 0; reflist->names = NULL; reflist->lengths = NULL; reflist->sequences = NULL; reflist->current = -1;
	if (strcmp(fastafile,"None") != 0)
	{
		if (read_fastaheader(fastafile,reflist) > 0) 
		{
			reflist->sequences = (char**)malloc(sizeof(char*)*reflist->ns);
			for (i=0;i<reflist->ns;i++)
			{
				reflist->sequences[i] = (char*)malloc(reflist->lengths[i]+1);
				if (i < 5) fprintf(stderr,"contig %s length %d\n",reflist->names[i],reflist->lengths[i]);
			}
			read_fasta(fastafile,reflist);
		}
	}
	//return 1;
	if (readsorted ==0 && bamfiles > 0)
	{
		for (i=0;i<bamfiles;i++) 
		{
			parse_bamfile_sorted(bamfilelist[i],&ht,chromvars,varlist,reflist);
		}
	}
	if (logfile != NULL) fclose(logfile); 

	// print variant-list for each chromosome as it is finished processing 
        int xor = pow(2,16)-1; int c0=0,c1=0; float ratio = 0.4; int flag = 0;

	if (VCFformat ==1) 
	{
		for (i=0;i<variants;i++)
		{
			if (varlist[i].type !=0) continue;
			if (varlist[i].genotype[0] == varlist[i].genotype[2]) continue;
			c0 = (varlist[i].A1>>16) + (varlist[i].A1 & xor); c1 = (varlist[i].A2>>16) + (varlist[i].A2 & xor); 
			if (varlist[i].depth < 8 || c1 < 4) continue; 
			ratio = (float)c1/varlist[i].depth;  flag = 0;
			if (varlist[i].depth < 20 && ( ratio < 0.3 || ratio > 0.7)) flag = 1; 
			if (ratio < 0.3 || ratio > 0.7) flag = 1; 
			if (flag ==1) continue; 
			fprintf(stdout,"#variant %d %s %s %d %d %s %s ref:%d:%d alt:%d:%d %d:%d:%0.3f\n",i,varlist[i].genotype,varlist[i].chrom,varlist[i].position,varlist[i].type,varlist[i].RA,varlist[i].AA,varlist[i].A1>>16,varlist[i].A1 & xor,varlist[i].A2>>16,varlist[i].A2 & xor,c0,c1,ratio);
			//fprintf(stdout,"#VCF %s\t%d\t%d\t%s\t%s\t.\tPASS\tREF=%d,ALT=%d\tGT:DP\t0/1:%d\n",varlist[i].chrom,varlist[i].position,i,varlist[i].RA,varlist[i].AA,c0,c1,varlist[i].depth);
		}
	}
	/*
		python code can read entire file and get list of het variants, store it in a hashtable... and then filter reads on the fly
        */


	return 0;
}


