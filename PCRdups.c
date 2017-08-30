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
//#define SOFT_CLIP_EXTEND  // if enabled, soft-clip portion of read is used for finding variants 

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


int MINQ = 13; // minimum base quality
int MIN_MQ = 30; // minimum read mapping quality
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
int OUTPUT_VCF = 0;
int OUTPUT_FR = 0;  // only output first in pair reads (flag = 64) that overlap het vars
int TREAT_SE = 0; // treat PE bam file as single-end, read1 if value = 1, read2 is output if value = 2
int FILTER_DUPS = 0; 

// last modified 02/18/2016, works for both PE and SE reads 

//int get_chrom_name(struct alignedread* read,HASHTABLE* ht,REFLIST* reflist);

#include "parsebamread.c"

void print_options();
int parse_bamfile_sorted(char* bamfile,HASHTABLE* ht,CHROMVARS* chromvars,VARIANT* varlist,REFLIST* reflist);

void print_options()
{
	fprintf(stderr,"\n PROGRAM TO estimate read duplicate cluster counts and extract reads overlapping heterozygous variants from sorted BAM files \n\n");
	fprintf(stderr,"./extract_duplicates [options] --bam reads.sorted.bam --VCF variants.vcf   > output.reads \n\n");
	fprintf(stderr,"=============== PROGRAM OPTIONS ======================================== \n\n");
	fprintf(stderr,"--mbq  : minimum base quality to consider a base, default 13\n");
	fprintf(stderr,"--mmq : minimum read mapping quality to consider a read, default 30\n");
	fprintf(stderr,"--treatSE 0/1: consider paired-end reads as Single End (first in pair), default 0\n");
	fprintf(stderr,"--filterdups 0/1: discard reads marked as Duplicates (cigar flag of 2048), default 0\n");
	//fprintf(stderr,"--variants : variant file with genotypes for a single individual in genotype format\n");
	//fprintf(stderr,"--indels 0/1 : extract reads spanning INDELS, default is 0, variants need to specified in VCF format to use this option\n");
	//fprintf(stderr,"--outputVCF : 0 is default (output reads and filtered variants), 1: output VCF file of filtered variants, 2: do not output filtered variants\n");
	fprintf(stderr,"\n ============================================================================ \n\n");
	//fprintf(stderr,"--out : output file for haplotype informative fragments (hairs)\n\n");
}

double ncr(int n,int r)
{
        if (r ==0 || n ==r) return 0;
        double ll =0; int i=0;
        if (2*r > n) r = n-r;
        for (i=0;i<r;i++)  ll += log10((double)(n-i)/(i+1));
        return ll;
}

double binomial_test(int R, int A, double e,double* pv0,double* pv1)
{
        double e1 = log10(e); double e2 = log10(1.0-e); int r=0;
        double ll = ncr(R,A) + A*e1 + (R-A)*e2; double pvlog = ll; *pv1 = ll;
        for (r=A+1;r<R+1;r++)
        {
                pvlog += log10(R-r+1) - log10(r) + e1-e2;
                if (pvlog > *pv1) *pv1 = pvlog + log10(1.0+pow(10,*pv1-pvlog));
                else *pv1 += log10(1.0+pow(10,pvlog-*pv1));
        }
        pvlog = R*e2; *pv0 = pvlog;
        for (r=1;r< A+1;r++)
        {
                pvlog += log10(R-r+1) - log10(r) + e1-e2;
                if (pvlog > *pv0) *pv0 = pvlog + log10(1.0+pow(10,*pv0-pvlog));
                else *pv0 += log10(1.0+pow(10,pvlog-*pv0));
        }
        return 1;
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
        int MAX_CSIZE = 2048;
        int dup_clusters[MAX_CSIZE]; for (i=0;i<MAX_CSIZE;i++) dup_clusters[i] =0; i=0; // # of clusters  


	samfile_t *fp;
	if ((fp = samopen(bamfile, "rb", 0)) == 0) { fprintf(stderr, "Fail to open BAM file %s\n", bamfile); return -1; }
	bam1_t *b = bam_init1();
	int flag = 0; int reads1=0; int r1 = 0; int strand = 0;

	while (samread(fp, b) >= 0)
	{
		reads1 +=1;
		flag = 0; 
		fetch_func(b, fp,read);  
		if (FILTER_DUPS ==1 && read->flag > 256) flag = 1; 
		if ((read->flag & (BAM_FUNMAP)) || read->mquality < MIN_MQ || (read->flag & 256)) flag = 1; 
		// also ignore non-primary alignments, BWA MEM and output cigar for read
		if ((read->flag & 1) && (read->IS ==0 || (read->flag & 8) || read->tid != read->mtid) )  flag = 1; 
		// we should filter single-end reads in a file with both single-end and paired-end reads ?? unlikely 
		if (FILTER_SE ==1 && (read->flag & 1) ==0) flag = 1; 

		if (flag ==0 && TREAT_SE ==1) // treat PE read as single-end read, read1 only 
		{
			r1 = read->flag & 64; if (r1 == 0 || read->IS ==0) flag =1;  
			read->IS = 0;  FILTER_SE = 0;  strand = read->flag & 16; read->flag = strand; read->mateposition = 0; read->mtid = -1;
		}
		reads+=1; if (reads%2000000 ==0) fprintf(stderr,"processed %d reads, useful fragments %d, filtered-reads %d current read->tid %d\n",reads,useful_reads,PCRdups_stats[0],read->tid);
		if (flag ==1) 
		{
			PCRdups_stats[0] +=1; 	free_readmemory(read); continue;
		}
		if ((read->flag & 2048)) PCRdups_stats[1] +=1; PCRdups_stats[2] +=1; 


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
			if (fragment.variants > 0 && (OUTPUT_FR == 0 || ((read->flag & 64) == 64) ) ) 
			{
				if ((read->flag & 64) == 64) r1 = 1; else r1 = 2;
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


int main (int argc, char *argv[])
{
	char samfile[2048]; char bamfile[2048]; char variantfile[2048]; char fastafile[2048];
	strcpy(samfile,"None"); strcpy(bamfile,"None"); strcpy(variantfile,"None"); strcpy(fastafile,"None");
	GROUPNAME = NULL;
	int readsorted = 0;
	char* sampleid = (char*)malloc(2048); sampleid[0] = '-'; sampleid[1] = '\0';
	int samplecol=10; // default if there is a single sample in the VCF file
	int i=0,j=0,variants=0,hetvariants=0;
	int bamfiles =0, vfile = 0;

	fprintf(stderr,"DEBUG: # of args %d \n",argc); for (i=0;i<argc;i++) fprintf(stderr,"argument %d %s \n",i,argv[i]); 

	logfile = NULL;
	for (i=1;i<argc;i+=2)
	{
		if (strcmp(argv[i],"--bam") ==0 || strcmp(argv[i],"--bamfile") ==0 )  { strcpy(bamfile,argv[i+1]); bamfiles++; } 
		else if (strcmp(argv[i],"--reffile") ==0 || strcmp(argv[i],"--ref") ==0)        strcpy(fastafile,argv[i+1]);
		else if (strcmp(argv[i],"--VCF") ==0 || strcmp(argv[i],"--vcf") ==0 )    {     strcpy(variantfile,argv[i+1]); VCFformat =1; vfile = 1;  }
		else if (strcmp(argv[i],"--sorted") ==0)       readsorted = atoi(argv[i+1]);
		else if (strcmp(argv[i],"--mbq") ==0 )       MINQ = atoi(argv[i+1]);
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
		else if (strcmp(argv[i],"--outputVCF")==0) OUTPUT_VCF = atoi(argv[i+1]);  
		else if (strcmp(argv[i],"--outputFR")==0) OUTPUT_FR = atoi(argv[i+1]);  
		else if (strcmp(argv[i],"--treatSE")==0) TREAT_SE = atoi(argv[i+1]);  
		else if (strcmp(argv[i],"--filterdups")==0) { FILTER_DUPS = atoi(argv[i+1]); fprintf(stderr,"read duplicates marked in BAM file will be ignored \n"); } 
		else if (strcmp(argv[i],"--variants") ==0)      {  strcpy(variantfile,argv[i+1]); vfile = 1; } // old variant format 
	}
	if (bamfiles > 0 && vfile > 0)
	{
		fprintf(stderr,"\n extracting reads covering SNPs from bamfile %s minQV %d minMQ %d maxIS %d \n\n",bamfile,MINQ,MIN_MQ,MAX_IS);
	}
	else
	{
		int flag1=0;
		strcpy(bamfile,"sample.bam"); strcpy(variantfile,"sample.VCF"); 
		FILE* file1 = fopen(bamfile,"r"); if (file1 != NULL) { flag1 +=1; fclose(file1); bamfiles +=1; } 
		file1 = fopen(variantfile,"r"); if (file1 != NULL) { flag1 +=1; fclose(file1); VCFformat =1; }
		fprintf(stderr,"flag %d \n",flag1);
		if (flag1 < 2) 
		{
			fprintf(stderr,"program options are missing, please specify at least one BAM file and a VCF file \n bamfiles %d  variantfile %s \n",bamfiles,variantfile); 
			print_options(); return -1;
		}
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

	VARIANTS = variants;  
	CHROMVARS* chromvars  = (CHROMVARS*)malloc(sizeof(CHROMVARS)*chromosomes);
	build_intervalmap(chromvars,chromosomes,varlist,VARIANTS);

	// read reference fasta file for INDELS//
	REFLIST* reflist = (REFLIST*)malloc(sizeof(REFLIST)); 
	reflist->ns = 0; reflist->names = NULL; reflist->lengths = NULL; reflist->sequences = NULL; reflist->current = -1;
	if (strcmp(fastafile,"None") != 0)
	{
		if (read_fastaheader(fastafile,reflist) > 0) 
		{
			reflist->sequences = (unsigned char**)malloc(sizeof(unsigned char*)*reflist->ns);
			for (i=0;i<reflist->ns;i++)
			{
				reflist->sequences[i] = (unsigned char*)malloc(reflist->lengths[i]+1);
				if (i < 5) fprintf(stderr,"contig %s length %d\n",reflist->names[i],reflist->lengths[i]);
			}
			read_fasta(fastafile,reflist);
		}
	}
	if (readsorted ==0 && bamfiles > 0)
	{
		parse_bamfile_sorted(bamfile,&ht,chromvars,varlist,reflist);
	}
	if (logfile != NULL) fclose(logfile); 

	// print variant-list for each chromosome as it is finished processing 
        int xor = pow(2,16)-1; int c0=0,c1=0; float ratio = 0.4; int flag = 0;
	double pv0=0,pv1=0,pv=0;

	if (VCFformat ==1) 
	{
		for (i=0;i<variants;i++)
		{
			if (varlist[i].type !=0) continue;
			if (varlist[i].genotype[0] == varlist[i].genotype[2]) continue;
			c0 = (varlist[i].A1>>16) + (varlist[i].A1 & xor); c1 = (varlist[i].A2>>16) + (varlist[i].A2 & xor); 
			if (OUTPUT_VCF ==3 && (ratio < 0.3 || ratio > 0.7) && varlist[i].depth >= 30) fprintf(stdout,"#variant-all %d %s %s %d %d %s %s ref:%d:%d alt:%d:%d %d:%d:%0.3f depth %d\n",i,varlist[i].genotype,varlist[i].chrom,varlist[i].position,varlist[i].type,varlist[i].RA,varlist[i].AA,varlist[i].A1>>16,varlist[i].A1 & xor,varlist[i].A2>>16,varlist[i].A2 & xor,c0,c1,ratio,varlist[i].depth);
			if (varlist[i].depth > 0) ratio = (float)c1/varlist[i].depth;  else ratio = 0.0; 

			if (OUTPUT_VCF ==0) 
			{
				binomial_test(c0+c1,c0,0.5,&pv0,&pv1); pv = pv0; if (pv1 < pv) pv = pv1; 
				flag = 0; if ( (ratio < 0.1 || ratio > 0.9) && varlist[i].depth >= 30) flag = 1; if (varlist[i].depth < 1) flag = 1;   
				if (flag ==1) 
				{
					if (varlist[i].depth > 0) fprintf(stdout,"#filtered %d %s %s %d %d %s %s ref:%d:%d alt:%d:%d %d:%d:%0.3f pv:%f\n",i,varlist[i].genotype,varlist[i].chrom,varlist[i].position,varlist[i].type,varlist[i].RA,varlist[i].AA,varlist[i].A1>>16,varlist[i].A1 & xor,varlist[i].A2>>16,varlist[i].A2 & xor,c0,c1,ratio,pv);
				}
				else fprintf(stdout,"#variant %d %s %s %d %d %s %s ref:%d:%d alt:%d:%d %d:%d:%0.3f\n",i,varlist[i].genotype,varlist[i].chrom,varlist[i].position,varlist[i].type,varlist[i].RA,varlist[i].AA,varlist[i].A1>>16,varlist[i].A1 & xor,varlist[i].A2>>16,varlist[i].A2 & xor,c0,c1,ratio);
			}

			flag = 0; if (c1 < 3 || varlist[i].depth < 8) flag = 1; if (ratio < 0.3 || ratio > 0.7) flag = 1; 
			if (flag ==0) 
			{
				if (OUTPUT_VCF == 4) fprintf(stdout,"#variant %d %s %s %d %d %s %s ref:%d:%d alt:%d:%d %d:%d:%0.3f pv:%f\n",i,varlist[i].genotype,varlist[i].chrom,varlist[i].position,varlist[i].type,varlist[i].RA,varlist[i].AA,varlist[i].A1>>16,varlist[i].A1 & xor,varlist[i].A2>>16,varlist[i].A2 & xor,c0,c1,ratio,pv);
			}
			if (flag ==1) continue; 
			if (OUTPUT_VCF ==1) fprintf(stdout,"#VCFchr%s\t%d\t%d\t%s\t%s\t.\tPASS\tREF=%d,ALT=%d\tGT:DP\t0/1:%d\n",varlist[i].chrom,varlist[i].position,i,varlist[i].RA,varlist[i].AA,c0,c1,varlist[i].depth);
		}
	}
	return 0;
}


