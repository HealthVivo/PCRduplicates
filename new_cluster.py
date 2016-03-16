#! /usr/bin/env python
# AUTHOR VIKAS BANSAL last edited dec 2014
import sys, os, glob, string, subprocess,time, math, random
from scipy import optimize# import minimize_scalar
import extra_functions

## latest code, last updated feb 17 2016 

#python code/cluster_PCRdups.py data/2_B_12.reads.new > clusters

## calculate false PCR duplication rate using output from PCRdups.c 
## also determine number of duplicates per fragment (poisson or other distribution), depends on fragment length, GC content... 
## probability that a pair of false PCR duplicates is real also depends on cluster size
## we are only looking at reads that cover a heterozygous SNP... limited sample but almost random since SNPs should be independnet of PCR duplication
## correct estimation of duplicates should also improve fit of count distribution (possion | exp | binomial) 
## examine relation between cluster size, GC content and insert size via regression
## we can stratify estimate of false duplicates by these parameters if we have enough data...

QV_THRESH = 0.01;
VERBOSE =1;
WEIGHTED =1;
Coeffs = [];  CSIZE = 2; 
OPTIMIZE = 1;
full_counts = []; ## cluster counts, total (0) 0000 matching (1)
SER = 0.01;
AB_counts=[0,0,0];

MAX_CLUSTER_SIZE = 100;

def fact(n):
	p = 1;
	for i in xrange(2,n+1): p *= i; 
	return p; 

# http://code.activestate.com/recipes/218332-generator-for-integer-partitions/
def partitions(n):
	if n ==0: yield []; return 
	for p in partitions(n-1): 
		yield [1] + p  
		if p and (len(p) < 2 or p[1]  > p[0]): 
			yield [p[0] + 1] + p[1:]

## pass Lambdas, counts_IND to this function 
def sorted_partitions(n,Lambda,counts_IND,c0,c1,c2,pcr_prev):
	## ignore first and last element of  list since they need to be estimated... we can calculate the sum of remaining elements... 
	it = partitions(n); list_part = [];
	for l in it: list_part.append(l);

	print >>sys.stderr,"----calling partiton function with n=",n,'C(....)=',c0,'C(00..00)=',c1,'C(00..001)=',c2,c1/c0,c2/c0;
	sum0 = 0.0; sum1 = 0.0; PDF = []; 
	sum_sm = 0.0; # 0000001 expression
	sum_delta = 0.0; # sum0-sum1; 
	for partition in list_part: 
		partition.sort(); l = len(partition); coeff3 = math.pow(2,l-1);
		if l == n or l == 1: PDF.append(0.0); continue; print >>sys.stderr, partition; continue; 
		prod = counts_IND[l-1]; 
		coeff2 = 1.0; N=0; coeff1 = fact(l); c = 1; 
		singles=0;
		for i in xrange(len(partition)): 
			if partition[i] == 1: singles +=1; 

		for i in xrange(1,len(partition)): 
			if partition[i] != partition[i-1]: coeff1 /= fact(c); c = 1; 
			else: c +=1;
		coeff1 /= fact(c);  

		for el in partition: 
			#for j in xrange(el,1,-1): prod *= Lambda[j-1]; print 'l_' + `j-1`,
			prod *= Lambda[el-1]; 
			N += el-1; coeff2 /= fact(el-1);
		coeff2 *= fact(N); 
		prod *= coeff1*coeff2; 
	
		## sum1 is wrong !!! or sum2 is wrong...  	
		sum1 += prod; sum_delta += prod*(coeff3-1);  sum_sm += prod*singles; sum0 += prod*coeff3;
		PDF.append(prod*coeff3);
		#print >>sys.stderr, '|',partition,coeff1*coeff2,coeff3,'|',prod

	flag =0;
	if sum1 > c1: print >>sys.stderr,"flag sum1 exceeds c1 count from data",c1-sum1; flag =1;
	# if ind_(n-1) == 0 -> ind_n should be zero, if pcr_(n-1) = 0, then pcr_n = 0
	print >>sys.stderr, 'n','ind_est_1',c1-sum1,c0-sum0,math.pow(2,n-1),'three-sums',sum0,sum1,sum_sm;
	if pcr_prev < 0.01: 
		pcr_n  = 0.0; ind_n = c0-sum0;
		## ind_n can be easily estimated from c0 only...
		ind_n /= math.pow(2,n-1); 
		if counts_IND[-1] < 0.01 or  ind_n < 0: ind_n = 0.0;
	else: 
		ind_n = c0-c1-sum_delta; ind_n /= math.pow(2,n-1)-1; 
		if counts_IND[-1] < 0.01 or  ind_n < 0: ind_n = 0.0;
		pcr_n = c1-sum1-ind_n; 
		if pcr_n < 0: pcr_n = 0.0

	PDF[-1] = pcr_n; PDF[0] = ind_n*math.pow(2,n-1);
	total = 0.0; unique = 0.0;
	for i in xrange(len(PDF)): total += PDF[i];
	#print >>sys.stderr, 'FINAL-res',PDF;
	for i in xrange(len(PDF)): PDF[i] /= total;  unique += len(list_part[i])*PDF[i]; 
	print >>sys.stderr, 'IND-n',ind_n,ind_n*math.pow(2,n-1)+sum0,'PCR-n',pcr_n,
	print >>sys.stderr, 'FINAL-res',PDF, unique,unique*c0,'--------------------\n'; 
	## return [ind_n,pcr_n] # two numbers... and also the weighted number of unique molecules for this cluster size 
	return [ind_n,pcr_n,unique,flag]; 


def calculate_values(clus2,clus3,duplicate_counts):

	for i in xrange(1,MAX_CLUSTER_SIZE):  ## scale counts to avoid filtering clusters for large 'k'
		if duplicate_counts[i][0] < 10: continue; 
		scalef = float(duplicate_counts[i][0])/duplicate_counts[i][1]; 
		#for j in xrange(1,4): duplicate_counts[i][j] = scalef*duplicate_counts[i][j]; 
		print >>sys.stderr, "dup-counts",i,duplicate_counts[i],float(duplicate_counts[i][1])/duplicate_counts[i][0]

	## from left off..
	total = duplicate_counts[1][0]; unique = duplicate_counts[1][0]; rate = 1.0; 
	for c in xrange(2,99): 
		if duplicate_counts[c][1] < 10 or duplicate_counts[c][2] < 5: break;
		total += duplicate_counts[c][0]*c;
		r = float(duplicate_counts[c][2])/(duplicate_counts[c][1]); fn = 1.0-math.pow(r,1.0/(c-1)); 
		ex = float(c-1)*fn*2 + 1; unique += ex*duplicate_counts[c][0];
		rate = 1.0-float(unique)/total;  print >>sys.stderr, "PCRrate...",c,rate,total,unique

		
	#### new approach implemented 03/02/16 
	total = duplicate_counts[1][0]; unique = duplicate_counts[1][0]; rate = 1.0; 
	counts_PCR = [duplicate_counts[1][1]]; counts_PCR.append(float(2*duplicate_counts[2][2]-duplicate_counts[2][1])); 
	counts_IND = [duplicate_counts[1][1]]; counts_IND.append(float(duplicate_counts[2][1]-duplicate_counts[2][2])); 	
	PCR_sum = counts_PCR[0] + counts_PCR[1]; 
	Lambda = [1.0,float(counts_PCR[1])/counts_PCR[0]]; 
	exp_unique = float(1.0*counts_PCR[1] + 4.0*counts_IND[1])/duplicate_counts[2][1]; 
	total += duplicate_counts[2][0]*2; unique += exp_unique*duplicate_counts[2][0];
	rate = 1.0-float(unique)/total;  print >>sys.stderr, "PCRrate-new...",2,rate,total,unique,exp_unique,counts_PCR[1],counts_IND[1],Lambda[1]

	flag = 0; 
	for c in xrange(3,100):
		if duplicate_counts[c][1] < 10 or duplicate_counts[c][2] < 5: break;
		total += duplicate_counts[c][0]*c;
		if flag ==0 and c < 20: 
			[ind_n,pcr_n,exp_unique,flag] = sorted_partitions(c,Lambda,counts_IND,duplicate_counts[c][1],duplicate_counts[c][2],duplicate_counts[c][3],counts_PCR[-1]); #sys.exit();
		if flag ==0: 
			unique += exp_unique*duplicate_counts[c][0]; 
			if pcr_n < 0.0001: Lambda.append(0.0); 
			else: Lambda.append(pcr_n/counts_PCR[0]);
			counts_PCR.append(pcr_n); counts_IND.append(ind_n); 
			rate = 1.0-float(unique)/total;  print >>sys.stderr, "PCRrate-new...",c,rate,total,unique,exp_unique,counts_PCR[-1],counts_IND[-1],Lambda[-1]
	
		else:
			r = float(duplicate_counts[c][2])/(duplicate_counts[c][1]); fn = 1.0-math.pow(r,1.0/(c-1)); 
			ex = float(c-1)*fn*2 + 1; unique += ex*duplicate_counts[c][0];
			rate = 1.0-float(unique)/total;  print >>sys.stderr, "PCRrate-simple...",c,rate,total,unique
			

	f1 = 0.0; f2=0.0;
	if clus2[0] + clus2[1] >= 10:	
		y1 = float(clus2[1])/(clus2[0]+clus2[1]+1); f1 = y1*2; print "Final results for",sys.argv[1], 'clusters-2',clus2[0],clus2[1],y1,'f2=',f1;
	if clus3[0] + clus3[1] >= 10: 
		y2 = 1.0 - math.sqrt(float(clus3[0])/(clus3[0]+clus3[1])); f2 = y2*2; print "Final results for",sys.argv[1], 'clusters-3',clus3[0],clus3[1],y2,'f3=',f2;
	return [f1,f2];

def estimate_allele_bias(read,counts):
	varlist = read[6]
	for i in xrange(len(varlist)):
		var = varlist[i].split(':'); 
		if len(var[2]) < 1: var[2] = ':'; 
		e1 = math.pow(10,-0.1*(ord(var[2][0])-33))/1;
		if (e1 < QV_THRESH and var[1] != '?'):
			if var[1] == '0': AB_counts[0] +=1; 
			if var[1] == '1': AB_counts[1] +=1; 
			if i ==0: AB_counts[2] +=1; counts[1][1] +=1; counts[1][2] += 1;

def llfunc(f):
	h = 0.5;
	a = float(AB_counts[0])/(AB_counts[0] + AB_counts[1]); b = 1.0-a;
	#a= 0.5; b = 1.0-a; 
	n = len(Coeffs); ll = 0.0; 
	for i in xrange(n):
		cvec = Coeffs[i];
		if (CSIZE ==2 or CSIZE == -1) and cvec[0] == 2: ll += math.log(cvec[1]*(1.0-f +(a*a+b*b)*f) + cvec[2]*f*2*a*b);
		if (CSIZE ==3 or CSIZE == -1) and cvec[0] == 3: ll += math.log(cvec[1]*(1.0-0.5*f)*(1.0-0.5*f) + cvec[2]*f*(1.0-0.25*f) );
		if (CSIZE ==4 or CSIZE == -1) and cvec[0] == 4: ll += math.log(cvec[1]*math.pow(1.0-0.5*f,3) + cvec[2]*(f*(1.0-f)*(1.0-f) + 1.5*f*f*(1.0-f) + f*f*f*0.5) + cvec[3]*(0.375*f*f*f + 0.75*f*f*(1.0-f) + 0.5*f*(1.0-f)*(1.0-f) ) );
	#print 'f',f,'ll',ll;
	return -1*ll; 

def process_duplicate_cluster_gen(readlist,s,csize,counts):
	
	varlist = []; 
	flag = 0; prev = []; mismatch =0; acounts = [0,0]; qarray = []
	for i in xrange(csize): 
		var = readlist[s+i][6][0].split(':'); 
		if len(var[2]) < 1: var[2] = ':'; 
		e = math.pow(10,-0.1*(ord(var[2][0])-33))
		if e >= QV_THRESH: flag = 1; 
		if var[1] == '?': flag = 1; 
		if len(prev) > 0 and var[0] != prev[0]: flag = 1; 
		if len(prev) > 0 and var[1] != prev[1] and var[0] == prev[0]: mismatch +=1; 
		prev = [var[0],var[1]];
		if var[1] == '0': acounts[0] +=1; 
		if var[1] == '1': acounts[1] +=1; 
		qarray.append([e,var[1]]); 
		#varlist.append([var[0],int(var[1]),var[2],e]); 
	if flag ==0: 
		if acounts[0] == csize or acounts[1] == csize: 
			if WEIGHTED ==0: counts[csize][2] +=1;  
			else:  
				qprod =1.0; 
				for i in xrange(csize): qprod *= 1.0-qarray[i][0]; 
				counts[csize][2] += qprod; 
			Coeffs.append([csize,1,0,0]);  
		elif acounts[0] == csize-1 or acounts[1] == csize-1: 
			Coeffs.append([csize,0,1,0]); 
			if WEIGHTED ==0: counts[csize][3] +=1; 
			if WEIGHTED ==1: 
				if csize ==2: counts[csize][2] += (1.0-qarray[0][0])*qarray[1][0]/3 + (1.0-qarray[1][0])*qarray[0][0]/3;
				else: 
					qprod =1.0; 
					for i in xrange(csize): 
						if qarray[i][1] == '1' and acounts[0] == csize-1: qprod *= qarray[i][0]/3; 
						elif qarray[i][1] == '0' and acounts[1] == csize-1: qprod *= qarray[i][0]/3; 
						else: qprod *= 1.0-qarray[i][0]; 
					counts[csize][2] += qprod; 
				counts[csize][3] +=1; 
				
		elif acounts[0] == csize-2 or acounts[1] == csize-2: Coeffs.append([csize,0,0,1]);
		counts[csize][1] += 1; ## total 
		#if mismatch ==0: counts[0] +=1; 
		#elif mismatch >= 1: counts[1] +=1; 
	else: pass;
	

## function to process cluster of reads with identical mapping coordinates and estimate # of true and false duplicates 
def process_duplicate_cluster(readlist,s,e,csize,clus2,clus3):

	## simple goal is to find two clusters and their relative sizes... 

	## BUG issue, if quality value == ':' -> splitting is incorrect...
	if csize ==2: ## most common size 
		list1 = readlist[s][6]; list2 = readlist[s+1][6]; 
		var1 = list1[0].split(':'); 
		if len(var1[2]) < 1: var1[2] = ':';
		var2 = list2[0].split(':'); 
		if len(var2[2]) < 1: var2[2] = ':';
#		if len(var1[2]) < 1 or len(var2[2]) < 1: print >>sys.stderr, 'ERROR',var1,var2,readlist[s],readlist[s+1]
		e1 = math.pow(10,-0.1*(ord(var1[2][0])-33))/1; e2 = math.pow(10,-0.1*(ord(var2[2][0])-33))/1;
		E1 = (1.0-e1)*(1.0-e2); E2 = (e1/3)*(1.0-e2) + (e2/3)*(1.0-e1);


		if var1[0] == var2[0] and e1 < QV_THRESH and e2 < QV_THRESH and var1[1] != '?' and var2[1] != '?':

			allele1 = int(var1[1]); allele2 = int(var2[1]); 
			#if random.random() < e1: allele1 = 1-allele1; 
			#if random.random() < e2: allele2 = 1-allele1; 
			## quality value filter ?? 
			if allele1 != allele2: 
				#if VERBOSE: print 'FALSE-DUP',E1,E2,var1[0],var1[1],var1[2],var2[0],var2[1],var2[2]; 
				if WEIGHTED ==0: clus2[1] += 1;
				else: 
					clus2[1] += E1; clus2[0] += E2; 
				if WEIGHTED ==0: Coeffs.append([2,0,1]);
				else: Coeffs.append([2,E2,E1]);
			if allele1 == allele2: 
				#if VERBOSE: print 'PCR-DUP',E2,E1,var1[0],var1[1],var1[2],var2[0],var2[1],var2[2]; 
				if WEIGHTED ==0: clus2[0] += 1;
				else: 
					clus2[0] += E1; clus2[1] += E2; 
				if WEIGHTED ==0: Coeffs.append([2,1,0]);
				else: Coeffs.append([2,E1,E2]);
				if allele1 ==0: clus2[2] +=1
			
			#if var1[1] != var2[1]: print 'FALSE-DUP',E1,E2,var1[0],var1[1],var1[2],var2[0],var2[1],var2[2]; clus2[1] +=E1; clus2[0] += E2;  
			#if var1[1] == var2[1]: print 'PCR-DUP',E2,E1,var1[0],var1[1],var1[2],var2[0],var2[1],var2[2]; clus2[0] +=E1; clus2[1] += E2;  


	if csize ==3: 
		list1 = readlist[s][6]; list2 = readlist[s+1][6]; list3 = readlist[s+2][6];
		var1 = list1[0].split(':'); 
		if len(var1[2]) < 1: var1[2] = ':'; # special case when quality value = ':' same as splitter
		var2 = list2[0].split(':'); 
		if len(var2[2]) < 1: var2[2] = ':'; 
		var3 = list3[0].split(':'); 
		if len(var3[2]) < 1: var3[2] = ':'; 
		e1 = math.pow(10,-0.1*(ord(var1[2][0])-33))/1; e2 = math.pow(10,-0.1*(ord(var2[2][0])-33))/1; e3 = math.pow(10,-0.1*(ord(var3[2][0])-33))/1
		vec = [e1/3,1.0-e1,e2/3,1.0-e2,e3/3,1.0-e3]; E1 = vec[1]*vec[3]*vec[5];

		mismatch = 0;
		## we should consider reads that cover multiple snps...
		if var1[0] == var2[0] and var2[0] == var3[0] and e1 < QV_THRESH and e2 < QV_THRESH and e3 < QV_THRESH and var1[1] != '?' and var2[1] != '?' and var3[1] != '?': ## first snp matches 
			allele1 = int(var1[1]); allele2 = int(var2[1]); allele3 = int(var3[1])

			if var1[1] == var2[1] and var2[1] == var3[1]: mismatch = 0; 
			#if allele1  == allele2 and allele2 == allele3: mismatch = 0; 
			elif var1[1] == var2[1]: mismatch = 3; 
			elif var2[1] == var3[1]: mismatch = 1; 
			elif var1[1] == var3[1]: mismatch = 2; 
			if mismatch >=1 and VERBOSE: print 'MIS-C3',var1,var2,var3; #clus3[1] += E1; clus3[0] += E2
			elif mismatch ==0 and VERBOSE: print 'MAT-C3',var1,var2,var3; #clus3[0] +=E1; clus3[1] += E2; 
			if mismatch >=1: 
				a = E1 + E1*(vec[0]/vec[1] + vec[2]/vec[3] + vec[4]/vec[5] - vec[mismatch*2-2]/vec[mismatch*2-1]); b= E1*vec[mismatch*2-2]/vec[mismatch*2-1];
				if WEIGHTED ==0: clus3[1] +=1; 
				else:  clus3[1] += a; clus3[0] += b; 
				if WEIGHTED ==0: Coeffs.append([3,0,1]);
				else: Coeffs.append([3,b,a]);
			elif mismatch ==0: 
				b = E1*(vec[0]/vec[1] + vec[2]/vec[3] + vec[4]/vec[5]);  
				if WEIGHTED ==0: clus3[0] +=1;  
				else: clus3[0] += E1; clus3[1] += b; 
				if WEIGHTED ==0: Coeffs.append([3,1,0]); 
				else: Coeffs.append([3,E1,b]);


	if VERBOSE and csize ==3:  
		j = s;
		if 'PCR' in readlist[j+1][0] and 'sampling' not in readlist[j+1][0] and 'PCR' in readlist[j+2][0] and 'sampling' not in readlist[j+2][0]: print 'c_000'; 
		elif 'PCR' not in readlist[j+1][0] and 'sampling' in readlist[j+1][0] and 'sampling' in readlist[j+2][0]  and 'PCR' not in readlist[j+2][0] and mismatch ==0: print 'c_0_0_0'; 
		elif 'PCR' in readlist[j+1][0] and 'sampling' in readlist[j+2][0] and mismatch ==0: print 'c_00_0'; 
		elif 'PCR' in readlist[j+2][0] and 'sampling' in readlist[j+1][0] and mismatch ==0: print 'c_00_0'; 
		for j in xrange(s,e): print readlist[j][0],
		if csize ==3:
			#if 'PCR' in readlist[j+1][0] and 'PCR' in readlist[j+2][0]: print '[000]',
			#if 'PCR' in readlist[j+1][0] and 'sampling' in readlist[j+2][0]: print '[000]',
			print 'mismatch',mismatch,  
			pass;
		print 'csize',csize,'\n', 


# process reads for each chromosome separately..
def process_readlist(readlist,duplicate_counts,clus2,clus3,chrom):

	clus2_sc = [0,0,0]; clus3_sc = [0,0,0]; # for single chromosome 
		
	reads = len(readlist);
	readlist.sort(key=lambda i: (i[0],i[4]));
	newreadlist = [];

	if reads > 1 and readlist[0][0] != readlist[1][0]: newreadlist.append(readlist[0]);
	i = 0; 
	while i < reads-1:
		if readlist[i][0] == readlist[i+1][0]: 
			## merge pair of reads and variant list into single list 
			newread = extra_functions.merge_readpair(readlist[i+1],readlist[i]); 
			if len(newread[6]) > 0: newreadlist.append(newread);  
			#print readlist[i-1],readlist[i],'\n';
			i +=2; 
		else: newreadlist.append(readlist[i]); i +=1; 


	newreadlist.sort(key=lambda i: (i[2],i[4],i[5])); ## sort readlist by position > IS >  flag:cigar | useful since for SE reads, we should require additional constraint on strand...
	reads = len(newreadlist);

	i = 1; csize = 1; start = 0; 
	while i < reads-1:
		flag1 = int(newreadlist[i][5].split(',')[0]); flag2 = int(newreadlist[i-1][5].split(',')[0]);  ## for single end reads if strand does not match
		if newreadlist[i][2] != newreadlist[i-1][2] or newreadlist[i][4] != newreadlist[i-1][4] or (flag1 + flag2 ==16): ## current read does not match previous cluster

			if csize ==1: estimate_allele_bias(newreadlist[i],duplicate_counts);
			if csize ==2 or csize == 3: process_duplicate_cluster(newreadlist,start,i,csize,clus2_sc,clus3_sc);
			if csize >=2 and csize < MAX_CLUSTER_SIZE: process_duplicate_cluster_gen(newreadlist,start,csize,duplicate_counts); 
			#elif csize ==5: process_duplicate_cluster_gen(newreadlist,start,csize,duplicate_counts); 
			if csize < MAX_CLUSTER_SIZE: duplicate_counts[csize][0] +=1; ## this is independent of sequencing errors, etc...
			csize = 0; start = i; 
		csize +=1; 
		i +=1;

	#for i in xrange(reads): print newreadlist[i]
	clus2[0] += clus2_sc[0]; clus2[1] += clus2_sc[1]; clus2[2] += clus2_sc[2]
	clus3[0] += clus3_sc[0]; clus3[1] += clus3_sc[1]; 

	if clus2_sc[0] + clus2_sc[1] > 0:
		y1 = float(clus2_sc[1])/(clus2_sc[0]+clus2_sc[1]); x1 = 2*y1; print >>sys.stderr, 'clus2-singlechrom',chrom,clus2_sc[0],clus2_sc[1],y1,x1,
		a0 = float(clus2[2])/(clus2[0]); a1 = 1.0-a0; a01 = 2*a0*a1; 
		y1 = float(clus2[1])/(clus2[0]+clus2[1]); x1 = 2*y1; print >>sys.stderr, 'clusters-2',clus2[0],clus2[1],y1,'f2=',x1,clus2[2],clus2[0]-clus2[2],a0,a1,a01
		print 'AB_counts',AB_counts[0],AB_counts[1],float(AB_counts[1])/(AB_counts[0]+AB_counts[1])

		#y2 = float(clus3_sc[1])/(2*clus3_sc[0]+2*clus3_sc[1]); x2 = y2*2; print >>sys.stderr, 'clus3',clus3_sc[0],clus3_sc[1],'f:',x2;
		#print >>sys.stderr, 'new-est',r0,r1,r0*(clus3[0]+clus3[1]),r1*(clus3[0]+clus3[1]); 
	if clus3_sc[0] + clus3_sc[1] > 0:
		y2 = 1.0 - math.sqrt(float(clus3_sc[0])/(clus3_sc[0]+clus3_sc[1])); x2 = y2*2; print >>sys.stderr, 'clus3',clus3_sc[0],clus3_sc[1],'f2:',x2,
		y2 = 1.0 - math.sqrt(float(clus3[0])/(clus3[0]+clus3[1])); x2 = y2*2; print >>sys.stderr, 'clus3',clus3[0],clus3[1],'f3=',x2;

	f4total = duplicate_counts[4][1];
	if f4total > 1:
		r = float(duplicate_counts[4][2])/(f4total); f4 = 2*(1.0-math.pow(r,1.0/3)); 
		#frac_0 = math.pow(1.0-f*0.5,3); 
		frac_1 = f4*(1.0-f4)*(1.0-f4) + 1.5*f4*f4*(1.0-f4) + f4*f4*f4*0.5; 
		print 'f4_counts',f4total,duplicate_counts[4][2],f4,'0001-counts',(frac_1)*float(f4total);


########################################################################################################
########################################## MAIN CODE ##########################################################
	


clus2 = [0,0,0]; clus3 = [0,0,0,0]; 

#### get list of heterozygous variants from PCRdups.c output file 02/22/16, initial pass through file
[het_vars,nhets] = extra_functions.get_variants(sys.argv[1]);


readlist = []; chrom = '-';
cluster_counts = [0]; new_counts = [0]; errorcounts= [0,0]; 
duplicate_counts = []; ## store how many times each read is duplicated
for c in xrange(MAX_CLUSTER_SIZE): duplicate_counts.append([0,0.0,0.0,0.0]);

filteredreads = 0;
File = open(sys.argv[1],'r');
for line in File: 
	read = line.strip().split();

	if read[0] == '#clusters': 
		for i in xrange(1,len(read)): cl = read[i].split(':'); cluster_counts.append(int(cl[1])); 
		continue;
	if read[0] == '#variant' or read[0] == '#VCF': continue; 
	if len(read) < 7: continue; 	

	if read[1] != chrom and len(readlist) > 0: 
		lr = len(readlist);
		print >>sys.stderr, "processing reads for chrom",chrom,'with',lr,'reads','filteredreads',filteredreads,errorcounts,float(errorcounts[1])*1.5/(errorcounts[0]+0.00001);
		#if '6' not in chrom and 'M' not in chrom and  'X' not in chrom and 'Y' not in chrom: process_readlist(readlist,duplicate_counts,clus2,clus3,chrom);  
		if 'M' not in chrom and  'X' not in chrom and 'Y' not in chrom: process_readlist(readlist,duplicate_counts,clus2,clus3,chrom);  
		for i in xrange(lr): readlist.pop();
		chrom = read[1]; 
		filteredreads =0;
	elif read[1] != chrom: chrom = read[1] 
	else: 
		## filter reads that do not overlap a variant in 'het_vars' (assuming nhets > 0) 
		if nhets > 0: 
			varlist = []; f= 0;
			for i in xrange(6,len(read)): 
				v = read[i].split(':'); 
				if (int(v[0])) in het_vars: 
					if v[1] == '?': errorcounts[1] +=1; 
					errorcounts[0] +=1; 
					varlist.append(read[i]); 
					f +=1; 
			if f > 0: 
				readlist.append([read[0],read[1],int(read[2]),int(read[3]),int(read[4]),(read[5]),varlist]);
			else: filteredreads +=1; 
		else: 
			readlist.append([read[0],read[1],int(read[2]),int(read[3]),int(read[4]),(read[5]),read[6:]]);

## final chromosome
if len(readlist) > 0 and 'X' not in chrom and 'Y' not in chrom: 
	lr = len(readlist);
	print >>sys.stderr, "processing reads for chrom",chrom,'with',lr,'reads'
	process_readlist(readlist,duplicate_counts,clus2,clus3,chrom);  
	for i in xrange(lr): readlist.pop();

###########################################################################################################
##### FINAL CALCULATIONS #######

## call function with clus2, clus3, duplicate_counts 

[f1,f2] = calculate_values(clus2,clus3,duplicate_counts); fest = f1; 
#CSIZE = 2; x_min = optimize.fminbound(llfunc,0.0001,0.9999,full_output=1,disp=1); print x_min[0]; fest = x_min[0];
#CSIZE = -1; x_min = optimize.fminbound(llfunc,0.0001,0.9999,full_output=1,disp=1); print x_min[0]; fjoint = x_min[0] #fest = x_min[0];
## calculate duplication rate; 
## print distribution of duplicate cluster sizes
total = 0; unique = 0; unique_corrected= [0.0,0.0];
for c in xrange(1,99): 
	## call function to estimate unique reads using 'f' estimated previously..
	#ex = expected_unique_reads(c-1,x1); 
	if duplicate_counts[c][1] < 10 or duplicate_counts[c][2] < 5: break;
	total += duplicate_counts[c][0]*c; unique += duplicate_counts[c][0];
	if c > 1: r = float(duplicate_counts[c][2])/(duplicate_counts[c][1]); fn = 1.0-math.pow(r,1.0/(c-1)); 
	else: fn = 0; 
	ex = float(c-1)*fn*2 + 1; unique_corrected[1] += ex*duplicate_counts[c][0];
	ex = float(c-1)*fest + 1; unique_corrected[0] += ex*duplicate_counts[c][0];
	print >>sys.stderr, c,duplicate_counts[c][0],ex,unique_corrected,'dup-rate',1.0-unique_corrected[1]/total,'C-total',duplicate_counts[c][1],'C-0000',duplicate_counts[c][2],fn
	#print >>sys.stderr, c,duplicate_counts[c],float(duplicate_counts[c])/duplicate_counts[c+1],ex,unique_corrected
	#print "Final results for",sys.argv[1], c,duplicate_counts[c],#,float(duplicate_counts[c])/duplicate_counts[c+1]
	#print ex,float(c-1)*fest+1;
print >>sys.stderr,"total reads-truncated",total,"unique",unique,"duplicate fraction",1.0-float(unique)/total,'PCR-dup-est_f1',1.0-unique_corrected[0]/total,'PCR-dup-est_fvar:',1.0-unique_corrected[1]/total,'singletons',AB_counts[2]


### this includes chrom 6 -> we should exclude ?? 
### analysis of total duplication clusters (not just that overlap SNPs) 
total = 0; unique = 0; unique_corrected= 0.0; 
for c in xrange(1,len(cluster_counts)): 
	new_counts.append(0.0);
	extra_functions.update_cluster_counts(new_counts,c,f1,cluster_counts[c]); 
	
for c in xrange(1,len(cluster_counts)): 
	total += cluster_counts[c]*c; unique += cluster_counts[c];
	#ex = expected_unique_reads(c-1,x1); 
	ex = float(c-1)*fest + 1; unique_corrected += ex*cluster_counts[c];
	if cluster_counts[c] <100: continue; 
	print 'FINAL-readduplicates',c,cluster_counts[c],ex,new_counts[c]

if total > 10: 
	print >>sys.stderr,"total reads",total,"unique",unique,"duplicate fraction",1.0-float(unique)/total,'PCR-dup-est:',1.0-unique_corrected/total;
	#print >>sys.stderr,'f2',f1,'f3',f2,'f_joint',fest;

if OPTIMIZE: 
	CSIZE = 2; x_min = optimize.fminbound(llfunc,0.0001,0.9999,full_output=1,disp=1); print x_min[0],'f2'
	CSIZE = 3; x_min = optimize.fminbound(llfunc,0.0001,0.9999,full_output=1,disp=1); print x_min[0],'f3'
	CSIZE = 4; x_min = optimize.fminbound(llfunc,0.0001,0.9999,full_output=1,disp=1); print x_min[0],'f4'
	#x_min = optimize.fminbound(llc4,0.0001,0.9999,full_output=1,disp=1); print x_min
