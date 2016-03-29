#! /usr/bin/env python
# AUTHOR VIKAS BANSAL last edited dec 2014
import sys, os, glob, string, subprocess,time, math, random
from scipy import optimize# import minimize_scalar

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
WEIGHTED =0;
OPTIMIZE = 0;
full_counts = []; ## cluster counts, total (0) 0000 matching (1)
SER = 0.01;
AB_counts=[0,0,0];

MAX_CLUSTER_SIZE = 100;
max_cluster_size = 0;

import extra_functions



## note that c2 is under-estimated... but we are not using c2 in the estimation, so it is okay 03/23/16
## pass Lambdas, counts_IND to this function 
def sorted_partitions(n,Lambda,counts_IND,c0,c1,c2,pcr_prev):
	## ignore first and last element of  list since they need to be estimated... we can calculate the sum of remaining elements... 
	it = extra_functions.partitions(n); list_part = [];
	for l in it: list_part.append(l);

	## coeff1 is for permutation of the elements in each partition, coeff3 is for 2^k-1 (label 0/1 for each element in partition) 
	## coeff2 is for  ??? not needed 03/16/16, set to 1.0 or remove 

	print >>sys.stderr,"\n----calling partiton function with n=",n,'C(....)=',c0,'C(00..00)=',c1,'C(00..001)=',c2#,c1/c0,c2/c0;
	sum0 = 0.0; sum1 = 0.0; PDF = []; 
	sum_0001 = 0.0; # 0000001 expression, we are not using it for calculation but can use to get alternate estimates... 
	sum_delta = 0.0; # sum0-sum1; 
	for partition in list_part: 
		partition.sort(); l = len(partition); coeff3 = math.pow(2,l-1);
		if l == n or l == 1: PDF.append(0.0); continue; print >>sys.stderr, partition; continue; 
		prod = counts_IND[l-1]; 
		coeff2 = 1.0; N=0; coeff1 = extra_functions.fact(l); c = 1; 
		singles=0;
		for i in xrange(len(partition)): 
			if partition[i] == 1: singles +=1; 

		for i in xrange(1,len(partition)): 
			if partition[i] != partition[i-1]: coeff1 /= extra_functions.fact(c); c = 1; 
			else: c +=1;
		coeff1 /= extra_functions.fact(c);  

		for el in partition: 
			#for j in xrange(el,1,-1): prod *= Lambda[j-1]; print 'l_' + `j-1`,
			prod *= Lambda[el-1]; 
			#N += el-1; coeff2 /= fact(el-1);
		#coeff2 *= fact(N); 
		prod *= coeff1*coeff2; 
	
		sum1 += prod; sum_delta += prod*(coeff3-1);  sum_0001 += prod*singles; sum0 += prod*coeff3;
		PDF.append(prod*coeff3);
		#print >>sys.stderr, '|',partition,'C1:',coeff1,'C2:',coeff2,'C3:2^k-1',coeff3,'|',prod,'sum1',sum1

	flag =0;
	if sum1 > c1: print >>sys.stderr,"flag sum1 exceeds c1 count from data",c1-sum1; flag =0;
	# if ind_(n-1) == 0 -> ind_n should be zero, if pcr_(n-1) = 0, then pcr_n = 0
	print >>sys.stderr, 'n','ind_est_1',c1-sum1,c0-sum0,math.pow(2,n-1),'three-sums',sum0,sum1,sum_0001;

	if c2 > sum_0001: 
		ind_n1 = (c2-sum_0001)/n; 
		pcr_n1 = c1-sum1-ind_n1;
		pcr_n2 = c0-sum0-ind_n1*math.pow(2,n-1);
		print >>sys.stderr,"alt estimate of ind_n is:",ind_n1,'pcr_n:',pcr_n1,pcr_n2 
	else: print >>sys.stderr,"ERROR"; 

	if pcr_prev < 0.01: 
		pcr_n  = 0.0; ind_n = c0-sum0;
		## ind_n can be easily estimated from c0 only...
		ind_n /= math.pow(2,n-1); 
		if counts_IND[-1] < 0.01 or  ind_n < 0: ind_n = 0.0;
	else: 
		ind_n = c0-c1-(sum0-sum1); ind_n /= math.pow(2,n-1)-1; 
		if counts_IND[-1] < 0.01 or  ind_n < 0: ind_n = 0.0;
		pcr_n = c1-sum1-ind_n; 
		if pcr_n < 0: pcr_n = 0.0
	
	#if c2 > sum_0001: ind_n = ind_n1; pcr_n = pcr_n;

	PDF[-1] = pcr_n; PDF[0] = ind_n*math.pow(2,n-1);
	total = 0.0; unique = 0.0;
	for i in xrange(len(PDF)): total += PDF[i];
	#print >>sys.stderr, 'FINAL-res',PDF;
	for i in xrange(len(PDF)): PDF[i] /= total;  unique += len(list_part[i])*PDF[i]; 
	if len(PDF) < 12: print >>sys.stderr, 'FINAL-res',PDF;
	print >>sys.stderr, 'IND-n',ind_n,'IND-n-scaled',ind_n*math.pow(2,n-1),'PCR-n',pcr_n,
	print >>sys.stderr, 'unique-exp:',round(unique,4),n,'PCR-dup-rate-forcluster:',1.0-unique/n,'--------'; 
	## return [ind_n,pcr_n] # two numbers... and also the weighted number of unique molecules for this cluster size 
	return [ind_n,pcr_n,unique,flag]; 


def calculate_values(clus2,clus3,duplicate_counts1,cluster_counts):

	duplicate_counts = []; ## store how many times each read is duplicated
        for c in xrange(MAX_CLUSTER_SIZE): duplicate_counts.append([0,0.0,0.0,0.0]);
	FL = 1;

        if max_cluster_size > 1 and FL ==1:
                for i in xrange(1,min(max_cluster_size,MAX_CLUSTER_SIZE)):  ## scale counts to avoid filtering clusters for large 'k'
                        if cluster_counts[i] < 100 or duplicate_counts1[i][0] < 100: continue;
                        #scalef = float(duplicate_counts[i][0])/duplicate_counts[i][1]; 
                        #for j in xrange(1,4): duplicate_counts[i][j] = scalef*duplicate_counts[i][j]; 
                        r0 = float(duplicate_counts1[i][2])/(duplicate_counts1[i][1]+1e-8); 
                        duplicate_counts[i][0] = cluster_counts[i]; duplicate_counts[i][1] = cluster_counts[i]; duplicate_counts[i][2] = int(r0*cluster_counts[i]); 
			if len(duplicate_counts1[i]) >=4: r1 = float(duplicate_counts1[i][3])/(duplicate_counts1[i][1]+1e-8); duplicate_counts[i][3] = int(r1*cluster_counts[i]);
                        print >>sys.stderr, "dup-counts",i,duplicate_counts1[i][0],duplicate_counts[i][0],round(float(duplicate_counts1[i][0])/duplicate_counts[i][0],4),round(r0,4),duplicate_counts1[i],duplicate_counts[i]
        else: ## for simulated data, cluster_counts is empty
                for i in xrange(1,MAX_CLUSTER_SIZE):
                        for j in xrange(min(4,len(duplicate_counts1[i]))): duplicate_counts[i][j] = duplicate_counts1[i][j];
                        if duplicate_counts1[i][0] >= 10: print >>sys.stderr, "dup-counts-simdata",i,duplicate_counts[i]
	print >>sys.stderr, "\n",


	#### new approach implemented 03/02/16 
	total = duplicate_counts[1][1]; unique = duplicate_counts[1][1]; rate = 1.0;  uniqlist = [1]; 
	if max_cluster_size > 1: total_f = cluster_counts[1]; unique_f = cluster_counts[1]; ## full data, all clusters
	counts_PCR = [duplicate_counts[1][1]]; counts_PCR.append(float(2*duplicate_counts[2][2]-duplicate_counts[2][1]));  # total = [00] + 2.[0][0], (00) = [00] + [0][0] 
	counts_IND = [duplicate_counts[1][1]]; counts_IND.append(float(duplicate_counts[2][1]-duplicate_counts[2][2])); # total - 00	
	PCR_sum = counts_PCR[0] + counts_PCR[1]; 
	Lambda = [1.0,float(counts_PCR[1])/counts_PCR[0]]; 
	exp_unique = float(1.0*counts_PCR[1] + 4.0*counts_IND[1])/duplicate_counts[2][1]; 
	total += duplicate_counts[2][1]*2; unique += exp_unique*duplicate_counts[2][1];
	if max_cluster_size > 2: total_f += cluster_counts[2]*2; unique_f += exp_unique*float(cluster_counts[2]); ## full data, all clusters
	rate = 1.0-float(unique)/total;  print >>sys.stderr, "PCRrate-new...",2,round(rate,4),total,unique,'PCR-rate-forcluster',1.0-exp_unique/2,exp_unique,'PCR-1',counts_PCR[1],'IND-1',counts_IND[1],counts_IND[1]*4,'lambda_2:',Lambda[1]
	uniqlist.append(round(exp_unique,4));

	flag = 0; 
	for c in xrange(3,100):
		if duplicate_counts[c][1] < 10 or duplicate_counts[c][2] < 5: break;
		total += duplicate_counts[c][1]*c; 
		if max_cluster_size > c: total_f += cluster_counts[c]*c;
		if flag ==0 and c < 20: 
			[ind_n,pcr_n,exp_unique,flag] = sorted_partitions(c,Lambda,counts_IND,duplicate_counts[c][1],duplicate_counts[c][2],duplicate_counts[c][3],counts_PCR[-1]); #sys.exit();
		print >>sys.stderr,"small-counts",duplicate_counts1[c]
		if flag ==0: 
			uniqlist.append(round(exp_unique,4));
			unique += exp_unique*duplicate_counts[c][1]; 
			if pcr_n < 0.0001: Lambda.append(0.0); 
			else: Lambda.append(pcr_n/counts_PCR[0]);
			counts_PCR.append(pcr_n); counts_IND.append(ind_n); 
			rate = 1.0-float(unique)/total;  
			print >>sys.stderr, "PCRrate-new...",c,round(rate,4),total,unique,exp_unique,counts_PCR[-1],counts_IND[-1],'lambda_'+`c`+':',Lambda[-1],'full-data',
			if max_cluster_size > c: 
				unique_f += exp_unique*float(cluster_counts[c]); rate_f = 1.0-float(unique_f)/total_f;  
				print >>sys.stderr, total_f,unique_f,rate_f;
			else: print >>sys.stderr, '\n',
	
		else:
			r = float(duplicate_counts[c][2])/(duplicate_counts[c][1]); fn = 1.0-math.pow(r,1.0/(c-1)); 
			ex = float(c-1)*fn*2 + 1; unique += ex*duplicate_counts[c][1]; 
			uniqlist.append(round(ex,4));
			rate = 1.0-float(unique)/total;  
			print >>sys.stderr, "data-no-fit PCRrate-simple...",c,rate,total,unique,
			if max_cluster_size > c: 
				unique_f += ex*float(cluster_counts[c]);  rate_f = 1.0-float(unique_f)/total_f;  
				print >>sys.stderr, total_f,unique_f,rate_f;
			else: print >>sys.stderr, '\n',
			
	print >>sys.stderr, '\nuniqlist:',uniqlist,'\n';
	for i in xrange(1,len(Lambda)):
		print >>sys.stderr,"L_"+`i`+':',Lambda[i]#(i+1)*Lambda[i+1]/Lambda[i];

	#extra_functions.simple_method_estimates(duplicate_counts,MAX_CLUSTER_SIZE);

	
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
		if acounts[0] == csize or acounts[1] == csize: ## 000000 or 11111 case  
			if WEIGHTED ==0: counts[csize][2] +=1;  
			else:  
				qprod =1.0; 
				for i in xrange(csize): qprod *= 1.0-qarray[i][0]; 
				counts[csize][2] += qprod; 
				for i in xrange(csize): counts[csize][3] += qprod*qarray[i][0]/(1.0-qarray[i][0]); 

		elif acounts[0] == csize-1 or acounts[1] == csize-1: ## 00001 or 11110 case 
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
		else:
			for j in xrange(2,csize/2+1): 
				if acounts[0] == csize-j or acounts[1] == csize-j: counts[csize][j+2] +=1;
		counts[csize][1] += 1; ## total 
	else: pass;
	

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

			if csize ==1: extra_functions.estimate_allele_bias(newreadlist[i],duplicate_counts,QV_THRESH,AB_counts);
			if csize ==2 or csize == 3: extra_functions.process_duplicate_cluster(newreadlist,start,i,csize,clus2_sc,clus3_sc,QV_THRESH,WEIGHTED);
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
		y1 = float(clus2_sc[1])/(clus2_sc[0]+clus2_sc[1]); x1 = 2*y1; print >>sys.stderr, 'clus2-singlechrom',chrom,clus2_sc[0],clus2_sc[1],y1;
		a0 = float(clus2[2])/(clus2[0]); a1 = 1.0-a0; a01 = 2*a0*a1; 
		y1 = float(clus2[1])/(clus2[0]+clus2[1]); x1 = 2*y1; print >>sys.stderr, 'clusters-2 00:',clus2[0],'01:',clus2[1],'f2=',1.0-y1,'A0:',a0,'dup-counts',duplicate_counts[2][1],duplicate_counts[2][2],float(duplicate_counts[2][2])/duplicate_counts[2][1],'\n';
		#clus2[2],clus2[0]-clus2[2],'A0',a0,'A1',a1,'\n'
		#print 'AB_counts',AB_counts[0],AB_counts[1],float(AB_counts[1])/(AB_counts[0]+AB_counts[1])

		#y2 = float(clus3_sc[1])/(2*clus3_sc[0]+2*clus3_sc[1]); x2 = y2*2; print >>sys.stderr, 'clus3',clus3_sc[0],clus3_sc[1],'f:',x2;
		#print >>sys.stderr, 'new-est',r0,r1,r0*(clus3[0]+clus3[1]),r1*(clus3[0]+clus3[1]); 
	if clus3_sc[0] + clus3_sc[1] > 0:
		y2 = 1.0 - math.sqrt(float(clus3_sc[0])/(clus3_sc[0]+clus3_sc[1])); x2 = y2*2; #print >>sys.stderr, 'clus3',clus3_sc[0],clus3_sc[1],'f2:',x2,
		y2 = 1.0 - math.sqrt(float(clus3[0])/(clus3[0]+clus3[1])); x2 = y2*2; #print >>sys.stderr, 'clus3',clus3[0],clus3[1],'f3=',x2;

	f4total = duplicate_counts[4][1];
	if f4total > 1:
		r = float(duplicate_counts[4][2])/(f4total); f4 = 2*(1.0-math.pow(r,1.0/3)); 
		#frac_0 = math.pow(1.0-f*0.5,3); 
		frac_1 = f4*(1.0-f4)*(1.0-f4) + 1.5*f4*f4*(1.0-f4) + f4*f4*f4*0.5; 
		#print 'f4_counts',f4total,duplicate_counts[4][2],f4,'0001-counts',(frac_1)*float(f4total);


########################################################################################################
########################################## MAIN CODE ##########################################################
	

clus2 = [0,0,0]; clus3 = [0,0,0,0]; 

#### get list of heterozygous variants from PCRdups.c output file 02/22/16, initial pass through file
[het_vars,nhets] = extra_functions.get_variants(sys.argv[1],"#variant");


readlist = []; chrom = '-';
cluster_counts = [0]; new_counts = [0]; errorcounts= [0,0]; 
duplicate_counts = []; ## store how many times each read is duplicated
for c in xrange(MAX_CLUSTER_SIZE): 
	duplicate_counts.append([0,0.0]); 
	for i in xrange(c/2 + 1): duplicate_counts[-1].append(0.0); 
	#duplicate_counts.append([0,0.0,0.0,0.0,0.0]);
## for i=2, total_uf, total, 00,01  | 000 001 | 0000 0001 0011 | 00000 00001 00011 | 000000 000001 000011 000111 size =  n/2 + 1 

filteredreads = 0;
File = open(sys.argv[1],'r');
for line in File: 
	read = line.strip().split();

	if read[0] == '#clusters': 
		for i in xrange(1,len(read)): 
			cl = read[i].split(':'); cluster_counts.append(int(cl[1])); max_cluster_size +=1; 
		continue;
	if read[0] == '#variant' or read[0] == '#VCF' or read[0] == '#filtered': continue; 
	if len(read) < 7: continue; 	

	if read[1] != chrom and len(readlist) > 0: 
		lr = len(readlist);
		print >>sys.stderr, "processing reads for chrom",chrom,'with',lr,'reads','filteredreads',filteredreads,errorcounts,float(errorcounts[1])*1.5/(errorcounts[0]+0.00001);
		#if '6' not in chrom and 'M' not in chrom and  'X' not in chrom and 'Y' not in chrom: process_readlist(readlist,duplicate_counts,clus2,clus3,chrom);  
		if 'M' not in chrom and  'X' not in chrom and 'Y' not in chrom: 
			if '16' in chrom or '6' not in chrom: process_readlist(readlist,duplicate_counts,clus2,clus3,chrom);  
			else: print >>sys.stderr, "#########not calling function on chrom###########",chrom
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

calculate_values(clus2,clus3,duplicate_counts,cluster_counts); 



## call function with clus2, clus3, duplicate_counts 
f1 = 0.0; f2=0.0;
if clus2[0] + clus2[1] >= 10:	
	y1 = float(clus2[1])/(clus2[0]+clus2[1]+1); f1 = y1*2; print "Final results for",sys.argv[1], 'clusters-2',clus2[0],clus2[1],y1,'f2=',f1;
if clus3[0] + clus3[1] >= 10: 
	y2 = 1.0 - math.sqrt(float(clus3[0])/(clus3[0]+clus3[1])); f2 = y2*2; print "Final results for",sys.argv[1], 'clusters-3',clus3[0],clus3[1],y2,'f3=',f2;
#return [f1,f2];
fest = f1; 

#CSIZE = 2; x_min = optimize.fminbound(llfunc,0.0001,0.9999,full_output=1,disp=1); print x_min[0]; fest = x_min[0];
#CSIZE = -1; x_min = optimize.fminbound(llfunc,0.0001,0.9999,full_output=1,disp=1); print x_min[0]; fjoint = x_min[0] #fest = x_min[0];
## calculate duplication rate; 
## print distribution of duplicate cluster sizes
total = 0; unique = 0; unique_corrected= [0.0,0.0];
print >>sys.stderr, "\n ------------simple procedure for estimating PCR-rate-----------";
for c in xrange(1,99): 
	## call function to estimate unique reads using 'f' estimated previously..
	#ex = expected_unique_reads(c-1,x1); 
	if duplicate_counts[c][1] < 10 or duplicate_counts[c][2] < 5: break;
	total += duplicate_counts[c][0]*c; unique += duplicate_counts[c][0];
	if c > 1: r = float(duplicate_counts[c][2])/(duplicate_counts[c][1]); fn = 1.0-math.pow(r,1.0/(c-1)); 
	else: fn = 0; 
	ex = float(c-1)*fn*2 + 1; unique_corrected[1] += ex*duplicate_counts[c][0];
	ex = float(c-1)*fest + 1; unique_corrected[0] += ex*duplicate_counts[c][0];
	print >>sys.stderr, c,duplicate_counts[c][0],round(ex,4),'PCR-forclus:',round(1.0-ex/c,4),'dup-rate',round(1.0-unique_corrected[1]/total,4),'C-total',duplicate_counts[c][1],'C-0000',duplicate_counts[c][2],round(fn,4)
	#print >>sys.stderr, c,duplicate_counts[c],float(duplicate_counts[c])/duplicate_counts[c+1],ex,unique_corrected
	#print "Final results for",sys.argv[1], c,duplicate_counts[c],#,float(duplicate_counts[c])/duplicate_counts[c+1]
	#print ex,float(c-1)*fest+1;
print >>sys.stderr,"\ntotal reads-truncated",total,"unique",unique,"duplicate fraction",1.0-float(unique)/total,'PCR-dup-est_f1',1.0-unique_corrected[0]/total,'PCR-dup-est_fvar:',1.0-unique_corrected[1]/total,'singletons',AB_counts[2]


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
