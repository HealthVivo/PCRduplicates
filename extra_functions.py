#! /usr/bin/env python
# AUTHOR VIKAS BANSAL last edited dec 2014
import sys, os, glob, string, subprocess,time, math, random
from scipy import stats
from numpy import random

VERBOSE =0;

## use multinomial distribution to sample with replacement from counts for each class (cluster = i), Bootstrap of PCR duplication rate estimate 
def generate_bsample(duplicate_counts,n):
	new_counts= []; new_counts.append(duplicate_counts[0]); new_counts.append(duplicate_counts[1]); 

	for i in xrange(2,n):
		N = duplicate_counts[i][1]; 
		if N < 2: new_counts.append(duplicate_counts[i]); continue; 
		pvals = []; 
		for j in xrange(2,len(duplicate_counts[i])): pvals.append(float(duplicate_counts[i][j])/N); 
		#print >>sys.stderr, "PVALS",pvals,duplicate_counts[i]
		S = random.multinomial(N,pvals); 
		new_counts.append(duplicate_counts[i]); 
		for j in xrange(2,len(duplicate_counts[i])): new_counts[-1][j] = S[j-2]; 
		print 'BOOTSTRAP-counts',i,new_counts[i]; 
	return new_counts;
	

## this function should be applied to clusters estimated from entire data | not just clusters that overlap het SNPS !! | makes sense 
# n = cluster size - 1 | p = f/2 
def expected_unique_reads(n,p):
        #(1-f/2)^k.1 + 2k(1-f/2)^(k-1).f/2 + ....
        pow1 = math.pow(1.0-p,n); sum1=0;
        for i in xrange(0,n+1):
                sum1 += pow1*(i+1);
                if i < n: pow1 /= 1.0-p; pow1 *= p; pow1 *= (n-i); pow1 /= i+1;
                #print 'excal',n,p,sum1,pow1;
        return sum1;

def update_cluster_counts(new_counts,n,p,csize):
        if n ==1: new_counts[1] += csize;
        if n ==2: new_counts[1] += 2.0*p*csize; new_counts[2] += (1.0-p)*csize;
        if n ==3:
                new_counts[1] += 3.0*p*p*csize;
                new_counts[3] += (1.0-p)*(1.0-p)*csize;
                new_counts[2] += 2.0*p*(1.0-p)*csize; new_counts[1] += 2.0*p*(1.0-p)*csize;
        if n ==4:
                new_counts[1] += 4.0*math.pow(p,3)*csize;
                new_counts[4] += math.pow(1.0-p,3)*csize;
                pr = 3*p*p*(1.0-p)*csize; new_counts[1] += 2*pr; new_counts[2] += pr;
                pr = p*(1.0-p)*(1.0-p)*csize; new_counts[1] += 2*pr; new_counts[3] += 2*pr; new_counts[2] += 2*pr;



# function to merge paired-end read into single read with list of overlapping variants, shared variants accounted for
def merge_readpair(read1,read2):
        newread = [read1[0],read1[1],read1[2],read1[3],read1[4],read1[5] + '|' + read2[5],[]]; flag = 0;
        var1 = read1[6]; var2 = read2[6]; ## should be non-empty lists by default
        varlist = [];
        for v in var1: varlist.append(v.split(":"));
        for v in var2: varlist.append(v.split(":"));
        varlist.sort(key=lambda i: int(i[0]));
        #print varlist,'PE read';
        i = 0;
        while i < len(varlist)-1:
                if varlist[i][0] != varlist[i+1][0]: newread[6].append(varlist[i][0] + ':' + varlist[i][1] + ':'+varlist[i][2]);
		else: ## duplicate coverage of variant  
			if (varlist[i+1][1] == '?' or varlist[i+1][1] == '*'): newread[6].append(varlist[i][0] + ':' + varlist[i][1] + ':'+varlist[i][2]);
			elif (varlist[i][1] == '?' or varlist[i][1] == '*'): newread[6].append(varlist[i+1][0] + ':' + varlist[i+1][1] + ':'+varlist[i+1][2]);
			elif varlist[i][1] != varlist[i+1][1]: newread[6].append(varlist[i][0] + ':' + '?' + ':'+varlist[i][2]);
			elif varlist[i][1] == varlist[i+1][1] and varlist[i][2] >= varlist[i+1][2]: newread[6].append(varlist[i][0] + ':' + varlist[i][1] + ':'+varlist[i][2]);
			elif varlist[i][1] == varlist[i+1][1] and varlist[i][2] < varlist[i+1][2]: newread[6].append(varlist[i][0] + ':' + varlist[i][1] + ':'+varlist[i+1][2]);

                if varlist[i][0] == varlist[i+1][0]: i +=2;
                else: i +=1;
        if i < len(varlist) and varlist[i][0] != varlist[i-1][0]: newread[6].append(varlist[i][0] + ':' + varlist[i][1] + ':'+varlist[i][2]);
        return newread;
        print newread;

# function to analyze heterozygous variants from output of PCRdups.c  
def get_variants(filename,filtertype):
	het_vars = {}; nhets=0;
	variants = {}; varinfo = {}; 
	File = open(filename,'r');  print >>sys.stderr, "reading input file to obtain list of het variants";
	for line in File:
		if line[0] == '#':
			if line[1] == 'c' and line[2] == 'l': print >>sys.stderr, line
			if (line[1] == 'v' and line[2] == 'a') or (line[1] == 'f' and line[2] == 'i'): 
				variant = line.strip().split(); varid = int(variant[1]); info = variant[3] + ':'+variant[4] + '-' + `(int(variant[4])+1)`; 
				varinfo[varid] = info;
			pass;
			#if line.strip().split()[0] == firstcol: v = line.split(); het_vars[int(v[1])] = 1; nhets +=1;
		else: 
			read = line.strip().split();
			if len(read) < 5: continue;
			try: 
				pos = int(read[2]); IS = int(read[4]); cigar  = read[5]; 
				for i in xrange(6,len(read)): 
					var = read[i].split(':'); varid = int(var[0]); allele = var[1]; qv = var[2];
				rn = random.random(); # append random number to list so that when we sort, we can pick first element from each duplicate cluster to get random sample 
				try: 
					variants[varid].append([pos,IS,rn,allele,qv,cigar]); 
				except KeyError: variants[varid] = [[pos,IS,rn,allele,qv,cigar]]; 
			except ValueError: pass; 
	File.close(); 

	filtered = 0; high = [0.0000001,0.0,0.0];
	for var in variants.iterkeys():
		readlist = variants[var]; readlist.sort(key=lambda i: (i[0],i[1],i[2])); l = len(readlist);
		## find # of unique reads in list and then allele counts for a random set of unique reads from list 
		unique = 1; a0 =0; a1 = 0; a0t = 0; a1t = 0; coverage=1;
		if readlist[0][3] == '0': a0 +=1; a0t +=1; 
		elif readlist[0][3] == '1': a1 +=1;  a1t +=1; 

		if filtertype == "nobias": ## randomly select two reads from list, make sure they are not the same and they don't have missing alleles.... 
			a0 = 0; a1 = 0; coverage =0;  
			#print 'var:',var,len(readlist),'clusters:',readlist;
			for i in xrange(l):
				#if ord(readlist[i][4][0])-33 < 20: continue; 
				if readlist[i][3] == '0': a0 +=1; 
				elif readlist[i][3] == '1': a1 +=1;  
			coverage = a0 + a1; 
			if coverage > 1 and random.random()*coverage*coverage < 2*a0*a1: het_vars[var] = 1; nhets +=1;
                        else: filtered +=1;
			continue; 

		if filtertype == "exome.nodups":
			
			csize = 1; a0 = 0; a1 = 0; coverage =0;  
			for i in xrange(1,l):
				if readlist[i-1][0] != readlist[i][0] or readlist[i-1][1] != readlist[i][1]: pass; 
				
		elif filtertype == "exome.allreads":
			csize = 1; 
			print 'var:',var,len(readlist),'clusters:',
			for i in xrange(1,l):
				if readlist[i-1][0] != readlist[i][0] or readlist[i-1][1] != readlist[i][1]: 
					unique +=1; 
					if csize ==1: 
						if readlist[i][3] == '0': a0 +=1; 
						elif readlist[i][3] == '1': a1 +=1;  
						coverage +=1; 
					elif csize ==4: print 'cs:'+ `csize` + ':' + readlist[i-1][3] + readlist[i-2][3] + readlist[i-3][3] + readlist[i-4][3],
					csize = 0;
				if readlist[i][3] == '0': a0 +=1; 
				elif readlist[i][3] == '1': a1 +=1;  
				coverage +=1;
				csize +=1; 
			print '| unique',unique,a0,a1;
			
		else: 
			for i in xrange(1,l):
				if readlist[i][3] == '0': a0t +=1; 
				elif readlist[i][3] == '1': a1t +=1;  
				
				if readlist[i-1][0] != readlist[i][0] or readlist[i-1][1] != readlist[i][1]: 
					unique +=1; 
					if readlist[i][3] == '0': a0 +=1; 
					elif readlist[i][3] == '1': a1 +=1;  
					coverage +=1;


		ratio = float(a0)/(a0+a1+1.0e-8); pval = 1.0;
		if (ratio < 0.4 or ratio > 0.6) and coverage >= 10: pval = stats.binom_test(a0,a0+a1,0.5); 
		
		if filtertype == "none": het_vars[var] = 1; nhets +=1; ## no filtering at all..
		elif "exome" in filtertype or filtertype == 'dna': 
			if pval >= 0.001 and a1 >= 2 and coverage >= 8: het_vars[var] = 1; nhets +=1;
			else: filtered +=1;  
		elif filtertype == "rna" or filtertype == "rnaout": 
			if (ratio < 0.1 or ratio > 0.9) and coverage >=10 and var in varinfo: print 'FALSE-het',varinfo[var],len(readlist),unique,a0,a1,round(ratio,4),pval; high[1] +=1; 
			if coverage >=10: high[0] +=1; 
			if coverage >=30: high[2] +=1; 
			if (ratio < 0.1 or ratio > 0.9) and coverage >=30: 
				filtered +=1; pass;  
			else: het_vars[var] = 1; nhets +=1;
		elif filtertype == "rnastrict": 
			if (ratio < 0.2 or ratio > 0.8) and coverage >=25: filtered +=1; pass;  
			else: het_vars[var] = 1; nhets +=1;

		"""
		if coverage >= 10: 
			print 'VARCOV',var,len(readlist),unique,a0,a1,round(ratio,4),pval,
			if var in varinfo: print varinfo[var],
			if pval < 0.001: print 'flag';
			else: print '\n',
		"""
		
		#,a0t,a1t#,readlist[0:min(20,len(readlist))],'\n'
	
	print >>sys.stderr, "identified",nhets,"het variants",'filtered',filtered,'high-cov',high[0],high[1],high[1]/high[0],high[2];
	if filtertype == "rnaout": sys.exit();
	return [het_vars,nhets];


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


### march 28 2016

def calculate_expected_unique(n,counts):
	# length of counts is int(n/2) + 1 
	EXP = 1.0;  k= n; 
	while k >= 2: 
		if k > 2: EXP += 2.0*counts[1]/(counts[1] + k*counts[0]+1.0e-8); 
		else: EXP += 2.0*counts[1]/(counts[1] + counts[0]+1.0e-8);
		#print >>sys.stderr, "EXP",EXP,k,counts;

		newcounts = []; 
		for i in xrange((k-1)/2 + 1): newcounts.append(0.0); 

		for i in xrange(k/2 + 1): 
			ones = i; zeros = k-i; 
			### if it is 0011 -> we need to flip it before updating special case, when counts are equal, otherwise zeroes < ones 
			if ones == zeros: 
				newcounts[ones-1] += counts[i]; 
			else: 
				## reduces ones by 1 and add (i/k)*counts[i] to newcounts [xx] ## reduces zeros by 1 and add (k-i/k)*counts[i] to newcounts [xx] 
				if i > 0: newcounts[i-1] += counts[i]*float(i)/k;  
				newcounts[i] += counts[i]*float(k-i)/k; 

		for i in xrange((k-1)/2 + 1): counts[i] = newcounts[i]; 
		k -=1;
	return EXP;

def simple_method_estimates(duplicate_counts,MAX_CLUSTER_SIZE):
	### simplest method for estimation #### have to do the estimation recursively for each cluster... rather than adding the EXP's 
	EXP_SUM = 0.0; total_reads =0.0; unique_reads =0.0; 
	for i in xrange(1,MAX_CLUSTER_SIZE):
		if duplicate_counts[i][0] <= 50: continue; 
		if i ==1: EXP = 1.0000; 
		elif i ==2:
			EXP = 2.0*duplicate_counts[i][3]/(duplicate_counts[i][3] + duplicate_counts[i][2]); 
		else: 
			EXP = 2.0*duplicate_counts[i][3]/(duplicate_counts[i][3] + i*duplicate_counts[i][2]+1.0e-8); 
		EXP_SUM += EXP; total_reads += i*duplicate_counts[i][0]; unique_reads += EXP_SUM*duplicate_counts[i][0]; PCRrate = 1.0-unique_reads/total_reads;
		print >>sys.stderr, "dup-counts-simdata",i,'EXP_i',round(EXP,4),round(EXP_SUM,4),'PCR-rate',round(PCRrate,4),duplicate_counts[i]


def estimate_allele_bias(read,counts,QV_THRESH,AB_counts):
	varlist = read[6]
	for i in xrange(len(varlist)):
		var = varlist[i].split(':'); 
		if len(var[2]) < 1: var[2] = ':'; 
		e1 = math.pow(10,-0.1*(ord(var[2][0])-33))/1;
		if (e1 < QV_THRESH and var[1] != '?'):
			if var[1] == '0': AB_counts[0] +=1; 
			if var[1] == '1': AB_counts[1] +=1; 
			if i ==0: AB_counts[2] +=1; counts[1][1] +=1; counts[1][2] += 1;

def process_duplicate_cluster(readlist,s,e,csize,clus2,clus3,QV_THRESH,WEIGHTED):

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
			if allele1 == allele2: 
				#if VERBOSE: print 'PCR-DUP',E2,E1,var1[0],var1[1],var1[2],var2[0],var2[1],var2[2]; 
				if WEIGHTED ==0: clus2[0] += 1;
				else: 
					clus2[0] += E1; clus2[1] += E2; 
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
			elif mismatch ==0: 
				b = E1*(vec[0]/vec[1] + vec[2]/vec[3] + vec[4]/vec[5]);  
				if WEIGHTED ==0: clus3[0] +=1;  
				else: clus3[0] += E1; clus3[1] += b; 


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

