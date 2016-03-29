#! /usr/bin/env python
# AUTHOR VIKAS BANSAL last edited dec 2014
import sys, os, glob, string, subprocess,time, math, random

VERBOSE =0;

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
                elif varlist[i][0] == varlist[i+1][0] and varlist[i][1] != varlist[i+1][1]: newread[6].append(varlist[i][0] + ':' + '?' + ':'+varlist[i][2]);
                elif varlist[i][0] == varlist[i+1][0] and varlist[i][1] == varlist[i+1][1] and varlist[i][2] >= varlist[i+1][2]: newread[6].append(varlist[i][0] + ':' + varlist[i][1] + ':'+varlist[i][2]);
                elif varlist[i][0] == varlist[i+1][0] and varlist[i][1] == varlist[i+1][1] and varlist[i][2] < varlist[i+1][2]: newread[6].append(varlist[i][0] + ':' + varlist[i][1] + ':'+varlist[i+1][2]);
                if varlist[i][0] == varlist[i+1][0]: i +=2;
                else: i +=1;
        if i < len(varlist) and varlist[i][0] != varlist[i-1][0]: newread[6].append(varlist[i][0] + ':' + varlist[i][1] + ':'+varlist[i][2]);
        return newread;
        print newread;


def get_variants(filename,firstcol):
	het_vars = {};nhets=0;
	File = open(filename,'r');  print >>sys.stderr, "reading input file to obtain list of het variants";
	for line in File:
		if line[0] == '#':
			if line.strip().split()[0] == firstcol: v = line.split(); het_vars[int(v[1])] = 1; nhets +=1;
		#if v[0] == '#variant': 
	File.close(); print >>sys.stderr, "identified",nhets,"het variants";
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

