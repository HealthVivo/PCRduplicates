#! /usr/bin/env python
# AUTHOR VIKAS BANSAL last edited dec 2014
import sys, os, glob, string, subprocess,time, math, random


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


def get_variants(filename):
	het_vars = {};nhets=0;
	File = open(filename,'r');  print >>sys.stderr, "reading input file to obtain list of het variants";
	for line in File:
		if line[0] == '#' and line[1] == 'v' and line[2] == 'a': v = line.split(); het_vars[int(v[1])] = 1; nhets +=1;

		#if v[0] == '#variant': 
	File.close(); print >>sys.stderr, "identified",nhets,"het variants";
	return [het_vars,nhets];


