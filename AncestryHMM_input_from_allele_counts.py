#!/usr/bin/python

'''
generate primary input file for "Ancestry_HMM" software.
A.A. Comeault
29 July 2019

This script requires (i.e. it is designed to only work with the following file formats!):
(1) allele frequencies estimated for both parental species (P1 and P2) 
	as generated from vcftools "--freq" tool.
(2) "sync" file for hybrid swarms sequenced as a pool of individuals.
	This script is designed to take the .sync file format for
	a pool of indiduals from an admixed population as generated with the 
	"mpileup2sync.pl" script distributed with the Popoolation2 set of scripts.
	*Was developed off of output from Popoolation2 v1201*

!!	NOTE that recombination rate is currently assumed to be uniform,
	and the genetic distance between markers is set to ave. of 1 recombination 
 	event every 20,000,000bp ~ 1cM == 200,000bp ~ 0.000005cM/bp (my_rec_rate parameter).

e.g. usage:
python AncestryHMM_input_from_allele_counts.py /path/and/af/file/for/P1_allele_freq.frq /path/and/af/file/for/P1_allele_freq.frq /path/and/af/file/for/hybrid_readcounts_sync.txt /path/to/output/directory/
'''

## ---- python packages and user-defined arguments ---- ##

import sys, glob, fileinput, re, math, os
from pprint import pprint
from subprocess import call
from collections import defaultdict
import numpy as np
from itertools import imap


if len(sys.argv[1:]) == 4:
    vals = sys.argv[1:]
    print "arguments successfully loaded!"
    print "P1-AFs:   ", vals[0]
    print "P2-AFs:   ", vals[1]
    print "hybrid allele counts:   ", vals[2]
    print "output directory:   ", vals[3]
else:
    print "\nUsage: <python AncestryHMM_input_from_allele_counts.py> [/full/path/to/parentalAFs/P1.freq] [/full/path/to/parentalAFs/P2.freq file] [/full/path/to/sync/file/hybrids.sync file] [/path/to/output/dir/]"


# ---- parse arguments ---- #

my_p1_afs       = str(vals[0])
my_p2_afs       = str(vals[1])
my_hyb_counts   = str(vals[2])
out_dir         = str(vals[3])
my_rec_rate     = 0.000005  ## !! HARD CODED uniform per base-pair recombination rate !! ##


## ---- get parental allele frequencies ---- ##
# sites are restricted to chromosomes (arms) with the following names:
## !! HARD CODED so needs to be changed to work with different genomes !! ##
parental_afs = {'2L': {},
				'2R': {},
				'3L': {},
				'3R': {},
				'X': {},
				'4': {}}

# p1 AFs
p1_af_file = open( my_p1_afs )

for i, line in enumerate(p1_af_file):
	if i > 0:
		line = line.rstrip('\n')
		line = re.split(r'\t', line)
		
		if line[0] in parental_afs and int(line[2]) == 2 and int(line[3]) > 1:
			parental_afs[line[0]][line[1]] = {}
			parental_afs[line[0]][line[1]][ re.split(r':', line[4])[0] ] = [ int(round( int(line[3]) * float(re.split(r':', line[4])[1]) )) ]
			parental_afs[line[0]][line[1]][ re.split(r':', line[5])[0] ] = [ int(round( int(line[3]) * float(re.split(r':', line[5])[1]) )) ]
			
		else:
			next

p1_af_file.close()

print "done with P1 allele counts"

# P2 AFs (only sites sites that overlap with P1):
p2_af_file = open( my_p2_afs )

for i, line in enumerate(p2_af_file):
	if i > 0:
		line = line.rstrip('\n')
		line = re.split(r'\t', line)
		
		if line[0] in parental_afs and line[1] in parental_afs[line[0]] and int(line[2]) == 2 and int(line[3]) > 1:
			parental_afs[line[0]][line[1]][ re.split(r':', line[4])[0] ].append( int(round( int(line[3]) * float(re.split(r':', line[4])[1]) )) )
			parental_afs[line[0]][line[1]][ re.split(r':', line[5])[0] ].append( int(round( int(line[3]) * float(re.split(r':', line[5])[1]) )) )
			
		else:
			next

p2_af_file.close()

print "done with P2 allele counts"



## ---- get read-counts for hybrid swarm (pooled-sequencing reads) and standardize to number of haplotypes in pool ---- ##
allele_map={ 'A':0,
			 'T':1,
			 'C':2,
			 'G':3,
			 'N':4,
			 'del':5 }

def get_allele_counts( read_counts, allele ):
	return int( read_counts[ allele_map[allele] ] )


hyb_sync_file = open(my_hyb_counts)

for i, line in enumerate(hyb_sync_file):
	line = line.rstrip('\n')
	line = re.split(r'\t', line)
	
	if line[0] in parental_afs and line[1] in parental_afs[line[0]]:
		
		for allele in sorted(parental_afs[line[0]][line[1]].keys()):
			parental_afs[line[0]][line[1]][ allele ].append( get_allele_counts( re.split(r':', line[3]), allele ) )
		
	else:
		next

hyb_sync_file.close()
		



## ---- write output ---- ##
out = open( out_dir + re.split('\.', os.path.basename( my_hyb_counts ))[0] + '_ancHMM.input', "a" )

print "writing input panel with allele counts to " + out_dir + re.split('\.', os.path.basename( my_hyb_counts ))[0] + "_ancHMM.input"


for chrom in sorted(parental_afs.keys()):
	
	my_poses = sorted(map( int, parental_afs[chrom].keys() ))
	pos_countr=-1
	
	for pos in map(str, my_poses):
		
		pos_countr+=1
		allele_counts=[0,0,0,0,0,0]
		countr=0
		
		for allele in sorted(parental_afs[chrom][pos].keys()):
			
			if countr == 0 and sum(parental_afs[chrom][pos][allele]) > 1 and len(parental_afs[chrom][pos][allele]) == 3:
				allele_counts[0]=parental_afs[chrom][pos][allele][0]
				allele_counts[2]=parental_afs[chrom][pos][allele][1]
				allele_counts[4]=parental_afs[chrom][pos][allele][2]
				
				countr+=1
				
			elif countr == 1 and sum(parental_afs[chrom][pos][allele]) > 1 and len(parental_afs[chrom][pos][allele]) == 3:
				allele_counts[1]=parental_afs[chrom][pos][allele][0]
				allele_counts[3]=parental_afs[chrom][pos][allele][1]
				allele_counts[5]=parental_afs[chrom][pos][allele][2]
				
				af_dif=abs( (float(allele_counts[0]) / sum(allele_counts[0:2])) - (float(allele_counts[2]) / sum(allele_counts[2:4])) )
				
				if af_dif > 0.5: ## !! HARD CODED - minimum difference in allele freqeucny between parental panels for a site to be considered !! ##
					if pos_countr > 0:
						cm_dist = float(my_poses[pos_countr] - my_poses[pos_countr-1])*float(my_rec_rate)
					else:
						cm_dist = 0
					out.write( '\t'.join([chrom, pos]) +'\t'+ '\t'.join( map(str, allele_counts[0:4]) ) +'\t'+ str(cm_dist) +'\t'+ '\t'.join( map(str, allele_counts[4:6]) ) +'\n' )
				
			else:
				next
			
	
out.close()

