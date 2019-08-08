#!/usr/bin/python
#python AncestryHMM_summarize_posteriors_by_window.pyy /path/to/directory/containing/posterior/files/ file_name 10000 0.5 100


# ---- overhead   ---- #
import sys, glob, fileinput, re, math, os, ntpath
from collections import defaultdict
import numpy as np

if len(sys.argv[1:]) == 6:
    vals = sys.argv[1:]
    print "\narguments successfully loaded:\npath_2_files = ", vals[0], "\npopn_prefix = ", vals[1], "\nwind_size = ", vals[2], "\nmin_posterior = ", vals[3], "\nmin_n_sites = ", vals[4], "\noutput2 = ", vals[5]
else:
    print "\nSCRIPT NOT CALLED CORRECTLY:\n e.g.: python AncestryHMM_summarize_posteriors_by_window.py path_2_files(string) popn_prefix(string) wind_size(int) min_posterior(float) min_n_sites(int) out_dir(str)" 


def path_leaf(path):
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)



# ---- user defined command line arguments   ---- #

path        = str(vals[0])   # path to directory containing .posterior files -- e.g. '/pine/scr/e/a/earleyej/ancestry_hmm/aaron_against_yakRef/'
prefix      = str(vals[1])   # population ID prefix (i.e. base file name) -- e.g. 'hyb_nacl_8'
wind_size   = int(vals[2])   # window size (in bp) -- e.g. int(10000)
threshold_1 = float(vals[3]) # minimum posterior prob for a site to be considered ancestry-informative -- e.g. float(0.5)
threshold_2 = int(vals[4])   # minimum number of sites within a window to infer ancestry -- e.g. int(100) 
out_dir     = str(vals[5])   # path to where output should be written

## !! HARD-CODED - varies depending on the ploidy specified when running ancestry HMM !! ##
ancestry = [0.000, 0.125, 0.250, 0.375, 0.500, 0.625, 0.750, 0.875, 1.000]


# ---- go through and summarize each individual's posterior file ---- #

my_inds = glob.glob(path+prefix+'*.posterior')

for ind in my_inds:
	
	my_ind = path_leaf(ind).split(".")[0]
	my_ind_post = open(ind, 'r')
	out = open(out_dir + my_ind +'_wind'+ str(wind_size) +'_minSNPs'+ str(threshold_2) + '.post.wind.data', "a")
	out.write( 'chr' +','+ 'wind_start' +','+ 'mean_ancestry' + '\n' )
	
	my_wind = 1
	wind_ancestry = {}
	
	for i, line in enumerate(my_ind_post):
		
		if re.match('^chrom', line):
			continue
		
		else:
			line = line.rstrip('\n')
			line = re.split(r'\t', line)
			
			my_chr = line[0]
			my_pos = int(line[1])
						
			if my_chr not in wind_ancestry:
				wind_ancestry[my_chr] = defaultdict(list)
			
			if my_pos < wind_size:
				my_wind = 1
			else:
				my_wind = (my_pos // wind_size) * wind_size
			
			
			post_probs = map(float, line[2:])
			m = max(post_probs)
			
			if m > threshold_1:
				site_ancestry = ancestry[ [i for i, j in enumerate(post_probs) if j == m][0] ]
				wind_ancestry[my_chr][my_wind].append( site_ancestry )
			else:
				continue
	
	
	for chr in wind_ancestry:
		for wind in wind_ancestry[chr]:
			
			if len( wind_ancestry[chr][wind] ) > threshold_2:
				mean_wind_anc = float( sum(wind_ancestry[chr][wind]) ) / float( len(wind_ancestry[chr][wind]) )
				out.write( chr +','+ str(wind) +','+ str(mean_wind_anc) + '\n' )
			
			else:
				out.write( chr +','+ str(wind) +','+ 'NA' + '\n' )
	
	
	my_ind_post.close()
	out.close()

