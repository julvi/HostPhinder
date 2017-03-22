#!/usr/bin/env python

# Takes sorted result file with hosts, decision measure and alpha and returns accn\tpredicted_host
# input example:
# ./get_host_alpha.py -a 2.797 -f ../HostPhinder/results/paropt/234_subpart/Y13918.fsa_species_pred_16mers_evalue0.05_host_sorted -c FALSE

# Julia Villarroel
# last modified: 150219

import argparse
import sys
from compiler.ast import flatten
import os
import re

parser= argparse.ArgumentParser(description='Get prediction and result from ALL table for best result.')

parser.add_argument('-f','--file', type=str, help="File that contains results \
with host SORTED by user decision")

parser.add_argument('-d','--decision', type=str, default='coverage', 
help='Value which should be chosen to find majority host: Score, z, frac_q, frac_d or coverage')

parser.add_argument('-a','--alpha', type=float, default=6.0,
 help='alpha value')

parser.add_argument('-c', '--clustering', choices=['TRUE', 'FALSE'],
 default='TRUE', help='performs homology clustering by default.\
 If set to FALSE it does not')

args = parser.parse_args()

get_col = {
	'Score' : 1,
	'z' : 3,
	'frac_q' : 5,
	'frac_d' : 6,
	'coverage' : 7
}

if not args.file:
	sys.stderr.write('WARNING: Please specify the sorted results file with hosts!\n')
	sys.exit(2)

#------------------------------------------------------------------------------
#	Read predictions
#------------------------------------------------------------------------------
pred = [pred.strip().split('\t') for pred in open(args.file) 
	if not pred.startswith('Template')]
#accnname = re.search('([0-9]\.?|[A-Z]+_|[A-Z])+\d(?=\.fsa)', args.file).group(0)
#accnname = re.search('([0-9]\.?|[A-Z]+_|[A-Z])+\d(?=_)', args.file).group(0)
#------------------------------------------------------------------------------
#	frac threshold
#------------------------------------------------------------------------------
# First hit value
first_value = float(pred[0][get_col[args.decision]])

# Remove hits that belong to the same clusters of higher hits
if args.clustering == 'TRUE':
    #---------------------------------------------------------------------------
    # Cluster similar phages
    #---------------------------------------------------------------------------
    links = open('skipped_16mers_0.7.list', 'r')

    # Use the kept genome as key and append the skipped ones in the corresponding
    # value list.
    clusters={}

    for line in links:
            line = line.split(' ')
            if line[10] in clusters:
                    clusters[line[10]].append(line[3])
            else:
                    clusters[line[10]] = [line[3]]

    # Convert the dictionary into a list of lists
    merged = [flatten(subtuple) for subtuple in clusters.iteritems()]
    
    groupset = set()	
    remove = []
    for sublist in pred:
	    accn = sublist[0].strip(' ')
	    # Save accn to remove because already in a group
	    for group1 in groupset:
		    if accn in group1:
			    remove.append(accn)
			    break
	    # Save the represented groups into a set
	    group = [group for group in merged if accn in group]
	    groupset.add(tuple(flatten(group)))	
    # Consider only prediction for accn that are not in remove	
    pred = [hit for hit in pred if not hit[0].strip(' ') in remove]

# Make host => [score, dec_value] dictionary
host_set= set()
for hit in pred:
	host_set.add(hit[-1])

hosts = {k: [0, 0] for k in host_set}
for hit in pred:
	score = (float(hit[get_col[args.decision]]) / first_value)**args.alpha
	hosts[hit[-1]][0] += score
	dec_value = hit[get_col[args.decision]]
	# Take the measure value of the highest hit
	if hosts[hit[-1]][1] == 0: 
		hosts[hit[-1]][1] = dec_value

# Take host with highest value
winner = sorted(host_set, key=lambda x: (hosts[x][0], hosts[x][1]),
	 reverse=True)[0]
#for host in winner:
 # print '%s\t%s\t%s' % (args.file, host, hosts[host][1])

##print '%s\t%s\t%s' % (accnname, winner, hosts[winner][1])
input = args.file

print '%s\t%s\t%s' % (input.replace("_pred_16mers_evalue0.05_host_sorted",""), winner, hosts[winner][1])
#print '%s\t%s' % (args.file, winner)
