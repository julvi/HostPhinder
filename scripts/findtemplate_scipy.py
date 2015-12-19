#!/usr/bin/env python

# On cge-s2 run with /tools/opt/anaconda/bin/python

#
# Import libraries
#
import argparse
import sys, time
import os
from math import sqrt, pow
import scipy
from scipy.stats import norm
from operator import itemgetter
import re
import cPickle as pickle

#
# Functions
#
# Construct the reverse complemant from a sequence
#
#
def reversecomplement(seq):
    '''Reverse complement'''
    comp = ''
    for s in seq:
        if s == 'A': comp = comp + 'T'
        elif s == 'T': comp = comp + 'A'
        elif s == 'C': comp = comp + 'G'
        elif s == 'G': comp = comp + 'C'
        else: comp = comp + s
    return comp[::-1]

#------------------------------------------------------------------------------
# Parse command line options
#------------------------------------------------------------------------------
parser= argparse.ArgumentParser(description='Look for overlapping k-mers\
 between query and templates in the database.')

parser.add_argument("-i", "--inputfile", type=str, help="Query fasta filename")
parser.add_argument("-t", "--templatefile", type=str, help="Database name\
--> makedatabase.py output")

parser.add_argument("-o", "--outputfile", type=str, help="Output filename")
parser.add_argument("-k", "--kmersize", type=int, default=16,
 help="Size of KMER")
parser.add_argument("-x", "--prefix", type=str, default='', help="prefix")

# The following option is not used as far as I can see
#parser.add_argument("-a", "--printall", help="Print matches to all templates\
# in templatefile unsorted") 
#parser.add_option("-a", "--printall", dest="printall",action="store_true", 
#help="Print matches to all templates in templatefile unsorted") 

parser.add_argument("-p", "--pickleinput", action="store_true",
 help="use pickle input on by default, option kept for backwords compatibility") 

# The first hit is the one to which most of query kmers mathch to
# with this option (not currently working) the matching k-mers are then eliminated
# from the search
parser.add_argument("-w", "--winnertakesitall", dest="wta", action="store_true",
 help="kmer hits are only assigned to most similar template")

parser.add_argument("-e", "--evalue", type=float, default=0.05,
 help="Maximum E-value")

args = parser.parse_args()
#
# set up prefix filtering
#
prefix = args.prefix
prefixlen = len(prefix)
#
#
#
evalue = args.evalue
# 
# Open files 
# 
t0 = time.time()
if args.inputfile != None:
  if args.inputfile == "--":
    inputfile = sys.stdin
  else:
    inputfile = open(args.inputfile,"r")
#
#args.pickleinput = True
if args.templatefile != None:
  if args.pickleinput == True:
    templatefile = open(args.templatefile+".p", "rb" )
    templatefile_lengths = open(args.templatefile+".len.p", "rb" )
    #try:
    templatefile_ulengths = open(args.templatefile+".ulen.p", "rb" )
    #except:
      # do nothing
     # pass
    templatefile_descriptions = open(args.templatefile+".desc.p", "rb" )
#  else:
#    templatefile = open(args.templatefile,"r")
else:
  sys.exit("No template file specified")
#
if args.outputfile != None:
  outputfile = open(args.outputfile,"w")
else:
  outputfile = sys.stdout

#
# Size of K-mer
#
kmersize = args.kmersize
#
# Make database of nmers
#
templates = {}
templateentries = {}
templateentries_tot = {}
templatebackground = {}
templatebackgroundtot = 0
Ntemplates =0
#seq = ""
#name= "None"
# This has to be made an option
oligolen=kmersize
#
# Read Template file
#
sys.stdout.write("%s\n" % ("# Reading database of templates\
 and ckecking for database integrity"))
#
# Harcode this to be true so I do not need to use the -p option
#
if args.pickleinput == True:
  templates = pickle.load(templatefile)
  templates_lengths = pickle.load(templatefile_lengths)
  try:
    templates_ulengths = pickle.load(templatefile_ulengths)
  except:
    	sys.stderr.write('No ulen.p file found for database')
	SystemExit() 
  templates_descriptions = pickle.load(templatefile_descriptions)

#
# Check for database integrity --> all template must have > 0 unique k-mers
#
if any(ulength == 0 for ulength in templates_ulengths.values()):
  sys.stderr.write('ERROR: database did not pass integrity check!\n')
  sys.exit(2)

#
#
#
template_tot_len = 0
template_tot_ulen = 0
Ntemplates = 0
#length added
#print templates_ulengths
for name in templates_lengths:
  template_tot_len  += templates_lengths[name]
  template_tot_ulen += templates_ulengths[name]
  Ntemplates += 1

#print template_tot_len
#
# Read inputfile
#
queryseq = []
queryseqsegments = []
consensusseq = []
queryname = []
querydesc = []
Nquerys=0
i=0
if args.inputfile != None:
  sys.stdout.write("%s\n" % ("# Reading inputfile"))
  for line in inputfile:
    fields=line.split()
    if len(line)>1:
      if fields[0][0] == ">":
	if (i>0):
	  queryseq[-1] = ''.join(queryseqsegments)
        del queryseqsegments
        queryseqsegments = []
	i=0
        queryseq.append("")
        consensusseq.append("")
        queryname.append(fields[0][1:]) 
        querydesc.append(re.sub(r"^[^\s]+\s","",line.strip()))
      elif fields[0][0] == "@":
        # Fasq file
	if (i>0):
	  queryseq[-1] = ''.join(queryseqsegments)
        del queryseqsegments
        queryseqsegments = []
	i=0
        queryseq.append("")
        consensusseq.append("")
        queryname.append(fields[0][1:]) 
        querydesc.append(re.sub(r"^[^\s]+\s","",line.strip()))
	try:
          line = inputfile.next()
          fields=line.split()
          queryseqsegments.append("")
          queryseqsegments[i] = fields[0]
	  i+=1
          line = inputfile.next()
          line = inputfile.next()
	except:
	  break
      else:
        queryseqsegments.append("")
        queryseqsegments[i] = fields[0]
	i+=1
  queryseq[-1] = ''.join(queryseqsegments)
del queryseqsegments
#
# Make list of unique k-mers in inputfile
#
sys.stdout.write("%s\n" % ("# Searching for matches in template file"))
queryindex = {}
qtotlen=0
querymers=0
uquerymers=0
for i in range(0, len(queryseq)):
  seqlen = len(queryseq[i])
  qtotlen += seqlen;
  #for qseq in [queryseq[i],reversecomplement(queryseq[i])]:
  for qseq in [queryseq[i]]:
    for j in range(0, seqlen-oligolen+1):
      submer = qseq[j:j+oligolen]
      if prefix == qseq[j:j+prefixlen]:
        if submer in queryindex:
          queryindex[submer] += 1
          querymers += 1 
        else:
          queryindex[submer] = 1
          querymers += 1 
          uquerymers += 1 
#
# Search for matches
#
sys.stdout.write("%s\n" % ("# Searching for matches of input in template"))
mincoverage = 1
Nhits=0
for submer in queryindex:
  if submer in templates:
    if queryindex[submer] >= mincoverage:
      matches = templates[submer].split(",")
      for match in matches:
        Nhits += 1
        if match in templateentries:
          templateentries[match] += 1
        else:
          templateentries[match] = 1
      umatches = list(set(matches))
      for match in umatches:
        if match in templateentries_tot:
          templateentries_tot[match] += queryindex[submer]
        else:
          templateentries_tot[match] = queryindex[submer]

#
# Print best scoring entries sorted
#
minscore = 0
etta = 0.001
sys.stdout.write("%s\n" % ("# Search statistics"))
sys.stdout.write("%s\n" % ("# Total number of hits: %s") % (Nhits))
sys.stdout.write("%s\n" % ("# Total number of kmers in templates : %s")
 % (template_tot_len))
#sys.stdout.write("%s\n" % ("# Minimum number of k-mer hits to report template: %s") % (minscore))
sys.stdout.write("%s\n" % ("# Maximum multiple testing corrected E-value\
 to report match : %s") % (evalue))
#sys.stdout.write("%s\n" % ("# Printing best matches"))
outputfile.write("Template\tScore\tExpected\tz\tp\tfrac_q\tfrac_d\tcoverage\
\tunique_Kmers_in_template\tunique_Kmers_in_query\tDescription\n")

def calc_values(Nhits, template_ulength, template_tot_ulen, score,
 etta, Ntemplates, uquerymers, templateentries, template_length):
  """Calculate expected, z, p, p_corr, frac_q, frac_d, coverage"""
  values = []
  expected = float(Nhits)*float(template_ulength)/float(template_tot_ulen)
  z = (score - expected)/sqrt(score + expected+etta)
  # convert Z-score to twosided p-value
  p = scipy.stats.norm.sf(z)*2
  p_corr = p*Ntemplates
  # frac_q = templateentries_tot[template]/(float(querymers)+etta)
  frac_q = score/float(uquerymers)+etta
  frac_d = score/(template_ulength+etta)
  # print score, querymers, uquerymers,templates_ulengths[template], templateentries_tot[template]
  coverage = 2*templateentries/float(template_length)
  values.extend((expected, z, p, p_corr, frac_q, frac_d, coverage))
  return values

if args.pickleinput == True:
  sortedlist= sorted(templateentries.items(), key = itemgetter(1), reverse=True)
  if not args.wta == True:   
    for template,score in sortedlist:
      if score > minscore:
        (expected, z, p, p_corr, frac_q, frac_d, coverage) = calc_values(Nhits,
         templates_ulengths[template], template_tot_ulen, score, etta,
         Ntemplates, uquerymers, templateentries_tot[template],
         templates_lengths[template])
	if p_corr <= evalue:
       		outputfile.write("%-12s\t%8s\t%8.3f\t%8.3f\t%4.1e\t%4.1e\t%4.1e\t%4.1e\t%8d\t%8d\t%s\n" % 
            		(template, score, expected, z, p_corr, frac_q, frac_d, coverage, templates_ulengths[template], uquerymers, templates_descriptions[template].strip()))
  else:
    templateentries2 = {}
    templateentries_tot2 = {}
    Nhits2=0
    for submer in queryindex:
      if submer in templates:
        matches = templates[submer].split(",")
	for match,score in sortedlist:
	  #print match,score	
          if match in matches:
	    #print "gotcha"
            Nhits2 += 1
            if match in templateentries2:
              templateentries2[match] += 1
	      templateentries_tot2[match] += queryindex[submer]
            else:
              templateentries2[match] = 1
	      templateentries_tot2[match] = queryindex[submer]
	    break
    sortedlist2= sorted(templateentries2.items(), key = itemgetter(1), reverse=True)
    for template,score in sortedlist2:
      #outputfile.write("%s %s\n" % (template,score))
      if score > minscore:
        (expected, z, p, p_corr, frac_q, frac_d, coverage) = calc_values(Nhits, 
         templates_ulengths[template], template_tot_ulen, score, etta, 
         Ntemplates, uquerymers, templateentries_tot[template], 
         templates_lengths[template])
	if p_corr <= evalue:
          outputfile.write("%-12s\t%8s\t%8.3f\t%8.3f\t%4.1e\t%4.1e\t%4.1e\t%4.1e\t%8d\t%8d\t%s\n" % 
            (template, score, expected, z, p_corr, frac_q, frac_d, coverage, templates_ulengths[template], uquerymers, templates_descriptions[template].strip()))
#
# Close files
#
t1 = time.time()
sys.stdout.write("# %s kmers (%s kmers / s). Total time used: %s sec" % ("{:,}".format(querymers), "{:,}".format(querymers / (t1-t0)),int(t1-t0)))
sys.stdout.flush()
sys.stdout.write("\n")
#sys.stderr.write("Done")
sys.stdout.write("# Closing files\n")
