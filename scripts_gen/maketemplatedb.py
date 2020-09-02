#!/usr/bin/env python

# Copyright (c) 2014, Ole Lund, Technical University of Denmark
# All rights reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#
# Import libraries
#
import sys
import time
import os
import argparse
from operator import itemgetter
from itertools import groupby
import cPickle as pickle


#################################################################
# FUNCTIONS:
#################################################################


def reversecomplement(seq):
    ''' Reverse complement '''
    comp = ''
    for s in seq:
        if s == 'A':
            comp = comp + 'T'
        elif s == 'T':
            comp = comp + 'A'
        elif s == 'C':
            comp = comp + 'G'
        elif s == 'G':
            comp = comp + 'C'
        else:
            comp = comp + s
    return comp[::-1]


def fasta_iter(fasta_name):
    ''' Given a fasta file. yield tuples of header, sequence '''
    fh = open(fasta_name)
    
    # Create a FASTA ITERATOR
    # groupby groups together items that have the same key
    # --> lines that have the same header
    # x[0] is a boolean (True for lines that start with '>' and False for other)
    # x[1] are the groups (header and sequence)
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))

    for header in faiter:
        # drop the ">"
        header = header.next()[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip().upper() for s in faiter.next())
        yield header, seq


def window(seq, kmersize, stepsize):
#        kmers_count = ((len(sequence) - kmersize) / stepsize) + 1
    for i in range(0, len(seq)-kmersize +1, stepsize):
            yield seq[i : i + kmersize]


def check_homology(inputseq, inputs):
    ''' Check homology between input sequence and stored sequence '''
    # Make list of unique k-mers in entry
    queryindex = {}
    uquerymers = 0
    seqlen = len(inputseq)
    for qseq in[inputseq, reversecomplement(inputseq)]:
        for submer in window(qseq, kmersize, 1):
            if submer.startswith(prefix):
                if submer in queryindex:
                    queryindex[submer] += 1
                else:
                    queryindex[submer] = 1
                    uquerymers += 1

    # Search for matches:
    mincoverage = 1
    templateentries = {}

    for submer in queryindex:
        if submer in inputs:
            if queryindex[submer] >= mincoverage:
                matches = inputs[submer].split(",")
                # get unique list:
                umatches = list(set(matches))
                for match in umatches:
                    # count matches:
                    if match in templateentries:
                        templateentries[match] += 1
                    else:
                        templateentries[match] = 1

    # calculate frac_q:
    frac_q = 0.0
    hitname = ""
    score = 0
    sortedlist = sorted(
        templateentries.items(), key=itemgetter(1), reverse=True)
    for template, score in sortedlist:
        # multiplication by 2 because of reverse complement
        frac_q = score / (float(uquerymers) + etta)
        hitname = template
        break

    del templateentries
    del queryindex

    return (hitname, frac_q, score)


def update_database(inputseq, inputname, filters):
    ''' Update the database with the new sequence '''
    # define global variables:
    global inputs, lengths, ulengths, descriptions, desc, Nstored, Nstored_old, Nustored, Nustored_old
    global kmer_count, kmersize, prefixlen, stepsize, t0, t1, printfreq
    global filterfilename, homthres

    # Start of database update
    for seq in[inputseq, reversecomplement(inputseq)]:

        expected_kmers = divmod((((len(seq) - kmersize) / stepsize) + 1),1)[0]
        kmer_count += expected_kmers

        for submer in window(seq, kmersize, stepsize):
            if submer.startswith(prefix):
                if (filterfilename != None and submer not in filters) or filterfilename == None:
                    Nstored += 1
                    if submer in inputs:
                        if (inputs[submer].find(inputname) == -1):
                            Nustored += 1
                        inputs[submer] = inputs[submer] + "," + inputname
                    else:
                        inputs[submer] = inputname
                        Nustored += 1

    # update nr of kmers:
    lengths[inputname] = Nstored - Nstored_old
    # update nr of unique kmers:
    ulengths[inputname] = Nustored - Nustored_old
    # update descriptions:
    descriptions[inputname] = desc
    # print inputname,lengths[inputname], Nstored, Nstored_old,"i: ",i,
    # len(inputseq)
    Nstored_old = Nstored
    Nustored_old = Nustored


def process_entry(inputname):
    ''' Decide weather to include or exclude entry '''
    global inputseq, inputs, homthres, filters

    sys.stdout.write("%s %s\n" % ("# Entry read", inputname))

    # check homology:
    if homthres != None:
        sys.stdout.write("%s\n" % ("# Checking for homology"))

        (hitname, frac_q, score) = check_homology(inputseq, inputs)
        sys.stdout.write("# Max frac_q similarity of %s to %s frac_q: %s Score: %s\n" % (
            inputname, hitname, frac_q, score))

    # decide to exclude / include entry:
    if (homthres != None and frac_q >= homthres):
        # exclude entry:
        sys.stdout.write("# Skipping entry: %s in database due to similarity to %s frac_q: %s\n" % (
            inputname, hitname, frac_q))
    if (homthres == None or (homthres != None and frac_q < homthres)):
        # include entry -> update database:
        sys.stdout.write("%s %s\n" % ("# Including entry: ", inputname))
        update_database(inputseq, inputname, filters)

##########################################################################
# DEFINE GLOBAL VARIABLES:
##########################################################################

global inputs, lengths, ulengths, descriptions, desc, Nstored, Nstored_old, Nustored, Nustored_old
global kmer_count, kmersize, prefixlen, stepsize, t0, t1, printfreq
global filterfilename, homthres, filters, organismlist


# Parse command line options
#
#
parser = argparse.ArgumentParser(description='Make kmer database')
parser.add_argument('-i', '--inputfile', dest="inputfilename",help="read from INFILE", metavar="INFILE")
parser.add_argument("-l", "--inputfilelist", dest="inputfilelist",help="read a list of fatsa file locations", metavar="INFILELIST")
parser.add_argument("-f", "--filterfile", dest="filterfilename",help="filter (ignore) K-mers present in FILTERFILE", metavar="FILTERFILE")
parser.add_argument("-o", "--outputfile", dest="outputfilename",help="write to OUTFILE", metavar="OUTFILE")
parser.add_argument("-k", "--kmersize", dest="kmersize", help="Size of KMER", metavar="KMERSIZE")
parser.add_argument("-t", "--homthres", dest="homthres",help="Threshold for homology reduction", metavar="HOMTHRES")
parser.add_argument("-s", "--stepsize", dest="stepsize",help="Size of step between K-mers", metavar="STEPSIZE")
parser.add_argument("-x", "--prefix", dest="prefix", help="type of prefix", metavar="PREFIX")
parser.add_argument("-a", "--templatefile", dest="templatefilename",help="add to database TEMFILE", metavar="TEMFILE")
parser.add_argument("-c", "--organismlistname", dest="organismlistname",help="provide organism list to replace IDs ORGLIST", metavar="ORGLIST")
args = parser.parse_args()


# Open file for input sequence with kmers to save in database:
if args.inputfilename != None:
    if args.inputfilename == "--":
        inputfile = sys.stdin
    else:
        #inputfile = open(args.inputfilename,"r")
        inputfile = args.inputfilename

# Open list of FASTA file locations:
if args.inputfilelist != None:

    if args.inputfilelist == "--":
        inputfilelist = sys.stdin
    else:
        inputfilelist = args.inputfilelist
    #    inputfilelist = open(args.inputfilelist, "r")
elif inputfile != "":
    inputfilelist = [inputfile]

# Open file to filter on (kmers not to save in database):
if args.filterfilename != None:
    filterfilename = args.filterfilename
else:
    filterfilename = None

# Open templatefile to add to already existing databse:
if args.templatefilename != None:
    templatefile = open(args.templatefilename + ".p", "rb")
    templatefile_lengths = open(args.templatefilename + ".len.p", "rb")
    templatefile_descriptions = open(args.templatefilename + ".desc.p", "rb")
    try:
        templatefile_ulengths = open(
            args.templatefilename + ".ulen.p", "rb")
    except:
        # do nothing
        pass

# Harcode this to be true so I do not need to use the -p option

if args.outputfilename != None:
    outputfile = open(args.outputfilename+".p", "wb")
    outputfile_lengths = open(args.outputfilename+".len.p", "wb")
    outputfile_ulengths = open(args.outputfilename+".ulen.p", "wb")
    outputfile_descriptions = open(args.outputfilename+".desc.p", "wb")
else:
    outputfile = sys.stdout


# Size of K-mer
if args.kmersize != None:
    kmersize = int(args.kmersize)
else:
    kmersize = 16


# Size of step when looking for K-mers in the sequence
if args.stepsize != None:
    stepsize = int(args.stepsize)
else:
    stepsize = 1

# Homology threshold for when to include entry in database
if args.homthres != None:
    homthres = float(args.homthres)
else:
    homthres = None

# Prefix to use fro filtering sequences
if args.prefix != None:
    prefix = args.prefix
    prefixlist = [prefix]
    prefixlen = len(prefixlist[0])
else:
    prefix = ''
    prefixlist = [prefix]
    prefixlen = len(prefixlist[0])


# get organism list:
if args.organismlistname != None:
    organismlist = open(args.organismlistname, "r")
else:
    organismlist = None

##################################################################
# INITIALIZE STATISTICS
##################################################################

# # of kmers
kmer_count = 0
# Start time to keep track of progress
t0 = time.time()
# Print progress
printfreq = 100000
# frequenct to save sorted list to db
dbsavefreq = 30000000

etta = 0.0001

###################################################################
# READ SEQUENCES FROM FILTERFILE AND SAVE KMERS
###################################################################

filters = {}
#Nfilters=0 # Unused variable
t1 = time.time()

if args.filterfilename != None:
    sys.stdout.write("%s\n" % ("# Reading filterfile"))
    fasta_generator = fasta_iter(args.filterfilename)

    for head_seq_tuple in fasta_generator:
        header, filterseq= head_seq_tuple
        header = header.split()
        filtername = header[0]

    for seq in [filterseq,reversecomplement(filterseq)]:

        kmer_count += ((len(seq) - kmersize) / stepsize) + 1

        for submer in window(seq, kmersize, kmerstep):
            if submer.startswith(prefix):
                if not submer in filters:
                    filters[submer] = filtername

            if kmer_count % printfreq == 0:
                t1 = time.time()
                sys.stdout.write("\r%s kmers (%s kmers / s)"
                 % ("{:,}".format(kmer_count),
                 "{:,}".format(kmer_count / (t1 - t0))))
                sys.stdout.flush()
            start += stepsize
    #
    # Print final statistics for filterfile
    #
    t1 = time.time()
    sys.stdout.write("\r%s kmers (%s kmers / s)"
     % ("{:,}".format(kmer_count), "{:,}".format(kmer_count / (t1-t0))))
    sys.stdout.flush()
    #sys.stdout.write("\n")


##########################################################################
# READ TEMPLATEFILE
##########################################################################

inputs = {} # dictionary kmer => templates having that kmer
lengths = {} # dictionary template_id => number of kmers
ulengths = {} # dictionary template_id => number of unique kmers
descriptions = {} # dictionary template_id => description

Nstored = 0
Nstored_old = Nstored
# count only one per sequence:
Nustored = 0
Nustored_old = Nustored

if args.templatefilename != None:
    sys.stdout.write("%s\n" % ("# Reading database of templates"))
    inputs = pickle.load(templatefile)
    lengths = pickle.load(templatefile_lengths)
    try:
        ulengths = pickle.load(templatefile_ulengths)
    except:
        sys.stderr.write('No ulen.p file found for database')
        SystemExit()
    descriptions = pickle.load(templatefile_descriptions)

    # Count number of k-mers and number of unique k-mers:
    for name in lengths:
        Nstored += 1
        if name in ulengths:
            Nustored += 1
    Nstored_old = Nstored
    Nustored_old = Nustored


##########################################################################
# READ ORGANISM LIST
##########################################################################


if args.organismlistname != None:
    sys.stdout.write("%s\n" % ("# Reading organism list"))

    organism = {}

    for l in organismlist:
        l = l.strip()
        l = l.split("\t")
        organism[l[0]] = l[1]


##########################################################################
#       PROCESS A LIST OF FASTA FILE LOCATIONS
##########################################################################

t1 = time.time()

sys.stdout.write("%s\n" % ("# Reading inputfile(s)"))

kmer_count = 0
desc = ''

for l in inputfilelist:   
    fasta_generator=fasta_iter(l)
    for head_seq_tuple in fasta_generator:
        header, inputseq = head_seq_tuple
        header = header.split()

        # get input name and translate if necessary:
        # and store description
        inputname = header[0]
        if organismlist != None:
            inputname = organism[inputname]
            desc = inputname
        else:
            desc = ' '.join(header[1:])

        # process new entry (homology check and include in
        # database):
        process_entry(inputname)


############################################################
# PRINT DATABASE OF KMERS
############################################################

pickle.dump(inputs, outputfile,2)
pickle.dump(lengths, outputfile_lengths,2)
pickle.dump(ulengths, outputfile_ulengths,2)
pickle.dump(descriptions, outputfile_descriptions,2)

###########################################################
# PRINT FINAL STATISTICS
###########################################################
t2 = time.time()
sys.stdout.write("\r%s kmers (%s kmers / s)"
 % ("{:,}".format(kmer_count), "{:,}".format(kmer_count / (t2 - t1))))
sys.stdout.flush()
sys.stdout.write("\n")
sys.stdout.write("# Total time used: %s s\n" % (t2 - t0))
#
# Done
#
