#!/bin/bash

# Author: Julia Villarroel

evalue=0.05
decision=coverage
outdir="HostPhinder/results"
kmersize=16

while getopts e:i:o:d:b:t:k: opt
do
   case "$opt" in
      e) evalue=$OPTARG;;
      i) input=$OPTARG;;
      o) outdir=$OPTARG;; # output folder
      d) decision=$OPTARG;;
      b) database=$OPTARG;;  
      t) taxonomy=$OPTARG;;
      k) kmersize=$OPTARG;;
   esac
done

#create: and define hash to get column after which should be sorted
declare -A meas2col=(
    ["Score"]="2"
    ["z"]="4"
    ["frac_q"]="6"
    ["frac_d"]="7"
    ["coverage"]="8"
)

col_dec=${meas2col["$decision"]}
workdir="/home/projects/pr_phage/people/juliav"
outdir="$workdir/$outdir"
#mkdir -p $outdir
#output_file="${outdir}_${kmersize}_${taxonomy}_first_${decision}_evalue$evalue"


if [ $taxonomy = species ]
then
    meta="$workdir/meta_1871_species_150506.tab"    
    database="$workdir/HostPhinder/database/accnwhost1871.list_step1_kmer16_thres1.0_db"
elif [ $taxonomy = genus ]
then
    meta="$workdir/meta_2196_genus_150506.tab"
    database="$workdir/HostPhinder/database/accnwhost2196.list_step1_kmer16_thres1.0_db"
fi
#------------------------ Find templates in database with match ---------------
findtempl_out="$outdir/${input##*/}_${taxonomy}_pred_${kmersize}mers_evalue$evalue"

python $workdir/HostPhinder/scripts/findtemplate_scipy.py -t $database -p\
 -k $kmersize -o $findtempl_out -i $input -e $evalue

# If a prediction has been made --> more than one line (header +..)
if [[ $(wc -l <$findtempl_out) -ge 2 ]]
then
#------------------------- Description --> Host column ------------------------
    # Replace description column with host column on the significant hits table
    predWhost=${findtempl_out}_host 
    rm $predWhost
    IFS=$'\n'       # make newlines the only separator
    for line in $(cat $findtempl_out)
    do
        if [[ $line == Template* ]]; then
            continue
        fi
        template=`echo "$line" | awk -F "\t" '{print $1}' | perl -pe 's/\s+$//'`
        pred=`grep $template $meta | awk -F "\t" '{print $NF}'` 
        echo -e "$line$pred" >> $predWhost
    done

#------------------------------ Sort and get results --------------------------
# Sort according to the column the user chooses
    sort_output=${predWhost}_sorted
    sort -rgk $col_dec $predWhost > $sort_output
    rm $predWhost

#------------------------------ Alpha method (=6) ------------------------------
    alpha_out=`$workdir/HostPhinder/scripts/get_host_alpha.py -f $sort_output -d coverage -a 6.0 -c FALSE`
    echo $alpha_out
else
    echo -e "$input\tNo_significant_match_found"

fi
