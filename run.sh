#!/bin/bash

# arguments
nuc_pos=$1
link_seq=$2
iter=$3
out_folder=$4

# check
if [ -z "$out_folder" ]; then
    echo "Missing arguments. Exiting..."
    exit
elif [ -d $out_folder ]; then
    echo "Output directory already exists. Exiting..."
    exit
fi

# start
ulimit -s 50000

mkdir $out_folder
mkdir $out_folder/output_pos $out_folder/output_cart $out_folder/output_nucl $out_folder/output_tables_helpar $out_folder/output_schnarp $out_folder/output_triad

# check overlap
./MC_Chromatin_olap $nuc_pos $link_seq $out_folder
cd $out_folder/output_schnarp
filename=vol_olap_000000.dat
line=$(head -n 1 $filename)

# run
if [ $line = "0" ]; then
    echo "This structure is not overlapping physically."
    cd -
    ./MC_Chromatin $nuc_pos $link_seq $iter $out_folder
else
    echo "The created structure physically overlaps. This function only works for non-overlapping starting structures. Change the input file for nucleosome positions so that the structure does not overlap. You can see the results for the created overlapping starting structure."
fi

