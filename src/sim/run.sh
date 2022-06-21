#!/bin/bash

nuc_pos=$1
link_seq=$2
iter=$3
helpars=$4
out_folder=$5

if [ -z "$out_folder" ]; then
    echo "Missing arguments. Exiting..."
    exit
elif [ -d $out_folder ]; then
    echo "Output directory already exists. Exiting..."
    exit
fi

ulimit -s 50000
ulimit -a

mkdir $out_folder
mkdir $out_folder/output_acc $out_folder/output_cart $out_folder/output_dnaflex $out_folder/output_nucl $out_folder/output_pos $out_folder/output_schnarp $out_folder/output_tables_helpar $out_folder/output_triad

./MC_Chromatin $nuc_pos $link_seq $iter $helpars $out_folder

