#!/bin/bash

##### files needed
#pfam_interpro_go_table.tsv
#go-basic.obo
#relevant combined_results_*_collapse_*.csv
#go_terms_all_levels.tsv

LABEL=$1
COLLAPSE=$2
TYPE=$3
number_of_tests_threshold=$4
folder_name=$5

#need to run map_all_genera_background_to_go.py first
#mkdir All
#python3 map_all_genera_background_to_go.py $folder_name $LABEL $COLLAPSE

mkdir $LABEL
python3 summarize_per_clade.py $LABEL $COLLAPSE $TYPE

python3 map_go_terms_per_clade.py $LABEL $COLLAPSE $TYPE $number_of_tests_threshold

mkdir results
mkdir results/all_clades
python3 hypergeometric_distribution_per_clade.py $LABEL $COLLAPSE $TYPE $number_of_tests_threshold

mkdir results/to_plot
python3 GO_terms_to_plot_per_clade.py $LABEL $COLLAPSE $TYPE $number_of_tests_threshold
