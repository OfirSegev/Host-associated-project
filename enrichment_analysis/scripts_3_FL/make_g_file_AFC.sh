#!/bin/bash

current_user=$(whoami)

if [ "$current_user" == "alexlevylab" ]; then
	source /sci/labs/asafle/alexlevylab/icore-data/miniconda3/etc/profile.d/conda.sh
	conda activate ofir_base_like
else
    if [ "$current_user" == "ofirse" ]; then
	source /sci/labs/asafle/ofirse/icore-data/mambaforge/etc/profile.d/conda.sh
	conda activate base
    else
        echo "Unknown user: $current_user"
    fi
fi


CLADE=$1 
SAMPLE=$2
COLLAPSE=$3
THRESHOLD=$4
LABEL=$5
AFC_MASTER_SHEET=$6
prefix=$(echo /sci/labs/asafle/alexlevylab/icore-data/alex/workspace/HA_vs_NHA/data_sources/final_genomes/predictions/scripts_3_FL)

traits="./files_for_scoary/clade_"$CLADE"_sample_"$SAMPLE"/t_traits__collapse_"$COLLAPSE"__clade_"$CLADE"_sample_"$SAMPLE".csv"
#traits="../HvN_collapse_G_distance_0.8_toy3/files_for_scoary/clade_"$CLADE"_sample_"$SAMPLE"/t_traits__collapse_"$COLLAPSE"__clade_"$CLADE"_sample_"$SAMPLE".csv"
echo $traits


awk -F"," ' NR>1 {print $1}' $traits > tmp_searcher__"$CLADE"__"$SAMPLE"

echo doing search
awk -F"," 'NR==FNR {search[$0]; next} ($3 in search) {print $3 "," $5}' tmp_searcher__"$CLADE"__"$SAMPLE" $AFC_MASTER_SHEET > relevant_gs__"$CLADE"__"$SAMPLE"
#awk -F"," '{print $3 "," $5}' $AFC_MASTER_SHEET | rg -F -w -f  tmp_searcher__"$CLADE"__"$SAMPLE" > relevant_gs__"$CLADE"__"$SAMPLE"
#awk -F"," '{print $3 "," $5}' /sci/labs/asafle/alexlevylab/icore-data/alex/workspace/HA_vs_NHA/data_sources/final_genomes/genome_data/new_full_protein_data/mmseqs_to_AFC_master_sheet.csv | rg -F -w -f  tmp_searcher__"$CLADE"__"$SAMPLE" > relevant_gs__"$CLADE"__"$SAMPLE"
#awk -F"," '{print $3 "," $6}' /sci/labs/asafle/alexlevylab/icore-data/alex/workspace/HA_vs_NHA/data_sources/final_genomes/genome_data/full_protein_scripts_and_data/mmseqs_to_AFC_master_sheet.csv | rg -F -w -f  tmp_searcher__"$CLADE"__"$SAMPLE" > relevant_gs__"$CLADE"__"$SAMPLE"


rm -f tmp_searcher__"$CLADE"__"$SAMPLE"

echo doing python
python3 $prefix/make_g_file_AFC_helper.py relevant_gs__"$CLADE"__"$SAMPLE" g_AFC_presence_absence__collapse_"$COLLAPSE"__clade_"$CLADE"_sample_"$SAMPLE".csv r_genomelist__collapse_"$COLLAPSE"__clade_"$CLADE"_sample_"$SAMPLE".csv $THRESHOLD g_AFC_full_counts__collapse_"$COLLAPSE"__clade_"$CLADE"_sample_"$SAMPLE".csv 
