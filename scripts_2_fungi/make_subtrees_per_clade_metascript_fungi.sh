#!/bin/bash

prefix=$(echo /sci/labs/asafle/alexlevylab/icore-data/alex/workspace/HA_vs_NHA/data_sources/final_genomes/predictions/scripts_2_fungi)
COLLAPSE=$1
LABEL=$2
FUNGI_TREE_PATH=$3

#make genome lists of the clades we want (per sample)
mkdir genome_lists_of_clades
job_ID2=($(sbatch --wrap="python3 $prefix/get_genome_list_by_clade_sample.py $COLLAPSE $LABEL" -o slurms/slurm_genome_lists_%J.out | grep -oE "[[:digit:]]+"))

#wait
sbatch --wait --dependency=afterok:$(echo $job_ID2) --wrap="echo finished lists" -o slurms/slurm_finished_genome_lists_%J.out

mkdir tree_files_sample_1/subtrees_for_scoary_collapse_$COLLAPSE
mkdir tree_files_sample_2
cp tree_files_sample_1/tree_sample1_collapse_${COLLAPSE}_clade_with_metadata.tsv tree_files_sample_2/tree_sample2_collapse_${COLLAPSE}_clade_with_metadata.tsv # assumes identical samples
mkdir tree_files_sample_2/subtrees_for_scoary_collapse_$COLLAPSE
mkdir tree_files_sample_3
cp tree_files_sample_1/tree_sample1_collapse_${COLLAPSE}_clade_with_metadata.tsv tree_files_sample_3/tree_sample3_collapse_${COLLAPSE}_clade_with_metadata.tsv # assumes identical samples
mkdir tree_files_sample_3/subtrees_for_scoary_collapse_$COLLAPSE

#job3 R script to create trees

filename="scored_clades_collapse_${COLLAPSE}.csv"

# Use tail to skip the first row, then use process substitution to feed the output of awk into the while loop
while IFS= read -r clade_id; do

	Rscript $prefix/cut_tree_by_genome_list.R $clade_id 1 $COLLAPSE $FUNGI_TREE_PATH
	Rscript $prefix/cut_tree_by_genome_list.R $clade_id 2 $COLLAPSE $FUNGI_TREE_PATH
	Rscript $prefix/cut_tree_by_genome_list.R $clade_id 3 $COLLAPSE $FUNGI_TREE_PATH

done < <(tail -n +2 "$filename" | awk -F ',' '{print $1}')

