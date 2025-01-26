#!/bin/bash

while IFS= read -r genome_list; do

	folder_name=$(echo $genome_list | cut -d'_' -f 3,4,5,6 | cut -d'.' -f 1)
	mkdir faa_folders/$folder_name
	sample=$(echo $genome_list | cut -d'_' -f 6 | cut -d'.' -f 1)

	while read genome_id;
		do
			ln /sci/labs/asafle/alexlevylab/icore-data/alex/workspace/HA_vs_NHA/data_sources/final_genomes/genome_data/faa/$genome_id.faa faa_folders/$folder_name

		done < genome_lists_of_clades/$genome_list

done < <( ls genome_lists_of_clades )














