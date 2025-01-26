#!/bin/bash

mkdir scoary_results
mkdir scoary_results/raw

COLLAPSE=$1
LABEL=$2

while IFS= read -r path_name; do

	folder_name=$(echo $path_name | cut -d'/' -f 2)

	ln $path_name/outputs_and_analysis/all_csvs/"$LABEL"_collapse_"$COLLAPSE"__"$folder_name"_combined.csv scoary_results/raw

done < <( ls -d files_for_scoary/* )
