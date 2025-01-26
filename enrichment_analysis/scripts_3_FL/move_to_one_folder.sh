#!/bin/bash

mkdir scoary_results_FL
mkdir scoary_results_FL/raw

COLLAPSE=$1
LABEL=$2

while IFS= read -r path_name; do

	folder_name=$(echo $path_name | cut -d'/' -f 2)

	ln $path_name/outputs_and_analysis_FL/all_csvs/"$LABEL"_collapse_"$COLLAPSE"__"$folder_name"_combined.csv scoary_results_FL/raw

done < <( ls -d files_for_scoary/* )
