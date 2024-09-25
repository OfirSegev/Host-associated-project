#!/bin/bash

mkdir all_clades_gff_pfams_tables
mkdir scaffold_correlation_matrixs
mkdir slurms
mkdir figures
mkdir mcl_inputs

sbatch --wrap="python3 0_genus_by_genus_scaffold_approach_new_coords_lenient.py" -c30 --time=24:0:0

./1_convert_corr_loop.sh

mkdir results_with_clusters_ids
#run 2 or 3
#./2_run_mcl_command.sh
python3 3_parse_mcl_into_cluster_ids_v2.py

python3 4_combine_results_new.py

#add informative columns
python3 5_further_operon_analysis.py

#find full operons in img to investigate
