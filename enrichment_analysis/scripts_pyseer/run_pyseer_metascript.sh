#!/bin/bash

current_user=$(whoami)

if [ "$current_user" == "alexlevylab" ]; then
	source /sci/labs/asafle/alexlevylab/icore-data/miniconda3/etc/profile.d/conda.sh
else
    if [ "$current_user" == "ofirse" ]; then
	source /sci/labs/asafle/ofirse/icore-data/mambaforge/etc/profile.d/conda.sh
    else
        echo "Unknown user: $current_user"
    fi
fi

conda activate /sci/labs/asafle/alexlevylab/icore-data/miniconda3/envs/pyseer_env/

COLLAPSE=$1
LABEL=$2
PFAM_THRESHOLD=$3

prefix=$(echo /sci/labs/asafle/alexlevylab/icore-data/alex/workspace/HA_vs_NHA/data_sources/final_genomes/predictions/scripts_pyseer)

python_path=$(echo /sci/labs/asafle/alexlevylab/icore-data/miniconda3/envs/pyseer_env/bin/python)

#make folder for pyseer
mkdir files_for_pyseer

#job id array
declare -a job_IDs

#create folder and run analysis on all the best clades with all samples
while read CLADE;
do

  #get clade size
  CLADE_SIZE=$(awk -F',' -v clade="$CLADE" '$1 == clade {print $2}' scored_clades_collapse_$COLLAPSE.csv)

  # Calculate memory in GB without the suffix
  MEM_GB=$(( (CLADE_SIZE / 130) - 13 ))

  # Ensure MEM_GB is at least 5
  if [ $MEM_GB -lt 5 ]; then
       MEM_GB=5
  fi

  # Append the G suffix to MEM_GB
  MEM_GB="${MEM_GB}G"

  #loop through 3 samples
  for SAMPLE in 1 2 3
    do
       #organize files
       mkdir files_for_pyseer/clade_"$CLADE"_sample_"$SAMPLE"

       #hard link all relevant files
       ln files_for_scoary/clade_"$CLADE"_sample_"$SAMPLE"/n_newicktree* files_for_pyseer/clade_"$CLADE"_sample_"$SAMPLE"
       sed 's/,/\t/g' files_for_scoary/clade_"$CLADE"_sample_"$SAMPLE"/t_traits* > files_for_pyseer/clade_"$CLADE"_sample_"$SAMPLE"/t_traits__collapse_"$COLLAPSE"__clade_"$CLADE"_sample_"$SAMPLE".tsv
       sed 's/,/\t/g' files_for_scoary/clade_"$CLADE"_sample_"$SAMPLE"/g_gene_presence_absence* > files_for_pyseer/clade_"$CLADE"_sample_"$SAMPLE"/g_gene_presence_absence__collapse_"$COLLAPSE"__clade_"$CLADE"_sample_"$SAMPLE".tsv
       #add line to g_gene_presence_absence file
       sed -i '1s/^/Pfam_id /' files_for_pyseer/clade_"$CLADE"_sample_"$SAMPLE"/g_gene_presence_absence__collapse_"$COLLAPSE"__clade_"$CLADE"_sample_"$SAMPLE".tsv

       cd files_for_pyseer/clade_"$CLADE"_sample_"$SAMPLE"
       mkdir outputs_and_analysis
       mkdir slurms
       if [ $CLADE_SIZE -lt 15000 ]
       then
       #run pyseer
        job_IDs+=($(sbatch --wrap="$python_path $prefix/make_files_and_run_pyseer.py $CLADE $SAMPLE $COLLAPSE $PFAM_THRESHOLD" -o slurms/slurm_make_files_and_run_pyseer_"$CLADE"_sample_"$SAMPLE"_%J.out -c3 --mem=$MEM_GB --time=4:0:0 | grep -oE "[[:digit:]]+"))
       else
        job_IDs+=($(sbatch --wrap="$python_path $prefix/make_files_and_run_pyseer.py $CLADE $SAMPLE $COLLAPSE $PFAM_THRESHOLD" -o slurms/slurm_make_files_and_run_pyseer_"$CLADE"_sample_"$SAMPLE"_%J.out -c3 --mem=500G --time=10:0:0 | grep -oE "[[:digit:]]+"))
       fi
       cd ../..
    done

done < list_of_best_clades.tsv

#wait for pyseer to finish
job_ID1=($(sbatch --wait --dependency=afterok:$(echo ${job_IDs[@]} | tr ' ' :) --wrap="echo pyseer_done" -o slurms/slurm_pyseer_done_%J.out | grep -oE "[[:digit:]]+"))

#move lmm pyseer results to one folder
mkdir lmm_results
ln files_for_pyseer/*/outputs_and_analysis/lmm*.assoc lmm_results

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

mkdir final_results

#fix no distances problem was removed from this script (can be found in backup)

#fdr for lmm 
job_ID2=($(sbatch --wrap="python3 $prefix/analyze_and_fix_output.py $COLLAPSE $LABEL lmm" -o slurms/slurm_analyze_lmm_output_%J.out --time=8:0:0 -J analyze_lmm | grep -oE "[[:digit:]]+"))

#script done
job_ID3=($(sbatch --wait --dependency=afterok:$(echo $job_ID2) --wrap="echo pyseer metascript done" -o slurms/slurm_pyseer_metascript_done_%J.out | grep -oE "[[:digit:]]+"))