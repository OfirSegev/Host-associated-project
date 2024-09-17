#!/bin/bash

#remember pyseer enviroment!!
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

prefix=$(echo /sci/labs/asafle/alexlevylab/icore-data/alex/workspace/HA_vs_NHA/data_sources/final_genomes/predictions/scripts_pyseer_FL)

python_path=$(echo /sci/labs/asafle/alexlevylab/icore-data/miniconda3/envs/pyseer_env/bin/python)


#job id array
declare -a job_IDs

#create folder and run analysis on all the best clades with all samples
if [ ! -d files_for_pyseer ]; then
    mkdir files_for_pyseer
fi


while read CLADE;
do

  #get clade size
  CLADE_SIZE=$(awk -F',' -v clade="$CLADE" '$1 == clade {print $2}' scored_clades_collapse_$COLLAPSE.csv)

  # Calculate memory in GB without the suffix
  MEM_GB=$(( (CLADE_SIZE / 300) - 13 ))

  # Ensure MEM_GB is at least 5
  if [ $MEM_GB -lt 5 ]; then
       MEM_GB=5
  fi

  # Append the G suffix to MEM_GB
  MEM_GB="${MEM_GB}G"


  #loop through 3 samples
  for SAMPLE in 1 2 3
    do
       if [ ! -d files_for_pyseer/clade_"$CLADE"_sample_"$SAMPLE" ]; then
            mkdir files_for_pyseer/clade_"$CLADE"_sample_"$SAMPLE" 
       fi

       #get g_AFC_presence_absence from scoary folder (trait file doesn't neeed to be changed)
       sed 's/,/\t/g' files_for_scoary/clade_"$CLADE"_sample_"$SAMPLE"/g_AFC_presence_absence* > files_for_pyseer/clade_"$CLADE"_sample_"$SAMPLE"/g_AFC_presence_absence__collapse_"$COLLAPSE"__clade_"$CLADE"_sample_"$SAMPLE".tsv
       #add line to g_AFC_presence_absence file
       sed -i '1s/^/Pfam_id /' files_for_pyseer/clade_"$CLADE"_sample_"$SAMPLE"/g_AFC_presence_absence__collapse_"$COLLAPSE"__clade_"$CLADE"_sample_"$SAMPLE".tsv

       cd files_for_pyseer/clade_"$CLADE"_sample_"$SAMPLE"
       mkdir outputs_and_analysis_FL


        if [ ! -d slurms ]; then
            mkdir slurms 
        fi

       if [ $CLADE_SIZE -lt 15000 ]
       then
       #run pyseer
            job_IDs+=($(sbatch --wrap="$python_path $prefix/make_files_and_run_pyseer.py $CLADE $SAMPLE $COLLAPSE $PFAM_THRESHOLD" -o slurms/slurm_make_files_and_run_pyseer_FL_"$CLADE"_sample_"$SAMPLE"_%J.out -c7 --mem=$MEM_GB --time=3:0:0 | grep -oE "[[:digit:]]+"))
       else
            job_IDs+=($(sbatch --wrap="$python_path $prefix/make_files_and_run_pyseer.py $CLADE $SAMPLE $COLLAPSE $PFAM_THRESHOLD" -o slurms/slurm_make_files_and_run_pyseer_FL_"$CLADE"_sample_"$SAMPLE"_%J.out -c7 --mem=300G --time=10:0:0 | grep -oE "[[:digit:]]+"))
        fi
       cd ../..
    done

done < list_of_best_clades.tsv

#wait for pyseer to finish
job_ID1=($(sbatch --wait --dependency=afterok:$(echo ${job_IDs[@]} | tr ' ' :) --wrap="echo pyseer_done_FL" -o slurms/slurm_pyseer_FL_done_%J.out | grep -oE "[[:digit:]]+"))

#move lmm pyseer results to one folder
mkdir lmm_results_FL
ln files_for_pyseer/*/outputs_and_analysis_FL/lmm*.assoc lmm_results_FL

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
#conda activate base

if [ ! -d final_results_FL ]; then
        mkdir final_results_FL
fi

#fix no distances problem was removed (can be found in backup)

#fdr for lmm 
job_ID2=($(sbatch --wrap="python3 $prefix/analyze_and_fix_output.py $COLLAPSE $LABEL lmm" -o slurms/slurm_analyze_lmm_output_FL_%J.out --time=8:0:0 -J analyze_lmm | grep -oE "[[:digit:]]+"))

#script done
job_ID3=($(sbatch --wait --dependency=afterok:$(echo $job_ID2) --wrap="echo pyseer FL metascript done" -o slurms/slurm_pyseer_metascript_done_FL_%J.out | grep -oE "[[:digit:]]+"))
