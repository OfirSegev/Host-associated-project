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

#conda activate /sci/labs/asafle/alexlevylab/icore-data/miniconda3/envs/scoary_env/
conda activate /sci/labs/asafle/ofirse/icore-data/miniconda3/envs/scoary_env/

COLLAPSE=$1
DISTANCE=$2
PFAM_THRESHOLD=$3
LABEL=$4
SAMPLE_1_PATH=$5
Pfam_E_value_Threshold=$6
HMM_coverage_threshold=$7
AFC_MASTER_SHEET=$8
#sample_percentage="0.8"
#num_of_bootstrap="10"
prefix=$(echo /sci/labs/asafle/alexlevylab/icore-data/alex/workspace/HA_vs_NHA/data_sources/final_genomes/predictions/scripts_3_FL)


#job id array
declare -a job_IDs

#run analysis on all the best clades with all samples - create scoary outputs
while read clade;
do

  job_IDs+=($(sbatch --wrap="bash $prefix/run_analysis_on_clade_and_sample.sh $clade 1 $COLLAPSE $PFAM_THRESHOLD $LABEL $Pfam_E_value_Threshold $HMM_coverage_threshold $AFC_MASTER_SHEET" -o slurms/slurm_run_analysis_FL_"$clade"_sample_1_%J.out --time=12:0:0 | grep -oE "[[:digit:]]+"))
  job_IDs+=($(sbatch --wrap="bash $prefix/run_analysis_on_clade_and_sample.sh $clade 2 $COLLAPSE $PFAM_THRESHOLD $LABEL $Pfam_E_value_Threshold $HMM_coverage_threshold $AFC_MASTER_SHEET" -o slurms/slurm_run_analysis_FL_"$clade"_sample_2_%J.out --time=12:0:0 | grep -oE "[[:digit:]]+"))
  job_IDs+=($(sbatch --wrap="bash $prefix/run_analysis_on_clade_and_sample.sh $clade 3 $COLLAPSE $PFAM_THRESHOLD $LABEL $Pfam_E_value_Threshold $HMM_coverage_threshold $AFC_MASTER_SHEET" -o slurms/slurm_run_analysis_FL_"$clade"_sample_3_%J.out --time=12:0:0 | grep -oE "[[:digit:]]+"))

done < list_of_best_clades.tsv

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

#move all combined files to one folder
job_ID1=($(sbatch --dependency=afterok:$(echo ${job_IDs[@]} | tr ' ' :) --wrap="bash $prefix/move_to_one_folder.sh $COLLAPSE $LABEL" -o slurms/slurm_move_to_one_folder_%J.out | grep -oE "[[:digit:]]+"))

sbatch --wait --dependency=afterok:$(echo $job_ID1) --wrap="echo move_to_one_folder is done" -o slurms/slurm_done_move_to_one_folder_%J.out

#combine 3 sample to one file for each clade
job_ID2=($(sbatch --wrap="bash $prefix/combine_3_samples.sh $COLLAPSE $LABEL" -o slurms/slurm_combine_3_samples_%J.out | grep -oE "[[:digit:]]+"))

#new fdr correction only full sampled tables (0)
job_ID3=($(sbatch --dependency=afterok:$(echo $job_ID2) --wrap="python3 $prefix/new_all_clades_fdr_correction.py $COLLAPSE $LABEL" -o slurms/slurm_new_all_clades_fdr_correction_%J.out | grep -oE "[[:digit:]]+"))

#the scripts ends when the job ends
sbatch --wait --dependency=afterok:$(echo $job_ID3) --wrap="rm -f relevant_gs__*; echo script3 FL is done" -o slurms/slurm_done_script3_%J.out

### exit code 1 if something fails
exit_code=$(seff $job_ID3 | grep State: | cut -d' ' -f 2 )
if [ $exit_code != "COMPLETED" ];
then
        exit 1
fi
