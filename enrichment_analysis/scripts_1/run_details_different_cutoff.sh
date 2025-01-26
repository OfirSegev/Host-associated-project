#!/bin/bash
# activate ete3 enriroment!!!
set -e

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

prefix=$(echo /sci/labs/asafle/alexlevylab/icore-data/alex/workspace/HA_vs_NHA/data_sources/final_genomes/predictions/scripts_1)
COLLAPSE=$1
DISTANCE=$2
LABEL=$3
SAMPLE_1_PATH=$4
MIN_CLADE_SIZE=$5
MIN_ENTROPY=$6


job_ID1=($(sbatch --wrap="python3 $prefix/devide_into_clades.py $COLLAPSE $DISTANCE $LABEL $SAMPLE_1_PATH" -o slurms/slurm_devide_into_clades_"$COLLAPSE"_%J.out -J devide_into_clades | grep -oE "[[:digit:]]+"))

job_ID2=($(sbatch --wrap="python3 $prefix/clade_entropy.py $COLLAPSE $LABEL" -o slurms/slurm_clade_entropy_"$COLLAPSE"_%J.out --dependency=afterok:$(echo $job_ID1) -J clade_entropy | grep -oE "[[:digit:]]+"))

job_ID3=($(sbatch --wrap="python3 $prefix/score_the_clades.py $COLLAPSE $MIN_CLADE_SIZE $MIN_ENTROPY" -o slurms/slurm_score_the_clades_entropy_"$COLLAPSE"_%J.out --dependency=afterok:$(echo $job_ID2) -J score_clades | grep -oE "[[:digit:]]+"))

sbatch --wait --wrap="echo run_details_different_cutoff.sh done" -o slurms/slurm_done_"$COLLAPSE"_%J.out --dependency=afterok:$(echo $job_ID3)
