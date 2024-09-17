#!/bin/bash
#export TMPDIR=/temps/$1
#remember to activate enviroment!!
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
#conda activate /sci/labs/asafle/alexlevylab/icore-data/miniconda3/envs/phylophlan/
conda activate /sci/labs/asafle/ofirse/icore-data/anaconda_files/envs/new_phylophlan/
#@ change diversity to low/medium
#@ change min_num_of_entries to 90% of the clade size
#@ change -i to the faa folder
#@ change -f full path
#@ change -o output directory
# decide on a name for the output file (find flag)
#@ need to get sample number for -o

# inputs to the script: 1=clade_x_sample_y, 2=clade size*90%, 3=collapse

sample=$(echo $1 | cut -d'_' -f 4)
clade=$(echo $1 | cut -d'_' -f 2)

COLLAPSE=$3

job_ID1=($(sbatch -c20 --mem=160g --time=18:0:0 -o slurms/slurm_phylophlan_"$1"_%J.out -J tree_"$clade"_"$sample" --wrap="phylophlan \
-i faa_folders/$1 \
-d phylophlan \
-f /sci/labs/asafle/alexlevylab/icore-data/phylophlan_db/supermatrix_aa_custom.cfg \
-o tree_files_sample_$sample/phylophlan_output/$1 \
--nproc 128 \
--min_num_entries $2 \
--diversity medium \
--submat pfasum60 \
--trim greedy \
--remove_fragmentary_entries \
--fragmentary_threshold 0.5 \
--subsample twentyfivepercent \
--scoring_function trident \
--not_variant_threshold 0.9 \
--gap_perc_threshold 0.7 \
--db_type a" \
| grep -oE "[[:digit:]]+"))


#copy tree to subtrees folder
job_ID2=($(sbatch --wait --dependency=afterok:$(echo $job_ID1) --wrap="cp tree_files_sample_$sample/phylophlan_output/$1/$1.tre tree_files_sample_$sample/subtrees_for_scoary_collapse_$COLLAPSE/tree_'$COLLAPSE'__cladeID_$clade.newick" -o slurms/slurm_copy_tree_"$1"_%J.out | grep -oE "[[:digit:]]+"))

#wait untill the jobs are done
sbatch --wait --dependency=afterok:$(echo $job_ID2) --wrap="echo copy_tree_done" -o slurms/slurm_copy_tree_done_"$1"_%J.out

