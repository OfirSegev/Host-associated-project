#!/bin/bash

prefix=$(echo /sci/labs/asafle/alexlevylab/icore-data/alex/workspace/HA_vs_NHA/data_sources/final_genomes/predictions/scripts_2)
COLLAPSE=$1
DISTANCE=$2
LABEL=$3
SAMPLE_1_PATH=$4

#create samples combined csv
job_ID1=($(sbatch --wrap="python3 $prefix/fuse_in_data.py $COLLAPSE $DISTANCE $LABEL $SAMPLE_1_PATH" -o slurms/slurm_fuse_in_data_%J.out | grep -oE "[[:digit:]]+"))

#make genome lists of the clades we want (per sample)
mkdir genome_lists_of_clades
job_ID2=($(sbatch --dependency=afterok:$(echo $job_ID1) --wrap="python3 $prefix/get_genome_list_by_clade_sample.py $COLLAPSE $LABEL" -o slurms/slurm_genome_lists_%J.out | grep -oE "[[:digit:]]+"))

#make a folder of faa for each list
mkdir faa_folders
job_ID3=($(sbatch --dependency=afterok:$(echo $job_ID2) --wrap="bash $prefix/make_faa_folder_for_clade.sh" -o slurms/slurm_faa_folders_%J.out -J make_faa | grep -oE "[[:digit:]]+"))

#wait
sbatch --wait --dependency=afterok:$(echo $job_ID3) --wrap="echo finished folders" -o slurms/slurm_finished_folders_%J.out

#make folders for phylophlan
### add phylophlan_output + subtrees_for_scoary for sample 1
mkdir tree_files_sample_1/phylophlan_output
mkdir tree_files_sample_1/subtrees_for_scoary_collapse_$COLLAPSE
mkdir tree_files_sample_2
mkdir tree_files_sample_2/phylophlan_output
mkdir tree_files_sample_2/subtrees_for_scoary_collapse_$COLLAPSE
mkdir tree_files_sample_3
mkdir tree_files_sample_3/phylophlan_output
mkdir tree_files_sample_3/subtrees_for_scoary_collapse_$COLLAPSE

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
conda activate /sci/labs/asafle/alexlevylab/icore-data/miniconda3/envs/phylophlan/

cp -r ../phylophlan_databases/ .

i=0

#run phlophlan on each clade_sample
while IFS= read -r clade_sample_list; do

	clade_size=$(ls $clade_sample_list | wc -l)
	num_entries=$(echo $clade_size*0.9 | bc | awk '{print int($1)}')
	clade_sample_phrase=$(echo $clade_sample_list | cut -d'/' -f 2)

	if [ $clade_size < 14000 ]; then
		bash $prefix/phylophlan_command.sh $clade_sample_phrase $num_entries $COLLAPSE &
	else
		bash $prefix/phylophlan_command_heavy_clades.sh $clade_sample_phrase $num_entries $COLLAPSE &
	fi
	pids[${i}]=$!
	let "i=i+1"

done < <( ls -d faa_folders/* )

#wait till all processes are done
for pid in ${pids[*]}; do
    wait $pid
done

# delete tmp files
sbatch --wait --wrap="bash $prefix/delete_tmp.sh" -o slurms/slurm_delete_tmp_%J.out --time=24:0:0 -J delete_tmp
