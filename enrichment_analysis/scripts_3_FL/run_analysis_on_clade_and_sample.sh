#!/bin/bash

#inputs 1 - clade number, 2 -sample number 3 - collapse , 4 - pfam abundance threshold, 5 - label type
# 8 - afc master sheet
#remember to turn on socary enrivroment!!!!!

CLADE=$1
SAMPLE=$2
COLLAPSE=$3
THRESHOLD=$4
LABEL=$5
Pfam_E_value_Threshold=$6
HMM_coverage_threshold=$7
AFC_MASTER_SHEET=$8
sample_percentage="0.7"
num_of_bootstrap="0"

prefix=$(echo /sci/labs/asafle/alexlevylab/icore-data/alex/workspace/HA_vs_NHA/data_sources/final_genomes/predictions/scripts_3_FL)

#create bootstraps lists
for i in $(seq 1 1 $num_of_bootstrap)
do
   python3 $prefix/bootstrap_sample.py $CLADE $SAMPLE $i $COLLAPSE $sample_percentage
done

#calculate clade size
clade_size=$(grep $CLADE, scored_clades_collapse_$COLLAPSE.csv | awk -F "," '{print $2}')

if [ $clade_size -lt 6000 ]
then
   #create files for scoary To Do -- make G files
   job_ID1=($(sbatch -c3 --time=5:0:0 --wrap="bash $prefix/make_g_file_AFC.sh $CLADE $SAMPLE $COLLAPSE $THRESHOLD $LABEL $AFC_MASTER_SHEET" -o slurm_FL_makefile_clade_"$CLADE"_sample_"$SAMPLE"_%J.out | grep -oE "[[:digit:]]+"))
else
   job_ID1=($(sbatch -c3 --mem=100G --time=5:0:0 --wrap="bash $prefix/make_g_file_AFC.sh $CLADE $SAMPLE $COLLAPSE $THRESHOLD $LABEL $AFC_MASTER_SHEET" -o slurm_FL_makefile_clade_"$CLADE"_sample_"$SAMPLE"_%J.out | grep -oE "[[:digit:]]+"))
fi

sbatch --wait --dependency=afterok:$(echo $job_ID1) --wrap="echo make_files_done" -o slurm_FL_makefile_done_clade_"$CLADE"_sample_"$SAMPLE"_%J.out

#wait
#organize files (needs to wait for finish)
#mkdir files_for_scoary
#mkdir files_for_scoary/clade_"$CLADE"_sample_"$SAMPLE"
mv *clade_"$CLADE"_sample_"$SAMPLE"* files_for_scoary/clade_"$CLADE"_sample_"$SAMPLE"

cd files_for_scoary/clade_"$CLADE"_sample_"$SAMPLE"
#mkdir outputs_and_analysis
#mkdir slurms

#create job id array
declare -a job_IDs
#make vriable names
clade_sample=$(echo clade_"$CLADE"_sample_"$SAMPLE")

#run scoary on all genomes

#differet resources for different clade sizes
if [ $clade_size -lt 6000 ]
then
   #24 cores as default
   job_IDs+=($(sbatch --wrap="$prefix/command_t_g_n.sh t_traits__collapse_'$COLLAPSE'__$clade_sample.csv g_AFC_presence_absence__collapse_'$COLLAPSE'__$clade_sample.csv n_newicktree_tree_'$COLLAPSE'__$clade_sample.newick $clade_sample $COLLAPSE $LABEL" -c18 --time=1:0:0 -o slurms/slurm_"$clade_sample"_all_genomes_FL_%J.out | grep -oE "[[:digit:]]+"))
else
   #memory needed: clade_size/40
   mem=$(expr $clade_size / 40)
   cpu=$(expr $clade_size / 320)
   job_IDs+=($(sbatch --wrap="$prefix/command_t_g_n.sh t_traits__collapse_'$COLLAPSE'__$clade_sample.csv g_AFC_presence_absence__collapse_'$COLLAPSE'__$clade_sample.csv n_newicktree_tree_'$COLLAPSE'__$clade_sample.newick $clade_sample $COLLAPSE $LABEL" -c$cpu --mem="$mem"G --time=3:0:0 -o slurms/slurm_"$clade_sample"_all_genomes_FL_%J.out | grep -oE "[[:digit:]]+"))
 fi


#run scoary on bootstraped
for i in $(seq 1 1 $num_of_bootstrap)
do
	echo $i
	job_IDs+=($(sbatch --wrap="$prefix/scoary_command_bootstrap.sh t_traits__collapse_'$COLLAPSE'__$clade_sample.csv g_AFC_presence_absence__collapse_'$COLLAPSE'__$clade_sample.csv n_newicktree_tree_'$COLLAPSE'__$clade_sample.newick $clade_sample $i bootstrap_'$i'_collapse_'$COLLAPSE'__$clade_sample.txt $COLLAPSE $LABEL"  -c18 --time=4:0:0 -o slurms/slurm_"$clade_sample"_bootstrap_"$i"_%J.out | grep -oE "[[:digit:]]+"))

done

#arrange output files inside one folder
# merge data to one table (remember scoary needs to finish)
job_ID2=($(sbatch --dependency=afterok:$(echo ${job_IDs[@]} | tr ' ' :) --wrap="$prefix/merge_scoary_output.sh $clade_sample $COLLAPSE $prefix $LABEL" -o slurms/slurm_merge_script_%J.out | grep -oE "[[:digit:]]+"))

sbatch --wait --dependency=afterok:$(echo $job_ID2) --wrap="echo merge scoary output done" -o slurm_merge_scoary_done_clade_"$CLADE"_sample_"$SAMPLE"_%J.out

