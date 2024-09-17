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

conda activate /sci/labs/asafle/ofirse/icore-data/mambaforge/envs/Evolink

export R_LIBS="/sci/labs/asafle/ofirse/icore-data/mambaforge/envs/Evolink/lib/R/library"
#export R_LIBS="/sci/labs/asafle/alexlevylab/icore-data/miniconda3/envs/Evolink/lib/R/library"

prefix=$(echo /sci/labs/asafle/alexlevylab/icore-data/alex/workspace/HA_vs_NHA/data_sources/final_genomes/predictions/scripts_evolink)

#inputs
COLLAPSE=$1
LABEL=$2
sample_threshold=$3

#make folder for evolink
mkdir files_for_evolink
mkdir evolink_results
mkdir evolink_results/presence_absence_tables
mkdir evolink_results/counts_tables

#job id array
declare -a job_IDs

#create folder and run analysis on all the best clades with all samples
while read CLADE;
do

  #loop through 3 samples
  for SAMPLE in 1 2 3
    do
       #organize files
       mkdir files_for_evolink/clade_"$CLADE"_sample_"$SAMPLE"

       INPUTFILE_T="t_traits__collapse_"$COLLAPSE"__clade_"$CLADE"_sample_"$SAMPLE".tsv"
       INPUTFILE_G="g_gene_presence_absence__collapse_"$COLLAPSE"__clade_"$CLADE"_sample_"$SAMPLE".tsv"
       INPUTFILE_G_COUNTS="g_gene_full_counts__collapse_"$COLLAPSE"__clade_"$CLADE"_sample_"$SAMPLE".tsv"
       NWK="n_newicktree_tree_"$COLLAPSE"__clade_"$CLADE"_sample_"$SAMPLE".newick"

       #Traits file
       sed  "1s/.*/Tip,Status/" files_for_scoary/clade_"$CLADE"_sample_"$SAMPLE"/t_traits* | tr "," "\t" > files_for_evolink/clade_"$CLADE"_sample_"$SAMPLE"/$INPUTFILE_T
       #Gene presence absence file
       sed 's/,/\t/g' files_for_scoary/clade_"$CLADE"_sample_"$SAMPLE"/g_gene_presence_absence* > files_for_evolink/clade_"$CLADE"_sample_"$SAMPLE"/$INPUTFILE_G
       #Gene counts file
       sed 's/,/\t/g' files_for_scoary/clade_"$CLADE"_sample_"$SAMPLE"/g_gene_full_counts* > files_for_evolink/clade_"$CLADE"_sample_"$SAMPLE"/$INPUTFILE_G_COUNTS
       #Tree file
       ln files_for_scoary/clade_"$CLADE"_sample_"$SAMPLE"/n_newicktree_tree* files_for_evolink/clade_"$CLADE"_sample_"$SAMPLE"

       cd files_for_evolink/clade_"$CLADE"_sample_"$SAMPLE"
       mkdir outputs_and_analysis
       mkdir slurms

       P_A_OUTPUT="evolink_output__presence_absence__collapse_"$COLLAPSE"__clade_"$CLADE"_sample_"$SAMPLE
       COUNT_OUTPUT="evolink_output__full_counts__collapse_"$COLLAPSE"__clade_"$CLADE"_sample_"$SAMPLE

       #run evolink presence absence
       job_ID1=($(sbatch -c3 --wrap="python3 /sci/labs/asafle/ofirse/icore-data/tools/Evolink/Evolink.py -g $INPUTFILE_G -t $INPUTFILE_T -n $NWK -o $P_A_OUTPUT" -o slurms/slurm_run_evolink_P_A_"$CLADE"_sample_"$SAMPLE"_%J.out | grep -oE "[[:digit:]]+"))
       #put in folder
       job_IDs+=($(sbatch --dependency=afterok:$(echo $job_ID1) --wrap="mv $P_A_OUTPUT/result.tsv ../../evolink_results/presence_absence_tables/evolink_presence_absence_table_result_"$CLADE"_sample_$SAMPLE.tsv" -o slurms/slurm_copy_results_evolink_P_A_"$CLADE"_sample_"$SAMPLE"_%J.out | grep -oE "[[:digit:]]+"))
       #run evolink gene counts
       job_ID3=($(sbatch -c3 --wrap="python3 /sci/labs/asafle/ofirse/icore-data/tools/Evolink/Evolink.py -g $INPUTFILE_G_COUNTS -t $INPUTFILE_T -n $NWK -o $COUNT_OUTPUT" -o slurms/slurm_run_evolink_gene_counts_"$CLADE"_sample_"$SAMPLE"_%J.out | grep -oE "[[:digit:]]+"))
       #put in folder
       job_IDs+=($(sbatch --dependency=afterok:$(echo $job_ID3) --wrap="mv $COUNT_OUTPUT/result.tsv ../../evolink_results/counts_tables/evolink_counts_tables_result_"$CLADE"_sample_$SAMPLE.tsv" -o slurms/slurm_copy_results_evolink_counts_"$CLADE"_sample_"$SAMPLE"_%J.out | grep -oE "[[:digit:]]+"))

#       echo ${job_IDs[@]}
       cd ../..
    done

done < list_of_best_clades.tsv

#combine all results presence_absence_tables
job_ID5=($(sbatch --dependency=afterok:$(echo ${job_IDs[@]} | tr ' ' :) --wrap="python3 $prefix/combine_to_one_list_evolink.py presence_absence_tables $COLLAPSE $LABEL" -o slurms/slurm_combine_evolink_presence_absence_%J.out | grep -oE "[[:digit:]]+"))

#combine all results counts_tables
job_ID6=($(sbatch --dependency=afterok:$(echo ${job_IDs[@]} | tr ' ' :) --wrap="python3 $prefix/combine_to_one_list_evolink.py counts_tables $COLLAPSE $LABEL" -o slurms/slurm_combine_evolink_counts_%J.out | grep -oE "[[:digit:]]+"))

#filter (sig, in 2 samples)
job_ID7=($(sbatch --dependency=afterok:$(echo $job_ID5),$(echo $job_ID6) --wrap="python3 $prefix/filter_evolink_results.py $COLLAPSE $LABEL $sample_threshold" -o slurms/slurm_filter_evolink_results_%J.out | grep -oE "[[:digit:]]+"))

#wait to finish
job_ID8=($(sbatch --wait --dependency=afterok:$(echo $job_ID7) --wrap="echo evolink metascript done" -o slurms/slurm_evolink_metascript_done_%J.out | grep -oE "[[:digit:]]+"))