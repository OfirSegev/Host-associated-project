#!/bin/bash
# Exit immediately if any command exits with a non-zero status
set -e

############################################################################
#global parameters
COLLAPSE=$1
DISTANCE="0.8" # just a relic, not relevant anymore
LABEL=$2
PFAM_THRESHOLD="0.97" #threshold for maximum shared pfam between all (1-THRESHOLD = minimum)
MIN_CLADE_SIZE="20" #minimum clade size to do the anlysis on
MIN_ENTROPY="0.25" #minimum label entropy per clade to do the anlysis on (only for larger than 100)
Pfam_E_value_Threshold="1e-5" # E_value threshold of hmm table
HMM_coverage_threshold="0.7" # minimum HMM coverage for pfam
sample_threshold="2" # minimum samples to pass
MIN_CLADES="3" # for pivot table
#bootstrap_threshold="10"

FULL_LENGTH_RUN=$3  #"both", "full_only", "pfam_only"
BACTERIAL_OR_FUNGI=$4  # 'B' or 'F' (fungi has scripts suitable for pfam only)
CHECKPOINT=$5 #3, none
CHECKPOINT_FOLDER=$6 

###
FUNGI_TREE_PATH="/sci/labs/asafle/alexlevylab/icore-data/alex/workspace/HA_vs_NHA/data_sources/final_genomes/fungi/mimi/iq_tree/right_label_trimmed.fasta.treefile.rooted"

### tests parameters ###
#scoary parameters
scoary_odds_ratio_enriched_threshold="1.5"
scoary_odds_ratio_depleted_threshold="0.67"
scoary_enrichment_pval_threshold="0.05"
scoary_tree_pval_threshold="0.05"
tree_cutoff_type="naive" #fdr/naive (scoary)
enrichment_cutoff_type="naive" #fdr/naive (scoary + fisher)
#lmm parameters
lmm_odds_ratio_enriched_threshold="1.05"
lmm_odds_ratio_depleted_threshold="0.95"
lmm_enrichment_pval_threshold="0.01"
lmm_tree_pval_threshold="0.05"
#fisher parameters
fisher_odds_ratio_enriched_threshold="2"
fisher_odds_ratio_depleted_threshold="0.5"
fisher_enrichment_pval_threshold="0.001"


############################################################################
#?Alphafold clusters master sheet path
AFC_MASTER_SHEET="/sci/labs/asafle/alexlevylab/icore-data/alex/workspace/HA_vs_NHA/data_sources/final_genomes/genome_data/new_full_protein_data/mmseqs_to_AFC_master_sheet.csv"
############################################################################
#sample needs to end with number_1.csv
#?Root vs Shoot
#SAMPLE_1_PATH="/sci/labs/asafle/alexlevylab/icore-data/alex/workspace/HA_vs_NHA/data_sources/final_genomes/samples/RvS_samples/new_final_genome_0.8_new_sample_RvS_number_1.csv"
#?insects
#SAMPLE_1_PATH="/sci/labs/asafle/alexlevylab/icore-data/alex/workspace/HA_vs_NHA/data_sources/final_genomes/insect_project/InsvAni_samples/new_final_genome_0.8_new_sample_InsvAni_number_1.csv"
#?mags
#SAMPLE_1_PATH="/sci/labs/asafle/alexlevylab/icore-data/alex/workspace/HA_vs_NHA/data_sources/final_genomes/samples/only_mags_HvN_samples/new_updated_final_genome_0.8_new_sample_HvN_number_1.csv"
#?toy 2 small genera
#SAMPLE_1_PATH="/sci/labs/asafle/alexlevylab/icore-data/alex/workspace/HA_vs_NHA/data_sources/final_genomes/samples/two_small_gen_toy/new_final_genome_0.8_new_sample_HvN_number_1.csv"
#?toy paraburkhulderia
#SAMPLE_1_PATH="/sci/labs/asafle/alexlevylab/icore-data/alex/workspace/HA_vs_NHA/data_sources/final_genomes/samples/toy_Paraburkholderia/new_final_genome_0.8_new_sample_HvN_number_1.csv"
#?toy Flavobacterium_Bradyrhizobium
#SAMPLE_1_PATH="/sci/labs/asafle/alexlevylab/icore-data/alex/workspace/HA_vs_NHA/data_sources/final_genomes/samples/toy_Flavobacterium_Bradyrhizobium/new_final_genome_0.8_new_sample_HvN_number_1.csv"
#?toy burkhulderia AvP
#SAMPLE_1_PATH="/sci/labs/asafle/alexlevylab/icore-data/alex/workspace/HA_vs_NHA/data_sources/final_genomes/samples/toy_Burkholderia_AvP/new_final_genome_0.8_new_sample_AvP_number_1.csv"
#?fungi HvN
#SAMPLE_1_PATH="/sci/labs/asafle/alexlevylab/icore-data/alex/workspace/HA_vs_NHA/data_sources/final_genomes/samples/HvN_fungi_samples/HvN_fungi_sample_number_1.csv"
#?fungi_AvP
#SAMPLE_1_PATH="/sci/labs/asafle/alexlevylab/icore-data/alex/workspace/HA_vs_NHA/data_sources/final_genomes/samples/AvP_fungi_samples/AvP_fungi_sample_number_1.csv"
############################################################################
#?udated AvP/HvN
#SAMPLE_1_PATH="/sci/labs/asafle/alexlevylab/icore-data/alex/workspace/HA_vs_NHA/data_sources/final_genomes/samples/new_samples/new_final_genome_0.8_new_sample_AvP_number_1.csv"
SAMPLE_1_PATH="/sci/labs/asafle/alexlevylab/icore-data/alex/workspace/HA_vs_NHA/data_sources/final_genomes/samples/new_samples/new_final_genome_0.8_new_sample_HvN_number_1.csv"
#?updated_mags
#SAMPLE_1_PATH="/sci/labs/asafle/alexlevylab/icore-data/alex/workspace/HA_vs_NHA/data_sources/final_genomes/samples/only_mags_29.7.24/new_final_genome_0.8_new_sample_HvN_number_1.csv"
############################################################################

if [ ! "$FULL_LENGTH_RUN" = "full_only" ] 
then
	echo in first part
	if [ "$CHECKPOINT" = "3" ]
	then
		echo in second part
		cd $CHECKPOINT_FOLDER

		job_ID3=($(sbatch --wrap="bash ../scripts_3/all_clades_sample_scoary.sh $COLLAPSE $DISTANCE $PFAM_THRESHOLD $LABEL $SAMPLE_1_PATH $Pfam_E_value_Threshold $HMM_coverage_threshold" -o slurms/slurm_all_clades_sample_scoary_metascript_%J.out --time=24:0:0 -J run_scoary| grep -oE "[[:digit:]]+"))

	else

		mkdir "$LABEL"_collapse_"$COLLAPSE"_distance_"$DISTANCE"
		cd "$LABEL"_collapse_"$COLLAPSE"_distance_"$DISTANCE"
		mkdir slurms
		mkdir tree_files_sample_1

		job_ID1=($(sbatch --wrap="bash ../scripts_1/run_details_different_cutoff.sh $COLLAPSE $DISTANCE $LABEL $SAMPLE_1_PATH $MIN_CLADE_SIZE $MIN_ENTROPY" -o slurms/slurm_run_details_different_cutoff_"$COLLAPSE"_%J.out --time=10:0:0 -J split_clades | grep -oE "[[:digit:]]+"))

		if [ "$BACTERIAL_OR_FUNGI" = "F" ]
		then
			job_ID2=($(sbatch --dependency=afterok:$(echo $job_ID1) --wrap="bash ../scripts_2_fungi/make_subtrees_per_clade_metascript_fungi.sh $COLLAPSE $LABEL $FUNGI_TREE_PATH" -o slurms/slurm_make_subtrees_per_clade_metascript_%J.out --time=40:0:0 -J subtrees | grep -oE "[[:digit:]]+"))
		else
			job_ID2=($(sbatch --dependency=afterok:$(echo $job_ID1) --wrap="bash ../scripts_2/make_subtrees_per_clade_metascript.sh $COLLAPSE $DISTANCE $LABEL $SAMPLE_1_PATH" -o slurms/slurm_make_subtrees_per_clade_metascript_%J.out --time=30:0:0 -J subtrees | grep -oE "[[:digit:]]+"))
		fi

		job_ID3=($(sbatch --dependency=afterok:$(echo $job_ID2) --wrap="bash ../scripts_3/all_clades_sample_scoary.sh $COLLAPSE $DISTANCE $PFAM_THRESHOLD $LABEL $SAMPLE_1_PATH $Pfam_E_value_Threshold $HMM_coverage_threshold" -o slurms/slurm_all_clades_sample_scoary_metascript_%J.out --time=24:0:0 -J run_scoary| grep -oE "[[:digit:]]+"))

	fi

echo near job4
job_ID4=($(sbatch --dependency=afterok:$(echo $job_ID3) --wrap="bash ../scripts_pyseer/run_pyseer_metascript.sh $COLLAPSE $LABEL $PFAM_THRESHOLD" -o slurms/slurm_run_pyseer_metascript_%J.out --time=10:0:0 -J run_pyseer | grep -oE "[[:digit:]]+"))

job_ID5=($(sbatch --dependency=afterok:$(echo $job_ID3) --wrap="bash ../scripts_evolink/run_evolink_metascript.sh $COLLAPSE $LABEL $sample_threshold" -o slurms/slurm_evolink_metascript_%J.out --time=4:0:0 -J evolink | grep -oE "[[:digit:]]+"))

job_ID6=($(sbatch --dependency=afterok:$(echo $job_ID4),$(echo $job_ID5) --wrap="bash ../scripts_all_results/all_results_metascript.sh $COLLAPSE $LABEL $sample_threshold $scoary_odds_ratio_enriched_threshold $scoary_odds_ratio_depleted_threshold $scoary_enrichment_pval_threshold $scoary_tree_pval_threshold $tree_cutoff_type $enrichment_cutoff_type $lmm_odds_ratio_enriched_threshold $lmm_odds_ratio_depleted_threshold $lmm_enrichment_pval_threshold $lmm_tree_pval_threshold $fisher_odds_ratio_enriched_threshold $fisher_odds_ratio_depleted_threshold $fisher_enrichment_pval_threshold $MIN_CLADES" -o slurms/slurm_all_results_metascript_%J.out --time=2:0:0 -J all_results | grep -oE "[[:digit:]]+"))

fi

#####################################
if [ "$FULL_LENGTH_RUN" = "both" -o "$FULL_LENGTH_RUN" = "full_only" ] 
then
	echo in full part

	if [ "$FULL_LENGTH_RUN" = "full_only" ] 
	then
		cd $CHECKPOINT_FOLDER
        
		job_ID3_FL=($(sbatch --wrap="bash ../scripts_3_FL/all_clades_sample_scoary.sh $COLLAPSE $DISTANCE $PFAM_THRESHOLD $LABEL $SAMPLE_1_PATH $Pfam_E_value_Threshold $HMM_coverage_threshold $AFC_MASTER_SHEET" -o slurms/slurm_all_clades_sample_scoary_FL_metascript_%J.out --time=20:0:0 -J run_scoary| grep -oE "[[:digit:]]+"))

	else
		job_ID3_FL=($(sbatch --dependency=afterok:$(echo $job_ID6) --wrap="bash ../scripts_3_FL/all_clades_sample_scoary.sh $COLLAPSE $DISTANCE $PFAM_THRESHOLD $LABEL $SAMPLE_1_PATH $Pfam_E_value_Threshold $HMM_coverage_threshold $AFC_MASTER_SHEET" -o slurms/slurm_all_clades_sample_scoary_FL_metascript_%J.out --time=20:0:0 -J run_scoary| grep -oE "[[:digit:]]+"))
	fi

	job_ID4_FL=($(sbatch --dependency=afterok:$(echo $job_ID3_FL) --wrap="bash ../scripts_pyseer_FL/run_pyseer_metascript.sh $COLLAPSE $LABEL $PFAM_THRESHOLD" -o slurms/slurm_run_pyseer_metascript_FL_%J.out --time=14:0:0 -J run_pyseer | grep -oE "[[:digit:]]+"))

	job_ID5_FL=($(sbatch --dependency=afterok:$(echo $job_ID3_FL) --wrap="bash ../scripts_evolink_FL/run_evolink_metascript.sh $COLLAPSE $LABEL $sample_threshold" -o slurms/slurm_evolink_metascript_FL_%J.out --time=4:0:0 -J evolink | grep -oE "[[:digit:]]+"))

	job_ID6_FL=($(sbatch --dependency=afterok:$(echo $job_ID4_FL),$(echo $job_ID5_FL) --wrap="bash ../scripts_all_results_FL/all_results_metascript.sh $COLLAPSE $LABEL $sample_threshold $scoary_odds_ratio_enriched_threshold $scoary_odds_ratio_depleted_threshold $scoary_enrichment_pval_threshold $scoary_tree_pval_threshold $tree_cutoff_type $enrichment_cutoff_type $lmm_odds_ratio_enriched_threshold $lmm_odds_ratio_depleted_threshold $lmm_enrichment_pval_threshold $lmm_tree_pval_threshold $fisher_odds_ratio_enriched_threshold $fisher_odds_ratio_depleted_threshold $fisher_enrichment_pval_threshold $MIN_CLADES" -o slurms/slurm_all_results_metascript_FL_%J.out --time=2:0:0 -J all_results | grep -oE "[[:digit:]]+"))

fi
