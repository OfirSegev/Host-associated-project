#!/bin/bash

prefix=$(echo /sci/labs/asafle/alexlevylab/icore-data/alex/workspace/HA_vs_NHA/data_sources/final_genomes/predictions/scripts_all_results_FL)

#global prameters
COLLAPSE=$1
LABEL=$2
sample_threshold=$3
#scoary parameters
scoary_odds_ratio_enriched_threshold=$4
scoary_odds_ratio_depleted_threshold=$5
scoary_enrichment_pval_threshold=$6
scoary_tree_pval_threshold=$7
tree_cutoff_type=$8 #fdr/naive (scoary)
enrichment_cutoff_type=$9 #fdr/naive (scoary + fisher)
#lmm parameters
lmm_odds_ratio_enriched_threshold=${10}
lmm_odds_ratio_depleted_threshold=${11}
lmm_enrichment_pval_threshold=${12}
lmm_tree_pval_threshold=${13}
#fisher parameters
fisher_odds_ratio_enriched_threshold=${14}
fisher_odds_ratio_depleted_threshold=${15}
fisher_enrichment_pval_threshold=${16}
#minimum clades for pivot table
MIN_CLADES=${17} 

#filter scoary, lmm and fisher and merge them
echo "running filter_and_merge_scoary_lmm_fisher_tables.py"
python3 $prefix/filter_and_merge_scoary_lmm_fisher_tables.py $COLLAPSE $LABEL $sample_threshold $scoary_odds_ratio_enriched_threshold $scoary_odds_ratio_depleted_threshold $scoary_enrichment_pval_threshold $scoary_tree_pval_threshold $tree_cutoff_type $enrichment_cutoff_type $lmm_odds_ratio_enriched_threshold $lmm_odds_ratio_depleted_threshold $lmm_enrichment_pval_threshold $lmm_tree_pval_threshold $fisher_odds_ratio_enriched_threshold $fisher_odds_ratio_depleted_threshold $fisher_enrichment_pval_threshold
#add evolink and make the final table
echo "running add_evolink_create_final_table.py"
python3 $prefix/add_evolink_create_final_table.py $COLLAPSE $LABEL
# make pivot, filter for only fisher
echo "running make_pivot_table.py"
python3 $prefix/make_pivot_table.py $COLLAPSE $LABEL $MIN_CLADES
# folder for exported tables
echo "running tables_to_export.py"
python3 $prefix/tables_to_export.py $COLLAPSE $LABEL

bash $prefix/print_slurms.sh
