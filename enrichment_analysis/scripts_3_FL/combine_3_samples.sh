#!/bin/bash

COLLAPSE=$1
LABEL=$2
#before
mkdir scoary_results_FL/raw_combined_samples

#loop throgh clade ids
while read CLADE;
do
  #combine to one file per clade
  less scoary_results_FL/raw/*clade_"$CLADE"_* | grep -v "Pfam_id,Clade_id,Sample_num,Bootstrap,Odds_ratio,Naive_p,Benjamini_H_p,Best_pairwise_comp_p,Worst_pairwise_comp_p,Empirical_p" >> scoary_results_FL/raw_combined_samples/"$LABEL"_collapse_"$COLLAPSE"__clade_"$CLADE"_combined.csv

  #add header to file
  sed -i '1s/^/Pfam_id,Clade_id,Sample_num,Bootstrap,Odds_ratio,Naive_p,Benjamini_H_p,Best_pairwise_comp_p,Worst_pairwise_comp_p,Empirical_p\n/' scoary_results_FL/raw_combined_samples/"$LABEL"_collapse_"$COLLAPSE"__clade_"$CLADE"_combined.csv

done < list_of_best_clades.tsv
