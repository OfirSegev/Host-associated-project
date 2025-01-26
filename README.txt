Here is presented the main script used for the identification of host-associated bacterial proteins and protein domains.

enrichment_analysis folder:
*run_all_analysis_megascript_new.sh contains the code to that utilizes all of the code required for the analysis. Parameters can be optimized within this script.
workflow:
* scripts_1 - division into taxonomic clades.
* scripts_2 - genetation of phylogenetic trees.
* scripts_3 - utilize the gene/trait metadata and apply Scoary and Fisher exact test
* scripts_pyseer - apply pyseer's lmm
* scripts_evolink - apply evolink presence/absence version and counts based version
* scripts_all_results - combine all results for exprotable tables.
#phylophlan databases are neeeded in order for the script to run smoothly

operons folder:
* pipleline for finding Pfam domians that co-occure together within gene clusters 

heatmap folder:
* Code used for visualization of the top enriched results spread across multiple taxa

GO_analysis:
* Code used for finding enriched GO terms in based on Pfam and AFCs results with visualization using R.
