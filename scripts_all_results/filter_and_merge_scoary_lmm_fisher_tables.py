import pandas as pd
import sys

#global prameters
COLLAPSE=sys.argv[1]
LABEL=sys.argv[2]
sample_threshold = float(sys.argv[3])
#scoary parameters
scoary_odds_ratio_enriched_threshold = float(sys.argv[4])
scoary_odds_ratio_depleted_threshold = float(sys.argv[5])
scoary_enrichment_pval_threshold = float(sys.argv[6])
scoary_tree_pval_threshold = float(sys.argv[7])
tree_cutoff_type = sys.argv[8] #fdr/naive (scoary)
enrichment_cutoff_type = sys.argv[9] #fdr/naive (scoary + fisher)
#lmm parameters
lmm_odds_ratio_enriched_threshold = float(sys.argv[10])
lmm_odds_ratio_depleted_threshold = float(sys.argv[11])
lmm_enrichment_pval_threshold = float(sys.argv[12])
lmm_tree_pval_threshold = float(sys.argv[13])
#fisher parameters
fisher_odds_ratio_enriched_threshold=float(sys.argv[14])
fisher_odds_ratio_depleted_threshold=float(sys.argv[15])
fisher_enrichment_pval_threshold=float(sys.argv[16])

#choose the right columns for filtering (either fdr corrected or native)
if (enrichment_cutoff_type == "fdr"):
	enrichment_cutoff_column = "Naive_pval_corrected"
else:
	enrichment_cutoff_column = "Benjamini_H_p"

if(tree_cutoff_type == "fdr"):
	tree_cutoff_column = "Worst_pval_corrected"
else:
	tree_cutoff_column = "Worst_pairwise_comp_p"


folder_location="final_results/"

#count samples
def count_samples(df):
	counts=len(df)
	return counts

#return final dataframe of enriched/depleted pfams
def analyze_data(df,flag,type):

	sample_counts_total=df.groupby(["Pfam_id","Clade_id"]).apply(count_samples)

	#filter table by samples
	filtered_by_sample_counts= sample_counts_total[sample_counts_total>=sample_threshold]
	print (filtered_by_sample_counts)
	
    #calculate average score of all samples
	averages = df.groupby(["Pfam_id","Clade_id"]).mean()
	print (averages)
	
	#keep relevant columns
	if( (type=="scoary") | (type=="fisher") ):
		averages = averages[["Odds_ratio",enrichment_cutoff_column]]
	elif(type=="lmm"):
		averages = averages[["Odds_ratio","q_value_chisq"]]

	#merge filtered by sample with average of odds_ratio, pvalue
	mer = pd.merge(averages, filtered_by_sample_counts.to_frame(), left_index=True, right_index=True, how="inner")
	mer.reset_index(inplace=True)
	#change column names to be unite
	mer.columns = ["pfam_id","clade_id","odds_ratio","p_value","num_of_samples"]
	print (mer)
	
	#add description of pfam
	enriched_pfams_with_info = mer.merge(pfam_description, on="pfam_id", how="left")
	print (enriched_pfams_with_info)
	
	if(flag):
		enriched_pfams_with_info.to_csv(folder_location+type+"_enriched_pfams_per_clade_"+LABEL+"_collapse_"+COLLAPSE+".csv",index=False)
	else:
		enriched_pfams_with_info.to_csv(folder_location+type+"_depleted_pfams_per_clade_"+LABEL+"_collapse_"+COLLAPSE+".csv",index=False)

	return enriched_pfams_with_info


#run output for a specific type (scoary,lmm,pyseer)
def run_main_output(df_enriched, df_depleted, type):

	if(len(df_enriched) == 0):
		output_enriched = pd.DataFrame(columns=['pfam_id', 'clade_id', 'odds_ratio', 'p_value', 'num_of_samples','description'])
		output_enriched.to_csv(folder_location+type+"_enriched_pfams_per_clade_"+LABEL+"_collapse_"+COLLAPSE+".csv",index=False)
	else:
		output_enriched = analyze_data(df_enriched,True,type)
	
	if(len(df_depleted) == 0):
		output_depleted = pd.DataFrame(columns=['pfam_id', 'clade_id', 'odds_ratio', 'p_value', 'num_of_samples','description'])
		output_depleted.to_csv(folder_location+type+"_depleted_pfams_per_clade_"+LABEL+"_collapse_"+COLLAPSE+".csv",index=False)
	else:
		output_depleted = analyze_data(df_depleted,False,type)

	#merge depleted/enriched
	output = pd.concat([output_enriched,output_depleted])
	#add type column
	output["type"] = type

	output.to_csv(folder_location+type+"_combined_results_pfams_per_clade_"+LABEL+"_collapse_"+COLLAPSE+".csv",index=False)

	return output

######### main ##########

#read pfam description table
pfam_description=pd.read_csv("/sci/labs/asafle/alexlevylab/icore-data/alex/workspace/HA_vs_NHA/misc_metadata/pfamID_to_description.tsv", sep="\t")

#read tables
df_scoary=pd.read_csv("all_clade_sample_collapse_"+COLLAPSE+"_scoary_output_"+LABEL+"_fdr_updated.csv")
df_lmm=pd.read_csv("all_clade_sample_collapse_"+COLLAPSE+"_lmm_output_"+LABEL+"_fdr_updated.csv")
#df_evolink_P_A
#df_evolink_counts

#table filtering by thresholds scoary
df_enriched_scoary=df_scoary.loc[(df_scoary[tree_cutoff_column]<=scoary_tree_pval_threshold) & (df_scoary[enrichment_cutoff_column]<=scoary_enrichment_pval_threshold) & (df_scoary["Odds_ratio"]>=scoary_odds_ratio_enriched_threshold)]
df_depleted_scoary=df_scoary.loc[(df_scoary[tree_cutoff_column]<=scoary_tree_pval_threshold) & (df_scoary[enrichment_cutoff_column]<=scoary_enrichment_pval_threshold) & (df_scoary["Odds_ratio"]<=scoary_odds_ratio_depleted_threshold)]

#table filtering by thresholds lmm
df_enriched_lmm=df_lmm[(df_lmm['Odds_ratio']>=lmm_odds_ratio_enriched_threshold) & (df_lmm['q_value_chisq'] <= lmm_enrichment_pval_threshold) & (df_lmm['q_value_lrt'] <= lmm_tree_pval_threshold)]
df_depleted_lmm=df_lmm[(df_lmm['Odds_ratio']<=lmm_odds_ratio_depleted_threshold) & (df_lmm['q_value_chisq'] <= lmm_enrichment_pval_threshold) & (df_lmm['q_value_lrt'] <= lmm_tree_pval_threshold)]

#table filtering by thresholds fisher
df_enriched_fisher=df_scoary.loc[(df_scoary[enrichment_cutoff_column]<=fisher_enrichment_pval_threshold) & (df_scoary["Odds_ratio"]>=fisher_odds_ratio_enriched_threshold)]
df_depleted_fisher=df_scoary.loc[(df_scoary[enrichment_cutoff_column]<=fisher_enrichment_pval_threshold) & (df_scoary["Odds_ratio"]<=fisher_odds_ratio_depleted_threshold)]

#run main function for each type
output_scoary = run_main_output(df_enriched_scoary, df_depleted_scoary, "scoary")
output_lmm = run_main_output(df_enriched_lmm, df_depleted_lmm, "lmm")
output_fisher = run_main_output(df_enriched_fisher, df_depleted_fisher, "fisher")

#merge all outputs to one
output_combined = pd.concat([output_scoary,output_lmm,output_fisher])
output_combined.to_csv("final_results/scoary_lmm_fisher_combined_results_pfams_per_clade_"+LABEL+"_collapse_"+COLLAPSE+".csv",index=False)




