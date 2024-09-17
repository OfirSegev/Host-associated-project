import pandas as pd
import sys

COLLAPSE=sys.argv[1]
LABEL=sys.argv[2]

df_scoary_lmm_fisher = pd.read_csv("final_results/scoary_lmm_fisher_combined_results_pfams_per_clade_"+LABEL+"_collapse_"+COLLAPSE+".csv")

#split enriched depleted
df_scoary_lmm_fisher_enriched = df_scoary_lmm_fisher[df_scoary_lmm_fisher["odds_ratio"] > 1]
df_scoary_lmm_fisher_depleted = df_scoary_lmm_fisher[df_scoary_lmm_fisher["odds_ratio"] < 1]
#evolink presence absence tables
df_evolink_P_A_enriched = pd.read_csv("evolink_results/filtered_enriched_"+LABEL+"_collapse_"+COLLAPSE+"__evolink_presence_absence_tables_combined.csv")
df_evolink_P_A_enriched["type"] = "evolink_P_A"
df_evolink_P_A_depleted = pd.read_csv("evolink_results/filtered_depleted_"+LABEL+"_collapse_"+COLLAPSE+"__evolink_presence_absence_tables_combined.csv")
df_evolink_P_A_depleted["type"] = "evolink_P_A"
#evolink counts tables
df_evolink_counts_enriched = pd.read_csv("evolink_results/filtered_enriched_"+LABEL+"_collapse_"+COLLAPSE+"__evolink_counts_tables_combined.csv")
df_evolink_counts_enriched["type"] = "evolink_counts"
df_evolink_counts_depleted = pd.read_csv("evolink_results/filtered_depleted_"+LABEL+"_collapse_"+COLLAPSE+"__evolink_counts_tables_combined.csv")
df_evolink_counts_depleted["type"] = "evolink_counts"

#relevant columns
columns = ["pfam_id","clade_id","type"]
#merge all tables
merged_enriched_tables = pd.concat([df_scoary_lmm_fisher_enriched[columns],df_evolink_P_A_enriched[columns],df_evolink_counts_enriched[columns]])
merged_depleted_tables = pd.concat([df_scoary_lmm_fisher_depleted[columns],df_evolink_P_A_depleted[columns],df_evolink_counts_depleted[columns]])
#convert format
merged_enriched_tables_converted = pd.get_dummies(merged_enriched_tables, columns=["type"])
merged_depleted_tables_converted = pd.get_dummies(merged_depleted_tables, columns=["type"])
#combine rows
merged_enriched_tables_converted_combined = merged_enriched_tables_converted.groupby(["pfam_id","clade_id"]).sum()
merged_depleted_tables_converted_combined = merged_depleted_tables_converted.groupby(["pfam_id","clade_id"]).sum()
#add counts column
merged_enriched_tables_converted_combined["tests_passed"] = merged_enriched_tables_converted_combined.sum(axis=1)
merged_depleted_tables_converted_combined["tests_passed"] = merged_depleted_tables_converted_combined.sum(axis=1)
#reset index
merged_enriched_tables_converted_combined.reset_index(inplace=True)
merged_depleted_tables_converted_combined.reset_index(inplace=True)
#add pfam description
pfam_description=pd.read_csv("/sci/labs/asafle/alexlevylab/icore-data/alex/workspace/HA_vs_NHA/misc_metadata/pfamID_to_description.tsv", sep="\t")
merged_enriched_tables_converted_combined = merged_enriched_tables_converted_combined.merge(pfam_description, left_on="pfam_id", right_on="pfam_id", how="left")
merged_depleted_tables_converted_combined = merged_depleted_tables_converted_combined.merge(pfam_description, left_on="pfam_id", right_on="pfam_id", how="left")
#change column names
merged_enriched_tables_converted_combined.rename(columns = {'type_evolink_P_A':'evolink_P_A','type_evolink_counts':'evolink_counts',
                                                            'type_scoary':'scoary','type_lmm':'lmm','type_fisher':'fisher'}, inplace=True)
merged_depleted_tables_converted_combined.rename(columns = {'type_evolink_P_A':'evolink_P_A','type_evolink_counts':'evolink_counts',
                                                            'type_scoary':'scoary','type_lmm':'lmm','type_fisher':'fisher'}, inplace=True)

print(merged_enriched_tables_converted_combined)

merged_enriched_tables_converted_combined.to_csv("final_results/final_table_all_tests_enriched_results_pfams_"+LABEL+"_collapse_"+COLLAPSE+".csv", index=False)
merged_depleted_tables_converted_combined.to_csv("final_results/final_table_all_tests_depleted_results_pfams_"+LABEL+"_collapse_"+COLLAPSE+".csv", index=False)
