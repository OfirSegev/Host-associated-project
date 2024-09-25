import pandas as pd
import sys

label = sys.argv[1]
collapse = sys.argv[2]
table_type = sys.argv[3]

if (label == "HvN"):
	if (table_type == "enriched"):
		enrichment="Host-associated"
	else:
		enrichment="Environmental"
elif (label == "AvP"):
	if (table_type == "enriched"):
		enrichment="Animal-associated"
	else:
		enrichment="Plant-associated"

df_results = pd.read_csv("FL_combined_results_"+label+"_collapse_"+collapse+".csv")

df_results = df_results[df_results["enrichment"] == enrichment]

uniprot_go_table = pd.read_csv("AFC_repID_to_GO_master_sheet__memIDCluster_within_PercentPresence_gte_0.75_updated_obsolete_alt.csv")

final_table = df_results.merge(uniprot_go_table, left_on = "pfam_id", right_on = "repID", how="left")
final_table = final_table.drop(["repID"], axis=1)

print (final_table)
final_table.to_csv(label+"/"+label+"_"+enrichment+"_results_per_clade_collapse_"+collapse+"_with_GO_terms.tsv", index=False, sep="\t")
