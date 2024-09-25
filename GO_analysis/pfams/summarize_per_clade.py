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

df_results = pd.read_csv("combined_results_"+label+"_collapse_"+collapse+".csv")

df_results = df_results[df_results["enrichment"] == enrichment]

#PF,description,sum,num_of_clades
interpro_go_table = pd.read_csv("pfam_interpro_go_table.tsv", sep="\t")
interpro_go_table = interpro_go_table.drop(["pfam_id","pfam_description"], axis=1)

final_table = df_results.merge(interpro_go_table, left_on = "pfam_id", right_on = "pfam_id_short", how="left")
final_table = final_table.drop(["pfam_id_short"], axis=1)

print (final_table)
final_table.to_csv(label+"/"+label+"_"+enrichment+"_results_per_clade_collapse_"+collapse+"_with_GO_terms.tsv", index=False, sep="\t")
