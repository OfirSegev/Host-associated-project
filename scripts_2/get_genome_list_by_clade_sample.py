import pandas as pd
import sys

COLLAPSE=sys.argv[1]
LABEL=sys.argv[2]

current_label = "label_"+LABEL

high_score_clades = pd.read_csv("scored_clades_collapse_"+COLLAPSE+".csv", usecols=["clade_id"], dtype="str")
high_score_clades = high_score_clades["clade_id"].tolist()

df=pd.read_csv("samples_combined_per_clade_"+COLLAPSE+".csv", dtype="str")

valid_labels=["0","1"]

def generate_genome_lists(clade,sample):

	genome_list = df[(df["clade_id"]==clade) & (df[current_label+"_"+sample].isin(valid_labels))]["genome_id_"+sample]
	genome_list.to_csv("genome_lists_of_clades/list_of_clade_"+clade+"_sample_"+sample+".tsv", index=False, header=False)

for clade in high_score_clades:

	generate_genome_lists(clade,"1")
	generate_genome_lists(clade,"2")
	generate_genome_lists(clade,"3")

