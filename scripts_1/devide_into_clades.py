import pandas as pd
import sys

COLLAPSE=sys.argv[1]
DISTANCE=sys.argv[2]
LABEL=sys.argv[3]
SAMPLE_1_PATH=sys.argv[4]

FOLDER = "tree_files_sample_1"

if COLLAPSE == "P":
	collapse_column = "phylum"
elif COLLAPSE == "C":
	collapse_column = "class"
elif COLLAPSE == "O":
	collapse_column = "order"
elif COLLAPSE == "F":
	collapse_column = "family"
elif COLLAPSE == "G":
	collapse_column = "genus"
elif COLLAPSE == "S":
	collapse_column = "species"
else:
	print ("unvalid label!!")


current_label = "label_"+LABEL

#read sample 1 csv
df_sample_1 = pd.read_csv(SAMPLE_1_PATH, dtype="str")

#add column of clade_id
clades = df_sample_1[collapse_column]
df_sample_1.insert(loc = 1, column = 'clade_id', value = clades)

#change clade names that have "_" in the name
df_sample_1['clade_id'] = df_sample_1['clade_id'].str.replace('_','-')

#output metadata df
valid_labels = ["0","1"]
metadata_df = df_sample_1[df_sample_1[current_label].isin(valid_labels)]
metadata_df.to_csv(FOLDER + "/" +"tree_sample1_collapse_"+COLLAPSE+"_clade_with_metadata.tsv", sep = '\t', index = False)

#calculate clade sizes
only_clades_df = df_sample_1[["genome_id","clade_id"]]
len = only_clades_df.groupby('clade_id').apply(lambda x: len(x))
l = len.sort_values(ascending = False)
l.to_frame().rename(columns = {0:'size'}).to_csv(FOLDER + "/" + "tree_collapse_"+ COLLAPSE + "_clade_sizes.tsv", sep = '\t')

