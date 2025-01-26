import pandas as pd
import sys

COLLAPSE=sys.argv[1]
SAMPLE=sys.argv[2]
DISTANCE=sys.argv[3]
LABEL=sys.argv[4]
SAMPLE_1_PATH=sys.argv[5]

current_label = "label_"+LABEL

FOLDER = "tree_files_sample_"+SAMPLE

all_samples=pd.read_csv("samples_combined_per_clade_"+COLLAPSE+".csv", dtype="str")

df_sample = all_samples[["genome_id_"+SAMPLE,"clade_id"]]
df_sample = df_sample.rename(columns={"genome_id_"+SAMPLE: "genome_id"})

#read csv of the current sample
data = pd.read_csv(SAMPLE_1_PATH.replace("1.csv",SAMPLE+".csv"), dtype="str")

mer = df_sample.merge(data, on = 'genome_id', how = 'left')

valid_labels=["0","1"]

mer = mer[mer[current_label].isin(valid_labels)]

mer.to_csv(FOLDER + "/" +"tree_sample"+SAMPLE+"_collapse_"+COLLAPSE+"_clade_with_metadata.tsv", sep = '\t', index = False)

