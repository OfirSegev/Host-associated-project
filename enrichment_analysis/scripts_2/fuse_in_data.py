import pandas as pd
import sys

COLLAPSE=sys.argv[1]
DISTANCE=sys.argv[2]
LABEL=sys.argv[3]
SAMPLE_1_PATH=sys.argv[4]

current_label = "label_"+LABEL

#read metadata csv of sample 1
mer = pd.read_csv("./tree_files_sample_1/tree_sample1_collapse_"+COLLAPSE+"_clade_with_metadata.tsv", sep = '\t', dtype="str")

#reduce columns of mer
mer = mer[['genome_id','cluster','clade_id',current_label,'species']]

#read other samples

sample2 = pd.read_csv(SAMPLE_1_PATH.replace("1.csv","2.csv"), usecols=["genome_id",current_label,"species","cluster"], dtype="str")

sample3 = pd.read_csv(SAMPLE_1_PATH.replace("1.csv","3.csv"), usecols=["genome_id",current_label,"species","cluster"], dtype="str")


mer = mer.merge(sample2, on = 'cluster', how = 'left')

mer = mer.merge(sample3, on = 'cluster', how = 'left')

#change column names
mer = mer.rename(columns={"genome_id_x":"genome_id_1", current_label+"_x":current_label+"_1", "species_x":"species_1"})
mer = mer.rename(columns={"genome_id_y":"genome_id_2", current_label+"_y":current_label+"_2", "species_y":"species_2"})
mer = mer.rename(columns={"genome_id":"genome_id_3", current_label:current_label+"_3", "species":"species_3"})

#rearrange columns

mer = mer[["cluster","clade_id","genome_id_1",current_label+"_1","species_1","genome_id_2",current_label+"_2","species_2","genome_id_3",current_label+"_3","species_3"]]

print (mer)

mer = mer[~mer["cluster"].isna()]

mer.to_csv("samples_combined_per_clade_"+COLLAPSE+".csv", index=False)



