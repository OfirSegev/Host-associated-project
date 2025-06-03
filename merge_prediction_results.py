import pandas as pd
import sys

LABEL = sys.argv[1]

#if not HvN , read the new csv

columns = ['metadata processed', 'genome_id', 'meta_not_stemmed', 'source',
       'Host (1) Env(0) unknown (-1)', 'pathogen(1) commensal(0)',
       'animal (1) plant(0)', 'human(1) mouse (0)',
       'GI tract(0) Oral(1) Skin (2)', 'Arabidopsis(1) Poplar(0)',
       'roots(1) shoots(0)', 'tresterial(1) aquatic(0)',
       'inverterbrate(1) non(-1)', 'Arthropod (1) non(-1)', 'Insect(1) non(-1)', 'Hemolymph(1) Other(0)']

if LABEL == "HvN":
	original = pd.read_csv("edited_mannual_labeling_sheet.csv", dtype="str", usecols=columns)
else:
	original = pd.read_csv("merged_prediction_results.csv", dtype="str")


pred = pd.read_csv("new_"+LABEL+"_labels_with_prediction_LR.csv", usecols = ["metadata processed","predicted_label_"+LABEL], dtype="str")

mer = original.merge(pred, how ="left")

print (original)
print (pred)
print (mer)

mer.to_csv("merged_prediction_results.csv", index=False)
