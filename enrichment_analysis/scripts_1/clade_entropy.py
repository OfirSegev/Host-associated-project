import pandas as pd
from skbio.diversity.alpha import shannon
from  sklearn.preprocessing import OrdinalEncoder
import numpy as np
import sys

#input collapse number

COLLAPSE=sys.argv[1]
LABEL=sys.argv[2]

df=pd.read_csv("tree_files_sample_1/tree_sample1_collapse_"+COLLAPSE+"_clade_with_metadata.tsv", sep="\t")

gb = df.groupby('clade_id')

#gb = gb.filter(lambda x: len(x) >= 100).groupby("clade_id")

current_label = "label_"+LABEL

for i in gb:
	cladeid = i[0]
	x = i[1]
	print("Clade " + str(cladeid))
	print(str(len(x)) + " members")
	print(x[current_label].value_counts())
	purity = shannon(x[current_label].value_counts().tolist())
	print("label ENTROPY: ", str(purity))
	print(x['phylum'].value_counts())
	print(x['genus'].value_counts())
	print("\n\n\n")

def entropy(inp):
	purity = shannon(inp[current_label].value_counts().tolist())
	return purity


a = gb.apply(lambda x: entropy(x))
print(a)

a.to_csv("tree_files_sample_1/sample1_"+COLLAPSE+"_clade_entropy.csv")
