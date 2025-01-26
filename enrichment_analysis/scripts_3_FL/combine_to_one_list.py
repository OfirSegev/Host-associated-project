import pandas as pd
#import dask.dataframe as dd
import concurrent.futures
import glob
import os
import sys

clade_sample=sys.argv[1]
COLLAPSE=sys.argv[2]
LABEL=sys.argv[3]

cols=["Pfam_id","Clade_id","Sample_num","Bootstrap","Odds_ratio","Naive_p","Benjamini_H_p","Best_pairwise_comp_p","Worst_pairwise_comp_p","Empirical_p"] # genome_names?
#output=pd.Dataframe(columns=cols)


file_list = glob.glob("*csv")
#"","SAMN05439496","SAMN05440341","Number_pos_present_in","Number_neg_present_in","Number_pos_not_present_in","Number_neg_not_present_in","Sensitivity","Specificity","Odds_ratio","Naive_p","Bonferroni_p","Benjamini_H_p","Max_Pairwise_comparisons","Max_supporting_pairs","Max_opposing_pairs","Best_pairwise_comp_p","Worst_pairwise_comp_p","Empirical_p"

def remove_prefix_suffix(text, prefix, suffix):
    if text.startswith(prefix):
        return text[len(prefix): len(text)-len(suffix)]
    print ("false file name")

def make_table(file):
	print (file)
	striped_name = remove_prefix_suffix(file, LABEL+"_collapse_"+COLLAPSE+"__", ".csv")
	clade = striped_name.split("_")[1]
	sample = striped_name.split("_")[3]
	if "bootstrap" in striped_name:
		bootstrap = striped_name.split("_")[5]
	else:
		bootstrap = "0"

	data=pd.read_csv(file, usecols=["Unnamed: 0","Odds_ratio","Naive_p","Benjamini_H_p","Best_pairwise_comp_p","Worst_pairwise_comp_p","Empirical_p"], dtype="str")

	data["Clade_id"] = clade
	data["Sample_num"] = sample
	data["Bootstrap"] = bootstrap

	return data


with concurrent.futures.ProcessPoolExecutor(max_workers=10) as executor:

	future = list(executor.map(make_table, file_list, chunksize=1))


df=pd.concat(future)
df=df.rename(columns={"Unnamed: 0": "Pfam_id"})
#sort columns in different order
df=df[["Pfam_id","Clade_id","Sample_num","Bootstrap","Odds_ratio","Naive_p","Benjamini_H_p","Best_pairwise_comp_p","Worst_pairwise_comp_p","Empirical_p"]]

df.to_csv(LABEL+"_collapse_"+COLLAPSE+"__"+clade_sample+"_combined.csv", index=False)

