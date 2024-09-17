import pandas as pd
import concurrent.futures
import glob
import sys

#evolink_presence_absence_table_result_"$CLADE"_sample_$SAMPLE.tsv
#evolink_counts_tables_result_"$CLADE"_sample_$SAMPLE.tsv
#evolink_results/counts_tables/evolink_counts_tables_result_"$CLADE"_sample_$SAMPLE.tsv

FOLDER=sys.argv[1]
COLLAPSE=sys.argv[2]
LABEL=sys.argv[3]

file_list = glob.glob("evolink_results/"+FOLDER+"/*tsv")
print (file_list)

#make table with the new columns
def make_table(file):

	file_name = remove_suffix(file,".tsv")
	clade = file_name.split("_")[-3]
	sample = file_name.split("_")[-1]

	data=pd.read_csv(file, sep="\t")

	data["Clade_id"] = clade
	data["Sample_num"] = sample
	print(data)
	return data

#revise name
def remove_suffix(text, suffix):
	text = text.split("/")[-1]
	if text.endswith(suffix):
		return text[: len(text)-len(suffix)]
	print ("false file name")


with concurrent.futures.ProcessPoolExecutor(max_workers=30) as executor:

	future = list(executor.map(make_table, file_list, chunksize=13))

#combine tables
df=pd.concat(future)
df=df.rename(columns={"orthoID": "Pfam_id"})
#sort columns in different order
df=df[["Pfam_id","Clade_id","Sample_num","Prevalence_index","Evolink_index","scores","significance"]]

df.to_csv("evolink_results/"+LABEL+"_collapse_"+COLLAPSE+"__evolink_"+FOLDER+"_combined.csv", index=False)
