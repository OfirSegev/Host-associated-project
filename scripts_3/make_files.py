import pandas as pd
pd.set_option("display.max_columns", None)
import concurrent.futures
import subprocess
import sys
import os

#inputs 1- clade number , 2- sample number , 3-collapse number , 4 - pfam abundance threshold, 5 - label type, 6 - e-value threshold, 7 - hmm coverage

prefix="/sci/labs/asafle/alexlevylab/icore-data/alex/workspace/HA_vs_NHA/data_sources/final_genomes/predictions/scripts_3/"

def make_traitfile(file, CLADE, NAME, LABEL):

	df = pd.read_csv(file, sep="\t", dtype="str")

	print(df)
	print(df.columns)

	df2 = df.loc[df['clade_id']== str(CLADE)]

	current_label = "label_"+LABEL

	t = df2[['genome_id', current_label]]

	out = t.rename(columns = {'genome_id' : "", current_label : LABEL})

	out.to_csv("t_traits__" + NAME + ".csv", index = False)

	return t

def onehotencoder(genome_id):
	print(genome_id)
	df=pd.read_csv("/sci/labs/asafle/alexlevylab/icore-data/alex/workspace/HA_vs_NHA/data_sources/final_genomes/genome_data/new_pfams/"+genome_id+".pfam", sep="\t", usecols=["Pfam_id","E-value","Coverage"])

	#filter for E-value threshold
	df = df[df["E-value"] <= Pfam_E_value_Threshold]
	#filter for hmm coverage
	df = df[df["Coverage"] >= HMM_coverage_threshold]

	df = df[["Pfam_id"]]

	df.insert(loc=0, column='genome_id', value=genome_id)
#	print(df)
	pv = df.pivot_table(index = 'Pfam_id', columns = 'genome_id', aggfunc = 'size')
#	print(pv)
	return pv

def get_tree(CLADE,SAMPLE,COLLAPSE):
	subprocess.Popen(['bash', prefix+'get_tree.sh', str(CLADE), str(SAMPLE), str(COLLAPSE)])
	os.wait()
'''
removes pfams rows with pfams that apear in most genomes (threshold)
'''


if __name__ == "__main__":

	SAMPLE = sys.argv[2]
	COLLAPSE = sys.argv[3]
	THRESHOLD = float(sys.argv[4])
	LABEL = sys.argv[5]
	file = "tree_files_sample_"+SAMPLE+"/tree_sample"+SAMPLE+"_collapse_"+COLLAPSE+"_clade_with_metadata.tsv"
	CLADE = sys.argv[1]
	Pfam_E_value_Threshold=float(sys.argv[6])
	HMM_coverage_threshold=float(sys.argv[7])

	NAME = "collapse_" + COLLAPSE + "__clade_" + str(CLADE) + "_sample_" + SAMPLE

	genome_list = make_traitfile(file, CLADE, NAME, LABEL)

	l = genome_list['genome_id'].tolist()

	futures = []

	with concurrent.futures.ProcessPoolExecutor(max_workers=30) as executor:
		output = executor.map(onehotencoder, l, chunksize=1)

	for i in output:
		futures.append(i)

	df_merged = pd.concat(futures, axis = 1, join = 'outer')
	df_merged.fillna(0, inplace = True) # full counts of pfam

	df_merged.to_csv("g_gene_full_counts__" + NAME + ".csv")

	#convert to just 1 or 0
	df_merged_ones = df_merged >=1
	df_merged_ones = df_merged_ones.astype(int)

	print(df_merged_ones)

	#add column of the pfam abundance
	df_merged_ones["score"] = df_merged_ones.apply(sum, axis=1) / len(df_merged_ones.columns)

	#removes pfams rows with pfams that apear in most genomes (threshold)
	df_filtered_presence_absence = df_merged_ones[df_merged_ones["score"] < THRESHOLD]

	#remove pfams rows with pfams that does not apear in most genomes (1-threshold)
	min_threshold = 1-THRESHOLD
	df_filtered_presence_absence = df_filtered_presence_absence[df_filtered_presence_absence["score"] > min_threshold]

	#remove the added column
	del df_filtered_presence_absence["score"]

	df_filtered_presence_absence.to_csv("g_gene_presence_absence__" + NAME + ".csv")

	get_tree(CLADE,SAMPLE,COLLAPSE)
