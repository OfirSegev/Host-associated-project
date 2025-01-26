import pandas as pd
import numpy as np
import glob
from statsmodels.stats.multitest import fdrcorrection
import sys

#inputs 1- collapse  ,2 - label , 3- type of pyseer analysis
COLLAPSE=sys.argv[1]
LABEL=sys.argv[2]
TYPE=sys.argv[3]
folder_location="final_results_FL/"

####import and organize files
files = glob.glob(TYPE+"_results_FL/*")

#print(files)

outputs_list = []

for file in files:
	collapse = str(file.split("_")[4])
	clade = str(file.split("_")[7])
	sample = int(file.split("_")[9].split(".")[0])

	try:
		df = pd.read_csv(file, sep = '\t')

		arr = np.hstack((np.full((len(df), 1), collapse), np.full((len(df), 1), clade), np.full((len(df), 1), sample)))

		ids = pd.DataFrame(columns = ['collapse', 'clade', 'sample'], data = arr, index = range(0, len(df)))

#	fdr correction per sample
####BENJAMINI HOCHBERG
        
		df['q_value_chisq'] = fdrcorrection(df['filter-pvalue'].to_list(), alpha = 0.05)[1]
		df['q_value_lrt'] = fdrcorrection(df['lrt-pvalue'].to_list(), alpha = 0.05)[1]
		df['odds_ratio'] =  np.exp(df['beta'])

		together = pd.concat([ids, df], axis = 1)

		outputs_list.append(together)

	except:
		print("empty file")
		print("collapse", "clade", "sample")
		print(collapse, clade, sample)

total = pd.concat(outputs_list)

print (total)

#raw results with added fdr and odds ratio
total.to_csv(folder_location+TYPE+"_total_results_fixed_raw.csv", index = False)

#keep only relevant columns
reduced=total[['variant','clade','sample','odds_ratio','q_value_chisq','q_value_lrt']]
#change columns names
reduced.columns=["Pfam_id","Clade_id","Sample_num","Odds_ratio","q_value_chisq","q_value_lrt"]

#reduced.to_csv(folder_location+"raw_all_clade_sample_collapse_"+COLLAPSE+"_"+TYPE+"_output_"+LABEL+"_fdr_updated.csv", index=False)
reduced.to_csv("all_clade_sample_collapse_"+COLLAPSE+"_"+TYPE+"_output_"+LABEL+"_fdr_updated_FL.csv", index=False)