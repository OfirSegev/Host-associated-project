import pandas as pd
from statsmodels.stats.multitest import fdrcorrection
import os
import sys

COLLAPSE=sys.argv[1]
LABEL=sys.argv[2]

table_list = os.listdir('scoary_results/raw_combined_samples')

first=True

for table_name in table_list:

	data=pd.read_csv("scoary_results/raw_combined_samples/"+table_name)
	#take only full sampled tables
	data=data[data["Bootstrap"]==0]

	if(first):
		df=data.copy()
		first=False
	else:
		df=pd.concat([df,data],ignore_index=True)


#TODO add correction for benjaminiv chi test
def update_table(df):

        naive_p=df["Naive_p"]
        best_pvals=df["Best_pairwise_comp_p"]
        worst_pvals=df["Worst_pairwise_comp_p"]

        corrected_array_naive = fdrcorrection(naive_p, alpha=0.05, method='indep', is_sorted=False)
        corrected_array_best = fdrcorrection(best_pvals, alpha=0.05, method='indep', is_sorted=False)
        corrected_array_worst = fdrcorrection(worst_pvals, alpha=0.05, method='indep', is_sorted=False)

        df["Naive_pval_corrected"]=corrected_array_naive[1]
        df["Best_pval_corrected"]=corrected_array_best[1]
        df["Worst_pval_corrected"]=corrected_array_worst[1]

        print (df)
        return df


output=update_table(df)

output.to_csv("all_clade_sample_collapse_"+COLLAPSE+"_scoary_output_"+LABEL+"_fdr_updated.csv", index=False)

