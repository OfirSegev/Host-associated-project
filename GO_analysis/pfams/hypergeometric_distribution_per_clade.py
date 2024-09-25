import pandas as pd
from scipy.stats import hypergeom
import sys
import statsmodels.stats.multitest as smm
import concurrent.futures

# k = success in sample
# n = success in population
# N = sample size
# M = population size

label = sys.argv[1]
collapse = sys.argv[2]
table_type = sys.argv[3] # enriched/depleted
number_of_tests_threshold = sys.argv[4] if len(sys.argv) > 4 else "2"

if (label == "HvN"):
	if (table_type == "enriched"):
		enrichment="Host-associated"
	else:
		enrichment="Environmental"
elif (label == "AvP"):
	if (table_type == "enriched"):
		enrichment="Animal-associated"
	else:
		enrichment="Plant-associated"
		

df_test = pd.read_csv(label+"/"+label+"_larger_than_"+number_of_tests_threshold+"_filtered_"+enrichment+"_all_go_terms_per_clade_collapse_"+collapse+".tsv", sep="\t")
df_all = pd.read_csv("All/all_clades_background_go_terms_"+label+"_"+collapse+".tsv", sep="\t")

clade_list = set(df_test["Clade_id"])

def hypergeom_test(row):
    return hypergeom.sf(k=(row["k"]-1), M=row["M"], n=row["n"], N=row["N"], loc=0)

def calculate_hypergeom(clade_id,category):

    df1 = df_test[(df_test["GO_category"] == category) & (df_test["Clade_id"] == clade_id)]
    df1 = df1.drop(["Clade_id","GO_category"], axis=1)
	
    df2 = df_all[(df_all["GO_category"] == category) & (df_all["Clade_id"] == clade_id)]
    df2 = df2.drop(["Clade_id","GO_category"], axis=1)
	
    #when there are no results in this category for this clade
    if((len(df1) == 0) | (len(df2) == 0)):
        return pd.DataFrame(columns=['Clade_id','GO_term','k','n','N','M','p_value','fdr_corrected', 'GO_category'])

    mer = df1.merge(df2, on = "GO_term", how="left")
    mer.rename(columns={"Occurrences_x": "k", "Occurrences_y": "n"}, inplace=True)

    mer["N"] = df1[df1["GO_term"] == category].iloc[0]["Occurrences"]
    mer["M"] = df2[df2["GO_term"] == category].iloc[0]["Occurrences"]

    mer["p_value"] = mer.apply(hypergeom_test, axis=1)
	
    # Add FDR corrected p-values to DataFrame as a new column
    fdrs = smm.fdrcorrection(mer["p_value"], alpha=0.05, method='indep', is_sorted=False)
    mer["fdr_corrected"] = fdrs[1]
	
    #add clade id and category columns
    mer.insert(loc=0, column='Clade_id', value=clade_id)
    mer["GO_category"] = category

    mer = mer.sort_values("p_value")
    print (mer)
	
    return mer

    #mer.to_csv("results/"+category+"/"+label+"_"+enrichment+"_larger_than_"+number_of_tests_threshold+"_"+category+"_hypergeometric_distribution.tsv", index=False, sep="\t")
	
def run_per_clade(clade_id):
	
    mf = calculate_hypergeom(clade_id,"molecular_function")
    bp = calculate_hypergeom(clade_id,"biological_process")
    cc = calculate_hypergeom(clade_id,"cellular_component")
	
    combined = pd.concat([mf,bp,cc])
    return combined


with concurrent.futures.ProcessPoolExecutor(max_workers=10) as executor:

        future = list(executor.map(run_per_clade, clade_list, chunksize=1))

all_outputs = pd.concat(future)
print(all_outputs)
all_outputs.to_csv("results/all_clades/"+label+"_"+enrichment+"_larger_than_"+number_of_tests_threshold+"_hypergeometric_distribution_all_clades_collapse_"+collapse+".tsv", index=False, sep="\t")

    #hypergeom.pmf(k, M, n, N, loc=0)
