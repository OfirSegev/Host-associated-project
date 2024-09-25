import pandas as pd
import numpy as np
import concurrent.futures

#parameters
base = "/sci/labs/asafle/alexlevylab/icore-data/alex/workspace/HA_vs_NHA/data_sources/final_genomes/predictions/"
run= "HvN_collapse_G_distance_0.8"
COLLAPSE="G"
LABEL="HvN"
ENRICHMENT="Host-associated"
TESTS_PASSED = 1

list_of_clades = pd.read_csv(base+"/"+run+"/list_of_best_clades.tsv", sep="\t", header=None)
list_of_clades = list_of_clades[0]

#add columns to each dataframe
def update_dfs(clade_id):

    print (clade_id)
    df = pd.read_csv("results_with_clusters_ids/"+clade_id+"_results_with_clusters_ids.tsv", sep = '\t')#, usecols = ['ID', "Cluster"])

    df['operon_id'] = clade_id + '_operon_' + df['Cluster'].astype(str)
    df['clade_id'] = clade_id
    return df


with concurrent.futures.ProcessPoolExecutor(max_workers=10) as executor:

        future = list(executor.map(update_dfs, list_of_clades, chunksize=1))


c = pd.concat(future)

##raw operons##
# this includes cross-genus hits, i.e. hits in X but also analyzed in Y
#a = c.groupby("operon_id").apply(lambda x: len(x) == 1).reset_index()
#d = c.merge(a, on = 'operon_id', how = 'left')
#d.loc[d[0] == True, 'operon_id'] = np.nan 
#print(d)
#d[['ID', 'operon_id', 'descriptions_collection']].to_csv("raw_operons_including_enriched_in_other_genera.csv", index = False)
###

print(c)
########

hits = pd.read_csv(base+"/"+run+"/tables_to_export/combined_results_"+LABEL+"_collapse_"+COLLAPSE+".csv")

hits = hits.loc[hits['enrichment'] == ENRICHMENT]

hits2 = hits.loc[hits['tests_passed'] >= TESTS_PASSED]

print (hits2)

#########

mer = hits2.merge(c, how = 'left', left_on = ['clade_id', 'pfam_id'], right_on = ['clade_id', 'Pfam_id'])

print (mer)

#remove singletons - not really "operons", they are singleton clusters
a = mer.groupby("operon_id").apply(lambda x: len(x) == 1).reset_index()
print(a[0].value_counts())

mer = mer.merge(a, on = 'operon_id', how = 'left')
mer.loc[mer[0] == True, 'operon_id'] = np.nan 

mer.drop(columns = ['Pfam_id', 'Cluster', 'Description'], inplace = True)
mer.sort_values(['clade_id', 'operon_id'], inplace = True)



new_order = [
    'pfam_id', 'clade_id', "operon_id", 'evolink_P_A', 'evolink_counts','fisher', 'lmm', 'scoary', 
    'tests_passed', 'enrichment', 'description', 
]

mer = mer[new_order]

print(mer)

mer.to_csv(ENRICHMENT+"_results_"+LABEL+"_collapse_"+COLLAPSE+"_with_operon_data.csv", index = False)

#remove the singletons
mer = mer.loc[~mer['operon_id'].isna()]

mer.to_csv(ENRICHMENT+"_results_"+LABEL+"_collapse_"+COLLAPSE+"_with_operon_data_filtered_for_singletons.csv", index = False)
