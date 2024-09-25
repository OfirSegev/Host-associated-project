import pandas as pd
import subprocess
pd.set_option('display.max_columns', None)
import glob
import concurrent.futures

#parameters
base = "/sci/labs/asafle/alexlevylab/icore-data/alex/workspace/HA_vs_NHA/data_sources/final_genomes/predictions/"
run= "HvN_collapse_G_distance_0.8"
COLLAPSE="G"
LABEL="HvN"
ENRICHMENT="Host-associated"
##################

list_of_clades = pd.read_csv(base+"/"+run+"/list_of_best_clades.tsv", sep="\t", header=None)
list_of_clades = list_of_clades[0]

def combine_gff(clade_id):
    d = pd.read_csv("all_clades_gff_pfams_tables/all_gff_with_pfams_"+clade_id+".csv", usecols = ['scaffold_ID','Start','End','Strand','gene_ID','genome_ID','Pfam_id'])
    d["clade_id"] = clade_id
    return d

with concurrent.futures.ProcessPoolExecutor(max_workers=10) as executor:

        future = list(executor.map(combine_gff, list_of_clades, chunksize=1))

#add gffs
gff_raw = pd.concat(future)

corr_data2 = pd.read_csv(ENRICHMENT+"_results_"+LABEL+"_collapse_"+COLLAPSE+"_with_operon_data_filtered_for_singletons.csv")
corr_data = corr_data2[['pfam_id', 'clade_id', 'operon_id','description']].copy()
print("Operons: " + str(corr_data2['operon_id'].nunique()))

#calculate operon max size
corr_data['size'] = corr_data.groupby('operon_id')['pfam_id'].transform('nunique')

#gff with data on operons
gff = corr_data.merge(gff_raw, left_on = ['pfam_id', 'clade_id'], right_on = ['Pfam_id', 'clade_id'], how = 'left')
gff.sort_values(['genome_ID', 'scaffold_ID',"Start"], inplace = True)

print(gff)

group = ['operon_id', 'genome_ID', 'scaffold_ID']
gb = gff.groupby(group).filter(lambda x: len(x) > 1)
print(gb)

gb['start_backshifted'] = gb.groupby(group)['Start'].shift(-1)
# if the distance is negative, this means the domians are in the same protein
gb['distance'] = gb['start_backshifted'] - gb['End']
gb['distance'].fillna(0, inplace = True)

# calculate mean distance and ignore the domains that are in the same protein
gb['mean_distance'] = gb.groupby(group)['distance'].transform(lambda x: x[x >= 0].mean())

#gb['mean_distance'] = gb.groupby(group)['distance'].transform('mean')
gb['direction_homogeneity'] = gb.groupby(group)['Strand'].transform(lambda x: x.value_counts(normalize = True).max())

#calculate number of unique pfams in the operon
gb['number_of_unique_hits'] = gb.groupby(group)['pfam_id'].transform('nunique')
#completeness = unique pfams/ max
gb["completeness"] = gb["number_of_unique_hits"]/gb["size"]

#remove unnessesary columns
gb = gb.drop(["Pfam_id","number_of_unique_hits","start_backshifted"], axis=1)
#sort the table
gb = gb.sort_values(by=group)

print(gb)
gb.to_csv(ENRICHMENT+"_"+LABEL+"_collapse_"+COLLAPSE+"_operon_parameters_with_genome_data.csv", index = False)

############

#summarize each operon
x = gb.groupby('operon_id')[['mean_distance', 'direction_homogeneity']].mean().reset_index()

print(x['mean_distance'].value_counts(bins = 10).sort_index())
print(x['direction_homogeneity'].value_counts(bins = 10).sort_index())

condition1 = x['mean_distance'] < 3000 
condition2 = x['direction_homogeneity'] > 0.65 

x_filtered = x.loc[(condition1) & (condition2)]
print(x_filtered)
final = corr_data2.loc[corr_data2['operon_id'].isin(x_filtered['operon_id'].unique())]

mean_tests_passed = final.groupby('operon_id')[['tests_passed']].mean().reset_index()
mean_tests_passed = mean_tests_passed.rename(columns={"tests_passed": "mean_tests_passed"})
#final with mean_tests_passed value
final = final.merge(mean_tests_passed, on="operon_id", how="left")
#final mean_distance, direction_homogeneity (from x_filtered)
final = final.merge(x_filtered, on="operon_id", how="left")

print("after Operons: " + str(final['operon_id'].nunique()))
print(final)
final.to_csv(ENRICHMENT+"_results_"+LABEL+"_collapse_"+COLLAPSE+"_with_operon_data_filtered_with_distance.csv", index = False)