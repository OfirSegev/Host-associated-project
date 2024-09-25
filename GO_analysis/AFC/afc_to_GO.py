import pandas as pd
#from matplotlib_venn import venn2
import matplotlib.pyplot as plt

#filename = "/sci/labs/asafle/alexlevylab/icore-data/alex/workspace/HA_vs_NHA/data_sources/final_genomes/genome_data/full_protein_scripts_and_data/uniprot_metadata/AFC_memberID_Descriptions_and_CrossRefs_raw_eggnogRaw.csv"
filename = "AFC_memberID_Descriptions_and_CrossRefs_only_GO.csv"

df = pd.read_csv(filename)
print(df)

df['cr_list'] = df['Cross References'].str.split(";")

# WE DROP NA HERE -- those that have no crossrefs
df.dropna(subset = "cr_list", inplace = True)

df['GO_term'] = df['cr_list'].apply(lambda x: [i for i in x if "GO!" in i])
#print(df)
df2 = df.explode('GO_term')
df2['GO_term'] = df2['GO_term'].str.replace("GO!", "")

df2.dropna(subset = 'GO_term', inplace = True)

print(df2)

m = pd.read_csv("/sci/labs/asafle/alexlevylab/icore-data/alex/workspace/HA_vs_NHA/data_sources/final_genomes/genome_data/new_full_protein_data/mmseqs_to_AFC_master_sheet.csv", usecols = ['repID', 'memID'])
print(m)

x = m.merge(df2, how = 'left', on = 'memID')

x.drop(columns = ['Description', 'Cross References', 'cr_list'], inplace = True)
x.drop_duplicates(inplace = True)
x.dropna(subset = 'GO_term', inplace = True)

print(x)

x.to_csv("Total_AFC_to_GO_based_on_uniprot.csv", index = False)
x.drop(columns = 'memID', inplace = True)

x.drop_duplicates(inplace = True)
x.to_csv("repID_to_GO_based_on_uniprot.csv", index = False)
