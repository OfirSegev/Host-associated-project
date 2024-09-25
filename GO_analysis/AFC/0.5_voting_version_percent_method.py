import pandas as pd
pd.set_option("display.max_columns", None)

df = pd.read_csv("Total_AFC_to_GO_based_on_uniprot.csv")
#df = pd.read_csv("../../../../genome_data/full_protein_scripts_and_data/uniprot_metadata/Total_AFC_to_cog_based_on_uniprot_wCategories.csv")
#df = df.loc[~df['list_version'].isna()]

print(df)
# calculate Z score of GO appearences

o = df.groupby("repID")['GO_term'].value_counts().reset_index()
p = df.groupby("repID")['memID'].nunique().reset_index().rename(columns = {"memID" : "memID_count"})

mer = o.merge(p, on = 'repID', how = 'left')

mer['percent_presence'] = mer['count'] / mer['memID_count'] 

print(mer)
print(mer['percent_presence'].value_counts())
##########

THRESHOLD = 0.75 

telem = mer.loc[(mer['percent_presence'] >= THRESHOLD)]

print(telem)

#out = telem[["repID", "GO_term"]].merge(df[['cog', 'list_version', '2']].drop_duplicates(), on = 'GO_term',  how = 'left')
#print(out)

telem.to_csv("AFC_repID_to_GO_master_sheet__memIDCluster_within_PercentPresence_gte_" + str(THRESHOLD)+".csv", index = False)
#out.to_csv("AFC_repID_to_GO_master_sheet__memIDCluster_within_PercentPresence_gte_" + str(THRESHOLD)+".csv", index = False)
