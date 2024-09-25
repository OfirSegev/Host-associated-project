import pandas as pd
import obonet

graph_with_obs = obonet.read_obo("go-basic.obo", ignore_obsolete=False)

#########################
#! remove and update obsolete

# Extract information for obsolete terms
obsolete_terms = []
for node, data in graph_with_obs.nodes(data=True):
    if data.get('is_obsolete') == 'true':
        term_info = {
            'id': node,
            'name': data.get('name'),
            'replaced_by': data.get('replaced_by', [None])[0]  # get the first replacement if it exists
        }
        obsolete_terms.append(term_info)

# Create a DataFrame
df_obsolete = pd.DataFrame(obsolete_terms, columns=['id', 'name', 'replaced_by'])

# Display the DataFrame
print(df_obsolete)
#df_obsolete.to_csv("obsolete_test.csv", index=False)

#read df to update
df = pd.read_csv("AFC_repID_to_GO_master_sheet__memIDCluster_within_PercentPresence_gte_0.75.csv")

# Merge the DataFrames to get the replacement information
merged_df = df.merge(df_obsolete, left_on='GO_term', right_on='id', how='left')

# Replace obsolete GO terms with new ones
merged_df['GO_term'] = merged_df.apply(lambda row: row['replaced_by'] if pd.notnull(row['replaced_by']) else row['GO_term'], axis=1)

# Remove rows where replaced_by is None
updated_df = merged_df[pd.notnull(merged_df['replaced_by']) | merged_df['id'].isnull()]

# Drop the columns from df_obsolete as they're no longer needed
updated_df = updated_df.drop(columns=['id', 'name', 'replaced_by'])

# Display the updated DataFrame
print(updated_df)

##################
#! remove and update alternative ids
alt_ids_to_change = ["GO:0016437","GO:0005887","GO:0043784","GO:0016021","GO:0004367","GO:0052928","GO:0034290"]
ids_to_change_to = ["GO:0004810","GO:0005886","GO:0008817","GO:0016020","GO:0047952","GO:0052927","GO:0140911"]

# Create a dictionary for the additional mapping
mapping_dict = dict(zip(alt_ids_to_change, ids_to_change_to))

# Replace the GO_term in updated_df based on the mapping_dict
updated_df['GO_term'] = updated_df['GO_term'].replace(mapping_dict)

# Display the updated DataFrame
print(updated_df)

#####################
updated_df.to_csv("AFC_repID_to_GO_master_sheet__memIDCluster_within_PercentPresence_gte_0.75_updated_obsolete_alt.csv", index=False)

