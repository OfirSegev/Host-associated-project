import pandas as pd

#parameters
base = "/sci/labs/asafle/alexlevylab/icore-data/alex/workspace/HA_vs_NHA/data_sources/final_genomes/predictions/"
run= "HvN_collapse_G_distance_0.8"
COLLAPSE="G"
LABEL="HvN"


def parse_mcl_clusters(clade_id):

    # Read the TSV file into a DataFrame
    data = pd.read_csv("mcl_inputs/scaf_correlation_clade_"+clade_id+"_edges_above_0.1.tsv", sep='\t')

    #remove negative correlations
    data = data[data["Correlation"] > 0]

    # Create a dictionary to store pairwise correlations
    correlation_dict = {}
    for index, row in data.iterrows():
        pfam1, pfam2, correlation = row['Gene1'], row['Gene2'], row['Correlation']
        correlation_dict.setdefault(pfam1, {})[pfam2] = correlation
        correlation_dict.setdefault(pfam2, {})[pfam1] = correlation

    # Initialize the clusters list
    clusters = []

    # Function to check if a new member can join the cluster
    def can_join_cluster(cluster, new_member, correlation_dict, threshold=0.0):
        for member in cluster:
            if correlation_dict.get(member, {}).get(new_member, 0) <= threshold:
                return False
        return True

    # Iterate through each pair to form clusters
    for index, row in data.iterrows():
        pfam1, pfam2 = row['Gene1'], row['Gene2']
        
        # Check if the pair can join any existing cluster
        added_to_cluster = False
        for cluster in clusters:
            if can_join_cluster(cluster, pfam1, correlation_dict) and can_join_cluster(cluster, pfam2, correlation_dict):
                cluster.add(pfam1)
                cluster.add(pfam2)
                added_to_cluster = True
                break
        
        # If the pair cannot join any existing cluster, create a new cluster
        if not added_to_cluster:
            clusters.append(set([pfam1, pfam2]))

    # Create a DataFrame for the clusters with numerical cluster numbers
    cluster_list = []
    for i, cluster in enumerate(clusters):
        for pfam_id in cluster:
            cluster_list.append({'Pfam_id': pfam_id, 'Cluster': i})

    cluster_df = pd.DataFrame(cluster_list)

    #add description
    output = cluster_df.merge(df_description, how="left")
    print (output)


    output.to_csv("results_with_clusters_ids/"+clade_id+"_results_with_clusters_ids.tsv", sep="\t", index=False)


##################

#read pfam descreption file
df_description = pd.read_csv("/sci/labs/asafle/alexlevylab/icore-data/alex/workspace/HA_vs_NHA/misc_metadata/pfamID_to_description.tsv", sep="\t")
df_description["pfam_id"] = df_description["pfam_id"].str.split(".").str[0]
df_description.columns = ["Pfam_id","Description"]

list_of_clades = pd.read_csv(base+"/"+run+"/list_of_best_clades.tsv", sep="\t", header=None)

for clade_id in list_of_clades[0]:
    print (clade_id)
    parse_mcl_clusters(clade_id)