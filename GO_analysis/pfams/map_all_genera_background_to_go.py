import pandas as pd
import networkx
import obonet
import sys
import concurrent.futures

path = "/sci/labs/asafle/alexlevylab/icore-data/alex/workspace/HA_vs_NHA/data_sources/final_genomes/predictions"
folder_name = sys.argv[1]
LABEL = sys.argv[2]
COLLAPSE = sys.argv[3]

df = pd.read_csv("pfam_interpro_go_table.tsv", sep="\t")
df_reduced = df.dropna(subset=["GO_id"])

graph = obonet.read_obo("go-basic.obo")

#mapping
id_to_name = {id_: data.get('name') for id_, data in graph.nodes(data=True)}
name_to_id = {data['name']: id_ for id_, data in graph.nodes(data=True) if 'name' in data}

def map_to_go(clade_id):
    #get a set of the pfams in the 3 samples
    sample1 = pd.read_csv(path+"/"+folder_name+"/files_for_scoary/clade_"+clade_id+"_sample_1/g_gene_presence_absence__collapse_"+COLLAPSE+"__clade_"+clade_id+"_sample_1.csv", usecols=[0], names=["pfam_id"]).dropna()
    sample2 = pd.read_csv(path+"/"+folder_name+"/files_for_scoary/clade_"+clade_id+"_sample_2/g_gene_presence_absence__collapse_"+COLLAPSE+"__clade_"+clade_id+"_sample_2.csv", usecols=[0], names=["pfam_id"]).dropna()
    sample3 = pd.read_csv(path+"/"+folder_name+"/files_for_scoary/clade_"+clade_id+"_sample_3/g_gene_presence_absence__collapse_"+COLLAPSE+"__clade_"+clade_id+"_sample_3.csv", usecols=[0], names=["pfam_id"]).dropna()
    pfams_in_genus = set(sample1["pfam_id"]).union(sample2["pfam_id"], sample3["pfam_id"])

    #filter for only genus with pfams
    clade_df = df_reduced[df_reduced["pfam_id"].isin(pfams_in_genus)]
    # Your empty DataFrame
    output = pd.DataFrame(columns=['Clade_id','GO_term', 'Occurrences', 'GO_category'])
    
    for GO_term in clade_df["GO_id"]:
        GO_list = sorted(id_to_name[superterm] for superterm in networkx.descendants(graph, GO_term))
        
        #add to the list the current go term
        GO_list.append(id_to_name[GO_term])

        #GO category column
        if("molecular_function" in GO_list):
            category = "molecular_function"
        elif("biological_process" in GO_list):
            category = "biological_process"
        else:
            category = "cellular_component"

        # Update existing counts or append new rows
        for string in GO_list:
            if string in output['GO_term'].values:
                output.loc[output['GO_term'] == string, 'Occurrences'] += 1                 
            else:
                output = pd.concat([output, pd.DataFrame([{'Clade_id': clade_id, 'GO_term': string, 'Occurrences': 1, 'GO_category': category}])], ignore_index=True)

    # Display the updated DataFrame
    print(output)

    output = output.sort_values("Occurrences", ascending = False)

    return output

  
#################
clade_list = pd.read_csv(path+"/"+folder_name+"/list_of_best_clades.tsv", sep="\t", header=None)[0].to_list()


with concurrent.futures.ProcessPoolExecutor(max_workers=50) as executor:

        future = list(executor.map(map_to_go, clade_list, chunksize=1))


all_outputs = pd.concat(future)

all_outputs.to_csv("All/all_clades_background_go_terms_"+LABEL+"_"+COLLAPSE+".tsv", index=False, sep="\t")