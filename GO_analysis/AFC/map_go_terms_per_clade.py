import pandas as pd
import networkx
import obonet
import sys
import concurrent.futures

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



df = pd.read_csv(label+"/"+label+"_"+enrichment+"_results_per_clade_collapse_"+collapse+"_with_GO_terms.tsv", sep="\t")
df_reduced = df[df["tests_passed"] >= float(number_of_tests_threshold)]

df_reduced = df_reduced.dropna(subset=["GO_term"])

#all unique clades
clade_list = set(df_reduced["clade_id"])
#clade_list = ["Achromobacter"]

graph = obonet.read_obo("go-basic.obo")

#mapping
id_to_name = {id_: data.get('name') for id_, data in graph.nodes(data=True)}
name_to_id = {data['name']: id_ for id_, data in graph.nodes(data=True) if 'name' in data}

def map_clade(clade_id):

    df_clade = df_reduced[df_reduced["clade_id"] == clade_id]
    print (df_clade)

    # Your empty DataFrame
    output = pd.DataFrame(columns=['Clade_id','GO_term', 'Occurrences', 'GO_category'])

    for GO_term in df_clade["GO_term"]:
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
        
    return output


with concurrent.futures.ProcessPoolExecutor(max_workers=10) as executor:

        future = list(executor.map(map_clade, clade_list, chunksize=1))

# Display the updated DataFrame
all_outputs = pd.concat(future)
print(all_outputs)

#all_outputs = all_outputs.sort_values("Occurrences", ascending = False)

all_outputs.to_csv(label+"/"+label+"_larger_than_"+number_of_tests_threshold+"_filtered_"+enrichment+"_all_go_terms_per_clade_collapse_"+collapse+".tsv", index=False, sep="\t")

'''
# Find edges to parent terms
node = name_to_id['pilus part']
for child, parent, key in graph.out_edges(node, keys=True):
    print(f'• {id_to_name[child]} ⟶ {key} ⟶ {id_to_name[parent]}')
'''
