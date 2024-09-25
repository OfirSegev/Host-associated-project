import pandas as pd
import sys

#infile = "scaf_correlation_matrix_genus_Gottfriedia.csv"
infile = sys.argv[1]
df = pd.read_csv(infile, index_col = 'Pfam_id')
CLADE = infile.split(".")[0].split("_")[-1] 

THRESHOLD = 0.1

edges = []
for i in range(df.shape[0]):
#    for j in range(i + 1, df.shape[1]):  
    for j in range(i, df.shape[1]):  
        if df.iloc[i, j] > THRESHOLD:  
            edges.append((df.index[i], df.columns[j], df.iloc[i, j]))

edges_df = pd.DataFrame(edges, columns=['Gene1', 'Gene2', 'Correlation'])

edges_df.to_csv('mcl_inputs/scaf_correlation_clade_'+CLADE+'_edges_above_'+ str(THRESHOLD)+ '.tsv', sep='\t', index=False)
