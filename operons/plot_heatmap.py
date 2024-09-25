import pandas as pd
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
from scipy.spatial.distance import squareform
import seaborn as sns
import matplotlib.pyplot as plt
import sys

#m = pd.read_csv("/sci/labs/asafle/alexlevylab/icore-data/alex/workspace/HA_vs_NHA/data_sources/final_genomes/genome_data/full_protein_scripts_and_data/uniprot_metadata/sample.tsv",sep = '\t', header = None, usecols = [0,3])

#infile = "scaf_correlation_matrix_genus_Aliivibrio.csv"
#infile = "scaffold_correlation_matrixs/scaf_correlation_matrix_clade_Paraburkholderia.csv"
infile = sys.argv[1] 
correlation_matrix = pd.read_csv(infile, index_col = 'Pfam_id')
CLADE = infile.split(".")[0].split("_")[-1] 

#descriptions to merge
descriptions = pd.read_csv("all_clades_gff_pfams_tables/all_gff_with_pfams_"+CLADE+".csv", usecols = ["Pfam_id","Description"])
descriptions = descriptions.drop_duplicates()

condensed_distance_matrix = squareform(1 - correlation_matrix)

Z = linkage(condensed_distance_matrix, 'average')
leaf_order = leaves_list(Z)

ordered_corr_matrix = correlation_matrix.iloc[leaf_order, :].iloc[:, leaf_order]

#add description column 
#mer = ordered_corr_matrix.merge(m, how = 'left', left_index = True, right_on = 0 )
mer = ordered_corr_matrix.merge(descriptions, how = 'left', left_index = True, right_on = "Pfam_id" )

mer.drop(columns = "Pfam_id", inplace = True)
mer.set_index("Description", inplace = True)
print(mer)

plt.figure(figsize=(90, 90))
sns.heatmap(mer, cmap='coolwarm', vmin=-1, vmax=1)

plt.savefig('figures/adjusted_correlation_matrix_sorted_50kwindow_clade_'+CLADE+'.png') 

