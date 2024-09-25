import pandas as pd

path="/sci/labs/asafle/alexlevylab/icore-data/alex/workspace/HA_vs_NHA/data_sources/final_genomes/updated_results_for_paper/"

label = "HvN"
collapse="G"
FL = True
enrichment = "Host-associated"

# Function to remove extra spaces inside the tuple-like string
def remove_extra_spaces(index_str):
    # Remove spaces after commas or within the parentheses
    return index_str.replace(", '  ", ", '").replace("( '", "('").replace(" ')", "')").strip()


if (FL):
    df = pd.read_csv(path+label+"_"+collapse+"/FL_pivot_table_"+enrichment+"_results_"+collapse+"_"+label+".csv", index_col=0)
    # Apply the function to the DataFrame's index
    df.index = df.index.map(remove_extra_spaces)
else:
    df = pd.read_csv(path+label+"_"+collapse+"/pivot_table_"+enrichment+"_results_"+collapse+"_"+label+".csv", index_col=0)


#remove specific clade (for HvN_G_FL_host)
#df = df.drop("Enterococcus", axis=1)
#remove specific clade (for HvN_G_FL_environment)
#df = df.drop("Pseudomonas", axis=1)
#remove specific clade (for HvN_F_pfams_host)
#df = df.drop("Desulfovibrionaceae", axis=1)
#remove specific clade (for HvN_G_pfams_host - 13 )
#df = df.drop(["Yersinia","Rhodanobacter","Corynebacterium"], axis=1)
#remove specific clade (for HvN_G_pfams_host - 13 )
df = df.drop(["Bacillus-A","Xanthomonas","Enterococcus"], axis=1)

top_rows = 100
final_number_of_clades = 12
number_of_rows = 30

top_df = df.head(top_rows)
print (top_df)

#value counts all clades scored within these lines
column_sums = (top_df**2).sum().sort_values(ascending=False)
#print(column_sums)

#filter for the top clades
top_clades = list(column_sums.head(final_number_of_clades).index)
print (top_clades)
top_df = top_df[top_clades]

#sort the table again
top_df["sum"] = top_df.sum(axis=1)
top_df = top_df.sort_values("sum", ascending=False)

#pick the best rows
top_df = top_df.head(number_of_rows)
print (top_df)
if(FL):
    top_df.to_csv("best_results_FL_"+label+"_"+collapse+"_"+enrichment+"_"+str(final_number_of_clades)+"_"+str(number_of_rows)+".csv")
else:
    top_df.to_csv("best_results_pfams_"+label+"_"+collapse+"_"+enrichment+"_"+str(final_number_of_clades)+"_"+str(number_of_rows)+".csv")