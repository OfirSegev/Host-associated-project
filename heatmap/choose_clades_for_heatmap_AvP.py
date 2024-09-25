import pandas as pd

path="/sci/labs/asafle/alexlevylab/icore-data/alex/workspace/HA_vs_NHA/data_sources/final_genomes/updated_results_for_paper/"

label = "AvP"
collapse="G"
FL = True
#enrichment = "Animal-associated"

# Function to remove extra spaces inside the tuple-like string
def remove_extra_spaces(index_str):
    # Remove spaces after commas or within the parentheses
    return index_str.replace(", '  ", ", '").replace("( '", "('").replace(" ')", "')").strip()

if (FL):
    df_animals = pd.read_csv(path+label+"_"+collapse+"/FL_pivot_table_Animal-associated_results_"+collapse+"_"+label+".csv", index_col=0)
    df_plants = pd.read_csv(path+label+"_"+collapse+"/FL_pivot_table_Plant-associated_results_"+collapse+"_"+label+".csv", index_col=0)
    #make plants negative values
    df_plants = df_plants.applymap(lambda x: -x if x != 0 else x)
    #combine dataframes
    df = df_plants.add(df_animals, fill_value=0)
    # Apply the function to the DataFrame's index
    df.index = df.index.map(remove_extra_spaces)
else:
    df_animals = pd.read_csv(path+label+"_"+collapse+"/pivot_table_Animal-associated_results_"+collapse+"_"+label+".csv", index_col=0)
    df_plants = pd.read_csv(path+label+"_"+collapse+"/pivot_table_Plant-associated_results_"+collapse+"_"+label+".csv", index_col=0)
    #make plants negative values
    df_plants = df_plants.applymap(lambda x: -x if x != 0 else x)
    #combine dataframes
    df = df_plants.add(df_animals, fill_value=0)


top_rows = 100
final_number_of_clades = 10
number_of_rows = 30

#drop clade (AvP_F_pfams)
#df = df.drop(["NBRC-103111","Acetobacteraceae"], axis=1)

#sort new dataframe
df["sum"] = (df**2).sum(axis=1)
df = df.sort_values("sum", ascending=False)

top_df = df.head(top_rows)
print (top_df)

#value counts all clades scored within these lines
column_sums = (top_df**2).sum().sort_values(ascending=False).drop("sum")
#print(column_sums)

#filter for the top clades
top_clades = list(column_sums.head(final_number_of_clades).index)
print (top_clades)
top_df = top_df[top_clades]

#sort the table again
#?used for AvP_F_Pfams
#top_df["sum"] = (top_df**2).sum(axis=1)
#?used for AvP_G_Pfams and AvP_G_AFC
top_df["sum"] = (top_df).sum(axis=1).abs()
#?used for AvP_F_AFC
#top_df["sum"] = (top_df**3).sum(axis=1).abs()
top_df = top_df.sort_values("sum", ascending=False)

#pick the best rows
top_df = top_df.head(number_of_rows)
print (top_df)
if(FL):
    top_df.to_csv("best_results_FL_"+label+"_"+collapse+"_both_"+str(final_number_of_clades)+"_"+str(number_of_rows)+".csv")
else:
    top_df.to_csv("best_results_pfams_"+label+"_"+collapse+"_both_"+str(final_number_of_clades)+"_"+str(number_of_rows)+".csv")