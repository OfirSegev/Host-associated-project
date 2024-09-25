import pandas as pd
import sys

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

all_levels = pd.read_csv("go_terms_all_levels.tsv", sep="\t")

#filters the dataframe contains all results of all clades
def filter_df(min_clades):
	df = pd.read_csv("results/all_clades/"+label+"_"+enrichment+"_larger_than_"+number_of_tests_threshold+"_hypergeometric_distribution_all_clades_collapse_"+collapse+".tsv", sep="\t")

	#filter by q-value
	df = df[df["fdr_corrected"] <= 0.05]

	# Group by 'GO_term' and calculate the mean for the grouped columns
	grouped_df = df.groupby('GO_term').agg({'k': 'mean', 'n': 'mean', 'N': 'mean', 'M': 'mean', 'p_value': 'mean', 'fdr_corrected': 'mean'}).reset_index()

	# Add a new column 'clades_count' that contains the count of grouped rows for each 'GO_term'
	grouped_df['clades_count'] = df.groupby('GO_term').size().values

	# include the 'GO_category' in the aggregation 
	grouped_df['GO_category'] = df.groupby('GO_term')['GO_category'].first().values

	#keep terms that apear in multiple clades
	filtered_df = grouped_df[grouped_df["clades_count"] > min_clades]

	return filtered_df

#returns True if the go term adds new information (the higher term has more in k than the childs)
def is_informative(data, GO_term):
	k = data[data["GO_term"] == GO_term].iloc[0]["k"]
	fdr = data[data["GO_term"] == GO_term].iloc[0]["fdr_corrected"]

	for child in data[data["parent"] == GO_term]["GO_term"]:
		#if child has the same amount of results as parent and if child has better p-value than parent
		if ((data[data["GO_term"] == child].iloc[0]["k"] == k) | (data[data["GO_term"] == child].iloc[0]["fdr_corrected"] < fdr)):
			return False
	return True


def best_to_plot_all(category):
	df = filter_df(0)
	df = df[df["GO_category"] == category]
	df = df.drop(["GO_category"], axis=1)

	final_set = set()

	all_levels_filtered = all_levels[all_levels["GO_category"] == category]
	
	mer = df.merge(all_levels_filtered, how="left", on="GO_term")

	for GO_term in set(mer["GO_term"]):

		if is_informative(mer,GO_term):
			final_set.add(GO_term)

	output = mer[mer["GO_term"].isin(final_set)]
#	output.to_csv("results/to_plot/"+label+"_"+enrichment+"_larger_than_"+number_of_tests_threshold+"_"+category+"_best_candidates_to_plot_all_clades_collapse_"+collapse+".tsv", sep="\t", index=False)
	return output

def run_best_to_plot_all_clades():

	output1 = best_to_plot_all("molecular_function")
	output2 = best_to_plot_all("biological_process")
	output3 = best_to_plot_all("cellular_component")
	merged_output = pd.concat([output1,output2,output3])
#	merged_output.to_csv("results/to_plot/"+label+"_"+enrichment+"_larger_than_"+number_of_tests_threshold+"_combined_best_candidates_to_plot_all_clades_collapse_"+collapse+".tsv", sep="\t", index=False)

def best_to_plot_per_clade(df,clade_id,category):

	df = df[(df["GO_category"] == category) & (df["Clade_id"] == clade_id)]
	df = df[df["fdr_corrected"] <= 0.05]
	df = df.drop(["GO_category"], axis=1)

	final_set = set()

	all_levels_filtered = all_levels[all_levels["GO_category"] == category]
	
	mer = df.merge(all_levels_filtered, how="left", on="GO_term")

	for GO_term in set(mer["GO_term"]):

		if is_informative(mer,GO_term):
			final_set.add(GO_term)

	output = mer[mer["GO_term"].isin(final_set)]
	#output.to_csv("results/to_plot/"+label+"_"+enrichment+"_larger_than_"+number_of_tests_threshold+"_"+category+"_best_candidates_to_plot_all_clades.tsv", sep="\t", index=False)
	return output

def run_best_to_plot_per_clade():
	df = pd.read_csv("results/all_clades/"+label+"_"+enrichment+"_larger_than_"+number_of_tests_threshold+"_hypergeometric_distribution_all_clades_collapse_"+collapse+".tsv", sep="\t")
	clade_list = set(df["Clade_id"])

	all_output = pd.DataFrame(columns=["Clade_id","GO_term","k","n","N","M","p_value","fdr_corrected","level","parent","GO_category"])

	for clade_id in clade_list:
		output1 = best_to_plot_per_clade(df,clade_id,"molecular_function")
		output2 = best_to_plot_per_clade(df,clade_id,"biological_process")
		output3 = best_to_plot_per_clade(df,clade_id,"cellular_component")
		merged_output = pd.concat([output1,output2,output3])

		all_output = pd.concat([all_output, merged_output])
	
	all_output = all_output.drop_duplicates(subset=["Clade_id", "GO_term"], keep='first')

	#add odds ratio
	if (table_type == "enriched"):
		all_output["odds_ratio"] = (all_output["k"]/all_output["N"]) / (all_output["n"]/all_output["M"])
	elif(table_type == "depleted"):
		all_output["odds_ratio"] = (all_output["n"]/all_output["M"]) / (all_output["k"]/all_output["N"])

	all_output.to_csv("results/to_plot/"+label+"_"+enrichment+"_larger_than_"+number_of_tests_threshold+"_combined_best_candidates_to_plot_per_clade_collapse_"+collapse+"_with_odds_ratio.tsv", sep="\t", index=False)


#run all clades (mean of all columns)
run_best_to_plot_all_clades()

#seperated for each clade
run_best_to_plot_per_clade()
