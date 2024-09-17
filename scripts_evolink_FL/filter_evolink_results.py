import pandas as pd
import sys

COLLAPSE=sys.argv[1]
LABEL=sys.argv[2]
sample_threshold=int(sys.argv[3])

def merge_samples(df):
    df["sample_counts"] = len(df)
    df["avg_Prevalence_index"] = df["Prevalence_index"].mean()
    df["avg_Evolink_index"] = df["Evolink_index"].mean()
    df["avg_scores"] = df["scores"].mean()
    output = df.iloc[0]
    return (output[["sample_counts","avg_Prevalence_index","avg_Evolink_index","avg_scores"]])

def filter_results(df, type):
    #split enriched/depleted
    if (type=="enriched"):
        df = df[df["Evolink_index"] > 0]
    elif (type=="depleted"):
        df = df[df["Evolink_index"] < 0]

    #get significant results only
    df = df[df["significance"] == "sig"]

    #group df by samples
    grouped = df.groupby(["Pfam_id","Clade_id"]).apply(merge_samples)
    grouped = grouped.reset_index()
    print(type)
    print (grouped)

    #filtered by sample counts
    if(len(grouped) == 0):
        return (pd.DataFrame(columns = ["pfam_id","clade_id","sample_counts","avg_Prevalence_index","avg_Evolink_index","avg_scores"]))
    filtered_grouped = grouped[grouped["sample_counts"] >= sample_threshold]

    #change column names
    filtered_grouped.rename(columns = {'Pfam_id':'pfam_id','Clade_id':'clade_id'}, inplace=True)
    
    return (filtered_grouped)

df_P_A = pd.read_csv("evolink_results_FL/"+LABEL+"_collapse_"+COLLAPSE+"__evolink_presence_absence_tables_combined.csv")
df_counts = pd.read_csv("evolink_results_FL/"+LABEL+"_collapse_"+COLLAPSE+"__evolink_counts_tables_combined.csv")

#filter results
df_P_A_enriched = filter_results(df_P_A, "enriched")
df_P_A_depleted = filter_results(df_P_A, "depleted")
df_counts_enriched = filter_results(df_counts, "enriched")
df_counts_depleted = filter_results(df_counts, "depleted")

#to csv outputs
df_P_A_enriched.to_csv("evolink_results_FL/filtered_enriched_"+LABEL+"_collapse_"+COLLAPSE+"__evolink_presence_absence_tables_combined.csv", index=False)
df_P_A_depleted.to_csv("evolink_results_FL/filtered_depleted_"+LABEL+"_collapse_"+COLLAPSE+"__evolink_presence_absence_tables_combined.csv", index=False)
df_counts_enriched.to_csv("evolink_results_FL/filtered_enriched_"+LABEL+"_collapse_"+COLLAPSE+"__evolink_counts_tables_combined.csv", index=False)
df_counts_depleted.to_csv("evolink_results_FL/filtered_depleted_"+LABEL+"_collapse_"+COLLAPSE+"__evolink_counts_tables_combined.csv", index=False)
