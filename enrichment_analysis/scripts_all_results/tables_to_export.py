import pandas as pd
import subprocess
import sys

COLLAPSE = sys.argv[1]
LABEL = sys.argv[2]

if (LABEL=="HvN"):
    enriched_phrase = "Host-associated"
    depleted_phrase = "Environmental"
elif (LABEL=="AvP"):
    enriched_phrase = "Animal-associated"
    depleted_phrase = "Plant-associated"
elif (LABEL=="RvS"):
    enriched_phrase = "Roots"
    depleted_phrase = "Shoots"
else:
    enriched_phrase = "Enriched"
    depleted_phrase = "Depleted"


subprocess.call(['mkdir','tables_to_export'])

# combine enriched, depleted table
df_enriched = pd.read_csv("final_results/final_table_all_tests_enriched_results_pfams_"+LABEL+"_collapse_"+COLLAPSE+".csv")
df_enriched["enrichment"] = enriched_phrase
df_depleted = pd.read_csv("final_results/final_table_all_tests_depleted_results_pfams_"+LABEL+"_collapse_"+COLLAPSE+".csv")
df_depleted["enrichment"] = depleted_phrase

df = pd.concat([df_enriched,df_depleted])

#remove unsure lmm results (lmm could be wrong and find the opposite enrichment)
df = df[~((df["lmm"] == 1) & (df["tests_passed"] == 1))]

# remove . at the end of the pfam_id
df["pfam_id"] = df["pfam_id"].str.split(".").str[0]
print (df)
df.to_csv("tables_to_export/combined_results_"+LABEL+"_collapse_"+COLLAPSE+".csv", index=False)

#get pivot table ready to export
def remove_sums(enrichment):
    df_pivot = pd.read_csv("final_results/pivot_table_"+enrichment+"_results_"+COLLAPSE+"_"+LABEL+"_at_least_0_clades.csv",index_col=0)
    output = df_pivot.tail(-1)
    output = output.drop("sum",axis=1)
    #remove . in pfam
    output.index = output.index.map(lambda x: ("('"+x.split("'")[1].split(".")[0]+"'"+x.split("'")[2]+"'"+x.split("'")[3]+"')"))

    if (enrichment == "enriched"):
        output.to_csv("tables_to_export/pivot_table_"+enriched_phrase+"_results_"+COLLAPSE+"_"+LABEL+".csv")
    else:
        output.to_csv("tables_to_export/pivot_table_"+depleted_phrase+"_results_"+COLLAPSE+"_"+LABEL+".csv")

remove_sums("enriched")
remove_sums("depleted")

# make 5 more tables (one for each test) (combined enriched and depleted)
#Fisher
