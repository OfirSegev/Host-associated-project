import pandas as pd
import sys

operon_id = sys.argv[1]
COLLAPSE="G"
LABEL="HvN"
ENRICHMENT="Host-associated"

df = pd.read_csv(ENRICHMENT+"_"+LABEL+"_collapse_"+COLLAPSE+"_operon_parameters_with_genome_data.csv")

#filter for operon_id and keep only img operons
df_img = df[(df["operon_id"] == operon_id) & (df['gene_ID'].str.startswith(('2', '6')))]

#sort by completeness
df_img = df_img.sort_values(by = ["Start"], ascending=True)
df_img = df_img.sort_values(by = ["completeness","genome_ID","scaffold_ID"], ascending=False)


#to print
best_genome_id = df_img["genome_ID"].iloc[0]
best_scaffold_id = df_img["scaffold_ID"].iloc[0]
operon_start = str(df_img[(df_img["genome_ID"] == best_genome_id) & (df_img["scaffold_ID"] == best_scaffold_id)]["Start"].iloc[0])
operon_end = str(df_img[(df_img["genome_ID"] == best_genome_id) & (df_img["scaffold_ID"] == best_scaffold_id)]["End"].iloc[-1])
first_gene = df_img["gene_ID"].iloc[0]

print (df_img)
print ("BEST HIT")
print ("genome_id: "+best_genome_id)
print ("scaffold: "+best_scaffold_id)
print ("operon start and end: "+operon_start+"  "+operon_end)
print ("first_gene: "+first_gene)

#to csv
df_img.to_csv("operons_info/"+operon_id+"_details.csv", index=False)