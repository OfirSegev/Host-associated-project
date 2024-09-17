import pandas as pd
import sys

#input collapse number

COLLAPSE=sys.argv[1]
MIN_CLADE_SIZE=float(sys.argv[2])
MIN_ENTROPY=float(sys.argv[3])

sizes=pd.read_csv("tree_files_sample_1/tree_collapse_"+COLLAPSE+"_clade_sizes.tsv", sep="\t")

entropy=pd.read_csv("tree_files_sample_1/sample1_"+COLLAPSE+"_clade_entropy.csv", names=["clade_id","entropy"])

sizes["clade_id"] = sizes["clade_id"].astype(str)

mer = sizes.merge(entropy, on = "clade_id", how="left")

#keep the best clades
mer = mer.loc[(mer["size"]>=MIN_CLADE_SIZE) & (mer["entropy"]>=MIN_ENTROPY)]
#remove very small uninfomative clades
mer = mer.loc[~((mer["size"]<=100) & (mer["entropy"]<=0.5))]
mer = mer.loc[~((mer["size"]<=50) & (mer["entropy"]<=0.8))]

#remove unclassified clade
mer = mer.loc[mer["clade_id"]!="unclassified"]

print (mer)

mer.to_csv("scored_clades_collapse_"+COLLAPSE+".csv", index=False)
