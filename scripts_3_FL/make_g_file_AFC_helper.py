import pandas as pd
import sys

#inputs
#1 - relevant gs clade_sample (gene clusters for each genomes)
#2 - AFC presence absence name
#3 - relvant genomelist name
#4 - threshold?
#5 - AFC full counts name

df = pd.read_csv(sys.argv[1], header = None, names = ['genome', 'repID'], dtype = 'str')
df['one'] = 1

f = pd.pivot_table(data = df, index = 'repID', columns = 'genome', values = 'one', aggfunc = 'sum')
f.index.name = None
f.columns.name = None
f.fillna(0, inplace = True)

g = f > 0
g = g.astype(int)
# g has all the data
print(g)
# filter out the very very common and the very very uncommon
df_merged_ones = g
df_merged_ones["score"] = df_merged_ones.apply(sum, axis=1) / len(df_merged_ones.columns)

THRESHOLD=float(sys.argv[4])

#removes pfams rows with pfams that apear in most genomes (threshold)
df_filtered_presence_absence = df_merged_ones[df_merged_ones["score"] < THRESHOLD]

#remove pfams rows with pfams that does not apear in most genomes (1-threshold)
min_threshold = 1-THRESHOLD
df_filtered_presence_absence = df_filtered_presence_absence[df_filtered_presence_absence["score"] > min_threshold]

#remove the added column
del df_filtered_presence_absence["score"]

df_filtered_presence_absence.to_csv(sys.argv[2])

r = pd.DataFrame(columns = g.columns)
r.to_csv(sys.argv[3], index = None)

############3

full_counts = f.loc[df_filtered_presence_absence.index]
print(full_counts)
full_counts.to_csv(sys.argv[5])

