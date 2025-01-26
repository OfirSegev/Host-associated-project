import pandas as pd
import sys

#inputs: 1 = clade_number, 2 = biological sample number ,  3 = number_of_bootstrap_sample , 4 = collapse_num , 5 = sample_percentage

clade=sys.argv[1]
bio_sample_num=sys.argv[2]
bootstrap_sample_num=sys.argv[3]
COLLAPSE=sys.argv[4]
sample_percentage=float(sys.argv[5])

df=pd.read_csv("samples_combined_per_clade_"+COLLAPSE+".csv", dtype="str")

clade_df=df[df["clade_id"]==clade]
sample_size=int(len(clade_df)*sample_percentage)
clade_sample = clade_df.sample(n=sample_size)

sampled_genome_list=clade_sample["genome_id_"+bio_sample_num].tolist()

count=0
with open('bootstrap_'+bootstrap_sample_num+'_collapse_'+COLLAPSE+'__clade_'+clade+'_sample_'+bio_sample_num+'.txt', 'w') as f:
	for item in sampled_genome_list:
		if count==0:
			f.write("%s" % item)
		else:
			f.write(",%s" % item)
		count=count+1
