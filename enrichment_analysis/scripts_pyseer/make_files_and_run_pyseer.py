#!/sci/labs/asafle/alexlevylab/icore-data/miniconda3/envs/pyseer_env/bin/python3

import pandas as pd
import subprocess
import sys

python_path="/sci/labs/asafle/alexlevylab/icore-data/miniconda3/envs/pyseer_env/bin/python"

#inputs 1- clade number , 2- sample number, 3- collapse , 4 - pfam abundance threshold

#needs to make a new tree for combination clade & sample (get_tree function)

def make_phylogeny(CLADE,SAMPLE,NAME,COLLAPSE):
	tree= "n_newicktree_tree_"+str(COLLAPSE)+"__clade_" + str(CLADE) + "_sample_" + str(SAMPLE) + ".newick"
	p = subprocess.Popen([python_path, prefix+'phylogeny_distance.py', "--lmm", tree], stdout=subprocess.PIPE)
	std, err = p.communicate()

	with open(NAME + "_phylogeny_similarity_lmm.tsv", "a+") as f:
		f.write(std.decode("utf-8"))
		f.close()

	p = subprocess.Popen([python_path, prefix+'phylogeny_distance.py', tree], stdout=subprocess.PIPE)
	std, err = p.communicate()

	with open(NAME + "_phylogeny_distances.tsv", "a+") as f:
		f.write(std.decode("utf-8"))
		f.close()

def run_pyseer(NAME):
	# commands $1 = traits , $2 = g (pfam), $3 = distance, $4 = name_prefix , $5 = max-af (pfam abundance threshold)
	#run pyseer with distances + without distances
	commands = ["t_traits__" + NAME + ".tsv",
			"g_gene_presence_absence__" + NAME + ".tsv",
			NAME + "_phylogeny_distances.tsv",
			NAME,
			THRESHOLD]

	#p = subprocess.Popen(['bash', prefix+'send_pyseer_job.sh'] + commands, stdout=subprocess.PIPE)

	#std, err = p.communicate()

	#print(std.decode("utf-8"))

	#run lmm with distances
	commands2 = ["t_traits__" + NAME + ".tsv",
			"g_gene_presence_absence__" + NAME + ".tsv",
			NAME + "_phylogeny_similarity_lmm.tsv",
			NAME,
			THRESHOLD]

	q = subprocess.Popen(['bash', prefix+'send_pyseer_job_lmm_version.sh'] + commands2, stdout=subprocess.PIPE)

	std, err = q.communicate()

	print(std.decode("utf-8"))


if __name__ == "__main__":

	COLLAPSE= sys.argv[3]
	SAMPLE = sys.argv[2]
	THRESHOLD = sys.argv[4]
	CLADE = sys.argv[1]
	prefix="/sci/labs/asafle/alexlevylab/icore-data/alex/workspace/HA_vs_NHA/data_sources/final_genomes/predictions/scripts_pyseer/"

	NAME = "collapse_" + str(COLLAPSE) + "__clade_" + str(CLADE) + "_sample_" + str(SAMPLE)

	# make matrix from tree
	make_phylogeny(CLADE,SAMPLE,NAME,COLLAPSE)

	run_pyseer(NAME)
