#!/bin/bash

#pyseer \
#--phenotypes  \ # HA or NHA (t file)
#--pres \ # presence absence of pfams (g file)
#--distances \ # seems to be for the fixed model
#--similarity  \ # they provide a custom python script to generate this kinship matrix. the custom script needs phylogeny. this seems to be for mixed model.
#--lmm \ #mixed model analysis...
#--lineage \ # lineage cluster also is a possibility, need to give poppunk clusters, though
#--lineage-file \
#--min-af 0.01 \
#--max-af 0.99 \
#--cpu 15 \
#--filter-pvalue 1E-8 > pyseer.assoc

#python scripts/phylogeny_distance.py --lmm core_genome.tree > phylogeny_similarity.tsv

#send lmm with distances

pyseer \
--lmm \
--phenotypes $1 \
--pres $2 \
--similarity $3 \
--cpu 15 \
--max-af $5 \
--filter-pvalue 1 > outputs_and_analysis/lmm_$4.assoc


#pyseer \
#--lmm \
#--phenotypes t_traits__collapse_0.18__clade_491_sample_1.tsv \
#--pres g_gene_presence_absence__collapse_0.18__clade_491_sample_1.tsv \
#--similarity collapse_0.18__clade_491_sample_1_phylogeny_similarity_lmm.tsv \
#--cpu 15 \
#--filter-pvalue 1E-4 > pyseer_lmm_test.assoc


