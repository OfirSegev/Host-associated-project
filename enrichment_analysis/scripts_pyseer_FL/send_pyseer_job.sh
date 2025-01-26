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

#send pyseer job with distances

pyseer \
--phenotypes $1 \
--pres $2 \
--distances $3 \
--cpu 15 \
--lineage \
--lineage-file $4_lineage_effects.out \
--max-af $5 \
--filter-pvalue 1 > outputs_and_analysis_FL/pyseer_$4.assoc

#send pyseer job without distances

pyseer \
--phenotypes $1 \
--pres $2 \
--no-distances \
--cpu 15 \
--max-af $5 \
--filter-pvalue 1 > outputs_and_analysis_FL/no_distances_pyseer_$4.assoc

