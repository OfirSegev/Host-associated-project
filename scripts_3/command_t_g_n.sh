#!/bin/bash

#REMEMBER:
#conda activate scoary_env
# $1=trait file, $2=gene absence/presence file,$3=newick tree, $4=clade_sample, $5=Collapse, $6=label

scoary \
-t $1 \
-g $2 \
-n $3 \
-s 2 \
-p 1.0 \
-e 10 \
--threads 50 \
-o outputs_and_analysis/"$6"_collapse_"$5"__"$4"_all_genomes #creates a directory instead of file

#-r
#-o outputs_and_analysis/sample1/host_associated_collapse_0.4_sample1__clade_$4.results.csv #creates a directory instead of file
