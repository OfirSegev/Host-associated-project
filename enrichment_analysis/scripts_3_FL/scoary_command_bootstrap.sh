#!/bin/bash

#REMEMBER:
#conda activate scoary_env
# $1=trait file, $2=gene absence/presence file,$3=newick tree, $4=clade_sample, $5=bootstrap, $6=genome names file, $7=Collapse, $8=label

scoary \
-t $1 \
-g $2 \
-n $3 \
-s 2 \
-p 1.0 \
-e 10 \
--threads 50 \
-o outputs_and_analysis_FL/"$8"_collapse_$7__$4_bootstrap_$5 \
-r $6

#-r
#-o outputs_and_analysis_FL/sample1/host_associated_collapse_0.4_sample1__clade_$4_bootstrap_$5.results.csv \ #creats a directory instead of file
