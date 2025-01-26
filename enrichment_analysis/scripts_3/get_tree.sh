#!/bin/bash

CLADE=$1
SAMPLE=$2
COLLAPSE=$3

cp tree_files_sample_$SAMPLE/subtrees_for_scoary_collapse_$COLLAPSE/tree_"$COLLAPSE"__cladeID_$CLADE.newick n_newicktree_tree_"$COLLAPSE"__clade_"$CLADE"_sample_"$SAMPLE".newick

echo "tree moved"
