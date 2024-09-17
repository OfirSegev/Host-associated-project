library(ape)
library(readr)

#inputs
args = commandArgs(trailingOnly=TRUE)
clade_id <- args[1]
sample <- args[2]
collapse <- args[3]
tree_path <- args[4]

# read fungi tree
tree <- read.tree(tree_path)

# read my df
df <- read_tsv(paste0("genome_lists_of_clades/list_of_clade_", clade_id, "_sample_", sample, ".tsv"), col_names = "genome_id")

# filter tree to keep only the relevant genome ids
filtered_tree <- drop.tip(tree, setdiff(tree$tip.label, df$genome_id))

# write out update represent_tree for iTol
write.tree(filtered_tree, file = paste0("tree_files_sample_", sample, "/subtrees_for_scoary_collapse_", collapse, "/tree_", collapse, "__cladeID_", clade_id, ".newick"))


