library(tidyverse)
library(ggplot2)
library(pheatmap)
library(dendextend)
library(readxl)
library(grid)

# Function to remove columns with all zeros
remove_columns_with_all_zeros <- function(data_frame) {
  # Identify columns where all values are 0
  zero_columns <- sapply(data_frame, function(col) all(col == 0))
  
  # Remove identified columns
  data_frame <- data_frame[, !zero_columns]
  
  return(data_frame)
}

#df = read.csv("pvalue_0.05_evalue_e-5_coverage_0.7/pivot_table_enriched_results_G_HvN_at_least_0_clades_colored.csv")
#name = "best_results_pfams_HvN_G_Host-associated_14_30.csv"
#name = "best_results_FL_HvN_G_Host-associated_10_30.csv"
name = "best_results_pfams_AvP_F_both_13_30.csv"
folder = "AvP_F"
df = read.csv(paste0(folder,"/",name))

#only for pfam_HvN_G_host
df$X[df$X == "('PF04809', 'HupH hydrogenase expression protein, C-terminal conserved region')"] <- "('PF04809', 'HupH hydrogenase expression protein)"
#only for FL_Environment_HvN_G
df$X[df$X == "('K6ULH7', 'Dihydrolipoyllysine-residue succinyltransferase component of 2-oxoglutarate dehydrogenase complex')"] <- "('K6ULH7', 'Component of 2-oxoglutarate dehydrogenase complex')"
#only for FL_Environment_HvN_F
df$X[df$X == "('A0A2G9U759', 'Glyceraldehyde 3-phosphate dehydrogenase domain protein')"] <- "('A0A2G9U759', 'Glyceraldehyde 3-phosphate dehydrogenase protein')"
#only for pfam_AvP_F
df$X[df$X == "('PF02275', 'Linear amide C-N hydrolases, choloylglycine hydrolase family')"] <- "('PF02275', 'Linear amide C-N hydrolases')"



#order of genera
#genera_phylum_info = read.csv("genera_phylum_info.csv")

#pfams to apear in heatmap
# 70ad47=green, 5b9bd5=blue, ffc000=yellow, 00ffff=cyan, a5a5a5=grey, ff00ff=purple, ed7d31=orange
#colors_to_plot=c("#70ad47","#5b9bd5","#ffc000","#00ffff")

#filtered_df <- df[df$Hex_Colors %in% colors_to_plot, ]
filtered_df <- df

# Set row names to the first column (Pfam_id)
row.names(filtered_df) <- filtered_df[, 1]
filtered_df <- filtered_df[, -1]

#maybe keep the color information
#filtered_df <- subset(filtered_df, select = -Hex_Colors)

#remove columns with sum less than x
#column_sums <- colSums(filtered_df)
#selected_columns <- column_sums > 2
#filtered_df <- filtered_df[, selected_columns]

#remove sum column
filtered_df <- subset(filtered_df, select = -c(sum))

#sort
#columns_order <- sort(colnames(filtered_df))
#columns_order <- c("Paraburkholderia","Bradyrhizobium","Sinorhizobium","Mesorhizobium",
#                   "Burkholderia","Vibrio","Pseudomonas.E","Acinetobacter","Escherichia","Rhodanobacter")

#filtered_df <- filtered_df[, columns_order]

#create data for heatmap
data <- as.matrix(filtered_df)

# Customize the color palette (choose your desired colors)
#for 5 tests with white
my_color_palette <- c("#FFFFFF", "#FFD3D3", "#FFAFAF", "#FF8C8C", "#FF6969", "#FF4646")
#for 4 tests with white
my_color_palette <- c("#FFFFFF", "#FFD3D3", "#FF8C8C", "#FF6969", "#FF4646")
#for 5 tests with blue
#my_color_palette <- c("#86A9CC", "#FFD3D3", "#FFAFAF", "#FF8C8C", "#FF6969", "#FF4646")
#for 4 tests with blue
# my_color_palette <- c("#86A9CC", "#FFD3D3", "#FF8C8C", "#FF6969", "#FF4646")
#for both enriched and depleted 11 pallete
my_color_palette <- c("#468AFF","#69A1FF","#8CB8FF","#AFCFFF","#D3E7FF","#FFFFFF", "#FFD3D3", "#FFAFAF", "#FF8C8C", "#FF6969", "#FF4646")
#for 5 enriched and 4 depleted AvP
my_color_palette <- c("#468AFF","#69A1FF","#AFCFFF","#D3E7FF","#FFFFFF", "#FFD3D3", "#FFAFAF", "#FF8C8C", "#FF6969", "#FF4646")
#for 5 enriched and 3 depleted AvP
#my_color_palette <- c("#468AFF","#69A1FF","#AFCFFF","#FFFFFF", "#FFD3D3", "#FFAFAF", "#FF8C8C", "#FF6969", "#FF4646")




#make col names italics
genus_names_italics <- lapply(colnames(data), function(x) bquote(italic(.(x))))

# Create the heatmap using pheatmap with your custom color palette
pheatmap(data, color = my_color_palette, cluster_rows = T,cluster_cols = T,
         labels_col = as.expression(genus_names_italics))

FL_pfams=strsplit(name, split = "_")[[1]][3]
label=strsplit(name, split = "_")[[1]][4]
collapse=strsplit(name, split = "_")[[1]][5]
enrichment=strsplit(name, split = "_")[[1]][6]


jpeg(paste0(folder,"/top_hits_heatmap_13_",FL_pfams,"_",label,"_",collapse,"_",enrichment,".jpeg"),width=8.5,height=6.2,units="in",res=2000)
pheatmap(data, color = my_color_palette, cluster_rows = T,cluster_cols = T,
         labels_col = as.expression(genus_names_italics))
dev.off()






