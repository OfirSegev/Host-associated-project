library(tidyverse)
library(ggrepel)

df_HvN_hosts_with_clades <- read_tsv("HvN_Host-associated_larger_than_2_combined_best_candidates_to_plot_per_clade_collapse_F_with_odds_ratio.tsv")
df_HvN_hosts_with_clades <- mutate(df_HvN_hosts_with_clades, dataframe_type = "Host-associated")

df_HvN_environment_with_clades <- read_tsv("HvN_Environmental_larger_than_2_combined_best_candidates_to_plot_per_clade_collapse_F_with_odds_ratio.tsv")
df_HvN_environment_with_clades <- mutate(df_HvN_environment_with_clades, dataframe_type = "Environmental")

df <- rbind(df_HvN_hosts_with_clades,df_HvN_environment_with_clades)

#remove cllular_componenets
df <- df[df$GO_category != "cellular_component", ]

df$diffexpressed <- NA
#color trait hosts
df$diffexpressed[df$GO_term == "nitrogen fixation"] <- "nitrogen fixation"
df$diffexpressed[df$GO_term == "carbohydrate metabolic process"] <- "carbohydrate metabolism"
df$diffexpressed[df$GO_term == "heme oxidation"] <- "heme oxidation"
df$diffexpressed[df$GO_term == "xylan catabolic process"] <- "xylan catabolism"
df$diffexpressed[df$GO_term == "hydrolase activity, hydrolyzing O-glycosyl compounds"] <- "hydrolyzing O-glycosyl"
#color traits en
df$diffexpressed[df$GO_term == "NADH dehydrogenase (ubiquinone) activity"] <- "NADH dehydrogenase"
df$diffexpressed[df$GO_term == "bacterial-type flagellum organization"] <- "flagellum organization"
df$diffexpressed[df$GO_term == "bacterial-type flagellum-dependent cell motility"] <- "flagellum-dependent motility"
df$diffexpressed[df$GO_term == "chemotaxis"] <- "chemotaxis"
df$diffexpressed[df$GO_term == "cytoskeletal motor activity"] <- "cytoskeletal motor activity"


#add labels to apear on plot
df$delabel <- NA
#df$delabel[df$diffexpressed != ""] <- paste0(df$Clade_id[df$diffexpressed != ""])
df$delabel[!is.na(df$diffexpressed)] <- paste0(df$Clade_id[!is.na(df$diffexpressed)], ":\n", df$diffexpressed[!is.na(df$diffexpressed)])
#df$delabel[!is.na(df$diffexpressed)] <- paste0(df$Clade_id[!is.na(df$diffexpressed)], ": ", df$diffexpressed[!is.na(df$diffexpressed)])

# Create a new column for alpha values
df$alpha_value <- ifelse(is.na(df$diffexpressed), 0.4, 1)


jpeg("volcano_plot_colored_hits_HvN_F_pfams.png",width=8,height=5.25,units="in",res=2000)
# plot adding up all layers we have seen so far
ggplot(data=df, aes(x=log10(odds_ratio), y=-log10(fdr_corrected), shape=GO_category, col=diffexpressed, label=delabel,alpha=alpha_value)) +
  geom_point(size=1) + 
  theme_minimal() +
  theme(legend.text = element_text(size = 7), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6))+
  labs(y= "-log10(q-value)", x = "log10(fold change)", shape = "GO category", color="GO term")+
  geom_text_repel(size = 2.3,max.overlaps = 100,lineheight =0.9, fontface="bold") +
  guides(shape = guide_legend(order = 1), color = guide_legend(order = 2))+
  #scale_color_manual(values=c("#696969","red","blue"))+
  geom_vline(xintercept=0, col="black", linetype="dotted")+
  scale_alpha_identity()
  #geom_vline(xintercept=c(-1, 1), col="red") +
  #geom_hline(yintercept=-log10(0.03), col="red")
dev.off()

#########################################################################

FL_df_HvN_hosts_with_clades <- read_tsv("FL_HvN_Host-associated_larger_than_2_combined_best_candidates_to_plot_per_clade_collapse_F_with_odds_ratio.tsv")
FL_df_HvN_hosts_with_clades <- mutate(FL_df_HvN_hosts_with_clades, dataframe_type = "Host-associated")

FL_df_HvN_environment_with_clades <- read_tsv("FL_HvN_Environmental_larger_than_2_combined_best_candidates_to_plot_per_clade_collapse_F_with_odds_ratio.tsv")
FL_df_HvN_environment_with_clades <- mutate(FL_df_HvN_environment_with_clades, dataframe_type = "Environmental")

FL_df <- rbind(FL_df_HvN_hosts_with_clades,FL_df_HvN_environment_with_clades)

#remove cllular_componenets
FL_df <- FL_df[FL_df$GO_category != "cellular_component", ]

FL_df$diffexpressed <- NA
#color trait hosts
FL_df$diffexpressed[FL_df$GO_term == "nitrogen fixation"] <- "nitrogen fixation"
FL_df$diffexpressed[FL_df$GO_term == "nitrogenase activity"] <- "nitrogenase activity"
#FL_df$diffexpressed[FL_df$GO_term == "DNA binding"] <- "DNA binding"
FL_df$diffexpressed[(FL_df$GO_term == "protein secretion by the type III secretion system") & (FL_df$dataframe_type == "Host-associated")] <- "type III secretion system"
FL_df$diffexpressed[FL_df$GO_term == "carbohydrate metabolic process"] <- "carbohydrate metabolism"
FL_df$diffexpressed[FL_df$GO_term == "hydrolase activity, hydrolyzing O-glycosyl compounds"] <- "hydrolyzing O-glycosyl"
FL_df$diffexpressed[FL_df$GO_term == "animal organ development"] <- "animal organ development"
FL_df$diffexpressed[FL_df$GO_term == "xylan catabolic process"] <- "xylan catabolism"
#FL_df$diffexpressed[FL_df$GO_term == "toxin activity"] <- "toxin activity"
FL_df$diffexpressed[FL_df$GO_term == "biological process involved in interaction with host"] <- "interaction with host"
FL_df$diffexpressed[FL_df$GO_term == "symbiont entry into host cell"] <- "entry into host cell"
FL_df$diffexpressed[FL_df$GO_term == "muscle organ development"] <- "muscle organ development"
FL_df$diffexpressed[FL_df$GO_term == "visual system development"] <- "visual system development"
#color trait env
FL_df$diffexpressed[FL_df$GO_term == "response to oxidative stress"] <- "response to oxidative stress"
FL_df$diffexpressed[FL_df$GO_term == "cytoskeletal motor activity"] <- "cytoskeletal motor activity"
FL_df$diffexpressed[FL_df$GO_term == "bacterial-type flagellum-dependent cell motility"] <- "flagellum-dependent motility"
FL_df$diffexpressed[FL_df$GO_term == "chemotaxis"] <- "chemotaxis"
FL_df$diffexpressed[FL_df$GO_term == "NADH dehydrogenase (ubiquinone) activity"] <- "NADH dehydrogenase"
FL_df$diffexpressed[FL_df$GO_term == "bacterial-type flagellum assembly"] <- "flagellum assembly"


#add labels to apear on plot
FL_df$delabel <- NA
FL_df$delabel[!is.na(FL_df$diffexpressed)] <- paste0(FL_df$Clade_id[!is.na(FL_df$diffexpressed)], ":\n", FL_df$diffexpressed[!is.na(FL_df$diffexpressed)])
#df$delabel[df$diffexpressed != "NO"] <- df$GO_term[df$diffexpressed != "NO"]

# Create a new column for alpha values
FL_df$alpha_value <- ifelse(is.na(FL_df$diffexpressed), 0.4, 1)

jpeg("volcano_plot_colored_hits_HvN_F_AFC.png",width=8,height=5.25,units="in",res=2000)
#jpeg("short_volcano_plot_colored_hits_HvN_G_AFC.png",width=7.8,height=5.3,units="in",res=2000)
# plot adding up all layers we have seen so far
ggplot(data=FL_df, aes(x=log10(odds_ratio), y=-log10(fdr_corrected), shape=GO_category, col=diffexpressed, label=delabel, alpha=alpha_value)) +
  geom_point(size=1) + 
  theme_minimal() +
  theme(legend.text = element_text(size = 7), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6))+
  labs(y= "-log10(q-value)", x = "log10(fold change)", shape = "GO category", color="GO term")+
  geom_text_repel(size = 2.3,max.overlaps = 100,lineheight =0.9, fontface="bold") +
  guides(shape = guide_legend(order = 1), color = guide_legend(order = 2))+
  #scale_color_manual(values=c("#696969","red","blue"))+
  geom_vline(xintercept=0, col="black", linetype="dotted")+
  scale_alpha_identity()
#geom_vline(xintercept=c(-1, 1), col="red") +
#geom_hline(yintercept=-log10(0.03), col="red")
dev.off()
###############################################################

df_AvP_animals_with_clades <- read_tsv("AvP_Animal-associated_larger_than_2_combined_best_candidates_to_plot_per_clade_collapse_F_with_odds_ratio.tsv")
df_AvP_animals_with_clades <- mutate(df_AvP_animals_with_clades, dataframe_type = "Animal-associated")

df_AvP_plants_with_clades <- read_tsv("AvP_Plant-associated_larger_than_2_combined_best_candidates_to_plot_per_clade_collapse_F_with_odds_ratio.tsv")
df_AvP_plants_with_clades <- mutate(df_AvP_plants_with_clades, dataframe_type = "Plant-associated")

df <- rbind(df_AvP_animals_with_clades,df_AvP_plants_with_clades)
#df <- df[order(df$fdr_corrected),]

#remove cllular_componenets
df <- df[df$GO_category != "cellular_component", ]

df$diffexpressed <- NA
#color trait hosts
df$diffexpressed[df$GO_term == "detoxification of mercury ion"] <- "detoxification of mercury"
df$diffexpressed[df$GO_term == "biotin transmembrane transporter activity"] <- "biotin transmembrane transport"
df$diffexpressed[df$GO_term == "toxin sequestering activity"] <- "toxin sequestering activity"
#color traits en
df$diffexpressed[df$GO_term == "hydrolase activity, hydrolyzing O-glycosyl compounds"] <- "hydrolyzing O-glycosyl"
df$diffexpressed[df$GO_term == "carbohydrate metabolic process"] <- "carbohydrate metabolism"
df$diffexpressed[df$GO_term == "cobalamin biosynthetic process"] <- "cobalamin biosynthesis"
df$diffexpressed[df$GO_term == "xylan catabolic process"] <- "xylan catabolism"
df$diffexpressed[df$GO_term == "nitrogen fixation"] <- "nitrogen fixation"
df$diffexpressed[df$GO_term == "antioxidant activity"] <- "antioxidant activity"


#add labels to apear on plot
df$delabel <- NA
#df$delabel[df$diffexpressed != ""] <- paste0(df$Clade_id[df$diffexpressed != ""])
df$delabel[!is.na(df$diffexpressed)] <- paste0(df$Clade_id[!is.na(df$diffexpressed)], ":\n", df$diffexpressed[!is.na(df$diffexpressed)])
#df$delabel[!is.na(df$diffexpressed)] <- paste0(df$Clade_id[!is.na(df$diffexpressed)], ": ", df$diffexpressed[!is.na(df$diffexpressed)])

# Create a new column for alpha values
df$alpha_value <- ifelse(is.na(df$diffexpressed), 0.4, 1)


jpeg("volcano_plot_colored_hits_AvP_F_pfams_test.jpeg",width=18,height=10.5,units="cm",res=2000)
# plot adding up all layers we have seen so far
ggplot(data=df, aes(x=log10(odds_ratio), y=-log10(fdr_corrected), shape=GO_category, col=diffexpressed, label=delabel,alpha=alpha_value)) +
  geom_point(size=1) + 
  theme_minimal() +
  theme(legend.title = element_text(size = 6.5), legend.text = element_text(size = 6), axis.text.x = element_text(size = 6), 
        axis.text.y = element_text(size = 6), legend.spacing = unit(0.1, 'cm'), legend.margin = margin(t = 0, b = 0))+
  labs(y= "-log10(q-value)", x = "log10(fold change)", shape = "GO category", color="GO term")+
  geom_text_repel(size = 2.3, max.overlaps = Inf, lineheight = 0.9, fontface = "bold", 
                  box.padding = 0.1, point.padding = 0.1, force = 2) +
  guides(shape = guide_legend(order = 1), color = guide_legend(order = 2))+
  #scale_color_manual(values=c("#696969","red","blue"))+
  geom_vline(xintercept=0, col="black", linetype="dotted")+
  #scale_color_brewer(palette = "Set3") +  # Use a different Brewer palette
  scale_color_hue(direction = -1,h.start=45) +  # Close to default ggplot2 palette
  scale_alpha_identity()
#geom_vline(xintercept=c(-1, 1), col="red") +
#geom_hline(yintercept=-log10(0.03), col="red")
dev.off()

###############################################################
FL_df_AvP_animals_with_clades <- read_tsv("FL_AvP_Animal-associated_larger_than_2_combined_best_candidates_to_plot_per_clade_collapse_F_with_odds_ratio.tsv")
FL_df_AvP_animals_with_clades <- mutate(FL_df_AvP_animals_with_clades, dataframe_type = "Animal-associated")

FL_df_AvP_plants_with_clades <- read_tsv("FL_AvP_Plant-associated_larger_than_2_combined_best_candidates_to_plot_per_clade_collapse_F_with_odds_ratio.tsv")
FL_df_AvP_plants_with_clades <- mutate(FL_df_AvP_plants_with_clades, dataframe_type = "Plant-associated")

FL_df <- rbind(FL_df_AvP_animals_with_clades,FL_df_AvP_plants_with_clades)
#df <- df[order(df$fdr_corrected),]

#remove cllular_componenets
FL_df <- FL_df[FL_df$GO_category != "cellular_component", ]

FL_df$diffexpressed <- NA
#color trait hosts
FL_df$diffexpressed[FL_df$GO_term == "animal organ development"] <- "animal organ development"
FL_df$diffexpressed[FL_df$GO_term == "nitrogen utilization"] <- "nitrogen utilization"
FL_df$diffexpressed[FL_df$GO_term == "mercury ion transmembrane transporter activity"] <- "mercury ion transport"
#FL_df$diffexpressed[FL_df$GO_term == "aromatic amino acid family biosynthetic process"] <- "aromatic amino acid biosynthetsis"
FL_df$diffexpressed[FL_df$GO_term == "DNA transposition"] <- "DNA transposition"
#color traits en
FL_df$diffexpressed[FL_df$GO_term == "nitrate assimilation"] <- "nitrate assimilation"
FL_df$diffexpressed[FL_df$GO_term == "carbohydrate metabolic process"] <- "carbohydrate metabolism"
#FL_df$diffexpressed[FL_df$GO_term == "carbohydrate catabolic process"] <- "carbohydrate catabolism"
FL_df$diffexpressed[FL_df$GO_term == "chemotaxis"] <- "chemotaxis"
FL_df$diffexpressed[FL_df$GO_term == "hydrolase activity, hydrolyzing O-glycosyl compounds"] <- "hydrolyzing O-glycosyl"
FL_df$diffexpressed[FL_df$GO_term == "potassium ion transport"] <- "potassium ion transport"


#add labels to apear on plot
FL_df$delabel <- NA
FL_df$delabel[!is.na(FL_df$diffexpressed)] <- paste0(FL_df$Clade_id[!is.na(FL_df$diffexpressed)], ":\n", FL_df$diffexpressed[!is.na(FL_df$diffexpressed)])
#df$delabel[df$diffexpressed != "NO"] <- df$GO_term[df$diffexpressed != "NO"]

# Create a new column for alpha values
FL_df$alpha_value <- ifelse(is.na(FL_df$diffexpressed), 0.4, 1)

jpeg("volcano_plot_colored_hits_AvP_F_AFC.jpeg",width=8,height=5.25,units="in",res=3000)
#jpeg("short_volcano_plot_colored_hits_HvN_G_AFC.png",width=7.8,height=5.3,units="in",res=2000)
# plot adding up all layers we have seen so far
ggplot(data=FL_df, aes(x=log10(odds_ratio), y=-log10(fdr_corrected), shape=GO_category, col=diffexpressed, label=delabel, alpha=alpha_value)) +
  geom_point(size=1) + 
  theme_minimal() +
  theme(legend.text = element_text(size = 7), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6))+
  labs(y= "-log10(q-value)", x = "log10(fold change)", shape = "GO category", color="GO term")+
  geom_text_repel(size = 2.3,max.overlaps = 100,lineheight =0.9, fontface="bold") +
  guides(shape = guide_legend(order = 1), color = guide_legend(order = 2))+
  #scale_color_manual(values=c("#696969","red","blue"))+
  geom_vline(xintercept=0, col="black", linetype="dotted")+
  scale_alpha_identity()
#geom_vline(xintercept=c(-1, 1), col="red") +
#geom_hline(yintercept=-log10(0.03), col="red")
dev.off()

