library(tidyverse)
library(ggrepel)

df_HvN_hosts_with_clades <- read_tsv("HvN_Host-associated_larger_than_2_combined_best_candidates_to_plot_per_clade_collapse_G_with_odds_ratio.tsv")
df_HvN_hosts_with_clades <- mutate(df_HvN_hosts_with_clades, dataframe_type = "Host-associated")

df_HvN_environment_with_clades <- read_tsv("HvN_Environmental_larger_than_2_combined_best_candidates_to_plot_per_clade_collapse_G_with_odds_ratio.tsv")
df_HvN_environment_with_clades <- mutate(df_HvN_environment_with_clades, dataframe_type = "Environmental")

df <- rbind(df_HvN_hosts_with_clades,df_HvN_environment_with_clades)

#remove cllular_componenets
df <- df[df$GO_category != "cellular_component", ]

df$diffexpressed <- NA
#color trait hosts
df$diffexpressed[df$GO_term == "nitrogen fixation"] <- "nitrogen fixation"
df$diffexpressed[df$GO_term == "trehalose biosynthetic process"] <- "trehalose biosynthetic"
df$diffexpressed[df$GO_term == "antibiotic catabolic process"] <- "antibiotic catabolism"
df$diffexpressed[df$GO_term == "detoxification of mercury ion"] <- "detoxification of mercury"
df$diffexpressed[df$GO_term == "toxin transport"] <- "toxin transport"
df$diffexpressed[df$GO_term == "metalloendopeptidase activity"] <- "metalloendopeptidase activity"
df$diffexpressed[df$GO_term == "cytolysis in another organism"] <- "cytolysis in another organism"
#color traits en
df$diffexpressed[df$GO_term == "spore germination"] <- "spore germination"
df$diffexpressed[df$GO_term == "biotin transmembrane transporter activity"] <- "biotin transport"
df$diffexpressed[df$GO_term == "alpha-mannosidase activity"] <- "alpha-mannosidase activity"
df$diffexpressed[df$GO_term == "mannose metabolic process"] <- "mannose metabolism"

'
host_traits_to_color = c("nitrogen fixation","trehalose biosynthetic process","antibiotic catabolic process",
                         "detoxification of mercury ion","toxin transport","metalloendopeptidase activity",
                         "cytolysis in another organism")

env_traits_to_color = c("spore germination","biotin transmembrane transporter activity",
                        "alpha-mannosidase activity","phosphoenolpyruvate-dependent sugar phosphotransferase system")
'

#add labels to apear on plot
df$delabel <- NA
#df$delabel[df$diffexpressed != ""] <- paste0(df$Clade_id[df$diffexpressed != ""])
df$delabel[!is.na(df$diffexpressed)] <- paste0(df$Clade_id[!is.na(df$diffexpressed)], ":\n", df$diffexpressed[!is.na(df$diffexpressed)])
#df$delabel[!is.na(df$diffexpressed)] <- paste0(df$Clade_id[!is.na(df$diffexpressed)], ": ", df$diffexpressed[!is.na(df$diffexpressed)])

# Create a new column for alpha values
df$alpha_value <- ifelse(is.na(df$diffexpressed), 0.4, 1)


jpeg("volcano_plot_colored_hits_HvN_G_pfams_new_dimentions.jpeg",width=18,height=10.5,units="cm",res=3000)
# plot adding up all layers we have seen so far
ggplot(data=df, aes(x=log10(odds_ratio), y=-log10(fdr_corrected), shape=GO_category, col=diffexpressed, label=delabel,alpha=alpha_value)) +
  geom_point(size=1) + 
  theme_minimal() +
  theme(legend.title = element_text(size = 6.5), legend.text = element_text(size = 6), axis.text.x = element_text(size = 6), 
        axis.text.y = element_text(size = 6), legend.spacing = unit(0.1, 'cm'), legend.margin = margin(t = 0, b = 0))+
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

FL_df_HvN_hosts_with_clades <- read_tsv("FL_HvN_Host-associated_larger_than_2_combined_best_candidates_to_plot_per_clade_collapse_G_with_odds_ratio.tsv")
FL_df_HvN_hosts_with_clades <- mutate(FL_df_HvN_hosts_with_clades, dataframe_type = "Host-associated")

FL_df_HvN_environment_with_clades <- read_tsv("FL_HvN_Environmental_larger_than_2_combined_best_candidates_to_plot_per_clade_collapse_G_with_odds_ratio.tsv")
FL_df_HvN_environment_with_clades <- mutate(FL_df_HvN_environment_with_clades, dataframe_type = "Environmental")

FL_df <- rbind(FL_df_HvN_hosts_with_clades,FL_df_HvN_environment_with_clades)

#remove cllular_componenets
FL_df <- FL_df[FL_df$GO_category != "cellular_component", ]

FL_df$diffexpressed <- NA
#color trait hosts
FL_df$diffexpressed[FL_df$GO_term == "nitrogen fixation"] <- "nitrogen fixation"
#FL_df$diffexpressed[FL_df$GO_term == "nitrogenase activity"] <- "nitrogenase activity"
FL_df$diffexpressed[FL_df$GO_term == "folic acid-containing compound biosynthetic process"] <- "folic acid biosynthesis"
FL_df$diffexpressed[FL_df$GO_term == "protein secretion by the type III secretion system"] <- "type III secretion system"
FL_df$diffexpressed[FL_df$GO_term == "cell adhesion involved in biofilm formation"] <- "cell adhesion"
FL_df$diffexpressed[FL_df$GO_term == "cell adhesion"] <- "cell adhesion"
#FL_df$diffexpressed[FL_df$GO_term == "metalloexopeptidase activity"] <- "metalloexopeptidase activity"
FL_df$diffexpressed[FL_df$GO_term == "response to antibiotic"] <- "response to antibiotic"
#color trait env
#FL_df$diffexpressed[FL_df$GO_term == "anion binding"] <- "anion binding"
FL_df$diffexpressed[FL_df$GO_term == "cytoskeletal motor activity"] <- "cytoskeletal motor activity"
FL_df$diffexpressed[FL_df$GO_term == "bacterial-type flagellum-dependent cell motility"] <- "flagellum-dependent motility"
FL_df$diffexpressed[FL_df$GO_term == "protein repair"] <- "protein repair"
#FL_df$diffexpressed[FL_df$GO_term == "nitrogen utilization"] <- "nitrogen utilization"
FL_df$diffexpressed[FL_df$GO_term == "sporulation"] <- "sporulation"

'
host_FL_traits_to_color = c("nitrogen fixation","nitrogenase activity","folic acid-containing compound biosynthetic process",
                         "iron-sulfur cluster binding","protein secretion by the type III secretion system",
                         "cell adhesion involved in biofilm formation", "cell adhesion", "metalloexopeptidase activity",
                         "response to antibiotic")

env_FL_traits_to_color = c("anion binding",
                           "nitronate monooxygenase activity","cytoskeletal motor activity",
                           "bacterial-type flagellum-dependent cell motility","protein repair",
                           "nitrogen utilization","sporulation")
'

#add labels to apear on plot
FL_df$delabel <- NA
FL_df$delabel[!is.na(FL_df$diffexpressed)] <- paste0(FL_df$Clade_id[!is.na(FL_df$diffexpressed)], ":\n", FL_df$diffexpressed[!is.na(FL_df$diffexpressed)])
#df$delabel[df$diffexpressed != "NO"] <- df$GO_term[df$diffexpressed != "NO"]

# Create a new column for alpha values
FL_df$alpha_value <- ifelse(is.na(FL_df$diffexpressed), 0.4, 1)

jpeg("volcano_plot_colored_hits_HvN_G_AFC_new_dimentions3.jpeg",width=18,height=10.5,units="cm",res=3000)
#jpeg("short_volcano_plot_colored_hits_HvN_G_AFC.png",width=7.8,height=5.3,units="in",res=2000)
# plot adding up all layers we have seen so far
ggplot(data=FL_df, aes(x=log10(odds_ratio), y=-log10(fdr_corrected), shape=GO_category, col=diffexpressed, label=delabel, alpha=alpha_value)) +
  geom_point(size=1) + 
  theme_minimal() +
  theme(legend.title = element_text(size = 6.5), legend.text = element_text(size = 6), axis.text.x = element_text(size = 6), 
        axis.text.y = element_text(size = 6), legend.spacing = unit(0.1, 'cm'), legend.margin = margin(t = 0, b = 0))+
  labs(y= "-log10(q-value)", x = "log10(fold change)", shape = "GO category", color="GO term")+
  geom_text_repel(size = 2.3, max.overlaps = Inf, lineheight = 0.9, fontface = "bold", 
                  box.padding = 0.3, point.padding = 0.1, force = 3) +
  guides(shape = guide_legend(order = 1), color = guide_legend(order = 2))+
  #scale_color_manual(values=c("#696969","red","blue"))+
  geom_vline(xintercept=0, col="black", linetype="dotted")+
  scale_alpha_identity()
#geom_vline(xintercept=c(-1, 1), col="red") +
#geom_hline(yintercept=-log10(0.03), col="red")
dev.off()
###############################################################

df_AvP_animals_with_clades <- read_tsv("AvP_Animal-associated_larger_than_2_combined_best_candidates_to_plot_per_clade_collapse_G_with_odds_ratio.tsv")
df_AvP_animals_with_clades <- mutate(df_AvP_animals_with_clades, dataframe_type = "Animal-associated")

df_AvP_plants_with_clades <- read_tsv("AvP_Plant-associated_larger_than_2_combined_best_candidates_to_plot_per_clade_collapse_G_with_odds_ratio.tsv")
df_AvP_plants_with_clades <- mutate(df_AvP_plants_with_clades, dataframe_type = "Plant-associated")

df <- rbind(df_AvP_animals_with_clades,df_AvP_plants_with_clades)
#df <- df[order(df$fdr_corrected),]

#remove cllular_componenets
df <- df[df$GO_category != "cellular_component", ]

# add a column of NAs
df$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
df$diffexpressed[df$odds_ratio > 5 & df$fdr_corrected < 0.04] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
df$diffexpressed[df$odds_ratio < 0.2 & df$fdr_corrected < 0.04] <- "DOWN"


#add labels to apear on plot
df$delabel <- NA
df$delabel[df$diffexpressed != "NO"] <- paste0(df$Clade_id[df$diffexpressed != "NO"],": ",df$GO_term[df$diffexpressed != "NO"])
#df$delabel[df$diffexpressed != "NO"] <- df$GO_term[df$diffexpressed != "NO"]


# plot adding up all layers we have seen so far
ggplot(data=df, aes(x=log10(odds_ratio), y=-log10(fdr_corrected), shape=GO_category, col=diffexpressed, label=delabel, )) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel(size = 2.3, max.overlaps = 100) +
  scale_color_manual(values=c("blue", "black", "red"))+
  geom_vline(xintercept=0, col="black", linetype="dotted")

###############################################################
FL_df_AvP_animals_with_clades <- read_tsv("FL_AvP_Animal-associated_larger_than_2_combined_best_candidates_to_plot_per_clade_collapse_G_with_odds_ratio.tsv")
FL_df_AvP_animals_with_clades <- mutate(FL_df_AvP_animals_with_clades, dataframe_type = "Animal-associated")

FL_df_AvP_plants_with_clades <- read_tsv("FL_AvP_Plant-associated_larger_than_2_combined_best_candidates_to_plot_per_clade_collapse_G_with_odds_ratio.tsv")
FL_df_AvP_plants_with_clades <- mutate(FL_df_AvP_plants_with_clades, dataframe_type = "Plant-associated")

FL_df <- rbind(FL_df_AvP_animals_with_clades,FL_df_AvP_plants_with_clades)
#df <- df[order(df$fdr_corrected),]

#remove cllular_componenets
FL_df <- FL_df[FL_df$GO_category != "cellular_component", ]

# add a column of NAs
FL_df$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
FL_df$diffexpressed[FL_df$odds_ratio > 5 & FL_df$fdr_corrected < 0.04] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
FL_df$diffexpressed[FL_df$odds_ratio < 0.2 & FL_df$fdr_corrected < 0.04] <- "DOWN"


#add labels to apear on plot
FL_df$delabel <- NA
FL_df$delabel[FL_df$diffexpressed != "NO"] <- paste0(FL_df$Clade_id[FL_df$diffexpressed != "NO"],": ",FL_df$GO_term[FL_df$diffexpressed != "NO"])
#df$delabel[df$diffexpressed != "NO"] <- df$GO_term[df$diffexpressed != "NO"]


# plot adding up all layers we have seen so far
ggplot(data=FL_df, aes(x=log10(odds_ratio), y=-log10(fdr_corrected), shape=GO_category, col=diffexpressed, label=delabel, )) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel(size = 2.3, max.overlaps = 100) +
  scale_color_manual(values=c("blue", "black", "red"))+
  geom_vline(xintercept=0, col="black", linetype="dotted")

