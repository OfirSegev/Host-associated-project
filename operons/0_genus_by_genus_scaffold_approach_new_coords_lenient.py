import pandas as pd
import concurrent.futures

#parameters
base = "/sci/labs/asafle/alexlevylab/icore-data/alex/workspace/HA_vs_NHA/data_sources/final_genomes/predictions/"
#run= "HvN_collapse_G_distance_0.8_pvalue_0.05_evalue_e-5_coverage_0.7_FL_run"
run= "HvN_collapse_G_distance_0.8"
TESTS_PASSED = 1
COLLAPSE="G"
LABEL="HvN"
ENRICHMENT="Host-associated"
window_size = 20000

#######
def fix_patric(gff):
    new_name = gff['scaffold_ID'] + "_" + gff['gene_ID'].split("_")[1]
    return new_name
########

#mapper = pd.read_csv("/sci/labs/asafle/alexlevylab/icore-data/alex/workspace/HA_vs_NHA/data_sources/final_genomes/genome_data/full_protein_scripts_and_data/mmseqs_to_AFC_master_sheet.csv", 
#usecols = ['lc_member', 'genome', 'repID'])
#print(mapper)

########

metadata = pd.read_csv(base + run + "/tree_files_sample_1/tree_sample1_collapse_"+COLLAPSE+"_clade_with_metadata.tsv", sep="\t", usecols = ['genome_id', 'clade_id','source', 'label_'+LABEL ])

#TODO for test
#metadata = metadata[metadata["clade_id"] == "Klebsiella"]

#metadata = pd.read_csv("/sci/labs/asafle/alexlevylab/icore-data/alex/workspace/HA_vs_NHA/data_sources/final_genomes/samples/HvN_samples/new_final_genome_0.8_new_sample_HvN_number_1.csv", usecols = ['genome_id', 'source', 'label_HvN', 'genus'])

#metadata = metadata.loc[~metadata['source'].isin(["poyet", "zou"])] #these sources have badly mismatched gff and faa files, no quick fix

#metadata['genus'] = metadata['genus'].str.replace("_", "-")

###########

#filter for PCOA improvment
'''
#filename = "/sci/labs/asafle/alexlevylab/icore-data/alex/workspace/HA_vs_NHA/data_sources/final_genomes/predictions/pcoa/table_outputs/HvN_collapse_G_distance_0.8_pvalue_0.05_evalue_e-5_coverage_0.7_FL_run__confusion_matrix_metrics_enriched_and_depleted_pcoa_liblin_C0.05_processed_deltaAIC_min2__deltaF1_gt0.tsv"
filename = "/sci/labs/asafle/alexlevylab/icore-data/alex/workspace/HA_vs_NHA/data_sources/final_genomes/predictions/pcoa/table_outputs/HvN_collapse_G_distance_0.8_pvalue_0.05_evalue_e-5_coverage_0.7_FL_run__confusion_matrix_metrics_pcoa_liblin_C0.05_distmet_canberra_processed_deltaAIC_min2__deltaF1_gt0.tsv"

pcoa_processed = pd.read_csv(filename, sep = '\t')
pcoa_processed_goods = pcoa_processed.loc[pcoa_processed['improvement'] == True]

metadata = metadata.loc[metadata['genus'].isin(pcoa_processed_goods['clade'].unique())]

#TEMP
#metadata = metadata.loc[metadata['genus'] == "Aliivibrio"]
#metadata = metadata.loc[metadata['genus'].isin(['Paenibacillus-A', 'Arthrobacter-I', 'Paenibacillus-C', 'Pseudomonas-A', 'Clostridium-F'])]
print(metadata)
#########3
'''


#hits = pd.read_csv(base + run + "/tables_to_export/combined_results_HvN_collapse_G_FILTEREDforPCoA.csv", usecols = ['pfam_id','clade_id', 'tests_passed', 'enrichment', 'evolink_P_A', 'evolink_counts', 'lmm', 'scoary'])
hits = pd.read_csv(base + run + "/tables_to_export/combined_results_"+LABEL+"_collapse_"+COLLAPSE+".csv", usecols = ['pfam_id','clade_id', 'tests_passed', 'enrichment', 'evolink_P_A', 'evolink_counts', 'lmm', 'scoary'])
hits = hits.loc[hits['enrichment'] == ENRICHMENT]

#list only informative clades
list_of_clades = hits['clade_id'].unique()


#pfs_of_interest = hits.loc[hits['tests_passed'] >= TESTS_PASSED]['pfam_id'].unique()

##########
def get_correlation(CLADE):
#for CLADE in list_of_clades:

    pfs_of_interest = hits.loc[(hits['tests_passed'] >= TESTS_PASSED) & (hits['clade_id'] == CLADE)]['pfam_id'].unique()
    print (pfs_of_interest)
    metadata2 = metadata.loc[metadata['clade_id'] == CLADE]

    gff_list = list()

    for genome_id in metadata2['genome_id'].unique():
        print(genome_id)
        gff3_columns = ["scaffold_ID","Source","Feature_Type","Start","End","Score","Strand","Phase","Attributes"]
        gff_df = pd.read_csv("/sci/labs/asafle/alexlevylab/icore-data/alex/workspace/HA_vs_NHA/data_sources/final_genomes/genome_data/gff/"+genome_id+".gff", 
        comment='#', 
        names=gff3_columns, 
        sep="\t",
        usecols = ['scaffold_ID', 'Start', 'End', 'Strand', "Attributes"])

        gff_df['gene_ID'] = gff_df['Attributes'].str.extract(r'ID=(.+?);')

        gff_df['genome_ID'] = genome_id

        # some gff and faa have slight mismatches in names, but it is fixable here
        #IMG sources all have matches already so no need to change
        if metadata2.loc[metadata2['genome_id'] == genome_id]['source'].values[0] in ["patric", "Garrido-Oter", "Garrido-Oter-new", "IMG_listed", "blackwell", "karasov", "tomas-white", "poyet","zou"]:
            gff_df['gene_ID'] = gff_df.apply(fix_patric, axis = 1)
        #print(gff_df) 
        
        gff_df['window_index'] = (gff_df['Start'] // window_size) + 1 
        gff_df['scaffold_ID'] = gff_df['scaffold_ID'].astype(str) + "_" + gff_df['window_index'].astype(str)
        print (gff_df)

        #read pfam file and merge
        pfam_df = pd.read_csv("/sci/labs/asafle/alexlevylab/icore-data/alex/workspace/HA_vs_NHA/data_sources/final_genomes/genome_data/new_pfams/"+genome_id+".pfam", sep="\t")
        #change column type to string for merging later
        pfam_df["Query_ID"] = pfam_df["Query_ID"].astype(str)
        #filter for evalue and coverage
        pfam_df = pfam_df[(pfam_df["Coverage"] >= 0.7) & (pfam_df["E-value"] <= 1e-5)]
        #remove "." from the end of pfam_id
        pfam_df["Pfam_id"] = pfam_df["Pfam_id"].str.split(".").str[0]
        #take olny relevant pfams - save gene_ids
        filtered_pfam_df = pfam_df[pfam_df["Pfam_id"].isin(pfs_of_interest)]
        #reduce pfam_df
        filtered_pfam_df = filtered_pfam_df[["Pfam_id","Query_ID","Description"]]
        print (filtered_pfam_df)
        #merge
        #gff_list = gff_df.merge(filtered_pfam_df, how = 'left', left_on = ['gene_ID'], right_on = ['Query_ID'])
        merged_gff_pfam = gff_df.merge(filtered_pfam_df, how = 'inner', left_on = ['gene_ID'], right_on = ['Query_ID'])
        print (merged_gff_pfam)

        merged_gff_pfam.drop(columns = ['Attributes','Query_ID'], inplace = True)
        gff_list.append(merged_gff_pfam)

        #gff_df.drop(columns = ['Attributes', 'Start'], inplace = True)
        #gff_list.append(gff_df)


    #filtered_gff_df = gff_df[gff_df["ID_value"].isin(gene_ids)]

    gff_df_all = pd.concat(gff_list)

    print(gff_df_all)

    #mer = gff_df_all.merge(mapper, how = 'left', left_on = ['gene_ID', 'genome_ID'], right_on = ['lc_member', 'genome'])

    #print(mer)

    #mer.to_csv("gff_and_repID_data_"+CLADE+".csv", index = False)

    gff_df_all.to_csv("all_clades_gff_pfams_tables/all_gff_with_pfams_"+CLADE+".csv", index=False)

################

    #mer_filtered = mer.loc[mer['repID'].isin(pfs_of_interest)]
    #print(mer_filtered)
###############
    #piv

    piv = pd.pivot_table(data = gff_df_all[['Pfam_id', 'genome_ID', 'scaffold_ID']], index = 'Pfam_id', columns = ['genome_ID', 'scaffold_ID'], aggfunc=len, fill_value=0)

    print(piv)

    scaf_correlation_matrix = piv.T.corr()

    print(scaf_correlation_matrix)

    scaf_correlation_matrix.to_csv("scaffold_correlation_matrixs/scaf_correlation_matrix_clade_" + CLADE + ".csv")

    """ 
    masked_scaf_correlation_matrix = scaf_correlation_matrix
    masked_scaf_correlation_matrix[masked_scaf_correlation_matrix < 0 ] = 0
    print(masked_scaf_correlation_matrix) 
    combined_data_final = raw_corr * masked_scaf_correlation_matrix

    print(combined_data_final)

    combined_data_final.to_csv("corr_matrix_adjusted_genus_" + GENUS + ".csv") 
    """

#make it concurrent
with concurrent.futures.ProcessPoolExecutor(max_workers=50) as executor:

        future = list(executor.map(get_correlation, list_of_clades, chunksize=1))