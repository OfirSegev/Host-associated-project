import pandas as pd
import sys

COLLAPSE=sys.argv[1]
LABEL=sys.argv[2]
MIN_CLADES=int(sys.argv[3])

# Define a function to calculate the weighted sum of a row
def weighted_sum(row):
    return sum(cell ** 2 for cell in row)

#filter pfam results that are only fisher
def filter_for_only_fisher(row,df):

    pfam_df = df[df["pfam_id"] == row["pfam_id"]]

    if(pfam_df["fisher"].sum() == pfam_df["tests_passed"].sum()):
        return False
    else:
        return True

def pivot(enrichment, min_clades):

    df=pd.read_csv("final_results_FL/final_table_all_tests_"+enrichment+"_results_FL_"+LABEL+"_collapse_"+COLLAPSE+".csv")

    #add number of clades column
    numclades = df.groupby("pfam_id").apply(lambda x: len(x['clade_id']))
    numclades = numclades.reset_index()
    df = df.merge(numclades, on = 'pfam_id', how = 'left').rename(columns={0:'number_of_clades'})

    print (df)

    #remove rows that are only fisher
    filtered_for_fisher_df = df[df.apply(filter_for_only_fisher, axis=1, args=(df,))]
    #filter for number of clades
    filtered_df = filtered_for_fisher_df[filtered_for_fisher_df["number_of_clades"] >= min_clades]

    p = pd.pivot_table(filtered_df, index = ['pfam_id','description'], columns = 'clade_id', values = 'tests_passed', fill_value = 0)
    #p['sum'] = p.sum(axis = 1)
    p['sum'] = p.apply(weighted_sum, axis=1)

    # Calculate the sum of each column
    sum_row = p.sum()
    # Convert the sum to a DataFrame with a single row
    sum_df = pd.DataFrame([sum_row], columns=p.columns, index=['sum_rows'])
    # Append the sum row to the original DataFrame
    p = pd.concat([p, sum_df])

    p = p.sort_values(by='sum_rows', axis=1, ascending = False)
    p = p.sort_values('sum', ascending = False)

    print (p)

    p.to_csv("final_results_FL/FL_pivot_table_"+enrichment+"_results_"+COLLAPSE+"_"+LABEL+"_at_least_"+str(min_clades)+"_clades.csv")


pivot("enriched",MIN_CLADES)
pivot("depleted",MIN_CLADES)
pivot("enriched",0)
pivot("depleted",0)