import pandas as pd
from dataframe_operations import *
from plot_operations import *
from helper_functions import get_sub_loci, get_top_performers
import numpy as np

def main():
    df_geno = pd.read_excel("data/Sequencing_data_Inholland.xlsx")
    clean_geno = filter_failed_homozygote(df_geno)
    print("loaded genotype info")

    clean_geno_copy = clean_geno.copy()

    df_pheno = pd.read_excel('data/Phenotype_data_PennxMB_population.xlsx')
    df_pheno = clean_phenotype_data(df_pheno)
    pheno_copy = df_pheno.copy()
    print("loaded phenotype info")

    merged, info = merge_geno_and_fenotype(clean_geno, df_pheno)
    print("merged geno and phenotype")

    result = group_by_location(merged)
    print("grouped by location")

    p_values = get_p_values(result)
    print("got p-values")

    final_merge = merge_pvalues_and_genotype_info(p_values,info)
    print("final merge compleated")

    preped = prep_df_for_plotting(final_merge)
    print("preped data for plotting")

    meta_data = ["Markers", "Chromosome", "location"]

    top = get_top_performers(pheno_copy, asc=True).index.to_list()
    print(top)

    parental = [
        "F1.1","F1.2","Parent1.1","Parent1.2","Parent2.1","Parent2.2",
    ]

    cols = meta_data + top + parental

    sub_set = clean_geno_copy[cols]

    significant = get_sub_loci(preped, threshold=6.5)

    sign_list = significant["Markers"].to_list()
    print(sign_list)

    filtered = sub_set[
        (clean_geno_copy["Chromosome"].str.contains("ch04")) &
        (clean_geno_copy["Markers"].isin(sign_list)) 
        ]

    filtered.to_csv("data/significant_markers_bottem.csv")



if __name__ == "__main__":
    main()
