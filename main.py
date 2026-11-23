import pandas as pd
from dataframe_operations import *
from plot_operations import *
from helper_functions import get_sub_loci, get_top_performers
import numpy as np

def main():
    df_geno = pd.read_excel("data/Sequencing_data_Inholland.xlsx")
    print(f"total length: {len(df_geno)}")
    clean_geno = filter_failed_homozygote(df_geno)
    print(f"length after filtering: {len(clean_geno)}")
    print(f"unique markers: {len(clean_geno["Markers"].unique())}")
    print("loaded genotype info")

    df_pheno = pd.read_excel('data/Phenotype_data_PennxMB_population.xlsx')
    df_pheno = clean_phenotype_data(df_pheno)
    print("loaded phenotype info")

    merged, info = merge_geno_and_phenotype(clean_geno, df_pheno)
    print("merged geno and phenotype")

    result = group_by_location(merged)
    print("grouped by location")

    p_values = get_p_values(result)
    print("got p-values")
    print(f"length after p-values {len(p_values)}")

    final_merge = merge_pvalues_and_genotype_info(p_values,info)
    print(f"length after final merge: {len(final_merge)}")
    print("final merge compleated")

    preped = prep_df_for_plotting(final_merge)
    print("preped data for plotting")
    print(f"length after prep: {len(preped)}")

    create_manhattan_plot(preped, "location_offset", "-logp")



if __name__ == "__main__":
    main()
