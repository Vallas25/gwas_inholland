import pandas as pd
from dataframe_operations import *
from plot_operations import *
from helper_functions import get_sub_loci, get_top_performers
import numpy as np

def main():
    df_geno = pd.read_excel("data/Sequencing_data_Inholland.xlsx")
    print(f"total length: {len(df_geno)}")
    clean_geno = filter_failed_homozygote(df_geno)
    clean_geno_copy = clean_geno.copy()
    print(f"length after filtering: {len(clean_geno)}")
    print(f"unique markers: {len(clean_geno["Markers"].unique())}")
    print("loaded genotype info")

    df_pheno = pd.read_excel('data/Phenotype_data_PennxMB_population.xlsx')
    df_pheno = clean_phenotype_data(df_pheno)
    pheno_copy = df_pheno.copy()
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

    meta_data = ["Markers", "Chromosome", "location"]

    top = get_top_performers(pheno_copy, asc=False).index.to_list()
    print(top)

    parental = [
        "F1.1","F1.2","Parent1.1","Parent1.2","Parent2.1","Parent2.2",
    ]

    cols = meta_data + top + parental

    sub_set = clean_geno_copy[cols]

    significant = get_sub_loci(preped,chr=2, threshold=3)

    sign_list = significant["Markers"].to_list()
    print(sign_list)

    filtered = sub_set[
        (sub_set["Chromosome"] == "SL4.0ch02") &
        (sub_set["Markers"].isin(sign_list)) 
        ]

    print(preped[(preped["chr_num"] == 2) & (preped["-logp"] > 4)])

    filtered.to_csv("data/ch02_top.csv")

if __name__ == "__main__":
    main()
