import pandas as pd
from dataframe_operations import *
from plot_operations import *
import numpy as np

def main():
    df_geno = pd.read_excel("data/Sequencing_data_Inholland.xlsx")
    print(df_geno.head())
    print("loaded genotype info")

    clean_geno = filter_failed_homozygote(df_geno)
    print(clean_geno.head())

    df_pheno = pd.read_excel('data/Phenotype_data_PennxMB_population.xlsx')
    df_pheno = clean_phenotype_data(df_pheno)

    print(df_pheno.head())
    print("loaded phenotype info")

    merged, info = merge_geno_and_fenotype(clean_geno, df_pheno)
    print(merged.head())
    print("merged geno and phenotype")

    result = group_by_location(merged)
    print(result[:6])
    print("grouped by location")

    p_values = get_p_values(result)
    print(p_values.head())
    print("got p-values")

    p_values["-logp"] = -np.log10(p_values["p_value"])
    final_merge = merge_pvalues_and_genotype_info(p_values,info)
    print("final merge compleated")

    create_manhattan_plot(final_merge, y_axis_name="-logp", x_axis_name="location")

if __name__ == "__main__":
    main()
