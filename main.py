import pandas as pd
from dataframe_operations import *
from plot_operations import *
import numpy as np

def main():
    df_geno = pd.read_excel("data/Sequencing_data_Inholland.xlsx")
    clean_geno = filter_failed_homozygote(df_geno)
    print("loaded genotype info")

    df_pheno = pd.read_excel('data/Phenotype_data_PennxMB_population.xlsx')
    df_pheno = clean_phenotype_data(df_pheno)
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

    create_manhattan_plot(preped, x_axis_name="location_offset", y_axis_name="-logp")


if __name__ == "__main__":
    main()
