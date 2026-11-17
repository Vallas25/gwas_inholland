import pandas as pd
from dataframe_operations import *
from plot_operations import *
import numpy as np

def main():
    df = pd.read_excel("data/Phenotype_data_PennxMB_population.xlsx")
    df = clean_phenotype_data(df)
    grouped = t_test_M(df)
    grouped["p_value"] = -np.log10(grouped["p_value"])
    create_manhattan_like_plot(grouped, "p_value")


if __name__ == "__main__":
    main()
