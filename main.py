import pandas as pd
from dataframe_operations import *

def main():
    df = pd.read_excel("data/Phenotype_data_PennxMB_population.xlsx")
    df = df.drop(["Unnamed: 10","Unnamed: 11","Unnamed: 12","cellcount/seconds/flowers*120(=seconds to measure complete sample)"],
                 axis=1)
    df = change_series_type(df, "#Pollen per flower")
    df = is_above_mean(df, "#Pollen per flower")
    df.to_csv("data/Phenotype_data.csv")


if __name__ == "__main__":
    main()
