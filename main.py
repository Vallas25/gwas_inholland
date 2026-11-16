import pandas as pd
from dataframe_operations import *

def main():
    df = pd.read_excel("data/Sequencing_data_Inholland.xlsx")
    df = filter_failed_hetrozygote(df)
    df.to_csv("data/get_average_test.csv")


if __name__ == "__main__":
    main()
